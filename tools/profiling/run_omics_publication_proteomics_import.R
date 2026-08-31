#!/usr/bin/env Rscript

proteomicsImportRunnerPath <- function() {
    file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
    if (length(file_arg)) {
        return(normalizePath(sub("^--file=", "", file_arg[[1L]]), mustWork = TRUE))
    }
    normalizePath(
        file.path(
            getwd(),
            "tools",
            "profiling",
            "run_omics_publication_proteomics_import.R"
        ),
        mustWork = TRUE
    )
}

.PROTEOMICS_IMPORT_REPO_ROOT <- normalizePath(
    file.path(dirname(proteomicsImportRunnerPath()), "..", ".."),
    mustWork = TRUE
)

for (source_file in c(
    "omics_publication_protocol.R",
    "omics_publication_comparators.R",
    "omics_publication_lock.R",
    "omics_publication_builds.R",
    "omics_publication_restore_reproducibility.R",
    "omics_publication_comparator_builds.R",
    "omics_publication_comparator_build_reproducibility.R",
    "omics_publication_comparator_evidence.R",
    "omics_publication_proteomics_contracts.R",
    "omics_publication_proteomics_truth.R",
    "omics_publication_workload_proteomics.R"
)) {
    source(
        file.path(
            .PROTEOMICS_IMPORT_REPO_ROOT,
            "tools",
            "profiling",
            source_file
        ),
        local = FALSE
    )
}

proteomicsImportRunnerArgs <- function(argv) {
    values <- list(
        contract = NULL,
        payload = NULL,
        truth = NULL,
        build_receipt = NULL,
        output = NULL
    )
    index <- 1L
    while (index <= length(argv)) {
        token <- argv[[index]]
        key <- gsub("-", "_", sub("^--", "", token), fixed = TRUE)
        if (!startsWith(token, "--") || !key %in% names(values) ||
            index == length(argv)) {
            stop("Comparator import arguments are invalid", call. = FALSE)
        }
        values[[key]] <- argv[[index + 1L]]
        index <- index + 2L
    }
    if (any(vapply(values, is.null, logical(1)))) {
        stop("Comparator import options are incomplete", call. = FALSE)
    }
    values
}

proteomicsImportRunnerLibraries <- function(path) {
    evidence <- publicationValidateComparatorBuildReceiptEvidence(path)
    receipt <- evidence$receipt
    list(
        receipt = receipt,
        package_library = dirname(receipt$build$installed_inventory$root),
        dependency_library = file.path(
            dirname(receipt$environment$restore_receipt$path),
            "library"
        )
    )
}

proteomicsImportRunnerMain <- function(argv = commandArgs(trailingOnly = TRUE)) {
    args <- proteomicsImportRunnerArgs(argv)
    if (file.exists(args$output)) {
        stop("Comparator import output already exists", call. = FALSE)
    }
    contract <- proteomicsPublicationReadContract(args$contract)
    if (!identical(
        publicationFileDigest(args$payload),
        contract$expected_digests$payload_sha256
    ) || !identical(
        publicationFileDigest(args$truth),
        contract$expected_digests$truth_sha256
    )) {
        proteomicsPublicationAbort("comparator import input differs")
    }
    libraries <- proteomicsImportRunnerLibraries(args$build_receipt)
    .libPaths(c(
        libraries$package_library,
        libraries$dependency_library,
        .Library
    ))
    package_path <- find.package(
        "MultiScholaR",
        lib.loc = libraries$package_library
    )
    result <- proteomicsPublicationRunImported(
        contract,
        args$payload,
        args$truth
    )
    evidence <- list(
        schema = "multischolar.omics_publication_proteomics_import_evidence",
        schema_version = "1.0.0",
        owner_ticket_id = "OMICS-ART-065",
        status = "passed",
        comparator_id = libraries$receipt$comparator_id,
        comparator_revision = libraries$receipt$revision,
        build_receipt_sha256 = publicationFileDigest(args$build_receipt),
        package_path_sha256 = digest::digest(
            normalizePath(package_path),
            algo = "sha256",
            serialize = FALSE
        ),
        contract_sha256 = publicationFileDigest(args$contract),
        payload_sha256 = publicationFileDigest(args$payload),
        truth_sha256 = publicationFileDigest(args$truth),
        workflow_evidence = result$workflow_evidence,
        candidate_loaded = FALSE,
        promotion_authority = FALSE,
        publication_authority = FALSE
    )
    publicationWriteJson(evidence, args$output)
    invisible(0L)
}

if (identical(environment(), globalenv()) && !interactive()) {
    status <- tryCatch(
        proteomicsImportRunnerMain(),
        error = \(error) {
            message(conditionMessage(error))
            1L
        }
    )
    quit(status = as.integer(status), save = "no")
}
