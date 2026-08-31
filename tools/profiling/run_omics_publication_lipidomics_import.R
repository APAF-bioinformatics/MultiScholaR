#!/usr/bin/env Rscript

lipidImportRunnerPath <- function() {
    file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
    if (length(file_arg)) {
        return(normalizePath(sub("^--file=", "", file_arg[[1L]]), mustWork = TRUE))
    }
    normalizePath(
        file.path(
            getwd(),
            "tools",
            "profiling",
            "run_omics_publication_lipidomics_import.R"
        ),
        mustWork = TRUE
    )
}

.LIPID_IMPORT_REPO_ROOT <- normalizePath(
    file.path(dirname(lipidImportRunnerPath()), "..", ".."),
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
    "omics_publication_lipidomics_contracts.R",
    "omics_publication_lipidomics_serializers.R",
    "omics_publication_lipidomics_truth.R",
    "omics_publication_lipidomics_sources.R",
    "omics_publication_lipidomics_governance.R",
    "omics_publication_workload_lipidomics.R"
)) {
    source(
        file.path(.LIPID_IMPORT_REPO_ROOT, "tools", "profiling", source_file),
        local = FALSE
    )
}

lipidImportRunnerArgs <- function(argv) {
    values <- list(
        contract = NULL,
        payload_root = NULL,
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
            stop("Lipidomics import arguments are invalid", call. = FALSE)
        }
        values[[key]] <- argv[[index + 1L]]
        index <- index + 2L
    }
    if (any(vapply(values, is.null, logical(1)))) {
        stop("Lipidomics import options are incomplete", call. = FALSE)
    }
    values
}

lipidImportRunnerLibraries <- function(path) {
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

lipidImportRunnerMain <- function(argv = commandArgs(trailingOnly = TRUE)) {
    args <- lipidImportRunnerArgs(argv)
    if (file.exists(args$output)) {
        stop("Lipidomics import output already exists", call. = FALSE)
    }
    contract <- lipidPublicationReadContract(args$contract)
    members <- unlist(contract$generator$output_members, use.names = FALSE)
    payload <- lipidPublicationPayloadBinding(file.path(args$payload_root, members))
    if (!identical(
        payload$payload_set_sha256,
        contract$expected_digests$payload_set_sha256
    ) || !identical(
        publicationFileDigest(args$truth),
        contract$expected_digests$truth_sha256
    )) {
        lipidPublicationAbort("comparator lipidomics input differs")
    }
    libraries <- lipidImportRunnerLibraries(args$build_receipt)
    .libPaths(c(
        libraries$package_library,
        libraries$dependency_library,
        .Library
    ))
    package_path <- find.package(
        "MultiScholaR",
        lib.loc = libraries$package_library
    )
    result <- lipidPublicationRunImported(
        contract,
        args$payload_root,
        args$truth
    )
    evidence <- list(
        schema = "multischolar.omics_publication_lipidomics_import_evidence",
        schema_version = "1.0.0",
        owner_ticket_id = "OMICS-ART-067",
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
        payload_set_sha256 = payload$payload_set_sha256,
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
        lipidImportRunnerMain(),
        error = function(error) {
            message(conditionMessage(error))
            1L
        }
    )
    quit(status = as.integer(status), save = "no")
}
