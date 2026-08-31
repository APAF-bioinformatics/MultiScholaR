#!/usr/bin/env Rscript

auxExecuteRunnerPath <- function() {
    file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
    if (length(file_arg)) {
        return(normalizePath(sub("^--file=", "", file_arg[[1L]]), mustWork = TRUE))
    }
    normalizePath(
        file.path(
            getwd(),
            "tools",
            "profiling",
            "run_omics_publication_auxiliary_execute.R"
        ),
        mustWork = TRUE
    )
}

.AUX_EXECUTE_REPO_ROOT <- normalizePath(
    file.path(dirname(auxExecuteRunnerPath()), "..", ".."),
    mustWork = TRUE
)

.AUX_EXECUTE_LIBRARY_CACHE <- new.env(parent = emptyenv())

for (source_file in c(
    "omics_publication_protocol.R",
    "omics_publication_comparators.R",
    "omics_publication_lock.R",
    "omics_publication_builds.R",
    "omics_publication_restore_reproducibility.R",
    "omics_publication_comparator_builds.R",
    "omics_publication_comparator_build_reproducibility.R",
    "omics_publication_comparator_evidence.R",
    "omics_publication_auxiliary_contracts.R",
    "omics_publication_auxiliary_model.R",
    "omics_publication_auxiliary_responses.R",
    "omics_publication_auxiliary_sources.R",
    "omics_publication_auxiliary_governance.R",
    "omics_publication_auxiliary_truth.R",
    "omics_publication_workload_auxiliary.R"
)) {
    source(
        file.path(.AUX_EXECUTE_REPO_ROOT, "tools", "profiling", source_file),
        local = FALSE
    )
}

auxExecuteRunnerArgs <- function(argv) {
    values <- list(
        contract = NULL,
        payload = NULL,
        truth = NULL,
        build_receipt = NULL,
        output_root = NULL,
        output = NULL
    )
    index <- 1L
    while (index <= length(argv)) {
        token <- argv[[index]]
        key <- gsub("-", "_", sub("^--", "", token), fixed = TRUE)
        if (!startsWith(token, "--") || !key %in% names(values) ||
            index == length(argv)) {
            stop("Auxiliary execute arguments are invalid", call. = FALSE)
        }
        values[[key]] <- argv[[index + 1L]]
        index <- index + 2L
    }
    if (any(vapply(values, is.null, logical(1)))) {
        stop("Auxiliary execute options are incomplete", call. = FALSE)
    }
    values
}

auxExecuteRunnerLibraries <- function(path) {
    normalized <- normalizePath(path, mustWork = TRUE)
    digest <- publicationFileDigest(normalized)
    key <- digest::digest(normalized, algo = "sha256", serialize = FALSE)
    if (exists(key, envir = .AUX_EXECUTE_LIBRARY_CACHE, inherits = FALSE)) {
        cached <- get(key, envir = .AUX_EXECUTE_LIBRARY_CACHE, inherits = FALSE)
        if (identical(cached$receipt_sha256, digest)) return(cached$libraries)
    }
    evidence <- publicationValidateComparatorBuildReceiptEvidence(path)
    receipt <- evidence$receipt
    libraries <- list(
        receipt = receipt,
        package_library = dirname(receipt$build$installed_inventory$root),
        dependency_library = file.path(
            dirname(receipt$environment$restore_receipt$path),
            "library"
        )
    )
    assign(
        key,
        list(receipt_sha256 = digest, libraries = libraries),
        envir = .AUX_EXECUTE_LIBRARY_CACHE
    )
    libraries
}

auxExecuteWarningCode <- function(message) {
    if (grepl("cross2", message, fixed = TRUE)) {
        return("purrr_cross2_deprecated")
    }
    if (grepl("already exists", message, fixed = TRUE)) {
        return("output_directory_exists")
    }
    if (grepl("deprecated|deprecation", message, ignore.case = TRUE)) {
        return("dependency_deprecation")
    }
    "unclassified_warning"
}

auxExecuteRunnerMain <- function(argv = commandArgs(trailingOnly = TRUE)) {
    args <- auxExecuteRunnerArgs(argv)
    if (file.exists(args$output) || file.exists(args$output_root) ||
        dir.exists(args$output_root)) {
        stop("Auxiliary execute output already exists", call. = FALSE)
    }
    contract <- auxPublicationReadContract(args$contract)
    payload <- auxPublicationLoadPayload(args$payload, contract)
    truth <- publicationReadJson(args$truth)
    auxPublicationValidateTruth(truth, contract, args$payload)
    if (!identical(
        publicationFileDigest(args$truth),
        contract$expected_digests$truth_sha256
    )) {
        auxPublicationAbort("auxiliary comparator truth binding differs")
    }
    libraries <- auxExecuteRunnerLibraries(args$build_receipt)
    .libPaths(c(
        libraries$package_library,
        libraries$dependency_library,
        .Library
    ))
    package_path <- find.package(
        "MultiScholaR",
        lib.loc = libraries$package_library
    )
    loadNamespace("MultiScholaR", lib.loc = libraries$package_library)
    warnings <- character()
    output_text <- capture.output(
        result <- withCallingHandlers(
            suppressMessages(auxPublicationRunPayload(
                contract,
                payload,
                args$output_root
            )),
            warning = function(warning) {
                warnings <<- c(warnings, conditionMessage(warning))
                invokeRestart("muffleWarning")
            }
        ),
        type = "output"
    )
    validation_error <- tryCatch({
        auxPublicationValidateRunResult(result, truth, contract)
        NULL
    }, error = function(error) error)
    passed <- is.null(validation_error)
    warning_codes <- unname(vapply(
        warnings,
        auxExecuteWarningCode,
        character(1)
    ))
    warning_code_counts <- as.list(table(factor(
        warning_codes,
        levels = sort(unique(warning_codes), method = "radix")
    )))
    evidence <- list(
        schema = "multischolar.omics_publication_auxiliary_execution_evidence",
        schema_version = .AUX_PUBLICATION_VERSION,
        owner_ticket_id = .AUX_PUBLICATION_OWNER,
        status = if (passed) "passed" else "compatibility_mismatch",
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
        workload_id = contract$workload_id,
        route_id = contract$route$route_id,
        workload_class = contract$workload_class,
        workflow_evidence = auxPublicationResultEvidence(result),
        validation = if (passed) {
            list(status = "passed")
        } else {
            list(
                status = "failed",
                condition_class = class(validation_error)[[1L]],
                message = conditionMessage(validation_error),
                expected_facts_sha256 = publicationObjectDigest(truth$facts)
            )
        },
        warning_count = as.integer(length(warnings)),
        warning_codes = as.list(warning_codes),
        warning_code_counts = warning_code_counts,
        warning_messages_sha256 = publicationObjectDigest(as.list(warnings)),
        console_line_count = as.integer(length(output_text)),
        network_allowed = FALSE,
        candidate_loaded = FALSE,
        promotion_authority = FALSE,
        publication_authority = FALSE
    )
    publicationWriteJson(evidence, args$output)
    cat(publicationFileDigest(args$output), "\n")
    invisible(if (passed) 0L else 2L)
}

if (identical(environment(), globalenv()) && !interactive()) {
    status <- tryCatch(
        auxExecuteRunnerMain(),
        error = function(error) {
            message(conditionMessage(error))
            1L
        }
    )
    quit(status = as.integer(status), save = "no")
}
