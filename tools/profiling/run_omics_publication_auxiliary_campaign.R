#!/usr/bin/env Rscript

auxCampaignRunnerPath <- function() {
    file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
    if (length(file_arg)) {
        return(normalizePath(sub("^--file=", "", file_arg[[1L]]), mustWork = TRUE))
    }
    normalizePath(
        file.path(
            getwd(),
            "tools",
            "profiling",
            "run_omics_publication_auxiliary_campaign.R"
        ),
        mustWork = TRUE
    )
}

.AUX_CAMPAIGN_REPO_ROOT <- normalizePath(
    file.path(dirname(auxCampaignRunnerPath()), "..", ".."),
    mustWork = TRUE
)

.aux_campaign_environment <- new.env(parent = globalenv())
sys.source(
    file.path(
        .AUX_CAMPAIGN_REPO_ROOT,
        "tools",
        "profiling",
        "run_omics_publication_auxiliary_execute.R"
    ),
    envir = .aux_campaign_environment
)

auxCampaignRunnerArgs <- function(argv) {
    values <- list(
        build_receipt = NULL,
        freeze_root = NULL,
        contract_root = NULL,
        truth_root = NULL,
        output_root = NULL,
        routes = NULL,
        classes = NULL,
        summary = NULL
    )
    index <- 1L
    while (index <= length(argv)) {
        token <- argv[[index]]
        key <- gsub("-", "_", sub("^--", "", token), fixed = TRUE)
        if (!startsWith(token, "--") || !key %in% names(values) ||
            index == length(argv)) {
            stop("Auxiliary campaign arguments are invalid", call. = FALSE)
        }
        values[[key]] <- argv[[index + 1L]]
        index <- index + 2L
    }
    if (any(vapply(values, is.null, logical(1)))) {
        stop("Auxiliary campaign options are incomplete", call. = FALSE)
    }
    values$routes <- strsplit(values$routes, ",", fixed = TRUE)[[1L]]
    values$classes <- strsplit(values$classes, ",", fixed = TRUE)[[1L]]
    values
}

auxCampaignWorkloadId <- function(route_id, workload_class) {
    paste("auxiliary", route_id, workload_class, "v1", sep = ".")
}

auxCampaignRunnerMain <- function(argv = commandArgs(trailingOnly = TRUE)) {
    args <- auxCampaignRunnerArgs(argv)
    if (file.exists(args$summary)) {
        stop("Auxiliary campaign summary already exists", call. = FALSE)
    }
    libraries <- .aux_campaign_environment$auxExecuteRunnerLibraries(
        args$build_receipt
    )
    .libPaths(c(
        libraries$package_library,
        libraries$dependency_library,
        .Library
    ))
    records <- list()
    for (workload_class in args$classes) {
        for (route_id in args$routes) {
            workload_id <- auxCampaignWorkloadId(route_id, workload_class)
            output_dir <- file.path(args$output_root, workload_id)
            evidence_path <- paste0(output_dir, ".json")
            if (file.exists(evidence_path)) {
                evidence <- publicationReadJson(evidence_path)
                status <- if (identical(evidence$status, "passed")) 0L else 2L
            } else {
                if (file.exists(output_dir) || dir.exists(output_dir)) {
                    stop("Auxiliary campaign contains a partial route", call. = FALSE)
                }
                status <- .aux_campaign_environment$auxExecuteRunnerMain(c(
                    "--contract",
                    file.path(args$contract_root, paste0(workload_id, ".json")),
                    "--payload",
                    file.path(args$freeze_root, workload_id, "payload.rds"),
                    "--truth",
                    file.path(args$truth_root, paste0(workload_id, ".json")),
                    "--build-receipt",
                    args$build_receipt,
                    "--output-root",
                    output_dir,
                    "--output",
                    evidence_path
                ))
                evidence <- publicationReadJson(evidence_path)
            }
            records[[length(records) + 1L]] <- list(
                workload_id = workload_id,
                route_id = route_id,
                workload_class = workload_class,
                status = evidence$status,
                process_status = as.integer(status),
                evidence_path = evidence_path,
                evidence_sha256 = publicationFileDigest(evidence_path),
                promotion_authority = FALSE
            )
        }
    }
    passed <- all(vapply(records, function(record) {
        identical(record$status, "passed") &&
            identical(record$process_status, 0L)
    }, logical(1)))
    summary <- list(
        schema = "multischolar.omics_publication_auxiliary_campaign",
        schema_version = "1.0.0",
        owner_ticket_id = "OMICS-ART-068",
        status = if (passed) "passed" else "compatibility_mismatch",
        comparator_id = libraries$receipt$comparator_id,
        comparator_revision = libraries$receipt$revision,
        build_receipt_sha256 =
            publicationFileDigest(args$build_receipt),
        records = records,
        candidate_loaded = FALSE,
        promotion_authority = FALSE,
        publication_authority = FALSE
    )
    publicationWriteJson(summary, args$summary)
    cat(publicationFileDigest(args$summary), "\n")
    invisible(if (passed) 0L else 2L)
}

if (identical(environment(), globalenv()) && !interactive()) {
    status <- tryCatch(
        auxCampaignRunnerMain(),
        error = function(error) {
            message(conditionMessage(error))
            1L
        }
    )
    quit(status = as.integer(status), save = "no")
}
