#!/usr/bin/env Rscript

proteomicsPreparationRunnerPath <- function() {
    file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
    if (length(file_arg)) {
        return(normalizePath(sub("^--file=", "", file_arg[[1L]]), mustWork = TRUE))
    }
    normalizePath(
        file.path(
            getwd(),
            "tools",
            "profiling",
            "run_omics_publication_workload_prepare.R"
        ),
        mustWork = TRUE
    )
}

.PROTEOMICS_PREPARATION_REPO_ROOT <- normalizePath(
    file.path(dirname(proteomicsPreparationRunnerPath()), "..", ".."),
    mustWork = TRUE
)

for (source_file in c(
    "omics_publication_protocol.R",
    "omics_publication_proteomics_contracts.R",
    "omics_publication_proteomics_model.R",
    "omics_publication_proteomics_serializers.R",
    "omics_publication_proteomics_truth.R",
    "omics_publication_proteomics_sources.R",
    "omics_publication_proteomics_governance.R",
    "omics_publication_workload_proteomics.R"
)) {
    source(
        file.path(
            .PROTEOMICS_PREPARATION_REPO_ROOT,
            "tools",
            "profiling",
            source_file
        ),
        local = FALSE
    )
}

proteomicsPreparationRunnerArgs <- function(argv) {
    values <- list(contract = NULL, output_root = NULL, summary = NULL)
    index <- 1L
    while (index <= length(argv)) {
        token <- argv[[index]]
        key <- gsub("-", "_", sub("^--", "", token), fixed = TRUE)
        if (!startsWith(token, "--") || !key %in% names(values) ||
            index == length(argv)) {
            stop(
                paste(
                    "Usage: --contract <json> --output-root <dir>",
                    "--summary <json>"
                ),
                call. = FALSE
            )
        }
        values[[key]] <- argv[[index + 1L]]
        index <- index + 2L
    }
    if (any(vapply(values, is.null, logical(1)))) {
        stop("Preparation runner options are incomplete", call. = FALSE)
    }
    values
}

proteomicsPreparationRunnerMain <- function(
    argv = commandArgs(trailingOnly = TRUE)
) {
    args <- proteomicsPreparationRunnerArgs(argv)
    if (file.exists(args$summary)) {
        stop("Preparation summary already exists", call. = FALSE)
    }
    contract <- proteomicsPublicationReadContract(args$contract)
    prepared <- if (identical(
        contract$workload_class,
        "fixture_correctness"
    )) {
        proteomicsPublicationPrepareFixture(contract, args$output_root)
    } else {
        proteomicsPublicationPrepareGenerated(
            contract,
            args$output_root,
            verify_expected = TRUE
        )
    }
    summary <- list(
        schema = "multischolar.omics_publication_preparation_summary",
        schema_version = "1.0.0",
        owner_ticket_id = "OMICS-ART-065",
        status = "passed",
        workload_id = contract$workload_id,
        contract_sha256 = publicationFileDigest(args$contract),
        payload_sha256 = prepared$payload$sha256,
        truth_sha256 = prepared$truth$sha256,
        preparation_sha256 = prepared$receipt$sha256,
        candidate_loaded = FALSE,
        publication_authority = FALSE
    )
    publicationWriteJson(summary, args$summary)
    invisible(0L)
}

if (identical(environment(), globalenv()) && !interactive()) {
    status <- tryCatch(
        proteomicsPreparationRunnerMain(),
        error = \(error) {
            message(conditionMessage(error))
            1L
        }
    )
    quit(status = as.integer(status), save = "no")
}
