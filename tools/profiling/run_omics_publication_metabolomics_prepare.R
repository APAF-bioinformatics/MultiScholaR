#!/usr/bin/env Rscript

metabPreparationRunnerPath <- function() {
    file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
    if (length(file_arg)) {
        return(normalizePath(sub("^--file=", "", file_arg[[1L]]), mustWork = TRUE))
    }
    normalizePath(
        file.path(
            getwd(),
            "tools",
            "profiling",
            "run_omics_publication_metabolomics_prepare.R"
        ),
        mustWork = TRUE
    )
}

.METAB_PREPARATION_REPO_ROOT <- normalizePath(
    file.path(dirname(metabPreparationRunnerPath()), "..", ".."),
    mustWork = TRUE
)

for (source_file in c(
    "omics_publication_protocol.R",
    "omics_publication_metabolomics_contracts.R",
    "omics_publication_metabolomics_model.R",
    "omics_publication_metabolomics_serializers.R",
    "omics_publication_metabolomics_truth.R",
    "omics_publication_metabolomics_sources.R",
    "omics_publication_metabolomics_governance.R",
    "omics_publication_metabolomics_negative.R",
    "omics_publication_metabolomics_fixtures.R",
    "omics_publication_metabolomics_freeze.R",
    "omics_publication_workload_metabolomics.R"
)) {
    source(
        file.path(.METAB_PREPARATION_REPO_ROOT, "tools", "profiling", source_file),
        local = FALSE
    )
}

metabPreparationRunnerArgs <- function(argv) {
    values <- list(
        contract = NULL,
        output_root = NULL,
        summary = NULL,
        verify_expected = "true"
    )
    index <- 1L
    while (index <= length(argv)) {
        token <- argv[[index]]
        key <- gsub("-", "_", sub("^--", "", token), fixed = TRUE)
        if (!startsWith(token, "--") || !key %in% names(values) ||
            index == length(argv)) {
            stop("Metabolomics preparation arguments are invalid", call. = FALSE)
        }
        values[[key]] <- argv[[index + 1L]]
        index <- index + 2L
    }
    if (any(vapply(values[c("contract", "output_root", "summary")],
        is.null, logical(1)))) {
        stop("Metabolomics preparation options are incomplete", call. = FALSE)
    }
    if (!values$verify_expected %in% c("true", "false")) {
        stop("verify-expected must be true or false", call. = FALSE)
    }
    values$verify_expected <- identical(values$verify_expected, "true")
    values
}

metabPreparationRunnerMain <- function(argv = commandArgs(trailingOnly = TRUE)) {
    args <- metabPreparationRunnerArgs(argv)
    if (file.exists(args$summary)) {
        stop("Metabolomics preparation summary already exists", call. = FALSE)
    }
    contract <- metabPublicationReadContract(args$contract)
    prepared <- if (identical(
        contract$workload_class,
        "fixture_correctness"
    )) {
        metabPublicationPrepareFixture(contract, args$output_root)
    } else {
        metabPublicationPrepareGenerated(
            contract,
            args$output_root,
            verify_expected = args$verify_expected
        )
    }
    summary <- list(
        schema = "multischolar.omics_publication_metabolomics_preparation_summary",
        schema_version = "1.0.0",
        owner_ticket_id = "OMICS-ART-066",
        status = "passed",
        workload_id = contract$workload_id,
        contract_sha256 = publicationFileDigest(args$contract),
        payload_set_sha256 = prepared$payload$payload_set_sha256,
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
        metabPreparationRunnerMain(),
        error = function(error) {
            message(conditionMessage(error))
            1L
        }
    )
    quit(status = as.integer(status), save = "no")
}
