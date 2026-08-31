#!/usr/bin/env Rscript

lipidPreparationRunnerPath <- function() {
    file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
    if (length(file_arg)) {
        return(normalizePath(sub("^--file=", "", file_arg[[1L]]), mustWork = TRUE))
    }
    normalizePath(
        file.path(
            getwd(),
            "tools",
            "profiling",
            "run_omics_publication_lipidomics_prepare.R"
        ),
        mustWork = TRUE
    )
}

.LIPID_PREPARATION_REPO_ROOT <- normalizePath(
    file.path(dirname(lipidPreparationRunnerPath()), "..", ".."),
    mustWork = TRUE
)

for (source_file in c(
    "omics_publication_protocol.R",
    "omics_publication_lipidomics_contracts.R",
    "omics_publication_lipidomics_model.R",
    "omics_publication_lipidomics_serializers.R",
    "omics_publication_lipidomics_truth.R",
    "omics_publication_lipidomics_sources.R",
    "omics_publication_lipidomics_governance.R",
    "omics_publication_lipidomics_negative.R",
    "omics_publication_lipidomics_fixtures.R",
    "omics_publication_lipidomics_freeze.R",
    "omics_publication_workload_lipidomics.R"
)) {
    source(
        file.path(.LIPID_PREPARATION_REPO_ROOT, "tools", "profiling", source_file),
        local = FALSE
    )
}

lipidPreparationRunnerArgs <- function(argv) {
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
            stop("Lipidomics preparation arguments are invalid", call. = FALSE)
        }
        values[[key]] <- argv[[index + 1L]]
        index <- index + 2L
    }
    if (any(vapply(values[c("contract", "output_root", "summary")],
        is.null, logical(1)))) {
        stop("Lipidomics preparation options are incomplete", call. = FALSE)
    }
    if (!values$verify_expected %in% c("true", "false")) {
        stop("verify-expected must be true or false", call. = FALSE)
    }
    values$verify_expected <- identical(values$verify_expected, "true")
    values
}

lipidPreparationRunnerMain <- function(argv = commandArgs(trailingOnly = TRUE)) {
    args <- lipidPreparationRunnerArgs(argv)
    if (file.exists(args$summary)) {
        stop("Lipidomics preparation summary already exists", call. = FALSE)
    }
    contract <- lipidPublicationReadContract(args$contract)
    prepared <- if (identical(
        contract$workload_class,
        "fixture_correctness"
    )) {
        lipidPublicationPrepareFixture(contract, args$output_root)
    } else {
        lipidPublicationPrepareGenerated(
            contract,
            args$output_root,
            verify_expected = args$verify_expected
        )
    }
    summary <- list(
        schema = "multischolar.omics_publication_lipidomics_preparation_summary",
        schema_version = "1.0.0",
        owner_ticket_id = "OMICS-ART-067",
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
        lipidPreparationRunnerMain(),
        error = function(error) {
            message(conditionMessage(error))
            1L
        }
    )
    quit(status = as.integer(status), save = "no")
}
