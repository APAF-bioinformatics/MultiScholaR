#!/usr/bin/env Rscript

auxPreparationRunnerPath <- function() {
    file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
    if (length(file_arg)) {
        return(normalizePath(sub("^--file=", "", file_arg[[1L]]), mustWork = TRUE))
    }
    normalizePath(
        file.path(
            getwd(),
            "tools",
            "profiling",
            "run_omics_publication_auxiliary_prepare.R"
        ),
        mustWork = TRUE
    )
}

.AUX_PREPARATION_REPO_ROOT <- normalizePath(
    file.path(dirname(auxPreparationRunnerPath()), "..", ".."),
    mustWork = TRUE
)

for (source_file in c(
    "omics_publication_protocol.R",
    "omics_publication_auxiliary_contracts.R",
    "omics_publication_auxiliary_model.R",
    "omics_publication_auxiliary_responses.R",
    "omics_publication_auxiliary_sources.R",
    "omics_publication_auxiliary_governance.R",
    "omics_publication_auxiliary_exclusions.R",
    "omics_publication_auxiliary_truth.R",
    "omics_publication_workload_auxiliary.R",
    "omics_publication_auxiliary_negative.R",
    "omics_publication_auxiliary_freeze.R"
)) {
    source(
        file.path(.AUX_PREPARATION_REPO_ROOT, "tools", "profiling", source_file),
        local = FALSE
    )
}

auxPreparationRunnerArgs <- function(argv) {
    values <- list(
        output_root = NULL,
        contract_root = NULL,
        truth_root = NULL,
        classes = "fixture_correctness,representative",
        summary = NULL
    )
    index <- 1L
    while (index <= length(argv)) {
        token <- argv[[index]]
        key <- gsub("-", "_", sub("^--", "", token), fixed = TRUE)
        if (!startsWith(token, "--") || !key %in% names(values) ||
            index == length(argv)) {
            stop("Auxiliary preparation arguments are invalid", call. = FALSE)
        }
        values[[key]] <- argv[[index + 1L]]
        index <- index + 2L
    }
    required <- c("output_root", "contract_root", "truth_root", "summary")
    if (any(vapply(values[required], is.null, logical(1)))) {
        stop("Auxiliary preparation options are incomplete", call. = FALSE)
    }
    values$classes <- strsplit(values$classes, ",", fixed = TRUE)[[1L]]
    values
}

auxPreparationRunnerMain <- function(argv = commandArgs(trailingOnly = TRUE)) {
    args <- auxPreparationRunnerArgs(argv)
    if (file.exists(args$summary)) {
        stop("Auxiliary preparation summary already exists", call. = FALSE)
    }
    result <- auxPublicationFreeze(
        args$output_root,
        args$contract_root,
        args$truth_root,
        args$classes
    )
    publicationWriteJson(result, args$summary)
    cat(publicationFileDigest(args$summary), "\n")
    invisible(0L)
}

if (identical(environment(), globalenv()) && !interactive()) {
    status <- tryCatch(
        auxPreparationRunnerMain(),
        error = function(error) {
            message(conditionMessage(error))
            1L
        }
    )
    quit(status = as.integer(status), save = "no")
}
