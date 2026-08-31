#!/usr/bin/env Rscript

metabCleanupRunnerPath <- function() {
    file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
    if (length(file_arg)) {
        return(normalizePath(sub("^--file=", "", file_arg[[1L]]), mustWork = TRUE))
    }
    normalizePath(
        file.path(
            getwd(),
            "tools",
            "profiling",
            "run_omics_publication_metabolomics_cleanup.R"
        ),
        mustWork = TRUE
    )
}

.METAB_CLEANUP_REPO_ROOT <- normalizePath(
    file.path(dirname(metabCleanupRunnerPath()), "..", ".."),
    mustWork = TRUE
)

for (source_file in c(
    "omics_publication_protocol.R",
    "omics_publication_metabolomics_contracts.R",
    "omics_publication_metabolomics_cleanup.R"
)) {
    source(
        file.path(
            .METAB_CLEANUP_REPO_ROOT,
            "tools",
            "profiling",
            source_file
        ),
        local = FALSE
    )
}

metabCleanupRunnerArgs <- function(argv) {
    values <- list(plan = NULL, output = NULL, execute = FALSE)
    index <- 1L
    while (index <= length(argv)) {
        key <- gsub("-", "_", sub("^--", "", argv[[index]]), fixed = TRUE)
        if (!key %in% names(values) || index == length(argv)) {
            stop("Cleanup runner arguments are invalid", call. = FALSE)
        }
        value <- argv[[index + 1L]]
        values[[key]] <- if (identical(key, "execute")) {
            identical(tolower(value), "true")
        } else {
            value
        }
        index <- index + 2L
    }
    if (is.null(values$plan) || is.null(values$output)) {
        stop("Cleanup plan and output are required", call. = FALSE)
    }
    values
}

metabCleanupRunnerMain <- function(argv = commandArgs(trailingOnly = TRUE)) {
    args <- metabCleanupRunnerArgs(argv)
    if (file.exists(args$output)) {
        stop("Cleanup output already exists", call. = FALSE)
    }
    plan <- publicationReadJson(args$plan)
    result <- metabPublicationRunCleanup(plan, execute = args$execute)
    publicationWriteJson(result, args$output)
    cat(publicationFileDigest(args$output), "\n")
    invisible(0L)
}

if (identical(environment(), globalenv()) && !interactive()) {
    status <- tryCatch(
        metabCleanupRunnerMain(),
        error = \(error) {
            message(conditionMessage(error))
            1L
        }
    )
    quit(status = as.integer(status), save = "no")
}
