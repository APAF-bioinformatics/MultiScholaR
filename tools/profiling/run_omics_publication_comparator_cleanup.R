#!/usr/bin/env Rscript

cleanupRunnerPath <- function() {
    file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
    if (length(file_arg)) {
        return(normalizePath(sub("^--file=", "", file_arg[[1L]]), mustWork = TRUE))
    }
    normalizePath(
        file.path(
            getwd(),
            "tools",
            "profiling",
            "run_omics_publication_comparator_cleanup.R"
        ),
        mustWork = TRUE
    )
}

.PUBLICATION_CLEANUP_ROOT <- normalizePath(
    file.path(dirname(cleanupRunnerPath()), "..", ".."),
    mustWork = TRUE
)

for (source_file in c(
    "omics_publication_protocol.R",
    "omics_publication_comparators.R",
    "omics_publication_builds.R",
    "omics_publication_comparator_cleanup.R"
)) {
    source(
        file.path(
            .PUBLICATION_CLEANUP_ROOT,
            "tools",
            "profiling",
            source_file
        ),
        local = FALSE
    )
}

cleanupRunnerArgs <- function(argv) {
    values <- list(plan = NULL, output = NULL, execute = FALSE)
    index <- 1L
    while (index <= length(argv)) {
        key <- sub("^--", "", argv[[index]])
        if (!key %in% names(values) || index == length(argv)) {
            stop("Usage: --plan <json> --output <json> --execute <true|false>")
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

cleanupRunnerMain <- function(argv = commandArgs(trailingOnly = TRUE)) {
    args <- cleanupRunnerArgs(argv)
    if (file.exists(args$output)) {
        stop("Cleanup output already exists", call. = FALSE)
    }
    plan <- publicationReadJson(args$plan)
    result <- publicationRunComparatorCleanup(plan, args$execute)
    publicationWriteJson(result, args$output)
    cat(publicationFileDigest(args$output), "\n")
    invisible(0L)
}

if (identical(environment(), globalenv()) && !interactive()) {
    status <- tryCatch(
        cleanupRunnerMain(),
        error = \(error) {
            message(conditionMessage(error))
            1L
        }
    )
    quit(status = as.integer(status), save = "no")
}
