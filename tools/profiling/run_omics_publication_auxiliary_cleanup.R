#!/usr/bin/env Rscript

auxCleanupRunnerPath <- function() {
    file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
    if (length(file_arg)) {
        return(normalizePath(sub("^--file=", "", file_arg[[1L]]), mustWork = TRUE))
    }
    normalizePath(
        file.path(
            getwd(),
            "tools",
            "profiling",
            "run_omics_publication_auxiliary_cleanup.R"
        ),
        mustWork = TRUE
    )
}

.AUX_CLEANUP_REPO_ROOT <- normalizePath(
    file.path(dirname(auxCleanupRunnerPath()), "..", ".."),
    mustWork = TRUE
)

for (source_file in c(
    "omics_publication_protocol.R",
    "omics_publication_auxiliary_contracts.R",
    "omics_publication_auxiliary_cleanup.R"
)) {
    source(
        file.path(.AUX_CLEANUP_REPO_ROOT, "tools", "profiling", source_file),
        local = FALSE
    )
}

auxCleanupRunnerArgs <- function(argv) {
    values <- list(plan = NULL, output = NULL, dry_run = "true")
    index <- 1L
    while (index <= length(argv)) {
        token <- argv[[index]]
        key <- gsub("-", "_", sub("^--", "", token), fixed = TRUE)
        if (!startsWith(token, "--") || !key %in% names(values) ||
            index == length(argv)) {
            stop("Auxiliary cleanup arguments are invalid", call. = FALSE)
        }
        values[[key]] <- argv[[index + 1L]]
        index <- index + 2L
    }
    if (is.null(values$plan) || is.null(values$output)) {
        stop("Auxiliary cleanup options are incomplete", call. = FALSE)
    }
    values$dry_run <- identical(tolower(values$dry_run), "true")
    values
}

auxCleanupRunnerMain <- function(argv = commandArgs(trailingOnly = TRUE)) {
    args <- auxCleanupRunnerArgs(argv)
    if (file.exists(args$output)) {
        stop("Auxiliary cleanup output already exists", call. = FALSE)
    }
    plan <- publicationReadJson(args$plan)
    result <- auxPublicationApplyCleanup(plan, args$dry_run)
    publicationWriteJson(result, args$output)
    cat(publicationFileDigest(args$output), "\n")
    invisible(0L)
}

if (identical(environment(), globalenv()) && !interactive()) {
    status <- tryCatch(
        auxCleanupRunnerMain(),
        error = function(error) {
            message(conditionMessage(error))
            1L
        }
    )
    quit(status = as.integer(status), save = "no")
}
