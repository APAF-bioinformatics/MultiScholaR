#!/usr/bin/env Rscript

publicationComparatorRunnerPath <- function() {
    file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
    if (length(file_arg)) {
        return(normalizePath(sub("^--file=", "", file_arg[[1L]]), mustWork = TRUE))
    }
    normalizePath(
        file.path(
            getwd(),
            "tools",
            "profiling",
            "run_omics_publication_comparator.R"
        ),
        mustWork = TRUE
    )
}

.PUBLICATION_COMPARATOR_RUNNER_ROOT <- normalizePath(
    file.path(dirname(publicationComparatorRunnerPath()), "..", ".."),
    mustWork = TRUE
)

for (source_file in c(
    "omics_publication_protocol.R",
    "omics_publication_comparators.R",
    "omics_publication_lock.R",
    "omics_publication_builds.R",
    "omics_publication_repository_inputs.R",
    "omics_publication_remote_installs.R",
    "omics_publication_restore_reproducibility.R",
    "omics_publication_comparator_builds.R",
    "omics_publication_comparator_build_reproducibility.R",
    "omics_publication_comparator_evidence.R",
    "omics_publication_comparator_envelopes.R"
)) {
    source(
        file.path(
            .PUBLICATION_COMPARATOR_RUNNER_ROOT,
            "tools",
            "profiling",
            source_file
        ),
        local = FALSE
    )
}

publicationComparatorRunnerArgs <- function() {
    list(
        comparator_id = NULL,
        backend = NULL,
        build_receipt = NULL,
        attempt = 1L,
        output = NULL,
        help = FALSE
    )
}

publicationComparatorRunnerParse <- function(argv) {
    args <- publicationComparatorRunnerArgs()
    index <- 1L
    while (index <= length(argv)) {
        token <- argv[[index]]
        if (identical(token, "--help")) {
            args$help <- TRUE
            index <- index + 1L
            next
        }
        if (!startsWith(token, "--")) {
            stop(paste("Unexpected argument:", token), call. = FALSE)
        }
        key <- gsub("-", "_", sub("^--", "", token), fixed = TRUE)
        if (!key %in% names(args)) {
            stop(paste("Unknown option:", token), call. = FALSE)
        }
        index <- index + 1L
        if (index > length(argv)) {
            stop(paste("Missing value for:", token), call. = FALSE)
        }
        value <- argv[[index]]
        args[[key]] <- if (identical(key, "attempt")) {
            as.integer(value)
        } else {
            value
        }
        index <- index + 1L
    }
    args
}

publicationComparatorRunnerUsage <- function() {
    cat(paste(
        "Usage: Rscript --vanilla tools/profiling/run_omics_publication_comparator.R",
        "  --comparator-id <id> --backend <backend>",
        "  [--build-receipt <json>] [--attempt <n>] --output <json>",
        sep = "\n"
    ))
}

publicationComparatorRunnerMain <- function(argv = commandArgs(trailingOnly = TRUE)) {
    args <- publicationComparatorRunnerParse(argv)
    if (isTRUE(args$help)) {
        publicationComparatorRunnerUsage()
        return(invisible(0L))
    }
    required <- c("comparator_id", "backend", "output")
    missing <- required[vapply(required, \(name) {
        is.null(args[[name]])
    }, logical(1))]
    if (length(missing)) {
        stop(paste("Missing comparator option:", missing[[1L]]), call. = FALSE)
    }
    if (file.exists(args$output)) {
        stop("Comparator runner output already exists", call. = FALSE)
    }
    result <- publicationRunComparatorTinyEnvelope(
        args$comparator_id,
        args$backend,
        args$build_receipt,
        attempt = args$attempt
    )
    publicationWriteJson(result, args$output)
    invisible(0L)
}

if (identical(environment(), globalenv()) && !interactive()) {
    status <- tryCatch(
        publicationComparatorRunnerMain(),
        error = \(error) {
            message(conditionMessage(error))
            1L
        }
    )
    quit(status = as.integer(status), save = "no")
}
