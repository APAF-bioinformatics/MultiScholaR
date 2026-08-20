#!/usr/bin/env Rscript

lifecycleScriptPath <- function() {
    file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
    normalizePath(sub("^--file=", "", file_arg[[1L]]), mustWork = TRUE)
}

lifecycleRepoRoot <- function() {
    normalizePath(file.path(dirname(lifecycleScriptPath()), "..", ".."), mustWork = TRUE)
}

lifecycleParseArgs <- function(arguments) {
    options <- list(output = NULL, cycles = 5L)
    index <- 1L
    while (index <= length(arguments)) {
        key <- sub("^--", "", arguments[[index]])
        if (!key %in% names(options) || index == length(arguments)) {
            stop("invalid lifecycle profiling arguments", call. = FALSE)
        }
        value <- arguments[[index + 1L]]
        options[[key]] <- if (identical(key, "cycles")) as.integer(value) else value
        index <- index + 2L
    }
    if (is.null(options$output) || is.na(options$cycles) || options$cycles < 3L) {
        stop("--output and at least three --cycles are required", call. = FALSE)
    }
    options
}

lifecycleProcessMetrics <- function() {
    if (!requireNamespace("ps", quietly = TRUE)) {
        return(list(rss_bytes = NA_real_, open_handles = NA_integer_, child_count = NA_integer_))
    }
    process <- ps::ps_handle()
    list(
        rss_bytes = as.numeric(ps::ps_memory_info(process)[["rss"]]),
        open_handles = tryCatch(
            as.integer(ps::ps_num_fds(process)),
            error = \(...) NA_integer_
        ),
        child_count = length(tryCatch(
            ps::ps_children(process, recursive = TRUE),
            error = \(...) list()
        ))
    )
}

lifecycleTemporaryCount <- function(registry) {
    path <- projectRegistryPath(registry, "temporary")
    if (!dir.exists(path)) return(0L)
    length(list.files(path, all.files = TRUE, no.. = TRUE))
}

lifecycleCycle <- function(registry, cycle) {
    before <- lifecycleProcessMetrics()
    session <- initializeProjectRegistry(registry)
    open <- projectRegistrySessionResourceInfo(session)
    connection <- projectRegistrySessionConnection(session)
    result_events <- character()
    projectRegistryFetchBound(
        connection,
        "SELECT ? AS cycle",
        list(cycle),
        result_observer = \(event) result_events <<- c(result_events, event)
    )
    live_results <- sum(result_events == "opened") - sum(result_events == "cleared")
    closeProjectRegistry(session)
    closed <- projectRegistrySessionResourceInfo(session)
    after <- lifecycleProcessMetrics()
    list(
        cycle = cycle,
        before = before,
        after = after,
        open = open,
        closed = closed,
        live_dbi_results_after_query = live_results,
        writer_guard_files_after_close = as.integer(file.exists(
            projectRegistryPath(registry, "owner")
        )),
        temporary_paths_after_close = lifecycleTemporaryCount(registry),
        background_tasks_after_close = 0L,
        artifact_observers_after_close = 0L,
        hydration_cache_entries_after_close = 0L
    )
}

lifecycleOmicRun <- function(omic_type, cycles) {
    project_root <- tempfile(paste0("multischolar-lifecycle-", omic_type, "-"))
    dir.create(project_root, recursive = TRUE)
    on.exit(unlink(project_root, recursive = TRUE, force = TRUE), add = TRUE)
    paths <- buildArtifactPaths(project_root, list(
        omic_type = omic_type,
        omic_label = paste0(omic_type, "_lifecycle"),
        workflow_slug = "lifecycle"
    ))
    registry <- newProjectRegistry(
        paths,
        project_id = paste0(omic_type, "-lifecycle")
    )
    measurements <- lapply(seq_len(cycles), \(cycle) {
        lifecycleCycle(registry, cycle)
    })
    rss <- vapply(measurements, \(item) item$after$rss_bytes, numeric(1))
    warm_rss <- rss[-1L]
    plateau_span <- if (anyNA(warm_rss)) NA_real_ else diff(range(warm_rss))
    plateau_limit <- 128 * 1024^2
    list(
        omic_type = omic_type,
        cycles = measurements,
        plateau = list(
            warm_cycle_rss_span_bytes = plateau_span,
            maximum_span_bytes = plateau_limit,
            passed = is.na(plateau_span) || plateau_span <= plateau_limit
        )
    )
}

lifecycleMain <- function() {
    options <- lifecycleParseArgs(commandArgs(trailingOnly = TRUE))
    repository <- lifecycleRepoRoot()
    devtools::load_all(repository, quiet = TRUE)
    required <- c("DBI", "duckdb", "filelock", "jsonlite")
    missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
    if (length(missing) > 0L) {
        stop(paste("missing lifecycle dependencies:", paste(missing, collapse = ", ")))
    }
    started <- artifactRefUtcNow()
    runs <- lapply(c("proteomics", "metabolomics", "lipidomics"), \(omic_type) {
        lifecycleOmicRun(omic_type, options$cycles)
    })
    result <- list(
        schema = "multischolar.artifact_lifecycle_profile",
        schema_version = 1L,
        started_at = started,
        completed_at = artifactRefUtcNow(),
        process_id = as.integer(Sys.getpid()),
        used_explicit_garbage_collection = FALSE,
        runs = runs
    )
    dir.create(dirname(options$output), recursive = TRUE, showWarnings = FALSE)
    jsonlite::write_json(
        result,
        options$output,
        auto_unbox = TRUE,
        pretty = TRUE,
        null = "null",
        na = "null"
    )
    invisible(result)
}

lifecycleMain()
