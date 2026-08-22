omicsResourceUtcNow <- function() {
    format(
        as.POSIXct(Sys.time(), tz = "UTC"),
        "%Y-%m-%dT%H:%M:%OS3Z",
        tz = "UTC"
    )
}

omicsResourceThreadEnvironment <- function() {
    as.list(Sys.getenv(
        c(
            "OMP_NUM_THREADS",
            "OPENBLAS_NUM_THREADS",
            "MKL_NUM_THREADS",
            "ARROW_NUM_THREADS",
            "DUCKDB_THREADS"
        ),
        unset = NA_character_
    ))
}

omicsResourceWriteJson <- function(value, path) {
    dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
    jsonlite::write_json(
        value,
        path,
        auto_unbox = TRUE,
        pretty = TRUE,
        null = "null",
        na = "null",
        digits = 17
    )
    invisible(path)
}

omicsResourceDisplayCommandArgs <- function(command_args, run_dir, repo_root) {
    replacements <- c(
        "<repo-root>" = normalizePath(repo_root, mustWork = FALSE),
        "<run-dir>" = normalizePath(run_dir, mustWork = FALSE)
    )
    vapply(command_args, \(argument) {
        displayed <- argument
        for (label in names(replacements)) {
            displayed <- gsub(
                replacements[[label]],
                label,
                displayed,
                fixed = TRUE
            )
        }
        displayed
    }, character(1), USE.NAMES = FALSE)
}

omicsResourceDirMetrics <- function(path, categories) {
    files <- if (dir.exists(path)) {
        list.files(
            path,
            recursive = TRUE,
            full.names = TRUE,
            all.files = TRUE,
            no.. = TRUE
        )
    } else {
        character()
    }
    files <- files[file.exists(files) & !dir.exists(files)]
    relative <- if (length(files)) {
        substring(files, nchar(normalizePath(path)) + 2L)
    } else {
        character()
    }
    sizes <- if (length(files)) as.numeric(file.info(files)$size) else numeric()
    category_names <- vapply(categories, `[[`, character(1), "category")
    bytes <- stats::setNames(as.list(rep(0, length(categories))), category_names)
    counts <- stats::setNames(as.list(rep(0L, length(categories))), category_names)
    assigned <- rep(FALSE, length(files))

    for (category in categories) {
        matches <- !assigned & grepl(category$pattern, relative, perl = TRUE)
        bytes[[category$category]] <- sum(sizes[matches], na.rm = TRUE)
        counts[[category$category]] <- sum(matches)
        assigned[matches] <- TRUE
    }
    list(
        total_bytes = sum(sizes, na.rm = TRUE),
        file_count = length(files),
        bytes = bytes,
        counts = counts
    )
}

omicsResourceProcessTreeMetrics <- function(pid) {
    if (!requireNamespace("ps", quietly = TRUE)) {
        stop("Package 'ps' is required for process-tree RSS measurement", call. = FALSE)
    }
    root <- tryCatch(ps::ps_handle(pid), error = \(...) NULL)
    running <- !is.null(root) && isTRUE(tryCatch(
        ps::ps_is_running(root),
        error = \(...) FALSE
    ))
    if (!running) {
        return(list(
            rss_bytes = 0,
            vms_bytes = 0,
            process_count = 0L,
            pids = list()
        ))
    }
    handles <- c(
        list(root),
        tryCatch(ps::ps_children(root, recursive = TRUE), error = \(...) list())
    )
    memory <- lapply(handles, \(handle) {
        tryCatch(
            ps::ps_memory_info(handle),
            error = \(...) c(rss = 0, vms = 0)
        )
    })
    list(
        rss_bytes = sum(vapply(
            memory,
            \(info) as.numeric(info[["rss"]]),
            numeric(1)
        )),
        vms_bytes = sum(vapply(
            memory,
            \(info) as.numeric(info[["vms"]]),
            numeric(1)
        )),
        process_count = length(handles),
        pids = as.list(vapply(handles, ps::ps_pid, integer(1)))
    )
}

omicsResourceCompactSamples <- function(samples, maximum_retained_samples) {
    sample_count <- length(samples)
    if (!is.finite(maximum_retained_samples) ||
        sample_count <= maximum_retained_samples) {
        return(list(samples = samples, truncated = FALSE))
    }
    maximum_retained_samples <- max(2L, as.integer(maximum_retained_samples))
    retained_indexes <- unique(as.integer(round(seq(
        1,
        sample_count,
        length.out = maximum_retained_samples
    ))))
    list(samples = samples[retained_indexes], truncated = TRUE)
}

omicsResourceMeasureSubprocess <- function(
    command,
    command_args,
    run_dir,
    execution,
    categories,
    repo_root,
    env = character(),
    capture_output = TRUE,
    timeout_ms = Inf,
    maximum_retained_samples = Inf,
    retention_marker = "retention-ready",
    require_retention_marker = FALSE
) {
    if (!requireNamespace("processx", quietly = TRUE)) {
        stop("Package 'processx' is required for fresh-process measurement", call. = FALSE)
    }
    dir.create(run_dir, recursive = TRUE, showWarnings = FALSE)
    stdout_path <- file.path(run_dir, "final", "stdout.log")
    stderr_path <- file.path(run_dir, "final", "stderr.log")
    dir.create(dirname(stdout_path), recursive = TRUE, showWarnings = FALSE)
    stdout_destination <- if (capture_output) stdout_path else nullfile()
    stderr_destination <- if (capture_output) stderr_path else nullfile()
    started_at <- omicsResourceUtcNow()
    started <- proc.time()[["elapsed"]]
    process <- processx::process$new(
        command,
        command_args,
        stdout = stdout_destination,
        stderr = stderr_destination,
        env = env,
        wd = repo_root,
        cleanup = TRUE,
        supervise = TRUE
    )
    samples <- list()
    peak_category_bytes <- stats::setNames(
        rep(0, length(categories)),
        vapply(categories, `[[`, character(1), "category")
    )
    retained <- NULL
    timed_out <- FALSE

    repeat {
        tree <- omicsResourceProcessTreeMetrics(process$get_pid())
        disk <- omicsResourceDirMetrics(run_dir, categories)
        for (category in names(peak_category_bytes)) {
            peak_category_bytes[[category]] <- max(
                peak_category_bytes[[category]],
                disk$bytes[[category]]
            )
        }
        sample <- list(
            elapsed_seconds = proc.time()[["elapsed"]] - started,
            tree_rss_bytes = tree$rss_bytes,
            tree_vms_bytes = tree$vms_bytes,
            process_count = tree$process_count,
            pids = tree$pids,
            disk_bytes = disk$total_bytes,
            file_count = disk$file_count,
            disk_category_bytes = disk$bytes
        )
        samples[[length(samples) + 1L]] <- sample
        marker_path <- file.path(run_dir, retention_marker)
        if (file.exists(marker_path) && tree$rss_bytes > 0) {
            retained <- sample
        }
        elapsed_ms <- 1000 * (proc.time()[["elapsed"]] - started)
        if (is.finite(timeout_ms) && elapsed_ms > timeout_ms && process$is_alive()) {
            timed_out <- TRUE
            process$kill()
        }
        if (!process$is_alive()) {
            break
        }
        Sys.sleep(as.numeric(execution$sampling_interval_ms) / 1000)
    }
    process$wait(timeout = 1000)
    final_disk <- omicsResourceDirMetrics(run_dir, categories)
    rss_values <- vapply(samples, `[[`, numeric(1), "tree_rss_bytes")
    vms_values <- vapply(samples, `[[`, numeric(1), "tree_vms_bytes")
    disk_values <- vapply(samples, `[[`, numeric(1), "disk_bytes")
    process_values <- vapply(samples, `[[`, integer(1), "process_count")
    artifact_categories <- intersect(
        c("committed", "staging_snapshot", "duckdb_spill"),
        names(peak_category_bytes)
    )
    artifact_disk_values <- vapply(samples, \(sample) {
        sum(unlist(
            sample$disk_category_bytes[artifact_categories],
            use.names = FALSE
        ))
    }, numeric(1))
    retention_marker_observed <- !is.null(retained)
    if (is.null(retained)) {
        retained <- samples[[length(samples)]]
    }
    compacted <- omicsResourceCompactSamples(
        samples,
        maximum_retained_samples
    )
    result <- list(
        status = if (!timed_out &&
            (!require_retention_marker || retention_marker_observed) &&
            identical(
            process$get_exit_status(),
            0L
        )) "passed" else "failed",
        exit_status = process$get_exit_status(),
        started_at = started_at,
        finished_at = omicsResourceUtcNow(),
        command = basename(command),
        arguments = as.list(omicsResourceDisplayCommandArgs(
            command_args,
            run_dir,
            repo_root
        )),
        stdout = if (capture_output) file.path("final", "stdout.log") else NULL,
        stderr = if (capture_output) file.path("final", "stderr.log") else NULL,
        output_retained = capture_output,
        metrics = list(
            elapsed_seconds = proc.time()[["elapsed"]] - started,
            peak_tree_rss_bytes = max(rss_values, na.rm = TRUE),
            retained_tree_rss_bytes = retained$tree_rss_bytes,
            peak_tree_vms_bytes = max(vms_values, na.rm = TRUE),
            maximum_process_count = max(process_values, na.rm = TRUE),
            peak_disk_bytes = max(disk_values, na.rm = TRUE),
            peak_artifact_disk_bytes = max(artifact_disk_values, na.rm = TRUE),
            peak_disk_category_bytes = as.list(peak_category_bytes),
            final_disk_bytes = final_disk$total_bytes,
            final_disk_category_bytes = final_disk$bytes,
            final_file_count = final_disk$file_count,
            sample_count = length(samples)
        ),
        samples = compacted$samples
    )
    if (timed_out) {
        result$failure <- "timeout"
    }
    if (require_retention_marker) {
        result$retention_marker_observed <- retention_marker_observed
        if (!retention_marker_observed && is.null(result$failure)) {
            result$failure <- "retention_marker_not_observed"
        }
    }
    if (compacted$truncated) {
        result$samples_truncated <- TRUE
        result$retained_sample_count <- length(compacted$samples)
    }
    result
}

omicsResourceProcSelfRss <- function() {
    status_path <- sprintf("/proc/%d/status", Sys.getpid())
    if (!file.exists(status_path)) {
        return(NA_real_)
    }
    line <- grep("^VmRSS:", readLines(status_path, warn = FALSE), value = TRUE)
    if (!length(line)) {
        return(NA_real_)
    }
    as.numeric(gsub("[^0-9]", "", line[[1L]])) * 1024
}

omicsResourceNativeMetrics <- function() {
    arrow_available <- requireNamespace("arrow", quietly = TRUE)
    duckdb_available <- requireNamespace("duckdb", quietly = TRUE)
    arrow_bytes <- if (arrow_available) {
        tryCatch(
            as.numeric(arrow::default_memory_pool()$bytes_allocated),
            error = \(...) NA_real_
        )
    } else {
        NA_real_
    }
    list(
        arrow = list(
            available = arrow_available,
            loaded = "arrow" %in% loadedNamespaces(),
            bytes_allocated = arrow_bytes
        ),
        duckdb = list(
            available = duckdb_available,
            loaded = "duckdb" %in% loadedNamespaces(),
            connections_opened = 0L
        )
    )
}

omicsResourceObservedNativeMetrics <- function() {
    loaded <- loadedNamespaces()
    arrow_loaded <- "arrow" %in% loaded
    duckdb_loaded <- "duckdb" %in% loaded
    arrow_bytes <- if (arrow_loaded) {
        tryCatch(
            as.numeric(arrow::default_memory_pool()$bytes_allocated),
            error = \(...) NA_real_
        )
    } else {
        NA_real_
    }
    list(
        arrow = list(
            available = nzchar(find.package("arrow", quiet = TRUE)),
            loaded = arrow_loaded,
            bytes_allocated = arrow_bytes
        ),
        duckdb = list(
            available = nzchar(find.package("duckdb", quiet = TRUE)),
            loaded = duckdb_loaded,
            connections_opened = 0L
        )
    )
}

omicsResourceAllocationMetrics <- function(path, sha256_file) {
    if (!file.exists(path)) {
        return(list(
            available = FALSE,
            allocated_bytes = NA_real_,
            allocation_records = 0L
        ))
    }
    lines <- readLines(path, warn = FALSE)
    bytes <- suppressWarnings(as.numeric(sub("^([0-9]+).*$", "\\1", lines)))
    list(
        available = TRUE,
        allocated_bytes = sum(bytes[is.finite(bytes)]),
        allocation_records = sum(is.finite(bytes)),
        profile_sha256 = sha256_file(path)
    )
}

omicsResourceInstallPackage <- function(repo_root, output_dir) {
    library_path <- tempfile("multischolar-measurement-library-")
    log_path <- file.path(output_dir, "harness", "package-install.log")
    dir.create(library_path, recursive = TRUE, showWarnings = FALSE)
    dir.create(dirname(log_path), recursive = TRUE, showWarnings = FALSE)
    status <- system2(
        file.path(R.home("bin"), "R"),
        c(
            "CMD",
            "INSTALL",
            "--no-test-load",
            paste0("--library=", shQuote(library_path)),
            shQuote(repo_root)
        ),
        stdout = log_path,
        stderr = log_path
    )
    installed <- file.path(library_path, "MultiScholaR")
    if (!identical(as.integer(status), 0L) || !dir.exists(installed)) {
        stop("measurement package installation failed", call. = FALSE)
    }
    normalizePath(library_path, mustWork = TRUE)
}

omicsResourceWorkingTreeDigest <- function(repo_root) {
    tracked_diff <- system2(
        "git",
        c(
            "-C",
            repo_root,
            "diff",
            "--no-ext-diff",
            "--binary",
            "HEAD",
            "--",
            "DESCRIPTION",
            "R",
            "tools/profiling"
        ),
        stdout = TRUE,
        stderr = FALSE
    )
    untracked <- system2(
        "git",
        c(
            "-C",
            repo_root,
            "ls-files",
            "--others",
            "--exclude-standard",
            "--",
            "DESCRIPTION",
            "R",
            "tools/profiling"
        ),
        stdout = TRUE,
        stderr = FALSE
    )
    untracked <- sort(untracked[nzchar(untracked)])
    untracked_digests <- vapply(
        file.path(repo_root, untracked),
        \(path) digest::digest(file = path, algo = "sha256", serialize = FALSE),
        character(1)
    )
    digest::digest(
        paste(
            c(
                paste(tracked_diff, collapse = "\n"),
                paste(
                    untracked,
                    untracked_digests,
                    sep = "=",
                    collapse = "\n"
                )
            ),
            collapse = "\n--untracked--\n"
        ),
        algo = "sha256",
        serialize = FALSE
    )
}

omicsResourceCodeRevision <- function(repo_root) {
    git_sha <- tryCatch(
        system2(
            "git",
            c("-C", repo_root, "rev-parse", "HEAD"),
            stdout = TRUE,
            stderr = FALSE
        )[[1L]],
        error = \(...) NA_character_
    )
    list(
        git_sha = git_sha,
        working_tree_sha256 = omicsResourceWorkingTreeDigest(repo_root)
    )
}

# Compatibility aliases retain the established DIA baseline API and schema.
baselineUtcNow <- omicsResourceUtcNow
baselineThreadEnvironment <- omicsResourceThreadEnvironment
baselineWriteJson <- omicsResourceWriteJson
baselineDirMetrics <- omicsResourceDirMetrics
baselineProcessTreeMetrics <- omicsResourceProcessTreeMetrics
baselineProcSelfRss <- omicsResourceProcSelfRss
baselineNativeMetrics <- omicsResourceNativeMetrics

baselineDisplayCommandArgs <- function(command_args, run_dir) {
    omicsResourceDisplayCommandArgs(
        command_args,
        run_dir,
        .BASELINE_REPO_ROOT
    )
}

baselineMeasureSubprocess <- function(
    command,
    command_args,
    run_dir,
    execution,
    categories,
    env = character(),
    capture_output = TRUE
) {
    omicsResourceMeasureSubprocess(
        command = command,
        command_args = command_args,
        run_dir = run_dir,
        execution = execution,
        categories = categories,
        repo_root = .BASELINE_REPO_ROOT,
        env = env,
        capture_output = capture_output
    )
}

baselineAllocationMetrics <- function(path) {
    omicsResourceAllocationMetrics(path, omicsParitySha256File)
}

baselineInstallMeasurementPackage <- function(output_dir) {
    omicsResourceInstallPackage(.BASELINE_REPO_ROOT, output_dir)
}

baselineWorkingTreeDigest <- function(repo_root = .BASELINE_REPO_ROOT) {
    omicsResourceWorkingTreeDigest(repo_root)
}

baselineCodeRevision <- function(repo_root = .BASELINE_REPO_ROOT) {
    omicsResourceCodeRevision(repo_root)
}
