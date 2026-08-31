publicationReadLinesSafe <- function(path) {
    if (!file.exists(path)) return(character())
    tryCatch(
        suppressWarnings(readLines(path, warn = FALSE)),
        error = \(...) character()
    )
}

publicationReadScalarFile <- function(path) {
    lines <- publicationReadLinesSafe(path)
    if (!length(lines)) return(NA_real_)
    value <- suppressWarnings(as.numeric(trimws(lines[[1L]])))
    if (length(value) != 1L || !is.finite(value)) NA_real_ else value
}

publicationMonotonicSeconds <- function() {
    lines <- publicationReadLinesSafe("/proc/uptime")
    if (!length(lines)) return(NA_real_)
    value <- suppressWarnings(as.numeric(strsplit(
        trimws(lines[[1L]]),
        "[[:space:]]+",
        perl = TRUE
    )[[1L]][[1L]]))
    if (publicationScalarNumber(value)) value else NA_real_
}

publicationParseNamedNumbers <- function(lines) {
    if (!length(lines)) return(list())
    pieces <- strsplit(trimws(lines), "[[:space:]]+", perl = TRUE)
    values <- lapply(pieces, \(piece) {
        if (length(piece) < 2L) return(NULL)
        value <- suppressWarnings(as.numeric(piece[[2L]]))
        if (!is.finite(value)) return(NULL)
        list(name = piece[[1L]], value = value)
    })
    values <- Filter(Negate(is.null), values)
    if (!length(values)) return(list())
    stats::setNames(
        lapply(values, `[[`, "value"),
        vapply(values, `[[`, character(1), "name")
    )
}

publicationParseIoStat <- function(lines) {
    totals <- c(rbytes = 0, wbytes = 0, rios = 0, wios = 0, dbytes = 0, dios = 0)
    if (!length(lines)) return(as.list(totals))
    for (line in lines) {
        fields <- strsplit(trimws(line), "[[:space:]]+", perl = TRUE)[[1L]]
        if (length(fields) < 2L) next
        for (field in fields[-1L]) {
            pair <- strsplit(field, "=", fixed = TRUE)[[1L]]
            if (length(pair) != 2L || !pair[[1L]] %in% names(totals)) next
            value <- suppressWarnings(as.numeric(pair[[2L]]))
            if (is.finite(value)) totals[[pair[[1L]]]] <- totals[[pair[[1L]]]] + value
        }
    }
    as.list(totals)
}

publicationParsePressure <- function(lines) {
    result <- list()
    for (line in lines) {
        fields <- strsplit(trimws(line), "[[:space:]]+", perl = TRUE)[[1L]]
        if (length(fields) < 2L) next
        values <- list()
        for (field in fields[-1L]) {
            pair <- strsplit(field, "=", fixed = TRUE)[[1L]]
            if (length(pair) != 2L) next
            value <- suppressWarnings(as.numeric(pair[[2L]]))
            if (is.finite(value)) values[[pair[[1L]]]] <- value
        }
        result[[fields[[1L]]]] <- values
    }
    result
}

publicationCgroupRoot <- function() {
    "/sys/fs/cgroup"
}

publicationCurrentCgroupRelativePath <- function(pid = "self") {
    lines <- publicationReadLinesSafe(file.path("/proc", as.character(pid), "cgroup"))
    match <- grep("^0::", lines, value = TRUE)
    if (length(match) != 1L) return(NULL)
    sub("^0::", "", match[[1L]])
}

publicationCurrentCgroupPath <- function(pid = "self") {
    relative <- publicationCurrentCgroupRelativePath(pid)
    if (is.null(relative)) return(NULL)
    file.path(publicationCgroupRoot(), sub("^/", "", relative))
}

publicationUserAppSlicePath <- function() {
    relative <- publicationCurrentCgroupRelativePath()
    if (is.null(relative)) return(NULL)
    pieces <- strsplit(sub("^/", "", relative), "/", fixed = TRUE)[[1L]]
    index <- which(pieces == "app.slice")
    if (length(index) != 1L) return(NULL)
    do.call(
        file.path,
        as.list(c(publicationCgroupRoot(), pieces[seq_len(index)]))
    )
}

publicationSystemdUnitName <- function(prefix = "multischolar-publication") {
    token <- substr(digest::digest(
        paste(Sys.getpid(), publicationUtcNow(), stats::runif(1L)),
        algo = "sha256",
        serialize = FALSE
    ), 1L, 12L)
    paste0(gsub("[^A-Za-z0-9-]", "-", prefix), "-", token)
}

publicationExpectedUnitCgroupPath <- function(unit_name) {
    app_slice <- publicationUserAppSlicePath()
    if (is.null(app_slice)) return(NULL)
    file.path(app_slice, paste0(unit_name, ".service"))
}

publicationSystemdUnitMainPid <- function(unit_name) {
    result <- processx::run(
        "systemctl",
        c(
            "--user", "show", paste0(unit_name, ".service"),
            "--property=MainPID", "--value"
        ),
        error_on_status = FALSE,
        echo = FALSE,
        timeout = 5000
    )
    if (result$status != 0L) return(NA_integer_)
    value <- suppressWarnings(as.integer(trimws(result$stdout)))
    if (length(value) == 1L && is.finite(value) && value > 0L) value else
        NA_integer_
}

publicationProcessParentPid <- function(pid) {
    lines <- publicationReadLinesSafe(file.path("/proc", pid, "stat"))
    if (length(lines) != 1L) return(NA_integer_)
    close <- regexpr("\\)[[:space:]][A-Z][[:space:]]", lines[[1L]], perl = TRUE)
    if (close[[1L]] < 0L) return(NA_integer_)
    suffix <- substring(lines[[1L]], close[[1L]] + attr(close, "match.length"))
    fields <- strsplit(trimws(suffix), "[[:space:]]+", perl = TRUE)[[1L]]
    value <- suppressWarnings(as.integer(fields[[1L]]))
    if (length(value) == 1L && is.finite(value)) value else NA_integer_
}

publicationPidOwnedByRoot <- function(pid, root_pid, maximum_depth = 128L) {
    pid <- as.integer(pid)
    root_pid <- as.integer(root_pid)
    if (!is.finite(pid) || !is.finite(root_pid)) return(NA)
    current <- pid
    visited <- integer()
    for (depth in seq_len(maximum_depth)) {
        if (identical(current, root_pid)) return(TRUE)
        if (current <= 1L || current %in% visited) return(FALSE)
        visited <- c(visited, current)
        current <- publicationProcessParentPid(current)
        if (!is.finite(current)) return(NA)
    }
    FALSE
}

publicationPidOwnership <- function(pids, root_pid) {
    pids <- unique(as.integer(pids))
    owned <- vapply(
        pids,
        publicationPidOwnedByRoot,
        logical(1),
        root_pid = root_pid
    )
    unexpected <- !is.na(owned) & !owned
    indeterminate <- is.na(owned)
    list(
        valid = !any(unexpected),
        root_pid = as.integer(root_pid),
        observed_pids = as.list(pids),
        unexpected_pids = as.list(pids[unexpected]),
        indeterminate_pids = as.list(pids[indeterminate])
    )
}

publicationReadSmapsRollup <- function(pid) {
    path <- file.path("/proc", as.character(pid), "smaps_rollup")
    lines <- publicationReadLinesSafe(path)
    publicationParseSmapsRollup(lines, pid)
}

publicationParseSmapsRollup <- function(lines, pid) {
    if (!length(lines)) {
        return(list(
            pid = as.integer(pid),
            available = FALSE,
            rss_bytes = NA_real_,
            pss_bytes = NA_real_,
            uss_bytes = NA_real_,
            missing_fields = list()
        ))
    }
    fields <- c("Rss", "Pss", "Private_Clean", "Private_Dirty", "Private_Hugetlb")
    values <- stats::setNames(rep(0, length(fields)), fields)
    observed <- stats::setNames(rep(FALSE, length(fields)), fields)
    for (field in fields) {
        match <- grep(paste0("^", field, ":"), lines, value = TRUE)
        if (!length(match)) next
        value <- suppressWarnings(as.numeric(gsub("[^0-9]", "", match[[1L]])))
        if (is.finite(value)) {
            values[[field]] <- value * 1024
            observed[[field]] <- TRUE
        }
    }
    if (!all(observed)) {
        return(list(
            pid = as.integer(pid),
            available = FALSE,
            rss_bytes = NA_real_,
            pss_bytes = NA_real_,
            uss_bytes = NA_real_,
            missing_fields = as.list(names(observed)[!observed])
        ))
    }
    list(
        pid = as.integer(pid),
        available = TRUE,
        rss_bytes = values[["Rss"]],
        pss_bytes = values[["Pss"]],
        uss_bytes = sum(values[c(
            "Private_Clean", "Private_Dirty", "Private_Hugetlb"
        )]),
        missing_fields = list()
    )
}

publicationCgroupPids <- function(cgroup_path) {
    values <- suppressWarnings(as.integer(publicationReadLinesSafe(
        file.path(cgroup_path, "cgroup.procs")
    )))
    values[is.finite(values)]
}

publicationAggregateSmaps <- function(pids) {
    records <- lapply(unique(as.integer(pids)), publicationReadSmapsRollup)
    available <- vapply(records, `[[`, logical(1), "available")
    if (!length(records) || !all(available)) {
        return(list(
            available = FALSE,
            rss_bytes = NA_real_,
            pss_bytes = NA_real_,
            uss_bytes = NA_real_,
            pids = as.list(unique(as.integer(pids))),
            unavailable_pids = if (length(records)) {
                as.list(vapply(records[!available], `[[`, integer(1), "pid"))
            } else {
                list()
            }
        ))
    }
    list(
        available = TRUE,
        rss_bytes = sum(vapply(records, `[[`, numeric(1), "rss_bytes")),
        pss_bytes = sum(vapply(records, `[[`, numeric(1), "pss_bytes")),
        uss_bytes = sum(vapply(records, `[[`, numeric(1), "uss_bytes")),
        pids = as.list(vapply(records, `[[`, integer(1), "pid")),
        unavailable_pids = list()
    )
}

publicationProcessResourceCounts <- function(pids) {
    pids <- unique(as.integer(pids))
    fd_counts <- vapply(pids, \(pid) {
        path <- file.path("/proc", pid, "fd")
        if (!dir.exists(path)) return(NA_real_)
        length(list.files(path, all.files = TRUE, no.. = TRUE))
    }, numeric(1))
    task_counts <- vapply(pids, \(pid) {
        path <- file.path("/proc", pid, "task")
        if (!dir.exists(path)) return(NA_real_)
        length(list.files(path, all.files = TRUE, no.. = TRUE))
    }, numeric(1))
    lock_lines <- publicationReadLinesSafe("/proc/locks")
    lock_count <- if (length(lock_lines) && length(pids)) {
        sum(vapply(lock_lines, \(line) {
            any(vapply(pids, \(pid) {
                grepl(paste0("[[:space:]]", pid, "[[:space:]]"), line)
            }, logical(1)))
        }, logical(1)))
    } else {
        0L
    }
    complete <- all(is.finite(fd_counts)) && all(is.finite(task_counts))
    list(
        available = complete,
        file_descriptors = if (complete) sum(fd_counts) else NA_real_,
        thread_tasks = if (complete) sum(task_counts) else NA_real_,
        kernel_locks = as.integer(lock_count),
        inaccessible_pids = as.list(pids[!is.finite(fd_counts) |
            !is.finite(task_counts)])
    )
}

publicationDefaultDiskCategories <- function() {
    list(
        list(category = "diagnostics", pattern = "(^|/)diagnostics(/|$)"),
        list(category = "staging_snapshot", pattern = "(^|/)(staging|snapshots)(/|$)"),
        list(category = "duckdb_spill", pattern = "(^|/)duckdb-spill(/|$)"),
        list(category = "committed", pattern = "(^|/)committed(/|$)"),
        list(category = "final", pattern = "(^|/)final(/|$)"),
        list(category = "harness", pattern = ".*")
    )
}

publicationAllocatedDirectoryBytes <- function(path) {
    output <- processx::run(
        "du",
        c("-s", "-B1", normalizePath(path, mustWork = TRUE)),
        error_on_status = FALSE,
        echo = FALSE
    )
    if (output$status != 0L) return(NA_real_)
    value <- suppressWarnings(as.numeric(strsplit(
        trimws(output$stdout),
        "[[:space:]]+",
        perl = TRUE
    )[[1L]][[1L]]))
    if (is.finite(value)) value else NA_real_
}

publicationDirectoryMetrics <- function(
    path,
    categories = publicationDefaultDiskCategories()
) {
    files <- list.files(
        path,
        recursive = TRUE,
        full.names = TRUE,
        all.files = TRUE,
        no.. = TRUE
    )
    files <- files[file.exists(files) & !dir.exists(files)]
    sizes <- if (length(files)) as.numeric(file.info(files)$size) else numeric()
    prefix <- paste0(normalizePath(path, mustWork = TRUE), .Platform$file.sep)
    relative <- if (length(files)) sub(prefix, "", normalizePath(files), fixed = TRUE) else
        character()
    category_names <- vapply(categories, `[[`, character(1), "category")
    logical_bytes <- stats::setNames(as.list(rep(0, length(categories))), category_names)
    file_counts <- stats::setNames(as.list(rep(0L, length(categories))), category_names)
    assigned <- rep(FALSE, length(files))
    for (category in categories) {
        matches <- !assigned & grepl(category$pattern, relative, perl = TRUE)
        logical_bytes[[category$category]] <- sum(sizes[matches], na.rm = TRUE)
        file_counts[[category$category]] <- sum(matches)
        assigned[matches] <- TRUE
    }
    list(
        logical_bytes = sum(sizes, na.rm = TRUE),
        allocated_bytes = publicationAllocatedDirectoryBytes(path),
        file_count = length(files),
        category_logical_bytes = logical_bytes,
        category_file_counts = file_counts
    )
}

publicationCgroupSnapshot <- function(cgroup_path) {
    if (is.null(cgroup_path) || !dir.exists(cgroup_path)) return(NULL)
    memory_stat <- publicationParseNamedNumbers(publicationReadLinesSafe(
        file.path(cgroup_path, "memory.stat")
    ))
    cpu_stat <- publicationParseNamedNumbers(publicationReadLinesSafe(
        file.path(cgroup_path, "cpu.stat")
    ))
    memory_events <- publicationParseNamedNumbers(publicationReadLinesSafe(
        file.path(cgroup_path, "memory.events")
    ))
    cgroup_events <- publicationParseNamedNumbers(publicationReadLinesSafe(
        file.path(cgroup_path, "cgroup.events")
    ))
    pids <- publicationCgroupPids(cgroup_path)
    list(
        memory = list(
            current_bytes = publicationReadScalarFile(file.path(
                cgroup_path,
                "memory.current"
            )),
            peak_bytes = publicationReadScalarFile(file.path(
                cgroup_path,
                "memory.peak"
            )),
            swap_current_bytes = publicationReadScalarFile(file.path(
                cgroup_path,
                "memory.swap.current"
            )),
            stat = memory_stat
        ),
        cpu = cpu_stat,
        io = publicationParseIoStat(publicationReadLinesSafe(file.path(
            cgroup_path,
            "io.stat"
        ))),
        pressure = list(
            memory = publicationParsePressure(publicationReadLinesSafe(file.path(
                cgroup_path,
                "memory.pressure"
            ))),
            cpu = publicationParsePressure(publicationReadLinesSafe(file.path(
                cgroup_path,
                "cpu.pressure"
            ))),
            io = publicationParsePressure(publicationReadLinesSafe(file.path(
                cgroup_path,
                "io.pressure"
            )))
        ),
        events = list(
            memory = memory_events,
            cgroup = cgroup_events
        ),
        smaps = publicationAggregateSmaps(pids),
        resources = publicationProcessResourceCounts(pids),
        process_count = length(pids),
        pids = as.list(pids)
    )
}

publicationSystemdRunArgs <- function(
    unit_name,
    command,
    command_args,
    env,
    extra_properties = character()
) {
    env_args <- if (length(env)) {
        paste0("--setenv=", names(env), "=", unname(env))
    } else {
        character()
    }
    c(
        "--user", "--wait", "--pipe",
        paste0("--unit=", unit_name),
        paste0("--working-directory=", .PUBLICATION_REPO_ROOT),
        "--property=MemoryAccounting=yes",
        "--property=IOAccounting=yes",
        "--property=TasksMax=512",
        "--property=KillMode=mixed",
        "--property=OOMPolicy=stop",
        extra_properties,
        env_args,
        command,
        command_args
    )
}

publicationStopSystemdUnit <- function(unit_name) {
    try(processx::run(
        "systemctl",
        c("--user", "stop", paste0(unit_name, ".service")),
        error_on_status = FALSE,
        echo = FALSE,
        timeout = 5000
    ), silent = TRUE)
    invisible(TRUE)
}

publicationSystemdUnitProperty <- function(unit_name, property) {
    result <- processx::run(
        "systemctl",
        c(
            "--user", "show", paste0(unit_name, ".service"),
            paste0("--property=", property), "--value"
        ),
        error_on_status = FALSE,
        echo = FALSE,
        timeout = 5000
    )
    if (result$status != 0L) return(NULL)
    trimws(result$stdout)
}

publicationSystemdUnitResult <- function(unit_name, systemd_run_exit_status) {
    result <- publicationSystemdUnitProperty(unit_name, "Result")
    main_code <- suppressWarnings(as.integer(publicationSystemdUnitProperty(
        unit_name,
        "ExecMainCode"
    )))
    main_status <- suppressWarnings(as.integer(publicationSystemdUnitProperty(
        unit_name,
        "ExecMainStatus"
    )))
    memory_peak <- suppressWarnings(as.numeric(publicationSystemdUnitProperty(
        unit_name,
        "MemoryPeak"
    )))
    property_evidence_available <- publicationScalarString(result) &&
        publicationScalarNumber(main_code) &&
        publicationScalarNumber(main_status)
    successful_wait_fallback <- !property_evidence_available &&
        identical(systemd_run_exit_status, 0L)
    list(
        available = property_evidence_available || successful_wait_fallback,
        source = if (property_evidence_available) {
            "systemd_unit_properties"
        } else if (successful_wait_fallback) {
            "systemd_run_wait_exit_status"
        } else {
            "unavailable"
        },
        result = if (successful_wait_fallback) "success" else result,
        exec_main_code = if (successful_wait_fallback) 1L else main_code,
        exec_main_status = if (successful_wait_fallback) 0L else main_status,
        memory_peak_bytes = memory_peak,
        oom_killed = identical(result, "oom-kill")
    )
}

publicationOwnedUnitCleanup <- function(unit_name, cgroup_path, timeout_seconds = 5) {
    started <- proc.time()[["elapsed"]]
    publicationStopSystemdUnit(unit_name)
    repeat {
        unit <- processx::run(
            "systemctl",
            c("--user", "is-active", paste0(unit_name, ".service")),
            error_on_status = FALSE,
            echo = FALSE,
            timeout = 5000
        )
        unit_active <- identical(trimws(unit$stdout), "active") ||
            identical(trimws(unit$stdout), "activating")
        cgroup_exists <- dir.exists(cgroup_path)
        if (!unit_active && !cgroup_exists) break
        if (proc.time()[["elapsed"]] - started >= timeout_seconds) break
        Sys.sleep(0.05)
    }
    try(processx::run(
        "systemctl",
        c("--user", "reset-failed", paste0(unit_name, ".service")),
        error_on_status = FALSE,
        echo = FALSE,
        timeout = 5000
    ), silent = TRUE)
    list(
        valid = !unit_active && !cgroup_exists,
        unit_active = unit_active,
        cgroup_exists = cgroup_exists,
        elapsed_seconds = proc.time()[["elapsed"]] - started
    )
}

publicationAcknowledgeRetention <- function(path) {
    connection <- fifo(path, open = "wb", blocking = TRUE)
    on.exit(close(connection), add = TRUE)
    writeBin(as.raw(1L), connection)
    flush(connection)
    TRUE
}

publicationCgroupExtraProperties <- function(execution) {
    if (is.null(execution$memory_max_bytes)) return(character())
    if (!publicationScalarNumber(execution$memory_max_bytes, positive = TRUE)) {
        publicationAbort(
            "Fault-injection memory limit is invalid",
            "multischolar_publication_cgroup_error"
        )
    }
    c(
        paste0("--property=MemoryMax=", execution$memory_max_bytes),
        "--property=MemorySwapMax=0"
    )
}

publicationPrepareCgroupMeasurement <- function(
    command,
    command_args,
    run_dir,
    execution,
    env,
    unit_name,
    require_certified_host,
    host_preflight
) {
    if (!publicationSystemdUserAvailable()) {
        publicationAbort(
            "Owned user systemd cgroup is unavailable",
            "multischolar_publication_cgroup_unavailable"
        )
    }
    if (isTRUE(require_certified_host) &&
        (is.null(host_preflight) || !isTRUE(host_preflight$certified))) {
        publicationAbort(
            "Publication workload requires a certified host preflight",
            "multischolar_publication_host_preflight_error"
        )
    }
    dir.create(run_dir, recursive = TRUE, showWarnings = FALSE)
    cgroup_path <- publicationExpectedUnitCgroupPath(unit_name)
    if (is.null(cgroup_path)) {
        publicationAbort(
            "Could not derive owned systemd cgroup path",
            "multischolar_publication_cgroup_unavailable"
        )
    }
    paths <- list(
        stdout = file.path(run_dir, "stdout.log"),
        stderr = file.path(run_dir, "stderr.log"),
        marker = file.path(run_dir, "retention-ready"),
        retention_state = file.path(run_dir, "retention-state.json"),
        acknowledgement = file.path(run_dir, "retention-sampled.fifo")
    )
    args <- publicationSystemdRunArgs(
        unit_name,
        command,
        command_args,
        env,
        publicationCgroupExtraProperties(execution)
    )
    started_at <- publicationUtcNow()
    started <- proc.time()[["elapsed"]]
    started_monotonic <- publicationMonotonicSeconds()
    process <- processx::process$new(
        "systemd-run",
        args,
        stdout = paths$stdout,
        stderr = paths$stderr,
        wd = .PUBLICATION_REPO_ROOT,
        cleanup = TRUE,
        supervise = TRUE
    )
    list(
        process = process,
        paths = paths,
        cgroup_path = cgroup_path,
        started_at = started_at,
        started = started,
        started_monotonic = started_monotonic
    )
}

publicationRunCgroupSampling <- function(
    prepared,
    run_dir,
    execution,
    unit_name,
    disk_categories,
    safety_check_fn
) {
    process <- prepared$process
    samples <- list()
    marker_elapsed <- NA_real_
    retention_state <- NULL
    cgroup_observed <- FALSE
    cgroup_lost <- FALSE
    expected_root_pid <- NA_integer_
    pid_ownership_valid <- TRUE
    unexpected_pids <- integer()
    indeterminate_pids <- integer()
    timed_out <- FALSE
    safety_aborted <- FALSE
    safety_reason <- NULL
    safety_trace <- list()
    safety_interval_seconds <- attr(
        safety_check_fn,
        "publication_safety_interval_seconds",
        exact = TRUE
    )
    safety_trace_required <- is.function(safety_check_fn) &&
        publicationScalarNumber(safety_interval_seconds, positive = TRUE)
    if (!publicationScalarNumber(safety_interval_seconds, positive = TRUE)) {
        safety_interval_seconds <- 1
    }
    retention_acknowledged <- FALSE
    disk_interval <- if (is.null(execution$disk_sampling_interval_ms)) {
        500
    } else {
        as.numeric(execution$disk_sampling_interval_ms)
    }
    boundary_tolerance <- if (is.null(
        execution$retained_boundary_tolerance_ms
    )) {
        0.1
    } else {
        execution$retained_boundary_tolerance_ms / 1000
    }
    last_disk <- publicationDirectoryMetrics(run_dir, disk_categories)
    last_disk_elapsed <- -Inf
    repeat {
        elapsed <- proc.time()[["elapsed"]] - prepared$started
        if (elapsed - last_disk_elapsed >= disk_interval / 1000) {
            last_disk <- publicationDirectoryMetrics(run_dir, disk_categories)
            last_disk_elapsed <- elapsed
        }
        sample <- publicationSampleCgroup(
            prepared$cgroup_path,
            elapsed,
            last_disk
        )
        if (!is.null(sample)) {
            cgroup_observed <- TRUE
            if (!is.finite(expected_root_pid)) {
                expected_root_pid <- publicationSystemdUnitMainPid(unit_name)
            }
            if (is.finite(expected_root_pid)) {
                ownership <- publicationPidOwnership(
                    unlist(sample$pids, use.names = FALSE),
                    expected_root_pid
                )
                if (!isTRUE(ownership$valid)) {
                    pid_ownership_valid <- FALSE
                    unexpected_pids <- unique(c(
                        unexpected_pids,
                        unlist(ownership$unexpected_pids, use.names = FALSE)
                    ))
                    safety_aborted <- TRUE
                    safety_reason <- "unexpected_pid_entered_owned_cgroup"
                    publicationStopSystemdUnit(unit_name)
                    process$kill()
                }
                indeterminate_pids <- unique(c(
                    indeterminate_pids,
                    unlist(ownership$indeterminate_pids, use.names = FALSE)
                ))
            }
            samples[[length(samples) + 1L]] <- sample
            if (is.function(safety_check_fn)) {
                safety <- safety_check_fn(sample, elapsed)
                if (is.list(safety) && isTRUE(safety$fresh)) {
                    safety$fresh <- NULL
                    safety_trace[[length(safety_trace) + 1L]] <- safety
                }
                if (!is.list(safety) || !isTRUE(safety$safe)) {
                    safety_aborted <- TRUE
                    safety_reason <- if (is.list(safety)) safety$reason else
                        "invalid_safety_result"
                    publicationStopSystemdUnit(unit_name)
                    process$kill()
                }
            }
        } else if (cgroup_observed && process$is_alive()) {
            dwell_complete <- is.finite(marker_elapsed) &&
                elapsed >= marker_elapsed + execution$retained_dwell_seconds -
                    boundary_tolerance
            if (!dwell_complete) {
                cgroup_lost <- TRUE
                safety_aborted <- TRUE
                safety_reason <- "owned_cgroup_disappeared"
                publicationStopSystemdUnit(unit_name)
                process$kill()
            }
        }
        if (!is.finite(marker_elapsed) && file.exists(prepared$paths$marker)) {
            retention_state <- publicationReadRetentionState(
                prepared$paths$retention_state
            )
            settled <- retention_state$settled_monotonic_seconds
            marker_elapsed <- if (publicationScalarNumber(settled) &&
                publicationScalarNumber(prepared$started_monotonic)) {
                settled - prepared$started_monotonic
            } else {
                NA_real_
            }
        }
        if (!retention_acknowledged && !is.null(sample) &&
            is.finite(marker_elapsed) && sample$elapsed_seconds >=
                marker_elapsed + execution$retained_dwell_seconds &&
            file.exists(prepared$paths$acknowledgement)) {
            acknowledged <- tryCatch(
                publicationAcknowledgeRetention(
                    prepared$paths$acknowledgement
                ),
                error = \(error) FALSE
            )
            if (isTRUE(acknowledged)) {
                retention_acknowledged <- TRUE
            } else {
                safety_aborted <- TRUE
                safety_reason <- "retention_acknowledgement_failed"
                publicationStopSystemdUnit(unit_name)
                process$kill()
            }
        }
        if (elapsed > execution$timeout_seconds && process$is_alive()) {
            timed_out <- TRUE
            publicationStopSystemdUnit(unit_name)
            process$kill()
        }
        if (!process$is_alive()) break
        Sys.sleep(execution$sampling_interval_ms / 1000)
    }
    list(
        samples = samples,
        marker_elapsed = marker_elapsed,
        retention_state = retention_state,
        cgroup_observed = cgroup_observed,
        cgroup_lost = cgroup_lost,
        expected_root_pid = expected_root_pid,
        pid_ownership_valid = pid_ownership_valid &&
            is.finite(expected_root_pid),
        unexpected_pids = unexpected_pids,
        indeterminate_pids = indeterminate_pids,
        retention_acknowledged = retention_acknowledged,
        timed_out = timed_out,
        safety_aborted = safety_aborted,
        safety_reason = safety_reason,
        safety_trace = safety_trace,
        safety_trace_required = safety_trace_required,
        safety_interval_seconds = safety_interval_seconds
    )
}

publicationFinalizeCgroupMeasurement <- function(
    prepared,
    sampled,
    run_dir,
    execution,
    unit_name,
    disk_categories
) {
    process <- prepared$process
    process$wait(timeout = 5000)
    total_elapsed <- proc.time()[["elapsed"]] - prepared$started
    final_disk <- publicationDirectoryMetrics(run_dir, disk_categories)
    resource_path <- file.path(run_dir, "worker-resources.json")
    worker_resources <- if (file.exists(resource_path)) {
        publicationReadRetentionState(resource_path)
    } else {
        NULL
    }
    unit_result <- publicationSystemdUnitResult(
        unit_name,
        process$get_exit_status()
    )
    cleanup <- publicationOwnedUnitCleanup(unit_name, prepared$cgroup_path)
    result <- publicationCgroupResult(
        process = process,
        samples = sampled$samples,
        marker_elapsed = sampled$marker_elapsed,
        retention_state = sampled$retention_state,
        execution = execution,
        timed_out = sampled$timed_out,
        cgroup_observed = sampled$cgroup_observed,
        cgroup_lost = sampled$cgroup_lost,
        expected_root_pid = sampled$expected_root_pid,
        pid_ownership_valid = sampled$pid_ownership_valid,
        unexpected_pids = sampled$unexpected_pids,
        indeterminate_pids = sampled$indeterminate_pids,
        retention_acknowledged = sampled$retention_acknowledged,
        unit_result = unit_result,
        worker_resource_ledger = worker_resources,
        cleanup = cleanup,
        unit_name = unit_name,
        cgroup_path = prepared$cgroup_path,
        stdout_path = prepared$paths$stdout,
        stderr_path = prepared$paths$stderr,
        started_at = prepared$started_at,
        total_elapsed = total_elapsed,
        final_disk = final_disk,
        safety_aborted = sampled$safety_aborted,
        safety_reason = sampled$safety_reason
    )
    result$safety_evidence <- if (isTRUE(sampled$safety_trace_required)) {
        publicationSafetyTraceEvidence(
            sampled$safety_trace,
            execution$timeout_seconds,
            sampled$safety_interval_seconds
        )
    } else {
        NULL
    }
    if (isTRUE(sampled$safety_trace_required) &&
        !isTRUE(tryCatch({
            publicationValidateSafetyTraceEvidence(result$safety_evidence)
            TRUE
        }, error = \(...) FALSE))) {
        result$status <- "failed"
        result$publication_certifiable <- FALSE
        result$safety_aborted <- TRUE
        result$safety_reason <- "safety_trace_invalid"
    }
    result
}

publicationMeasureCgroupSubprocess <- function(
    command,
    command_args,
    run_dir,
    execution = list(
        sampling_interval_ms = 20,
        timeout_seconds = 60,
        retained_dwell_seconds = 5,
        retained_window_seconds = 2,
        maximum_boundary_bracket_seconds = 0.5
    ),
    env = character(),
    unit_name = publicationSystemdUnitName(),
    disk_categories = publicationDefaultDiskCategories(),
    require_certified_host = FALSE,
    host_preflight = NULL,
    safety_check_fn = NULL
) {
    prepared <- publicationPrepareCgroupMeasurement(
        command,
        command_args,
        run_dir,
        execution,
        env,
        unit_name,
        require_certified_host,
        host_preflight
    )
    process <- prepared$process
    on.exit({
        if (process$is_alive()) process$kill()
        publicationStopSystemdUnit(unit_name)
    }, add = TRUE)
    sampled <- publicationRunCgroupSampling(
        prepared,
        run_dir,
        execution,
        unit_name,
        disk_categories,
        safety_check_fn
    )
    publicationFinalizeCgroupMeasurement(
        prepared,
        sampled,
        run_dir,
        execution,
        unit_name,
        disk_categories
    )
}
