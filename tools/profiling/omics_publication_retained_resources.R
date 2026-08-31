publicationReadRetentionState <- function(path) {
    if (!file.exists(path)) return(NULL)
    tryCatch(jsonlite::read_json(path, simplifyVector = TRUE), error = \(...) NULL)
}

publicationWorkerResourceNames <- function() {
    c(
        "arrow_pool_bytes", "duckdb_memory_bytes", "duckdb_spill_bytes",
        "duckdb_connections", "duckdb_results", "duckdb_prepared_statements",
        "temporary_paths", "cache_entries", "active_tasks", "observers",
        "native_resources"
    )
}

publicationEmptyWorkerResources <- function() {
    stats::setNames(
        as.list(rep(0, length(publicationWorkerResourceNames()))),
        publicationWorkerResourceNames()
    )
}

publicationWorkerResourceLedgerValid <- function(ledger) {
    if (!is.list(ledger) || !setequal(names(ledger), c(
        "schema", "schema_version", "high_water", "retained", "terminal"
    )) || !identical(
        ledger$schema,
        "multischolar.omics_publication_worker_resources"
    ) || !identical(ledger$schema_version, "1.0.0")) {
        return(FALSE)
    }
    expected <- publicationWorkerResourceNames()
    states <- c("high_water", "retained", "terminal")
    if (any(!vapply(states, \(state) {
        is.list(ledger[[state]]) && setequal(names(ledger[[state]]), expected) &&
            all(vapply(ledger[[state]], \(value) {
                publicationScalarNumber(value) && value >= 0
            }, logical(1)))
    }, logical(1)))) {
        return(FALSE)
    }
    high_water <- unlist(ledger$high_water, use.names = TRUE)
    retained <- unlist(ledger$retained, use.names = TRUE)
    terminal <- unlist(ledger$terminal, use.names = TRUE)
    all(high_water >= retained) && all(retained == 0) && all(terminal == 0)
}

publicationRetentionStateValid <- function(state) {
    required <- c(
        "active_tasks", "open_queries", "open_writers", "open_leases",
        "open_connections", "open_results", "active_child_processes",
        "arrow_pool_bytes", "duckdb_memory_bytes", "duckdb_spill_bytes",
        "duckdb_prepared_statements", "temporary_paths", "cache_entries",
        "observers", "native_resources",
        "retained_dwell_seconds", "retention_acknowledgement",
        "settled_monotonic_seconds"
    )
    resources <- setdiff(required, c(
        "settled_monotonic_seconds", "retained_dwell_seconds",
        "retention_acknowledgement"
    ))
    is.list(state) && all(required %in% names(state)) &&
        publicationScalarNumber(state$settled_monotonic_seconds, positive = TRUE) &&
        identical(as.numeric(state$retained_dwell_seconds), 5) &&
        identical(state$retention_acknowledgement, "fifo_v1") &&
        all(vapply(state[resources], \(value) {
            is.numeric(value) && length(value) == 1L && is.finite(value) && value == 0
        }, logical(1)))
}

publicationSampleCgroup <- function(
    cgroup_path,
    elapsed_seconds,
    disk
) {
    snapshot <- publicationCgroupSnapshot(cgroup_path)
    if (is.null(snapshot)) return(NULL)
    snapshot$elapsed_seconds <- elapsed_seconds
    snapshot$disk <- disk
    snapshot
}

publicationMetricMaximum <- function(samples, accessor) {
    values <- vapply(samples, \(sample) {
        value <- accessor(sample)
        if (length(value) != 1L || !is.finite(value)) NA_real_ else value
    }, numeric(1))
    if (!any(is.finite(values))) return(NA_real_)
    max(values, na.rm = TRUE)
}

publicationRetainedSamples <- function(
    samples,
    marker_elapsed,
    dwell_seconds,
    retained_window_seconds
) {
    lower <- marker_elapsed + dwell_seconds - retained_window_seconds
    upper <- marker_elapsed + dwell_seconds
    Filter(
        \(sample) sample$elapsed_seconds >= lower &&
            sample$elapsed_seconds <= upper,
        samples
    )
}

publicationRetainedMedian <- function(samples, accessor) {
    values <- vapply(samples, \(sample) {
        value <- accessor(sample)
        if (length(value) != 1L || !is.finite(value)) NA_real_ else value
    }, numeric(1))
    values <- values[is.finite(values)]
    if (!length(values)) return(NA_real_)
    stats::median(values)
}

publicationMetricDelta <- function(samples, accessor) {
    if (length(samples) < 2L) return(NA_real_)
    values <- vapply(samples, \(sample) {
        value <- accessor(sample)
        if (length(value) != 1L || !is.finite(value)) NA_real_ else value
    }, numeric(1))
    if (any(!is.finite(values))) return(NA_real_)
    values[[length(values)]] - values[[1L]]
}

publicationMetricSlope <- function(samples, accessor) {
    if (length(samples) < 2L) return(NA_real_)
    elapsed <- vapply(samples, `[[`, numeric(1), "elapsed_seconds")
    values <- vapply(samples, \(sample) {
        value <- accessor(sample)
        if (length(value) != 1L || !is.finite(value)) NA_real_ else value
    }, numeric(1))
    centered <- elapsed - mean(elapsed)
    denominator <- sum(centered^2)
    if (any(!is.finite(values)) || !is.finite(denominator) || denominator <= 0) {
        return(NA_real_)
    }
    sum(centered * (values - mean(values))) / denominator
}

publicationRetainedTrace <- function(samples) {
    lapply(samples, \(sample) {
        list(
            elapsed_seconds = sample$elapsed_seconds,
            charged_memory_bytes = sample$memory$current_bytes,
            anonymous_memory_bytes = sample$memory$stat$anon,
            file_memory_bytes = sample$memory$stat$file,
            kernel_memory_bytes = sample$memory$stat$kernel,
            slab_memory_bytes = sample$memory$stat$slab,
            swap_bytes = sample$memory$swap_current_bytes,
            pss_bytes = sample$smaps$pss_bytes,
            uss_bytes = sample$smaps$uss_bytes,
            rss_bytes = sample$smaps$rss_bytes,
            cpu_usage_usec = sample$cpu$usage_usec,
            io_read_bytes = sample$io$rbytes,
            io_write_bytes = sample$io$wbytes,
            process_count = sample$process_count,
            file_descriptors = sample$resources$file_descriptors,
            thread_tasks = sample$resources$thread_tasks,
            kernel_locks = sample$resources$kernel_locks
        )
    })
}

publicationRetainedDiagnostics <- function(samples) {
    trace <- publicationRetainedTrace(samples)
    elapsed <- vapply(samples, `[[`, numeric(1), "elapsed_seconds")
    cpu_activity <- publicationMetricDelta(
        samples,
        \(sample) sample$cpu$usage_usec
    ) / 1e6
    io_read_activity <- publicationMetricDelta(
        samples,
        \(sample) sample$io$rbytes
    )
    io_write_activity <- publicationMetricDelta(
        samples,
        \(sample) sample$io$wbytes
    )
    list(
        sample_count = length(samples),
        observed_window_seconds = if (length(elapsed) > 1L) {
            max(elapsed) - min(elapsed)
        } else {
            NA_real_
        },
        charged_memory_slope_bytes_per_second = publicationMetricSlope(
            samples,
            \(sample) sample$memory$current_bytes
        ),
        anonymous_memory_slope_bytes_per_second = publicationMetricSlope(
            samples,
            \(sample) sample$memory$stat$anon
        ),
        cpu_activity_seconds = cpu_activity,
        io_read_activity_bytes = io_read_activity,
        io_write_activity_bytes = io_write_activity,
        background_activity_observed = {
            activity <- c(
                cpu_activity,
                io_read_activity,
                io_write_activity
            )
            any(!is.finite(activity)) || any(activity != 0)
        },
        maximum_swap_bytes = publicationMetricMaximum(
            samples,
            \(sample) sample$memory$swap_current_bytes
        ),
        maximum_process_count = publicationMetricMaximum(
            samples,
            \(sample) sample$process_count
        ),
        trace_sha256 = publicationObjectDigest(trace),
        trace = trace
    )
}

publicationRetainedWindowEvidence <- function(
    retained_samples,
    all_samples,
    marker_elapsed,
    dwell_seconds,
    retained_window_seconds,
    boundary_tolerance_seconds,
    maximum_boundary_bracket_seconds
) {
    expected_start <- marker_elapsed + dwell_seconds - retained_window_seconds
    expected_stop <- marker_elapsed + dwell_seconds
    elapsed <- vapply(retained_samples, `[[`, numeric(1), "elapsed_seconds")
    all_elapsed <- vapply(all_samples, `[[`, numeric(1), "elapsed_seconds")
    pre_boundary <- all_elapsed[all_elapsed <= expected_start]
    post_boundary <- all_elapsed[all_elapsed >= expected_stop]
    lower_bracket_seconds <- if (length(elapsed) && length(pre_boundary)) {
        min(elapsed) - max(pre_boundary)
    } else {
        NA_real_
    }
    upper_bracket_seconds <- if (length(elapsed) && length(post_boundary)) {
        min(post_boundary) - max(elapsed)
    } else {
        NA_real_
    }
    lower_bracket_valid <- is.finite(lower_bracket_seconds) &&
        lower_bracket_seconds <= maximum_boundary_bracket_seconds
    upper_bracket_valid <- is.finite(upper_bracket_seconds) &&
        upper_bracket_seconds <= maximum_boundary_bracket_seconds
    valid <- is.finite(marker_elapsed) &&
        publicationScalarNumber(boundary_tolerance_seconds, positive = TRUE) &&
        publicationScalarNumber(
            maximum_boundary_bracket_seconds,
            positive = TRUE
        ) &&
        length(elapsed) >= 2L && all(is.finite(elapsed)) &&
        (min(elapsed) <= expected_start + boundary_tolerance_seconds ||
            lower_bracket_valid) &&
        (max(elapsed) >= expected_stop - boundary_tolerance_seconds ||
            upper_bracket_valid) &&
        all(elapsed >= expected_start) && all(elapsed <= expected_stop)
    list(
        valid = valid,
        expected_start_seconds = expected_start,
        expected_stop_seconds = expected_stop,
        observed_start_seconds = if (length(elapsed)) min(elapsed) else NA_real_,
        observed_stop_seconds = if (length(elapsed)) max(elapsed) else NA_real_,
        pre_window_sample_seconds = if (length(pre_boundary)) {
            max(pre_boundary)
        } else {
            NA_real_
        },
        post_boundary_sample_seconds = if (length(post_boundary)) {
            min(post_boundary)
        } else {
            NA_real_
        },
        lower_boundary_bracketed = lower_bracket_valid,
        lower_boundary_bracket_seconds = lower_bracket_seconds,
        upper_boundary_bracketed = upper_bracket_valid,
        upper_boundary_bracket_seconds = upper_bracket_seconds,
        boundary_tolerance_seconds = boundary_tolerance_seconds,
        maximum_boundary_bracket_seconds = maximum_boundary_bracket_seconds
    )
}

publicationCgroupResult <- function(
    process,
    samples,
    marker_elapsed,
    retention_state,
    execution,
    timed_out,
    cgroup_observed,
    cgroup_lost,
    expected_root_pid,
    pid_ownership_valid,
    unexpected_pids,
    indeterminate_pids,
    retention_acknowledged,
    unit_result,
    worker_resource_ledger,
    cleanup,
    unit_name,
    cgroup_path,
    stdout_path,
    stderr_path,
    started_at,
    total_elapsed,
    final_disk,
    safety_aborted,
    safety_reason
) {
    retained <- publicationRetainedSamples(
        samples,
        marker_elapsed,
        execution$retained_dwell_seconds,
        execution$retained_window_seconds
    )
    boundary_tolerance_seconds <- if (is.null(
        execution$retained_boundary_tolerance_ms
    )) {
        0.1
    } else {
        execution$retained_boundary_tolerance_ms / 1000
    }
    retained_window <- publicationRetainedWindowEvidence(
        retained,
        samples,
        marker_elapsed,
        execution$retained_dwell_seconds,
        execution$retained_window_seconds,
        boundary_tolerance_seconds,
        execution$maximum_boundary_bracket_seconds
    )
    retained_diagnostics <- publicationRetainedDiagnostics(retained)
    swap_peak <- publicationMetricMaximum(
        samples,
        \(sample) sample$memory$swap_current_bytes
    )
    oom_count <- publicationMetricMaximum(
        samples,
        \(sample) as.numeric(sample$events$memory$oom)
    )
    oom_kill_count <- publicationMetricMaximum(
        samples,
        \(sample) as.numeric(sample$events$memory$oom_kill)
    )
    retained_valid <- is.finite(marker_elapsed) &&
        publicationRetentionStateValid(retention_state) &&
        isTRUE(retained_window$valid) &&
        !isTRUE(retained_diagnostics$background_activity_observed) &&
        identical(retained_diagnostics$maximum_process_count, 1) &&
        is.finite(swap_peak) && swap_peak == 0
    metrics <- publicationCgroupMetrics(
        samples,
        retained,
        total_elapsed,
        final_disk
    )
    metrics_valid <- publicationCgroupMetricsValid(metrics)
    terminal_valid <- !timed_out && !safety_aborted && cgroup_observed &&
        !cgroup_lost && isTRUE(pid_ownership_valid) &&
        !length(unexpected_pids) && isTRUE(unit_result$available) &&
        isTRUE(retention_acknowledged) &&
        identical(unit_result$result, "success") &&
        !isTRUE(unit_result$oom_killed) &&
        publicationWorkerResourceLedgerValid(worker_resource_ledger) &&
        isTRUE(cleanup$valid) &&
        retained_valid && metrics_valid && is.finite(oom_count) &&
        oom_count == 0 && is.finite(oom_kill_count) && oom_kill_count == 0 &&
        identical(process$get_exit_status(), 0L)
    list(
        schema = "multischolar.omics_publication_cgroup_measurement",
        schema_version = "1.0.0",
        metric_contract = publicationMetricContractBinding(),
        status = if (terminal_valid) "passed" else "failed",
        unit_name = unit_name,
        cgroup_path = cgroup_path,
        cgroup_observed = cgroup_observed,
        cgroup_lost = cgroup_lost,
        expected_root_pid = expected_root_pid,
        pid_ownership_valid = pid_ownership_valid,
        unexpected_pids = as.list(unexpected_pids),
        indeterminate_pids = as.list(indeterminate_pids),
        retention_acknowledged = retention_acknowledged,
        unit_result = unit_result,
        worker_resource_ledger = worker_resource_ledger,
        worker_resource_ledger_valid = publicationWorkerResourceLedgerValid(
            worker_resource_ledger
        ),
        cleanup = cleanup,
        publication_certifiable = terminal_valid,
        started_at = started_at,
        finished_at = publicationUtcNow(),
        exit_status = process$get_exit_status(),
        timed_out = timed_out,
        safety_aborted = safety_aborted,
        safety_reason = safety_reason,
        retention_marker_elapsed_seconds = marker_elapsed,
        retention_state = retention_state,
        retention_state_valid = publicationRetentionStateValid(retention_state),
        retained_window = retained_window,
        retained_diagnostics = retained_diagnostics,
        metrics_valid = metrics_valid,
        retained_sample_count = length(retained),
        metrics = metrics,
        stdout_path = stdout_path,
        stderr_path = stderr_path,
        samples = samples
    )
}

publicationCgroupMetricsValid <- function(metrics) {
    positive <- c(
        "elapsed_seconds", "peak_charged_memory_bytes",
        "retained_charged_memory_bytes", "peak_anonymous_memory_bytes",
        "retained_anonymous_memory_bytes", "peak_pss_bytes",
        "retained_pss_bytes", "peak_uss_bytes", "retained_uss_bytes",
        "peak_rss_bytes", "retained_rss_bytes", "cpu_usage_seconds",
        "maximum_process_count", "peak_file_descriptors",
        "retained_file_descriptors", "peak_thread_tasks",
        "retained_thread_tasks", "peak_logical_disk_bytes",
        "peak_allocated_disk_bytes", "final_logical_disk_bytes",
        "final_allocated_disk_bytes", "final_file_count", "sample_count"
    )
    nonnegative <- c(
        "peak_file_memory_bytes", "retained_file_memory_bytes",
        "peak_kernel_memory_bytes", "retained_kernel_memory_bytes",
        "peak_slab_memory_bytes", "retained_slab_memory_bytes",
        "peak_sampled_charged_memory_bytes",
        "peak_sampled_accounted_memory_bytes",
        "retained_accounted_memory_bytes",
        "peak_reconciliation_absolute_residual_bytes",
        "retained_reconciliation_absolute_residual_bytes", "cpu_user_seconds",
        "cpu_system_seconds", "peak_kernel_locks", "retained_kernel_locks"
    )
    if (!is.list(metrics) || !all(c(positive, nonnegative, "peak_swap_bytes") %in%
        names(metrics)) || !is.list(metrics$io) ||
        !identical(sort(names(metrics$io)), sort(c(
            "rbytes", "wbytes", "rios", "wios", "dbytes", "dios"
        ))) || !is.list(metrics$memory_events) ||
        !all(c("oom", "oom_kill") %in% names(metrics$memory_events))) {
        return(FALSE)
    }
    finite_positive <- all(vapply(metrics[positive], \(value) {
        publicationScalarNumber(value, positive = TRUE)
    }, logical(1)))
    finite_nonnegative <- all(vapply(metrics[nonnegative], \(value) {
        publicationScalarNumber(value) && value >= 0
    }, logical(1)))
    finite_positive && publicationScalarNumber(metrics$peak_swap_bytes) &&
        metrics$peak_swap_bytes == 0 && finite_nonnegative &&
        all(vapply(metrics$io, publicationScalarNumber, logical(1))) &&
        all(vapply(metrics$memory_events, \(value) {
            publicationScalarNumber(value) && value >= 0
        }, logical(1))) && metrics$memory_events$oom == 0 &&
        metrics$memory_events$oom_kill == 0
}

publicationMemoryStatValue <- function(sample, name) {
    as.numeric(sample$memory$stat[[name]])
}

publicationAccountedMemoryBytes <- function(sample) {
    sum(as.numeric(unlist(
        sample$memory$stat[c("anon", "file", "kernel")],
        use.names = FALSE
    )))
}

publicationMemoryReconciliationResidual <- function(sample) {
    abs(sample$memory$current_bytes - publicationAccountedMemoryBytes(sample))
}

publicationCgroupChargedMemoryMetrics <- function(samples, retained) {
    statMax <- \(name) publicationMetricMaximum(
        samples,
        \(sample) publicationMemoryStatValue(sample, name)
    )
    statMedian <- \(name) publicationRetainedMedian(
        retained,
        \(sample) publicationMemoryStatValue(sample, name)
    )
    list(
        peak_charged_memory_bytes = publicationMetricMaximum(
            samples,
            \(sample) sample$memory$peak_bytes
        ),
        peak_sampled_charged_memory_bytes = publicationMetricMaximum(
            samples,
            \(sample) sample$memory$current_bytes
        ),
        retained_charged_memory_bytes = publicationRetainedMedian(
            retained,
            \(sample) sample$memory$current_bytes
        ),
        peak_anonymous_memory_bytes = statMax("anon"),
        retained_anonymous_memory_bytes = statMedian("anon"),
        peak_file_memory_bytes = statMax("file"),
        retained_file_memory_bytes = statMedian("file"),
        peak_kernel_memory_bytes = statMax("kernel"),
        retained_kernel_memory_bytes = statMedian("kernel"),
        peak_slab_memory_bytes = statMax("slab"),
        retained_slab_memory_bytes = statMedian("slab"),
        peak_sampled_accounted_memory_bytes = publicationMetricMaximum(
            samples,
            publicationAccountedMemoryBytes
        ),
        retained_accounted_memory_bytes = publicationRetainedMedian(
            retained,
            publicationAccountedMemoryBytes
        ),
        peak_reconciliation_absolute_residual_bytes = publicationMetricMaximum(
            samples,
            publicationMemoryReconciliationResidual
        ),
        retained_reconciliation_absolute_residual_bytes =
            publicationRetainedMedian(
                retained,
                publicationMemoryReconciliationResidual
            ),
        peak_swap_bytes = publicationMetricMaximum(
            samples,
            \(sample) sample$memory$swap_current_bytes
        )
    )
}

publicationCgroupProcessMemoryMetrics <- function(samples, retained) {
    list(
        peak_pss_bytes = publicationMetricMaximum(
            samples,
            \(sample) sample$smaps$pss_bytes
        ),
        retained_pss_bytes = publicationRetainedMedian(
            retained,
            \(sample) sample$smaps$pss_bytes
        ),
        peak_uss_bytes = publicationMetricMaximum(
            samples,
            \(sample) sample$smaps$uss_bytes
        ),
        retained_uss_bytes = publicationRetainedMedian(
            retained,
            \(sample) sample$smaps$uss_bytes
        ),
        peak_rss_bytes = publicationMetricMaximum(
            samples,
            \(sample) sample$smaps$rss_bytes
        ),
        retained_rss_bytes = publicationRetainedMedian(
            retained,
            \(sample) sample$smaps$rss_bytes
        )
    )
}

publicationCgroupCpuIoMetrics <- function(samples) {
    cpuMax <- \(name) publicationMetricMaximum(
        samples,
        \(sample) as.numeric(sample$cpu[[name]])
    ) / 1e6
    ioMax <- \(name) publicationMetricMaximum(
        samples,
        \(sample) as.numeric(sample$io[[name]])
    )
    pressureMax <- \(resource, level) publicationMetricMaximum(
        samples,
        \(sample) as.numeric(sample$pressure[[resource]][[level]]$total)
    )
    eventMax <- \(name) publicationMetricMaximum(
        samples,
        \(sample) as.numeric(sample$events$memory[[name]])
    )
    io_names <- c("rbytes", "wbytes", "rios", "wios", "dbytes", "dios")
    event_names <- c("low", "high", "max", "oom", "oom_kill", "oom_group_kill")
    list(
        cpu_usage_seconds = cpuMax("usage_usec"),
        cpu_user_seconds = cpuMax("user_usec"),
        cpu_system_seconds = cpuMax("system_usec"),
        io = stats::setNames(lapply(io_names, ioMax), io_names),
        pressure_total_usec = list(
            memory_some = pressureMax("memory", "some"),
            memory_full = pressureMax("memory", "full"),
            cpu_some = pressureMax("cpu", "some"),
            cpu_full = pressureMax("cpu", "full"),
            io_some = pressureMax("io", "some"),
            io_full = pressureMax("io", "full")
        ),
        memory_events = stats::setNames(
            lapply(event_names, eventMax),
            event_names
        )
    )
}

publicationCgroupLifecycleMetrics <- function(samples, retained) {
    list(
        maximum_process_count = publicationMetricMaximum(
            samples,
            \(sample) sample$process_count
        ),
        peak_file_descriptors = publicationMetricMaximum(
            samples,
            \(sample) sample$resources$file_descriptors
        ),
        retained_file_descriptors = publicationRetainedMedian(
            retained,
            \(sample) sample$resources$file_descriptors
        ),
        peak_thread_tasks = publicationMetricMaximum(
            samples,
            \(sample) sample$resources$thread_tasks
        ),
        retained_thread_tasks = publicationRetainedMedian(
            retained,
            \(sample) sample$resources$thread_tasks
        ),
        peak_kernel_locks = publicationMetricMaximum(
            samples,
            \(sample) sample$resources$kernel_locks
        ),
        retained_kernel_locks = publicationRetainedMedian(
            retained,
            \(sample) sample$resources$kernel_locks
        )
    )
}

publicationCgroupDiskMetrics <- function(samples, final_disk) {
    category_names <- names(samples[[1L]]$disk$category_logical_bytes)
    categoryMax <- \(name) publicationMetricMaximum(
        samples,
        \(sample) as.numeric(sample$disk$category_logical_bytes[[name]])
    )
    list(
        peak_logical_disk_bytes = publicationMetricMaximum(
            samples,
            \(sample) sample$disk$logical_bytes
        ),
        peak_allocated_disk_bytes = publicationMetricMaximum(
            samples,
            \(sample) sample$disk$allocated_bytes
        ),
        peak_disk_category_logical_bytes = stats::setNames(
            lapply(category_names, categoryMax),
            category_names
        ),
        final_logical_disk_bytes = as.numeric(final_disk$logical_bytes),
        final_allocated_disk_bytes = as.numeric(final_disk$allocated_bytes),
        final_file_count = as.numeric(final_disk$file_count),
        sample_count = length(samples)
    )
}

publicationCgroupMetrics <- function(
    samples,
    retained,
    elapsed_seconds,
    final_disk
) {
    if (!length(samples)) return(list())
    c(
        list(elapsed_seconds = elapsed_seconds),
        publicationCgroupChargedMemoryMetrics(samples, retained),
        publicationCgroupProcessMemoryMetrics(samples, retained),
        publicationCgroupCpuIoMetrics(samples),
        publicationCgroupLifecycleMetrics(samples, retained),
        publicationCgroupDiskMetrics(samples, final_disk)
    )
}
