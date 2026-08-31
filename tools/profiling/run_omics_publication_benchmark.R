#!/usr/bin/env Rscript

publicationRunnerScriptPath <- function() {
    file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
    if (length(file_arg)) {
        return(normalizePath(sub("^--file=", "", file_arg[[1L]]), mustWork = TRUE))
    }
    normalizePath(
        file.path(getwd(), "tools", "profiling", "run_omics_publication_benchmark.R"),
        mustWork = TRUE
    )
}

.PUBLICATION_RUNNER_REPO_ROOT <- normalizePath(
    file.path(dirname(publicationRunnerScriptPath()), "..", ".."),
    mustWork = TRUE
)

for (source_file in c(
    "omics_publication_protocol.R",
    "omics_publication_measure_spec.R",
    "omics_publication_linux_resources.R",
    "omics_publication_retained_resources.R",
    "omics_publication_host_safety.R",
    "omics_publication_schedule.R",
    "omics_publication_campaign_state.R",
    "omics_publication_statistics.R"
)) {
    source(
        file.path(
            .PUBLICATION_RUNNER_REPO_ROOT,
            "tools",
            "profiling",
            source_file
        ),
        local = FALSE
    )
}

publicationRunnerDefaultArgs <- function() {
    list(
        worker = NULL,
        run_dir = NULL,
        dwell_seconds = 5,
        null_calibration = FALSE,
        runs = 30L,
        output = NULL,
        measure_spec = NULL,
        host_envelope = FALSE,
        help = FALSE
    )
}

publicationRunnerParseBool <- function(value) {
    tolower(as.character(value)) %in% c("1", "true", "yes", "y")
}

publicationRunnerParseArgs <- function(argv) {
    args <- publicationRunnerDefaultArgs()
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
        parts <- strsplit(sub("^--", "", token), "=", fixed = TRUE)[[1L]]
        key <- gsub("-", "_", parts[[1L]], fixed = TRUE)
        if (!key %in% names(args)) stop(paste("Unknown option:", token), call. = FALSE)
        if (length(parts) == 2L) {
            value <- parts[[2L]]
        } else {
            index <- index + 1L
            if (index > length(argv)) stop(paste("Missing value for:", token), call. = FALSE)
            value <- argv[[index]]
        }
        if (key %in% c("null_calibration", "host_envelope")) {
            args[[key]] <- publicationRunnerParseBool(value)
        } else if (identical(key, "runs")) {
            args[[key]] <- as.integer(value)
        } else if (identical(key, "dwell_seconds")) {
            args[[key]] <- as.numeric(value)
        } else {
            args[[key]] <- value
        }
        index <- index + 1L
    }
    args
}

publicationRunnerUsage <- function() {
    cat(paste(
        "Usage: Rscript --vanilla tools/profiling/run_omics_publication_benchmark.R [options]",
        "",
        "  --host-envelope true --output <json>",
        "  --null-calibration true --runs 30 --output <json>",
        "  --measure-spec <json> --output <json>",
        "  --worker null --run-dir <path> --dwell-seconds 5",
        sep = "\n"
    ))
}

publicationPhaseCounterSnapshot <- function() {
    cgroup_path <- publicationCurrentCgroupPath()
    if (is.null(cgroup_path) || !dir.exists(cgroup_path)) {
        stop("Worker phase cgroup is unavailable", call. = FALSE)
    }
    list(
        cpu = publicationParseNamedNumbers(publicationReadLinesSafe(file.path(
            cgroup_path,
            "cpu.stat"
        ))),
        io = publicationParseIoStat(publicationReadLinesSafe(file.path(
            cgroup_path,
            "io.stat"
        )))
    )
}

publicationWritePhaseBoundary <- function(run_dir, marker) {
    evidence <- list(
        marker = marker,
        pid = Sys.getpid(),
        counters = publicationPhaseCounterSnapshot()
    )
    publicationWriteJson(evidence, file.path(run_dir, paste0(marker, ".json")))
    evidence
}

publicationNullWorker <- function(run_dir, dwell_seconds) {
    if (!publicationScalarNumber(dwell_seconds, positive = TRUE)) {
        stop("Null worker dwell must be positive", call. = FALSE)
    }
    dir.create(run_dir, recursive = TRUE, showWarnings = FALSE)
    acknowledgement_path <- file.path(run_dir, "retention-sampled.fifo")
    fifo <- processx::run(
        "mkfifo",
        acknowledgement_path,
        error_on_status = FALSE,
        echo = FALSE
    )
    if (fifo$status != 0L) {
        stop("Could not create retention acknowledgement FIFO", call. = FALSE)
    }
    started <- proc.time()[["elapsed"]]
    cpu_started <- sum(proc.time()[c("user.self", "sys.self")])
    publicationWritePhaseBoundary(run_dir, "measured_worker_start")
    payload <- numeric(1024L * 1024L)
    payload[] <- seq_along(payload) %% 17L
    checksum <- sum(payload)
    phase_elapsed <- proc.time()[["elapsed"]] - started
    phase_cpu <- sum(proc.time()[c("user.self", "sys.self")]) - cpu_started
    publicationWritePhaseBoundary(run_dir, "measured_worker_stop")
    publicationWriteJson(
        list(
            schema = "multischolar.omics_publication_null_worker",
            schema_version = "1.0.0",
            status = "passed",
            pid = Sys.getpid(),
            phase_elapsed_seconds = phase_elapsed,
            phase_cpu_seconds = phase_cpu,
            work_unit_id = "null_operation",
            work_count = 1,
            checksum = checksum
        ),
        file.path(run_dir, "worker-result.json")
    )
    empty_resources <- publicationEmptyWorkerResources()
    publicationWriteJson(
        list(
            schema = "multischolar.omics_publication_worker_resources",
            schema_version = "1.0.0",
            high_water = empty_resources,
            retained = empty_resources,
            terminal = empty_resources
        ),
        file.path(run_dir, "worker-resources.json")
    )
    publicationWriteJson(
        list(
            active_tasks = 0L,
            open_queries = 0L,
            open_writers = 0L,
            open_leases = 0L,
            open_connections = 0L,
            open_results = 0L,
            active_child_processes = 0L,
            arrow_pool_bytes = 0L,
            duckdb_memory_bytes = 0L,
            duckdb_spill_bytes = 0L,
            duckdb_prepared_statements = 0L,
            temporary_paths = 0L,
            cache_entries = 0L,
            observers = 0L,
            native_resources = 0L,
            retained_dwell_seconds = dwell_seconds,
            retention_acknowledgement = "fifo_v1",
            settled_monotonic_seconds = publicationMonotonicSeconds()
        ),
        file.path(run_dir, "retention-state.json")
    )
    file.create(file.path(run_dir, "retention-ready"))
    acknowledgement <- fifo(
        acknowledgement_path,
        open = "rb",
        blocking = TRUE
    )
    on.exit(close(acknowledgement), add = TRUE)
    token <- readBin(acknowledgement, what = "raw", n = 1L)
    if (!identical(token, as.raw(1L))) {
        stop("Retention acknowledgement is invalid", call. = FALSE)
    }
    publicationWriteJson(
        list(marker = "worker_exit", pid = Sys.getpid()),
        file.path(run_dir, "worker-exit.json")
    )
    invisible(0L)
}

publicationFaultWorker <- function(mode, run_dir, dwell_seconds) {
    if (!identical(
        Sys.getenv("MULTISCHOLAR_PUBLICATION_FAULT_INJECTION"),
        "true"
    )) {
        stop("Publication fault worker is disabled", call. = FALSE)
    }
    if (identical(mode, "crash")) {
        stop("Injected worker crash", call. = FALSE)
    }
    if (identical(mode, "timeout")) {
        Sys.sleep(dwell_seconds + 30)
        return(invisible(0L))
    }
    if (identical(mode, "child_leak")) {
        child <- processx::process$new("sleep", "30", cleanup = FALSE)
        on.exit(if (child$is_alive()) child$kill(), add = TRUE)
        return(publicationNullWorker(run_dir, dwell_seconds))
    }
    if (identical(mode, "oom")) {
        payload <- vector("list", 128L)
        for (index in seq_along(payload)) {
            payload[[index]] <- raw(16L * 1024L * 1024L)
            payload[[index]][] <- as.raw(index %% 256L)
        }
        stop("Injected OOM limit did not terminate worker", call. = FALSE)
    }
    stop("Unknown publication fault mode", call. = FALSE)
}

publicationReadPhaseEvidence <- function(run_dir) {
    start <- publicationReadJson(file.path(
        run_dir,
        "measured_worker_start.json"
    ))
    stop <- publicationReadJson(file.path(
        run_dir,
        "measured_worker_stop.json"
    ))
    worker <- publicationReadJson(file.path(run_dir, "worker-result.json"))
    counter_names <- c("usage_usec", "user_usec", "system_usec")
    counters_valid <- all(counter_names %in% names(start$counters$cpu)) &&
        all(counter_names %in% names(stop$counters$cpu)) &&
        all(vapply(c(start$counters$cpu[counter_names],
            stop$counters$cpu[counter_names]), publicationScalarNumber, logical(1))) &&
        all(vapply(c(start$counters$io, stop$counters$io),
            publicationScalarNumber, logical(1)))
    valid <- identical(start$marker, "measured_worker_start") &&
        identical(stop$marker, "measured_worker_stop") &&
        identical(as.integer(start$pid), as.integer(worker$pid)) &&
        identical(as.integer(stop$pid), as.integer(worker$pid)) &&
        identical(worker$status, "passed") &&
        publicationScalarNumber(worker$phase_elapsed_seconds, positive = TRUE) &&
        publicationScalarNumber(worker$phase_cpu_seconds, positive = TRUE) &&
        publicationScalarNumber(worker$work_count, positive = TRUE) &&
        publicationScalarString(worker$work_unit_id) && counters_valid
    if (!valid) stop("Worker phase evidence is invalid", call. = FALSE)
    list(start = start, stop = stop, worker = worker)
}

publicationAttachPhaseEvidence <- function(measurement, run_dir) {
    phase <- publicationReadPhaseEvidence(run_dir)
    worker <- phase$worker
    cpu_delta <- stats::setNames(lapply(
        c("usage_usec", "user_usec", "system_usec"),
        \(name) as.numeric(phase$stop$counters$cpu[[name]]) -
            as.numeric(phase$start$counters$cpu[[name]])
    ), c("usage_usec", "user_usec", "system_usec"))
    io_delta <- stats::setNames(lapply(names(phase$start$counters$io), \(name) {
        as.numeric(phase$stop$counters$io[[name]]) -
            as.numeric(phase$start$counters$io[[name]])
    }), names(phase$start$counters$io))
    if (any(unlist(c(cpu_delta, io_delta), use.names = FALSE) < 0)) {
        stop("Worker phase counters decreased", call. = FALSE)
    }
    phase_cgroup_cpu_seconds <- cpu_delta$usage_usec / 1e6
    if (!publicationScalarNumber(phase_cgroup_cpu_seconds, positive = TRUE)) {
        stop("Worker phase cgroup CPU is invalid", call. = FALSE)
    }
    measurement$phase_evidence <- list(
        valid = TRUE,
        worker_pid = worker$pid,
        work_unit_id = worker$work_unit_id,
        work_count = worker$work_count,
        phase_elapsed_seconds = worker$phase_elapsed_seconds,
        phase_cgroup_cpu_seconds = phase_cgroup_cpu_seconds,
        phase_process_self_cpu_seconds = worker$phase_cpu_seconds,
        phase_cgroup_user_seconds = cpu_delta$user_usec / 1e6,
        phase_cgroup_system_seconds = cpu_delta$system_usec / 1e6,
        phase_cgroup_io = io_delta
    )
    measurement$metrics$cgroup_lifecycle_elapsed_seconds <-
        measurement$metrics$elapsed_seconds
    measurement$metrics$elapsed_seconds <- worker$phase_elapsed_seconds
    measurement$metrics$phase_cpu_seconds <- phase_cgroup_cpu_seconds
    measurement$metrics$phase_process_self_cpu_seconds <- worker$phase_cpu_seconds
    measurement$metrics$phase_io <- io_delta
    measurement$metrics$primary_work_units_per_wall_second <-
        worker$work_count / worker$phase_elapsed_seconds
    measurement$metrics$primary_work_units_per_cpu_second <-
        worker$work_count / phase_cgroup_cpu_seconds
    measurement
}

publicationNullExecution <- function(dwell_seconds) {
    list(
        sampling_interval_ms = 20,
        disk_sampling_interval_ms = 500,
        timeout_seconds = dwell_seconds + 30,
        retained_dwell_seconds = dwell_seconds,
        retained_window_seconds = min(2, dwell_seconds / 2),
        retention_acknowledgement = "fifo_v1",
        maximum_boundary_bracket_seconds = 0.5,
        retained_boundary_tolerance_ms = 100
    )
}

publicationNullThreadEnvironment <- function() {
    c(
        OMP_NUM_THREADS = "1",
        OPENBLAS_NUM_THREADS = "1",
        MKL_NUM_THREADS = "1",
        ARROW_NUM_THREADS = "1",
        DUCKDB_THREADS = "1",
        TZ = "UTC"
    )
}

publicationMeasureNullWorker <- function(run_dir, dwell_seconds) {
    publicationMeasureCgroupSubprocess(
        command = file.path(R.home("bin"), "Rscript"),
        command_args = c(
            "--vanilla",
            publicationRunnerScriptPath(),
            "--worker", "null",
            "--run-dir", run_dir,
            "--dwell-seconds", as.character(dwell_seconds)
        ),
        run_dir = run_dir,
        execution = publicationNullExecution(dwell_seconds),
        env = publicationNullThreadEnvironment(),
        unit_name = publicationSystemdUnitName("multischolar-null")
    )
}

publicationAttachNullPhase <- function(measurement, run_dir) {
    if (!identical(measurement$status, "passed")) {
        measurement$phase_evidence <- list(
            valid = FALSE,
            reason = "measurement_failed_before_valid_phase"
        )
        return(measurement)
    }
    tryCatch(
        publicationAttachPhaseEvidence(measurement, run_dir),
        error = \(error) {
            measurement$status <- "failed"
            measurement$publication_certifiable <- FALSE
            measurement$phase_evidence <- list(
                valid = FALSE,
                reason = conditionMessage(error)
            )
            measurement
        }
    )
}

publicationNullWorkerResult <- function(run_dir) {
    path <- file.path(run_dir, "worker-result.json")
    if (file.exists(path)) return(publicationReadJson(path))
    list(
        schema = "multischolar.omics_publication_null_worker",
        schema_version = "1.0.0",
        status = "failed",
        reason = "worker_result_missing"
    )
}

publicationCompactNullMeasurement <- function(
    measurement,
    worker,
    raw_path
) {
    diagnostics <- measurement$retained_diagnostics
    diagnostics$trace <- NULL
    list(
        run_id = measurement$run_id,
        status = measurement$status,
        metric_contract = measurement$metric_contract,
        exit_status = measurement$exit_status,
        timed_out = measurement$timed_out,
        safety_aborted = measurement$safety_aborted,
        safety_reason = measurement$safety_reason,
        cgroup_observed = measurement$cgroup_observed,
        cgroup_lost = measurement$cgroup_lost,
        expected_root_pid = measurement$expected_root_pid,
        pid_ownership_valid = measurement$pid_ownership_valid,
        unexpected_pids = measurement$unexpected_pids,
        indeterminate_pids = measurement$indeterminate_pids,
        retention_acknowledged = measurement$retention_acknowledged,
        unit_result = measurement$unit_result,
        worker_resource_ledger = measurement$worker_resource_ledger,
        worker_resource_ledger_valid = measurement$worker_resource_ledger_valid,
        cleanup = measurement$cleanup,
        publication_certifiable = measurement$publication_certifiable,
        retention_state = measurement$retention_state,
        retention_state_valid = measurement$retention_state_valid,
        retained_window = measurement$retained_window,
        retained_diagnostics = diagnostics,
        metrics_valid = measurement$metrics_valid,
        retained_sample_count = measurement$retained_sample_count,
        metrics = measurement$metrics,
        phase_evidence = measurement$phase_evidence,
        worker = worker,
        raw_measurement_sha256 = publicationFileDigest(raw_path)
    )
}

publicationNullCalibrationRun <- function(root, index, dwell_seconds) {
    run_dir <- file.path(root, sprintf("run-%03d", index))
    measurement <- publicationMeasureNullWorker(run_dir, dwell_seconds)
    measurement <- publicationAttachNullPhase(measurement, run_dir)
    measurement$run_id <- sprintf("null-%03d", index)
    measurement$metrics$harness_elapsed_seconds <- if (identical(
        measurement$status,
        "passed"
    )) {
        measurement$metrics$cgroup_lifecycle_elapsed_seconds - dwell_seconds -
            measurement$metrics$elapsed_seconds
    } else {
        NA_real_
    }
    raw_path <- file.path(run_dir, "measurement.json")
    publicationWriteJson(measurement, raw_path)
    publicationCompactNullMeasurement(
        measurement,
        publicationNullWorkerResult(run_dir),
        raw_path
    )
}

publicationNullCalibrationMethod <- function(dwell_seconds) {
    list(
        owned_cgroup_v2 = TRUE,
        sampling_interval_ms = 20,
        disk_sampling_interval_ms = 500,
        retained_dwell_seconds = dwell_seconds,
        retained_window_seconds = min(2, dwell_seconds / 2),
        retention_acknowledgement = "fifo_v1",
        maximum_boundary_bracket_seconds = 0.5,
        retained_boundary_tolerance_ms = 100,
        candidate_loaded = FALSE,
        worker_thread_environment = list(
            OMP_NUM_THREADS = 1L,
            OPENBLAS_NUM_THREADS = 1L,
            MKL_NUM_THREADS = 1L,
            ARROW_NUM_THREADS = 1L,
            DUCKDB_THREADS = 1L
        )
    )
}

publicationNullRecordsPassed <- function(records) {
    all(vapply(records, \(record) {
        identical(record$status, "passed") &&
            isTRUE(record$publication_certifiable) &&
            identical(record$exit_status, 0L) && !isTRUE(record$timed_out) &&
            !isTRUE(record$safety_aborted) && !isTRUE(record$cgroup_lost) &&
            isTRUE(record$pid_ownership_valid) && !length(record$unexpected_pids) &&
            isTRUE(record$retention_acknowledged) &&
            identical(record$unit_result$result, "success") &&
            !isTRUE(record$unit_result$oom_killed) &&
            isTRUE(record$worker_resource_ledger_valid) &&
            isTRUE(record$cleanup$valid)
    }, logical(1)))
}

publicationNullCalibration <- function(output, runs, dwell_seconds) {
    if (!is.numeric(runs) || length(runs) != 1L || runs < 2L) {
        stop("Null calibration runs must be at least two", call. = FALSE)
    }
    output <- normalizePath(output, mustWork = FALSE)
    root <- file.path(
        dirname(output),
        paste0(tools::file_path_sans_ext(basename(output)), "-runs")
    )
    if (file.exists(output) || dir.exists(root)) {
        stop("Null calibration output or raw root already exists", call. = FALSE)
    }
    dir.create(root, recursive = TRUE, showWarnings = FALSE)
    records <- lapply(
        seq_len(runs),
        \(index) publicationNullCalibrationRun(root, index, dwell_seconds)
    )
    host <- publicationHostEnvelope(dirname(output))
    host$cgroup$current_path <- "<redacted-owned-user-cgroup>"
    result <- list(
        schema = "multischolar.omics_publication_null_calibration",
        schema_version = "1.0.0",
        calibration_id = "multischolar.null_cgroup.2026-08-27.v6",
        generated_at = publicationUtcNow(),
        candidate_loaded = FALSE,
        implementation_sources = publicationKernelSourceBindings(),
        governance_bindings = publicationNullGovernanceBindings(),
        method = publicationNullCalibrationMethod(dwell_seconds),
        host = host,
        runs = records,
        status = if (publicationNullRecordsPassed(records)) "passed" else "failed"
    )
    publicationWriteJson(result, output)
    result
}

publicationMeasureSpec <- function(spec_path, output) {
    spec <- publicationReadJson(spec_path)
    publicationValidateMeasureSpec(spec)
    if (file.exists(output) || dir.exists(spec$run_dir)) {
        stop("Measure output or run directory already exists", call. = FALSE)
    }
    safety_monitor <- publicationRuntimeSafetyMonitor(
        spec$safety_limits,
        dirname(spec$run_dir)
    )
    result <- publicationMeasureCgroupSubprocess(
        command = spec$command,
        command_args = unlist(spec$arguments, use.names = FALSE),
        run_dir = spec$run_dir,
        execution = lapply(spec$execution, as.numeric),
        env = unlist(spec$environment, use.names = TRUE),
        unit_name = publicationSystemdUnitName("multischolar-workload"),
        require_certified_host = TRUE,
        host_preflight = spec$host_preflight,
        safety_check_fn = safety_monitor
    )
    if (identical(result$status, "passed")) {
        result <- tryCatch(
            publicationAttachPhaseEvidence(result, spec$run_dir),
            error = \(error) {
                result$status <- "failed"
                result$publication_certifiable <- FALSE
                result$phase_evidence <- list(
                    valid = FALSE,
                    reason = conditionMessage(error)
                )
                result
            }
        )
    }
    work_valid <- identical(result$status, "passed") &&
        identical(
            result$phase_evidence$work_unit_id,
            spec$work$primary_work_unit_id
        ) && result$phase_evidence$work_count ==
            spec$work$expected_primary_work_count
    if (!work_valid) {
        result$status <- "failed"
        result$publication_certifiable <- FALSE
    }
    result$measure_spec_binding <- list(
        measure_spec_id = spec$measure_spec_id,
        path = spec_path,
        sha256 = publicationFileDigest(spec_path)
    )
    result$protocol_binding <- spec$protocol_binding
    result$estimand_binding <- spec$estimand_binding
    result$schedule_binding <- spec$schedule_binding
    result$source_binding <- spec$source_binding
    result$candidate_binding <- spec$candidate_binding
    result$cache_evidence <- spec$cache_evidence
    result$work_binding <- spec$work
    result$work_binding_valid <- work_valid
    publicationWriteJson(result, output)
    result
}

publicationRunnerMain <- function(argv = commandArgs(trailingOnly = TRUE)) {
    args <- publicationRunnerParseArgs(argv)
    if (isTRUE(args$help)) {
        publicationRunnerUsage()
        return(invisible(0L))
    }
    if (!is.null(args$worker)) {
        if (is.null(args$run_dir)) {
            stop("Worker --run-dir is required", call. = FALSE)
        }
        if (identical(args$worker, "null")) {
            return(publicationNullWorker(args$run_dir, args$dwell_seconds))
        }
        return(publicationFaultWorker(
            args$worker,
            args$run_dir,
            args$dwell_seconds
        ))
    }
    if (isTRUE(args$host_envelope)) {
        if (is.null(args$output)) stop("--output is required", call. = FALSE)
        publicationWriteJson(publicationHostEnvelope(), args$output)
        return(invisible(0L))
    }
    if (isTRUE(args$null_calibration)) {
        if (is.null(args$output)) stop("--output is required", call. = FALSE)
        result <- publicationNullCalibration(args$output, args$runs, args$dwell_seconds)
        return(invisible(if (identical(result$status, "passed")) 0L else 1L))
    }
    if (!is.null(args$measure_spec)) {
        if (is.null(args$output)) stop("--output is required", call. = FALSE)
        result <- publicationMeasureSpec(args$measure_spec, args$output)
        return(invisible(if (identical(result$status, "passed")) 0L else 1L))
    }
    publicationRunnerUsage()
    invisible(2L)
}

if (identical(environment(), globalenv()) && !interactive()) {
    status <- tryCatch(
        publicationRunnerMain(),
        error = \(error) {
            message(conditionMessage(error))
            1L
        }
    )
    quit(status = as.integer(status), save = "no")
}
