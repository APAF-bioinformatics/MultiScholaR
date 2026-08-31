proteomicsPublicationPilotExecution <- function() {
    list(
        sampling_interval_ms = 20,
        disk_sampling_interval_ms = 500,
        timeout_seconds = 1800,
        retained_dwell_seconds = 5,
        retained_window_seconds = 2,
        retention_acknowledgement = "fifo_v1",
        maximum_boundary_bracket_seconds = 0.5,
        retained_boundary_tolerance_ms = 100
    )
}

proteomicsPublicationPilotThreadEnvironment <- function(
    home,
    temp,
    package_library,
    dependency_library,
    site_library
) {
    library_path <- paste(
        c(package_library, dependency_library),
        collapse = .Platform$path.sep
    )
    c(
        OMP_NUM_THREADS = "1",
        OPENBLAS_NUM_THREADS = "1",
        MKL_NUM_THREADS = "1",
        ARROW_NUM_THREADS = "1",
        DUCKDB_THREADS = "1",
        TZ = "UTC",
        LANG = "C.UTF-8",
        LC_ALL = "C.UTF-8",
        R_PROFILE_USER = "",
        R_ENVIRON_USER = "",
        R_LIBS = library_path,
        R_LIBS_USER = library_path,
        R_LIBS_SITE = site_library,
        HOME = home,
        TMPDIR = temp
    )
}

proteomicsPublicationPilotPhaseSnapshot <- function() {
    cgroup_path <- publicationCurrentCgroupPath()
    if (is.null(cgroup_path) || !dir.exists(cgroup_path)) {
        proteomicsPublicationAbort("pilot worker cgroup is unavailable")
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

proteomicsPublicationWritePilotBoundary <- function(run_dir, marker) {
    record <- list(
        marker = marker,
        pid = Sys.getpid(),
        counters = proteomicsPublicationPilotPhaseSnapshot()
    )
    publicationWriteJson(record, file.path(run_dir, paste0(marker, ".json")))
    invisible(record)
}

proteomicsPublicationPilotImporterName <- function(format) {
    switch(
        format,
        diann = "importDIANNData",
        maxquant = "importMaxQuantData",
        fragpipe = "importFragPipeData",
        pd_tmt = "importProteomeDiscovererTMTData",
        proteomicsPublicationAbort("pilot format is unsupported")
    )
}

proteomicsPublicationPilotResourceLedger <- function() {
    empty <- publicationEmptyWorkerResources()
    list(
        schema = "multischolar.omics_publication_worker_resources",
        schema_version = "1.0.0",
        high_water = empty,
        retained = empty,
        terminal = empty
    )
}

proteomicsPublicationPilotRetentionState <- function(dwell_seconds) {
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
    )
}

proteomicsPublicationPilotWorker <- function(
    contract_path,
    payload_path,
    truth_path,
    run_dir,
    package_library,
    dependency_library,
    historical_revision,
    dwell_seconds = 5
) {
    dir.create(run_dir, recursive = TRUE, showWarnings = FALSE)
    contract <- proteomicsPublicationReadContract(contract_path)
    if (!identical(contract$workload_class, "operational_heavy") ||
        isTRUE(contract$claim_scope$performance_authority)) {
        proteomicsPublicationAbort("pilot contract is not calibration-only")
    }
    if (!identical(
        publicationFileDigest(payload_path),
        contract$expected_digests$payload_sha256
    ) || !identical(
        publicationFileDigest(truth_path),
        contract$expected_digests$truth_sha256
    )) {
        proteomicsPublicationAbort("pilot payload or truth differs")
    }
    .libPaths(c(package_library, dependency_library, .Library))
    namespace <- loadNamespace("MultiScholaR")
    package_path <- find.package("MultiScholaR", lib.loc = package_library)
    if (!identical(normalizePath(package_path), normalizePath(file.path(
        package_library,
        "MultiScholaR"
    )))) {
        proteomicsPublicationAbort("historical pilot package path differs")
    }
    importer <- get(
        proteomicsPublicationPilotImporterName(
            contract$capability$input_format
        ),
        envir = namespace,
        inherits = FALSE
    )
    acknowledgement_path <- file.path(run_dir, "retention-sampled.fifo")
    fifo <- processx::run(
        "mkfifo",
        acknowledgement_path,
        error_on_status = FALSE,
        echo = FALSE
    )
    if (fifo$status != 0L) proteomicsPublicationAbort("pilot FIFO creation failed")
    started <- proc.time()[["elapsed"]]
    cpu_started <- sum(proc.time()[c("user.self", "sys.self")])
    proteomicsPublicationWritePilotBoundary(run_dir, "measured_worker_start")
    imported <- suppressMessages(suppressWarnings(importer(payload_path)))
    phase_elapsed <- proc.time()[["elapsed"]] - started
    phase_cpu <- sum(proc.time()[c("user.self", "sys.self")]) - cpu_started
    proteomicsPublicationWritePilotBoundary(run_dir, "measured_worker_stop")
    data <- imported$data
    publicationWriteJson(
        list(
            schema = "multischolar.omics_publication_proteomics_pilot_worker",
            schema_version = "1.0.0",
            status = "passed",
            pid = Sys.getpid(),
            comparator_role = "historical_janitor",
            comparator_revision = historical_revision,
            package_path_sha256 = digest::digest(
                normalizePath(package_path),
                algo = "sha256",
                serialize = FALSE
            ),
            candidate_loaded = FALSE,
            phase_elapsed_seconds = phase_elapsed,
            phase_cpu_seconds = phase_cpu,
            work_unit_id = "validated_input_byte",
            work_count = as.numeric(file.info(payload_path)$size),
            imported_rows = as.numeric(nrow(data)),
            imported_columns = as.integer(ncol(data)),
            imported_data_type = imported$data_type
        ),
        file.path(run_dir, "worker-result.json")
    )
    publicationWriteJson(
        proteomicsPublicationPilotResourceLedger(),
        file.path(run_dir, "worker-resources.json")
    )
    publicationWriteJson(
        proteomicsPublicationPilotRetentionState(dwell_seconds),
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
        proteomicsPublicationAbort("pilot retention acknowledgement differs")
    }
    rm(data, imported)
    gc(verbose = FALSE)
    publicationWriteJson(
        list(marker = "worker_exit", pid = Sys.getpid()),
        file.path(run_dir, "worker-exit.json")
    )
    invisible(0L)
}

proteomicsPublicationReadPilotPhase <- function(run_dir) {
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
    valid <- identical(start$marker, "measured_worker_start") &&
        identical(stop$marker, "measured_worker_stop") &&
        identical(as.integer(start$pid), as.integer(worker$pid)) &&
        identical(as.integer(stop$pid), as.integer(worker$pid)) &&
        identical(worker$status, "passed") &&
        all(counter_names %in% names(start$counters$cpu)) &&
        all(counter_names %in% names(stop$counters$cpu))
    if (!valid) proteomicsPublicationAbort("pilot phase evidence differs")
    list(start = start, stop = stop, worker = worker)
}

proteomicsPublicationAttachPilotPhase <- function(measurement, run_dir) {
    if (!identical(measurement$status, "passed")) return(measurement)
    phase <- proteomicsPublicationReadPilotPhase(run_dir)
    cpu_names <- c("usage_usec", "user_usec", "system_usec")
    cpu <- stats::setNames(lapply(cpu_names, \(name) {
        as.numeric(phase$stop$counters$cpu[[name]]) -
            as.numeric(phase$start$counters$cpu[[name]])
    }), cpu_names)
    io <- stats::setNames(lapply(names(phase$start$counters$io), \(name) {
        as.numeric(phase$stop$counters$io[[name]]) -
            as.numeric(phase$start$counters$io[[name]])
    }), names(phase$start$counters$io))
    if (any(unlist(c(cpu, io), use.names = FALSE) < 0)) {
        proteomicsPublicationAbort("pilot phase counters decreased")
    }
    worker <- phase$worker
    phase_cpu <- cpu$usage_usec / 1e6
    measurement$phase_evidence <- list(
        valid = TRUE,
        worker_pid = worker$pid,
        work_unit_id = worker$work_unit_id,
        work_count = worker$work_count,
        phase_elapsed_seconds = worker$phase_elapsed_seconds,
        phase_cgroup_cpu_seconds = phase_cpu,
        phase_process_self_cpu_seconds = worker$phase_cpu_seconds,
        phase_cgroup_user_seconds = cpu$user_usec / 1e6,
        phase_cgroup_system_seconds = cpu$system_usec / 1e6,
        phase_cgroup_io = io
    )
    measurement$metrics$cgroup_lifecycle_elapsed_seconds <-
        measurement$metrics$elapsed_seconds
    measurement$metrics$elapsed_seconds <- worker$phase_elapsed_seconds
    measurement$metrics$phase_cpu_seconds <- phase_cpu
    measurement$metrics$phase_process_self_cpu_seconds <- worker$phase_cpu_seconds
    measurement$metrics$phase_io <- io
    measurement$metrics$primary_work_units_per_wall_second <-
        worker$work_count / worker$phase_elapsed_seconds
    measurement$metrics$primary_work_units_per_cpu_second <-
        worker$work_count / phase_cpu
    measurement$worker <- worker
    measurement
}

proteomicsPublicationPilotSafetyLimits <- function(host, preflight) {
    publicationDynamicRuntimeSafetyLimits(
        host,
        preflight,
        maximum_run_allocated_disk_bytes = 53687091200
    )
}

proteomicsPublicationPilotQualification <- function(measurement) {
    publicationPilotQualificationDecision(measurement)
}

proteomicsPublicationPilotStatus <- function(measurement, qualification) {
    publicationPilotTerminalStatus(measurement, qualification)
}

proteomicsPublicationPilotRecord <- function(
    contract_path,
    payload_path,
    truth_path,
    build_receipt_path,
    host,
    preflight,
    measurement,
    raw_path
) {
    contract <- proteomicsPublicationReadContract(contract_path)
    qualification <- proteomicsPublicationPilotQualification(measurement)
    list(
        schema = "multischolar.omics_publication_proteomics_pilot",
        schema_version = "1.0.0",
        owner_ticket_id = "OMICS-ART-065",
        status = proteomicsPublicationPilotStatus(measurement, qualification),
        comparator_role = "historical_janitor",
        comparator_revision = publicationComparatorRevisions()[[
            "historical_janitor"
        ]],
        build_receipt = list(
            path = build_receipt_path,
            sha256 = publicationFileDigest(build_receipt_path)
        ),
        contract = list(
            path = contract_path,
            sha256 = publicationFileDigest(contract_path),
            workload_id = contract$workload_id
        ),
        payload = list(
            sha256 = publicationFileDigest(payload_path),
            size_bytes = as.numeric(file.info(payload_path)$size)
        ),
        truth = list(sha256 = publicationFileDigest(truth_path)),
        host_preflight = list(
            sha256 = publicationObjectDigest(preflight),
            certified = preflight$certified
        ),
        host_scope = list(
            kernel = host$os$release,
            machine = host$os$machine,
            logical_cores = host$cpu$logical_cores
        ),
        measurement = list(
            raw_sha256 = publicationFileDigest(raw_path),
            status = measurement$status,
            publication_certifiable = measurement$publication_certifiable,
            timed_out = measurement$timed_out,
            safety_aborted = measurement$safety_aborted,
            safety_reason = measurement$safety_reason,
            phase_evidence = measurement$phase_evidence,
            safety_evidence = measurement$safety_evidence,
            metrics = measurement$metrics,
            cleanup = measurement$cleanup,
            worker = measurement$worker
        ),
        qualification = qualification,
        candidate_loaded = FALSE,
        candidate_inspected = FALSE,
        promotion_authority = FALSE,
        publication_authority = FALSE
    )
}
