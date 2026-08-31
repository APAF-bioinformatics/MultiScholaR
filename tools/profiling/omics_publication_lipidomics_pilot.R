lipidPublicationPilotExecution <- function() {
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

lipidPublicationValidatePilotContract <- function(contract) {
    publicationRequireNames(
        contract,
        lipidPublicationWorkloadFields(),
        "Lipidomics pilot contract"
    )
    valid <- identical(
        contract$schema,
        "multischolar.omics_publication_lipidomics_pilot_contract"
    ) && identical(contract$schema_version, .LIPID_PUBLICATION_VERSION) &&
        identical(contract$owner_ticket_id, .LIPID_PUBLICATION_OWNER) &&
        identical(contract$status, "frozen_pilot_calibration") &&
        identical(contract$workload_class, "operational_heavy") &&
        identical(contract$assay_profile$profile_id, "mixed_lc") &&
        !isTRUE(contract$claim_scope$scientific_authority) &&
        !isTRUE(contract$claim_scope$performance_authority) &&
        !isTRUE(contract$claim_scope$cross_project_authority) &&
        !isTRUE(contract$claim_scope$promotion_authority) &&
        !isTRUE(contract$publication_authority)
    lipidPublicationValidateCapability(contract$capability)
    lipidPublicationValidateAssayProfile(
        contract$assay_profile,
        contract$capability$input_format
    )
    lipidPublicationValidateDimensions(contract$dimensions, contract$assay_profile)
    lipidPublicationValidateScale(contract)
    for (field in c(
        "parameter_authority", "source_authority", "split_authority",
        "mapping_authority", "support_boundary"
    )) {
        lipidPublicationValidateBinding(contract[[field]], field)
    }
    lipidPublicationValidateGenerator(contract$generator, contract)
    lipidPublicationValidateRng(contract$rng, contract$workload_class)
    lipidPublicationValidateTruthContract(contract$truth_contract, contract)
    lipidPublicationValidateExecution(contract$execution, contract)
    lipidPublicationValidatePrivacy(contract$privacy, contract)
    publicationRequireNames(contract$claim_scope, c(
        "evidence_class", "verified_stages", "scientific_authority",
        "performance_authority", "cross_project_authority",
        "vendor_authority", "gc_vendor_authority",
        "three_file_workflow_authority",
        "promotion_authority", "limitations"
    ), "Lipidomics pilot claim scope")
    publicationRequireNames(contract$expected_digests, c(
        "payload_set_sha256", "truth_sha256"
    ), "Lipidomics pilot digests")
    lapply(contract$expected_digests, lipidPublicationRequireDigest, "Pilot")
    splits <- publicationReadJson(contract$split_authority$path)
    lipidPublicationValidateSplits(splits)
    assignment <- Filter(function(value) identical(
        value$assignment_id,
        contract$workload_id
    ), splits$pilot_calibration_assignments)
    valid <- valid && length(assignment) == 1L && identical(
        assignment[[1L]]$seed,
        contract$rng$seed
    ) && identical(
        contract$claim_scope$evidence_class,
        "pilot_calibration"
    ) && identical(
        contract$claim_scope$verified_stages,
        contract$capability$verified_stages
    ) && !isTRUE(contract$claim_scope$vendor_authority) &&
        !isTRUE(contract$claim_scope$gc_vendor_authority) &&
        !isTRUE(contract$claim_scope$three_file_workflow_authority) &&
        length(contract$claim_scope$limitations) > 0L
    if (!valid) lipidPublicationAbort("lipidomics pilot contract differs")
    invisible(contract)
}

lipidPublicationPreparePilotGenerated <- function(
    contract,
    output_root,
    verify_expected = TRUE
) {
    lipidPublicationValidatePilotContract(contract)
    if (file.exists(output_root) || dir.exists(output_root)) {
        lipidPublicationAbort("lipidomics pilot preparation root exists")
    }
    dir.create(output_root, recursive = TRUE, showWarnings = FALSE)
    parameters <- publicationReadJson(contract$parameter_authority$path)
    lipidPublicationValidateParameters(parameters)
    plan <- lipidPublicationModelPlan(
        contract$dimensions$assay_feature_counts,
        contract$dimensions$sample_count,
        parameters,
        contract$rng$streams,
        contract$generator$chunk_rows
    )
    state <- lipidPublicationTruthState(plan)
    payload <- lipidPublicationSerialize(
        plan,
        contract$capability$input_format,
        contract$assay_profile$profile_id,
        file.path(output_root, "payload"),
        contract$generator$chunk_rows,
        observer = function(block, index) {
            lipidPublicationObserveTruth(state, block, index)
        }
    )
    truth_record <- lipidPublicationFinalizeTruth(
        state,
        plan,
        contract,
        payload
    )
    lipidPublicationValidateTruth(truth_record, contract)
    truth <- lipidPublicationWriteTruth(
        truth_record,
        file.path(output_root, contract$generator$truth_filename)
    )
    receipt <- lipidPublicationPreparationReceipt(contract, payload, truth)
    lipidPublicationValidatePreparation(receipt, contract)
    receipt_path <- file.path(output_root, "preparation-receipt.json")
    publicationWriteJson(receipt, receipt_path)
    if (isTRUE(verify_expected) && (!identical(
        payload$payload_set_sha256,
        contract$expected_digests$payload_set_sha256
    ) || !identical(
        truth$sha256,
        contract$expected_digests$truth_sha256
    ))) {
        lipidPublicationAbort("lipidomics pilot preparation digest differs")
    }
    list(
        payload = payload,
        truth = truth,
        receipt = list(
            path = receipt_path,
            sha256 = publicationFileDigest(receipt_path),
            record = receipt
        )
    )
}

lipidPublicationBuildPilotContract <- function(
    route,
    aggregate_feature_count,
    sample_count,
    grid_id
) {
    if (!route %in% names(lipidPublicationCapabilities()) ||
        !publicationScalarString(grid_id)) {
        lipidPublicationAbort("lipidomics pilot identity is invalid")
    }
    paths <- lipidPublicationAuthorityPaths()
    splits <- publicationReadJson(paths$splits)
    assignments <- Filter(function(value) {
        identical(value$route, route) &&
            identical(value$profile_id, "mixed_lc")
    }, splits$pilot_calibration_assignments)
    if (length(assignments) != 1L) {
        lipidPublicationAbort("lipidomics pilot split is absent")
    }
    assignment <- assignments[[1L]]
    counts <- lipidPublicationAssayCounts("mixed_lc", aggregate_feature_count)
    member_count <- 2L
    dimensions <- list(
        aggregate_feature_count = as.integer(aggregate_feature_count),
        assay_feature_counts = counts,
        sample_count = as.integer(sample_count),
        assay_count = 2L,
        quantity_count = as.numeric(aggregate_feature_count * sample_count),
        payload_member_count = member_count
    )
    generator <- lipidPublicationGeneratedGenerator(route, "mixed_lc")
    pilot_path <- "tools/profiling/omics_publication_lipidomics_pilot.R"
    generator$source_bindings[[length(generator$source_bindings) + 1L]] <-
        lipidPublicationAuthorityBinding(pilot_path)
    contract <- list(
        schema = "multischolar.omics_publication_lipidomics_pilot_contract",
        schema_version = .LIPID_PUBLICATION_VERSION,
        contract_id = paste(
            "lipidomics", route, "mixed_lc", "pilot", grid_id, "v1", sep = "."
        ),
        owner_ticket_id = .LIPID_PUBLICATION_OWNER,
        status = "frozen_pilot_calibration",
        workload_id = assignment$assignment_id,
        workload_class = "operational_heavy",
        capability = lipidPublicationExpectedCapability(route),
        assay_profile = list(
            profile_id = "mixed_lc",
            assays = lipidPublicationAssayProfiles()$mixed_lc$assays,
            payload_mode = "bundle",
            member_count = member_count
        ),
        dimensions = dimensions,
        model_profile_id = paste(
            assignment$generator_parameter_profile_id,
            grid_id,
            sep = "."
        ),
        parameter_authority = lipidPublicationAuthorityBinding(paths$parameters),
        source_authority = lipidPublicationAuthorityBinding(paths$sources),
        split_authority = lipidPublicationAuthorityBinding(paths$splits),
        mapping_authority = lipidPublicationAuthorityBinding(paths$mappings),
        contract_mapping_id = paste(
            "lipidomics.mapping",
            route,
            "mixed_lc",
            sep = "."
        ),
        support_boundary = lipidPublicationAuthorityBinding(paths$boundaries),
        generator = generator,
        rng = lipidPublicationRng(assignment, "operational_heavy"),
        truth_contract = lipidPublicationTruthContract(
            route,
            "mixed_lc",
            "operational_heavy"
        ),
        execution = lipidPublicationExecution(dimensions, "operational_heavy"),
        privacy = lipidPublicationPrivacy(route, "operational_heavy"),
        claim_scope = list(
            evidence_class = "pilot_calibration",
            verified_stages =
                lipidPublicationCapabilities()[[route]]$verified_stages,
            scientific_authority = FALSE,
            performance_authority = FALSE,
            cross_project_authority = FALSE,
            vendor_authority = FALSE,
            gc_vendor_authority = FALSE,
            three_file_workflow_authority = FALSE,
            promotion_authority = FALSE,
            limitations = list(
                "Historical-only dimension qualification; no candidate access.",
                "Pilot observations cannot become confirmatory evidence."
            )
        ),
        expected_digests = list(
            payload_set_sha256 = strrep("0", 64L),
            truth_sha256 = strrep("0", 64L)
        ),
        publication_authority = FALSE
    )
    lipidPublicationValidatePilotContract(contract)
    contract
}

lipidPublicationPilotThreadEnvironment <- function(
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

lipidPublicationPilotPhaseSnapshot <- function() {
    cgroup_path <- publicationCurrentCgroupPath()
    if (is.null(cgroup_path) || !dir.exists(cgroup_path)) {
        lipidPublicationAbort("pilot worker cgroup is unavailable")
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

lipidPublicationWritePilotBoundary <- function(run_dir, marker) {
    record <- list(
        marker = marker,
        pid = Sys.getpid(),
        counters = lipidPublicationPilotPhaseSnapshot()
    )
    publicationWriteJson(record, file.path(run_dir, paste0(marker, ".json")))
    invisible(record)
}

lipidPublicationPilotResourceLedger <- function() {
    empty <- publicationEmptyWorkerResources()
    list(
        schema = "multischolar.omics_publication_worker_resources",
        schema_version = "1.0.0",
        high_water = empty,
        retained = empty,
        terminal = empty
    )
}

lipidPublicationPilotRetentionState <- function(dwell_seconds) {
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

lipidPublicationPilotWorker <- function(
    contract_path,
    payload_root,
    truth_path,
    run_dir,
    package_library,
    dependency_library,
    historical_revision,
    dwell_seconds = 5
) {
    dir.create(run_dir, recursive = TRUE, showWarnings = FALSE)
    contract <- publicationReadJson(contract_path)
    lipidPublicationValidatePilotContract(contract)
    member_names <- unlist(contract$generator$output_members, use.names = FALSE)
    member_paths <- file.path(payload_root, member_names)
    payload <- lipidPublicationPayloadBinding(member_paths)
    if (!identical(
        payload$payload_set_sha256,
        contract$expected_digests$payload_set_sha256
    ) || !identical(
        publicationFileDigest(truth_path),
        contract$expected_digests$truth_sha256
    )) {
        lipidPublicationAbort("pilot payload or truth differs")
    }
    .libPaths(c(package_library, dependency_library, .Library))
    namespace <- loadNamespace("MultiScholaR")
    package_path <- find.package("MultiScholaR", lib.loc = package_library)
    if (!identical(normalizePath(package_path), normalizePath(file.path(
        package_library,
        "MultiScholaR"
    )))) {
        lipidPublicationAbort("historical pilot package path differs")
    }
    importer <- get(
        if (identical(
            contract$capability$input_format,
            "lipidsearch"
        )) "importLipidSearchData" else "importLipidMSDIALData",
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
    if (fifo$status != 0L) lipidPublicationAbort("pilot FIFO creation failed")
    started <- proc.time()[["elapsed"]]
    cpu_started <- sum(proc.time()[c("user.self", "sys.self")])
    lipidPublicationWritePilotBoundary(run_dir, "measured_worker_start")
    imported <- lapply(seq_along(member_paths), function(index) {
        arguments <- list(filepath = member_paths[[index]])
        if (identical(contract$capability$input_format, "custom")) {
            arguments$lipid_id_column <- "lipid_id"
            arguments$annotation_column <- "annotation"
        }
        list(
            filename = member_names[[index]],
            assay = lipidPublicationMemberAssay(member_names[[index]], contract),
            imported = suppressMessages(suppressWarnings(do.call(
                importer,
                arguments
            )))
        )
    })
    phase_elapsed <- proc.time()[["elapsed"]] - started
    phase_cpu <- sum(proc.time()[c("user.self", "sys.self")]) - cpu_started
    lipidPublicationWritePilotBoundary(run_dir, "measured_worker_stop")
    truth <- publicationReadJson(truth_path)
    lipidPublicationValidateTruth(truth, contract)
    summary <- lipidPublicationImportSummary(contract, imported)
    lipidPublicationValidateImportSummary(summary, truth)
    direction <- lipidPublicationDirectionEvidence(summary, truth)
    if (!isTRUE(direction$valid)) {
        lipidPublicationAbort("pilot lipidomics direction differs")
    }
    payload_bytes <- sum(file.info(member_paths)$size)
    publicationWriteJson(
        list(
            schema = "multischolar.omics_publication_lipidomics_pilot_worker",
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
            work_count = as.numeric(payload_bytes),
            imported_rows = as.numeric(summary$aggregate_feature_count),
            imported_columns = as.integer(
                contract$dimensions$sample_count + 9L
            ),
            imported_data_type = "lipid",
            member_count = as.integer(length(imported)),
            truth_valid = TRUE
        ),
        file.path(run_dir, "worker-result.json")
    )
    publicationWriteJson(
        lipidPublicationPilotResourceLedger(),
        file.path(run_dir, "worker-resources.json")
    )
    publicationWriteJson(
        lipidPublicationPilotRetentionState(dwell_seconds),
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
        lipidPublicationAbort("pilot retention acknowledgement differs")
    }
    rm(imported, summary, direction)
    gc(verbose = FALSE)
    publicationWriteJson(
        list(marker = "worker_exit", pid = Sys.getpid()),
        file.path(run_dir, "worker-exit.json")
    )
    invisible(0L)
}

lipidPublicationReadPilotPhase <- function(run_dir) {
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
    if (!valid) lipidPublicationAbort("pilot phase evidence differs")
    list(start = start, stop = stop, worker = worker)
}

lipidPublicationAttachPilotPhase <- function(measurement, run_dir) {
    if (!identical(measurement$status, "passed")) return(measurement)
    phase <- lipidPublicationReadPilotPhase(run_dir)
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
        lipidPublicationAbort("pilot phase counters decreased")
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

lipidPublicationPilotSafetyLimits <- function(host, preflight) {
    publicationDynamicRuntimeSafetyLimits(
        host,
        preflight,
        maximum_run_allocated_disk_bytes = 53687091200
    )
}

lipidPublicationPilotQualification <- function(measurement) {
    publicationPilotQualificationDecision(measurement)
}

lipidPublicationPilotStatus <- function(measurement, qualification) {
    publicationPilotTerminalStatus(measurement, qualification)
}

lipidPublicationPilotRecord <- function(
    contract_path,
    payload_root,
    truth_path,
    build_receipt_path,
    host,
    preflight,
    measurement,
    raw_path
) {
    contract <- publicationReadJson(contract_path)
    lipidPublicationValidatePilotContract(contract)
    members <- file.path(
        payload_root,
        unlist(contract$generator$output_members, use.names = FALSE)
    )
    payload <- lipidPublicationPayloadBinding(members)
    qualification <- lipidPublicationPilotQualification(measurement)
    list(
        schema = "multischolar.omics_publication_lipidomics_pilot",
        schema_version = "1.0.0",
        owner_ticket_id = "OMICS-ART-067",
        status = lipidPublicationPilotStatus(measurement, qualification),
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
        payload = payload,
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
