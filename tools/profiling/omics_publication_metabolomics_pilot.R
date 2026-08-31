metabPublicationPilotExecution <- function() {
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

metabPublicationValidatePilotContract <- function(contract) {
    publicationRequireNames(
        contract,
        metabPublicationWorkloadFields(),
        "Metabolomics pilot contract"
    )
    valid <- identical(
        contract$schema,
        "multischolar.omics_publication_metabolomics_pilot_contract"
    ) && identical(contract$schema_version, .METAB_PUBLICATION_VERSION) &&
        identical(contract$owner_ticket_id, .METAB_PUBLICATION_OWNER) &&
        identical(contract$status, "frozen_pilot_calibration") &&
        identical(contract$workload_class, "operational_heavy") &&
        identical(contract$assay_profile$profile_id, "mixed") &&
        !isTRUE(contract$claim_scope$scientific_authority) &&
        !isTRUE(contract$claim_scope$performance_authority) &&
        !isTRUE(contract$claim_scope$cross_project_authority) &&
        !isTRUE(contract$claim_scope$promotion_authority) &&
        !isTRUE(contract$publication_authority)
    metabPublicationValidateCapability(contract$capability)
    metabPublicationValidateAssayProfile(
        contract$assay_profile,
        contract$capability$input_format
    )
    metabPublicationValidateDimensions(contract$dimensions, contract$assay_profile)
    metabPublicationValidateScale(contract)
    for (field in c(
        "parameter_authority", "source_authority", "split_authority",
        "mapping_authority", "support_boundary"
    )) {
        metabPublicationValidateBinding(contract[[field]], field)
    }
    metabPublicationValidateGenerator(contract$generator, contract)
    metabPublicationValidateRng(contract$rng, contract$workload_class)
    metabPublicationValidateTruthContract(contract$truth_contract, contract)
    metabPublicationValidateExecution(contract$execution, contract)
    metabPublicationValidatePrivacy(contract$privacy, contract)
    publicationRequireNames(contract$claim_scope, c(
        "evidence_class", "verified_stages", "scientific_authority",
        "performance_authority", "cross_project_authority",
        "promotion_authority", "limitations"
    ), "Metabolomics pilot claim scope")
    publicationRequireNames(contract$expected_digests, c(
        "payload_set_sha256", "truth_sha256"
    ), "Metabolomics pilot digests")
    lapply(contract$expected_digests, metabPublicationRequireDigest, "Pilot")
    splits <- publicationReadJson(contract$split_authority$path)
    metabPublicationValidateSplits(splits)
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
    ) && length(contract$claim_scope$limitations) > 0L
    if (!valid) metabPublicationAbort("metabolomics pilot contract differs")
    invisible(contract)
}

metabPublicationPreparePilotGenerated <- function(
    contract,
    output_root,
    verify_expected = TRUE
) {
    metabPublicationValidatePilotContract(contract)
    if (file.exists(output_root) || dir.exists(output_root)) {
        metabPublicationAbort("metabolomics pilot preparation root exists")
    }
    dir.create(output_root, recursive = TRUE, showWarnings = FALSE)
    parameters <- publicationReadJson(contract$parameter_authority$path)
    metabPublicationValidateParameters(parameters)
    plan <- metabPublicationModelPlan(
        contract$dimensions$assay_feature_counts,
        contract$dimensions$sample_count,
        parameters,
        contract$rng$streams,
        contract$generator$chunk_rows
    )
    state <- metabPublicationTruthState(plan)
    payload <- metabPublicationSerialize(
        plan,
        contract$capability$input_format,
        contract$assay_profile$profile_id,
        file.path(output_root, "payload"),
        contract$generator$chunk_rows,
        observer = function(block, index) {
            metabPublicationObserveTruth(state, block, index)
        }
    )
    truth_record <- metabPublicationFinalizeTruth(
        state,
        plan,
        contract,
        payload
    )
    metabPublicationValidateTruth(truth_record, contract)
    truth <- metabPublicationWriteTruth(
        truth_record,
        file.path(output_root, contract$generator$truth_filename)
    )
    receipt <- metabPublicationPreparationReceipt(contract, payload, truth)
    metabPublicationValidatePreparation(receipt, contract)
    receipt_path <- file.path(output_root, "preparation-receipt.json")
    publicationWriteJson(receipt, receipt_path)
    if (isTRUE(verify_expected) && (!identical(
        payload$payload_set_sha256,
        contract$expected_digests$payload_set_sha256
    ) || !identical(
        truth$sha256,
        contract$expected_digests$truth_sha256
    ))) {
        metabPublicationAbort("metabolomics pilot preparation digest differs")
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

metabPublicationBuildPilotContract <- function(
    route,
    aggregate_feature_count,
    sample_count,
    grid_id
) {
    if (!route %in% names(metabPublicationCapabilities()) ||
        !publicationScalarString(grid_id)) {
        metabPublicationAbort("metabolomics pilot identity is invalid")
    }
    paths <- metabPublicationAuthorityPaths()
    splits <- publicationReadJson(paths$splits)
    assignments <- Filter(function(value) {
        identical(value$route, route) &&
            identical(value$profile_id, "mixed")
    }, splits$pilot_calibration_assignments)
    if (length(assignments) != 1L) {
        metabPublicationAbort("metabolomics pilot split is absent")
    }
    assignment <- assignments[[1L]]
    counts <- metabPublicationAssayCounts("mixed", aggregate_feature_count)
    member_count <- if (identical(route, "msdial")) 3L else 1L
    dimensions <- list(
        aggregate_feature_count = as.integer(aggregate_feature_count),
        assay_feature_counts = counts,
        sample_count = as.integer(sample_count),
        assay_count = 3L,
        quantity_count = as.numeric(aggregate_feature_count * sample_count),
        payload_member_count = member_count
    )
    generator <- metabPublicationGeneratedGenerator(route, "mixed")
    pilot_path <- "tools/profiling/omics_publication_metabolomics_pilot.R"
    generator$source_bindings[[length(generator$source_bindings) + 1L]] <-
        metabPublicationAuthorityBinding(pilot_path)
    contract <- list(
        schema = "multischolar.omics_publication_metabolomics_pilot_contract",
        schema_version = .METAB_PUBLICATION_VERSION,
        contract_id = paste(
            "metabolomics", route, "mixed", "pilot", grid_id, "v1", sep = "."
        ),
        owner_ticket_id = .METAB_PUBLICATION_OWNER,
        status = "frozen_pilot_calibration",
        workload_id = assignment$assignment_id,
        workload_class = "operational_heavy",
        capability = metabPublicationExpectedCapability(route),
        assay_profile = list(
            profile_id = "mixed",
            assays = metabPublicationAssayProfiles()$mixed$assays,
            payload_mode = if (member_count == 3L) "bundle" else "single",
            member_count = member_count
        ),
        dimensions = dimensions,
        model_profile_id = paste(
            assignment$generator_parameter_profile_id,
            grid_id,
            sep = "."
        ),
        parameter_authority = metabPublicationAuthorityBinding(paths$parameters),
        source_authority = metabPublicationAuthorityBinding(paths$sources),
        split_authority = metabPublicationAuthorityBinding(paths$splits),
        mapping_authority = metabPublicationAuthorityBinding(paths$mappings),
        contract_mapping_id = paste(
            "metabolomics.mapping",
            route,
            "mixed",
            sep = "."
        ),
        support_boundary = metabPublicationAuthorityBinding(paths$boundaries),
        generator = generator,
        rng = metabPublicationRng(assignment, "operational_heavy"),
        truth_contract = metabPublicationTruthContract("operational_heavy"),
        execution = metabPublicationExecution(dimensions, "operational_heavy"),
        privacy = metabPublicationPrivacy("operational_heavy"),
        claim_scope = list(
            evidence_class = "pilot_calibration",
            verified_stages =
                metabPublicationCapabilities()[[route]]$verified_stages,
            scientific_authority = FALSE,
            performance_authority = FALSE,
            cross_project_authority = FALSE,
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
    metabPublicationValidatePilotContract(contract)
    contract
}

metabPublicationPilotThreadEnvironment <- function(
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

metabPublicationPilotPhaseSnapshot <- function() {
    cgroup_path <- publicationCurrentCgroupPath()
    if (is.null(cgroup_path) || !dir.exists(cgroup_path)) {
        metabPublicationAbort("pilot worker cgroup is unavailable")
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

metabPublicationWritePilotBoundary <- function(run_dir, marker) {
    record <- list(
        marker = marker,
        pid = Sys.getpid(),
        counters = metabPublicationPilotPhaseSnapshot()
    )
    publicationWriteJson(record, file.path(run_dir, paste0(marker, ".json")))
    invisible(record)
}

metabPublicationPilotResourceLedger <- function() {
    empty <- publicationEmptyWorkerResources()
    list(
        schema = "multischolar.omics_publication_worker_resources",
        schema_version = "1.0.0",
        high_water = empty,
        retained = empty,
        terminal = empty
    )
}

metabPublicationPilotRetentionState <- function(dwell_seconds) {
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

metabPublicationPilotWorker <- function(
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
    metabPublicationValidatePilotContract(contract)
    member_names <- unlist(contract$generator$output_members, use.names = FALSE)
    member_paths <- file.path(payload_root, member_names)
    payload <- metabPublicationPayloadBinding(member_paths)
    if (!identical(
        payload$payload_set_sha256,
        contract$expected_digests$payload_set_sha256
    ) || !identical(
        publicationFileDigest(truth_path),
        contract$expected_digests$truth_sha256
    )) {
        metabPublicationAbort("pilot payload or truth differs")
    }
    .libPaths(c(package_library, dependency_library, .Library))
    namespace <- loadNamespace("MultiScholaR")
    package_path <- find.package("MultiScholaR", lib.loc = package_library)
    if (!identical(normalizePath(package_path), normalizePath(file.path(
        package_library,
        "MultiScholaR"
    )))) {
        metabPublicationAbort("historical pilot package path differs")
    }
    importer <- get(
        "importMetabMSDIALData",
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
    if (fifo$status != 0L) metabPublicationAbort("pilot FIFO creation failed")
    started <- proc.time()[["elapsed"]]
    cpu_started <- sum(proc.time()[c("user.self", "sys.self")])
    metabPublicationWritePilotBoundary(run_dir, "measured_worker_start")
    imported <- lapply(seq_along(member_paths), function(index) {
        arguments <- list(filepath = member_paths[[index]])
        if (identical(contract$capability$input_format, "custom")) {
            arguments$metabolite_id_column <- "feature_id"
            arguments$annotation_column <- "annotation"
        }
        list(
            filename = member_names[[index]],
            assay = metabPublicationMemberAssay(member_names[[index]], contract),
            imported = suppressMessages(suppressWarnings(do.call(
                importer,
                arguments
            )))
        )
    })
    phase_elapsed <- proc.time()[["elapsed"]] - started
    phase_cpu <- sum(proc.time()[c("user.self", "sys.self")]) - cpu_started
    metabPublicationWritePilotBoundary(run_dir, "measured_worker_stop")
    truth <- publicationReadJson(truth_path)
    metabPublicationValidateTruth(truth, contract)
    summary <- metabPublicationImportSummary(contract, imported)
    metabPublicationValidateImportSummary(summary, truth)
    direction <- metabPublicationDirectionEvidence(summary, truth)
    if (!isTRUE(direction$valid)) {
        metabPublicationAbort("pilot metabolomics direction differs")
    }
    payload_bytes <- sum(file.info(member_paths)$size)
    publicationWriteJson(
        list(
            schema = "multischolar.omics_publication_metabolomics_pilot_worker",
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
                contract$dimensions$sample_count + 5L
            ),
            imported_data_type = "metabolite",
            member_count = as.integer(length(imported)),
            truth_valid = TRUE
        ),
        file.path(run_dir, "worker-result.json")
    )
    publicationWriteJson(
        metabPublicationPilotResourceLedger(),
        file.path(run_dir, "worker-resources.json")
    )
    publicationWriteJson(
        metabPublicationPilotRetentionState(dwell_seconds),
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
        metabPublicationAbort("pilot retention acknowledgement differs")
    }
    rm(imported, summary, direction)
    gc(verbose = FALSE)
    publicationWriteJson(
        list(marker = "worker_exit", pid = Sys.getpid()),
        file.path(run_dir, "worker-exit.json")
    )
    invisible(0L)
}

metabPublicationReadPilotPhase <- function(run_dir) {
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
    if (!valid) metabPublicationAbort("pilot phase evidence differs")
    list(start = start, stop = stop, worker = worker)
}

metabPublicationAttachPilotPhase <- function(measurement, run_dir) {
    if (!identical(measurement$status, "passed")) return(measurement)
    phase <- metabPublicationReadPilotPhase(run_dir)
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
        metabPublicationAbort("pilot phase counters decreased")
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

metabPublicationPilotSafetyLimits <- function(host, preflight) {
    publicationDynamicRuntimeSafetyLimits(
        host,
        preflight,
        maximum_run_allocated_disk_bytes = 53687091200
    )
}

metabPublicationPilotQualification <- function(measurement) {
    publicationPilotQualificationDecision(measurement)
}

metabPublicationPilotStatus <- function(measurement, qualification) {
    publicationPilotTerminalStatus(measurement, qualification)
}

metabPublicationPilotRecord <- function(
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
    metabPublicationValidatePilotContract(contract)
    members <- file.path(
        payload_root,
        unlist(contract$generator$output_members, use.names = FALSE)
    )
    payload <- metabPublicationPayloadBinding(members)
    qualification <- metabPublicationPilotQualification(measurement)
    list(
        schema = "multischolar.omics_publication_metabolomics_pilot",
        schema_version = "1.0.0",
        owner_ticket_id = "OMICS-ART-066",
        status = metabPublicationPilotStatus(measurement, qualification),
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
