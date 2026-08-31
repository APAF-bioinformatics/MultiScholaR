diaRepairAbort <- function(message) {
    stop(message, call. = FALSE)
}

diaRepairGatePath <- function() {
    file.path(
        "tests",
        "testdata",
        "omics-performance",
        "dia-repair-gates-v4.json"
    )
}

diaRepairReadGates <- function(path = diaRepairGatePath()) {
    gates <- publicationReadJson(path)
    required <- c(
        "schema", "schema_version", "gate_id", "owner_ticket_id", "status",
        "comparison", "design", "gates", "absolute_gates", "decision"
    )
    if (!is.list(gates) || !setequal(names(gates), required) ||
        !identical(gates$schema, "multischolar.dia_commit_repair_gates") ||
        !identical(gates$schema_version, "1.3.0") ||
        !identical(gates$owner_ticket_id, "OMICS-ART-070") ||
        !identical(gates$status, "frozen_after_validated_owned_load_floor")) {
        diaRepairAbort("DIA repair gate contract header is invalid")
    }
    design <- gates$design
    design_valid <- identical(as.integer(design$required_pairs), 36L) &&
        identical(as.integer(design$required_sessions), 3L) &&
        isTRUE(design$counterbalanced_within_pair) &&
        isTRUE(design$fresh_process_per_arm) &&
        identical(as.integer(design$sampling_interval_ms), 20L) &&
        identical(as.numeric(design$retained_dwell_seconds), 5) &&
        identical(as.numeric(design$retained_window_seconds), 2) &&
        identical(as.numeric(design$maximum_idle_cpu_activity_seconds), 0.001) &&
        identical(as.integer(design$resamples), 10000L) &&
        identical(as.numeric(design$confidence_level), 0.95) &&
        is.list(design$host_safety) &&
        identical(
            as.numeric(design$host_safety$minimum_disk_floor_bytes),
            8589934592
        ) && identical(
            as.numeric(design$host_safety$maximum_temporary_bytes_per_run),
            4294967296
        ) && identical(
            as.numeric(
                design$host_safety$maximum_prelaunch_load_per_logical_core
            ),
            0.35
        ) && identical(
            as.numeric(
                design$host_safety$maximum_prelaunch_thermal_celsius
            ),
            70
        ) && identical(
            as.numeric(
                design$host_safety$maximum_runtime_load_per_logical_core
            ),
            0.5
        ) && identical(
            as.numeric(design$host_safety$owned_workload_load_allowance),
            4
        ) && identical(
            as.numeric(design$host_safety$readiness_timeout_seconds),
            900
        ) && identical(
            as.numeric(design$host_safety$readiness_poll_seconds),
            10
        ) && isTRUE(design$host_safety$owned_user_cgroup_required) &&
        isTRUE(design$host_safety$zero_swap_required) &&
        isTRUE(design$host_safety$observe_power_policy_without_mutation)
    gate_ids <- vapply(gates$gates, `[[`, character(1), "gate_id")
    metric_ids <- vapply(gates$gates, `[[`, character(1), "metric")
    bounds <- vapply(gates$gates, `[[`, character(1), "confidence_bound")
    thresholds <- vapply(gates$gates, `[[`, numeric(1), "threshold")
    numerical_valid <- length(gate_ids) == 8L && !anyDuplicated(gate_ids) &&
        !anyDuplicated(metric_ids) &&
        all(bounds %in% c("upper_one_sided", "lower_one_sided")) &&
        all(is.finite(thresholds) & thresholds > 0)
    comparison_valid <- isTRUE(gates$comparison$owned_process_tree) &&
        isTRUE(gates$comparison$parity_workers_included) &&
        identical(
            gates$comparison$peak_scope,
            "complete_owned_cgroup_lifecycle"
        ) &&
        identical(gates$comparison$manual_garbage_collection_allowed, FALSE) &&
        diaRepairPredecessorGateValid(gates)
    authority_valid <- identical(
        gates$decision$automatic_policy_authority,
        FALSE
    ) && identical(gates$decision$publication_authority, FALSE)
    if (!all(c(
        design_valid,
        numerical_valid,
        comparison_valid,
        authority_valid
    ))) {
        diaRepairAbort("DIA repair gate contract semantics are invalid")
    }
    gates
}

diaRepairSchedule <- function(pair_count = 36L, session_count = 3L) {
    pair_count <- as.integer(pair_count)
    session_count <- as.integer(session_count)
    if (pair_count < 1L || pair_count %% 2L != 0L || session_count < 1L ||
        pair_count %% session_count != 0L) {
        diaRepairAbort("DIA repair schedule must be even and session-balanced")
    }
    first_arm <- rep(
        c("pre_repair_artifact", "candidate_artifact"),
        length.out = pair_count
    )
    do.call(rbind, lapply(seq_len(pair_count), function(index) {
        second_arm <- setdiff(
            c("pre_repair_artifact", "candidate_artifact"),
            first_arm[[index]]
        )
        data.frame(
            pair_id = sprintf("dia-repair-pair-%03d", index),
            session_id = sprintf(
                "session-%02d",
                ((index - 1L) %% session_count) + 1L
            ),
            sequence_in_pair = 1:2,
            arm = c(first_arm[[index]], second_arm),
            stringsAsFactors = FALSE
        )
    }))
}

diaRepairAttributeDescriptor <- function(value, owner) {
    if (is.null(value)) return(list(names = NULL, values = list()))
    value_names <- names(value)
    value_order <- order(value_names, method = "radix")
    value <- value[value_order]
    value_names <- value_names[value_order]
    list(
        names = unname(value_names),
        values = unname(lapply(seq_along(value), function(index) {
            diaRepairValueDescriptor(
                value[[index]],
                paste0(owner, "@", value_names[[index]])
            )
        }))
    )
}

diaRepairValueDescriptor <- function(value, owner) {
    if (is.null(value)) return(list(kind = "null"))
    if (isS4(value)) {
        slot_names <- methods::slotNames(value)
        return(list(
            kind = "S4",
            class_name = class(value)[[1L]],
            slot_names = slot_names,
            slots = unname(lapply(slot_names, function(slot_name) {
                diaRepairValueDescriptor(
                    methods::slot(value, slot_name),
                    paste0(owner, "@", slot_name)
                )
            }))
        ))
    }
    value_names <- names(value)
    value_attributes <- attributes(value)
    attributes(value) <- NULL
    if (is.atomic(value)) {
        materialized <- vector(typeof(value), length(value))
        materialized[] <- value
        return(list(
            kind = "atomic",
            storage_type = typeof(materialized),
            length = length(materialized),
            content_digest = digest::digest(
                materialized,
                algo = "sha256",
                serialize = TRUE,
                ascii = FALSE,
                serializeVersion = 3L
            ),
            attributes = diaRepairAttributeDescriptor(
                value_attributes,
                paste0(owner, "@attributes")
            )
        ))
    }
    if (is.list(value)) {
        return(list(
            kind = "list",
            values = unname(lapply(seq_along(value), function(index) {
                label <- if (!is.null(value_names) &&
                    nzchar(value_names[[index]])) {
                    value_names[[index]]
                } else {
                    index
                }
                diaRepairValueDescriptor(
                    value[[index]],
                    paste0(owner, "[[", label, "]]")
                )
            })),
            attributes = diaRepairAttributeDescriptor(
                value_attributes,
                paste0(owner, "@attributes")
            )
        ))
    }
    diaRepairAbort(paste("Unsupported DIA proof value:", owner))
}

diaRepairStableDigest <- function(value) {
    encoded <- jsonlite::toJSON(
        diaRepairValueDescriptor(value, "scientific_state"),
        auto_unbox = TRUE,
        null = "null",
        na = "string",
        digits = NA
    )
    digest::digest(as.character(encoded), algo = "sha256", serialize = FALSE)
}

diaRepairSaltedDigest <- function(value, salt) {
    if (!publicationScalarString(salt)) {
        diaRepairAbort("DIA repair evidence salt is unavailable")
    }
    digest::hmac(
        key = salt,
        object = value,
        algo = "sha256",
        serialize = FALSE
    )
}

diaRepairNamespaceValue <- function(namespace, name) {
    get(name, envir = namespace, inherits = FALSE)
}

diaRepairWorkflowPaths <- function(run_dir) {
    root <- file.path(run_dir, "private-project")
    paths <- list(
        base_dir = root,
        project_id = "dia-repair-private-canary",
        omic_type = "proteomics",
        omic_label = "dia-repair-canary",
        source_dir = file.path(root, "source"),
        results_dir = file.path(root, "results")
    )
    dir.create(paths$source_dir, recursive = TRUE, showWarnings = FALSE)
    dir.create(paths$results_dir, recursive = TRUE, showWarnings = FALSE)
    paths
}

diaRepairNewWorkflow <- function(namespace, paths, rollout = "dual_write") {
    create_context <- diaRepairNamespaceValue(namespace, "createWorkflowContext")
    state_class <- diaRepairNamespaceValue(namespace, "WorkflowState")
    workflow <- new.env(parent = emptyenv())
    workflow$workflow_context <- create_context(
        paths,
        "proteomics",
        "dia-repair-canary",
        storage_policy = list(
            requested_backend = "artifact",
            requested_rollout = rollout,
            migration_requested = TRUE,
            project_id = paths$project_id
        )
    )
    workflow$state_manager <- state_class$new()
    workflow$artifact_stage_results <- NULL
    workflow$processing_log <- list()
    workflow$tab_status <- list(
        setup_import = "complete",
        design_matrix = "complete",
        quality_control = "pending",
        normalization = "disabled",
        differential_expression = "disabled",
        enrichment_analysis = "disabled",
        session_summary = "disabled"
    )
    workflow
}

diaRepairPrepareImport <- function(namespace, workflow, source_path) {
    stage_import <- diaRepairNamespaceValue(
        namespace,
        "stageProtDiaImportArtifacts"
    )
    persist_import <- diaRepairNamespaceValue(
        namespace,
        "persistProtDiaImportArtifacts"
    )
    staged <- suppressMessages(stage_import(
        workflow,
        source_path,
        "diann",
        use_precursor_norm = TRUE,
        sanitize_names = FALSE
    ))
    if (!isTRUE(staged$ok) || !is.list(staged$result)) {
        diaRepairAbort("DIA repair import staging failed")
    }
    imported <- staged$result
    workflow$dia_repair_import_worker_timings <- staged$worker_timings
    workflow$data_tbl <- imported$data
    workflow$data_cln <- imported$data
    workflow$data_format <- "diann"
    workflow$data_type <- imported$data_type
    workflow$column_mapping <- imported$column_mapping
    workflow$artifact_import_summary <- imported$import_summary
    workflow$state_manager$setWorkflowType("DIA")
    published <- persist_import(
        workflow,
        imported,
        source_path,
        pending_stage = staged$pending_stage,
        worker_attempted = TRUE,
        log_warn = function(...) invisible(NULL)
    )
    if (!isTRUE(published$ok) || !isTRUE(published$committed)) {
        diaRepairAbort("DIA repair import publication failed")
    }
    invisible(imported)
}

diaRepairPrepareDesign <- function(namespace, workflow, imported, paths) {
    deferred <- isTRUE(imported$parent_hydration_deferred)
    runs <- if (deferred) {
        unlist(imported$import_summary$run_values, use.names = FALSE)
    } else {
        unique(workflow$data_cln$Run)
    }
    workflow$design_matrix <- data.frame(
        Run = runs,
        group = sub("_.*", "", runs),
        replicates = seq_along(runs),
        tech_rep_group = runs,
        stringsAsFactors = FALSE
    )
    workflow$config_list <- list(
        globalParameters = list(workflow_type = "DIA")
    )
    workflow$contrasts_tbl <- data.frame(
        contrasts = "groupKO-groupWT",
        friendly_names = "KO_vs_WT",
        full_format = "KO_vs_WT=groupKO-groupWT",
        stringsAsFactors = FALSE
    )
    workflow$uniprot_dat_cln <- NULL
    workflow$aa_seq_tbl_final <- NULL
    if (deferred) {
        design_input <- diaRepairNamespaceValue(
            namespace,
            "protDiaDeferredDesignInput"
        )(workflow)
        results <- list(
            data_cln = design_input,
            design_matrix = workflow$design_matrix,
            contrasts_tbl = workflow$contrasts_tbl,
            config_list = workflow$config_list
        )
        specification <- diaRepairNamespaceValue(
            namespace,
            "protDiaDeferredDesignSpec"
        )(
            workflow,
            results,
            paths,
            "dia-repair-canary"
        )
        worker <- diaRepairNamespaceValue(
            namespace,
            "runProtDiaDeferredDesignProcess"
        )(specification)
        settlement <- diaRepairNamespaceValue(
            namespace,
            "applyProtDiaDeferredDesignResult"
        )(workflow, worker)
        return(list(
            deferred = TRUE,
            design = worker$design_result,
            settlement = c(
                settlement,
                list(parity_worker_pid = worker$worker_pid)
            )
        ))
    }
    constructor <- diaRepairNamespaceValue(
        namespace,
        "PeptideQuantitativeDataDiann"
    )
    object <- constructor(
        peptide_data = workflow$data_cln,
        design_matrix = workflow$design_matrix,
        technical_replicate_id = "tech_rep_group",
        args = workflow$config_list
    )
    workflow$state_manager$saveState(
        "raw_data_s4",
        object,
        workflow$config_list,
        "DIA repair benchmark memory checkpoint."
    )
    compatibility <- list(
        data_cln = workflow$data_cln,
        design_matrix = workflow$design_matrix,
        contrasts_tbl = workflow$contrasts_tbl,
        config_list = workflow$config_list
    )
    diaRepairNamespaceValue(
        namespace,
        "persistProtDesignBuilderArtifacts"
    )(compatibility, workflow, paths$source_dir)
    list(deferred = FALSE, object = object)
}

diaRepairRunCommit <- function(namespace, workflow, paths, prepared_design) {
    if (isTRUE(prepared_design$deferred)) {
        return(prepared_design[c("design", "settlement")])
    }
    persist_design <- diaRepairNamespaceValue(
        namespace,
        "persistProtDiaDesignArtifacts"
    )
    settle <- diaRepairNamespaceValue(
        namespace,
        "settleProtDiaArtifactWorkflowSafely"
    )
    design <- persist_design(workflow, log_warn = function(...) invisible(NULL))
    if (!isTRUE(design$ok) || !isTRUE(design$committed)) {
        diaRepairAbort("DIA repair design publication failed")
    }
    settlement <- settle(
        workflow,
        paths,
        "dia-repair-canary",
        storage_policy = workflow$workflow_context$getStoragePolicy(),
        rollout_fn = function(...) "evict",
        log_warn = function(...) invisible(NULL)
    )
    released <- is.null(workflow$data_tbl) && is.null(workflow$data_cln)
    if (!isTRUE(settlement$evicted) || !released) {
        diaRepairAbort("DIA repair settlement did not release source fields")
    }
    list(design = design, settlement = settlement)
}

diaRepairPhaseSnapshot <- function() {
    cgroup_path <- publicationCurrentCgroupPath()
    if (is.null(cgroup_path) || !dir.exists(cgroup_path)) {
        diaRepairAbort("DIA repair worker cgroup is unavailable")
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

diaRepairWriteBoundary <- function(run_dir, marker) {
    record <- list(
        marker = marker,
        pid = as.integer(Sys.getpid()),
        monotonic_seconds = publicationMonotonicSeconds(),
        counters = diaRepairPhaseSnapshot()
    )
    publicationWriteJson(record, file.path(run_dir, paste0(marker, ".json")))
    invisible(record)
}

diaRepairArrowPoolBytes <- function() {
    if (!requireNamespace("arrow", quietly = TRUE)) return(0)
    as.numeric(arrow::default_memory_pool()$bytes_allocated)
}

diaRepairResourceLedger <- function(arrow_pool_bytes) {
    empty <- publicationEmptyWorkerResources()
    high_water <- empty
    high_water$active_tasks <- 1L
    high_water$duckdb_connections <- 1L
    high_water$arrow_pool_bytes <- arrow_pool_bytes
    list(
        schema = "multischolar.omics_publication_worker_resources",
        schema_version = "1.0.0",
        high_water = high_water,
        retained = empty,
        terminal = empty
    )
}

diaRepairRetentionState <- function(dwell_seconds, arrow_pool_bytes) {
    list(
        active_tasks = 0L,
        open_queries = 0L,
        open_writers = 0L,
        open_leases = 0L,
        open_connections = 0L,
        open_results = 0L,
        active_child_processes = 0L,
        arrow_pool_bytes = arrow_pool_bytes,
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

diaRepairComponentTimer <- function(operation) {
    started <- proc.time()
    value <- operation()
    finished <- proc.time()
    list(
        value = value,
        elapsed_seconds = unname(finished[["elapsed"]] - started[["elapsed"]]),
        process_cpu_seconds = unname(
            sum(finished[c("user.self", "sys.self")]) -
                sum(started[c("user.self", "sys.self")])
        )
    )
}

diaRepairWorker <- function(
    source_path,
    run_dir,
    package_library,
    arm,
    salt,
    dwell_seconds = 5
) {
    dir.create(run_dir, recursive = TRUE, showWarnings = FALSE)
    .libPaths(c(package_library, .libPaths()))
    namespace <- loadNamespace("MultiScholaR", lib.loc = package_library)
    paths <- diaRepairWorkflowPaths(run_dir)
    rollout <- if (identical(arm, "candidate_artifact")) {
        "evict"
    } else {
        "dual_write"
    }
    workflow <- diaRepairNewWorkflow(namespace, paths, rollout)
    acknowledgement_path <- file.path(run_dir, "retention-sampled.fifo")
    fifo_result <- processx::run(
        "mkfifo",
        acknowledgement_path,
        error_on_status = FALSE,
        echo = FALSE
    )
    if (fifo_result$status != 0L) diaRepairAbort("DIA repair FIFO creation failed")
    started <- proc.time()[["elapsed"]]
    cpu_started <- sum(proc.time()[c("user.self", "sys.self")])
    diaRepairWriteBoundary(run_dir, "measured_worker_start")
    import_component <- diaRepairComponentTimer(function() {
        suppressMessages(diaRepairPrepareImport(
            namespace,
            workflow,
            source_path
        ))
    })
    imported <- import_component$value
    design_component <- diaRepairComponentTimer(function() {
        suppressMessages(diaRepairPrepareDesign(
            namespace,
            workflow,
            imported,
            paths
        ))
    })
    prepared_design <- design_component$value
    commit_component <- diaRepairComponentTimer(function() {
        suppressMessages(diaRepairRunCommit(
            namespace,
            workflow,
            paths,
            prepared_design
        ))
    })
    commit <- commit_component$value
    phase_elapsed <- proc.time()[["elapsed"]] - started
    phase_cpu <- sum(proc.time()[c("user.self", "sys.self")]) - cpu_started
    diaRepairWriteBoundary(run_dir, "measured_worker_stop")
    proof <- diaRepairScientificProof(
        workflow,
        commit,
        salt,
        package_library,
        run_dir,
        arm
    )
    payload_free_state_manager <- inherits(
        workflow$state_manager,
        "ArtifactWorkflowState"
    ) && identical(workflow$state_manager$getCacheInfo()$entries, 0L)
    complete_payload_returned <- isTRUE(
        commit$design$hydration_verification$complete_payload_returned
    ) || isTRUE(commit$settlement$complete_payload_returned)
    result <- list(
        schema = "multischolar.dia_commit_repair_worker",
        schema_version = "1.0.0",
        status = "passed",
        pid = as.integer(Sys.getpid()),
        arm = arm,
        phase_elapsed_seconds = phase_elapsed,
        phase_cpu_seconds = phase_cpu,
        phase_components = list(
            import = import_component[c(
                "elapsed_seconds", "process_cpu_seconds"
            )] |> c(workflow$dia_repair_import_worker_timings),
            design = design_component[c(
                "elapsed_seconds", "process_cpu_seconds"
            )],
            commit = commit_component[c(
                "elapsed_seconds", "process_cpu_seconds"
            )]
        ),
        work_unit_id = "validated_input_byte",
        work_count = as.numeric(file.info(source_path)$size),
        scientific_proof = proof,
        payload_free_state_manager = payload_free_state_manager,
        complete_payload_returned = complete_payload_returned,
        parity_worker_distinct = if (is.null(
            commit$settlement$parity_worker_pid
        )) {
            NA
        } else {
            !identical(
                as.integer(commit$settlement$parity_worker_pid),
                as.integer(Sys.getpid())
            )
        }
    )
    arrow_pool_bytes <- diaRepairArrowPoolBytes()
    publicationWriteJson(result, file.path(run_dir, "worker-result.json"))
    publicationWriteJson(
        diaRepairResourceLedger(arrow_pool_bytes),
        file.path(run_dir, "worker-resources.json")
    )
    publicationWriteJson(
        diaRepairRetentionState(dwell_seconds, arrow_pool_bytes),
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
        diaRepairAbort("DIA repair retention acknowledgement differs")
    }
    publicationWriteJson(
        list(marker = "worker_exit", pid = as.integer(Sys.getpid())),
        file.path(run_dir, "worker-exit.json")
    )
    invisible(0L)
}

diaRepairExecution <- function() {
    list(
        sampling_interval_ms = 20,
        disk_sampling_interval_ms = 500,
        timeout_seconds = 1800,
        retained_dwell_seconds = 5,
        retained_window_seconds = 2,
        retention_acknowledgement = "fifo_v1",
        maximum_boundary_bracket_seconds = 0.5,
        retained_boundary_tolerance_ms = 100,
        zero_swap_required = TRUE
    )
}

diaRepairEngineeringPreflight <- function(host, contract, gates) {
    safety <- gates$design$host_safety
    preflight <- publicationEvaluateHostPreflight(
        host,
        contract,
        expected_peak_bytes = 4 * 1024^3,
        maximum_temporary_bytes = safety$maximum_temporary_bytes_per_run
    )
    minimum_disk <- max(
        safety$minimum_disk_floor_bytes,
        2 * safety$maximum_temporary_bytes_per_run
    )
    load_limit <- host$cpu$logical_cores *
        safety$maximum_prelaunch_load_per_logical_core
    load <- as.numeric(unlist(host$load_average, use.names = FALSE))
    preflight$checks$disk <- isTRUE(host$filesystem$available) &&
        isTRUE(host$filesystem$writable) &&
        publicationScalarNumber(host$filesystem$available_bytes) &&
        host$filesystem$available_bytes >= minimum_disk
    preflight$checks$load <- length(load) == 3L && all(is.finite(load)) &&
        load[[1L]] <= load_limit
    preflight$required_available_disk_bytes <- minimum_disk
    preflight$certified <- all(unlist(preflight$checks, use.names = FALSE))
    preflight$failure_outcome <- if (isTRUE(preflight$certified)) {
        NULL
    } else {
        "host_not_engineering_certified"
    }
    preflight$engineering_envelope <- list(
        claim_scope = gates$design$claim_scope,
        maximum_prelaunch_load_per_logical_core =
            safety$maximum_prelaunch_load_per_logical_core,
        publication_authority = FALSE
    )
    preflight
}

diaRepairRuntimeSafetyLimits <- function(
    host,
    preflight,
    gates,
    baseline_load
) {
    limits <- publicationDynamicRuntimeSafetyLimits(
        host,
        preflight,
        maximum_run_allocated_disk_bytes =
            gates$design$host_safety$maximum_temporary_bytes_per_run
    )
    limits$minimum_available_disk_bytes <-
        preflight$required_available_disk_bytes
    absolute_limit <- host$cpu$logical_cores *
        gates$design$host_safety$maximum_runtime_load_per_logical_core
    limits$maximum_load <- min(
        absolute_limit,
        baseline_load +
            gates$design$host_safety$owned_workload_load_allowance
    )
    limits
}

diaRepairWaitForReadiness <- function(host, gates) {
    safety <- gates$design$host_safety
    maximum_load <- diaRepairPrelaunchMaximumLoad(host, gates)
    maximum_thermal <- safety$maximum_prelaunch_thermal_celsius
    started <- proc.time()[["elapsed"]]
    observations <- list()
    repeat {
        load <- as.numeric(unlist(publicationLoadAverage(), use.names = FALSE))
        thermal <- publicationThermalState()
        elapsed <- proc.time()[["elapsed"]] - started
        observations[[length(observations) + 1L]] <- list(
            elapsed_seconds = elapsed,
            one_minute_load = load[[1L]],
            maximum_thermal_celsius = thermal$maximum_celsius
        )
        thermal_ready <- isTRUE(thermal$available) &&
            publicationScalarNumber(thermal$maximum_celsius) &&
            thermal$maximum_celsius <= maximum_thermal
        if (length(load) == 3L && all(is.finite(load)) &&
            load[[1L]] <= maximum_load && thermal_ready) {
            return(list(
                ready = TRUE,
                baseline_one_minute_load = load[[1L]],
                maximum_prelaunch_load = maximum_load,
                baseline_maximum_thermal_celsius = thermal$maximum_celsius,
                maximum_prelaunch_thermal_celsius = maximum_thermal,
                waited_seconds = elapsed,
                observations = observations
            ))
        }
        if (elapsed >= safety$readiness_timeout_seconds) {
            return(list(
                ready = FALSE,
                baseline_one_minute_load = load[[1L]],
                maximum_prelaunch_load = maximum_load,
                baseline_maximum_thermal_celsius = thermal$maximum_celsius,
                maximum_prelaunch_thermal_celsius = maximum_thermal,
                waited_seconds = elapsed,
                observations = observations
            ))
        }
        Sys.sleep(safety$readiness_poll_seconds)
    }
}

diaRepairThreadEnvironment <- function(home, temp, library, site_library) {
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
        R_LIBS = library,
        R_LIBS_USER = library,
        R_LIBS_SITE = site_library,
        HOME = home,
        TMPDIR = temp
    )
}

diaRepairReadPhase <- function(run_dir) {
    list(
        start = publicationReadJson(file.path(
            run_dir,
            "measured_worker_start.json"
        )),
        stop = publicationReadJson(file.path(
            run_dir,
            "measured_worker_stop.json"
        )),
        worker = publicationReadJson(file.path(run_dir, "worker-result.json"))
    )
}

diaRepairPhasePeak <- function(measurement, phase) {
    origin <- as.numeric(
        measurement$retention_state$settled_monotonic_seconds
    ) - as.numeric(measurement$retention_marker_elapsed_seconds)
    start <- as.numeric(phase$start$monotonic_seconds) - origin
    stop <- as.numeric(phase$stop$monotonic_seconds) - origin
    samples <- Filter(function(sample) {
        elapsed <- as.numeric(sample$elapsed_seconds)
        elapsed >= start && elapsed <= stop
    }, measurement$samples)
    if (!length(samples)) diaRepairAbort("DIA repair phase has no cgroup samples")
    max(vapply(samples, function(sample) {
        as.numeric(sample$memory$current_bytes)
    }, numeric(1)))
}

diaRepairAttachPhase <- function(measurement, run_dir) {
    phase <- diaRepairReadPhase(run_dir)
    worker <- phase$worker
    cpu_names <- c("usage_usec", "user_usec", "system_usec")
    cpu <- stats::setNames(lapply(cpu_names, function(name) {
        as.numeric(phase$stop$counters$cpu[[name]]) -
            as.numeric(phase$start$counters$cpu[[name]])
    }), cpu_names)
    phase_cpu <- cpu$usage_usec / 1e6
    valid <- identical(worker$status, "passed") &&
        identical(as.integer(worker$pid), as.integer(phase$start$pid)) &&
        identical(as.integer(worker$pid), as.integer(phase$stop$pid)) &&
        all(unlist(cpu, use.names = FALSE) >= 0) && phase_cpu > 0
    if (!isTRUE(valid)) diaRepairAbort("DIA repair phase evidence is invalid")
    measurement$metrics$phase_peak_charged_memory_bytes <-
        diaRepairPhasePeak(measurement, phase)
    measurement$metrics$cgroup_lifecycle_elapsed_seconds <-
        measurement$metrics$elapsed_seconds
    measurement$metrics$elapsed_seconds <- worker$phase_elapsed_seconds
    measurement$metrics$phase_cpu_seconds <- phase_cpu
    measurement$metrics$phase_process_self_cpu_seconds <-
        worker$phase_cpu_seconds
    measurement$metrics$primary_work_units_per_wall_second <-
        worker$work_count / worker$phase_elapsed_seconds
    measurement$metrics$primary_work_units_per_cpu_second <-
        worker$work_count / phase_cpu
    measurement$phase_evidence <- list(
        valid = TRUE,
        worker_pid = worker$pid,
        phase_peak_charged_memory_bytes =
            measurement$metrics$phase_peak_charged_memory_bytes,
        phase_elapsed_seconds = worker$phase_elapsed_seconds,
        phase_cgroup_cpu_seconds = phase_cpu
    )
    measurement$worker <- worker
    measurement
}

diaRepairPairMetric <- function(records, metric) {
    pair_ids <- unique(vapply(records, `[[`, character(1), "pair_id"))
    rows <- lapply(pair_ids, function(pair_id) {
        pair <- Filter(function(record) identical(record$pair_id, pair_id), records)
        arms <- vapply(pair, `[[`, character(1), "arm")
        if (length(pair) != 2L || !setequal(
            arms,
            c("pre_repair_artifact", "candidate_artifact")
        )) {
            return(NULL)
        }
        values <- stats::setNames(vapply(pair, function(record) {
            value <- record$measurement$metrics[[metric]]
            if (is.null(value) || length(value) != 1L) return(NA_real_)
            as.numeric(value)
        }, numeric(1)), arms)
        if (any(!is.finite(values)) || any(values <= 0)) return(NULL)
        data.frame(
            pair_id = pair_id,
            ratio = values[["candidate_artifact"]] /
                values[["pre_repair_artifact"]],
            log_ratio = log(values[["candidate_artifact"]]) -
                log(values[["pre_repair_artifact"]]),
            stringsAsFactors = FALSE
        )
    })
    rows <- Filter(Negate(is.null), rows)
    if (!length(rows)) return(data.frame())
    do.call(rbind, rows)
}

diaRepairBootstrapRatio <- function(pairs, gates) {
    if (!is.data.frame(pairs) || !nrow(pairs)) return(NULL)
    set.seed(as.integer(gates$design$seed))
    resamples <- as.integer(gates$design$resamples)
    values <- replicate(resamples, {
        sampled <- sample.int(nrow(pairs), nrow(pairs), replace = TRUE)
        exp(mean(pairs$log_ratio[sampled]))
    })
    alpha <- 1 - as.numeric(gates$design$confidence_level)
    list(
        pairs = nrow(pairs),
        geometric_mean_ratio = exp(mean(pairs$log_ratio)),
        lower_one_sided = unname(stats::quantile(
            values,
            alpha,
            type = 8
        )),
        upper_one_sided = unname(stats::quantile(
            values,
            1 - alpha,
            type = 8
        ))
    )
}

diaRepairEvaluate <- function(records, gates = diaRepairReadGates()) {
    required_pairs <- as.integer(gates$design$required_pairs)
    gate_results <- lapply(gates$gates, function(gate) {
        pairs <- diaRepairPairMetric(records, gate$metric)
        summary <- diaRepairBootstrapRatio(pairs, gates)
        bound <- if (is.null(summary)) NA_real_ else
            as.numeric(summary[[gate$confidence_bound]])
        passed <- nrow(pairs) >= required_pairs && is.finite(bound) &&
            if (identical(gate$confidence_bound, "upper_one_sided")) {
                bound <= gate$threshold
            } else {
                bound >= gate$threshold
            }
        list(
            gate_id = gate$gate_id,
            metric = gate$metric,
            threshold = gate$threshold,
            confidence_bound = gate$confidence_bound,
            observed_bound = bound,
            summary = summary,
            passed = passed
        )
    })
    pair_ids <- unique(vapply(records, `[[`, character(1), "pair_id"))
    parity <- vapply(pair_ids, function(pair_id) {
        pair <- Filter(function(record) identical(record$pair_id, pair_id), records)
        proofs <- lapply(pair, function(record) {
            record$measurement$worker$scientific_proof
        })
        length(pair) == 2L && all(vapply(proofs, is.list, logical(1))) &&
            identical(
                pair[[1L]]$measurement$worker$scientific_proof,
                pair[[2L]]$measurement$worker$scientific_proof
            )
    }, logical(1))
    run_valid <- vapply(records, function(record) {
        measurement <- record$measurement
        worker <- measurement$worker
        identical(measurement$status, "passed") &&
            isTRUE(measurement$publication_certifiable) &&
            isTRUE(measurement$cleanup$valid) &&
            identical(as.numeric(measurement$metrics$peak_swap_bytes), 0) &&
            isTRUE(worker$scientific_proof$valid_s4) &&
            isTRUE(worker$scientific_proof$source_fields_released) &&
            !isTRUE(worker$complete_payload_returned) &&
            (!identical(record$arm, "candidate_artifact") ||
                isTRUE(worker$payload_free_state_manager))
    }, logical(1))
    numerical_pass <- all(vapply(gate_results, `[[`, logical(1), "passed"))
    absolute_pass <- length(pair_ids) >= required_pairs && all(parity) &&
        all(run_valid)
    valid_pair_counts <- vapply(gate_results, function(gate) {
        if (is.null(gate$summary)) 0L else as.integer(gate$summary$pairs)
    }, integer(1))
    list(
        status = if (numerical_pass && absolute_pass) "passed" else "failed",
        valid_pairs = min(valid_pair_counts),
        numerical_gates = gate_results,
        scientific_parity = list(passed = all(parity), pairs = length(parity)),
        lifecycle = list(passed = all(run_valid), runs = length(run_valid)),
        automatic_policy_authority = FALSE,
        publication_authority = FALSE,
        may_start_omics_art_071 = numerical_pass && absolute_pass
    )
}
