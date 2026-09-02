#!/usr/bin/env Rscript

metabPilotScriptPath <- function() {
    file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
    if (length(file_arg)) {
        return(normalizePath(sub("^--file=", "", file_arg[[1L]])))
    }
    normalizePath(file.path(
        getwd(),
        "tools",
        "profiling",
        "run_metab_engineering_pilot.R"
    ))
}

.METAB_PILOT_ROOT <- normalizePath(
    file.path(dirname(metabPilotScriptPath()), "..", "..")
)

for (source_file in c(
    "omics_publication_protocol.R",
    "omics_publication_host_safety.R",
    "omics_publication_linux_resources.R",
    "omics_publication_retained_resources.R"
)) {
    source(file.path(.METAB_PILOT_ROOT, "tools", "profiling", source_file))
}

`%||%` <- function(left, right) {
    if (is.null(left)) right else left
}

metabPilotArgs <- function(argv) {
    values <- list(
        mode = "campaign",
        arm = NULL,
        source_a = NULL,
        source_b = NULL,
        package_library = NULL,
        old_library = NULL,
        candidate_library = NULL,
        run_dir = NULL,
        output_root = NULL,
        result = NULL,
        repetitions = 3L,
        dwell_seconds = 5
    )
    index <- 1L
    while (index <= length(argv)) {
        token <- argv[[index]]
        key <- gsub("-", "_", sub("^--", "", token), fixed = TRUE)
        if (!startsWith(token, "--") || !key %in% names(values) ||
            index == length(argv)) {
            stop("metabolomics pilot arguments are invalid", call. = FALSE)
        }
        value <- argv[[index + 1L]]
        values[[key]] <- if (identical(key, "repetitions")) {
            as.integer(value)
        } else if (identical(key, "dwell_seconds")) {
            as.numeric(value)
        } else {
            value
        }
        index <- index + 2L
    }
    values
}

metabPilotRequire <- function(args, fields) {
    missing <- fields[vapply(fields, \(field) {
        is.null(args[[field]]) || !nzchar(args[[field]])
    }, logical(1))]
    if (length(missing)) {
        stop(paste("metabolomics pilot option is missing:", missing[[1L]]),
            call. = FALSE)
    }
    invisible(args)
}

metabPilotValue <- function(namespace, name) {
    get(name, envir = namespace, inherits = FALSE)
}

metabPilotPaths <- function(run_dir) {
    root <- file.path(run_dir, "project")
    paths <- list(
        base_dir = root,
        project_id = "metab-pilot",
        omic_type = "metabolomics",
        omic_label = "metab-pilot",
        source_dir = file.path(root, "source"),
        results_dir = file.path(root, "results")
    )
    dir.create(paths$source_dir, recursive = TRUE, showWarnings = FALSE)
    dir.create(paths$results_dir, recursive = TRUE, showWarnings = FALSE)
    paths
}

metabPilotWorkflow <- function(namespace, paths) {
    workflow <- new.env(parent = emptyenv())
    workflow$workflow_context <- metabPilotValue(namespace, "createWorkflowContext")(
        paths,
        "metabolomics",
        "metab-pilot",
        storage_policy = list(
            requested_backend = "artifact",
            requested_rollout = "evict",
            migration_requested = TRUE,
            project_id = paths$project_id
        )
    )
    workflow$state_manager <- metabPilotValue(namespace, "WorkflowState")$new()
    workflow$artifact_stage_results <- list()
    workflow$processing_log <- list()
    workflow$tab_status <- list(
        setup_import = "complete",
        design_matrix = "complete",
        quality_control = "disabled",
        normalization = "disabled"
    )
    workflow
}

metabPilotMapping <- function() {
    list(
        metabolite_id_col = "Feature.Name",
        annotation_col = "Annotation",
        sample_columns = c(
            sprintf("CTRL_%02d", 1:6),
            sprintf("TREAT_%02d", 1:6)
        ),
        is_pattern = NA_character_
    )
}

metabPilotMemoryPayload <- function(namespace, args, mapping) {
    importer <- metabPilotValue(namespace, "importMetabMSDIALData")
    assays <- list(
        LCMS_Pos = suppressMessages(importer(args$source_a))$data,
        GCMS = suppressMessages(importer(args$source_b))$data
    )
    payload <- list(
        assayList = assays,
        sampleCols = mapping$sample_columns,
        sampleNamesSanitized = FALSE,
        dataFormat = "custom",
        columnMapping = mapping,
        processingLog = list(
            evidence_class = "generated_operational_heavy",
            n_assays = 2L
        ),
        sourceFiles = list(
            assay1 = args$source_a,
            assay2 = args$source_b
        ),
        assaySourceFiles = list(
            LCMS_Pos = args$source_a,
            GCMS = args$source_b
        )
    )
    list(
        payload = metabPilotValue(
            namespace,
            "coerceMetabImportWorkflowPayloadSamples"
        )(payload),
        selection_payload = assays
    )
}

metabPilotImport <- function(namespace, workflow, args) {
    mapping <- metabPilotMapping()
    candidate <- identical(args$arm, "candidate_artifact")
    staged <- NULL
    memory_import <- metabPilotMemoryPayload(namespace, args, mapping)
    payload <- memory_import$payload
    if (candidate) {
        memory_import$selection_payload <- NULL
        metabPilotValue(namespace, "artifactReleaseTransientMemory")(full = TRUE)
    }
    metabPilotValue(namespace, "applyMetabImportWorkflowPayload")(
        workflow,
        payload,
        logInfo = \(...) invisible(NULL)
    )
    persist_import <- metabPilotValue(namespace, "persistMetabImportArtifacts")
    persist_args <- list(
        workflow_data = workflow,
        workflow_payload = payload,
        log_warn = \(...) invisible(NULL)
    )
    if ("pending_stage" %in% names(formals(persist_import))) {
        persist_args$pending_stage <- NULL
        persist_args$worker_attempted <- FALSE
    }
    persisted <- do.call(persist_import, persist_args)
    if (!isTRUE(persisted$ok) || !isTRUE(persisted$committed)) {
        stop("metabolomics pilot import did not commit", call. = FALSE)
    }
    list(
        payload = payload,
        persisted = persisted,
        staged = staged,
        selection_payload = if (candidate) NULL else {
            memory_import$selection_payload
        }
    )
}

metabPilotDesign <- function(namespace, workflow, payload) {
    mapping <- payload$columnMapping
    metadata <- c(mapping$metabolite_id_col, mapping$annotation_col)
    workflow$data_cln <- lapply(payload$assayList, \(assay) {
        assay[c(metadata, mapping$sample_columns)]
    })
    samples <- mapping$sample_columns
    groups <- rep(c("CTRL", "TREAT"), each = length(samples) / 2L)
    replicates <- ave(seq_along(samples), groups, FUN = seq_along)
    workflow$design_matrix <- data.frame(
        Run = samples,
        group = groups,
        batch = rep(c("B1", "B2"), length.out = length(samples)),
        replicates = as.integer(replicates),
        tech_rep_group = paste(groups, replicates, sep = "_"),
        stringsAsFactors = FALSE
    )
    workflow$contrasts_tbl <- data.frame(
        contrasts = "groupTREAT-groupCTRL",
        friendly_names = "TREAT_vs_CTRL",
        full_format = "TREAT_vs_CTRL=groupTREAT-groupCTRL",
        stringsAsFactors = FALSE
    )
    workflow$config_list <- list(
        globalParameters = list(workflow_type = "metabolomics_standard"),
        deAnalysisParameters = list(formula_string = "~ 0 + group")
    )
    object <- metabPilotValue(namespace, "createMetaboliteAssayData")(
        workflow$data_cln,
        workflow$design_matrix,
        metabolite_id_column = mapping$metabolite_id_col,
        annotation_id_column = mapping$annotation_col,
        sample_id = "Run",
        group_id = "group",
        technical_replicate_id = "tech_rep_group",
        args = workflow$config_list
    )
    workflow$state_manager$setWorkflowType("metabolomics_standard")
    workflow$state_manager$saveState(
        "metab_raw_data_s4",
        object,
        workflow$config_list,
        "metabolomics engineering pilot"
    )
    persisted <- metabPilotValue(namespace, "persistMetabDesignArtifacts")(
        workflow,
        log_warn = \(...) invisible(NULL)
    )
    if (!isTRUE(persisted$ok) || !isTRUE(persisted$committed)) {
        stop("metabolomics pilot design did not commit", call. = FALSE)
    }
    list(object = object, persisted = persisted)
}

metabPilotStateDigest <- function(object, digest_fn) {
    slots <- methods::slotNames(object)
    digests <- vapply(slots, \(slot_name) {
        digest_fn(methods::slot(object, slot_name))
    }, character(1))
    list(
        digest = digest_fn(list(
            class = class(object),
            slots = slots,
            digests = digests
        )),
        slot_digests = as.list(digests)
    )
}

metabPilotWorker <- function(args) {
    metabPilotRequire(args, c(
        "arm", "source_a", "source_b", "package_library", "run_dir"
    ))
    .libPaths(c(args$package_library, .libPaths()))
    namespace <- loadNamespace("MultiScholaR", lib.loc = args$package_library)
    acknowledgement_path <- file.path(args$run_dir, "retention-sampled.fifo")
    dir.create(args$run_dir, recursive = TRUE, showWarnings = FALSE)
    fifo_result <- processx::run(
        "mkfifo",
        acknowledgement_path,
        error_on_status = FALSE,
        echo = FALSE
    )
    if (fifo_result$status != 0L) stop("metabolomics pilot FIFO creation failed")
    paths <- metabPilotPaths(args$run_dir)
    workflow <- metabPilotWorkflow(namespace, paths)
    imported <- metabPilotImport(namespace, workflow, args)
    design <- metabPilotDesign(namespace, workflow, imported$payload)
    state_evidence <- metabPilotStateDigest(
        design$object,
        metabPilotValue(namespace, "artifactSemanticDigest")
    )
    settlement <- design$persisted$settlement
    if (!is.list(settlement) || !isTRUE(settlement$evicted)) {
        settlement <- metabPilotValue(
            namespace,
            "settleMetabArtifactWorkflowSafely"
        )(
            workflow,
            rollout_fn = \(...) "evict",
            log_warn = \(...) invisible(NULL)
        )
    }
    if (!isTRUE(settlement$evicted)) {
        stop("metabolomics pilot settlement failed", call. = FALSE)
    }
    imported$payload <- NULL
    imported$staged <- NULL
    design$object <- NULL
    result <- list(
        status = "passed",
        arm = args$arm,
        state_digest = state_evidence$digest,
        state_slot_digests = state_evidence$slot_digests,
        source_fields_released = is.null(workflow$data_tbl) &&
            is.null(workflow$data_cln),
        cache_entries = workflow$state_manager$getCacheInfo()$entries,
        import_process_evidence = imported$persisted$process_evidence,
        state_manager_replaced = isTRUE(settlement$state_manager_replaced),
        complete_payload_returned = isTRUE(settlement$complete_payload_returned),
        historical_selection_payload_retained =
            !is.null(imported$selection_payload)
    )
    result_path <- file.path(args$run_dir, "final", "worker-result.json")
    dir.create(dirname(result_path), recursive = TRUE, showWarnings = FALSE)
    jsonlite::write_json(
        result,
        result_path,
        auto_unbox = TRUE,
        pretty = TRUE,
        null = "null"
    )
    arrow_pool_bytes <- as.numeric(arrow::default_memory_pool()$bytes_allocated)
    resources <- publicationEmptyWorkerResources()
    resources$arrow_pool_bytes <- arrow_pool_bytes
    publicationWriteJson(list(
        schema = "multischolar.omics_publication_worker_resources",
        schema_version = "1.0.0",
        high_water = resources,
        retained = resources,
        terminal = publicationEmptyWorkerResources()
    ), file.path(args$run_dir, "worker-resources.json"))
    retention_state <- c(list(
        active_tasks = 0,
        open_queries = 0,
        open_writers = 0,
        open_leases = 0,
        open_connections = 0,
        open_results = 0,
        active_child_processes = 0,
        arrow_pool_bytes = arrow_pool_bytes,
        duckdb_memory_bytes = 0,
        duckdb_spill_bytes = 0,
        duckdb_prepared_statements = 0,
        temporary_paths = 0,
        cache_entries = 0,
        observers = 0,
        native_resources = 0,
        retained_dwell_seconds = args$dwell_seconds,
        retention_acknowledgement = "fifo_v1",
        settled_monotonic_seconds = publicationMonotonicSeconds()
    ))
    publicationWriteJson(
        retention_state,
        file.path(args$run_dir, "retention-state.json")
    )
    writeLines("ready", file.path(args$run_dir, "retention-ready"))
    acknowledgement <- fifo(acknowledgement_path, open = "rb", blocking = TRUE)
    on.exit(close(acknowledgement), add = TRUE)
    token <- readBin(acknowledgement, what = "raw", n = 1L)
    if (!identical(token, as.raw(1L))) {
        stop("metabolomics pilot acknowledgement differs", call. = FALSE)
    }
    invisible(0L)
}

metabPilotCategories <- function() {
    list(
        list(category = "diagnostics", pattern = "(^|/)diagnostics(/|$)"),
        list(category = "staging_snapshot", pattern = "(^|/)project(/|$)"),
        list(category = "duckdb_spill", pattern = "(^|/)duckdb-spill(/|$)"),
        list(category = "committed", pattern = "(^|/)project(/|$)"),
        list(category = "final", pattern = "(^|/)final(/|$)"),
        list(category = "harness", pattern = ".*")
    )
}

metabPilotMeasure <- function(args, arm, repetition) {
    library <- if (identical(arm, "pre_repair_artifact")) {
        args$old_library
    } else {
        args$candidate_library
    }
    run_dir <- file.path(
        args$output_root,
        sprintf("pair-%02d", repetition),
        arm
    )
    worker_args <- c(
        "--vanilla", metabPilotScriptPath(),
        "--mode", "worker",
        "--arm", arm,
        "--source-a", normalizePath(args$source_a),
        "--source-b", normalizePath(args$source_b),
        "--package-library", normalizePath(library),
        "--run-dir", normalizePath(run_dir, mustWork = FALSE),
        "--dwell-seconds", as.character(args$dwell_seconds)
    )
    measurement <- publicationMeasureCgroupSubprocess(
        command = file.path(R.home("bin"), "Rscript"),
        command_args = worker_args,
        run_dir = run_dir,
        execution = list(
            sampling_interval_ms = 20,
            timeout_seconds = 900,
            retained_dwell_seconds = args$dwell_seconds,
            retained_window_seconds = 2,
            maximum_boundary_bracket_seconds = 0.5,
            retention_acknowledgement = "fifo_v1",
            maximum_idle_cpu_activity_seconds = 0.001,
            zero_swap_required = TRUE
        ),
        unit_name = publicationSystemdUnitName("multischolar-metab-pilot")
    )
    worker_path <- file.path(run_dir, "final", "worker-result.json")
    worker <- if (file.exists(worker_path)) {
        jsonlite::read_json(worker_path, simplifyVector = FALSE)
    } else {
        NULL
    }
    list(
        pair = repetition,
        arm = arm,
        measurement = measurement,
        worker = worker
    )
}

metabPilotRatio <- function(records, metric) {
    pairs <- split(records, vapply(records, `[[`, integer(1), "pair"))
    ratios <- vapply(pairs, \(pair) {
        values <- stats::setNames(vapply(pair, \(record) {
            as.numeric(record$measurement$metrics[[metric]])
        }, numeric(1)), vapply(pair, `[[`, character(1), "arm"))
        values[["candidate_artifact"]] / values[["pre_repair_artifact"]]
    }, numeric(1))
    list(
        metric = metric,
        ratios = as.list(ratios),
        median_ratio = unname(stats::median(ratios))
    )
}

metabPilotCampaign <- function(args) {
    metabPilotRequire(args, c(
        "source_a", "source_b", "old_library", "candidate_library",
        "output_root", "result"
    ))
    if (dir.exists(args$output_root) || file.exists(args$result)) {
        stop("metabolomics pilot output already exists", call. = FALSE)
    }
    dir.create(args$output_root, recursive = TRUE, showWarnings = FALSE)
    records <- list()
    for (pair in seq_len(args$repetitions)) {
        order <- if (pair %% 2L) {
            c("pre_repair_artifact", "candidate_artifact")
        } else {
            c("candidate_artifact", "pre_repair_artifact")
        }
        for (arm in order) {
            records[[length(records) + 1L]] <- metabPilotMeasure(
                args,
                arm,
                pair
            )
        }
    }
    valid <- all(vapply(records, \(record) {
        identical(as.integer(record$measurement$exit_status), 0L) &&
            isTRUE(record$measurement$cgroup_observed) &&
            isTRUE(record$measurement$metrics_valid) &&
            isTRUE(record$measurement$cleanup$valid) &&
            isTRUE(record$measurement$retention_acknowledged) &&
            !isTRUE(record$measurement$safety_aborted) &&
            identical(record$worker$status, "passed") &&
            isTRUE(record$worker$source_fields_released)
    }, logical(1)))
    digests <- unique(vapply(records, \(record) {
        if (is.list(record$worker)) {
            record$worker$state_digest %||% NA_character_
        } else {
            NA_character_
        }
    }, character(1)))
    digests <- digests[!is.na(digests)]
    result <- list(
        schema = "multischolar.metabolomics_engineering_pilot",
        schema_version = "1.0.0",
        owner_ticket_id = "OMICS-ART-075",
        status = if (valid && length(digests) == 1L) "passed" else "failed",
        source_bytes = sum(file.info(c(args$source_a, args$source_b))$size),
        repetitions = args$repetitions,
        records = records,
        ratios = if (valid) lapply(c(
            "peak_charged_memory_bytes",
            "retained_charged_memory_bytes",
            "elapsed_seconds",
            "peak_logical_disk_bytes",
            "final_logical_disk_bytes"
        ), \(metric) metabPilotRatio(records, metric)) else list(),
        scientific_parity = length(digests) == 1L,
        publication_certifiable = all(vapply(records, \(record) {
            isTRUE(record$measurement$publication_certifiable)
        }, logical(1))),
        promotion_authority = FALSE,
        publication_authority = FALSE
    )
    jsonlite::write_json(
        result,
        args$result,
        auto_unbox = TRUE,
        pretty = TRUE,
        null = "null",
        digits = NA
    )
    invisible(if (identical(result$status, "passed")) 0L else 1L)
}

if (identical(environment(), globalenv()) && !interactive()) {
    args <- metabPilotArgs(commandArgs(trailingOnly = TRUE))
    status <- tryCatch(
        if (identical(args$mode, "worker")) {
            metabPilotWorker(args)
        } else {
            metabPilotCampaign(args)
        },
        error = \(error) {
            message(conditionMessage(error))
            1L
        }
    )
    quit(status = as.integer(status), save = "no")
}
