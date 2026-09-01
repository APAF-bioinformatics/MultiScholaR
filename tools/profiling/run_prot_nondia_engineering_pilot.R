#!/usr/bin/env Rscript

protNonDiaPilotScriptPath <- function() {
    file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
    if (length(file_arg)) {
        return(normalizePath(sub("^--file=", "", file_arg[[1L]])))
    }
    normalizePath(file.path(
        getwd(),
        "tools",
        "profiling",
        "run_prot_nondia_engineering_pilot.R"
    ))
}

.PROT_NONDIA_PILOT_ROOT <- normalizePath(
    file.path(dirname(protNonDiaPilotScriptPath()), "..", "..")
)

for (source_file in c(
    "omics_publication_protocol.R",
    "omics_publication_host_safety.R",
    "omics_publication_linux_resources.R",
    "omics_publication_retained_resources.R"
)) {
    source(file.path(
        .PROT_NONDIA_PILOT_ROOT,
        "tools",
        "profiling",
        source_file
    ))
}

protNonDiaPilotArgs <- function(argv) {
    values <- list(
        mode = "campaign",
        arm = NULL,
        format = NULL,
        source = NULL,
        package_library = NULL,
        old_library = NULL,
        candidate_library = NULL,
        run_dir = NULL,
        output_root = NULL,
        result = NULL,
        repetitions = 2L,
        dwell_seconds = 5
    )
    index <- 1L
    while (index <= length(argv)) {
        token <- argv[[index]]
        key <- gsub("-", "_", sub("^--", "", token), fixed = TRUE)
        if (!startsWith(token, "--") || !key %in% names(values) ||
            index == length(argv)) {
            stop("non-DIA pilot arguments are invalid", call. = FALSE)
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

protNonDiaPilotRequire <- function(args, fields) {
    missing <- fields[vapply(fields, \(field) {
        is.null(args[[field]]) || !nzchar(args[[field]])
    }, logical(1))]
    if (length(missing)) {
        stop(paste("non-DIA pilot option is missing:", missing[[1L]]), call. = FALSE)
    }
    invisible(args)
}

protNonDiaPilotNamespaceValue <- function(namespace, name) {
    get(name, envir = namespace, inherits = FALSE)
}

protNonDiaPilotImporter <- function(namespace, format) {
    name <- switch(format,
        maxquant = "importMaxQuantData",
        fragpipe = "importFragPipeData",
        stop("non-DIA pilot format is unsupported", call. = FALSE)
    )
    protNonDiaPilotNamespaceValue(namespace, name)
}

protNonDiaPilotPaths <- function(run_dir, format) {
    root <- file.path(run_dir, "project")
    paths <- list(
        base_dir = root,
        project_id = paste0("nondia-pilot-", format),
        omic_type = "proteomics",
        omic_label = "nondia-pilot",
        source_dir = file.path(root, "source"),
        results_dir = file.path(root, "results")
    )
    dir.create(paths$source_dir, recursive = TRUE, showWarnings = FALSE)
    dir.create(paths$results_dir, recursive = TRUE, showWarnings = FALSE)
    paths
}

protNonDiaPilotWorkflow <- function(namespace, paths) {
    workflow <- new.env(parent = emptyenv())
    workflow$workflow_context <- protNonDiaPilotNamespaceValue(
        namespace,
        "createWorkflowContext"
    )(
        paths,
        "proteomics",
        "nondia-pilot",
        storage_policy = list(
            requested_backend = "artifact",
            requested_rollout = "dual_write",
            migration_requested = TRUE,
            project_id = paths$project_id
        )
    )
    workflow$state_manager <- protNonDiaPilotNamespaceValue(
        namespace,
        "WorkflowState"
    )$new()
    workflow$artifact_stage_results <- NULL
    workflow$tab_status <- list(
        setup_import = "complete",
        design_matrix = "complete"
    )
    workflow
}

protNonDiaPilotImport <- function(namespace, workflow, args) {
    candidate <- identical(args$arm, "candidate_artifact")
    staged <- if (candidate) {
        suppressMessages(protNonDiaPilotNamespaceValue(
            namespace,
            "stageProtNonDiaImportArtifacts"
        )(
            workflow,
            args$source,
            args$format
        ))
    } else {
        NULL
    }
    imported <- if (is.list(staged) && isTRUE(staged$ok)) {
        staged$result
    } else {
        suppressMessages(protNonDiaPilotImporter(namespace, args$format)(
            args$source
        ))
    }
    workflow$data_tbl <- imported$data
    workflow$data_cln <- imported$data
    workflow$data_format <- args$format
    workflow$data_type <- imported$data_type
    workflow$column_mapping <- imported$column_mapping
    workflow$state_manager$setWorkflowType("LFQ")
    persist_import <- protNonDiaPilotNamespaceValue(
        namespace,
        "persistProtNonDiaImportArtifacts"
    )
    persisted <- if (candidate) {
        persist_import(
            workflow,
            imported,
            args$source,
            pending_stage = staged$pending_stage,
            worker_attempted = TRUE
        )
    } else {
        persist_import(workflow, imported, args$source)
    }
    if (!isTRUE(persisted$ok) || !isTRUE(persisted$committed)) {
        stop("non-DIA pilot import did not commit", call. = FALSE)
    }
    list(imported = imported, persisted = persisted, staged = staged)
}

protNonDiaPilotDesign <- function(namespace, workflow) {
    mapping <- workflow$column_mapping
    runs <- unique(as.character(workflow$data_cln[[mapping$run_col]]))
    group <- rep(c("A", "B"), length.out = length(runs))
    workflow$design_matrix <- data.frame(
        Run = runs,
        group = group,
        replicates = ave(seq_along(runs), group, FUN = seq_along),
        stringsAsFactors = FALSE
    )
    workflow$config_list <- list(
        globalParameters = list(workflow_type = "LFQ")
    )
    workflow$contrasts_tbl <- data.frame(
        contrasts = "groupB-groupA",
        friendly_names = "B_vs_A",
        full_format = "B_vs_A=groupB-groupA",
        stringsAsFactors = FALSE
    )
    proteins <- unique(as.character(
        workflow$data_cln[[mapping$protein_col]]
    ))
    workflow$uniprot_dat_cln <- data.frame(
        Protein.Ids = proteins,
        Gene = paste0("GENE", seq_along(proteins)),
        stringsAsFactors = FALSE
    )
    workflow$aa_seq_tbl_final <- data.frame(
        accession = proteins,
        sequence = rep("PEPTIDE", length(proteins)),
        stringsAsFactors = FALSE
    )
    suppressMessages(protNonDiaPilotNamespaceValue(
        namespace,
        "buildProtDesignStateCheckpoint"
    )(
        workflow,
        "LFQ",
        "non-DIA engineering pilot",
        validateColumnMapping = TRUE
    ))
    object <- workflow$state_manager$getState("protein_s4_initial")
    persisted <- protNonDiaPilotNamespaceValue(
        namespace,
        "persistProtNonDiaDesignArtifacts"
    )(workflow)
    if (!isTRUE(persisted$ok) || !isTRUE(persisted$committed)) {
        stop("non-DIA pilot design did not commit", call. = FALSE)
    }
    list(object = object, persisted = persisted)
}

protNonDiaPilotStateDigest <- function(object) {
    slots <- methods::slotNames(object)
    digests <- vapply(slots, \(slot_name) {
        value <- methods::slot(object, slot_name)
        digest::digest(
            if (is.data.frame(value)) as.data.frame(value) else value,
            algo = "sha256",
            serialize = TRUE,
            ascii = FALSE,
            serializeVersion = 3L
        )
    }, character(1))
    digest::digest(
        list(class = class(object), slots = slots, digests = digests),
        algo = "sha256",
        serialize = TRUE,
        serializeVersion = 3L
    )
}

protNonDiaPilotWorker <- function(args) {
    protNonDiaPilotRequire(args, c(
        "arm", "format", "source", "package_library", "run_dir"
    ))
    .libPaths(c(args$package_library, .libPaths()))
    namespace <- loadNamespace(
        "MultiScholaR",
        lib.loc = args$package_library
    )
    acknowledgement_path <- file.path(
        args$run_dir,
        "retention-sampled.fifo"
    )
    dir.create(args$run_dir, recursive = TRUE, showWarnings = FALSE)
    fifo_result <- processx::run(
        "mkfifo",
        acknowledgement_path,
        error_on_status = FALSE,
        echo = FALSE
    )
    if (fifo_result$status != 0L) {
        stop("non-DIA pilot FIFO creation failed", call. = FALSE)
    }
    paths <- protNonDiaPilotPaths(args$run_dir, args$format)
    workflow <- protNonDiaPilotWorkflow(namespace, paths)
    imported <- protNonDiaPilotImport(namespace, workflow, args)
    design <- protNonDiaPilotDesign(namespace, workflow)
    state_digest <- protNonDiaPilotStateDigest(design$object)
    capability_id <- paste0(
        "proteomics.", args$format, ".protein.lfq.v1"
    )
    descriptor <- protNonDiaPilotNamespaceValue(
        namespace,
        "protNonDiaReadthroughDescriptor"
    )(capability_id)
    settlement <- protNonDiaPilotNamespaceValue(
        namespace,
        "settleProtNonDiaArtifactWorkflowSafely"
    )(
        workflow,
        paths,
        "nondia-pilot",
        descriptor,
        rollout_fn = \(...) "evict"
    )
    if (!isTRUE(settlement$evicted)) {
        stop("non-DIA pilot settlement failed", call. = FALSE)
    }
    imported$imported <- NULL
    imported$staged <- NULL
    design$object <- NULL
    result <- list(
        status = "passed",
        arm = args$arm,
        format = args$format,
        state_digest = state_digest,
        source_fields_released = is.null(workflow$data_tbl) &&
            is.null(workflow$data_cln),
        cache_entries = if (inherits(
            workflow$state_manager,
            "ArtifactWorkflowState"
        )) {
            workflow$state_manager$getCacheInfo()$entries
        } else {
            NA_integer_
        },
        import_process_evidence = imported$persisted$process_evidence,
        design_process_evidence = design$persisted$process_evidence,
        state_manager_replaced = isTRUE(settlement$state_manager_replaced),
        complete_payload_returned = isTRUE(
            settlement$complete_payload_returned
        )
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
    arrow_pool_bytes <- if (requireNamespace("arrow", quietly = TRUE)) {
        as.numeric(arrow::default_memory_pool()$bytes_allocated)
    } else {
        0
    }
    resources <- publicationEmptyWorkerResources()
    resources$arrow_pool_bytes <- arrow_pool_bytes
    ledger <- list(
        schema = "multischolar.omics_publication_worker_resources",
        schema_version = "1.0.0",
        high_water = resources,
        retained = resources,
        terminal = publicationEmptyWorkerResources()
    )
    publicationWriteJson(
        ledger,
        file.path(args$run_dir, "worker-resources.json")
    )
    retention_state <- c(
        list(
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
        )
    )
    publicationWriteJson(
        retention_state,
        file.path(args$run_dir, "retention-state.json")
    )
    writeLines("ready", file.path(args$run_dir, "retention-ready"))
    acknowledgement <- fifo(
        acknowledgement_path,
        open = "rb",
        blocking = TRUE
    )
    on.exit(close(acknowledgement), add = TRUE)
    token <- readBin(acknowledgement, what = "raw", n = 1L)
    if (!identical(token, as.raw(1L))) {
        stop("non-DIA pilot retention acknowledgement differs", call. = FALSE)
    }
    invisible(0L)
}

protNonDiaPilotCategories <- function() {
    list(
        list(category = "diagnostics", pattern = "(^|/)diagnostics(/|$)"),
        list(category = "staging_snapshot", pattern = "(^|/)project(/|$)"),
        list(category = "duckdb_spill", pattern = "(^|/)duckdb-spill(/|$)"),
        list(category = "committed", pattern = "(^|/)project(/|$)"),
        list(category = "final", pattern = "(^|/)final(/|$)"),
        list(category = "harness", pattern = ".*")
    )
}

protNonDiaPilotMeasure <- function(args, arm, repetition) {
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
        "--vanilla", protNonDiaPilotScriptPath(),
        "--mode", "worker",
        "--arm", arm,
        "--format", args$format,
        "--source", normalizePath(args$source),
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
        unit_name = publicationSystemdUnitName("multischolar-nondia-pilot")
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

protNonDiaPilotRatio <- function(records, metric) {
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

protNonDiaPilotCampaign <- function(args) {
    protNonDiaPilotRequire(args, c(
        "format", "source", "old_library", "candidate_library",
        "output_root", "result"
    ))
    if (dir.exists(args$output_root) || file.exists(args$result)) {
        stop("non-DIA pilot output already exists", call. = FALSE)
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
            records[[length(records) + 1L]] <- protNonDiaPilotMeasure(
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
        schema = "multischolar.prot_nondia_engineering_pilot",
        schema_version = "1.0.0",
        owner_ticket_id = "OMICS-ART-073",
        status = if (valid && length(digests) == 1L) "passed" else "failed",
        format = args$format,
        source_bytes = unname(as.numeric(file.info(args$source)$size)),
        repetitions = args$repetitions,
        records = records,
        ratios = if (isTRUE(valid)) lapply(c(
            "peak_charged_memory_bytes",
            "retained_charged_memory_bytes",
            "elapsed_seconds",
            "peak_logical_disk_bytes",
            "final_logical_disk_bytes"
        ), \(metric) protNonDiaPilotRatio(records, metric)) else list(),
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

`%||%` <- function(left, right) {
    if (is.null(left)) right else left
}

if (identical(environment(), globalenv()) && !interactive()) {
    args <- protNonDiaPilotArgs(commandArgs(trailingOnly = TRUE))
    if (identical(args$mode, "worker")) {
        status <- tryCatch(
            protNonDiaPilotWorker(args),
            error = \(error) {
                message(conditionMessage(error))
                1L
            }
        )
    } else {
        status <- tryCatch(
            protNonDiaPilotCampaign(args),
            error = \(error) {
                message(conditionMessage(error))
                1L
            }
        )
    }
    quit(status = as.integer(status), save = "no")
}
