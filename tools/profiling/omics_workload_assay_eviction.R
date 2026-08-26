omicsAssayEvictionReadContract <- function(repo_root, contract_path, omic) {
    environment <- new.env(parent = .GlobalEnv)
    sys.source(
        file.path(repo_root, "tools", "profiling", "omics_workload_contract.R"),
        envir = environment
    )
    adapter_name <- paste0("omics_workload_", omic, ".R")
    adapter <- environment$omicsWorkloadLoadAdapter(
        file.path(repo_root, "tools", "profiling", adapter_name),
        environment$omicsWorkloadReadContract(contract_path)
    )
    contract <- environment$omicsWorkloadReadContract(contract_path)
    list(environment = environment, adapter = adapter, contract = contract)
}

omicsAssayEvictionPaths <- function(root, project_id, omic) {
    list(
        base_dir = root,
        project_id = project_id,
        omic_type = omic,
        omic_label = paste0(omic, "-study"),
        source_dir = file.path(root, "source"),
        results_dir = file.path(root, "results")
    )
}

omicsAssayEvictionWorkflow <- function(paths, omic, backend = "memory") {
    dir.create(paths$source_dir, recursive = TRUE, showWarnings = FALSE)
    dir.create(paths$results_dir, recursive = TRUE, showWarnings = FALSE)
    workflow <- new.env(parent = emptyenv())
    workflow$workflow_context <- createWorkflowContext(
        paths,
        omic,
        paste0(omic, "-study"),
        storage_policy = list(
            requested_backend = backend,
            requested_rollout = "dual_write",
            migration_requested = identical(backend, "artifact"),
            project_id = paths$project_id
        )
    )
    workflow$state_manager <- WorkflowState$new()
    workflow$artifact_stage_results <- list()
    workflow$processing_log <- list()
    workflow$tab_status <- list(
        setup_import = "pending",
        design_matrix = "disabled",
        quality_control = "disabled",
        normalization = "disabled"
    )
    workflow
}

omicsAssayEvictionInput <- function(payload_path, truth_path, omic) {
    data <- utils::read.delim(
        payload_path,
        check.names = FALSE,
        stringsAsFactors = FALSE,
        na.strings = "NA"
    )
    truth <- jsonlite::read_json(truth_path, simplifyVector = FALSE)
    assays <- unlist(truth$assays, use.names = FALSE)
    samples <- unlist(truth$sample_ids, use.names = FALSE)
    values <- lapply(assays, \(assay_name) {
        assay <- data[data$assay == assay_name, , drop = FALSE]
        assay$assay <- NULL
        row.names(assay) <- NULL
        assay
    })
    names(values) <- assays
    groups <- unlist(truth$group_assignments, use.names = FALSE)
    batches <- unlist(truth$batch_assignments, use.names = FALSE)
    replicates <- ave(seq_along(samples), groups, FUN = seq_along)
    design <- data.frame(
        Run = samples,
        group = groups,
        batch = batches,
        replicates = as.integer(replicates),
        tech_rep_group = paste(groups, replicates, sep = "_"),
        stringsAsFactors = FALSE
    )
    levels <- unique(groups)
    contrast <- paste0("group", levels[[2L]], "-group", levels[[1L]])
    friendly <- paste0(levels[[2L]], "_vs_", levels[[1L]])
    contrasts <- data.frame(
        contrasts = contrast,
        friendly_names = friendly,
        full_format = paste0(friendly, "=", contrast),
        stringsAsFactors = FALSE
    )
    if (identical(omic, "metabolomics")) {
        id_column <- "feature_id"
        annotation_column <- "annotation"
        workflow_type <- "metabolomics_standard"
        data_format <- "custom"
    } else {
        id_column <- "lipid_id"
        annotation_column <- "lipid_class"
        workflow_type <- "lipidomics_standard"
        data_format <- "lipidsearch"
    }
    cleaned <- lapply(values, \(assay) {
        assay[unique(c(id_column, annotation_column, samples))]
    })
    mapping <- if (identical(omic, "metabolomics")) {
        list(
            metabolite_id_col = id_column,
            annotation_col = annotation_column,
            sample_columns = samples,
            is_pattern = NA_character_
        )
    } else {
        list(
            lipid_id_col = id_column,
            annotation_col = annotation_column,
            sample_columns = samples,
            is_pattern = NA_character_
        )
    }
    source_files <- stats::setNames(
        as.list(rep(payload_path, length(assays))),
        assays
    )
    list(
        payload = list(
            assayList = values,
            sampleCols = samples,
            sampleNamesSanitized = FALSE,
            dataFormat = data_format,
            columnMapping = mapping,
            processingLog = list(evidence_class = "generated_scaling"),
            assaySourceFiles = source_files,
            sourceFiles = source_files
        ),
        cleaned = cleaned,
        design = design,
        contrasts = contrasts,
        config = list(
            globalParameters = list(workflow_type = workflow_type),
            deAnalysisParameters = list(formula_string = "~ 0 + group")
        ),
        id_column = id_column,
        annotation_column = annotation_column,
        source_bytes = as.numeric(utils::object.size(values)) +
            as.numeric(utils::object.size(cleaned))
    )
}

omicsAssayEvictionApplyImport <- function(workflow, input, omic) {
    payload <- input$payload
    if (identical(omic, "metabolomics")) {
        applyMetabImportWorkflowPayload(
            workflow,
            payload,
            logInfo = \(...) invisible(NULL)
        )
    } else {
        applyLipidImportResultToWorkflow(
            workflowData = workflow,
            assayList = payload$assayList,
            dataFormat = payload$dataFormat,
            lipidIdCol = payload$columnMapping$lipid_id_col,
            annotationCol = payload$columnMapping$annotation_col,
            sampleColumns = payload$columnMapping$sample_columns,
            isPattern = payload$columnMapping$is_pattern,
            logInfo = \(...) invisible(NULL)
        )
        workflow$processing_log$setup_import <- payload$processingLog
    }
    invisible(workflow)
}

omicsAssayEvictionCreateState <- function(workflow, input, omic) {
    workflow$data_cln <- input$cleaned
    workflow$design_matrix <- input$design
    workflow$contrasts_tbl <- input$contrasts
    workflow$config_list <- input$config
    mapping <- input$payload$columnMapping
    if (identical(omic, "metabolomics")) {
        object <- createMetaboliteAssayData(
            workflow$data_cln,
            workflow$design_matrix,
            metabolite_id_column = input$id_column,
            annotation_id_column = input$annotation_column,
            sample_id = "Run",
            group_id = "group",
            technical_replicate_id = "tech_rep_group",
            args = workflow$config_list
        )
        state_name <- "metab_raw_data_s4"
        workflow_type <- "metabolomics_standard"
    } else {
        object <- createLipidomicsAssayData(
            workflow$data_cln,
            workflow$design_matrix,
            lipid_id_column = input$id_column,
            annotation_id_column = input$annotation_column,
            sample_id = "Run",
            group_id = "group",
            technical_replicate_id = "tech_rep_group",
            args = workflow$config_list
        )
        state_name <- "lipid_raw_data_s4"
        workflow_type <- "lipidomics_standard"
    }
    workflow$column_mapping <- mapping
    workflow$state_manager$setWorkflowType(workflow_type)
    workflow$state_manager$saveState(
        state_name,
        object,
        workflow$config_list,
        paste(omic, "eviction workload")
    )
    invisible(workflow)
}

omicsAssayEvictionBuildWorkflow <- function(
    payload_path,
    truth_path,
    paths,
    omic,
    backend
) {
    input <- omicsAssayEvictionInput(payload_path, truth_path, omic)
    workflow <- omicsAssayEvictionWorkflow(paths, omic, backend)
    omicsAssayEvictionApplyImport(workflow, input, omic)
    if (identical(backend, "artifact")) {
        import <- if (identical(omic, "metabolomics")) {
            persistMetabImportArtifacts(
                workflow,
                input$payload,
                log_warn = \(...) invisible(NULL)
            )
        } else {
            persistLipidImportArtifacts(
                workflow,
                input$payload,
                log_warn = \(...) invisible(NULL)
            )
        }
        if (!isTRUE(import$ok)) {
            stop(paste(
                "artifact import preparation failed:",
                import$error_message %||% import$reason %||% "unknown error"
            ))
        }
    }
    omicsAssayEvictionCreateState(workflow, input, omic)
    if (identical(backend, "artifact")) {
        design <- if (identical(omic, "metabolomics")) {
            persistMetabDesignArtifacts(
                workflow,
                log_warn = \(...) invisible(NULL)
            )
        } else {
            persistLipidDesignArtifacts(
                workflow,
                log_warn = \(...) invisible(NULL)
            )
        }
        if (!isTRUE(design$ok)) {
            stop(paste(
                "artifact design preparation failed:",
                design$error_message %||% design$reason %||% "unknown error"
            ))
        }
    }
    list(workflow = workflow, source_bytes = input$source_bytes)
}

omicsAssayEvictionPrepare <- function(
    repo_root,
    contract_path,
    run_dir,
    backend,
    omic
) {
    devtools::load_all(repo_root, quiet = TRUE)
    loaded <- omicsAssayEvictionReadContract(repo_root, contract_path, omic)
    contract <- loaded$contract
    do.call(RNGkind, as.list(unlist(contract$rng$kind, use.names = FALSE)))
    set.seed(as.integer(contract$rng$seed))
    prepared <- loaded$adapter$prepare(list(contract = contract, run_dir = run_dir))
    project_root <- file.path(run_dir, "project")
    project_id <- paste0("eviction-", omic)
    paths <- omicsAssayEvictionPaths(project_root, project_id, omic)
    source_bytes <- omicsAssayEvictionInput(
        prepared$payload_path,
        prepared$truth_path,
        omic
    )$source_bytes
    if (identical(backend, "artifact")) {
        built <- omicsAssayEvictionBuildWorkflow(
            prepared$payload_path,
            prepared$truth_path,
            paths,
            omic,
            "artifact"
        )
        source_bytes <- built$source_bytes
    }
    result <- list(
        payload_path = normalizePath(prepared$payload_path, mustWork = TRUE),
        truth_path = normalizePath(prepared$truth_path, mustWork = TRUE),
        paths = paths,
        source_bytes = source_bytes
    )
    jsonlite::write_json(
        result,
        file.path(run_dir, "prepared.json"),
        auto_unbox = TRUE,
        pretty = TRUE,
        digits = 17
    )
    invisible(result)
}

omicsAssayEvictionWarmRuntime <- function() {
    root <- tempfile("eviction-warm-")
    dir.create(root)
    on.exit(unlink(root, recursive = TRUE, force = TRUE), add = TRUE)
    encoded <- encodeArtifactTable(
        data.frame(id = 1L, value = 1, stringsAsFactors = FALSE),
        owner = "eviction warmup"
    )
    path <- file.path(root, "warm.parquet")
    do.call(
        arrow::write_parquet,
        c(
            list(x = encoded$payload, sink = path),
            artifactParquetWriteArgs(encoded)
        )
    )
    invisible(arrow::read_parquet(path, as_data_frame = TRUE))
}

omicsAssayEvictionWarmMemory <- function(omic) {
    samples <- sprintf("WARM_S%02d", 1:6)
    groups <- rep(c("CTRL", "TREAT"), each = 3L)
    table <- data.frame(
        feature_id = c("warm_1", "warm_2"),
        annotation = c("warm", "warm"),
        matrix(seq_len(12L), nrow = 2L),
        stringsAsFactors = FALSE,
        check.names = FALSE
    )
    names(table)[-(1:2)] <- samples
    if (identical(omic, "lipidomics")) {
        names(table)[1:2] <- c("lipid_id", "lipid_class")
    }
    id_column <- names(table)[[1L]]
    annotation_column <- names(table)[[2L]]
    design <- data.frame(
        Run = samples,
        group = groups,
        batch = rep(c("B1", "B2"), length.out = 6L),
        replicates = rep(1:3, 2L),
        tech_rep_group = paste(groups, rep(1:3, 2L), sep = "_"),
        stringsAsFactors = FALSE
    )
    config <- list(
        globalParameters = list(workflow_type = if (identical(
            omic,
            "metabolomics"
        )) {
            "metabolomics_standard"
        } else {
            "lipidomics_standard"
        }),
        deAnalysisParameters = list(formula_string = "~ 0 + group")
    )
    if (identical(omic, "metabolomics")) {
        object <- createMetaboliteAssayData(
            list(LCMS_Pos = table),
            design,
            metabolite_id_column = id_column,
            annotation_id_column = annotation_column,
            sample_id = "Run",
            group_id = "group",
            technical_replicate_id = "tech_rep_group",
            args = config
        )
    } else {
        object <- createLipidomicsAssayData(
            list(LCMS_Pos = table),
            design,
            lipid_id_column = id_column,
            annotation_id_column = annotation_column,
            sample_id = "Run",
            group_id = "group",
            technical_replicate_id = "tech_rep_group",
            args = config
        )
    }
    manager <- WorkflowState$new()
    manager$saveState("warm_state", object, config, "eviction warmup")
    invisible(TRUE)
}

omicsAssayEvictionWarmBackend <- function(prepared, omic, backend) {
    started <- proc.time()[["elapsed"]]
    if (identical(backend, "artifact")) {
        settled <- omicsAssayEvictionResumeSettled(prepared, omic)
        settled$workflow$state_manager$close()
    } else {
        omicsAssayEvictionWarmMemory(omic)
    }
    invisible(gc(full = TRUE))
    artifactReleaseProcessAllocator()
    proc.time()[["elapsed"]] - started
}

omicsAssayEvictionBlankWorkflow <- function(paths, omic) {
    workflow <- new.env(parent = emptyenv())
    workflow$workflow_context <- createWorkflowContext(
        paths,
        omic,
        paste0(omic, "-study"),
        storage_policy = list()
    )
    workflow$state_manager <- WorkflowState$new()
    workflow$artifact_stage_results <- list()
    workflow$tab_status <- list(
        setup_import = "pending",
        design_matrix = "disabled",
        quality_control = "disabled",
        normalization = "disabled"
    )
    workflow
}

omicsAssayEvictionResumeSettled <- function(prepared, omic) {
    rss <- list(start = omicsAssayEvictionRss())
    workflow <- omicsAssayEvictionBlankWorkflow(prepared$paths, omic)
    rss$blank <- omicsAssayEvictionRss()
    if (identical(omic, "metabolomics")) {
        context <- createMetabResumeContext(
            prepared$paths,
            "metabolomics-study"
        )$context
        rss$context <- omicsAssayEvictionRss()
        bundle <- hydrateMetabResumeBundle(
            context,
            retain_source_payloads = FALSE
        )
        rss$bundle <- omicsAssayEvictionRss()
        applyMetabResumeBundle(workflow, bundle)
        rss$applied <- omicsAssayEvictionRss()
        result <- settleMetabArtifactWorkflowSafely(
            workflow,
            rollout_fn = \(...) "evict",
            log_warn = \(...) invisible(NULL)
        )
    } else {
        context <- createLipidResumeContext(
            prepared$paths,
            "lipidomics-study"
        )$context
        rss$context <- omicsAssayEvictionRss()
        bundle <- hydrateLipidResumeBundle(
            context,
            retain_source_payloads = FALSE
        )
        rss$bundle <- omicsAssayEvictionRss()
        applyLipidResumeBundle(workflow, bundle)
        rss$applied <- omicsAssayEvictionRss()
        result <- settleLipidArtifactWorkflowSafely(
            workflow,
            rollout_fn = \(...) "evict",
            log_warn = \(...) invisible(NULL)
        )
    }
    if (!isTRUE(result$evicted) ||
        !isTRUE(result$source_hydration_avoided) ||
        !identical(workflow$state_manager$getCacheInfo()$entries, 0L)) {
        stop("artifact workflow did not settle payload-free")
    }
    rss$settled <- omicsAssayEvictionRss()
    list(workflow = workflow, settlement = result, rss_checkpoints = rss)
}

omicsAssayEvictionRss <- function() {
    as.numeric(ps::ps_memory_info(ps::ps_handle())[["rss"]])
}

omicsAssayEvictionArtifactBytes <- function(paths) {
    files <- list.files(paths$base_dir, recursive = TRUE, full.names = TRUE)
    files <- files[file.exists(files) & !dir.exists(files)]
    if (!length(files)) return(0)
    sum(as.numeric(file.info(files)$size), na.rm = TRUE)
}

omicsAssayEvictionMain <- function(args, omic, schema, script_path, repo_root) {
    if (identical(args[[1L]], "--prepare")) {
        if (length(args) != 4L) stop("invalid eviction preparation arguments")
        omicsAssayEvictionPrepare(
            repo_root,
            normalizePath(args[[2L]], mustWork = TRUE),
            normalizePath(args[[3L]], mustWork = TRUE),
            match.arg(args[[4L]], c("memory", "artifact")),
            omic
        )
        return(invisible(TRUE))
    }
    if (length(args) != 3L) {
        stop("usage: <contract.json> <memory|artifact> <output.json>")
    }
    contract_path <- normalizePath(args[[1L]], mustWork = TRUE)
    backend <- match.arg(args[[2L]], c("memory", "artifact"))
    output_path <- args[[3L]]
    run_dir <- tempfile(paste0(omic, "-eviction-worker-"))
    dir.create(run_dir, recursive = TRUE)
    on.exit(unlink(run_dir, recursive = TRUE, force = TRUE), add = TRUE)
    preparation_started <- proc.time()[["elapsed"]]
    preparation <- processx::run(
        file.path(R.home("bin"), "Rscript"),
        c(
            "--vanilla", script_path, "--prepare", contract_path, run_dir,
            backend
        ),
        wd = repo_root,
        error_on_status = TRUE,
        echo = FALSE
    )
    preparation_elapsed <- proc.time()[["elapsed"]] - preparation_started
    if (!identical(preparation$status, 0L)) {
        stop("eviction preparation process failed")
    }
    prepared <- jsonlite::read_json(
        file.path(run_dir, "prepared.json"),
        simplifyVector = FALSE
    )
    loaded <- omicsAssayEvictionReadContract(repo_root, contract_path, omic)
    contract <- loaded$contract
    payload_digest <- digest::digest(
        file = prepared$payload_path,
        algo = "sha256",
        serialize = FALSE
    )
    truth_digest <- digest::digest(
        file = prepared$truth_path,
        algo = "sha256",
        serialize = FALSE
    )
    if (!identical(payload_digest, contract$expected_digests$payload_sha256) ||
        !identical(truth_digest, contract$expected_digests$truth_sha256)) {
        stop("frozen eviction workload digest mismatch")
    }
    devtools::load_all(repo_root, quiet = TRUE)
    omicsAssayEvictionWarmRuntime()
    lifecycle_warm_elapsed <- omicsAssayEvictionWarmBackend(
        prepared,
        omic,
        backend
    )
    invisible(gc(full = TRUE, reset = TRUE))
    Sys.sleep(0.25)
    rss_baseline <- omicsAssayEvictionRss()
    workload_started <- proc.time()[["elapsed"]]
    if (identical(backend, "memory")) {
        memory_root <- file.path(run_dir, "memory-project")
        paths <- omicsAssayEvictionPaths(
            memory_root,
            paste0("eviction-", omic, "-memory"),
            omic
        )
        built <- omicsAssayEvictionBuildWorkflow(
            prepared$payload_path,
            prepared$truth_path,
            paths,
            omic,
            "memory"
        )
        workflow <- built$workflow
        source_bytes <- built$source_bytes
        state <- workflow$state_manager$getState()
        retained_graph_bytes <- source_bytes +
            as.numeric(utils::object.size(state))
        cache_entries <- 0L
        artifact_bytes <- 0
        settlement_reason <- "memory_source_graph_retained"
    } else {
        settled <- omicsAssayEvictionResumeSettled(prepared, omic)
        workflow <- settled$workflow
        source_bytes <- as.numeric(prepared$source_bytes)
        retained_graph_bytes <- sum(vapply(
            c("data_tbl", "data_cln"),
            \(name) as.numeric(utils::object.size(workflow[[name]])),
            numeric(1)
        ))
        cache_entries <- workflow$state_manager$getCacheInfo()$entries
        artifact_bytes <- omicsAssayEvictionArtifactBytes(prepared$paths)
        settlement_reason <- settled$settlement$reason
        rss_checkpoints <- settled$rss_checkpoints
    }
    workload_elapsed <- proc.time()[["elapsed"]] - workload_started
    invisible(gc(full = TRUE))
    allocator_pages_released <- if (identical(backend, "artifact")) {
        isTRUE(artifactReleaseProcessAllocator())
    } else {
        FALSE
    }
    Sys.sleep(0.25)
    rss_settled <- omicsAssayEvictionRss()
    workload_rss <- max(0, rss_settled - rss_baseline)
    result <- list(
        schema = schema,
        schema_version = "2.0.0",
        workload_id = contract$workload_id,
        contract_digest = loaded$environment$omicsWorkloadDigest(contract),
        privacy_classification = contract$privacy$classification,
        backend = backend,
        payload_sha256 = payload_digest,
        truth_sha256 = truth_digest,
        feature_count = as.integer(contract$dimensions$feature_count),
        sample_count = as.integer(contract$dimensions$sample_count),
        assay_count = as.integer(contract$dimensions$assay_count),
        source_bytes = source_bytes,
        artifact_bytes = artifact_bytes,
        preparation_elapsed_seconds = preparation_elapsed,
        lifecycle_warm_elapsed_seconds = lifecycle_warm_elapsed,
        workload_elapsed_seconds = workload_elapsed,
        rss_warm_baseline = rss_baseline,
        rss_settled = rss_settled,
        workload_retained_rss_bytes = workload_rss,
        rss_before_eviction = rss_baseline,
        rss_after_eviction = rss_settled,
        retained_source_bytes = if (identical(backend, "artifact")) {
            0
        } else {
            source_bytes
        },
        retained_workflow_graph_bytes = retained_graph_bytes,
        post_settlement_cache_entries = cache_entries,
        allocator_pages_released = allocator_pages_released,
        settlement_reason = settlement_reason,
        rss_checkpoints = if (identical(backend, "artifact")) {
            rss_checkpoints
        } else {
            list(start = rss_baseline, settled = rss_settled)
        },
        rss_measurement = paste(
            "fresh process workload-attributable retained RSS above an identical",
            "package and Arrow warm baseline; absolute RSS retained separately"
        ),
        generated_claim_boundary = paste(
            "resource and eviction evidence only;",
            "not parser or biological validation"
        )
    )
    jsonlite::write_json(
        result,
        output_path,
        auto_unbox = TRUE,
        pretty = TRUE,
        digits = 17
    )
    close <- tryCatch(workflow$state_manager$close, error = \(...) NULL)
    if (is.function(close)) try(close(), silent = TRUE)
    invisible(result)
}
