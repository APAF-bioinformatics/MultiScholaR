omics014SkipDependencies <- function() {
    for (package in c("DBI", "duckdb", "filelock", "arrow")) {
        testthat::skip_if_not_installed(package)
    }
}

omics014Rss <- function() {
    if (!requireNamespace("ps", quietly = TRUE)) return(NA_real_)
    as.numeric(ps::ps_memory_info(ps::ps_handle())[["rss"]])
}

omics014Capability <- function() {
    candidates <- Filter(\(capability) {
        identical(capability$identity$input_format, "diann") &&
            identical(capability$identity$data_level, "peptide")
    }, workflowCapabilityCatalogue())
    stopifnot(length(candidates) == 1L)
    capability <- candidates[[1L]]
    capability$artifact_eligible <- TRUE
    capability$auto_eligible <- TRUE
    capability$maximum_artifact_rollout <- "dual_write"
    capability
}

omics014Context <- function(project_root, project_id = "omics-014-project") {
    context <- createWorkflowContext(
        experiment_paths = list(
            base_dir = project_root,
            omic_label = "proteomics_study"
        ),
        omic_type = "proteomics",
        storage_policy = list(
            requested_backend = "artifact",
            requested_rollout = "dual_write",
            migration_requested = TRUE,
            project_id = project_id
        )
    )
    bindWorkflowContextFromImport(
        context,
        workflow_type = "DIA",
        input_format = "diann",
        data_level = "peptide",
        capabilities = list(omics014Capability())
    )
    context
}

omics014Manager <- function(context) {
    manager <- newProtDiaArtifactStateManager(context)
    manager$setWorkflowType("DIA")
    manager
}

omics014Protein <- function() {
    values <- module_ci_prot_norm_values()
    table <- data.frame(
        Protein.Group = rownames(values),
        as.data.frame(values, check.names = FALSE),
        check.names = FALSE,
        stringsAsFactors = FALSE
    )
    rownames(table) <- NULL
    design <- module_ci_prot_norm_design(colnames(values))
    args <- module_ci_prot_norm_args(rownames(values))
    args$peptide_qc_audit <- list(
        records = list(list(record_id = "peptide-audit-final")),
        current_record_id = "peptide-audit-final"
    )
    args$protein_qc_audit <- list(
        records = list(list(record_id = "protein-audit-final")),
        current_record_id = "protein-audit-final"
    )
    ProteinQuantitativeData(
        protein_quant_table = table,
        protein_id_column = "Protein.Group",
        protein_id_table = data.frame(
            Protein.Group = rownames(values),
            Protein.Ids = paste0("ACCESSION_", seq_len(nrow(values))),
            stringsAsFactors = FALSE
        ),
        design_matrix = design,
        sample_id = "Run",
        group_id = "group",
        technical_replicate_id = "replicates",
        args = args
    )
}

omics014Workflow <- function(manager, context = NULL) {
    workflow <- new.env(parent = emptyenv())
    workflow$state_manager <- manager
    workflow$config_list <- list(globalParameters = list(workflow_type = "DIA"))
    workflow$tab_status <- list(
        normalization = "pending",
        differential_expression = "disabled"
    )
    workflow$state_update_trigger <- NULL
    workflow$protein_counts <- list()
    workflow$artifact_stage_results <- list()
    workflow$normalization_state_refs <- list()
    workflow$ruv_normalised_for_da_analysis_obj <- NULL
    workflow$final_for_da_ref <- NULL
    if (!is.null(context)) workflow$workflow_context <- context
    workflow
}

omics014NormData <- function(manager) {
    norm_data <- new.env(parent = emptyenv())
    norm_data$state_manager <- manager
    norm_data$state_refs <- list()
    norm_data$normalized_protein_obj <- NULL
    norm_data$ruv_normalized_obj <- NULL
    norm_data$correlation_filtered_obj <- NULL
    norm_data$normalization_complete <- FALSE
    norm_data$ruv_complete <- FALSE
    norm_data$correlation_filtering_complete <- FALSE
    norm_data
}

omics014SaveInitial <- function(manager, object) {
    manager$saveState(
        "protein_replicate_filtered",
        object,
        list(stage = "protein_replicate_filtered"),
        "Saved protein-QC parent",
        audit_metadata = list(record_id = "initial:protein_replicate_filtered")
    )
    invisible(object)
}

omics014Manifest <- function(context, manager, state_name) {
    row <- manager$states[[state_name]]
    store <- newArtifactStore(context$getPaths(), context$getIdentity()$project_id)
    artifactWorkflowStateReadManifest(store, row$manifest_relative_path)
}

omics014Metadata <- function(context, manager, state_name) {
    manifest <- omics014Manifest(context, manager, state_name)
    artifactWorkflowStateUnserializeMetadata(
        manifest$data$metadata_json,
        "OMICS-ART-014 test metadata"
    )
}

omics014ExpectScientificEqual <- function(expected, actual) {
    for (slot_name in setdiff(methods::slotNames(expected), "args")) {
        expect_identical(
            methods::slot(actual, slot_name),
            methods::slot(expected, slot_name),
            info = slot_name
        )
    }
    expected_args <- expected@args
    actual_args <- actual@args
    expected_args$normalization_artifact_audit <- NULL
    actual_args$normalization_artifact_audit <- NULL
    expected_args$normalization_artifact_parameters <- NULL
    actual_args$normalization_artifact_parameters <- NULL
    expect_identical(actual_args, expected_args)
    expect_identical(methods::validObject(actual, test = TRUE), TRUE)
}

omics014FailAt <- function(target_stage) {
    force(target_stage)
    \(stage, context) {
        if (identical(stage, target_stage)) {
            rlang::abort(
                paste("injected OMICS-ART-014 failure at", stage),
                class = "multischolar_test_artifact_failure"
            )
        }
        invisible(context)
    }
}

omics014SaveNormalized <- function(
    workflow,
    norm_data,
    before,
    normalized,
    method,
    ruv_mode = "manual"
) {
    parameters <- list(
        normalization_method = method,
        ruv_mode = ruv_mode,
        ruv_grouping_variable = "group",
        skip_reason = if (identical(ruv_mode, "skip")) "dataset constraints" else NULL,
        normalized_matrix_file = "normalized_protein_matrix_pre_ruv.tsv"
    )
    saved <- saveProtNormState(
        workflow,
        workflow$state_manager,
        before,
        normalized,
        "normalization",
        "normalized",
        parameters,
        "Applied between-samples normalization",
        parameters,
        transformation_type = "normalization"
    )
    settleProtNormArtifactState(
        workflow,
        norm_data,
        "normalization",
        "normalized",
        saved
    )
}

omics014RunChain <- function(
    manager,
    workflow,
    norm_data,
    threshold,
    ruv_mode = "manual"
) {
    before <- manager$getState()
    methods <- module_ci_prot_norm_s4_methods()
    normalized <- methods$normalise(before, normalisation_method = "none")
    normalized <- omics014SaveNormalized(
        workflow,
        norm_data,
        before,
        normalized,
        "none",
        ruv_mode
    )
    workflow$normalization_context <- list(
        normalization_method = "none",
        ruv_mode = ruv_mode,
        ruv_grouping_variable = "group",
        ruv_parameters = list()
    )
    controls <- module_ci_prot_norm_control_index(
        normalized@protein_quant_table$Protein.Group
    )
    if (!identical(ruv_mode, "skip")) {
        corrected <- methods$ruviii(
            normalized,
            ruv_grouping_variable = "group",
            ruv_number_k = 2L,
            ctrl = controls
        )
        optimization_results <- data.frame(
            percentage = seq_len(300L),
            separation_score = seq_len(300L) / 1000,
            best_k = rep(2L, 300L),
            composite_score = rev(seq_len(300L)) / 1000,
            status = rep("ok", 300L),
            stringsAsFactors = FALSE
        )
        parameters <- list(
            normalization_method = "none",
            ruv_mode = ruv_mode,
            ruv_grouping_variable = "group",
            ruv_k = 2L,
            percentage_as_neg_ctrl = 40,
            control_genes_index = list(
                owner = "ProteinQuantitativeData@args$ruvIII_C_Varying$ctrl",
                selected = as.integer(sum(controls)),
                digest = .peptideQcDigest(controls)
            ),
            optimization_results = optimization_results,
            separation_metric = "max_difference",
            k_penalty_weight = 0.5,
            adaptive_k_penalty = TRUE,
            ruv_result_file = "ruv_optimization_results.RDS"
        )
        corrected <- saveProtNormState(
            workflow,
            manager,
            normalized,
            corrected,
            "ruv_correction",
            "ruv_corrected",
            protDiaNormAuditParameters(parameters),
            "Applied RUV-III correction",
            parameters,
            transformation_type = "batch_correction"
        )
        settleProtNormArtifactState(
            workflow,
            norm_data,
            "ruv_correction",
            "ruv_corrected",
            corrected
        )
        workflow$ruv_optimization_result <- list(
            best_k = 2L,
            best_percentage = 40,
            ruv_skipped = FALSE,
            optimization_results = optimization_results,
            separation_metric_used = "max_difference"
        )
        correlation_input <- corrected
    } else {
        workflow$ruv_optimization_result <- buildProtNormSkippedRuvResult()
        correlation_input <- normalized
    }
    norm_data$ruv_normalized_obj <- correlation_input
    pairs <- module_ci_prot_norm_correlation_pairs("pass_all")
    if (identical(threshold, 1)) {
        pairs$pearson_correlation <- 1
    }
    filtered <- filterSamplesByProteinCorrelationThreshold(
        correlation_input,
        pearson_correlation_per_pair = pairs,
        min_pearson_correlation_threshold = threshold
    )
    norm_data$correlation_filtered_obj <- filtered
    finalizeProtNormCorrelationWorkflowState(
        filtered,
        workflow,
        normData = norm_data,
        correlationThreshold = threshold,
        skipped = FALSE,
        timeFn = \() as.POSIXct("2026-08-21 12:00:00", tz = "UTC"),
        messageFn = \(...) invisible(NULL),
        catFn = \(...) invisible(NULL)
    )
    final <- manager$getState("correlation_filtered")
    omics014ExpectScientificEqual(filtered, final)
    final
}

test_that("all protein normalization engines retain memory and artifact parity", {
    omics014SkipDependencies()
    methods <- module_ci_prot_norm_s4_methods()
    for (method in c("none", "cyclicloess", "quantile", "scale")) {
        before <- omics014Protein()
        expected <- methods$normalise(before, normalisation_method = method)

        memory_manager <- newWorkflowState()
        memory_manager$setWorkflowType("DIA")
        omics014SaveInitial(memory_manager, before)
        memory_workflow <- omics014Workflow(memory_manager)
        memory_norm <- omics014NormData(memory_manager)
        memory <- omics014SaveNormalized(
            memory_workflow,
            memory_norm,
            before,
            expected,
            method
        )

        context <- omics014Context(
            withr::local_tempdir(),
            paste0("omics-014-normalization-", method)
        )
        artifact_manager <- omics014Manager(context)
        omics014SaveInitial(artifact_manager, before)
        artifact_workflow <- omics014Workflow(artifact_manager, context)
        artifact_norm <- omics014NormData(artifact_manager)
        artifact <- omics014SaveNormalized(
            artifact_workflow,
            artifact_norm,
            before,
            expected,
            method
        )

        omics014ExpectScientificEqual(memory, artifact)
        expect_identical(resolveProteinQuantIdentityColumn(artifact), "Protein.Group")
        record <- protDiaNormCurrentRecord(artifact)
        expect_identical(record$stage_id, "normalization")
        expect_identical(record$resolved_parameters$normalization_method, method)
        expect_identical(record$protein_qc_record_id, "protein-audit-final")
        expect_identical(artifact_manager$getCacheInfo()$entries, 1L)
        artifact_manager$close()
    }
})

test_that("RUV and correlation generations externalize args and settle S4 owners", {
    omics014SkipDependencies()
    for (threshold in c(0, 0.75, 1)) {
        rss_before <- omics014Rss()
        before <- omics014Protein()
        memory_manager <- newWorkflowState()
        memory_manager$setWorkflowType("DIA")
        omics014SaveInitial(memory_manager, before)
        memory_workflow <- omics014Workflow(memory_manager)
        memory_final <- omics014RunChain(
            memory_manager,
            memory_workflow,
            omics014NormData(memory_manager),
            threshold,
            ruv_mode = if (threshold == 0.75) "automatic" else "manual"
        )
        rss_after_memory <- omics014Rss()

        context <- omics014Context(
            withr::local_tempdir(),
            paste0("omics-014-chain-", gsub("\\.", "-", threshold))
        )
        manager <- omics014Manager(context)
        omics014SaveInitial(manager, before)
        workflow <- omics014Workflow(manager, context)
        norm_data <- omics014NormData(manager)
        final <- omics014RunChain(
            manager,
            workflow,
            norm_data,
            threshold,
            ruv_mode = if (threshold == 0.75) "automatic" else "manual"
        )
        rss_after_artifact <- omics014Rss()

        omics014ExpectScientificEqual(memory_final, final)
        rss <- c(rss_before, rss_after_memory, rss_after_artifact)
        expect_true(all(is.na(rss) | (is.finite(rss) & rss > 0)))
        expect_identical(memory_manager$getHistory(), manager$getHistory())
        expect_identical(
            manager$getHistory(),
            c(
                "initial",
                "protein_replicate_filtered",
                "normalized",
                "ruv_corrected",
                "correlation_filtered"
            )
        )
        expect_identical(manager$getCacheInfo()$entries, 1L)
        for (stage_id in .PROT_DIA_NORM_STAGES) {
            metrics <- workflow$artifact_stage_results[[stage_id]]$metrics
            expect_gt(metrics$hydrated_object_bytes, 0)
            expect_gt(metrics$generation_artifact_bytes, 0)
            expect_gte(metrics$generation_artifact_count, 1L)
        }
        expect_null(norm_data$normalized_protein_obj)
        expect_null(norm_data$ruv_normalized_obj)
        expect_null(norm_data$correlation_filtered_obj)
        expect_null(workflow$ruv_normalised_for_da_analysis_obj)
        expect_identical(workflow$final_for_da_ref$state_name, "correlation_filtered")
        expect_s4_class(final, "ProteinQuantitativeData")
        expect_identical(final@protein_id_column, "Protein.Group")
        expect_identical(
            final@args$filterSamplesByProteinCorrelationThreshold$
                min_pearson_correlation_threshold,
            threshold
        )
        expect_identical(
            methods::validObject(final, test = TRUE),
            TRUE
        )

        ruv_metadata <- omics014Metadata(context, manager, "ruv_corrected")
        ruv_owners <- vapply(ruv_metadata$payloads, `[[`, character(1), "owner")
        expect_true(any(grepl("optimization_results", ruv_owners, fixed = TRUE)))
        correlation_metadata <- omics014Metadata(
            context,
            manager,
            "correlation_filtered"
        )
        correlation_owners <- vapply(
            correlation_metadata$payloads,
            `[[`,
            character(1),
            "owner"
        )
        expect_true(any(grepl(
            "pearson_correlation_per_pair",
            correlation_owners,
            fixed = TRUE
        )))
        expect_identical(
            final@args$normalization_artifact_parameters$ruv_correction$
                optimization_results,
            workflow$ruv_optimization_result$optimization_results
        )
        stages <- vapply(
            final@args$normalization_artifact_audit$records,
            `[[`,
            character(1),
            "stage_id"
        )
        expect_identical(stages, c(
            "normalization",
            "ruv_correction",
            "correlation_filter"
        ))
        manager$close()
    }
})

test_that("RUV and correlation skips preserve exact values and audit decisions", {
    omics014SkipDependencies()
    context <- omics014Context(withr::local_tempdir(), "omics-014-skip")
    manager <- omics014Manager(context)
    before <- omics014Protein()
    omics014SaveInitial(manager, before)
    workflow <- omics014Workflow(manager, context)
    norm_data <- omics014NormData(manager)
    methods <- module_ci_prot_norm_s4_methods()
    normalized <- methods$normalise(before, normalisation_method = "none")
    normalized <- omics014SaveNormalized(
        workflow,
        norm_data,
        before,
        normalized,
        "none",
        "skip"
    )
    workflow$normalization_context <- list(
        ruv_mode = "skip",
        ruv_grouping_variable = "group"
    )
    workflow$ruv_optimization_result <- buildProtNormSkippedRuvResult()
    norm_data$correlation_filtered_obj <- normalized
    finalizeProtNormCorrelationWorkflowState(
        normalized,
        workflow,
        normData = norm_data,
        correlationThreshold = 0,
        skipped = TRUE,
        timeFn = \() as.POSIXct("2026-08-21 12:00:00", tz = "UTC"),
        messageFn = \(...) invisible(NULL),
        catFn = \(...) invisible(NULL)
    )
    final <- manager$getState()
    omics014ExpectScientificEqual(normalized, final)
    expect_identical(
        manager$getHistory(),
        c(
            "initial",
            "protein_replicate_filtered",
            "normalized",
            "correlation_filtered"
        )
    )
    records <- final@args$normalization_artifact_audit$records
    expect_identical(vapply(records, `[[`, character(1), "status"), c(
        "applied",
        "skipped"
    ))
    expect_identical(
        records[[1L]]$resolved_parameters$skip_reason,
        "dataset constraints"
    )
    expect_identical(
        records[[2L]]$resolved_parameters$output_files,
        c(
            "normalised_results_cln_with_replicates.tsv",
            "normalised_results_cln_with_replicates.RDS"
        )
    )
    manager$close()
})

test_that("artifact failures, reverts, globals, and stage switches are isolated", {
    omics014SkipDependencies()
    context <- omics014Context(withr::local_tempdir(), "omics-014-rollback")
    manager <- omics014Manager(context)
    before <- omics014Protein()
    omics014SaveInitial(manager, before)
    workflow <- omics014Workflow(manager, context)
    norm_data <- omics014NormData(manager)
    methods <- module_ci_prot_norm_s4_methods()
    normalized <- methods$normalise(before, normalisation_method = "none")
    normalized <- omics014SaveNormalized(
        workflow,
        norm_data,
        before,
        normalized,
        "none"
    )
    parent_generation <- manager$getCurrentGenerationId()
    controls <- module_ci_prot_norm_control_index(
        normalized@protein_quant_table$Protein.Group
    )
    corrected <- methods$ruviii(
        normalized,
        ruv_grouping_variable = "group",
        ruv_number_k = 2L,
        ctrl = controls
    )
    parameters <- list(
        normalization_method = "none",
        ruv_mode = "manual",
        ruv_grouping_variable = "group",
        ruv_k = 2L,
        control_genes_index = controls
    )
    expect_error(
        saveProtNormState(
            workflow,
            manager,
            normalized,
            corrected,
            "ruv_correction",
            "ruv_corrected",
            protDiaNormAuditParameters(parameters),
            "Failed RUV correction",
            parameters,
            failure_injector = omics014FailAt("before_state_registry_commit")
        ),
        class = "multischolar_test_artifact_failure"
    )
    expect_identical(manager$getCurrentGenerationId(), parent_generation)
    expect_identical(manager$getCurrentStateName(), "normalized")
    expect_identical(manager$getState(), normalized)
    corrected <- saveProtNormState(
        workflow,
        manager,
        normalized,
        corrected,
        "ruv_correction",
        "ruv_corrected",
        protDiaNormAuditParameters(parameters),
        "Applied RUV correction",
        parameters
    )
    expect_identical(manager$getState(), corrected)
    expect_identical(revertProtDiaNormState(workflow, "normalized"), normalized)
    expect_identical(manager$getCurrentStateName(), "normalized")

    old_global <- get0("config_list", envir = .GlobalEnv, inherits = FALSE)
    had_global <- exists("config_list", envir = .GlobalEnv, inherits = FALSE)
    withr::defer({
        if (had_global) {
            assign("config_list", old_global, envir = .GlobalEnv)
        } else if (exists("config_list", envir = .GlobalEnv, inherits = FALSE)) {
            rm("config_list", envir = .GlobalEnv)
        }
    })
    hostile <- list(normaliseBetweenSamples = list(method = "hostile"))
    assign("config_list", hostile, envir = .GlobalEnv)
    workflow$config_list <- list(normaliseBetweenSamples = list(method = "none"))
    observed_method <- NULL
    runProtNormBetweenSamplesStep(
        currentS4 = before,
        normMethod = "quantile",
        normData = norm_data,
        workflowData = workflow,
        proteinQcDir = NULL,
        normaliseBetweenSamplesFn = \(object, normalisation_method) {
            observed_method <<- normalisation_method
            object
        },
        captureCheckpointFn = \(...) invisible(NULL),
        messageFn = \(...) invisible(NULL)
    )
    expect_identical(observed_method, "quantile")
    expect_identical(workflow$config_list$normaliseBetweenSamples$method, "quantile")
    expect_identical(get("config_list", envir = .GlobalEnv), hostile)
    manager$close()

    for (stage_id in .PROT_DIA_NORM_STAGES) {
        stage_context <- omics014Context(
            withr::local_tempdir(),
            paste0("omics-014-disabled-", stage_id)
        )
        artifact_manager <- omics014Manager(stage_context)
        memory_manager <- newWorkflowState()
        memory_manager$setWorkflowType("DIA")
        parent <- omics014Protein()
        omics014SaveInitial(artifact_manager, parent)
        omics014SaveInitial(memory_manager, parent)
        stage_workflow <- omics014Workflow(memory_manager, stage_context)
        option <- stats::setNames(
            list("disabled"),
            paste0("multischolar.prot_dia.normalization_artifacts.", stage_id)
        )
        withr::with_options(option, {
            child <- saveProtNormState(
                stage_workflow,
                memory_manager,
                parent,
                parent,
                stage_id,
                paste0(stage_id, "_memory_only"),
                list(stage = stage_id),
                "Memory-only stage"
            )
            expect_identical(child@protein_quant_table, parent@protein_quant_table)
            expect_identical(
                artifact_manager$getCurrentStateName(),
                "protein_replicate_filtered"
            )
            reads <- 0L
            resolver_manager <- new.env(parent = emptyenv())
            resolver_manager$getState <- \(...) {
                reads <<- reads + 1L
                parent
            }
            resolver_workflow <- new.env(parent = emptyenv())
            resolver_workflow$state_manager <- resolver_manager
            expect_null(resolveProtNormStateObject(
                workflow_data = resolver_workflow,
                state_names = "protein_replicate_filtered",
                stage_id = stage_id
            ))
            expect_identical(reads, 0L)
        })
        artifact_manager$close()
    }
})

test_that("normalization QC products preserve current filenames and plot getters", {
    object <- omics014Protein()
    paths <- list()
    saved <- list()
    make_plot <- \(...) ggplot2::ggplot() + ggplot2::geom_blank()
    save_plot <- function(qcDir, filename, plotObject, width, height, ...) {
        saved[[filename]] <<- list(
            class = class(plotObject),
            width = width,
            height = height
        )
        file.path(qcDir, filename)
    }
    paths <- generateProtNormPostNormalizationQcArtifacts(
        object,
        qcDir = tempdir(),
        aesthetics = list(color_var = "group", shape_var = "batch"),
        qcPlotPaths = paths,
        messageFn = \(...) invisible(NULL),
        gcFn = \(...) invisible(NULL),
        plotPcaFn = make_plot,
        plotRleFn = make_plot,
        buildDensityFn = make_plot,
        buildCorrelationFn = make_plot,
        saveArtifactFn = save_plot
    )
    expect_identical(
        names(saved),
        c(
            "post_norm_pca.png",
            "post_norm_rle.png",
            "post_norm_density.png",
            "post_norm_correlation.png"
        )
    )
    expect_true(all(vapply(saved, \(item) "ggplot2::ggplot" %in% item$class, logical(1))))
    expect_identical(basename(paths$post_normalization$pca), "post_norm_pca.png")
    plot <- module_ci_prot_norm_cancor_plot()
    result <- list(
        best_cancor_plot = plot,
        optimization_results = data.frame(percentage = 1:2)
    )
    expect_identical(getProtNormRuvCanonicalCorrelationPlot(result), plot)
    prepared <- prepareProtNormOptimizationResultsTable(result)
    expect_true(prepared$hasResults)
    expect_identical(prepared$results, result$optimization_results)
    expect_null(prepared$bestPercentage)
})

test_that("normalization artifact helper is collated once before every caller", {
    description <- read.dcf(testthat::test_path("..", "..", "DESCRIPTION"))
    collate <- strsplit(description[[1L, "Collate"]], "[[:space:]]+")[[1L]]
    collate <- gsub("^'|'$", "", collate[nzchar(collate)])
    helper <- "mod_prot_norm_artifact_helpers.R"
    callers <- c(
        "mod_prot_norm_qc_support_helpers.R",
        "mod_prot_norm_workflow_helpers.R",
        "mod_prot_norm_ruv_helpers.R",
        "mod_prot_norm_correlation_helpers.R",
        "mod_prot_norm_session_helpers.R",
        "mod_prot_norm_observer_helpers.R"
    )
    expect_identical(sum(collate == helper), 1L)
    expect_true(all(match(helper, collate) < match(callers, collate)))
})
