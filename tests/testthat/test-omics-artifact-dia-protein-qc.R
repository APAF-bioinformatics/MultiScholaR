omics013SkipDependencies <- function() {
    for (package in c("DBI", "duckdb", "filelock", "arrow")) {
        testthat::skip_if_not_installed(package)
    }
}

omics013Capability <- function() {
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

omics013Context <- function(project_root, project_id = "omics-013-project") {
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
        capabilities = list(omics013Capability())
    )
    context
}

omics013Manager <- function(context) {
    manager <- newProtDiaArtifactStateManager(context)
    manager$setWorkflowType("DIA")
    manager
}

omics013Peptide <- function() {
    design <- data.frame(
        Run = c("S1", "S2", "S3", "S4"),
        group = factor(c("control", "control", "case", "case")),
        replicates = c("R1", "R2", "R1", "R2"),
        stringsAsFactors = FALSE
    )
    data <- data.frame(
        Protein.Group = rep(c("PG1", "PG2"), each = 4L),
        Protein.Ids = rep(c("P1;P1-2", "P2"), each = 4L),
        Stripped.Sequence = rep(c("PEPTIDE_A", "PEPTIDE_B"), each = 4L),
        Modified.Sequence = rep(c("_PEPTIDE_A_", "_PEPTIDE_B_"), each = 4L),
        Precursor.Id = rep(c("precursor-a", "precursor-b"), each = 4L),
        Precursor.Charge = 2L,
        Run = rep(design$Run, 2L),
        Q.Value = 0.001,
        Global.Q.Value = 0.002,
        Global.PG.Q.Value = 0.003,
        Proteotypic = 1L,
        Precursor.Quantity = c(10, NA, 30, 40, 50, 60, 70, 80),
        Precursor.Normalised = c(1, NA, 3, 4, 5, 6, 7, 8),
        Peptide.Imputed = c(1, NA, 3, 4, 5, 6, 7, 8),
        identification_peptide_count = 2L,
        identification_peptidoform_count = 2L,
        stringsAsFactors = FALSE
    )
    object <- PeptideQuantitativeDataDiann(
        data,
        design,
        sample_id = "Run",
        group_id = "group",
        technical_replicate_id = "replicates",
        args = list(source = "omics-art-013")
    )
    object@norm_quantity_column <- "Peptide.Imputed"
    object <- calcPeptideMatrix(object)
    object <- .ensurePeptideFeatureKeyMap(object, record_migration = FALSE)
    object@args$peptide_qc_audit <- list(
        records = list(list(record_id = "peptide-audit-final")),
        current_record_id = "peptide-audit-final"
    )
    object
}

omics013Protein <- function(peptide = omics013Peptide(), duplicate = FALSE) {
    groups <- if (duplicate) c("PG1", "PG1", "PG1", "PG2") else c("PG1", "PG2")
    values <- if (duplicate) {
        rbind(c(1, 2, 3, 4), c(3, 4, 5, 6), c(100, 200, 300, 400), c(8, 9, 10, 11))
    } else {
        rbind(c(2, 3, 4, 5), c(6, 7, 8, 9))
    }
    table <- data.frame(
        Protein.Group = groups,
        as.data.frame(values, check.names = FALSE),
        check.names = FALSE,
        stringsAsFactors = FALSE
    )
    names(table) <- c("Protein.Group", as.character(peptide@design_matrix$Run))
    ProteinQuantitativeData(
        protein_quant_table = table,
        protein_id_column = "Protein.Group",
        protein_id_table = .proteinIdTableFromPeptideLineage(peptide, groups),
        design_matrix = peptide@design_matrix,
        sample_id = "Run",
        group_id = "group",
        technical_replicate_id = "replicates",
        args = peptide@args
    )
}

omics013Workflow <- function(manager, context = NULL) {
    workflow <- new.env(parent = emptyenv())
    workflow$state_manager <- manager
    workflow$config_list <- list(globalParameters = list(workflow_type = "DIA"))
    workflow$qc_params <- list()
    if (!is.null(context)) workflow$workflow_context <- context
    workflow
}

omics013SaveInitial <- function(manager, name, object) {
    manager$saveState(
        name,
        object,
        list(stage = name),
        paste("Saved", name),
        audit_metadata = list(record_id = paste0("initial:", name))
    )
    invisible(object)
}

omics013Manifest <- function(context, manager, state_name = NULL) {
    if (is.null(state_name)) state_name <- manager$getCurrentStateName()
    row <- manager$states[[state_name]]
    store <- newArtifactStore(context$getPaths(), context$getIdentity()$project_id)
    artifactWorkflowStateReadManifest(store, row$manifest_relative_path)
}

omics013Metadata <- function(context, manager, state_name = NULL) {
    manifest <- omics013Manifest(context, manager, state_name)
    artifactWorkflowStateUnserializeMetadata(
        manifest$data$metadata_json,
        "OMICS-ART-013 test metadata"
    )
}

omics013ExpectScientificEqual <- function(expected, actual) {
    expect_identical(actual@protein_quant_table, expected@protein_quant_table)
    expect_identical(actual@protein_id_column, expected@protein_id_column)
    expect_identical(actual@protein_id_table, expected@protein_id_table)
    expect_identical(actual@design_matrix, expected@design_matrix)
    expect_identical(actual@sample_id, expected@sample_id)
    expect_identical(actual@group_id, expected@group_id)
    expect_identical(methods::validObject(actual, test = TRUE), TRUE)
}

omics013FailAt <- function(target_stage) {
    force(target_stage)
    \(stage, context) {
        if (identical(stage, target_stage)) {
            rlang::abort(
                paste("injected OMICS-ART-013 failure at", stage),
                class = "multischolar_test_artifact_failure"
            )
        }
        invisible(context)
    }
}

test_that("protein-QC artifacts preserve lineage and semantic representations", {
    omics013SkipDependencies()
    context <- omics013Context(withr::local_tempdir())
    manager <- omics013Manager(context)
    withr::defer(manager$close())
    workflow <- omics013Workflow(manager, context)
    peptide <- omics013Peptide()
    omics013SaveInitial(manager, "imputed", peptide)

    protein <- omics013Protein(peptide)
    current <- saveProtProteinQcState(
        workflow,
        manager,
        peptide,
        protein,
        "protein_rollup",
        "protein_s4_created",
        list(rollup_method = "iq_maxlfq"),
        "IQ MaxLFQ rollup",
        audit_parameters = list(
            rollup_method = "iq_maxlfq",
            normalization = "none",
            q_value_compatibility_value = 0.0009,
            pg_q_value_compatibility_value = 0.009,
            missing_intensity_compatibility_value = 0
        ),
        transformation_type = "aggregation",
        now = as.POSIXct("2026-08-21", tz = "UTC")
    )
    expect_null(omics013Metadata(context, manager)$derivation)
    expect_identical(resolveProteinQuantIdentityColumn(current), "Protein.Group")
    expect_true(all(c("Protein.Group", "Protein.Ids") %in%
        names(current@protein_id_table)))
    expect_identical(
        current@args$peptide_feature_key_map,
        peptide@args$peptide_feature_key_map
    )

    cleaned <- current
    cleaned@protein_quant_table$S1[[1L]] <- 20
    current <- saveProtProteinQcState(
        workflow,
        manager,
        current,
        cleaned,
        "protein_accession_cleanup",
        "protein_accession_cleaned",
        list(aggregation_method = "mean"),
        "Accession cleanup",
        transformation_type = "aggregation"
    )
    expect_null(omics013Metadata(context, manager)$derivation)

    filtered <- current
    filtered@protein_quant_table <- filtered@protein_quant_table[-1L, , drop = FALSE]
    filtered@protein_id_table <- filtered@protein_id_table[
        filtered@protein_id_table$Protein.Group == "PG2",
        ,
        drop = FALSE
    ]
    current <- saveProtProteinQcState(
        workflow,
        manager,
        current,
        filtered,
        "protein_intensity_filter",
        "protein_intensity_filtered",
        list(strict_mode = TRUE),
        "Intensity filter"
    )
    expect_identical(
        omics013Metadata(context, manager)$derivation$representation,
        "row_selection"
    )

    replicated <- saveProtProteinQcState(
        workflow,
        manager,
        current,
        current,
        "protein_replicate_filter",
        "protein_replicate_filtered",
        list(grouping_variable = "group"),
        "Replicate filter"
    )
    metadata <- omics013Metadata(context, manager)
    expect_identical(metadata$derivation$representation, "row_selection")
    expect_identical(manager$getCacheInfo()$entries, 1L)
    expect_identical(manager$getState(), replicated)
    expect_identical(methods::validObject(replicated, test = TRUE), TRUE)

    all_failed <- replicated
    all_failed@protein_quant_table <- all_failed@protein_quant_table[0L, , drop = FALSE]
    all_failed@protein_id_table <- all_failed@protein_id_table[0L, , drop = FALSE]
    all_failed <- saveProtProteinQcState(
        workflow,
        manager,
        replicated,
        all_failed,
        "protein_intensity_filter",
        "protein_all_failed",
        list(strict_mode = TRUE),
        "All-fail intensity filter"
    )
    all_fail_metadata <- omics013Metadata(context, manager)
    expect_identical(all_fail_metadata$derivation$selected_count, 0L)
    expect_identical(all_fail_metadata$derivation$rejected_count, 1L)
    expect_identical(manager$getCacheInfo()$entries, 1L)

    audit <- all_failed@args$protein_qc_audit
    expect_identical(length(audit$records), 5L)
    expect_identical(
        vapply(audit$records, `[[`, character(1), "stage_id"),
        c(
            "protein_rollup",
            "protein_accession_cleanup",
            "protein_intensity_filter",
            "protein_replicate_filter",
            "protein_intensity_filter"
        )
    )
    expect_identical(
        audit$records[[1L]]$peptide_audit_record_id,
        "peptide-audit-final"
    )
    expect_identical(
        audit$records[[5L]]$protein_id_table_digest,
        .peptideQcDigest(all_failed@protein_id_table)
    )
    expect_identical(
        audit$records[[5L]]$parent_record_id,
        audit$records[[4L]]$record_id
    )
    result <- workflow$artifact_stage_results$protein_intensity_filter
    expect_gt(result$metrics$new_artifact_bytes, 0)
    expect_gt(result$metrics$hydrated_object_bytes, 0)
    expect_false(result$compaction$enabled)
})

test_that("duplicate mean median and max always materialize exact values", {
    omics013SkipDependencies()
    expected <- c(mean = 104 / 3, median = 3, max = 100)
    for (method in names(expected)) {
        root <- withr::local_tempdir()
        context <- omics013Context(root, paste0("omics-013-", method))
        manager <- omics013Manager(context)
        workflow <- omics013Workflow(manager, context)
        parent <- omics013Protein(duplicate = TRUE)
        omics013SaveInitial(manager, "protein_s4_created", parent)

        result <- runProteinDuplicateRemovalStep(
            workflow,
            aggregationMethod = method,
            logInfoFn = \(...) invisible(NULL)
        )
        pg1 <- result$deduplicatedS4@protein_quant_table$Protein.Group == "PG1"
        expect_equal(
            result$deduplicatedS4@protein_quant_table$S1[pg1],
            unname(expected[[method]]),
            info = method
        )
        expect_null(omics013Metadata(context, manager)$derivation, info = method)
        expect_identical(manager$getState(), result$deduplicatedS4)
        manager$close()
    }
})

test_that("IQ and limpa rollups retain scientific parity in artifact mode", {
    omics013SkipDependencies()
    peptide <- omics013Peptide()
    paths <- list(
        peptide_qc_dir = withr::local_tempdir(),
        protein_qc_dir = withr::local_tempdir(),
        source_dir = withr::local_tempdir()
    )
    iq_observed <- new.env(parent = emptyenv())
    process_iq <- function(
        input_filename,
        output_filename,
        sample_id,
        primary_id,
        secondary_id,
        intensity_col,
        filter_double_less,
        normalization
    ) {
        input <- readr::read_tsv(input_filename, show_col_types = FALSE)
        iq_observed$input <- input
        iq_observed$parameters <- list(
            sample_id = sample_id,
            primary_id = primary_id,
            secondary_id = secondary_id,
            intensity_col = intensity_col,
            filter_double_less = filter_double_less,
            normalization = normalization
        )
        output <- data.frame(
            Protein.Group = c("PG1", "PG2"),
            S_001 = c(2, 6),
            S_002 = c(0, 7),
            S_003 = c(4, 8),
            S_004 = c(5, 9),
            check.names = FALSE
        )
        readr::write_tsv(output, output_filename)
    }
    run_iq <- function(manager, workflow) {
        omics013SaveInitial(manager, "imputed", peptide)
        runProteinIqRollupApplyStep(
            workflow,
            experimentPaths = paths,
            processLongFormatFn = process_iq,
            readTsvFn = \(file, .name_repair) {
                readr::read_tsv(file, show_col_types = FALSE)
            },
            captureCheckpointFn = \(...) invisible(NULL),
            showNotificationFn = \(...) invisible(NULL),
            logInfoFn = \(...) invisible(NULL),
            logWarnFn = \(...) invisible(NULL),
            sleepFn = \(...) invisible(NULL)
        )$proteinObj
    }

    withr::local_options(multischolar.prot_dia.protein_qc_artifacts = "disabled")
    memory_manager <- newWorkflowState()
    memory_manager$setWorkflowType("DIA")
    memory <- run_iq(memory_manager, omics013Workflow(memory_manager))
    context <- omics013Context(withr::local_tempdir(), "omics-013-iq")
    artifact_manager <- omics013Manager(context)
    withr::defer(artifact_manager$close())
    artifact <- run_iq(
        artifact_manager,
        omics013Workflow(artifact_manager, context)
    )
    omics013ExpectScientificEqual(memory, artifact)
    expect_identical(artifact@protein_id_column, "Protein.Group")
    expect_identical(unique(artifact@protein_id_table$Protein.Ids), c("P1;P1-2", "P2"))
    expect_true(all(iq_observed$input$Q.Value == 0.0009))
    expect_true(all(iq_observed$input$PG.Q.Value == 0.009))
    expect_false(anyNA(iq_observed$input$Peptide.Imputed))
    expect_identical(iq_observed$parameters$normalization, "none")
    iq_audit <- protDiaProteinQcCurrentRecord(artifact)
    expect_identical(iq_audit$resolved_parameters$rollup_method, "iq_maxlfq")
    expect_true(nzchar(iq_audit$software$scientific_engine$version))

    limpa_fn <- function(peptide_object, verbose) {
        object <- omics013Protein(peptide_object)
        object@args$proteinMissingValueImputationLimpa <- list(dpc_slope = 0.8)
        object@args$limpa_dpc_quant_results <- list(
            dpc_parameters_used = c(0.2, 0.8),
            quantification_method = "limpa_dpc_quant",
            dpc_method = "limpa_dpc_quant",
            missing_percentage_peptides = 12.5,
            missing_percentage_proteins = 0,
            total_proteins_quantified = 2L,
            total_peptides_used = 8L
        )
        object
    }
    run_limpa <- function(manager, workflow) {
        omics013SaveInitial(manager, "imputed", peptide)
        runProteinLimpaRollupApplyStep(
            workflow,
            experimentPaths = paths,
            proteinLimpaFn = limpa_fn,
            requireNamespaceFn = \(...) TRUE,
            captureCheckpointFn = \(...) invisible(NULL),
            logInfoFn = \(...) invisible(NULL)
        )$proteinObj
    }
    limpa_memory_manager <- newWorkflowState()
    limpa_memory_manager$setWorkflowType("DIA")
    limpa_memory <- run_limpa(
        limpa_memory_manager,
        omics013Workflow(limpa_memory_manager)
    )
    limpa_context <- omics013Context(withr::local_tempdir(), "omics-013-limpa")
    limpa_artifact_manager <- omics013Manager(limpa_context)
    withr::defer(limpa_artifact_manager$close())
    limpa_artifact <- run_limpa(
        limpa_artifact_manager,
        omics013Workflow(limpa_artifact_manager, limpa_context)
    )
    omics013ExpectScientificEqual(limpa_memory, limpa_artifact)
    limpa_audit <- protDiaProteinQcCurrentRecord(limpa_artifact)
    expect_identical(limpa_audit$resolved_parameters$dpc_parameters_used, c(0.2, 0.8))
    expect_identical(limpa_audit$software$scientific_engine$name, "limpa")
})

test_that("protein artifact failures and reverts preserve the active generation", {
    omics013SkipDependencies()
    context <- omics013Context(withr::local_tempdir())
    manager <- omics013Manager(context)
    withr::defer(manager$close())
    workflow <- omics013Workflow(manager, context)
    parent <- omics013Protein()
    omics013SaveInitial(manager, "protein_s4_created", parent)
    parent_generation <- manager$getCurrentGenerationId()
    child <- parent
    child@protein_quant_table <- child@protein_quant_table[-1L, , drop = FALSE]
    child@protein_id_table <- child@protein_id_table[
        child@protein_id_table$Protein.Group == "PG2",
        ,
        drop = FALSE
    ]

    for (failure_stage in c(
        "after_audit_creation",
        "before_state_registry_commit"
    )) {
        expect_error(
            saveProtProteinQcState(
                workflow,
                manager,
                parent,
                child,
                "protein_intensity_filter",
                paste0("failed_", failure_stage),
                list(strict_mode = TRUE),
                "Failed filter",
                failure_injector = omics013FailAt(failure_stage)
            ),
            class = "multischolar_test_artifact_failure"
        )
        expect_identical(manager$getCurrentGenerationId(), parent_generation)
        expect_identical(manager$getCurrentStateName(), "protein_s4_created")
        expect_identical(manager$getState(), parent)
    }

    saved <- saveProtProteinQcState(
        workflow,
        manager,
        parent,
        child,
        "protein_intensity_filter",
        "protein_intensity_filtered",
        list(strict_mode = TRUE),
        "Intensity filter"
    )
    expect_identical(manager$getState(), saved)
    reverted <- revertProtDiaProteinQcState(workflow, "protein_s4_created")
    expect_identical(reverted, parent)
    manager$close()
    resumed <- omics013Manager(context)
    withr::defer(resumed$close())
    expect_identical(resumed$getCurrentStateName(), "protein_s4_created")
    expect_identical(resumed$getState(), parent)
})

test_that("artifact protein-QC context remains workflow-session scoped", {
    omics013SkipDependencies()
    context <- omics013Context(withr::local_tempdir())
    manager <- omics013Manager(context)
    withr::defer(manager$close())
    workflow <- omics013Workflow(manager, context)
    paths <- list(
        source_dir = withr::local_tempdir(),
        peptide_qc_dir = withr::local_tempdir(),
        protein_qc_dir = withr::local_tempdir(),
        publication_graphs_dir = withr::local_tempdir(),
        time_dir = withr::local_tempdir()
    )
    for (name in c("filtering_progress", "project_dirs", "experiment_paths")) {
        if (exists(name, envir = .GlobalEnv, inherits = FALSE)) {
            old <- get(name, envir = .GlobalEnv, inherits = FALSE)
            withr::defer(assign(name, old, envir = .GlobalEnv))
            rm(list = name, envir = .GlobalEnv)
        }
    }
    initializeProtDiaProteinQcSessionContext(
        workflow,
        paths,
        "proteomics",
        "study"
    )
    protein <- omics013Protein()
    omics013SaveInitial(manager, "protein_s4_created", protein)
    workflow$aa_seq_tbl_final <- data.frame(database_id = c("P1", "P2"))
    cleanup <- runProteinAccessionCleanupStep(
        workflow,
        delimiter = ";",
        aggregationMethod = "mean",
        chooseBestProteinAccessionFn = \(theObject, ...) theObject,
        logInfoFn = \(...) invisible(NULL),
        logWarnFn = \(...) invisible(NULL),
        saveRdsFn = \(...) invisible(NULL),
        existsFn = \(...) stop("artifact cleanup consulted the global fallback")
    )
    expect_true(cleanup$cleanupApplied)
    observed <- NULL
    update_fn <- function(
        data,
        step_name,
        omic_type,
        experiment_label,
        return_grid,
        overwrite,
        project_dirs,
        progress_env
    ) {
        observed <<- list(
            project_dirs = project_dirs,
            progress_env = progress_env,
            step_name = step_name
        )
        "plot-grid"
    }
    result <- protDiaProteinQcUpdateFiltering(
        update_fn,
        workflow,
        data.frame(Protein.Group = "PG1", S1 = 1),
        "protein_s4_created",
        "proteomics",
        "study",
        TRUE,
        TRUE
    )
    expect_identical(result, "plot-grid")
    expect_identical(observed$project_dirs$proteomics_study, paths)
    expect_identical(observed$progress_env, workflow$protein_qc_context$progress_env)
    expect_identical(
        resolveLimpaQcSaveDir("peptide_qc", observed$project_dirs, FALSE),
        paths$peptide_qc_dir
    )
    expect_true("project_dirs" %in% names(formals(generateLimpaQCPlots)))
    expect_false(exists("filtering_progress", envir = .GlobalEnv, inherits = FALSE))
    expect_false(exists("project_dirs", envir = .GlobalEnv, inherits = FALSE))
    expect_false(exists("experiment_paths", envir = .GlobalEnv, inherits = FALSE))
})
