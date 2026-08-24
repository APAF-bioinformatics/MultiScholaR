nondiaArtifact025SkipDependencies <- function() {
    for (package in c("arrow", "DBI", "duckdb", "filelock", "ps")) {
        testthat::skip_if_not_installed(package)
    }
}

nondiaArtifact025RepoPath <- function(...) {
    file.path(normalizePath(testthat::test_path("..", ".."), mustWork = TRUE), ...)
}

nondiaArtifact025Scenarios <- function() {
    manifest <- jsonlite::read_json(
        nondiaArtifact025RepoPath(
            "tests", "testdata", "omics-parity", "proteomics", "manifest.json"
        ),
        simplifyVector = FALSE
    )
    manifest$fixture_scenarios
}

nondiaArtifact025Importer <- function(format) {
    switch(
        format,
        maxquant = importMaxQuantData,
        fragpipe = importFragPipeData,
        pd_tmt = importProteomeDiscovererTMTData
    )
}

nondiaArtifact025WorkflowType <- function(format) {
    if (identical(format, "pd_tmt")) "TMT" else "LFQ"
}

nondiaArtifact025CapabilityId <- function(format) {
    paste0(
        "proteomics.",
        switch(format, maxquant = "maxquant", fragpipe = "fragpipe", pd_tmt = "pd_tmt"),
        ".protein.",
        if (identical(format, "pd_tmt")) "tmt" else "lfq",
        ".v1"
    )
}

nondiaArtifact025VendorObject <- function(scenario) {
    imported <- suppressMessages(
        nondiaArtifact025Importer(scenario$format)(
            nondiaArtifact025RepoPath(scenario$fixture_path)
        )
    )
    data <- as.data.frame(imported$data)
    mapping <- imported$column_mapping
    runs <- unique(as.character(data[[mapping$run_col]]))
    groups <- ifelse(grepl("KO", runs, fixed = TRUE), "KO", "WT")
    replicates <- ave(seq_along(runs), groups, FUN = seq_along)
    workflow_type <- nondiaArtifact025WorkflowType(scenario$format)
    manager <- new.env(parent = emptyenv())
    manager$saveState <- \(state_name, s4_data_object, ...) {
        manager$object <- s4_data_object
        invisible(state_name)
    }
    workflow <- list2env(list(
        data_cln = data,
        column_mapping = mapping,
        design_matrix = data.frame(
            Run = runs,
            group = groups,
            replicates = paste0("R", replicates),
            stringsAsFactors = FALSE
        ),
        config_list = list(
            globalParameters = list(workflow_type = workflow_type)
        ),
        state_manager = manager
    ), parent = emptyenv())
    suppressMessages(buildProtDesignStateCheckpoint(
        workflow,
        workflow_type,
        "non-DIA normalization fixture",
        validateColumnMapping = TRUE
    ))
    manager$object
}

nondiaArtifact025Context <- function(
    root,
    descriptor,
    project_id
) {
    context <- createWorkflowContext(
        list(
            base_dir = root,
            project_id = project_id,
            omic_label = "normalization-study"
        ),
        "proteomics",
        "normalization-study",
        storage_policy = list(
            requested_backend = "artifact",
            requested_rollout = "dual_write",
            migration_requested = TRUE,
            project_id = project_id
        )
    )
    identity <- descriptor$identity
    bindWorkflowContextFromImport(
        context,
        workflow_type = identity$workflow_type,
        input_format = identity$input_format,
        data_level = identity$data_level,
        descriptor_catalogue = artifactWorkflowDescriptorCatalogue()
    )
    context
}

nondiaArtifact025Manager <- function(context, descriptor, initial) {
    manager <- newWorkflowState(
        workflow_context = context,
        workflow_descriptor = descriptor,
        descriptor_catalogue = artifactWorkflowDescriptorCatalogue(),
        codec_catalogue = artifactS4CodecCatalogue()
    )
    manager$setWorkflowType(descriptor$identity$workflow_type)
    manager$saveState(
        "protein_s4_initial",
        initial,
        initial@args,
        "Initial non-DIA protein state",
        audit_metadata = list(record_id = "initial:protein_s4_initial")
    )
    manager
}

nondiaArtifact025Workflow <- function(manager, context, workflow_type) {
    workflow <- new.env(parent = emptyenv())
    workflow$state_manager <- manager
    workflow$workflow_context <- context
    workflow$config_list <- list(
        globalParameters = list(workflow_type = workflow_type)
    )
    workflow$artifact_stage_results <- list()
    workflow$normalization_state_refs <- list()
    workflow$ruv_normalised_for_da_analysis_obj <- NULL
    workflow$final_for_da_ref <- NULL
    workflow$protein_counts <- list()
    workflow$tab_status <- list(
        normalization = "pending",
        differential_expression = "disabled"
    )
    workflow$normalization_context <- list(
        normalization_method = "none",
        ruv_mode = "manual",
        ruv_grouping_variable = "group",
        ruv_parameters = list()
    )
    workflow
}

nondiaArtifact025NormData <- function(manager) {
    norm_data <- new.env(parent = emptyenv())
    norm_data$state_manager <- manager
    norm_data$state_refs <- list()
    norm_data$normalized_protein_obj <- NULL
    norm_data$ruv_normalized_obj <- NULL
    norm_data$correlation_filtered_obj <- NULL
    norm_data
}

nondiaArtifact025ExpectScientificEqual <- function(expected, actual) {
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

nondiaArtifact025SaveNormalized <- function(
    workflow,
    norm_data,
    before,
    normalized,
    method,
    ruv_mode = "manual",
    failure_injector = NULL
) {
    parameters <- list(
        normalization_method = method,
        ruv_mode = ruv_mode,
        ruv_grouping_variable = "group",
        skip_reason = if (identical(ruv_mode, "skip")) {
            "dataset constraints"
        } else {
            NULL
        },
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
        transformation_type = "normalization",
        failure_injector = failure_injector
    )
    settleProtNormArtifactState(
        workflow,
        norm_data,
        "normalization",
        "normalized",
        saved
    )
}

nondiaArtifact025Metadata <- function(context, manager, state_name) {
    row <- manager$states[[state_name]]
    store <- newArtifactStore(context$getPaths(), context$getIdentity()$project_id)
    manifest <- artifactWorkflowStateReadManifest(
        store,
        row$manifest_relative_path
    )
    artifactWorkflowStateUnserializeMetadata(
        manifest$data$metadata_json,
        "OMICS-ART-025 test metadata"
    )
}

nondiaArtifact025ModuleObject <- function(descriptor) {
    object <- module_ci_prot_norm_object()
    object@args$globalParameters$workflow_type <-
        descriptor$identity$workflow_type
    object
}

nondiaArtifact025RunChain <- function(
    manager,
    workflow,
    norm_data,
    descriptor,
    ruv_mode = "manual",
    correlation_skipped = FALSE
) {
    methods <- module_ci_prot_norm_s4_methods()
    before <- manager$getState()
    normalized <- methods$normalise(before, normalisation_method = "none")
    normalized <- nondiaArtifact025SaveNormalized(
        workflow,
        norm_data,
        before,
        normalized,
        "none",
        ruv_mode
    )
    controls <- module_ci_prot_norm_control_index(
        normalized@protein_quant_table$Protein.Ids
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
            control_genes_index = controls,
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
            optimization_results = optimization_results
        )
        correlation_input <- corrected
    } else {
        workflow$ruv_optimization_result <- buildProtNormSkippedRuvResult()
        correlation_input <- normalized
    }
    norm_data$ruv_normalized_obj <- correlation_input
    pairs <- module_ci_prot_norm_correlation_pairs("pass_all")
    filtered <- filterSamplesByProteinCorrelationThreshold(
        correlation_input,
        pearson_correlation_per_pair = pairs,
        min_pearson_correlation_threshold = 0
    )
    norm_data$correlation_filtered_obj <- filtered
    finalizeProtNormCorrelationWorkflowState(
        filtered,
        workflow,
        normData = norm_data,
        correlationThreshold = 0,
        skipped = correlation_skipped,
        timeFn = \() as.POSIXct("2026-08-24 00:00:00", tz = "UTC"),
        messageFn = \(...) invisible(NULL),
        catFn = \(...) invisible(NULL)
    )
    manager$getState("correlation_filtered")
}

test_that("all non-DIA normalization engines retain exact artifact parity", {
    nondiaArtifact025SkipDependencies()
    methods <- module_ci_prot_norm_s4_methods()
    for (scenario in nondiaArtifact025Scenarios()) {
        descriptor <- protNonDiaReadthroughDescriptor(
            nondiaArtifact025CapabilityId(scenario$format)
        )
        for (method in c("none", "cyclicloess", "quantile", "scale")) {
            before <- nondiaArtifact025VendorObject(scenario)
            expected <- methods$normalise(
                before,
                normalisation_method = method
            )
            root <- withr::local_tempdir(
                pattern = paste0("nondia-025-", scenario$format, "-", method, "-")
            )
            context <- nondiaArtifact025Context(
                root,
                descriptor,
                paste0("nondia-025-", scenario$format, "-", method)
            )
            manager <- nondiaArtifact025Manager(context, descriptor, before)
            workflow <- nondiaArtifact025Workflow(
                manager,
                context,
                descriptor$identity$workflow_type
            )
            norm_data <- nondiaArtifact025NormData(manager)
            rss_before <- as.numeric(ps::ps_memory_info(ps::ps_handle())[["rss"]])
            actual <- nondiaArtifact025SaveNormalized(
                workflow,
                norm_data,
                before,
                expected,
                method
            )
            rss_after <- as.numeric(ps::ps_memory_info(ps::ps_handle())[["rss"]])
            nondiaArtifact025ExpectScientificEqual(expected, actual)
            nondiaArtifact025ExpectScientificEqual(
                actual,
                manager$getState("normalized")
            )
            record <- protNonDiaNormCurrentRecord(actual)
            expect_identical(record$capability_id, descriptor$descriptor_id)
            expect_identical(record$stage_id, "normalization")
            expect_identical(
                record$resolved_parameters$normalization_method,
                method
            )
            expect_identical(
                workflow$normalization_state_refs$normalization$state_name,
                "normalized"
            )
            expect_true(all(is.finite(c(rss_before, rss_after))))
            manager$close()
        }
    }
})

test_that("non-DIA reachable normalization state maps and refs stay exact", {
    nondiaArtifact025SkipDependencies()
    for (descriptor in protNonDiaReadthroughDescriptors()) {
        for (ruv_mode in c("manual", "automatic", "skip")) {
            root <- withr::local_tempdir(
                pattern = paste0("nondia-025-chain-", descriptor$identity$input_format, "-")
            )
            before <- nondiaArtifact025ModuleObject(descriptor)
            context <- nondiaArtifact025Context(
                root,
                descriptor,
                paste0("nondia-025-chain-", descriptor$identity$input_format, "-", ruv_mode)
            )
            manager <- nondiaArtifact025Manager(context, descriptor, before)
            workflow <- nondiaArtifact025Workflow(
                manager,
                context,
                descriptor$identity$workflow_type
            )
            norm_data <- nondiaArtifact025NormData(manager)
            final <- nondiaArtifact025RunChain(
                manager,
                workflow,
                norm_data,
                descriptor,
                ruv_mode,
                correlation_skipped = identical(ruv_mode, "skip")
            )
            expected_history <- if (identical(ruv_mode, "skip")) {
                c("initial", "protein_s4_initial", "normalized", "correlation_filtered")
            } else {
                c(
                    "initial", "protein_s4_initial", "normalized",
                    "ruv_corrected", "correlation_filtered"
                )
            }
            expect_identical(manager$getHistory(), expected_history)
            expect_false("protein_qc" %in% manager$getHistory())
            expect_false("peptide_qc" %in% manager$getHistory())
            stages <- vapply(
                final@args$normalization_artifact_audit$records,
                `[[`,
                character(1),
                "stage_id"
            )
            expect_identical(
                stages,
                if (identical(ruv_mode, "skip")) {
                    c("normalization", "correlation_filter")
                } else {
                    .PROT_DIA_NORM_STAGES
                }
            )
            expect_identical(
                workflow$final_for_da_ref$state_name,
                "correlation_filtered"
            )
            expect_identical(
                workflow$final_for_da_ref$capability_id,
                descriptor$descriptor_id
            )
            expect_null(norm_data$normalized_protein_obj)
            expect_null(norm_data$ruv_normalized_obj)
            metadata <- if (!identical(ruv_mode, "skip")) {
                nondiaArtifact025Metadata(context, manager, "ruv_corrected")
            } else {
                list(payloads = list())
            }
            if (!identical(ruv_mode, "skip")) {
                owners <- vapply(metadata$payloads, `[[`, character(1), "owner")
                expect_true(any(grepl(
                    "optimization_results",
                    owners,
                    fixed = TRUE
                )))
            }
            correlation_metadata <- nondiaArtifact025Metadata(
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
            current_record <- protNonDiaNormCurrentRecord(final)
            expect_identical(
                current_record$resolved_parameters$output_files,
                if (identical(ruv_mode, "skip")) {
                    c(
                        "normalised_results_cln_with_replicates.tsv",
                        "normalised_results_cln_with_replicates.RDS"
                    )
                } else {
                    c(
                        "ruv_normalised_results_cln_with_replicates.tsv",
                        "ruv_normalised_results_cln_with_replicates.RDS"
                    )
                }
            )
            expect_identical(
                revertProtNormState(workflow, "normalized")@protein_quant_table,
                manager$getState("normalized")@protein_quant_table
            )
            manager$close()
        }
    }
})

test_that("non-DIA normalization failures, switches, and globals fail closed", {
    nondiaArtifact025SkipDependencies()
    descriptor <- protNonDiaReadthroughDescriptors()[[1L]]
    before <- nondiaArtifact025ModuleObject(descriptor)
    root <- withr::local_tempdir(pattern = "nondia-025-failure-")
    context <- nondiaArtifact025Context(root, descriptor, "nondia-025-failure")
    artifact_manager <- nondiaArtifact025Manager(context, descriptor, before)
    parent_generation <- artifact_manager$getCurrentGenerationId()
    artifact_manager$close()
    artifact_manager <- NULL
    memory_manager <- WorkflowState$new()
    memory_manager$setWorkflowType(descriptor$identity$workflow_type)
    memory_manager$saveState(
        "protein_s4_initial",
        before,
        before@args,
        "Initial memory state"
    )
    workflow <- nondiaArtifact025Workflow(
        memory_manager,
        context,
        descriptor$identity$workflow_type
    )
    norm_data <- nondiaArtifact025NormData(memory_manager)
    normalized <- module_ci_prot_norm_s4_methods()$normalise(
        before,
        normalisation_method = "none"
    )
    expect_error(
        nondiaArtifact025SaveNormalized(
            workflow,
            norm_data,
            before,
            normalized,
            "none",
            failure_injector = \(stage, ...) {
                if (identical(stage, "before_state_registry_commit")) {
                    rlang::abort(
                        "injected normalization failure",
                        class = "multischolar_test_artifact_failure"
                    )
                }
            }
        ),
        class = "multischolar_test_artifact_failure"
    )
    artifact_manager <- newWorkflowState(
        workflow_context = context,
        workflow_descriptor = descriptor,
        descriptor_catalogue = artifactWorkflowDescriptorCatalogue(),
        codec_catalogue = artifactS4CodecCatalogue()
    )
    expect_identical(artifact_manager$getCurrentGenerationId(), parent_generation)
    expect_identical(memory_manager$getCurrentStateName(), "protein_s4_initial")

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
    expect_identical(
        workflow$config_list$normaliseBetweenSamples$method,
        "quantile"
    )
    expect_identical(get("config_list", envir = .GlobalEnv), hostile)
    updateProtNormRuvAuditTrail(
        ruvK = 2L,
        controlGenesIndex = c(TRUE, FALSE),
        percentageAsNegCtrl = 40,
        modeLabel = "manual",
        workflowData = workflow,
        updateRuvParametersFn = \(config, k, controls, percentage) {
            config$ruv <- list(
                k = k,
                controls = controls,
                percentage = percentage
            )
            config
        },
        messageFn = \(...) invisible(NULL)
    )
    expect_identical(workflow$config_list$ruv$k, 2L)
    expect_identical(get("config_list", envir = .GlobalEnv), hostile)

    withr::with_options(list(
        multischolar.prot_nondia.normalization_artifacts.normalization =
            "disabled"
    ), {
        child <- nondiaArtifact025SaveNormalized(
            workflow,
            norm_data,
            before,
            normalized,
            "none"
        )
        nondiaArtifact025ExpectScientificEqual(normalized, child)
        expect_identical(
            artifact_manager$getCurrentGenerationId(),
            parent_generation
        )
    })

    memory_only <- WorkflowState$new()
    memory_only$setWorkflowType("LFQ")
    memory_only$saveState("protein_s4_initial", before, before@args, "memory")
    memory_workflow <- new.env(parent = emptyenv())
    memory_workflow$state_manager <- memory_only
    memory_workflow$config_list <- before@args
    memory_workflow$artifact_stage_results <- list()
    saved_memory <- saveProtNormState(
        memory_workflow,
        memory_only,
        before,
        normalized,
        "normalization",
        "normalized",
        list(method = "none"),
        "memory normalization"
    )
    expect_identical(saved_memory, normalized)
    expect_identical(memory_only$getState(), normalized)
    artifact_manager$close()
})

test_that("moved normalized states hydrate exact and compatibility failure rolls back", {
    nondiaArtifact025SkipDependencies()
    for (descriptor in protNonDiaReadthroughDescriptors()) {
        original <- withr::local_tempdir(
            pattern = paste0("nondia-025-move-", descriptor$identity$input_format, "-")
        )
        moved <- withr::local_tempdir(
            pattern = paste0("nondia-025-moved-", descriptor$identity$input_format, "-")
        )
        project_id <- paste0("nondia-025-move-", descriptor$identity$input_format)
        before <- nondiaArtifact025ModuleObject(descriptor)
        context <- nondiaArtifact025Context(original, descriptor, project_id)
        manager <- nondiaArtifact025Manager(context, descriptor, before)
        workflow <- nondiaArtifact025Workflow(
            manager,
            context,
            descriptor$identity$workflow_type
        )
        normalized <- module_ci_prot_norm_s4_methods()$normalise(
            before,
            normalisation_method = "none"
        )
        normalized <- nondiaArtifact025SaveNormalized(
            workflow,
            nondiaArtifact025NormData(manager),
            before,
            normalized,
            "none"
        )
        manager$close()
        entries <- list.files(
            original,
            all.files = TRUE,
            no.. = TRUE,
            full.names = TRUE
        )
        expect_true(all(file.copy(entries, moved, recursive = TRUE)))
        unlink(original, recursive = TRUE, force = TRUE)
        moved_context <- nondiaArtifact025Context(
            moved,
            descriptor,
            project_id
        )
        moved_manager <- newWorkflowState(
            workflow_context = moved_context,
            workflow_descriptor = descriptor,
            descriptor_catalogue = artifactWorkflowDescriptorCatalogue(),
            codec_catalogue = artifactS4CodecCatalogue()
        )
        nondiaArtifact025ExpectScientificEqual(
            normalized,
            moved_manager$getState("normalized")
        )
        expect_identical(
            moved_manager$getCurrentStateName(),
            "normalized"
        )
        moved_manager$close()

        rollback_root <- withr::local_tempdir(
            pattern = paste0("nondia-025-compat-", descriptor$identity$input_format, "-")
        )
        rollback_context <- nondiaArtifact025Context(
            rollback_root,
            descriptor,
            paste0(project_id, "-compat")
        )
        setup_manager <- nondiaArtifact025Manager(
            rollback_context,
            descriptor,
            before
        )
        parent_generation <- setup_manager$getCurrentGenerationId()
        setup_manager$close()
        memory_manager <- WorkflowState$new()
        memory_manager$setWorkflowType(descriptor$identity$workflow_type)
        memory_manager$saveState(
            "protein_s4_initial",
            before,
            before@args,
            "memory parent"
        )
        failing_manager <- new.env(parent = emptyenv())
        failing_manager$getWorkflowType <- memory_manager$getWorkflowType
        failing_manager$getCurrentStateName <- memory_manager$getCurrentStateName
        failing_manager$getState <- memory_manager$getState
        failing_manager$hasState <- memory_manager$hasState
        failing_manager$revertToState <- memory_manager$revertToState
        failing_manager$saveState <- \(...) stop("compatibility save failure")
        rollback_workflow <- nondiaArtifact025Workflow(
            failing_manager,
            rollback_context,
            descriptor$identity$workflow_type
        )
        expect_error(
            saveProtNormState(
                rollback_workflow,
                failing_manager,
                before,
                normalized,
                "normalization",
                "normalized",
                list(normalization_method = "none"),
                "compatibility failure"
            ),
            "compatibility save failure",
            fixed = TRUE
        )
        inspection <- newWorkflowState(
            workflow_context = rollback_context,
            workflow_descriptor = descriptor,
            descriptor_catalogue = artifactWorkflowDescriptorCatalogue(),
            codec_catalogue = artifactS4CodecCatalogue()
        )
        expect_identical(inspection$getCurrentGenerationId(), parent_generation)
        expect_identical(inspection$getCurrentStateName(), "protein_s4_initial")
        inspection$close()
    }
})
