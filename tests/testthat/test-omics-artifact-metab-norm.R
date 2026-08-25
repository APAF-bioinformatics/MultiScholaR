metab034SkipDependencies <- function() {
    for (package in c("DBI", "duckdb", "filelock", "arrow", "limma")) {
        testthat::skip_if_not_installed(package)
    }
}

metab034Context <- function(root, project_id) {
    paths <- list(
        base_dir = root,
        project_id = project_id,
        omic_type = "metabolomics",
        omic_label = "metabolomics-study",
        source_dir = file.path(root, "source"),
        results_dir = file.path(root, "results")
    )
    dir.create(paths$source_dir, recursive = TRUE, showWarnings = FALSE)
    dir.create(paths$results_dir, recursive = TRUE, showWarnings = FALSE)
    capability <- workflowCapabilityCatalogue()[[
        "metabolomics.custom.metabolite.standard.v1"
    ]]
    capability$artifact_eligible <- TRUE
    capability$auto_eligible <- FALSE
    capability$maximum_artifact_rollout <- "dual_write"
    context <- createWorkflowContext(
        paths,
        "metabolomics",
        "metabolomics-study",
        storage_policy = list(
            requested_backend = "artifact",
            requested_rollout = "dual_write",
            migration_requested = TRUE,
            project_id = project_id
        )
    )
    bindWorkflowContextFromImport(
        context,
        workflow_type = "metabolomics_standard",
        input_format = "custom",
        data_level = "metabolite",
        capabilities = list(capability)
    )
    context
}

metab034Manager <- function(context) {
    manager <- ArtifactWorkflowState$new(
        workflow_context = context,
        dehydrate_fn = dehydrateMetabolomicsS4Artifact,
        validate_bundle_fn = validateMetabolomicsS4Bundle,
        hydrate_fn = hydrateMetabolomicsS4Artifact,
        descriptor_contract = artifactStageDescriptorContract(
            artifactMetabolomicsWorkflowDescriptor()
        )
    )
    manager$setWorkflowType("metabolomics_standard")
    manager
}

metab034Workflow <- function(manager, context) {
    workflow <- new.env(parent = emptyenv())
    workflow$state_manager <- manager
    workflow$workflow_context <- context
    workflow$config_list <- list(
        globalParameters = list(workflow_type = "metabolomics_standard"),
        normalization = list(mode = "artifact-ticket")
    )
    workflow$artifact_stage_results <- list()
    workflow$processing_log <- list()
    workflow$tab_status <- list(
        quality_control = "complete",
        normalization = "pending",
        differential_analysis = "locked"
    )
    workflow$data_format <- "custom"
    workflow$data_type <- "metabolite"
    workflow
}

metab034SaveInitial <- function(manager, object) {
    manager$saveState(
        "metab_qc_complete",
        object,
        list(stage = "qc"),
        "metabolomics normalization parent"
    )
    invisible(object)
}

metab034Manifest <- function(context, manager, state_name = NULL) {
    if (is.null(state_name)) state_name <- manager$getCurrentStateName()
    row <- manager$states[[state_name]]
    store <- newArtifactStore(
        context$getPaths(),
        context$getIdentity()$project_id
    )
    artifactWorkflowStateReadManifest(store, row$manifest_relative_path)
}

metab034ScientificSlotsEqual <- function(expected, actual) {
    all(vapply(methods::slotNames(expected), function(slot_name) {
        identical(
            methods::slot(expected, slot_name),
            methods::slot(actual, slot_name)
        )
    }, logical(1)))
}

metab034ReviewedObject <- function() {
    repo <- normalizePath(testthat::test_path("..", ".."), mustWork = TRUE)
    assays <- list(
        GCMS = utils::read.delim(
            file.path(
                repo, "tests", "testdata", "e2e", "metab_combined",
                "combined_gcms_features.tsv"
            ),
            check.names = FALSE,
            stringsAsFactors = FALSE
        ),
        LCMS_Pos = utils::read.delim(
            file.path(
                repo, "tests", "testdata", "e2e", "metab_combined",
                "combined_lcms_features.tsv"
            ),
            check.names = FALSE,
            stringsAsFactors = FALSE
        )
    )
    samples <- grep("^(WT|KO)_", names(assays[[1L]]), value = TRUE)
    createMetaboliteAssayData(
        assays,
        data.frame(
            Run = samples,
            group = sub("_.*$", "", samples),
            batch = rep(c("B1", "B2"), length.out = length(samples)),
            tech_rep_group = samples,
            stringsAsFactors = FALSE
        ),
        metabolite_id_column = "Feature.Name",
        annotation_id_column = "Feature.Name",
        sample_id = "Run",
        group_id = "group",
        technical_replicate_id = "tech_rep_group",
        internal_standard_regex = "^IS_",
        args = list(normalisation_method = "quantile")
    )
}

test_that("metabolomics normalization methods persist exact materialized states", {
    metab034SkipDependencies()
    for (method in c("cyclicloess", "quantile", "scale")) {
        root <- withr::local_tempdir(pattern = paste0("metab034-", method, "-"))
        context <- metab034Context(root, paste0("omics-art-034-", method))
        manager <- metab034Manager(context)
        withr::defer(manager$close())
        workflow <- metab034Workflow(manager, context)
        before <- module_ci_metab_norm_object(
            layout = "combined",
            positive_only = TRUE,
            args = module_ci_metab_norm_args(method)
        )
        metab034SaveInitial(manager, before)
        expected <- suppressMessages(normaliseBetweenSamples(
            before,
            normalisation_method = method
        ))
        saved <- saveMetabNormState(
            workflow,
            before,
            expected,
            "between_sample_normalization",
            "metab_normalized",
            workflow$config_list,
            paste("between-sample", method),
            parameters = list(method = method),
            transformation_type = "normalization"
        )

        expect_true(metab034ScientificSlotsEqual(expected, saved), info = method)
        expect_identical(manager$getState(), expected, info = method)
        expect_identical(
            workflow$artifact_stage_results$between_sample_normalization$representation,
            "materialized",
            info = method
        )
        expect_identical(names(saved@metabolite_data), c("LCMS_Pos", "GCMS"))
        expect_identical(saved@design_matrix, before@design_matrix)
    }
})

test_that("normalization none and skip states write metadata-only generations", {
    metab034SkipDependencies()
    context <- metab034Context(withr::local_tempdir(), "omics-art-034-noops")
    manager <- metab034Manager(context)
    withr::defer(manager$close())
    workflow <- metab034Workflow(manager, context)
    before <- module_ci_metab_norm_object(positive_only = TRUE)
    metab034SaveInitial(manager, before)
    norm_data <- module_ci_metab_norm_data(before)

    state <- runMetabNormBetweenSampleNormalizationStep(
        before,
        "none",
        workflow,
        norm_data
    )
    expect_identical(state$currentS4, before)
    expect_identical(
        workflow$artifact_stage_results$between_sample_normalization$representation,
        "state_reference"
    )
    expect_identical(
        workflow$artifact_stage_results$between_sample_normalization$metrics$new_artifact_bytes,
        0
    )
    skipped <- runMetabNormRuvProgressApplyShell(
        currentS4 = state$currentS4,
        totalSteps = 6,
        ruvMode = "skip",
        autoPercentageMin = 5,
        autoPercentageMax = 10,
        ruvGroupingVariable = "tech_rep_group",
        separationMetric = "max_difference",
        kPenaltyWeight = 0.5,
        adaptiveKPenalty = TRUE,
        manualK = 2,
        manualPercentage = 10,
        experimentPaths = list(),
        groupingVariable = "group",
        shapeVariable = NULL,
        workflowData = workflow,
        normData = norm_data,
        incProgressFn = function(...) invisible(NULL)
    )
    expect_identical(skipped$currentS4, before)
    expect_identical(
        workflow$artifact_stage_results$ruv_skip$representation,
        "state_reference"
    )
    expect_identical(manager$getCurrentStateName(), "metab_normalized")
})

test_that("ITSD and log transforms preserve parameters and boundary behavior", {
    metab034SkipDependencies()
    context <- metab034Context(withr::local_tempdir(), "omics-art-034-itsd-log")
    manager <- metab034Manager(context)
    withr::defer(manager$close())
    workflow <- metab034Workflow(manager, context)
    before <- module_ci_metab_norm_object(
        layout = "combined",
        positive_only = FALSE
    )
    metab034SaveInitial(manager, before)
    norm_data <- module_ci_metab_norm_data(before)
    itsd <- suppressMessages(runMetabNormItsdNormalizationStep(
        before,
        "median",
        workflowData = workflow,
        normData = norm_data
    ))
    itsd_audit <- manager$getStateAudit("metab_itsd_norm")

    expect_identical(itsd_audit$parameters$aggregation, "median")
    expect_identical(
        itsd_audit$parameters$feature_selection$source,
        "object_internal_standard_regex"
    )
    expect_true(all(vapply(
        itsd$currentS4@metabolite_data,
        function(assay) !any(grepl("^IS_", assay$database_identifier)),
        logical(1)
    )))
    expect_null(norm_data$post_itsd_obj)
    expect_identical(
        resolveMetabNormIntermediate(workflow, norm_data, "post_itsd_obj"),
        itsd$currentS4
    )
    log_state <- suppressMessages(runMetabNormLog2TransformationStep(
        itsd$currentS4,
        1,
        workflow,
        norm_data
    ))
    log_audit <- manager$getStateAudit("metab_log2")
    expect_identical(log_audit$parameters$offset, 1)
    expect_identical(log_audit$parameters$non_finite_policy, "preserve_as_NA")
    expect_true(all(vapply(
        log_state$currentS4@metabolite_data,
        function(assay) {
            samples <- intersect(names(assay), before@design_matrix$Run)
            values <- unlist(assay[samples], use.names = FALSE)
            all(is.na(values) | is.finite(values)) && all(values >= 0, na.rm = TRUE)
        },
        logical(1)
    )))
    expect_null(norm_data$post_log2_obj)
    expect_identical(
        resolveMetabNormIntermediate(workflow, norm_data, "post_log2_obj"),
        log_state$currentS4
    )
    expect_identical(manager$getState(), log_state$currentS4)
})

test_that("RUV correction retains per-assay controls and large args by ref", {
    metab034SkipDependencies()
    context <- metab034Context(withr::local_tempdir(), "omics-art-034-ruv")
    manager <- metab034Manager(context)
    withr::defer(manager$close())
    workflow <- metab034Workflow(manager, context)
    before <- module_ci_metab_norm_object(positive_only = TRUE)
    metab034SaveInitial(manager, before)
    norm_data <- module_ci_metab_norm_data(before)
    controls <- module_ci_metab_norm_ctrl(before)
    best_k <- stats::setNames(rep(list(2L), length(controls)), names(controls))
    optimization <- list(
        ruvParams = list(mode = "automatic", percentage_min = 5, percentage_max = 20),
        ruvResults = module_ci_metab_norm_ruv_results(before, "automatic")
    )
    corrected <- runMetabNormRuvCorrectionStep(
        before,
        "tech_rep_group",
        best_k,
        controls,
        workflow,
        norm_data,
        ruvIII_C_VaryingFn = function(
            theObject,
            ruv_grouping_variable,
            ruv_number_k,
            ctrl
        ) {
            theObject@args$ruvIII_C_Varying <- list(
                grouping = ruv_grouping_variable,
                k = ruv_number_k,
                ctrl = ctrl,
                optimizer_trace = rep(seq_len(10000L), 20L)
            )
            theObject
        },
        ruvMode = "automatic",
        optimizationState = optimization
    )
    audit <- manager$getStateAudit("metab_ruv_corrected")
    manifest <- metab034Manifest(context, manager)

    expect_identical(audit$parameters$mode, "automatic")
    expect_identical(audit$parameters$optimizer$status, "recorded")
    expect_identical(
        vapply(audit$parameters$assays, `[[`, integer(1), "control_count"),
        rep(6L, length(controls))
    )
    expect_gt(audit$after_args$object_bytes, 512 * 1024)
    expect_gt(length(manifest$data$artifact_refs), length(before@metabolite_data))
    expect_identical(manager$getState(), corrected$currentS4)
    expect_identical(
        manager$getState()@args$ruvIII_C_Varying$optimizer_trace,
        corrected$currentS4@args$ruvIII_C_Varying$optimizer_trace
    )
})

test_that("correlation filter and completion preserve exact DA handoff", {
    metab034SkipDependencies()
    context <- metab034Context(withr::local_tempdir(), "omics-art-034-correlation")
    manager <- metab034Manager(context)
    withr::defer(manager$close())
    workflow <- metab034Workflow(manager, context)
    before <- module_ci_metab_norm_object(positive_only = TRUE)
    metab034SaveInitial(manager, before)
    norm_data <- module_ci_metab_norm_data(before)
    result <- dispatchMetabNormApplyCorrelation(
        workflow,
        norm_data,
        observerState = list(
            currentS4 = before,
            threshold = 0.75,
            groupingVariable = "tech_rep_group",
            notificationId = "artifact-correlation"
        ),
        showNotificationFn = function(...) invisible(NULL),
        removeNotificationFn = function(...) invisible(NULL),
        reqFn = function(value) value,
        calculateCorrelationsFn = function(...) {
            module_ci_metab_norm_correlation_pairs(
                "fail_one",
                names(before@metabolite_data)
            )
        }
    )
    current <- manager$getState()
    expect_identical(result$status, "success")
    expect_identical(current@design_matrix$Run, c("S1", "S2"))
    expect_true(all(vapply(current@metabolite_data, function(assay) {
        all(c("S1", "S2") %in% names(assay)) &&
            !any(c("S3", "S4") %in% names(assay))
    }, logical(1))))
    expect_identical(
        manager$getStateAudit("metab_correlation_filtered")$parameters$threshold,
        0.75
    )
    expect_identical(workflow$tab_status$normalization, "complete")
    expect_identical(methods::validObject(current, test = TRUE), TRUE)
})

test_that("reviewed mixed fixtures retain numerical and DA-entry parity", {
    metab034SkipDependencies()
    context <- metab034Context(withr::local_tempdir(), "omics-art-034-reviewed")
    manager <- metab034Manager(context)
    withr::defer(manager$close())
    workflow <- metab034Workflow(manager, context)
    before <- metab034ReviewedObject()
    metab034SaveInitial(manager, before)
    expected <- suppressMessages(normaliseBetweenSamples(
        before,
        normalisation_method = "quantile"
    ))
    saved <- saveMetabNormState(
        workflow,
        before,
        expected,
        "between_sample_normalization",
        "metab_normalized",
        workflow$config_list,
        "reviewed quantile normalization",
        parameters = list(
            method = "quantile",
            evidence_class = "independently_reviewed_fixture"
        ),
        transformation_type = "normalization"
    )
    completed <- saveMetabNormState(
        workflow,
        saved,
        saved,
        "correlation_skip",
        "metab_norm_complete",
        workflow$config_list,
        "reviewed normalization complete",
        parameters = list(reason = "oracle_completion"),
        status = "skipped",
        transformation_type = "no_op"
    )

    expect_identical(completed, expected)
    expect_identical(manager$getState(), expected)
    expect_identical(names(completed@metabolite_data), c("GCMS", "LCMS_Pos"))
    expect_identical(completed@design_matrix, before@design_matrix)
    expect_identical(methods::validObject(completed, test = TRUE), TRUE)
    expect_identical(
        workflow$artifact_stage_results$correlation_skip$representation,
        "state_reference"
    )
    expect_identical(
        workflow$artifact_stage_results$correlation_skip$metrics$new_artifact_bytes,
        0
    )
    manager$close()
    manager <- metab034Manager(context)
    expect_identical(manager$getState(), expected)
    expect_identical(manager$getCurrentStateName(), "metab_norm_complete")
    expect_identical(
        manager$revertToState("metab_normalized"),
        expected
    )
})

test_that("normalization reset restores the captured pre-normalization state", {
    metab034SkipDependencies()
    context <- metab034Context(withr::local_tempdir(), "omics-art-034-reset")
    manager <- metab034Manager(context)
    withr::defer(manager$close())
    workflow <- metab034Workflow(manager, context)
    before <- module_ci_metab_norm_object(positive_only = TRUE)
    metab034SaveInitial(manager, before)
    transformed <- logTransformAssays(before, offset = 1)
    saveMetabNormState(
        workflow,
        before,
        transformed,
        "log2_transformation",
        "metab_log2",
        workflow$config_list,
        "reset parent",
        parameters = list(offset = 1),
        transformation_type = "transformation"
    )
    norm_data <- module_ci_metab_norm_data(before)
    norm_data$post_filter_obj <- before
    result <- runMetabNormResetNormalizationObserverShell(
        workflow,
        norm_data,
        showNotificationFn = function(...) invisible(NULL),
        reqFn = function(value) value
    )

    expect_identical(result$status, "success")
    expect_identical(manager$getCurrentStateName(), "metab_reset")
    expect_identical(manager$getState(), before)
    expect_identical(
        manager$getStateAudit("metab_reset")$stage_id,
        "normalization_reset"
    )
    expect_false(norm_data$normalization_complete)
    expect_false(norm_data$ruv_complete)
    expect_false(norm_data$correlation_filtering_complete)
})

test_that("frozen mixed CI workload reaches supported normalization states", {
    metab034SkipDependencies()
    repo <- normalizePath(testthat::test_path("..", ".."), mustWork = TRUE)
    adapter_env <- new.env(parent = globalenv())
    sys.source(
        file.path(repo, "tools", "profiling", "omics_workload_metabolomics.R"),
        envir = adapter_env
    )
    contract <- jsonlite::read_json(
        file.path(
            repo, "tests", "testdata", "omics-parity", "metabolomics",
            "workloads", "mixed-public-ci-v1.json"
        ),
        simplifyVector = FALSE
    )
    rng_kind <- unlist(contract$rng$kind, use.names = FALSE)
    do.call(RNGkind, as.list(rng_kind))
    set.seed(as.integer(contract$rng$seed))
    root <- withr::local_tempdir()
    prepared <- adapter_env$metabWorkloadPrepare(list(
        contract = contract,
        run_dir = root
    ))
    expect_identical(
        digest::digest(
            prepared$payload_path,
            algo = "sha256",
            serialize = FALSE,
            file = TRUE
        ),
        contract$expected_digests$payload_sha256
    )
    expect_identical(
        digest::digest(
            prepared$truth_path,
            algo = "sha256",
            serialize = FALSE,
            file = TRUE
        ),
        contract$expected_digests$truth_sha256
    )
    generated <- adapter_env$metabWorkloadRead(prepared$payload_path)
    truth <- jsonlite::read_json(prepared$truth_path, simplifyVector = TRUE)
    samples <- unlist(truth$sample_ids, use.names = FALSE)
    assay_names <- unlist(truth$assays, use.names = FALSE)
    assays <- lapply(assay_names, function(assay_name) {
        value <- generated[generated$assay == assay_name, , drop = FALSE]
        value$assay <- NULL
        rownames(value) <- NULL
        value
    })
    names(assays) <- assay_names
    object <- createMetaboliteAssayData(
        assays,
        data.frame(
            Run = samples,
            group = unlist(truth$group_assignments, use.names = FALSE),
            batch = unlist(truth$batch_assignments, use.names = FALSE),
            tech_rep_group = rep(
                paste0("pair_", seq_len(length(samples) / 2L)),
                each = 2L
            ),
            stringsAsFactors = FALSE
        ),
        metabolite_id_column = "feature_id",
        annotation_id_column = "annotation",
        sample_id = "Run",
        group_id = "group",
        technical_replicate_id = "tech_rep_group",
        internal_standard_regex = "^unused_generated_regex$",
        args = list(normalisation_method = "quantile")
    )
    resolved <- resolveMetabDuplicateAssayData(
        object@metabolite_data,
        object@metabolite_id_column,
        logWarnFn = function(...) invisible(NULL)
    )
    object@metabolite_data <- resolved$resolvedAssayList
    context <- metab034Context(root, "omics-art-034-frozen-ci")
    manager <- metab034Manager(context)
    withr::defer(manager$close())
    workflow <- metab034Workflow(manager, context)
    metab034SaveInitial(manager, object)
    norm_data <- module_ci_metab_norm_data(object)
    itsd_ids <- lapply(object@metabolite_data, function(assay) {
        assay$feature_id[assay$internal_standard %in% TRUE]
    })
    current <- suppressMessages(runMetabNormItsdNormalizationStep(
        object,
        "median",
        itsdFeatureIds = itsd_ids,
        workflowData = workflow,
        normData = norm_data
    ))$currentS4
    current <- suppressMessages(runMetabNormLog2TransformationStep(
        current,
        1,
        workflow,
        norm_data
    ))$currentS4
    current <- suppressMessages(runMetabNormBetweenSampleNormalizationStep(
        current,
        "quantile",
        workflow,
        norm_data
    ))$currentS4
    controls <- lapply(current@metabolite_data, function(assay) {
        value <- seq_len(nrow(assay)) %% 10L == 0L
        names(value) <- assay$feature_id
        value
    })
    best_k <- stats::setNames(rep(list(1L), length(controls)), names(controls))
    current <- runMetabNormRuvCorrectionStep(
        current,
        "tech_rep_group",
        best_k,
        controls,
        workflow,
        norm_data,
        ruvIII_C_VaryingFn = function(
            theObject,
            ruv_grouping_variable,
            ruv_number_k,
            ctrl
        ) {
            theObject@args$ruvIII_C_Varying <- list(
                grouping = ruv_grouping_variable,
                k = ruv_number_k,
                ctrl = ctrl,
                evidence_class = "generated_scaling"
            )
            theObject
        },
        ruvMode = "manual"
    )$currentS4
    completed <- saveMetabNormState(
        workflow,
        current,
        current,
        "correlation_skip",
        "metab_norm_complete",
        workflow$config_list,
        "generated normalization complete",
        parameters = list(
            evidence_class = "generated_scaling",
            claim_boundary = paste(
                "contract invariants and resource behavior only;",
                "not parser or biological validation"
            )
        ),
        status = "skipped",
        transformation_type = "no_op"
    )

    expect_identical(manager$getState(), completed)
    expect_identical(names(completed@metabolite_data), assay_names)
    expect_identical(completed@design_matrix$Run, samples)
    expect_identical(methods::validObject(completed, test = TRUE), TRUE)
    expect_true(all(c(
        "itsd_normalization", "log2_transformation",
        "between_sample_normalization", "ruv_correction", "correlation_skip"
    ) %in% names(workflow$artifact_stage_results)))
    expect_identical(
        workflow$artifact_stage_results$correlation_skip$metrics$new_artifact_bytes,
        0
    )
})

test_that("normalization artifact mode disables independently", {
    before <- module_ci_metab_norm_object(positive_only = TRUE)
    manager <- WorkflowState$new()
    manager$setWorkflowType("metabolomics_standard")
    metab034SaveInitial(manager, before)
    workflow <- new.env(parent = emptyenv())
    workflow$state_manager <- manager
    workflow$config_list <- list()
    withr::local_options(
        multischolar.metabolomics.norm_artifacts = "disabled"
    )
    after <- logTransformAssays(before, offset = 1)
    saved <- saveMetabNormState(
        workflow,
        before,
        after,
        "log2_transformation",
        "metab_log2",
        workflow$config_list,
        "memory-only log2",
        parameters = list(offset = 1),
        transformation_type = "transformation"
    )

    expect_identical(manager$getState(), saved)
    expect_null(workflow$artifact_stage_results)
})
