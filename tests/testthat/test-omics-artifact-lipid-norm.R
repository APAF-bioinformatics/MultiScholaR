lipid043SkipDependencies <- function() {
    for (package in c("DBI", "duckdb", "filelock", "arrow", "limma")) {
        testthat::skip_if_not_installed(package)
    }
}

lipid043Context <- function(root, project_id) {
    paths <- list(
        base_dir = root,
        project_id = project_id,
        omic_type = "lipidomics",
        omic_label = "lipidomics-study",
        source_dir = file.path(root, "source"),
        results_dir = file.path(root, "results")
    )
    dir.create(paths$source_dir, recursive = TRUE, showWarnings = FALSE)
    dir.create(paths$results_dir, recursive = TRUE, showWarnings = FALSE)
    capability <- workflowCapabilityCatalogue()[[
        "lipidomics.lipidsearch.lipid.standard.v1"
    ]]
    capability$artifact_eligible <- TRUE
    capability$auto_eligible <- FALSE
    capability$maximum_artifact_rollout <- "dual_write"
    capability$explicit_maximum_artifact_rollout <- "dual_write"
    context <- createWorkflowContext(
        paths,
        "lipidomics",
        "lipidomics-study",
        storage_policy = list(
            requested_backend = "artifact",
            requested_rollout = "dual_write",
            migration_requested = TRUE,
            project_id = project_id
        )
    )
    bindWorkflowContextFromImport(
        context,
        workflow_type = "lipidomics_standard",
        input_format = "lipidsearch",
        data_level = "lipid",
        capabilities = list(capability)
    )
    context
}

lipid043Manager <- function(context) {
    manager <- ArtifactWorkflowState$new(
        workflow_context = context,
        dehydrate_fn = dehydrateLipidomicsS4Artifact,
        validate_bundle_fn = validateLipidomicsS4Bundle,
        hydrate_fn = hydrateLipidomicsS4Artifact,
        descriptor_contract = artifactStageDescriptorContract(
            artifactLipidomicsWorkflowDescriptor()
        )
    )
    manager$setWorkflowType("lipidomics_standard")
    manager
}

lipid043Workflow <- function(manager, context) {
    workflow <- new.env(parent = emptyenv())
    workflow$state_manager <- manager
    workflow$workflow_context <- context
    workflow$config_list <- list(
        globalParameters = list(workflow_type = "lipidomics_standard"),
        normalization = list(mode = "artifact-ticket")
    )
    workflow$artifact_stage_results <- list()
    workflow$processing_log <- list()
    workflow$tab_status <- list(
        quality_control = "complete",
        normalization = "pending",
        differential_analysis = "locked"
    )
    workflow$data_format <- "lipidsearch"
    workflow$data_type <- "lipid"
    workflow
}

lipid043SaveInitial <- function(manager, object) {
    manager$saveState(
        "lipid_qc_complete",
        object,
        list(stage = "qc"),
        "lipidomics normalization parent"
    )
    invisible(object)
}

lipid043Manifest <- function(context, manager, state_name = NULL) {
    if (is.null(state_name)) state_name <- manager$getCurrentStateName()
    row <- manager$states[[state_name]]
    store <- newArtifactStore(
        context$getPaths(),
        context$getIdentity()$project_id
    )
    artifactWorkflowStateReadManifest(store, row$manifest_relative_path)
}

lipid043ScientificSlotsEqual <- function(expected, actual) {
    all(vapply(methods::slotNames(expected), function(slot_name) {
        identical(
            methods::slot(expected, slot_name),
            methods::slot(actual, slot_name)
        )
    }, logical(1)))
}

lipid043ReviewedObject <- function() {
    object <- lipid042ReviewedObject()
    object@args <- list(normalisation_method = "quantile")
    object
}

test_that("lipidomics normalization methods persist exact materialized states", {
    lipid043SkipDependencies()
    for (method in c("cyclicloess", "quantile", "scale")) {
        root <- withr::local_tempdir(pattern = paste0("lipid043-", method, "-"))
        context <- lipid043Context(root, paste0("omics-art-043-", method))
        manager <- lipid043Manager(context)
        withr::defer(manager$close())
        workflow <- lipid043Workflow(manager, context)
        before <- module_ci_lipid_norm_object(
            layout = "combined",
            positive_only = TRUE,
            args = module_ci_lipid_norm_args(method)
        )
        lipid043SaveInitial(manager, before)
        expected <- suppressMessages(normaliseBetweenSamples(
            before,
            normalisation_method = method
        ))
        saved <- saveLipidNormState(
            workflow,
            before,
            expected,
            "between_sample_normalization",
            "lipid_normalized",
            workflow$config_list,
            paste("between-sample", method),
            parameters = list(method = method),
            transformation_type = "normalization"
        )

        expect_true(lipid043ScientificSlotsEqual(expected, saved), info = method)
        expect_identical(manager$getState(), expected, info = method)
        expect_identical(
            workflow$artifact_stage_results$between_sample_normalization$representation,
            "materialized",
            info = method
        )
        expect_identical(names(saved@lipid_data), c("LCMS_Pos", "GCMS"))
        expect_identical(saved@design_matrix, before@design_matrix)
    }
})

test_that("normalization none and skip states write metadata-only generations", {
    lipid043SkipDependencies()
    context <- lipid043Context(withr::local_tempdir(), "omics-art-043-noops")
    manager <- lipid043Manager(context)
    withr::defer(manager$close())
    workflow <- lipid043Workflow(manager, context)
    before <- module_ci_lipid_norm_object(positive_only = TRUE)
    lipid043SaveInitial(manager, before)
    norm_data <- module_ci_lipid_norm_data(before)

    current <- persistLipidNormBetweenSampleState(
        workflow,
        before,
        before,
        "none",
        "between-sample none"
    )
    expect_identical(current, before)
    expect_identical(
        workflow$artifact_stage_results$between_sample_normalization$representation,
        "state_reference"
    )
    expect_identical(
        workflow$artifact_stage_results$between_sample_normalization$metrics$new_artifact_bytes,
        0
    )
    skipped <- recordLipidNormNoOp(
        workflow,
        current,
        "ruv_skip",
        list(
            mode = "skip",
            grouping_variable = "tech_rep_group",
            reason = "user_selected_skip"
        )
    )
    expect_identical(skipped, before)
    expect_identical(
        workflow$artifact_stage_results$ruv_skip$representation,
        "state_reference"
    )
    expect_identical(manager$getCurrentStateName(), "lipid_normalized")
})

test_that("ITSD and log transforms preserve parameters and boundary behavior", {
    lipid043SkipDependencies()
    context <- lipid043Context(withr::local_tempdir(), "omics-art-043-itsd-log")
    manager <- lipid043Manager(context)
    withr::defer(manager$close())
    workflow <- lipid043Workflow(manager, context)
    before <- module_ci_lipid_norm_object(
        layout = "combined",
        positive_only = FALSE
    )
    lipid043SaveInitial(manager, before)
    norm_data <- module_ci_lipid_norm_data(before)
    normalized_itsd <- suppressMessages(normaliseUntransformedData(
        before,
        method = "ITSD",
        itsd_aggregation = "median"
    ))
    current <- persistLipidNormItsdState(
        workflow,
        norm_data,
        before,
        normalized_itsd,
        "median",
        NULL,
        "ITSD median"
    )
    itsd_audit <- manager$getStateAudit("lipid_itsd_norm")

    expect_identical(itsd_audit$parameters$aggregation, "median")
    expect_identical(
        itsd_audit$parameters$feature_selection$source,
        "object_internal_standard_regex"
    )
    expect_true(all(vapply(
        current@lipid_data,
        function(assay) !any(grepl(
            "^IS_",
            assay[[current@lipid_id_column]]
        )),
        logical(1)
    )))
    expect_null(norm_data$post_itsd_obj)
    expect_identical(
        resolveLipidNormIntermediate(workflow, norm_data, "post_itsd_obj"),
        current
    )
    transformed <- suppressMessages(logTransformAssays(current, offset = 1))
    log_state <- persistLipidNormLog2State(
        workflow,
        norm_data,
        current,
        transformed,
        1,
        "log2 offset 1"
    )
    log_audit <- manager$getStateAudit("lipid_log2")
    expect_identical(log_audit$parameters$offset, 1)
    expect_identical(log_audit$parameters$non_finite_policy, "preserve_as_NA")
    expect_true(all(vapply(
        log_state@lipid_data,
        function(assay) {
            samples <- intersect(names(assay), before@design_matrix$Run)
            values <- unlist(assay[samples], use.names = FALSE)
            all(is.na(values) | is.finite(values)) && all(values >= 0, na.rm = TRUE)
        },
        logical(1)
    )))
    expect_null(norm_data$post_log2_obj)
    expect_identical(
        resolveLipidNormIntermediate(workflow, norm_data, "post_log2_obj"),
        log_state
    )
    expect_identical(manager$getState(), log_state)
})

test_that("RUV correction retains per-assay controls and large args by ref", {
    lipid043SkipDependencies()
    context <- lipid043Context(withr::local_tempdir(), "omics-art-043-ruv")
    manager <- lipid043Manager(context)
    withr::defer(manager$close())
    workflow <- lipid043Workflow(manager, context)
    before <- module_ci_lipid_norm_object(positive_only = TRUE)
    lipid043SaveInitial(manager, before)
    norm_data <- module_ci_lipid_norm_data(before)
    controls <- module_ci_lipid_norm_ctrl(before)
    best_k <- stats::setNames(rep(list(2L), length(controls)), names(controls))
    optimization <- list(
        ruvParams = list(mode = "automatic", percentage_min = 5, percentage_max = 20),
        ruvResults = module_ci_lipid_norm_ruv_results(before, "automatic")
    )
    after <- before
    after@args$ruvIII_C_Varying <- list(
        grouping = "tech_rep_group",
        k = best_k,
        ctrl = controls,
        optimizer_trace = rep(seq_len(10000L), 20L)
    )
    corrected <- persistLipidNormRuvState(
        workflow,
        before,
        after,
        "tech_rep_group",
        best_k,
        controls,
        "automatic",
        optimization,
        "RUV correction"
    )
    audit <- manager$getStateAudit("lipid_ruv_corrected")
    manifest <- lipid043Manifest(context, manager)

    expect_identical(audit$parameters$mode, "automatic")
    expect_identical(audit$parameters$optimizer$status, "recorded")
    expect_identical(
        vapply(audit$parameters$assays, `[[`, integer(1), "control_count"),
        rep(6L, length(controls))
    )
    expect_gt(audit$after_args$object_bytes, 512 * 1024)
    expect_gt(length(manifest$data$artifact_refs), length(before@lipid_data))
    expect_identical(manager$getState(), corrected)
    expect_identical(
        manager$getState()@args$ruvIII_C_Varying$optimizer_trace,
        corrected@args$ruvIII_C_Varying$optimizer_trace
    )
})

test_that("correlation filter and completion preserve exact DA handoff", {
    lipid043SkipDependencies()
    context <- lipid043Context(withr::local_tempdir(), "omics-art-043-correlation")
    manager <- lipid043Manager(context)
    withr::defer(manager$close())
    workflow <- lipid043Workflow(manager, context)
    before <- module_ci_lipid_norm_object(positive_only = TRUE)
    lipid043SaveInitial(manager, before)
    norm_data <- module_ci_lipid_norm_data(before)
    norm_data$post_norm_obj <- before
    norm_data$normalization_complete <- TRUE
    result <- handleLipidNormApplyCorrelationFilter(
        input = list(
            min_pearson_correlation_threshold = 0.75,
            ruv_grouping_variable = "tech_rep_group"
        ),
        workflowData = workflow,
        normData = norm_data,
        addLog = function(...) invisible(NULL),
        showNotificationFn = function(...) invisible(NULL),
        removeNotificationFn = function(...) invisible(NULL),
        reqFn = function(value) value,
        pearsonCorForSamplePairsFn = function(...) {
            module_ci_lipid_norm_correlation_pairs(
                "fail_one",
                names(before@lipid_data)
            )
        },
        filterSamplesByLipidCorrelationThresholdFn = function(
            theObject,
            pearson_correlation_per_pair,
            min_pearson_correlation_threshold
        ) {
            filtered <- theObject
            filtered@design_matrix <- filtered@design_matrix[
                filtered@design_matrix$Run %in% c("S1", "S2"),
                ,
                drop = FALSE
            ]
            filtered@lipid_data <- lapply(filtered@lipid_data, function(assay) {
                assay[setdiff(names(assay), c("S3", "S4"))]
            })
            methods::validObject(filtered)
            filtered
        },
        logInfoFn = function(...) invisible(NULL),
        logErrorFn = function(message) stop(message)
    )
    current <- manager$getState()
    expect_true(result)
    expect_identical(current@design_matrix$Run, c("S1", "S2"))
    expect_true(all(vapply(current@lipid_data, function(assay) {
        all(c("S1", "S2") %in% names(assay)) &&
            !any(c("S3", "S4") %in% names(assay))
    }, logical(1))))
    expect_identical(
        manager$getStateAudit("lipid_correlation_filtered")$parameters$threshold,
        0.75
    )
    expect_identical(workflow$tab_status$normalization, "complete")
    expect_identical(methods::validObject(current, test = TRUE), TRUE)
})

test_that("reviewed mixed fixtures retain numerical and DA-entry parity", {
    lipid043SkipDependencies()
    context <- lipid043Context(withr::local_tempdir(), "omics-art-043-reviewed")
    manager <- lipid043Manager(context)
    withr::defer(manager$close())
    workflow <- lipid043Workflow(manager, context)
    before <- lipid043ReviewedObject()
    lipid043SaveInitial(manager, before)
    expected <- suppressMessages(normaliseBetweenSamples(
        before,
        normalisation_method = "quantile"
    ))
    saved <- saveLipidNormState(
        workflow,
        before,
        expected,
        "between_sample_normalization",
        "lipid_normalized",
        workflow$config_list,
        "reviewed quantile normalization",
        parameters = list(
            method = "quantile",
            evidence_class = "independently_reviewed_fixture"
        ),
        transformation_type = "normalization"
    )
    completed <- saveLipidNormState(
        workflow,
        saved,
        saved,
        "correlation_skip",
        "lipid_norm_complete",
        workflow$config_list,
        "reviewed normalization complete",
        parameters = list(reason = "oracle_completion"),
        status = "skipped",
        transformation_type = "no_op"
    )

    expect_identical(completed, expected)
    expect_identical(manager$getState(), expected)
    expect_identical(
        names(completed@lipid_data),
        c("LCMS_Pos", "LCMS_Neg", "GCMS")
    )
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
    manager <- lipid043Manager(context)
    expect_identical(manager$getState(), expected)
    expect_identical(manager$getCurrentStateName(), "lipid_norm_complete")
    expect_identical(
        manager$revertToState("lipid_normalized"),
        expected
    )
})

test_that("normalization reset restores the captured pre-normalization state", {
    lipid043SkipDependencies()
    context <- lipid043Context(withr::local_tempdir(), "omics-art-043-reset")
    manager <- lipid043Manager(context)
    withr::defer(manager$close())
    workflow <- lipid043Workflow(manager, context)
    before <- module_ci_lipid_norm_object(positive_only = TRUE)
    lipid043SaveInitial(manager, before)
    transformed <- logTransformAssays(before, offset = 1)
    saveLipidNormState(
        workflow,
        before,
        transformed,
        "log2_transformation",
        "lipid_log2",
        workflow$config_list,
        "reset parent",
        parameters = list(offset = 1),
        transformation_type = "transformation"
    )
    norm_data <- module_ci_lipid_norm_data(before)
    norm_data$post_filter_obj <- before
    result <- handleLipidNormResetNormalization(
        workflowData = workflow,
        normData = norm_data,
        addLog = function(...) invisible(NULL),
        showNotificationFn = function(...) invisible(NULL),
        reqFn = function(value) value
    )

    expect_true(result)
    expect_identical(manager$getCurrentStateName(), "lipid_reset")
    expect_identical(manager$getState(), before)
    expect_identical(
        manager$getStateAudit("lipid_reset")$stage_id,
        "normalization_reset"
    )
    expect_false(norm_data$normalization_complete)
    expect_false(norm_data$ruv_complete)
    expect_false(norm_data$correlation_filtering_complete)
})

test_that("frozen mixed CI workload reaches supported normalization states", {
    lipid043SkipDependencies()
    repo <- normalizePath(testthat::test_path("..", ".."), mustWork = TRUE)
    adapter_env <- new.env(parent = globalenv())
    sys.source(
        file.path(repo, "tools", "profiling", "omics_workload_lipidomics.R"),
        envir = adapter_env
    )
    contract <- jsonlite::read_json(
        file.path(
            repo, "tests", "testdata", "omics-parity", "lipidomics",
            "workloads", "mixed-public-ci-v1.json"
        ),
        simplifyVector = FALSE
    )
    rng_kind <- unlist(contract$rng$kind, use.names = FALSE)
    do.call(RNGkind, as.list(rng_kind))
    set.seed(as.integer(contract$rng$seed))
    root <- withr::local_tempdir()
    prepared <- adapter_env$lipidWorkloadPrepare(list(
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
    generated <- utils::read.delim(
        prepared$payload_path,
        check.names = FALSE,
        stringsAsFactors = FALSE
    )
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
    object <- createLipidomicsAssayData(
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
        lipid_id_column = "lipid_id",
        annotation_id_column = "lipid_class",
        sample_id = "Run",
        group_id = "group",
        technical_replicate_id = "tech_rep_group",
        internal_standard_regex = "^unused_generated_regex$",
        args = list(normalisation_method = "quantile")
    )
    resolved <- lipid042ResolveDuplicates(object)
    object@lipid_data <- resolved$resolvedAssayList
    context <- lipid043Context(root, "omics-art-043-frozen-ci")
    manager <- lipid043Manager(context)
    withr::defer(manager$close())
    workflow <- lipid043Workflow(manager, context)
    lipid043SaveInitial(manager, object)
    norm_data <- module_ci_lipid_norm_data(object)
    itsd_ids <- lapply(object@lipid_data, function(assay) {
        assay$lipid_id[assay$internal_standard %in% TRUE]
    })
    normalized_itsd <- suppressMessages(normaliseUntransformedData(
        object,
        method = "ITSD",
        itsd_aggregation = "median",
        itsd_feature_ids = itsd_ids
    ))
    current <- persistLipidNormItsdState(
        workflow,
        norm_data,
        object,
        normalized_itsd,
        "median",
        itsd_ids,
        "generated ITSD normalization"
    )
    transformed <- suppressMessages(logTransformAssays(current, offset = 1))
    current <- persistLipidNormLog2State(
        workflow,
        norm_data,
        current,
        transformed,
        1,
        "generated log2"
    )
    normalized <- suppressMessages(normaliseBetweenSamples(
        current,
        normalisation_method = "quantile"
    ))
    current <- persistLipidNormBetweenSampleState(
        workflow,
        current,
        normalized,
        "quantile",
        "generated quantile normalization"
    )
    controls <- lapply(current@lipid_data, function(assay) {
        value <- seq_len(nrow(assay)) %% 10L == 0L
        names(value) <- assay$lipid_id
        value
    })
    best_k <- stats::setNames(rep(list(1L), length(controls)), names(controls))
    corrected <- current
    corrected@args$ruvIII_C_Varying <- list(
        grouping = "tech_rep_group",
        k = best_k,
        ctrl = controls,
        evidence_class = "generated_scaling"
    )
    current <- persistLipidNormRuvState(
        workflow,
        current,
        corrected,
        "tech_rep_group",
        best_k,
        controls,
        "manual",
        NULL,
        "generated RUV correction"
    )
    completed <- saveLipidNormState(
        workflow,
        current,
        current,
        "correlation_skip",
        "lipid_norm_complete",
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
    expect_identical(names(completed@lipid_data), assay_names)
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

test_that("production normalization pipeline keeps redundant S4 states by ref", {
    lipid043SkipDependencies()
    context <- lipid043Context(
        withr::local_tempdir(),
        "omics-art-043-production-pipeline"
    )
    manager <- lipid043Manager(context)
    withr::defer(manager$close())
    workflow <- lipid043Workflow(manager, context)
    before <- module_ci_lipid_norm_object(positive_only = TRUE)
    lipid043SaveInitial(manager, before)
    norm_data <- module_ci_lipid_norm_data(before)
    result <- handleLipidNormRunNormalization(
        input = module_ci_lipid_norm_input(
            norm_method = "none",
            ruv_mode = "skip",
            apply_itsd = FALSE,
            log_offset = 1
        ),
        workflowData = workflow,
        experimentPaths = list(lipid_qc_dir = withr::local_tempdir()),
        omicType = "lipidomics",
        normData = norm_data,
        addLog = function(...) invisible(NULL),
        getPlotAestheticsFn = function() {
            list(color_var = "group", shape_var = NULL)
        },
        reqFn = function(value) {
            if (is.null(value)) stop("required value missing")
            invisible(value)
        },
        withProgressFn = function(message, value, expr) force(expr),
        incProgressFn = function(...) invisible(NULL),
        generateLipidQcPlotsFn = function(...) invisible(NULL),
        generateCompositeFromFilesFn = function(...) NULL,
        savePlotFn = function(...) invisible(NULL),
        dirExistsFn = function(...) FALSE,
        showNotificationFn = function(...) invisible(NULL),
        logWarnFn = function(...) invisible(NULL),
        logErrorFn = function(message) stop(message)
    )

    expect_true(result)
    expect_null(norm_data$post_filter_obj)
    expect_null(norm_data$post_itsd_obj)
    expect_null(norm_data$post_log2_obj)
    expect_s4_class(norm_data$post_norm_obj, "LipidomicsAssayData")
    expect_s4_class(norm_data$ruv_corrected_obj, "LipidomicsAssayData")
    expect_identical(
        names(norm_data$artifact_state_refs),
        c("post_filter_obj", "post_log2_obj")
    )
    expect_identical(
        workflow$artifact_stage_results$itsd_skip$representation,
        "state_reference"
    )
    expect_identical(
        workflow$artifact_stage_results$between_sample_normalization$representation,
        "state_reference"
    )
    expect_identical(
        workflow$artifact_stage_results$ruv_skip$representation,
        "state_reference"
    )
    expect_identical(manager$getCurrentStateName(), "lipid_normalized")
    expect_identical(methods::validObject(manager$getState(), test = TRUE), TRUE)
})

test_that("normalization artifact mode disables independently", {
    before <- module_ci_lipid_norm_object(positive_only = TRUE)
    manager <- WorkflowState$new()
    manager$setWorkflowType("lipidomics_standard")
    lipid043SaveInitial(manager, before)
    workflow <- new.env(parent = emptyenv())
    workflow$state_manager <- manager
    workflow$config_list <- list()
    withr::local_options(
        multischolar.lipidomics.norm_artifacts = "disabled"
    )
    after <- logTransformAssays(before, offset = 1)
    saved <- saveLipidNormState(
        workflow,
        before,
        after,
        "log2_transformation",
        "lipid_log2",
        workflow$config_list,
        "memory-only log2",
        parameters = list(offset = 1),
        transformation_type = "transformation"
    )

    expect_identical(manager$getState(), saved)
    expect_null(workflow$artifact_stage_results)
})
