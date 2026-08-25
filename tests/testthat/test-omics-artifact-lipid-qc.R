test_that("lipidomics intensity recipes bind each exact parent assay", {
    lipid042SkipDependencies()
    context <- lipid042Context(withr::local_tempdir())
    manager <- lipid042Manager(context)
    withr::defer(manager$close())
    workflow <- lipid042Workflow(manager, context)
    parent <- lipid042Object(include_duplicates = FALSE)
    lipid042SaveInitial(manager, parent)
    parent_manifest <- lipid042Manifest(context, manager)
    child <- lipid042Filter(parent)

    saved <- saveLipidQcState(
        workflow,
        parent,
        child,
        stage_id = "intensity_filter",
        state_name = "lipid_intensity_filtered",
        config_object = workflow$config_list,
        description = "ticket intensity filter",
        parameters = list(
            intensity_cutoff_percentile = 10,
            proportion_below_cutoff = 0.5
        ),
        transformation_type = "filter",
        now = as.POSIXct("2026-08-25", tz = "UTC")
    )
    hydrated <- manager$getState()
    manifest <- lipid042Manifest(context, manager)
    metadata <- lipid042Metadata(context, manager)
    derivation <- metadata$derivation

    expect_identical(hydrated, saved)
    expect_identical(methods::validObject(hydrated, test = TRUE), TRUE)
    expect_identical(derivation$representation, "assay_selection")
    expect_identical(
        vapply(derivation$assays, `[[`, character(1), "assay_name"),
        names(parent@lipid_data)
    )
    expect_identical(
        derivation$parent_generation_id,
        parent_manifest$generation_id
    )
    expect_identical(
        derivation$parent_semantic_digest,
        parent_manifest$data$semantic_digest
    )
    for (entry in derivation$assays) {
        parent_ref <- parent_manifest$data$artifact_refs[[
            entry$parent_payload_key
        ]]
        expect_identical(entry$parent_artifact_id, parent_ref$artifact_id)
        recipe <- lipid042ReadRef(
            context,
            manifest$data$artifact_refs[[entry$recipe_ref_name]]
        )
        expect_identical(sum(recipe$selection_status == "selected"), 3L)
        expect_true(all(nzchar(recipe$failure_reason[
            recipe$selection_status == "rejected"
        ])))
        expect_true(all(entry$key_columns %in% names(recipe)))
    }
    expect_identical(
        workflow$artifact_stage_results$intensity_filter$representation,
        "assay_selection"
    )
    expect_true(workflow$artifact_stage_results$intensity_filter$source_payloads_retained)
    expect_false(workflow$artifact_stage_results$intensity_filter$eviction_performed)
})

test_that("lipidomics QC materializes duplicate resolution and records policy", {
    lipid042SkipDependencies()
    context <- lipid042Context(
        withr::local_tempdir(),
        "omics-art-042-duplicates"
    )
    manager <- lipid042Manager(context)
    withr::defer(manager$close())
    workflow <- lipid042Workflow(manager, context)
    parent <- lipid042Object(include_duplicates = TRUE, include_itsd = FALSE)
    lipid042SaveInitial(manager, parent)
    resolution <- lipid042ResolveDuplicates(parent)
    child <- parent
    child@lipid_data <- resolution$resolvedAssayList

    saved <- saveLipidQcState(
        workflow,
        parent,
        child,
        stage_id = "duplicate_resolution",
        state_name = "lipid_duplicates_resolved",
        config_object = workflow$config_list,
        description = "ticket duplicate resolution",
        parameters = lipidQcDuplicateProvenance(
            parent,
            child,
            resolution$statsList
        ),
        transformation_type = "materialization",
        now = as.POSIXct("2026-08-25", tz = "UTC")
    )
    record <- manager$getStateAudit()

    expect_identical(manager$getState(), saved)
    expect_identical(
        workflow$artifact_stage_results$duplicate_resolution$representation,
        "materialized"
    )
    expect_identical(record$parameters$method, "highest_mean_intensity")
    expect_identical(record$parameters$aggregation, "none")
    expect_identical(record$parameters$tie_break, "first_stable_input_row")
    expect_true(all(vapply(
        record$parameters$assays,
        \(assay) length(assay$value_columns) > 0L &&
            nzchar(assay$value_digest) &&
            nzchar(assay$annotation_digest) &&
            assay$removed_rows > 0L,
        logical(1)
    )))
    expect_true(all(vapply(
        saved@lipid_data,
        function(assay) anyDuplicated(assay[[saved@lipid_id_column]]) == 0L,
        logical(1)
    )))
})

test_that("lipidomics assay recipes chain through recipe-backed parents", {
    lipid042SkipDependencies()
    context <- lipid042Context(withr::local_tempdir(), "omics-art-042-chain")
    manager <- lipid042Manager(context)
    withr::defer(manager$close())
    workflow <- lipid042Workflow(manager, context)
    parent <- lipid042Object(include_duplicates = FALSE)
    lipid042SaveInitial(manager, parent)
    first <- lipid042Filter(parent)
    first <- saveLipidQcState(
        workflow,
        parent,
        first,
        "intensity_filter",
        "lipid_intensity_filtered",
        workflow$config_list,
        "first filter",
        parameters = list(percentile = 10, proportion = 0.5),
        transformation_type = "filter"
    )
    first_manifest <- lipid042Manifest(context, manager)
    first_metadata <- lipid042Metadata(context, manager)
    second <- first
    second@lipid_data <- lapply(first@lipid_data, function(assay) {
        value <- assay[c(2L, 1L), , drop = FALSE]
        rownames(value) <- NULL
        value
    })
    second <- saveLipidQcState(
        workflow,
        first,
        second,
        "intensity_filter",
        "lipid_intensity_refiltered",
        workflow$config_list,
        "second filter",
        parameters = list(percentile = 20, proportion = 0.5),
        transformation_type = "filter"
    )
    second_metadata <- lipid042Metadata(context, manager)

    expect_identical(manager$getState(), second)
    expect_identical(second_metadata$derivation$representation, "assay_selection")
    for (index in seq_along(second_metadata$derivation$assays)) {
        entry <- second_metadata$derivation$assays[[index]]
        parent_entry <- first_metadata$derivation$assays[[index]]
        expect_identical(entry$parent_ref_name, parent_entry$recipe_ref_name)
        expect_identical(
            entry$parent_artifact_id,
            first_manifest$data$artifact_refs[[parent_entry$recipe_ref_name]]$artifact_id
        )
    }
})

test_that("lipidomics internal-standard state retains assay provenance", {
    lipid042SkipDependencies()
    context <- lipid042Context(withr::local_tempdir(), "omics-art-042-itsd")
    manager <- lipid042Manager(context)
    withr::defer(manager$close())
    workflow <- lipid042Workflow(manager, context)
    parent <- lipid042Object(include_duplicates = FALSE, include_itsd = TRUE)
    lipid042SaveInitial(manager, parent)
    analysis <- lipid042ItsdAnalysis(parent, "^IS_")

    result <- recordLipidQcItsdAnalysis(
        workflow,
        parent,
        analysis,
        "^IS_"
    )
    record <- manager$getStateAudit()
    provenance <- result$provenance

    expect_true(result$enabled)
    expect_identical(provenance$requested_regex, "^IS_")
    expect_identical(provenance$regex_source, "user_input")
    expect_identical(provenance$assay_scope, names(parent@lipid_data))
    expect_identical(provenance$total_matches, analysis$nIsTotal)
    expect_true(all(vapply(
        provenance$assays,
        function(assay) assay$match_count > 0L &&
            length(assay$matched_ids) == assay$match_count,
        logical(1)
    )))
    expect_identical(record$parameters$resolved_regex, "^IS_")
    expect_identical(
        workflow$artifact_stage_results$internal_standard_analysis$representation,
        "state_reference"
    )
    manifest <- lipid042Manifest(context, manager)
    metadata <- lipid042Metadata(context, manager)
    expect_length(manifest$data$artifact_refs, 0L)
    expect_identical(
        workflow$artifact_stage_results$internal_standard_analysis$metrics$new_artifact_bytes,
        0
    )
    tampered <- metadata
    tampered$parent_manifest_digest <- paste(rep("0", 64L), collapse = "")
    store <- newArtifactStore(
        context$getPaths(),
        context$getIdentity()$project_id
    )
    expect_error(
        artifactWorkflowStateHydrateStateReference(
            store,
            manifest,
            tampered,
            hydrateLipidomicsS4Artifact,
            character()
        ),
        class = "multischolar_artifact_state_reference_parent_mismatch"
    )
})

test_that("lipidomics artifact hydration preserves QC plot data", {
    lipid042SkipDependencies()
    context <- lipid042Context(withr::local_tempdir(), "omics-art-042-plots")
    manager <- lipid042Manager(context)
    withr::defer(manager$close())
    workflow <- lipid042Workflow(manager, context)
    parent <- lipid042Object(
        include_duplicates = FALSE,
        include_itsd = FALSE,
        layout = "gc"
    )
    lipid042SaveInitial(manager, parent)
    child <- lipid042Filter(parent)
    saved <- saveLipidQcState(
        workflow,
        parent,
        child,
        "intensity_filter",
        "lipid_intensity_filtered",
        workflow$config_list,
        "plot parity",
        parameters = list(percentile = 10, proportion = 0.5),
        transformation_type = "filter"
    )
    hydrated <- manager$getState()

    expect_identical(
        lipid042PlotBuildData(saved, plotPca, grouping_variable = "group"),
        lipid042PlotBuildData(hydrated, plotPca, grouping_variable = "group")
    )
    expect_identical(
        lipid042PlotBuildData(saved, plotRle, grouping_variable = "group"),
        lipid042PlotBuildData(hydrated, plotRle, grouping_variable = "group")
    )
    expect_identical(
        lipid042PlotBuildData(saved, plotDensity, grouping_variable = "group"),
        lipid042PlotBuildData(hydrated, plotDensity, grouping_variable = "group")
    )
    expect_identical(
        lipid042PlotBuildData(
            saved,
            plotPearson,
            tech_rep_remove_regex = "",
            correlation_group = "tech_rep_group"
        ),
        lipid042PlotBuildData(
            hydrated,
            plotPearson,
            tech_rep_remove_regex = "",
            correlation_group = "tech_rep_group"
        )
    )
})

test_that("reviewed mixed fixtures retain scientific state and explicit IS skip", {
    lipid042SkipDependencies()
    context <- lipid042Context(withr::local_tempdir(), "omics-art-042-reviewed")
    manager <- lipid042Manager(context)
    withr::defer(manager$close())
    workflow <- lipid042Workflow(manager, context)
    parent <- lipid042ReviewedObject()
    lipid042SaveInitial(manager, parent)
    expected <- lipidIntensityFiltering(
        parent,
        lipids_intensity_cutoff_percentile = 10,
        lipids_proportion_of_samples_below_cutoff = 0.5
    )
    current <- saveLipidQcState(
        workflow,
        parent,
        expected,
        "intensity_filter",
        "lipid_intensity_filtered",
        workflow$config_list,
        "reviewed fixture intensity",
        parameters = list(percentile = 10, proportion = 0.5),
        transformation_type = "filter"
    )
    error <- tryCatch(
        {
            lipid042ItsdAnalysis(current, "^IS_", require_matches = TRUE)
            NULL
        },
        error = identity
    )
    expect_s3_class(error, "error")
    skipped <- recordLipidQcItsdFailure(
        workflow,
        current,
        "^IS_",
        error
    )
    current <- manager$getState()
    finalized <- saveLipidQcState(
        workflow,
        current,
        current,
        "qc_finalization",
        "lipid_qc_complete",
        workflow$config_list,
        "reviewed fixture QC complete",
        parameters = list(evidence_class = "independently_reviewed_fixture"),
        status = "complete",
        transformation_type = "finalization"
    )

    scientific_slots <- setdiff(methods::slotNames(expected), "args")
    expect_true(all(vapply(scientific_slots, function(slot_name) {
        identical(
            methods::slot(finalized, slot_name),
            methods::slot(expected, slot_name)
        )
    }, logical(1))))
    expect_identical(
        names(finalized@lipid_data),
        c("LCMS_Pos", "LCMS_Neg", "GCMS")
    )
    expect_identical(
        skipped$provenance$assay_scope,
        c("LCMS_Pos", "LCMS_Neg", "GCMS")
    )
    expect_true(all(vapply(
        skipped$provenance$assays,
        function(assay) identical(assay$status, "no_match") &&
            identical(assay$match_count, 0L),
        logical(1)
    )))
    expect_match(
        skipped$provenance$failure_reason,
        "No internal standards found",
        fixed = TRUE
    )
    skip_record <- manager$getStateAudit("lipid_internal_standards_skipped")
    expect_identical(skip_record$status, "skipped")
    expect_identical(
        skip_record$parameters$failure_reason,
        skipped$provenance$failure_reason
    )
    expect_identical(methods::validObject(finalized, test = TRUE), TRUE)
})

test_that("lipidomics QC resumes and reverts exact mixed assay generations", {
    lipid042SkipDependencies()
    root <- withr::local_tempdir()
    context <- lipid042Context(root, "omics-art-042-resume")
    manager <- lipid042Manager(context)
    workflow <- lipid042Workflow(manager, context)
    parent <- lipid042Object(include_duplicates = FALSE)
    lipid042SaveInitial(manager, parent)
    child <- lipid042Filter(parent)
    saved <- saveLipidQcState(
        workflow,
        parent,
        child,
        "intensity_filter",
        "lipid_intensity_filtered",
        workflow$config_list,
        "resume intensity",
        parameters = list(percentile = 10, proportion = 0.5),
        transformation_type = "filter"
    )
    manager$close()

    resumed <- lipid042Manager(context)
    withr::defer(resumed$close())
    expect_identical(resumed$getState(), saved)
    expect_identical(
        resumed$getHistory(),
        c("initial", "lipid_raw_data_s4", "lipid_intensity_filtered")
    )
    expect_identical(resumed$revertToState("lipid_raw_data_s4"), parent)
    expect_identical(resumed$getState(), parent)
})

test_that("fresh lipidomics dual-write keeps memory and artifact QC aligned", {
    lipid042SkipDependencies()
    root <- withr::local_tempdir()
    context <- lipid042Context(root, "omics-art-042-dual")
    artifact_manager <- lipid042Manager(context)
    parent <- lipid042Object(include_duplicates = FALSE)
    lipid042SaveInitial(artifact_manager, parent)
    artifact_manager$close()
    memory_manager <- WorkflowState$new()
    memory_manager$setWorkflowType("lipidomics_standard")
    memory_manager$saveState(
        "lipid_raw_data_s4",
        parent,
        list(stage = "design"),
        "lipidomics QC parent"
    )
    workflow <- lipid042Workflow(memory_manager, context)
    child <- lipid042Filter(parent)
    saved <- saveLipidQcState(
        workflow,
        parent,
        child,
        "intensity_filter",
        "lipid_intensity_filtered",
        workflow$config_list,
        "dual manager intensity",
        parameters = list(percentile = 10, proportion = 0.5),
        transformation_type = "filter"
    )
    artifact_manager <- lipid042Manager(context)
    withr::defer(artifact_manager$close())

    expect_identical(memory_manager$getState(), saved)
    expect_identical(artifact_manager$getState(), saved)
    expect_identical(
        memory_manager$getCurrentStateName(),
        artifact_manager$getCurrentStateName()
    )
    artifact_manager$close()
    reverted <- revertLipidQcState(workflow, "lipid_raw_data_s4")
    expect_identical(reverted, parent)
    expect_identical(memory_manager$getState(), parent)
    artifact_manager <- lipid042Manager(context)
    expect_identical(artifact_manager$getState(), parent)
})

test_that("lipidomics QC commit failure leaves its parent generation active", {
    lipid042SkipDependencies()
    context <- lipid042Context(withr::local_tempdir(), "omics-art-042-failure")
    manager <- lipid042Manager(context)
    withr::defer(manager$close())
    workflow <- lipid042Workflow(manager, context)
    parent <- lipid042Object(include_duplicates = FALSE)
    lipid042SaveInitial(manager, parent)
    parent_generation <- manager$getCurrentGenerationId()
    child <- lipid042Filter(parent)
    injector <- function(stage, context) {
        if (identical(stage, "before_state_registry_commit")) {
            rlang::abort(
                "injected OMICS-ART-042 failure",
                class = "multischolar_test_artifact_failure"
            )
        }
        invisible(context)
    }

    expect_error(
        saveLipidQcState(
            workflow,
            parent,
            child,
            "intensity_filter",
            "lipid_intensity_filtered",
            workflow$config_list,
            "failure injection",
            parameters = list(percentile = 10, proportion = 0.5),
            transformation_type = "filter",
            failure_injector = injector
        ),
        class = "multischolar_test_artifact_failure"
    )
    expect_identical(manager$getCurrentGenerationId(), parent_generation)
    expect_identical(manager$getCurrentStateName(), "lipid_raw_data_s4")
    expect_identical(manager$getState(), parent)
    expect_null(workflow$artifact_stage_results$intensity_filter)
})

test_that("lipidomics artifact progress is session-owned and released", {
    had_global <- exists(
        "filtering_progress_lipidomics",
        envir = .GlobalEnv,
        inherits = FALSE
    )
    previous <- if (had_global) {
        get("filtering_progress_lipidomics", envir = .GlobalEnv)
    } else {
        NULL
    }
    withr::defer({
        if (had_global) {
            assign(
                "filtering_progress_lipidomics",
                previous,
                envir = .GlobalEnv
            )
        } else if (exists(
            "filtering_progress_lipidomics",
            envir = .GlobalEnv,
            inherits = FALSE
        )) {
            rm("filtering_progress_lipidomics", envir = .GlobalEnv)
        }
    })
    if (!had_global && exists(
        "filtering_progress_lipidomics",
        envir = .GlobalEnv,
        inherits = FALSE
    )) {
        rm("filtering_progress_lipidomics", envir = .GlobalEnv)
    }
    owner <- new.env(parent = emptyenv())
    progress <- getFilteringProgressLipidomics(owner)
    expect_s4_class(progress, "FilteringProgressLipidomics")
    expect_identical(owner$filtering_progress_lipidomics, progress)
    if (!had_global) {
        expect_false(exists(
            "filtering_progress_lipidomics",
            envir = .GlobalEnv,
            inherits = FALSE
        ))
    }
    releaseFilteringProgressLipidomics(owner)
    expect_null(owner$filtering_progress_lipidomics)
})

test_that("artifact QC tracking receives session-owned paths and progress", {
    workflow <- new.env(parent = emptyenv())
    workflow$state_manager <- structure(
        new.env(parent = emptyenv()),
        class = "ArtifactWorkflowState"
    )
    output_dir <- withr::local_tempdir()
    initializeLipidQcProjectContext(
        workflow,
        list(time_dir = output_dir),
        "lipidomics",
        "session-context"
    )
    captured <- new.env(parent = emptyenv())
    result <- invokeLipidQcTracking(
        function(
            theObject,
            omics_type,
            time_dir,
            progress_owner
        ) {
            captured$args <- as.list(environment())
            "tracked"
        },
        list(theObject = "s4-token"),
        workflow
    )

    expect_identical(result, "tracked")
    expect_identical(captured$args$theObject, "s4-token")
    expect_identical(captured$args$omics_type, "lipidomics")
    expect_identical(captured$args$time_dir, output_dir)
    expect_identical(captured$args$progress_owner, workflow)
})

test_that("artifact QC session cleanup releases progress and project context", {
    workflow <- new.env(parent = emptyenv())
    workflow$state_manager <- structure(
        new.env(parent = emptyenv()),
        class = "ArtifactWorkflowState"
    )
    initializeLipidQcProjectContext(
        workflow,
        list(time_dir = withr::local_tempdir()),
        "lipidomics",
        "cleanup"
    )
    cleanup <- NULL
    session <- list(onSessionEnded = \(callback) cleanup <<- callback)

    expect_true(initializeLipidQcProgressOwnership(workflow, session))
    expect_s4_class(
        workflow$filtering_progress_lipidomics,
        "FilteringProgressLipidomics"
    )
    expect_true(is.list(workflow$lipid_qc_project_context))
    cleanup()
    expect_null(workflow$filtering_progress_lipidomics)
    expect_null(workflow$lipid_qc_project_context)
    expect_null(workflow$lipid_qc_itsd_provenance)
})

test_that("frozen mixed CI workload reaches every lipidomics QC state", {
    lipid042SkipDependencies()
    repo <- normalizePath(testthat::test_path("..", ".."), mustWork = TRUE)
    adapter_env <- new.env(parent = globalenv())
    sys.source(
        file.path(repo, "tools", "profiling", "omics_workload_lipidomics.R"),
        envir = adapter_env
    )
    contract <- jsonlite::read_json(
        file.path(
            repo,
            "tests",
            "testdata",
            "omics-parity",
            "lipidomics",
            "workloads",
            "mixed-public-ci-v1.json"
        ),
        simplifyVector = FALSE
    )
    root <- withr::local_tempdir()
    do.call(RNGkind, as.list(unlist(contract$rng$kind, use.names = FALSE)))
    set.seed(as.integer(contract$rng$seed))
    prepared <- adapter_env$lipidWorkloadPrepare(list(
        contract = contract,
        run_dir = root
    ))
    generated <- utils::read.delim(
        prepared$payload_path,
        check.names = FALSE,
        stringsAsFactors = FALSE
    )
    truth <- jsonlite::read_json(prepared$truth_path, simplifyVector = TRUE)
    sample_columns <- unlist(truth$sample_ids, use.names = FALSE)
    assay_names <- unlist(truth$assays, use.names = FALSE)
    assays <- lapply(assay_names, function(assay_name) {
        value <- generated[generated$assay == assay_name, , drop = FALSE]
        value$assay <- NULL
        rownames(value) <- NULL
        value
    })
    names(assays) <- assay_names
    design <- data.frame(
        Run = sample_columns,
        group = unlist(truth$group_assignments, use.names = FALSE),
        batch = unlist(truth$batch_assignments, use.names = FALSE),
        tech_rep_group = paste0("rep_", seq_along(sample_columns)),
        stringsAsFactors = FALSE
    )
    object <- createLipidomicsAssayData(
        assays,
        design,
        lipid_id_column = "lipid_id",
        annotation_id_column = "lipid_class",
        sample_id = "Run",
        group_id = "group",
        technical_replicate_id = "tech_rep_group",
        internal_standard_regex = "^IS_",
        args = list(evidence_class = "generated_scaling")
    )
    context <- lipid042Context(root, "omics-art-042-frozen-ci")
    manager <- lipid042Manager(context)
    withr::defer(manager$close())
    workflow <- lipid042Workflow(manager, context)
    lipid042SaveInitial(manager, object)

    filtered <- lipidIntensityFiltering(
        object,
        lipids_intensity_cutoff_percentile = 10,
        lipids_proportion_of_samples_below_cutoff = 0.5
    )
    current <- saveLipidQcState(
        workflow,
        object,
        filtered,
        "intensity_filter",
        "lipid_intensity_filtered",
        workflow$config_list,
        "frozen workload intensity filter",
        parameters = list(percentile = 10, proportion = 0.5),
        transformation_type = "filter"
    )
    resolved <- lipid042ResolveDuplicates(current)
    after_duplicates <- current
    after_duplicates@lipid_data <- resolved$resolvedAssayList
    current <- saveLipidQcState(
        workflow,
        current,
        after_duplicates,
        "duplicate_resolution",
        "lipid_duplicates_resolved",
        workflow$config_list,
        "frozen workload duplicate resolution",
        parameters = lipidQcDuplicateProvenance(
            current,
            after_duplicates,
            resolved$statsList
        ),
        transformation_type = "materialization"
    )
    analysis <- lipid042ItsdAnalysis(current, "^IS_")
    recordLipidQcItsdAnalysis(
        workflow,
        current,
        analysis,
        "^IS_"
    )
    current <- manager$getState()
    finalized <- saveLipidQcState(
        workflow,
        current,
        current,
        "qc_finalization",
        "lipid_qc_complete",
        workflow$config_list,
        "frozen workload QC complete",
        parameters = list(
            evidence_class = "generated_scaling",
            claim_boundary = paste(
                "declared invariants and resource behavior only;",
                "not parser or biological validation"
            )
        ),
        status = "complete",
        transformation_type = "finalization"
    )

    expect_identical(manager$getState(), finalized)
    expect_identical(methods::validObject(finalized, test = TRUE), TRUE)
    expect_identical(names(finalized@lipid_data), assay_names)
    expect_identical(finalized@design_matrix, design)
    expect_true(all(c(
        "intensity_filter", "duplicate_resolution",
        "internal_standard_analysis", "qc_finalization"
    ) %in% names(workflow$artifact_stage_results)))
    expect_identical(
        workflow$artifact_stage_results$internal_standard_analysis$metrics$new_artifact_bytes,
        0
    )
    expect_identical(
        workflow$artifact_stage_results$qc_finalization$metrics$new_artifact_bytes,
        0
    )
    expect_identical(
        manager$getStateAudit("lipid_qc_complete")$parameters$evidence_class,
        "generated_scaling"
    )
})

test_that("production intensity seam commits an assay-selection generation", {
    lipid042SkipDependencies()
    context <- lipid042Context(
        withr::local_tempdir(),
        "omics-art-042-intensity-seam"
    )
    manager <- lipid042Manager(context)
    withr::defer(manager$close())
    workflow <- lipid042Workflow(manager, context)
    parent <- lipid042Object(include_duplicates = FALSE)
    lipid042SaveInitial(manager, parent)
    output <- new.env(parent = emptyenv())

    registerLipidIntensityApplyFilterObserver(
        input = list(
            apply_filter = 1L,
            intensity_cutoff_percentile = 10L,
            proportion_below_cutoff = 0.5
        ),
        workflowData = workflow,
        output = output,
        filterStatsVal = \(value) invisible(value),
        filterPlotVal = \(value) invisible(value),
        omicType = "lipidomics",
        observeEventFn = \(event, handler) force(handler),
        reqFn = \(value) {
            if (is.null(value)) stop("required value missing")
            invisible(value)
        },
        showNotificationFn = \(...) invisible(NULL),
        removeNotificationFn = \(...) invisible(NULL),
        renderTextFn = \(value) value,
        logInfoFn = \(...) invisible(NULL),
        logWarnFn = \(...) invisible(NULL),
        logErrorFn = \(...) invisible(NULL),
        lipidIntensityFilteringFn = \(...) lipid042Filter(parent),
        updateLipidFilteringFn = \(...) NULL
    )

    expect_identical(manager$getCurrentStateName(), "lipid_intensity_filtered")
    expect_identical(manager$getState(), lipid042Filter(parent))
    expect_identical(
        workflow$artifact_stage_results$intensity_filter$representation,
        "assay_selection"
    )
})

test_that("production duplicate and finalization seams retain exact S4 state", {
    lipid042SkipDependencies()
    context <- lipid042Context(
        withr::local_tempdir(),
        "omics-art-042-finalization-seam"
    )
    manager <- lipid042Manager(context)
    withr::defer(manager$close())
    workflow <- lipid042Workflow(manager, context)
    parent <- lipid042Object(include_duplicates = TRUE, include_itsd = FALSE)
    lipid042SaveInitial(manager, parent)

    resolved <- handleLipidDuplicateResolution(
        workflow,
        "lipidomics",
        reqFn = \(value) {
            if (is.null(value)) stop("required value missing")
            invisible(value)
        },
        updateFilteringFn = \(...) NULL,
        logWarnFn = \(...) invisible(NULL)
    )
    expect_identical(manager$getState(), resolved$currentS4)
    expect_identical(
        workflow$artifact_stage_results$duplicate_resolution$representation,
        "materialized"
    )

    output <- new.env(parent = emptyenv())
    runLipidQcS4FinalizeWorkflow(
        workflowData = workflow,
        omicType = "lipidomics",
        filterPlot = \(value) invisible(value),
        output = output,
        updateTrackingPlotFn = \(
            currentS4,
            omicType,
            setFilterPlotFn,
            workflowData
        ) {
            setFilterPlotFn(NULL)
        },
        reportFinalizeSuccessFn = \(...) invisible(NULL),
        reportFinalizeErrorFn = \(error) stop(error)
    )

    expect_identical(manager$getCurrentStateName(), "lipid_qc_complete")
    expect_identical(manager$getState(), resolved$currentS4)
    expect_identical(methods::validObject(manager$getState(), test = TRUE), TRUE)
    expect_identical(
        workflow$artifact_stage_results$qc_finalization$representation,
        "state_reference"
    )
    expect_identical(
        workflow$artifact_stage_results$qc_finalization$metrics$new_artifact_bytes,
        0
    )
})
