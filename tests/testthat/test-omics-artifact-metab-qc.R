metab033SkipDependencies <- function() {
    for (package in c("DBI", "duckdb", "filelock", "arrow")) {
        testthat::skip_if_not_installed(package)
    }
}

metab033Capability <- function() {
    capabilities <- workflowCapabilityCatalogue()
    capability <- capabilities[[
        "metabolomics.custom.metabolite.standard.v1"
    ]]
    capability$artifact_eligible <- TRUE
    capability$auto_eligible <- FALSE
    capability$maximum_artifact_rollout <- "dual_write"
    capability$explicit_maximum_artifact_rollout <- "dual_write"
    capability
}

metab033Context <- function(root, project_id = "omics-art-033") {
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
        capabilities = list(metab033Capability())
    )
    context
}

metab033Manager <- function(context) {
    descriptor <- artifactMetabolomicsWorkflowDescriptor()
    manager <- ArtifactWorkflowState$new(
        workflow_context = context,
        dehydrate_fn = dehydrateMetabolomicsS4Artifact,
        validate_bundle_fn = validateMetabolomicsS4Bundle,
        hydrate_fn = hydrateMetabolomicsS4Artifact,
        descriptor_contract = artifactStageDescriptorContract(descriptor)
    )
    manager$setWorkflowType("metabolomics_standard")
    manager
}

metab033Object <- function(
    include_duplicates = FALSE,
    include_itsd = TRUE,
    layout = "combined"
) {
    module_ci_metab_qc_object(
        layout = layout,
        include_duplicates = include_duplicates,
        include_itsd = include_itsd
    )
}

metab033Workflow <- function(manager, context) {
    workflow <- new.env(parent = emptyenv())
    workflow$state_manager <- manager
    workflow$workflow_context <- context
    workflow$config_list <- list(
        globalParameters = list(workflow_type = "metabolomics_standard"),
        qc = list(mode = "artifact-ticket")
    )
    workflow$artifact_stage_results <- list()
    workflow$processing_log <- list()
    workflow$tab_status <- list(
        setup_import = "complete",
        design_matrix = "complete",
        quality_control = "pending",
        normalization = "disabled"
    )
    workflow$data_format <- "custom"
    workflow$data_type <- "metabolite"
    workflow
}

metab033SaveInitial <- function(manager, object) {
    manager$saveState(
        "metab_raw_data_s4",
        object,
        list(stage = "design"),
        "metabolomics QC parent"
    )
    invisible(object)
}

metab033Manifest <- function(context, manager, state_name = NULL) {
    if (is.null(state_name)) state_name <- manager$getCurrentStateName()
    row <- manager$states[[state_name]]
    store <- newArtifactStore(
        context$getPaths(),
        context$getIdentity()$project_id
    )
    artifactWorkflowStateReadManifest(store, row$manifest_relative_path)
}

metab033Metadata <- function(context, manager, state_name = NULL) {
    manifest <- metab033Manifest(context, manager, state_name)
    artifactWorkflowStateUnserializeMetadata(
        manifest$data$metadata_json,
        "OMICS-ART-033 metadata"
    )
}

metab033ReadRef <- function(context, ref) {
    store <- newArtifactStore(
        context$getPaths(),
        context$getIdentity()$project_id
    )
    sidecar <- artifactStoreReadSidecar(
        store,
        artifactStoreManagedPaths(
            store,
            ref$logical_key,
            ref$artifact_id
        )$sidecar,
        validate_payload = TRUE
    )
    payload <- arrow::read_parquet(
        artifactStoreResolveFile(store, ref$relative_path, must_exist = TRUE),
        as_data_frame = FALSE
    )
    decodeArtifactRectangular(payload, sidecar$codec_metadata)
}

metab033Filter <- function(object) {
    filtered <- object
    filtered@metabolite_data <- lapply(
        object@metabolite_data,
        function(assay) {
            value <- assay[c(4L, 1L, 3L), , drop = FALSE]
            rownames(value) <- NULL
            value
        }
    )
    methods::validObject(filtered)
    filtered
}

metab033PlotBuildData <- function(object, plot_fn, ...) {
    plots <- suppressWarnings(suppressMessages(plot_fn(object, ...)))
    lapply(plots, function(plot) ggplot2::ggplot_build(plot)$data)
}

metab033ReviewedObject <- function() {
    repo <- normalizePath(testthat::test_path("..", ".."), mustWork = TRUE)
    assays <- list(
        GCMS = utils::read.delim(
            file.path(
                repo,
                "tests",
                "testdata",
                "e2e",
                "metab_combined",
                "combined_gcms_features.tsv"
            ),
            check.names = FALSE,
            stringsAsFactors = FALSE
        ),
        LCMS_Pos = utils::read.delim(
            file.path(
                repo,
                "tests",
                "testdata",
                "e2e",
                "metab_combined",
                "combined_lcms_features.tsv"
            ),
            check.names = FALSE,
            stringsAsFactors = FALSE
        )
    )
    samples <- grep("^(WT|KO)_", names(assays[[1L]]), value = TRUE)
    design <- data.frame(
        Run = samples,
        group = sub("_.*$", "", samples),
        batch = rep(c("B1", "B2"), length.out = length(samples)),
        tech_rep_group = samples,
        stringsAsFactors = FALSE
    )
    createMetaboliteAssayData(
        assays,
        design,
        metabolite_id_column = "Feature.Name",
        annotation_id_column = "Feature.Name",
        sample_id = "Run",
        group_id = "group",
        technical_replicate_id = "tech_rep_group",
        internal_standard_regex = "^IS_",
        args = list(evidence_class = "independently_reviewed_fixture")
    )
}

test_that("metabolomics intensity recipes bind each exact parent assay", {
    metab033SkipDependencies()
    context <- metab033Context(withr::local_tempdir())
    manager <- metab033Manager(context)
    withr::defer(manager$close())
    workflow <- metab033Workflow(manager, context)
    parent <- metab033Object(include_duplicates = FALSE)
    metab033SaveInitial(manager, parent)
    parent_manifest <- metab033Manifest(context, manager)
    child <- metab033Filter(parent)

    saved <- saveMetabQcState(
        workflow,
        parent,
        child,
        stage_id = "intensity_filter",
        state_name = "metab_intensity_filtered",
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
    manifest <- metab033Manifest(context, manager)
    metadata <- metab033Metadata(context, manager)
    derivation <- metadata$derivation

    expect_identical(hydrated, saved)
    expect_identical(methods::validObject(hydrated, test = TRUE), TRUE)
    expect_identical(derivation$representation, "assay_selection")
    expect_identical(
        vapply(derivation$assays, `[[`, character(1), "assay_name"),
        names(parent@metabolite_data)
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
        recipe <- metab033ReadRef(
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

test_that("metabolomics QC materializes duplicate resolution and records policy", {
    metab033SkipDependencies()
    context <- metab033Context(
        withr::local_tempdir(),
        "omics-art-033-duplicates"
    )
    manager <- metab033Manager(context)
    withr::defer(manager$close())
    workflow <- metab033Workflow(manager, context)
    parent <- metab033Object(include_duplicates = TRUE, include_itsd = FALSE)
    metab033SaveInitial(manager, parent)
    resolution <- resolveMetabDuplicateAssayData(
        parent@metabolite_data,
        parent@metabolite_id_column,
        logWarnFn = function(...) invisible(NULL)
    )
    child <- parent
    child@metabolite_data <- resolution$resolvedAssayList

    saved <- saveMetabQcState(
        workflow,
        parent,
        child,
        stage_id = "duplicate_resolution",
        state_name = "metab_duplicates_resolved",
        config_object = workflow$config_list,
        description = "ticket duplicate resolution",
        parameters = list(
            method = "highest_mean_intensity",
            representative = "maximum_row_mean",
            tie_break = "first_stable_input_feature",
            aggregation = "none",
            per_assay = resolution$statsList
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
    expect_identical(record$parameters$tie_break, "first_stable_input_feature")
    expect_true(all(vapply(
        saved@metabolite_data,
        function(assay) anyDuplicated(assay[[saved@metabolite_id_column]]) == 0L,
        logical(1)
    )))
})

test_that("metabolomics assay recipes chain through recipe-backed parents", {
    metab033SkipDependencies()
    context <- metab033Context(withr::local_tempdir(), "omics-art-033-chain")
    manager <- metab033Manager(context)
    withr::defer(manager$close())
    workflow <- metab033Workflow(manager, context)
    parent <- metab033Object(include_duplicates = FALSE)
    metab033SaveInitial(manager, parent)
    first <- metab033Filter(parent)
    first <- saveMetabQcState(
        workflow,
        parent,
        first,
        "intensity_filter",
        "metab_intensity_filtered",
        workflow$config_list,
        "first filter",
        parameters = list(percentile = 10, proportion = 0.5),
        transformation_type = "filter"
    )
    first_manifest <- metab033Manifest(context, manager)
    first_metadata <- metab033Metadata(context, manager)
    second <- first
    second@metabolite_data <- lapply(first@metabolite_data, function(assay) {
        value <- assay[c(2L, 1L), , drop = FALSE]
        rownames(value) <- NULL
        value
    })
    second <- saveMetabQcState(
        workflow,
        first,
        second,
        "intensity_filter",
        "metab_intensity_refiltered",
        workflow$config_list,
        "second filter",
        parameters = list(percentile = 20, proportion = 0.5),
        transformation_type = "filter"
    )
    second_metadata <- metab033Metadata(context, manager)

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

test_that("metabolomics internal-standard state retains assay provenance", {
    metab033SkipDependencies()
    context <- metab033Context(withr::local_tempdir(), "omics-art-033-itsd")
    manager <- metab033Manager(context)
    withr::defer(manager$close())
    workflow <- metab033Workflow(manager, context)
    parent <- metab033Object(include_duplicates = FALSE, include_itsd = TRUE)
    metab033SaveInitial(manager, parent)
    analysis <- analyzeMetabQcItsdData(
        parent,
        inputPattern = "^IS_",
        logInfoFn = function(...) invisible(NULL)
    )

    result <- recordMetabQcItsdAnalysis(
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
    expect_identical(provenance$assay_scope, names(parent@metabolite_data))
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
    manifest <- metab033Manifest(context, manager)
    metadata <- metab033Metadata(context, manager)
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
            hydrateMetabolomicsS4Artifact,
            character()
        ),
        class = "multischolar_artifact_state_reference_parent_mismatch"
    )
})

test_that("metabolomics artifact hydration preserves QC plot data", {
    metab033SkipDependencies()
    context <- metab033Context(withr::local_tempdir(), "omics-art-033-plots")
    manager <- metab033Manager(context)
    withr::defer(manager$close())
    workflow <- metab033Workflow(manager, context)
    parent <- metab033Object(
        include_duplicates = FALSE,
        include_itsd = FALSE,
        layout = "gc"
    )
    metab033SaveInitial(manager, parent)
    child <- metab033Filter(parent)
    saved <- saveMetabQcState(
        workflow,
        parent,
        child,
        "intensity_filter",
        "metab_intensity_filtered",
        workflow$config_list,
        "plot parity",
        parameters = list(percentile = 10, proportion = 0.5),
        transformation_type = "filter"
    )
    hydrated <- manager$getState()

    expect_identical(
        metab033PlotBuildData(saved, plotPca, grouping_variable = "group"),
        metab033PlotBuildData(hydrated, plotPca, grouping_variable = "group")
    )
    expect_identical(
        metab033PlotBuildData(saved, plotRle, grouping_variable = "group"),
        metab033PlotBuildData(hydrated, plotRle, grouping_variable = "group")
    )
    expect_identical(
        metab033PlotBuildData(saved, plotDensity, grouping_variable = "group"),
        metab033PlotBuildData(hydrated, plotDensity, grouping_variable = "group")
    )
    expect_identical(
        metab033PlotBuildData(
            saved,
            plotPearson,
            tech_rep_remove_regex = "",
            correlation_group = "tech_rep_group"
        ),
        metab033PlotBuildData(
            hydrated,
            plotPearson,
            tech_rep_remove_regex = "",
            correlation_group = "tech_rep_group"
        )
    )
})

test_that("reviewed mixed fixtures retain scientific state and explicit IS skip", {
    metab033SkipDependencies()
    context <- metab033Context(withr::local_tempdir(), "omics-art-033-reviewed")
    manager <- metab033Manager(context)
    withr::defer(manager$close())
    workflow <- metab033Workflow(manager, context)
    parent <- metab033ReviewedObject()
    metab033SaveInitial(manager, parent)
    expected <- metaboliteIntensityFiltering(
        parent,
        metabolites_intensity_cutoff_percentile = 10,
        metabolites_proportion_of_samples_below_cutoff = 0.5
    )
    current <- saveMetabQcState(
        workflow,
        parent,
        expected,
        "intensity_filter",
        "metab_intensity_filtered",
        workflow$config_list,
        "reviewed fixture intensity",
        parameters = list(percentile = 10, proportion = 0.5),
        transformation_type = "filter"
    )
    error <- tryCatch(
        {
            analyzeMetabQcItsdData(
                current,
                inputPattern = "^IS_",
                logInfoFn = function(...) invisible(NULL)
            )
            NULL
        },
        error = identity
    )
    expect_s3_class(error, "error")
    skipped <- recordMetabQcItsdFailure(
        workflow,
        current,
        "^IS_",
        error
    )
    current <- manager$getState()
    finalized <- saveMetabQcState(
        workflow,
        current,
        current,
        "qc_finalization",
        "metab_qc_complete",
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
    expect_identical(names(finalized@metabolite_data), c("GCMS", "LCMS_Pos"))
    expect_identical(skipped$provenance$assay_scope, c("GCMS", "LCMS_Pos"))
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
    skip_record <- manager$getStateAudit("metab_internal_standards_skipped")
    expect_identical(skip_record$status, "skipped")
    expect_identical(
        skip_record$parameters$failure_reason,
        skipped$provenance$failure_reason
    )
    expect_identical(methods::validObject(finalized, test = TRUE), TRUE)
})

test_that("metabolomics QC resumes and reverts exact mixed assay generations", {
    metab033SkipDependencies()
    root <- withr::local_tempdir()
    context <- metab033Context(root, "omics-art-033-resume")
    manager <- metab033Manager(context)
    workflow <- metab033Workflow(manager, context)
    parent <- metab033Object(include_duplicates = FALSE)
    metab033SaveInitial(manager, parent)
    child <- metab033Filter(parent)
    saved <- saveMetabQcState(
        workflow,
        parent,
        child,
        "intensity_filter",
        "metab_intensity_filtered",
        workflow$config_list,
        "resume intensity",
        parameters = list(percentile = 10, proportion = 0.5),
        transformation_type = "filter"
    )
    manager$close()

    resumed <- metab033Manager(context)
    withr::defer(resumed$close())
    expect_identical(resumed$getState(), saved)
    expect_identical(
        resumed$getHistory(),
        c("initial", "metab_raw_data_s4", "metab_intensity_filtered")
    )
    expect_identical(resumed$revertToState("metab_raw_data_s4"), parent)
    expect_identical(resumed$getState(), parent)
})

test_that("fresh metabolomics dual-write keeps memory and artifact QC aligned", {
    metab033SkipDependencies()
    root <- withr::local_tempdir()
    context <- metab033Context(root, "omics-art-033-dual")
    artifact_manager <- metab033Manager(context)
    parent <- metab033Object(include_duplicates = FALSE)
    metab033SaveInitial(artifact_manager, parent)
    artifact_manager$close()
    memory_manager <- WorkflowState$new()
    memory_manager$setWorkflowType("metabolomics_standard")
    memory_manager$saveState(
        "metab_raw_data_s4",
        parent,
        list(stage = "design"),
        "metabolomics QC parent"
    )
    workflow <- metab033Workflow(memory_manager, context)
    child <- metab033Filter(parent)
    saved <- saveMetabQcState(
        workflow,
        parent,
        child,
        "intensity_filter",
        "metab_intensity_filtered",
        workflow$config_list,
        "dual manager intensity",
        parameters = list(percentile = 10, proportion = 0.5),
        transformation_type = "filter"
    )
    artifact_manager <- metab033Manager(context)
    withr::defer(artifact_manager$close())

    expect_identical(memory_manager$getState(), saved)
    expect_identical(artifact_manager$getState(), saved)
    expect_identical(
        memory_manager$getCurrentStateName(),
        artifact_manager$getCurrentStateName()
    )
    artifact_manager$close()
    reverted <- revertMetabQcState(workflow, "metab_raw_data_s4")
    expect_identical(reverted, parent)
    expect_identical(memory_manager$getState(), parent)
    artifact_manager <- metab033Manager(context)
    expect_identical(artifact_manager$getState(), parent)
})

test_that("metabolomics QC commit failure leaves its parent generation active", {
    metab033SkipDependencies()
    context <- metab033Context(withr::local_tempdir(), "omics-art-033-failure")
    manager <- metab033Manager(context)
    withr::defer(manager$close())
    workflow <- metab033Workflow(manager, context)
    parent <- metab033Object(include_duplicates = FALSE)
    metab033SaveInitial(manager, parent)
    parent_generation <- manager$getCurrentGenerationId()
    child <- metab033Filter(parent)
    injector <- function(stage, context) {
        if (identical(stage, "before_state_registry_commit")) {
            rlang::abort(
                "injected OMICS-ART-033 failure",
                class = "multischolar_test_artifact_failure"
            )
        }
        invisible(context)
    }

    expect_error(
        saveMetabQcState(
            workflow,
            parent,
            child,
            "intensity_filter",
            "metab_intensity_filtered",
            workflow$config_list,
            "failure injection",
            parameters = list(percentile = 10, proportion = 0.5),
            transformation_type = "filter",
            failure_injector = injector
        ),
        class = "multischolar_test_artifact_failure"
    )
    expect_identical(manager$getCurrentGenerationId(), parent_generation)
    expect_identical(manager$getCurrentStateName(), "metab_raw_data_s4")
    expect_identical(manager$getState(), parent)
    expect_null(workflow$artifact_stage_results$intensity_filter)
})

test_that("metabolomics artifact progress is session-owned and released", {
    had_global <- exists(
        "filtering_progress_metabolomics",
        envir = .GlobalEnv,
        inherits = FALSE
    )
    previous <- if (had_global) {
        get("filtering_progress_metabolomics", envir = .GlobalEnv)
    } else {
        NULL
    }
    withr::defer({
        if (had_global) {
            assign(
                "filtering_progress_metabolomics",
                previous,
                envir = .GlobalEnv
            )
        } else if (exists(
            "filtering_progress_metabolomics",
            envir = .GlobalEnv,
            inherits = FALSE
        )) {
            rm("filtering_progress_metabolomics", envir = .GlobalEnv)
        }
    })
    if (!had_global && exists(
        "filtering_progress_metabolomics",
        envir = .GlobalEnv,
        inherits = FALSE
    )) {
        rm("filtering_progress_metabolomics", envir = .GlobalEnv)
    }
    owner <- new.env(parent = emptyenv())
    progress <- getFilteringProgressMetabolomics(owner)
    expect_s4_class(progress, "FilteringProgressMetabolomics")
    expect_identical(owner$filtering_progress_metabolomics, progress)
    if (!had_global) {
        expect_false(exists(
            "filtering_progress_metabolomics",
            envir = .GlobalEnv,
            inherits = FALSE
        ))
    }
    releaseFilteringProgressMetabolomics(owner)
    expect_null(owner$filtering_progress_metabolomics)
})

test_that("artifact QC tracking receives session-owned paths and progress", {
    workflow <- new.env(parent = emptyenv())
    workflow$state_manager <- structure(
        new.env(parent = emptyenv()),
        class = "ArtifactWorkflowState"
    )
    output_dir <- withr::local_tempdir()
    initializeMetabQcProjectContext(
        workflow,
        list(time_dir = output_dir),
        "metabolomics",
        "session-context"
    )
    captured <- new.env(parent = emptyenv())
    result <- invokeMetabQcTracking(
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
    expect_identical(captured$args$omics_type, "metabolomics")
    expect_identical(captured$args$time_dir, output_dir)
    expect_identical(captured$args$progress_owner, workflow)
})

test_that("frozen mixed CI workload reaches every metabolomics QC state", {
    metab033SkipDependencies()
    repo <- normalizePath(testthat::test_path("..", ".."), mustWork = TRUE)
    adapter_env <- new.env(parent = globalenv())
    sys.source(
        file.path(repo, "tools", "profiling", "omics_workload_metabolomics.R"),
        envir = adapter_env
    )
    contract <- jsonlite::read_json(
        file.path(
            repo,
            "tests",
            "testdata",
            "omics-parity",
            "metabolomics",
            "workloads",
            "mixed-public-ci-v1.json"
        ),
        simplifyVector = FALSE
    )
    root <- withr::local_tempdir()
    prepared <- adapter_env$metabWorkloadPrepare(list(
        contract = contract,
        run_dir = root
    ))
    generated <- adapter_env$metabWorkloadRead(prepared$payload_path)
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
    object <- createMetaboliteAssayData(
        assays,
        design,
        metabolite_id_column = "feature_id",
        annotation_id_column = "annotation",
        sample_id = "Run",
        group_id = "group",
        technical_replicate_id = "tech_rep_group",
        internal_standard_regex = "^IS_",
        args = list(evidence_class = "generated_scaling")
    )
    context <- metab033Context(root, "omics-art-033-frozen-ci")
    manager <- metab033Manager(context)
    withr::defer(manager$close())
    workflow <- metab033Workflow(manager, context)
    metab033SaveInitial(manager, object)

    filtered <- metaboliteIntensityFiltering(
        object,
        metabolites_intensity_cutoff_percentile = 10,
        metabolites_proportion_of_samples_below_cutoff = 0.5
    )
    current <- saveMetabQcState(
        workflow,
        object,
        filtered,
        "intensity_filter",
        "metab_intensity_filtered",
        workflow$config_list,
        "frozen workload intensity filter",
        parameters = list(percentile = 10, proportion = 0.5),
        transformation_type = "filter"
    )
    resolved <- resolveMetabDuplicateAssayData(
        current@metabolite_data,
        current@metabolite_id_column,
        logWarnFn = function(...) invisible(NULL)
    )
    after_duplicates <- current
    after_duplicates@metabolite_data <- resolved$resolvedAssayList
    current <- saveMetabQcState(
        workflow,
        current,
        after_duplicates,
        "duplicate_resolution",
        "metab_duplicates_resolved",
        workflow$config_list,
        "frozen workload duplicate resolution",
        parameters = list(
            method = "highest_mean_intensity",
            aggregation = "none",
            per_assay = resolved$statsList
        ),
        transformation_type = "materialization"
    )
    metrics <- lapply(current@metabolite_data, function(assay) {
        selected <- assay$internal_standard %in% TRUE
        data.frame(
            is_id = assay$feature_id[selected],
            mean_intensity = rowMeans(
                assay[selected, sample_columns, drop = FALSE],
                na.rm = TRUE
            ),
            cv = 0,
            stringsAsFactors = FALSE
        )
    })
    analysis <- list(
        metrics = do.call(rbind, Map(function(value, assay_name) {
            value$assay <- assay_name
            value
        }, metrics, names(metrics))),
        metricsByAssay = metrics,
        longData = NULL,
        resultText = "generated invariant evidence",
        nIsTotal = sum(vapply(metrics, nrow, integer(1))),
        pattern = "internal_standard == TRUE"
    )
    recordMetabQcItsdAnalysis(
        workflow,
        current,
        analysis,
        "internal_standard == TRUE"
    )
    current <- manager$getState()
    finalized <- saveMetabQcState(
        workflow,
        current,
        current,
        "qc_finalization",
        "metab_qc_complete",
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
    expect_identical(names(finalized@metabolite_data), assay_names)
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
        manager$getStateAudit("metab_qc_complete")$parameters$evidence_class,
        "generated_scaling"
    )
})
