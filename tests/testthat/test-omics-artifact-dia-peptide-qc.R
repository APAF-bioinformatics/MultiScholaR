omics012SkipDependencies <- function() {
    for (package in c("DBI", "duckdb", "filelock", "arrow")) {
        testthat::skip_if_not_installed(package)
    }
}

omics012Capability <- function() {
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

omics012Context <- function(project_root) {
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
            project_id = "omics-012-project"
        )
    )
    bindWorkflowContextFromImport(
        context,
        workflow_type = "DIA",
        input_format = "diann",
        data_level = "peptide",
        capabilities = list(omics012Capability())
    )
    context
}

omics012Peptide <- function(rows = 8L, duplicate_key = FALSE) {
    runs <- rep(c("S1", "S2", "S3", "S4"), length.out = rows)
    peptide_data <- data.frame(
        Protein.Group = paste0("PG", rep(seq_len(4L), each = 2L))[seq_len(rows)],
        Protein.Ids = paste0("P", rep(seq_len(4L), each = 2L))[seq_len(rows)],
        Stripped.Sequence = paste0("PEPTIDE_", seq_len(rows)),
        Modified.Sequence = paste0("_PEPTIDE_", seq_len(rows), "_"),
        Precursor.Id = paste0("precursor-", seq_len(rows)),
        Precursor.Charge = rep(c(2L, 3L), length.out = rows),
        Run = runs,
        Q.Value = seq(0.001, 0.008, length.out = rows),
        Global.Q.Value = seq(0.002, 0.009, length.out = rows),
        Global.PG.Q.Value = seq(0.003, 0.01, length.out = rows),
        Proteotypic = rep(1L, rows),
        Precursor.Quantity = seq(100, by = 10, length.out = rows),
        Precursor.Normalised = seq(10, by = 1, length.out = rows),
        identification_peptide_count = rep(2L, rows),
        identification_peptidoform_count = rep(2L, rows),
        stringsAsFactors = FALSE
    )
    if (isTRUE(duplicate_key) && rows >= 2L) {
        key_columns <- c(
            "Protein.Group", "Protein.Ids", "Stripped.Sequence",
            "Modified.Sequence", "Precursor.Id", "Precursor.Charge", "Run"
        )
        peptide_data[2L, key_columns] <- peptide_data[1L, key_columns]
    }
    design <- data.frame(
        Run = c("S1", "S2", "S3", "S4"),
        group = factor(c("control", "control", "case", "case")),
        replicates = c("R1", "R1", "R2", "R2"),
        stringsAsFactors = FALSE
    )
    object <- PeptideQuantitativeDataDiann(
        peptide_data,
        design,
        sample_id = "Run",
        group_id = "group",
        technical_replicate_id = "replicates",
        args = list(source = "omics-art-012")
    )
    object@peptide_matrix <- matrix(
        seq_len(16L),
        nrow = 4L,
        dimnames = list(paste0("feature-", seq_len(4L)), design$Run)
    )
    .ensurePeptideFeatureKeyMap(object, record_migration = FALSE)
}

omics012ExpectExact <- function(expected, actual) {
    expect_identical(class(actual), class(expected))
    for (slot_name in methods::slotNames(expected)) {
        expect_identical(
            methods::slot(actual, slot_name),
            methods::slot(expected, slot_name),
            info = slot_name
        )
    }
    expect_identical(methods::validObject(actual, test = TRUE), TRUE)
}

omics012Manager <- function(context) {
    manager <- newProtDiaArtifactStateManager(context)
    manager$setWorkflowType("DIA")
    manager
}

omics012Save <- function(
    manager,
    before,
    after,
    stage_id,
    state_name,
    transformation_type = "filter",
    status = "applied",
    failure_injector = NULL,
    workflow_data = NULL
) {
    if (is.null(workflow_data)) {
        workflow_data <- new.env(parent = emptyenv())
        workflow_data$state_manager <- manager
    }
    .savePeptideQcState(
        state_manager = manager,
        workflow_data = workflow_data,
        before = before,
        after = after,
        stage_id = stage_id,
        state_name = state_name,
        config_object = list(stage_id = stage_id),
        description = paste("Applied", stage_id),
        audit_parameters = list(threshold = 0.01, stage_id = stage_id),
        now = as.POSIXct("2026-08-20 00:00:00", tz = "UTC"),
        status = status,
        decision_reason = if (identical(status, "skipped")) {
            "no_technical_replicates"
        } else {
            NA_character_
        },
        transformation_type = transformation_type,
        failure_injector = failure_injector
    )
}

omics012Manifest <- function(context, manager, state_name = NULL) {
    if (is.null(state_name)) state_name <- manager$getCurrentStateName()
    row <- manager$states[[state_name]]
    store <- newArtifactStore(context$getPaths(), context$getIdentity()$project_id)
    artifactWorkflowStateReadManifest(store, row$manifest_relative_path)
}

omics012Metadata <- function(context, manager, state_name = NULL) {
    manifest <- omics012Manifest(context, manager, state_name)
    artifactWorkflowStateUnserializeMetadata(
        manifest$data$metadata_json,
        "OMICS-ART-012 test metadata"
    )
}

omics012ReadPayload <- function(context, manifest, reference_name) {
    store <- newArtifactStore(context$getPaths(), context$getIdentity()$project_id)
    ref <- manifest$data$artifact_refs[[reference_name]]
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

omics012FailAt <- function(target_stage) {
    force(target_stage)
    \(stage, context) {
        if (identical(stage, target_stage)) {
            rlang::abort(
                paste("injected OMICS-ART-012 failure at", stage),
                class = "multischolar_test_artifact_failure"
            )
        }
        invisible(context)
    }
}

test_that("DIA peptide-QC selection recipes bind exact immutable lineage", {
    omics012SkipDependencies()
    context <- omics012Context(withr::local_tempdir())
    manager <- omics012Manager(context)
    withr::defer(manager$close())
    parent <- omics012Peptide()
    manager$saveState("raw_data_s4", parent, list(stage = "import"), "Imported")
    parent_manifest <- omics012Manifest(context, manager)

    child <- parent
    child@peptide_data <- child@peptide_data[c(4L, 1L, 7L, 3L), , drop = FALSE]
    child@peptide_data$confidence_class <- factor(
        c("pass", "review", "pass", "review"),
        levels = c("review", "pass")
    )
    workflow_data <- new.env(parent = emptyenv())
    workflow_data$state_manager <- manager
    audited <- omics012Save(
        manager,
        parent,
        child,
        "qvalue_filter",
        "qvalue_filtered",
        workflow_data = workflow_data
    )
    hydrated <- manager$getState()
    manifest <- omics012Manifest(context, manager)
    metadata <- omics012Metadata(context, manager)
    derivation <- metadata$derivation
    recipe <- omics012ReadPayload(
        context,
        manifest,
        .ARTIFACT_ROW_SELECTION_RECIPE_REF
    )

    omics012ExpectExact(audited, hydrated)
    expect_identical(derivation$representation, "row_selection")
    expect_identical(derivation$method, "qvalue_filter")
    expect_identical(derivation$selected_count, 4L)
    expect_identical(derivation$rejected_count, 4L)
    expect_identical(derivation$patched_columns, "confidence_class")
    expect_identical(derivation$normalized_parameters$threshold, 0.01)
    expect_identical(derivation$software$name, "MultiScholaR")
    expect_identical(nrow(recipe), nrow(parent@peptide_data))
    expect_identical(sum(recipe$selection_status == "selected"), 4L)
    expect_identical(sum(recipe$selection_status == "rejected"), 4L)
    expect_true(all(nzchar(recipe$failure_reason[
        recipe$selection_status == "rejected"
    ])))
    expect_identical(
        derivation$parent_generation_id,
        parent_manifest$generation_id
    )
    expect_identical(
        derivation$lineage$parent_content_identity,
        parent_manifest$data$semantic_digest
    )
    expect_identical(
        derivation$lineage$audit_record_id,
        audited@args$peptide_qc_audit$current_record_id
    )
    result <- manager$getEvents()
    expect_true(any(vapply(result, \(event) {
        identical(event$event_type, "state_saved")
    }, logical(1))))
    expect_identical(manager$getCacheInfo()$entries, 1L)
    artifact_result <- workflow_data$artifact_stage_results$qvalue_filter
    expect_identical(artifact_result$representation, "row_selection")
    expect_gt(artifact_result$metrics$new_artifact_bytes, 0)
    expect_gt(artifact_result$metrics$hydrated_object_bytes, 0)

    recipe_artifact_id <- manifest$data$artifact_refs[[
        .ARTIFACT_ROW_SELECTION_RECIPE_REF
    ]]$artifact_id
    parent_artifact_count <- length(parent_manifest$data$artifact_refs)
    manager$close()
    registry <- newProjectRegistry(
        context$getPaths(),
        context$getIdentity()$project_id
    )
    registry_session <- openProjectRegistryReadOnly(registry)
    dependencies <- projectRegistryQuery(
        registry_session,
        "dependencies",
        filters = list(
            workflow_id = context$getIdentity()$workflow_id,
            artifact_id = recipe_artifact_id
        )
    )
    closeProjectRegistry(registry_session)
    expect_identical(nrow(dependencies), parent_artifact_count)
    expect_true(all(
        dependencies$relationship_type == "selection_parent_generation"
    ))
    resumed <- omics012Manager(context)
    withr::defer(resumed$close())
    omics012ExpectExact(audited, resumed$getState())
})

test_that("all peptide-QC stages use selection or materialization by semantics", {
    omics012SkipDependencies()
    context <- omics012Context(withr::local_tempdir())
    manager <- omics012Manager(context)
    withr::defer(manager$close())
    current <- omics012Peptide()
    manager$saveState("raw_data_s4", current, list(stage = "import"), "Imported")

    current <- omics012Save(
        manager,
        current,
        current,
        "qvalue_filter",
        "qvalue_no_op",
        status = "no_op"
    )
    expect_identical(
        omics012Metadata(context, manager)$derivation$selected_count,
        8L
    )

    filter_stages <- c(
        "qvalue_filter",
        "intensity_filter",
        "sample_filter",
        "protein_evidence_filter",
        "replicate_filter"
    )
    filter_states <- c(
        "qvalue_filtered",
        "intensity_filtered",
        "sample_filtered",
        "protein_peptide_filtered",
        "replicate_filtered"
    )
    for (index in seq_along(filter_stages)) {
        child <- current
        child@peptide_data <- child@peptide_data[-1L, , drop = FALSE]
        current <- omics012Save(
            manager,
            current,
            child,
            filter_stages[[index]],
            filter_states[[index]]
        )
        expect_identical(
            omics012Metadata(context, manager)$derivation$representation,
            "row_selection",
            info = filter_stages[[index]]
        )
        omics012ExpectExact(current, manager$getState())
    }

    rolled_up <- current
    rolled_up@peptide_data$Precursor.Normalised <-
        rolled_up@peptide_data$Precursor.Normalised * 2
    current <- omics012Save(
        manager,
        current,
        rolled_up,
        "precursor_rollup",
        "precursor_rollup",
        transformation_type = "aggregation"
    )
    expect_null(omics012Metadata(context, manager)$derivation)

    imputed <- current
    imputed@peptide_data$Precursor.Normalised[[1L]] <- 42
    current <- omics012Save(
        manager,
        current,
        imputed,
        "imputation",
        "imputed",
        transformation_type = "imputation"
    )
    expect_null(omics012Metadata(context, manager)$derivation)

    skipped <- omics012Save(
        manager,
        current,
        current,
        "imputation",
        "imputation_skipped",
        transformation_type = "imputation",
        status = "skipped"
    )
    expect_identical(
        omics012Metadata(context, manager)$derivation$representation,
        "row_selection"
    )
    omics012ExpectExact(skipped, manager$getState())
    expect_identical(manager$getCacheInfo()$entries, 1L)

    expect_true(methods::is(manager$revertToState("qvalue_filtered"),
        "PeptideQuantitativeData"
    ))
    expect_identical(manager$getCurrentStateName(), "qvalue_filtered")
    expect_true(methods::is(manager$revertToState("intensity_filtered"),
        "PeptideQuantitativeData"
    ))
})

test_that("generic row selections represent all-pass and all-fail boundaries", {
    parent <- omics012Peptide()
    key_columns <- protDiaPeptideQcStableKeyColumns(parent)
    parent_keys <- artifactWorkflowStateEntityKeys(
        parent@peptide_data,
        key_columns
    )
    base_hint <- list(
        slot_name = "peptide_data",
        key_columns = key_columns,
        method = "synthetic_filter",
        normalized_parameters = list(threshold = 1L),
        software = list(name = "test", version = "1", source = "testthat"),
        lineage = list(
            audit_record_id = "audit-1",
            state_name = "selected",
            parent_state = "parent",
            parent_record_id = "audit-0"
        )
    )
    all_pass <- do.call(newArtifactRowSelectionHint, c(
        base_hint,
        list(rejected_reasons = character())
    ))
    pass_plan <- artifactWorkflowStateSelectionPlan(parent, parent, all_pass)
    expect_identical(pass_plan$selected_count, 8L)
    expect_identical(pass_plan$rejected_count, 0L)

    child <- parent
    child@peptide_data <- child@peptide_data[0L, , drop = FALSE]
    all_fail <- do.call(newArtifactRowSelectionHint, c(
        base_hint,
        list(rejected_reasons = stats::setNames(
            rep("failed_synthetic_contract", length(parent_keys)),
            parent_keys
        ))
    ))
    fail_plan <- artifactWorkflowStateSelectionPlan(parent, child, all_fail)
    expect_identical(fail_plan$selected_count, 0L)
    expect_identical(fail_plan$rejected_count, 8L)
    expect_true(all(!is.na(fail_plan$recipe$failure_reason)))
})

test_that("ambiguous peptide identity falls back to exact materialization", {
    omics012SkipDependencies()
    context <- omics012Context(withr::local_tempdir())
    manager <- omics012Manager(context)
    withr::defer(manager$close())
    parent <- omics012Peptide(duplicate_key = TRUE)
    manager$saveState("raw_data_s4", parent, list(stage = "import"), "Imported")
    child <- parent
    child@peptide_data <- child@peptide_data[-nrow(child@peptide_data), , drop = FALSE]

    audited <- omics012Save(
        manager,
        parent,
        child,
        "sample_filter",
        "sample_filtered"
    )

    expect_null(omics012Metadata(context, manager)$derivation)
    omics012ExpectExact(audited, manager$getState())
})

test_that("artifact failures never advance the active peptide generation", {
    omics012SkipDependencies()
    context <- omics012Context(withr::local_tempdir())
    manager <- omics012Manager(context)
    withr::defer(manager$close())
    parent <- omics012Peptide()
    manager$saveState("raw_data_s4", parent, list(stage = "import"), "Imported")
    parent_generation <- manager$getCurrentGenerationId()
    child <- parent
    child@peptide_data <- child@peptide_data[-1L, , drop = FALSE]

    for (failure_stage in c(
        "after_audit_creation",
        "before_state_registry_commit"
    )) {
        expect_error(
            omics012Save(
                manager,
                parent,
                child,
                "sample_filter",
                paste0("failed_", failure_stage),
                failure_injector = omics012FailAt(failure_stage)
            ),
            class = "multischolar_test_artifact_failure",
            info = failure_stage
        )
        expect_identical(manager$getCurrentGenerationId(), parent_generation)
        expect_identical(manager$getCurrentStateName(), "raw_data_s4")
        omics012ExpectExact(parent, manager$getState())
    }
})

test_that("compatibility save failure restores the artifact parent", {
    omics012SkipDependencies()
    context <- omics012Context(withr::local_tempdir())
    parent <- omics012Peptide()
    artifact_manager <- omics012Manager(context)
    artifact_manager$saveState(
        "raw_data_s4",
        parent,
        list(stage = "import"),
        "Imported"
    )
    parent_generation <- artifact_manager$getCurrentGenerationId()
    artifact_manager$close()
    failing_manager <- list(
        getCurrentStateName = \() "raw_data_s4",
        saveState = \(...) stop("injected compatibility save failure")
    )
    workflow_data <- new.env(parent = emptyenv())
    workflow_data$state_manager <- failing_manager
    workflow_data$workflow_context <- context
    child <- parent
    child@peptide_data <- child@peptide_data[-1L, , drop = FALSE]

    expect_error(
        omics012Save(
            failing_manager,
            parent,
            child,
            "sample_filter",
            "sample_filtered",
            workflow_data = workflow_data
        ),
        "injected compatibility save failure"
    )
    reopened <- omics012Manager(context)
    withr::defer(reopened$close())
    expect_identical(reopened$getCurrentGenerationId(), parent_generation)
    expect_identical(reopened$getCurrentStateName(), "raw_data_s4")
    expect_true(reopened$hasState("sample_filtered"))
    omics012ExpectExact(parent, reopened$getState())
})

test_that("reactive dual-write and revert keep memory and artifact states exact", {
    omics012SkipDependencies()
    context <- omics012Context(withr::local_tempdir())
    parent <- omics012Peptide()
    artifact_manager <- omics012Manager(context)
    artifact_manager$saveState(
        "raw_data_s4",
        parent,
        list(stage = "import"),
        "Imported"
    )
    artifact_manager$close()
    memory_manager <- newWorkflowState()
    memory_manager$setWorkflowType("DIA")
    memory_manager$saveState(
        "raw_data_s4",
        parent,
        list(stage = "import"),
        "Imported"
    )
    workflow_data <- shiny::reactiveValues(
        state_manager = memory_manager,
        workflow_context = context
    )
    child <- parent
    child@peptide_data <- child@peptide_data[-1L, , drop = FALSE]

    audited <- shiny::isolate(omics012Save(
        memory_manager,
        parent,
        child,
        "sample_filter",
        "sample_filtered",
        workflow_data = workflow_data
    ))
    expect_identical(memory_manager$getCurrentStateName(), "sample_filtered")
    omics012ExpectExact(audited, memory_manager$getState())
    reopened <- omics012Manager(context)
    omics012ExpectExact(audited, reopened$getState())
    reopened$close()

    reverted <- shiny::isolate(revertProtDiaPeptideQcState(
        workflow_data,
        "raw_data_s4"
    ))
    omics012ExpectExact(parent, reverted)
    expect_identical(memory_manager$getCurrentStateName(), "raw_data_s4")
    reopened <- omics012Manager(context)
    withr::defer(reopened$close())
    expect_identical(reopened$getCurrentStateName(), "raw_data_s4")
    omics012ExpectExact(parent, reopened$getState())
})
