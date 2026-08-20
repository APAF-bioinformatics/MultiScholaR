artifactStateTestCapability <- function() {
    candidates <- Filter(function(capability) {
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

artifactStateTestContext <- function(project_root) {
    capability <- artifactStateTestCapability()
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
            project_id = "artifact-state-project"
        )
    )
    bindWorkflowContextFromImport(
        context,
        workflow_type = "DIA",
        input_format = "diann",
        data_level = "peptide",
        capabilities = list(capability)
    )
    context
}

artifactStateTestPeptide <- function(offset = 0) {
    peptide_data <- data.frame(
        Protein.Group = rep(c("PG1", "PG2"), each = 2L),
        Protein.Ids = rep(c("P1", "P2"), each = 2L),
        Stripped.Sequence = rep(c("PEPTIDE_A", "PEPTIDE_B"), each = 2L),
        Run = rep(c("S1", "S2"), 2L),
        Q.Value = c(0.001, 0.002, 0.003, 0.004),
        Global.Q.Value = c(0.005, 0.006, 0.007, 0.008),
        Proteotypic = c(1L, 1L, 0L, 0L),
        Precursor.Quantity = c(100, 200, 300, 400) + offset,
        Precursor.Normalised = c(10, 20, 30, 40) + offset,
        stringsAsFactors = FALSE
    )
    design <- data.frame(
        Run = c("S1", "S2"),
        group = factor(c("control", "case"), levels = c("case", "control")),
        replicates = c("R1", "R1"),
        stringsAsFactors = FALSE
    )
    object <- PeptideQuantitativeDataDiann(
        peptide_data,
        design,
        sample_id = "Run",
        group_id = "group",
        technical_replicate_id = "replicates",
        args = list(branch = "artifact-state-test")
    )
    object@peptide_matrix <- matrix(
        c(10, 30, 20, 40) + offset,
        nrow = 2L,
        dimnames = list(c("PG1%PEPTIDE_A", "PG2%PEPTIDE_B"), c("S1", "S2"))
    )
    .ensurePeptideFeatureKeyMap(object, record_migration = FALSE)
}

expectArtifactStateExact <- function(before, after) {
    expect_identical(class(after), class(before))
    for (slot_name in methods::slotNames(before)) {
        expect_identical(
            methods::slot(after, slot_name),
            methods::slot(before, slot_name),
            info = slot_name
        )
    }
    expect_identical(methods::validObject(after, test = TRUE), TRUE)
}

artifactStateSkipDependencies <- function() {
    for (package in c("DBI", "duckdb", "filelock", "arrow")) {
        testthat::skip_if_not_installed(package)
    }
}

test_that("WorkflowState factory preserves inert memory delegation", {
    artifactStateSkipDependencies()
    project_root <- withr::local_tempdir()
    context <- createWorkflowContext(
        experiment_paths = list(base_dir = project_root),
        omic_type = "proteomics",
        storage_policy = list(requested_backend = "memory")
    )
    bindWorkflowContextFromImport(
        context,
        workflow_type = "DIA",
        input_format = "diann",
        data_level = "peptide"
    )

    manager <- newWorkflowState(
        audit_enabled = FALSE,
        workflow_context = context
    )

    expect_s3_class(manager, "WorkflowState")
    expect_false(manager$isAuditEnabled())
    expect_identical(manager$getCurrentStateName(), "initial")
    expect_false(dir.exists(file.path(project_root, "state")))
    expect_s3_class(newWorkflowState(), "WorkflowState")
})

test_that("artifact states persist exact refs and reuse unchanged payloads", {
    artifactStateSkipDependencies()
    project_root <- withr::local_tempdir()
    manager <- newWorkflowState(
        workflow_context = artifactStateTestContext(project_root)
    )
    withr::defer(manager$close())
    object <- artifactStateTestPeptide()

    expect_s3_class(manager, "ArtifactWorkflowState")
    expect_identical(manager$getCurrentStateName(), "initial")
    expect_identical(manager$setWorkflowType("DIA"), "DIA")
    expect_identical(
        manager$saveState(
            "imported",
            object,
            list(stage = "import", threshold = 0.01),
            "Imported DIA peptides"
        ),
        "imported"
    )
    expectArtifactStateExact(object, manager$getState())
    expect_identical(
        manager$getStateConfig(),
        list(stage = "import", threshold = 0.01)
    )
    first_refs <- manager$states$imported$artifact_refs
    first_generation <- manager$getCurrentGenerationId()
    result <- manager$commitState(
        "metadata_updated",
        object,
        list(stage = "import", threshold = 0.02),
        "Updated metadata only",
        action_id = "metadata-only-action",
        expected_parent_generation_id = first_generation
    )

    expect_identical(result$status, "accepted")
    expect_true(result$artifacts_reused)
    expect_identical(
        manager$states$metadata_updated$artifact_refs,
        first_refs
    )
    expect_identical(manager$getCacheInfo()$entries, 1L)
    expect_false(any(vapply(manager$states, isS4, logical(1))))
    snapshot <- manager$exportState()
    expect_identical(snapshot$backend, "artifact")
    expect_false("states" %in% names(snapshot))
    expect_false(grepl(
        project_root,
        paste(capture.output(str(snapshot)), collapse = ""),
        fixed = TRUE
    ))

    manifest_path <- file.path(
        project_root,
        manager$states$metadata_updated$manifest_relative_path
    )
    manifest_text <- paste(readLines(manifest_path, warn = FALSE), collapse = "\n")
    expect_false(grepl(project_root, manifest_text, fixed = TRUE))
    expect_false(grepl('"data"[[:space:]]*:[[:space:]]*"', manifest_text))
})

test_that("row-selection hints are additive to materialized WorkflowState", {
    artifactStateSkipDependencies()
    project_root <- withr::local_tempdir()
    manager <- newWorkflowState(
        workflow_context = artifactStateTestContext(project_root)
    )
    withr::defer(manager$close())
    manager$setWorkflowType("DIA")
    parent <- artifactStateTestPeptide()
    manager$saveState("imported", parent, list(stage = "import"), "Imported")
    child <- parent
    child@peptide_data <- child@peptide_data[c(4L, 1L), , drop = FALSE]
    key_columns <- c("Protein.Group", "Stripped.Sequence", "Run")
    parent_keys <- artifactWorkflowStateEntityKeys(parent@peptide_data, key_columns)
    child_keys <- artifactWorkflowStateEntityKeys(child@peptide_data, key_columns)
    rejected_keys <- parent_keys[!parent_keys %in% child_keys]
    hint <- newArtifactRowSelectionHint(
        slot_name = "peptide_data",
        key_columns = key_columns,
        method = "workflow_state_test_filter",
        normalized_parameters = list(enabled = TRUE),
        software = list(
            name = "MultiScholaR",
            version = "test",
            source = "testthat"
        ),
        lineage = list(
            audit_record_id = "audit-test-1",
            state_name = "selected",
            parent_state = "imported",
            parent_record_id = "audit-test-0"
        ),
        rejected_reasons = stats::setNames(
            rep("test_filter", length(rejected_keys)),
            rejected_keys
        )
    )
    selected <- manager$commitState(
        "selected",
        child,
        list(stage = "filter"),
        "Selected rows",
        persistence_hint = hint,
        expected_parent_generation_id = manager$getCurrentGenerationId()
    )
    expect_identical(selected$representation, "row_selection")
    expectArtifactStateExact(child, manager$getState())

    changed <- child
    changed@peptide_data$Precursor.Normalised <-
        changed@peptide_data$Precursor.Normalised + 1
    materialized <- manager$commitState(
        "changed",
        changed,
        list(stage = "transform"),
        "Changed values",
        expected_parent_generation_id = manager$getCurrentGenerationId()
    )
    expect_identical(materialized$representation, "materialized")
    expectArtifactStateExact(changed, manager$getState())
})

test_that("artifact actions are idempotent and parent compare-and-set is strict", {
    artifactStateSkipDependencies()
    project_root <- withr::local_tempdir()
    manager <- newWorkflowState(
        workflow_context = artifactStateTestContext(project_root)
    )
    withr::defer(manager$close())
    manager$setWorkflowType("DIA")
    initial_generation <- manager$getCurrentGenerationId()
    object <- artifactStateTestPeptide()

    accepted <- manager$commitState(
        "first",
        object,
        list(step = 1L),
        "First action",
        action_id = "stable-action",
        expected_parent_generation_id = initial_generation
    )
    payloads_after_accept <- list.files(
        project_root,
        pattern = "[.]parquet$",
        recursive = TRUE
    )
    duplicate <- manager$commitState(
        "ignored-duplicate",
        artifactStateTestPeptide(100),
        list(step = 999L),
        "Duplicate action",
        action_id = "stable-action",
        expected_parent_generation_id = initial_generation
    )
    rejected <- manager$commitState(
        "late",
        artifactStateTestPeptide(200),
        list(step = 2L),
        "Out-of-order completion",
        action_id = "late-action",
        expected_parent_generation_id = initial_generation
    )

    expect_identical(accepted$status, "accepted")
    expect_true(duplicate$idempotent)
    expect_identical(duplicate$generation_id, accepted$generation_id)
    expect_identical(rejected$status, "rejected_stale_parent")
    expect_identical(manager$getCurrentGenerationId(), accepted$generation_id)
    expect_false(manager$hasState("late"))
    expect_identical(
        list.files(project_root, pattern = "[.]parquet$", recursive = TRUE),
        payloads_after_accept
    )
    expect_true(any(vapply(manager$getEvents(), function(event) {
        identical(event$event_type, "state_save_rejected")
    }, logical(1))))
})

test_that("revert preserves evidence, marks descendants stale, and permits branching", {
    artifactStateSkipDependencies()
    project_root <- withr::local_tempdir()
    manager <- newWorkflowState(
        workflow_context = artifactStateTestContext(project_root)
    )
    withr::defer(manager$close())
    manager$setWorkflowType("DIA")
    first_object <- artifactStateTestPeptide()
    second_object <- artifactStateTestPeptide(10)
    third_object <- artifactStateTestPeptide(20)
    manager$saveState("first", first_object, list(step = 1L), "First")
    first_generation <- manager$getCurrentGenerationId()
    manager$saveState("second", second_object, list(step = 2L), "Second")
    second_refs <- manager$states$second$artifact_refs
    second_paths <- vapply(second_refs, `[[`, character(1), "relative_path")

    expectArtifactStateExact(first_object, manager$revertToState("first"))
    expect_identical(manager$getCurrentGenerationId(), first_generation)
    expect_identical(manager$states$second$status, "stale")
    expect_true(all(file.exists(file.path(project_root, second_paths))))
    cache_before <- manager$getCacheInfo()
    expectArtifactStateExact(second_object, manager$getState("second"))
    expect_identical(manager$getCacheInfo()$entries, 1L)
    expect_identical(
        manager$getCacheInfo()$generation_id,
        cache_before$generation_id
    )

    manager$saveState("third", third_object, list(step = 3L), "Branch")
    expect_identical(manager$getHistory(), c("initial", "first", "third"))
    manager$saveState("first", third_object, list(step = 4L), "Repeated name")
    expectArtifactStateExact(third_object, manager$getState("first"))
    expect_identical(tail(manager$getHistory(), 1L), "first")
    expect_true(any(vapply(manager$getEvents(), function(event) {
        identical(event$event_type, "state_reverted")
    }, logical(1))))
})

test_that("unsafe metadata and future restores cannot advance artifact state", {
    artifactStateSkipDependencies()
    project_root <- withr::local_tempdir()
    manager <- newWorkflowState(
        workflow_context = artifactStateTestContext(project_root)
    )
    withr::defer(manager$close())
    manager$setWorkflowType("DIA")
    object <- artifactStateTestPeptide()
    baseline_generation <- manager$getCurrentGenerationId()
    baseline_payloads <- list.files(
        project_root,
        pattern = "[.]parquet$",
        recursive = TRUE
    )

    expect_error(
        manager$saveState(
            "unsafe",
            object,
            list(cache = new.env(parent = emptyenv())),
            "Unsafe"
        ),
        class = "multischolar_unsafe_artifact_state_metadata"
    )
    expect_error(
        manager$saveState(
            "absolute",
            object,
            list(input_path = project_root),
            "Absolute path"
        ),
        class = "multischolar_absolute_path_in_artifact_state"
    )
    expect_identical(manager$getCurrentGenerationId(), baseline_generation)
    expect_identical(
        list.files(project_root, pattern = "[.]parquet$", recursive = TRUE),
        baseline_payloads
    )

    snapshot <- manager$exportState()
    future <- snapshot
    future$schema_version <- ARTIFACT_WORKFLOW_STATE_VERSION + 1L
    expect_error(
        manager$restoreState(future),
        class = "multischolar_future_artifact_state_version"
    )
    incompatible <- snapshot
    incompatible$current_generation_id <- artifactOpaqueId("gen")
    expect_error(
        manager$restoreState(incompatible),
        class = "multischolar_unsupported_artifact_state_restore"
    )
    expect_identical(manager$getCurrentGenerationId(), baseline_generation)
})

test_that("artifact state hydrates exact current data in a fresh process", {
    artifactStateSkipDependencies()
    skip_if(Sys.which("Rscript") == "")
    project_root <- withr::local_tempdir()
    manager <- newWorkflowState(
        workflow_context = artifactStateTestContext(project_root)
    )
    manager$setWorkflowType("DIA")
    manager$saveState(
        "fresh_process",
        artifactStateTestPeptide(),
        list(stage = "resume"),
        "Fresh-process state"
    )
    expected_generation <- manager$getCurrentGenerationId()
    expect_true(manager$close())

    script <- tempfile(fileext = ".R")
    output <- tempfile(fileext = ".json")
    log <- tempfile(fileext = ".log")
    repo_root <- normalizePath(testthat::test_path("..", ".."), mustWork = TRUE)
    writeLines(c(
        "args <- commandArgs(trailingOnly = TRUE)",
        "devtools::load_all(args[[1L]], quiet = TRUE)",
        "caps <- workflowCapabilityCatalogue()",
        "caps <- Filter(function(x) x$identity$input_format == 'diann' &&",
        "    x$identity$data_level == 'peptide', caps)",
        "cap <- caps[[1L]]",
        "cap$artifact_eligible <- TRUE",
        "cap$auto_eligible <- TRUE",
        "cap$maximum_artifact_rollout <- 'dual_write'",
        "context <- createWorkflowContext(",
        "    list(base_dir = args[[2L]], omic_label = 'proteomics_study'),",
        "    'proteomics',",
        "    storage_policy = list(",
        "        requested_backend = 'artifact',",
        "        requested_rollout = 'dual_write',",
        "        migration_requested = TRUE,",
        "        project_id = 'artifact-state-project'",
        "    )",
        ")",
        "bindWorkflowContextFromImport(",
        "    context, 'DIA', 'diann', 'peptide', capabilities = list(cap)",
        ")",
        "manager <- newWorkflowState(workflow_context = context)",
        "object <- manager$getState()",
        "result <- list(",
        "    state = manager$getCurrentStateName(),",
        "    generation = manager$getCurrentGenerationId(),",
        "    class = class(object)[[1L]],",
        "    valid = identical(methods::validObject(object, test = TRUE), TRUE),",
        "    cache_entries = manager$getCacheInfo()$entries,",
        "    hydration_count = manager$getCacheInfo()$hydration_count",
        ")",
        "jsonlite::write_json(result, args[[3L]], auto_unbox = TRUE)",
        "manager$close()"
    ), script)
    status <- system2(
        Sys.which("Rscript"),
        c("--vanilla", script, repo_root, project_root, output),
        stdout = log,
        stderr = log
    )
    if (!identical(as.integer(status), 0L)) {
        testthat::fail(paste(readLines(log, warn = FALSE), collapse = "\n"))
    }
    result <- jsonlite::read_json(output, simplifyVector = TRUE)

    expect_identical(result$state, "fresh_process")
    expect_identical(result$generation, expected_generation)
    expect_identical(result$class, "PeptideQuantitativeData")
    expect_true(result$valid)
    expect_identical(result$cache_entries, 1L)
    expect_identical(result$hydration_count, 1L)
})
