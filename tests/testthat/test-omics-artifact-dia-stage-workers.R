diaStageWorkerSkipDependencies <- function() {
    for (package in c(
        "arrow", "DBI", "duckdb", "filelock", "processx", "pkgload"
    )) {
        testthat::skip_if_not_installed(package)
    }
}

diaStageWorkerWorkflow <- function(root, backend = "artifact") {
    paths <- list(
        base_dir = root,
        project_id = paste0("dia-stage-worker-", backend),
        omic_type = "proteomics",
        omic_label = "dia-stage-worker",
        source_dir = file.path(root, "source"),
        results_dir = file.path(root, "results")
    )
    dir.create(paths$source_dir, recursive = TRUE)
    dir.create(paths$results_dir, recursive = TRUE)
    workflow <- new.env(parent = emptyenv())
    workflow$workflow_context <- createWorkflowContext(
        paths,
        "proteomics",
        "dia-stage-worker",
        storage_policy = list(
            requested_backend = backend,
            requested_rollout = "dual_write",
            migration_requested = identical(backend, "artifact"),
            project_id = paths$project_id
        )
    )
    workflow$state_manager <- WorkflowState$new()
    workflow$artifact_stage_results <- NULL
    workflow
}

diaStageWorkerSource <- function() {
    testthat::test_path("..", "testdata", "e2e", "prot_dia", "report.tsv")
}

diaStageWorkerStore <- function(workflow) {
    context <- workflow$workflow_context
    identity <- context$getIdentity()
    newArtifactStore(context$getPaths(), identity$project_id)
}

diaStageWorkerRegistryArtifacts <- function(workflow) {
    registry <- projectRegistryForContext(workflow$workflow_context)
    session <- openProjectRegistryReadOnly(registry)
    on.exit(closeProjectRegistry(session), add = TRUE)
    projectRegistryQuery(
        session,
        "artifacts",
        filters = list(
            workflow_id = workflow$workflow_context$getIdentity()$workflow_id
        )
    )
}

expectDiaStageWorkerGenerationAbsent <- function(store) {
    expect_length(artifactStoreSidecarPaths(store), 0L)
    expect_length(artifactStoreIntentPaths(store), 0L)
    expect_length(
        list.files(
            store$project_root,
            pattern = "[.]parquet$",
            recursive = TRUE,
            full.names = TRUE
        ),
        0L
    )
}

test_that("DIA-NN import stages write, verify, hydrate, then publish exactly", {
    diaStageWorkerSkipDependencies()
    root <- withr::local_tempdir()
    workflow <- diaStageWorkerWorkflow(root)
    source <- diaStageWorkerSource()
    expected <- suppressMessages(importDIANNData(source))

    staged <- stageProtDiaImportArtifacts(
        workflow_data = workflow,
        source_path = source,
        format = "diann",
        use_precursor_norm = TRUE,
        sanitize_names = FALSE
    )
    withr::defer({
        if (!isTRUE(staged$published)) {
            try(discardProtDiaArtifactWorkerStage(
                staged$pending_stage$stage
            ), silent = TRUE)
        }
    })

    expect_true(staged$enabled)
    expect_true(staged$attempted)
    expect_true(staged$ok)
    expect_identical(staged$result, expected)
    expect_silent(artifactResourceDataOnly(staged$pending_stage))
    process_ids <- unlist(staged$process_evidence, use.names = FALSE)
    expect_length(unique(process_ids), 3L)
    expect_identical(staged$process_evidence$parent_pid, as.integer(Sys.getpid()))
    expect_identical(nrow(diaStageWorkerRegistryArtifacts(workflow)), 0L)

    workflow$data_tbl <- staged$result$data
    output <- persistProtDiaImportArtifacts(
        workflow_data = workflow,
        data_import_result = staged$result,
        source_path = source,
        pending_stage = staged$pending_stage,
        worker_attempted = TRUE
    )
    staged$published <- isTRUE(output$committed)
    expect_true(output$ok)
    expect_true(output$committed)
    artifacts <- diaStageWorkerRegistryArtifacts(workflow)
    expect_identical(nrow(artifacts), 1L)
    expect_identical(artifacts$artifact_id, output$refs$canonical_data$artifact_id)
})

test_that("memory mode never starts a DIA-NN artifact stage process", {
    root <- withr::local_tempdir()
    workflow <- diaStageWorkerWorkflow(root, backend = "memory")
    result <- stageProtDiaImportArtifacts(
        workflow_data = workflow,
        source_path = diaStageWorkerSource(),
        format = "diann"
    )

    expect_false(result$enabled)
    expect_false(result$attempted)
    expect_true(result$ok)
    expect_false(dir.exists(file.path(root, "state")))
})

test_that("worker hard exits leave no visible or orphaned generation", {
    diaStageWorkerSkipDependencies()
    root <- withr::local_tempdir()
    workflow <- diaStageWorkerWorkflow(root)
    result <- stageProtDiaImportArtifactsSafely(
        workflow_data = workflow,
        source_path = diaStageWorkerSource(),
        format = "diann",
        writer_failure_stage = "hard_exit_after_write",
        log_warn = \(...) invisible(NULL)
    )

    expect_true(result$attempted)
    expect_false(result$ok)
    expectDiaStageWorkerGenerationAbsent(diaStageWorkerStore(workflow))
    expect_identical(nrow(diaStageWorkerRegistryArtifacts(workflow)), 0L)
})

test_that("worker timeouts fall back without partial artifact state", {
    diaStageWorkerSkipDependencies()
    root <- withr::local_tempdir()
    workflow <- diaStageWorkerWorkflow(root)
    result <- stageProtDiaImportArtifactsSafely(
        workflow_data = workflow,
        source_path = diaStageWorkerSource(),
        format = "diann",
        timeout_ms = 1000L,
        writer_failure_stage = "timeout_before_import",
        log_warn = \(...) invisible(NULL)
    )

    expect_true(result$attempted)
    expect_false(result$ok)
    expectDiaStageWorkerGenerationAbsent(diaStageWorkerStore(workflow))
    expect_identical(nrow(diaStageWorkerRegistryArtifacts(workflow)), 0L)
})
