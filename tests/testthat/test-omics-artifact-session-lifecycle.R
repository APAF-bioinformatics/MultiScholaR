lifecycleSkipDependencies <- function(include_arrow = FALSE) {
    packages <- c("DBI", "duckdb", "filelock")
    if (isTRUE(include_arrow)) packages <- c(packages, "arrow")
    for (package in packages) testthat::skip_if_not_installed(package)
}

lifecycleRegistryFixture <- function(project_root = NULL) {
    if (is.null(project_root)) {
        project_root <- tempfile("artifact-lifecycle-registry-")
        dir.create(project_root, recursive = TRUE)
    }
    paths <- buildArtifactPaths(project_root, list(
        omic_type = "proteomics",
        omic_label = "proteomics_lifecycle",
        workflow_slug = "lifecycle"
    ))
    list(
        root = project_root,
        paths = paths,
        registry = newProjectRegistry(paths, "lifecycle-project")
    )
}

lifecycleArtifactCapability <- function() {
    matches <- Filter(\(capability) {
        identical(capability$identity$input_format, "diann") &&
            identical(capability$identity$data_level, "peptide")
    }, workflowCapabilityCatalogue())
    stopifnot(length(matches) == 1L)
    capability <- matches[[1L]]
    capability$artifact_eligible <- TRUE
    capability$auto_eligible <- TRUE
    capability$maximum_artifact_rollout <- "dual_write"
    capability
}

lifecycleArtifactContext <- function(project_root, project_id) {
    dir.create(project_root, recursive = TRUE, showWarnings = FALSE)
    context <- createWorkflowContext(
        experiment_paths = list(
            base_dir = project_root,
            omic_label = "proteomics_lifecycle"
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
        capabilities = list(lifecycleArtifactCapability())
    )
    context
}

lifecycleMemoryContext <- function(project_root) {
    dir.create(project_root, recursive = TRUE, showWarnings = FALSE)
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
    context
}

lifecyclePeptideObject <- function() {
    peptide_data <- data.frame(
        Protein.Group = rep(c("PG1", "PG2"), each = 2L),
        Protein.Ids = rep(c("P1", "P2"), each = 2L),
        Stripped.Sequence = rep(c("PEPTIDE_A", "PEPTIDE_B"), each = 2L),
        Run = rep(c("S1", "S2"), 2L),
        Q.Value = c(0.001, 0.002, 0.003, 0.004),
        Global.Q.Value = c(0.005, 0.006, 0.007, 0.008),
        Proteotypic = c(1L, 1L, 0L, 0L),
        Precursor.Quantity = c(100, 200, 300, 400),
        Precursor.Normalised = c(10, 20, 30, 40),
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
        args = list(branch = "artifact-lifecycle-test")
    )
    object@peptide_matrix <- matrix(
        c(10, 30, 20, 40),
        nrow = 2L,
        dimnames = list(c("PG1%PEPTIDE_A", "PG2%PEPTIDE_B"), c("S1", "S2"))
    )
    .ensurePeptideFeatureKeyMap(object, record_migration = FALSE)
}

lifecycleTempEntryCount <- function(registry) {
    path <- projectRegistryPath(registry, "temporary")
    if (!dir.exists(path)) return(0L)
    length(list.files(path, all.files = TRUE, no.. = TRUE))
}

test_that("artifact resource ownership is complete and data only", {
    ownership <- validateArtifactResourceOwnership()
    expect_setequal(
        ownership$resource_type,
        c(
            "registry_connection", "dbi_result", "writer_guard",
            "hydration_cache", "arrow_query_object", "duckdb_query_handle",
            "project_temporary_path", "background_task", "artifact_observer",
            "artifact_workflow_state"
        )
    )
    expect_true(all(ownership$process_bound[c(1L, 2L, 3L, 5L, 6L, 8L, 9L, 10L)]))
    expect_silent(artifactResourceDataOnly(ownership))
    expect_error(
        artifactResourceDataOnly(new.env(parent = emptyenv())),
        class = "multischolar_process_bound_worker_payload"
    )
    expect_error(
        artifactResourceDataOnly(ArtifactResourceScope$new("unsafe-worker")),
        class = "multischolar_process_bound_worker_payload"
    )
})

test_that("resource scopes unwind in reverse order and report every failure", {
    order <- character()
    scope <- ArtifactResourceScope$new("lifecycle-test")
    scope$register(
        "first",
        "artifact_observer",
        \() order <<- c(order, "first")
    )
    scope$register(
        "failing",
        "background_task",
        \() {
            order <<- c(order, "failing")
            rlang::abort("injected cleanup failure", class = "lifecycle_cleanup_test")
        }
    )
    scope$register(
        "last",
        "project_temporary_path",
        \() order <<- c(order, "last")
    )

    expect_error(
        scope$close("test_close"),
        class = "multischolar_artifact_resource_cleanup_failed"
    )
    expect_identical(order, c("last", "failing", "first"))
    info <- scope$getInfo()
    expect_true(info$closed)
    expect_identical(info$resource_count, 0L)
    expect_length(info$cleanup_errors, 1L)
    expect_false(scope$close("second_close"))
})

test_that("resource event evidence is bounded", {
    scope <- ArtifactResourceScope$new("bounded-events")
    for (index in seq_len(60L)) {
        resource_id <- paste0("observer-", index)
        scope$register(resource_id, "artifact_observer", \() invisible(TRUE))
        scope$release(resource_id)
    }
    expect_length(scope$getInfo()$events, .artifactResourceEventLimit)
    expect_true(scope$close())
})

test_that("background tasks are cancelled and joined without unsafe payloads", {
    events <- character()
    task <- new.env(parent = emptyenv())
    scope <- ArtifactResourceScope$new("task-owner")
    registerArtifactBackgroundTask(
        scope,
        "task-1",
        task,
        cancel_fn = \(task) events <<- c(events, "cancel"),
        join_fn = \(task) events <<- c(events, "join"),
        worker_payload = list(values = 1:3)
    )
    expect_true(scope$close("cancelled"))
    expect_identical(events, c("cancel", "join"))

    unsafe <- ArtifactResourceScope$new("unsafe-task")
    expect_error(
        registerArtifactBackgroundTask(
            unsafe,
            "task-2",
            task,
            cancel_fn = \(task) invisible(task),
            join_fn = \(task) invisible(task),
            worker_payload = task
        ),
        class = "multischolar_process_bound_worker_payload"
    )
    expect_identical(unsafe$getInfo()$resource_count, 0L)
    unsafe$close()
})

test_that("registry connections, results, writers, and temp paths close explicitly", {
    lifecycleSkipDependencies()
    fixture <- lifecycleRegistryFixture()
    session <- initializeProjectRegistry(fixture$registry)
    open <- projectRegistrySessionResourceInfo(session)
    expect_identical(open$connection_count, 1L)
    expect_identical(open$writer_guard_count, 1L)
    expect_identical(open$temporary_path_count, 1L)

    connection <- projectRegistrySessionConnection(session)
    result_events <- character()
    expect_identical(
        projectRegistryFetchBound(
            connection,
            "SELECT ? AS value",
            list(1L),
            result_observer = \(event) result_events <<- c(result_events, event)
        )$value,
        1L
    )
    expect_error(
        projectRegistryFetchBound(
            connection,
            "SELECT ? + ?",
            list(1L),
            result_observer = \(event) result_events <<- c(result_events, event)
        )
    )
    expect_identical(
        result_events,
        c("opened", "cleared", "opened", "cleared")
    )

    expect_true(closeProjectRegistry(session))
    closed <- projectRegistrySessionResourceInfo(session)
    expect_identical(closed$connection_count, 0L)
    expect_identical(closed$writer_guard_count, 0L)
    expect_identical(closed$temporary_path_count, 0L)
    expect_false(file.exists(projectRegistryPath(fixture$registry, "owner")))
    expect_identical(lifecycleTempEntryCount(fixture$registry), 0L)
    expect_false(closeProjectRegistry(session))
})

test_that("bounded query handles clean temp paths after query errors", {
    lifecycleSkipDependencies()
    fixture <- lifecycleRegistryFixture()
    store <- newArtifactStore(fixture$paths, "lifecycle-project")
    handle <- artifactQueryConnect(store, normalizeProjectRegistryPolicy())
    temporary_path <- handle$temporary_path
    expect_true(dir.exists(temporary_path))
    expect_error(
        projectRegistryFetchBound(
            handle$connection,
            "SELECT ? + ?",
            list(1L)
        )
    )
    expect_true(artifactQueryDisconnect(handle))
    expect_false(dir.exists(temporary_path))
})

test_that("failed registry startup releases writer and temporary resources", {
    lifecycleSkipDependencies()
    fixture <- lifecycleRegistryFixture()
    fail_startup <- \(stage, migration) {
        expect_identical(stage, "after_schema")
        expect_identical(migration$version, 1L)
        rlang::abort("injected startup failure", class = "lifecycle_startup_test")
    }
    expect_error(
        initializeProjectRegistry(fixture$registry, failure_injector = fail_startup),
        class = "lifecycle_startup_test"
    )
    expect_false(file.exists(projectRegistryPath(fixture$registry, "owner")))
    expect_identical(lifecycleTempEntryCount(fixture$registry), 0L)
})

test_that("artifact WorkflowState bounds cache and clears strong references", {
    lifecycleSkipDependencies(include_arrow = TRUE)
    root <- withr::local_tempdir()
    manager <- newWorkflowState(
        workflow_context = lifecycleArtifactContext(root, "lifecycle-cache")
    )
    withr::defer(try(manager$close(), silent = TRUE))
    manager$saveState(
        "imported",
        lifecyclePeptideObject(),
        list(stage = "import"),
        "Lifecycle cache"
    )
    first <- manager$getState()
    second <- manager$getState()
    expect_identical(first, second)
    stop_observing <- manager$observeTransitions(\(event) invisible(event))
    before <- manager$getResourceInfo()
    expect_identical(before$hydration_cache_entries, 1L)
    expect_identical(before$observer_count, 1L)
    expect_true(before$registry_connection)
    expect_true(before$writer_guard)

    expect_true(manager$close())
    after <- manager$getResourceInfo()
    expect_true(after$closed)
    expect_false(after$registry_connection)
    expect_false(after$writer_guard)
    expect_identical(after$hydration_cache_entries, 0L)
    expect_identical(after$observer_count, 0L)
    expect_false(manager$close())
    expect_silent(stop_observing())
})

test_that("cleanup errors release writer ownership without advancing state", {
    lifecycleSkipDependencies(include_arrow = TRUE)
    root <- withr::local_tempdir()
    context <- lifecycleArtifactContext(root, "lifecycle-cleanup-error")
    manager <- newWorkflowState(workflow_context = context)
    generation <- manager$getCurrentGenerationId()
    owner_path <- projectRegistryPath(
        projectRegistryForContext(context),
        "owner",
        must_exist = TRUE
    )
    testthat::local_mocked_bindings(
        artifactCleanupTemporaryPath = \(...) {
            rlang::abort(
                "injected temporary cleanup failure",
                class = "multischolar_artifact_temporary_cleanup_failed"
            )
        },
        .package = "MultiScholaR"
    )

    expect_error(
        manager$close(),
        class = "multischolar_artifact_temporary_cleanup_failed"
    )
    expect_identical(manager$getCurrentGenerationId(), generation)
    expect_true(manager$getResourceInfo()$closed)
    expect_false(file.exists(owner_path))
    expect_false(manager$close())
})

test_that("persistence errors preserve lineage and remain explicitly closable", {
    lifecycleSkipDependencies(include_arrow = TRUE)
    root <- withr::local_tempdir()
    context <- lifecycleArtifactContext(root, "lifecycle-persistence-error")
    manager <- newWorkflowState(workflow_context = context)
    withr::defer(try(manager$close(), silent = TRUE))
    generation <- manager$getCurrentGenerationId()
    payloads <- list.files(root, pattern = "[.]parquet$", recursive = TRUE)
    testthat::local_mocked_bindings(
        artifactStoreWriteParquet = \(...) {
            rlang::abort(
                "injected persistence failure",
                class = "lifecycle_persistence_test"
            )
        },
        .package = "MultiScholaR"
    )

    expect_error(
        manager$saveState(
            "failed",
            lifecyclePeptideObject(),
            list(stage = "import"),
            "Injected failure"
        ),
        class = "lifecycle_persistence_test"
    )
    expect_identical(manager$getCurrentGenerationId(), generation)
    expect_identical(
        list.files(root, pattern = "[.]parquet$", recursive = TRUE),
        payloads
    )
    expect_true(manager$getResourceInfo()$registry_connection)
    expect_true(manager$getResourceInfo()$writer_guard)
    expect_true(manager$close())
    expect_true(manager$getResourceInfo()$closed)
})

test_that("forked processes reject inherited artifact resources", {
    lifecycleSkipDependencies()
    skip_on_os("windows")
    fixture <- lifecycleRegistryFixture()
    session <- initializeProjectRegistry(fixture$registry)
    withr::defer(try(closeProjectRegistry(session), silent = TRUE))
    scope <- ArtifactResourceScope$new("fork-owner")
    withr::defer(try(scope$close(), silent = TRUE))
    child <- parallel::mcparallel({
        c(
            registry = tryCatch(
                projectRegistrySessionResourceInfo(session),
                error = \(error) class(error)[[1L]]
            ),
            scope = tryCatch(
                scope$getInfo(),
                error = \(error) class(error)[[1L]]
            )
        )
    })
    result <- parallel::mccollect(child)[[1L]]
    expect_identical(
        unname(result),
        rep("multischolar_cross_process_artifact_resource", 2L)
    )
    expect_true(projectRegistrySessionResourceInfo(session)$connection_count == 1L)
    expect_true(closeProjectRegistry(session))
    expect_true(scope$close())
})

test_that("session lifecycle closes replaced managers and module resources", {
    lifecycleSkipDependencies(include_arrow = TRUE)
    root <- withr::local_tempdir()
    first <- newWorkflowState(
        workflow_context = lifecycleArtifactContext(
            file.path(root, "first"),
            "lifecycle-first"
        )
    )
    second <- newWorkflowState(
        workflow_context = lifecycleArtifactContext(
            file.path(root, "second"),
            "lifecycle-second"
        )
    )
    session <- shiny::MockShinySession$new()
    workflow_data <- shiny::reactiveValues(
        workflow_context = lifecycleMemoryContext(file.path(root, "memory")),
        state_manager = first
    )
    scope <- shiny::withReactiveDomain(session, {
        shiny::isolate(registerWorkflowArtifactLifecycle(workflow_data, session))
    })
    session$flushReact()
    expect_true(scope$getInfo()$resource_count >= 2L)

    workflow_data$state_manager <- second
    session$flushReact()
    expect_true(first$getResourceInfo()$closed)
    expect_false(second$getResourceInfo()$closed)

    expect_true(closeWorkflowArtifactLifecycle(workflow_data, "module_removed"))
    expect_true(second$getResourceInfo()$closed)
    expect_true(scope$isClosed())
    expect_null(shiny::isolate(workflow_data$state_manager))
    session$close()
})

test_that("memory lifecycle creates no artifact-native resources", {
    root <- withr::local_tempdir()
    session <- shiny::MockShinySession$new()
    workflow_data <- shiny::reactiveValues(
        workflow_context = lifecycleMemoryContext(root),
        state_manager = newWorkflowState()
    )
    scope <- shiny::withReactiveDomain(session, {
        shiny::isolate(registerWorkflowArtifactLifecycle(workflow_data, session))
    })
    session$flushReact()
    info <- scope$getInfo()
    expect_false("artifact_workflow_state" %in% names(info$resources_by_type))
    expect_false(dir.exists(file.path(root, "state")))
    session$close()
    expect_true(scope$isClosed())
})

test_that("fresh-process lifecycle cycles plateau without explicit collection", {
    lifecycleSkipDependencies()
    testthat::skip_if_not_installed("processx")
    testthat::skip_if_not_installed("ps")
    output <- tempfile(fileext = ".json")
    log <- tempfile(fileext = ".log")
    script <- testthat::test_path(
        "..",
        "..",
        "tools",
        "profiling",
        "run_omics_artifact_lifecycle.R"
    )
    status <- system2(
        Sys.which("Rscript"),
        c("--vanilla", script, "--output", output, "--cycles", "4"),
        stdout = log,
        stderr = log
    )
    if (!identical(as.integer(status), 0L)) {
        testthat::fail(paste(readLines(log, warn = FALSE), collapse = "\n"))
    }
    result <- jsonlite::read_json(output, simplifyVector = FALSE)
    expect_false(result$used_explicit_garbage_collection)
    expect_setequal(
        vapply(result$runs, `[[`, character(1), "omic_type"),
        c("proteomics", "metabolomics", "lipidomics")
    )
    for (run in result$runs) {
        expect_true(run$plateau$passed, info = run$omic_type)
        for (cycle in run$cycles) {
            expect_identical(cycle$closed$connection_count, 0L)
            expect_identical(cycle$closed$writer_guard_count, 0L)
            expect_identical(cycle$closed$temporary_path_count, 0L)
            expect_identical(cycle$live_dbi_results_after_query, 0L)
            expect_identical(cycle$writer_guard_files_after_close, 0L)
            expect_identical(cycle$temporary_paths_after_close, 0L)
            expect_identical(cycle$background_tasks_after_close, 0L)
            expect_identical(cycle$artifact_observers_after_close, 0L)
            expect_identical(cycle$hydration_cache_entries_after_close, 0L)
            expect_identical(cycle$after$child_count, 0L)
        }
    }
})
