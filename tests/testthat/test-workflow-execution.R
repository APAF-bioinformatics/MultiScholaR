# Installed-only WFX-003 lifecycle suites. No implementation files are sourced.

makeExecutionTestOwner <- function(application, label) {
    root <- file.path(application$config$root, label)
    dir.create(root)
    store <- newWorkflowEntityStore(root, label, "DIA")
    manager <- newDIAStateAdapter(store, WorkflowState$new(audit_enabled = FALSE))
    manager$.context$reference_states <- TRUE
    object <- makeWorkflowCheckpointPeptide(readWorkflowCheckpointRecipe("small", workflow_test_root))
    manager$saveState("protein_peptide_filtered", object, list(), "execution fixture")
    rm(object)
    owner <- newWorkflowExecutionController(manager, application)
    owner
}

awaitExecutionTest <- function(owner, until_started = FALSE, seconds = 70) {
    deadline <- proc.time()[["elapsed"]] + seconds
    repeat {
        if (until_started && !is.null(owner$job) &&
            file.exists(file.path(owner$job$paths$request_dir, "worker-started.rds")))
            return(list(status = "started"))
        result <- pollWorkflowAction(owner)
        if (!result$status %in% c("admitting", "running")) return(result)
        if (proc.time()[["elapsed"]] > deadline) stop("bounded lifecycle wait exceeded")
        httpuv::service(0.01)
        Sys.sleep(0.02)
    }
}

test_that("qualification decisions match frozen cases and reject malformed bindings", {
    skip_if(!exists("workflow_test_mode") || !workflow_test_mode %in% c("protocol", "resources"), "installed lifecycle suite only")
    manifest <- jsonlite::read_json(file.path(workflow_test_root,
        "tests/testdata/workflow-execution/qualification.json"))
    keys <- c("operation", "input", "parameters", "runtime", "boundary", "policy", "outputs")
    bindings <- stats::setNames(as.list(rep(strrep("a", 64), length(keys))), keys)
    for (fixture in manifest$decision_fixtures$cases) {
        evidence <- fixture
        evidence$schema <- 1L
        evidence$bindings <- bindings
        evidence$open_holds <- as.character(unlist(fixture$open_holds, use.names = FALSE))
        expect_identical(evaluateWorkflowQualification(evidence, bindings)$status,
            fixture$expected, info = fixture$id)
    }
    evidence <- manifest$decision_fixtures$cases[[1L]]
    evidence$schema <- 1L
    evidence$bindings <- bindings
    evidence$open_holds <- character()
    changed <- evidence
    changed$bindings$runtime <- strrep("b", 64)
    expect_identical(evaluateWorkflowQualification(changed, bindings)$status, "held")
    changed <- evidence
    changed$science <- "true"
    expect_identical(evaluateWorkflowQualification(changed, bindings)$status, "held")
    changed <- evidence
    changed$pairs[[1]]$candidate$peak_bytes <- 17 * 1024^3
    expect_identical(evaluateWorkflowQualification(changed, bindings)$status, "rejected")
    expect_identical(evaluateWorkflowQualification(evidence, list())$status, "held")
})

test_that("installed controller commits once, preserves output-only retry and releases leases", {
    skip_if(!exists("workflow_test_mode") || !workflow_test_mode %in% c("protocol", "resources"), "installed lifecycle suite only")
    app <- workflow_test_application
    expect_identical(environment(runWorkflowAction), asNamespace("MultiScholaR"))
    owner <- makeExecutionTestOwner(app, "protocol")
    on.exit(closeDIAStateAdapter(owner$manager), add = TRUE)
    store <- owner$manager$.context$store
    other_store <- newWorkflowEntityStore(dirname(store$directory), "protocol", "DIA")
    expect_error(acquireWorkflowMutationLease(other_store), class = "workflow_execution_busy")
    foreign <- callr::r(function(root) {
        ns <- asNamespace("MultiScholaR")
        store <- get("newWorkflowEntityStore", ns)(root, "protocol", "DIA")
        tryCatch({get("acquireWorkflowMutationLease", ns)(store); FALSE},
            workflow_execution_busy = function(error) TRUE)
    }, args = list(dirname(store$directory)))
    expect_true(foreign)
    input <- list(state = getDIAStateReference(owner$manager))
    id <- runWorkflowAction("wfx_fixture", input, list(mode = "identity"), owner)
    paths <- owner$job$paths
    withr::local_options(multischolar.workflow_execution_boundary = function(stage) {
        if (stage == "publication") stop("intentional publication failure")
    })
    done <- awaitExecutionTest(owner)
    expect_identical(done$status, "publication_pending")
    accepted <- readRDS(file.path(paths$request_dir, "accepted-result.rds"))
    worker_hash <- digest::digest(file = file.path(paths$request_dir, "worker-started.rds"), algo = "sha256")
    options(multischolar.workflow_execution_boundary = NULL)
    resumed <- resumeWorkflowAction(owner, id)
    expect_identical(resumed$state_reference, accepted$state_reference)
    duplicate <- acceptWorkflowActionResult(owner, readRDS(file.path(paths$request_dir, "worker-result.rds")))
    expect_identical(duplicate$status, "duplicate")
    expect_true(duplicate$delivery_current)
    expect_identical(digest::digest(file = file.path(paths$request_dir, "worker-started.rds"), algo = "sha256"), worker_hash)
    expect_true(done$cleanup$worker_stopped && done$cleanup$monitor_stopped)
    expect_identical(store$open_resources, 0L)
    for (mode in c("fail", "forbidden_commit")) {
        before <- workflowExecutionHead(store)
        runWorkflowAction("wfx_fixture", list(state = getDIAStateReference(owner$manager)), list(mode = mode), owner)
        done <- awaitExecutionTest(owner)
        expect_identical(done$status, "worker_failed")
        expect_identical(workflowExecutionHead(store), before)
        expect_identical(store$open_resources, 0L)
    }
    runWorkflowAction("wfx_fixture", list(state = getDIAStateReference(owner$manager)),
        list(mode = "sleep", seconds = 20), owner)
    started <- awaitExecutionTest(owner, until_started = TRUE)
    expect_identical(started$status, "started")
    old <- owner$job$request$token
    previous_message <- readRDS(file.path(paths$request_dir, "worker-result.rds"))
    duplicate_while_busy <- acceptWorkflowActionResult(owner, previous_message)
    expect_identical(duplicate_while_busy$status, "duplicate")
    expect_false(duplicate_while_busy$delivery_current)
    expect_identical(owner$job$request$token, old)
    before <- workflowExecutionHead(store)
    cancelWorkflowAction(owner, "reverted")
    expect_identical(workflowExecutionHead(store), before)
    expect_gt(owner$control$generation, old$generation)
    expect_error(acceptWorkflowActionResult(owner, list(token = old)), class = "workflow_execution_stale")
    expect_null(app$active)
    expect_false(acceptWorkflowActionResult(owner, previous_message)$delivery_current)
    invisible(owner$manager$revertToState("protein_peptide_filtered"))
    expect_false(acceptWorkflowActionResult(owner, previous_message)$delivery_current)
    closeWorkflowEntityStore(other_store)
    cancelWorkflowAction(owner, "project_replaced")
    expect_false(acceptWorkflowActionResult(owner, previous_message)$delivery_current)
    closeWorkflowExecutionController(owner)
    expect_error(acceptWorkflowActionResult(owner, previous_message), class = "workflow_execution_stale")
    expect_error(resumeWorkflowAction(owner, id), class = "workflow_execution_stale")
    lock <- acquireWorkflowMutationLease(store)
    expect_false(is.null(lock))
    filelock::unlock(lock)
})

test_that("actual worker descendants are contained while the GUI controller stays usable", {
    skip_if(!exists("workflow_test_mode") || workflow_test_mode != "resources", "resource suite only")
    app <- workflow_test_application
    owner <- makeExecutionTestOwner(app, "resources")
    on.exit(closeDIAStateAdapter(owner$manager), add = TRUE)
    # A live HTTP event loop in the controller process remains outside W.
    port <- httpuv::randomPort()
    server <- httpuv::startServer("127.0.0.1", port, list(call = function(request)
        list(status = 200L, headers = list("Content-Type" = "text/plain"), body = "controller-alive")))
    on.exit(httpuv::stopServer(server), add = TRUE)
    shiny::testServer(function(input, output, session) {
        output$heartbeat <- shiny::renderText(input$ping)
    }, {
    session$setInputs(ping = "before")
    expect_identical(output$heartbeat, "before")
    before_pid <- Sys.getpid()
    before_head <- workflowExecutionHead(owner$manager$.context$store)
    id <- runWorkflowAction("wfx_fixture", list(state = getDIAStateReference(owner$manager)),
        list(mode = "allocate", bytes = 2 * 1024^3), owner)
    paths <- owner$job$paths
    worker_unit <- owner$job$worker_unit
    done <- awaitExecutionTest(owner)
    expect_true(file.exists(file.path(paths$request_dir, "worker-started.rds")))
    allocation <- readRDS(file.path(paths$request_dir, "allocation-started.rds"))
    started <- readRDS(file.path(paths$request_dir, "worker-started.rds"))
    expect_false(identical(allocation$pid, started$pid))
    expect_true(any(grepl(worker_unit, allocation$cgroup, fixed = TRUE)))
    session$setInputs(ping = "after")
    expect_identical(output$heartbeat, "after")
    expect_true(done$status %in% c("memory_guard", "nonworker_reservation"))
    expect_identical(Sys.getpid(), before_pid)
    expect_identical(workflowExecutionHead(owner$manager$.context$store), before_head)
    expect_lt(done$peak_bytes, app$config$policy$memory_max + 64 * 1024^2)
    expect_true(done$cleanup$worker_stopped && done$cleanup$monitor_stopped)
    client <- callr::r_bg(function(port) {
        readLines(paste0("http://127.0.0.1:", port), warn = FALSE)
    }, args = list(port))
    deadline <- proc.time()[["elapsed"]] + 10
    while (client$is_alive() && proc.time()[["elapsed"]] < deadline) httpuv::service(0.01)
    expect_false(client$is_alive())
    expect_identical(client$get_result(), "controller-alive")
    # A freed transient must remain visible, and a fresh window must reset below it.
    gc()
    first <- openWorkflowPeakWindow(app)
    values <- rep.int(1, 96 * 1024^2 / 8)
    rm(values)
    gc()
    current <- workflowMemoryCounter(file.path(app$cgroup, "memory.current"))
    high <- closeWorkflowPeakWindow(first)
    second <- openWorkflowPeakWindow(app)
    low <- closeWorkflowPeakWindow(second)
    valid_peak <- function(high, low, current) high > low + 32 * 1024^2 &&
        high > current + 32 * 1024^2
    expect_true(valid_peak(high, low, current))
    expect_false(valid_peak(high, high, current)) # lifetime-only reader
    expect_false(valid_peak(current, low, current)) # current-only reader
    later <- openWorkflowPeakWindow(app)
    values <- rep.int(1, 192 * 1024^2 / 8)
    rm(values)
    gc()
    later_peak <- closeWorkflowPeakWindow(later)
    expect_gt(later_peak, high + 32 * 1024^2)
    expect_lt(low, later_peak - 32 * 1024^2)
    saveRDS(list(high = high, low = low, current = current, later = later_peak),
        file.path(app$config$root, "peak-capability.rds"))
    expect_null(app$active)
    expect_identical(owner$manager$.context$store$open_resources, 0L)
    })
})

test_that("controller death on each side of head publication recovers one checkpoint", {
    skip_if(!exists("workflow_test_mode") || !workflow_test_mode %in% c("protocol", "resources"), "installed lifecycle suite only")
    app <- workflow_test_application
    for (stage in c("before_head_publication", "after_head_publication")) {
        label <- if (startsWith(stage, "before")) "crashbefore" else "crashafter"
        owner <- makeExecutionTestOwner(app, label)
        project <- dirname(owner$manager$.context$store$directory)
        closeDIAStateAdapter(owner$manager)
        config <- app$config
        config$root <- file.path(app$config$root, paste0(label, "-application"))
        dir.create(config$root)
        moved <- file.path(config$root, "project")
        expect_true(file.rename(project, moved))
        project <- moved
        config$stem <- paste0("multischolarwfx", substr(digest::digest(list(label, tempfile()), algo = "sha256"), 1, 16))
        config$slice <- paste0(config$stem, ".slice")
        config$control_group <- createWorkflowApplicationSlice(config$stem, config$policy)
        config$mode <- "fault_case"
        config$controller_unit <- paste0(config$stem, "-controller.service")
        config$fault_stage <- stage
        config$fault_project <- project
        config$fault_marker <- file.path(config$root, paste0(label, "-marker.rds"))
        path <- file.path(config$root, paste0(label, ".rds"))
        writeWorkflowExecutionMessage(path, config)
        workflowLaunchService(config, config$controller_unit, "application", path,
            digest::digest(file = path, algo = "sha256"), binds = app$config$controller_unit, deadline = 90L)
        deadline <- proc.time()[["elapsed"]] + 60
        while (!file.exists(config$fault_marker) && proc.time()[["elapsed"]] < deadline) Sys.sleep(0.05)
        expect_true(file.exists(config$fault_marker))
        marker <- readRDS(config$fault_marker)
        workflowSystemd(c("kill", "--signal=KILL", config$controller_unit))
        workflowSystemd(c("stop", config$controller_unit), required = FALSE)
        store <- newWorkflowEntityStore(project, label, "DIA")
        manager <- newDIAStateAdapter(store, reopen = TRUE)
        recovered <- newWorkflowExecutionController(manager, app)
        paths <- workflowExportPaths(store$directory, marker$request_id)
        journal <- readWorkflowExecutionJournal(paths)
        expect_identical(journal$status, "science_committed")
        accepted <- readRDS(file.path(paths$request_dir, "accepted-result.rds"))
        expect_identical(getDIAStateReference(manager), accepted$state_reference)
        started <- digest::digest(file = file.path(paths$request_dir, "worker-started.rds"), algo = "sha256")
        result <- resumeWorkflowAction(recovered, marker$request_id)
        expect_identical(result$state_reference, accepted$state_reference)
        expect_identical(digest::digest(file = file.path(paths$request_dir, "worker-started.rds"), algo = "sha256"), started)
        expect_identical(readWorkflowExecutionJournal(paths)$status, "complete")
        closeDIAStateAdapter(manager)
        workflowSystemd(c("stop", config$slice), required = FALSE)
    }
})

test_that("project replacement and monitor loss discard worker delivery and release resources", {
    skip_if(!exists("workflow_test_mode") || workflow_test_mode != "resources", "resource suite only")
    app <- workflow_test_application
    owner <- makeExecutionTestOwner(app, "replacement")
    on.exit(closeDIAStateAdapter(owner$manager), add = TRUE)
    runWorkflowAction("wfx_fixture", list(state = getDIAStateReference(owner$manager)),
        list(mode = "sleep", seconds = 20), owner)
    expect_identical(awaitExecutionTest(owner, until_started = TRUE)$status, "started")
    token <- owner$job$request$token
    old_head <- workflowExecutionHead(owner$manager$.context$store)
    cancelWorkflowAction(owner, "project_replaced")
    expect_error(acceptWorkflowActionResult(owner, list(token = token)), class = "workflow_execution_stale")
    expect_identical(workflowExecutionHead(owner$manager$.context$store), old_head)
    runWorkflowAction("wfx_fixture", list(state = getDIAStateReference(owner$manager)),
        list(mode = "sleep", seconds = 20), owner)
    expect_identical(awaitExecutionTest(owner, until_started = TRUE)$status, "started")
    worker <- owner$job$worker_unit
    workflowSystemd(c("kill", "--signal=KILL", owner$job$monitor_unit))
    result <- awaitExecutionTest(owner)
    expect_identical(result$status, "monitor_lost")
    expect_true(workflowServiceEmpty(app, worker))
    expect_null(app$active)
})

test_that("affected saved-entity regressions preserve existing boundaries", {
    skip_if(!exists("workflow_test_mode") || workflow_test_mode != "resources", "final installed suite only")
    library <- workflow_test_application$config$library
    result <- callr::r(function(library, root) {
        .libPaths(c(library, .libPaths()))
        ns <- asNamespace("MultiScholaR")
        stopifnot(identical(normalizePath(getNamespaceInfo(ns, "path")),
            normalizePath(file.path(library, "MultiScholaR"))))
        env <- new.env(parent = ns)
        sys.source(file.path(root, "tests/testthat/helper-workflow-checkpoint-fixtures.R"), env)
        names <- paste0("test-workflow-saved-", c("resources", "recovery", "identity", "projection", "exports"), ".R")
        results <- lapply(names, function(name) {
            results <- testthat::test_file(file.path(root, "tests/testthat", name), env = env,
                load_helpers = FALSE, package = "MultiScholaR", load_package = "installed",
                reporter = "summary", stop_on_failure = FALSE)
            results <- as.data.frame(results)
            list(file = name, passed = sum(results$passed), failed = sum(results$failed),
                errors = sum(results$error), details = results)
        })
        list(results = results)
    }, args = list(library, workflow_test_root), timeout = 360)
    saveRDS(result, file.path(workflow_test_application$config$root, "affected-regressions.rds"))
    for (item in result$results) {
        expect_identical(item$failed, 0L, info = item$file)
        expect_identical(item$errors, 0L, info = item$file)
    }
})

test_that("capability and request failures refuse before dispatch", {
    skip_if(!exists("workflow_test_mode"), "installed lifecycle suite only")
    app <- workflow_test_application
    fake <- WorkflowExecutionApplication$new()
    fake$config <- app$config
    fake$cgroup <- app$cgroup
    fake$pid <- Sys.getpid()
    fake$config$controller_unit <- paste0(app$config$stem, "-missing.service")
    expect_error(checkWorkflowApplication(fake), class = "workflow_execution_unsupported")
    fake$config <- app$config
    fake$config$policy$memory_max <- fake$config$policy$memory_max + 1
    expect_error(checkWorkflowApplication(fake), class = "workflow_execution_unsupported")
    expect_error(testthat::with_mocked_bindings(checkWorkflowApplication(app),
        workflowServiceProperties = function(unit) list(available = FALSE), .package = "MultiScholaR"),
        class = "workflow_execution_unsupported")
    path <- file.path(app$config$root, "wrong-request.rds")
    saveRDS(list(schema = 2L, build = workflowExecutionBuildIdentity()), path)
    expect_error(readWorkflowExecutionRequest(c(path, digest::digest(file = path, algo = "sha256"))),
        class = "workflow_execution_invalid")
    saveRDS(list(schema = 1L, build = strrep("0", 64)), path)
    expect_error(readWorkflowExecutionRequest(c(path, digest::digest(file = path, algo = "sha256"))),
        class = "workflow_execution_invalid")
    expect_null(app$active)
    expect_identical(workflowExecutionPolicy()$threads, 8L)
    expect_identical(workflowExecutionPolicy()$memory_max, 16 * 1024^3)
    expect_identical(workflowExecutionPolicy()$startup_headroom, 4 * 1024^3)
    expect_identical(workflowExecutionPolicy()$host_floor, 2 * 1024^3)
})


test_that("failed teardown retains ownership until a verified retry", {
    skip_if(!exists("workflow_test_mode") || workflow_test_mode != "resources", "resource suite only")
    app <- workflow_test_application
    owner <- makeExecutionTestOwner(app, "cleanup-repair")
    on.exit(closeDIAStateAdapter(owner$manager), add = TRUE)
    runWorkflowAction("wfx_fixture", list(state = getDIAStateReference(owner$manager)),
        list(mode = "sleep", seconds = 20), owner)
    expect_identical(awaitExecutionTest(owner, until_started = TRUE)$status, "started")
    id <- owner$control$active_request
    testthat::with_mocked_bindings({
        expect_error(closeWorkflowExecutionController(owner), class = "workflow_execution_cleanup_failed")
        expect_false(owner$closed)
        expect_false(is.null(owner$lease))
        expect_identical(owner$manager$.context$store$execution_owner, owner)
        expect_identical(owner$control$active_request, id)
        expect_identical(app$active, id)
        expect_error(runWorkflowAction("wfx_fixture",
            list(state = getDIAStateReference(owner$manager)), list(mode = "identity"), owner),
            class = "workflow_execution_busy")
    }, workflowServiceEmpty = function(...) FALSE, .package = "MultiScholaR")
    closeWorkflowExecutionController(owner)
    expect_true(owner$closed)
    expect_null(owner$job)
    expect_null(owner$lease)
    expect_null(app$active)
    expect_silent(closeWorkflowExecutionController(owner))
})

test_that("replacement project paths are untouched during teardown", {
    skip_if(!exists("workflow_test_mode") || workflow_test_mode != "resources", "resource suite only")
    app <- workflow_test_application
    owner <- makeExecutionTestOwner(app, "replacement-repair")
    on.exit(closeDIAStateAdapter(owner$manager), add = TRUE)
    runWorkflowAction("wfx_fixture", list(state = getDIAStateReference(owner$manager)),
        list(mode = "sleep", seconds = 20), owner)
    expect_identical(awaitExecutionTest(owner, until_started = TRUE)$status, "started")
    job <- owner$job
    store <- owner$manager$.context$store
    original <- store$directory
    moved <- paste0(original, "-original")
    expect_true(file.rename(original, moved))
    replacement <- newWorkflowEntityStore(dirname(original), "replacement-other", "DIA")
    on.exit(closeWorkflowEntityStore(replacement), add = TRUE)
    dir.create(job$paths$request_dir, recursive = TRUE)
    for (role in c("worker", "monitor")) {
        path <- file.path(job$paths$request_dir, paste0(role, "-scratch"))
        dir.create(path)
        writeLines("replacement-owned", file.path(path, "keep"))
    }
    hashes <- function() {
        files <- sort(list.files(original, recursive = TRUE, full.names = TRUE, all.files = TRUE))
        files <- files[!dir.exists(files)]
        stats::setNames(vapply(files, function(path) digest::digest(file = path, algo = "sha256"),
            character(1)), files)
    }
    before <- hashes()
    done <- pollWorkflowAction(owner)
    expect_identical(done$status, "project_replaced")
    expect_identical(hashes(), before)
    expect_true(done$cleanup$worker_stopped && done$cleanup$monitor_stopped)
    expect_true(file.exists(file.path(app$config$root, paste0(job$request$token$request_id, "-completion.rds"))))
    expect_false(done$measurement_valid)
    expect_null(app$active)
})

test_that("monitor launch refusals preserve receipts and allow safe readmission", {
    skip_if(!exists("workflow_test_mode") || workflow_test_mode != "resources", "resource suite only")
    app <- workflow_test_application
    owner <- makeExecutionTestOwner(app, "launch-repair")
    on.exit(closeDIAStateAdapter(owner$manager), add = TRUE)
    launch <- workflowLaunchService
    for (created in c(FALSE, TRUE)) {
        error <- testthat::with_mocked_bindings(tryCatch(
            runWorkflowAction("wfx_fixture", list(state = getDIAStateReference(owner$manager)),
                list(mode = "identity"), owner), workflow_execution_refused = identity),
            workflowLaunchService = function(...) {
                if (created) launch(...)
                cli::cli_abort("injected monitor launch acknowledgement failure",
                    class = "workflow_execution_unsupported")
            }, .package = "MultiScholaR")
        expect_s3_class(error, "workflow_execution_refused")
        paths <- workflowExportPaths(owner$manager$.context$store$directory, error$request_id)
        receipt <- readRDS(file.path(paths$request_dir, "completion.rds"))
        expect_identical(receipt$status, "admission_refused")
        expect_true(receipt$cleanup$worker_stopped && receipt$cleanup$monitor_stopped)
        expect_identical(readWorkflowExecutionJournal(paths)$status, "admission_refused")
        expect_null(owner$control$active_request)
        expect_null(owner$job)
        expect_null(app$active)
        expect_false(dir.exists(file.path(paths$request_dir, "monitor-scratch")))
    }
    runWorkflowAction("wfx_fixture", list(state = getDIAStateReference(owner$manager)), list(mode = "identity"), owner)
    expect_identical(awaitExecutionTest(owner)$status, "complete")
})
