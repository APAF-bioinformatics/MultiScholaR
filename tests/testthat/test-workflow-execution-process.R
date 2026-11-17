# KIT-131: the managed worker without systemd (the plain `processx` launcher), forced here by its configuration
# (`launcher = "process"`) and by the option `multischolar.workflow_launcher`. Layer 1 drives the process side of the
# unit functions with trivial R scripts and needs no project. Layer 2 drives the real protocol with the diagnostic
# fixture through a plain application whose controller is this process.
#
# A worker is `Rscript` loading the package. Against the installed package the launcher's own command line runs. From
# source (`load_all`) no installed build exists for a child to load, so `processTestRuntime()` replaces the entry
# expression of each child with a script that loads this source tree and then calls the same entry point, and fixes
# the build identity, as the other source tests do; the command line the launcher built is asserted before it is
# replaced. Each test starts at most one worker at a time.

processTestRscript <- function() {
    file.path(R.home("bin"), if (identical(.Platform$OS.type, "windows")) "Rscript.exe" else "Rscript")
}

processTestLibraries <- function() c(R_LIBS = paste(.libPaths(), collapse = .Platform$path.sep))

#' The source tree when the package is loaded from source, `"-"` when it is installed (an empty argument does not
#' survive every command line).
processTestSourceRoot <- function() {
    if (nzchar(system.file("R", "MultiScholaR.rdb", package = "MultiScholaR"))) return("-")
    normalizePath(testthat::test_path("..", ".."), winslash = "/")
}

processTestUnit <- function(label) paste0("multischolarwfxtest", Sys.getpid(), "-", label, ".service")

#' Start one unit running an R script, stopped and forgotten when the calling test ends.
processTestStart <- function(root, label, lines, arguments = character(), deadline = 60L, env = parent.frame()) {
    script <- file.path(root, paste0(label, ".R"))
    writeLines(lines, script)
    unit <- processTestUnit(label)
    withr::defer({
        workflowProcessStop(unit)
        workflowProcessForget(unit)
    }, envir = env)
    workflowProcessStart(unit, processTestRscript(), c("--vanilla", script, arguments), processTestLibraries(),
        root, file.path(root, paste0(label, ".log")), deadline)
    unit
}

processTestAwait <- function(config, unit, seconds = 60) {
    deadline <- proc.time()[["elapsed"]] + seconds
    repeat {
        state <- workflowUnitState(config, unit)
        if (!isTRUE(state$available) || !state$SubState %in% c("running", "start")) return(state)
        if (proc.time()[["elapsed"]] > deadline) stop("bounded unit wait exceeded")
        Sys.sleep(0.05)
    }
}

processTestFile <- function(path, seconds = 180) {
    deadline <- proc.time()[["elapsed"]] + seconds
    while (!file.exists(path) && proc.time()[["elapsed"]] < deadline) Sys.sleep(0.05)
    file.exists(path)
}

processTestAlive <- function(pid) {
    tryCatch({
        handle <- ps::ps_handle(as.integer(pid))
        ps::ps_is_running(handle) && !identical(ps::ps_status(handle), "zombie")
    }, error = function(error) FALSE)
}

processTestGone <- function(pids, seconds = 5) {
    deadline <- proc.time()[["elapsed"]] + seconds
    repeat {
        alive <- vapply(pids, processTestAlive, logical(1))
        if (!any(alive) || proc.time()[["elapsed"]] > deadline) return(!any(alive))
        Sys.sleep(0.05)
    }
}

processTestApplication <- function(root) {
    application <- WorkflowExecutionApplication$new()
    application$config <- list(launcher = "process", root = root, mode = "gui",
        library = dirname(system.file(package = "MultiScholaR")))
    application$pid <- Sys.getpid()
    application
}

# --- Layer 1: the process side of the unit functions ---

test_that("L1.1 a plain unit reports its exit status, and a command that cannot start is refused", {
    root <- withr::local_tempdir()
    config <- processTestApplication(root)$config
    succeeded <- processTestStart(root, "succeeds", "quit(status = 0L)")
    running <- workflowUnitState(config, succeeded)
    expect_true(running$available)
    state <- processTestAwait(config, succeeded)
    expect_identical(state$SubState, "exited")
    expect_identical(state$Result, "success")
    expect_identical(state$MainPID, running$MainPID)
    # The handle keeps the status until the unit is forgotten: a second read gives the same answer.
    expect_identical(workflowUnitState(config, succeeded)$Result, "success")
    failed <- processTestStart(root, "fails", "stop('intentional unit failure')")
    state <- processTestAwait(config, failed)
    expect_identical(state$SubState, "failed")
    expect_identical(state$Result, "exit-code")
    missing <- processTestUnit("missing")
    expect_error(workflowProcessStart(missing, file.path(root, "no-such-command"), character(), character(), root,
        file.path(root, "missing.log"), 60L), class = "workflow_execution_unsupported")
    expect_false(workflowUnitState(config, missing)$available)
    # One name is one unit.
    expect_error(workflowProcessStart(succeeded, processTestRscript(), c("--vanilla", "-e", "1"), character(), root,
        file.path(root, "again.log"), 60L), class = "workflow_execution_invalid")
})

test_that("L1.2 stopping a unit ends everything it started, and the unit is then empty", {
    root <- withr::local_tempdir()
    application <- processTestApplication(root)
    record <- file.path(root, "pids.rds")
    unit <- processTestStart(root, "tree", c(
        "record <- commandArgs(TRUE)[[1L]]",
        "rscript <- file.path(R.home('bin'), if (identical(.Platform$OS.type, 'windows')) 'Rscript.exe' else 'Rscript')",
        # The grandchild leaves the unit's process group (processx starts it in its own session).
        "child <- processx::process$new(rscript, c('--vanilla', '-e', 'Sys.sleep(120)'), cleanup = FALSE)",
        "saveRDS(list(unit = Sys.getpid(), child = child$get_pid()), paste0(record, '.tmp'))",
        "file.rename(paste0(record, '.tmp'), record)",
        "Sys.sleep(120)"), arguments = record)
    expect_true(processTestFile(record))
    pids <- readRDS(record)
    launched <- workflowUnitState(application$config, unit)$MainPID
    expect_true(processTestAlive(pids$unit) && processTestAlive(pids$child))
    expect_false(workflowServiceEmpty(application, unit))
    # Whether the R session is the launched process or its child is the platform's business (on windows-latest with
    # R 4.6.1, run 1, it was the launched `Rscript.exe` itself): the stop below must end both either way.
    started <- proc.time()[["elapsed"]]
    workflowUnitStop(application$config, unit)
    expect_lt(proc.time()[["elapsed"]] - started, 3)
    expect_true(processTestGone(c(pids$unit, pids$child, launched)))
    expect_true(workflowServiceEmpty(application, unit))
    state <- workflowUnitState(application$config, unit)
    expect_identical(state$SubState, "failed")
    expect_false(identical(state$Result, "success"))
    expect_silent(workflowUnitStop(application$config, unit))
})

test_that("L1.3 a unit that outlives its deadline is stopped by the state read and reported as a timeout", {
    root <- withr::local_tempdir()
    application <- processTestApplication(root)
    record <- file.path(root, "pid.rds")
    started <- proc.time()[["elapsed"]]
    unit <- processTestStart(root, "deadline", c("saveRDS(Sys.getpid(), commandArgs(TRUE)[[1L]])", "Sys.sleep(30)"),
        arguments = record, deadline = 2L)
    expect_identical(workflowUnitState(application$config, unit)$SubState, "running")
    state <- processTestAwait(application$config, unit)
    elapsed <- proc.time()[["elapsed"]] - started
    expect_identical(state$Result, "timeout")
    expect_identical(state$SubState, "failed")
    expect_gte(elapsed, 2)
    expect_lt(elapsed, 6)
    expect_true(processTestGone(c(readRDS(record), state$MainPID)))
    expect_true(workflowServiceEmpty(application, unit))
    expect_identical(workflowUnitState(application$config, unit)$Result, "timeout")
})

test_that("L1.4 a unit does not outlive a controller that is killed", {
    root <- withr::local_tempdir()
    sleeper <- file.path(root, "sleeper.R")
    writeLines("Sys.sleep(120)", sleeper)
    record <- file.path(root, "started.rds")
    script <- file.path(root, "controller.R")
    writeLines(c(
        "arguments <- commandArgs(TRUE)",
        "if (!identical(arguments[[1L]], '-')) suppressMessages(pkgload::load_all(arguments[[1L]], quiet = TRUE, helpers = FALSE, attach_testthat = FALSE))",
        "ns <- asNamespace('MultiScholaR')",
        "rscript <- file.path(R.home('bin'), if (identical(.Platform$OS.type, 'windows')) 'Rscript.exe' else 'Rscript')",
        "ns$workflowProcessStart('lifetime.service', rscript, c('--vanilla', arguments[[3L]]), character(), arguments[[2L]],",
        "    file.path(arguments[[2L]], 'unit.log'), 120L)",
        "saveRDS(list(controller = Sys.getpid(), unit = ns$workflowProcessState('lifetime.service')$MainPID), paste0(arguments[[4L]], '.tmp'))",
        "file.rename(paste0(arguments[[4L]], '.tmp'), arguments[[4L]])",
        "Sys.sleep(120)"), script)
    controller <- processx::process$new(processTestRscript(), c("--vanilla", script, processTestSourceRoot(), root, sleeper, record),
        env = c("current", processTestLibraries()), stdout = file.path(root, "controller.log"), stderr = "2>&1", cleanup = TRUE)
    withr::defer(controller$kill())
    expect_true(processTestFile(record))
    pids <- readRDS(record)
    withr::defer(if (processTestAlive(pids$unit)) ps::ps_kill(ps::ps_handle(as.integer(pids$unit))))
    expect_true(processTestAlive(pids$unit))
    # A hard kill of the controller alone: no exit handler runs, and the unit is in its own process group.
    ps::ps_kill(ps::ps_handle(as.integer(pids$controller)))
    expect_true(processTestGone(pids$controller))
    expect_true(processTestGone(pids$unit))
})

test_that("L1.5 a worker unit gets the contract's environment, its scratch directory as working directory, and its own log", {
    root <- withr::local_tempdir()
    application <- processTestApplication(root)
    request_dir <- file.path(root, "request")
    dir.create(request_dir)
    path <- file.path(request_dir, "request.rds")
    saveRDS(list(schema = 1L), path)
    hash <- digest::digest(file = path, algo = "sha256")
    record <- file.path(root, "seen.rds")
    script <- file.path(root, "report.R")
    writeLines(c(
        "message('unit standard error line')",
        "names <- c('MULTISCHOLAR_WORKFLOW_WORKER', 'TMPDIR', 'OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS',",
        "    'DUCKPLYR_FALLBACK_AUTOUPLOAD', 'DUCKPLYR_FALLBACK_COLLECT', 'R_LIBS')",
        "saveRDS(list(environment = Sys.getenv(names), directory = getwd(), arguments = commandArgs(TRUE)), commandArgs(TRUE)[[1L]])"),
        script)
    unit <- processTestUnit("worker")
    withr::defer({
        workflowProcessStop(unit)
        workflowProcessForget(unit)
    })
    start <- workflowProcessStart
    built <- NULL
    testthat::local_mocked_bindings(workflowProcessStart = function(unit, command, args, ...) {
        built <<- list(command = command, args = args)
        start(unit, command, c("--vanilla", script, record), ...)
    }, .package = "MultiScholaR")
    before <- names(Sys.getenv())
    set.seed(131)
    expected_draw <- stats::runif(1)
    set.seed(131)
    workflowLaunchService(application$config, unit, "worker", path, hash, cap = 2 * 1024^3,
        binds = "ignored.service", deadline = 60L)
    # Starting a unit leaves the controller's random stream and environment as they were.
    expect_identical(stats::runif(1), expected_draw)
    expect_identical(names(Sys.getenv()), before)
    expect_identical(basename(built$command), basename(processTestRscript()))
    expect_identical(built$args, c("--vanilla", "-e", "MultiScholaR:::workflowWorkerMain(commandArgs(TRUE))", path, hash))
    expect_identical(processTestAwait(application$config, unit)$Result, "success")
    seen <- readRDS(record)
    scratch <- file.path(request_dir, "worker-scratch")
    expect_identical(normalizePath(seen$directory, winslash = "/"), normalizePath(scratch, winslash = "/"))
    environment <- as.list(seen$environment)
    expect_identical(environment$MULTISCHOLAR_WORKFLOW_WORKER, "1")
    expect_identical(normalizePath(environment$TMPDIR, winslash = "/"), normalizePath(scratch, winslash = "/"))
    for (name in c("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS")) expect_identical(environment[[name]], "1")
    for (name in c("DUCKPLYR_FALLBACK_AUTOUPLOAD", "DUCKPLYR_FALLBACK_COLLECT")) expect_identical(environment[[name]], "0")
    expect_identical(environment$R_LIBS, paste(unique(c(application$config$library, .libPaths())), collapse = .Platform$path.sep))
    log <- file.path(root, paste0("worker-", unit, ".log"))
    expect_true(any(grepl("unit standard error line", readLines(log, warn = FALSE), fixed = TRUE)))
    # The plain launcher starts workers only: the application runs in the launching process and there is no monitor.
    for (role in c("application", "monitor")) expect_error(workflowLaunchService(application$config,
        processTestUnit(role), role, path, hash), class = "workflow_execution_unsupported")
})

test_that("L1.6 a unit that was never started is unavailable, empty, and silent to stop", {
    application <- processTestApplication(withr::local_tempdir())
    unit <- processTestUnit("never")
    expect_identical(workflowUnitState(application$config, unit), list(available = FALSE))
    expect_silent(workflowUnitStop(application$config, unit))
    expect_true(workflowServiceEmpty(application, unit))
})

test_that("L1.7 the option chooses the launcher; without it only a system with systemctl gets systemd", {
    expect_identical(workflowLauncherChoice("process", "Linux", "/usr/bin/systemctl"), "process")
    expect_identical(workflowLauncherChoice("systemd", "Windows", ""), "systemd")
    expect_identical(workflowLauncherChoice(NULL, "Linux", "/usr/bin/systemctl"), "systemd")
    expect_identical(workflowLauncherChoice(NULL, "Linux", ""), "process")
    expect_identical(workflowLauncherChoice(NULL, "Windows", ""), "process")
    expect_identical(workflowLauncherChoice(NULL, "Darwin", ""), "process")
    for (invalid in list("plain", c("process", "systemd"), NA_character_, TRUE)) expect_error(
        workflowLauncherChoice(invalid, "Linux", "/usr/bin/systemctl"), class = "workflow_execution_invalid")
    # The configuration carries the choice: an absent launcher is systemd, as every stored configuration is today.
    expect_false(workflowPlainLauncher(list(mode = "gui")))
    expect_false(workflowPlainLauncher(NULL))
    expect_true(workflowPlainLauncher(list(launcher = "process")))
})

# --- D4: the state owner's paths where normalizePath() answers with backslashes (Windows) ---

test_that("D4 the state owner's containment checks hold where normalizePath() defaults to backslashes", {
    root <- withr::local_tempdir()
    real <- normalizePath
    # Windows' normalizePath() joins with "\\" unless asked for "/"; the store's checks compare with "/"-joined paths.
    windows <- function(path, winslash = "\\", mustWork = NA) gsub("/", winslash, real(path, "/", mustWork), fixed = TRUE)
    store <- testthat::with_mocked_bindings(newWorkflowEntityStore(root, "d4", "DIA"), normalizePath = windows,
        .package = "base")
    withr::defer(closeWorkflowEntityStore(store))
    payload <- file.path(store$directory, "entities", "d4", "payload.rds")
    dir.create(dirname(payload), recursive = TRUE)
    saveRDS(1L, payload)
    testthat::with_mocked_bindings({
        expect_true(checkWorkflowEntityStore(store))
        expect_identical(checkSavedEntityPath(store, "entities/d4/payload.rds"), "entities/d4/payload.rds")
        expect_true(syncSavedWorkflowPaths(store, c(payload, store$directory)))
        expect_error(syncSavedWorkflowPaths(store, root), NA)
        expect_error(syncSavedWorkflowPaths(store, tempdir()), class = "saved_entity_path")
        expect_identical(hashWorkflowOwnedFile(store, payload), digest::digest(file = payload, algo = "sha256"))
        paths <- workflowExportPaths(root, "d4-request")
        expect_true(startsWith(paths$request_dir, paste0(paths$root, "/")))
        expect_true(syncWorkflowExportPaths(paths$request_dir, paths$root))
    }, normalizePath = windows, .package = "base")
})

test_that("D4 reading a staged table's metadata leaves no handle or mapping on the file (Windows cannot rename around one)", {
    # Windows refuses to rename a directory while a file inside it is open or mapped: the staged table's commit
    # (run 1 on windows-latest: "cannot rename file ... reason 'Access is denied'"). Linux shows the same leak in /proc.
    skip_if_not(dir.exists("/proc/self/fd") && file.exists("/proc/self/maps"))
    root <- withr::local_tempdir()
    path <- file.path(root, "payload.parquet")
    writeSavedWorkflowParquet(data.frame(a = 1:3, b = c("x", "y", "z")), ".rownames", path, metadata_version = 3L)
    target <- normalizePath(path)
    held <- function() {
        descriptors <- list.files("/proc/self/fd", full.names = TRUE)
        # A descriptor can close between the listing and its read (the listing's own): its link reads as NA.
        any(grepl(target, readLines("/proc/self/maps", warn = FALSE), fixed = TRUE)) ||
            target %in% vapply(descriptors, Sys.readlink, character(1))
    }
    invisible(gc())
    expect_false(held())
    expect_identical(readSavedWorkflowTableMetadata(path)$rows, 3L)
    expect_false(held())
})

# --- Layer 2: the protocol under the plain launcher, with the diagnostic fixture ---

#' From source, make every child of this test load this source tree (see the head of the file). Against the installed
#' package nothing is replaced.
#' @return The script a child controller sources first, or `"-"` against the installed package.
processTestRuntime <- function(root, env = parent.frame()) {
    source_root <- processTestSourceRoot()
    if (identical(source_root, "-")) return("-")
    testthat::skip_if_not_installed("pkgload")
    preamble <- file.path(root, "source-runtime.R")
    worker <- file.path(root, "source-worker.R")
    entry <- c("--vanilla", "-e", "MultiScholaR:::workflowWorkerMain(commandArgs(TRUE))")
    writeLines(c(
        sprintf("suppressMessages(pkgload::load_all(%s, quiet = TRUE, helpers = FALSE, attach_testthat = FALSE))", deparse(source_root)),
        "utils::assignInNamespace('workflowExecutionBuildIdentity', function() strrep('f', 64), 'MultiScholaR')",
        "local({",
        "    start <- get('workflowProcessStart', asNamespace('MultiScholaR'))",
        sprintf("    worker <- %s", deparse(worker)),
        "    utils::assignInNamespace('workflowProcessStart', function(unit, command, args, ...)",
        "        start(unit, command, c('--vanilla', worker, args[-(1:3)]), ...), 'MultiScholaR')",
        "})"), preamble)
    writeLines(c(sprintf("source(%s)", deparse(preamble)), entry[[3L]]), worker)
    start <- workflowProcessStart
    testthat::local_mocked_bindings(
        workflowExecutionBuildIdentity = function() strrep("f", 64),
        workflowProcessStart = function(unit, command, args, ...) {
            stopifnot(identical(args[1:3], entry))
            start(unit, command, c("--vanilla", worker, args[-(1:3)]), ...)
        }, .package = "MultiScholaR", .env = env)
    preamble
}

#' A plain managed application whose controller is this process, with the diagnostic policy the fixture needs.
processTestManagedApplication <- function(root) {
    directory <- tempfile("application-", tmpdir = root)
    dir.create(directory)
    application <- WorkflowExecutionApplication$new()
    application$config <- workflowProcessApplicationConfig(directory, diagnostic = TRUE, port = 0L)
    application$pid <- Sys.getpid()
    application
}

processTestProject <- function(root, label) {
    project <- file.path(root, label)
    dir.create(project)
    store <- newWorkflowEntityStore(project, label, "DIA")
    manager <- newDIAStateAdapter(store, WorkflowState$new(audit_enabled = FALSE))
    manager$.context$reference_states <- TRUE
    object <- makeWorkflowCheckpointPeptide(readWorkflowCheckpointRecipe("small", testthat::test_path("..", "..")))
    manager$saveState("protein_peptide_filtered", object, list(), "execution fixture")
    manager
}

processTestAction <- function(owner, until_started = FALSE, seconds = 180) {
    deadline <- proc.time()[["elapsed"]] + seconds
    repeat {
        if (until_started && !is.null(owner$job) &&
            file.exists(file.path(owner$job$paths$request_dir, "worker-started.rds"))) return(list(status = "started"))
        result <- pollWorkflowAction(owner)
        if (!result$status %in% c("admitting", "running")) return(result)
        if (proc.time()[["elapsed"]] > deadline) stop("bounded lifecycle wait exceeded")
        Sys.sleep(0.05)
    }
}

processTestFixture <- function(owner, parameters, reservation = NULL) {
    runWorkflowAction("wfx_fixture", list(state = getDIAStateReference(owner$manager)), parameters, owner,
        reservation = reservation)
}

test_that("L2.1 a GUI action runs on a plain worker: one commit, the result read back, a receipt without measurements", {
    skip_on_cran()
    root <- withr::local_tempdir()
    processTestRuntime(root)
    application <- processTestManagedApplication(root)
    expect_true(checkWorkflowApplication(application))
    manager <- processTestProject(root, "identity")
    withr::defer(closeDIAStateAdapter(manager))
    owner <- newWorkflowExecutionController(manager, application)
    withr::defer(closeWorkflowExecutionController(owner))
    store <- manager$.context$store
    before <- workflowExecutionHead(store)
    states <- length(manager$getHistory())
    reservation <- beginWorkflowGUIAction(application, manager, "wfx_fixture", list(experiment_label = "test",
        experiment_paths = list(source_dir = root), omic_type = "proteomics"))
    expect_null(reservation$peak)
    id <- processTestFixture(owner, list(mode = "identity"), reservation)
    paths <- owner$job$paths
    unit <- owner$job$worker_unit
    done <- processTestAction(owner)
    expect_identical(done$status, "complete")
    # One commit: the head moved once, by one state.
    expect_false(identical(workflowExecutionHead(store)$entity_id, before$entity_id))
    expect_identical(length(manager$getHistory()), states + 1L)
    # The result was computed in another process and read back from its file.
    started <- readRDS(file.path(paths$request_dir, "worker-started.rds"))
    expect_false(identical(started$pid, Sys.getpid()))
    handoff <- readRDS(file.path(paths$request_dir, "worker-result.rds"))
    expect_identical(handoff$status, "complete")
    expect_identical(handoff$token$request_id, id)
    accepted <- resumeWorkflowAction(owner, id)
    expect_identical(accepted$state_reference, getDIAStateReference(manager))
    expect_identical(readWorkflowExecutionJournal(paths)$status, "complete")
    # L2.7: the admitted cap is a number handed to the step; nothing enforces it.
    admitted <- readRDS(file.path(paths$request_dir, "admitted.rds"))
    expect_true(is.finite(admitted$worker_cap) && admitted$worker_cap > 0)
    expect_identical(readRDS(file.path(paths$request_dir, "request.rds"))$worker_cap, admitted$worker_cap)
    # The receipt: the worker stopped, and no measurement is claimed.
    expect_identical(done$worker$Result, "success")
    expect_true(done$cleanup$worker_stopped)
    expect_true(is.na(done$cleanup$application_monitor_active))
    expect_false(done$measurement_valid)
    expect_true(is.na(done$peak_bytes) && is.na(done$whole_service_peak_bytes))
    expect_true(all(is.na(done$accounting_delta)))
    expect_null(done$monitor)
    expect_identical(readRDS(file.path(paths$request_dir, "completion.rds"))$status, "complete")
    expect_false(dir.exists(file.path(paths$request_dir, "worker-scratch")))
    expect_false(workflowUnitState(application$config, unit)$available)
    expect_null(owner$job)
    # The GUI action closes on the browser's acknowledgement, without waiting for a monitor.
    closing <- proc.time()[["elapsed"]]
    receipt <- finishWorkflowGUIAction(reservation, "complete", list(request_id = reservation$id,
        content_id = "content", plot_ready = TRUE, stable_seconds = 0.2))
    expect_lt(proc.time()[["elapsed"]] - closing, 1.5)
    expect_identical(receipt$status, "complete")
    expect_identical(receipt$worker_request_id, id)
    expect_false(receipt$measurement_valid)
    expect_true(is.na(receipt$peak_bytes) && is.na(receipt$whole_service_peak_bytes))
    expect_null(application$active)
    expect_identical(store$open_resources, 0L)
})

test_that("L2.2 a failing worker and a worker that tries to commit end as worker_failed, and the head stays", {
    skip_on_cran()
    root <- withr::local_tempdir()
    processTestRuntime(root)
    application <- processTestManagedApplication(root)
    manager <- processTestProject(root, "failure")
    withr::defer(closeDIAStateAdapter(manager))
    owner <- newWorkflowExecutionController(manager, application)
    withr::defer(closeWorkflowExecutionController(owner))
    store <- manager$.context$store
    before <- workflowExecutionHead(store)
    expected <- c(fail = "intentional bounded worker failure", forbidden_commit = "workers cannot mutate workflow heads")
    for (mode in names(expected)) {
        processTestFixture(owner, list(mode = mode))
        done <- processTestAction(owner)
        expect_identical(done$status, "worker_failed", info = mode)
        expect_match(done$error, expected[[mode]], fixed = TRUE, info = mode)
        expect_identical(done$worker$Result, "exit-code", info = mode)
        expect_true(done$cleanup$worker_stopped, info = mode)
        expect_identical(workflowExecutionHead(store), before, info = mode)
        expect_null(application$active)
        expect_identical(store$open_resources, 0L)
    }
})

test_that("L2.3 a cancellation stops the worker and leaves the head, the history and the action lock as they were", {
    skip_on_cran()
    root <- withr::local_tempdir()
    processTestRuntime(root)
    application <- processTestManagedApplication(root)
    manager <- processTestProject(root, "cancel")
    withr::defer(closeDIAStateAdapter(manager))
    owner <- newWorkflowExecutionController(manager, application)
    withr::defer(closeWorkflowExecutionController(owner))
    store <- manager$.context$store
    before <- workflowExecutionHead(store)
    history <- manager$getHistory()
    processTestFixture(owner, list(mode = "sleep", seconds = 60))
    paths <- owner$job$paths
    token <- owner$job$request$token
    expect_identical(processTestAction(owner, until_started = TRUE)$status, "started")
    worker <- readRDS(file.path(paths$request_dir, "worker-started.rds"))$pid
    expect_true(processTestAlive(worker))
    started <- proc.time()[["elapsed"]]
    cancelWorkflowAction(owner, "reverted")
    expect_lt(proc.time()[["elapsed"]] - started, 4)
    expect_true(processTestGone(worker))
    expect_identical(workflowExecutionHead(store), before)
    expect_identical(manager$getHistory(), history)
    expect_gt(owner$control$generation, token$generation)
    expect_error(acceptWorkflowActionResult(owner, list(token = token)), class = "workflow_execution_stale")
    receipt <- readRDS(file.path(paths$request_dir, "completion.rds"))
    expect_identical(receipt$status, "reverted")
    expect_true(receipt$cleanup$worker_stopped)
    expect_false(dir.exists(file.path(paths$request_dir, "worker-scratch")))
    expect_null(owner$job)
    expect_null(application$active)
    lock <- filelock::lock(file.path(application$config$root, ".heavy-action.lock"), timeout = 0)
    expect_false(is.null(lock))
    filelock::unlock(lock)
})

test_that("L2.4 a worker past the per-step deadline ends as worker_failed", {
    skip_on_cran()
    root <- withr::local_tempdir()
    processTestRuntime(root)
    testthat::local_mocked_bindings(workflowActionWorkerDeadline = function(operation, policy) 2L, .package = "MultiScholaR")
    application <- processTestManagedApplication(root)
    manager <- processTestProject(root, "deadline")
    withr::defer(closeDIAStateAdapter(manager))
    owner <- newWorkflowExecutionController(manager, application)
    withr::defer(closeWorkflowExecutionController(owner))
    before <- workflowExecutionHead(manager$.context$store)
    processTestFixture(owner, list(mode = "sleep", seconds = 60))
    started <- proc.time()[["elapsed"]]
    done <- processTestAction(owner)
    elapsed <- proc.time()[["elapsed"]] - started
    expect_identical(done$status, "worker_failed")
    expect_identical(done$worker$Result, "timeout")
    expect_gte(elapsed, 2)
    expect_lt(elapsed, 8)
    expect_true(done$cleanup$worker_stopped)
    expect_identical(workflowExecutionHead(manager$.context$store), before)
    expect_null(application$active)
})

test_that("L2.5 closing the application stops a running worker and releases the project", {
    skip_on_cran()
    root <- withr::local_tempdir()
    processTestRuntime(root)
    application <- processTestManagedApplication(root)
    manager <- processTestProject(root, "shutdown")
    withr::defer(closeDIAStateAdapter(manager))
    owner <- newWorkflowExecutionController(manager, application)
    withr::defer(closeWorkflowExecutionController(owner))
    store <- manager$.context$store
    before <- workflowExecutionHead(store)
    processTestFixture(owner, list(mode = "sleep", seconds = 60))
    paths <- owner$job$paths
    expect_identical(processTestAction(owner, until_started = TRUE)$status, "started")
    worker <- readRDS(file.path(paths$request_dir, "worker-started.rds"))$pid
    closeWorkflowManagedApplication(application)
    expect_true(processTestGone(worker))
    expect_true(application$closed)
    expect_true(owner$closed)
    expect_identical(readRDS(file.path(paths$request_dir, "completion.rds"))$status, "controller_closed")
    expect_identical(workflowExecutionHead(store), before)
    expect_error(checkWorkflowApplication(application), class = "workflow_execution_unsupported")
    lease <- acquireWorkflowMutationLease(store)
    expect_false(is.null(lease))
    filelock::unlock(lease)
})

test_that("L2.6 after a controller is killed its worker goes, and a new controller recovers the project and runs an action", {
    skip_on_cran()
    root <- withr::local_tempdir()
    preamble <- processTestRuntime(root)
    manager <- processTestProject(root, "lost")
    store <- manager$.context$store
    project <- dirname(store$directory)
    before <- workflowExecutionHead(store)
    closeDIAStateAdapter(manager)
    lost_root <- file.path(root, "lost-application")
    dir.create(lost_root)
    record <- file.path(root, "lost.rds")
    script <- file.path(root, "lost-controller.R")
    writeLines(c(
        "arguments <- commandArgs(TRUE)",
        "if (!identical(arguments[[1L]], '-')) source(arguments[[1L]])",
        "ns <- asNamespace('MultiScholaR')",
        "application <- ns$WorkflowExecutionApplication$new()",
        "application$config <- ns$workflowProcessApplicationConfig(arguments[[2L]], TRUE, 0L)",
        "application$pid <- Sys.getpid()",
        "store <- ns$newWorkflowEntityStore(arguments[[3L]], 'lost', 'DIA')",
        "manager <- ns$newDIAStateAdapter(store, reopen = TRUE)",
        "owner <- ns$newWorkflowExecutionController(manager, application)",
        "id <- ns$runWorkflowAction('wfx_fixture', list(state = ns$getDIAStateReference(manager)), list(mode = 'sleep', seconds = 60), owner)",
        "started <- file.path(owner$job$paths$request_dir, 'worker-started.rds')",
        "repeat {",
        "    ns$pollWorkflowAction(owner)",
        "    if (file.exists(started)) break",
        "    Sys.sleep(0.05)",
        "}",
        "saveRDS(list(controller = Sys.getpid(), worker = readRDS(started)$pid, request_id = id), paste0(arguments[[4L]], '.tmp'))",
        "file.rename(paste0(arguments[[4L]], '.tmp'), arguments[[4L]])",
        "repeat {",
        "    ns$pollWorkflowAction(owner)",
        "    Sys.sleep(0.05)",
        "}"), script)
    controller <- processx::process$new(processTestRscript(), c("--vanilla", script, preamble, lost_root, project, record),
        env = c("current", processTestLibraries()), stdout = file.path(root, "lost-controller.log"), stderr = "2>&1", cleanup = TRUE)
    withr::defer(controller$kill())
    expect_true(processTestFile(record, seconds = 240))
    pids <- readRDS(record)
    withr::defer(if (processTestAlive(pids$worker)) ps::ps_kill(ps::ps_handle(as.integer(pids$worker))))
    expect_true(processTestAlive(pids$worker))
    ps::ps_kill(ps::ps_handle(as.integer(pids$controller)))
    expect_true(processTestGone(pids$controller))
    expect_true(processTestGone(pids$worker))
    application <- processTestManagedApplication(root)
    store <- newWorkflowEntityStore(project, "lost", "DIA")
    manager <- newDIAStateAdapter(store, reopen = TRUE)
    withr::defer(closeDIAStateAdapter(manager))
    owner <- newWorkflowExecutionController(manager, application)
    withr::defer(closeWorkflowExecutionController(owner))
    expect_identical(workflowExecutionHead(store), before)
    expect_identical(readWorkflowExecutionJournal(workflowExportPaths(store$directory, pids$request_id))$status, "controller_lost")
    expect_null(owner$control$active_request)
    processTestFixture(owner, list(mode = "identity"))
    expect_identical(processTestAction(owner)$status, "complete")
    expect_false(identical(workflowExecutionHead(store)$entity_id, before$entity_id))
})

test_that("L2.8 the plain launch is the GUI application in the calling process; suites, a supervisor and a journey record are refused", {
    root <- withr::local_tempdir()
    withr::local_options(multischolar.workflow_launcher = "process")
    arguments <- NULL
    testthat::local_mocked_bindings(workflowExecutionBuildIdentity = function() strrep("f", 64),
        workflowHostAvailable = function() 64 * 1024^3,
        workflowSystemd = function(...) stop("the plain launch must not reach systemd"),
        workflowManagedMain = function(arguments) arguments <<- arguments, .package = "MultiScholaR")
    for (mode in c("protocol", "resources", "adapters")) expect_error(
        launchWorkflowApplication(root, diagnostic = TRUE, mode = mode, test_root = root), class = "workflow_execution_unsupported")
    expect_error(launchWorkflowApplication(root, supervisor_unit = "multischolar-dia003-browser-smoke.abcdef.service"),
        class = "workflow_execution_unsupported")
    expect_error(launchWorkflowApplication(file.path(root, "absent")), class = "workflow_execution_unsupported")
    journey <- file.path(root, "journey")
    dir.create(journey)
    saveRDS(list(schema = 1L, programme = "kit-journey-v1", deadline = 600), file.path(journey, "gui-member.rds"))
    expect_error(launchWorkflowApplication(journey), class = "workflow_execution_unsupported")
    expect_null(arguments)
    expect_false(file.exists(file.path(root, "application.rds")))
    config <- launchWorkflowApplication(root, port = 4321L)
    expect_identical(config$launcher, "process")
    expect_identical(config$mode, "gui")
    expect_identical(config$port, 4321L)
    expect_identical(config$policy, workflowExecutionPolicy(FALSE))
    expect_null(config$slice)
    expect_identical(arguments, c(file.path(config$root, "application.rds"),
        digest::digest(file = file.path(config$root, "application.rds"), algo = "sha256")))
    expect_identical(readWorkflowExecutionRequest(arguments), config)
    expect_error(launchWorkflowApplication(root), "already used")
})

test_that("L2.8 the plain controller is the calling process: no control group, checked by class, owner and launcher, closed at exit", {
    root <- withr::local_tempdir()
    withr::local_options(multischolar.workflow_launcher = "process")
    served <- NULL
    testthat::local_mocked_bindings(workflowExecutionBuildIdentity = function() strrep("f", 64),
        workflowHostAvailable = function() 64 * 1024^3,
        workflowSystemd = function(...) stop("the plain launch must not reach systemd"),
        run_app = function(...) list(...)$workflow_execution, .package = "MultiScholaR")
    testthat::local_mocked_bindings(runApp = function(appDir, ...) {
        served <<- list(application = appDir, checked = checkWorkflowApplication(appDir), closed = appDir$closed, options = list(...))
        invisible(NULL)
    }, .package = "shiny")
    config <- launchWorkflowApplication(root, port = 4322L)
    expect_s3_class(served$application, "WorkflowExecutionApplication")
    expect_identical(served$application$config, config)
    expect_null(served$application$cgroup)
    expect_identical(served$application$pid, Sys.getpid())
    expect_true(served$checked)
    expect_false(served$closed)
    expect_identical(served$options$port, 4322L)
    expect_true(served$application$closed)
    # Another process, or a closed application, is not a managed application.
    other <- WorkflowExecutionApplication$new()
    other$config <- config
    other$pid <- Sys.getpid() + 1L
    expect_error(checkWorkflowApplication(other), class = "workflow_execution_unsupported")
})

test_that("L2.8 without the option a system with systemctl takes the systemd launch, as before", {
    skip_if_not(identical(Sys.info()[["sysname"]], "Linux") && nzchar(Sys.which("systemctl")))
    root <- withr::local_tempdir()
    reached <- function(args, ...) cli::cli_abort("the systemd launch was reached", class = "process_test_systemd")
    testthat::local_mocked_bindings(workflowHostAvailable = function() 64 * 1024^3, workflowSystemd = reached,
        workflowManagedMain = function(arguments) stop("the plain launch was reached"), .package = "MultiScholaR")
    withr::with_options(list(multischolar.workflow_launcher = NULL),
        expect_error(launchWorkflowApplication(root), class = "process_test_systemd"))
    withr::with_options(list(multischolar.workflow_launcher = "systemd"),
        expect_error(launchWorkflowApplication(root), class = "process_test_systemd"))
    expect_false(file.exists(file.path(root, "application.rds")))
})
