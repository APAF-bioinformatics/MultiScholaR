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

#' The source tree when the package is loaded from source, `""` when it is installed.
processTestSourceRoot <- function() {
    if (nzchar(system.file("R", "MultiScholaR.rdb", package = "MultiScholaR"))) return("")
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
    # `Rscript.exe` is a front end on Windows: the R session is its child, so a stop must reach the tree.
    if (identical(.Platform$OS.type, "windows")) expect_false(identical(pids$unit, launched)) else
        expect_identical(pids$unit, launched)
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
        "if (nzchar(arguments[[1L]])) suppressMessages(pkgload::load_all(arguments[[1L]], quiet = TRUE, helpers = FALSE, attach_testthat = FALSE))",
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
