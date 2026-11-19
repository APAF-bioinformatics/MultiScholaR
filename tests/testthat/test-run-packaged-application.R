# KIT-134 P3, the package side: the packaged application's entry (`runPackagedApplication`), its launcher order, and
# the environment a systemd unit started from the bundle receives (`workflowLaunchService` with `MSR_BUNDLE`).
# Launchers, the trial unit, the shell watchdog and shiny are mocked; nothing is started.

packagedRun <- function(probe, envir = parent.frame()) {
    root <- withr::local_tempdir(.local_envir = envir)
    withr::local_envvar(MSR_LAUNCHER_PROBE = probe, .local_envir = envir)
    list(root = root, arguments = c("43210", root, "999"))
}

launcherLines <- function(output) {
    lines <- grep("^MSR-LAUNCHER ", output, value = TRUE)
    lapply(sub("^MSR-LAUNCHER ", "", lines), jsonlite::fromJSON)
}

test_that("the probe from the shell is read; a missing or unreadable probe means the plain launcher", {
    expect_identical(readPackagedLauncherProbe('{"launcher":"systemd","reason":"every host requirement holds"}'),
        list(launcher = "systemd", reason = "every host requirement holds"))
    expect_identical(readPackagedLauncherProbe('{"launcher":"process","reason":"kernel 6.8 is older than 6.12","facts":{}}')$launcher,
        "process")
    for (text in c("", "not json", '{"launcher":"docker"}', "[1, 2]")) {
        probe <- readPackagedLauncherProbe(text)
        expect_identical(probe$launcher, "process", info = text)
        expect_true(nzchar(probe$reason), info = text)
    }
})

test_that("invalid entry arguments are refused before anything starts", {
    root <- withr::local_tempdir()
    for (arguments in list(c("43210", root), c("port", root, "999"), c("43210", file.path(root, "absent"), "999"),
        c("43210", root, "pid"))) {
        expect_error(runPackagedApplication(arguments), class = "workflow_execution_invalid")
    }
})

test_that("systemd from the probe: the systemd launch under the option, one MSR-LAUNCHER line, then the hold", {
    run <- packagedRun('{"launcher":"systemd","reason":"systemd worker: every host requirement holds"}')
    calls <- character()
    testthat::local_mocked_bindings(
        watchPackagedShell = function(shell) calls <<- c(calls, paste("watch", shell)),
        launchWorkflowApplication = function(root, mode = "gui", port = 3838L, ...) {
            calls <<- c(calls, paste("launch", getOption("multischolar.workflow_launcher"), basename(root), mode, port))
            list(slice = "multischolarwfx0123456789abcdef.slice")
        },
        holdPackagedSystemd = function(config, shell) calls <<- c(calls, paste("hold", config$slice, shell)),
        runPackagedTrialUnit = function(dir) stop("the trial must not run after a systemd launch"),
        loadDependencies = function(...) stop("loadDependencies must never be called"),
        .package = "MultiScholaR")
    output <- utils::capture.output(runPackagedApplication(run$arguments))
    expect_identical(calls, c("watch 999", "launch systemd systemd gui 43210",
        "hold multischolarwfx0123456789abcdef.slice 999"))
    lines <- launcherLines(output)
    expect_length(lines, 1L)
    expect_identical(lines[[1L]]$launcher, "systemd")
    expect_identical(lines[[1L]]$slice, "multischolarwfx0123456789abcdef.slice")
    expect_null(getOption("multischolar.workflow_launcher"))
})

test_that("a systemd launch that is unsupported or refused falls to the plain launcher after its trial", {
    for (class in c("workflow_execution_unsupported", "workflow_execution_refused")) {
        run <- packagedRun('{"launcher":"systemd","reason":"systemd worker: every host requirement holds"}')
        calls <- character()
        testthat::local_mocked_bindings(
            watchPackagedShell = function(shell) invisible(NULL),
            launchWorkflowApplication = function(root, mode = "gui", port = 3838L, ...) {
                launcher <- getOption("multischolar.workflow_launcher")
                calls <<- c(calls, paste("launch", launcher, basename(root)))
                if (identical(launcher, "systemd")) cli::cli_abort("aggregate memory controls are unavailable", class = class)
                list(launcher = "process")
            },
            runPackagedTrialUnit = function(dir) { calls <<- c(calls, paste("trial", basename(dir))); TRUE },
            holdPackagedSystemd = function(config, shell) stop("no hold without a systemd launch"),
            loadDependencies = function(...) stop("loadDependencies must never be called"),
            .package = "MultiScholaR")
        output <- utils::capture.output(runPackagedApplication(run$arguments))
        expect_identical(calls, c("launch systemd systemd", "trial trial", "launch process process"), info = class)
        lines <- launcherLines(output)
        expect_length(lines, 1L)
        expect_identical(lines[[1L]]$launcher, "process", info = class)
        expect_match(lines[[1L]]$reason, "aggregate memory controls are unavailable", fixed = TRUE, info = class)
    }
})

test_that("process from the probe: the trial unit, then the plain launch under the option; systemd is not tried", {
    run <- packagedRun('{"launcher":"process","reason":"kernel 6.8 is older than 6.12"}')
    calls <- character()
    testthat::local_mocked_bindings(
        watchPackagedShell = function(shell) invisible(NULL),
        launchWorkflowApplication = function(root, mode = "gui", port = 3838L, ...) {
            calls <<- c(calls, paste("launch", getOption("multischolar.workflow_launcher"), basename(root), port))
            list(launcher = "process")
        },
        runPackagedTrialUnit = function(dir) { calls <<- c(calls, "trial"); TRUE },
        loadDependencies = function(...) stop("loadDependencies must never be called"),
        .package = "MultiScholaR")
    output <- utils::capture.output(runPackagedApplication(run$arguments))
    expect_identical(calls, c("trial", "launch process process 43210"))
    lines <- launcherLines(output)
    expect_length(lines, 1L)
    expect_identical(lines[[1L]]$launcher, "process")
    expect_identical(lines[[1L]]$reason, "kernel 6.8 is older than 6.12")
})

test_that("a failing trial unit runs the application in this process, said once with the reason", {
    for (trial in list(function(dir) "the trial unit ended with exit-code (exit 1)",
        function(dir) cli::cli_abort("managed unit could not start: no such file", class = "workflow_execution_unsupported"))) {
        run <- packagedRun('{"launcher":"process","reason":"win32: no systemd"}')
        served <- NULL
        testthat::local_mocked_bindings(
            watchPackagedShell = function(shell) invisible(NULL),
            launchWorkflowApplication = function(...) stop("no managed launch after a failed trial"),
            runPackagedTrialUnit = trial,
            run_app = function(...) "application",
            loadDependencies = function(...) stop("loadDependencies must never be called"),
            .package = "MultiScholaR")
        testthat::local_mocked_bindings(runApp = function(appDir, ...) {
            served <<- c(list(application = appDir, directory = getwd()), list(...))
            invisible(NULL)
        }, .package = "shiny")
        output <- utils::capture.output(runPackagedApplication(run$arguments))
        expect_identical(served$application, "application")
        expect_identical(served$host, "127.0.0.1")
        expect_identical(served$port, 43210L)
        expect_false(served$launch.browser)
        expect_identical(normalizePath(served$directory), normalizePath(run$root))
        lines <- launcherLines(output)
        expect_length(lines, 1L)
        expect_identical(lines[[1L]]$launcher, "in-process")
        expect_match(lines[[1L]]$reason, "win32: no systemd", fixed = TRUE)
        expect_match(lines[[1L]]$reason, "trial", fixed = TRUE)
    }
})

test_that("the systemd hold stops the slice when the shell or the controller is gone", {
    config <- list(slice = "multischolarwfx0123456789abcdef.slice",
        controller_unit = "multischolarwfx0123456789abcdef-controller.service")
    for (case in list(list(alive = FALSE, state = "active"), list(alive = TRUE, state = "inactive"))) {
        stopped <- character()
        testthat::local_mocked_bindings(
            packagedShellAlive = function(shell) case$alive,
            workflowUnitState = function(config, unit) list(available = TRUE, ActiveState = case$state),
            workflowSystemd = function(args, required = TRUE, ...) stopped <<- c(stopped, paste(args, collapse = " ")),
            .package = "MultiScholaR")
        holdPackagedSystemd(config, 999L)
        expect_identical(stopped, "stop multischolarwfx0123456789abcdef.slice")
    }
})

serviceArguments <- function(bundle = NULL, root = withr::local_tempdir(.local_envir = envir), envir = parent.frame()) {
    path <- file.path(root, "request.rds")
    config <- list(root = root, slice = "multischolarwfx0123456789abcdef.slice", library = "/library", mode = "gui",
        policy = list(diagnostic = FALSE))
    captured <- NULL
    testthat::local_mocked_bindings(workflowSystemd = function(args, required = TRUE, command = "systemctl") {
        captured <<- args
        invisible(NULL)
    }, .package = "MultiScholaR", .env = envir)
    if (is.null(bundle)) withr::local_envvar(MSR_BUNDLE = NA, .local_envir = envir) else
        withr::local_envvar(bundle, .local_envir = envir)
    workflowLaunchService(config, "multischolarwfx0123456789abcdef-w1.service", "worker", path, "hash",
        cap = 2 * 1024^3, deadline = 60L)
    list(args = captured, root = root, path = path)
}

test_that("without MSR_BUNDLE a systemd unit's arguments are exactly as before", {
    run <- serviceArguments()
    expected <- c("--unit", "multischolarwfx0123456789abcdef-w1.service", "--slice", "multischolarwfx0123456789abcdef.slice",
        "--property=Type=exec", "--property=KillMode=control-group", "--property=TimeoutStopSec=1",
        "--property=MemoryAccounting=yes", "--property=CPUAccounting=yes", "--property=IOAccounting=yes",
        "--property=MemorySwapMax=0", "--property=RuntimeMaxSec=60",
        paste0("--property=StandardOutput=append:", file.path(run$root, "worker-request.rds.log")),
        "--property=StandardError=inherit",
        paste0("--setenv=R_LIBS=", paste(unique(c("/library", .libPaths())), collapse = .Platform$path.sep)),
        "--setenv=OMP_NUM_THREADS=1", "--setenv=OPENBLAS_NUM_THREADS=1", "--setenv=MKL_NUM_THREADS=1",
        "--setenv=DUCKPLYR_FALLBACK_AUTOUPLOAD=0", "--setenv=DUCKPLYR_FALLBACK_COLLECT=0",
        "--remain-after-exit", "--property=MemoryMax=2147483648", "--property=OOMPolicy=kill",
        paste0("--setenv=TMPDIR=", file.path(run$root, "worker-scratch")), "--setenv=MULTISCHOLAR_WORKFLOW_WORKER=1",
        file.path(R.home("bin"), "Rscript"), "--vanilla", "-e", "MultiScholaR:::workflowWorkerMain(commandArgs(TRUE))",
        run$path, "hash")
    expect_identical(run$args, expected)
})

test_that("inside the bundle a systemd unit also gets the bundle's environment, where set, before the command", {
    bundle <- c(MSR_BUNDLE = "/opt/msr", RHOME = "/opt/msr/R", LD_LIBRARY_PATH = "/opt/msr/syslib",
        R_LIBS_USER = "/opt/msr/R/no-user-library", R_LIBS_SITE = "/opt/msr/R/no-site-library",
        RETICULATE_PYTHON = "/opt/msr/python/bin/python3", PYTHONNOUSERSITE = "1", RSTUDIO_PANDOC = "/opt/msr/pandoc",
        PATH = "/usr/bin:/bin", LANG = "C.UTF-8", LC_CTYPE = "C.UTF-8", LC_ALL = NA)
    root <- withr::local_tempdir()
    plain <- serviceArguments(root = root)$args
    inside <- serviceArguments(bundle, root = root)$args
    added <- c("--setenv=RHOME=/opt/msr/R", "--setenv=LD_LIBRARY_PATH=/opt/msr/syslib",
        "--setenv=R_LIBS_USER=/opt/msr/R/no-user-library", "--setenv=R_LIBS_SITE=/opt/msr/R/no-site-library",
        "--setenv=RETICULATE_PYTHON=/opt/msr/python/bin/python3", "--setenv=PYTHONNOUSERSITE=1",
        "--setenv=RSTUDIO_PANDOC=/opt/msr/pandoc", "--setenv=MSR_BUNDLE=/opt/msr", "--setenv=PATH=/usr/bin:/bin",
        "--setenv=LANG=C.UTF-8", "--setenv=LC_CTYPE=C.UTF-8")
    # The same vector with the bundle's values inserted just before the command (the bundled R/bin/Rscript).
    expect_identical(inside, append(plain, added, after = match(file.path(R.home("bin"), "Rscript"), plain) - 1L))
    expect_false(any(startsWith(inside, "--setenv=LC_ALL=")))
})
