# MultiScholaR: Interactive Multi-Omics Analysis
# Copyright (C) 2024-2026 Ignatius Pang, William Klare, and APAF-bioinformatics
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

#' Run the Shiny Application
#'
#' @param test_mode Logical. When TRUE, sets the \code{multischolar.test_mode}
#'   option before launching, enabling test-only UI elements (stable selectors,
#'   state digest output) and disabling native file dialogs.
#' @param ... arguments to pass to golem_opts.
#' See `?golem::get_golem_options` for more details.
#' @inheritParams shiny::shinyApp
#'
#' @export
#' @importFrom shiny shinyApp addResourcePath
#' @importFrom golem with_golem_options
run_app <- function(
  onStart = NULL,
  options = list(),
  enableBookmarking = NULL,
  uiPattern = "/",
  test_mode = FALSE,
  ...
) {
  if (isTRUE(test_mode)) {
    options(multischolar.test_mode = TRUE) # KIT-064 (roar BUG-22): intentional for the process's lifetime; the harness's test mode
  }
  # Set up resource path for static files in inst/shiny/www before UI renders
  # Handle both development (pkgload) and installed package scenarios
  shiny_resource_path <- system.file("shiny/www", package = "MultiScholaR")
  if (shiny_resource_path == "" || !dir.exists(shiny_resource_path)) {
    # Development mode: try to find inst/shiny/www relative to common locations
    # Check current directory and parent directories
    possible_paths <- c(
      file.path(getwd(), "inst", "shiny", "www"),
      file.path(dirname(getwd()), "inst", "shiny", "www"),
      file.path(".", "inst", "shiny", "www")
    )
    for (path in possible_paths) {
      if (dir.exists(path)) {
        shiny_resource_path <- normalizePath(path)
        break
      }
    }
  }
  if (dir.exists(shiny_resource_path)) {
    shiny::addResourcePath("www", shiny_resource_path)
  }
  
  execution <- list(...)$workflow_execution
  application_server <- app_server
  if (!is.null(execution)) {
    checkWorkflowApplication(execution)
    application_server <- workflowManagedServer(execution)
  }

  # Combine with user-provided onStart if any
  combined_onStart <- function() {
    if (!is.null(onStart)) {
      onStart()
    }
  }
  
  with_golem_options(
    app = shiny::shinyApp(
      ui = app_ui,
      server = application_server,
      onStart = combined_onStart,
      options = options, 
      enableBookmarking = enableBookmarking, 
      uiPattern = uiPattern
    ), 
    golem_opts = list(...)
  )
}

#' The packaged application's entry (KIT-134 P3): what the bundle's shell starts, as
#' `Rscript --vanilla entry.R <port> <run folder> <shell pid>`. The shell's start-up probe (`MSR_LAUNCHER_PROBE`)
#' names the first launcher to try; the order is the systemd worker, then the plain launcher after a trial unit, then
#' the application in this process. The choice is passed to the managed launch as the option
#' `multischolar.workflow_launcher`. One line `MSR-LAUNCHER <json>` on standard output says which launcher runs and why.
#' Nothing is installed or downloaded (`loadDependencies()` is never called).
#' @param arguments Character vector: the loopback port, the run folder, the shell's process id.
#' @return NULL invisibly, once the application has stopped (systemd: once the shell or the controller is gone).
#' @noRd
runPackagedApplication <- function(arguments) {
    if (!is.character(arguments) || length(arguments) != 3L)
        cli::cli_abort("the packaged entry needs a port, a run folder and the shell's pid",
            class = "workflow_execution_invalid")
    port <- suppressWarnings(as.integer(arguments[[1L]]))
    shell <- suppressWarnings(as.integer(arguments[[3L]]))
    root <- arguments[[2L]]
    if (is.na(port) || is.na(shell) || !dir.exists(root))
        cli::cli_abort("invalid packaged entry arguments: {.val {arguments}}", class = "workflow_execution_invalid")
    probe <- readPackagedLauncherProbe(Sys.getenv("MSR_LAUNCHER_PROBE"))
    reason <- probe$reason
    watchPackagedShell(shell)
    if (identical(probe$launcher, "systemd")) {
        config <- tryCatch(withr::with_options(list(multischolar.workflow_launcher = "systemd"),
            launchWorkflowApplication(packagedLauncherFolder(root, "systemd"), mode = "gui", port = port)),
            workflow_execution_unsupported = function(error) error, workflow_execution_refused = function(error) error)
        if (!inherits(config, "condition")) {
            announcePackagedLauncher("systemd", reason, slice = config$slice)
            return(holdPackagedSystemd(config, shell))
        }
        reason <- paste("the systemd worker could not start:", conditionMessage(config))
    }
    trial <- tryCatch(runPackagedTrialUnit(packagedLauncherFolder(root, "trial")),
        workflow_execution_unsupported = function(error) conditionMessage(error))
    if (isTRUE(trial)) {
        announcePackagedLauncher("process", reason)
        withr::local_options(multischolar.workflow_launcher = "process")
        launchWorkflowApplication(packagedLauncherFolder(root, "process"), mode = "gui", port = port)
        return(invisible(NULL))
    }
    announcePackagedLauncher("in-process", paste0(reason, "; the plain launcher's trial failed: ", trial))
    withr::local_dir(root)
    shiny::runApp(run_app(), host = "127.0.0.1", port = port, launch.browser = FALSE)
    invisible(NULL)
}

#' Read the shell's start-up probe (KIT-134 P3, decision U1: the host facts are read by the shell).
#' @param text The value of `MSR_LAUNCHER_PROBE`: JSON with `launcher` (`"systemd"` or `"process"`) and `reason`.
#' @return `list(launcher, reason)`; a missing or unreadable probe gives the plain launcher, with the reason.
#' @noRd
readPackagedLauncherProbe <- function(text) {
    if (!nzchar(text)) return(list(launcher = "process", reason = "no start-up probe from the shell"))
    probe <- tryCatch(jsonlite::fromJSON(text, simplifyVector = FALSE), error = function(error) NULL)
    if (!is.list(probe) || is.null(names(probe)) || !isTRUE(probe$launcher %in% c("systemd", "process")))
        return(list(launcher = "process", reason = "the start-up probe from the shell is unreadable"))
    list(launcher = probe$launcher, reason = paste(as.character(probe$reason %||% ""), collapse = " "))
}

#' A new folder for one launcher's managed application inside the run folder.
#' @param root The run folder.
#' @param name The launcher's folder name.
#' @return The folder's path.
#' @noRd
packagedLauncherFolder <- function(root, name) {
    folder <- file.path(root, name)
    dir.create(folder, showWarnings = FALSE)
    folder
}

#' Say which launcher runs: one line for the shell on standard output (its protocol, not a diagnostic), and the log.
#' @param launcher `"systemd"`, `"process"` or `"in-process"`.
#' @param reason Why.
#' @param ... Further fields (the systemd slice, which the shell stops when its window closes).
#' @return NULL invisibly.
#' @noRd
announcePackagedLauncher <- function(launcher, reason, ...) {
    line <- jsonlite::toJSON(c(list(launcher = launcher, reason = reason), list(...)), auto_unbox = TRUE)
    cat("MSR-LAUNCHER ", line, "\n", sep = "")
    flush(stdout())
    # The JSON's braces are not glue expressions.
    logger::log_info(logger::skip_formatter(paste("MSR-LAUNCHER", line)))
    invisible(NULL)
}

#' Whether the shell's process is still running.
#' @param shell The shell's process id.
#' @return TRUE or FALSE.
#' @noRd
packagedShellAlive <- function(shell) {
    tryCatch(ps::ps_is_running(ps::ps_handle(shell)), error = function(error) FALSE)
}

#' End R when the shell is gone (KIT-130: R must not outlive a shell that was killed). It runs from R's event loop
#' every 2 s, so it cannot fire while a step blocks the loop (audit F9).
#' @param shell The shell's process id.
#' @return NULL invisibly.
#' @noRd
watchPackagedShell <- function(shell) {
    check <- function() {
        if (!packagedShellAlive(shell)) {
            logger::log_info(sprintf("the shell (pid %d) is gone: stopping R", shell))
            quit(save = "no", status = 0L)
        }
        later::later(check, 2)
    }
    later::later(check, 2)
    invisible(NULL)
}

#' Start and stop one trivial unit through the plain launcher, so a host that cannot start a child is found at start,
#' and check the host headroom the plain launch requires.
#' @param folder The trial's folder.
#' @return TRUE, or why the trial failed; `workflow_execution_unsupported` when the unit cannot start.
#' @noRd
runPackagedTrialUnit <- function(folder) {
    unit <- paste0("multischolar-trial-", Sys.getpid())
    rscript <- file.path(R.home("bin"), if (identical(.Platform$OS.type, "windows")) "Rscript.exe" else "Rscript")
    workflowProcessStart(unit, rscript, c("--vanilla", "-e", "quit(status = 0)"), character(), folder,
        file.path(folder, "trial.log"), 30)
    on.exit({
        workflowProcessStop(unit)
        workflowProcessForget(unit)
    }, add = TRUE)
    repeat {
        state <- workflowProcessState(unit)
        if (!identical(state$SubState, "running")) break
        Sys.sleep(0.1)
    }
    if (!identical(state$Result, "success"))
        return(sprintf("the trial unit ended with %s (exit %s)", state$Result, state$ExitStatus))
    needed <- workflowExecutionPolicy()$startup_headroom
    available <- workflowHostAvailable()
    if (available < needed)
        return(sprintf("%.1f GiB of memory available, %.0f GiB needed", available / 1024^3, needed / 1024^3))
    TRUE
}

#' Under the systemd worker the entry stays alive and stops the application's slice when the shell or the controller
#' unit is gone, so nothing of the application outlives the shell.
#' @param config The managed application's configuration.
#' @param shell The shell's process id.
#' @return NULL invisibly.
#' @noRd
holdPackagedSystemd <- function(config, shell) {
    on.exit(workflowSystemd(c("stop", config$slice), required = FALSE), add = TRUE)
    repeat {
        if (!packagedShellAlive(shell)) break
        state <- workflowUnitState(config, config$controller_unit)
        if (!isTRUE(state$available) || !isTRUE(state$ActiveState %in% c("active", "activating"))) break
        Sys.sleep(2)
    }
    invisible(NULL)
}

