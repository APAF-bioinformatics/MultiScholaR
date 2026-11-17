# KIT-131 step 7 (scratch driver, not a repository test; tests/ are not edited): one GUI-entry lane with every managed
# step on a plain worker (launcher "process"), against the installed package, compared with its frozen oracle.
#
# Reuses tests/testthat/helper-gui-entry.R unchanged. Three things differ from the in-process lane:
#  1. the module's session carries a plain managed application (`session$userData$workflow_execution`, what
#     `run_app(workflow_execution = app)` gives a session under `launchWorkflowProcessApplication`); the application is
#     built as `workflowManagedMain` builds it for the plain launcher (config, pid, check, close on exit);
#  2. after every input the lane sets (each click), the driver waits until the omic's runner has no active action, advancing
#     the mock session's timers (the 50 ms poll) and applying the browser double's queued updates, as a user waits for a
#     step before the next click; before each capture it reads the outputs once and waits again, so displays the worker
#     draws (`workflowOmicsDisplay`) are delivered before the capture reads them;
#  3. the comparison reports the managed step's own working notification ("Applying <omic> step ...") separately: the
#     in-process path has none. Everything else must be identical to the oracle.
# Usage (cwd = repository root, R_LIBS with the installed scratch library first):
#   Rscript lane_worker.R <lane_id> <out_dir>
arguments <- commandArgs(TRUE)
lane_id <- arguments[[1L]]
out_dir <- arguments[[2L]]
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
.libPaths(c(Sys.getenv("KIT131_LIBRARY"), .libPaths()))
Sys.setenv(R_LIBS = paste(.libPaths(), collapse = .Platform$path.sep), NOT_CRAN = "true")
options(multischolar.workflow_launcher = "process")
ns <- asNamespace("MultiScholaR")
cat("MultiScholaR from", getNamespaceInfo(ns, "path"), "\n")
stopifnot(normalizePath(getNamespaceInfo(ns, "path")) == normalizePath(file.path(Sys.getenv("KIT131_LIBRARY"), "MultiScholaR")))
suppressMessages(library(shiny))
env <- new.env(parent = ns)
for (helper in c("helper-workflow-omics-residency.R", "helper-e2e-enrichment-doubles.R", "helper-gui-entry.R"))
    sys.source(file.path("tests/testthat", helper), env)

# The plain application, as workflowManagedMain builds it for the plain launcher.
app_root <- file.path(out_dir, "application")
dir.create(app_root)
app <- ns$WorkflowExecutionApplication$new()
app$config <- ns$workflowProcessApplicationConfig(app_root, FALSE, 0L)
stopifnot(identical(app$config$launcher, "process"))
app$pid <- Sys.getpid()
stopifnot(ns$checkWorkflowApplication(app))

driver <- new.env()
driver$waits <- 0L
driver$worker_seconds <- 0
original_spec <- env$guiEntryServerSpec
env$guiEntryServerSpec <- function(omic) {
    spec <- original_spec(omic)
    module <- spec$server
    driver$omic <- omic
    spec$server <- function(id, ...) {
        domain <- shiny::getDefaultReactiveDomain()
        domain$userData$workflow_execution <- app
        module(id, ...)
    }
    spec
}
original_double <- env$guiEntryInstallBrowserDouble
env$guiEntryInstallBrowserDouble <- function(log, env = parent.frame()) {
    driver$log <- log
    original_double(log, env)
}
# Wait until the omic's runner is idle: timers advanced in the mock session, real time for the worker, queued updates
# applied as the lane's own `apply_updates` applies them.
drain <- function(api, ceiling = 900) {
    root <- api$session$rootScope()
    started <- proc.time()[["elapsed"]]
    seen_active <- FALSE
    repeat {
        runner <- api$session$userData$workflow_omics_actions[[driver$omic]]
        if (is.null(runner)) stop("no managed runner on the session: the lane would run in process")
        active <- shiny::isolate(runner$step$active())
        if (!is.null(shiny::isolate(runner$step$pending()))) stop("a managed publication is pending")
        if (is.null(active)) break
        seen_active <- TRUE
        if (proc.time()[["elapsed"]] - started > ceiling) stop("bounded managed wait exceeded")
        Sys.sleep(0.05)
        root$elapse(50)
        for (i in seq_len(50L)) {
            if (!length(driver$log$updates)) break
            pending <- driver$log$updates; driver$log$updates <- list()
            values <- stats::setNames(lapply(pending, `[[`, "value"), vapply(pending, `[[`, "", "id"))
            values <- values[!duplicated(names(values), fromLast = TRUE)]
            do.call(root$setInputs, values)
        }
    }
    if (seen_active) {
        driver$waits <- driver$waits + 1L
        driver$worker_seconds <- driver$worker_seconds + (proc.time()[["elapsed"]] - started)
    }
    invisible(seen_active)
}
original_steps <- env$guiEntrySteps
env$guiEntrySteps <- function(lane) {
    lapply(original_steps(lane), function(step) {
        run <- step$run
        step$run <- function(api) {
            managed <- api
            managed$set <- function(values) { api$set(values); drain(api) }
            managed$upload <- function(id, path) { api$upload(id, path); drain(api) }
            managed$defaults <- function(ui_fn, id) { api$defaults(ui_fn, id); drain(api) }
            run(managed)
            drain(api)
            # The displays the capture will read: asked for now, delivered before the capture.
            for (round in seq_len(10L)) {
                api$session$flushReact()
                env$guiEntryOutputs(api$session$rootScope(), character())
                if (!drain(api)) break
            }
            if (!nzchar(api$session$userData$workflow_execution$config$launcher %||% "")) stop("lost the application")
        }
        step
    })
}

started <- Sys.time()
record <- tryCatch(env$guiEntryRunLane(lane_id), error = function(error) error)
ns$closeWorkflowManagedApplication(app)
if (inherits(record, "error")) {
    cat("KIT-131 LANE", lane_id, "ERROR:", conditionMessage(record), "\n")
    quit(status = 2L)
}
receipts <- list.files(app_root, "receipt|journal", recursive = TRUE)
requests <- list.files(app_root, recursive = TRUE)
steps <- env$guiEntryComparable(record)
env$guiEntryWriteJson(list(lane = lane_id, steps = steps), file.path(out_dir, paste0(lane_id, "-worker.json")))
oracle <- jsonlite::read_json(env$guiEntryOraclePath(lane_id), simplifyVector = FALSE)
working <- paste0("^Applying ", driver$omic, " step \\.\\.\\.$")
strip <- function(step) {
    notes <- step$notifications
    keep <- vapply(notes, function(note) !(identical(note$type, "default") && grepl(working, note$message)), logical(1))
    step$notifications <- notes[keep]
    list(step = step, working = sum(!keep))
}
identical_raw <- 0L; identical_stripped <- 0L; working_total <- 0L
for (id in names(oracle$steps)) {
    mine <- steps[[id]]
    if (identical(mine, oracle$steps[[id]])) identical_raw <- identical_raw + 1L
    stripped <- strip(mine)
    working_total <- working_total + stripped$working
    if (identical(stripped$step, oracle$steps[[id]])) {
        identical_stripped <- identical_stripped + 1L
    } else {
        fields <- union(names(stripped$step), names(oracle$steps[[id]]))
        differ <- fields[!vapply(fields, function(f) identical(stripped$step[[f]], oracle$steps[[id]][[f]]), logical(1))]
        cat("KIT-131 LANE DIFF", lane_id, "step", id, "fields:", paste(differ, collapse = ", "), "\n")
    }
}
cat(sprintf("KIT-131 LANE %s: steps %d (oracle %d, order %s); identical %d raw, %d without the working notification (%d removed); managed waits %d, %.1f s waiting; %d files under the application root; lane %.1f s\n",
    lane_id, length(steps), length(oracle$steps), identical(names(steps), names(oracle$steps)), identical_raw, identical_stripped,
    working_total, driver$waits, driver$worker_seconds, length(requests), record$seconds))
