# KIT-131 step 7 on Windows (build-time verification, not package code; nothing in R/ or tests/ is changed): one
# GUI-entry lane with every managed step on a plain worker (launcher "process"), with the installed application's R,
# compared with the same lane in process on the same runner (design A1) and with its frozen oracle.
#
# The driver is docs/qa/workflow-execution/kit/kit-131/step8/lane_worker.R (step 7 on Linux), adapted:
#  - the package is the installed application's (its R's own library); the test tools are a separate library after it;
#  - an expectation inside a step is counted and muffled, not fatal (the lane runs outside test_that, as on Linux);
#  - the record is compared with the in-process record of the same job (run-lanes.R) and with the oracle.
# Unchanged from the Linux driver: the module's session carries a plain managed application, built as workflowManagedMain
# builds it for the plain launcher; after every input the driver waits until the omic's runner is idle; the managed
# step's own working notification ("Applying <omic> step ...") is set aside before comparing (the in-process path has none).
#
# Usage: Rscript --vanilla lane-worker.R <checkout> <test-tool library> <output directory> <lane id> [<in-process record>]
# Prints KIT-131 lines. Exit 0 when the lane ran to its end on workers (whatever its differences), 2 when it did not.

arguments <- commandArgs(trailingOnly = TRUE)
if (length(arguments) < 4L) stop("usage: Rscript --vanilla lane-worker.R <checkout> <test-tool library> <output directory> <lane id> [<in-process record>]")
checkout <- normalizePath(arguments[[1L]], mustWork = TRUE)
tools_lib <- normalizePath(arguments[[2L]], mustWork = TRUE)
out_dir <- arguments[[3L]]
lane_id <- arguments[[4L]]
inprocess_path <- if (length(arguments) >= 5L) arguments[[5L]] else ""
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
out_dir <- normalizePath(out_dir)
setwd(checkout)
.libPaths(c(.libPaths(), tools_lib))
# Workers find the same libraries.
Sys.setenv(R_LIBS = paste(.libPaths(), collapse = .Platform$path.sep), NOT_CRAN = "true")
options(multischolar.workflow_launcher = "process")
suppressPackageStartupMessages(library(MultiScholaR))
ns <- asNamespace("MultiScholaR")
cat(sprintf("KIT-131 SETUP: MultiScholaR %s from %s; R %s; %s\n", utils::packageVersion("MultiScholaR"),
    getNamespaceInfo(ns, "path"), getRversion(), Sys.info()[["sysname"]]))
suppressMessages(library(shiny))
env <- new.env(parent = ns)
for (helper in c("helper-workflow-omics-residency.R", "helper-e2e-enrichment-doubles.R", "helper-gui-entry.R"))
    sys.source(file.path("tests", "testthat", helper), env)

app_root <- file.path(out_dir, "application")
dir.create(app_root, showWarnings = FALSE)
app <- ns$WorkflowExecutionApplication$new()
app$config <- ns$workflowProcessApplicationConfig(app_root, FALSE, 0L)
stopifnot(identical(app$config$launcher, "process"))
app$pid <- Sys.getpid()
stopifnot(ns$checkWorkflowApplication(app))
cat(sprintf("KIT-131 SETUP: plain application, launcher %s, controller pid %d, root %s\n", app$config$launcher, app$pid, app_root))

driver <- new.env()
driver$waits <- 0L
driver$worker_seconds <- 0
driver$pids <- integer()

# Diagnostic (KIT-131, the silent loss on Windows workers): every worker the plain launcher starts is recorded (unit,
# request, operation, its TMPDIR); before the first worker of each lipid norm-tab publication (qc_figures, composite)
# and before the first worker at all, a probe (worker-probe.R) is started by the same launcher function with the same
# request, so it has that worker's environment, TMPDIR and working directory.
driver$launches <- list()
driver$probing <- FALSE
driver$probed <- character()
Sys.setenv(KIT131_PROBE = file.path(checkout, "packaging", "verify", "worker-probe.R"), KIT131_CHECKOUT = checkout,
    KIT131_OUT = out_dir, KIT131_SHORT = tempdir())
.kit131Launch <- function(config, unit, path) {
    if (isTRUE(driver$probing)) return(invisible())
    operation <- tryCatch(readRDS(path)$operation, error = function(error) "?")
    scratch <- file.path(dirname(path), "worker-scratch")
    driver$launches[[length(driver$launches) + 1L]] <- list(unit = unit, path = path, operation = operation, scratch = scratch)
    cat(sprintf("KIT-131 LAUNCH %d %s: unit %s; TMPDIR (worker-scratch) %d characters: %s\n", length(driver$launches), operation,
        unit, nchar(scratch), gsub("\\\\", "/", scratch)))
    kind <- if (grepl("_qc_figures$", operation)) "qc_figures" else if (grepl("_composite_figure$", operation)) "composite" else
        if (!length(driver$probed)) "first" else ""
    if (!nzchar(kind) || kind %in% driver$probed) return(invisible())
    driver$probed <- c(driver$probed, kind)
    driver$probing <- TRUE
    on.exit(driver$probing <- FALSE)
    probe <- paste0("kit131-probe-", length(driver$probed))
    result <- tryCatch({
        ns$workflowLaunchProcess(config, probe, "worker", "source(Sys.getenv(\"KIT131_PROBE\"))", path, "probe", 600L)
        ns$.workflowProcessUnits[[probe]]$process$wait(600000)
        status <- ns$.workflowProcessUnits[[probe]]$process$get_exit_status()
        ns$workflowProcessForget(probe)
        paste("exit", status)
    }, error = function(error) paste("error:", conditionMessage(error)))
    log <- file.path(config$root, paste0("worker-", probe, ".log"))
    cat(sprintf("KIT-131 PROBE %s for %s (%s); its log:\n", probe, operation, result))
    if (file.exists(log)) for (line in readLines(log, warn = FALSE))
        cat("KIT-131 PROBE |", substr(gsub("[\r\t]+", " ", line), 1L, 700L), "\n")
    invisible()
}
suppressMessages(trace("workflowLaunchProcess", tracer = quote(.kit131Launch(config, unit, path)), where = ns, print = FALSE))
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
# The R processes that are children of this controller while a managed step runs (the plain workers).
workers <- function() {
    children <- tryCatch(ps::ps_children(ps::ps_handle(), recursive = TRUE), error = function(error) list())
    vapply(children, function(p) tryCatch(ps::ps_pid(p), error = function(error) NA_integer_), integer(1))
}
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
        driver$pids <- unique(c(driver$pids, workers()))
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
# Outside test_that, as the Linux driver ran it; an expectation inside a step is counted and muffled, not fatal.
expectations <- list()
driver$error <- tryCatch(withCallingHandlers({
    driver$record <- env$guiEntryRunLane(lane_id)
    NULL
}, expectation = function(x) {
    expectations[[length(expectations) + 1L]] <<- x
    invokeRestart("muffle_expectation")
}), error = function(error) error)
closed <- tryCatch({ ns$closeWorkflowManagedApplication(app); "closed" }, error = function(error) paste("close failed:", conditionMessage(error)))
cat(sprintf("KIT-131 APPLICATION %s; worker processes still running under the controller: %d\n", closed, length(workers())))
# Diagnostic: each worker's log; in full for the norm tab's publications (qc_figures, composite) and the norm steps,
# else its warnings and errors only; the request's publication folder and result where they are still there.
for (i in seq_along(driver$launches)) {
    launch <- driver$launches[[i]]
    log <- file.path(app$config$root, paste0("worker-", launch$unit, ".log"))
    lines <- if (file.exists(log)) readLines(log, warn = FALSE) else character()
    result <- tryCatch(readRDS(file.path(dirname(launch$path), "worker-result.rds")), error = function(error) NULL)
    published <- list.files(file.path(dirname(launch$path), "p"), recursive = TRUE, full.names = TRUE)
    cat(sprintf("KIT-131 WORKER %d %s: log %d lines; result %s%s; request publication folder %s\n", i, launch$operation, length(lines),
        if (is.null(result)) "not there" else result$status, if (!is.null(result$error)) paste0(" (", substr(result$error, 1L, 400L), ")") else "",
        if (length(published)) paste(length(published), "files, longest", max(nchar(published))) else "empty or gone"))
    full <- grepl("qc_figures$|composite_figure$|_norm", launch$operation)
    shown <- if (full) lines else grep("WARN|ERROR|rror|arning|cannot|failed", lines, value = TRUE)
    for (line in utils::head(shown, if (full) 120L else 12L))
        cat(sprintf("KIT-131 WLOG %d %s | %s\n", i, launch$operation, substr(gsub("[\r\t]+", " ", line), 1L, 600L)))
}
if (!is.null(driver$error)) cat("KIT-131 PROBLEM run error:", substr(gsub("[\r\n\t]+", " ", conditionMessage(driver$error)), 1L, 900L), "\n")
classes <- vapply(expectations, function(x) class(x)[[1L]], character(1))
cat(sprintf("KIT-131 LANE %s expectations: pass %d fail %d error %d skip %d\n", lane_id, sum(classes == "expectation_success"),
    sum(classes == "expectation_failure"), sum(classes == "expectation_error"), sum(classes == "expectation_skip")))
for (x in utils::head(expectations[!classes %in% c("expectation_success", "expectation_skip")], 6L))
    cat("KIT-131 PROBLEM", class(x)[[1L]], ":", substr(gsub("[\r\n\t]+", " ", conditionMessage(x)), 1L, 700L), "\n")
if (is.null(driver$record)) {
    cat(sprintf("KIT-131 LANE %s: did not run to its end on workers (managed waits %d)\n", lane_id, driver$waits))
    quit(save = "no", status = 2L)
}
steps <- env$guiEntryComparable(driver$record)
env$guiEntryWriteJson(list(lane = lane_id, steps = steps), file.path(out_dir, paste0(lane_id, "-worker.json")))

working <- paste0("^Applying ", driver$omic, " step \\.\\.\\.$")
strip <- function(step) {
    notes <- step$notifications
    keep <- vapply(notes, function(note) !(identical(note$type, "default") && grepl(working, note$message)), logical(1))
    step$notifications <- notes[keep]
    list(step = step, working = sum(!keep))
}
short <- function(x, n = 160L) {
    x <- gsub("[\t\r\n]+", " ", paste(as.character(unlist(x)), collapse = " "))
    if (nchar(x) > n) paste0(substr(x, 1L, n), "...") else x
}
# Each step of ours against theirs: identical raw, identical without the working notification, else the differing fields
# and, for named lists, the differing keys (at most 6 per field).
compare <- function(label, theirs) {
    raw <- 0L; stripped_same <- 0L; working_total <- 0L
    for (id in union(names(theirs), names(steps))) {
        mine <- steps[[id]]
        other <- theirs[[id]]
        if (is.null(mine) || is.null(other)) {
            cat(sprintf("KIT-131 DIFF %s %s step %s: %s\n", label, lane_id, id, if (is.null(mine)) "absent on workers" else "absent there"))
            next
        }
        if (identical(mine, other)) raw <- raw + 1L
        s <- strip(mine)
        working_total <- working_total + s$working
        if (identical(s$step, other)) { stripped_same <- stripped_same + 1L; next }
        # Diagnostic: the norm steps' difference in full (every key, longer values), one line per key.
        whole <- label == "in-process" && id %in% c("norm_refusal", "norm")
        for (field in union(names(s$step), names(other))) {
            a <- s$step[[field]]; b <- other[[field]]
            if (identical(a, b)) next
            keyed <- is.list(a) && is.list(b) && !is.null(names(a)) && !is.null(names(b))
            keys <- if (keyed) Filter(function(k) !identical(a[[k]], b[[k]]), union(names(a), names(b))) else character()
            if (whole && length(keys)) {
                cat(sprintf("KIT-131 DIFF %s %s step %s field %s: %d keys\n", label, lane_id, id, field, length(keys)))
                for (k in keys) cat(sprintf("KIT-131 FULLDIFF %s field %s key %s || workers: %s || in-process: %s\n", id, field, k,
                    short(a[[k]], 650L), short(b[[k]], 650L)))
                next
            }
            if (whole) {
                cat(sprintf("KIT-131 FULLDIFF %s field %s || workers: %s || in-process: %s\n", id, field, short(a, 650L), short(b, 650L)))
                next
            }
            cat(sprintf("KIT-131 DIFF %s %s step %s field %s%s\n", label, lane_id, id, field,
                if (length(keys)) paste0(": ", length(keys), " keys: ", paste(vapply(utils::head(keys, 6L), function(k)
                    sprintf("%s [workers %s | %s %s]", short(k, 120L), short(a[[k]], 60L), label, short(b[[k]], 60L)), ""), collapse = "; "))
                else sprintf(" [workers %s | %s %s]", short(a, 200L), label, short(b, 200L))))
        }
    }
    cat(sprintf("KIT-131 LANE %s against %s: steps %d (there %d, same order %s); identical %d raw, %d without the working notification (%d removed)\n",
        lane_id, label, length(steps), length(theirs), identical(names(steps), names(theirs)), raw, stripped_same, working_total))
}
oracle <- jsonlite::read_json(env$guiEntryOraclePath(lane_id), simplifyVector = FALSE)
compare("oracle", oracle$steps)
if (nzchar(inprocess_path) && file.exists(inprocess_path)) {
    inprocess <- jsonlite::read_json(inprocess_path, simplifyVector = FALSE)
    compare("in-process", inprocess$steps)
} else {
    cat(sprintf("KIT-131 LANE %s against in-process: no in-process record at '%s'\n", lane_id, inprocess_path))
}
files <- list.files(app_root, recursive = TRUE)
cat(sprintf("KIT-131 LANE %s: managed waits %d, %.1f s waiting; distinct worker processes seen %d; %d files under the application root (receipts or journals %d); lane %.1f s\n",
    lane_id, driver$waits, driver$worker_seconds, length(driver$pids), length(files), sum(grepl("receipt|journal", files)),
    as.numeric(difftime(Sys.time(), started, units = "secs"))))
quit(save = "no", status = 0L)
