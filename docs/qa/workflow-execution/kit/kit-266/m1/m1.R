# KIT-266 M1 (controller memory of the worker branch). Usage (run-tests scripts job, package root):
#   Rscript docs/qa/workflow-execution/kit/kit-266/m1/m1.R
# The limpa lane is driven on the managed route with the test double, except that the worker role runs in a CHILD Rscript
# (`MultiScholaR:::workflowWorkerMain`, MULTISCHOLAR_WORKFLOW_WORKER=1, the production launcher's command shape), so this
# process is the controller only. VmHWM is reset ("5" to /proc/self/clear_refs) right before each limpa_composite action's
# child starts; the rise is VmHWM after that action's delivery minus VmRSS before it, read from this process. Controls:
# the reset (VmHWM falls to VmRSS); the path (worker path taken: files written from the store, no fallback, PNG size);
# the positive control (G0's in-process build in a fresh Rscript, memory-build.R, must read 0.3 GiB or more).
out <- normalizePath("kit266-m1", mustWork = FALSE)
dir.create(out, recursive = TRUE, showWarnings = FALSE)
say <- function(...) {
    line <- paste0(..., collapse = "")
    cat(line, "\n", sep = "")
    cat(line, "\n", sep = "", file = file.path(out, "m1.txt"), append = TRUE)
}
status <- function() {
    lines <- readLines("/proc/self/status")
    lines[grepl("^(VmPeak|VmSize|VmHWM|VmRSS):", lines)]
}
kb <- function(field) as.numeric(gsub("[^0-9]", "", status()[startsWith(status(), paste0(field, ":"))]))
GiB <- 1024^2
suppressMessages(pkgload::load_all(".", quiet = TRUE, export_all = TRUE, helpers = TRUE))
rows <- list()
path <- new.env(parent = emptyenv())
path$written <- logical()
path$logged <- character()
lines <- function(...) say(...)

testthat::test_that("KIT-266 M1", {
    testthat::skip_if_not(workflowOwnerCapable())
    fixture <- protNormW1LimpaOwner()
    log <- new.env(parent = emptyenv())
    protNormW3Notices(log, fixture$root)
    probe <- new.env()
    controller <- omicsManagedChainMocks(list(owner = fixture$manager, paths = fixture$paths), probe)
    child <- function(request_path) {
        hash <- digest::digest(file = request_path, algo = "sha256")
        expression <- paste0("ns <- asNamespace('MultiScholaR'); unlockBinding('workflowExecutionBuildIdentity', ns); ",
            "assign('workflowExecutionBuildIdentity', function() strrep('f', 64), ns); MultiScholaR:::workflowWorkerMain(commandArgs(TRUE))")
        system2(file.path(R.home("bin"), "Rscript"), c("--vanilla", "-e", shQuote(expression), request_path, hash),
            env = c("MULTISCHOLAR_WORKFLOW_WORKER=1", paste0("R_LIBS=", paste(.libPaths(), collapse = .Platform$path.sep))),
            stdout = file.path(out, "child.log"), stderr = file.path(out, "child.log"))
    }
    # protNormW2Dispatch with the worker role in a child process.
    dispatch <- function(operation, parameters, input, context = NULL) {
        store <- fixture$manager$.context$store
        id <- paste0("wfx-", substr(digest::digest(list(tempfile(), operation), algo = "sha256"), 1L, 32L))
        paths <- workflowExportPaths(store$directory, id)
        token <- list(request_id = id, incarnation = store$identity$incarnation, expected_head = workflowExecutionHead(store)$entity_id,
            generation = controller$control$generation, controller_id = controller$control$controller_id)
        writeWorkflowExecutionMessage(file.path(paths$request_dir, "publication-context.rds"),
            list(experiment_label = "core", experiment_paths = fixture$paths, omic_type = "proteomics"))
        writeWorkflowExportJournal(paths, list(schema = "multischolar.workflow_action.v1", request_id = id, status = "admitting",
            expected_head = token$expected_head, token_hash = digest::digest(token, algo = "sha256")))
        request <- list(schema = 1L, build = workflowExecutionBuildIdentity(), operation = operation,
            input_references = list(state = input), parameters = parameters, token = token, owner = store$identity,
            project_root = dirname(store$directory), directory = paths$request_dir, config = list(mode = "adapters"),
            worker_cap = 12 * 1024^3)
        request_path <- file.path(paths$request_dir, "request.rds")
        writeWorkflowExecutionMessage(request_path, request)
        composite <- identical(parameters$figures, "limpa_composite")
        if (composite) {
            writeLines("5", "/proc/self/clear_refs")
            path$rss <- kb("VmRSS")
            path$hwm_after_reset <- kb("VmHWM")
            path$mode <- parameters$ruv_mode
        }
        child(request_path)
        message <- readRDS(file.path(paths$request_dir, "worker-result.rds"))
        if (!identical(message$status, "complete")) stop(message$error, call. = FALSE)
        controller$job <- list(request = list(token = token, operation = operation, config = list(mode = "adapters"),
            parameters = parameters), paths = paths, completion_context = validateWorkflowOmicsCompletionContext(operation, context))
        receipt <- acceptWorkflowActionResult(controller, message[c("schema", "status", "token", "result")])
        accepted <- resumeWorkflowAction(controller, id)
        c(readManagedOmicsResult(fixture$manager, accepted, NULL), list(acceptance = receipt$status, result = message$result))
    }
    original_deliver <- deliverDiaLimpaCompositeFigure
    original_write <- writeDIACompositeFigureFiles
    old <- logger::log_appender(function(lines) path$logged <- c(path$logged, lines))
    withr::defer(logger::log_appender(old))
    testthat::local_mocked_bindings(
        runWorkflowAction = function(operation, input_references, parameters, owner, reservation, completion_context = NULL) {
            probe$actions <- c(probe$actions, list(list(operation = operation, input = input_references$state$entity_id)))
            out_value <- dispatch(operation, parameters, input_references$state, completion_context)
            probe$handoff <- list(accepted = list(operation = operation, reference = out_value$reference,
                state_reference = out_value$stateReference, display = out_value$display))
            paste0("worker-", operation)
        },
        pollWorkflowAction = function(owner) list(status = "complete"),
        writeDIACompositeFigureFiles = function(...) {
            got <- original_write(...)
            path$written <- c(path$written, got)
            got
        },
        deliverDiaLimpaCompositeFigure = function(...) {
            result <- original_deliver(...)
            hwm <- kb("VmHWM")
            rows[[length(rows) + 1L]] <<- list(mode = path$mode, rss_before = path$rss, hwm_after_reset = path$hwm_after_reset,
                hwm_after = hwm, rise = (hwm - path$rss) / GiB)
            result
        }, .package = "MultiScholaR")
    click <- protNormW3Clicker()
    steps <- list(pass_through = do.call(click, c(list("run_normalization"), protNormW3Inputs("covariates"))),
        correlation = click("apply_correlation_filter"), ruv = click("apply_ruv_covariates"))
    selected <- shiny::reactiveVal("import")
    shiny::testServer(mod_prot_norm_server, args = list(workflow_data = fixture$workflow, experiment_paths = fixture$paths,
        omic_type = "proteomics", experiment_label = "core", selected_tab = selected), {
        session$userData$dia_capacity_manager <- fixture$manager
        for (case in names(steps)) {
            steps[[case]](session)
            for (i in 1:3) {
                session$elapse(100)
                session$flushReact()
            }
            sizes <- file.path(fixture$paths$protein_qc_dir, "proteomics_composite_QC_figure.png")
            png <- if (file.exists(sizes)) readBin(readBin(sizes, "raw", 24L)[17:24], "integer", n = 2L, size = 4L, endian = "big") else c(NA, NA)
            say("STEP ", case, " composite PNG ", png[[1L]], " x ", png[[2L]], " px; actions ", length(probe$actions))
        }
    })
    dir.create(file.path(out, "images"), showWarnings = FALSE)
    file.copy(list.files(fixture$paths$protein_qc_dir, pattern = "[.]png$", full.names = TRUE), file.path(out, "images"))
})

say("CONTROL reset: after clear_refs VmHWM is within 5 percent of VmRSS: ",
    if (length(rows)) all(vapply(rows, function(r) r$hwm_after_reset <= 1.05 * r$rss_before, NA)) else NA)
for (r in rows) say(sprintf("M1_RISE mode=%s VmRSS_before_kB=%.0f VmHWM_after_reset_kB=%.0f VmHWM_after_kB=%.0f rise_GiB=%.3f verdict=%s",
    r$mode, r$rss_before, r$hwm_after_reset, r$hwm_after, r$rise, if (r$rise < 0.3) "under 0.3 GiB" else "AT OR OVER 0.3 GiB"))
say("CONTROL path: writeDIACompositeFigureFiles returned ", paste(path$written, collapse = " "), "; fallback logged: ",
    any(grepl("not found saved", path$logged, fixed = TRUE)))
say("rows measured: ", length(rows))
for (mode in c("two", "three")) {
    result <- system2(file.path(R.home("bin"), "Rscript"), c("docs/qa/workflow-execution/kit/kit-266/m1/memory-build.R", mode,
        file.path(out, "images")), stdout = TRUE, stderr = TRUE)
    writeLines(result, file.path(out, paste0("positive-control-", mode, ".txt")))
    say(grep("MEMORY_RISE", result, value = TRUE))
}
