# KIT-252 (KIT-233 finding A3): a managed step that ends without an error of its own tells the user in plain words what
# happened, whether its result was saved and what to do, instead of an internal code ("managed action ended:
# memory_guard"); the code stays in the log by an explicit line. Every sentence is a DRAFT for the owner's wording pass
# (VAL-41, row QC-8; docs/qa/workflow-execution/kit/kit-252/wording.md). Written before the fix: on the tests-only
# commit the E rows fail at the base value "managed action ended: <status>" (the code), and G1 to G3 pass on both
# commits. The E rows drive the step's poll observer (`installWorkflowManagedStep`, R/utils_workflow_execution.R) with
# the poll's receipt mocked; E7, G2 the proteomics normalisation tab; G3 the lipidomics normalisation chain.

# ------------------------------------------------------------------------------------------------
# The drafts (ticket section 7), each after "Error applying <step>: ".
# ------------------------------------------------------------------------------------------------

kit252Text <- list(
    memory_limit = paste("the step needed more memory than the application allows for one step (at most 16 GB, less",
        "when the computer has less free memory), so it was stopped before it finished. Its result was not saved.",
        "Closing other programs may help; if it stops again, this data set is too large for this step on this computer."),
    low_memory = paste("the computer's free memory fell below 2 GB, so the step was stopped to protect the computer.",
        "Its result was not saved. Close other programs and press the button again."),
    time_limit = paste("the step reached its time limit of 30 minutes and was stopped. Its result was not saved.",
        "Running it again on the same data will stop at the same limit."),
    stopped = paste("the step's process stopped before it finished (for example, the computer ran out of memory or the",
        "process was closed). Its result was not saved. Press the button again; if it stops again, the application's log",
        "has the details."),
    not_started = paste("the step could not be started (too little free memory, or its process could not start).",
        "Nothing was run. Close other programs and press the button again."),
    publication_pending = paste("the result was saved, but its tables and figures could not be written to the project",
        "folder. Press the button again to write them; the step is not run again."),
    changed = paste("the table or the project changed while the step ran, so its result was not used.",
        "Press the button again."),
    internal = paste("the step's result could not be received (an internal check failed). Its result was not saved.",
        "Press the button again; if it happens again, the application's log file has the details."),
    not_closed = paste("the step's process could not be closed cleanly. Its result may have been saved; check this tab",
        "before running the step again, and restart the application if this happens again.")
)

# ------------------------------------------------------------------------------------------------
# The harness.
# ------------------------------------------------------------------------------------------------

#' One managed step (label "kit252", operation `operation`) on a stand-in project, pressed once, whose poll returns
#' `receipt` (the form of test-workflow-execution.R's KIT-115 case). The controller carries the server's policy, so a
#' step without a limit of its own has 1,800 s. Returns the error notices shown, the messages the step's `failed` hears
#' and the error log lines.
kit252Press <- function(receipt, operation = "protein_replicate_filter") {
    seen <- new.env(parent = emptyenv())
    seen$notes <- character()
    seen$errors <- character()
    seen$log <- character()
    manager <- structure(new.env(), class = c("DIAStateAdapter", "environment"))
    manager$.context <- list(store = list(execution_owner = NULL))
    workflow_data <- shiny::reactiveValues(state_manager = manager)
    ui <- list(label = "kit252", notification = "kit252_working",
        failed = function(error) seen$errors <- c(seen$errors, conditionMessage(error)))
    testthat::local_mocked_bindings(
        workflowSessionApplication = function(session) list(),
        beginWorkflowGUIAction = function(...) {
            reservation <- new.env()
            reservation$closed <- FALSE
            reservation
        },
        finishWorkflowGUIAction = function(reservation, status, ...) reservation$closed <- TRUE,
        usesWorkflowReferenceStates = function(manager) TRUE,
        newWorkflowExecutionController = function(manager, application) list(manager = manager, job = NULL,
            application = list(config = list(policy = workflowExecutionPolicy()))),
        runWorkflowAction = function(...) "wfx-kit252",
        workflowManagedCompletionContext = function(...) NULL,
        getWorkflowStateReference = function(manager) "reference",
        pollWorkflowAction = function(owner) receipt,
        .package = "MultiScholaR")
    testthat::local_mocked_bindings(
        showNotification = function(ui, action = NULL, duration = 5, closeButton = TRUE, id = NULL,
            type = c("default", "message", "warning", "error"), session = NULL) {
            if (identical(type[[1L]], "error")) seen$notes <- c(seen$notes, paste(as.character(ui), collapse = " "))
            invisible(id %||% "notification")
        },
        removeNotification = function(id, session = NULL) invisible(id),
        .package = "shiny")
    testthat::local_mocked_bindings(
        log_error = function(...) {
            seen$log <- c(seen$log, paste0(...))
            invisible(NULL)
        },
        .package = "logger")
    holder <- new.env(parent = emptyenv())
    shiny::testServer(function(input, output, session) {
        holder$step <- installWorkflowManagedStep(session, input, workflow_data, operation, NULL, ui, function(result) NULL)
    }, {
        shiny::isolate(holder$step$apply(list(a = 1)))
        session$flushReact()
    })
    seen
}

#' The notice the user sees for `receipt`, and what the step's `failed` hears (the same text without the prefix).
kit252Expect <- function(receipt, text) {
    seen <- kit252Press(receipt)
    expect_identical(seen$notes, paste("Error applying kit252:", text))
    expect_identical(seen$errors, text)
    invisible(seen)
}

#' Copy of test-prot-norm-core.R's `protNormAuditManaged` (:1722-1730): the proteomics normalisation tab on the managed
#' chain mocks, its notices in `log$notes`.
kit252ProtNormManaged <- function(fixture, probe, drive) {
    log <- new.env(parent = emptyenv())
    protNormW3Notices(log, fixture$root)
    protNormW3Managed(fixture, probe)
    protNormCoreQuiet(utils::capture.output(shiny::testServer(mod_prot_norm_server, args = list(
        workflow_data = fixture$workflow, experiment_paths = fixture$paths, omic_type = "proteomics",
        experiment_label = "core"), drive(session, session$getReturned(), log, output))))
    invisible(log)
}

#' The proteomics normalisation's first action ends with `receipt` (the form of test-prot-norm-core.R:1839-1858).
#' Returns the notices, the error log lines, the dispatched operation and the head afterwards.
kit252ProtNormPress <- function(receipt) {
    fixture <- protNormCoreOwner()
    seen <- new.env(parent = emptyenv())
    seen$log <- character()
    testthat::local_mocked_bindings(
        log_error = function(...) {
            seen$log <- c(seen$log, paste0(...))
            invisible(NULL)
        },
        .package = "logger")
    kit252ProtNormManaged(fixture, new.env(), function(session, norm, log, output) {
        do.call(session$setInputs, protNormW3Inputs("manual"))
        testthat::with_mocked_bindings({
            session$setInputs(run_normalization = 1)
            session$flushReact()
        }, runWorkflowAction = function(operation, ...) {
            seen$operation <- operation
            "wfx-not-run"
        }, pollWorkflowAction = function(owner) receipt, .package = "MultiScholaR")
        seen$notes <- log$notes
    })
    seen$head <- fixture$manager$current_state
    seen
}

#' The lipidomics normalisation chain on the managed step (the KIT-078 tab server and chain mocks): QC as the KIT-078
#' journey, then "Run Normalization Pipeline", whose RUV action ends with `memory_guard` before it delivers (the earlier
#' actions, ITSD, log2 and between-sample normalisation, are committed as they are delivered). Returns the notices, the
#' head and the history.
kit252LipidChain <- function(root) {
    fixture <- lipidDeclarationsOwner(file.path(root, "lipid"))
    withr::defer(closeDIAStateAdapter(fixture$owner))
    binding <- omicsDeclarationBindings()[["lipidomics"]]
    ruv <- paste0(binding$prefix, "_ruv")
    probe <- new.env()
    log <- new.env(parent = emptyenv())
    protNormW3Notices(log, root)
    omicsManagedChainMocks(fixture, probe)
    chain_run <- get("runWorkflowAction", envir = asNamespace("MultiScholaR"))
    testthat::local_mocked_bindings(
        runWorkflowAction = function(operation, ...) {
            probe$last <- operation
            # The RUV worker is stopped before it delivers: nothing of it is accepted.
            if (identical(operation, ruv)) return(paste0("worker-", operation))
            chain_run(operation, ...)
        },
        pollWorkflowAction = function(owner) list(status = if (identical(probe$last, ruv)) "memory_guard" else "complete",
            error = NULL),
        .package = "MultiScholaR")
    assays <- names(methods::slot(fixture$object, binding$data_slot))
    id_column <- methods::slot(fixture$object, binding$id_slot)
    ids <- stats::setNames(lapply(assays, function(assay) methods::slot(fixture$object, binding$data_slot)[[assay]][[id_column]]), assays)
    selections <- list(assays = assays, candidates = ids, rows = stats::setNames(lapply(assays, function(assay) 1:2), assays))
    suppressWarnings(suppressMessages(utils::capture.output(shiny::testServer(
        omicsTabServer(fixture, "lipidomics", probe, TRUE), {
            do.call(session$setInputs, omicsTabInputs())
            session$setInputs(`intensity-intensity_cutoff_percentile` = 1, `intensity-proportion_below_cutoff` = 0.5,
                `intensity-apply_filter` = 1)
            session$flushReact()
            session$setInputs(`duplicates-resolve_duplicates` = 1)
            session$flushReact()
            commitWorkflowOmicsStoreStep(fixture$owner, "lipid_qc_finalize", "lipid_qc_complete", list(),
                "QC processing complete - ready for normalization")
            norm <- session$userData$probe_norm
            norm$itsd_selections <- selections$rows
            norm$itsd_candidate_ids <- selections$candidates
            norm$assay_names <- selections$assays
            session$setInputs(run_normalization = 1)
            session$flushReact()
        }))))
    list(notes = log$notes, head = fixture$owner$current_state, history = unlist(fixture$owner$getHistory()),
        operations = probe$operations, prefix = binding$prefix)
}

# ------------------------------------------------------------------------------------------------
# E rows: each status that reaches the user through the poll, and an unknown one.
# ------------------------------------------------------------------------------------------------

test_that("KIT-252 E1a: memory_guard (the application's memory limit for one step) is told as the memory limit", {
    seen <- kit252Expect(list(status = "memory_guard", error = NULL), kit252Text$memory_limit)
    expect_false(any(grepl("memory_guard", seen$notes, fixed = TRUE)))
})

test_that("KIT-252 E1b: host_headroom (the computer's free memory below 2 GB) is told as low memory", {
    kit252Expect(list(status = "host_headroom", error = NULL), kit252Text$low_memory)
})

test_that("KIT-252 E2a: a worker stopped at its time limit (1,800 s) is told as the 30-minute time limit", {
    for (error in list(NULL, "")) {
        kit252Expect(list(status = "worker_failed", error = error, worker = list(Result = "timeout")), kit252Text$time_limit)
    }
})

test_that("KIT-252 E2b: a worker that stopped otherwise, with no error of its own, is told as stopped", {
    for (worker in list(list(Result = "signal"), list(Result = "exit-code"), NULL)) {
        for (error in list(NULL, "")) {
            kit252Expect(list(status = "worker_failed", error = error, worker = worker), kit252Text$stopped)
        }
    }
})

test_that("KIT-252 E3: admission_refused (too little memory to start, or a launch failure) and monitor_start_timeout are told as not started", {
    kit252Expect(list(status = "admission_refused", error = NULL), kit252Text$not_started)
    kit252Expect(list(status = "admission_refused", error = "the worker unit could not be started"), kit252Text$not_started)
    kit252Expect(list(status = "monitor_start_timeout", error = NULL), kit252Text$not_started)
})

test_that("KIT-252 E4: stale and project_replaced are told as changed meanwhile", {
    for (status in c("stale", "project_replaced")) kit252Expect(list(status = status, error = NULL), kit252Text$changed)
})

test_that("KIT-252 E5: monitor_lost, invalid_result, idle, the six guard reasons and an unknown status are told as internal", {
    statuses <- c("monitor_lost", "invalid_result", "idle", "nonworker_reservation", "invalid_worker_registration",
        "invalid_monitor_snapshot", "invalid_guard_injection", "injected_guard_stop", "guard_payload_invalid", "zzz_unknown")
    for (status in statuses) kit252Expect(list(status = status, error = NULL), kit252Text$internal)
})

test_that("KIT-252 E5b: cleanup_failed and duplicate (either can follow a saved result) never say the result was not saved", {
    for (status in c("cleanup_failed", "duplicate")) {
        seen <- kit252Expect(list(status = status, error = NULL), kit252Text$not_closed)
        expect_false(any(grepl("not saved", seen$notes, fixed = TRUE)))
    }
})

test_that("KIT-252 E6: publication_pending (the result committed, its files not written) says it was saved", {
    seen <- kit252Expect(list(status = "publication_pending", error = NULL), kit252Text$publication_pending)
    expect_false(any(grepl("not saved", seen$notes, fixed = TRUE)))
})

test_that("KIT-252 E7: through the proteomics normalisation tab, memory_guard is the memory sentence and its code is in the log", {
    skip_if_not(workflowOwnerCapable())
    seen <- kit252ProtNormPress(list(status = "memory_guard", error = NULL))
    expect_identical(grep("^\\[error", seen$notes, value = TRUE), c(
        paste("[error 15] Error applying proteomics step:", kit252Text$memory_limit),
        paste("[error 10] Error in normalization:", kit252Text$memory_limit)))
    expect_true(sprintf("managed action %s ended: memory_guard", seen$operation) %in% seen$log,
        info = paste(seen$log, collapse = " | "))
})

# ------------------------------------------------------------------------------------------------
# G rows: what must not change.
# ------------------------------------------------------------------------------------------------

test_that("KIT-252 G1: a worker that failed with its own error still shows that error (GUI audit A111-03)", {
    kit252Expect(list(status = "worker_failed", error = "No valid percentage found", worker = list(Result = "exit-code")),
        "No valid percentage found")
})

test_that("KIT-252 G2: a single action's failure leaves the head where it was", {
    skip_if_not(workflowOwnerCapable())
    seen <- kit252ProtNormPress(list(status = "memory_guard", error = NULL))
    expect_identical(seen$head, "protein_replicate_filtered")
})

test_that("KIT-252 G3: a lipidomics normalisation whose RUV action ends memory_guard keeps the normalised head; no notice says 'Nothing was saved'", {
    skip_if_not(workflowOwnerCapable())
    chain <- kit252LipidChain(withr::local_tempdir())
    prefix <- chain$prefix
    expect_true(paste0(prefix, "_ruv") %in% chain$operations, info = paste(chain$operations, collapse = " "))
    expect_identical(chain$head, paste0(prefix, "_normalized"))
    expect_true(all(paste0(prefix, c("_log2", "_normalized")) %in% chain$history))
    expect_false(paste0(prefix, "_ruv_corrected") %in% chain$history)
    expect_length(grep("Nothing was saved", chain$notes, fixed = TRUE), 0L)
    expect_length(grep("^\\[error 15\\] Error applying lipidomics step: ", chain$notes), 1L)
})

test_that("KIT-252 E1a in a chain: the lipidomics RUV action's memory_guard is told as the memory limit", {
    skip_if_not(workflowOwnerCapable())
    chain <- kit252LipidChain(withr::local_tempdir())
    expect_identical(grep("^\\[error 15\\]", chain$notes, value = TRUE),
        paste("[error 15] Error applying lipidomics step:", kit252Text$memory_limit))
})
