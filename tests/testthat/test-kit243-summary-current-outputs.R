# KIT-243 (VAL-212 (1), (2); the coordinator's ruling on the plan audit's P02): each differential analysis first removes the
# figures and per-run files its publication writes by name that an earlier analysis left (Save Heatmap's folder kept);
# while that publication has not completed, Copy and the Report refuse (`workflow_da_publication_pending`); and they
# refuse while the last "Save Workflow Arguments" does not describe the current history and DA (`workflow_summary_save_stale`).
# Written before the fix: on the tests-only commit each E row runs as before, so it fails at the base value the ticket
# states (section 10); the G rows pass on both commits. The new functions are called only where they exist (head), after
# the base assertions. Lines `KIT-243 <row>` are the evidence the ticket's records copy (outcome.md).

library(testthat)

# ------------------------------------------------------------------------------------------------
# Copied unchanged from tests/testthat/test-kit212-da-current.R at v0.9-staging 16112917 (KIT-207 E-7's rule): lines 16-666
# (KIT-194's, KIT-207's and KIT-212's helpers). Then, from tests/testthat/test-workflow-saved-da.R at 16112917: lines
# 664-685 (savedDaW1Owner), 742-746 (savedDaW1Contrasts), 1209-1272 (savedDaW2Dispatch, savedDaW2Front, savedDaW2Managed),
# 1287-1299 (savedDaW2Dispatcher) and 1316-1320 (savedDaW2Mocks).
# ------------------------------------------------------------------------------------------------

k194Omic <- function(class) c(lipid = "lipidomics", metab = "metabolomics")[[class]]

#' Copy of `kit172PreNorm` (test-kit172-downstream.R:127-137; A11): an owner at the QC-complete state (intensity filter,
#' QC finalise), as the normalisation tab finds it.
k194PreNorm <- function(root, values, class) {
    fixture <- kit172Owner(root, kit172Object(values, class), class)
    binding <- kit172Binding(class)
    kit172Quiet({
        runWorkflowOmicsReferenceIntensityFilter(binding, fixture$owner, 1, 0.5, fixture$workflow$config_list)
        commitWorkflowOmicsStoreStep(fixture$owner, paste0(binding$prefix, "_qc_finalize"), paste0(binding$prefix, "_qc_complete"),
            fixture$workflow$config_list, "QC complete")
    })
    fixture
}

#' Copy of `kit172InProcessRunner` (test-kit172-downstream.R:140-160; A11): the managed actions run in this process, each
#' commit through its declaration's wrapper, its result delivered as the worker's; a product delivers an empty display.
k194InProcessRunner <- function(manager) {
    list(progress = function(...) list(inc = function(...) invisible(NULL), close = function() invisible(NULL)),
        run = function(actions, done = function() NULL, failed = function(error) stop(error)) {
            outcome <- tryCatch({
                for (action in actions) {
                    if (is.function(action$before)) action$before()
                    if (is.null(action$context)) {
                        action$deliver(list(display = list()))
                    } else {
                        step <- list(operation = action$operation, parameters = action$parameters, state_name = action$context$state_name,
                            description = action$context$description)
                        result <- runWorkflowOmicsReferenceStep(manager, step, action$context$config)
                        action$deliver(list(stateReference = result$reference, value = result$value))
                    }
                }
                NULL
            }, error = function(e) e)
            if (inherits(outcome, "error")) failed(outcome) else done()
        })
}

#' The recording runner of the managed cases (B04): `step$active()` is NULL while idle and an action while `busy`; a
#' press while `busy` or `hold` is queued, and `start()` runs what was queued; `actions` given as a function are resolved
#' when they start (`R/utils_workflow_omics_actions.R:43`), a refusal then reaching `failed`; the actions run as
#' `k194InProcessRunner` runs them, on the bus's owner of the moment. `operations` records every action run.
k194Runner <- function(workflow, busy = FALSE, hold = FALSE) {
    rec <- new.env(parent = emptyenv())
    rec$operations <- character()
    rec$queued <- list()
    rec$busy <- busy
    rec$hold <- hold
    execute <- function(actions, done, failed) {
        resolved <- tryCatch(if (is.function(actions)) actions() else actions, error = function(error) error)
        if (inherits(resolved, "error")) {
            if (is.function(failed)) failed(resolved) else stop(resolved)
            return(invisible(NULL))
        }
        rec$operations <- c(rec$operations, vapply(resolved, function(action) action$operation, ""))
        k194InProcessRunner(workflow$state_manager)$run(resolved, done = done %||% function() NULL,
            failed = failed %||% function(error) stop(error))
        invisible(NULL)
    }
    rec$start <- function() {
        queued <- rec$queued
        rec$queued <- list()
        rec$hold <- FALSE
        for (waiting in queued) execute(waiting$actions, waiting$done, waiting$failed)
        invisible(NULL)
    }
    runner <- new.env(parent = emptyenv())
    runner$step <- list(active = function() if (isTRUE(rec$busy)) "a running action")
    runner$progress <- function(...) list(set = function(...) invisible(NULL), inc = function(...) invisible(NULL),
        close = function(...) invisible(NULL))
    runner$run <- function(actions, done = NULL, failed = NULL, dropped = NULL) {
        if (isTRUE(rec$busy) || isTRUE(rec$hold)) {
            rec$queued[[length(rec$queued) + 1L]] <- list(actions = actions, done = done, failed = failed)
            return(invisible(NULL))
        }
        execute(actions, done, failed)
    }
    rec$runner <- runner
    rec
}

#' A lipid or metab project at QC complete on F172 (seed 172001), with the omic's GUI configuration (its QC figures
#' drawn by nothing) and the tab's log.
k194Fixture <- function(root, class) {
    fixture <- k194PreNorm(root, kit172Values(172001L)$values, class)
    fixture$class <- class
    fixture$omic <- k194Omic(class)
    fixture$binding <- kit172Binding(class)
    fixture$prefix <- fixture$binding$prefix
    config <- omicGuiConfig(fixture$omic)
    config$binding$binding$qc_plots <- function(...) invisible(NULL)
    fixture$config <- config
    log <- new.env(parent = emptyenv())
    log$lines <- character()
    fixture$log <- log
    fixture$addLog <- function(entry) log$lines <- c(log$lines, entry)
    fixture
}

k194Close <- function(fixture) try(closeDIAStateAdapter(fixture$workflow$state_manager), silent = TRUE)

#' Close the bus's owner of the moment (a load replaces it) when the calling test ends; `fixture` is forced here, so a
#' loop's fixtures are each closed.
k194Defer <- function(fixture, envir = parent.frame()) {
    force(fixture)
    withr::defer(k194Close(fixture), envir = envir)
}

k194Owner <- function(fixture) fixture$workflow$state_manager
k194History <- function(fixture) as.character(unlist(k194Owner(fixture)$getHistory()))
k194Head <- function(fixture) k194Owner(fixture)$current_state
k194Ref <- function(fixture, state = NULL) getWorkflowStateReference(k194Owner(fixture), state)

#' Shiny's notifications recorded as "[type] text"; removal, progress and selector updates quiet.
k194Notes <- function(env = parent.frame()) {
    seen <- new.env(parent = emptyenv())
    seen$notes <- character()
    testthat::local_mocked_bindings(
        showNotification = function(ui, ..., type = "default") {
            seen$notes <- c(seen$notes, paste0("[", type, "] ", paste(as.character(ui), collapse = " ")))
            invisible(NULL)
        },
        removeNotification = function(...) invisible(NULL),
        withProgress = function(expr, ...) expr,
        incProgress = function(...) invisible(NULL),
        updateSelectInput = function(...) invisible(NULL),
        .package = "shiny", .env = env)
    seen
}

k194Noted <- function(seen, text) any(grepl(text, seen$notes, fixed = TRUE))

#' The tab's inputs: ITSD off unless asked (by pattern), offset 0, `method`, RUV skip unless `ruv` (the legacy list).
k194Input <- function(method = "quantile", itsd = FALSE, ruv = NULL, threshold = 0.5) {
    list(apply_itsd = itsd, itsd_aggregation = "median", log_offset = 0, norm_method = method,
        ruv_mode = if (is.null(ruv)) "skip" else ruv$ruv_mode, auto_percentage_min = 1, auto_percentage_max = 20,
        ruv_grouping_variable = "group", separation_metric = "max_difference", k_penalty_weight = 0.5,
        adaptive_k_penalty = TRUE, ruv_k = ruv$manual_k %||% 1L, ruv_percentage = ruv$manual_percentage %||% 20,
        min_pearson_correlation_threshold = threshold)
}

#' "Run Normalization Pipeline": the frame's `runOmicNormPipeline`, in process (`runner` NULL) or managed.
k194Run <- function(fixture, nd, input = k194Input(), runner = NULL) {
    testthat::local_mocked_bindings(workflowOmicsManagedActions = function(omic_type) runner, .package = "MultiScholaR")
    kit172Quiet(shiny::isolate(runOmicNormPipeline(input, fixture$workflow, list(), nd, fixture$addLog,
        function(...) list(color_var = "group", shape_var = NULL), fixture$config)))
}

#' The run a case builds on completed (an instrument check, VAL-54): head the normalised or RUV state, flag set.
k194ExpectRun <- function(fixture, nd, info = "") {
    expect_true(k194Head(fixture) %in% paste0(fixture$prefix, c("_normalized", "_ruv_corrected")),
        info = paste("run precondition:", info, k194Head(fixture), paste(fixture$log$lines, collapse = " | ")))
    expect_true(isTRUE(shiny::isolate(nd$normalization_complete)), info = paste("run precondition:", info))
}

k194Apply <- function(fixture, nd, runner = NULL, threshold = 0.5) {
    testthat::local_mocked_bindings(workflowOmicsManagedActions = function(omic_type) runner, .package = "MultiScholaR")
    kit172Quiet(shiny::isolate(applyOmicNormCorrelation(k194Input(threshold = threshold), fixture$workflow, nd, fixture$addLog,
        fixture$config)))
}

k194Skip <- function(fixture, nd) kit172Quiet(shiny::isolate(skipOmicNormCorrelation(fixture$workflow, nd, fixture$addLog, fixture$config)))

k194Reset <- function(fixture, nd) kit172Quiet(shiny::isolate(resetOmicNorm(fixture$workflow, nd, fixture$addLog, fixture$config)))

#' Revert: the Intensity Filter sub-tab's real observer under `testServer`, pressed until the head is `until` (or
#' `times` presses).
k194Revert <- function(fixture, until = NULL, times = 10L) {
    kit172Quiet(utils::capture.output(shiny::testServer(mod_omic_qc_intensity_server, args = list(workflow_data = fixture$workflow,
        omic_type = fixture$omic, experiment_label = "kit194", config = fixture$config), {
        press <- 0L
        while (press < times && (is.null(until) || !identical(fixture$workflow$state_manager$current_state, until))) {
            press <- press + 1L
            session$setInputs(revert_filter = press)
        }
    })))
    k194Head(fixture)
}

#' A state's LCMS_Pos matrix (features by the design's runs present).
k194Matrix <- function(fixture, reference) {
    object <- readWorkflowOmicsState(fixture$binding, k194Owner(fixture), reference)
    table <- as.data.frame(methods::slot(object, fixture$binding$data_slot)$LCMS_Pos)
    id <- methods::slot(object, fixture$binding$id_slot)
    runs <- intersect(kit172Design()$Run, colnames(table))
    matrix(as.matrix(table[, runs, drop = FALSE]), ncol = length(runs), dimnames = list(table[[id]], runs))
}

k194Max <- function(fixture, reference) max(k194Matrix(fixture, reference), na.rm = TRUE)

#' A session file as an older build exported it (the fields `collectWorkflowOmicsNormSession` writes), holding `object`
#' under the state name `state`, in its own directory with the DA's file name.
k194LegacySession <- function(fixture, object, state, dir) {
    dir.create(dir, recursive = TRUE, showWarnings = FALSE)
    session <- list(r6_current_state_name = state, current_s4_object = object, contrasts_tbl = fixture$workflow$contrasts_tbl,
        design_matrix = object@design_matrix, config_list = list(), itsd_selections = list(), ruv_optimization_results = list(),
        correlation_results = NULL, assay_names = names(methods::slot(object, fixture$binding$data_slot)),
        export_timestamp = as.POSIXct("2026-10-28 12:00:00", tz = "UTC"), omic_type = fixture$omic, experiment_label = "legacy",
        normalization_method = "quantile", ruv_mode = "skip", itsd_applied = FALSE, itsd_aggregation = NA, log_offset = 0,
        correlation_threshold = 0.5, ruv_grouping_variable = "group", feature_counts = NULL, qc_params = list(),
        normalization_complete = TRUE, ruv_complete = FALSE, correlation_filtering_complete = TRUE)
    saveRDS(session, file.path(dir, fixture$config$da$session_load$file))
    dir
}

#' "Load Filtered Session" of the DA tab (the core's load, in the module's session, B04); the bus's owner is replaced by
#' the adopted one. Returns the DA tab's values.
k194Load <- function(fixture, dir) {
    tab <- new.env(parent = emptyenv())
    out <- kit172Quiet(shiny::withReactiveDomain(fixture$session, loadWorkflowOmicsDaSession(fixture$binding, fixture$workflow, tab,
        list(source_dir = dir), fixture$config$da$session_load)))
    expect_identical(out$status, "success", info = paste("load precondition:", out$message))
    tab
}

#' A legacy session of F172's table loaded under `state`.
k194LoadLegacy <- function(fixture, state, root) {
    # A self-contained object, as an older build's file holds one (a state read from the owner is reference-backed).
    object <- kit172Object(kit172Values(172001L)$values, fixture$class)
    k194Load(fixture, k194LegacySession(fixture, object, state, file.path(root, paste0("legacy-", state))))
}

#' trx's "Run Normalization Pipeline" (`runTrxNormalisation`, as test-trx-norm-core.R drives it) with the notifications
#' recorded in `seen$notes` and the tab's log in `nd`.
k194TrxRun <- function(project, nd, inputs, seen, runner = NULL) {
    testthat::local_mocked_bindings(workflowOmicsManagedActions = function(omic_type) runner, .package = "MultiScholaR")
    notify <- function(ui, ..., type = "default", duration = 5) {
        seen$notes <- c(seen$notes, paste0("[", type, "] ", paste(as.character(ui), collapse = " ")))
        invisible(NULL)
    }
    set_status <- function(tab, value) {
        status <- project$workflow$tab_status
        status[[tab]] <- value
        project$workflow$tab_status <- status
    }
    withr::local_seed(20261124L)
    kit172Quiet(utils::capture.output(shiny::isolate(runTrxNormalisation(list(managed = runTrxNormManaged, in_process = runTrxNormInProcess),
        transcriptomicsDeclarationBinding, project$workflow$state_manager, trxNormGridParameters(inputs), nd, project$workflow$config_list,
        project$paths, omicNormAddLog(nd), notify, set_status))))
    invisible(NULL)
}

#' A fresh project's directories for the omic (the Summary's `projectDirs`).
k207Dirs <- function(root, omic) {
    paths <- list(source_dir = file.path(root, "scripts", omic), integration_dir = file.path(root, "integration"), base_dir = root)
    dir.create(paths$source_dir, recursive = TRUE, showWarnings = FALSE)
    stats::setNames(list(paths), omic)
}

k207IntegrationDir <- function(paths) paths$integration_dir %||% file.path(paths$base_dir, "integration")

#' "Save Workflow Arguments": the core's `saveWorkflowSummaryParameters` under the omic's group, in process (`runner`
#' NULL) or managed, with the real study-parameter and record writers; the notifications ("[type] text"), `done`,
#' `failed`, the in-process RDS writes and the console lines recorded.
k207Save <- function(workflow, omic, dirs, runner = NULL) {
    rec <- new.env(parent = emptyenv())
    rec$notes <- character()
    rec$done <- 0L
    rec$failed <- list()
    rec$rds <- list()
    rec$cat <- character()
    testthat::local_mocked_bindings(workflowOmicsManagedActions = function(omic_type) runner, .package = "MultiScholaR")
    kit172Quiet(utils::capture.output(saveWorkflowSummaryParameters(omicGuiConfig(omic)$summary, workflow, dirs, omic, "k207", "",
        notify = function(message, type, duration) rec$notes <- c(rec$notes, paste0("[", type, "] ", message)),
        done = function(saved) rec$done <- rec$done + 1L,
        failed = function(error) rec$failed[[length(rec$failed) + 1L]] <- error,
        saveRDSFn = function(object, path) {
            rec$rds[[length(rec$rds) + 1L]] <- list(object = object, path = path)
            invisible(NULL)
        },
        catFn = function(...) rec$cat <- c(rec$cat, paste(..., collapse = " ")))))
    rec
}

#' The reference DA on the head (as `kit172Journey`); it commits no state.
k207Da <- function(fx) {
    kit172Quiet(runWorkflowOmicsReferenceDaAnalysis(fx$binding, k194Owner(fx), k194Ref(fx),
        list(formula_string = "~ 0 + group", da_q_val_thresh = 0.05, treat_lfc_cutoff = 0), kit172Contrasts()))
    invisible(fx)
}

#' A lipid or metab fixture after "Run Normalization Pipeline" (`input`), Apply (unless `apply` is FALSE) and the DA
#' (when `da`); closed when the calling test ends.
k207Normalised <- function(root, class, input = k194Input(), apply = TRUE, da = FALSE, env = parent.frame()) {
    fx <- k194Fixture(root, class)
    k194Defer(fx, envir = env)
    nd <- omicNormReactiveState()
    k194Run(fx, nd, input)
    k194ExpectRun(fx, nd, paste("KIT-207", class, basename(root)))
    if (apply) expect_true(isTRUE(k194Apply(fx, nd)), info = paste("apply precondition", class, basename(root)))
    if (da) k207Da(fx)
    list(fx = fx, nd = nd)
}

#' A sequence's head, checked as the row's precondition (an instrument check, VAL-54).
k207ExpectHead <- function(fx, state, info) expect_identical(k194Head(fx), state, info = paste("head precondition:", info))

# ------------------------------------------------------------------------------------------------
# KIT-212's harness.
# ------------------------------------------------------------------------------------------------

# The texts of section 7, as the ticket drafts them (VAL-41: the owner's wording pass may change them).
k212Texts <- list(
    stale = paste("The differential analysis was run on an earlier table, not the current one: the table changed",
        "after the analysis ran (a new normalisation run, a Revert, or the correlation filter applied again). Run the",
        "differential analysis again, then press 'Save Workflow Arguments' before copying or generating the report."),
    unmade_lead = "holds differential analysis files that no differential analysis of this session made",
    replay = paste("This record's differential analysis was run on a table other than the record's final table, so a replay",
        "would test a different table. Open the session, run the differential analysis on the final table, save, and replay",
        "the new record."))

k212Settings <- list(formula_string = "~ 0 + group", da_q_val_thresh = 0.05, treat_lfc_cutoff = 0)

#' The DA status and the unmade files, where the new functions exist (head), else "<absent>".
k212Status <- function(manager) if (exists("workflowDaResultStatus", mode = "function")) workflowDaResultStatus(manager) else "<absent>"
k212Unmade <- function(manager, paths, omic) {
    if (!exists("workflowDaFilesUnmade", mode = "function")) return("<absent>")
    basename(workflowDaFilesUnmade(manager, paths, omic))
}

#' A fresh project's directories for the omic, with the DA folders (the Summary's `projectDirs`).
k212Dirs <- function(root, omic) {
    paths <- list(source_dir = file.path(root, "scripts", omic), integration_dir = file.path(root, "integration"), base_dir = root,
        da_output_dir = file.path(root, "da_output"), publication_graphs_dir = file.path(root, "publication_graphs"),
        results_summary_dir = file.path(root, "results_summary"))
    for (path in paths) dir.create(path, recursive = TRUE, showWarnings = FALSE)
    stats::setNames(list(paths), omic)
}

#' The `project_dirs` global as it was before the calling test (Copy assigns it when absent).
k212KeepGlobal <- function(env = parent.frame()) {
    had <- exists("project_dirs", envir = .GlobalEnv, inherits = FALSE)
    old <- if (had) get("project_dirs", envir = .GlobalEnv, inherits = FALSE)
    withr::defer(if (had) assign("project_dirs", old, envir = .GlobalEnv) else
        if (exists("project_dirs", envir = .GlobalEnv, inherits = FALSE)) rm("project_dirs", envir = .GlobalEnv), envir = env)
}

#' "Copy to Publication Directory": the core's `copyWorkflowSummaryResults` under the omic's group with a recording
#' `copyFn`; the notifications ("[type] text"), the copy calls and the value.
k212Copy <- function(workflow, omic, dirs, argsSaved = FALSE) {
    rec <- new.env(parent = emptyenv())
    rec$notes <- character()
    rec$copies <- 0L
    k212KeepGlobal()
    testthat::local_mocked_bindings(workflowOmicsManagedActions = function(omic_type) NULL, .package = "MultiScholaR")
    kit172Quiet(utils::capture.output(rec$value <- copyWorkflowSummaryResults(omicGuiConfig(omic)$summary, workflow, dirs, omic,
        "k212", "d", argsSaved = argsSaved, notify = function(message, type, duration) rec$notes <- c(rec$notes,
            paste0("[", type, "] ", message)), copyFn = function(...) {
            rec$copies <- rec$copies + 1L
            list()
        }, catFn = function(...) NULL)))
    rec
}

#' "Generate Report": the core's `renderWorkflowSummaryReport` with the template in the project, a recording render
#' and record bundle; the notifications and the calls.
k212Report <- function(workflow, omic, dirs) {
    config <- omicGuiConfig(omic)$summary
    template <- if (is.function(config$report$template)) "k212.rmd" else config$report$template
    if (is.function(config$report$template)) config$report$template <- function(...) list(templateFilename = template)
    scripts <- file.path(dirs[[omic]]$base_dir, "scripts", omic)
    dir.create(scripts, recursive = TRUE, showWarnings = FALSE)
    writeLines("KIT-212 template", file.path(scripts, template))
    rec <- new.env(parent = emptyenv())
    rec$notes <- character()
    rec$render <- 0L
    rec$bundle <- 0L
    report <- file.path(dirs[[omic]]$base_dir, "k212-report.html")
    testthat::local_mocked_bindings(workflowOmicsManagedActions = function(omic_type) NULL, .package = "MultiScholaR")
    kit172Quiet(utils::capture.output(renderWorkflowSummaryReport(config, workflow, dirs, omic, "k212",
        notify = function(message, type, duration) rec$notes <- c(rec$notes, paste0("[", type, "] ", message)),
        progress = function(...) NULL, renderFn = function(omic_type, experiment_label, rmd_filename) {
            rec$render <- rec$render + 1L
            writeLines("KIT-212 report", report)
            report
        }, exportRecordFn = function(...) {
            rec$bundle <- rec$bundle + 1L
            invisible(NULL)
        }, catFn = function(...) NULL)))
    rec
}

k212Noted <- function(rec, text) any(grepl(text, rec$notes, fixed = TRUE))
k212Line <- function(row, label, rec, extra = "") cat(sprintf("KIT-212 %s %s\tnotes=%s%s\n", row, label,
    paste(rec$notes, collapse = " | "), extra))

#' A refused Save (E1's head expectation): the notice, `failed` once with the class, nothing written, `done` not called.
k212ExpectSaveRefused <- function(rec, dirs, omic, info) {
    source_dir <- dirs[[omic]]$source_dir
    info <- paste(info, "| notes:", paste(rec$notes, collapse = " | "))
    expect_true(paste("[error] The session could not be saved:", k212Texts$stale) %in% rec$notes, info = info)
    expect_identical(rec$done, 0L, info = info)
    expect_true(length(rec$failed) == 1L && inherits(rec$failed[[1L]], "workflow_summary_da_stale"), info = info)
    expect_false(file.exists(file.path(source_dir, "study_parameters.txt")), info = info)
    expect_false(file.exists(file.path(source_dir, "analysis_record.json")), info = info)
}

#' A refused Copy: the notice with `text`, nothing copied, no fallback study parameters.
k212ExpectCopyRefused <- function(rec, dirs, omic, text, info) {
    info <- paste(info, "| notes:", paste(rec$notes, collapse = " | "))
    expect_true(any(startsWith(rec$notes, "[error] Copy error:") & grepl(text, rec$notes, fixed = TRUE)), info = info)
    expect_identical(rec$copies, 0L, info = info)
    expect_false(file.exists(file.path(dirs[[omic]]$source_dir, "study_parameters.txt")), info = info)
}

#' A refused Report: the notice with `text`, no render, no bundle.
k212ExpectReportRefused <- function(rec, text, info) {
    info <- paste(info, "| notes:", paste(rec$notes, collapse = " | "))
    expect_true(any(startsWith(rec$notes, "[error] Report generation failed:") & grepl(text, rec$notes, fixed = TRUE)), info = info)
    expect_identical(rec$render, 0L, info = info)
    expect_identical(rec$bundle, 0L, info = info)
}

#' Q1 (lipid, metab): run, Apply, DA, Reset, run 2, Apply, no DA; head `<prefix>_correlation_filtered` of run 2.
k212Q1 <- function(root, class, env = parent.frame()) {
    q <- k207Normalised(root, class, da = TRUE, env = env)
    k194Reset(q$fx, q$nd)
    k194Run(q$fx, q$nd)
    k194ExpectRun(q$fx, q$nd, paste("KIT-212 Q1 run 2", class))
    expect_true(isTRUE(k194Apply(q$fx, q$nd)), info = paste("Q1 Apply 2 precondition", class))
    k207ExpectHead(q$fx, paste0(q$fx$prefix, "_correlation_filtered"), paste("Q1", class))
    q
}

#' The publication of the DA in process (the binding's `<prefix>_da_publication`) into `paths`.
k212Publish <- function(binding, manager, paths) kit172Quiet(runWorkflowOmicsPublicationInProcess(binding, manager,
    paste0(binding$prefix, "_da_publication"), list(), paths))

#' A record of the owner's head, as Save writes it (`writeWorkflowAnalysisRecord`, unchanged by the ticket), replayed on a
#' copy of the store cut to the design state (KIT-194's `k194Replay`, test-omics-norm-state-guards.R:1162-1183, with the
#' replay called directly): the record and the replay's error, or NULL.
k212Replay <- function(fx, root) {
    owner <- k194Owner(fx)
    dir.create(root, recursive = TRUE, showWarnings = FALSE)
    json <- file.path(root, "analysis_record.json")
    writeWorkflowAnalysisRecord(owner, json)
    record <- jsonlite::read_json(json, simplifyVector = FALSE)
    copy <- file.path(root, "copy")
    dir.create(copy)
    stopifnot(system2("cp", c("-a", file.path(fx$source_dir, "."), copy)) == 0L)
    harness <- new.env(parent = globalenv())
    sys.source(file.path(testthat::test_path("..", ".."), "tools", "profiling", "workflow_kit.R"), harness)
    api <- function(name) get(name, asNamespace("MultiScholaR"))
    store_root <- dirname(dirname(list.files(fx$source_dir, "^store[.]json$", recursive = TRUE, all.files = TRUE, full.names = TRUE)[[1L]]))
    copy_root <- file.path(copy, substring(store_root, nchar(fx$source_dir) + 2L))
    prepared <- harness$kitPrepareReplayProject(copy_root, c("initial", paste0(fx$prefix, "_raw_data_s4")), api)
    inputs <- stats::setNames(c(file.path(fx$source_dir, "import_sources", "LCMS_Pos_export.txt"), file.path(fx$source_dir, "design_matrix.tab")),
        c(paste0(fx$binding$design_input_prefix, 1L), "design_file"))
    scratch <- file.path(root, "scratch")
    dir.create(scratch)
    error <- tryCatch({
        kit172Quiet(replayWorkflowAnalysisRecord(record, prepared, inputs, fx$paths, scratch))
        NULL
    }, error = function(e) e)
    da <- record$analyses$differential_abundance
    cat(sprintf("KIT-212 replay %s\tsource_state=%s\tsteps=%s\terror=%s\n", basename(root), format(da$source_state %||% "<null>"),
        paste(vapply(record$steps, function(step) as.character(step$state), ""), collapse = ","),
        if (is.null(error)) "none" else paste(class(error)[[1L]], conditionMessage(error))))
    list(record = record, error = error)
}

#' The declarations journey of T10 (test-omics-norm-state-guards.R:1190-1202) up to the correlation filter: log2 (offset 0),
#' quantile normalisation, the correlation filter (group, 0.5).
k212Journey <- function(fx) {
    owner <- k194Owner(fx)
    config <- fx$workflow$config_list
    kit172Quiet({
        runWorkflowOmicsReferenceLog2(fx$binding, owner, 0, config)
        runWorkflowOmicsReferenceBetweenSampleNormalisation(fx$binding, owner, "quantile", config)
        runWorkflowOmicsReferenceCorrelation(fx$binding, owner, "group", 0.5, config)
    })
    invisible(fx)
}

#' The proteomics native owner of KIT-206 (`k206Owner`, test-kit206-prot-norm-export-guard.R:35-39) normalised (RUV
#' skipped) and correlation-filtered at 0.5 (`k206Click`, `k206Correlation`, :42-55), with the DA folders added to its paths.
k212ProtContrasts <- data.frame(contrasts = "grouptreated-groupcontrol", friendly_names = "treated_vs_control",
    full_format = "treated_vs_control=grouptreated-groupcontrol", stringsAsFactors = FALSE)
k212Prot <- function(envir = parent.frame()) {
    fixture <- protNormCoreOwner(envir)
    fixture$workflow$contrasts_tbl <- k212ProtContrasts
    fixture$paths$da_output_dir <- file.path(fixture$root, "results", "proteomics", "da_output")
    fixture$paths$pathway_dir <- file.path(fixture$root, "results", "proteomics", "pathway_enrichment")
    for (path in fixture$paths) dir.create(path, recursive = TRUE, showWarnings = FALSE)
    norm <- new.env(parent = emptyenv())
    settings <- list(norm_method = "quantile", ruv_mode = "skip", auto_percentage_min = NULL, auto_percentage_max = NULL,
        separation_metric = NULL, k_penalty_weight = NULL, adaptive_k_penalty = NULL, ruv_percentage = 50, ruv_k = 1,
        color_variable = "group", shape_variable = NULL)
    args <- list(settings = settings, workflow = fixture$workflow, norm_data = norm, paths = fixture$paths, omic_type = "proteomics",
        experiment_label = "core", post_qc = function(reference) invisible(NULL), ruv_qc = function(reference) invisible(NULL),
        grouping = function() "group", pre_qc = NULL)
    protNormCoreQuiet(do.call(runDIAReferenceNormalization, args[intersect(names(args), names(formals(runDIAReferenceNormalization)))]))
    protNormCoreQuiet(runDIAReferenceCorrelation(fixture$workflow, norm, fixture$paths, 0.5, FALSE, "proteomics", "core",
        function() "group", function() list(color_var = "group", shape_var = NULL)))
    norm$correlation_filtering_complete <- TRUE
    fixture$norm <- norm
    fixture
}

#' The proteomics DA in process (`runDIAReferenceDaAnalysis`, as helper-go-propagation.R:149-153 calls it), its
#' publication when `publish` (`produceDIADaPublication` into the fixture's DA folders, the pages not saved as assets).
k212ProtDa <- function(fixture, publish = FALSE) {
    manager <- fixture$manager
    protNormCoreQuiet(withDIAStateScope(manager, function() {
        da_data <- new.env(parent = emptyenv())
        da_data$current_s4_object <- getDIAStateReference(manager)
        da_data$contrasts_available <- k212ProtContrasts$contrasts
        manifest <- runDIAReferenceDaAnalysis(k212Settings, da_data, fixture$workflow, fixture$paths, contrasts = k212ProtContrasts,
            publish = FALSE)
        if (publish) produceDIADaPublication(manager, manifest, data.frame(), 0.05,
            fixture$paths[c("da_output_dir", "publication_graphs_dir")], save_assets = FALSE)
        manifest
    }))
}

#' The Summary's directories for a proteomics fixture (its paths, a fresh source folder, the integration folder).
k212ProtDirs <- function(fixture, name) {
    paths <- fixture$paths
    paths$source_dir <- file.path(fixture$root, name, "scripts", "proteomics")
    paths$integration_dir <- file.path(fixture$root, name, "integration")
    paths$base_dir <- file.path(fixture$root, name)
    dir.create(paths$source_dir, recursive = TRUE, showWarnings = FALSE)
    list(proteomics = paths)
}

#' A completed enrichment on the head (`k206Enrichment`, test-kit206-prot-norm-export-guard.R:90-96), raised on error.
k212Enrichment <- function(fixture) {
    results <- methods::new("EnrichmentResults", contrasts = tibble::tibble(contrasts = "grouptreated-groupcontrol"),
        enrichment_data = list(`grouptreated-groupcontrol` = list(up = data.frame(term_id = "GO:0008150", p_value = 0.01))),
        enrichment_plots = list(), enrichment_plotly = list(), enrichment_summaries = list())
    protNormCoreQuiet(fixture$manager$saveState("enrichment_completed", results, list(q_cutoff = 0.05), "fixture enrichment"))
}

#' "Export Filtered Data Session" of the proteomics tab: `k206Export` (test-kit206-prot-norm-export-guard.R:125-145) with the
#' real export workflow; its notices and errors.
k212Export <- function(fixture, norm, mode = "skip") {
    log <- new.env(parent = emptyenv())
    log$notes <- character()
    log$errors <- list()
    note <- function(ui, type = "default", duration = 5, ...) {
        log$notes <- c(log$notes, sprintf("[%s %s] %s", type, format(duration), paste(as.character(ui), collapse = " ")))
        invisible(NULL)
    }
    quiet <- function(...) invisible(NULL)
    run <- function(workflowData, normData, input, sourceDir) runProtNormExportSessionWorkflow(
        workflowData, normData, input, sourceDir, withProgressFn = function(message, value, expr) force(expr),
        incProgressFn = function(...) NULL)
    log$value <- protNormCoreQuiet(withDIAStateScope(fixture$manager, function() runProtNormExportObserver(
        list(norm_method = "quantile", ruv_mode = mode), fixture$workflow, norm, fixture$paths,
        runExportSessionWorkflowFn = run, showNotificationFn = note, messageFn = quiet,
        handleExportErrorFn = function(e) {
            log$errors <- c(log$errors, list(e))
            handleProtNormExportError(e, showNotificationFn = note, messageFn = quiet)
        })))
    expect_length(log$errors, 0L)
    as.list(log)
}

#' The archive reopened into a new owner and project (`k206Reopen`, test-kit206-prot-norm-export-guard.R:178-199).
k212Reopen <- function(archive, envir = parent.frame()) {
    root <- withr::local_tempdir(.local_envir = envir)
    paths <- list(source_dir = file.path(root, "scripts", "proteomics"),
        protein_qc_dir = file.path(root, "results", "proteomics", "protein_qc"),
        peptide_qc_dir = file.path(root, "results", "proteomics", "peptide_qc"))
    paths$publication_graphs_dir <- file.path(root, "results", "proteomics", "publication_graphs")
    paths$time_dir <- file.path(paths$publication_graphs_dir, "filtering_qc")
    for (path in unlist(paths)) dir.create(path, recursive = TRUE, showWarnings = FALSE)
    session <- new.env(parent = emptyenv())
    session$userData <- new.env(parent = emptyenv())
    session$onFlush <- session$onFlushed <- session$onSessionEnded <- function(...) function() invisible(NULL)
    workflow <- new.env(parent = emptyenv())
    workflow$state_manager <- WorkflowState$new(audit_enabled = FALSE)
    workflow$qc_params <- list()
    imported <- protNormCoreQuiet(openDIANativeSession(archive, paths$source_dir))
    protNormCoreQuiet(adoptDIANativeSession(imported, workflow, session, paths$source_dir, paths))
    withr::defer(closeDIAStateAdapter(workflow$state_manager), envir = envir)
    norm <- new.env(parent = emptyenv())
    testthat::local_mocked_bindings(updateSelectInput = function(...) invisible(NULL), .package = "shiny")
    protNormCoreQuiet(restoreProtNormTabState(NULL, workflow, norm))
    list(manager = workflow$state_manager, workflow = workflow, paths = paths, root = root, norm = norm)
}

#' The Save of KIT-206 E11b (test-kit206-prot-norm-export-guard.R:641-653): the writers of the parameters and the record
#' replaced, the notices, `done` and `failed` recorded.
k212ProtSave <- function(fixture, dirs) {
    seen <- new.env(parent = emptyenv())
    seen$notes <- character()
    seen$failed <- list()
    seen$done <- list()
    protNormCoreQuiet(saveWorkflowSummaryParameters(omicGuiConfig("proteomics")$summary, fixture$workflow, dirs,
        "proteomics", "k212", "d", notify = function(message, type = NULL, duration = NULL) {
            seen$notes <- c(seen$notes, as.character(message))
            invisible(NULL)
        }, done = function(saved) seen$done <- c(seen$done, list(saved)),
        failed = function(error) seen$failed <- c(seen$failed, list(error)),
        createWorkflowArgsFn = function(...) "study_parameters.txt", writeAnalysisRecordFn = function(...) NULL,
        catFn = function(...) NULL))
    seen
}

#' The module-CI proteomics owner of W1's `prot_owner` (test-omic-summary-module.R:824-843): a reference-state session
#' owner holding MCI-009's object as `correlation_filtered`, no DA; with the module-CI project and its publication inputs.
k212ProtW1Owner <- function(env = parent.frame()) {
    config <- module_ci_prot_summary_branch_config("DIA")
    bus <- module_ci_prot_summary_workflow(config)
    project_dirs <- module_ci_prot_summary_paths()
    paths <- project_dirs$proteomics
    withr::defer(unlink(paths$base_dir, recursive = TRUE), envir = env)
    module_ci_prot_summary_write_publication_inputs(project_dirs)
    design <- module_ci_prot_da_design()
    rownames(design) <- NULL
    quant <- module_ci_prot_da_table()
    rownames(quant) <- NULL
    object <- module_ci_prot_da_object(workflow_type = config$workflow_type, use_limpa = config$use_limpa,
        report_template = config$report_template, design = design, protein_quant_table = quant,
        args = module_ci_prot_summary_args(config))
    store_root <- file.path(paths$base_dir, "store")
    dir.create(store_root, showWarnings = FALSE)
    manager <- newDIAStateAdapter(newWorkflowEntityStore(store_root, "kit212-prot", config$workflow_type),
        WorkflowState$new(audit_enabled = FALSE))
    withr::defer(closeDIAStateAdapter(manager), envir = env)
    manager$.context$reference_states <- TRUE
    manager$.context$experiment_paths <- paths
    manager$saveState("correlation_filtered", object, list(), "MCI-009 fixture state")
    bus$state_manager <- manager
    list(bus = bus, project_dirs = project_dirs, manager = manager)
}

#' The DA tab of lipid mounted on a `shiny::reactiveValues` bus over the fixture's owner (a plain-environment bus cannot
#' fire the tab's observer). `code(tab, run_da, fire, bus, session)` runs inside `testServer`: `tab` the tab's values,
#' `run_da()` presses Run DA (in process), `fire()` sets `state_update_trigger` and flushes.
k212Tab <- function(fx, code) {
    bus <- shiny::reactiveValues(state_manager = k194Owner(fx), contrasts_tbl = kit172Contrasts(), design_matrix = kit172Design(),
        config_list = fx$workflow$config_list, tab_status = list(normalization = "complete", differential_analysis = "pending"),
        state_update_trigger = NULL)
    testthat::local_mocked_bindings(workflowOmicsManagedActions = function(omic_type) NULL, .package = "MultiScholaR")
    kit172Quiet(utils::capture.output(shiny::testServer(mod_lipid_da_server, args = list(workflow_data = bus, experiment_paths = list(),
        omic_type = "lipidomics", experiment_label = "k212"), {
        run_da <- function(n = 1L) {
            session$setInputs(formula_string = "~ 0 + group", da_q_val_thresh = 0.05, treat_lfc_cutoff = 0)
            session$setInputs(run_da_analysis = n)
            session$flushReact()
        }
        fire <- function() {
            bus$state_update_trigger <- Sys.time()
            session$flushReact()
        }
        code(session$returned, run_da, fire, bus, session)
    })))
    invisible(bus)
}


#' A DIA (`iq` route) owner on the "small" recipe's protein state, accepted at `correlation_filtered`, with the paths
#' the DA publication writes to, the tab's `da_data` and a plain-environment workflow. The owner is closed at the end.
savedDaW1Owner <- function(project = "kit112-w1", envir = parent.frame(), before_save = NULL) {
    root <- normalizePath(withr::local_tempdir(.local_envir = envir))
    paths <- list(source_dir = file.path(root, "source"), da_output_dir = file.path(root, "da"),
        publication_graphs_dir = file.path(root, "graphs"))
    for (path in paths) dir.create(path)
    manager <- newDIAStateAdapter(newWorkflowEntityStore(root, project, "DIA"), WorkflowState$new(audit_enabled = FALSE))
    manager$.context$reference_states <- TRUE
    manager$setWorkflowType("DIA")
    manager$.context$experiment_paths <- paths
    # GUI audit A112-01: a test names the provenance the state is committed under (none pending is `legacy`).
    if (is.function(before_save)) before_save(manager)
    manager$saveState("correlation_filtered", makeWorkflowCheckpointProtein(readWorkflowCheckpointRecipe("small")),
        list(), "fixture")
    workflow <- list2env(list(state_manager = manager, tab_status = list()))
    withr::defer(closeDIAStateAdapter(workflow$state_manager), envir = envir)
    da <- new.env(parent = emptyenv())
    da$current_s4_object <- getDIAStateReference(manager)
    list(root = root, project = project, paths = paths, workflow = workflow, da = da)
}

savedDaW1Contrasts <- function(n = 1L) {
    data.frame(contrasts = c("treated.vs.control", "control.vs.treated")[seq_len(n)],
        full_format = c("treated.vs.control=grouptreated-groupcontrol", "control.vs.treated=groupcontrol-grouptreated")[seq_len(n)],
        friendly_names = c("treated.vs.control", "control.vs.treated")[seq_len(n)], stringsAsFactors = FALSE)
}

#' One managed action in this process: the request (with its named input) through the worker's entry point, the
#' acceptance and what the tab reads. A worker failure is raised with the worker's message. Each delivered value's size
#' and names are kept in `log`.
savedDaW2Dispatch <- function(fixture, controller, operation, parameters, input, log = NULL) {
    store <- fixture$workflow$state_manager$.context$store
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
    path <- file.path(paths$request_dir, "request.rds")
    writeWorkflowExecutionMessage(path, request)
    withr::with_envvar(c(MULTISCHOLAR_WORKFLOW_WORKER = "1"), suppressWarnings(suppressMessages(utils::capture.output(try(
        workflowWorkerMain(c(path, digest::digest(file = path, algo = "sha256"))), silent = TRUE)))))
    message <- readRDS(file.path(paths$request_dir, "worker-result.rds"))
    if (!identical(message$status, "complete")) stop(message$error, call. = FALSE)
    if (is.environment(log)) log$values[[length(log$values) + 1L]] <- list(operation = operation,
        bytes = as.numeric(object.size(message$result)), names = names(message$result), status = message$result$status,
        base_fit = !is.null(message$result$display$base_fit), mutable = hasMutableWorkflowEntityOwner(message$result))
    controller$job <- list(request = list(token = token, operation = operation, config = list(mode = "adapters"),
        parameters = parameters), paths = paths, completion_context = validateWorkflowOmicsCompletionContext(operation, NULL))
    receipt <- acceptWorkflowActionResult(controller, message[c("schema", "status", "token", "result")])
    accepted <- resumeWorkflowAction(controller, id)
    c(readManagedOmicsResult(fixture$workflow$state_manager, accepted, NULL), list(acceptance = receipt$status, accepted = accepted))
}

#' The managed DA click's front: the request half in the controller, each contrast on the worker (the route in the
#' request), the base fit adopted only when delivered (R1), the accept half in the controller. `today = TRUE` is
#' today's managed front instead: the whole runner in the controller without its publication.
savedDaW2Front <- function(fixture, dispatch, settings, contrasts, today = FALSE) {
    manager <- fixture$workflow$state_manager
    annotations <- fixture$workflow$uniprot_dat_cln
    if (today) return(runDIAReferenceDaAnalysis(settings, fixture$da, fixture$workflow, fixture$paths, contrasts, annotations,
        publish = FALSE))
    request <- withDIAStateScope(manager, function() diaDaRequest(manager, settings, fixture$da, fixture$workflow, contrasts,
        annotations))
    route <- withDIAStateScope(manager, function() diaActiveRoute(manager))
    results <- list()
    for (step in request$steps) {
        delivered <- dispatch("DA-contrast-output", c(step, list(route = route)), request$source)
        results[[step$name]] <- delivered$value
        if (!is.null(delivered$display$base_fit)) manager$.context$da_base_fit_reference <- delivered$display$base_fit
    }
    withDIAStateScope(manager, function() acceptDIADaManifest(manager, request$manifest, results))
}

#' One managed DA click: the front, then the publication step on the worker and its delivery to the tab.
savedDaW2Managed <- function(fixture, dispatch, settings, contrasts, today = FALSE) {
    manager <- fixture$workflow$state_manager
    manifest <- savedDaW2Front(fixture, dispatch, settings, contrasts, today)
    # GUI audit A112-06: the annotation table travels as the run plugin sends it, a saved display table (none: NULL).
    published <- dispatch("da_publication", list(annotations = saveDIADisplayTable(manager, fixture$workflow$uniprot_dat_cln)),
        getDIAStateReference(manager))
    delivered <- withDIAStateScope(manager, function() readManagedDaResult(manager, published$accepted, NULL))
    completeDIADaDisplay(fixture$da, fixture$workflow, delivered$manifest)
    manifest
}

#' A dispatcher on the fixture's current owner, with a controller per owner (a reopen gets its own).
savedDaW2Dispatcher <- function(fixture, log) {
    owner <- NULL
    controller <- NULL
    function(operation, parameters, input) {
        manager <- fixture$workflow$state_manager
        if (!identical(manager, owner)) {
            owner <<- manager
            controller <<- omicsWorkerController(list(owner = manager))
        }
        savedDaW2Dispatch(fixture, controller, operation, parameters, input, log)
    }
}

savedDaW2Mocks <- function(env = parent.frame()) {
    testthat::local_mocked_bindings(workflowExecutionBuildIdentity = function() strrep("f", 64),
        finishWorkflowExecutionJob = function(owner, status) { owner$job <- NULL; list(status = status) },
        .package = "MultiScholaR", .env = env)
}

# ------------------------------------------------------------------------------------------------
# KIT-243's harness.
# ------------------------------------------------------------------------------------------------

# The texts of section 7, as the ticket drafts them (VAL-41: the owner's wording pass may change them).
k243Texts <- list(
    stale = paste("The saved workflow parameters do not describe the current analysis: a step or the differential analysis was",
        "run, reverted or changed after the analysis record (analysis_record.json) was last saved, or it was saved by another",
        "version of MultiScholaR. Press 'Save Workflow Arguments' again, then copy the results or generate the report. If this",
        "message returns after saving, the analysis record could not be written: see the warning shown when saving."),
    pending = paste("The files of the last differential analysis were not all written to the project folder, so it may still",
        "hold files of an earlier analysis. Run the differential analysis again, then copy the results or generate the report."),
    removal = "An earlier file of the differential analysis could not be removed before this analysis's files were written:")

k243MarkerName <- ".multischolar_da_publication_pending"
k243Marker <- function(paths) file.path(paths$publication_graphs_dir, k243MarkerName)

# "Earlier figures" are planted as these bytes under the publication's own names (KIT-212's planted tables, section 4).
k243Bytes <- charToRaw("KIT-243 planted: a file of an earlier differential analysis")
k243Plant <- function(files) {
    for (file in files) {
        dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
        writeBin(k243Bytes, file)
    }
    invisible(files)
}
#' Whether a file still holds the planted bytes.
k243Planted <- function(file) isTRUE(file.exists(file) && !dir.exists(file) &&
    identical(readBin(file, "raw", n = length(k243Bytes) + 1L), k243Bytes))
k243PlantedOf <- function(files) basename(Filter(k243Planted, files))

#' The last Save's status, where the new function exists (head), else "<absent>".
k243Status <- function(manager, sourceDir) {
    if (!exists("workflowSummarySaveStatus", mode = "function")) return("<absent>")
    workflowSummarySaveStatus(manager, sourceDir)
}

#' The names the store's last DA's publication writes (`workflowDaMadeFiles`, KIT-212), in the owner's scope.
k243Made <- function(manager, paths, omic) withWorkflowStateScope(manager, function() workflowDaMadeFiles(manager, paths, omic))

#' The DA tab's run (`executeOmicDaRun`, what Run DA calls) on the head: in process (`runner` NULL) or managed; the
#' notices are recorded by the calling test's `k194Notes()`.
k243Da <- function(workflow, paths, da, binding, contrasts, settings = k212Settings, runner = NULL) {
    testthat::local_mocked_bindings(workflowOmicsManagedActions = function(omic_type) runner, .package = "MultiScholaR")
    tab <- list2env(list(current_state_name = NULL, analysis_complete = FALSE))
    outcome <- kit172Quiet(executeOmicDaRun(tab, workflow, NULL, paths, settings, contrasts, da, binding))
    c(outcome, list(tab = tab))
}

#' A managed runner for lipidomics, metabolomics and transcriptomics that runs each requested action here as the worker
#' would and delivers it (KIT-194's `k194Runner` delivers an empty display for a product, so a managed DA through it would
#' accept no contrast): a contrast through its declaration on the head (as `runWorkflowOmicsReferenceDaAnalysis`), the tab
#' results through the worker's handler, the publication through the in-process publication into `paths` (what the
#' worker stages and `publishWorkflowOmicsResult` copies). `fail` names an operation that fails as the worker's would.
k243Runner <- function(manager, binding, paths, fail = NULL) {
    rec <- new.env(parent = emptyenv())
    rec$operations <- character()
    worker <- function(action) {
        operation <- action$operation
        head <- getWorkflowStateReference(manager)
        if (identical(operation, paste0(binding$prefix, "_da_contrast"))) return(withWorkflowStateScope(manager, function() {
            chosen <- workflowRequestOperation(manager, workflowOperationEntry(operation), head, action$parameters)
            list(value = runSavedOperation(manager, chosen$declaration, head, chosen$parameters)$value)
        }))
        if (identical(operation, paste0(binding$prefix, "_da_tab_results"))) return(withWorkflowStateScope(manager, function()
            computeWorkflowOmicsDaTabResultsAction(manager, list(operation = operation, input_references = list(state = head)), NULL)))
        if (identical(operation, paste0(binding$prefix, "_da_publication")))
            return(list(display = runWorkflowOmicsPublicationInProcess(binding, manager, operation, list(), paths)))
        stop("KIT-243 runner: unexpected operation ", operation, call. = FALSE)
    }
    runner <- new.env(parent = emptyenv())
    runner$step <- list(active = function() NULL, pending = function() NULL)
    runner$progress <- function(...) list(set = function(...) invisible(NULL), inc = function(...) invisible(NULL),
        close = function(...) invisible(NULL))
    runner$run <- function(actions, done = NULL, failed = NULL, dropped = NULL) {
        outcome <- tryCatch({
            if (is.function(actions)) actions <- actions()
            for (action in actions) {
                rec$operations <- c(rec$operations, action$operation)
                if (identical(action$operation, fail)) stop("KIT-243: the worker's ", fail, " failed", call. = FALSE)
                action$deliver(worker(action))
            }
            NULL
        }, error = function(error) error)
        if (inherits(outcome, "error")) {
            if (is.function(failed)) failed(outcome) else stop(outcome)
        } else if (is.function(done)) done()
        invisible(NULL)
    }
    rec$runner <- runner
    rec
}

#' A refused Copy (no Save check): the notice with `text`, nothing copied.
k243ExpectCopyRefused <- function(rec, text, info) {
    info <- paste(info, "| notes:", paste(rec$notes, collapse = " | "))
    expect_true(any(startsWith(rec$notes, "[error] Copy error:") & grepl(text, rec$notes, fixed = TRUE)), info = info)
    expect_identical(rec$copies, 0L, info = info)
}

#' The proteomics Save with the real writers (P03): `k212ProtSave` without its two writer overrides.
k243ProtSave <- function(fixture, dirs) {
    seen <- new.env(parent = emptyenv())
    seen$notes <- character()
    seen$failed <- list()
    seen$done <- list()
    protNormCoreQuiet(saveWorkflowSummaryParameters(omicGuiConfig("proteomics")$summary, fixture$workflow, dirs,
        "proteomics", "k243", "d", notify = function(message, type = NULL, duration = NULL) {
            seen$notes <- c(seen$notes, as.character(message))
            invisible(NULL)
        }, done = function(saved) seen$done <- c(seen$done, list(saved)),
        failed = function(error) seen$failed <- c(seen$failed, list(error)), catFn = function(...) NULL))
    seen
}

#' `k212ProtDa` with the settings named (KIT-206's fixture; the publication, when asked, straight into its folders).
k243ProtDa <- function(fixture, settings = k212Settings, publish = FALSE) {
    manager <- fixture$manager
    protNormCoreQuiet(withDIAStateScope(manager, function() {
        da_data <- new.env(parent = emptyenv())
        da_data$current_s4_object <- getDIAStateReference(manager)
        da_data$contrasts_available <- k212ProtContrasts$contrasts
        manifest <- runDIAReferenceDaAnalysis(settings, da_data, fixture$workflow, fixture$paths, contrasts = k212ProtContrasts,
            publish = FALSE)
        if (publish) produceDIADaPublication(manager, manifest, data.frame(), settings$da_q_val_thresh,
            fixture$paths[c("da_output_dir", "publication_graphs_dir")], save_assets = FALSE)
        manifest
    }))
}

#' The proteomics DA click in process with its publication (`runDIAReferenceDaAnalysis`, publish = TRUE: S1), on
#' KIT-112's `savedDaW1Owner`; the error, or NULL.
k243ProtClick <- function(fixture, qvalue = 0.05) {
    contrasts <- savedDaW1Contrasts(1L)
    fixture$workflow$contrasts_tbl <- contrasts
    fixture$da$contrasts_available <- contrasts$contrasts
    settings <- list(formula_string = "~ 0 + group", da_q_val_thresh = qvalue, treat_lfc_cutoff = 0)
    tryCatch({
        suppressWarnings(suppressMessages(utils::capture.output(runDIAReferenceDaAnalysis(settings, fixture$da, fixture$workflow,
            fixture$paths, contrasts, NULL))))
        NULL
    }, error = function(error) error)
}

#' The proteomics DA click on the managed path through the worker (`savedDaW2Managed`: each contrast and the publication
#' on the worker's entry point, the publication through `publishWorkflowDaResult`: S2 (a)); the error, or NULL.
k243ProtWorkerClick <- function(fixture, dispatch, qvalue = 0.05) {
    contrasts <- savedDaW1Contrasts(1L)
    fixture$workflow$contrasts_tbl <- contrasts
    fixture$da$contrasts_available <- contrasts$contrasts
    settings <- list(formula_string = "~ 0 + group", da_q_val_thresh = qvalue, treat_lfc_cutoff = 0)
    tryCatch({
        suppressWarnings(suppressMessages(utils::capture.output(savedDaW2Managed(fixture, dispatch, settings, contrasts))))
        NULL
    }, error = function(error) error)
}

#' The Summary's directories for a `savedDaW1Owner` fixture.
k243ProtFixtureDirs <- function(fixture) list(proteomics = c(fixture$paths, list(base_dir = fixture$root,
    integration_dir = file.path(fixture$root, "integration"))))

#' The proteomics DA tab's Run DA on a managed session (R1-04): the handler (`da_server_run_analysis_handler`) mounted with
#' the managed application present, the contrasts run here as the worker would on a runner set on the session, and the
#' publication step (`installWorkflowManagedStep`) replaced by one whose `apply` records the call and never completes
#' (the worker's publication failed; `publishWorkflowDaResult` is never called).
k243ProtManagedPress <- function(fixture, qvalue) {
    manager <- fixture$workflow$state_manager
    rec <- new.env(parent = emptyenv())
    rec$applied <- 0L
    rec$contrasts <- 0L
    step <- list(apply = function(...) {
        rec$applied <- rec$applied + 1L
        invisible(NULL)
    }, active = function() NULL, pending = function() NULL)
    runner <- new.env(parent = emptyenv())
    runner$step <- list(active = function() NULL, pending = function() NULL)
    runner$run <- function(actions, done = NULL, failed = NULL, dropped = NULL) {
        outcome <- tryCatch({
            for (action in actions) {
                delivered <- withDIAStateScope(manager, function() {
                    parameters <- action$parameters
                    route <- parameters$route
                    parameters$route <- NULL
                    computed <- computeDIAReferenceDaContrast(manager, action$input, parameters, route)
                    list(value = computed$run$value, display = list(base_fit = computed$base_fit))
                })
                rec$contrasts <- rec$contrasts + 1L
                action$deliver(delivered)
            }
            NULL
        }, error = function(error) error)
        if (inherits(outcome, "error")) failed(outcome) else done()
        invisible(NULL)
    }
    testthat::local_mocked_bindings(installWorkflowManagedStep = function(...) step,
        workflowSessionApplication = function(session) "kit243-managed-application", .package = "MultiScholaR")
    contrasts <- savedDaW1Contrasts(1L)
    fixture$workflow$contrasts_tbl <- contrasts
    fixture$da$contrasts_available <- contrasts$contrasts
    server <- function(input, output, session) {
        session$userData$workflow_omics_actions <- list(proteomics = runner)
        session$userData$dia_capacity_manager <- manager
        da_server_run_analysis_handler(input, output, session, session$ns, fixture$da, fixture$workflow, fixture$paths)
    }
    kit172Quiet(utils::capture.output(shiny::testServer(server, {
        session$setInputs(formula_string = "~ 0 + group", da_q_val_thresh = qvalue, treat_lfc_cutoff = 0)
        session$setInputs(run_da_analysis = 1L)
    })))
    rec
}

k243Line <- function(row, label, ...) cat(sprintf("KIT-243 %s %s\t%s\n", row, label, paste(..., sep = "\t")))

# ------------------------------------------------------------------------------------------------
# Part A, E1 to E3 (failing first) with G8 (lipid): lipid and metabolomics heatmaps, in process (S3) and managed (S4).
# ------------------------------------------------------------------------------------------------

for (case in list(list(row = "E1", class = "lipid", managed = FALSE), list(row = "E2", class = "lipid", managed = TRUE),
    list(row = "E3", class = "metab", managed = FALSE))) {
    test_that(paste("KIT-243", case$row, case$class, if (case$managed) "managed" else "in process",
        "a second DA removes the heatmaps an earlier DA left under its names", if (case$class == "lipid") "(G8: Save Heatmap kept)"), {
        skip_if_not(workflowOwnerCapable())
        root <- withr::local_tempdir()
        seen <- k194Notes()
        omic <- k194Omic(case$class)
        q <- k207Normalised(file.path(root, "fixture"), case$class)
        owner <- k194Owner(q$fx)
        paths <- k212Dirs(file.path(root, "project"), omic)[[omic]]
        run <- function() {
            runner <- if (case$managed) k243Runner(owner, q$fx$binding, paths)$runner
            k243Da(q$fx$workflow, paths, q$fx$config$da, q$fx$binding, kit172Contrasts(), runner = runner)
        }
        first <- run()
        expect_identical(first$status, if (case$managed) "dispatched" else "success", info = paste(case$row, "DA1 precondition",
            first$errorMessage %||% "", paste(seen$notes, collapse = " | ")))
        heatmaps <- grep("/Heatmaps/", k243Made(owner, paths, omic), value = TRUE, fixed = TRUE)
        drawn <- basename(heatmaps[file.exists(heatmaps)])
        expect_gt(length(heatmaps), 1L)
        k243Plant(heatmaps)
        # G8: heatmaps the user saved with Save Heatmap (`Heatmap/`, one letter less).
        saved <- file.path(paths$publication_graphs_dir, "Heatmap", paste0(omicDaHeatmapSavePrefix("lipid_", "groupB-groupA"),
            c("_heatmap.png", "_heatmap.pdf", "_heatmap_methods.md", "_clusters.csv")))
        if (identical(case$class, "lipid")) k243Plant(saved)
        second <- run()
        left <- k243PlantedOf(heatmaps)
        candidates <- k243PlantedOf(workflowDaCandidateFiles(paths, omic))
        k243Line(case$row, case$class, sprintf("drawn_by_DA1=%s", paste(drawn, collapse = ",")), sprintf("planted=%d", length(heatmaps)),
            sprintf("left=%s", paste(left, collapse = ",")), sprintf("marker=%s", file.exists(k243Marker(paths))),
            sprintf("notes=%s", paste(seen$notes, collapse = " | ")))
        expect_identical(second$status, first$status, info = paste(case$row, "DA2 precondition", second$errorMessage %||% ""))
        expect_identical(left, character(), info = paste(case$row, "planted heatmaps left"))
        expect_identical(candidates, character(), info = paste(case$row, "planted candidates"))
        expect_false(file.exists(k243Marker(paths)), info = paste(case$row, "marker after a completed publication"))
        if (identical(case$class, "lipid")) {
            expect_identical(k243PlantedOf(saved), basename(saved), info = paste(case$row, "G8 Save Heatmap files"))
        }
    })
}

# ------------------------------------------------------------------------------------------------
# Part A, E4, E4b, E4m (failing first) with G8 (proteomics): the count barplot and the interactive page-or-note pair.
# ------------------------------------------------------------------------------------------------

test_that("KIT-243 E4, G8: proteomics in process, a DA drawing no count barplot removes the earlier one; Save Heatmap kept", {
    skip_if_not_installed("Glimma")
    skip_if_not(workflowOwnerCapable())
    fixture <- savedDaW1Owner("kit243-e4")
    graphs <- fixture$paths$publication_graphs_dir
    expect_null(k243ProtClick(fixture), info = "E4 DA1 precondition")
    barplot <- file.path(graphs, "NumSigDaMolecules", paste0("da_proteins_num_sig_da_molecules.", c("png", "pdf")))
    drawn <- file.exists(barplot)
    k243Plant(barplot)
    saved <- file.path(graphs, "Heatmap", paste0(omicDaHeatmapSavePrefix("prot_", "treated.vs.control"),
        c("_heatmap.png", "_heatmap.pdf", "_heatmap_methods.md", "_clusters.csv")))
    k243Plant(saved)
    # DA2 at a threshold nothing passes: no protein significant, so no barplot is drawn.
    error <- k243ProtClick(fixture, qvalue = 1e-300)
    left <- k243PlantedOf(barplot)
    k243Line("E4", "prot", sprintf("barplot_drawn_by_DA1=%s", paste(drawn, collapse = ",")), sprintf("left=%s", paste(left, collapse = ",")),
        sprintf("error=%s", if (is.null(error)) "none" else conditionMessage(error)), sprintf("marker=%s", file.exists(k243Marker(fixture$paths))))
    expect_null(error, info = "E4 DA2 precondition")
    expect_identical(left, character(), info = "E4 planted barplot left")
    expect_false(file.exists(k243Marker(fixture$paths)), info = "E4 marker after a completed publication")
    expect_identical(k243PlantedOf(saved), basename(saved), info = "G8 Save Heatmap files")
})

#' E4b and E4m: the interactive volcano's page or its not-available note, whichever DA1 did not write, planted; DA2 is
#' the identical request (so the publication's `produce` is skipped and the staged files copied again; R1-01).
k243PageOrNote <- function(fixture, click, row) {
    interactive <- file.path(fixture$paths$publication_graphs_dir, "Interactive_Volcano_Plots")
    expect_null(click(), info = paste(row, "DA1 precondition"))
    page <- file.path(interactive, diaInteractiveVolcanoFileName("treated.vs.control"))
    note <- diaInteractiveVolcanoUnavailablePath(interactive, "treated.vs.control")
    written <- c(page = file.exists(page), note = file.exists(note))
    expect_identical(sum(written), 1L, info = paste(row, "DA1 wrote one of page and note"))
    other <- if (written[["page"]]) note else page
    k243Plant(other)
    error <- click()
    k243Line(row, "prot", sprintf("DA1_wrote=%s", paste(names(written)[written], collapse = ",")), sprintf("planted=%s", basename(other)),
        sprintf("left=%s", k243Planted(other)), sprintf("error=%s", if (is.null(error)) "none" else conditionMessage(error)),
        sprintf("marker=%s", file.exists(k243Marker(fixture$paths))))
    expect_null(error, info = paste(row, "DA2 precondition"))
    expect_false(k243Planted(other), info = paste(row, "planted", basename(other)))
    expect_false(file.exists(k243Marker(fixture$paths)), info = paste(row, "marker after a completed publication"))
}

test_that("KIT-243 E4b: proteomics in process, the interactive page or note DA1 did not write is removed by DA2 (S1)", {
    skip_if_not_installed("Glimma")
    skip_if_not(workflowOwnerCapable())
    fixture <- savedDaW1Owner("kit243-e4b")
    k243PageOrNote(fixture, function() k243ProtClick(fixture), "E4b")
})

test_that("KIT-243 E4m: proteomics managed, the interactive page or note DA1 did not write is removed by DA2 (S2 (a))", {
    skip_if_not_installed("Glimma")
    skip_if_not(workflowOwnerCapable())
    savedDaW2Mocks()
    fixture <- savedDaW1Owner("kit243-e4m")
    dispatch <- savedDaW2Dispatcher(fixture, NULL)
    k243PageOrNote(fixture, function() k243ProtWorkerClick(fixture, dispatch), "E4m")
})

# ------------------------------------------------------------------------------------------------
# Part B, E5 to E10 (failing first) with G2: Save, then a change, then Copy or Report without Save.
# ------------------------------------------------------------------------------------------------

for (case in list(list(row = "E5, E6, G2", class = "lipid"), list(row = "E8", class = "metab"))) {
    test_that(paste("KIT-243", case$row, ":", case$class, "DA, Save, a DA at another q: Copy and Report refuse until Save; then proceed"), {
        skip_if_not(workflowOwnerCapable())
        root <- withr::local_tempdir()
        seen <- k194Notes()
        omic <- k194Omic(case$class)
        q <- k207Normalised(file.path(root, "fixture"), case$class)
        owner <- k194Owner(q$fx)
        dirs <- k212Dirs(file.path(root, "project"), omic)
        da <- function(qvalue) k243Da(q$fx$workflow, dirs[[omic]], q$fx$config$da, q$fx$binding, kit172Contrasts(),
            settings = list(formula_string = "~ 0 + group", da_q_val_thresh = qvalue, treat_lfc_cutoff = 0))
        expect_identical(da(0.05)$status, "success", info = paste(case$row, "DA1 precondition"))
        rec <- k207Save(q$fx$workflow, omic, dirs)
        expect_identical(rec$done, 1L, info = paste(case$row, "Save precondition", paste(rec$notes, collapse = " | ")))
        expect_identical(da(0.1)$status, "success", info = paste(case$row, "DA2 precondition"))
        status <- k243Status(owner, dirs[[omic]]$source_dir)
        copy <- k212Copy(q$fx$workflow, omic, dirs, argsSaved = TRUE)
        report <- k212Report(q$fx$workflow, omic, dirs)
        k243Line(case$row, case$class, sprintf("status=%s", status), sprintf("copies=%d", copy$copies), sprintf("render=%d", report$render),
            sprintf("copy_notes=%s", paste(copy$notes, collapse = " | ")), sprintf("report_notes=%s", paste(report$notes, collapse = " | ")))
        k243ExpectCopyRefused(copy, k243Texts$stale, paste(case$row, "Copy"))
        expect_false(isTRUE(copy$value$fallbackWritten), info = paste(case$row, "fallback"))
        k212ExpectReportRefused(report, k243Texts$stale, paste(case$row, "Report"))
        if (exists("workflowSummarySaveStatus", mode = "function")) expect_identical(status, "changed", info = case$row)
        # G2: Save again, then Copy and Report proceed.
        rec <- k207Save(q$fx$workflow, omic, dirs)
        expect_identical(rec$done, 1L, info = paste("G2", case$class, "Save", paste(rec$notes, collapse = " | ")))
        copy <- k212Copy(q$fx$workflow, omic, dirs, argsSaved = TRUE)
        expect_identical(copy$copies, 1L, info = paste("G2", case$class, "Copy", paste(copy$notes, collapse = " | ")))
        report <- k212Report(q$fx$workflow, omic, dirs)
        expect_identical(report$render, 1L, info = paste("G2", case$class, "Report", paste(report$notes, collapse = " | ")))
        if (exists("workflowSummarySaveStatus", mode = "function")) {
            expect_identical(k243Status(owner, dirs[[omic]]$source_dir), "current", info = paste("G2", case$class))
        }
    })
}

test_that("KIT-243 E7: lipid, run 1, Apply, DA, Save, Reset, run 2, Apply, DA: Copy refuses until Save", {
    skip_if_not(workflowOwnerCapable())
    root <- withr::local_tempdir()
    seen <- k194Notes()
    q <- k207Normalised(file.path(root, "fixture"), "lipid")
    dirs <- k212Dirs(file.path(root, "project"), "lipidomics")
    da <- function() k243Da(q$fx$workflow, dirs$lipidomics, q$fx$config$da, q$fx$binding, kit172Contrasts())
    expect_identical(da()$status, "success", info = "E7 DA1 precondition")
    rec <- k207Save(q$fx$workflow, "lipidomics", dirs)
    expect_identical(rec$done, 1L, info = "E7 Save precondition")
    k194Reset(q$fx, q$nd)
    k194Run(q$fx, q$nd)
    k194ExpectRun(q$fx, q$nd, "KIT-243 E7 run 2")
    expect_true(isTRUE(k194Apply(q$fx, q$nd)), info = "E7 Apply 2 precondition")
    expect_identical(da()$status, "success", info = "E7 DA2 precondition")
    status <- k243Status(k194Owner(q$fx), dirs$lipidomics$source_dir)
    copy <- k212Copy(q$fx$workflow, "lipidomics", dirs, argsSaved = TRUE)
    k243Line("E7", "lipid", sprintf("status=%s", status), sprintf("copies=%d", copy$copies), sprintf("notes=%s", paste(copy$notes, collapse = " | ")))
    k243ExpectCopyRefused(copy, k243Texts$stale, "E7 Copy")
    if (exists("workflowSummarySaveStatus", mode = "function")) expect_identical(status, "changed", info = "E7")
})

test_that("KIT-243 E9: transcriptomics, DA, Save, a DA at another q: Copy and Report refuse", {
    skip_if_not_installed("edgeR")
    seen <- k194Notes()
    project <- module_ci_trx_project(module_ci_trx_corpus())
    owner <- project$owner
    k194TrxRun(project, omicNormReactiveState(), list(ruv_mode = "manual", ruv_k = 2L, ruv_percentage = 10), seen)
    expect_identical(owner$current_state, "trx_ruvs", info = "E9 run precondition")
    da <- function(qvalue) k243Da(project$workflow, project$paths, omicGuiConfig("transcriptomics")$da, transcriptomicsDeclarationBinding,
        project$workflow$contrasts_tbl, settings = list(formula_string = "~ 0 + group", da_q_val_thresh = qvalue, treat_lfc_cutoff = 0))
    expect_identical(da(0.05)$status, "success", info = paste("E9 DA1 precondition", paste(seen$notes, collapse = " | ")))
    rec <- k207Save(project$workflow, "transcriptomics", project$project_dirs)
    expect_identical(rec$done, 1L, info = paste("E9 Save precondition", paste(rec$notes, collapse = " | ")))
    expect_identical(da(0.1)$status, "success", info = "E9 DA2 precondition")
    status <- k243Status(owner, project$project_dirs$transcriptomics$source_dir)
    copy <- k212Copy(project$workflow, "transcriptomics", project$project_dirs, argsSaved = TRUE)
    report <- k212Report(project$workflow, "transcriptomics", project$project_dirs)
    k243Line("E9", "trx", sprintf("status=%s", status), sprintf("copies=%d", copy$copies), sprintf("render=%d", report$render),
        sprintf("copy_notes=%s", paste(copy$notes, collapse = " | ")))
    k243ExpectCopyRefused(copy, k243Texts$stale, "E9 Copy")
    k212ExpectReportRefused(report, k243Texts$stale, "E9 Report")
    if (exists("workflowSummarySaveStatus", mode = "function")) expect_identical(status, "changed", info = "E9")
})

test_that("KIT-243 E10: proteomics, DA, Save (real writers), a DA at another q: Copy and Report refuse", {
    skip_if_not(workflowOwnerCapable())
    fixture <- k212Prot()
    k212Export(fixture, fixture$norm)
    k243ProtDa(fixture, publish = TRUE)
    dirs <- k212ProtDirs(fixture, "e10")
    seen <- k243ProtSave(fixture, dirs)
    expect_length(seen$failed, 0L)
    expect_true(file.exists(file.path(dirs$proteomics$source_dir, "analysis_record.json")), info = paste("E10 record precondition",
        paste(seen$notes, collapse = " | ")))
    k243ProtDa(fixture, settings = list(formula_string = "~ 0 + group", da_q_val_thresh = 0.1, treat_lfc_cutoff = 0), publish = TRUE)
    status <- k243Status(fixture$manager, dirs$proteomics$source_dir)
    copy <- k212Copy(fixture$workflow, "proteomics", dirs, argsSaved = TRUE)
    report <- k212Report(fixture$workflow, "proteomics", dirs)
    k243Line("E10", "prot", sprintf("status=%s", status), sprintf("copies=%d", copy$copies), sprintf("render=%d", report$render),
        sprintf("copy_notes=%s", paste(copy$notes, collapse = " | ")))
    k243ExpectCopyRefused(copy, k243Texts$stale, "E10 Copy")
    k212ExpectReportRefused(report, k243Texts$stale, "E10 Report")
    if (exists("workflowSummarySaveStatus", mode = "function")) expect_identical(status, "changed", info = "E10")
})

# ------------------------------------------------------------------------------------------------
# Part D, E12 to E14 (failing first) with G14: a removal or a publication that did not complete.
# ------------------------------------------------------------------------------------------------

test_that("KIT-243 E12: proteomics in process, an earlier file that cannot be removed: the DA says so; Copy and Report refuse", {
    skip_if_not_installed("Glimma")
    skip_if_not(workflowOwnerCapable())
    seen <- k194Notes()
    fixture <- savedDaW1Owner("kit243-e12")
    expect_null(k243ProtClick(fixture), info = "E12 DA1 precondition")
    # R1-03: a non-empty directory under the barplot's name, which DA2 does not write (no barplot); a directory is never
    # deleted.
    blocked <- file.path(fixture$paths$publication_graphs_dir, "NumSigDaMolecules", "da_proteins_num_sig_da_molecules.png")
    unlink(blocked)
    dir.create(blocked, recursive = TRUE)
    writeLines("KIT-243", file.path(blocked, "inside.txt"))
    error <- k243ProtClick(fixture, qvalue = 1e-300)
    dirs <- k243ProtFixtureDirs(fixture)
    copy <- k212Copy(fixture$workflow, "proteomics", dirs)
    report <- k212Report(fixture$workflow, "proteomics", dirs)
    k243Line("E12", "prot", sprintf("error=%s", if (is.null(error)) "none" else paste(class(error)[[1L]], conditionMessage(error))),
        sprintf("marker=%s", file.exists(k243Marker(fixture$paths))), sprintf("copies=%d", copy$copies), sprintf("render=%d", report$render),
        sprintf("copy_notes=%s", paste(copy$notes, collapse = " | ")))
    expect_s3_class(error, "workflow_figure_unwritten")
    expect_true(!is.null(error) && grepl(k243Texts$removal, conditionMessage(error), fixed = TRUE) &&
        grepl("da_proteins_num_sig_da_molecules.png", conditionMessage(error), fixed = TRUE), info = "E12 the removal's text")
    k243ExpectCopyRefused(copy, k243Texts$pending, "E12 Copy")
    expect_false(file.exists(file.path(dirs$proteomics$source_dir, "study_parameters.txt")), info = "E12 no fallback file")
    k212ExpectReportRefused(report, k243Texts$pending, "E12 Report")
})

test_that("KIT-243 E13, G14: lipid in process, a publication that fails after the removal: Copy and Report refuse until the DA runs again", {
    skip_if_not(workflowOwnerCapable())
    root <- withr::local_tempdir()
    seen <- k194Notes()
    q <- k207Normalised(file.path(root, "fixture"), "lipid")
    owner <- k194Owner(q$fx)
    dirs <- k212Dirs(file.path(root, "project"), "lipidomics")
    paths <- dirs$lipidomics
    da <- function() k243Da(q$fx$workflow, paths, q$fx$config$da, q$fx$binding, kit172Contrasts())
    expect_identical(da()$status, "success", info = "E13 DA1 precondition")
    tables <- list.files(paths$da_output_dir, "_long_annot[.]xlsx$")
    heatmaps <- grep("/Heatmaps/", k243Made(owner, paths, "lipidomics"), value = TRUE, fixed = TRUE)
    k243Plant(heatmaps)
    # R1-04: the publication's handler stages a file, then stops (resolved by name from the package, so mocked there).
    second <- testthat::with_mocked_bindings(da(), computeWorkflowOmicsDaPublicationAction = function(manager, request, variant) {
        writeLines("staged", file.path(omicsStagedPublicationDir(request, "da_output_dir"), "kit243-staged.txt"))
        stop("KIT-243: the publication failed after staging", call. = FALSE)
    }, .package = "MultiScholaR")
    copy <- k212Copy(q$fx$workflow, "lipidomics", dirs)
    report <- k212Report(q$fx$workflow, "lipidomics", dirs)
    left <- k243PlantedOf(heatmaps)
    k243Line("E13", "lipid", sprintf("status=%s", second$status), sprintf("left=%d", length(left)),
        sprintf("tables=%s", paste(list.files(paths$da_output_dir, "_long_annot[.]xlsx$"), collapse = ",")),
        sprintf("marker=%s", file.exists(k243Marker(paths))), sprintf("copies=%d", copy$copies), sprintf("render=%d", report$render),
        sprintf("notes=%s", paste(seen$notes, collapse = " | ")))
    expect_identical(second$status, "success", info = "E13 DA2 (the analysis completes; its files do not)")
    expect_true(k194Noted(seen, "[warning] Warning: Could not save results to disk: KIT-243: the publication failed after staging"),
        info = paste("E13 the tab's warning", paste(seen$notes, collapse = " | ")))
    expect_true(length(tables) > 0L && all(file.exists(file.path(paths$da_output_dir, tables))), info = "E13 the earlier tables kept")
    expect_identical(left, character(), info = "E13 planted heatmaps left")
    k243ExpectCopyRefused(copy, k243Texts$pending, "E13 Copy")
    k212ExpectReportRefused(report, k243Texts$pending, "E13 Report")
    # G14: the DA run again with its publication, Save, Copy.
    expect_identical(da()$status, "success", info = "G14 DA3 precondition")
    rec <- k207Save(q$fx$workflow, "lipidomics", dirs)
    expect_identical(rec$done, 1L, info = paste("G14 Save", paste(rec$notes, collapse = " | ")))
    copy <- k212Copy(q$fx$workflow, "lipidomics", dirs, argsSaved = TRUE)
    k243Line("G14", "lipid", sprintf("marker=%s", file.exists(k243Marker(paths))), sprintf("copies=%d", copy$copies),
        sprintf("notes=%s", paste(copy$notes, collapse = " | ")))
    expect_false(file.exists(k243Marker(paths)), info = "G14 marker")
    expect_identical(copy$copies, 1L, info = paste("G14 Copy", paste(copy$notes, collapse = " | ")))
})

test_that("KIT-243 E14: proteomics managed, the worker's publication fails after the accept: earlier figures gone, Copy and Report refuse", {
    skip_if_not_installed("Glimma")
    skip_if_not(workflowOwnerCapable())
    seen <- k194Notes()
    fixture <- savedDaW1Owner("kit243-e14")
    expect_null(k243ProtClick(fixture), info = "E14 DA1 precondition")
    graphs <- fixture$paths$publication_graphs_dir
    figures <- setdiff(list.files(graphs, recursive = TRUE, full.names = TRUE), list.files(file.path(graphs, "Heatmap"), full.names = TRUE))
    expect_gt(length(figures), 0L)
    press <- k243ProtManagedPress(fixture, qvalue = 1e-300)
    left <- basename(figures[file.exists(figures)])
    dirs <- k243ProtFixtureDirs(fixture)
    copy <- k212Copy(fixture$workflow, "proteomics", dirs)
    report <- k212Report(fixture$workflow, "proteomics", dirs)
    k243Line("E14", "prot", sprintf("contrasts=%d", press$contrasts), sprintf("publication_requested=%d", press$applied),
        sprintf("figures=%d", length(figures)), sprintf("left=%s", paste(left, collapse = ",")), sprintf("marker=%s", file.exists(k243Marker(fixture$paths))),
        sprintf("copies=%d", copy$copies), sprintf("render=%d", report$render), sprintf("notes=%s", paste(seen$notes, collapse = " | ")))
    expect_identical(press$contrasts, 1L, info = "E14 precondition: the contrast ran")
    expect_identical(press$applied, 1L, info = "E14 precondition: the publication was requested")
    expect_identical(left, character(), info = "E14 earlier figures left")
    k243ExpectCopyRefused(copy, k243Texts$pending, "E14 Copy")
    k212ExpectReportRefused(report, k243Texts$pending, "E14 Report")
})

# ------------------------------------------------------------------------------------------------
# G rows: allowed, and unchanged (both commits; the new functions are called only where they exist).
# ------------------------------------------------------------------------------------------------

#' DA, Save, Copy and Report at once: Save done, Copy copies, Report renders, no marker; on head the Save "current".
k243ExpectDefault <- function(row, workflow, omic, dirs, manager, save = NULL) {
    if (is.null(save)) {
        rec <- k207Save(workflow, omic, dirs)
        expect_identical(rec$done, 1L, info = paste(row, omic, "Save", paste(rec$notes, collapse = " | ")))
    } else {
        expect_length(save$failed, 0L)
        expect_length(save$done, 1L)
    }
    status <- k243Status(manager, dirs[[omic]]$source_dir)
    copy <- k212Copy(workflow, omic, dirs, argsSaved = TRUE)
    report <- k212Report(workflow, omic, dirs)
    marker <- file.exists(k243Marker(dirs[[omic]]))
    k243Line(row, omic, sprintf("status=%s", status), sprintf("record=%s", file.exists(file.path(dirs[[omic]]$source_dir, "analysis_record.json"))),
        sprintf("copies=%d", copy$copies), sprintf("render=%d", report$render), sprintf("marker=%s", marker),
        sprintf("notes=%s", paste(c(copy$notes, report$notes), collapse = " | ")))
    expect_identical(copy$copies, 1L, info = paste(row, omic, "Copy", paste(copy$notes, collapse = " | ")))
    expect_identical(report$render, 1L, info = paste(row, omic, "Report", paste(report$notes, collapse = " | ")))
    expect_false(marker, info = paste(row, omic, "marker"))
    if (exists("workflowSummarySaveStatus", mode = "function")) expect_identical(status, "current", info = paste(row, omic))
}

for (case in list(list(class = "lipid", managed = FALSE), list(class = "metab", managed = FALSE), list(class = "lipid", managed = TRUE))) {
    test_that(paste("KIT-243 G1:", case$class, if (case$managed) "managed" else "in process", "DA, Save, Copy, Report at once"), {
        skip_if_not(workflowOwnerCapable())
        root <- withr::local_tempdir()
        seen <- k194Notes()
        omic <- k194Omic(case$class)
        q <- k207Normalised(file.path(root, "fixture"), case$class)
        dirs <- k212Dirs(file.path(root, "project"), omic)
        runner <- if (case$managed) k243Runner(k194Owner(q$fx), q$fx$binding, dirs[[omic]])
        out <- k243Da(q$fx$workflow, dirs[[omic]], q$fx$config$da, q$fx$binding, kit172Contrasts(), runner = runner$runner)
        expect_identical(out$status, if (case$managed) "dispatched" else "success", info = paste("G1 DA precondition", out$errorMessage %||% ""))
        if (case$managed) expect_true(paste0(q$fx$prefix, "_da_publication") %in% runner$operations, info = "G1 managed publication")
        k243ExpectDefault(paste("G1", if (case$managed) "managed" else "in-process"), q$fx$workflow, omic, dirs, k194Owner(q$fx))
    })
}

test_that("KIT-243 G1, G10: transcriptomics in process; a second DA rewrites every figure and removes nothing else", {
    skip_if_not_installed("edgeR")
    seen <- k194Notes()
    project <- module_ci_trx_project(module_ci_trx_corpus())
    owner <- project$owner
    k194TrxRun(project, omicNormReactiveState(), list(ruv_mode = "manual", ruv_k = 2L, ruv_percentage = 10), seen)
    expect_identical(owner$current_state, "trx_ruvs", info = "G1 trx run precondition")
    da <- function() k243Da(project$workflow, project$paths, omicGuiConfig("transcriptomics")$da, transcriptomicsDeclarationBinding,
        project$workflow$contrasts_tbl)
    expect_identical(da()$status, "success", info = paste("G1 trx DA precondition", paste(seen$notes, collapse = " | ")))
    k243ExpectDefault("G1 in-process", project$workflow, "transcriptomics", project$project_dirs, owner)
    # G10: every figure name the DA writes planted, and one file no DA writes; DA2.
    graphs <- project$paths$publication_graphs_dir
    figures <- Filter(function(path) startsWith(path, paste0(graphs, "/")), k243Made(owner, project$paths, "transcriptomics"))
    other <- file.path(graphs, "Volcano_Plots", "kit243_not_a_da_figure.png")
    k243Plant(c(figures, other))
    expect_identical(da()$status, "success", info = "G10 DA2 precondition")
    k243Line("G10", "trx", sprintf("figures=%d", length(figures)), sprintf("planted_left=%s", paste(k243PlantedOf(figures), collapse = ",")),
        sprintf("other=%s", k243Planted(other)))
    expect_gt(length(figures), 0L)
    expect_identical(k243PlantedOf(figures), character(), info = "G10 every figure rewritten")
    expect_true(k243Planted(other), info = "G10 the other file kept")
})

test_that("KIT-243 G1: proteomics in process (S1) and managed (S2 (a)), DA, Save (real writers), Copy, Report at once", {
    skip_if_not_installed("Glimma")
    skip_if_not(workflowOwnerCapable())
    savedDaW2Mocks()
    for (managed in c(FALSE, TRUE)) {
        fixture <- savedDaW1Owner(if (managed) "kit243-g1-managed" else "kit243-g1")
        error <- if (managed) k243ProtWorkerClick(fixture, savedDaW2Dispatcher(fixture, NULL)) else k243ProtClick(fixture)
        expect_null(error, info = paste("G1 prot DA precondition", managed))
        dirs <- k243ProtFixtureDirs(fixture)
        save <- k243ProtSave(fixture, dirs)
        k243ExpectDefault(paste("G1", if (managed) "managed" else "in-process"), fixture$workflow, "proteomics", dirs,
            fixture$workflow$state_manager, save = save)
    }
})

test_that("KIT-243 G3: lipid, a fresh project, DA, Copy with no Save ever: the fallback parameters and the copy, as today", {
    skip_if_not(workflowOwnerCapable())
    root <- withr::local_tempdir()
    seen <- k194Notes()
    q <- k207Normalised(file.path(root, "fixture"), "lipid")
    dirs <- k212Dirs(file.path(root, "project"), "lipidomics")
    expect_identical(k243Da(q$fx$workflow, dirs$lipidomics, q$fx$config$da, q$fx$binding, kit172Contrasts())$status, "success")
    status <- k243Status(k194Owner(q$fx), dirs$lipidomics$source_dir)
    copy <- k212Copy(q$fx$workflow, "lipidomics", dirs)
    k243Line("G3", "lipid", sprintf("status=%s", status), sprintf("copies=%d", copy$copies), sprintf("fallback=%s", copy$value$fallbackWritten))
    expect_identical(copy$copies, 1L, info = paste("G3 Copy", paste(copy$notes, collapse = " | ")))
    expect_true(isTRUE(copy$value$fallbackWritten), info = "G3 fallback")
    if (exists("workflowSummarySaveStatus", mode = "function")) expect_identical(status, "none", info = "G3")
})

test_that("KIT-243 G4: the module-CI in-memory proteomics manager still copies (KIT-212 G15)", {
    skip_if_not(workflowOwnerCapable())
    fake <- module_ci_prot_summary_workflow(module_ci_prot_summary_branch_config("DIA"))
    dirs <- module_ci_prot_summary_paths()
    withr::defer(unlink(dirs$proteomics$base_dir, recursive = TRUE))
    module_ci_prot_summary_write_publication_inputs(dirs)
    rec <- k212Copy(fake, "proteomics", dirs, argsSaved = TRUE)
    k243Line("G4", "prot", sprintf("status=%s", k243Status(fake$state_manager, dirs$proteomics$source_dir)), sprintf("copies=%d", rec$copies))
    expect_identical(rec$copies, 1L, info = paste("G4", paste(rec$notes, collapse = " | ")))
    if (exists("workflowSummarySaveStatus", mode = "function")) {
        expect_identical(k243Status(fake$state_manager, dirs$proteomics$source_dir), "none", info = "G4")
    }
})

test_that("KIT-243 G5, G6 (proteomics), G6b: the archive reopened; an identical DA re-run in the session; and after a reopen", {
    skip_if_not(workflowOwnerCapable())
    # G5: DA, the export, Save (real writers), the archive reopened, Copy.
    fixture <- k212Prot()
    k243ProtDa(fixture)
    k212Export(fixture, fixture$norm)
    dirs <- k212ProtDirs(fixture, "g5")
    save <- k243ProtSave(fixture, dirs)
    expect_identical(vapply(save$failed, conditionMessage, ""), character(), info = "G5 Save")
    archives <- sort(list.files(fixture$paths$source_dir, "^filtered_session_.*[.]msdia[.]tar$", full.names = TRUE))
    expect_true(length(archives) >= 1L, info = "G5 precondition: an archive")
    reopened <- k212Reopen(archives[[length(archives)]])
    status <- k243Status(reopened$manager, dirs$proteomics$source_dir)
    copy <- k212Copy(reopened$workflow, "proteomics", dirs, argsSaved = TRUE)
    k243Line("G5", "prot", sprintf("status=%s", status), sprintf("copies=%d", copy$copies), sprintf("notes=%s", paste(copy$notes, collapse = " | ")))
    expect_identical(copy$copies, 1L, info = paste("G5 Copy", paste(copy$notes, collapse = " | ")))
    if (exists("workflowSummarySaveStatus", mode = "function")) expect_identical(status, "current", info = "G5")
    # G6b (an observation, not a stop): the identical DA run again on the reopened owner, then Copy.
    k243ProtDa(list(manager = reopened$manager, workflow = reopened$workflow, paths = reopened$paths))
    status <- k243Status(reopened$manager, dirs$proteomics$source_dir)
    copy <- k212Copy(reopened$workflow, "proteomics", dirs, argsSaved = TRUE)
    k243Line("G6b", "prot", sprintf("status=%s", status), sprintf("copies=%d", copy$copies), sprintf("notes=%s", paste(copy$notes, collapse = " | ")))
    # G6: DA, Save, the identical DA in the same session, Copy.
    fixture <- k212Prot()
    # Instrument repair (R1 45590602596, VAL-54): the proteomics Save follows the tab's export, as in every KIT-212 Save.
    k212Export(fixture, fixture$norm)
    k243ProtDa(fixture)
    dirs <- k212ProtDirs(fixture, "g6")
    save <- k243ProtSave(fixture, dirs)
    expect_identical(vapply(save$failed, conditionMessage, ""), character(), info = "G6 prot Save")
    k243ProtDa(fixture)
    status <- k243Status(fixture$manager, dirs$proteomics$source_dir)
    copy <- k212Copy(fixture$workflow, "proteomics", dirs, argsSaved = TRUE)
    k243Line("G6", "prot", sprintf("status=%s", status), sprintf("copies=%d", copy$copies))
    expect_identical(copy$copies, 1L, info = paste("G6 prot Copy", paste(copy$notes, collapse = " | ")))
    if (exists("workflowSummarySaveStatus", mode = "function")) expect_identical(status, "current", info = "G6 prot")
})

test_that("KIT-243 G6: lipid, DA, Save, the identical DA in the same session: Copy proceeds", {
    skip_if_not(workflowOwnerCapable())
    root <- withr::local_tempdir()
    seen <- k194Notes()
    q <- k207Normalised(file.path(root, "fixture"), "lipid")
    dirs <- k212Dirs(file.path(root, "project"), "lipidomics")
    da <- function() k243Da(q$fx$workflow, dirs$lipidomics, q$fx$config$da, q$fx$binding, kit172Contrasts())
    expect_identical(da()$status, "success")
    rec <- k207Save(q$fx$workflow, "lipidomics", dirs)
    expect_identical(rec$done, 1L)
    expect_identical(da()$status, "success")
    status <- k243Status(k194Owner(q$fx), dirs$lipidomics$source_dir)
    copy <- k212Copy(q$fx$workflow, "lipidomics", dirs, argsSaved = TRUE)
    k243Line("G6", "lipid", sprintf("status=%s", status), sprintf("copies=%d", copy$copies))
    expect_identical(copy$copies, 1L, info = paste("G6 lipid Copy", paste(copy$notes, collapse = " | ")))
    if (exists("workflowSummarySaveStatus", mode = "function")) expect_identical(status, "current", info = "G6 lipid")
})

test_that("KIT-243 G9: lipid, a DA of A and B, then of A only: Copy refuses naming B's table; B's figures not removed", {
    skip_if_not(workflowOwnerCapable())
    root <- withr::local_tempdir()
    seen <- k194Notes()
    q <- k207Normalised(file.path(root, "fixture"), "lipid")
    owner <- k194Owner(q$fx)
    dirs <- k212Dirs(file.path(root, "project"), "lipidomics")
    both <- data.frame(contrasts = c("groupB-groupA", "groupA-groupB"), friendly_names = c("B_vs_A", "A_vs_B"),
        full_format = c("B_vs_A=groupB-groupA", "A_vs_B=groupA-groupB"), stringsAsFactors = FALSE)
    da <- function(contrasts) k243Da(q$fx$workflow, dirs$lipidomics, q$fx$config$da, q$fx$binding, contrasts)
    expect_identical(da(both)$status, "success")
    b_heatmaps <- grep("groupA-groupB_heatmap[.](png|pdf)$", k243Made(owner, dirs$lipidomics, "lipidomics"), value = TRUE)
    k243Plant(b_heatmaps)
    expect_identical(da(both[1L, , drop = FALSE])$status, "success")
    copy <- k212Copy(q$fx$workflow, "lipidomics", dirs)
    table_b <- "de_posmode_lipids_groupA-groupB_long_annot.xlsx"
    k243Line("G9", "lipid", sprintf("b_heatmaps=%s", paste(k243PlantedOf(b_heatmaps), collapse = ",")), sprintf("copies=%d", copy$copies),
        sprintf("notes=%s", paste(copy$notes, collapse = " | ")))
    expect_gt(length(b_heatmaps), 0L)
    expect_identical(k243PlantedOf(b_heatmaps), basename(b_heatmaps), info = "G9 B's figures kept")
    k243ExpectCopyRefused(copy, table_b, "G9 Copy")
})

test_that("KIT-243 G11: replay after a published DA writes no marker", {
    skip_if_not(workflowOwnerCapable())
    root <- withr::local_tempdir()
    seen <- k194Notes()
    fx <- k194Fixture(file.path(root, "project"), "lipid")
    k194Defer(fx)
    k212Journey(fx)
    paths <- k212Dirs(file.path(root, "publication"), "lipidomics")$lipidomics
    expect_identical(k243Da(fx$workflow, paths, fx$config$da, fx$binding, kit172Contrasts())$status, "success")
    out <- k212Replay(fx, file.path(root, "replay"))
    markers <- list.files(root, "^[.]multischolar_da_publication_pending$", recursive = TRUE, all.files = TRUE)
    k243Line("G11", "lipid", sprintf("error=%s", if (is.null(out$error)) "none" else conditionMessage(out$error)),
        sprintf("markers=%s", paste(markers, collapse = ",")))
    expect_null(out$error, info = "G11 replay")
    expect_identical(markers, character(), info = "G11 markers")
})
