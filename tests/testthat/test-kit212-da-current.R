# KIT-212 (VAL-172): the Summary's Save, Report and Copy refuse, and the lipid, metab and trx DA tab clears, when the
# differential analysis was not run on the current table, or when the DA folders hold DA files no DA of this session
# made; replay refuses a record whose DA was run on another table. Written before the fix: on the tests-only commit
# each E row runs as before, so it fails at the value the ticket states (section 10); the G rows pass on both commits.
# The new functions are called directly only where they exist, after the base assertions. Lines `KIT-212 <row>` are the
# evidence the ticket's records copy (outcome.md).

library(testthat)

# ------------------------------------------------------------------------------------------------
# Copied unchanged from tests/testthat/test-kit207-summary-head-guard.R at 32766889 (KIT-207 E-7): lines 16-256 (KIT-194's
# helpers, themselves copied from test-omics-norm-state-guards.R), 280-309 (k207Dirs, k207IntegrationDir, k207Save) and
# 353-373 (k207Da, k207Normalised, k207ExpectHead). That file is edited only at E14 (the fix commit).
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

# ------------------------------------------------------------------------------------------------
# E1 to E4: lipidomics and metabolomics, the table changed after the DA (failing first).
# ------------------------------------------------------------------------------------------------

for (class in c("lipid", "metab")) {
    test_that(paste("KIT-212 E1, E2, E3:", class, "Q1 (a second run without a DA): Save, Copy and Report refuse"), {
        skip_if_not(workflowOwnerCapable())
        root <- withr::local_tempdir()
        seen <- k194Notes()
        omic <- k194Omic(class)
        q <- k212Q1(file.path(root, "q1"), class)
        owner <- k194Owner(q$fx)
        cat(sprintf("KIT-212 Q1 %s status=%s\n", class, k212Status(owner)))
        # E1: Save, an integration input planted first.
        dirs <- k212Dirs(file.path(root, "e1-project"), omic)
        input <- file.path(k207IntegrationDir(dirs[[omic]]), omic)
        dir.create(input, recursive = TRUE, showWarnings = FALSE)
        writeLines("planted", file.path(input, "k212-planted.txt"))
        rec <- k207Save(q$fx$workflow, omic, dirs)
        params <- file.path(dirs[[omic]]$source_dir, "study_parameters.txt")
        record <- file.path(dirs[[omic]]$source_dir, "analysis_record.json")
        source_state <- if (file.exists(record)) format(jsonlite::read_json(record)$analyses$differential_abundance$source_state %||% "<null>") else "<no record>"
        cat(sprintf("KIT-212 E1 %s\tdone=%d\tfailed=%s\tparams=%s\tsource_state=%s\tplanted=%s\tnotes=%s\n", class, rec$done,
            paste(vapply(rec$failed, function(e) class(e)[[1L]], ""), collapse = ","), file.exists(params), source_state,
            file.exists(file.path(input, "k212-planted.txt")), paste(rec$notes, collapse = " | ")))
        k212ExpectSaveRefused(rec, dirs, omic, paste("E1", class))
        expect_false(dir.exists(input), info = paste("E1", class, "the planted input"))
        # E2: Copy (the parameters not saved, so the base writes its fallback file).
        dirs <- k212Dirs(file.path(root, "e2-project"), omic)
        rec <- k212Copy(q$fx$workflow, omic, dirs)
        k212Line("E2", class, rec, sprintf("\tcopies=%d", rec$copies))
        k212ExpectCopyRefused(rec, dirs, omic, k212Texts$stale, paste("E2", class))
        # E3: Report.
        dirs <- k212Dirs(file.path(root, "e3-project"), omic)
        rec <- k212Report(q$fx$workflow, omic, dirs)
        k212Line("E3", class, rec, sprintf("\trender=%d\tbundle=%d", rec$render, rec$bundle))
        k212ExpectReportRefused(rec, k212Texts$stale, paste("E3", class))
        if (exists("workflowDaResultStatus", mode = "function")) expect_identical(workflowDaResultStatus(owner), "stale")
    })

    test_that(paste("KIT-212 E4:", class, "Q2 (one Revert) and Q3 (Apply again): Save refuses"), {
        skip_if_not(workflowOwnerCapable())
        root <- withr::local_tempdir()
        seen <- k194Notes()
        omic <- k194Omic(class)
        prefix <- kit172Binding(class)$prefix
        # Q2: DA at the correlation result, one Revert (the head the run's `_normalized`).
        q <- k207Normalised(file.path(root, "q2"), class, da = TRUE)
        k194Revert(q$fx, times = 1L)
        k207ExpectHead(q$fx, paste0(prefix, "_normalized"), paste("Q2", class))
        dirs <- k212Dirs(file.path(root, "e4-q2"), omic)
        rec <- k207Save(q$fx$workflow, omic, dirs)
        k212Line("E4-Q2", class, rec, sprintf("\tdone=%d", rec$done))
        k212ExpectSaveRefused(rec, dirs, omic, paste("E4 Q2", class))
        # Q3: DA, Apply again with another threshold.
        q <- k207Normalised(file.path(root, "q3"), class, da = TRUE)
        expect_true(isTRUE(k194Apply(q$fx, q$nd, threshold = 0.6)), info = paste("Q3 Apply precondition", class))
        k207ExpectHead(q$fx, paste0(prefix, "_correlation_filtered"), paste("Q3", class))
        dirs <- k212Dirs(file.path(root, "e4-q3"), omic)
        rec <- k207Save(q$fx$workflow, omic, dirs)
        k212Line("E4-Q3", class, rec, sprintf("\tdone=%d", rec$done))
        k212ExpectSaveRefused(rec, dirs, omic, paste("E4 Q3", class))
    })
}

# ------------------------------------------------------------------------------------------------
# E5: transcriptomics; E6: proteomics; E7: a QC Revert (failing first).
# ------------------------------------------------------------------------------------------------

test_that("KIT-212 E5: transcriptomics, a second run without a DA: Save and Copy refuse", {
    skip_if_not_installed("edgeR")
    seen <- new.env(parent = emptyenv())
    seen$notes <- character()
    project <- module_ci_trx_project(module_ci_trx_corpus())
    owner <- project$owner
    nd <- omicNormReactiveState()
    k194TrxRun(project, nd, list(ruv_mode = "manual", ruv_k = 2L, ruv_percentage = 10), seen)
    expect_identical(owner$current_state, "trx_ruvs", info = paste("E5 run precondition:", paste(seen$notes, collapse = " | ")))
    kit172Quiet(runTrxReferenceDaAnalysis(owner, getWorkflowStateReference(owner), list(da_q_val_thresh = 0.05, treat_lfc_cutoff = 0),
        trxDaContrastEntries(project$workflow$contrasts_tbl)))
    k194TrxRun(project, nd, list(ruv_mode = "manual", ruv_k = 1L, ruv_percentage = 10), seen)
    expect_identical(owner$current_state, "trx_ruvs", info = "E5 second run precondition")
    dirs <- k212Dirs(file.path(project$paths$base_dir, "k212-e5-save"), "transcriptomics")
    rec <- k207Save(project$workflow, "transcriptomics", dirs)
    k212Line("E5-save", "trx", rec, sprintf("\tdone=%d", rec$done))
    k212ExpectSaveRefused(rec, dirs, "transcriptomics", "E5 Save")
    dirs <- k212Dirs(file.path(project$paths$base_dir, "k212-e5-copy"), "transcriptomics")
    rec <- k212Copy(project$workflow, "transcriptomics", dirs)
    k212Line("E5-copy", "trx", rec, sprintf("\tcopies=%d", rec$copies))
    k212ExpectCopyRefused(rec, dirs, "transcriptomics", k212Texts$stale, "E5 Copy")
})

test_that("KIT-212 E6: proteomics Q6 (a QC Revert after the DA): Copy, Save and Report refuse", {
    skip_if_not(workflowOwnerCapable())
    fixture <- k212Prot()
    manager <- fixture$manager
    k212ProtDa(fixture)
    protNormCoreQuiet(withDIAStateScope(manager, function() selectDIAStateReference(manager, "protein_replicate_filtered")))
    expect_identical(manager$current_state, "protein_replicate_filtered", info = "E6 head precondition")
    cat(sprintf("KIT-212 E6 status=%s\n",
        k212Status(manager)))
    dirs <- k212ProtDirs(fixture, "e6-copy")
    rec <- k212Copy(fixture$workflow, "proteomics", dirs)
    k212Line("E6-copy", "prot", rec, sprintf("\tcopies=%d", rec$copies))
    k212ExpectCopyRefused(rec, dirs, "proteomics", k212Texts$stale, "E6 Copy")
    dirs <- k212ProtDirs(fixture, "e6-save")
    rec <- k207Save(fixture$workflow, "proteomics", dirs)
    k212Line("E6-save", "prot", rec, sprintf("\tdone=%d\tfailed=%s\tparams=%s", rec$done,
        paste(vapply(rec$failed, function(e) class(e)[[1L]], ""), collapse = ","),
        file.exists(file.path(dirs$proteomics$source_dir, "study_parameters.txt"))))
    info <- paste("E6 Save", paste(rec$notes, collapse = " | "))
    expect_true(length(rec$failed) == 1L && inherits(rec$failed[[1L]], c("workflow_summary_da_stale", "omic_norm_correlation_input_stale")),
        info = info)
    expect_false(file.exists(file.path(dirs$proteomics$source_dir, "study_parameters.txt")), info = info)
    dirs <- k212ProtDirs(fixture, "e6-report")
    rec <- k212Report(fixture$workflow, "proteomics", dirs)
    k212Line("E6-report", "prot", rec, sprintf("\trender=%d", rec$render))
    expect_identical(rec$render, 0L, info = paste("E6 Report", paste(rec$notes, collapse = " | ")))
    expect_true(any(startsWith(rec$notes, "[error] Report generation failed:")), info = paste("E6 Report", paste(rec$notes, collapse = " | ")))
})

test_that("KIT-212 E7: lipid, the DA, then Quality Control reverted to QC complete: Copy refuses", {
    skip_if_not(workflowOwnerCapable())
    root <- withr::local_tempdir()
    seen <- k194Notes()
    q <- k207Normalised(file.path(root, "e7"), "lipid", da = TRUE)
    expect_identical(k194Revert(q$fx, "lipid_qc_complete"), "lipid_qc_complete", info = "E7 precondition")
    dirs <- k212Dirs(file.path(root, "e7-project"), "lipidomics")
    rec <- k212Copy(q$fx$workflow, "lipidomics", dirs)
    k212Line("E7", "lipid", rec, sprintf("\tcopies=%d", rec$copies))
    k212ExpectCopyRefused(rec, dirs, "lipidomics", k212Texts$stale, "E7")
})

# ------------------------------------------------------------------------------------------------
# E8, E10 (failing first) and G14: the lipid DA tab on a reactive bus.
# ------------------------------------------------------------------------------------------------

test_that("KIT-212 E8: lipid, the tab's results are cleared when a second run moves the head after its DA", {
    skip_if_not(workflowOwnerCapable())
    root <- withr::local_tempdir()
    seen <- k194Notes()
    q <- k207Normalised(file.path(root, "e8"), "lipid")
    out <- new.env(parent = emptyenv())
    k212Tab(q$fx, function(tab, run_da, fire, bus, session) {
        run_da()
        out$after_da <- !is.null(shiny::isolate(tab$da_results_list))
        k194Reset(q$fx, q$nd)
        k194Run(q$fx, q$nd)
        k194Apply(q$fx, q$nd)
        fire()
        out$results <- !is.null(shiny::isolate(tab$da_results_list))
        out$complete <- shiny::isolate(tab$analysis_complete)
        out$heatmap <- !is.null(shiny::isolate(tab$current_heatmap_plot))
        out$clusters <- !is.null(shiny::isolate(tab$current_row_clusters)) || !is.null(shiny::isolate(tab$current_col_clusters))
        out$status <- shiny::isolate(bus$tab_status$differential_analysis)
    })
    cat(sprintf("KIT-212 E8 lipid\tafter_da=%s\tresults=%s\tcomplete=%s\theatmap=%s\tstatus=%s\thead=%s\n", out$after_da,
        out$results, format(out$complete), out$heatmap, format(out$status), k194Head(q$fx)))
    expect_true(isTRUE(out$after_da), info = "E8 precondition: the DA's results shown")
    expect_false(isTRUE(out$results), info = "E8 results")
    expect_identical(out$complete, FALSE, info = "E8 analysis_complete")
    expect_false(isTRUE(out$heatmap), info = "E8 heatmap")
    expect_false(isTRUE(out$clusters), info = "E8 clusters")
    expect_identical(out$status, "pending", info = "E8 status")
})

test_that("KIT-212 E10: lipid, a loaded state name is forgotten once the head moves, and the DA runs on the head", {
    skip_if_not(workflowOwnerCapable())
    root <- withr::local_tempdir()
    seen <- k194Notes()
    q <- k207Normalised(file.path(root, "e10"), "lipid")
    owner <- k194Owner(q$fx)
    loaded <- "lipid_correlation_filtered"
    out <- new.env(parent = emptyenv())
    k212Tab(q$fx, function(tab, run_da, fire, bus, session) {
        # Load Filtered Session of the owner's own export names the head (`tab$current_state_name`,
        # R/utils_workflow_omics_da.R:109); its effect on the tab is set here, the export and its reopening not driven.
        tab$current_state_name <- loaded
        protNormCoreQuiet(selectWorkflowStateReference(owner, "lipid_normalized"))
        fire()
        out$name <- shiny::isolate(tab$current_state_name)
        run_da()
    })
    manifest <- withWorkflowStateScope(owner, function() resolveWorkflowStateValue(owner$.context, owner$.context$snapshot$da_results))
    source <- names(owner$.context$snapshot$states)[vapply(owner$.context$snapshot$states, function(reference)
        identical(reference$entity_id, manifest$source$entity_id), logical(1))]
    dirs <- k207Dirs(file.path(root, "e10-project"), "lipidomics")
    rec <- k207Save(q$fx$workflow, "lipidomics", dirs)
    cat(sprintf("KIT-212 E10 lipid\thead=%s\tname_after=%s\tda_source=%s\tsave_done=%d\tnotes=%s\n", k194Head(q$fx),
        format(out$name %||% "<null>"), paste(source, collapse = ","), rec$done, paste(rec$notes, collapse = " | ")))
    expect_identical(k194Head(q$fx), "lipid_normalized", info = "E10 head precondition")
    expect_null(out$name, info = "E10 the loaded name")
    expect_identical(manifest$source$entity_id, getWorkflowStateReference(owner, "lipid_normalized")$entity_id, info = "E10 the DA's source")
    expect_identical(rec$done, 1L, info = paste("E10 Save", paste(rec$notes, collapse = " | ")))
})

test_that("KIT-212 G14: lipid, the tab with the owner attached keeps the results its own DA published", {
    skip_if_not(workflowOwnerCapable())
    root <- withr::local_tempdir()
    seen <- k194Notes()
    q <- k207Normalised(file.path(root, "g14"), "lipid")
    out <- new.env(parent = emptyenv())
    k212Tab(q$fx, function(tab, run_da, fire, bus, session) {
        attachDIAStateSession(k194Owner(q$fx), session, bus, omic_type = "lipidomics")
        run_da()
        session$flushReact()
        out$fired <- !is.null(shiny::isolate(bus$state_update_trigger))
        out$results <- !is.null(shiny::isolate(tab$da_results_list))
        out$status <- shiny::isolate(bus$tab_status$differential_analysis)
    })
    cat(sprintf("KIT-212 G14 lipid\tfired=%s\tresults=%s\tstatus=%s\n", out$fired, out$results, format(out$status)))
    expect_true(isTRUE(out$fired), info = "G14 precondition: the DA's publication fired the trigger")
    expect_true(isTRUE(out$results), info = "G14 results")
    expect_identical(out$status, "complete", info = "G14 status")
})

# ------------------------------------------------------------------------------------------------
# E9, E9b, E9c (failing first), G11, G12: replay (lipid; the declarations journey of KIT-194 T10).
# ------------------------------------------------------------------------------------------------

test_that("KIT-212 E9, E9b, E9c: replay refuses a record whose DA was run on another table than its final one", {
    skip_if_not(workflowOwnerCapable())
    root <- withr::local_tempdir()
    # E9: a second run after a Revert to QC complete, no DA (Q1's record shape: the DA's source name rebound, its
    # `source_state` NA; Reset has no replay entry, R/utils_workflow_omics_declarations.R:2110-2145).
    fx <- k194Fixture(file.path(root, "e9", "project"), "lipid")
    k194Defer(fx)
    owner <- k194Owner(fx)
    k212Journey(fx)
    k207Da(fx)
    protNormCoreQuiet(selectWorkflowStateReference(owner, "lipid_qc_complete"))
    k212Journey(fx)
    out <- k212Replay(fx, file.path(root, "e9"))
    recorded <- out$record$analyses$differential_abundance$source_state
    expect_true(is.null(recorded) || isTRUE(is.na(recorded)) || identical(recorded, "NA"), info = "E9 precondition: source_state NA")
    expect_s3_class(out$error, "workflow_replay_invalid")
    expect_identical(if (is.null(out$error)) NA_character_ else conditionMessage(out$error), k212Texts$replay)
    # E9b: one Revert after the DA (`source_state` named, outside the steps).
    fx <- k194Fixture(file.path(root, "e9b", "project"), "lipid")
    k194Defer(fx)
    owner <- k194Owner(fx)
    k212Journey(fx)
    k207Da(fx)
    protNormCoreQuiet(selectWorkflowStateReference(owner, "lipid_normalized"))
    out <- k212Replay(fx, file.path(root, "e9b"))
    expect_identical(out$record$analyses$differential_abundance$source_state, "lipid_correlation_filtered", info = "E9b precondition")
    expect_s3_class(out$error, "workflow_replay_invalid")
    # E9c: the DA on the log2 state, then the normalisation (the source among the steps, another table at the head).
    fx <- k194Fixture(file.path(root, "e9c", "project"), "lipid")
    k194Defer(fx)
    owner <- k194Owner(fx)
    kit172Quiet(runWorkflowOmicsReferenceLog2(fx$binding, owner, 0, fx$workflow$config_list))
    k207Da(fx)
    kit172Quiet(runWorkflowOmicsReferenceBetweenSampleNormalisation(fx$binding, owner, "quantile", fx$workflow$config_list))
    out <- k212Replay(fx, file.path(root, "e9c"))
    expect_identical(out$record$analyses$differential_abundance$source_state, "lipid_log2", info = "E9c precondition")
    expect_s3_class(out$error, "workflow_replay_invalid")
})

test_that("KIT-212 G11, G12: after a Revert past the DA's source (and an identical re-apply) Save refuses; the DA re-run saves and replays", {
    skip_if_not(workflowOwnerCapable())
    root <- withr::local_tempdir()
    seen <- k194Notes()
    for (route in c("G11", "G12")) {
        fx <- k194Fixture(file.path(root, route, "project"), "lipid")
        k194Defer(fx)
        owner <- k194Owner(fx)
        k212Journey(fx)
        k207Da(fx)
        before <- withWorkflowStateScope(owner, function() readWorkflowStateMetadata(owner)$data$components)
        protNormCoreQuiet(selectWorkflowStateReference(owner, "lipid_normalized"))
        # G11 (L7): the correlation filter applied again with the same threshold; G12 (L8): the head stays `_normalized`.
        if (identical(route, "G11")) kit172Quiet(runWorkflowOmicsReferenceCorrelation(fx$binding, owner, "group", 0.5, fx$workflow$config_list))
        same <- identical(before, withWorkflowStateScope(owner, function() readWorkflowStateMetadata(owner)$data$components))
        dirs <- k212Dirs(file.path(root, paste0(route, "-first")), "lipidomics")
        rec <- k207Save(fx$workflow, "lipidomics", dirs)
        cat(sprintf("KIT-212 %s lipid\thead=%s\tsame_table=%s\tfirst_save_done=%d\tnotes=%s\n", route, k194Head(fx), same, rec$done,
            paste(rec$notes, collapse = " | ")))
        k212ExpectSaveRefused(rec, dirs, "lipidomics", paste(route, "first Save"))
        k207Da(fx)
        dirs <- k212Dirs(file.path(root, paste0(route, "-second")), "lipidomics")
        rec <- k207Save(fx$workflow, "lipidomics", dirs)
        expect_identical(rec$done, 1L, info = paste(route, "second Save", paste(rec$notes, collapse = " | ")))
        out <- k212Replay(fx, file.path(root, paste0(route, "-replay")))
        expect_null(out$error, info = paste(route, "replay"))
    }
})

# ------------------------------------------------------------------------------------------------
# E11 to E13: DA files no DA of this session made (failing first).
# ------------------------------------------------------------------------------------------------

test_that("KIT-212 E11: lipid, a new session's owner on a project holding an earlier session's DA tables: Copy and Report refuse", {
    skip_if_not(workflowOwnerCapable())
    workflow <- module_ci_lipid_summary_workflow()
    dirs <- module_ci_lipid_summary_paths()
    withr::defer(unlink(dirs$lipidomics$base_dir, recursive = TRUE))
    module_ci_lipid_summary_write_publication_inputs(dirs)
    table <- "de_posmode_lipids_groupKO-groupWT_long_annot.xlsx"
    rec <- k212Copy(workflow, "lipidomics", dirs)
    k212Line("E11-copy", "lipid", rec, sprintf("\tcopies=%d\tunmade=%s", rec$copies,
        paste(k212Unmade(workflow$state_manager, dirs$lipidomics, "lipidomics"), collapse = ",")))
    info <- paste("E11 Copy", paste(rec$notes, collapse = " | "))
    expect_true(any(startsWith(rec$notes, "[error] Copy error:") & grepl(k212Texts$unmade_lead, rec$notes, fixed = TRUE) &
        grepl(table, rec$notes, fixed = TRUE)), info = info)
    expect_identical(rec$copies, 0L, info = info)
    rec <- k212Report(workflow, "lipidomics", dirs)
    k212Line("E11-report", "lipid", rec, sprintf("\trender=%d", rec$render))
    k212ExpectReportRefused(rec, k212Texts$unmade_lead, "E11 Report")
})

test_that("KIT-212 E12: lipid, a DA of contrasts A and B, then of A only: Copy refuses and names B's table", {
    skip_if_not(workflowOwnerCapable())
    root <- withr::local_tempdir()
    seen <- k194Notes()
    q <- k207Normalised(file.path(root, "e12"), "lipid")
    owner <- k194Owner(q$fx)
    dirs <- k212Dirs(file.path(root, "e12-project"), "lipidomics")
    both <- data.frame(contrasts = c("groupB-groupA", "groupA-groupB"), friendly_names = c("B_vs_A", "A_vs_B"),
        full_format = c("B_vs_A=groupB-groupA", "A_vs_B=groupA-groupB"), stringsAsFactors = FALSE)
    kit172Quiet(runWorkflowOmicsReferenceDaAnalysis(q$fx$binding, owner, k194Ref(q$fx), k212Settings, both))
    k212Publish(q$fx$binding, owner, dirs$lipidomics)
    kit172Quiet(runWorkflowOmicsReferenceDaAnalysis(q$fx$binding, owner, k194Ref(q$fx), k212Settings, both[1L, , drop = FALSE]))
    k212Publish(q$fx$binding, owner, dirs$lipidomics)
    table_b <- "de_posmode_lipids_groupA-groupB_long_annot.xlsx"
    on_disk <- list.files(dirs$lipidomics$da_output_dir)
    figures <- list.files(file.path(dirs$lipidomics$publication_graphs_dir, "Heatmaps"))
    rec <- k212Copy(q$fx$workflow, "lipidomics", dirs)
    unmade <- k212Unmade(owner, dirs$lipidomics, "lipidomics")
    k212Line("E12", "lipid", rec, sprintf("\tcopies=%d\tda_files=%s\theatmaps=%s\tunmade=%s", rec$copies, paste(on_disk, collapse = ","),
        paste(figures, collapse = ","), paste(unmade, collapse = ",")))
    expect_true(table_b %in% on_disk, info = "E12 precondition: B's table from the first publication")
    info <- paste("E12 Copy", paste(rec$notes, collapse = " | "))
    expect_true(any(startsWith(rec$notes, "[error] Copy error:") & grepl(table_b, rec$notes, fixed = TRUE)), info = info)
    expect_identical(rec$copies, 0L, info = info)
    if (exists("workflowDaFilesUnmade", mode = "function")) {
        # The unmade files are B's: its table, and its heatmaps where the publication drew them (lipid Copy collects
        # `Heatmaps`, not `Volcano_Plots`).
        expect_setequal(unmade, c(table_b, grep("groupA-groupB_heatmap[.](png|pdf)$", figures, value = TRUE)))
    }
})

test_that("KIT-212 E13: proteomics, a session owner without a DA on a project holding DA files: Copy and Report refuse", {
    skip_if_not(workflowOwnerCapable())
    w <- k212ProtW1Owner()
    rec <- k212Copy(w$bus, "proteomics", w$project_dirs)
    k212Line("E13-copy", "prot", rec, sprintf("\tcopies=%d\tunmade=%s", rec$copies,
        paste(k212Unmade(w$manager, w$project_dirs$proteomics, "proteomics"), collapse = ",")))
    info <- paste("E13 Copy", paste(rec$notes, collapse = " | "))
    expect_true(any(startsWith(rec$notes, "[error] Copy error:") & grepl(k212Texts$unmade_lead, rec$notes, fixed = TRUE) &
        grepl("da_groupB-groupA_long_annot.xlsx", rec$notes, fixed = TRUE)), info = info)
    expect_identical(rec$copies, 0L, info = info)
    rec <- k212Report(w$bus, "proteomics", w$project_dirs)
    k212Line("E13-report", "prot", rec, sprintf("\trender=%d", rec$render))
    k212ExpectReportRefused(rec, k212Texts$unmade_lead, "E13 Report")
})

# ------------------------------------------------------------------------------------------------
# G rows: allowed, and unchanged (both commits; the new functions are called only where they exist).
# ------------------------------------------------------------------------------------------------

#' The default path's Summary on both commits: Save done, Copy copies, Report renders; on head the DA "current" and no
#' unmade file.
k212ExpectDefault <- function(row, workflow, omic, dirs, manager, save = TRUE) {
    if (save) {
        rec <- k207Save(workflow, omic, dirs)
        expect_identical(rec$done, 1L, info = paste(row, omic, "Save", paste(rec$notes, collapse = " | ")))
    }
    rec <- k212Copy(workflow, omic, dirs, argsSaved = TRUE)
    unmade <- k212Unmade(manager, dirs[[omic]], omic)
    files <- c(list.files(dirs[[omic]]$da_output_dir), list.files(dirs[[omic]]$publication_graphs_dir, recursive = TRUE))
    cat(sprintf("KIT-212 %s %s\tstatus=%s\tcopies=%d\tunmade=%s\tfiles=%s\tnotes=%s\n", row, omic, k212Status(manager), rec$copies,
        paste(unmade, collapse = ","), paste(files, collapse = ","), paste(rec$notes, collapse = " | ")))
    expect_identical(rec$copies, 1L, info = paste(row, omic, "Copy", paste(rec$notes, collapse = " | ")))
    rec <- k212Report(workflow, omic, dirs)
    expect_identical(rec$render, 1L, info = paste(row, omic, "Report", paste(rec$notes, collapse = " | ")))
    if (exists("workflowDaResultStatus", mode = "function")) {
        expect_identical(workflowDaResultStatus(manager), "current", info = paste(row, omic))
        expect_identical(unmade, character(), info = paste(row, omic, "unmade"))
    }
}

for (class in c("lipid", "metab")) {
    test_that(paste("KIT-212 G1, G2, G3, G5, G6:", class, "the default path, no DA, a DA re-run, Skip after the DA, a loaded session's DA"), {
        skip_if_not(workflowOwnerCapable())
        root <- withr::local_tempdir()
        seen <- k194Notes()
        omic <- k194Omic(class)
        prefix <- kit172Binding(class)$prefix
        # G1: run, Apply, DA and its publication (the writers' own names pinned against the files they wrote), Save, Copy,
        # Report.
        q <- k207Normalised(file.path(root, "g1"), class, da = TRUE)
        dirs <- k212Dirs(file.path(root, "g1-project"), omic)
        k212Publish(q$fx$binding, k194Owner(q$fx), dirs[[omic]])
        expect_gt(length(list.files(dirs[[omic]]$da_output_dir, "_long_annot[.]xlsx$")), 0L)
        k212ExpectDefault("G1", q$fx$workflow, omic, dirs, k194Owner(q$fx))
        # G2: no DA and no DA file.
        q <- k207Normalised(file.path(root, "g2"), class)
        dirs <- k212Dirs(file.path(root, "g2-project"), omic)
        rec <- k207Save(q$fx$workflow, omic, dirs)
        expect_identical(rec$done, 1L, info = paste("G2", class, paste(rec$notes, collapse = " | ")))
        rec <- k212Copy(q$fx$workflow, omic, dirs, argsSaved = TRUE)
        expect_identical(rec$copies, 1L, info = paste("G2", class, paste(rec$notes, collapse = " | ")))
        if (exists("workflowDaResultStatus", mode = "function")) expect_identical(workflowDaResultStatus(k194Owner(q$fx)), "none")
        # G3: Q1, then the DA run again.
        q <- k212Q1(file.path(root, "g3"), class)
        k207Da(q$fx)
        k212ExpectDefault("G3", q$fx$workflow, omic, k212Dirs(file.path(root, "g3-project"), omic), k194Owner(q$fx))
        # G5: the DA at the correlation result, then Skip's `_norm_complete` on the same table.
        q <- k207Normalised(file.path(root, "g5"), class, da = TRUE)
        kit172Quiet(commitWorkflowOmicsStoreStep(k194Owner(q$fx), paste0(prefix, "_norm_complete"), paste0(prefix, "_norm_complete"),
            q$fx$workflow$config_list, "Correlation filter skipped"))
        k207ExpectHead(q$fx, paste0(prefix, "_norm_complete"), paste("G5", class))
        rec <- k207Save(q$fx$workflow, omic, k212Dirs(file.path(root, "g5-project"), omic))
        expect_identical(rec$done, 1L, info = paste("G5", class, paste(rec$notes, collapse = " | ")))
        if (exists("workflowDaResultStatus", mode = "function")) expect_identical(workflowDaResultStatus(k194Owner(q$fx)), "current")
        # G6: a session loaded at `_correlation_filtered`, its DA in process at the head, Save.
        fx <- k194Fixture(file.path(root, "g6"), class)
        k194Defer(fx)
        k194LoadLegacy(fx, paste0(prefix, "_correlation_filtered"), file.path(root, "g6-files"))
        k207Da(fx)
        rec <- k207Save(fx$workflow, omic, k212Dirs(file.path(root, "g6-project"), omic))
        expect_identical(rec$done, 1L, info = paste("G6", class, paste(rec$notes, collapse = " | ")))
        if (exists("workflowDaResultStatus", mode = "function")) expect_identical(workflowDaResultStatus(k194Owner(fx)), "current")
    })

    test_that(paste("KIT-212 G8, G9:", class, "the owner's own export reopened saves; the managed Save after Q1 queues nothing"), {
        skip_if_not(workflowOwnerCapable())
        root <- withr::local_tempdir()
        seen <- k194Notes()
        omic <- k194Omic(class)
        prefix <- kit172Binding(class)$prefix
        # G8: run, Apply, DA, the owner's export (the tab's "Export Normalized Session"), reopened by Load Filtered Session.
        q <- k207Normalised(file.path(root, "g8"), class, da = TRUE)
        export_dir <- file.path(root, "g8-export")
        dir.create(export_dir)
        kit172Quiet(shiny::isolate(exportOmicNormSession(k194Input(), q$fx$workflow, list(source_dir = export_dir), "k212", q$nd,
            q$fx$addLog, q$fx$config)))
        k194Load(q$fx, export_dir)
        rec <- k207Save(q$fx$workflow, omic, k212Dirs(file.path(root, "g8-project"), omic))
        cat(sprintf("KIT-212 G8 %s\thead=%s\tstatus=%s\tdone=%d\n", class, k194Head(q$fx), k212Status(k194Owner(q$fx)), rec$done))
        expect_identical(rec$done, 1L, info = paste("G8", class, paste(rec$notes, collapse = " | ")))
        if (exists("workflowDaResultStatus", mode = "function")) expect_identical(workflowDaResultStatus(k194Owner(q$fx)), "current")
        # G9: managed, Q1's Save queues nothing (head; the base queues the product); the default path queues the product.
        q <- k212Q1(file.path(root, "g9"), class)
        runner <- k194Runner(q$fx$workflow, hold = TRUE)
        rec <- k207Save(q$fx$workflow, omic, k212Dirs(file.path(root, "g9-q1"), omic), runner = runner$runner)
        cat(sprintf("KIT-212 G9 %s Q1\tqueued=%d\tnotes=%s\n", class, length(runner$queued), paste(rec$notes, collapse = " | ")))
        expect_identical(length(runner$queued), 0L, info = paste("G9 Q1", class))
        q <- k207Normalised(file.path(root, "g9-default"), class, da = TRUE)
        runner <- k194Runner(q$fx$workflow, hold = TRUE)
        rec <- k207Save(q$fx$workflow, omic, k212Dirs(file.path(root, "g9-default-project"), omic), runner = runner$runner)
        expect_identical(length(runner$queued), 1L, info = paste("G9 default", class, paste(rec$notes, collapse = " | ")))
    })
}

test_that("KIT-212 G1: transcriptomics, the default path with the DE publication", {
    skip_if_not_installed("edgeR")
    seen <- new.env(parent = emptyenv())
    seen$notes <- character()
    project <- module_ci_trx_project(module_ci_trx_corpus())
    owner <- project$owner
    k194TrxRun(project, omicNormReactiveState(), list(ruv_mode = "manual", ruv_k = 2L, ruv_percentage = 10), seen)
    expect_identical(owner$current_state, "trx_ruvs", info = "G1 trx run precondition")
    kit172Quiet(runTrxReferenceDaAnalysis(owner, getWorkflowStateReference(owner), list(da_q_val_thresh = 0.05, treat_lfc_cutoff = 0),
        trxDaContrastEntries(project$workflow$contrasts_tbl)))
    k212Publish(transcriptomicsDeclarationBinding, owner, project$paths)
    expect_gt(length(list.files(project$paths$da_output_dir, "^de_results_")), 0L)
    k212ExpectDefault("G1", project$workflow, "transcriptomics", project$project_dirs, owner)
})

test_that("KIT-212 G1, G4, G4b, G8: proteomics, the default path; the DA then enrichment (twice); the archive reopened", {
    skip_if_not(workflowOwnerCapable())
    skip_if_not_installed("Glimma")
    # G1: normalise, correlation, export, the DA with its publication, Save, Copy, Report.
    fixture <- k212Prot()
    k212Export(fixture, fixture$norm)
    k212ProtDa(fixture, publish = TRUE)
    dirs <- k212ProtDirs(fixture, "g1")
    seen <- k212ProtSave(fixture, dirs)
    expect_length(seen$failed, 0L)
    expect_length(seen$done, 1L)
    k212ExpectDefault("G1", fixture$workflow, "proteomics", dirs, fixture$manager, save = FALSE)
    # G4, G4b: the DA, an enrichment, Save; a second enrichment (it names itself as its previous state).
    fixture <- k212Prot()
    k212Export(fixture, fixture$norm)
    k212ProtDa(fixture)
    k212Enrichment(fixture)
    seen <- k212ProtSave(fixture, k212ProtDirs(fixture, "g4"))
    cat(sprintf("KIT-212 G4 prot\thead=%s\tstatus=%s\tfailed=%d\tnotes=%s\n", fixture$manager$current_state, k212Status(fixture$manager),
        length(seen$failed), paste(seen$notes, collapse = " | ")))
    expect_length(seen$failed, 0L)
    expect_length(seen$done, 1L)
    if (exists("workflowDaResultStatus", mode = "function")) expect_identical(workflowDaResultStatus(fixture$manager), "current")
    k212Enrichment(fixture)
    if (exists("workflowDaResultStatus", mode = "function")) expect_identical(workflowDaResultStatus(fixture$manager), "current", info = "G4b")
    # G8: the DA, the export, the archive reopened into a new owner, its export, Save.
    fixture <- k212Prot()
    k212ProtDa(fixture)
    k212Export(fixture, fixture$norm)
    archives <- sort(list.files(fixture$paths$source_dir, "^filtered_session_.*[.]msdia[.]tar$", full.names = TRUE))
    expect_true(length(archives) >= 1L, info = "G8 precondition: an archive")
    reopened <- k212Reopen(archives[[length(archives)]])
    status <- k212Status(reopened$manager)
    k212Export(reopened, reopened$norm)
    dirs <- list(proteomics = c(reopened$paths, list(integration_dir = file.path(reopened$root, "integration"), base_dir = reopened$root)))
    seen <- k212ProtSave(reopened, dirs)
    cat(sprintf("KIT-212 G8 prot\thead=%s\tstatus_after_reopen=%s\tfailed=%d\tnotes=%s\n", reopened$manager$current_state, status,
        length(seen$failed), paste(seen$notes, collapse = " | ")))
    expect_length(seen$failed, 0L)
    expect_length(seen$done, 1L)
    if (exists("workflowDaResultStatus", mode = "function")) expect_identical(status, "current", info = "G8 after reopen")
})

test_that("KIT-212 G7, G15: the W1 fixtures' status and unmade files; the in-memory proteomics manager still copies", {
    skip_if_not(workflowOwnerCapable())
    rows <- list()
    lipid <- module_ci_lipid_summary_workflow()
    metab <- module_ci_metab_summary_workflow()
    for (case in list(list("lipid", lipid, "lipidomics", module_ci_lipid_summary_paths, module_ci_lipid_summary_write_publication_inputs),
        list("metab", metab, "metabolomics", module_ci_metab_summary_paths, module_ci_metab_summary_write_publication_inputs))) {
        for (inputs in c("complete", "optional_missing", "none")) {
            dirs <- case[[4L]]()
            withr::defer(unlink(dirs[[case[[3L]]]]$base_dir, recursive = TRUE))
            if (!identical(inputs, "none")) case[[5L]](dirs, include_optional = identical(inputs, "complete"))
            unmade <- k212Unmade(case[[2L]]$state_manager, dirs[[case[[3L]]]], case[[3L]])
            rows[[length(rows) + 1L]] <- c(case[[1L]], inputs, k212Status(case[[2L]]$state_manager), length(unmade))
            if (exists("workflowDaFilesUnmade", mode = "function")) {
                expect_identical(length(unmade) > 0L, !identical(inputs, "none"), info = paste("G7", case[[1L]], inputs))
            }
        }
    }
    fake <- module_ci_prot_summary_workflow(module_ci_prot_summary_branch_config("DIA"))
    dirs <- module_ci_prot_summary_paths()
    withr::defer(unlink(dirs$proteomics$base_dir, recursive = TRUE))
    module_ci_prot_summary_write_publication_inputs(dirs)
    rows[[length(rows) + 1L]] <- c("prot_route", "complete", k212Status(fake$state_manager),
        length(k212Unmade(fake$state_manager, dirs$proteomics, "proteomics")))
    owner <- k212ProtW1Owner()
    rows[[length(rows) + 1L]] <- c("prot_owner", "complete", k212Status(owner$manager),
        length(k212Unmade(owner$manager, owner$project_dirs$proteomics, "proteomics")))
    for (row in rows) cat(sprintf("KIT-212 G7\t%s\n", paste(row, collapse = "\t")))
    # G15: the in-memory manager with the fixture's DA files: Copy proceeds.
    rec <- k212Copy(fake, "proteomics", dirs, argsSaved = TRUE)
    k212Line("G15", "prot", rec, sprintf("\tcopies=%d", rec$copies))
    expect_identical(rec$copies, 1L, info = paste("G15", paste(rec$notes, collapse = " | ")))
})
