# KIT-207 (VAL-105 (2), VAL-148): the lipidomics, metabolomics and transcriptomics Summary refuses to write its final
# object, study parameters, analysis record, integration file, report bundle and session state from a head that is not
# the normalisation's result: lipidomics and metabolomics by the DA's head rule of KIT-194 (without its precondition on
# the normalisation status), transcriptomics by its DA's reader (a normalised state). Written before the fix: on the
# tests-only commit each refused case (E1 to E9, E11) writes as before, so it fails at the values the ticket states
# (section 10); the G rows pass on both commits (G7 and G8 call the new guard only where it exists). Lines
# `KIT-207 <row>` are the evidence the ticket's records copy (outcome.md: the base values, 'G1 digests', G8's premise).

library(testthat)

# ------------------------------------------------------------------------------------------------
# Copied unchanged from tests/testthat/test-omics-norm-state-guards.R at 3d4a585d (KIT-194; KIT-207 E-7): lines 18-124,
# 131-134, 138-202, 219-260 and 303-321. That file is not edited.
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

# ------------------------------------------------------------------------------------------------
# KIT-207's harness.
# ------------------------------------------------------------------------------------------------

# The four texts of section 7, as the ticket drafts them (VAL-41: the owner's wording pass may change them).
k207Texts <- list(
    summary = paste("The current table is not the result of the last normalisation: Quality Control was reverted, or the",
        "normalisation was reset, after it ran. Run the normalisation (and the correlation filter or its Skip) again, or,",
        "if you have only just pressed Reset, press Revert once in a Quality Control sub-tab; then press this button again."),
    summary_trx = paste("The current table is not the result of the last normalisation: Quality Control was reverted, or the",
        "normalisation was reset, after it ran. Run the normalisation again, then press this button again."),
    summary_none = paste("No normalisation has been run in this session, so there is no final table to save. Run the",
        "normalisation (and the correlation filter or its Skip), then press this button again."),
    summary_loaded = paste("This session file was saved after Quality Control was reverted past the normalisation, so its",
        "table (\"%s\") was not normalised. Import the data again, run quality control and the normalisation, then press",
        "this button again."))

k207Refused <- function(text) paste("The session could not be saved:", text)

k207Classes <- c(lipid = "LipidomicsAssayData", metab = "MetaboliteAssayData")

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

#' The LCMS_Pos maximum of an object as `k194Max` reads a state's.
k207ObjectMax <- function(object, binding) {
    table <- as.data.frame(methods::slot(object, binding$data_slot)$LCMS_Pos)
    runs <- intersect(kit172Design()$Run, colnames(table))
    max(as.matrix(table[, runs, drop = FALSE]), na.rm = TRUE)
}

#' What the save handed the RDS writer: per call the object's class (and its LCMS_Pos maximum), or NULL.
k207RdsText <- function(rec, binding) {
    if (!length(rec$rds)) return("<none>")
    paste(vapply(rec$rds, function(call) {
        object <- call$object
        if (is.null(object)) return("NULL")
        if (!isS4(object) || is.null(binding)) return(class(object)[[1L]])
        sprintf("%s max=%s", class(object)[[1L]], format(k207ObjectMax(object, binding), digits = 12))
    }, ""), collapse = "; ")
}

#' One evidence line per refused row (both commits): the head, the state the save retrieved, what it wrote.
k207Line <- function(row, label, head, rec, binding = NULL, extra = "") {
    retrieved <- grep("Retrieved DATA S4 object|final transcriptomics state", rec$cat, value = TRUE)
    cat(sprintf("KIT-207 %s %s\thead=%s\tretrieved=%s\trds=%s\tdone=%d\tfailed=%s\tnotes=%s%s\n", row, label, head,
        if (length(retrieved)) trimws(retrieved[[1L]]) else "<none>", k207RdsText(rec, binding), rec$done,
        paste(vapply(rec$failed, function(e) class(e)[[1L]], ""), collapse = ","), paste(rec$notes, collapse = " | "), extra))
}

#' The refusal of a Save (E1's head expectation): the notification, nothing handed to the RDS writer, no parameters or
#' record in the source directory, `done` not called, `failed` once with the class, no state retrieved.
k207ExpectRefused <- function(rec, dirs, omic, text, info, binding = NULL) {
    source_dir <- dirs[[omic]]$source_dir
    info <- paste(info, "| notes:", paste(rec$notes, collapse = " | "), "| cat:", paste(trimws(rec$cat), collapse = " | "),
        "| rds:", k207RdsText(rec, binding))
    expect_true(paste0("[error] ", k207Refused(text)) %in% rec$notes, info = info)
    expect_identical(length(rec$rds), 0L, info = info)
    expect_false(file.exists(file.path(source_dir, "study_parameters.txt")), info = info)
    expect_false(file.exists(file.path(source_dir, "analysis_record.json")), info = info)
    expect_identical(rec$done, 0L, info = info)
    expect_true(length(rec$failed) == 1L && inherits(rec$failed[[1L]], "omic_norm_correlation_input_stale"), info = info)
    expect_false(any(grepl("Retrieved DATA S4 object", rec$cat, fixed = TRUE)), info = info)
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

#' The study parameters' SHA-256 without the lines that differ by run or by build (the timestamps, the version block's
#' commit, branch and source) and with the temporary root masked.
k207ParamsSha <- function(file, root) {
    lines <- readLines(file, warn = FALSE)
    lines <- lines[!grepl("Timestamp:|^(Repository|Branch|Commit|Package Version|Source):", lines)]
    for (prefix in unique(c(normalizePath(root, mustWork = FALSE), root))) lines <- gsub(prefix, "<root>", lines, fixed = TRUE)
    digest::digest(paste(lines, collapse = "\n"), algo = "sha256", serialize = FALSE)
}

# ------------------------------------------------------------------------------------------------
# E1 to E5, E7, E8: lipidomics and metabolomics, refused (written on the base).
# ------------------------------------------------------------------------------------------------

for (class in c("lipid", "metab")) {
    test_that(paste("KIT-207 E1, E2, E3:", class, "a Save after a Revert to QC complete, into the normalisation, or earlier is refused"), {
        skip_if_not(workflowOwnerCapable())
        root <- withr::local_tempdir()
        # The run, Apply and Revert notify and show progress (KIT-194's doubles).
        seen <- k194Notes()
        omic <- k194Omic(class)
        binding <- kit172Binding(class)
        prefix <- binding$prefix
        # E1 (Q1): run, Apply, DA, Revert to `<prefix>_qc_complete`.
        q <- k207Normalised(file.path(root, "e1"), class, da = TRUE)
        expect_identical(k194Revert(q$fx, paste0(prefix, "_qc_complete")), paste0(prefix, "_qc_complete"), info = "E1 precondition")
        qc_max <- k194Max(q$fx, k194Ref(q$fx, paste0(prefix, "_qc_complete")))
        dirs <- k207Dirs(file.path(root, "e1-project"), omic)
        rec <- k207Save(q$fx$workflow, omic, dirs)
        k207Line("E1", class, k194Head(q$fx), rec, binding, sprintf("\tqc_complete_max=%s", format(qc_max, digits = 12)))
        k207ExpectRefused(rec, dirs, omic, k207Texts$summary, paste("E1", class), binding)
        # E2 (Q2): run, Apply, Revert into the normalisation (head `<prefix>_log2`).
        q <- k207Normalised(file.path(root, "e2"), class)
        expect_identical(k194Revert(q$fx, paste0(prefix, "_log2")), paste0(prefix, "_log2"), info = "E2 precondition")
        dirs <- k207Dirs(file.path(root, "e2-project"), omic)
        rec <- k207Save(q$fx$workflow, omic, dirs)
        k207Line("E2", class, k194Head(q$fx), rec, binding)
        k207ExpectRefused(rec, dirs, omic, k207Texts$summary, paste("E2", class), binding)
        # E3 (Q3): run, Apply, Revert to `<prefix>_intensity_filtered`.
        q <- k207Normalised(file.path(root, "e3"), class)
        expect_identical(k194Revert(q$fx, paste0(prefix, "_intensity_filtered")), paste0(prefix, "_intensity_filtered"),
            info = "E3 precondition")
        dirs <- k207Dirs(file.path(root, "e3-project"), omic)
        rec <- k207Save(q$fx$workflow, omic, dirs)
        k207Line("E3", class, k194Head(q$fx), rec, binding)
        k207ExpectRefused(rec, dirs, omic, k207Texts$summary, paste("E3", class), binding)
    })

    test_that(paste("KIT-207 E4, E5:", class, "a Save after Reset, and on a session loaded under a QC state's name, is refused"), {
        skip_if_not(workflowOwnerCapable())
        root <- withr::local_tempdir()
        # The run, Apply and Revert notify and show progress (KIT-194's doubles).
        seen <- k194Notes()
        omic <- k194Omic(class)
        binding <- kit172Binding(class)
        prefix <- binding$prefix
        # E4 (Q4): run, Apply, DA, Reset.
        q <- k207Normalised(file.path(root, "e4"), class, da = TRUE)
        k194Reset(q$fx, q$nd)
        k207ExpectHead(q$fx, paste0(prefix, "_reset"), paste("E4", class))
        dirs <- k207Dirs(file.path(root, "e4-project"), omic)
        rec <- k207Save(q$fx$workflow, omic, dirs)
        k207Line("E4", class, k194Head(q$fx), rec, binding)
        k207ExpectRefused(rec, dirs, omic, k207Texts$summary, paste("E4", class), binding)
        # E5 (Q5): a session file named `<prefix>_qc_complete` loaded by "Load Filtered Session".
        fx <- k194Fixture(file.path(root, "e5"), class)
        k194Defer(fx)
        k194LoadLegacy(fx, paste0(prefix, "_qc_complete"), root)
        k207ExpectHead(fx, paste0(prefix, "_qc_complete"), paste("E5", class))
        dirs <- k207Dirs(file.path(root, "e5-project"), omic)
        rec <- k207Save(fx$workflow, omic, dirs)
        k207Line("E5", class, k194Head(fx), rec, binding)
        k207ExpectRefused(rec, dirs, omic, sprintf(k207Texts$summary_loaded, paste0(prefix, "_qc_complete")), paste("E5", class),
            binding)
    })

    test_that(paste("KIT-207 E7:", class, "a managed Save after a Revert to QC complete queues nothing"), {
        skip_if_not(workflowOwnerCapable())
        root <- withr::local_tempdir()
        # The run, Apply and Revert notify and show progress (KIT-194's doubles).
        seen <- k194Notes()
        omic <- k194Omic(class)
        binding <- kit172Binding(class)
        prefix <- binding$prefix
        q <- k207Normalised(file.path(root, "e7"), class, da = TRUE)
        expect_identical(k194Revert(q$fx, paste0(prefix, "_qc_complete")), paste0(prefix, "_qc_complete"), info = "E7 precondition")
        runner <- k194Runner(q$fx$workflow, hold = TRUE)
        dirs <- k207Dirs(file.path(root, "e7-project"), omic)
        rec <- k207Save(q$fx$workflow, omic, dirs, runner = runner$runner)
        queued <- if (length(runner$queued)) runner$queued[[1L]]$actions[[1L]]$parameters$state_name %||% "<unnamed>" else "<none>"
        k207Line("E7", class, k194Head(q$fx), rec, NULL, sprintf("\tqueued_state=%s", queued))
        info <- paste("E7", class, "queued state", queued)
        expect_identical(length(runner$queued), 0L, info = info)
        expect_identical(length(runner$operations), 0L, info = info)
        k207ExpectRefused(rec, dirs, omic, k207Texts$summary, paste("E7", class))
    })

    test_that(paste("KIT-207 E8:", class, "Generate Report and Export Session State refuse after a Revert to QC complete"), {
        skip_if_not(workflowOwnerCapable())
        root <- withr::local_tempdir()
        # The run, Apply and Revert notify and show progress (KIT-194's doubles).
        seen <- k194Notes()
        omic <- k194Omic(class)
        prefix <- kit172Binding(class)$prefix
        config <- omicGuiConfig(omic)$summary
        q <- k207Normalised(file.path(root, "e8"), class, da = TRUE)
        expect_identical(k194Revert(q$fx, paste0(prefix, "_qc_complete")), paste0(prefix, "_qc_complete"), info = "E8 precondition")
        dirs <- k207Dirs(file.path(root, "e8-project"), omic)
        # The project's scripts/<omic>/ holds the template, so the report goes straight to the render.
        writeLines("KIT-207 template", file.path(dirs[[omic]]$base_dir, "scripts", omic, config$report$template))
        calls <- new.env(parent = emptyenv())
        calls$render <- 0L
        calls$bundle <- 0L
        calls$notes <- character()
        notify <- function(message, type, duration) calls$notes <- c(calls$notes, paste0("[", type, "] ", message))
        kit172Quiet(utils::capture.output(renderWorkflowSummaryReport(config, q$fx$workflow, dirs, omic, "k207", notify = notify,
            progress = function(...) NULL, renderFn = function(omic_type, experiment_label, rmd_filename) {
                calls$render <- calls$render + 1L
                path <- file.path(root, "e8-report.html")
                writeLines("KIT-207 report", path)
                path
            }, exportRecordFn = function(...) {
                calls$bundle <- calls$bundle + 1L
                invisible(NULL)
            }, catFn = function(...) NULL)))
        cat(sprintf("KIT-207 E8 %s report\trender=%d\tbundle=%d\tnotes=%s\n", class, calls$render, calls$bundle,
            paste(calls$notes, collapse = " | ")))
        info <- paste("E8 report", class, paste(calls$notes, collapse = " | "))
        expect_true(paste("[error] Report generation failed:", k207Texts$summary) %in% calls$notes, info = info)
        expect_identical(calls$render, 0L, info = info)
        expect_identical(calls$bundle, 0L, info = info)
        # Export Session State.
        calls$notes <- character()
        kit172Quiet(utils::capture.output(writeWorkflowSummarySessionState(config, q$fx$workflow, dirs, omic, "k207", "",
            values = list(workflow_args_saved = FALSE, files_copied = FALSE, report_generated = FALSE, report_path = NULL),
            notify = notify)))
        files <- list.files(dirs[[omic]]$source_dir, "^session_state_.*[.]RDS$", full.names = TRUE)
        counts <- if (length(files)) paste(format(unlist(readRDS(files[[1L]])$feature_counts)), collapse = ",") else "<none>"
        cat(sprintf("KIT-207 E8 %s session_state\tfiles=%d\tfeature_counts=%s\tnotes=%s\n", class, length(files), counts,
            paste(calls$notes, collapse = " | ")))
        info <- paste("E8 session state", class, "files", length(files), "feature_counts", counts, paste(calls$notes, collapse = " | "))
        expect_true(paste("[error] Export failed:", k207Texts$summary) %in% calls$notes, info = info)
        expect_identical(length(files), 0L, info = info)
    })
}

# ------------------------------------------------------------------------------------------------
# E6, E9: transcriptomics, refused (written on the base).
# ------------------------------------------------------------------------------------------------

test_that("KIT-207 E6: transcriptomics, a Save after a Revert to trx_qc_complete is refused and writes no integration file", {
    skip_if_not_installed("edgeR")
    seen <- new.env(parent = emptyenv())
    seen$notes <- character()
    project <- module_ci_trx_project(module_ci_trx_corpus())
    owner <- project$owner
    nd <- omicNormReactiveState()
    k194TrxRun(project, nd, list(ruv_mode = "manual", ruv_k = 2L, ruv_percentage = 10), seen)
    expect_identical(owner$current_state, "trx_ruvs", info = paste("E6 run precondition:", paste(seen$notes, collapse = " | ")))
    kit172Quiet(runTrxReferenceRuviiiExport(owner, "group"))
    expect_false(is.null(owner$.context$snapshot$exports$trx_ruviii_export), info = "E6 export precondition")
    selectWorkflowStateReference(owner, "trx_qc_complete")
    expect_identical(owner$current_state, "trx_qc_complete", info = "E6 head precondition")
    dirs <- project$project_dirs
    tsv <- file.path(k207IntegrationDir(dirs$transcriptomics), "transcriptomics_normalized_data.tsv")
    rec <- k207Save(project$workflow, "transcriptomics", dirs)
    k207Line("E6", "trx", owner$current_state, rec, NULL, sprintf("\ttsv=%s", file.exists(tsv)))
    k207ExpectRefused(rec, dirs, "transcriptomics", k207Texts$summary_trx, "E6")
    expect_false(file.exists(tsv), info = "E6 integration file")
})

test_that("KIT-207 E9: transcriptomics, a second run without a new RUVIII export refuses the earlier export and says so", {
    skip_if_not_installed("edgeR")
    seen <- new.env(parent = emptyenv())
    seen$notes <- character()
    project <- module_ci_trx_project(module_ci_trx_corpus())
    owner <- project$owner
    nd <- omicNormReactiveState()
    k194TrxRun(project, nd, list(ruv_mode = "manual", ruv_k = 2L, ruv_percentage = 10), seen)
    expect_identical(owner$current_state, "trx_ruvs", info = paste("E9 run precondition:", paste(seen$notes, collapse = " | ")))
    kit172Quiet(runTrxReferenceRuviiiExport(owner, "group"))
    # A second "Run Normalization Pipeline" (D1 (b): a new trx_ruvs from the remembered start), k 1; no new export.
    k194TrxRun(project, nd, list(ruv_mode = "manual", ruv_k = 1L, ruv_percentage = 10), seen)
    expect_identical(owner$current_state, "trx_ruvs", info = paste("E9 second run precondition:", paste(seen$notes, collapse = " | ")))
    export <- owner$.context$snapshot$exports$trx_ruviii_export
    expect_false(identical(export$source$entity_id, getWorkflowStateReference(owner, "trx_ruvs")$entity_id),
        info = "E9 precondition: the accepted export is the first run's")
    dirs <- project$project_dirs
    tsv <- file.path(k207IntegrationDir(dirs$transcriptomics), "transcriptomics_normalized_data.tsv")
    rec <- k207Save(project$workflow, "transcriptomics", dirs)
    k207Line("E9", "trx", owner$current_state, rec, NULL, sprintf("\ttsv=%s", file.exists(tsv)))
    info <- paste("E9", paste(rec$notes, collapse = " | "), "| tsv", file.exists(tsv))
    expect_true(paste("[error]", k207Refused(paste("the accepted RUVIII export was made from an earlier RUV step: run the export",
        "again, then save again"))) %in% rec$notes, info = info)
    expect_true(length(rec$failed) == 1L && inherits(rec$failed[[1L]], "integration_input_invalid"), info = info)
    expect_false(file.exists(tsv), info = info)
    expect_true(file.exists(file.path(dirs$transcriptomics$source_dir, "study_parameters.txt")), info = info)
})

# ------------------------------------------------------------------------------------------------
# E11: the tab mounted (lipidomics).
# ------------------------------------------------------------------------------------------------

test_that("KIT-207 E11: lipid, the Summary tab mounted: a refused Save, then Copy stops at the missing study parameters", {
    skip_if_not(workflowOwnerCapable())
    root <- withr::local_tempdir()
    seen <- k194Notes()
    q <- k207Normalised(file.path(root, "e11"), "lipid", da = TRUE)
    fx <- q$fx
    expect_identical(k194Revert(fx, "lipid_qc_complete"), "lipid_qc_complete", info = "E11 precondition")
    fx$workflow$design_matrix <- kit172Design()
    withr::local_envvar(MULTISCHOLAR_TEST_MODE = "true")
    invisible(utils::capture.output(project_dirs <- suppressMessages(setupDirectories(base_dir = file.path(root, "project", "kit207"),
        omic_types = "lipidomics", label = NULL, force = TRUE, reuse_existing = TRUE))))
    had <- exists("project_dirs", envir = .GlobalEnv, inherits = FALSE)
    old <- if (had) get("project_dirs", envir = .GlobalEnv, inherits = FALSE)
    withr::defer(if (had) assign("project_dirs", old, envir = .GlobalEnv) else
        if (exists("project_dirs", envir = .GlobalEnv, inherits = FALSE)) rm("project_dirs", envir = .GlobalEnv))
    assign("project_dirs", project_dirs, envir = .GlobalEnv)
    testthat::local_mocked_bindings(workflowOmicsManagedActions = function(omic_type) NULL, .package = "MultiScholaR")
    kit172Quiet(utils::capture.output(shiny::testServer(mod_lipid_summary_server, args = list(project_dirs = project_dirs,
        omic_type = "lipidomics", experiment_label = NULL, workflow_data = fx$workflow), {
        session$setInputs(experiment_label = "k207", description = "d")
        session$setInputs(save_workflow_args = 1L)
        session$setInputs(copy_to_publication = 1L)
    })))
    params <- file.path(project_dirs$lipidomics$source_dir, "study_parameters.txt")
    cat(sprintf("KIT-207 E11 lipid\tstudy_parameters=%s\tnotes=%s\n", file.exists(params), paste(seen$notes, collapse = " | ")))
    info <- paste("E11", paste(seen$notes, collapse = " | "))
    expect_true(k194Noted(seen, paste("[error]", k207Refused(k207Texts$summary))), info = info)
    expect_true(k194Noted(seen, "Copy error: Required publication artifacts missing: Study Parameters"), info = info)
    expect_false(file.exists(params), info = info)
})

# ------------------------------------------------------------------------------------------------
# G1 to G9: allowed, and unchanged (both commits; G7 and G8 call the new guard only where it exists).
# ------------------------------------------------------------------------------------------------

for (class in c("lipid", "metab")) {
    test_that(paste("KIT-207 G1, G2, G9:", class, "the default path, Skip, and Reset then one Revert save as before"), {
        skip_if_not(workflowOwnerCapable())
        root <- withr::local_tempdir()
        # The run, Apply and Revert notify and show progress (KIT-194's doubles).
        seen <- k194Notes()
        omic <- k194Omic(class)
        binding <- kit172Binding(class)
        prefix <- binding$prefix
        # G1: run, Apply, DA, Save.
        q <- k207Normalised(file.path(root, "g1"), class, da = TRUE)
        k207ExpectHead(q$fx, paste0(prefix, "_correlation_filtered"), paste("G1", class))
        dirs <- k207Dirs(file.path(root, "g1-project"), omic)
        rec <- k207Save(q$fx$workflow, omic, dirs)
        info <- paste("G1", class, paste(rec$notes, collapse = " | "), "| rds:", k207RdsText(rec, binding))
        expect_identical(rec$done, 1L, info = info)
        expect_identical(length(rec$rds), 1L, info = info)
        object <- if (length(rec$rds)) rec$rds[[1L]]$object
        expect_true(methods::is(object, k207Classes[[class]]), info = info)
        if (methods::is(object, k207Classes[[class]])) {
            expect_identical(k207ObjectMax(object, binding), k194Max(q$fx, k194Ref(q$fx, paste0(prefix, "_correlation_filtered"))),
                info = info)
        }
        params <- file.path(dirs[[omic]]$source_dir, "study_parameters.txt")
        record <- file.path(dirs[[omic]]$source_dir, "analysis_record.json")
        expect_true(file.exists(params), info = info)
        expect_true(file.exists(record), info = info)
        if (file.exists(params) && file.exists(record) && !is.null(object)) {
            steps <- vapply(jsonlite::read_json(record, simplifyVector = FALSE)$steps, function(step) as.character(step$state), "")
            cat(sprintf("KIT-207 G1 %s object=%s params=%s record_steps=%s\n", omic, digest::digest(object, algo = "sha256"),
                k207ParamsSha(params, root), paste(steps, collapse = ",")))
            cat(sprintf("KIT-207 G1-assays %s assays=%s\n", omic,
                digest::digest(methods::slot(object, binding$data_slot), algo = "sha256")))
        }
        # G2: Skip at the result, Save (`<prefix>_norm_complete`).
        q <- k207Normalised(file.path(root, "g2"), class, apply = FALSE)
        expect_true(isTRUE(k194Skip(q$fx, q$nd)), info = paste("G2 Skip precondition", class))
        k207ExpectHead(q$fx, paste0(prefix, "_norm_complete"), paste("G2", class))
        rec <- k207Save(q$fx$workflow, omic, k207Dirs(file.path(root, "g2-project"), omic))
        info <- paste("G2", class, paste(rec$notes, collapse = " | "), paste(rec$cat, collapse = " | "))
        expect_identical(rec$done, 1L, info = info)
        expect_true(any(grepl(sprintf("Retrieved DATA S4 object from state '%s_norm_complete'", prefix), rec$cat, fixed = TRUE)),
            info = info)
        # G9: run, Apply, DA, Reset, one Revert in a QC sub-tab (the run's `<prefix>_correlation_filtered`), Save.
        q <- k207Normalised(file.path(root, "g9"), class, da = TRUE)
        k194Reset(q$fx, q$nd)
        k207ExpectHead(q$fx, paste0(prefix, "_reset"), paste("G9 Reset", class))
        k194Revert(q$fx, times = 1L)
        k207ExpectHead(q$fx, paste0(prefix, "_correlation_filtered"), paste("G9 Revert", class))
        rec <- k207Save(q$fx$workflow, omic, k207Dirs(file.path(root, "g9-project"), omic))
        expect_identical(rec$done, 1L, info = paste("G9", class, paste(rec$notes, collapse = " | ")))
    })

    test_that(paste("KIT-207 G3, G4, G6:", class, "one Revert inside the run, loaded normalised files and the managed path save as before"), {
        skip_if_not(workflowOwnerCapable())
        root <- withr::local_tempdir()
        # The run, Apply and Revert notify and show progress (KIT-194's doubles).
        seen <- k194Notes()
        omic <- k194Omic(class)
        prefix <- kit172Binding(class)$prefix
        # G3 (lipid; KIT-194 runs RUV on F172 for lipidomics only): run with RUV, Apply, one Revert to `_ruv_corrected`.
        if (identical(class, "lipid")) {
            q <- k207Normalised(file.path(root, "g3"), class, input = k194Input(ruv = kit172LegacyRuv("manual")))
            k207ExpectHead(q$fx, paste0(prefix, "_correlation_filtered"), "G3 Apply")
            k194Revert(q$fx, times = 1L)
            k207ExpectHead(q$fx, paste0(prefix, "_ruv_corrected"), "G3 Revert")
            rec <- k207Save(q$fx$workflow, omic, k207Dirs(file.path(root, "g3-project"), omic))
            info <- paste("G3", paste(rec$notes, collapse = " | "), paste(rec$cat, collapse = " | "))
            expect_identical(rec$done, 1L, info = info)
            expect_true(any(grepl(sprintf("Retrieved DATA S4 object from state '%s_ruv_corrected'", prefix), rec$cat, fixed = TRUE)),
                info = info)
        }
        # G4: session files loaded under `<prefix>_correlation_filtered` and `loaded_for_de`.
        for (state in c(paste0(prefix, "_correlation_filtered"), "loaded_for_de")) {
            fx <- k194Fixture(file.path(root, paste0("g4-", state)), class)
            k194Defer(fx)
            k194LoadLegacy(fx, state, file.path(root, paste0("g4-files-", state)))
            k207ExpectHead(fx, state, paste("G4", class))
            rec <- k207Save(fx$workflow, omic, k207Dirs(file.path(root, paste0("g4-project-", state)), omic))
            expect_identical(rec$done, 1L, info = paste("G4", class, state, paste(rec$notes, collapse = " | ")))
        }
        # G6: the managed default path queues the product with the correlation state.
        q <- k207Normalised(file.path(root, "g6"), class)
        runner <- k194Runner(q$fx$workflow, hold = TRUE)
        rec <- k207Save(q$fx$workflow, omic, k207Dirs(file.path(root, "g6-project"), omic), runner = runner$runner)
        info <- paste("G6", class, paste(rec$notes, collapse = " | "))
        expect_identical(length(runner$queued), 1L, info = info)
        if (length(runner$queued)) {
            expect_identical(runner$queued[[1L]]$actions[[1L]]$parameters$state_name, paste0(prefix, "_correlation_filtered"), info = info)
        }
        expect_identical(length(rec$failed), 0L, info = info)
    })
}

test_that("KIT-207 G5: transcriptomics after its run and DA, and after a run without RUV, saves as before", {
    skip_if_not_installed("edgeR")
    seen <- new.env(parent = emptyenv())
    seen$notes <- character()
    project <- module_ci_trx_project(module_ci_trx_corpus())
    owner <- project$owner
    k194TrxRun(project, omicNormReactiveState(), list(ruv_mode = "manual", ruv_k = 2L, ruv_percentage = 10), seen)
    expect_identical(owner$current_state, "trx_ruvs", info = paste("G5 run precondition:", paste(seen$notes, collapse = " | ")))
    kit172Quiet(runTrxReferenceDaAnalysis(owner, getWorkflowStateReference(owner), list(da_q_val_thresh = 0.05, treat_lfc_cutoff = 0),
        trxDaContrastEntries(project$workflow$contrasts_tbl)))
    rec <- k207Save(project$workflow, "transcriptomics", project$project_dirs)
    expect_identical(rec$done, 1L, info = paste("G5 after DA", paste(rec$notes, collapse = " | ")))
    expect_identical(length(rec$failed), 0L, info = paste("G5 after DA", paste(rec$notes, collapse = " | ")))
    # RUV design none (no factor): the TMM logCPM, no export needed.
    none <- module_ci_trx_project(module_ci_trx_corpus())
    kit172Quiet(runTrxReferenceRuvs(none$owner, NULL, list(), design_parameters = trxRuvDesignParameters("none")))
    expect_identical(none$owner$current_state, "trx_ruvs", info = "G5 none precondition")
    rec <- k207Save(none$workflow, "transcriptomics", none$project_dirs)
    expect_identical(rec$done, 1L, info = paste("G5 none", paste(rec$notes, collapse = " | ")))
    expect_identical(length(rec$failed), 0L, info = paste("G5 none", paste(rec$notes, collapse = " | ")))
})

test_that("KIT-207 G7: the guard passes proteomics and a manager that is not a session owner (head only)", {
    skip_if_not(workflowOwnerCapable())
    if (exists("assertWorkflowSummaryHeadFinal", mode = "function")) {
        root <- withr::local_tempdir()
        fx <- k194Fixture(file.path(root, "g7"), "lipid")
        k194Defer(fx)
        expect_true(assertWorkflowSummaryHeadFinal(fx$workflow, "proteomics"))
        expect_error(assertWorkflowSummaryHeadFinal(fx$workflow, "lipidomics"), class = "omic_norm_correlation_input_stale")
        expect_true(assertWorkflowSummaryHeadFinal(list(state_manager = new.env()), "lipidomics"))
        expect_true(assertWorkflowSummaryHeadFinal(list(state_manager = NULL), "transcriptomics"))
    }
})

test_that("KIT-207 G8: the module-CI Summary owners pass the guard, with the default configuration and with RUV skipped (head only)", {
    skip_if_not(workflowOwnerCapable())
    if (exists("assertWorkflowSummaryHeadFinal", mode = "function")) {
        cases <- list(lipidomics = list(default = function() module_ci_lipid_summary_workflow(),
            skip = function() module_ci_lipid_summary_workflow(module_ci_lipid_summary_config(ruv_mode = "skip"))),
            metabolomics = list(default = function() module_ci_metab_summary_workflow(),
            skip = function() module_ci_metab_summary_workflow(module_ci_metab_summary_config(ruv_mode = "skip"))))
        for (omic in names(cases)) {
            for (variant in names(cases[[omic]])) {
                workflow <- cases[[omic]][[variant]]()
                owner <- workflow$state_manager
                prefix <- omicGuiConfig(omic)$binding$binding$prefix
                same <- tryCatch(withWorkflowStateScope(owner, function() identical(
                    readWorkflowStateMetadata(owner, getWorkflowStateReference(owner))$data$components,
                    readWorkflowStateMetadata(owner, getWorkflowStateReference(owner, paste0(prefix, "_normalized")))$data$components)),
                    error = function(e) paste("unreadable:", conditionMessage(e)))
                cat(sprintf("KIT-207 G8 %s %s\thead=%s\tcomponents_identical_to_normalized=%s\n", omic, variant, owner$current_state, same))
                expect_true(assertWorkflowSummaryHeadFinal(workflow, omic), info = paste(omic, variant, owner$current_state))
            }
        }
    }
})
