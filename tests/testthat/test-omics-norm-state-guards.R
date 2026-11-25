# KIT-194 (VAL-99 (1), VAL-105): the normalisation state guards of lipidomics, metabolomics and transcriptomics.
# Quality-control steps refuse a normalised table; "Run Normalization Pipeline" refuses a normalised input and, pressed
# again at the tab's own result, runs again from the table the previous run started from (D1 (b)); for lipidomics and
# metabolomics "Apply Correlation Filter", its Skip and "Export Normalized Data" refuse unless the head is the last run's
# result (or a correlation-step state made from it), and the DA's head fallback takes only the normalisation's states
# (D2), with a loaded session named after a QC, reset, ITSD or log2 state refused (D4 (b1)). Written before the fix: on
# the tests-only commit each refused case acts, so T1 to T8 fail at the values the ticket states (section 10); T9, T10
# and the "both commits" cases pass on both. Every case drives the application's entry points; only F1 and the classes
# block reach the functions this ticket adds, through `get0`. Lines `KIT-194 effect`, `KIT-194 F1`, `KIT-194 F2` and
# `KIT-194 F3` are the evidence the ticket's records copy (effect.tsv, feasibility.md).

library(testthat)

# ------------------------------------------------------------------------------------------------
# The harness.
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

k194DeferOwner <- function(owner, envir = parent.frame()) {
    force(owner)
    withr::defer(try(closeDIAStateAdapter(owner), silent = TRUE), envir = envir)
}

k194Owner <- function(fixture) fixture$workflow$state_manager
k194History <- function(fixture) as.character(unlist(k194Owner(fixture)$getHistory()))
k194Head <- function(fixture) k194Owner(fixture)$current_state
k194Ref <- function(fixture, state = NULL) getWorkflowStateReference(k194Owner(fixture), state)
k194Snapshot <- function(fixture) list(history = k194History(fixture), head = k194Head(fixture),
    entity = k194Ref(fixture)$entity_id)

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

#' One QC press under `testServer` of its sub-tab's server: "Apply Filter" (`apply`), "Resolve Duplicates" (`resolve`),
#' "Finalize QC" (`finalize`), in process (`runner` NULL) or managed.
k194QcPress <- function(fixture, button, runner = NULL) {
    testthat::local_mocked_bindings(workflowOmicsManagedActions = function(omic_type) runner, .package = "MultiScholaR")
    server <- switch(button, apply = mod_omic_qc_intensity_server, resolve = mod_omic_qc_duplicates_server,
        finalize = mod_omic_qc_finalize_server)
    inputs <- switch(button, apply = list(intensity_cutoff_percentile = 1, proportion_below_cutoff = 0.5, apply_filter = 1L),
        resolve = list(resolve_duplicates = 1L), finalize = list(finalize_qc = 1L))
    kit172Quiet(utils::capture.output(shiny::testServer(server, args = list(workflow_data = fixture$workflow,
        omic_type = fixture$omic, experiment_label = "kit194", config = fixture$config), {
        do.call(session$setInputs, inputs)
    })))
    invisible(NULL)
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

#' "Run DA" with no state loaded (or the loaded `tab`): the frame's `runOmicDaAnalysis` with the run recorded.
k194Da <- function(fixture, tab = new.env(parent = emptyenv()), workflow = fixture$workflow, binding = fixture$binding,
    da = fixture$config$da) {
    calls <- new.env(parent = emptyenv())
    calls$runs <- 0L
    da$preflight <- NULL
    testthat::local_mocked_bindings(executeOmicDaRun = function(tab, workflowData, ...) {
        calls$runs <- calls$runs + 1L
        calls$state <- tab$current_state_name %||% workflowData$state_manager$current_state
        list(status = "success")
    }, workflowOmicsManagedActions = function(omic_type) NULL, .package = "MultiScholaR")
    kit172Quiet(shiny::isolate(runOmicDaAnalysis(list(da_q_val_thresh = 0.05, treat_lfc_cutoff = 0), NULL, workflow, tab, list(),
        da, binding)))
    as.list(calls)
}

k194SetNormStatus <- function(fixture, value) {
    status <- fixture$workflow$tab_status
    status$normalization <- value
    fixture$workflow$tab_status <- status
}

# The 7.6 texts, as the ticket proposes them (VAL-41).
k194QcLead <- "quality-control steps cannot be applied to a table that has already been normalised."
k194RunLead <- "this table has already been normalised, so it cannot be normalised again."
k194RunUnknownLead <- "this tab does not know which table the earlier normalisation started from, so it cannot run it again."
k194NoStep <- "This table was loaded from a saved normalised session, or no quality-control result before normalisation is left in this session"
k194Busy <- "A step is still running. Press \"Run Normalization Pipeline\" again when it has finished."
k194Stale <- "the current table is not the one the last normalisation produced"
k194DaText <- "The current table is not the result of the last normalisation"
k194Step <- function(prefix) sprintf("\"Finalize QC\" (shown as %s_qc_complete)", prefix)

omicsSessionInputK194 <- function() {
    list(norm_method = "quantile", ruv_mode = "skip", apply_itsd = FALSE, itsd_aggregation = NA, log_offset = 0,
        min_pearson_correlation_threshold = 0.5, ruv_grouping_variable = "group")
}

# ------------------------------------------------------------------------------------------------
# Transcriptomics' harness.
# ------------------------------------------------------------------------------------------------

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

k194TrxLog <- function(nd) shiny::isolate(nd$normalization_log)

k194TrxW <- function(owner, state = "trx_ruvs") {
    design <- readTrxWorkflowState(owner, getWorkflowStateReference(owner, state))@design_matrix
    design[grepl("^W_", colnames(design))]
}

k194TrxPress <- function(project, button, config = omicGuiConfig("transcriptomics"), runner = NULL) {
    testthat::local_mocked_bindings(workflowOmicsManagedActions = function(omic_type) runner, .package = "MultiScholaR")
    kit172Quiet(utils::capture.output(if (identical(button, "apply")) {
        shiny::testServer(mod_trx_qc_filter_server, args = list(workflow_data = project$workflow, omic_type = "transcriptomics",
            experiment_label = "kit194"), session$setInputs(apply_filter = 1L))
    } else {
        shiny::testServer(mod_omic_qc_finalize_server, args = list(workflow_data = project$workflow, omic_type = "transcriptomics",
            experiment_label = "kit194", config = config), session$setInputs(finalize_qc = 1L))
    }))
    invisible(NULL)
}

# ------------------------------------------------------------------------------------------------
# F1: the lineage predicate, the adoption marker and the derivation (written here; base has none).
# ------------------------------------------------------------------------------------------------

k194NormStates <- function(binding) {
    if (identical(binding$omic_type, "transcriptomics")) return("trx_ruvs")
    paste0(binding$prefix, c("_itsd_norm", "_log2", "_normalized", "_ruv_corrected", "_correlation_filtered", "_norm_complete"))
}

#' The predicate of 7.1, with the clause that fired and the labels it read.
k194LineageDetail <- function(binding, manager, reference) {
    trx <- identical(binding$omic_type, "transcriptomics")
    labels <- if (trx) c("trx_ruvs", "owner:trx-ruvs") else c(k194NormStates(binding),
        paste0("owner:", binding$prefix, c("-itsd-normalised", "-log2", "-normalised", "-ruv-corrected")))
    products <- if (trx) "trx_ruvs_science" else
        paste0(binding$prefix, c("_itsd_normalisation", "_log2", "_between_sample_normalisation", "_ruv"), "_science")
    spec <- tryCatch(readWorkflowStateMetadata(manager, reference)$data, error = function(e) NULL)
    if (is.null(spec)) return(list(normalised = FALSE, clause = "unreadable", labels = ""))
    seen <- vapply(spec$dependencies, function(d) if (is.list(d) && is.character(d$identity$parameters)) d$identity$parameters[[1L]] else NA_character_, "")
    first <- any(names(spec$products) %in% products)
    second <- any(seen %in% labels)
    list(normalised = first || second, clause = paste(c(if (first) "products", if (second) "dependencies"), collapse = "+"),
        labels = paste(sort(unique(stats::na.omit(seen))), collapse = ","))
}

k194Derives <- function(binding, manager, reference) {
    normalised <- getWorkflowStateReference(manager, paste0(binding$prefix, "_normalized"))
    if (is.null(normalised) || !isSavedWorkflowReference(reference)) return(FALSE)
    if (identical(reference$entity_id, normalised$entity_id)) return(TRUE)
    spec <- readWorkflowStateMetadata(manager, reference)$data
    any(vapply(spec$dependencies, function(d) is.list(d) && identical(d$entity_id, normalised$entity_id), NA)) ||
        identical(spec$components, readWorkflowStateMetadata(manager, normalised)$data$components)
}

#' The package's function when the fix has it, else the copy above.
k194Fn <- function(name, copy) get0(name, envir = asNamespace("MultiScholaR"), inherits = FALSE) %||% copy

k194F1 <- function(binding, manager, reference, expected, omic, path, state) {
    detail <- k194LineageDetail(binding, manager, reference)
    verdict <- k194Fn("omicNormalisedLineage", function(binding, manager, reference) k194LineageDetail(binding, manager, reference)$normalised)(
        binding, manager, reference)
    cat(sprintf("KIT-194 F1\t%s\t%s\t%s\texpected=%s\tverdict=%s\tclause=%s\tlabels=%s\n", omic, path, state, expected, verdict,
        detail$clause, detail$labels))
    expect_identical(isTRUE(verdict), expected, info = paste("F1", omic, path, state, detail$clause, detail$labels))
}

k194F1Derives <- function(binding, manager, reference, expected, omic, path, state) {
    verdict <- k194Fn("omicNormDerivesFromCurrentRun", k194Derives)(binding, manager, reference)
    cat(sprintf("KIT-194 F1-derivation\t%s\t%s\t%s\texpected=%s\tverdict=%s\n", omic, path, state, expected, verdict))
    expect_identical(isTRUE(verdict), expected, info = paste("F1 derivation", omic, path, state))
}

k194F1Adopted <- function(manager, omic, path) {
    marked <- vapply(names(manager$.context$snapshot$step_provenance), function(state)
        isTRUE(manager$.context$snapshot$step_provenance[[state]]$adopted), NA)
    cat(sprintf("KIT-194 F1-adopted\t%s\t%s\tmarked=%s\n", omic, path, paste(names(marked)[marked], collapse = ",")))
    expect_false(any(marked), info = paste("F1 adoption marker", omic, path))
}

for (class in c("lipid", "metab")) {
    test_that(paste("KIT-194 F1:", class, "the lineage holds for every normalised state and for none before normalisation (in process)"), {
        skip_if_not(workflowOwnerCapable())
        root <- withr::local_tempdir()
        seen <- k194Notes()
        fx <- k194Fixture(file.path(root, "a"), class)
        k194Defer(fx)
        binding <- fx$binding
        prefix <- fx$prefix
        pre <- lapply(stats::setNames(nm = paste0(prefix, c("_raw_data_s4", "_intensity_filtered", "_qc_complete"))),
            function(state) k194Ref(fx, state))
        nd <- omicNormReactiveState()
        ruv <- if (identical(class, "lipid")) kit172LegacyRuv("manual")
        k194Run(fx, nd, k194Input(method = "cyclicloess", itsd = TRUE, ruv = ruv))
        k194ExpectRun(fx, nd, "F1 run 1")
        # Skip at the result, then Revert to the result, Apply, Apply again and Skip at the correlation state.
        k194Skip(fx, nd)
        k194F1Derives(binding, k194Owner(fx), k194Ref(fx), TRUE, class, "in_process", "norm_complete (Skip at the result)")
        k194Revert(fx, times = 1L)
        k194Apply(fx, nd)
        first_correlation <- k194Ref(fx, paste0(prefix, "_correlation_filtered"))
        k194F1Derives(binding, k194Owner(fx), first_correlation, TRUE, class, "in_process", "correlation_filtered (applied)")
        k194Apply(fx, nd)
        k194F1Derives(binding, k194Owner(fx), k194Ref(fx), TRUE, class, "in_process", "correlation_filtered (applied again)")
        k194Skip(fx, nd)
        k194F1Derives(binding, k194Owner(fx), k194Ref(fx), TRUE, class, "in_process", "norm_complete (Skip at correlation)")
        k194F1Derives(binding, k194Owner(fx), k194Ref(fx, paste0(prefix, "_normalized")), TRUE, class, "in_process", "normalized")
        if (!is.null(ruv)) k194F1Derives(binding, k194Owner(fx), k194Ref(fx, paste0(prefix, "_ruv_corrected")), TRUE, class, "in_process", "ruv_corrected")
        states <- paste0(prefix, c("_itsd_norm", "_log2", "_normalized", if (!is.null(ruv)) "_ruv_corrected", "_correlation_filtered", "_norm_complete"))
        for (state in states) k194F1(binding, k194Owner(fx), k194Ref(fx, state), TRUE, class, "in_process", state)
        for (state in names(pre)) k194F1(binding, k194Owner(fx), pre[[state]], FALSE, class, "in_process", paste(state, "(before normalisation)"))
        k194F1Adopted(k194Owner(fx), class, "in_process")
        # Reset of the un-normalised start: `<prefix>_reset` is not normalised.
        first_ruv <- k194Ref(fx, paste0(prefix, "_ruv_corrected"))
        k194Reset(fx, nd)
        k194F1(binding, k194Owner(fx), k194Ref(fx, paste0(prefix, "_reset")), FALSE, class, "in_process", paste0(prefix, "_reset"))
        # A second run after Reset: the first run's correlation and RUV states no longer derive from the current run (A01).
        k194Run(fx, nd, k194Input(method = "quantile", ruv = ruv))
        k194ExpectRun(fx, nd, "F1 run 2")
        k194F1Derives(binding, k194Owner(fx), first_correlation, FALSE, class, "in_process", "first run's correlation_filtered after Reset and a second run")
        if (!is.null(ruv)) k194F1Derives(binding, k194Owner(fx), first_ruv, FALSE, class, "in_process", "first run's ruv_corrected after Reset and a second run")
        # QC states made by the wrappers on a normalised head (the scripted path; old sessions).
        selectWorkflowStateReference(k194Owner(fx), paste0(prefix, "_normalized"))
        config <- fx$workflow$config_list
        kit172Quiet(runWorkflowOmicsReferenceIntensityFilter(binding, k194Owner(fx), 1, 0.5, config))
        k194F1(binding, k194Owner(fx), k194Ref(fx), TRUE, class, "in_process", "intensity filter on a normalised head")
        kit172Quiet(runWorkflowOmicsReferenceDuplicateResolution(binding, k194Owner(fx), config))
        k194F1(binding, k194Owner(fx), k194Ref(fx), TRUE, class, "in_process", "duplicate resolution on a normalised head")
        kit172Quiet(commitWorkflowOmicsStoreStep(k194Owner(fx), paste0(prefix, "_qc_finalize"), paste0(prefix, "_qc_complete"), config, "QC"))
        k194F1(binding, k194Owner(fx), k194Ref(fx), TRUE, class, "in_process", "Finalize QC on a normalised head")
        # J1's correlation state (a correlation on the QC table after a Revert, by the wrapper): not derived.
        fx2 <- k194Fixture(file.path(root, "j1"), class)
        k194Defer(fx2)
        nd2 <- omicNormReactiveState()
        k194Run(fx2, nd2)
        k194ExpectRun(fx2, nd2, "F1 J1 run")
        selectWorkflowStateReference(k194Owner(fx2), paste0(prefix, "_qc_complete"))
        kit172Quiet(runWorkflowOmicsReferenceCorrelation(binding, k194Owner(fx2), "group", 0.5, fx2$workflow$config_list))
        k194F1Derives(binding, k194Owner(fx2), k194Ref(fx2), FALSE, class, "in_process", "J1 correlation_filtered")
        k194F1(binding, k194Owner(fx2), k194Ref(fx2), FALSE, class, "in_process", "J1 correlation_filtered")
    })
}

test_that("KIT-194 F1: the lineage and the derivation on the managed worker path (the KIT-078 tabs' journey)", {
    skip_if_not(workflowOwnerCapable())
    root <- withr::local_tempdir()
    for (omic in c("lipidomics", "metabolomics")) {
        lipid <- identical(omic, "lipidomics")
        binding <- omicsDeclarationBindings()[[omic]]
        fixture <- if (lipid) lipidDeclarationsOwner(file.path(root, omic)) else metabDeclarationsOwner(file.path(root, omic))
        k194DeferOwner(fixture$owner)
        probe <- new.env()
        local({
            omicsManagedChainMocks(fixture, probe)
            assays <- names(methods::slot(fixture$object, binding$data_slot))
            id_column <- methods::slot(fixture$object, binding$id_slot)
            ids <- stats::setNames(lapply(assays, function(assay) methods::slot(fixture$object, binding$data_slot)[[assay]][[id_column]]), assays)
            selections <- list(assays = assays, candidates = ids, rows = stats::setNames(lapply(assays, function(assay) 1:2), assays))
            suppressWarnings(suppressMessages(utils::capture.output(shiny::testServer(omicsTabServer(fixture, omic, probe, TRUE), {
                session$userData$probe_owner <- fixture$owner
                omicsTabJourney(session, omic, selections)
            }))))
        })
        owner <- fixture$owner
        prefix <- binding$prefix
        history <- as.character(unlist(owner$getHistory()))
        states <- paste0(prefix, c("_itsd_norm", "_log2", "_normalized", "_ruv_corrected", "_correlation_filtered"))
        expect_true(all(states %in% history), info = paste("worker journey precondition:", omic, paste(history, collapse = " ")))
        for (state in intersect(states, history)) {
            k194F1(binding, owner, getWorkflowStateReference(owner, state), TRUE, omic, "worker", state)
        }
        for (state in intersect(paste0(prefix, c("_raw_data_s4", "_intensity_filtered", "_duplicates_resolved", "_qc_complete")), history)) {
            k194F1(binding, owner, getWorkflowStateReference(owner, state), FALSE, omic, "worker", state)
        }
        for (state in intersect(paste0(prefix, c("_ruv_corrected", "_correlation_filtered")), history)) {
            k194F1Derives(binding, owner, getWorkflowStateReference(owner, state), TRUE, omic, "worker", state)
        }
        k194F1Adopted(owner, omic, "worker")
    }
})

test_that("KIT-194 F1: transcriptomics' RUVs state is normalised, its QC states and Reset are not, and a filter on it is", {
    skip_if_not_installed("edgeR")
    project <- module_ci_trx_project(module_ci_trx_corpus())
    owner <- project$owner
    binding <- transcriptomicsDeclarationBinding
    pre <- lapply(stats::setNames(nm = c("trx_raw_data_s4", "trx_expression_filtered", "trx_qc_complete")),
        function(state) getWorkflowStateReference(owner, state))
    nd <- omicNormReactiveState()
    seen <- new.env(parent = emptyenv())
    k194TrxRun(project, nd, list(ruv_mode = "manual", ruv_k = 2L, ruv_percentage = 10), seen)
    expect_identical(owner$current_state, "trx_ruvs", info = paste("run precondition:", paste(seen$notes, collapse = " | ")))
    k194F1(binding, owner, getWorkflowStateReference(owner, "trx_ruvs"), TRUE, "transcriptomics", "in_process", "trx_ruvs")
    for (state in names(pre)) k194F1(binding, owner, pre[[state]], FALSE, "transcriptomics", "in_process", state)
    k194F1Adopted(owner, "transcriptomics", "in_process")
    kit172Quiet(shiny::isolate(resetTrxNormalisation(owner, nd, project$workflow$config_list, function(...) NULL,
        function(...) NULL, function(...) NULL)))
    k194F1(binding, owner, getWorkflowStateReference(owner, "trx_reset"), FALSE, "transcriptomics", "in_process", "trx_reset")
    selectWorkflowStateReference(owner, "trx_ruvs")
    kit172Quiet(runTrxReferenceExpressionFilter(owner, list()))
    k194F1(binding, owner, getWorkflowStateReference(owner), TRUE, "transcriptomics", "in_process", "expression filter on trx_ruvs")
})

# ------------------------------------------------------------------------------------------------
# T1, T2: the QC guard.
# ------------------------------------------------------------------------------------------------

for (class in c("lipid", "metab")) {
    test_that(paste("KIT-194 T1:", class, "QC steps refuse a normalised table, in process and managed, and run on a QC table"), {
        skip_if_not(workflowOwnerCapable())
        root <- withr::local_tempdir()
        seen <- k194Notes()
        for (route in c("in_process", "managed")) {
            fx <- k194Fixture(file.path(root, route), class)
            k194Defer(fx)
            prefix <- fx$prefix
            rec <- k194Runner(fx$workflow)
            runner <- if (identical(route, "managed")) rec$runner
            nd <- omicNormReactiveState()
            ruv <- if (identical(class, "lipid")) kit172LegacyRuv("manual")
            k194Run(fx, nd, k194Input(ruv = ruv))
            k194ExpectRun(fx, nd, paste("T1", route))
            k194Apply(fx, nd)
            # KIT-254 (set-up only): each correlation-step head is met in its own sequence (a Skip at the correlation
            # state now returns to the run's result first); the heads and every refusal expectation are unchanged.
            result <- paste0(prefix, if (!is.null(ruv)) "_ruv_corrected" else "_normalized")
            expect_true(all(paste0(prefix, c("_correlation_filtered", if (!is.null(ruv)) "_ruv_corrected", "_normalized", "_log2")) %in% k194History(fx)),
                info = paste("T1 precondition (after the Apply):", paste(k194History(fx), collapse = " ")))
            refused <- function(head, label) {
                for (button in c("apply", "resolve", "finalize")) {
                    if (!is.null(head)) selectWorkflowStateReference(k194Owner(fx), head)
                    before <- k194Snapshot(fx)
                    seen$notes <- character()
                    rec$operations <- character()
                    k194QcPress(fx, button, runner)
                    info <- paste("T1", route, label, button, paste(seen$notes, collapse = " | "))
                    expect_identical(k194Snapshot(fx), before, info = info)
                    expect_true(k194Noted(seen, k194QcLead) && k194Noted(seen, k194Step(prefix)), info = info)
                    if (identical(route, "managed")) expect_length(rec$operations, 0L)
                }
            }
            refused(paste0(prefix, "_correlation_filtered"), paste0(prefix, "_correlation_filtered"))
            selectWorkflowStateReference(k194Owner(fx), result)
            k194Skip(fx, nd)
            expect_true(paste0(prefix, "_norm_complete") %in% k194History(fx),
                info = paste("T1 precondition (after the Skip):", paste(k194History(fx), collapse = " ")))
            for (head in paste0(prefix, c("_norm_complete", if (!is.null(ruv)) "_ruv_corrected", "_normalized"))) refused(head, head)
            # A QC state made by the wrapper on the normalised table (an old session's sequence, the scripted path).
            selectWorkflowStateReference(k194Owner(fx), paste0(prefix, "_normalized"))
            kit172Quiet(runWorkflowOmicsReferenceIntensityFilter(fx$binding, k194Owner(fx), 1, 0.5, fx$workflow$config_list))
            refused(NULL, "intensity filter on the normalised table")
            # A head inside the normalisation (one Revert there): `<prefix>_log2`, which precedes the states above.
            refused(paste0(prefix, "_log2"), paste0(prefix, "_log2"))
            # A QC table: the steps run (both commits).
            selectWorkflowStateReference(k194Owner(fx), paste0(prefix, "_qc_complete"))
            qc_reference <- k194Ref(fx)
            for (step in list(c("apply", "_intensity_filtered"), c("resolve", "_duplicates_resolved"), c("finalize", "_qc_complete"))) {
                seen$notes <- character()
                k194QcPress(fx, step[[1L]], runner)
                expect_identical(k194Head(fx), paste0(prefix, step[[2L]]), info = paste("T1 allowed", route, step[[1L]], paste(seen$notes, collapse = " | ")))
            }
            kit172Quiet(commitWorkflowOmicsStoreStep(k194Owner(fx), paste0(prefix, "_reset"), paste0(prefix, "_reset"),
                fx$workflow$config_list, "Reset to pre-normalization state", source = qc_reference))
            k194QcPress(fx, "apply", runner)
            expect_identical(k194Head(fx), paste0(prefix, "_intensity_filtered"), info = paste("T1 allowed after Reset", route))
            if (identical(route, "managed")) {
                # A02: a press queued behind a running action is judged when it starts. The press is made at a QC table,
                # a normalisation then commits, and the queued filter starts on its result.
                rec$hold <- TRUE
                seen$notes <- character()
                k194QcPress(fx, "apply", runner)
                expect_length(rec$queued, 1L)
                k194Run(fx, nd, k194Input())
                k194ExpectRun(fx, nd, "T1 A02")
                before <- k194Snapshot(fx)
                tryCatch(kit172Quiet(rec$start()), error = function(e) seen$notes <- c(seen$notes, paste("[start error]", conditionMessage(e))))
                info <- paste("T1 A02", paste(seen$notes, collapse = " | "))
                expect_identical(k194Snapshot(fx), before, info = info)
                expect_true(k194Noted(seen, k194QcLead), info = info)
            }
        }
        # A08: a session loaded by Load Filtered Session (adopted, no lineage) under `loaded_for_de` refuses QC.
        fx <- k194Fixture(file.path(root, "loaded"), class)
        k194Defer(fx)
        k194LoadLegacy(fx, "loaded_for_de", root)
        before <- k194Snapshot(fx)
        seen$notes <- character()
        k194QcPress(fx, "apply")
        info <- paste("T1 A08", paste(seen$notes, collapse = " | "))
        expect_identical(k194Snapshot(fx), before, info = info)
        expect_true(k194Noted(seen, k194QcLead) && k194Noted(seen, k194NoStep), info = info)
    })
}

test_that("KIT-194 T2: transcriptomics' Expression Filter and Finalize QC refuse trx_ruvs, and run on trx_qc_complete", {
    skip_if_not_installed("edgeR")
    seen <- k194Notes()
    project <- module_ci_trx_project(module_ci_trx_corpus())
    owner <- project$owner
    nd <- omicNormReactiveState()
    k194TrxRun(project, nd, list(ruv_mode = "manual", ruv_k = 2L, ruv_percentage = 10), seen)
    expect_identical(owner$current_state, "trx_ruvs", info = paste("run precondition:", paste(seen$notes, collapse = " | ")))
    for (button in c("apply", "finalize")) {
        selectWorkflowStateReference(owner, "trx_ruvs")
        before <- list(history = as.character(unlist(owner$getHistory())), head = owner$current_state)
        seen$notes <- character()
        k194TrxPress(project, button)
        info <- paste("T2", button, paste(seen$notes, collapse = " | "), "head", owner$current_state)
        expect_identical(list(history = as.character(unlist(owner$getHistory())), head = owner$current_state), before, info = info)
        expect_true(k194Noted(seen, k194QcLead) && k194Noted(seen, k194Step("trx")), info = info)
        expect_true(button != "apply" || k194Noted(seen, "Expression Filter sub-tab"), info = info)
        expect_false(k194Noted(seen, "omic_qc_input_normalised"), info = info)
        if (identical(owner$current_state, "trx_expression_filtered")) {
            cat(sprintf("KIT-194 T2 base\ttrx_expression_filtered committed on trx_ruvs with W columns: %s\n",
                paste(names(k194TrxW(owner, "trx_expression_filtered")), collapse = ",")))
        }
    }
    selectWorkflowStateReference(owner, "trx_qc_complete")
    k194TrxPress(project, "apply")
    expect_identical(owner$current_state, "trx_expression_filtered", info = paste("T2 allowed", paste(seen$notes, collapse = " | ")))
    # Managed (A02, B02): a press queued at a QC table, judged when it starts on the normalised table.
    selectWorkflowStateReference(owner, "trx_qc_complete")
    rec <- k194Runner(project$workflow, hold = TRUE)
    k194TrxPress(project, "apply", runner = rec$runner)
    expect_length(rec$queued, 1L)
    k194TrxRun(project, omicNormReactiveState(), list(ruv_mode = "manual", ruv_k = 2L, ruv_percentage = 10), seen)
    expect_identical(owner$current_state, "trx_ruvs", info = "T2 managed precondition")
    before <- list(history = as.character(unlist(owner$getHistory())), head = owner$current_state)
    seen$notes <- character()
    tryCatch(kit172Quiet(rec$start()), error = function(e) seen$notes <- c(seen$notes, paste("[start error]", conditionMessage(e))))
    info <- paste("T2 managed", paste(seen$notes, collapse = " | "))
    expect_identical(list(history = as.character(unlist(owner$getHistory())), head = owner$current_state), before, info = info)
    expect_true(k194Noted(seen, k194QcLead), info = info)
})

# ------------------------------------------------------------------------------------------------
# T3, T4, T5: the run rule.
# ------------------------------------------------------------------------------------------------

#' Two runs of the same settings; the second pressed at `at` (the first run's result, or after `between`).
k194TwoRuns <- function(fx, nd, input, runner = NULL, between = NULL) {
    k194Run(fx, nd, input, runner)
    k194ExpectRun(fx, nd, "first run")
    first <- list(start = shiny::isolate(nd$post_filter_obj), history = k194History(fx),
        normalised = k194Matrix(fx, paste0(fx$prefix, "_normalized")), log = length(fx$log$lines))
    if (is.function(between)) between()
    k194Run(fx, nd, input, runner)
    second <- list(start = shiny::isolate(nd$post_filter_obj), history = k194History(fx),
        normalised = k194Matrix(fx, paste0(fx$prefix, "_normalized")), log = fx$log$lines[-seq_len(first$log)])
    list(first = first, second = second)
}

k194ExpectRerun <- function(fx, runs, info, history = runs$first$history) {
    prefix <- fx$prefix
    line <- sprintf("'%s_qc_complete'", prefix)
    expect_true(any(grepl("runs again from", runs$second$log, fixed = TRUE) & grepl(line, runs$second$log, fixed = TRUE)),
        info = paste(info, paste(runs$second$log, collapse = " | ")))
    expect_identical(runs$second$start$entity_id, runs$first$start$entity_id, info = info)
    expect_identical(runs$second$history, history, info = info)
    expect_identical(runs$second$normalised, runs$first$normalised, info = info)
}

for (class in c("lipid", "metab")) {
    test_that(paste("KIT-194 T3 (a, b):", class, "a second run at the result runs again from the table the first started from"), {
        skip_if_not(workflowOwnerCapable())
        root <- withr::local_tempdir()
        seen <- k194Notes()
        for (route in c("in_process", "managed")) {
            fx <- k194Fixture(file.path(root, route), class)
            k194Defer(fx)
            rec <- k194Runner(fx$workflow)
            runner <- if (identical(route, "managed")) rec$runner
            runs <- k194TwoRuns(fx, omicNormReactiveState(), k194Input(), runner)
            cat(sprintf("KIT-194 effect\t%s\tF-c\t%s\tfirst_normalised_max\t%.6g\tsecond_normalised_max\t%.6g\tsecond_start\t%s\n", class, route,
                max(runs$first$normalised, na.rm = TRUE), max(runs$second$normalised, na.rm = TRUE),
                paste(names(Filter(function(r) identical(r$entity_id, runs$second$start$entity_id), k194Owner(fx)$.context$snapshot$states)), collapse = "|")))
            k194ExpectRerun(fx, runs, paste("T3 (a)", class, route))
        }
        # (b) Run, Apply, run again.
        fx <- k194Fixture(file.path(root, "b"), class)
        k194Defer(fx)
        nd <- omicNormReactiveState()
        runs <- k194TwoRuns(fx, nd, k194Input(), between = function() {
            k194Apply(fx, nd)
            expect_identical(k194Head(fx), paste0(fx$prefix, "_correlation_filtered"), info = "T3 (b) precondition")
        })
        k194ExpectRerun(fx, runs, paste("T3 (b)", class))
    })

    test_that(paste("KIT-194 T3 (c, d):", class, "ITSD by pattern runs again from the start; a fresh tab at the result refuses"), {
        skip_if_not(workflowOwnerCapable())
        root <- withr::local_tempdir()
        seen <- k194Notes()
        fx <- k194Fixture(file.path(root, "c"), class)
        k194Defer(fx)
        runs <- k194TwoRuns(fx, omicNormReactiveState(), k194Input(itsd = TRUE))
        errors <- grep("^ERROR", runs$second$log, value = TRUE)
        cat(sprintf("KIT-194 F3\t%s\tsecond run with ITSD by pattern: error=%s\tsecond_normalised_max=%.6g\tfirst_normalised_max=%.6g\n", class,
            if (length(errors)) errors[[1L]] else "none", max(runs$second$normalised, na.rm = TRUE), max(runs$first$normalised, na.rm = TRUE)))
        k194ExpectRerun(fx, runs, paste("T3 (c)", class))
        # (d) A reopened tab (its state new) at the head `<prefix>_normalized`.
        fx <- k194Fixture(file.path(root, "d"), class)
        k194Defer(fx)
        nd <- omicNormReactiveState()
        k194Run(fx, nd)
        k194ExpectRun(fx, nd, "T3 (d)")
        before <- k194Snapshot(fx)
        seen$notes <- character()
        k194Run(fx, omicNormReactiveState())
        info <- paste("T3 (d)", class, paste(seen$notes, collapse = " | "))
        expect_identical(k194Snapshot(fx), before, info = info)
        expect_true(k194Noted(seen, k194RunUnknownLead) && k194Noted(seen, k194Step(fx$prefix)), info = info)
    })

    test_that(paste("KIT-194 T3 (e, f, g):", class, "a run is refused while the runner is busy, and on a loaded session"), {
        skip_if_not(workflowOwnerCapable())
        root <- withr::local_tempdir()
        seen <- k194Notes()
        # (e) A02: managed, the runner busy.
        fx <- k194Fixture(file.path(root, "e"), class)
        k194Defer(fx)
        rec <- k194Runner(fx$workflow, busy = TRUE)
        before <- k194Snapshot(fx)
        k194Run(fx, omicNormReactiveState(), k194Input(), rec$runner)
        info <- paste("T3 (e)", class, paste(seen$notes, collapse = " | "))
        expect_length(rec$queued, 0L)
        expect_identical(k194Snapshot(fx), before, info = info)
        expect_true(k194Noted(seen, k194Busy), info = info)
        # (f) A07, (g) A08: Load Filtered Session of a legacy file, then Run.
        for (state in c(paste0(kit172Binding(class)$prefix, "_correlation_filtered"), "loaded_for_de")) {
            fx <- k194Fixture(file.path(root, state), class)
            k194Defer(fx)
            k194LoadLegacy(fx, state, root)
            before <- k194Snapshot(fx)
            seen$notes <- character()
            k194Run(fx, omicNormReactiveState())
            info <- paste("T3 (f, g)", class, state, paste(seen$notes, collapse = " | "))
            expect_identical(k194Snapshot(fx), before, info = info)
            expect_true(k194Noted(seen, k194NoStep), info = info)
        }
    })
}

test_that("KIT-194 T4: transcriptomics' second run (k 2, then k 1) runs again from trx_qc_complete; a fresh tab and a busy runner refuse", {
    skip_if_not_installed("edgeR")
    seen <- new.env(parent = emptyenv())
    seen$notes <- character()
    project <- module_ci_trx_project(module_ci_trx_corpus())
    owner <- project$owner
    nd <- omicNormReactiveState()
    k194TrxRun(project, nd, list(ruv_mode = "manual", ruv_k = 2L, ruv_percentage = 10), seen)
    expect_identical(owner$current_state, "trx_ruvs", info = paste("run precondition:", paste(seen$notes, collapse = " | ")))
    first_start <- shiny::isolate(nd$post_filter_obj)
    lines <- length(k194TrxLog(nd))
    seen$notes <- character()
    k194TrxRun(project, nd, list(ruv_mode = "manual", ruv_k = 1L, ruv_percentage = 10), seen)
    log <- k194TrxLog(nd)[-seq_len(lines)]
    errors <- grep("[[]error", seen$notes, value = TRUE)
    second_w <- tryCatch(names(k194TrxW(owner)), error = function(e) paste("unreadable:", conditionMessage(e)))
    cat(sprintf("KIT-194 F2\tsecond legacy-path run k 2 then k 1: error=%s\tW columns=%s\n",
        if (length(errors)) errors[[1L]] else "none", paste(second_w, collapse = ",")))
    info <- paste("T4", paste(seen$notes, collapse = " | "), paste(log, collapse = " | "))
    expect_length(errors, 0L)
    expect_true(any(grepl("runs again from 'trx_qc_complete'", log, fixed = TRUE)), info = info)
    expect_identical(shiny::isolate(nd$post_filter_obj)$entity_id, first_start$entity_id, info = info)
    expect_identical(second_w, "W_1", info = info)
    # Equal to a single k 1 run from trx_qc_complete.
    single <- module_ci_trx_project(module_ci_trx_corpus())
    k194TrxRun(single, omicNormReactiveState(), list(ruv_mode = "manual", ruv_k = 1L, ruv_percentage = 10), new.env())
    if (identical(second_w, "W_1")) expect_equal(k194TrxW(owner), k194TrxW(single$owner), info = info)
    # (d) A fresh tab at trx_ruvs.
    before <- list(history = as.character(unlist(owner$getHistory())), head = owner$current_state)
    seen$notes <- character()
    k194TrxRun(project, omicNormReactiveState(), list(ruv_mode = "manual", ruv_k = 1L, ruv_percentage = 10), seen)
    info <- paste("T4 (d)", paste(seen$notes, collapse = " | "))
    expect_identical(list(history = as.character(unlist(owner$getHistory())), head = owner$current_state), before, info = info)
    expect_true(k194Noted(seen, k194RunUnknownLead) && k194Noted(seen, k194Step("trx")), info = info)
    # (e) A02: the runner busy.
    rec <- k194Runner(project$workflow, busy = TRUE)
    seen$notes <- character()
    k194TrxRun(project, nd, list(ruv_mode = "manual", ruv_k = 1L, ruv_percentage = 10), seen, runner = rec$runner)
    info <- paste("T4 (e)", paste(seen$notes, collapse = " | "))
    expect_length(rec$queued, 0L)
    expect_identical(list(history = as.character(unlist(owner$getHistory())), head = owner$current_state), before, info = info)
    expect_true(k194Noted(seen, k194Busy), info = info)
})

test_that("KIT-194 T5: a QC state made on a normalised table (an old session) is not normalised again", {
    skip_if_not(workflowOwnerCapable())
    root <- withr::local_tempdir()
    seen <- k194Notes()
    for (class in c("lipid", "metab")) {
        fx <- k194Fixture(file.path(root, class), class)
        k194Defer(fx)
        nd <- omicNormReactiveState()
        k194Run(fx, nd)
        k194ExpectRun(fx, nd, "T5")
        kit172Quiet(runWorkflowOmicsReferenceIntensityFilter(fx$binding, k194Owner(fx), 1, 0.5, fx$workflow$config_list))
        before <- k194Snapshot(fx)
        seen$notes <- character()
        k194Run(fx, nd)
        info <- paste("T5", class, paste(seen$notes, collapse = " | "))
        expect_identical(k194Snapshot(fx), before, info = info)
        expect_true(k194Noted(seen, k194RunLead) && k194Noted(seen, k194Step(fx$prefix)), info = info)
    }
    skip_if_not_installed("edgeR")
    project <- module_ci_trx_project(module_ci_trx_corpus())
    owner <- project$owner
    nd <- omicNormReactiveState()
    trx_seen <- new.env(parent = emptyenv())
    k194TrxRun(project, nd, list(ruv_mode = "manual", ruv_k = 2L, ruv_percentage = 10), trx_seen)
    expect_identical(owner$current_state, "trx_ruvs", info = "T5 trx run precondition")
    kit172Quiet(runTrxReferenceExpressionFilter(owner, list()))
    before <- list(history = as.character(unlist(owner$getHistory())), head = owner$current_state)
    trx_seen$notes <- character()
    k194TrxRun(project, nd, list(ruv_mode = "manual", ruv_k = 1L, ruv_percentage = 10), trx_seen)
    info <- paste("T5 trx", paste(trx_seen$notes, collapse = " | "))
    expect_identical(list(history = as.character(unlist(owner$getHistory())), head = owner$current_state), before, info = info)
    expect_true(k194Noted(trx_seen, k194RunLead) && k194Noted(trx_seen, k194Step("trx")), info = info)
})

# ------------------------------------------------------------------------------------------------
# T6, T7: the correlation, Skip and export guard (lipidomics, metabolomics).
# ------------------------------------------------------------------------------------------------

for (class in c("lipid", "metab")) {
    test_that(paste("KIT-194 T6 (a, b, c):", class, "after a Revert to the QC table, Apply, Skip and the managed Apply are refused"), {
        skip_if_not(workflowOwnerCapable())
        root <- withr::local_tempdir()
        seen <- k194Notes()
        fx <- k194Fixture(file.path(root, "a"), class)
        k194Defer(fx)
        prefix <- fx$prefix
        nd <- omicNormReactiveState()
        k194Run(fx, nd, k194Input(method = "cyclicloess", itsd = TRUE))
        k194ExpectRun(fx, nd, "T6")
        normalised_max <- k194Max(fx, paste0(prefix, "_normalized"))
        qc_matrix <- k194Matrix(fx, paste0(prefix, "_qc_complete"))
        expect_identical(k194Revert(fx, paste0(prefix, "_qc_complete")), paste0(prefix, "_qc_complete"))
        status <- fx$workflow$tab_status$normalization
        flags <- shiny::isolate(list(nd$normalization_complete, nd$correlation_filtering_complete))
        # (a) Apply.
        before <- k194Snapshot(fx)
        seen$notes <- character()
        k194Apply(fx, nd)
        info <- paste("T6 (a)", class, paste(seen$notes, collapse = " | "))
        if (identical(k194Head(fx), paste0(prefix, "_correlation_filtered"))) {
            committed <- k194Matrix(fx, NULL)
            cat(sprintf("KIT-194 effect\t%s\tF-j J1\tDA-input max\t%.6g\tnormalised max\t%.6g\tequal to qc_complete on kept samples\t%s\n",
                class, max(committed, na.rm = TRUE), normalised_max,
                identical(unname(committed), unname(qc_matrix[rownames(committed), colnames(committed), drop = FALSE]))))
            selectWorkflowStateReference(k194Owner(fx), paste0(prefix, "_qc_complete"))
        }
        expect_identical(k194Snapshot(fx)[c("history", "head")], before[c("history", "head")], info = info)
        expect_true(k194Noted(seen, paste("Error:", k194Stale)), info = info)
        expect_identical(fx$workflow$tab_status$normalization, status, info = info)
        expect_identical(shiny::isolate(list(nd$normalization_complete, nd$correlation_filtering_complete)), flags, info = info)
        # (b) Skip.
        before <- k194Snapshot(fx)
        seen$notes <- character()
        value <- k194Skip(fx, nd)
        info <- paste("T6 (b)", class, paste(seen$notes, collapse = " | "))
        expect_false(isTRUE(value), info = info)
        expect_identical(k194Snapshot(fx)[c("history", "head")], before[c("history", "head")], info = info)
        expect_true(k194Noted(seen, k194Stale), info = info)
        if (!identical(k194Head(fx), before$head)) selectWorkflowStateReference(k194Owner(fx), before$head)
        # (c) Managed Apply: nothing dispatched.
        rec <- k194Runner(fx$workflow)
        seen$notes <- character()
        k194Apply(fx, nd, runner = rec$runner)
        info <- paste("T6 (c)", class, paste(rec$operations, collapse = " "), paste(seen$notes, collapse = " | "))
        expect_length(rec$operations, 0L)
        expect_true(k194Noted(seen, k194Stale), info = info)
    })

    test_that(paste("KIT-194 T6 (d, e, f, g):", class, "the correlation step at the result runs; off it (J6, A01) it is refused"), {
        skip_if_not(workflowOwnerCapable())
        root <- withr::local_tempdir()
        seen <- k194Notes()
        # (f) After an in-process and a managed run the tab's final is the head's saved reference (both commits).
        for (route in c("in_process", "managed")) {
            fx <- k194Fixture(file.path(root, paste0("f-", route)), class)
            k194Defer(fx)
            rec <- k194Runner(fx$workflow)
            nd <- omicNormReactiveState()
            k194Run(fx, nd, k194Input(), if (identical(route, "managed")) rec$runner)
            k194ExpectRun(fx, nd, paste("T6 (f)", route))
            expect_true(isSavedWorkflowReference(shiny::isolate(nd$ruv_corrected_obj)), info = route)
            expect_true(isSavedWorkflowReference(shiny::isolate(nd$post_norm_obj)), info = route)
            expect_identical(shiny::isolate(nd$ruv_corrected_obj)$entity_id, k194Ref(fx)$entity_id, info = route)
        }
        # (e) Unaffected (both commits): Apply at the result, Apply again and Skip at the correlation state, Revert to the
        # result and Apply, Revert to the result and Skip.
        prefix <- fx$prefix
        expect_true(isTRUE(k194Apply(fx, nd)), info = paste("T6 (e) Apply at the result", paste(seen$notes, collapse = " | ")))
        expect_identical(k194Head(fx), paste0(prefix, "_correlation_filtered"))
        expect_true(isTRUE(k194Apply(fx, nd)), info = "T6 (e) Apply again")
        expect_identical(k194Head(fx), paste0(prefix, "_correlation_filtered"))
        expect_true(isTRUE(k194Skip(fx, nd)), info = "T6 (e) Skip at the correlation state")
        expect_identical(k194Head(fx), paste0(prefix, "_norm_complete"))
        expect_identical(k194Revert(fx, paste0(prefix, "_normalized")), paste0(prefix, "_normalized"))
        expect_true(isTRUE(k194Apply(fx, nd)), info = "T6 (e) Revert to the result, Apply")
        expect_identical(k194Head(fx), paste0(prefix, "_correlation_filtered"))
        expect_identical(k194Revert(fx, paste0(prefix, "_normalized")), paste0(prefix, "_normalized"))
        expect_true(isTRUE(k194Skip(fx, nd)), info = "T6 (e) Skip at the result")
        expect_identical(k194Head(fx), paste0(prefix, "_norm_complete"))
        expect_false(k194Noted(seen, "[error]"), info = paste("T6 (e)", paste(seen$notes, collapse = " | ")))
        # (d) J6 (lipid): a RUV run, one Revert (head `_normalized`), Apply.
        if (identical(class, "lipid")) {
            fx <- k194Fixture(file.path(root, "d"), class)
            k194Defer(fx)
            nd <- omicNormReactiveState()
            k194Run(fx, nd, k194Input(ruv = kit172LegacyRuv("manual")))
            k194ExpectRun(fx, nd, "T6 (d)")
            expect_identical(k194Head(fx), paste0(prefix, "_ruv_corrected"), info = "T6 (d) RUV precondition")
            expect_identical(k194Revert(fx, times = 1L), paste0(prefix, "_normalized"))
            before <- k194Snapshot(fx)
            seen$notes <- character()
            k194Apply(fx, nd)
            info <- paste("T6 (d)", paste(seen$notes, collapse = " | "))
            expect_identical(k194Snapshot(fx)[c("history", "head")], before[c("history", "head")], info = info)
            expect_true(k194Noted(seen, k194Stale), info = info)
        }
        # (g) A01: run, Apply, Reset, a second run (cyclic loess), one Revert: the first run's correlation state is the head.
        fx <- k194Fixture(file.path(root, "g"), class)
        k194Defer(fx)
        nd <- omicNormReactiveState()
        k194Run(fx, nd, k194Input())
        k194ExpectRun(fx, nd, "T6 (g) run 1")
        k194Apply(fx, nd)
        k194Reset(fx, nd)
        k194Run(fx, nd, k194Input(method = "cyclicloess"))
        k194ExpectRun(fx, nd, "T6 (g) run 2")
        second <- k194Matrix(fx, paste0(prefix, "_normalized"))
        expect_identical(k194Revert(fx, times = 1L), paste0(prefix, "_correlation_filtered"), info = "T6 (g) precondition")
        for (action in c("skip", "apply")) {
            before <- k194Snapshot(fx)
            seen$notes <- character()
            if (identical(action, "skip")) k194Skip(fx, nd) else k194Apply(fx, nd)
            info <- paste("T6 (g)", class, action, paste(seen$notes, collapse = " | "))
            if (identical(action, "skip") && identical(k194Head(fx), paste0(prefix, "_norm_complete"))) {
                committed <- k194Matrix(fx, NULL)
                kept <- intersect(colnames(committed), colnames(second))
                rows <- intersect(rownames(committed), rownames(second))
                cat(sprintf("KIT-194 T6 (g) base\t%s\tnorm_complete equals the second run's normalised table on kept samples: %s\n", class,
                    identical(unname(committed[rows, kept]), unname(second[rows, kept]))))
                selectWorkflowStateReference(k194Owner(fx), before$head)
            }
            expect_identical(k194Snapshot(fx)[c("history", "head")], before[c("history", "head")], info = info)
            expect_true(k194Noted(seen, k194Stale), info = info)
        }
        # T8 case (A01): the DA on that head, status complete.
        k194SetNormStatus(fx, "complete")
        seen$notes <- character()
        da <- k194Da(fx)
        info <- paste("T8 (A01)", class, paste(seen$notes, collapse = " | "))
        expect_identical(da$runs, 0L, info = info)
        expect_true(k194Noted(seen, k194DaText), info = info)
    })

    test_that(paste("KIT-194 T7:", class, "Export Normalized Data refuses a reverted head and writes at the correlation state"), {
        skip_if_not(workflowOwnerCapable())
        root <- withr::local_tempdir()
        seen <- k194Notes()
        fx <- k194Fixture(file.path(root, "a"), class)
        k194Defer(fx)
        prefix <- fx$prefix
        nd <- omicNormReactiveState()
        k194Run(fx, nd)
        k194ExpectRun(fx, nd, "T7")
        k194Apply(fx, nd)
        export_dir <- file.path(root, "export")
        dir.create(export_dir)
        latest <- file.path(export_dir, paste0(prefix, "_filtered_session_data_latest.rds"))
        export <- function() kit172Quiet(shiny::isolate(exportOmicNormSession(omicsSessionInputK194(), fx$workflow,
            list(source_dir = export_dir), "kit194", nd, fx$addLog, fx$config)))
        expect_identical(k194Revert(fx, paste0(prefix, "_qc_complete")), paste0(prefix, "_qc_complete"))
        seen$notes <- character()
        export()
        info <- paste("T7", class, paste(seen$notes, collapse = " | "))
        written <- file.exists(latest)
        if (written) {
            cat(sprintf("KIT-194 T7 base\t%s\tthe session file names %s as its current state\n", class, readRDS(latest)$r6_current_state_name))
            unlink(export_dir, recursive = TRUE)
            dir.create(export_dir)
        }
        expect_false(written, info = info)
        expect_true(k194Noted(seen, paste("Export error:", k194Stale)), info = info)
        # At the correlation state the export is written (both commits): a run from the QC table, then Apply.
        k194Run(fx, nd)
        k194ExpectRun(fx, nd, "T7 second run")
        k194Apply(fx, nd)
        expect_identical(k194Head(fx), paste0(prefix, "_correlation_filtered"))
        seen$notes <- character()
        export()
        expect_true(file.exists(latest), info = paste("T7 at the correlation state", paste(seen$notes, collapse = " | ")))
    })
}

# ------------------------------------------------------------------------------------------------
# T8: the DA guard (lipidomics, metabolomics; D2, D4 (b1)).
# ------------------------------------------------------------------------------------------------

for (class in c("lipid", "metab")) {
    test_that(paste("KIT-194 T8:", class, "the DA's head fallback takes only the normalisation's states; loaded QC-named files are refused"), {
        skip_if_not(workflowOwnerCapable())
        root <- withr::local_tempdir()
        seen <- k194Notes()
        fx <- k194Fixture(file.path(root, "a"), class)
        k194Defer(fx)
        prefix <- fx$prefix
        nd <- omicNormReactiveState()
        k194Run(fx, nd, k194Input(itsd = TRUE))
        k194ExpectRun(fx, nd, "T8")
        k194Apply(fx, nd)
        k194SetNormStatus(fx, "complete")
        da <- function(head, expected) {
            if (!identical(k194Head(fx), head)) selectWorkflowStateReference(k194Owner(fx), head)
            seen$notes <- character()
            out <- k194Da(fx)
            info <- paste("T8", class, head, paste(seen$notes, collapse = " | "))
            expect_identical(out$runs, expected, info = info)
            if (identical(expected, 0L)) expect_true(k194Noted(seen, k194DaText), info = info)
        }
        # Accepted (both commits), latest first so that each selection keeps the next in the history. KIT-254 (set-up
        # only): the correlation state first, then the Skip pressed at the run's result for `_norm_complete`.
        da(paste0(prefix, "_correlation_filtered"), 1L)
        selectWorkflowStateReference(k194Owner(fx), paste0(prefix, "_normalized"))
        k194Skip(fx, nd)
        for (head in paste0(prefix, c("_norm_complete", "_normalized"))) da(head, 1L)
        # Refused heads: a state inside the normalisation, and the QC table after a Revert.
        for (head in paste0(prefix, c("_log2", "_itsd_norm", "_qc_complete"))) da(head, 0L)
        # After Reset (no Revert needed): `<prefix>_reset`.
        fx2 <- k194Fixture(file.path(root, "reset"), class)
        k194Defer(fx2)
        nd2 <- omicNormReactiveState()
        k194Run(fx2, nd2)
        k194ExpectRun(fx2, nd2, "T8 reset")
        k194Apply(fx2, nd2)
        k194Reset(fx2, nd2)
        k194SetNormStatus(fx2, "complete")
        expect_identical(k194Head(fx2), paste0(prefix, "_reset"))
        seen$notes <- character()
        out <- k194Da(fx2)
        expect_identical(out$runs, 0L, info = paste("T8 reset", paste(seen$notes, collapse = " | ")))
        expect_true(k194Noted(seen, k194DaText))
        # Lipid's RUV state.
        if (identical(class, "lipid")) {
            fx3 <- k194Fixture(file.path(root, "ruv"), class)
            k194Defer(fx3)
            nd3 <- omicNormReactiveState()
            k194Run(fx3, nd3, k194Input(ruv = kit172LegacyRuv("manual")))
            k194ExpectRun(fx3, nd3, "T8 ruv")
            k194SetNormStatus(fx3, "complete")
            expect_identical(k194Head(fx3), paste0(prefix, "_ruv_corrected"))
            expect_identical(k194Da(fx3)$runs, 1L)
        }
        # A bus without an owner, and transcriptomics: as before (both commits).
        bare <- list(contrasts_tbl = kit172Contrasts(), tab_status = list(normalization = "complete"))
        testthat::local_mocked_bindings(omicDaReadSource = function(...) list(class = fx$binding$class), .package = "MultiScholaR")
        expect_identical(k194Da(fx, workflow = bare)$runs, 1L)
        trx_config <- omicGuiConfig("transcriptomics")
        testthat::local_mocked_bindings(omicDaReadSource = function(...) list(class = trx_config$da$source$class), .package = "MultiScholaR")
        expect_identical(k194Da(fx, workflow = bare, binding = transcriptomicsDeclarationBinding, da = trx_config$da)$runs, 1L)
    })

    test_that(paste("KIT-194 T8 (D4 (b1)):", class, "a session loaded under a QC state's name is refused at Run DA"), {
        skip_if_not(workflowOwnerCapable())
        root <- withr::local_tempdir()
        seen <- k194Notes()
        prefix <- kit172Binding(class)$prefix
        for (case in list(c(paste0(prefix, "_qc_complete"), "0"), c("loaded_for_de", "1"), c(paste0(prefix, "_correlation_filtered"), "1"))) {
            fx <- k194Fixture(file.path(root, case[[1L]]), class)
            k194Defer(fx)
            tab <- k194LoadLegacy(fx, case[[1L]], root)
            seen$notes <- character()
            out <- k194Da(fx, tab = tab)
            info <- paste("T8 (b1)", class, case[[1L]], paste(seen$notes, collapse = " | "))
            expect_identical(out$runs, as.integer(case[[2L]]), info = info)
            if (identical(case[[2L]], "0")) expect_true(k194Noted(seen, k194DaText), info = info)
        }
    })
}

# ------------------------------------------------------------------------------------------------
# T9: the default path (both commits).
# ------------------------------------------------------------------------------------------------

for (class in c("lipid", "metab")) {
    test_that(paste("KIT-194 T9:", class, "QC, a run, Apply (or Skip) and the DA run as before"), {
        skip_if_not(workflowOwnerCapable())
        root <- withr::local_tempdir()
        seen <- k194Notes()
        for (finish in c("apply", "skip")) {
            fx <- k194Fixture(file.path(root, finish), class)
            k194Defer(fx)
            prefix <- fx$prefix
            # The QC presses from the design state.
            selectWorkflowStateReference(k194Owner(fx), paste0(prefix, "_raw_data_s4"))
            k194QcPress(fx, "apply")
            k194QcPress(fx, "finalize")
            expect_identical(k194History(fx), c("initial", paste0(prefix, c("_raw_data_s4", "_intensity_filtered", "_qc_complete"))))
            nd <- omicNormReactiveState()
            k194Run(fx, nd, k194Input(method = "cyclicloess", itsd = TRUE))
            k194ExpectRun(fx, nd, "T9")
            if (identical(finish, "apply")) expect_true(isTRUE(k194Apply(fx, nd))) else expect_true(isTRUE(k194Skip(fx, nd)))
            last <- paste0(prefix, if (identical(finish, "apply")) "_correlation_filtered" else "_norm_complete")
            expect_identical(k194History(fx), c("initial", paste0(prefix, c("_raw_data_s4", "_intensity_filtered", "_qc_complete",
                "_itsd_norm", "_log2", "_normalized")), last))
            expect_identical(fx$workflow$tab_status$normalization, "complete")
            out <- k194Da(fx)
            expect_identical(out$runs, 1L)
            expect_identical(out$state, last)
            expect_false(k194Noted(seen, "[error]"), info = paste("T9", class, finish, paste(seen$notes, collapse = " | ")))
        }
    })
}

test_that("KIT-194 T9: transcriptomics' QC, run and DA run as before", {
    skip_if_not_installed("edgeR")
    seen <- k194Notes()
    project <- module_ci_trx_project(module_ci_trx_corpus())
    owner <- project$owner
    selectWorkflowStateReference(owner, "trx_raw_data_s4")
    k194TrxPress(project, "apply")
    k194TrxPress(project, "finalize")
    expect_identical(as.character(unlist(owner$getHistory())), c("initial", "trx_raw_data_s4", "trx_expression_filtered", "trx_qc_complete"))
    nd <- omicNormReactiveState()
    k194TrxRun(project, nd, list(ruv_mode = "manual", ruv_k = 2L, ruv_percentage = 10), seen)
    expect_identical(owner$current_state, "trx_ruvs")
    expect_false(k194Noted(seen, "[error]"), info = paste(seen$notes, collapse = " | "))
    trx_config <- omicGuiConfig("transcriptomics")
    calls <- new.env(parent = emptyenv())
    calls$runs <- 0L
    da <- trx_config$da
    da$preflight <- NULL
    testthat::local_mocked_bindings(executeOmicDaRun = function(...) {
        calls$runs <- calls$runs + 1L
        list(status = "success")
    }, .package = "MultiScholaR")
    kit172Quiet(shiny::isolate(runOmicDaAnalysis(list(da_q_val_thresh = 0.05, treat_lfc_cutoff = 0), NULL, project$workflow,
        new.env(parent = emptyenv()), list(), da, transcriptomicsDeclarationBinding)))
    expect_identical(calls$runs, 1L)
})

# ------------------------------------------------------------------------------------------------
# T10: replay untouched (both commits).
# ------------------------------------------------------------------------------------------------

#' The record of a scripted journey replayed on a copy of the store cut to the design state (as
#' test-workflow-lipid-declarations.R replays one).
k194Replay <- function(fx, root) {
    owner <- k194Owner(fx)
    json <- file.path(root, "analysis_record.json")
    writeWorkflowAnalysisRecord(owner, json)
    record <- jsonlite::read_json(json, simplifyVector = FALSE)
    history <- k194History(fx)
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
    replay <- kit172Quiet(harness$kitReplayJourney(record, prepared, store_root, inputs, fx$paths, scratch, api))
    list(record = record, history = history, replay = replay)
}

for (class in c("lipid", "metab")) {
    test_that(paste("KIT-194 T10:", class, "a record of J1's sequence and one of the default path replay as recorded"), {
        skip_if_not(workflowOwnerCapable())
        root <- withr::local_tempdir()
        for (journey in c("j1", "default")) {
            fx <- k194Fixture(file.path(root, journey, "project"), class)
            k194Defer(fx)
            owner <- k194Owner(fx)
            prefix <- fx$prefix
            config <- fx$workflow$config_list
            kit172Quiet({
                runWorkflowOmicsReferenceLog2(fx$binding, owner, 0, config)
                runWorkflowOmicsReferenceBetweenSampleNormalisation(fx$binding, owner, "quantile", config)
                if (identical(journey, "j1")) selectWorkflowStateReference(owner, paste0(prefix, "_qc_complete"))
                runWorkflowOmicsReferenceCorrelation(fx$binding, owner, "group", 0.5, config)
                runWorkflowOmicsReferenceDaAnalysis(fx$binding, owner, getWorkflowStateReference(owner),
                    list(formula_string = "~ 0 + group", da_q_val_thresh = 0.05, treat_lfc_cutoff = 0), kit172Contrasts())
            })
            out <- k194Replay(fx, file.path(root, journey))
            operations <- vapply(out$record$steps, function(step) as.character(step$operation %||% ""), "")
            info <- paste("T10", class, journey, paste(operations, collapse = " "))
            expect_identical(out$replay$states, out$history, info = info)
            # The store also holds states a Revert left out of the active history (J1's log2 and normalised states); the
            # record and the replay are compared on the active history's states.
            mismatches <- as.character(out$replay$mismatches)
            expect_identical(mismatches[sub("/.*$", "", mismatches) %in% out$history], character(), info = paste(info, paste(mismatches, collapse = " ")))
            expect_true(any(vapply(out$replay$per_step, function(x) identical(x$operation, paste0(prefix, "_da_contrast")), NA)), info = info)
            # J1's record shows no log2 or normalisation step: the DA tested the QC table, as recorded.
            expect_identical(any(paste0(prefix, c("_log2", "_between_sample_normalisation")) %in% operations), identical(journey, "default"), info = info)
        }
    })
}

# ------------------------------------------------------------------------------------------------
# The refusals' classes (7.8).
# ------------------------------------------------------------------------------------------------

test_that("KIT-194: the guards raise their classes", {
    skip_if_not(workflowOwnerCapable())
    root <- withr::local_tempdir()
    seen <- k194Notes()
    run_start <- get0("omicNormRunStart", envir = asNamespace("MultiScholaR"), inherits = FALSE)
    qc_guard <- get0("assertOmicQcHeadNotNormalised", envir = asNamespace("MultiScholaR"), inherits = FALSE)
    result_guard <- get0("assertOmicNormHeadIsResult", envir = asNamespace("MultiScholaR"), inherits = FALSE)
    expect_true(is.function(run_start) && is.function(qc_guard) && is.function(result_guard))
    if (!(is.function(run_start) && is.function(qc_guard) && is.function(result_guard))) return(invisible(NULL))
    fx <- k194Fixture(file.path(root, "a"), "lipid")
    k194Defer(fx)
    nd <- omicNormReactiveState()
    k194Run(fx, nd)
    k194ExpectRun(fx, nd, "classes")
    binding <- fx$config$binding$binding
    expect_error(qc_guard(binding, k194Owner(fx)), class = "omic_qc_input_normalised")
    expect_error(shiny::isolate(run_start(binding, k194Owner(fx), omicNormReactiveState())), class = "omic_norm_input_normalised")
    rec <- k194Runner(fx$workflow, busy = TRUE)
    testthat::local_mocked_bindings(workflowOmicsManagedActions = function(omic_type) rec$runner, .package = "MultiScholaR")
    expect_error(shiny::isolate(run_start(binding, k194Owner(fx), nd)), class = "omic_norm_runner_busy")
    selectWorkflowStateReference(k194Owner(fx), paste0(fx$prefix, "_log2"))
    expect_error(shiny::isolate(result_guard(binding, k194Owner(fx), nd)), class = "omic_norm_correlation_input_stale")
    # No owner: every guard is a no-op.
    expect_silent(qc_guard(binding, "owner"))
    expect_true(isTRUE(shiny::isolate(result_guard(binding, "owner", nd))))
})
