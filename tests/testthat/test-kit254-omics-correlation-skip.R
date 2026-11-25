# KIT-254 (AUDIT-KIT-248-plan P03, U1; AUDIT-KIT-254-plan): lipidomics and metabolomics, "Skip Filtering & Proceed to
# DE" pressed after the correlation filter was applied, and the filter applied again, start from the table the filter
# first received (the run's result). Written before the fix: on the tests-only commit each E row runs as base does and
# fails at the value G0 recorded (docs/qa/workflow-execution/kit/kit-254/g0-base-values.tsv, run 45599722106); the G rows
# pass on both commits. Every row asserts the DESIGN's sample list (R1-07 (a)), never an assay's columns. Lines
# `KIT-254 <row>` are the evidence the ticket's records copy.

library(testthat)

# ------------------------------------------------------------------------------------------------
# Copied from tests/testthat/test-omics-norm-state-guards.R at e0747fc1 (k194PreNorm :22-31, k194Fixture :101-115,
# k194Close, k194Defer, k194Owner, k194History, k194Head, k194Ref :117-134, k194Notes :139-155, k194Input :158-164,
# k194Run :167-171, k194Apply :180-184, k194Skip :186), renamed kit254*; and from tests/testthat/test-kit212-da-current.R
# at e0747fc1 (k207Save :270-288, k207Da :291-295, k212Dirs :337-343, k212Texts$stale :319-321, k212ExpectSaveRefused
# :403-411, k212Tab :647-666, k212Replay :447-475, extended to read the replayed head's design), renamed kit254*.
# ------------------------------------------------------------------------------------------------

kit254Omic <- function(class) c(lipid = "lipidomics", metab = "metabolomics")[[class]]

kit254PreNorm <- function(root, values, class) {
    fixture <- kit172Owner(root, kit172Object(values, class), class)
    binding <- kit172Binding(class)
    kit172Quiet({
        runWorkflowOmicsReferenceIntensityFilter(binding, fixture$owner, 1, 0.5, fixture$workflow$config_list)
        commitWorkflowOmicsStoreStep(fixture$owner, paste0(binding$prefix, "_qc_finalize"), paste0(binding$prefix, "_qc_complete"),
            fixture$workflow$config_list, "QC complete")
    })
    fixture
}

#' A lipid or metab project at QC complete on F172 (seed 172001; `values` replaces it), with the omic's GUI configuration
#' (its QC figures drawn by nothing) and the tab's log.
kit254Fixture <- function(root, class, values = kit172Values(172001L)$values) {
    fixture <- kit254PreNorm(root, values, class)
    fixture$class <- class
    fixture$omic <- kit254Omic(class)
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

kit254Close <- function(fixture) try(closeDIAStateAdapter(fixture$workflow$state_manager), silent = TRUE)
kit254Defer <- function(fixture, envir = parent.frame()) {
    force(fixture)
    withr::defer(kit254Close(fixture), envir = envir)
}
kit254Owner <- function(fixture) fixture$workflow$state_manager
kit254History <- function(fixture) as.character(unlist(kit254Owner(fixture)$getHistory()))
kit254Head <- function(fixture) kit254Owner(fixture)$current_state
kit254Ref <- function(fixture, state = NULL) getWorkflowStateReference(kit254Owner(fixture), state)

kit254Notes <- function(env = parent.frame()) {
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

kit254Input <- function(method = "quantile", itsd = FALSE, ruv = NULL, threshold = 0.5) {
    list(apply_itsd = itsd, itsd_aggregation = "median", log_offset = 0, norm_method = method,
        ruv_mode = if (is.null(ruv)) "skip" else ruv$ruv_mode, auto_percentage_min = 1, auto_percentage_max = 20,
        ruv_grouping_variable = "group", separation_metric = "max_difference", k_penalty_weight = 0.5,
        adaptive_k_penalty = TRUE, ruv_k = ruv$manual_k %||% 1L, ruv_percentage = ruv$manual_percentage %||% 20,
        min_pearson_correlation_threshold = threshold)
}

kit254Run <- function(fixture, nd, input = kit254Input(), runner = NULL) {
    testthat::local_mocked_bindings(workflowOmicsManagedActions = function(omic_type) runner, .package = "MultiScholaR")
    kit172Quiet(shiny::isolate(runOmicNormPipeline(input, fixture$workflow, list(), nd, fixture$addLog,
        function(...) list(color_var = "group", shape_var = NULL), fixture$config)))
}

kit254Apply <- function(fixture, nd, runner = NULL, threshold = 0.5, input = kit254Input(threshold = threshold)) {
    testthat::local_mocked_bindings(workflowOmicsManagedActions = function(omic_type) runner, .package = "MultiScholaR")
    kit172Quiet(shiny::isolate(applyOmicNormCorrelation(input, fixture$workflow, nd, fixture$addLog, fixture$config)))
}

kit254Skip <- function(fixture, nd) kit172Quiet(shiny::isolate(skipOmicNormCorrelation(fixture$workflow, nd, fixture$addLog, fixture$config)))

kit254Save <- function(workflow, omic, dirs, runner = NULL) {
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

kit254Da <- function(fx) {
    kit172Quiet(runWorkflowOmicsReferenceDaAnalysis(fx$binding, kit254Owner(fx), kit254Ref(fx),
        list(formula_string = "~ 0 + group", da_q_val_thresh = 0.05, treat_lfc_cutoff = 0), kit172Contrasts()))
    invisible(fx)
}

kit254Dirs <- function(root, omic) {
    paths <- list(source_dir = file.path(root, "scripts", omic), integration_dir = file.path(root, "integration"), base_dir = root,
        da_output_dir = file.path(root, "da_output"), publication_graphs_dir = file.path(root, "publication_graphs"),
        results_summary_dir = file.path(root, "results_summary"))
    for (path in paths) dir.create(path, recursive = TRUE, showWarnings = FALSE)
    stats::setNames(list(paths), omic)
}

kit254Stale <- paste("The differential analysis was run on an earlier table, not the current one: the table changed",
    "after the analysis ran (a new normalisation run, a Revert, or the correlation filter applied again). Run the",
    "differential analysis again, then press 'Save Workflow Arguments' before copying or generating the report.")

kit254ExpectSaveRefused <- function(rec, dirs, omic, info) {
    source_dir <- dirs[[omic]]$source_dir
    info <- paste(info, "| notes:", paste(rec$notes, collapse = " | "))
    expect_true(paste("[error] The session could not be saved:", kit254Stale) %in% rec$notes, info = info)
    expect_identical(rec$done, 0L, info = info)
    expect_true(length(rec$failed) == 1L && inherits(rec$failed[[1L]], "workflow_summary_da_stale"), info = info)
    expect_false(file.exists(file.path(source_dir, "study_parameters.txt")), info = info)
    expect_false(file.exists(file.path(source_dir, "analysis_record.json")), info = info)
}

#' The DA tab of lipid on a `shiny::reactiveValues` bus over the fixture's owner; `code(tab, run_da, fire, bus, session)`.
kit254Tab <- function(fx, code) {
    bus <- shiny::reactiveValues(state_manager = kit254Owner(fx), contrasts_tbl = kit172Contrasts(), design_matrix = kit172Design(),
        config_list = fx$workflow$config_list, tab_status = list(normalization = "complete", differential_analysis = "pending"),
        state_update_trigger = NULL)
    testthat::local_mocked_bindings(workflowOmicsManagedActions = function(omic_type) NULL, .package = "MultiScholaR")
    kit172Quiet(utils::capture.output(shiny::testServer(mod_lipid_da_server, args = list(workflow_data = bus, experiment_paths = list(),
        omic_type = "lipidomics", experiment_label = "k254"), {
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

#' The owner's record (as Save writes it) replayed on a copy of the store cut to the design state; the replay's error
#' or value, and the replayed head with its design's samples (read from the replay store after the replay closed it).
kit254Replay <- function(fx, root) {
    owner <- kit254Owner(fx)
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
    # Each assay input matched to the record's role by its SHA-256 (the roles' order is the record's, not the assays').
    sources <- list.files(file.path(fx$source_dir, "import_sources"), full.names = TRUE)
    digests <- vapply(sources, function(path) digest::digest(file = path, algo = "sha256"), character(1))
    roles <- Filter(function(input) startsWith(input$role, fx$binding$design_input_prefix), record$inputs)
    inputs <- c(stats::setNames(unname(sources[match(vapply(roles, function(input) input$sha256, ""), digests)]),
        vapply(roles, function(input) input$role, "")), design_file = file.path(fx$source_dir, "design_matrix.tab"))
    scratch <- file.path(root, "scratch")
    dir.create(scratch)
    value <- NULL
    warned <- character()
    warned_design <- character()
    # KIT-254 (OD-1 (a), VAL-247): the replay's statement is also raised as a warning of its class; it is kept here.
    error <- tryCatch({
        value <- kit172Quiet(withCallingHandlers(replayWorkflowAnalysisRecord(record, prepared, inputs, fx$paths, scratch),
            workflow_replay_correlation_skip = function(w) {
                warned <<- c(warned, conditionMessage(w))
                invokeRestart("muffleWarning")
            }, workflow_replay_replicate_key_collision = function(w) {
                warned_design <<- c(warned_design, conditionMessage(w))
                invokeRestart("muffleWarning")
            }))
        NULL
    }, error = function(e) e)
    head <- NULL
    runs <- NULL
    if (is.null(error)) {
        identity <- readSavedStoreIdentity(file.path(prepared, ".workflow-entities", "store.json"))
        manager <- newDIAStateAdapter(newWorkflowEntityStore(prepared, identity$project_id, identity$workflow_id), reopen = TRUE)
        head <- manager$current_state
        runs <- kit254RunsOf(fx$binding, manager, NULL)
        closeDIAStateAdapter(manager)
    }
    cat(sprintf("KIT-254 replay %s\tsteps=%s\terror=%s\thead=%s\tremoved=%s\n", basename(root),
        paste(vapply(record$steps, function(step) as.character(step$state), ""), collapse = ","),
        if (is.null(error)) "none" else paste(class(error)[[1L]], conditionMessage(error)), format(head %||% "<none>"),
        paste(kit254Removed(runs %||% character()), collapse = ",")))
    list(record = record, error = error, value = value, head = head, runs = runs, warned = warned, warned_design = warned_design)
}

# ------------------------------------------------------------------------------------------------
# KIT-254's harness.
# ------------------------------------------------------------------------------------------------

kit254RunsOf <- function(binding, manager, reference) {
    meta <- readWorkflowOmicsMetadata(binding, manager, reference)
    if (is.null(meta)) return(character())
    sort(as.character(as.data.frame(meta$design_matrix)[["Run"]]))
}
#' The design's samples of a state (the head when NULL).
kit254Runs <- function(fx, reference = NULL) {
    if (is.character(reference)) reference <- kit254Ref(fx, reference)
    kit254RunsOf(fx$binding, kit254Owner(fx), reference)
}
kit254All <- sort(kit172Design()$Run)
kit254Removed <- function(runs) setdiff(kit254All, runs)
kit254Final <- function(nd) shiny::isolate(nd$ruv_corrected_obj %||% nd$post_norm_obj)
kit254Components <- function(fx, reference) readWorkflowStateMetadata(kit254Owner(fx), reference)$data$components
kit254Steps <- function(record) vapply(record$steps, function(step) as.character(step$operation), "")
kit254Fields <- function(nd) shiny::isolate(list(complete = nd$correlation_filtering_complete, obj = nd$correlation_filtered_obj,
    results = nd$correlation_results))
kit254Line <- function(row, label, ...) cat(sprintf("KIT-254 %s %s\t%s\n", row, label, paste(..., sep = "\t")))

#' The DA's source (a saved reference) as the owner recorded it.
kit254DaSource <- function(fx) {
    owner <- kit254Owner(fx)
    withWorkflowStateScope(owner, function() {
        da <- owner$.context$snapshot$da_results
        if (!isSavedWorkflowReference(da)) return(NULL)
        resolveWorkflowStateValue(owner$.context, da)$source
    })
}

#' The source of every in-process correlation: the head's entity id when the wrapper is called, in `probe$sources`.
kit254SourceProbe <- function(probe = new.env(parent = emptyenv()), env = parent.frame()) {
    wrapper <- get("runWorkflowOmicsReferenceCorrelation", envir = asNamespace("MultiScholaR"))
    probe$sources <- character()
    testthat::local_mocked_bindings(runWorkflowOmicsReferenceCorrelation = function(binding, manager, ...) {
        probe$sources <- c(probe$sources, getWorkflowStateReference(manager)$entity_id)
        wrapper(binding, manager, ...)
    }, .package = "MultiScholaR", .env = env)
    probe
}

#' A fixture after "Run Normalization Pipeline" in process (`ruv` the legacy list, or NULL to skip RUV); closed when the
#' calling test ends. The run is the row's precondition (an instrument check, VAL-54).
kit254Ran <- function(root, class, ruv = NULL, env = parent.frame()) {
    fx <- kit254Fixture(root, class)
    kit254Defer(fx, envir = env)
    nd <- omicNormReactiveState()
    kit254Run(fx, nd, kit254Input(ruv = ruv))
    expect_true(kit254Head(fx) %in% paste0(fx$prefix, c("_normalized", "_ruv_corrected")),
        info = paste("run precondition:", class, basename(root), kit254Head(fx), paste(fx$log$lines, collapse = " | ")))
    list(fx = fx, nd = nd)
}

#' A press's precondition (an instrument check, VAL-54): the head the press must reach on both commits.
kit254ExpectHead <- function(fx, state, info) expect_identical(kit254Head(fx), paste0(fx$prefix, state),
    info = paste("head precondition:", info, paste(fx$log$lines, collapse = " | ")))

# ---- the managed harness: KIT-194 F1's worker journey (omicsManagedChainMocks, the real runner) with the tab's run,
# Apply and Skip observers (as omicsTabServer, helper-workflow-omics-worker.R:193-238, mounts them, the Skip added) ----

kit254Server <- function(fx, probe) {
    force(fx); force(probe)
    function(input, output, session) {
        workflow <- fx$workflow
        runner <- installWorkflowOmicsActions(session, workflow, fx$omic, "kit254", fx$paths)
        runner$progress <- function(message) list(set = function(...) invisible(NULL), inc = function(...) invisible(NULL),
            close = function(...) invisible(NULL))
        probe$runner <- runner
        nd <- omicNormReactiveState()
        nd$assay_names <- "LCMS_Pos"
        probe$nd <- nd
        aesthetics <- function(...) list(color_var = "group", shape_var = NULL)
        config <- omicGuiConfig(fx$omic)
        shiny::observeEvent(input$run_normalization, runOmicNormPipeline(input, workflow, fx$paths, nd, fx$addLog, aesthetics, config))
        shiny::observeEvent(input$apply_correlation_filter, applyOmicNormCorrelation(input, workflow, nd, fx$addLog, config))
        shiny::observeEvent(input$skip_correlation_filter, skipOmicNormCorrelation(workflow, nd, fx$addLog, config))
    }
}

kit254Settle <- function(session, probe, rounds = 40L) {
    for (i in seq_len(rounds)) {
        if (is.null(shiny::isolate(probe$runner$step$active()))) break
        session$elapse(100)
    }
    is.null(shiny::isolate(probe$runner$step$active()))
}

#' One managed sequence on a fresh F172 owner at QC complete (closed when the calling test ends): inside the server the
#' tab's inputs are set and the run pressed, then `steps(session, press, probe, fx)`; `press(button, extra)` sets
#' `extra`, then the button, and flushes. The worker's input head per action is in `probe$inputs` (entity ids, with the
#' operation as name); `probe$hold` keeps the worker's poll at "running". Returns `list(fx, nd, probe)`, `nd` the tab's
#' values as a list read before the session ended.
kit254Managed <- function(root, class, ruv, steps, env = parent.frame()) {
    fx <- kit254Fixture(root, class)
    kit254Defer(fx, envir = env)
    probe <- new.env()
    probe$hold <- FALSE
    probe$inputs <- character()
    kit254ManagedRun(fx, probe, ruv, steps)
    list(fx = fx, nd = probe$values, probe = probe)
}

kit254ManagedRun <- function(fx, probe, ruv, steps) {
    omicsManagedChainMocks(fx, probe)
    dispatch <- get("runWorkflowAction", envir = asNamespace("MultiScholaR"))
    testthat::local_mocked_bindings(
        runWorkflowAction = function(operation, input_references, parameters, owner, reservation, completion_context = NULL) {
            probe$inputs <- c(probe$inputs, stats::setNames(getWorkflowStateReference(fx$owner)$entity_id, operation))
            dispatch(operation, input_references, parameters, owner, reservation, completion_context)
        },
        pollWorkflowAction = function(owner) if (isTRUE(probe$hold)) list(status = "running") else list(status = "complete"),
        .package = "MultiScholaR")
    suppressWarnings(suppressMessages(utils::capture.output(shiny::testServer(kit254Server(fx, probe), {
        session$userData$probe_owner <- fx$owner
        presses <- new.env(parent = emptyenv())
        press <- function(button, extra = list()) {
            presses[[button]] <- (presses[[button]] %||% 0L) + 1L
            if (length(extra)) do.call(session$setInputs, extra)
            do.call(session$setInputs, stats::setNames(list(presses[[button]]), button))
            session$flushReact()
        }
        do.call(session$setInputs, kit254Input(ruv = ruv))
        press("run_normalization")
        probe$settled_run <- kit254Settle(session, probe)
        steps(session, press, probe, fx)
        # The tab's values while the session lives (a module's reactive values cannot be read after it ends).
        probe$values <- shiny::isolate(shiny::reactiveValuesToList(probe$nd))
    }))))
    expect_true(isTRUE(probe$settled_run) && kit254Head(fx) %in% paste0(fx$prefix, c("_normalized", "_ruv_corrected", "_correlation_filtered",
        "_norm_complete")), info = paste("managed run precondition:", kit254Head(fx), paste(fx$log$lines, collapse = " | ")))
    invisible(probe)
}

#' The worker's input heads of the managed pairs actions, in order.
kit254PairsInputs <- function(fx, probe) unname(probe$inputs[names(probe$inputs) == paste0(fx$prefix, "_pearson_pairs")])

# ------------------------------------------------------------------------------------------------
# The configurations and G0's thresholds (run 45599722106, g0-base-values.tsv; the rule of ticket section 10).
# ------------------------------------------------------------------------------------------------

kit254Configs <- list(
    lipid = list(class = "lipid", ruv = NULL),
    metab = list(class = "metab", ruv = NULL),
    lipid_ruv = list(class = "lipid", ruv = kit172LegacyRuv("manual")))

# Lipid with RUV (AUDIT-KIT-254-impl I01): after KIT-224 (staging c4505c7f) G0's rule finds no listed threshold that
# removes some samples and keeps two per group (G0b, run 45628705535: none up to 0.95, ten at 0.99), so stop rule 2
# applies: t_hi 0.5, t_lo 0.3, nothing removed; E7 keeps its state and record assertions, its sample difference is not
# shown on this fixture (was 0.99 / 0.95 at e0747fc1, G0 45599722106: B_4).
kit254T <- list(lipid = c(hi = 0.95, lo = 0.9), metab = c(hi = 0.95, lo = 0.9), lipid_ruv = c(hi = 0.5, lo = 0.3))
# G0's design samples (OBSERVED, both omics in process; C3 and II0): Apply(t_hi) alone removes A_5 (B_4 with RUV);
# Apply(t_lo) alone removes none.
kit254HiRuns <- list(lipid = setdiff(kit254All, "A_5"), metab = setdiff(kit254All, "A_5"), lipid_ruv = kit254All)
kit254LoRuns <- list(lipid = kit254All, metab = kit254All, lipid_ruv = kit254All)

# ------------------------------------------------------------------------------------------------
# E1 to E4, E7: Skip after Apply (failing first).
# ------------------------------------------------------------------------------------------------

#' E1 to E4 on a sequence that ended with Apply(t_hi), Skip: the final state, the record and replay, the DA, the Save.
kit254ExpectSkipAfterApply <- function(fx, nd, label, root) {
    prefix <- fx$prefix
    kit254ExpectHead(fx, "_norm_complete", label)
    final <- kit254Final(nd)
    expect_true(isSavedWorkflowReference(final), info = paste("final precondition:", label))
    runs <- kit254Runs(fx)
    same <- identical(kit254Components(fx, kit254Ref(fx)), kit254Components(fx, final))
    history <- kit254History(fx)
    kit254Line("E1", label, paste("history", paste(history, collapse = " ")), paste("removed", paste(kit254Removed(runs), collapse = ",")),
        paste("components_equal", same))
    # E1: the final state is the run's result, every sample.
    expect_false(paste0(prefix, "_correlation_filtered") %in% history, info = paste("E1", label))
    expect_true(same, info = paste("E1 components", label))
    expect_identical(runs, kit254All, info = paste("E1 design", label))
    # E4: the record has no selection; the replay's head has every sample.
    out <- kit254Replay(fx, file.path(root, "replay"))
    expect_null(out$error, info = paste("E4 replay precondition", label))
    expect_false(paste0(prefix, "_pearson_selection") %in% kit254Steps(out$record), info = paste("E4 record", label))
    expect_identical(out$runs, kit254All, info = paste("E4 replay design", label))
    # G10 (OD-1 (a), VAL-247; both commits of the follow-up): a record made after the fix carries no statement.
    expect_null(out$value$correlation_statement)
    expect_length(out$warned, 0L)
    # E2: the DA's source has every sample.
    kit254Da(fx)
    source <- kit254DaSource(fx)
    expect_true(isSavedWorkflowReference(source), info = paste("E2 DA precondition", label))
    da_runs <- kit254RunsOf(fx$binding, kit254Owner(fx), source)
    kit254Line("E2", label, paste("da_removed", paste(kit254Removed(da_runs), collapse = ",")))
    expect_identical(da_runs, kit254All, info = paste("E2", label))
    # E3: the Save's integration object and study parameters describe every sample, no threshold.
    dirs <- kit254Dirs(file.path(root, "save"), fx$omic)
    rec <- kit254Save(fx$workflow, fx$omic, dirs)
    expect_identical(rec$done, 1L, info = paste("E3 Save precondition", label, paste(rec$notes, collapse = " | ")))
    objects <- Filter(function(call) isS4(call$object) && "design_matrix" %in% methods::slotNames(call$object), rec$rds)
    counts <- vapply(objects, function(call) nrow(methods::slot(call$object, "design_matrix")), integer(1))
    params <- file.path(dirs[[fx$omic]]$source_dir, "study_parameters.txt")
    lines <- if (file.exists(params)) readLines(params, warn = FALSE) else character()
    threshold <- grep("min_pearson_correlation_threshold", lines, value = TRUE, fixed = TRUE)
    kit254Line("E3", label, paste("objects", paste(counts, collapse = ",")), paste("threshold_lines", paste(trimws(threshold), collapse = " | ")))
    expect_true(length(counts) > 0L && file.exists(params), info = paste("E3 precondition: the integration object and parameters written", label))
    expect_true(all(counts == length(kit254All)), info = paste("E3 integration object", label))
    expect_length(threshold, 0L)
}

for (name in names(kit254Configs)) {
    config <- kit254Configs[[name]]
    rows <- if (is.null(config$ruv)) "E1 to E4" else "E7 (E1 to E4 with RUV)"
    test_that(paste("KIT-254", rows, name, "in process: a Skip after an Apply commits the run's result"), {
        skip_if_not(workflowOwnerCapable())
        root <- withr::local_tempdir()
        seen <- kit254Notes()
        q <- kit254Ran(file.path(root, "a"), config$class, config$ruv)
        expect_true(isTRUE(kit254Apply(q$fx, q$nd, threshold = kit254T[[name]][["hi"]])), info = paste("Apply precondition", name))
        kit254ExpectHead(q$fx, "_correlation_filtered", paste(name, "Apply"))
        # AUDIT-KIT-254-impl I01: the first Apply removed the sample G0 recorded (a moved threshold fails here, loudly).
        expect_identical(kit254Runs(q$fx), kit254HiRuns[[name]], info = paste("Apply(t_hi) precondition", name))
        expect_true(isTRUE(kit254Skip(q$fx, q$nd)), info = paste("Skip precondition", name, paste(seen$notes, collapse = " | ")))
        kit254ExpectSkipAfterApply(q$fx, q$nd, paste(name, "in_process"), root)
    })
    # The managed route with RUV is not run: in G0 the worker-journey harness dispatched no RUV action for these inputs
    # (run 45599722106, case I lipid_ruv managed; outcome.md "Not checked").
    if (is.null(config$ruv)) test_that(paste("KIT-254", rows, name, "managed: a Skip after an Apply commits the run's result"), {
        skip_if_not(workflowOwnerCapable())
        root <- withr::local_tempdir()
        seen <- kit254Notes()
        m <- kit254Managed(file.path(root, "m"), config$class, config$ruv, function(session, press, probe, fx) {
            press("apply_correlation_filter", list(min_pearson_correlation_threshold = kit254T[[name]][["hi"]]))
            probe$settled <- kit254Settle(session, probe)
            probe$applied <- fx$owner$current_state
            press("skip_correlation_filter")
        })
        expect_true(isTRUE(m$probe$settled) && identical(m$probe$applied, paste0(m$fx$prefix, "_correlation_filtered")),
            info = paste("managed Apply precondition", name, paste(m$fx$log$lines, collapse = " | ")))
        kit254ExpectSkipAfterApply(m$fx, m$nd, paste(name, "managed"), root)
    })
}

# ------------------------------------------------------------------------------------------------
# E5, E6, E7, E12: Apply again (failing first).
# ------------------------------------------------------------------------------------------------

#' E5, E6 on a sequence that ended with Apply(t_hi), Apply(t_lo): `source` the second correlation's input entity id.
kit254ExpectReapplied <- function(fx, nd, source, name, label, root) {
    kit254ExpectHead(fx, "_correlation_filtered", label)
    final <- kit254Final(nd)
    runs <- kit254Runs(fx)
    kit254Line("E5", label, paste("source_is_result", identical(source, final$entity_id)), paste("removed", paste(kit254Removed(runs), collapse = ",")))
    expect_identical(source, final$entity_id, info = paste("E5 source", label))
    expect_identical(runs, kit254LoRuns[[name]], info = paste("E5 design", label))
    out <- kit254Replay(fx, file.path(root, "replay"))
    expect_null(out$error, info = paste("E6 replay precondition", label))
    expect_identical(out$runs, runs, info = paste("E6 replay equals live", label))
}

for (name in names(kit254Configs)) {
    config <- kit254Configs[[name]]
    rows <- if (is.null(config$ruv)) "E5, E6" else "E7 (E5, E6 with RUV)"
    test_that(paste("KIT-254", rows, name, "in process: a re-applied filter reads the run's result"), {
        skip_if_not(workflowOwnerCapable())
        root <- withr::local_tempdir()
        seen <- kit254Notes()
        probe <- kit254SourceProbe()
        q <- kit254Ran(file.path(root, "a"), config$class, config$ruv)
        expect_true(isTRUE(kit254Apply(q$fx, q$nd, threshold = kit254T[[name]][["hi"]])), info = paste("Apply precondition", name))
        # AUDIT-KIT-254-impl I01: the first Apply removed the sample G0 recorded.
        expect_identical(kit254Runs(q$fx), kit254HiRuns[[name]], info = paste("Apply(t_hi) precondition", name))
        expect_true(isTRUE(kit254Apply(q$fx, q$nd, threshold = kit254T[[name]][["lo"]])), info = paste("Apply again precondition", name))
        expect_length(probe$sources, 2L)
        kit254ExpectReapplied(q$fx, q$nd, probe$sources[[2L]], name, paste(name, "in_process"), root)
    })
    if (is.null(config$ruv)) test_that(paste("KIT-254", rows, name, "managed: a re-applied filter reads the run's result"), {
        skip_if_not(workflowOwnerCapable())
        root <- withr::local_tempdir()
        seen <- kit254Notes()
        m <- kit254Managed(file.path(root, "m"), config$class, config$ruv, function(session, press, probe, fx) {
            press("apply_correlation_filter", list(min_pearson_correlation_threshold = kit254T[[name]][["hi"]]))
            probe$settled <- kit254Settle(session, probe)
            press("apply_correlation_filter", list(min_pearson_correlation_threshold = kit254T[[name]][["lo"]]))
            probe$settled_again <- kit254Settle(session, probe)
        })
        pairs <- kit254PairsInputs(m$fx, m$probe)
        expect_true(isTRUE(m$probe$settled) && isTRUE(m$probe$settled_again) && length(pairs) == 2L,
            info = paste("managed precondition", name, paste(names(m$probe$inputs), collapse = " ")))
        kit254ExpectReapplied(m$fx, m$nd, pairs[[2L]], name, paste(name, "managed"), root)
    })
}

for (name in c("lipid", "metab")) {
    test_that(paste("KIT-254 E12:", name, "managed: an Apply queued behind a running Apply reads the run's result"), {
        skip_if_not(workflowOwnerCapable())
        root <- withr::local_tempdir()
        seen <- kit254Notes()
        config <- kit254Configs[[name]]
        m <- kit254Managed(file.path(root, "m"), config$class, NULL, function(session, press, probe, fx) {
            probe$hold <- TRUE
            press("apply_correlation_filter", list(min_pearson_correlation_threshold = kit254T[[name]][["hi"]]))
            probe$held <- !is.null(shiny::isolate(probe$runner$step$active()))
            press("apply_correlation_filter", list(min_pearson_correlation_threshold = kit254T[[name]][["lo"]]))
            probe$queued <- length(probe$inputs)
            probe$hold <- FALSE
            probe$settled <- kit254Settle(session, probe)
        })
        pairs <- kit254PairsInputs(m$fx, m$probe)
        expect_true(isTRUE(m$probe$held) && isTRUE(m$probe$settled) && length(pairs) == 2L,
            info = paste("E12 precondition: the first press held, the second queued", name, paste(names(m$probe$inputs), collapse = " ")))
        kit254ExpectHead(m$fx, "_correlation_filtered", paste("E12", name))
        final <- kit254Final(m$nd)
        runs <- kit254Runs(m$fx)
        kit254Line("E12", name, paste("source_is_result", identical(pairs[[2L]], final$entity_id)), paste("removed", paste(kit254Removed(runs), collapse = ",")))
        expect_identical(pairs[[2L]], final$entity_id, info = paste("E12 source", name))
        expect_identical(runs, kit254LoRuns[[name]], info = paste("E12 design", name))
    })
}

# ------------------------------------------------------------------------------------------------
# E9: two assays, the cascade (failing first).
# ------------------------------------------------------------------------------------------------

#' Ticket 2.3's example on F172's group A (a = A_1, b = A_2): every within-group pair of `runs` (the input design's
#' samples). LCMS_Pos: every pair with A_1 at 0.1, every other at 0.99. LCMS_Neg: A_1-A_2 and the pairs among A_3 to
#' A_6 at 0.99, the other A pairs at 0.1; group B at 0.99.
kit254E9Pairs <- function(runs, assay) {
    groups <- split(sort(runs), sub("_.*$", "", sort(runs)))
    rows <- do.call(rbind, lapply(groups, function(g) {
        if (length(g) < 2L) return(NULL)
        pairs <- utils::combn(g, 2L)
        data.frame(Run.x = pairs[1L, ], Run.y = pairs[2L, ], stringsAsFactors = FALSE)
    }))
    key <- paste(rows$Run.x, rows$Run.y)
    pass <- if (identical(assay, "LCMS_Pos")) rows$Run.x != "A_1" & rows$Run.y != "A_1" else
        !startsWith(rows$Run.x, "A_") | key %in% c("A_1 A_2", "A_3 A_4", "A_3 A_5", "A_3 A_6", "A_4 A_5", "A_4 A_6", "A_5 A_6")
    rows$pearson_correlation <- ifelse(pass, 0.99, 0.1)
    rownames(rows) <- NULL
    rows
}

for (class in c("lipid", "metab")) {
    test_that(paste("KIT-254 E9:", class, "two assays: applying the filter again at the same threshold removes no further sample"), {
        skip_if_not(workflowOwnerCapable())
        root <- withr::local_tempdir()
        seen <- kit254Notes()
        values <- list(LCMS_Pos = kit172Values(172001L)$values, LCMS_Neg = kit172Values(172002L)$values)
        fx <- kit254Fixture(file.path(root, "a"), class, values)
        kit254Defer(fx)
        nd <- omicNormReactiveState()
        kit254Run(fx, nd)
        kit254ExpectHead(fx, "_normalized", paste("E9 run", class))
        # As MCI-021.4 binds the kernel (test-module-ci-lipid-norm.R:222-224): the pairs of the input design's samples.
        testthat::local_mocked_bindings(pearsonCorForSamplePairs = function(theObject, ...) {
            runs <- as.character(methods::slot(theObject, "design_matrix")$Run)
            assays <- names(methods::slot(theObject, fx$binding$data_slot))
            stats::setNames(lapply(assays, function(assay) kit254E9Pairs(runs, assay)), assays)
        }, .package = "MultiScholaR")
        expect_true(isTRUE(kit254Apply(fx, nd, threshold = 0.5)), info = paste("E9 first Apply precondition", class))
        expect_identical(kit254Removed(kit254Runs(fx)), "A_1", info = paste("E9 first Apply precondition", class))
        second <- kit254Apply(fx, nd, threshold = 0.5)
        runs <- kit254Runs(fx)
        kit254Line("E9", class, paste("head", kit254Head(fx)), paste("removed", paste(kit254Removed(runs), collapse = ",")),
            paste("log", paste(grep("ERROR", fx$log$lines, value = TRUE), collapse = " | ")))
        kit254ExpectHead(fx, "_correlation_filtered", paste("E9", class))
        expect_true(isTRUE(second), info = paste("E9 second Apply completes", class, paste(seen$notes, collapse = " | ")))
        expect_identical(kit254Removed(runs), "A_1", info = paste("E9", class))
    })
}

# ------------------------------------------------------------------------------------------------
# E10, E11: a failure after the cut (failing first).
# ------------------------------------------------------------------------------------------------

#' The tab and the head after a failure that followed the cut: the run's result, the three fields at their start.
kit254ExpectCleared <- function(fx, nd, row, label) {
    final <- kit254Final(nd)
    fields <- kit254Fields(nd)
    kit254Line(row, label, paste("head", kit254Head(fx)), paste("complete", format(fields$complete)),
        paste("obj", format(!is.null(fields$obj))), paste("results", length(fields$results)))
    expect_identical(kit254Ref(fx)$entity_id, final$entity_id, info = paste(row, label, "head"))
    expect_false(paste0(fx$prefix, "_correlation_filtered") %in% kit254History(fx), info = paste(row, label, "history"))
    expect_identical(fields$complete, FALSE, info = paste(row, label, "complete"))
    expect_null(fields$obj, info = paste(row, label, "filtered object"))
    expect_length(fields$results, 0L)
}

for (class in c("lipid", "metab")) {
    test_that(paste("KIT-254 E10, E11:", class, "a failure after the return to the run's result leaves the tab on that table"), {
        skip_if_not(workflowOwnerCapable())
        root <- withr::local_tempdir()
        seen <- kit254Notes()
        # E10: an Apply whose grouping column is not in the design (an error in the pairs, after the cut).
        q <- kit254Ran(file.path(root, "e10"), class)
        expect_true(isTRUE(kit254Apply(q$fx, q$nd, threshold = kit254T[[class]][["hi"]])), info = paste("E10 Apply precondition", class))
        input <- kit254Input(threshold = kit254T[[class]][["lo"]])
        input$ruv_grouping_variable <- "kit254_absent"
        seen$notes <- character()
        expect_false(isTRUE(kit254Apply(q$fx, q$nd, input = input)), info = paste("E10 failing Apply precondition", class))
        expect_true(any(startsWith(seen$notes, "[error] Error:")), info = paste("E10 notice precondition", paste(seen$notes, collapse = " | ")))
        kit254ExpectCleared(q$fx, q$nd, "E10", class)
        # E11: a Skip whose store commit is refused (F105-02's mock, test-omic-norm-module.R:516-519).
        q <- kit254Ran(file.path(root, "e11"), class)
        expect_true(isTRUE(kit254Apply(q$fx, q$nd, threshold = kit254T[[class]][["hi"]])), info = paste("E11 Apply precondition", class))
        local({
            testthat::local_mocked_bindings(commitWorkflowReferenceState = function(...) stop("kit254: the store commit is refused"),
                .package = "MultiScholaR")
            expect_false(isTRUE(kit254Skip(q$fx, q$nd)), info = paste("E11 refused Skip precondition", class))
        })
        kit254ExpectCleared(q$fx, q$nd, "E11", class)
    })
}

# ------------------------------------------------------------------------------------------------
# E13, E14 (failing first).
# ------------------------------------------------------------------------------------------------

for (class in c("lipid", "metab")) {
    test_that(paste("KIT-254 E13, E14:", class, "Apply, Skip, Apply; and a DA between Apply and Skip becomes stale"), {
        skip_if_not(workflowOwnerCapable())
        root <- withr::local_tempdir()
        seen <- kit254Notes()
        # E13: run, Apply(t_hi), Skip, Apply(t_lo).
        q <- kit254Ran(file.path(root, "e13"), class)
        expect_true(isTRUE(kit254Apply(q$fx, q$nd, threshold = kit254T[[class]][["hi"]])), info = paste("E13 Apply precondition", class))
        expect_true(isTRUE(kit254Skip(q$fx, q$nd)), info = paste("E13 Skip precondition", class))
        expect_true(isTRUE(kit254Apply(q$fx, q$nd, threshold = kit254T[[class]][["lo"]])), info = paste("E13 Apply again precondition", class))
        kit254ExpectHead(q$fx, "_correlation_filtered", paste("E13", class))
        runs <- kit254Runs(q$fx)
        history <- kit254History(q$fx)
        out <- kit254Replay(q$fx, file.path(root, "e13-replay"))
        kit254Line("E13", class, paste("history", paste(history, collapse = " ")), paste("removed", paste(kit254Removed(runs), collapse = ",")),
            paste("replay_removed", paste(kit254Removed(out$runs %||% character()), collapse = ",")))
        expect_false(paste0(q$fx$prefix, "_norm_complete") %in% history, info = paste("E13 history", class))
        expect_identical(runs, kit254LoRuns[[class]], info = paste("E13 design", class))
        expect_null(out$error, info = paste("E13 replay precondition", class))
        expect_identical(out$runs, runs, info = paste("E13 replay equals live", class))
        # E14: run, Apply(t_hi), DA, Skip, Save: refused as stale.
        q <- kit254Ran(file.path(root, "e14"), class)
        expect_true(isTRUE(kit254Apply(q$fx, q$nd, threshold = kit254T[[class]][["hi"]])), info = paste("E14 Apply precondition", class))
        kit254Da(q$fx)
        expect_true(isTRUE(kit254Skip(q$fx, q$nd)), info = paste("E14 Skip precondition", class))
        status <- workflowDaResultStatus(kit254Owner(q$fx))
        dirs <- kit254Dirs(file.path(root, "e14-save"), kit254Omic(class))
        rec <- kit254Save(q$fx$workflow, kit254Omic(class), dirs)
        kit254Line("E14", class, paste("status", status), paste("done", rec$done), paste("notes", paste(rec$notes, collapse = " | ")))
        expect_identical(status, "stale", info = paste("E14 DA status", class))
        kit254ExpectSaveRefused(rec, dirs, kit254Omic(class), paste("E14", class))
    })
}

test_that("KIT-254 E14: lipid, the DA tab's results are cleared when a Skip after its DA returns to the run's result", {
    skip_if_not(workflowOwnerCapable())
    root <- withr::local_tempdir()
    seen <- kit254Notes()
    q <- kit254Ran(file.path(root, "tab"), "lipid")
    expect_true(isTRUE(kit254Apply(q$fx, q$nd, threshold = kit254T$lipid[["hi"]])), info = "E14 tab Apply precondition")
    out <- new.env(parent = emptyenv())
    kit254Tab(q$fx, function(tab, run_da, fire, bus, session) {
        run_da()
        out$after_da <- !is.null(shiny::isolate(tab$da_results_list))
        kit254Skip(q$fx, q$nd)
        fire()
        out$results <- !is.null(shiny::isolate(tab$da_results_list))
        out$complete <- shiny::isolate(tab$analysis_complete)
    })
    kit254Line("E14-tab", "lipid", paste("after_da", out$after_da), paste("results", out$results), paste("complete", format(out$complete)))
    expect_true(isTRUE(out$after_da), info = "E14 tab precondition: the DA's results shown")
    expect_false(isTRUE(out$results), info = "E14 tab results")
    expect_identical(out$complete, FALSE, info = "E14 tab analysis_complete")
})

# ------------------------------------------------------------------------------------------------
# G1, G2, G3, E8 (both commits).
# ------------------------------------------------------------------------------------------------

for (class in c("lipid", "metab")) {
    test_that(paste("KIT-254 G1, G2, G3, E8:", class, "the first correlation action after a run, a higher threshold again, a scripted record"), {
        skip_if_not(workflowOwnerCapable())
        root <- withr::local_tempdir()
        seen <- kit254Notes()
        hi <- kit254T[[class]][["hi"]]
        lo <- kit254T[[class]][["lo"]]
        # G1: run, Skip.
        q <- kit254Ran(file.path(root, "g1"), class)
        expect_true(isTRUE(kit254Skip(q$fx, q$nd)))
        kit254ExpectHead(q$fx, "_norm_complete", paste("G1", class))
        expect_identical(kit254Components(q$fx, kit254Ref(q$fx)), kit254Components(q$fx, kit254Final(q$nd)), info = paste("G1", class))
        out <- kit254Replay(q$fx, file.path(root, "g1-replay"))
        expect_null(out$error, info = paste("G1 replay", class))
        expect_false(paste0(q$fx$prefix, "_pearson_selection") %in% kit254Steps(out$record), info = paste("G1 record", class))
        expect_identical(out$runs, kit254All, info = paste("G1 replay design", class))
        expect_null(out$value$correlation_statement)
        # G2: run, Apply(t_hi).
        q <- kit254Ran(file.path(root, "g2"), class)
        expect_true(isTRUE(kit254Apply(q$fx, q$nd, threshold = hi)))
        kit254ExpectHead(q$fx, "_correlation_filtered", paste("G2", class))
        expect_identical(kit254Runs(q$fx), kit254HiRuns[[class]], info = paste("G2 design", class))
        out <- kit254Replay(q$fx, file.path(root, "g2-replay"))
        expect_null(out$error, info = paste("G2 replay", class))
        expect_identical(sum(kit254Steps(out$record) == paste0(q$fx$prefix, "_pearson_selection")), 1L, info = paste("G2 record", class))
        expect_identical(out$runs, kit254HiRuns[[class]], info = paste("G2 replay design", class))
        # G9 (OD-1 (a); both commits of the follow-up): a record without the pattern (Apply only) carries no statement.
        expect_null(out$value$correlation_statement)
        expect_length(out$warned, 0L)
        # G3: run, Apply(t_lo), Apply(t_hi) = Apply(t_hi) alone.
        q <- kit254Ran(file.path(root, "g3"), class)
        expect_true(isTRUE(kit254Apply(q$fx, q$nd, threshold = lo)))
        expect_true(isTRUE(kit254Apply(q$fx, q$nd, threshold = hi)))
        expect_identical(kit254Runs(q$fx), kit254HiRuns[[class]], info = paste("G3 design", class))
        # E8 (a G row under OD-1 (a0)): a record of base's order made by the wrappers, the store step on the filtered
        # head (as KIT-212 G5), replays as recorded, with no statement.
        q <- kit254Ran(file.path(root, "e8"), class)
        owner <- kit254Owner(q$fx)
        config <- q$fx$workflow$config_list
        kit172Quiet(runWorkflowOmicsReferenceCorrelation(q$fx$binding, owner, "group", hi, config))
        kit172Quiet(commitWorkflowOmicsStoreStep(owner, paste0(q$fx$prefix, "_norm_complete"), paste0(q$fx$prefix, "_norm_complete"), config,
            "Normalization complete (correlation filtering skipped)"))
        live <- kit254Runs(q$fx)
        expect_identical(live, kit254HiRuns[[class]], info = paste("E8 precondition", class))
        out <- kit254Replay(q$fx, file.path(root, "e8-replay"))
        expect_null(out$error, info = paste("E8 replay", class))
        expect_null(out$value$design_statement)
        expect_true(paste0(q$fx$prefix, "_pearson_selection") %in% kit254Steps(out$record), info = paste("E8 record", class))
        expect_identical(out$runs, live, info = paste("E8 replay design", class))
        # E8, statement form (OD-1 (a), VAL-247; failing first in R1b): the replay states that the record shows the filter
        # applied and then skipped, that it reproduces the analysis as computed, and names per assay the samples left out.
        statement <- out$value$correlation_statement
        kit254Line("E8-statement", class, paste("statement", format(statement %||% "<none>")), paste("warned", length(out$warned)))
        expect_true(is.character(statement) && length(statement) == 1L, info = paste("E8 statement", class))
        expect_true(isTRUE(grepl("leaves out A_5 (below the threshold in LCMS_Pos)", statement, fixed = TRUE)),
            info = paste("E8 statement names A_5, left out of the design, failing in LCMS_Pos", class))
        expect_false(isTRUE(grepl("applied more than once", statement, fixed = TRUE)), info = paste("E8 one statement only", class))
        expect_true(isTRUE(grepl(format(hi), statement, fixed = TRUE)), info = paste("E8 statement names the threshold", class))
        expect_identical(out$warned, statement %||% "<none>", info = paste("E8 warning", class))
    })
}

# ------------------------------------------------------------------------------------------------
# G7, G8: the managed runner (both commits).
# ------------------------------------------------------------------------------------------------

for (class in c("lipid", "metab")) {
    test_that(paste("KIT-254 G7:", class, "managed: an Apply pressed while the final PCA is drawn completes; nothing is dropped"), {
        skip_if_not(workflowOwnerCapable())
        root <- withr::local_tempdir()
        seen <- kit254Notes()
        m <- kit254Managed(file.path(root, "m"), class, NULL, function(session, press, probe, fx) {
            probe$hold <- TRUE
            shiny::isolate(workflowOmicsFinalPca(fx$binding, probe$runner, fx$owner, kit254Final(probe$nd), "group", NULL))
            probe$display_running <- !is.null(shiny::isolate(probe$runner$step$active()))
            press("apply_correlation_filter", list(min_pearson_correlation_threshold = kit254T[[class]][["hi"]]))
            probe$hold <- FALSE
            probe$settled <- kit254Settle(session, probe)
        })
        kit254Line("G7", class, paste("statuses", paste(m$probe$statuses, collapse = ",")), paste("head", kit254Head(m$fx)))
        expect_true(isTRUE(m$probe$display_running) && isTRUE(m$probe$settled), info = paste("G7 precondition", class))
        expect_true(all(m$probe$statuses == "delivered"), info = paste("G7 statuses", paste(m$probe$statuses, collapse = ",")))
        kit254ExpectHead(m$fx, "_correlation_filtered", paste("G7", class))
        expect_identical(kit254Runs(m$fx), kit254HiRuns[[class]], info = paste("G7 design", class))
    })
}

test_that("KIT-254 G8: a managed Apply pressed while a publication retry is pending does not move the head", {
    skip_if_not(workflowOwnerCapable())
    root <- withr::local_tempdir()
    seen <- kit254Notes()
    for (class in c("lipid", "metab")) {
        q <- kit254Ran(file.path(root, class), class)
        expect_true(isTRUE(kit254Apply(q$fx, q$nd, threshold = kit254T[[class]][["hi"]])), info = paste("G8 Apply precondition", class))
        before <- kit254Ref(q$fx)$entity_id
        # The real runner's retry branch (R/utils_workflow_omics_actions.R:75-82): the actions are resolved at the press
        # only to compare operations, and a press that is not the retry waits its turn.
        runner <- new.env(parent = emptyenv())
        runner$retrying <- function() TRUE
        runner$progress <- function(...) list(set = function(...) invisible(NULL), inc = function(...) invisible(NULL), close = function(...) invisible(NULL))
        runner$run <- function(actions, done = NULL, failed = NULL, dropped = NULL) {
            runner$resolved <- if (is.function(actions)) actions() else actions
            invisible(NULL)
        }
        kit254Apply(q$fx, q$nd, runner = runner, threshold = kit254T[[class]][["lo"]])
        kit254Line("G8", class, paste("head", kit254Head(q$fx)), paste("unchanged", identical(kit254Ref(q$fx)$entity_id, before)))
        expect_true(length(runner$resolved) == 1L, info = paste("G8 precondition: the press reached the runner", class))
        kit254ExpectHead(q$fx, "_correlation_filtered", paste("G8", class))
        expect_identical(kit254Ref(q$fx)$entity_id, before, info = paste("G8", class))
    }
})


test_that("KIT-254 G11 (OD-1 (a), VAL-247): a proteomics record carries no correlation statement", {
    # The statement is built from the record's steps only for lipidomics and metabolomics prefixes; a proteomics record
    # with its correlation step followed by another step gets none (the builder is absent before the follow-up).
    builder <- get0("workflowReplayCorrelationStatement", envir = asNamespace("MultiScholaR"), inherits = FALSE)
    record <- list(steps = list(
        list(state = "normalized", operation = "normalization-ruv", parameters = list()),
        list(state = "correlation_filtered", operation = "pearson-selection",
            parameters = list(list(name = "threshold", value = 0.5)), input_entity_ids = list("a"), result_entity_id = "b"),
        list(state = "norm_complete", operation = "norm_complete", parameters = list())))
    expect_null(if (is.function(builder)) builder(record, NULL) else NULL)
})

# ------------------------------------------------------------------------------------------------
# OD-1 (a), VAL-247, and AUDIT-KIT-254-impl point 6 (I08, I09): the replay statements on base-form records. Failing first
# where a statement is expected (the statement code is absent on the tests-only commit); the guards pass on both.
# ------------------------------------------------------------------------------------------------

#' A record made with the wrappers in a given order on the run's result (as E8, KIT-212 G5): each element a threshold
#' (the selection, on the head) or "skip" (the store step on the head).
kit254Scripted <- function(fx, sequence) {
    owner <- kit254Owner(fx)
    config <- fx$workflow$config_list
    for (step in sequence) {
        if (identical(step, "skip")) {
            kit172Quiet(commitWorkflowOmicsStoreStep(owner, paste0(fx$prefix, "_norm_complete"), paste0(fx$prefix, "_norm_complete"),
                config, "Normalization complete (correlation filtering skipped)"))
        } else {
            kit172Quiet(runWorkflowOmicsReferenceCorrelation(fx$binding, owner, "group", step, config))
        }
    }
    invisible(fx)
}

kit254First <- "repeats both as recorded"
kit254Second <- "applied more than once"

for (class in c("lipid", "metab")) {
    test_that(paste("KIT-254 E15, E16, G12:", class, "a filter applied twice, and Apply, Skip, Apply, get the second statement only; no sample left, none"), {
        skip_if_not(workflowOwnerCapable())
        root <- withr::local_tempdir()
        seen <- kit254Notes()
        hi <- kit254T[[class]][["hi"]]
        lo <- kit254T[[class]][["lo"]]
        for (case in list(list(row = "E15", sequence = list(hi, lo)), list(row = "E16", sequence = list(hi, "skip", lo)))) {
            q <- kit254Ran(file.path(root, case$row), class)
            kit254Scripted(q$fx, case$sequence)
            out <- kit254Replay(q$fx, file.path(root, paste0(case$row, "-replay")))
            statement <- out$value$correlation_statement
            kit254Line(case$row, class, paste("statement", format(statement %||% "<none>")))
            expect_null(out$error, info = paste(case$row, "replay precondition", class))
            expect_true(isTRUE(grepl(kit254Second, statement, fixed = TRUE)) && isTRUE(grepl(paste0("(", format(lo), ")"), statement, fixed = TRUE)),
                info = paste(case$row, "the second statement with the recorded threshold", class))
            expect_false(isTRUE(grepl(kit254First, statement, fixed = TRUE)), info = paste(case$row, "never 'repeats both as recorded'", class))
            expect_identical(out$warned, statement %||% "<none>", info = paste(case$row, "warning", class))
        }
        # G12: the filter applied (t_lo, removing no sample) and then skipped: no statement.
        q <- kit254Ran(file.path(root, "g12"), class)
        kit254Scripted(q$fx, list(lo, "skip"))
        out <- kit254Replay(q$fx, file.path(root, "g12-replay"))
        expect_null(out$error, info = paste("G12 replay precondition", class))
        expect_identical(out$runs, kit254All, info = paste("G12 precondition: no sample left the design", class))
        expect_null(out$value$correlation_statement)
        expect_length(out$warned, 0L)
    })
}

for (class in c("lipid", "metab")) {
    test_that(paste("KIT-254 E17:", class, "two assays: the statement names the sample left out of the design, with the assay it failed in"), {
        skip_if_not(workflowOwnerCapable())
        root <- withr::local_tempdir()
        seen <- kit254Notes()
        values <- list(LCMS_Pos = kit172Values(172001L)$values, LCMS_Neg = kit172Values(172002L)$values)
        fx <- kit254Fixture(file.path(root, "a"), class, values)
        kit254Defer(fx)
        nd <- omicNormReactiveState()
        kit254Run(fx, nd)
        kit254ExpectHead(fx, "_normalized", paste("E17 run", class))
        testthat::local_mocked_bindings(pearsonCorForSamplePairs = function(theObject, ...) {
            runs <- as.character(methods::slot(theObject, "design_matrix")$Run)
            assays <- names(methods::slot(theObject, fx$binding$data_slot))
            stats::setNames(lapply(assays, function(assay) kit254E9Pairs(runs, assay)), assays)
        }, .package = "MultiScholaR")
        kit254Scripted(fx, list(0.5, "skip"))
        expect_identical(kit254Removed(kit254Runs(fx)), "A_1", info = paste("E17 precondition", class))
        out <- kit254Replay(fx, file.path(root, "replay"))
        statement <- out$value$correlation_statement
        kit254Line("E17", class, paste("statement", format(statement %||% "<none>")))
        expect_null(out$error, info = paste("E17 replay precondition", class))
        expect_true(isTRUE(grepl("leaves out A_1 (below the threshold in LCMS_Pos)", statement, fixed = TRUE)), info = paste("E17", class))
    })
}

test_that("KIT-254 E18, G13: with a merged design too, both statements and both warnings; the replayed states are unchanged", {
    skip_if_not(workflowOwnerCapable())
    root <- withr::local_tempdir()
    seen <- kit254Notes()
    q <- kit254Ran(file.path(root, "a"), "lipid")
    kit254Scripted(q$fx, list(kit254T$lipid[["hi"]], "skip"))
    kit254Da(q$fx)
    # G13 (both commits): the replay's states and head design are the same with the statement code bound to nothing.
    builder <- get0("workflowReplayCorrelationStatement", envir = asNamespace("MultiScholaR"), inherits = FALSE)
    plain <- local({
        if (is.function(builder)) testthat::local_mocked_bindings(workflowReplayCorrelationStatement = function(...) NULL, .package = "MultiScholaR")
        kit254Replay(q$fx, file.path(root, "plain"))
    })
    stated <- kit254Replay(q$fx, file.path(root, "stated"))
    expect_true(is.null(plain$error) && is.null(stated$error), info = "G13 replay precondition")
    expect_identical(stated$value$states, plain$value$states)
    expect_identical(stated$runs, plain$runs)
    # E18: KIT-222's recorder fed by a merged design during the replay (its consumers' check, bound to one collision).
    testthat::local_mocked_bindings(findOmicsDesignReplicateKeyCollisions = function(...) data.frame(key = "A_1", rule = "D1",
        runs = "A_1, A_2", stringsAsFactors = FALSE), .package = "MultiScholaR")
    out <- kit254Replay(q$fx, file.path(root, "both"))
    kit254Line("E18", "lipid", paste("design", format(out$value$design_statement %||% "<none>")),
        paste("correlation", format(out$value$correlation_statement %||% "<none>")))
    expect_null(out$error, info = "E18 replay precondition")
    expect_true(is.character(out$value$design_statement) && length(out$warned_design) == 1L, info = "E18 precondition: KIT-222's statement")
    expect_true(is.character(out$value$correlation_statement), info = "E18 the correlation statement beside it")
    expect_length(out$warned, 1L)
})
