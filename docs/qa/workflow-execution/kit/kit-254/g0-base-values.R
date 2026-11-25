# KIT-254 gate G0 (ticket r1.2, section 10): base values, recorded with no assertion, of the lipidomics and
# metabolomics correlation filter's Skip after an Apply and of a re-applied filter, at v0.9-staging e0747fc1.
# Observation only: branch kit-254-g0 carries this script and nothing under R/; it is never merged.
#
# Owners: F172 (kit172Values(172001L), one assay LCMS_Pos, groups A and B of 6) at QC complete (KIT-194's k194Fixture,
# copied below). Configurations: lipid with RUV skipped ("lipid"), metab with RUV skipped ("metab"), lipid with RUV
# applied (kit172LegacyRuv("manual"): k 2, 5 %; "lipid_ruv"). In process: the frame's runOmicNormPipeline,
# applyOmicNormCorrelation and skipOmicNormCorrelation (KIT-194's k194Run, k194Apply, k194Skip). Managed: the same owner
# (built from the F172 object, R1-01) under KIT-194 F1's worker journey (omicsManagedChainMocks; the real
# installWorkflowOmicsActions runner; the worker's entry point in this process), the tab's observers in a testServer.
# The file-local helpers of test-omics-norm-state-guards.R, test-kit207-summary-head-guard.R and test-kit212-da-current.R
# are COPIED here (k194*, k207Da, k207Save, k212Dirs, k212Replay); no test file is sourced.
#
# THRESHOLD RULE (fixed here before any value is read; ticket section 10). Per configuration, in process, a fresh
# owner per threshold: run, then Apply at each of 0, 0.3, 0.5, 0.7, 0.8, 0.9, 0.95, 0.99; record the design's samples
# kept. t_hi = the smallest listed threshold that removes at least one sample and keeps at least two per group; t_lo =
# the largest listed threshold below t_hi that removes strictly fewer samples than t_hi. If t_hi does not exist, stop
# rule 2 applies to that configuration: its cases run at t_hi = 0.5, t_lo = 0.3 and their sample differences are
# recorded as not shown on the fixture. If t_lo does not exist, (ii)'s cases run with t_lo = the smallest listed
# threshold (0) and are flagged. Both routes use the in-process t_hi and t_lo (R1-01: same F172 data).
# Supplementary (chooses nothing): each sample's best within-group partner correlation from the Apply at 0.
#
# PHASE 0 (R1-01; lifecycle check of the managed instrument): lipid and metab, managed: run; Apply(t); then, in a fresh
# owner, run; Apply(t); Skip (t = t_hi of the configuration). Passes when every action completes (heads
# <prefix>_correlation_filtered, then <prefix>_norm_complete) and the managed Apply's design sample list equals the
# in-process Apply's at the same t. The managed rows of the other cases are recorded only if phase 0 passes for that
# omic (else an instrument defect, VAL-54).
#
# Cases (fresh owner each): I = run, Apply(t_hi), Skip [(i); E1 to E4, E7]; II = run, Apply(t_hi), Apply(t_lo) [(ii);
# E5, E6, E7]; II0 = run, Apply(t_lo) alone; E13 = run, Apply(t_hi), Skip, Apply(t_lo); E14 = run, Apply(t_hi), DA,
# Skip, Summary Save; E10 = run, Apply(t_hi), Apply with a grouping column absent from the design; E11 = run, Apply(t_hi),
# Skip with commitWorkflowReferenceState refused (F105-02's mock); E9 = a two-assay owner (kit172Object with LCMS_Pos and
# LCMS_Neg), pearsonCorForSamplePairs bound to 2.3's constructed table restricted to the input design's samples, Apply
# at 0.5 twice; E12 (managed only) = run, Apply(t_hi) pressed, Apply(t_lo) pressed while the first is on the worker
# (pollWorkflowAction held at "running"), then released; C1 = run, Skip; C2 = run, Apply(t_lo), Apply(t_hi); C3 = run,
# Apply(t_hi) alone. I and II also managed (phase 0 permitting); E10, E11, E13, E14 lipid and metab only (RUV skipped).
#
# PREDICTIONS (base code, the ticket's reading; written before any value is read; also written to g0-predictions.tsv):
# see k254Predictions below.
#
# Usage (run-tests scripts lane, from the package root):
#   Rscript docs/qa/workflow-execution/kit/kit-254/g0-base-values.R <out>

args <- commandArgs(trailingOnly = TRUE)
k254Out <- if (length(args)) args[[1L]] else "kit254-g0"
dir.create(k254Out, recursive = TRUE, showWarnings = FALSE)

k254Predictions <- data.frame(stringsAsFactors = FALSE, rbind(
    c("SCAN", "t_hi, t_lo", "not predicted (stop rule 2 may apply); F172's within-group correlations are high"),
    c("PHASE0", "managed Apply; Apply then Skip", "complete: heads <prefix>_correlation_filtered, <prefix>_norm_complete; managed Apply design = in-process Apply design at the same t"),
    c("I", "head / history", "<prefix>_norm_complete; <prefix>_correlation_filtered in the active history before it"),
    c("I", "head design", "the t_hi Apply's kept samples (filtered)"),
    c("I", "head components equal the run result's", "FALSE (they are the Apply's)"),
    c("I", "summary original / after / removed", "N / N / 0 (both from the unfiltered run result)"),
    c("I", "final PCA reference samples", "N (the unfiltered run result)"),
    c("I", "DA source entity / samples", "the head (<prefix>_norm_complete); filtered count; status current"),
    c("I", "Save: integration object samples / study parameters", "filtered count; a pearson selection block with min_pearson_correlation_threshold = t_hi"),
    c("I", "record steps", "... <prefix>_correlation_filtered (<prefix>_pearson_selection at t_hi, preceded by pairs), <prefix>_norm_complete"),
    c("I", "replay", "completes, no statement; replay head design = the filtered samples (agrees with live)"),
    c("II", "second Apply source entity", "the t_hi state's entity (the head), not the run result"),
    c("II", "second Apply design", "t_hi's kept samples (one assay: max(t_hi, t_lo) = t_hi)"),
    c("II0", "design", "all samples (t_lo removes none by the rule)"),
    c("II", "record correlation step", "one <prefix>_pearson_selection at t_lo; input_entity_ids = result_entity_id"),
    c("II", "replay design", "t_lo on the run result: all samples; differs from live when t_hi removed a sample"),
    c("E13", "head / history / design / replay", "<prefix>_correlation_filtered; <prefix>_norm_complete after the head; t_hi's samples; replay all samples (differs)"),
    c("E14", "Save", "done 1, no failure; DA status current; filtered outputs"),
    c("E10", "head / fields", "the t_hi state; correlation_filtering_complete TRUE, correlation_filtered_obj = t_hi entity; an error notice"),
    c("E11", "head / fields", "the t_hi state; lipid: flags put back to t_hi's; metab: unchanged (t_hi's)"),
    c("E9", "design after two Applies", "A_1 and A_2 removed (cascade); one Apply removes A_1 only; or an error reading the second press's input (F2), recorded"),
    c("E12", "second press source / design", "the t_hi state (queued press reads the head); t_hi's samples; the queue is held (instrument)"),
    c("C1", "design / record / replay", "all samples; no pearson selection step; replay all samples"),
    c("C2 vs C3", "design", "equal (one assay)")))
names(k254Predictions) <- c("case", "observable", "predicted")
utils::write.table(k254Predictions, file.path(k254Out, "g0-predictions.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

suppressMessages(pkgload::load_all(".", quiet = TRUE, export_all = TRUE, helpers = TRUE))
library(testthat)

# ---- recording ----

k254Rows <- new.env(parent = emptyenv())
k254Rows$rows <- list()
k254Text <- function(x) gsub("[\t\r\n]+", " ", gsub("[0-9a-f]{64}", "<id>", x))
k254Record <- function(case, config, route, observable, value) {
    value <- k254Text(paste(as.character(unlist(value)), collapse = ","))
    if (!length(value) || !nzchar(value)) value <- "<none>"
    k254Rows$rows[[length(k254Rows$rows) + 1L]] <- data.frame(case = case, config = config, route = route,
        observable = observable, value = value, stringsAsFactors = FALSE)
    cat(sprintf("KIT-254 G0\t%s\t%s\t%s\t%s\t%s\n", case, config, route, observable, value))
}
k254Write <- function() {
    rows <- do.call(rbind, k254Rows$rows)
    utils::write.table(rows, file.path(k254Out, "g0-base-values.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
}
k254Quiet <- function(code) {
    value <- NULL
    utils::capture.output(value <- suppressWarnings(suppressMessages(code)))
    value
}
k254Error <- function(e) paste0("error [", paste(class(e), collapse = "/"), "] ", conditionMessage(e))

# ---- copied unchanged from tests/testthat/test-omics-norm-state-guards.R at e0747fc1 (lines 18-31, 101-188) ----

k194Omic <- function(class) c(lipid = "lipidomics", metab = "metabolomics")[[class]]

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
k194Owner <- function(fixture) fixture$workflow$state_manager
k194History <- function(fixture) as.character(unlist(k194Owner(fixture)$getHistory()))
k194Head <- function(fixture) k194Owner(fixture)$current_state
k194Ref <- function(fixture, state = NULL) getWorkflowStateReference(k194Owner(fixture), state)

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

k194Input <- function(method = "quantile", itsd = FALSE, ruv = NULL, threshold = 0.5) {
    list(apply_itsd = itsd, itsd_aggregation = "median", log_offset = 0, norm_method = method,
        ruv_mode = if (is.null(ruv)) "skip" else ruv$ruv_mode, auto_percentage_min = 1, auto_percentage_max = 20,
        ruv_grouping_variable = "group", separation_metric = "max_difference", k_penalty_weight = 0.5,
        adaptive_k_penalty = TRUE, ruv_k = ruv$manual_k %||% 1L, ruv_percentage = ruv$manual_percentage %||% 20,
        min_pearson_correlation_threshold = threshold)
}

k194Run <- function(fixture, nd, input = k194Input(), runner = NULL) {
    testthat::local_mocked_bindings(workflowOmicsManagedActions = function(omic_type) runner, .package = "MultiScholaR")
    kit172Quiet(shiny::isolate(runOmicNormPipeline(input, fixture$workflow, list(), nd, fixture$addLog,
        function(...) list(color_var = "group", shape_var = NULL), fixture$config)))
}

k194Apply <- function(fixture, nd, runner = NULL, threshold = 0.5) {
    testthat::local_mocked_bindings(workflowOmicsManagedActions = function(omic_type) runner, .package = "MultiScholaR")
    kit172Quiet(shiny::isolate(applyOmicNormCorrelation(k194Input(threshold = threshold), fixture$workflow, nd, fixture$addLog,
        fixture$config)))
}

k194Skip <- function(fixture, nd) kit172Quiet(shiny::isolate(skipOmicNormCorrelation(fixture$workflow, nd, fixture$addLog, fixture$config)))

# ---- copied unchanged from tests/testthat/test-kit212-da-current.R at e0747fc1 (lines 270-295, 337-343, 447-475;
# k212Replay extended at its end to return the prepared store and the replay's value, and the replay's head read) ----

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

k207Da <- function(fx) {
    kit172Quiet(runWorkflowOmicsReferenceDaAnalysis(fx$binding, k194Owner(fx), k194Ref(fx),
        list(formula_string = "~ 0 + group", da_q_val_thresh = 0.05, treat_lfc_cutoff = 0), kit172Contrasts()))
    invisible(fx)
}

k212Dirs <- function(root, omic) {
    paths <- list(source_dir = file.path(root, "scripts", omic), integration_dir = file.path(root, "integration"), base_dir = root,
        da_output_dir = file.path(root, "da_output"), publication_graphs_dir = file.path(root, "publication_graphs"),
        results_summary_dir = file.path(root, "results_summary"))
    for (path in paths) dir.create(path, recursive = TRUE, showWarnings = FALSE)
    stats::setNames(list(paths), omic)
}

k254Replay <- function(fx, root) {
    owner <- k194Owner(fx)
    dir.create(root, recursive = TRUE, showWarnings = FALSE)
    json <- file.path(root, "analysis_record.json")
    writeWorkflowAnalysisRecord(owner, json)
    record <- jsonlite::read_json(json, simplifyVector = FALSE)
    copy <- file.path(root, "copy")
    dir.create(copy)
    stopifnot(system2("cp", c("-a", file.path(fx$source_dir, "."), copy)) == 0L)
    harness <- new.env(parent = globalenv())
    sys.source(file.path("tools", "profiling", "workflow_kit.R"), harness)
    api <- function(name) get(name, asNamespace("MultiScholaR"))
    store_root <- dirname(dirname(list.files(fx$source_dir, "^store[.]json$", recursive = TRUE, all.files = TRUE, full.names = TRUE)[[1L]]))
    copy_root <- file.path(copy, substring(store_root, nchar(fx$source_dir) + 2L))
    prepared <- harness$kitPrepareReplayProject(copy_root, c("initial", paste0(fx$prefix, "_raw_data_s4")), api)
    inputs <- stats::setNames(c(file.path(fx$source_dir, "import_sources", "LCMS_Pos_export.txt"), file.path(fx$source_dir, "design_matrix.tab")),
        c(paste0(fx$binding$design_input_prefix, 1L), "design_file"))
    scratch <- file.path(root, "scratch")
    dir.create(scratch)
    value <- NULL
    error <- tryCatch({
        value <- kit172Quiet(replayWorkflowAnalysisRecord(record, prepared, inputs, fx$paths, scratch))
        NULL
    }, error = function(e) e)
    # KIT-254: the replayed head's design samples, read from the replay store after the replay closed it.
    runs <- NULL
    if (is.null(error)) {
        identity <- readSavedStoreIdentity(file.path(prepared, ".workflow-entities", "store.json"))
        manager <- newDIAStateAdapter(newWorkflowEntityStore(prepared, identity$project_id, identity$workflow_id), reopen = TRUE)
        runs <- tryCatch(k254RunsOf(fx$binding, manager, NULL), error = function(e) paste("<unreadable>", conditionMessage(e)))
        head <- manager$current_state
        closeDIAStateAdapter(manager)
        runs <- c(head = head, runs = paste(runs, collapse = ","))
    }
    list(record = record, error = error, value = value, replayed = runs)
}

# ---- KIT-254 readers ----

k254RunsOf <- function(binding, manager, reference) {
    meta <- readWorkflowOmicsMetadata(binding, manager, reference)
    if (is.null(meta)) return("<no state>")
    sort(as.character(as.data.frame(meta$design_matrix)[["Run"]]))
}
k254Runs <- function(fx, reference = NULL) {
    if (is.character(reference)) reference <- k194Ref(fx, reference)
    k254RunsOf(fx$binding, k194Owner(fx), reference)
}
k254All <- sort(kit172Design()$Run)
k254Removed <- function(runs) setdiff(k254All, runs)
k254Components <- function(fx, reference) readWorkflowStateMetadata(k194Owner(fx), reference)$data$components
k254Final <- function(nd) shiny::isolate(nd$ruv_corrected_obj %||% nd$post_norm_obj)
k254Entity <- function(reference) if (isSavedWorkflowReference(reference)) substr(reference$entity_id, 1L, 12L) else "<not a reference>"

#' The summary's three numbers as R/mod_omic_norm_panels.R:1314-1319 computes them, or the text it would show.
k254Summary <- function(fx, nd) {
    values <- shiny::isolate(shiny::reactiveValuesToList(nd))
    if (!isTRUE(values$correlation_filtering_complete)) return("Apply correlation filter to see results...")
    if (is.null(values$correlation_results) || !length(values$correlation_results)) return("No correlation results available.")
    original <- values$ruv_corrected_obj %||% values$post_norm_obj
    filtered <- values$correlation_filtered_obj
    if (is.null(original) || is.null(filtered)) return("<no sample lines>")
    read <- fx$config$binding$read_metadata
    o <- nrow(read(k194Owner(fx), original)$design_matrix)
    f <- nrow(read(k194Owner(fx), filtered)$design_matrix)
    sprintf("original=%d;after=%d;removed=%d", o, f, o - f)
}

k254PcaSamples <- function(fx, nd) {
    values <- shiny::isolate(shiny::reactiveValuesToList(nd))
    reference <- values$correlation_filtered_obj %||% values$ruv_corrected_obj %||% values$post_norm_obj
    if (is.null(reference)) return("<none>")
    length(k254Runs(fx, reference))
}

k254Fields <- function(nd) {
    values <- shiny::isolate(shiny::reactiveValuesToList(nd))
    sprintf("complete=%s;filtered_obj=%s;results=%s", format(values$correlation_filtering_complete),
        if (is.null(values$correlation_filtered_obj)) "NULL" else k254Entity(values$correlation_filtered_obj),
        paste(names(values$correlation_results) %||% "<empty>", collapse = "+"))
}

k254Steps <- function(record) paste(vapply(record$steps, function(step) paste0(step$state, ":", step$operation), ""), collapse = " ")

k254CorrelationStep <- function(record, prefix) {
    steps <- Filter(function(step) identical(step$operation, paste0(prefix, "_pearson_selection")), record$steps)
    if (!length(steps)) return("<no pearson selection step>")
    vapply(steps, function(step) {
        threshold <- Filter(function(p) identical(p$name, "min_pearson_correlation_threshold"), step$parameters)
        threshold <- if (length(threshold)) paste(unlist(threshold[[1L]]$value), collapse = ",") else "<none>"
        input <- if (length(step$input_entity_ids)) substr(step$input_entity_ids[[1L]], 1L, 12L) else "<none>"
        sprintf("state=%s;threshold=%s;input=%s;result=%s;input_is_result=%s;preceded_by=%s", step$state, threshold, input,
            substr(step$result_entity_id, 1L, 12L), identical(input, substr(step$result_entity_id, 1L, 12L)),
            paste(vapply(step$preceded_by %||% list(), function(p) p$operation, ""), collapse = "+"))
    }, "")
}

k254DaSource <- function(fx) {
    owner <- k194Owner(fx)
    withWorkflowStateScope(owner, function() {
        da <- owner$.context$snapshot$da_results
        if (!isSavedWorkflowReference(da)) return("<no DA>")
        source <- resolveWorkflowStateValue(owner$.context, da)$source
        name <- Find(function(state) identical(owner$.context$snapshot$states[[state]]$entity_id, source$entity_id),
            names(owner$.context$snapshot$states))
        sprintf("entity=%s;state=%s", substr(source$entity_id, 1L, 12L), name %||% "<not in states>")
    })
}

k254RecordRun <- function(case, config, route, fx, nd, root, replay = TRUE) {
    prefix <- fx$prefix
    k254Record(case, config, route, "head", k194Head(fx))
    k254Record(case, config, route, "head entity", k254Entity(k194Ref(fx)))
    k254Record(case, config, route, "history", paste(sub(paste0("^", prefix, "_"), "", k194History(fx)), collapse = " "))
    runs <- k254Runs(fx)
    k254Record(case, config, route, "head design n", length(runs))
    k254Record(case, config, route, "head design removed", paste(k254Removed(runs), collapse = ",") %||% "")
    final <- k254Final(nd)
    k254Record(case, config, route, "final (run result) entity", k254Entity(final))
    if (isSavedWorkflowReference(final)) {
        k254Record(case, config, route, "head components equal the run result's",
            identical(k254Components(fx, k194Ref(fx)), k254Components(fx, final)))
    }
    k254Record(case, config, route, "fields", k254Fields(nd))
    k254Record(case, config, route, "summary", k254Summary(fx, nd))
    k254Record(case, config, route, "final PCA samples", k254PcaSamples(fx, nd))
    if (!replay) return(invisible(NULL))
    out <- k254Replay(fx, file.path(root, "replay"))
    k254Record(case, config, route, "record steps", k254Steps(out$record))
    k254Record(case, config, route, "record correlation step", k254CorrelationStep(out$record, prefix))
    k254Record(case, config, route, "replay", if (is.null(out$error)) "completed" else k254Error(out$error))
    k254Record(case, config, route, "replay statement", out$value$design_statement %||% "<none>")
    if (!is.null(out$replayed)) {
        k254Record(case, config, route, "replay head", out$replayed[["head"]])
        k254Record(case, config, route, "replay design removed",
            paste(k254Removed(strsplit(out$replayed[["runs"]], ",", fixed = TRUE)[[1L]]), collapse = ","))
        k254Record(case, config, route, "replay equals live design", identical(out$replayed[["runs"]], paste(runs, collapse = ",")))
    }
    invisible(NULL)
}

# ---- configurations, fixtures ----

k254Configs <- list(
    lipid = list(class = "lipid", ruv = NULL),
    metab = list(class = "metab", ruv = NULL),
    lipid_ruv = list(class = "lipid", ruv = kit172LegacyRuv("manual")))
k254Root <- tempfile("kit254-g0-")
dir.create(k254Root)
k254Count <- new.env(parent = emptyenv())
k254Count$n <- 0L
k254NewRoot <- function(label) {
    k254Count$n <- k254Count$n + 1L
    file.path(k254Root, sprintf("%03d-%s", k254Count$n, label))
}

#' A fixture after "Run Normalization Pipeline" in process; the run's head and flag recorded as the case's precondition.
k254InProcess <- function(config, label) {
    fx <- k194Fixture(file.path(k254NewRoot(label), "project"), config$class)
    nd <- omicNormReactiveState()
    k194Run(fx, nd, k194Input(ruv = config$ruv))
    list(fx = fx, nd = nd, ran = k194Head(fx) %in% paste0(fx$prefix, c("_normalized", "_ruv_corrected")) &&
        isTRUE(shiny::isolate(nd$normalization_complete)))
}

k254Case <- function(case, name, route, code) {
    started <- proc.time()[["elapsed"]]
    result <- tryCatch(code(), error = function(e) {
        k254Record(case, name, route, "INSTRUMENT ERROR", k254Error(e))
        NULL
    })
    k254Record(case, name, route, "seconds", round(proc.time()[["elapsed"]] - started, 1))
    k254Write()
    result
}

# ---- SCAN ----

k254Thresholds <- c(0, 0.3, 0.5, 0.7, 0.8, 0.9, 0.95, 0.99)
k254Scan <- list()
k254BestPartner <- list()
for (name in names(k254Configs)) {
    config <- k254Configs[[name]]
    removed <- list()
    for (t in k254Thresholds) {
        k254Case("SCAN", name, "in_process", function() {
            seen <- k194Notes()
            run <- k254InProcess(config, paste0("scan-", name))
            on.exit(k194Close(run$fx), add = TRUE)
            k254Record("SCAN", name, "in_process", sprintf("run completed (t=%s)", t), run$ran)
            k194Apply(run$fx, run$nd, threshold = t)
            ok <- identical(k194Head(run$fx), paste0(run$fx$prefix, "_correlation_filtered"))
            runs <- if (ok) k254Runs(run$fx) else NULL
            k254Record("SCAN", name, "in_process", sprintf("t=%s head", t), k194Head(run$fx))
            k254Record("SCAN", name, "in_process", sprintf("t=%s removed", t), if (ok) paste(k254Removed(runs), collapse = ",") else "<no commit>")
            k254Record("SCAN", name, "in_process", sprintf("t=%s notes", t), paste(seen$notes, collapse = " | "))
            removed[[as.character(t)]] <<- if (ok) k254Removed(runs) else NA
            if (t == 0 && ok) {
                pairs <- shiny::isolate(run$nd$correlation_results)
                rows <- do.call(rbind, lapply(names(pairs), function(assay) {
                    p <- as.data.frame(pairs[[assay]])
                    long <- rbind(data.frame(sample = p$Run.x, partner = p$Run.y, r = p$pearson_correlation),
                        data.frame(sample = p$Run.y, partner = p$Run.x, r = p$pearson_correlation))
                    best <- do.call(rbind, lapply(split(long, long$sample), function(s) s[which.max(s$r), , drop = FALSE]))
                    data.frame(config = name, assay = assay, sample = best$sample, best_partner = best$partner,
                        best_r = sprintf("%.6f", best$r), stringsAsFactors = FALSE)
                }))
                k254BestPartner[[name]] <<- rows
            }
        })
    }
    k254Scan[[name]] <- removed
}
k254Choose <- function(removed) {
    groups <- stats::setNames(kit172Design()$group, kit172Design()$Run)
    keeps_two <- function(r) all(table(factor(groups[setdiff(k254All, r)], levels = c("A", "B"))) >= 2L)
    n <- vapply(removed, function(r) if (length(r) == 1L && is.na(r)) NA_integer_ else length(r), integer(1))
    ts <- as.numeric(names(removed))
    hi <- which(!is.na(n) & n >= 1L & vapply(removed, function(r) !(length(r) == 1L && is.na(r)) && keeps_two(r), logical(1)))
    if (!length(hi)) return(list(t_hi = 0.5, t_lo = 0.3, stop_rule_2 = TRUE, t_lo_flag = FALSE))
    t_hi <- min(ts[hi])
    lo <- which(ts < t_hi & !is.na(n) & n < n[ts == t_hi])
    if (!length(lo)) return(list(t_hi = t_hi, t_lo = min(ts), stop_rule_2 = FALSE, t_lo_flag = TRUE))
    list(t_hi = t_hi, t_lo = max(ts[lo]), stop_rule_2 = FALSE, t_lo_flag = FALSE)
}
k254T <- lapply(k254Scan, k254Choose)
for (name in names(k254T)) {
    k254Record("SCAN", name, "in_process", "t_hi", k254T[[name]]$t_hi)
    k254Record("SCAN", name, "in_process", "t_lo", k254T[[name]]$t_lo)
    k254Record("SCAN", name, "in_process", "stop rule 2", k254T[[name]]$stop_rule_2)
    k254Record("SCAN", name, "in_process", "t_lo flagged (none below t_hi removes fewer)", k254T[[name]]$t_lo_flag)
}
if (length(k254BestPartner)) utils::write.table(do.call(rbind, k254BestPartner), file.path(k254Out, "g0-best-partner.tsv"),
    sep = "\t", quote = FALSE, row.names = FALSE)
k254Write()

# ---- the managed harness (KIT-194 F1's worker journey, with an Apply and a Skip observer) ----

#' The tab's run, Apply and Skip observers on a managed session (the real runner), as omicsTabServer
#' (helper-workflow-omics-worker.R:193-238) mounts them, with the Skip added.
k254Server <- function(fx, probe) {
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

#' One managed case: a fresh F172 owner at QC complete under the worker-journey mocks. Inside the server the tab's
#' inputs are set and the run pressed, then `steps(session, press, probe, fx)`; `press(button, extra)` sets `extra`, then
#' the button, and flushes. After the server ends, `after(fx, nd, probe)` reads the owner (still open). The managed
#' actions' input heads are recorded in `probe$inputs`; `probe$hold` keeps the worker's poll at "running" while TRUE.
k254Managed <- function(config, label, steps, after) {
    fx <- k194Fixture(file.path(k254NewRoot(label), "project"), config$class)
    on.exit(k194Close(fx), add = TRUE)
    probe <- new.env()
    probe$hold <- FALSE
    probe$inputs <- character()
    omicsManagedChainMocks(fx, probe)
    dispatch <- get("runWorkflowAction", envir = asNamespace("MultiScholaR"))
    testthat::local_mocked_bindings(
        runWorkflowAction = function(operation, input_references, parameters, owner, reservation, completion_context = NULL) {
            probe$inputs <- c(probe$inputs, sprintf("%s<-%s(%s)", operation, k254Entity(getWorkflowStateReference(fx$owner)),
                fx$owner$current_state))
            dispatch(operation, input_references, parameters, owner, reservation, completion_context)
        },
        pollWorkflowAction = function(owner) if (isTRUE(probe$hold)) list(status = "running") else list(status = "complete"),
        .package = "MultiScholaR")
    suppressWarnings(suppressMessages(utils::capture.output(shiny::testServer(k254Server(fx, probe), {
        session$userData$probe_owner <- fx$owner
        presses <- new.env(parent = emptyenv())
        press <- function(button, extra = list()) {
            presses[[button]] <- (presses[[button]] %||% 0L) + 1L
            if (length(extra)) do.call(session$setInputs, extra)
            do.call(session$setInputs, stats::setNames(list(presses[[button]]), button))
            session$flushReact()
        }
        do.call(session$setInputs, k194Input(ruv = config$ruv))
        press("run_normalization")
        probe$settled_run <- k254Settle(session, probe)
        probe$ran <- isTRUE(fx$owner$current_state %in% paste0(fx$prefix, c("_normalized", "_ruv_corrected")))
        steps(session, press, probe, fx)
    }))))
    after(fx, probe$nd, probe)
}

k254Settle <- function(session, probe, rounds = 40L) {
    for (i in seq_len(rounds)) {
        if (is.null(shiny::isolate(probe$runner$step$active()))) break
        session$elapse(100)
    }
    is.null(shiny::isolate(probe$runner$step$active()))
}

k254ManagedCommon <- function(case, name, fx, probe) {
    k254Record(case, name, "managed", "run completed", isTRUE(probe$ran) && isTRUE(probe$settled_run))
    k254Record(case, name, "managed", "worker inputs", paste(probe$inputs, collapse = " ; "))
    k254Record(case, name, "managed", "statuses", paste(probe$statuses, collapse = ","))
    k254Record(case, name, "managed", "errors in log", paste(grep("ERROR", fx$log$lines, value = TRUE), collapse = " | "))
}

#' The source of every in-process correlation (the head when the wrapper is called), recorded in `probe$sources`.
k254SourceProbe <- function(probe, env = parent.frame()) {
    wrapper <- get("runWorkflowOmicsReferenceCorrelation", envir = asNamespace("MultiScholaR"))
    probe$sources <- character()
    testthat::local_mocked_bindings(runWorkflowOmicsReferenceCorrelation = function(binding, manager, ...) {
        probe$sources <- c(probe$sources, sprintf("%s(%s)", k254Entity(getWorkflowStateReference(manager)), manager$current_state))
        wrapper(binding, manager, ...)
    }, .package = "MultiScholaR", .env = env)
    invisible(probe)
}

k254RdsSamples <- function(rec) {
    if (!length(rec$rds)) return("<none>")
    paste(vapply(rec$rds, function(call) {
        object <- call$object
        if (isS4(object) && "design_matrix" %in% methods::slotNames(object)) {
            sprintf("%s:%d samples", basename(call$path), nrow(methods::slot(object, "design_matrix")))
        } else sprintf("%s:%s", basename(call$path), class(object)[[1L]])
    }, ""), collapse = "; ")
}

#' The DA on the head, then the Summary's Save, recorded.
k254DaAndSave <- function(case, name, route, fx) {
    k207Da(fx)
    k254Record(case, name, route, "DA source", k254DaSource(fx))
    k254Record(case, name, route, "DA source design n", length(k254Runs(fx)))
    k254Record(case, name, route, "DA status", workflowDaResultStatus(k194Owner(fx)))
    dirs <- k212Dirs(k254NewRoot(paste0(tolower(case), "-save")), fx$omic)
    rec <- k207Save(fx$workflow, fx$omic, dirs)
    k254Record(case, name, route, "Save done / failed", sprintf("done=%d;failed=%s", rec$done,
        paste(vapply(rec$failed, function(e) class(e)[[1L]], ""), collapse = ",")))
    k254Record(case, name, route, "Save notes", paste(rec$notes, collapse = " | "))
    k254Record(case, name, route, "Save integration object", k254RdsSamples(rec))
    params <- file.path(dirs[[fx$omic]]$source_dir, "study_parameters.txt")
    lines <- if (file.exists(params)) grep("pearson|correlation", readLines(params, warn = FALSE), value = TRUE, ignore.case = TRUE) else "<no file>"
    k254Record(case, name, route, "study parameters correlation lines", paste(trimws(lines), collapse = " | "))
}

#' E11's Skip with the store commit refused (F105-02's mock, test-omic-norm-module.R:516-519).
k254SkipRefused <- function(fx, nd) {
    testthat::local_mocked_bindings(commitWorkflowReferenceState = function(...) stop("kit254: the store commit is refused"),
        .package = "MultiScholaR")
    k194Skip(fx, nd)
}

#' E9's owner: F172 with two assays (LCMS_Pos from seed 172001, LCMS_Neg from seed 172002), at QC complete, with
#' k194Fixture's fields.
k254Fixture2 <- function(root, class) {
    values <- list(LCMS_Pos = kit172Values(172001L)$values, LCMS_Neg = kit172Values(172002L)$values)
    fixture <- k194PreNorm(root, values, class)
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

#' Ticket 2.3's example on F172's group A (a = A_1, b = A_2): every within-group pair of `runs` (the input design's
#' samples). LCMS_Pos: every pair with A_1 at 0.1, every other at 0.99. LCMS_Neg: A_1-A_2 and the pairs among A_3 to
#' A_6 at 0.99, the other A pairs at 0.1; group B at 0.99.
k254E9Pairs <- function(runs, assay) {
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

#' `pearsonCorForSamplePairs` bound in the namespace to E9's table (as MCI-021.4 binds it, test-module-ci-lipid-norm.R:222-224).
k254E9Mock <- function(env = parent.frame()) {
    testthat::local_mocked_bindings(pearsonCorForSamplePairs = function(theObject, ...) {
        runs <- as.character(methods::slot(theObject, "design_matrix")$Run)
        slot <- if (methods::is(theObject, "LipidomicsAssayData")) "lipid_data" else "metabolite_data"
        assays <- names(methods::slot(theObject, slot))
        stats::setNames(lapply(assays, function(assay) k254E9Pairs(runs, assay)), assays)
    }, .package = "MultiScholaR", .env = env)
}

k254Applied <- function(case, name, route, fx, label) {
    k254Record(case, name, route, paste(label, "head"), k194Head(fx))
    k254Record(case, name, route, paste(label, "entity"), k254Entity(k194Ref(fx)))
    k254Record(case, name, route, paste(label, "removed"), paste(k254Removed(k254Runs(fx)), collapse = ","))
}

# ---- PHASE 0 ----

k254Phase0 <- list()
for (name in c("lipid", "metab")) {
    config <- k254Configs[[name]]
    t <- k254T[[name]]$t_hi
    expected <- k254Case("PHASE0", name, "in_process", function() {
        seen <- k194Notes()
        run <- k254InProcess(config, paste0("p0-ip-", name))
        on.exit(k194Close(run$fx), add = TRUE)
        k194Apply(run$fx, run$nd, threshold = t)
        k254Applied("PHASE0", name, "in_process", run$fx, sprintf("Apply(%s)", t))
        k254Runs(run$fx)
    })
    apply_ok <- k254Case("PHASE0", name, "managed", function() {
        seen <- k194Notes()
        k254Managed(config, paste0("p0-apply-", name), function(session, press, probe, fx) {
            press("apply_correlation_filter", list(min_pearson_correlation_threshold = t))
            probe$settled <- k254Settle(session, probe)
        }, function(fx, nd, probe) {
            k254ManagedCommon("PHASE0", name, fx, probe)
            k254Applied("PHASE0", name, "managed", fx, sprintf("Apply(%s)", t))
            runs <- k254Runs(fx)
            ok <- isTRUE(probe$settled) && identical(k194Head(fx), paste0(fx$prefix, "_correlation_filtered")) && identical(runs, expected)
            k254Record("PHASE0", name, "managed", "Apply completed with the in-process sample list", ok)
            ok
        })
    })
    skip_ok <- k254Case("PHASE0", name, "managed", function() {
        seen <- k194Notes()
        k254Managed(config, paste0("p0-skip-", name), function(session, press, probe, fx) {
            press("apply_correlation_filter", list(min_pearson_correlation_threshold = t))
            probe$settled <- k254Settle(session, probe)
            press("skip_correlation_filter")
            probe$settled_skip <- k254Settle(session, probe)
        }, function(fx, nd, probe) {
            k254ManagedCommon("PHASE0", name, fx, probe)
            k254Applied("PHASE0", name, "managed", fx, "Apply then Skip")
            ok <- isTRUE(probe$settled) && isTRUE(probe$settled_skip) && identical(k194Head(fx), paste0(fx$prefix, "_norm_complete"))
            k254Record("PHASE0", name, "managed", "Apply then Skip completed", ok)
            ok
        })
    })
    k254Phase0[[name]] <- isTRUE(apply_ok) && isTRUE(skip_ok)
    k254Record("PHASE0", name, "managed", "PHASE 0 PASSED", k254Phase0[[name]])
}
k254Phase0$lipid_ruv <- k254Phase0$lipid

# ---- the cases ----

for (name in names(k254Configs)) {
    config <- k254Configs[[name]]
    t_hi <- k254T[[name]]$t_hi
    t_lo <- k254T[[name]]$t_lo
    ruv_skipped <- is.null(config$ruv)

    # I: run, Apply(t_hi), Skip (in process).
    k254Case("I", name, "in_process", function() {
        seen <- k194Notes()
        run <- k254InProcess(config, paste0("i-", name))
        on.exit(k194Close(run$fx), add = TRUE)
        k254Record("I", name, "in_process", "run completed", run$ran)
        k194Apply(run$fx, run$nd, threshold = t_hi)
        k254Applied("I", name, "in_process", run$fx, sprintf("Apply(%s)", t_hi))
        k194Skip(run$fx, run$nd)
        k254RecordRun("I", name, "in_process", run$fx, run$nd, k254NewRoot(paste0("i-replay-", name)))
        k254DaAndSave("I", name, "in_process", run$fx)
        k254Record("I", name, "in_process", "notes", paste(seen$notes, collapse = " | "))
    })

    # II: run, Apply(t_hi), Apply(t_lo) (in process); II0: run, Apply(t_lo) alone.
    k254Case("II", name, "in_process", function() {
        seen <- k194Notes()
        probe <- k254SourceProbe(new.env())
        run <- k254InProcess(config, paste0("ii-", name))
        on.exit(k194Close(run$fx), add = TRUE)
        k254Record("II", name, "in_process", "run completed", run$ran)
        k254Record("II", name, "in_process", "run result entity", k254Entity(k254Final(run$nd)))
        k194Apply(run$fx, run$nd, threshold = t_hi)
        k254Applied("II", name, "in_process", run$fx, sprintf("Apply(%s)", t_hi))
        k194Apply(run$fx, run$nd, threshold = t_lo)
        k254Applied("II", name, "in_process", run$fx, sprintf("then Apply(%s)", t_lo))
        k254Record("II", name, "in_process", "correlation sources (entity(head) per press)", paste(probe$sources, collapse = " ; "))
        k254RecordRun("II", name, "in_process", run$fx, run$nd, k254NewRoot(paste0("ii-replay-", name)))
    })
    k254Case("II0", name, "in_process", function() {
        seen <- k194Notes()
        run <- k254InProcess(config, paste0("ii0-", name))
        on.exit(k194Close(run$fx), add = TRUE)
        k194Apply(run$fx, run$nd, threshold = t_lo)
        k254Applied("II0", name, "in_process", run$fx, sprintf("Apply(%s) alone", t_lo))
    })

    # Controls: C1 run, Skip; C2 run, Apply(t_lo), Apply(t_hi); C3 run, Apply(t_hi) alone.
    k254Case("C1", name, "in_process", function() {
        seen <- k194Notes()
        run <- k254InProcess(config, paste0("c1-", name))
        on.exit(k194Close(run$fx), add = TRUE)
        k194Skip(run$fx, run$nd)
        k254RecordRun("C1", name, "in_process", run$fx, run$nd, k254NewRoot(paste0("c1-replay-", name)))
    })
    k254Case("C2", name, "in_process", function() {
        seen <- k194Notes()
        run <- k254InProcess(config, paste0("c2-", name))
        on.exit(k194Close(run$fx), add = TRUE)
        k194Apply(run$fx, run$nd, threshold = t_lo)
        k194Apply(run$fx, run$nd, threshold = t_hi)
        k254Applied("C2", name, "in_process", run$fx, sprintf("Apply(%s) then Apply(%s)", t_lo, t_hi))
    })
    k254Case("C3", name, "in_process", function() {
        seen <- k194Notes()
        run <- k254InProcess(config, paste0("c3-", name))
        on.exit(k194Close(run$fx), add = TRUE)
        k194Apply(run$fx, run$nd, threshold = t_hi)
        k254Applied("C3", name, "in_process", run$fx, sprintf("Apply(%s) alone", t_hi))
    })

    if (ruv_skipped) {
        # E13: run, Apply(t_hi), Skip, Apply(t_lo).
        k254Case("E13", name, "in_process", function() {
            seen <- k194Notes()
            run <- k254InProcess(config, paste0("e13-", name))
            on.exit(k194Close(run$fx), add = TRUE)
            k194Apply(run$fx, run$nd, threshold = t_hi)
            k194Skip(run$fx, run$nd)
            k194Apply(run$fx, run$nd, threshold = t_lo)
            k254RecordRun("E13", name, "in_process", run$fx, run$nd, k254NewRoot(paste0("e13-replay-", name)))
        })
        # E14: run, Apply(t_hi), DA, Skip, Save.
        k254Case("E14", name, "in_process", function() {
            seen <- k194Notes()
            run <- k254InProcess(config, paste0("e14-", name))
            on.exit(k194Close(run$fx), add = TRUE)
            k194Apply(run$fx, run$nd, threshold = t_hi)
            k207Da(run$fx)
            k254Record("E14", name, "in_process", "DA source (before Skip)", k254DaSource(run$fx))
            k194Skip(run$fx, run$nd)
            k254Record("E14", name, "in_process", "head after Skip", k194Head(run$fx))
            k254Record("E14", name, "in_process", "DA status after Skip", workflowDaResultStatus(k194Owner(run$fx)))
            dirs <- k212Dirs(k254NewRoot(paste0("e14-save-", name)), run$fx$omic)
            rec <- k207Save(run$fx$workflow, run$fx$omic, dirs)
            k254Record("E14", name, "in_process", "Save done / failed", sprintf("done=%d;failed=%s", rec$done,
                paste(vapply(rec$failed, function(e) class(e)[[1L]], ""), collapse = ",")))
            k254Record("E14", name, "in_process", "Save notes", paste(rec$notes, collapse = " | "))
            k254Record("E14", name, "in_process", "Save integration object", k254RdsSamples(rec))
        })
        # E10: run, Apply(t_hi), Apply with a grouping column absent from the design.
        k254Case("E10", name, "in_process", function() {
            seen <- k194Notes()
            run <- k254InProcess(config, paste0("e10-", name))
            on.exit(k194Close(run$fx), add = TRUE)
            k194Apply(run$fx, run$nd, threshold = t_hi)
            k254Applied("E10", name, "in_process", run$fx, sprintf("Apply(%s)", t_hi))
            k254Record("E10", name, "in_process", "fields after Apply", k254Fields(run$nd))
            seen$notes <- character()
            input <- k194Input(threshold = t_lo)
            input$ruv_grouping_variable <- "kit254_absent"
            testthat::local_mocked_bindings(workflowOmicsManagedActions = function(omic_type) NULL, .package = "MultiScholaR")
            value <- kit172Quiet(shiny::isolate(applyOmicNormCorrelation(input, run$fx$workflow, run$nd, run$fx$addLog, run$fx$config)))
            k254Record("E10", name, "in_process", "failed Apply value", format(value))
            k254RecordRun("E10", name, "in_process", run$fx, run$nd, NULL, replay = FALSE)
            k254Record("E10", name, "in_process", "notes", paste(seen$notes, collapse = " | "))
        })
        # E11: run, Apply(t_hi), Skip with the store commit refused (F105-02's mock).
        k254Case("E11", name, "in_process", function() {
            seen <- k194Notes()
            run <- k254InProcess(config, paste0("e11-", name))
            on.exit(k194Close(run$fx), add = TRUE)
            k194Apply(run$fx, run$nd, threshold = t_hi)
            k254Record("E11", name, "in_process", "fields after Apply", k254Fields(run$nd))
            seen$notes <- character()
            value <- k254SkipRefused(run$fx, run$nd)
            k254Record("E11", name, "in_process", "refused Skip value", format(value))
            k254RecordRun("E11", name, "in_process", run$fx, run$nd, NULL, replay = FALSE)
            k254Record("E11", name, "in_process", "notes", paste(seen$notes, collapse = " | "))
        })
        # E9: two assays, the constructed pairs, Apply(0.5) twice.
        k254Case("E9", name, "in_process", function() {
            seen <- k194Notes()
            fx <- k254Fixture2(file.path(k254NewRoot(paste0("e9-", name)), "project"), config$class)
            on.exit(k194Close(fx), add = TRUE)
            nd <- omicNormReactiveState()
            k194Run(fx, nd, k194Input())
            k254Record("E9", name, "in_process", "run head", k194Head(fx))
            k254E9Mock()
            k194Apply(fx, nd, threshold = 0.5)
            k254Applied("E9", name, "in_process", fx, "first Apply(0.5)")
            k194Apply(fx, nd, threshold = 0.5)
            k254Applied("E9", name, "in_process", fx, "second Apply(0.5)")
            k254Record("E9", name, "in_process", "errors in log", paste(grep("ERROR", fx$log$lines, value = TRUE), collapse = " | "))
            k254Record("E9", name, "in_process", "notes", paste(seen$notes, collapse = " | "))
        })
    }

    # Managed I and II (phase 0 permitting); E12 (RUV skipped only).
    if (!isTRUE(k254Phase0[[name]])) {
        k254Record("I/II/E12", name, "managed", "not run", "phase 0 did not pass for this omic (instrument; VAL-54)")
        next
    }
    k254Case("I", name, "managed", function() {
        seen <- k194Notes()
        k254Managed(config, paste0("mi-", name), function(session, press, probe, fx) {
            press("apply_correlation_filter", list(min_pearson_correlation_threshold = t_hi))
            probe$settled <- k254Settle(session, probe)
            press("skip_correlation_filter")
        }, function(fx, nd, probe) {
            k254ManagedCommon("I", name, fx, probe)
            k254RecordRun("I", name, "managed", fx, nd, k254NewRoot(paste0("mi-replay-", name)))
            k254DaAndSave("I", name, "managed", fx)
        })
    })
    k254Case("II", name, "managed", function() {
        seen <- k194Notes()
        k254Managed(config, paste0("mii-", name), function(session, press, probe, fx) {
            press("apply_correlation_filter", list(min_pearson_correlation_threshold = t_hi))
            probe$settled <- k254Settle(session, probe)
            probe$first <- c(head = fx$owner$current_state, entity = k254Entity(getWorkflowStateReference(fx$owner)))
            press("apply_correlation_filter", list(min_pearson_correlation_threshold = t_lo))
            probe$settled_second <- k254Settle(session, probe)
        }, function(fx, nd, probe) {
            k254ManagedCommon("II", name, fx, probe)
            k254Record("II", name, "managed", "after Apply(t_hi)", paste(probe$first, collapse = " "))
            k254RecordRun("II", name, "managed", fx, nd, k254NewRoot(paste0("mii-replay-", name)))
        })
    })
    if (ruv_skipped) k254Case("E12", name, "managed", function() {
        seen <- k194Notes()
        k254Managed(config, paste0("e12-", name), function(session, press, probe, fx) {
            probe$hold <- TRUE
            press("apply_correlation_filter", list(min_pearson_correlation_threshold = t_hi))
            probe$first_active <- !is.null(shiny::isolate(probe$runner$step$active()))
            probe$inputs_after_first <- length(probe$inputs)
            press("apply_correlation_filter", list(min_pearson_correlation_threshold = t_lo))
            probe$inputs_after_second <- length(probe$inputs)
            probe$hold <- FALSE
            probe$settled <- k254Settle(session, probe)
        }, function(fx, nd, probe) {
            k254ManagedCommon("E12", name, fx, probe)
            k254Record("E12", name, "managed", "first press held on the worker", isTRUE(probe$first_active))
            k254Record("E12", name, "managed", "worker actions after first / second press",
                sprintf("%s / %s", probe$inputs_after_first, probe$inputs_after_second))
            k254Record("E12", name, "managed", "settled", isTRUE(probe$settled))
            k254Applied("E12", name, "managed", fx, "after both presses")
            k254Record("E12", name, "managed", "fields", k254Fields(nd))
        })
    })
}

k254Write()
cat("KIT-254 G0 done\n")
