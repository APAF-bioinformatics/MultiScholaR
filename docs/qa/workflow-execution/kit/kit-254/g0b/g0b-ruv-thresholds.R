# KIT-254 correction pass, AUDIT-KIT-254-impl I01 step 2: lipid with RUV's thresholds re-derived on the merged tree
# (kit-254 at 993c454b = staging c4505c7f, with KIT-224's legacy RUV-III-C feature levels, plus KIT-254), by G0's rule
# unchanged (k254Choose below, copied with the rest of G0's SCAN from g0-base-values.R at 671aba2c). Observation only;
# scripts-only branch kit-254-g0b, never merged. R1c (run 45624689427) showed G0's lipid-with-RUV value (Apply(0.99)
# removes B_4 only) no longer holds on the merged tree.
#
# PREDICTIONS (before any value is read): t_hi and t_lo for lipid with RUV are not predicted; lipid and metab are not
# run (R1c's precondition passed for them). The configuration's Apply(t_hi) and Apply(t_lo) design samples are
# recorded (C3, II0) so the test constants follow from this run alone. If no listed threshold qualifies, stop rule 2
# applies to E7 and the coordinator decides.
#
args <- commandArgs(trailingOnly = TRUE)
k254Out <- if (length(args)) args[[1L]] else "kit254-g0b"
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
utils::write.table(k254Predictions[0, ], file.path(k254Out, "g0b-unused.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

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
    utils::write.table(rows, file.path(k254Out, "g0b-values.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
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
if (length(k254BestPartner)) utils::write.table(do.call(rbind, k254BestPartner), file.path(k254Out, "g0b-best-partner.tsv"),
    sep = "\t", quote = FALSE, row.names = FALSE)
k254Write()


# C3 and II0 at the derived thresholds (fresh owners).
for (name in names(k254Configs)) {
    config <- k254Configs[[name]]
    for (case in list(c("C3", "hi"), c("II0", "lo"))) {
        t <- k254T[[name]][[paste0("t_", case[[2L]])]]
        k254Case(case[[1L]], name, "in_process", function() {
            seen <- k194Notes()
            run <- k254InProcess(config, paste0(tolower(case[[1L]]), "-", name))
            on.exit(k194Close(run$fx), add = TRUE)
            k254Record(case[[1L]], name, "in_process", "run head", k194Head(run$fx))
            k194Apply(run$fx, run$nd, threshold = t)
            k254Record(case[[1L]], name, "in_process", sprintf("Apply(%s) head", t), k194Head(run$fx))
            k254Record(case[[1L]], name, "in_process", sprintf("Apply(%s) removed", t), paste(k254Removed(k254Runs(run$fx)), collapse = ","))
        })
    }
}
k254Write()
cat("KIT-254 G0b done\n")
