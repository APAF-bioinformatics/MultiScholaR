# KIT-255 gate G0 (ticket r1.2, section 10): observe, at the base commit, what a lipidomics or metabolomics
# correlation filter leaves when a sample fails in one assay and passes in another, and what each DA route and the
# integration input then do. Observation only: branch kit-255-g0 carries this script and nothing under R/; it is never
# merged. Nothing here judges a case: every observable is recorded; the comparison with the predictions below is made
# by reading the tables.
#
# What KIT-242's C4b (run 45594264982) already showed is cited, not repeated: the kernel's split on an emptied column,
# the lipid declaration-route refusal and the metabolomics declaration-route completion. K1 (C4b's object) is used here
# only as input to o2, o4 and o5, with one line on A2's normalisation (P02).
#
# Objects (ticket section 10):
#   K1  C4b's object (KIT-242 tools/kit242/g0.R at 7ea50bc1: g0Draw(4) with A2's standards of B_6 set missing; builder
#       COPIED below, seed 242004), through QC complete, ITSD (median), log2 (offset 0), cyclic loess, correlation 0.5.
#   K2  kit172Object with LCMS_Pos = kit172Values(255001L)$values, LCMS_Neg = kit172Values(255003L)$values; B_6's
#       observed LCMS_Neg cells replaced by 2^rnorm(k, m, s), m and s the mean and SD of LCMS_Neg's finite log2 values,
#       drawn under set.seed(255002, "Mersenne-Twister", "Inversion", "Rejection"); missing cells stay missing. On the
#       declaration route: QC complete (a store step, values unchanged), then correlation 0.5 (the design state's values).
#   K0  K2 without the replacement (control).
#
# PHASE 0 (lifecycle check of the managed instrument; AUDIT-KIT-254-plan R1-01; ticket stop rule (iv)): per omic, K2
#   on a two-assay owner at QC complete (kit172Owner = lipidDeclarationsOwner / metabDeclarationsOwner, intensity
#   filter 1 / 0.5, QC complete, as KIT-254's k194PreNorm), the tab's run (log2, quantile, RUV skipped) and one Apply at
#   0.5 through KIT-194 F1's worker journey (omicsManagedChainMocks; the real installWorkflowOmicsActions runner; the
#   worker's entry point in this process; the observers as omicsTabServer mounts them; KIT-254 G0's harness, copied).
#   Passes when the run and the Apply complete (head <prefix>_correlation_filtered) and the managed Apply's design
#   sample list equals the in-process Apply's on the same object. The managed harness is used by nothing else here.
#
# Items (lipid and metab unless stated):
#   o1 (K2)      design samples, each assay's sample columns.
#   o2 (K1, K2)  omicDaPreflight (R/mod_omic_da.R:398-407) with da = omicGuiLipidDa() / omicGuiMetabDa(),
#                kit172Contrasts(), settings$formula_string "~ 0 + group", current = list(reference =
#                getWorkflowStateReference(owner, "<prefix>_correlation_filtered")), the omic's binding, running = FALSE,
#                workflowOmicsManagedActions bound to NULL (in process), inside a MockShinySession, showNotification
#                recorded.
#   o3 (K2)      the DA through the worker (omicsWorkerRun(fixture, "<prefix>_da_contrast", ...), as KIT-242 G0 C1).
#   o4 (K1, K2, metab)  runMetabolitesDA and the declaration route (runWorkflowOmicsReferenceDaAnalysis) on the
#                correlation state, and on the same object with B_6's column removed by hand from every assay (a fresh
#                owner holding that object); per assay and contrast `identical`.
#   o5 (K1, K2)  writeOmicsIntegrationInput at the correlation head; the A1 / LCMS_Pos view's feature table and the
#                samples table (omicsIntegrationInputSource, the writer's own source).
#   K0           design samples, assay columns, both DAs (declaration route).
#   Stop rule (i): K2's B_6 best within-group correlation in each assay is recorded on the declaration route and in
#                phase 0's in-process run (pass in LCMS_Pos, fail in LCMS_Neg, else an instrument defect).
#
# PREDICTIONS: written below (k255Predictions) and to predictions.tsv before the package is loaded or any value read.
#
# Usage (run-tests scripts lane, from the package root):  Rscript tools/kit255/g0.R <out dir>

args <- commandArgs(trailingOnly = TRUE)
k255Out <- if (length(args)) args[[1L]] else "kit255-g0"
dir.create(k255Out, recursive = TRUE, showWarnings = FALSE)

k255Predictions <- data.frame(stringsAsFactors = FALSE, rbind(
    c("PHASE0", "K2", "lipid,metab", "managed run and Apply(0.5)", "complete; head <prefix>_correlation_filtered"),
    c("PHASE0", "K2", "lipid,metab", "managed design = in-process design", "TRUE"),
    c("PHASE0", "K2", "lipid,metab", "B_6 (after log2, quantile)", "out of the design and LCMS_Neg, in LCMS_Pos (both routes)"),
    c("STOP-i", "K2", "lipid,metab", "B_6 best within-group correlation", "LCMS_Pos >= 0.5 (passes); LCMS_Neg < 0.5 (fails); declaration route and phase 0 in process"),
    c("o1", "K2", "lipid,metab", "design samples", "11; B_6 absent"),
    c("o1", "K2", "lipid,metab", "LCMS_Pos sample columns", "12; B_6 present"),
    c("o1", "K2", "lipid,metab", "LCMS_Neg sample columns", "11; B_6 absent"),
    c("K1", "K1", "lipid,metab", "A2 between-sample warning", "present: 'Assay 'A2': At least one sample column contains only NA values. Skipping normalization.'"),
    c("K1", "K1", "lipid,metab", "A2 normalised equals its log2 input", "TRUE"),
    c("K1", "K1", "lipid,metab", "correlation state", "design 11, B_6 absent; A1 12 columns with B_6; A2 11 columns without (as C4b)"),
    c("o2", "K1", "lipid,metab", "omicDaPreflight value", "status=error; stage=preflight; errorMessage=Assay 'A1' has sample column(s) missing from design: B_6"),
    c("o2", "K2", "lipid,metab", "omicDaPreflight value", "status=error; stage=preflight; errorMessage=Assay 'LCMS_Pos' has sample column(s) missing from design: B_6"),
    c("o2", "K1,K2", "lipid,metab", "notification", "one [error] notification with the errorMessage"),
    c("o3", "K2", "lipid", "worker DA handoff", "carries 'Assay 'LCMS_Pos' has sample column(s) missing from design: B_6'; no assay table"),
    c("o3", "K2", "metab", "worker DA handoff", "complete; tables for LCMS_Pos and LCMS_Neg; no error"),
    c("o4", "K1,K2", "metab", "runMetabolitesDA: correlation state vs B_6 removed by hand from every assay", "identical per assay and contrast (PREDICTED BY READING: both select columns by design name)"),
    c("o4", "K1,K2", "metab", "declaration DA: correlation state vs B_6 removed by hand from every assay", "identical per assay and contrast (PREDICTED BY READING)"),
    c("o5", "K1,K2", "lipid,metab", "writer", "writes (TRUE)"),
    c("o5", "K1,K2", "lipid,metab", "A1 / LCMS_Pos view feature table has a column B_6", "TRUE"),
    c("o5", "K1,K2", "lipid,metab", "samples table rows", "11"),
    c("o5", "K1,K2", "lipid,metab", "every matrix's sample columns", "11"),
    c("K0", "K0", "lipid,metab", "correlation state", "design 12; both assays 12 columns (nothing removed)"),
    c("K0", "K0", "lipid,metab", "declaration DA", "completes; tables for LCMS_Pos and LCMS_Neg")))
names(k255Predictions) <- c("item", "object", "class", "observable", "predicted")
utils::write.table(k255Predictions, file.path(k255Out, "predictions.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

k255Git <- function(...) suppressWarnings(system2("git", c(...), stdout = TRUE, stderr = TRUE))
k255Provenance <- c(sprintf("head\t%s", k255Git("rev-parse", "HEAD")[[1L]]), "staging_base\t7327607fbed47c01a45bda9998cfb17c4999455a",
    sprintf("R\t%s", R.version.string), sprintf("started\t%s", format(Sys.time(), tz = "UTC")))

suppressMessages(pkgload::load_all(".", quiet = TRUE, export_all = TRUE, helpers = TRUE))
library(testthat)

# ---- recording ----

k255Rows <- new.env(parent = emptyenv())
k255Rows$obs <- list()
k255Rows$cond <- list()
k255Text <- function(x) gsub("[\t\r\n]+", " ", gsub("[0-9a-f]{64}", "<id>", x))
k255Record <- function(item, object, class, route, observable, value) {
    value <- k255Text(paste(as.character(unlist(value)), collapse = ","))
    if (!length(value) || !nzchar(value)) value <- "<none>"
    k255Rows$obs[[length(k255Rows$obs) + 1L]] <- data.frame(item = item, object = object, class = class, route = route,
        observable = observable, value = value, stringsAsFactors = FALSE)
    cat(sprintf("KIT-255 G0\t%s\t%s\t%s\t%s\t%s\t%s\n", item, object, class, route, observable, value))
}
k255Cond <- function(object, class, step, kind, text) {
    k255Rows$cond[[length(k255Rows$cond) + 1L]] <- data.frame(object = object, class = class, step = step, kind = kind,
        text = k255Text(text), stringsAsFactors = FALSE)
}
k255Write <- function() {
    if (length(k255Rows$obs)) utils::write.table(do.call(rbind, k255Rows$obs), file.path(k255Out, "observations.tsv"), sep = "\t",
        quote = FALSE, row.names = FALSE)
    if (length(k255Rows$cond)) utils::write.table(do.call(rbind, k255Rows$cond), file.path(k255Out, "conditions.tsv"), sep = "\t",
        quote = FALSE, row.names = FALSE)
}
k255Error <- function(e) paste0("error [", paste(class(e), collapse = "/"), "] ", conditionMessage(e))

#' Evaluate `expr`, keeping its warnings and an error (returned as the value), the messages muffled.
k255Capture <- function(expr) {
    warns <- character()
    value <- withCallingHandlers(tryCatch(expr, error = function(e) e),
        warning = function(w) {
            warns <<- c(warns, conditionMessage(w))
            invokeRestart("muffleWarning")
        },
        message = function(m) invokeRestart("muffleMessage"))
    list(value = value, failed = inherits(value, "error"), warnings = warns)
}

k255Case <- function(item, object, class, route, code) {
    started <- proc.time()[["elapsed"]]
    result <- tryCatch(code(), error = function(e) {
        k255Record(item, object, class, route, "INSTRUMENT ERROR", k255Error(e))
        NULL
    })
    k255Record(item, object, class, route, "seconds", round(proc.time()[["elapsed"]] - started, 1))
    k255Write()
    result
}

k255Root <- tempfile("kit255-g0-")
dir.create(k255Root)
k255Count <- new.env(parent = emptyenv())
k255Count$n <- 0L
k255NewRoot <- function(label) {
    k255Count$n <- k255Count$n + 1L
    file.path(k255Root, sprintf("%03d-%s", k255Count$n, label))
}

RUNS <- kit172Design()$Run
S <- "B_6"
DA_SETTINGS <- list(formula_string = "~ 0 + group", da_q_val_thresh = 0.05, treat_lfc_cutoff = 0)
k255Close <- function(owner) try(closeDIAStateAdapter(owner), silent = TRUE)

# ---- K1: copied unchanged from tools/kit242/g0.R at 7ea50bc1 (lines 25, 91-134; g0Matrix 142-152, g0Equal 155-160) ----

SEED <- function(k) set.seed(242000L + k, kind = "Mersenne-Twister", normal.kind = "Inversion", sample.kind = "Rejection")

g0Draw <- function(k, n = c(A1 = 40L, A2 = 40L), n_is = c(A1 = 5L, A2 = 5L)) {
  SEED(k)
  loading <- stats::rnorm(length(RUNS), 0, 0.5)
  out <- list()
  for (a in names(n)) {
    level <- stats::rnorm(n[[a]], 20, 1.5)
    noise <- matrix(stats::rnorm(n[[a]] * length(RUNS), 0, 0.25), n[[a]])
    y <- outer(level, rep(1, length(RUNS))) + outer(rep(1, n[[a]]), loading) + noise
    ids <- sprintf("%s_F%03d", a, seq_len(n[[a]]))
    if (n_is[[a]] > 0L) {
      s <- matrix(16 + stats::rnorm(n_is[[a]] * length(RUNS), 0, 0.1), n_is[[a]]) + outer(rep(1, n_is[[a]]), loading)
      y <- rbind(y, s)
      ids <- c(ids, sprintf("IS_%s_%d", a, seq_len(n_is[[a]])))
    }
    m <- 2^y
    dimnames(m) <- list(ids, RUNS)
    out[[a]] <- m
  }
  out
}

g0Table <- function(m, class, annotation = NULL) {
  ids <- rownames(m)
  ann <- if (is.null(annotation)) ifelse(startsWith(ids, "IS_"), ids, "PC") else annotation
  table <- if (identical(class, "lipid")) data.frame(LipidName = ids, LipidClass = ann, stringsAsFactors = FALSE) else
    data.frame(Metabolite = ids, Ontology = ann, stringsAsFactors = FALSE)
  cbind(table, stats::setNames(as.data.frame(unname(m)), colnames(m)))
}

g0Object <- function(values, class, annotations = list()) {
  assays <- lapply(stats::setNames(names(values), names(values)), function(a) g0Table(values[[a]], class, annotations[[a]]))
  design <- kit172Design()
  args <- list(deAnalysisParameters = list(formula_string = "~ 0 + group", de_q_val_thresh = 0.05, treat_lfc_cutoff = 0))
  if (identical(class, "lipid")) {
    methods::new("LipidomicsAssayData", lipid_data = assays, design_matrix = design, lipid_id_column = "LipidName",
      annotation_id_column = "LipidClass", sample_id = "Run", group_id = "group", technical_replicate_id = "tech_rep_group",
      database_identifier_type = "Unknown", internal_standard_regex = "^IS_", args = args)
  } else {
    methods::new("MetaboliteAssayData", metabolite_data = assays, design_matrix = design, metabolite_id_column = "Metabolite",
      annotation_id_column = "Ontology", sample_id = "Run", group_id = "group", technical_replicate_id = "tech_rep_group",
      database_identifier_type = "Unknown", internal_standard_regex = "^IS_", args = args)
  }
}

g0Matrix <- function(object, binding, assay) {
  if (is.null(object)) return(NULL)
  table <- methods::slot(object, binding$data_slot)[[assay]]
  if (is.null(table)) return(NULL)
  table <- as.data.frame(table, check.names = FALSE, stringsAsFactors = FALSE)
  runs <- intersect(RUNS, names(table))
  m <- as.matrix(table[, runs, drop = FALSE])
  storage.mode(m) <- "double"
  rownames(m) <- as.character(table[[methods::slot(object, binding$id_slot)]])
  m[order(rownames(m)), , drop = FALSE]
}

g0Equal <- function(a, b) {
  if (is.null(a) || is.null(b)) return(NA)
  rows <- intersect(rownames(a), rownames(b))
  cols <- intersect(colnames(a), colnames(b))
  isTRUE(all.equal(unname(a[rows, cols, drop = FALSE]), unname(b[rows, cols, drop = FALSE]), tolerance = 0))
}

#' K1: C4b's object (KIT-242 caseC4b, lines 470-472).
k255K1Object <- function(class) {
    v <- g0Draw(4L)
    v$A2[startsWith(rownames(v$A2), "IS_"), S] <- NA_real_
    g0Object(v, class)
}

# ---- K2, K0 ----

#' K2's values (ticket section 10); `replace = FALSE` gives K0.
k255K2Values <- function(replace = TRUE) {
    pos <- kit172Values(255001L)$values
    neg <- kit172Values(255003L)$values
    if (replace) {
        logs <- log2(neg)
        logs <- logs[is.finite(logs)]
        m <- mean(logs)
        s <- stats::sd(logs)
        set.seed(255002L, kind = "Mersenne-Twister", normal.kind = "Inversion", sample.kind = "Rejection")
        observed <- which(!is.na(neg[, S]))
        neg[observed, S] <- 2^stats::rnorm(length(observed), m, s)
    }
    list(LCMS_Pos = pos, LCMS_Neg = neg)
}

# ---- reading ----

k255State <- function(fx, binding, state) readWorkflowOmicsState(binding, fx$owner, state)
k255DesignRuns <- function(object) as.character(as.data.frame(object@design_matrix)[[object@sample_id]])

#' o1's observables on a state object: the design's samples and each assay's sample columns.
k255Layout <- function(item, object_name, class, route, object, binding) {
    if (is.null(object)) {
        k255Record(item, object_name, class, route, "state", "absent")
        return(invisible(NULL))
    }
    design <- k255DesignRuns(object)
    k255Record(item, object_name, class, route, "design samples", sprintf("%d; %s %s", length(design), S,
        if (S %in% design) "present" else "absent"))
    for (assay in names(methods::slot(object, binding$data_slot))) {
        cols <- intersect(RUNS, names(as.data.frame(methods::slot(object, binding$data_slot)[[assay]], check.names = FALSE)))
        k255Record(item, object_name, class, route, paste(assay, "sample columns"), sprintf("%d; %s %s", length(cols), S,
            if (S %in% cols) "present" else "absent"))
    }
    invisible(NULL)
}

#' Per assay, the best within-group correlation of `sample` in the pairs tables (stop rule (i)).
k255Best <- function(tables, sample = S) {
    vapply(names(tables), function(assay) {
        t <- as.data.frame(tables[[assay]])
        x <- grep("[.]x$", names(t), value = TRUE)[[1L]]
        y <- grep("[.]y$", names(t), value = TRUE)[[1L]]
        r <- t$pearson_correlation[as.character(t[[x]]) == sample | as.character(t[[y]]) == sample]
        r <- r[is.finite(r)]
        if (length(r)) max(r) else NA_real_
    }, numeric(1))
}
k255RecordBest <- function(item, object_name, class, route, tables) {
    best <- k255Best(tables)
    for (assay in names(best)) {
        k255Record(item, object_name, class, route, sprintf("%s best within-group correlation, %s", S, assay),
            sprintf("%s (%s at 0.5)", if (is.na(best[[assay]])) "NA" else sprintf("%.6f", best[[assay]]),
                if (isTRUE(best[[assay]] >= 0.5)) "passes" else "fails"))
    }
    ok <- isTRUE(best[["LCMS_Pos"]] >= 0.5) && !isTRUE(best[["LCMS_Neg"]] >= 0.5)
    k255Record(item, object_name, class, route, "stop rule (i): fails in exactly LCMS_Neg", ok)
    invisible(ok)
}

#' Find every element named `name` in a nested list (KIT-242 G0's findNamed).
k255FindNamed <- function(x, name) {
    found <- list()
    walk <- function(y) {
        if (!is.list(y)) return(invisible())
        if (!is.null(names(y)) && name %in% names(y)) found[[length(found) + 1L]] <<- y[[name]]
        for (z in y) walk(z)
    }
    walk(x)
    found
}

#' The declaration DA on `source`: per contrast the assays with a table, their rows and the errors (KIT-242's g0Da).
k255Da <- function(item, object_name, class, route, fx, binding, source) {
    da <- k255Capture(runWorkflowOmicsReferenceDaAnalysis(binding, fx$owner, source, DA_SETTINGS, kit172Contrasts()))
    if (da$failed) {
        k255Record(item, object_name, class, route, "DA", paste("failed:", conditionMessage(da$value)))
        return(invisible(NULL))
    }
    manifest <- da$value$manifest
    for (contrast in names(manifest$results)) {
        entry <- manifest$results[[contrast]]
        k255Record(item, object_name, class, route, paste("DA", contrast, "assays with a table"), names(entry$assays))
        k255Record(item, object_name, class, route, paste("DA", contrast, "rows per assay"),
            vapply(entry$assays, function(e) as.integer(e$rows), 1L))
        k255Record(item, object_name, class, route, paste("DA", contrast, "errors"), if (length(entry$errors)) entry$errors else "none")
    }
    invisible(da$value)
}

# ---- the declaration-route objects to the correlation state ----

#' A fresh owner holding `object`, taken to the correlation state at 0.5. K1: QC complete, ITSD (median), log2 (0),
#' cyclic loess, correlation (C4b's steps). K2, K0: QC complete, correlation.
k255Build <- function(object_name, class, object) {
    binding <- kit172Binding(class)
    prefix <- binding$prefix
    fx <- kit172Quiet(kit172Owner(file.path(k255NewRoot(paste0(tolower(object_name), "-", class)), "project"), object, class))
    owner <- fx$owner
    config <- fx$workflow$config_list
    steps <- list(qc_finalize = function() commitWorkflowOmicsStoreStep(owner, paste0(prefix, "_qc_finalize"),
        paste0(prefix, "_qc_complete"), config, "QC complete"))
    if (identical(object_name, "K1")) {
        steps <- c(steps, list(
            itsd = function() runWorkflowOmicsReferenceItsdNormalisation(binding, owner, "median", NULL, config),
            log2 = function() runWorkflowOmicsReferenceLog2(binding, owner, 0, config),
            between_sample = function() runWorkflowOmicsReferenceBetweenSampleNormalisation(binding, owner, "cyclicloess", config)))
    }
    steps <- c(steps, list(correlation = function() runWorkflowOmicsReferenceCorrelation(binding, owner, "group", 0.5, config)))
    out <- list(fx = fx, binding = binding, prefix = prefix, ok = TRUE)
    for (step in names(steps)) {
        run <- k255Capture(steps[[step]]())
        for (w in run$warnings) k255Cond(object_name, class, step, "warning", w)
        if (run$failed) {
            k255Cond(object_name, class, step, "error", conditionMessage(run$value))
            k255Record("BUILD", object_name, class, "declaration", paste(step, "result"), paste("failed:", conditionMessage(run$value)))
            out$ok <- FALSE
            break
        }
        out[[step]] <- run
    }
    k255Record("BUILD", object_name, class, "declaration", "head", owner$current_state)
    out$ok <- out$ok && identical(owner$current_state, paste0(prefix, "_correlation_filtered"))
    out
}

# ---- the managed harness: copied from KIT-254 G0 (docs/qa/workflow-execution/kit/kit-254/g0-base-values.R at
# 671aba2c, lines 104-175 and 454-517; k194* copied there unchanged from test-omics-norm-state-guards.R at e0747fc1),
# with the fixture's values given (two assays) and the server's assay names set to both assays ----

k194Omic <- function(class) c(lipid = "lipidomics", metab = "metabolomics")[[class]]

k255PreNorm <- function(root, values, class) {
    fixture <- kit172Owner(root, kit172Object(values, class), class)
    binding <- kit172Binding(class)
    kit172Quiet({
        runWorkflowOmicsReferenceIntensityFilter(binding, fixture$owner, 1, 0.5, fixture$workflow$config_list)
        commitWorkflowOmicsStoreStep(fixture$owner, paste0(binding$prefix, "_qc_finalize"), paste0(binding$prefix, "_qc_complete"),
            fixture$workflow$config_list, "QC complete")
    })
    fixture
}

k255Fixture <- function(root, values, class) {
    fixture <- k255PreNorm(root, values, class)
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

k255Server <- function(fx, probe) {
    force(fx); force(probe)
    function(input, output, session) {
        workflow <- fx$workflow
        runner <- installWorkflowOmicsActions(session, workflow, fx$omic, "kit255", fx$paths)
        runner$progress <- function(message) list(set = function(...) invisible(NULL), inc = function(...) invisible(NULL),
            close = function(...) invisible(NULL))
        probe$runner <- runner
        nd <- omicNormReactiveState()
        nd$assay_names <- c("LCMS_Pos", "LCMS_Neg")
        probe$nd <- nd
        aesthetics <- function(...) list(color_var = "group", shape_var = NULL)
        config <- omicGuiConfig(fx$omic)
        shiny::observeEvent(input$run_normalization, runOmicNormPipeline(input, workflow, fx$paths, nd, fx$addLog, aesthetics, config))
        shiny::observeEvent(input$apply_correlation_filter, applyOmicNormCorrelation(input, workflow, nd, fx$addLog, config))
    }
}

k255Settle <- function(session, probe, rounds = 40L) {
    for (i in seq_len(rounds)) {
        if (is.null(shiny::isolate(probe$runner$step$active()))) break
        session$elapse(100)
    }
    is.null(shiny::isolate(probe$runner$step$active()))
}

#' Phase 0's managed arm: K2's owner at QC complete under the worker-journey mocks; run, then Apply(0.5). After the
#' server ends, `after(fx, probe)` reads the owner (still open; no reactive value is read after the server).
k255Managed <- function(class, label, values, after) {
    fx <- k255Fixture(file.path(k255NewRoot(label), "project"), values, class)
    on.exit(k255Close(fx$owner), add = TRUE)
    probe <- new.env()
    omicsManagedChainMocks(fx, probe)
    suppressWarnings(suppressMessages(utils::capture.output(shiny::testServer(k255Server(fx, probe), {
        session$userData$probe_owner <- fx$owner
        do.call(session$setInputs, k194Input())
        session$setInputs(run_normalization = 1L)
        session$flushReact()
        probe$settled_run <- k255Settle(session, probe)
        probe$ran <- isTRUE(fx$owner$current_state %in% paste0(fx$prefix, c("_normalized", "_ruv_corrected")))
        session$setInputs(min_pearson_correlation_threshold = 0.5)
        session$setInputs(apply_correlation_filter = 1L)
        session$flushReact()
        probe$settled_apply <- k255Settle(session, probe)
    }))))
    after(fx, probe)
}

k255Runs <- function(fx) {
    meta <- readWorkflowOmicsMetadata(fx$binding, fx$owner, NULL)
    if (is.null(meta)) return("<no state>")
    sort(as.character(as.data.frame(meta$design_matrix)[["Run"]]))
}

# ---- PHASE 0 ----

k255Phase0 <- list()
for (class in c("lipid", "metab")) {
    expected <- k255Case("PHASE0", "K2", class, "in_process", function() {
        seen <- k194Notes()
        fx <- k255Fixture(file.path(k255NewRoot(paste0("p0-ip-", class)), "project"), k255K2Values(), class)
        on.exit(k255Close(fx$owner), add = TRUE)
        nd <- omicNormReactiveState()
        k194Run(fx, nd)
        k255Record("PHASE0", "K2", class, "in_process", "run head", fx$owner$current_state)
        pairs <- NULL
        trace <- get("runWorkflowOmicsReferenceCorrelation", envir = asNamespace("MultiScholaR"))
        testthat::local_mocked_bindings(runWorkflowOmicsReferenceCorrelation = function(...) {
            result <- trace(...)
            pairs <<- result$correlations
            result
        }, .package = "MultiScholaR")
        k194Apply(fx, nd, threshold = 0.5)
        k255Record("PHASE0", "K2", class, "in_process", "Apply head", fx$owner$current_state)
        # Recorders only: an error here is recorded and does not change the phase 0 comparison.
        tryCatch({
            if (!is.null(pairs)) k255RecordBest("STOP-i", "K2", class, "phase 0 in process (log2, quantile)", pairs) else
                k255Record("STOP-i", "K2", class, "phase 0 in process (log2, quantile)", "pairs", "<not captured>")
            k255Layout("PHASE0", "K2", class, "in_process", k255State(fx, fx$binding, NULL), fx$binding)
        }, error = function(e) k255Record("PHASE0", "K2", class, "in_process", "RECORDER ERROR", k255Error(e)))
        k255Record("PHASE0", "K2", class, "in_process", "notes", paste(seen$notes, collapse = " | "))
        if (identical(fx$owner$current_state, paste0(fx$prefix, "_correlation_filtered"))) k255Runs(fx) else NULL
    })
    managed <- k255Case("PHASE0", "K2", class, "managed", function() {
        seen <- k194Notes()
        k255Managed(class, paste0("p0-managed-", class), k255K2Values(), function(fx, probe) {
            k255Record("PHASE0", "K2", class, "managed", "run completed", isTRUE(probe$ran) && isTRUE(probe$settled_run))
            k255Record("PHASE0", "K2", class, "managed", "Apply settled", isTRUE(probe$settled_apply))
            k255Record("PHASE0", "K2", class, "managed", "operations", paste(probe$operations, collapse = ","))
            k255Record("PHASE0", "K2", class, "managed", "statuses", paste(probe$statuses, collapse = ","))
            k255Record("PHASE0", "K2", class, "managed", "head", fx$owner$current_state)
            k255Record("PHASE0", "K2", class, "managed", "errors in log", paste(grep("ERROR", fx$log$lines, value = TRUE), collapse = " | "))
            k255Record("PHASE0", "K2", class, "managed", "notes", paste(seen$notes, collapse = " | "))
            runs <- k255Runs(fx)
            tryCatch(k255Layout("PHASE0", "K2", class, "managed", k255State(fx, fx$binding, NULL), fx$binding),
                error = function(e) k255Record("PHASE0", "K2", class, "managed", "RECORDER ERROR", k255Error(e)))
            ok <- isTRUE(probe$ran) && isTRUE(probe$settled_run) && isTRUE(probe$settled_apply) &&
                identical(fx$owner$current_state, paste0(fx$prefix, "_correlation_filtered")) && !is.null(expected) &&
                identical(runs, expected)
            k255Record("PHASE0", "K2", class, "managed", "managed design = in-process design", identical(runs, expected))
            ok
        })
    })
    k255Phase0[[class]] <- isTRUE(managed)
    k255Record("PHASE0", "K2", class, "managed", "PHASE 0 PASSED", k255Phase0[[class]])
}
k255Write()

# ---- the declaration-route objects ----

k255Objects <- list()
for (class in c("lipid", "metab")) {
    binding <- kit172Binding(class)
    for (name in c("K1", "K2", "K0")) {
        k255Objects[[paste(name, class)]] <- k255Case("BUILD", name, class, "declaration", function() {
            object <- switch(name, K1 = k255K1Object(class), K2 = kit172Object(k255K2Values(), class),
                K0 = kit172Object(k255K2Values(replace = FALSE), class))
            built <- k255Build(name, class, object)
            corr <- built$correlation$value
            if (!is.null(corr$correlations)) {
                if (identical(name, "K2")) k255RecordBest("STOP-i", "K2", class, "declaration (design-state values)", corr$correlations)
                if (!identical(name, "K2")) {
                    best <- k255Best(corr$correlations)
                    k255Record("BUILD", name, class, "declaration", "B_6 best within-group correlation per assay",
                        paste(names(best), ifelse(is.na(best), "NA", sprintf("%.6f", best)), sep = "=", collapse = "; "))
                }
            }
            if (identical(name, "K1")) {
                warned <- any(grepl("Assay 'A2'.*Skipping normalization", built$between_sample$warnings))
                k255Record("K1", "K1", class, "declaration", "A2 between-sample warning",
                    if (warned) grep("Assay 'A2'", built$between_sample$warnings, value = TRUE) else "absent")
                log2 <- k255State(built$fx, binding, paste0(binding$prefix, "_log2"))
                norm <- k255State(built$fx, binding, paste0(binding$prefix, "_normalized"))
                feat <- function(m) if (is.null(m)) NULL else m[!startsWith(rownames(m), "IS_"), , drop = FALSE]
                k255Record("K1", "K1", class, "declaration", "A2 normalised equals its log2 input (non-IS rows; tol 0)",
                    g0Equal(feat(g0Matrix(log2, binding, "A2")), feat(g0Matrix(norm, binding, "A2"))))
                k255Record("K1", "K1", class, "declaration", "A1 normalised equals its log2 input (non-IS rows; tol 0)",
                    g0Equal(feat(g0Matrix(log2, binding, "A1")), feat(g0Matrix(norm, binding, "A1"))))
            }
            item <- switch(name, K1 = "K1", K2 = "o1", K0 = "K0")
            if (built$ok) k255Layout(item, name, class, "declaration", k255State(built$fx, binding, NULL), binding)
            built
        })
    }
}

# ---- o2: the tab's in-process preflight through the controller ----

k255Preflight <- function(built, class) {
    binding <- built$binding
    da <- if (identical(class, "lipid")) omicGuiLipidDa() else omicGuiMetabDa()
    seen <- k194Notes()
    testthat::local_mocked_bindings(workflowOmicsManagedActions = function(omic_type) NULL, .package = "MultiScholaR")
    session <- shiny::MockShinySession$new()
    current <- list(reference = getWorkflowStateReference(built$fx$owner, paste0(built$prefix, "_correlation_filtered")))
    value <- shiny::withReactiveDomain(session, omicDaPreflight(built$fx$workflow, current, kit172Contrasts(),
        list(formula_string = "~ 0 + group"), da, binding, running = FALSE))
    list(value = value, notes = seen$notes)
}

for (class in c("lipid", "metab")) {
    for (name in c("K1", "K2")) {
        built <- k255Objects[[paste(name, class)]]
        if (is.null(built) || !isTRUE(built$ok)) {
            k255Record("o2", name, class, "tab in process", "not run", "the object did not reach the correlation state")
            next
        }
        k255Case("o2", name, class, "tab in process", function() {
            out <- k255Preflight(built, class)
            v <- out$value
            k255Record("o2", name, class, "tab in process", "omicDaPreflight value",
                if (is.null(v)) "NULL (the run may go on)" else sprintf("status=%s; stage=%s; errorMessage=%s", v$status, v$stage, v$errorMessage))
            k255Record("o2", name, class, "tab in process", "notifications", if (length(out$notes)) paste(out$notes, collapse = " | ") else "none")
        })
    }
}

# ---- o3: the DA through the worker (K2) ----

for (class in c("lipid", "metab")) {
    built <- k255Objects[[paste("K2", class)]]
    if (is.null(built) || !isTRUE(built$ok)) {
        k255Record("o3", "K2", class, "worker", "not run", "the object did not reach the correlation state")
        next
    }
    k255Case("o3", "K2", class, "worker", function() {
        testthat::local_mocked_bindings(workflowExecutionBuildIdentity = function() strrep("f", 64), .package = "MultiScholaR")
        step <- omicsDaContrastSteps(built$binding, DA_SETTINGS, kit172Contrasts())[[1L]]
        k255Record("o3", "K2", class, "worker", "source head", built$fx$owner$current_state)
        handoff <- k255Capture(omicsWorkerRun(built$fx, step$operation, step$parameters))
        if (handoff$failed) {
            k255Record("o3", "K2", class, "worker", "handoff", paste("omicsWorkerRun failed:", conditionMessage(handoff$value)))
            return(invisible(NULL))
        }
        h <- handoff$value
        strings <- tryCatch(unlist(rapply(h, function(x) x, classes = "character", how = "unlist"), use.names = FALSE),
            error = function(e) character())
        k255Record("o3", "K2", class, "worker", "handoff status", h$status %||% "<none>")
        k255Record("o3", "K2", class, "worker", "handoff error", h$error %||% "none")
        k255Record("o3", "K2", class, "worker", "handoff texts naming 'missing from design'",
            unique(grep("missing from design", c(strings, h$error %||% character()), value = TRUE)) %||% "none")
        assays <- Filter(function(a) is.list(a) && length(names(a)), k255FindNamed(h, "assays"))
        k255Record("o3", "K2", class, "worker", "assay tables (names of each 'assays' list found)",
            if (length(assays)) vapply(assays, function(a) paste(names(a), collapse = "+"), "") else "none")
        rows <- unlist(lapply(assays, function(a) vapply(a, function(e) if (is.list(e) && !is.null(e$rows)) as.character(e$rows) else "-", "")))
        k255Record("o3", "K2", class, "worker", "rows per assay", if (length(rows)) paste(names(rows), rows, sep = "=") else "none")
        errors <- unlist(k255FindNamed(h, "errors"))
        k255Record("o3", "K2", class, "worker", "errors", if (length(errors)) errors else "none")
    })
}

# ---- o4: metabolomics, the correlation state against the same object with B_6 removed by hand ----

k255Tables <- function(results) results$per_assay_results
k255Compare <- function(object_name, route, a, b) {
    assays <- union(names(a), names(b))
    k255Record("o4", object_name, "metab", route, "assays (correlation state | removed by hand)",
        sprintf("%s | %s", paste(names(a), collapse = "+"), paste(names(b), collapse = "+")))
    for (assay in assays) {
        x <- a[[assay]]
        y <- b[[assay]]
        k255Record("o4", object_name, "metab", route, paste(assay, "rows (correlation state | removed by hand)"),
            sprintf("%s | %s", if (is.null(x)) "absent" else nrow(x), if (is.null(y)) "absent" else nrow(y)))
        k255Record("o4", object_name, "metab", route, paste(assay, "identical"), identical(x, y))
        if (!identical(x, y)) {
            k255Record("o4", object_name, "metab", route, paste(assay, "all.equal"),
                paste(utils::head(as.character(all.equal(x, y)), 5L), collapse = " | "))
        }
    }
    all(vapply(assays, function(assay) identical(a[[assay]], b[[assay]]), logical(1)))
}

for (name in c("K1", "K2")) {
    built <- k255Objects[[paste(name, "metab")]]
    if (is.null(built) || !isTRUE(built$ok)) {
        k255Record("o4", name, "metab", "-", "not run", "the object did not reach the correlation state")
        next
    }
    k255Case("o4", name, "metab", "both", function() {
        binding <- built$binding
        reference <- getWorkflowStateReference(built$fx$owner, "metab_correlation_filtered")
        corr <- readWorkflowOmicsState(binding, built$fx$owner, reference)
        removed <- corr
        removed@metabolite_data <- lapply(corr@metabolite_data, function(t) t[, setdiff(names(t), S), drop = FALSE])
        k255Record("o4", name, "metab", "both", "B_6 columns (correlation state | removed by hand)", sprintf("%s | %s",
            paste(names(corr@metabolite_data)[vapply(corr@metabolite_data, function(t) S %in% names(t), logical(1))], collapse = "+"),
            paste(names(removed@metabolite_data)[vapply(removed@metabolite_data, function(t) S %in% names(t), logical(1))], collapse = "+")))
        # runMetabolitesDA (exported).
        a <- k255Capture(runMetabolitesDA(corr, kit172Contrasts(), "~ 0 + group", 0.05, 0))
        b <- k255Capture(runMetabolitesDA(removed, kit172Contrasts(), "~ 0 + group", 0.05, 0))
        if (a$failed || b$failed) {
            k255Record("o4", name, "metab", "runMetabolitesDA", "failed", sprintf("%s | %s",
                if (a$failed) conditionMessage(a$value) else "ok", if (b$failed) conditionMessage(b$value) else "ok"))
        } else {
            k255Record("o4", name, "metab", "runMetabolitesDA", "ALL IDENTICAL", k255Compare(name, "runMetabolitesDA",
                k255Tables(a$value), k255Tables(b$value)))
        }
        # The declaration route: the correlation state on its owner; the removed object as a fresh owner's design state.
        da <- k255Capture(runWorkflowOmicsReferenceDaAnalysis(binding, built$fx$owner, reference, DA_SETTINGS, kit172Contrasts()))
        other <- kit172Quiet(kit172Owner(file.path(k255NewRoot(paste0("o4-", tolower(name))), "project"), removed, "metab"))
        on.exit(k255Close(other$owner), add = TRUE)
        db <- k255Capture(runWorkflowOmicsReferenceDaAnalysis(binding, other$owner, getWorkflowStateReference(other$owner),
            DA_SETTINGS, kit172Contrasts()))
        if (da$failed || db$failed) {
            k255Record("o4", name, "metab", "declaration", "failed", sprintf("%s | %s",
                if (da$failed) conditionMessage(da$value) else "ok", if (db$failed) conditionMessage(db$value) else "ok"))
        } else {
            k255Record("o4", name, "metab", "declaration", "ALL IDENTICAL", k255Compare(name, "declaration",
                k255Tables(da$value$results), k255Tables(db$value$results)))
        }
    })
}

# ---- o5: the integration input at the correlation head ----

for (class in c("lipid", "metab")) {
    for (name in c("K1", "K2")) {
        built <- k255Objects[[paste(name, class)]]
        if (is.null(built) || !isTRUE(built$ok)) {
            k255Record("o5", name, class, "integration", "not run", "the object did not reach the correlation state")
            next
        }
        k255Case("o5", name, class, "integration", function() {
            binding <- built$binding
            owner <- built$fx$owner
            dir <- file.path(k255NewRoot(paste0("o5-", tolower(name), "-", class)), "integration", binding$omic_type)
            dir.create(dir, recursive = TRUE, showWarnings = FALSE)
            written <- k255Capture(writeOmicsIntegrationInput(binding, owner, dir))
            k255Record("o5", name, class, "integration", "writer", if (written$failed) paste("failed:", conditionMessage(written$value)) else
                isTRUE(written$value))
            k255Record("o5", name, class, "integration", "files", sort(list.files(dir, recursive = TRUE)))
            source <- withWorkflowStateScope(owner, function() omicsIntegrationInputSource(binding, owner))
            assay <- if (identical(name, "K1")) "A1" else "LCMS_Pos"
            view <- integrationInputViewName(binding$omic_type, assay)
            features <- source$features[[view]]
            k255Record("o5", name, class, "integration", sprintf("%s view feature table has a column %s", assay, S), S %in% names(features))
            k255Record("o5", name, class, "integration", sprintf("%s view feature table columns", assay), names(features))
            k255Record("o5", name, class, "integration", "samples table rows", nrow(source$samples))
            k255Record("o5", name, class, "integration", "samples table has B_6", S %in% as.character(source$samples$Run))
            k255Record("o5", name, class, "integration", "matrix sample columns per view",
                paste(names(source$matrices), vapply(source$matrices, function(m) ncol(m) - 1L, 1L), sep = "=", collapse = "; "))
        })
    }
}

# ---- K0 and the remaining controls: the declaration DA on the correlation state ----

for (class in c("lipid", "metab")) {
    built <- k255Objects[[paste("K0", class)]]
    if (is.null(built) || !isTRUE(built$ok)) {
        k255Record("K0", "K0", class, "declaration", "not run", "the object did not reach the correlation state")
        next
    }
    k255Case("K0", "K0", class, "declaration", function() {
        k255Da("K0", "K0", class, "declaration", built$fx, built$binding, getWorkflowStateReference(built$fx$owner))
    })
}

for (key in names(k255Objects)) if (!is.null(k255Objects[[key]]$fx)) k255Close(k255Objects[[key]]$fx$owner)
k255Write()
writeLines(c(k255Provenance, sprintf("phase0_lipid\t%s", isTRUE(k255Phase0$lipid)), sprintf("phase0_metab\t%s", isTRUE(k255Phase0$metab)),
    sprintf("finished\t%s", format(Sys.time(), tz = "UTC"))), file.path(k255Out, "provenance.txt"))
cat("KIT-255 G0 done\n")
