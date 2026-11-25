# KIT-255 gate G1 (ticket r2.1 section 10.1): observe, at the base commit, what the objects the new rule makes do in
# every consumer, and what an object with RUV applied does today. Observation only: branch kit-255-g1 carries this
# script and nothing under R/; it is never merged. Nothing here judges a case: every observable is recorded; the
# comparison with the predictions below is made by reading the tables (ticket 10.1's decision rule).
#
# Objects (ticket 10.1; K2 as G0's, KIT-255-g0 tools/kit255/g0.R, copied never sourced):
#   K2          LCMS_Pos = kit172Values(255001L), LCMS_Neg = kit172Values(255003L) (10 % missing cells), B_6's observed
#               LCMS_Neg cells replaced by 2^rnorm(k, m, s) (set.seed(255002); m, s the mean and SD of LCMS_Neg's finite
#               log2 values): B_6 passes in LCMS_Pos and fails in LCMS_Neg.
#   K2-b        K2 with B_6's LCMS_Neg column NA_real_ and B_6 kept in the design (the end state of the rule, by hand).
#   K2-a        K2 with B_6 out of the design and of every assay (the earlier (a) rule's end state, by hand).
#   K3-b        K2-b with every group-B sample (B_1 to B_6) NA_real_ in LCMS_Neg.
#   K2-RUVnew   P-P2-08: the same draw WITHOUT missing cells (kit172Values(.., missing = FALSE)), B_6 replaced in LCMS_Neg
#               over all its rows; QC complete, log2 (offset 0), quantile, the new RUV design (groups_all_features, user,
#               k = 1), correlation 0.5 (declaration route; and the tab in process and managed).
#   K2-RUVleg   the same values; QC complete, log2, quantile, legacy RUV-III-C (manual, k = 1, 10 %), correlation 0.5.
# The RUV objects are built without missing cells because legacy RUV-III-C's search drops an assay and the new design may
# give k = 0 on the 10 % missing cells (KIT-242 C4, C12; plan audit P2-08).
#
# What is hand-built and how (declared before the run): K2-b, K2-a and K3-b are plain S4 objects given to the exported DA
# kernels (runLipidsDA, runMetabolitesDA), to the declaration route on a fresh owner, to the tab's preflight, to the plots
# and, on an owner holding QC complete and norm complete, to the integration writer and the model's alignment. The RUV-
# applied (b) states (P-e, P-f RUV arm) are the uncorrected values with the failing column NA, given to the same exported
# kernel with the factors a real runOmicsRuvFactors call estimated BEFORE the blanking (as the real filter blanks after
# RUV); the declared RUV DA's own wrapper is exercised in P-b on the real state. P-k hands the writer the real RUV-
# applied owner at its ruv_corrected head with the failing column blanked in the final read and in the read of the
# normalised matrix that entered RUV (mocked readWorkflowOmicsState, by state reference).
#
# PREDICTIONS: written below (k255Predictions) and to predictions.tsv before the package is loaded or any value read.
#
# Usage (run-tests scripts lane, from the package root):  Rscript tools/kit255/g1.R <out dir>

args <- commandArgs(trailingOnly = TRUE)
k255Out <- if (length(args)) args[[1L]] else "kit255-g1"
dir.create(k255Out, recursive = TRUE, showWarnings = FALSE)

k255Predictions <- data.frame(stringsAsFactors = FALSE, rbind(
    c("STOP-i", "K2-b,K3-b", "lipid,metab", "construction", "K2-b: B_6 all NA in LCMS_Neg only (LCMS_Pos untouched); K3-b: B_1 to B_6 all NA in LCMS_Neg"),
    c("STOP-i", "K2-RUVnew,K2-RUVleg", "lipid,metab", "B_6 best within-group correlation after RUV, declaration route", "LCMS_Pos >= 0.5 (passes); LCMS_Neg < 0.5 (fails)"),
    c("STOP-i", "K2-RUVnew,K2-RUVleg", "lipid,metab", "RUV applied: k per assay at the ruv_corrected head", "k = 1 in LCMS_Pos and LCMS_Neg (new design: status ok; legacy: applied)"),
    c("P-a", "K2-RUVnew", "lipid,metab", "Apply 0.5, declaration: design samples", "11; B_6 absent"),
    c("P-a", "K2-RUVnew", "lipid,metab", "Apply 0.5, declaration: LCMS_Pos / LCMS_Neg sample columns", "12 with B_6 / 11 without B_6"),
    c("P-a", "K2-RUVnew", "lipid,metab", "args$ruv_uncorrected_assays$LCMS_Neg: B_6 finite values, and copy before = copy after the press", "all 520 finite; the copy identical before and after"),
    c("P-a", "K2-RUVleg", "lipid,metab", "design samples; LCMS_Pos / LCMS_Neg sample columns", "11 without B_6; 12 with B_6 / 11 without"),
    c("P-a", "K2-RUVleg", "lipid,metab", "args$ruv_uncorrected_assays", "absent"),
    c("P-a", "K2-RUVnew", "lipid,metab", "tab in process then managed (new design, quantile): design, columns, copy, managed = in process", "design 11 without B_6; LCMS_Pos 12 with B_6; LCMS_Neg 11; copy's B_6 finite in all rows; managed design = in-process design"),
    c("P-b", "K2-RUVnew", "lipid", "declared RUV DA on the filtered state", "no assay table; errors name 'missing from design: B_6' (LCMS_Pos)"),
    c("P-b", "K2-RUVnew", "metab", "declared RUV DA on the filtered state", "completes; tables for LCMS_Pos and LCMS_Neg; no errors; B_6 is not in the design"),
    c("P-b", "K2-RUVleg", "lipid", "legacy DA on the filtered state", "no assay table; errors or a failure naming 'missing from design: B_6'"),
    c("P-b", "K2-RUVleg", "metab", "legacy DA on the filtered state", "completes; tables for LCMS_Pos and LCMS_Neg"),
    c("P-c", "K2-RUVnew,K2-RUVleg", "lipid,metab", "writeOmicsIntegrationInput on the filtered state", "writes; samples 11; every matrix 11 columns; corrected matrices 11 columns; the LCMS_Pos feature table has a column B_6; the uncorrected matrix of the corrected view has 11 columns without B_6; correction present"),
    c("P-d", "K2-b", "lipid,metab", "exported DA kernel on K2-b and K2-a", "both complete; tables for LCMS_Pos and LCMS_Neg"),
    c("P-d", "K2-b", "lipid,metab", "LCMS_Neg table of K2-b against K2-a: common numeric non-run columns differing above 1e-10", "none (equal to 1e-10)"),
    c("P-d", "K2-b", "lipid,metab", "LCMS_Pos table of K2-b against K2-a", "differs (the 12-sample fit includes B_6)"),
    c("P-d", "K2-b", "lipid,metab", "B_6 in the long tables: per-sample column, finite cells in LCMS_Neg / LCMS_Pos", "if the table carries per-sample columns: B_6 all NA in LCMS_Neg, finite in LCMS_Pos"),
    c("P-d", "K2-b,K2-a", "lipid,metab", "per-feature residual df of lmFit(~ 0 + group) on LCMS_Neg", "K2-b = K2-a feature by feature (df = observed - rank)"),
    c("P-d", "K2-b", "lipid,metab", "features whose group-B observed count falls from >= 3 to < 3 by the blanking", "a count is reported (information)"),
    c("P-d", "K2-b", "lipid,metab", "declaration route DA on a fresh owner holding K2-b", "completes; tables for LCMS_Pos and LCMS_Neg; no errors"),
    c("P-d", "K2-b", "lipid,metab", "tab in-process preflight (omicDaPreflight) on K2-b", "NULL (the run may go on)"),
    c("P-e", "K2-b-RUV", "lipid,metab", "exported kernel with factors on the uncorrected values, B_6 blank in LCMS_Neg, against B_6 removed (same factors)", "both complete; W has a B_6 row in LCMS_Neg; LCMS_Neg equal to 1e-10; LCMS_Pos differs"),
    c("P-f", "K3-b", "lipid,metab", "non-RUV kernel DA with the contrast groupB-groupA", "completes; a LCMS_Neg table with NA statistics for the contrast (no abort of both assays); LCMS_Pos tested"),
    c("P-f", "K3-b-RUV", "lipid,metab", "RUV-factors kernel DA, group B blank in LCMS_Neg", "completes; group B listed as not observed in LCMS_Neg; LCMS_Pos tested"),
    c("P-g", "K2-b", "lipid,metab", "plotPca, plotDensity", "LCMS_Neg shows 11 samples (B_6 omitted), LCMS_Pos 12"),
    c("P-g", "K2-b", "lipid,metab", "plotRle, static volcano, Glimma volcano (LCMS_Neg)", "build without error"),
    c("P-g", "K2-b", "lipid,metab", "DA heatmap LCMS_Neg (all features, row scaling)", "builds; B_6's column drawn as zeros (every cell 0)"),
    c("P-g", "K2-b", "lipid,metab", "DA heatmap Combined", "builds; B_6's LCMS_Neg cells drawn as zeros, its LCMS_Pos cells not"),
    c("P-h", "K2-b", "lipid,metab", "writer on an owner at norm complete", "writes; samples table 12 rows; every matrix 12 columns; LCMS_Neg B_6 numeric all NA; no column B_6 in either feature table"),
    c("P-h", "K2-b", "lipid,metab", "model alignment, whole and partitioned, missing_samples keep_missing and drop_incomplete", "no refusal; every view counts 12 samples measured (B_6 counted in LCMS_Neg); drop_incomplete keeps B_6; the two paths agree"),
    c("P-i", "K2-b", "lipid,metab", "design rows; samples counted per assay; the summary's removed count", "12; 12 and 12; Removed: 0 (arithmetic of the summary on the metadata reader's design counts)"),
    c("P-j", "K2", "lipid,metab", "project saved by the base through the declared journey (QC complete, correlation 0.5), analysis record, replay once", "design 11 without B_6; LCMS_Pos holds B_6; replay mismatches empty; tarball and sha written"),
    c("P-k", "K2-b-RUVnew,K2-b-RUVleg", "lipid,metab", "writer with the failing column blank in the final read and in the normalised matrix that entered RUV", "LCMS_Neg removed_factors status not_determined, 0 features measured in every run; warning of class integration_input_removed_factors_undetermined; input written without the factors; LCMS_Pos determined")))
names(k255Predictions) <- c("item", "object", "class", "observable", "predicted")
utils::write.table(k255Predictions, file.path(k255Out, "predictions.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

k255Git <- function(...) suppressWarnings(system2("git", c(...), stdout = TRUE, stderr = TRUE))
k255Provenance <- c(sprintf("head\t%s", k255Git("rev-parse", "HEAD")[[1L]]), sprintf("R\t%s", R.version.string),
    sprintf("started\t%s", format(Sys.time(), tz = "UTC")))

suppressMessages(pkgload::load_all(".", quiet = TRUE, export_all = TRUE, helpers = TRUE))
library(testthat)

# ---- recording (G0's) ----

k255Rows <- new.env(parent = emptyenv())
k255Rows$obs <- list()
k255Text <- function(x) gsub("[\t\r\n]+", " ", gsub("[0-9a-f]{64}", "<id>", x))
k255Record <- function(item, object, class, route, observable, value) {
    value <- k255Text(paste(as.character(unlist(value)), collapse = ","))
    if (!length(value) || !nzchar(value)) value <- "<none>"
    k255Rows$obs[[length(k255Rows$obs) + 1L]] <- data.frame(item = item, object = object, class = class, route = route,
        observable = observable, value = value, stringsAsFactors = FALSE)
    cat(sprintf("KIT-255 G1\t%s\t%s\t%s\t%s\t%s\t%s\n", item, object, class, route, observable, value))
}
k255Write <- function() {
    if (length(k255Rows$obs)) utils::write.table(do.call(rbind, k255Rows$obs), file.path(k255Out, "observations.tsv"), sep = "\t",
        quote = FALSE, row.names = FALSE)
}
k255Val <- function(x) if (is.null(x)) "NULL" else paste(as.character(unlist(x)), collapse = ",")
k255Error <- function(e) paste0("error [", paste(class(e), collapse = "/"), "] ", conditionMessage(e))

#' Evaluate `expr`, keeping its warnings (text and classes) and an error (returned as the value), the messages muffled.
k255Capture <- function(expr) {
    warns <- character()
    classes <- character()
    value <- withCallingHandlers(tryCatch(expr, error = function(e) e),
        warning = function(w) {
            warns <<- c(warns, conditionMessage(w))
            classes <<- c(classes, class(w))
            invokeRestart("muffleWarning")
        },
        message = function(m) invokeRestart("muffleMessage"))
    list(value = value, failed = inherits(value, "error"), warnings = warns, classes = classes)
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

k255Root <- tempfile("kit255-g1-")
dir.create(k255Root)
k255Count <- new.env(parent = emptyenv())
k255Count$n <- 0L
k255NewRoot <- function(label) {
    k255Count$n <- k255Count$n + 1L
    file.path(k255Root, sprintf("%03d-%s", k255Count$n, label))
}

RUNS <- kit172Design()$Run
S <- "B_6"
GROUP_B <- sprintf("B_%d", 1:6)
DA_SETTINGS <- list(formula_string = "~ 0 + group", da_q_val_thresh = 0.05, treat_lfc_cutoff = 0)
k255Close <- function(owner) try(closeDIAStateAdapter(owner), silent = TRUE)
k255Gate <- function(code) withr::with_options(list(MultiScholaR.ruv_design_panel = TRUE), code)

# ---- values and objects ----

#' K2's values (G0's draw; `missing = FALSE` for the RUV objects, B_6's replacement then over every row).
k255K2Values <- function(missing = TRUE) {
    pos <- kit172Values(255001L, missing = missing)$values
    neg <- kit172Values(255003L, missing = missing)$values
    logs <- log2(neg)
    logs <- logs[is.finite(logs)]
    m <- mean(logs)
    s <- stats::sd(logs)
    set.seed(255002L, kind = "Mersenne-Twister", normal.kind = "Inversion", sample.kind = "Rejection")
    observed <- which(!is.na(neg[, S]))
    neg[observed, S] <- 2^stats::rnorm(length(observed), m, s)
    list(LCMS_Pos = pos, LCMS_Neg = neg)
}

k255Blank <- function(values, assay, runs) {
    values[[assay]][, runs] <- NA_real_
    values
}
k255Without <- function(values, runs) lapply(values, function(m) m[, setdiff(colnames(m), runs), drop = FALSE])
k255Log <- function(values) lapply(values, log2)

#' kit172Object with the design replaced (K2-a has no B_6 in the design).
k255Obj <- function(values, class, design = NULL) {
    object <- kit172Object(values, class)
    if (!is.null(design)) object@design_matrix <- design
    object
}
k255DesignWithout <- function(runs) {
    design <- kit172Design()
    design[!(design$Run %in% runs), , drop = FALSE]
}

k255Slot <- function(class) if (identical(class, "lipid")) "lipid_data" else "metabolite_data"
k255DesignRuns <- function(object) as.character(as.data.frame(object@design_matrix)[[object@sample_id]])

#' The design's samples and each assay's sample columns of an object (G0's o1 observables).
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

#' STOP-i for the hand-built objects: B_6 (and group B) NA exactly where intended.
k255CheckConstruction <- function(item, name, class, values) {
    neg_b6 <- all(is.na(values$LCMS_Neg[, S]))
    pos_b6 <- any(!is.na(values$LCMS_Pos[, S]))
    pos_other <- any(!is.na(values$LCMS_Neg[, setdiff(RUNS, S)]))
    k255Record(item, name, class, "construction", "B_6 all NA in LCMS_Neg; finite cells in LCMS_Pos; other LCMS_Neg samples hold values",
        sprintf("%s; %s; %s", neg_b6, pos_b6, pos_other))
}

# ---- reading ----

k255State <- function(fx, binding, state) readWorkflowOmicsState(binding, fx$owner, state)

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

#' The declaration DA on `source`: per contrast the assays with a table, their rows and the errors (G0's k255Da).
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
        k255Record(item, object_name, class, route, paste("DA", contrast, "errors"), if (length(entry$errors)) unlist(entry$errors) else "none")
    }
    invisible(da$value)
}

#' The exported DA kernel on a plain object (optionally with factors).
k255Kernel <- function(class, object, factors = NULL) {
    run <- if (identical(class, "lipid")) runLipidsDA else runMetabolitesDA
    if (is.null(factors)) k255Capture(run(object, kit172Contrasts(), "~ 0 + group", 0.05, 0)) else
        k255Capture(run(object, kit172Contrasts(), "~ 0 + group", 0.05, 0, ruv_factors = factors))
}

k255RecordKernel <- function(item, object_name, class, route, run) {
    if (run$failed) {
        k255Record(item, object_name, class, route, "kernel DA", paste("failed:", conditionMessage(run$value)))
        return(invisible(FALSE))
    }
    tables <- run$value$per_assay_results
    k255Record(item, object_name, class, route, "kernel DA assays with a table", names(tables))
    k255Record(item, object_name, class, route, "kernel DA rows per assay", paste(names(tables), vapply(tables, nrow, 1L), sep = "=", collapse = "; "))
    if (length(run$warnings)) k255Record(item, object_name, class, route, "kernel DA warnings", unique(run$warnings))
    invisible(TRUE)
}

k255IdColumn <- function(table) intersect(c("lipid_id", "metabolite_id", "feature_id"), names(table))[[1L]]

#' Compare two long tables by their common numeric columns that are not run columns, to 1e-10.
k255CompareTables <- function(item, object_name, class, route, label, x, y) {
    if (is.null(x) || is.null(y)) {
        k255Record(item, object_name, class, route, label, "a table is absent")
        return(invisible(NA))
    }
    x <- as.data.frame(x)
    y <- as.data.frame(y)
    id <- k255IdColumn(x)
    key <- function(t) if ("comparison" %in% names(t)) paste(as.character(t[[id]]), as.character(t$comparison)) else as.character(t[[id]])
    x <- x[order(key(x)), , drop = FALSE]
    y <- y[order(key(y)), , drop = FALSE]
    if (!identical(as.character(x[[id]]), as.character(y[[id]]))) {
        k255Record(item, object_name, class, route, label, sprintf("the tables list different features (%d vs %d rows)", nrow(x), nrow(y)))
        return(invisible(FALSE))
    }
    common <- setdiff(intersect(names(x), names(y)), RUNS)
    numeric <- common[vapply(common, function(n) is.numeric(x[[n]]) && is.numeric(y[[n]]), logical(1))]
    diffs <- vapply(numeric, function(n) {
        both_na <- is.na(x[[n]]) & is.na(y[[n]])
        d <- abs(x[[n]] - y[[n]])[!both_na]
        if (anyNA(d)) Inf else if (length(d)) max(d) else 0
    }, numeric(1))
    k255Record(item, object_name, class, route, paste(label, "numeric columns compared"), length(numeric))
    k255Record(item, object_name, class, route, paste(label, "columns differing above 1e-10"),
        if (any(diffs > 1e-10)) paste(names(diffs)[diffs > 1e-10], signif(diffs[diffs > 1e-10], 3), sep = "=") else "none")
    k255Record(item, object_name, class, route, paste(label, "only in the first / only in the second"),
        sprintf("%s / %s", paste(setdiff(names(x), names(y)), collapse = "+"), paste(setdiff(names(y), names(x)), collapse = "+")))
    invisible(!any(diffs > 1e-10))
}

#' P-d's per-sample columns: whether the long table carries B_6, and its finite cells.
k255RecordBCells <- function(item, object_name, class, route, tables) {
    for (assay in names(tables)) {
        t <- as.data.frame(tables[[assay]])
        k255Record(item, object_name, class, route, paste(assay, "long table: B_6 column; finite cells; rows"),
            if (S %in% names(t)) sprintf("present; %d; %d", sum(is.finite(t[[S]])), nrow(t)) else sprintf("no per-sample column; -; %d", nrow(t)))
    }
}

#' Per-feature residual df of lmFit(~ 0 + group) on an assay of a plain object, as a table of counts.
k255ResidualDf <- function(object, class, assay) {
    table <- as.data.frame(methods::slot(object, k255Slot(class))[[assay]], check.names = FALSE)
    design <- as.data.frame(object@design_matrix)
    runs <- intersect(as.character(design$Run), names(table))
    m <- as.matrix(table[, runs, drop = FALSE])
    storage.mode(m) <- "double"
    mm <- stats::model.matrix(~ 0 + group, design[match(runs, as.character(design$Run)), , drop = FALSE])
    fit <- limma::lmFit(m, mm)
    stats::setNames(fit$df.residual, as.character(table[[1L]]))
}

# ==== object builders on the declaration route (K2-RUVnew, K2-RUVleg) ====

RUV_NEW <- list(ruv_design = "groups_all_features", ruv_grouping_variable = "group", ruv_k = 1L, ruv_k_rule = "user", ruv_k_seed = NULL,
    ruv_control_ids = NULL)
RUV_LEG <- utils::modifyList(kit172LegacyRuv("manual"), list(manual_k = 1L, manual_percentage = 10))

#' A fresh owner holding the object of `values` (no missing cells), taken to the ruv_corrected head: QC complete, log2 (0),
#' quantile, RUV. The correlation step is run separately (`k255Correlate`).
k255BuildRuv <- function(kind, class) {
    binding <- kit172Binding(class)
    prefix <- binding$prefix
    object_name <- paste0("K2-RUV", kind)
    fx <- kit172Quiet(kit172Owner(file.path(k255NewRoot(paste0(tolower(object_name), "-", class)), "project"),
        kit172Object(k255K2Values(missing = FALSE), class), class))
    owner <- fx$owner
    config <- fx$workflow$config_list
    steps <- list(
        qc_finalize = function() commitWorkflowOmicsStoreStep(owner, paste0(prefix, "_qc_finalize"), paste0(prefix, "_qc_complete"), config, "QC complete"),
        log2 = function() runWorkflowOmicsReferenceLog2(binding, owner, 0, config),
        between_sample = function() runWorkflowOmicsReferenceBetweenSampleNormalisation(binding, owner, "quantile", config),
        ruv = function() runWorkflowOmicsReferenceRuv(binding, owner, if (identical(kind, "new")) RUV_NEW else RUV_LEG, config))
    out <- list(fx = fx, binding = binding, prefix = prefix, name = object_name, ok = TRUE)
    for (step in names(steps)) {
        run <- k255Capture(k255Gate(steps[[step]]()))
        if (length(run$warnings)) k255Record("BUILD", object_name, class, "declaration", paste(step, "warnings"), unique(run$warnings))
        if (run$failed) {
            k255Record("BUILD", object_name, class, "declaration", paste(step, "result"), paste("failed:", conditionMessage(run$value)))
            out$ok <- FALSE
            break
        }
        out[[step]] <- run
    }
    k255Record("BUILD", object_name, class, "declaration", "head", owner$current_state)
    out$ok <- out$ok && identical(owner$current_state, paste0(prefix, "_ruv_corrected"))
    if (out$ok) {
        ruv <- out$ruv$value$value
        k255Record("STOP-i", object_name, class, "declaration", "RUV best_k per assay", unlist(ruv$best_k))
        k255Record("STOP-i", object_name, class, "declaration", "RUV applied per assay", unlist(ruv$ruv_applied))
        k255Record("STOP-i", object_name, class, "declaration", "RUV not-applied reasons", unlist(ruv$ruv_not_applied_reason))
    }
    out
}

#' The correlation step (0.5) on a built owner; records B_6's best correlation (stop rule (i)).
k255Correlate <- function(built, class) {
    binding <- built$binding
    run <- k255Capture(k255Gate(runWorkflowOmicsReferenceCorrelation(binding, built$fx$owner, "group", 0.5,
        built$fx$workflow$config_list)))
    if (length(run$warnings)) k255Record("BUILD", built$name, class, "declaration", "correlation warnings", unique(run$warnings))
    if (run$failed) {
        k255Record("BUILD", built$name, class, "declaration", "correlation result", paste("failed:", conditionMessage(run$value)))
        return(FALSE)
    }
    k255RecordBest("STOP-i", built$name, class, "declaration (after RUV)", run$value$correlations)
    k255Record("BUILD", built$name, class, "declaration", "head", built$fx$owner$current_state)
    identical(built$fx$owner$current_state, paste0(built$prefix, "_correlation_filtered"))
}

# ==== integration writer recording ====

#' Write the integration input of `owner` and record what the writer and the source hold; returns the written folder.
k255WriteInput <- function(item, object_name, class, route, binding, owner, view_for_features = "LCMS_Pos") {
    dir <- file.path(k255NewRoot(paste0(item, "-", tolower(gsub("[^A-Za-z0-9]", "", object_name)), "-", class)), "integration", binding$omic_type)
    dir.create(dir, recursive = TRUE, showWarnings = FALSE)
    written <- k255Capture(writeOmicsIntegrationInput(binding, owner, dir))
    k255Record(item, object_name, class, route, "writer", if (written$failed) paste("failed:", conditionMessage(written$value)) else isTRUE(written$value))
    k255Record(item, object_name, class, route, "writer warning classes", if (length(written$classes)) unique(written$classes) else "none")
    k255Record(item, object_name, class, route, "writer warnings", if (length(written$warnings)) unique(written$warnings) else "none")
    if (written$failed) return(NULL)
    k255Record(item, object_name, class, route, "files", sort(list.files(dir, recursive = TRUE)))
    source <- withWorkflowStateScope(owner, function() omicsIntegrationInputSource(binding, owner))
    if (is.null(source)) {
        k255Record(item, object_name, class, route, "source", "NULL")
        return(dir)
    }
    k255Record(item, object_name, class, route, "samples table rows", nrow(source$samples))
    k255Record(item, object_name, class, route, "samples table has B_6", S %in% as.character(source$samples$Run))
    for (view in names(source$matrices)) {
        m <- source$matrices[[view]]
        k255Record(item, object_name, class, route, paste(view, "matrix sample columns; B_6 column"),
            sprintf("%d; %s", ncol(m) - 1L, if (S %in% names(m)) sprintf("present, numeric %s, %d of %d NA", is.numeric(m[[S]]), sum(is.na(m[[S]])), nrow(m)) else "absent"))
        features <- source$features[[view]]
        k255Record(item, object_name, class, route, paste(view, "feature table has a column", S), S %in% names(features))
    }
    for (view in names(source$corrected)) {
        m <- source$corrected[[view]]
        k255Record(item, object_name, class, route, paste(view, "corrected matrix sample columns; B_6 column"),
            sprintf("%d; %s", ncol(m) - 1L, if (S %in% names(m)) sprintf("present, %d of %d NA", sum(is.na(m[[S]])), nrow(m)) else "absent"))
    }
    k255Record(item, object_name, class, route, "manifest correction",
        substr(tryCatch(as.character(jsonlite::toJSON(source$manifest$correction, auto_unbox = TRUE, null = "null")),
            error = function(e) paste("unreadable:", conditionMessage(e))), 1L, 500L))
    manifest <- tryCatch(jsonlite::read_json(file.path(dir, "integration_input.json")), error = function(e) NULL)
    for (entry in manifest$views) {
        removed <- entry$removed_factors
        if (!is.null(removed)) {
            k255Record(item, object_name, class, route, paste(entry$name, "removed_factors"),
                sprintf("status=%s; features=%s; k_declared=%s; k=%s; reason=%s", k255Val(removed$status), k255Val(removed$features),
                    k255Val(removed$k_declared), k255Val(removed$k), k255Val(removed$reason)))
        }
    }
    dir
}

# ==== a staged owner (QC complete, norm complete) holding a plain object ====

k255Staged <- function(class, object, label) {
    binding <- kit172Binding(class)
    fx <- kit172Quiet(kit172Owner(file.path(k255NewRoot(label), "project"), object, class))
    config <- fx$workflow$config_list
    kit172Quiet({
        commitWorkflowOmicsStoreStep(fx$owner, paste0(binding$prefix, "_qc_finalize"), paste0(binding$prefix, "_qc_complete"), config, "QC complete")
        commitWorkflowOmicsStoreStep(fx$owner, paste0(binding$prefix, "_norm_complete"), paste0(binding$prefix, "_norm_complete"), config,
            "Normalization complete")
    })
    list(fx = fx, binding = binding, prefix = binding$prefix, ok = identical(fx$owner$current_state, paste0(binding$prefix, "_norm_complete")))
}

# ==== the tab's preflight (G0's o2, on a head reference) ====

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

k255Preflight <- function(fx, binding, class) {
    da <- if (identical(class, "lipid")) omicGuiLipidDa() else omicGuiMetabDa()
    seen <- k194Notes()
    testthat::local_mocked_bindings(workflowOmicsManagedActions = function(omic_type) NULL, .package = "MultiScholaR")
    session <- shiny::MockShinySession$new()
    current <- list(reference = getWorkflowStateReference(fx$owner))
    value <- shiny::withReactiveDomain(session, omicDaPreflight(fx$workflow, current, kit172Contrasts(),
        list(formula_string = "~ 0 + group"), da, binding, running = FALSE))
    list(value = value, notes = seen$notes)
}

# ==== the managed harness: G0's (KIT-254 G0's), the tab's inputs extended with the new RUV design ====

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

k194InputNew <- function(threshold = 0.5) {
    list(apply_itsd = FALSE, itsd_aggregation = "median", log_offset = 0, norm_method = "quantile", ruv_mode = "skip",
        auto_percentage_min = 1, auto_percentage_max = 20, ruv_grouping_variable = "group", separation_metric = "max_difference",
        k_penalty_weight = 0.5, adaptive_k_penalty = TRUE, ruv_k = 1L, ruv_percentage = 20, min_pearson_correlation_threshold = threshold,
        ruv_design = "groups_all_features", ruv_k_rule = "user", ruv_k_user = 1, ruv_k_seed = RUV_FACTOR_COUNT_SEED)
}

k194Run <- function(fixture, nd, input) {
    testthat::local_mocked_bindings(workflowOmicsManagedActions = function(omic_type) NULL, .package = "MultiScholaR")
    kit172Quiet(shiny::isolate(runOmicNormPipeline(input, fixture$workflow, list(), nd, fixture$addLog,
        function(...) list(color_var = "group", shape_var = NULL), fixture$config)))
}

k194Apply <- function(fixture, nd, input) {
    testthat::local_mocked_bindings(workflowOmicsManagedActions = function(omic_type) NULL, .package = "MultiScholaR")
    kit172Quiet(shiny::isolate(applyOmicNormCorrelation(input, fixture$workflow, nd, fixture$addLog, fixture$config)))
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

#' The managed arm: the owner at QC complete under the worker-journey mocks; run (new RUV design), then Apply(0.5).
k255Managed <- function(class, label, values, after) {
    fx <- k255Fixture(file.path(k255NewRoot(label), "project"), values, class)
    on.exit(k255Close(fx$owner), add = TRUE)
    probe <- new.env()
    omicsManagedChainMocks(fx, probe)
    suppressWarnings(suppressMessages(utils::capture.output(shiny::testServer(k255Server(fx, probe), {
        session$userData$probe_owner <- fx$owner
        do.call(session$setInputs, k194InputNew())
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

#' B_6's finite cells in the uncorrected copy of an object's LCMS_Neg, or "absent".
k255CopyCells <- function(object) {
    copy <- object@args$ruv_uncorrected_assays
    if (is.null(copy)) return("absent")
    sprintf("%d of %d finite", sum(is.finite(as.numeric(as.data.frame(copy$LCMS_Neg)[[S]]))), nrow(copy$LCMS_Neg))
}

# ======================================================================================================================
# THE OBSERVATIONS
# ======================================================================================================================

# ---- construction check of the hand-built objects (stop rule (i), by construction) ----

for (class in c("lipid", "metab")) {
    k255Case("STOP-i", "K2-b,K3-b", class, "construction", function() {
        k255CheckConstruction("STOP-i", "K2-b", class, k255Blank(k255K2Values(), "LCMS_Neg", S))
        k255CheckConstruction("STOP-i", "K3-b", class, k255Blank(k255K2Values(), "LCMS_Neg", GROUP_B))
    })
}

# ---- P-d, P-f (non-RUV), P-g, P-i: the hand-built K2-b, K2-a, K3-b ----

k255Plain <- list()
for (class in c("lipid", "metab")) {
    binding <- kit172Binding(class)
    values <- k255K2Values()
    kb <- k255Obj(k255Blank(values, "LCMS_Neg", S), class)
    ka <- k255Obj(k255Without(values, S), class, k255DesignWithout(S))
    k3 <- k255Obj(k255Blank(values, "LCMS_Neg", GROUP_B), class)

    # P-d: the exported kernel on K2-b and K2-a; the tables; B_6's cells; the residual df; the features below the minimum.
    k255Case("P-d", "K2-b", class, "kernel", function() {
        db <- k255Kernel(class, kb)
        da <- k255Kernel(class, ka)
        okb <- k255RecordKernel("P-d", "K2-b", class, "kernel", db)
        oka <- k255RecordKernel("P-d", "K2-a", class, "kernel", da)
        if (okb && oka) {
            tb <- db$value$per_assay_results
            ta <- da$value$per_assay_results
            k255CompareTables("P-d", "K2-b", class, "kernel", "LCMS_Neg (K2-b vs K2-a)", tb$LCMS_Neg, ta$LCMS_Neg)
            k255CompareTables("P-d", "K2-b", class, "kernel", "LCMS_Pos (K2-b vs K2-a)", tb$LCMS_Pos, ta$LCMS_Pos)
            k255RecordBCells("P-d", "K2-b", class, "kernel", tb)
        }
        dfb <- k255ResidualDf(kb, class, "LCMS_Neg")
        dfa <- k255ResidualDf(ka, class, "LCMS_Neg")
        k255Record("P-d", "K2-b", class, "kernel", "LCMS_Neg residual df table (K2-b)", paste(names(table(dfb)), as.integer(table(dfb)), sep = ":", collapse = " "))
        k255Record("P-d", "K2-a", class, "kernel", "LCMS_Neg residual df table (K2-a)", paste(names(table(dfa)), as.integer(table(dfa)), sep = ":", collapse = " "))
        k255Record("P-d", "K2-b", class, "kernel", "residual df equal feature by feature (K2-b vs K2-a)", isTRUE(all.equal(unname(dfb), unname(dfa))))
        before <- as.matrix(as.data.frame(methods::slot(k255Obj(values, class), k255Slot(class))$LCMS_Neg)[, GROUP_B])
        after <- as.matrix(as.data.frame(methods::slot(kb, k255Slot(class))$LCMS_Neg)[, GROUP_B])
        storage.mode(before) <- "double"
        storage.mode(after) <- "double"
        k255Record("P-d", "K2-b", class, "kernel", "features with group-B observed count >= 3 before and < 3 after the blanking (of 520)",
            sum(rowSums(!is.na(before)) >= 3L & rowSums(!is.na(after)) < 3L))
    })

    # P-d: the declaration route on a fresh owner; the tab's preflight.
    k255Case("P-d", "K2-b", class, "declaration", function() {
        fx <- kit172Quiet(kit172Owner(file.path(k255NewRoot(paste0("pd-kb-", class)), "project"), kb, class))
        on.exit(k255Close(fx$owner), add = TRUE)
        k255Da("P-d", "K2-b", class, "declaration", fx, binding, getWorkflowStateReference(fx$owner))
        pre <- k255Preflight(fx, binding, class)
        k255Record("P-d", "K2-b", class, "tab in process", "omicDaPreflight value",
            if (is.null(pre$value)) "NULL (the run may go on)" else sprintf("status=%s; stage=%s; errorMessage=%s", pre$value$status, pre$value$stage, pre$value$errorMessage))
        k255Record("P-d", "K2-b", class, "tab in process", "notifications", if (length(pre$notes)) paste(pre$notes, collapse = " | ") else "none")
    })

    # P-f: K3-b, non-RUV.
    k255Case("P-f", "K3-b", class, "kernel", function() {
        run <- k255Kernel(class, k3)
        if (k255RecordKernel("P-f", "K3-b", class, "kernel", run)) {
            t <- as.data.frame(run$value$per_assay_results$LCMS_Neg)
            k255Record("P-f", "K3-b", class, "kernel", "LCMS_Neg table: rows; rows with finite logFC; rows with finite P.Value / raw p",
                sprintf("%d; %d; %d", nrow(t), if ("logFC" %in% names(t)) sum(is.finite(t$logFC)) else -1L,
                    sum(vapply(intersect(c("P.Value", "raw_pvalue", "pvalue"), names(t)), function(n) sum(is.finite(t[[n]])), 1L))))
            k255Record("P-f", "K3-b", class, "kernel", "LCMS_Neg table columns", names(t))
            k255Record("P-f", "K3-b", class, "kernel", "result fields", names(run$value))
        }
    })

    # P-e, P-f (RUV arm): the factors are estimated by a real RUV call on the object WITHOUT missing cells and logged values
    # (as the pipeline's RUV step sees them), then the exported kernel tests the uncorrected values with the failing column blank.
    k255Case("P-e", "K2-b-RUV", class, "kernel", function() {
        full <- k255Log(k255K2Values(missing = FALSE))
        ruv <- k255Capture(runOmicsRuvFactors(k255Obj(full, class), ruv_design = "groups_all_features", ruv_k = 1L, ruv_k_rule = "user",
            ruv_grouping_variable = "group"))
        if (ruv$failed) {
            k255Record("P-e", "K2-b-RUV", class, "kernel", "runOmicsRuvFactors", paste("failed:", conditionMessage(ruv$value)))
            return(invisible(NULL))
        }
        factors <- ruv$value$factors
        k255Record("P-e", "K2-b-RUV", class, "kernel", "factors: status; k; W rows; B_6 in W (per assay)",
            paste(names(factors), vapply(factors, function(f) sprintf("%s; %s; %s; %s", f$status, f$k %||% NA, if (is.null(f$W)) 0L else nrow(f$W),
                if (is.null(f$W)) NA else S %in% rownames(f$W)), ""), sep = "=", collapse = " | "))
        b <- k255Kernel(class, k255Obj(k255Blank(full, "LCMS_Neg", S), class), factors)
        a <- k255Kernel(class, k255Obj(k255Without(full, S), class, k255DesignWithout(S)), factors)
        okb <- k255RecordKernel("P-e", "K2-b-RUV", class, "kernel", b)
        oka <- k255RecordKernel("P-e", "K2-a-RUV", class, "kernel", a)
        if (okb && oka) {
            k255CompareTables("P-e", "K2-b-RUV", class, "kernel", "LCMS_Neg (K2-b vs K2-a, same factors)", b$value$per_assay_results$LCMS_Neg,
                a$value$per_assay_results$LCMS_Neg)
            k255CompareTables("P-e", "K2-b-RUV", class, "kernel", "LCMS_Pos (K2-b vs K2-a, same factors)", b$value$per_assay_results$LCMS_Pos,
                a$value$per_assay_results$LCMS_Pos)
        }
        # P-f, RUV arm: group B blank in LCMS_Neg.
        k3 <- k255Kernel(class, k255Obj(k255Blank(full, "LCMS_Neg", GROUP_B), class), factors)
        if (k255RecordKernel("P-f", "K3-b-RUV", class, "kernel", k3)) {
            k255Record("P-f", "K3-b-RUV", class, "kernel", "ruv_not_observed", utils::capture.output(str(k3$value$ruv_not_observed, give.attr = FALSE)))
            k255Record("P-f", "K3-b-RUV", class, "kernel", "ruv_not_tested (names)", names(k3$value$ruv_not_tested) %||% "none")
        }
    })

    # P-g: figures on K2-b.
    k255Case("P-g", "K2-b", class, "figures", function() {
        summarise <- function(p) {
            if (is.null(p)) return("NULL")
            if (inherits(p, "ggplot")) p <- list(p)
            paste(vapply(seq_along(p), function(i) {
                el <- p[[i]]
                nm <- names(p)[[i]] %||% as.character(i)
                if (is.null(el)) return(paste0(nm, ": NULL"))
                d <- tryCatch(as.data.frame(el$data), error = function(e) NULL)
                if (is.null(d)) return(paste0(nm, ": ", class(el)[[1L]]))
                col <- names(d)[vapply(d, function(x) any(as.character(x) %in% RUNS), logical(1))][1L]
                if (is.na(col)) return(sprintf("%s: %d rows, no run column", nm, nrow(d)))
                sprintf("%s: %d distinct samples, B_6 %s", nm, length(unique(as.character(d[[col]]))), S %in% as.character(d[[col]]))
            }, ""), collapse = " | ")
        }
        plot <- function(label, expr) {
            run <- k255Capture(expr)
            k255Record("P-g", "K2-b", class, "figures", label, if (run$failed) paste("failed:", conditionMessage(run$value)) else summarise(run$value))
        }
        plot("plotPca", plotPca(kb, grouping_variable = "group"))
        plot("plotDensity", plotDensity(kb, grouping_variable = "group"))
        plot("plotRle", plotRle(kb, grouping_variable = "group"))
        da <- k255Kernel(class, kb)
        if (da$failed) {
            k255Record("P-g", "K2-b", class, "figures", "DA for the figures", paste("failed:", conditionMessage(da$value)))
            return(invisible(NULL))
        }
        heat <- if (identical(class, "lipid")) generateLipidDAHeatmap else generateMetabDAHeatmap
        volcano_static <- if (identical(class, "lipid")) generateLipidDAVolcanoStatic else generateMetabDAVolcanoStatic
        glimma <- if (identical(class, "lipid")) generateLipidDAVolcanoPlotGlimma else generateMetabDAVolcanoPlotGlimma
        for (view in c("LCMS_Neg", "Combined")) {
            h <- k255Capture(heat(da$value, selected_contrast = "groupB-groupA", selected_assay = view, top_n = 50, da_q_val_thresh = 2))
            if (h$failed) {
                k255Record("P-g", "K2-b", class, "figures", paste("heatmap", view), paste("failed:", conditionMessage(h$value)))
            } else if (is.null(h$value)) {
                k255Record("P-g", "K2-b", class, "figures", paste("heatmap", view), "NULL (no heatmap)")
            } else {
                m <- tryCatch(h$value$plot@matrix, error = function(e) NULL)
                if (is.null(m)) {
                    k255Record("P-g", "K2-b", class, "figures", paste("heatmap", view), "built; the matrix is not readable")
                } else {
                    k255Record("P-g", "K2-b", class, "figures", paste("heatmap", view, "matrix dim; columns; B_6 column zero cells / non-zero cells"),
                        sprintf("%d x %d; %s; %s", nrow(m), ncol(m), paste(colnames(m), collapse = "+"),
                            if (S %in% colnames(m)) sprintf("%d / %d", sum(m[, S] == 0), sum(m[, S] != 0)) else "no B_6 column"))
                }
            }
        }
        st <- k255Capture(volcano_static(da$value, selected_contrast = "groupB-groupA", selected_assay = "LCMS_Neg"))
        k255Record("P-g", "K2-b", class, "figures", "static volcano LCMS_Neg", if (st$failed) paste("failed:", conditionMessage(st$value)) else class(st$value)[[1L]])
        gl <- k255Capture(glimma(da$value, selected_contrast = "groupB-groupA", selected_assay = "LCMS_Neg"))
        k255Record("P-g", "K2-b", class, "figures", "Glimma volcano LCMS_Neg", if (gl$failed) paste("failed:", conditionMessage(gl$value)) else if (is.null(gl$value)) "NULL" else class(gl$value)[[1L]])
    })

    # P-h and P-i: the writer, the model's alignment and the counts on an owner at norm complete holding K2-b.
    k255Case("P-h", "K2-b", class, "integration", function() {
        staged <- k255Staged(class, kb, paste0("ph-kb-", class))
        on.exit(k255Close(staged$fx$owner), add = TRUE)
        k255Record("P-h", "K2-b", class, "integration", "owner at norm complete", staged$ok)
        meta <- readWorkflowOmicsMetadata(binding, staged$fx$owner, NULL)
        k255Record("P-i", "K2-b", class, "tab summary", "design rows (metadata reader); Removed (original - filtered, both the same state)",
            sprintf("%d; %d", nrow(meta$design_matrix), nrow(meta$design_matrix) - nrow(meta$design_matrix)))
        object <- k255State(staged$fx, binding, NULL)
        design_runs <- k255DesignRuns(object)
        k255Record("P-i", "K2-b", class, "session counts", "samples counted per assay (intersect(columns, design))",
            vapply(methods::slot(object, binding$data_slot), function(t) length(intersect(colnames(t), design_runs)), 1L))
        dir <- k255WriteInput("P-h", "K2-b", class, "integration", binding, staged$fx$owner)
        if (is.null(dir)) return(invisible(NULL))
        read <- k255Capture(integrationViewsFromInputs(dir))
        if (read$failed) {
            k255Record("P-h", "K2-b", class, "model", "integrationViewsFromInputs", paste("failed:", conditionMessage(read$value)))
            return(invisible(NULL))
        }
        object <- read$value$object
        for (ms in c("keep_missing", "drop_incomplete")) {
            parameters <- k255Capture(integrationResolveMofaParameters(object, list(missing_samples = ms)))
            if (parameters$failed) {
                k255Record("P-h", "K2-b", class, "model", paste(ms, "parameters"), paste("failed:", conditionMessage(parameters$value)))
                next
            }
            whole <- k255Capture(integrationAlignViews(object, parameters$value))
            paged <- k255Capture(integrationAlignViewsPartitioned(object, parameters$value, integrationObjectViewPages(object)))
            for (route in c("whole", "partitioned")) {
                run <- if (identical(route, "whole")) whole else paged
                if (run$failed) {
                    k255Record("P-h", "K2-b", class, paste("model", route), paste(ms, "alignment"), paste("failed:", conditionMessage(run$value)))
                    next
                }
                v <- run$value
                k255Record("P-h", "K2-b", class, paste("model", route), paste(ms, ": shared samples; kept samples; columns per view; B_6 in a view"),
                    sprintf("%s; %s; %s; %s", v$shared_samples, NROW(v$samples), paste(names(v$views), vapply(v$views, ncol, 1L), sep = "=", collapse = " "),
                        paste(names(v$views), vapply(v$views, function(m) S %in% colnames(m), logical(1)), sep = "=", collapse = " ")))
                k255Record("P-h", "K2-b", class, paste("model", route), paste(ms, ": view summary"),
                    substr(tryCatch(as.character(jsonlite::toJSON(v$view_summary, auto_unbox = TRUE)), error = function(e) "unreadable"), 1L, 600L))
                k255Record("P-h", "K2-b", class, paste("model", route), paste(ms, ": warning classes"), if (length(run$classes)) unique(run$classes) else "none")
            }
            if (!whole$failed && !paged$failed) {
                k255Record("P-h", "K2-b", class, "model", paste(ms, ": whole and partitioned views identical"), identical(whole$value$views, paged$value$views))
            }
        }
    })
}

# ---- the RUV-applied objects: P-k (before the correlation step), then P-a, P-b, P-c ----

k255Real <- get("readWorkflowOmicsState", envir = asNamespace("MultiScholaR"))
options(MultiScholaR.ruv_design_panel = TRUE)

for (class in c("lipid", "metab")) {
    for (kind in c("new", "leg")) {
        built <- k255Case("BUILD", paste0("K2-RUV", kind), class, "declaration", function() k255BuildRuv(kind, class))
        if (is.null(built) || !isTRUE(built$ok)) {
            k255Record("BUILD", paste0("K2-RUV", kind), class, "declaration", "not run", "the object did not reach the ruv_corrected head")
            next
        }
        name <- built$name
        binding <- built$binding
        owner <- built$fx$owner

        # The copy and the layout before the press.
        before <- k255State(built$fx, binding, NULL)
        k255Record("P-a", name, class, "declaration (before the press)", "args$ruv_uncorrected_assays: B_6 cells in LCMS_Neg", k255CopyCells(before))
        copy_before <- before@args$ruv_uncorrected_assays

        # P-k: the writer on the pre-correlation head, the failing column blank in the final read and in the normalised read.
        k255Case("P-k", paste0(name, "-b"), class, "integration", function() {
            ruv_id <- getWorkflowStateReference(owner, paste0(built$prefix, "_ruv_corrected"))$entity_id
            norm_id <- getWorkflowStateReference(owner, paste0(built$prefix, "_normalized"))$entity_id
            blank <- function(object, copy) {
                assays <- methods::slot(object, binding$data_slot)
                assays$LCMS_Neg[[S]] <- NA_real_
                methods::slot(object, binding$data_slot) <- assays
                if (copy && !is.null(object@args$ruv_uncorrected_assays)) object@args$ruv_uncorrected_assays$LCMS_Neg[[S]] <- NA_real_
                object
            }
            testthat::local_mocked_bindings(readWorkflowOmicsState = function(binding, manager, reference = NULL) {
                object <- k255Real(binding, manager, reference)
                id <- if (is.character(reference)) getWorkflowStateReference(manager, reference)$entity_id else reference$entity_id
                if (identical(id, ruv_id)) blank(object, copy = TRUE) else if (identical(id, norm_id)) blank(object, copy = FALSE) else object
            }, .package = "MultiScholaR")
            k255WriteInput("P-k", paste0(name, "-b"), class, "integration (blanked by hand)", binding, owner)
        })

        # The correlation step, then P-a, P-b, P-c.
        if (!k255Correlate(built, class)) {
            k255Record("P-a", name, class, "declaration", "not run", "the correlation step did not complete")
            next
        }
        k255Case("P-a", name, class, "declaration", function() {
            after <- k255State(built$fx, binding, NULL)
            k255Layout("P-a", name, class, "declaration", after, binding)
            k255Record("P-a", name, class, "declaration", "args$ruv_uncorrected_assays: B_6 cells in LCMS_Neg", k255CopyCells(after))
            k255Record("P-a", name, class, "declaration", "copy identical before and after the press", identical(copy_before, after@args$ruv_uncorrected_assays))
        })
        k255Case("P-b", name, class, "declaration", function() {
            k255Da("P-b", name, class, "declaration", built$fx, binding, getWorkflowStateReference(owner))
        })
        k255Case("P-c", name, class, "integration", function() k255WriteInput("P-c", name, class, "integration", binding, owner))
        k255Close(owner)
    }
}

# ---- P-a, the tab: in process then managed, new design ----

for (class in c("lipid", "metab")) {
    in_process <- k255Case("P-a", "K2-RUVnew", class, "tab in process", function() {
        seen <- k194Notes()
        fx <- k255Fixture(file.path(k255NewRoot(paste0("pa-ip-", class)), "project"), k255K2Values(missing = FALSE), class)
        on.exit(k255Close(fx$owner), add = TRUE)
        nd <- omicNormReactiveState()
        k255Gate(k194Run(fx, nd, k194InputNew()))
        k255Record("P-a", "K2-RUVnew", class, "tab in process", "run head", fx$owner$current_state)
        pairs <- NULL
        trace <- get("runWorkflowOmicsReferenceCorrelation", envir = asNamespace("MultiScholaR"))
        testthat::local_mocked_bindings(runWorkflowOmicsReferenceCorrelation = function(...) {
            result <- trace(...)
            pairs <<- result$correlations
            result
        }, .package = "MultiScholaR")
        k255Gate(k194Apply(fx, nd, k194InputNew()))
        k255Record("P-a", "K2-RUVnew", class, "tab in process", "Apply head", fx$owner$current_state)
        tryCatch({
            if (!is.null(pairs)) k255RecordBest("STOP-i", "K2-RUVnew", class, "tab in process (intensity filter, log2, quantile, new RUV)", pairs)
            k255Layout("P-a", "K2-RUVnew", class, "tab in process", k255State(fx, fx$binding, NULL), fx$binding)
            k255Record("P-a", "K2-RUVnew", class, "tab in process", "args$ruv_uncorrected_assays: B_6 cells in LCMS_Neg", k255CopyCells(k255State(fx, fx$binding, NULL)))
        }, error = function(e) k255Record("P-a", "K2-RUVnew", class, "tab in process", "RECORDER ERROR", k255Error(e)))
        k255Record("P-a", "K2-RUVnew", class, "tab in process", "notes", paste(seen$notes, collapse = " | "))
        if (identical(fx$owner$current_state, paste0(fx$prefix, "_correlation_filtered"))) k255Runs(fx) else NULL
    })
    k255Case("P-a", "K2-RUVnew", class, "managed", function() {
        seen <- k194Notes()
        k255Gate(k255Managed(class, paste0("pa-managed-", class), k255K2Values(missing = FALSE), function(fx, probe) {
            k255Record("P-a", "K2-RUVnew", class, "managed", "run completed", isTRUE(probe$ran) && isTRUE(probe$settled_run))
            k255Record("P-a", "K2-RUVnew", class, "managed", "Apply settled", isTRUE(probe$settled_apply))
            k255Record("P-a", "K2-RUVnew", class, "managed", "operations", paste(probe$operations, collapse = ","))
            k255Record("P-a", "K2-RUVnew", class, "managed", "statuses", paste(probe$statuses, collapse = ","))
            k255Record("P-a", "K2-RUVnew", class, "managed", "head", fx$owner$current_state)
            k255Record("P-a", "K2-RUVnew", class, "managed", "errors in log", paste(grep("ERROR", fx$log$lines, value = TRUE), collapse = " | "))
            k255Record("P-a", "K2-RUVnew", class, "managed", "notes", paste(seen$notes, collapse = " | "))
            runs <- k255Runs(fx)
            tryCatch({
                k255Layout("P-a", "K2-RUVnew", class, "managed", k255State(fx, fx$binding, NULL), fx$binding)
                k255Record("P-a", "K2-RUVnew", class, "managed", "args$ruv_uncorrected_assays: B_6 cells in LCMS_Neg", k255CopyCells(k255State(fx, fx$binding, NULL)))
            }, error = function(e) k255Record("P-a", "K2-RUVnew", class, "managed", "RECORDER ERROR", k255Error(e)))
            k255Record("P-a", "K2-RUVnew", class, "managed", "managed design = in-process design", identical(runs, in_process))
        }))
    })
}

options(MultiScholaR.ruv_design_panel = NULL)

# ---- P-j: the older-form record: the base writes K2's one-assay-failure project, the record, one replay ----

k255Projects <- file.path(k255Root, "projects")
dir.create(k255Projects)
harness <- new.env(parent = globalenv())
sys.source(file.path("tools", "profiling", "workflow_kit.R"), harness)
k255Api <- function(name) get(name, asNamespace("MultiScholaR"))

for (class in c("lipid", "metab")) {
    k255Case("P-j", "K2", class, "project and replay", function() {
        binding <- kit172Binding(class)
        prefix <- binding$prefix
        root <- file.path(k255Projects, class)
        fx <- kit172Quiet(kit172Owner(root, kit172Object(k255K2Values(), class), class))
        owner <- fx$owner
        config <- fx$workflow$config_list
        kit172Quiet({
            commitWorkflowOmicsStoreStep(owner, paste0(prefix, "_qc_finalize"), paste0(prefix, "_qc_complete"), config, "QC complete")
            runWorkflowOmicsReferenceCorrelation(binding, owner, "group", 0.5, config)
        })
        k255Record("P-j", "K2", class, "project", "head", owner$current_state)
        tryCatch(k255Layout("P-j", "K2", class, "project", readWorkflowOmicsState(binding, owner, NULL), binding),
            error = function(e) k255Record("P-j", "K2", class, "project", "RECORDER ERROR", k255Error(e)))
        json <- file.path(root, "analysis_record.json")
        writeWorkflowAnalysisRecord(owner, json)
        k255Record("P-j", "K2", class, "project", "analysis record written (bytes)", file.size(json))
        record <- jsonlite::read_json(json, simplifyVector = FALSE)
        k255Record("P-j", "K2", class, "project", "record step operations", vapply(record$steps %||% list(), function(s) as.character(s$operation %||% "?"), ""))
        k255Close(owner)
        source_dir <- fx$source_dir
        store_root <- dirname(dirname(list.files(source_dir, "^store[.]json$", recursive = TRUE, all.files = TRUE, full.names = TRUE)[[1L]]))
        names_in_order <- names(methods::slot(kit172Object(k255K2Values(), class), binding$data_slot))
        orders <- list(names_in_order, rev(names_in_order))
        replayed <- FALSE
        for (i in seq_along(orders)) {
            if (replayed) break
            copy <- tempfile("kit255-replay-")
            dir.create(copy)
            stopifnot(system2("cp", c("-a", file.path(source_dir, "."), copy)) == 0L)
            copy_root <- file.path(copy, substring(store_root, nchar(source_dir) + 2L))
            prepared <- harness$kitPrepareReplayProject(copy_root, c("initial", paste0(prefix, "_raw_data_s4")), k255Api)
            inputs <- c(stats::setNames(file.path(source_dir, "import_sources", paste0(orders[[i]], "_export.txt")), paste0(prefix, "_assay_", seq_along(orders[[i]]))),
                design_file = file.path(source_dir, "design_matrix.tab"))
            scratch <- tempfile("kit255-scratch-")
            dir.create(scratch)
            replay <- k255Capture(kit172Quiet(harness$kitReplayJourney(record, prepared, store_root, inputs, fx$paths, scratch, k255Api)))
            if (replay$failed) {
                k255Record("P-j", "K2", class, "replay", sprintf("replay attempt %d (assay order %s)", i, paste(orders[[i]], collapse = ",")),
                    paste("failed:", conditionMessage(replay$value)))
            } else {
                k255Record("P-j", "K2", class, "replay", sprintf("replay attempt %d (assay order %s): states; mismatches", i, paste(orders[[i]], collapse = ",")),
                    sprintf("%s; %d (%s)", paste(replay$value$states, collapse = "+"), length(replay$value$mismatches), paste(utils::head(as.character(replay$value$mismatches), 5L), collapse = " | ")))
                replayed <- TRUE
            }
        }
        k255Record("P-j", "K2", class, "replay", "replayed", replayed)
    })
}
k255Case("P-j", "K2", "lipid,metab", "tarball", function() {
    archive <- file.path(k255Root, "base-record.tar.gz")
    old <- setwd(k255Projects)
    on.exit(setwd(old), add = TRUE)
    utils::tar(archive, files = c("lipid", "metab"), compression = "gzip", tar = "internal")
    setwd(old)
    writeLines(jsonlite::base64_enc(readBin(archive, "raw", file.size(archive))), file.path(k255Out, "base-record.tar.gz.b64.txt"))
    writeLines(c(k255Git("rev-parse", "HEAD")[[1L]], digest::digest(file = archive, algo = "sha256")), file.path(k255Out, "base-record.sha.txt"))
    k255Record("P-j", "K2", "lipid,metab", "tarball", "bytes (tar.gz)", file.size(archive))
})

k255Write()
writeLines(c(k255Provenance, sprintf("finished\t%s", format(Sys.time(), tz = "UTC"))), file.path(k255Out, "provenance.txt"))
cat("KIT-255 G1 done\n")
