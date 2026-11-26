# KIT-281 (VAL-354 (2), VAL-365 Q5): the scoring of the calibration gate's rows (ticket sections 5, 6, 9, 10). Evidence code under
# docs/qa: no package code, no package loaded.
#
# Usage, from the package root (any R with jsonlite and digest):
#   Rscript docs/qa/workflow-execution/kit/kit-281/score-k281.R <rows root> <out dir>
# <rows root> holds one directory per job (sim-k281.R's <cell>.tsv.gz, metrics.tsv, range.tsv, provenance-sim.json, and
# determinism.tsv where it ran); a cell's ranges are joined into one cell of replicates before resampling.
# Adapted from docs/qa/workflow-execution/kit/kit-211/score-k211.R (SHA-256 aefc66baa7dc42a86f5f69247911e51b8701200b1e7d95851bec676ae1110d54);
# the blocks between "verbatim" and "end verbatim" markers are KIT-197's lines unchanged and the functions named in script-diff.md
# are KIT-211's text unchanged (both checked by e0-k281.py); every other change is marked KIT-281 or is the new decision section.
#
# Decision (ticket section 5): the deciding quantity is rho_g = alpha_B1(H_B1 without H_B) / alpha_B(H_B) pooled over the twelve
# settings (A to D at vd5, vd10, vd25); the deciding limit is the bootstrap 0.975 quantile (type 7) of rhog.B1 over the 10,000
# resamples; Calibrated needs it at most 1.10 and the guards G1a, G1b, G2 to pass; the jackknife-normal limits are reported and
# decide nothing. A resample whose statistic is not finite is dropped and counted (boot_nonfinite). Unit of resampling: a cell's
# replicate with all its masks and arms. Bootstrap: B = 10,000, set.seed(280000001, Mersenne-Twister, Inversion, Rejection) once,
# then for cells A, B, C, D in that order sample.int(R, B * R, replace = TRUE) filled by row into a B x R index matrix.
# Environment: K211_B_BOOT (default 10000; the self-test's constructed tables read 1000, recorded in checks.tsv), K281_RANGE
# ("first,last": a declared replicate range for the self-test and lifecycle scorings; without it the range is rule G3's table,
# g3-rule.tsv, and the scoring is a counted one), K281_G3_TABLE (a path to a G3 table; the self-test's counted-mode cases),
# K281_GUARD_PROBE=1 (the self-test's trigger of the partial-match guard).
# Script lessons of KIT-188 (script audit 2): every table is read as text and converted by name with a check that no value turns
# into NA; `$` prefix matches (R spells the data-frame message "Partial match of": the text is lowered) are an instrument stop, never
# an ordinary failure, and no `$` is evaluated inside suppressWarnings (the scorer uses no suppressWarnings); provenance and a
# progress file are written before the first unit; the scorer stops on rows outside the declared range and on any smoke flag.

options(warnPartialMatchDollar = TRUE)
# Last net only (the guard itself is a LOCAL calling handler inside the tryCatch of the driver below, so its stop() is caught and recorded by stopRecord): an error that
# escaped every handler writes the Stop outcome before the non-zero exit. No global calling handler is used.
options(error = function() {
    try({
        if (!is.na(ctx[["out"]])) {
            kind_text <- paste("uncaught error:", substr(geterrmessage(), 1L, 300L))
            writeLines(c("outcome\tkind\treason\tcounted", paste("Stop", "instrument", gsub("[\t\n]", " ", kind_text), as.character(ctx[["counted"]]), sep = "\t")),
                file.path(ctx[["out"]], "outcome.tsv"))
        }
    }, silent = TRUE)
    quit(save = "no", status = 1L)
})

CELLS_DECIDING <- c("A", "B", "C", "D")
CELLS_ALL <- CELLS_DECIDING # KIT-281: no Av, no Ag
VD <- c("vd5", "vd10", "vd25")
CELL_MASKS <- list(A = VD, B = VD, C = VD, D = VD) # KIT-281
ARMS_MAIN <- c("S", "B", "B1", "N") # KIT-281
CANDIDATES <- c("B", "B1") # KIT-281
STRATA <- c("one_treated", "one_control", "both", "other_only", "none")
POOL <- list(A = VD, B = VD, C = VD, D = VD) # KIT-281: the twelve settings
B_BOOT <- as.integer(Sys.getenv("K211_B_BOOT", "10000"))
SEED <- 280000001L # KIT-281
RHO_MARGIN <- 1.10 # KIT-281 (VAL-365 Q5; fixed before any data of the gate exist)
FDP_LIMIT <- 0.05 # G1a: the nominal level
GUARD_MARGIN <- 0.0025 # G1b: KIT-211's guard margin
POWER_MARGIN <- 0.005 # G2: KIT-211's "no material cost"
PC_MARGIN <- 1.20 # the positive control (reported, decides nothing)
ALPHA_PRIMARY <- 0.025 # one-sided 97.5 %
ALPHA_SETTING <- 0.05 / 12 # the setting labels (reported)
ALPHA_F_POOL <- 0.05 # F5' is one falsifier: falsifier5()'s halving gives a two-sided 95 % interval
THIN <- c("Phi.S", "Phi.B", "Phi.B1", "rhog.B1", "B1.m6_l_H") # KIT-281
AG_STATS <- c("pow", "fdr_whole", "powbh", "fdr_bh", "tested_per_r", "kept_per_r") # KIT-281: the whole-list rows, by arm
G3_TABLE <- Sys.getenv("K281_G3_TABLE", "docs/qa/workflow-execution/kit/kit-281/g3-rule.tsv")
ROW_SUBSETS <- c("added", "kept_both", "whole", "whole0", "common", "g_added")
ROW_COLUMNS <- c("cell", "mask", "r", "arm", "subset", "class", "stratum", "ell", "band", "ell0", "n", "n_tested", "n_p05", "n_call",
    "n_bh", "sum_err", "sum_s_err", "sum_sgn_err")
ROW_NUMERIC <- c("r", "ell", "band", "ell0", "n", "n_tested", "n_p05", "n_call", "n_bh", "sum_err", "sum_s_err", "sum_sgn_err")
METRIC_COLUMNS <- c("cell", "mask", "r", "arm", "kept", "added", "calls", "t_r", "mono_ok", "qmethod", "threshold", "viol", "e10_dev",
    "e3_pass", "n_planted", "n_down_all", "n1_band_all", "n_onoff_all", "calls_both_s", "calls_new_s", "gained", "lost", "seconds")
METRIC_NUMERIC <- c("r", "kept", "added", "calls", "t_r", "threshold", "viol", "e10_dev", "n_planted", "n_down_all", "n1_band_all",
    "n_onoff_all", "calls_both_s", "calls_new_s", "gained", "lost")
METRIC_LOGICAL <- c("mono_ok", "e3_pass")

k211Abort <- function(class, message) stop(structure(class = c(class, "error", "condition"),
    list(message = message, call = NULL)))

# ==== verbatim: score-k197.R:44-49 at ef63c918 ====
num <- function(x) ifelse(is.na(x), "NA", sprintf("%.17g", x))

writeTsv <- function(df, path) {
    for (column in names(df)) if (is.double(df[[column]])) df[[column]] <- num(df[[column]])
    utils::write.table(df, path, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
}
# ==== end verbatim ====

# ==== verbatim: score-k197.R:219-223 at ef63c918 ====
jackSe <- function(theta) {
    n <- length(theta)
    if (n < 2L || any(!is.finite(theta))) return(NA_real_)
    sqrt((n - 1) / n * sum((theta - mean(theta))^2))
}
# ==== end verbatim ====

# ==== verbatim: score-k197.R:225-225 at ef63c918 ====
quant <- function(boot, p) apply(boot, 2L, stats::quantile, probs = p, type = 7, na.rm = TRUE, names = FALSE)
# ==== end verbatim ====

# ---- input (KIT-281: ADAPTED, script lesson 1: every column is read as text and converted by name) -----------------------------
# KIT-211's readInput let read.delim guess the classes of the numeric columns; here every column is character, `na.strings` is empty
# (a cell, mask or arm labelled "NA", "T" or "F" is a label and stays text; the vocabulary check below stops on it) and each numeric
# column is converted through toNum(), which stops when any value is not a number (no value turns into NA silently).
readText <- function(path) {
    connection <- if (grepl("\\.gz$", path)) gzfile(path) else path
    utils::read.delim(connection, colClasses = "character", na.strings = character(0), check.names = FALSE, stringsAsFactors = FALSE)
}

toNum <- function(df, columns, what) {
    for (column in columns) {
        text <- df[[column]]
        good <- grepl("^-?[0-9]+(\\.[0-9]+)?([eE][-+]?[0-9]+)?$", text)
        if (!all(good)) k211Abort("k211_instrument", sprintf("%s: column %s holds a value that is not a number ('%s'; the text is never coerced to NA)",
            what, column, text[!good][[1L]]))
        df[[column]] <- as.numeric(text)
    }
    df
}

toLogical <- function(df, columns, what) {
    for (column in columns) {
        text <- df[[column]]
        value <- ifelse(text %in% c("TRUE", "T"), TRUE, ifelse(text %in% c("FALSE", "F"), FALSE, NA))
        if (anyNA(value)) k211Abort("k211_instrument", sprintf("%s: column %s holds a value that is not TRUE, FALSE, T or F ('%s')",
            what, column, text[is.na(value)][[1L]]))
        df[[column]] <- value
    }
    df
}

needColumns <- function(df, columns, what) {
    missing <- setdiff(columns, names(df))
    if (length(missing)) k211Abort("k211_instrument", sprintf("%s: the columns %s are absent", what, paste(missing, collapse = ", ")))
}

needLabels <- function(values, allowed, what) {
    bad <- setdiff(unique(values), allowed)
    if (length(bad)) k211Abort("k211_instrument", sprintf("%s: unexpected label '%s' (allowed: %s)", what, bad[[1L]], paste(allowed, collapse = ", ")))
}

# The declared replicate range (ticket sections 7, 8; RANGE-1): rule G3's table for a counted scoring, K281_RANGE for the self-test's
# and the lifecycle's mechanics scorings. Returns the range by cell, its source and the table's text hash.
declaredRange <- function() {
    given <- Sys.getenv("K281_RANGE")
    if (nzchar(given)) {
        parts <- strsplit(given, ",", fixed = TRUE)[[1L]]
        if (length(parts) != 2L || !all(grepl("^[0-9]+$", parts))) k211Abort("k211_instrument", "K281_RANGE must read '<first>,<last>'")
        return(list(range = lapply(stats::setNames(CELLS_ALL, CELLS_ALL), \(cellname) seq.int(as.integer(parts[[1L]]), as.integer(parts[[2L]]))),
            source = "env", sha = NA_character_))
    }
    if (!file.exists(G3_TABLE)) k211Abort("k211_instrument", sprintf("no declared replicate range: %s is absent and K281_RANGE is not set", G3_TABLE))
    g3 <- readText(G3_TABLE)
    needColumns(g3, c("cell", "r_cell", "first_r", "last_r"), "the G3 table")
    g3 <- toNum(g3, c("r_cell", "first_r", "last_r"), "the G3 table")
    if (!setequal(g3[["cell"]], CELLS_ALL) || anyDuplicated(g3[["cell"]]) || any(g3[["last_r"]] - g3[["first_r"]] + 1 != g3[["r_cell"]]) ||
        length(unique(g3[["r_cell"]])) != 1L) k211Abort("k211_instrument", "the G3 table does not hold one R for cells A to D with last r - first r + 1 = R")
    list(range = lapply(stats::setNames(CELLS_ALL, CELLS_ALL), \(cellname) seq.int(g3[["first_r"]][g3[["cell"]] == cellname],
        g3[["last_r"]][g3[["cell"]] == cellname])), source = "g3-table", sha = digest::digest(file = G3_TABLE, algo = "sha256"))
}

# Every job directory under <rows root>: a cell's ranges are joined into one set of replicates before resampling.
readInput <- function(root, declared) {
    dirs <- sort(list.dirs(root, recursive = FALSE))
    dirs <- dirs[file.exists(file.path(dirs, "metrics.tsv"))]
    if (!length(dirs)) k211Abort("k211_instrument", "no job directory with metrics.tsv under the rows root")
    rows <- list()
    metrics <- list()
    ranges <- list()
    determinism <- list()
    for (dir in dirs) {
        what <- basename(dir)
        # a stopped or smoke job leaves its record in the job directory: the scorer refuses it
        for (stop_file in c("instrument-stop.txt", "errors.txt")) if (file.exists(file.path(dir, stop_file))) k211Abort("k211_instrument",
            sprintf("%s: the job directory holds %s (a partial-match guard or a task failure stopped the job; an instrument defect)", what, stop_file))
        provenance_file <- file.path(dir, "provenance-sim.json")
        if (!file.exists(provenance_file)) k211Abort("k211_instrument", sprintf("%s: provenance-sim.json is absent (the job did not complete)", what))
        # the provenance is read as text (no JSON parser is trusted with a partial match): the smoke flag must read false
        provenance_text <- paste(readLines(provenance_file, warn = FALSE), collapse = "\n")
        if (!grepl("\"smoke\"[[:space:]]*:[[:space:]]*false", provenance_text)) k211Abort("k211_instrument", sprintf(
            "%s: the provenance carries a smoke flag, or none (an instrument check is never scored)", what))
        if (!grepl("\"ticket\"[[:space:]]*:[[:space:]]*\"KIT-281\"", provenance_text) || !grepl("\"status\"[[:space:]]*:[[:space:]]*\"complete\"", provenance_text))
            k211Abort("k211_instrument", sprintf("%s: the provenance is not KIT-281's complete record", what))
        range <- readText(file.path(dir, "range.tsv"))
        needColumns(range, c("cell", "first", "last", "expected", "complete"), paste(what, "range.tsv"))
        range <- toNum(range, c("first", "last", "expected", "complete"), paste(what, "range.tsv"))
        needLabels(range[["cell"]], CELLS_ALL, paste(what, "range.tsv cell"))
        range[["dir"]] <- what
        ranges[[length(ranges) + 1L]] <- range
        m <- readText(file.path(dir, "metrics.tsv"))
        needColumns(m, METRIC_COLUMNS, paste(what, "metrics.tsv"))
        m <- toLogical(toNum(m, METRIC_NUMERIC, paste(what, "metrics.tsv")), METRIC_LOGICAL, paste(what, "metrics.tsv"))
        needLabels(m[["cell"]], CELLS_ALL, paste(what, "metrics.tsv cell"))
        needLabels(m[["mask"]], VD, paste(what, "metrics.tsv mask"))
        needLabels(m[["arm"]], ARMS_MAIN, paste(what, "metrics.tsv arm"))
        metrics[[length(metrics) + 1L]] <- m
        for (cellname in range[["cell"]]) {
            file <- file.path(dir, paste0(cellname, ".tsv.gz"))
            if (!file.exists(file)) k211Abort("k211_instrument", sprintf("%s: %s.tsv.gz is absent", what, cellname))
            d <- readText(file)
            needColumns(d, ROW_COLUMNS, paste(what, cellname))
            d <- toNum(d, ROW_NUMERIC, paste(what, cellname))
            needLabels(d[["cell"]], CELLS_ALL, paste(what, cellname, "cell"))
            needLabels(d[["mask"]], VD, paste(what, cellname, "mask"))
            needLabels(d[["arm"]], ARMS_MAIN, paste(what, cellname, "arm"))
            needLabels(d[["subset"]], ROW_SUBSETS, paste(what, cellname, "subset"))
            needLabels(d[["class"]], c("null", "up", "down"), paste(what, cellname, "class"))
            needLabels(d[["stratum"]], c(STRATA, "."), paste(what, cellname, "stratum"))
            rows[[cellname]] <- c(rows[[cellname]], list(d))
        }
        det_file <- file.path(dir, "determinism.tsv")
        if (file.exists(det_file)) {
            det <- readText(det_file)
            needColumns(det, c("cell", "r", "rows_identical", "metrics_identical"), paste(what, "determinism.tsv"))
            determinism[[length(determinism) + 1L]] <- toLogical(det, c("rows_identical", "metrics_identical"), paste(what, "determinism.tsv"))
        }
    }
    metrics <- do.call(rbind, metrics)
    rows <- lapply(rows, \(x) do.call(rbind, x))
    # RANGE-1: no scored row lies outside the declared range of its cell
    for (cellname in names(rows)) {
        if (!all(rows[[cellname]][["r"]] %in% declared[["range"]][[cellname]]))
            k211Abort("k211_instrument", sprintf("cell %s: a row lies outside the declared range %d to %d", cellname,
                min(declared[["range"]][[cellname]]), max(declared[["range"]][[cellname]])))
    }
    for (cellname in CELLS_ALL) {
        r_metrics <- metrics[["r"]][metrics[["cell"]] == cellname]
        if (!all(r_metrics %in% declared[["range"]][[cellname]])) k211Abort("k211_instrument", sprintf(
            "cell %s: a metrics row lies outside the declared range %d to %d", cellname, min(declared[["range"]][[cellname]]), max(declared[["range"]][[cellname]])))
    }
    list(rows = rows, metrics = metrics, ranges = do.call(rbind, ranges), determinism = if (length(determinism)) do.call(rbind, determinism))
}


# ---- atoms ----------------------------------------------------------------------------------------------------------
# The sums of one (cell, mask, arm) per replicate: an R x K matrix (rows in the order of `replicates`).
armAtoms <- function(d, m, replicates) {
    R <- length(replicates)
    key <- factor(d$r, levels = replicates)
    total <- function(keep, column) {
        if (!any(keep)) return(numeric(R))
        v <- tapply(d[[column]][keep], key[keep], sum)
        v[is.na(v)] <- 0
        as.numeric(v)
    }
    mm <- m[match(replicates, m$r), , drop = FALSE]
    is_null <- d$class == "null"
    is_planted <- d$class %in% c("up", "down")
    one_sided <- d$stratum %in% c("one_treated", "one_control")
    sub <- d$subset
    atoms <- list(one = rep(1, R))
    for (s in c("added", "kept_both", "whole")) {
        in_s <- sub == s
        atoms[[paste0("n0_", s)]] <- total(in_s & is_null, "n_tested")
        atoms[[paste0("p05_", s)]] <- total(in_s & is_null, "n_p05")
        atoms[[paste0("V_", s)]] <- total(in_s & is_null, "n_call")
        atoms[[paste0("n1_", s)]] <- total(in_s & is_planted, "n_tested")
        atoms[[paste0("T_", s)]] <- total(in_s & is_planted, "n_call")
        atoms[[paste0("err_", s)]] <- total(in_s & is_planted, "sum_err")
        atoms[[paste0("sgn_", s)]] <- total(in_s & is_planted, "sum_sgn_err")
    }
    for (l in 0:1) {
        atoms[[paste0("n0_added_l", l)]] <- total(sub == "added" & is_null & d$ell == l, "n_tested")
        atoms[[paste0("p05_added_l", l)]] <- total(sub == "added" & is_null & d$ell == l, "n_p05")
        atoms[[paste0("n0_whole_l", l)]] <- total(sub == "whole" & is_null & d$ell == l, "n_tested")
        atoms[[paste0("p05_whole_l", l)]] <- total(sub == "whole" & is_null & d$ell == l, "n_p05")
        atoms[[paste0("n0_whole0_l", l)]] <- total(sub == "whole0" & is_null & d$ell0 == l, "n_tested")
        atoms[[paste0("p05_whole0_l", l)]] <- total(sub == "whole0" & is_null & d$ell0 == l, "n_p05")
        keep <- sub == "added" & is_null & one_sided & d$ell == l
        atoms[[paste0("m6sum_", l)]] <- total(keep, "sum_s_err")
        atoms[[paste0("m6n_", l)]] <- total(keep, "n_tested")
        keep <- sub == "whole" & is_null & one_sided & d$ell == l
        atoms[[paste0("m6Hsum_", l)]] <- total(keep, "sum_s_err")
        atoms[[paste0("m6Hn_", l)]] <- total(keep, "n_tested")
    }
    atoms$n0_common_l1 <- total(sub == "common" & is_null & d$ell == 1L, "n_tested")
    atoms$p05_common_l1 <- total(sub == "common" & is_null & d$ell == 1L, "n_p05")
    for (st in STRATA) for (l in 0:1) {
        keep <- sub == "added" & is_null & d$stratum == st & d$ell == l
        atoms[[paste0("sn0_", st, "_", l)]] <- total(keep, "n_tested")
        atoms[[paste0("sp05_", st, "_", l)]] <- total(keep, "n_p05")
    }
    for (direction in c("up", "down")) {
        keep <- sub == "whole" & d$class == direction
        atoms[[paste0("n_", direction)]] <- total(keep, "n_tested")
        atoms[[paste0("T_", direction)]] <- total(keep, "n_call")
        atoms[[paste0("bh_", direction)]] <- total(keep, "n_bh")
        atoms[[paste0("sgn_", direction)]] <- total(keep, "sum_sgn_err")
    }
    band <- sub == "whole" & d$band == 1L
    atoms$n0b <- total(band & is_null, "n_tested")
    atoms$p05b <- total(band & is_null, "n_p05")
    atoms$Tb <- total(band & is_planted, "n_call")
    atoms$n1b <- total(band & is_planted, "n_tested")
    atoms$sgnb <- total(band & is_planted, "sum_sgn_err")
    atoms$tested_whole <- total(sub == "whole", "n_tested")
    atoms$calls_whole <- total(sub == "whole", "n_call")
    atoms$kept_whole <- total(sub == "whole", "n")
    atoms$on_kept <- total(sub == "onoff", "n")
    atoms$on_call <- total(sub == "onoff", "n_call")
    atoms$on_bh <- total(sub == "onoff", "n_bh")
    atoms$planted <- mm$n_planted
    atoms$down_all <- mm$n_down_all
    atoms$n1b_all <- mm$n1_band_all
    atoms$onoff_all <- mm$n_onoff_all
    atoms$calls_both_s <- mm$calls_both_s
    atoms$calls_new_s <- mm$calls_new_s
    atoms$gained <- mm$gained
    atoms$lost <- mm$lost
    atoms$n0t_added <- atoms$n0_added * mm$t_r
    atoms$n0t_kept_both <- atoms$n0_kept_both * mm$t_r
    # KIT-281: the gate's atoms. g_added = arm B1's proteins not in H_B (null rows by the flag ell); the BH calls of the whole list.
    in_g <- sub == "g_added"
    atoms$n0_gadd <- total(in_g & is_null, "n_tested")
    atoms$p05_gadd <- total(in_g & is_null, "n_p05")
    atoms$n0_gadd_l0 <- total(in_g & is_null & d$ell == 0L, "n_tested")
    atoms$p05_gadd_l0 <- total(in_g & is_null & d$ell == 0L, "n_p05")
    atoms$n0_gadd_l1 <- total(in_g & is_null & d$ell == 1L, "n_tested")
    atoms$p05_gadd_l1 <- total(in_g & is_null & d$ell == 1L, "n_p05")
    atoms$bhV_whole <- total(sub == "whole" & is_null, "n_bh")
    atoms$bhT_whole <- total(sub == "whole" & is_planted, "n_bh")
    calls <- atoms$V_whole + atoms$T_whole
    atoms$fdp_sum <- ifelse(calls > 0, atoms$V_whole / calls, 0)
    atoms$fdp_n <- rep(1, R)
    do.call(cbind, atoms)
}

# ---- estimands (6.3), on a matrix of summed atoms (one row per resample, jackknife unit or the estimate) ----------------
armStats <- function(S, arm) {
    g <- function(name) S[, paste0(arm, ".", name)]
    a_o <- g("p05_kept_both") / g("n0_kept_both")
    a_l <- g("p05_added_l1") / g("n0_added_l1")
    fdr <- function(s) g(paste0("V_", s)) / (g(paste0("V_", s)) + g(paste0("T_", s)))
    cal <- function(s) g(paste0("n0t_", s)) / (g(paste0("n0t_", s)) + g(paste0("T_", s)))
    fdr_d <- fdr("added")
    fdr_o <- fdr("kept_both")
    cal_d <- cal("added")
    cal_o <- cal("kept_both")
    others <- g("planted") - g("down_all")
    out <- list(
        alpha_gadd = g("p05_gadd") / g("n0_gadd"), alpha_gadd_l0 = g("p05_gadd_l0") / g("n0_gadd_l0"), # KIT-281
        alpha_gadd_l1 = g("p05_gadd_l1") / g("n0_gadd_l1"), alpha_whole = g("p05_whole") / g("n0_whole"),
        fdr_bh = g("bhV_whole") / (g("bhV_whole") + g("bhT_whole")), gadd_per_r = g("n0_gadd") / g("one"),
        alpha_O = a_o, alpha_added = g("p05_added") / g("n0_added"), alpha_L = a_l, rho_l = a_l / a_o,
        rho_O = (g("p05_added") / g("n0_added")) / a_o,
        rho_l_common = (g("p05_common_l1") / g("n0_common_l1")) / a_o,
        Rc = (g("p05_whole_l1") / g("n0_whole_l1")) / (g("p05_whole_l0") / g("n0_whole_l0")),
        rho_l0 = (g("p05_whole0_l1") / g("n0_whole0_l1")) / (g("p05_whole0_l0") / g("n0_whole0_l0")),
        n0_L_per_r = g("n0_added_l1") / g("one"), pow = g("T_whole") / g("planted"), fdr_whole = fdr("whole"),
        fdp_mean = g("fdp_sum") / g("fdp_n"), kept_per_r = g("kept_whole") / g("one"),
        tested_per_r = g("tested_whole") / g("one"), calls_per_r = g("calls_whole") / g("one"),
        untested_per_r = (g("kept_whole") - g("tested_whole")) / g("one"),
        m6_l_added = g("m6sum_1") / g("m6n_1"), m6_na_added = g("m6sum_0") / g("m6n_0"),
        m6_l_H = g("m6Hsum_1") / g("m6Hn_1"), m6_na_H = g("m6Hsum_0") / g("m6Hn_0"),
        M5_whole = g("err_whole") / g("n1_whole"), A_whole = g("sgn_whole") / g("n1_whole"),
        M5_added = g("err_added") / g("n1_added"), A_added = g("sgn_added") / g("n1_added"),
        A_up = g("sgn_up") / g("n_up"), A_down = g("sgn_down") / g("n_down"),
        E_cal = (fdr_d - cal_d) - (fdr_o - cal_o), E_comp = cal_d - cal_o,
        alpha_band = g("p05b") / g("n0b"), pow_band = g("Tb") / g("n1b_all"), A_band = g("sgnb") / g("n1b"),
        pow_up = g("T_up") / others, pow_down = g("T_down") / g("down_all"),
        powbh = (g("bh_up") + g("bh_down")) / g("planted"), powbh_up = g("bh_up") / others,
        powbh_down = g("bh_down") / g("down_all"), onoff_n_per_r = g("onoff_all") / g("one"),
        onoff_removed = 1 - g("on_kept") / g("onoff_all"), onoff_pow = g("on_call") / g("onoff_all"),
        onoff_powbh = g("on_bh") / g("onoff_all"), calls_both_s_per_r = g("calls_both_s") / g("one"),
        new_calls_per_r = g("calls_new_s") / g("one"), gained_per_r = g("gained") / g("one"),
        lost_per_r = g("lost") / g("one"))
    for (st in STRATA) for (l in 0:1) {
        out[[paste0("st_rho_l_", st, "_", l)]] <- (g(paste0("sp05_", st, "_", l)) / g(paste0("sn0_", st, "_", l))) / a_o
        out[[paste0("st_n0_per_r_", st, "_", l)]] <- g(paste0("sn0_", st, "_", l)) / g("one")
    }
    out <- do.call(cbind, out)
    colnames(out) <- paste0(arm, ".", colnames(out))
    out
}

# KIT-281 (ADAPTED): KIT-211's :232-233 (C0.alpha_band, C0.rho_l0) and :240-245 (the Bp1 pairing) are removed; the gate's statistics
# are new. rhog* = B1's null p < 0.05 rate among the proteins H_B1 adds, over B's rate in H_B (the deciding quantity and its split by
# the flag ell); dfdrbh = BH false share of B1 minus B's (G1b); dpowbh = BH power of B minus B1's (G2); Phi as KIT-211 defined it.
computeStats <- function(S, arms) {
    st <- do.call(cbind, lapply(arms, \(arm) armStats(S, arm)))
    v <- function(name) st[, name]
    out <- list()
    if ("N" %in% arms) {
        for (x in c("S", CANDIDATES)) {
            out[[paste0("Phi.", x)]] <- v(paste0(x, ".rho_l")) / v("N.rho_l")
            out[[paste0("PhiH.", x)]] <- v(paste0(x, ".Rc")) / v("N.Rc")
            out[[paste0("drho.", x)]] <- v("S.rho_l") - v(paste0(x, ".rho_l"))
        }
        for (x in CANDIDATES) out[[paste0("share.", x)]] <- (v("S.rho_l") - v(paste0(x, ".rho_l"))) /
            (v("S.rho_l") - v("N.rho_l"))
    }
    for (x in setdiff(arms, "S")) {
        out[[paste0("dpow.", x)]] <- v("S.pow") - v(paste0(x, ".pow"))
        out[[paste0("dfdr.", x)]] <- v(paste0(x, ".fdr_whole")) - v("S.fdr_whole")
        out[[paste0("callshare.", x)]] <- v(paste0(x, ".calls_both_s_per_r")) / v("S.calls_per_r")
    }
    out[["rhog.B1"]] <- v("B1.alpha_gadd") / v("B.alpha_whole")
    out[["rhog_l0.B1"]] <- v("B1.alpha_gadd_l0") / v("B.alpha_whole")
    out[["rhog_l1.B1"]] <- v("B1.alpha_gadd_l1") / v("B.alpha_whole")
    out[["dfdrbh.B1"]] <- v("B1.fdr_bh") - v("B.fdr_bh")
    out[["dpowbh.B1"]] <- v("B.powbh") - v("B1.powbh")
    out[["rhoratio.B1"]] <- v("B1.rho_l") / v("B.rho_l")
    cbind(st, do.call(cbind, out))
}


# ---- resampling ----------------------------------------------------------------------------------------------------
# B, the seed and the draw order are declared (ticket section 5): after set.seed(280000001, ...) one B x R index matrix per cell in the
# order A, B, C, D, shared by every mask, statistic and arm of the cell. (KIT-281: SEED and CELLS_ALL above are the only changes.)
bootstrapWeights <- function(replicate_counts) {
    set.seed(SEED, kind = "Mersenne-Twister", normal.kind = "Inversion", sample.kind = "Rejection")
    cells <- CELLS_ALL[CELLS_ALL %in% names(replicate_counts)]
    lapply(stats::setNames(cells, cells), \(cellname) {
        R <- replicate_counts[[cellname]]
        idx <- matrix(sample.int(R, B_BOOT * R, replace = TRUE), nrow = B_BOOT, byrow = TRUE)
        t(apply(idx, 1L, tabulate, nbins = R))
    })
}

# Estimate, bootstrap and jackknife draws of every statistic for sums given per cell: `parts` is a list by cell of
# R x K matrices (each the sum over the masks that the setting or pool takes from that cell). A resample whose statistic
# is not finite (no tested unchanged protein in a subgroup) holds NA there: the limits are taken over the finite
# resamples and the count of the others is reported (boot_nonfinite).
resampleStats <- function(parts, weights, arms) {
    cells <- names(parts)
    total <- Reduce(`+`, lapply(parts, colSums))
    estimate <- computeStats(matrix(total, nrow = 1L, dimnames = list(NULL, names(total))), arms)
    boot <- computeStats(Reduce(`+`, lapply(cells, \(cellname) weights[[cellname]] %*% parts[[cellname]])), arms)
    jack <- computeStats(do.call(rbind, lapply(cells, \(cellname) matrix(total, nrow(parts[[cellname]]),
        length(total), byrow = TRUE, dimnames = list(NULL, names(total))) - parts[[cellname]])), arms)
    nonfinite <- colSums(!is.finite(boot))
    boot[!is.finite(boot)] <- NA_real_
    list(estimate = estimate, boot = boot, jack = jack, R = sum(vapply(parts, nrow, integer(1L))), nonfinite = nonfinite)
}

summariseStats <- function(x, scope, cellname, mask) {
    stats_names <- colnames(x$boot)
    data.frame(scope = scope, cell = cellname, mask = mask, stat = stats_names, R = x$R,
        estimate = as.numeric(x$estimate[1L, ]), boot_se = apply(x$boot, 2L, stats::sd, na.rm = TRUE),
        jack_se = apply(x$jack, 2L, jackSe), lo95 = quant(x$boot, 0.025), hi95 = quant(x$boot, 0.975),
        boot_nonfinite = as.numeric(x$nonfinite), stringsAsFactors = FALSE, row.names = NULL)
}

q1 <- function(boot, p) stats::quantile(boot, p, type = 7, na.rm = TRUE, names = FALSE)

# ---- the rules (6.4) -----------------------------------------------------------------------------------------------
# A one-sided state: `lower_rule` and `upper_rule` read the limits. Limits are one-sided at `alpha`: the lower limit is the
# alpha quantile, the upper limit the 1 - alpha quantile; the jackknife-normal limits beside are estimate -/+ z SE.
limitsOf <- function(x, stat, alpha) {
    boot <- x$boot[, stat]
    estimate <- unname(x$estimate[1L, stat])
    se <- jackSe(x$jack[, stat])
    z <- stats::qnorm(1 - alpha)
    c(estimate = estimate, lower = q1(boot, alpha), upper = q1(boot, 1 - alpha), jack_se = se,
        lower_jack = estimate - z * se, upper_jack = estimate + z * se)
}

fmt <- function(x, digits = 3) ifelse(is.na(x), "NA", formatC(x, digits = digits, format = "f"))

mdTable <- function(df) {
    cell <- function(x) gsub("\\|", "/", as.character(x))
    c(paste0("| ", paste(names(df), collapse = " | "), " |"), paste0("|", paste(rep("---", ncol(df)), collapse = "|"), "|"),
        vapply(seq_len(nrow(df)), \(i) paste0("| ", paste(vapply(df[i, , drop = TRUE], cell, character(1L)), collapse = " | "), " |"),
            character(1L)))
}

# ---- the scoring (KIT-281: new code; replaces KIT-211's :336-593) -----------------------------------------------------------------
ctx <- new.env()
ctx[["checks"]] <- list()
ctx[["counted"]] <- NA
ctx[["out"]] <- NA_character_

addCheck <- function(check, value, pass, detail) {
    ctx[["checks"]][[length(ctx[["checks"]]) + 1L]] <- data.frame(check = check, value = as.character(value), pass = pass, detail = detail,
        stringsAsFactors = FALSE)
}

progress <- function(unit, status) cat(sprintf("%s\t%s\n", unit, status), file = file.path(ctx[["out"]], "progress-score.tsv"), append = TRUE)

# A stop is an outcome of its own (ticket section 10, row 1): outcome.tsv says Stop with its kind and reason; the checks table so far is kept.
stopRecord <- function(out, kind, reason) {
    dir.create(out, recursive = TRUE, showWarnings = FALSE)
    if (length(ctx[["checks"]])) writeTsv(do.call(rbind, ctx[["checks"]]), file.path(out, "checks.tsv"))
    writeTsv(data.frame(outcome = "Stop", kind = kind, reason = reason, counted = as.character(ctx[["counted"]]), stringsAsFactors = FALSE),
        file.path(out, "outcome.tsv"))
    cat(sprintf("KIT-281 score: STOP (%s): %s\n", kind, reason))
    "stop"
}

fileSha <- function(path) digest::digest(file = path, algo = "sha256")

settingLabel <- function(lower, upper) if (isTRUE(upper <= RHO_MARGIN)) "within margin" else if (isTRUE(lower > RHO_MARGIN)) "exceeds margin" else
    "not determined"

# The deciding limit decides: a state is "pass" when the bootstrap upper limit is at most the limit, as written (no tolerance); a
# limit that straddles it, or one that is not finite, is "fail". The jackknife-normal state is printed beside and decides nothing.
decisionRow <- function(pool, id, quantity, stat, limit) {
    l <- limitsOf(pool, stat, ALPHA_PRIMARY)
    data.frame(id = id, quantity = quantity, stat = stat, R = pool[["R"]], estimate = l[["estimate"]], upper = l[["upper"]], limit = limit,
        state = if (isTRUE(l[["upper"]] <= limit)) "pass" else "fail", jack_se = l[["jack_se"]], upper_jack = l[["upper_jack"]],
        state_jack = if (isTRUE(l[["upper_jack"]] <= limit)) "pass" else "fail", boot_nonfinite = unname(pool[["nonfinite"]][[stat]]),
        stringsAsFactors = FALSE)
}

# F5' (ticket section 9): the fold-change bias of the unchanged proteins of H_B1 that hold a value at or below the threshold in a compared
# group and have missing values in one compared group only (KIT-211's m6_l on the kept set, signed to the side of the missing value);
# fires when the pooled two-sided 95 % interval over the twelve settings excludes 0.
falsifier5 <- function(pool) {
    stat <- "B1.m6_l_H"
    boot <- pool[["boot"]][, stat]
    estimate <- unname(pool[["estimate"]][1L, stat])
    se <- jackSe(pool[["jack"]][, stat])
    lower <- q1(boot, ALPHA_F_POOL / 2)
    upper <- q1(boot, 1 - ALPHA_F_POOL / 2)
    z <- stats::qnorm(1 - ALPHA_F_POOL / 2)
    data.frame(falsifier = "F5'", stat = stat, R = pool[["R"]], estimate = estimate, lower = lower, upper = upper, jack_se = se,
        lower_jack = estimate - z * se, upper_jack = estimate + z * se, fires = isTRUE(lower > 0 || upper < 0),
        fires_jack = isTRUE(estimate - z * se > 0 || estimate + z * se < 0), boot_nonfinite = unname(pool[["nonfinite"]][[stat]]),
        stringsAsFactors = FALSE)
}

scoreMain <- function(rows_root, out) {
    ctx[["out"]] <- out
    writeLines("unit\tstatus", file.path(out, "progress-score.tsv"))
    declared <- declaredRange()
    counted <- identical(declared[["source"]], "g3-table") && B_BOOT == 10000L && !identical(Sys.getenv("K281_GUARD_PROBE"), "1")
    ctx[["counted"]] <- counted
    script_file <- normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[[1L]]))
    input_files <- sort(list.files(rows_root, "\\.(tsv|gz|json|txt)$", recursive = TRUE, full.names = TRUE))
    writeProvenance <- function(status, outcome) jsonlite::write_json(list(ticket = "KIT-281", status = status, outcome = outcome,
        script_sha256 = fileSha(script_file), counted = counted, range_source = declared[["source"]],
        declared_range = lapply(declared[["range"]], range), g3_table_sha256 = declared[["sha"]],
        inputs_sha256 = as.list(stats::setNames(vapply(input_files, fileSha, character(1L)), sub(paste0(rows_root, "/"), "", input_files, fixed = TRUE))),
        bootstrap = list(B = B_BOOT, seed = SEED, kind = "Mersenne-Twister", normal_kind = "Inversion", sample_kind = "Rejection",
            unit = "cell replicate, all masks and arms", quantile_type = 7L), rho_margin = RHO_MARGIN, fdp_limit = FDP_LIMIT,
        guard_margin = GUARD_MARGIN, power_margin = POWER_MARGIN, pc_margin = PC_MARGIN,
        alpha = list(primary = ALPHA_PRIMARY, setting = ALPHA_SETTING, falsifier_pool = ALPHA_F_POOL), r_version = R.version.string),
        file.path(out, "provenance-score.json"), auto_unbox = TRUE, pretty = TRUE, digits = NA)
    writeProvenance("started", NA_character_) # before the first unit (script lesson 3)
    if (identical(Sys.getenv("K281_GUARD_PROBE"), "1")) { # the self-test's trigger: a deliberate `$` prefix match, which the guard must stop
        probe <- list(abcd = 1)
        invisible(probe$ab)
    }
    input <- readInput(rows_root, declared)
    metrics <- input[["metrics"]]
    ranges <- input[["ranges"]]
    ranges[["incomplete"]] <- ranges[["expected"]] != ranges[["complete"]]
    writeTsv(ranges, file.path(out, "ranges.tsv"))
    cells <- CELLS_ALL[CELLS_ALL %in% unique(metrics[["cell"]])]
    replicates <- lapply(stats::setNames(cells, cells), \(cellname) sort(unique(metrics[["r"]][metrics[["cell"]] == cellname])))
    writeTsv(data.frame(cell = cells, R = vapply(replicates, length, integer(1L)), first = vapply(replicates, min, numeric(1L)),
        last = vapply(replicates, max, numeric(1L)), stringsAsFactors = FALSE), file.path(out, "replicates.tsv"))
    # ---- the checks, read before any outcome ---------------------------------------------------------------------------------
    absent <- lapply(stats::setNames(cells, cells), \(cellname) setdiff(declared[["range"]][[cellname]], replicates[[cellname]]))
    r_each <- vapply(replicates, length, integer(1L))
    expected <- sum(r_each * 12L)
    duplicates <- anyDuplicated(metrics[c("cell", "mask", "r", "arm")]) > 0L
    arms_masks_ok <- all(vapply(cells, \(cellname) {
        m <- metrics[metrics[["cell"]] == cellname, , drop = FALSE]
        setequal(unique(m[["arm"]]), ARMS_MAIN) && setequal(unique(m[["mask"]]), CELL_MASKS[[cellname]])
    }, logical(1L)))
    complete <- setequal(cells, CELLS_ALL) && !any(lengths(absent) > 0L) && length(unique(r_each)) == 1L && nrow(metrics) == expected && !duplicates &&
        arms_masks_ok
    addCheck(check = "complete", value = nrow(metrics), pass = complete, detail = sprintf("%d rows of %d expected; cells %s; replicates absent from the declared range: %d; duplicates %s; every arm and mask present %s",
        nrow(metrics), expected, paste(cells, collapse = ""), sum(lengths(absent)), duplicates, arms_masks_ok))
    if (!complete) k211Abort("k211_instrument", sprintf("the rows are incomplete: a replicate is absent from the declared range, a cell, arm or mask is missing, or a row is duplicated (%d rows of %d expected; %d replicate(s) absent)",
        nrow(metrics), expected, sum(lengths(absent))))
    for (cellname in cells) if (anyDuplicated(input[["rows"]][[cellname]][c("mask", "r", "arm", "subset", "class", "stratum", "ell", "band", "ell0")]) > 0L)
        k211Abort("k211_instrument", sprintf("cell %s: a duplicated row key (a replicate joined twice)", cellname))
    addCheck(check = "ranges_incomplete", value = sum(ranges[["incomplete"]]), pass = !any(ranges[["incomplete"]]), detail = "ranges whose complete count is below the expected count")
    if (any(ranges[["incomplete"]])) k211Abort("k211_instrument", "a range is incomplete (a job did not complete; it is continued once, not scored around)")
    # Coordinator's ruling (S03): a Benjamini-Hochberg fallback of qMethod() is deterministic and the decision reads no q-value (raw p < 0.05 and BH on
    # raw p), so the count of rows whose recorded q-method is not qvalue is REPORTED (per cell, mask and arm in qmethod-fallback.tsv, in the checks table and
    # in outcome.tsv), never a stop. (The way a silent switch could arise, a guard inside package code, is kept out of the simulation.)
    n_fallback <- sum(metrics[["qmethod"]] != "qvalue")
    fallback_by <- stats::aggregate(list(n_rows = rep(1, nrow(metrics)), n_fallback = as.numeric(metrics[["qmethod"]] != "qvalue")),
        by = list(cell = metrics[["cell"]], mask = metrics[["mask"]], arm = metrics[["arm"]]), FUN = sum)
    writeTsv(fallback_by, file.path(out, "qmethod-fallback.tsv"))
    addCheck(check = "qmethod", value = n_fallback, pass = n_fallback == 0L,
        detail = "metrics rows whose recorded q-method is not qvalue (KIT-211: 0); reported, not a stop: the decision reads no q-value")
    det <- input[["determinism"]]
    det_fail <- if (is.null(det)) 0L else sum(!det[["rows_identical"]]) + sum(!det[["metrics_identical"]])
    addCheck(check = "determinism", value = if (is.null(det)) 0L else nrow(det), pass = !is.null(det) && det_fail == 0L,
        detail = "cell A, first five replicates of a range run twice, rows and metrics identical (required in a counted scoring)")
    if (counted && is.null(det)) k211Abort("k211_instrument", "no determinism table in a counted scoring")
    f1 <- sum(metrics[["viol"]])
    f2 <- sum(metrics[["e10_dev"]] > 0) + sum(!metrics[["e3_pass"]]) + det_fail
    addCheck(check = "F1", value = f1, pass = f1 == 0, detail = "identity violations of F1 (O in H_S, H_S in H_B, H_B in H_B1), summed")
    addCheck(check = "F2", value = f2, pass = f2 == 0, detail = "E3 failures (B1's DA input), logFC identities (E10) not exact, determinism differences, counted")
    addCheck(check = "t_r_monotone", value = sum(!metrics[["mono_ok"]]), pass = all(metrics[["mono_ok"]]), detail = "replicate-masks where the count of tested p <= t_r differs from the calls (reported, not a stop)")
    addCheck(check = "bootstrap_B", value = B_BOOT, pass = B_BOOT == 10000L, detail = "B is 10000 in a counted scoring; another value is a mechanics run (self-test or lifecycle)")
    addCheck(check = "range_declared", value = declared[["source"]], pass = identical(declared[["source"]], "g3-table"), detail = "the declared range is rule G3's table in a counted scoring; K281_RANGE is a mechanics run")
    if (f1 > 0 || f2 > 0) k211Abort("k211_falsifier", if (f1 > 0) "F1: a kept-set identity fails" else "F2: B1's input, an identity of logFC or the determinism comparison fails")
    # ---- atoms per cell and mask, then the bootstrap weights --------------------------------------------------------------------
    rows <- input[["rows"]]
    row_index <- lapply(stats::setNames(cells, cells), \(cellname) split(seq_len(nrow(rows[[cellname]])),
        paste(rows[[cellname]][["mask"]], rows[[cellname]][["arm"]])))
    cellAtoms <- function(cellname, mask) {
        mats <- lapply(ARMS_MAIN, \(arm) {
            a <- armAtoms(rows[[cellname]][row_index[[cellname]][[paste(mask, arm)]], , drop = FALSE],
                metrics[metrics[["cell"]] == cellname & metrics[["mask"]] == mask & metrics[["arm"]] == arm, , drop = FALSE],
                replicates[[cellname]])
            colnames(a) <- paste0(arm, ".", colnames(a))
            a
        })
        progress(unit = paste("atoms", cellname, mask), status = "done")
        do.call(cbind, mats)
    }
    atoms <- lapply(stats::setNames(cells, cells), \(cellname) lapply(stats::setNames(CELL_MASKS[[cellname]], CELL_MASKS[[cellname]]),
        \(mask) cellAtoms(cellname, mask)))
    weights <- bootstrapWeights(lapply(replicates, length))
    thin <- function(x) list(estimate = x[["estimate"]][, THIN, drop = FALSE], boot = x[["boot"]][, THIN, drop = FALSE],
        jack = x[["jack"]][, THIN, drop = FALSE], R = x[["R"]], nonfinite = x[["nonfinite"]][THIN])
    settings <- list()
    setting_tables <- list()
    for (cellname in cells) for (mask in CELL_MASKS[[cellname]]) {
        x <- resampleStats(stats::setNames(list(atoms[[cellname]][[mask]]), cellname), weights, ARMS_MAIN)
        setting_tables[[length(setting_tables) + 1L]] <- summariseStats(x, scope = "setting", cellname = cellname, mask = mask)
        settings[[paste(cellname, mask, sep = ":")]] <- thin(x)
        progress(unit = paste("setting", cellname, mask), status = "done")
    }
    design_stats <- lapply(stats::setNames(cells, cells), \(cellname) {
        x <- resampleStats(stats::setNames(list(Reduce(`+`, atoms[[cellname]])), cellname), weights, ARMS_MAIN)
        progress(unit = paste("design", cellname), status = "done")
        x
    })
    pool <- resampleStats(lapply(stats::setNames(names(POOL), names(POOL)), \(cellname) Reduce(`+`, lapply(POOL[[cellname]],
        \(mask) atoms[[cellname]][[mask]]))), weights, ARMS_MAIN)
    progress(unit = "pool twelve settings", status = "done")
    table_settings <- do.call(rbind, setting_tables)
    table_designs <- do.call(rbind, lapply(cells, \(cellname) summariseStats(design_stats[[cellname]], scope = "design", cellname = cellname, mask = "vd")))
    table_pool <- summariseStats(pool, scope = "pool", cellname = "twelve", mask = "vd")
    writeTsv(table_settings, file.path(out, "settings.tsv"))
    writeTsv(table_designs, file.path(out, "designs.tsv"))
    writeTsv(table_pool, file.path(out, "pooled.tsv"))
    # ---- the deciding quantity, the guards, F5', the reported items ----------------------------------------------------------------
    decision <- rbind(
        decisionRow(pool, id = "rho_g", quantity = "alpha_B1(H_B1 without H_B) / alpha_B(H_B), unchanged proteins", stat = "rhog.B1", limit = RHO_MARGIN),
        decisionRow(pool, id = "G1a", quantity = "FDP_BH(B1) = false BH calls / all BH calls at 0.05", stat = "B1.fdr_bh", limit = FDP_LIMIT),
        decisionRow(pool, id = "G1b", quantity = "FDP_BH(B1) minus FDP_BH(B), paired", stat = "dfdrbh.B1", limit = GUARD_MARGIN),
        decisionRow(pool, id = "G2", quantity = "BH power(B) minus BH power(B1), paired, over the 600 changed proteins", stat = "dpowbh.B1", limit = POWER_MARGIN))
    writeTsv(decision, file.path(out, "decision.tsv"))
    # Coordinator's ruling (S02): a NON-FINITE ESTIMATE of rho_g or of a guard is an instrument stop (nothing is claimed; the default stays 2), not
    # "Not calibrated". Resample-level non-finite values stay dropped and counted (boot_nonfinite).
    not_formed <- decision[["id"]][!is.finite(decision[["estimate"]])]
    if (length(not_formed)) k211Abort("k211_instrument", sprintf("the deciding quantity or a guard could not be formed (not finite): %s", paste(not_formed, collapse = ", ")))
    f5 <- falsifier5(pool)
    writeTsv(f5, file.path(out, "falsifier.tsv"))
    pc_limits <- limitsOf(pool, "Phi.S", ALPHA_PRIMARY)
    pc <- isTRUE(pc_limits[["lower"]] > PC_MARGIN)
    writeTsv(data.frame(stat = "Phi.S", R = pool[["R"]], estimate = pc_limits[["estimate"]], lower_975 = pc_limits[["lower"]], margin = PC_MARGIN,
        passes = pc, phi_B = unname(pool[["estimate"]][1L, "Phi.B"]), phi_B1 = unname(pool[["estimate"]][1L, "Phi.B1"]), stringsAsFactors = FALSE),
        file.path(out, "positive-control.tsv"))
    setting_labels <- do.call(rbind, lapply(names(settings), \(key) {
        parts <- strsplit(key, ":", fixed = TRUE)[[1L]]
        l <- limitsOf(settings[[key]], "rhog.B1", ALPHA_SETTING)
        data.frame(cell = parts[[1L]], mask = parts[[2L]], R = settings[[key]][["R"]], estimate = l[["estimate"]], lower = l[["lower"]],
            upper = l[["upper"]], label = settingLabel(l[["lower"]], l[["upper"]]), label_jack = settingLabel(l[["lower_jack"]], l[["upper_jack"]]),
            stringsAsFactors = FALSE)
    }))
    writeTsv(setting_labels, file.path(out, "labels-settings.tsv"))
    rho_pass <- identical(decision[["state"]][decision[["id"]] == "rho_g"], "pass")
    guards_pass <- all(decision[["state"]][decision[["id"]] != "rho_g"] == "pass")
    outcome <- if (f5[["fires"]]) "Withheld" else if (rho_pass && guards_pass) "Calibrated" else "Not calibrated"
    consequence <- switch(outcome, Withheld = "no claim on calibration; the default stays 2; the table returns to the owner",
        Calibrated = "Part C: the tab default becomes 1 (only after KIT-280 is on staging)",
        "the default stays 2; the result, with the per-design and by-flag tables, returns to the owner")
    writeTsv(data.frame(outcome = outcome, kind = "outcome", reason = consequence, counted = as.character(counted), bootstrap_B = B_BOOT,
        range_source = declared[["source"]], rho_g = decision[["estimate"]][decision[["id"]] == "rho_g"],
        rho_g_upper = decision[["upper"]][decision[["id"]] == "rho_g"], rho_g_state = decision[["state"]][decision[["id"]] == "rho_g"],
        g1a_state = decision[["state"]][decision[["id"]] == "G1a"], g1b_state = decision[["state"]][decision[["id"]] == "G1b"],
        g2_state = decision[["state"]][decision[["id"]] == "G2"], f5_fires = f5[["fires"]], pc_pass = pc, f1 = f1, f2 = f2,
        n_qmethod_fallback = n_fallback, notes = sprintf("%d metrics row(s) with a q-method other than qvalue (reported; the decision reads no q-value)", n_fallback),
        stringsAsFactors = FALSE), file.path(out, "outcome.tsv"))
    writeTsv(do.call(rbind, ctx[["checks"]]), file.path(out, "checks.tsv"))
    # ---- the decision table (generated; the owner edits the wording, VAL-41) ---------------------------------------------------------
    ivOf <- function(x, stat, digits = 3) sprintf("%s [%s, %s]", fmt(x[["estimate"]][1L, stat], digits), fmt(q1(x[["boot"]][, stat], 0.025), digits),
        fmt(q1(x[["boot"]][, stat], 0.975), digits))
    units <- c(design_stats, list(pooled = pool))
    per_design <- data.frame(Design = c(cells, "twelve settings"),
        `tested added null proteins per masked data set` = unname(vapply(units, \(x) fmt(x[["estimate"]][1L, "B1.gadd_per_r"], 1), character(1L))),
        `rho_g (two-sided 95 % limits)` = unname(vapply(units, ivOf, character(1L), stat = "rhog.B1")),
        `rho_g, flagged proteins (ell = 1)` = unname(vapply(units, ivOf, character(1L), stat = "rhog_l1.B1")),
        `rho_g, other proteins (ell = 0)` = unname(vapply(units, ivOf, character(1L), stat = "rhog_l0.B1")),
        `Phi, B (against no filter)` = unname(vapply(units, ivOf, character(1L), stat = "Phi.B")),
        `Phi, B1` = unname(vapply(units, ivOf, character(1L), stat = "Phi.B1")),
        `Phi, S (positive control)` = unname(vapply(units, ivOf, character(1L), stat = "Phi.S")), check.names = FALSE, stringsAsFactors = FALSE)
    arm_labels <- c(S = "S (shipped)", B = "B (2 groups)", B1 = "B1 (1 group)", N = "N (no filter)")
    writeTsv(table_pool[table_pool[["stat"]] %in% as.vector(outer(names(arm_labels), AG_STATS, paste, sep = ".")), , drop = FALSE], file.path(out, "whole-list.tsv"))
    whole_list <- data.frame(Arm = unname(arm_labels), kept = vapply(names(arm_labels), \(a) fmt(pool[["estimate"]][1L, paste0(a, ".kept_per_r")], 1), character(1L)),
        tested = vapply(names(arm_labels), \(a) fmt(pool[["estimate"]][1L, paste0(a, ".tested_per_r")], 1), character(1L)),
        `power by q < 0.05` = vapply(names(arm_labels), \(a) ivOf(pool, paste0(a, ".pow"), 4), character(1L)),
        `power by BH` = vapply(names(arm_labels), \(a) ivOf(pool, paste0(a, ".powbh"), 4), character(1L)),
        `false share by q < 0.05` = vapply(names(arm_labels), \(a) ivOf(pool, paste0(a, ".fdr_whole"), 4), character(1L)),
        `false share by BH` = vapply(names(arm_labels), \(a) ivOf(pool, paste0(a, ".fdr_bh"), 4), character(1L)), check.names = FALSE,
        stringsAsFactors = FALSE, row.names = NULL)
    deciding <- data.frame(Item = decision[["id"]], Quantity = decision[["quantity"]], Estimate = fmt(decision[["estimate"]], 4),
        `Upper limit (bootstrap, one-sided 97.5 %; decides)` = fmt(decision[["upper"]], 4), Limit = fmt(decision[["limit"]], 4), State = decision[["state"]],
        `Upper limit (jackknife-normal; reported)` = fmt(decision[["upper_jack"]], 4), `State (jackknife)` = decision[["state_jack"]],
        `Resamples dropped (not finite)` = decision[["boot_nonfinite"]], check.names = FALSE, stringsAsFactors = FALSE)
    design_estimates <- vapply(design_stats, \(x) x[["estimate"]][1L, "rhog.B1"], numeric(1L))
    claim <- sprintf("In simulated data, among unchanged proteins, the proteins that 1 group adds reached p < 0.05 %s times as often as the proteins kept at 2 groups, pooled over the twelve settings (the prespecified analysis; one-sided 97.5 %% upper limit %s; margin %.2f); per design the ratio was %s to %s.",
        fmt(decision[["estimate"]][decision[["id"]] == "rho_g"], 2), fmt(decision[["upper"]][decision[["id"]] == "rho_g"], 2), RHO_MARGIN,
        fmt(min(design_estimates), 2), fmt(max(design_estimates), 2))
    writeLines(c("# KIT-281 decision table (generated by score-k281.R; the owner edits the wording, VAL-41)", "",
        sprintf("Outcome: **%s**. Consequence: %s.", outcome, consequence),
        sprintf("Replicates per cell: %s. B = %d. Counted scoring: %s (declared range from %s).", paste(sprintf("%s %d", cells, r_each), collapse = ", "), B_BOOT,
            if (counted) "yes" else "NO (a mechanics run)", declared[["source"]]),
        sprintf("Positive control (S's one-sided 97.5 %% lower limit of Phi above %.2f; reported, decides and stops nothing): %s.", PC_MARGIN,
            if (pc) "passes" else "fails"),
        sprintf("Metrics rows whose recorded q-method is not qvalue: %d (reported; the decision reads raw p < 0.05 and BH on raw p).", n_fallback),
        sprintf("F5' (the pooled two-sided 95 %% interval of the fold-change bias of H_B1's flagged one-sided proteins excludes 0): %s.", if (f5[["fires"]]) "FIRES" else "does not fire"), "",
        "## The deciding quantity and the guards (the twelve settings pooled; ratios of sums)", "", mdTable(deciding), "",
        "## Per design, split by the flag, Phi and the positive control (reported, deciding nothing)", "", mdTable(per_design), "",
        "## Whole list, pooled", "", mdTable(whole_list), "", "## Sentence (the claim carries the per-design range; the sentence that 1 group leaves no excess is never written)", "", claim),
        file.path(out, "decision-table.md"))
    writeProvenance("complete", outcome)
    cat(sprintf("KIT-281 score: outcome %s\n", outcome))
    "ok"
}

# ---- driver --------------------------------------------------------------------------------------------------------------------------
args <- commandArgs(TRUE)
if (length(args) != 2L) k211Abort("k211_usage", "usage: score-k281.R <rows root> <out dir>")
rows_root_arg <- args[[1L]]
out_arg <- args[[2L]]
dir.create(out_arg, recursive = TRUE, showWarnings = FALSE)
# The handler is a calling handler inside the tryCatch, so the error it raises is caught below: a `$` prefix match is an instrument
# stop that ends the scorer with a recorded outcome, never an ordinary failure (script lesson 2).
status <- tryCatch(withCallingHandlers(scoreMain(rows_root_arg, out_arg), warning = function(w) {
    if (grepl("partial match of", tolower(conditionMessage(w)), fixed = TRUE))
        k211Abort("k211_instrument", paste("partial-match guard:", conditionMessage(w)))
}), k211_falsifier = function(e) stopRecord(out = out_arg, kind = "falsifier", reason = conditionMessage(e)),
k211_instrument = function(e) stopRecord(out = out_arg, kind = "instrument", reason = conditionMessage(e)),
error = function(e) stopRecord(out = out_arg, kind = "instrument", reason = paste("unclassified error:", conditionMessage(e))))
quit(status = if (identical(status, "ok")) 0L else 1L)
