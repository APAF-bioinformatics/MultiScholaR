# KIT-211 (VAL-142 (2), VAL-171, VAL-315): the scoring of the simulation's rows (ticket sections 6.3 to 6.9, 18.3). Evidence
# code under docs/qa: no package code, no package loaded.
#
# Usage, from the package root (any R with jsonlite and digest):
#   Rscript docs/qa/workflow-execution/kit/kit-211/score-k211.R <rows root> <out dir>
# <rows root> holds one directory per job (sim-k211.R's <cell>.tsv.gz, metrics.tsv, range.tsv); a cell's ranges (201 to
# 300 and 301 to 400) are joined into one cell of replicates before resampling (6.5); a range left incomplete leaves its
# cell scored on the replicates that completed, R stated beside every figure (range tables written).
# Extends docs/qa/workflow-execution/kit/kit-197/score-k197.R (SHA-256 f546947985617c5a6c0e79edd2d2f4a062db2ed0bea922b9927af6a1a9e9e426):
# the blocks between "verbatim" and "end verbatim" markers are KIT-197's lines unchanged (checked by e0-k211.py).
# Ratios are ratios of sums over the replicates of a setting (or of the ten deciding settings for the pool). Unit of
# resampling: a cell's replicate, with all its masks and arms (one index set per cell). Bootstrap: B = 10,000,
# set.seed(211000001, Mersenne-Twister, Inversion, Rejection) once, then for cells A, B, C, D, Av, then Ag in that order
# sample.int(R, B * R, replace = TRUE) filled by row into a B x R index matrix; percentile limits by stats::quantile(type = 7).
# Jackknife: delete one replicate of one cell with all its masks; for a setting the units are its cell's R replicates, for
# a pool the replicates of its cells. K211_B_BOOT (default 10000) is read only by the lifecycle's constructed tables, whose
# tables are degenerate (every replicate alike), and is recorded in checks.tsv (bootstrap_B).

CELLS_DECIDING <- c("A", "B", "C", "D")
CELLS_ALL <- c(CELLS_DECIDING, "Av", "Ag")
VD <- c("vd5", "vd10", "vd25")
MCAR <- c("mcar5", "mcar10", "mcar25")
CELL_MASKS <- list(A = c(VD, MCAR), B = c(VD, MCAR), C = c(VD, MCAR), D = c(VD, MCAR), Av = VD, Ag = VD)
ARMS_MAIN <- c("S", "A", "B", "Bp", "N", "C0")
ARMS_AG <- c("S", "A", "B", "Bp", "B1", "Bp1")
CANDIDATES <- c("A", "B", "Bp")
STRATA <- c("one_treated", "one_control", "both", "other_only", "none")
POOL <- list(A = VD, B = VD, C = VD, D = "vd25")
B_BOOT <- as.integer(Sys.getenv("K211_B_BOOT", "10000"))
SEED <- 211000001L
MARGIN <- 1.20
POWER_MARGIN <- 0.005
GUARD_MARGIN <- 0.0025
ALPHA_PRIMARY <- 0.025 / 3
ALPHA_PC <- 0.025
ALPHA_SETTING <- 0.05 / 72
ALPHA_F_POOL <- 0.05 / 4
ALPHA_F_SETTING <- 0.05 / 48
THIN <- c("Phi.S", "Phi.A", "Phi.B", "Phi.Bp", "S.m6_na_added", "A.m6_l_H", "B.m6_l_H", "Bp.m6_l_H")
AG_STATS <- c("pow", "fdr_whole", "powbh", "pow_down", "pow_up", "powbh_down", "powbh_up", "onoff_removed", "onoff_pow",
    "onoff_powbh", "kept_per_r", "tested_per_r", "calls_per_r")

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

# ---- input ---------------------------------------------------------------------------------------------------------
# Every job directory under <rows root> (A.tsv.gz ..., metrics.tsv, range.tsv): a cell's ranges are joined into one set
# of replicates (6.5); a range left incomplete leaves its cell scored on the replicates that completed.
readInput <- function(root) {
    dirs <- sort(list.dirs(root, recursive = FALSE))
    dirs <- dirs[file.exists(file.path(dirs, "metrics.tsv"))]
    if (!length(dirs)) k211Abort("k211_instrument", "no job directory with metrics.tsv under the rows root")
    row_classes <- c(cell = "character", mask = "character", arm = "character", subset = "character",
        class = "character", stratum = "character")
    metric_classes <- c(cell = "character", mask = "character", arm = "character", qmethod = "character",
        mono_ok = "logical", e3_pass = "logical")
    rows <- list()
    metrics <- list()
    ranges <- list()
    for (dir in dirs) {
        range <- utils::read.delim(file.path(dir, "range.tsv"), colClasses = c(cell = "character"))
        range$dir <- basename(dir)
        ranges[[length(ranges) + 1L]] <- range
        metrics[[length(metrics) + 1L]] <- utils::read.delim(file.path(dir, "metrics.tsv"), colClasses = metric_classes)
        for (cellname in range$cell) rows[[cellname]] <- c(rows[[cellname]], list(utils::read.delim(file.path(dir,
            paste0(cellname, ".tsv.gz")), colClasses = row_classes)))
    }
    list(rows = lapply(rows, \(x) do.call(rbind, x)), metrics = do.call(rbind, metrics), ranges = do.call(rbind, ranges))
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
        for (x in setdiff(arms, "C0")) out[[paste0("bandratio.", x)]] <- v(paste0(x, ".alpha_band")) / v("C0.alpha_band")
        out[["rho_l0.C0"]] <- v("C0.rho_l0")
    }
    for (x in setdiff(arms, "S")) {
        out[[paste0("dpow.", x)]] <- v("S.pow") - v(paste0(x, ".pow"))
        out[[paste0("dfdr.", x)]] <- v(paste0(x, ".fdr_whole")) - v("S.fdr_whole")
        out[[paste0("callshare.", x)]] <- v(paste0(x, ".calls_both_s_per_r")) / v("S.calls_per_r")
    }
    if ("B1" %in% arms) {
        for (s in AG_STATS) {
            for (x in setdiff(arms, "S")) out[[paste0("dS.", x, ".", s)]] <- v(paste0(x, ".", s)) - v(paste0("S.", s))
            for (b in c("B", "Bp")) out[[paste0("d1.", b, ".", s)]] <- v(paste0(b, "1.", s)) - v(paste0(b, ".", s))
        }
    }
    cbind(st, do.call(cbind, out))
}

# ---- resampling ----------------------------------------------------------------------------------------------------
# B, the seed and the draw order are declared (6.5, 7.1): after set.seed(211000001, ...) one B x R index matrix per cell
# in the order A, B, C, D, Av, then Ag (the cells present), shared by every mask, statistic and arm of the cell.
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

primaryState <- function(l, margin, null_value) {
    state <- function(lower, upper) if (isTRUE(upper <= margin)) "removes" else if (isTRUE(lower > null_value)) "excess remains" else
        "not determined"
    c(state = state(l[["lower"]], l[["upper"]]), state_jack = state(l[["lower_jack"]], l[["upper_jack"]]),
        note = if (isTRUE(l[["upper"]] <= margin) && isTRUE(l[["lower"]] > null_value)) "non-zero, within the margin" else "")
}

costState <- function(l, margin, high, low) {
    state <- function(lower, upper) if (isTRUE(lower > margin)) high else if (isTRUE(upper <= margin)) low else "not determined"
    c(state = state(l[["lower"]], l[["upper"]]), state_jack = state(l[["lower_jack"]], l[["upper_jack"]]))
}

classOf <- function(primary, power, guard, withheld) {
    if (withheld) return("U")
    if (identical(primary, "excess remains")) return("R")
    if (identical(primary, "removes")) return(if (identical(power, "no material cost") && identical(guard, "not raised")) "E" else "C")
    "U"
}

# The first row of 6.8 that holds. `pc` is TRUE when the positive control passes; `classes` the three candidates'.
outcomeOf <- function(stop_reason, pc, classes) {
    if (nzchar(stop_reason)) return("Stop")
    if (!pc) return("N0")
    e <- sum(classes == "E")
    if (e >= 2L) "K2" else if (e == 1L) "K1" else if (any(classes == "C")) "K3" else if (all(classes == "R")) "K4" else "K5"
}

fmt <- function(x, digits = 3) ifelse(is.na(x), "NA", formatC(x, digits = digits, format = "f"))

mdTable <- function(df) {
    cell <- function(x) gsub("\\|", "/", as.character(x))
    c(paste0("| ", paste(names(df), collapse = " | "), " |"), paste0("|", paste(rep("---", ncol(df)), collapse = "|"), "|"),
        vapply(seq_len(nrow(df)), \(i) paste0("| ", paste(vapply(df[i, , drop = TRUE], cell, character(1L)), collapse = " | "), " |"),
            character(1L)))
}

# ---- main ----------------------------------------------------------------------------------------------------------
args <- commandArgs(TRUE)
if (length(args) != 2L) k211Abort("k211_usage", "usage: score-k211.R <rows root> <out dir>")
rows_root <- args[[1L]]
out <- args[[2L]]
dir.create(out, recursive = TRUE, showWarnings = FALSE)
input <- readInput(rows_root)
metrics <- input$metrics
cells <- CELLS_ALL[CELLS_ALL %in% unique(metrics$cell)]
if (!all(CELLS_DECIDING %in% cells)) k211Abort("k211_instrument", "cells A to D are required")
armsOf <- function(cellname) if (identical(cellname, "Ag")) ARMS_AG else ARMS_MAIN
replicates <- lapply(stats::setNames(cells, cells), \(cellname) sort(unique(metrics$r[metrics$cell == cellname])))
expected <- sum(vapply(cells, \(cellname) length(replicates[[cellname]]) * length(CELL_MASKS[[cellname]]) *
    length(armsOf(cellname)), numeric(1L)))
arms_masks_ok <- all(vapply(cells, \(cellname) {
    m <- metrics[metrics$cell == cellname, , drop = FALSE]
    setequal(unique(m$arm), armsOf(cellname)) && setequal(unique(m$mask), CELL_MASKS[[cellname]])
}, logical(1L)))
ranges <- input$ranges
ranges$incomplete <- ranges$expected != ranges$complete
writeTsv(ranges, file.path(out, "ranges.tsv"))
writeTsv(data.frame(cell = cells, R = vapply(replicates, length, integer(1L)), stringsAsFactors = FALSE),
    file.path(out, "replicates.tsv"))
f1 <- sum(metrics$viol)
f2 <- sum(metrics$e10_dev > 0) + sum(!metrics$e3_pass)
checks <- data.frame(check = c("complete", "ranges_incomplete", "F1", "F2", "t_r_monotone", "bh_fallback", "bootstrap_B"),
    value = c(nrow(metrics), sum(ranges$incomplete), f1, f2, sum(!metrics$mono_ok), sum(metrics$qmethod != "qvalue"), B_BOOT),
    pass = c(nrow(metrics) == expected && !anyDuplicated(metrics[c("cell", "mask", "r", "arm")]) && arms_masks_ok, TRUE,
        f1 == 0, f2 == 0, all(metrics$mono_ok), TRUE, B_BOOT == 10000L),
    detail = c(sprintf("%d rows of %d expected", nrow(metrics), expected),
        "ranges incomplete after their one continuation (a cell is scored on the replicates that completed; R beside every figure)",
        "identity violations of F1, summed", "E3 failures and logFC identities (E10) not exact, counted",
        "replicate-masks where the count of tested p <= t_r differs from the calls",
        "replicate-masks whose q-values fell back to BH (counted, not a stop)",
        "B is 10000 in a counted scoring; another value is a lifecycle mechanics run"), stringsAsFactors = FALSE)
det_files <- file.path(sort(list.dirs(rows_root, recursive = FALSE)), "determinism.tsv")
det_files <- det_files[file.exists(det_files)]
for (f in det_files) {
    det <- utils::read.delim(f, stringsAsFactors = FALSE)
    checks <- rbind(checks, data.frame(check = "determinism", value = nrow(det), pass = all(det$rows_identical) &&
        all(det$metrics_identical), detail = "cell A, first five replicates of the range run twice", stringsAsFactors = FALSE))
}
if (!checks$pass[checks$check == "complete"]) {
    writeTsv(checks, file.path(out, "checks.tsv"))
    k211Abort("k211_instrument", "the rows are incomplete")
}
if (f1 > 0 || f2 > 0) {
    writeTsv(checks, file.path(out, "checks.tsv"))
    writeTsv(data.frame(outcome = "Stop", reason = if (f1 > 0) "F1: a kept-set identity fails" else
        "F2: A's input or an identity of logFC fails", stringsAsFactors = FALSE), file.path(out, "outcome.tsv"))
    k211Abort("k211_falsifier", if (f1 > 0) "F1: a kept-set identity fails" else "F2: A's input or an identity of logFC fails")
}

# Atoms per cell and mask, then the bootstrap weights.
row_index <- lapply(stats::setNames(cells, cells), \(cellname) split(seq_len(nrow(input$rows[[cellname]])),
    paste(input$rows[[cellname]]$mask, input$rows[[cellname]]$arm)))
cellAtoms <- function(cellname, mask) {
    mats <- lapply(armsOf(cellname), \(arm) {
        a <- armAtoms(input$rows[[cellname]][row_index[[cellname]][[paste(mask, arm)]], , drop = FALSE],
            metrics[metrics$cell == cellname & metrics$mask == mask & metrics$arm == arm, , drop = FALSE],
            replicates[[cellname]])
        colnames(a) <- paste0(arm, ".", colnames(a))
        a
    })
    do.call(cbind, mats)
}
atoms <- lapply(stats::setNames(cells, cells), \(cellname) lapply(stats::setNames(CELL_MASKS[[cellname]],
    CELL_MASKS[[cellname]]), \(mask) cellAtoms(cellname, mask)))
weights <- bootstrapWeights(lapply(replicates, length))
thin <- function(x) list(estimate = x$estimate[, THIN, drop = FALSE], boot = x$boot[, THIN, drop = FALSE],
    jack = x$jack[, THIN, drop = FALSE], R = x$R, nonfinite = x$nonfinite[THIN])

settings <- list()
setting_tables <- list()
for (cellname in cells) for (mask in CELL_MASKS[[cellname]]) {
    x <- resampleStats(stats::setNames(list(atoms[[cellname]][[mask]]), cellname), weights, armsOf(cellname))
    setting_tables[[length(setting_tables) + 1L]] <- summariseStats(x, if (cellname %in% CELLS_DECIDING) "setting" else
        "reported", cellname, mask)
    if (cellname %in% CELLS_DECIDING) settings[[paste(cellname, mask, sep = ":")]] <- thin(x)
}
poolParts <- function(masks_of) lapply(stats::setNames(names(POOL), names(POOL)), \(cellname)
    Reduce(`+`, lapply(masks_of(POOL[[cellname]]), \(mask) atoms[[cellname]][[mask]])))
pool_vd <- resampleStats(poolParts(identity), weights, ARMS_MAIN)
pool_mcar <- resampleStats(poolParts(\(m) sub("^vd", "mcar", m)), weights, ARMS_MAIN)
pool_24 <- resampleStats(lapply(stats::setNames(CELLS_DECIDING, CELLS_DECIDING), \(cellname) Reduce(`+`,
    atoms[[cellname]])), weights, ARMS_MAIN)
table_settings <- do.call(rbind, setting_tables)
table_pools <- rbind(summariseStats(pool_vd, "pool", "ten", "vd"), summariseStats(pool_mcar, "pool", "ten", "mcar"),
    summariseStats(pool_24, "pool", "all24", "all"))
writeTsv(table_settings, file.path(out, "settings.tsv"))
writeTsv(table_pools, file.path(out, "pooled.tsv"))

# ---- falsifiers (6.7): the pool of the 24 settings at 0.05/4, the per-setting Bonferroni intervals beside ----------------
falsifier <- function(id, stat) {
    boot <- pool_24$boot[, stat]
    estimate <- unname(pool_24$estimate[1L, stat])
    se <- jackSe(pool_24$jack[, stat])
    lower <- q1(boot, ALPHA_F_POOL / 2)
    upper <- q1(boot, 1 - ALPHA_F_POOL / 2)
    z <- stats::qnorm(1 - ALPHA_F_POOL / 2)
    per_setting <- vapply(names(settings), \(key) {
        b <- settings[[key]]$boot[, stat]
        c(q1(b, ALPHA_F_SETTING), q1(b, 1 - ALPHA_F_SETTING))
    }, numeric(2L))
    data.frame(falsifier = id, stat = stat, R = pool_24$R, estimate = estimate, lower = lower, upper = upper,
        jack_se = se, lower_jack = estimate - z * se, upper_jack = estimate + z * se,
        fires = isTRUE(lower > 0 || upper < 0), fires_jack = isTRUE(estimate - z * se > 0 || estimate + z * se < 0),
        settings_excluding_zero = sum(per_setting[1L, ] > 0 | per_setting[2L, ] < 0, na.rm = TRUE),
        stringsAsFactors = FALSE)
}
fals <- rbind(falsifier("F3", "A.m6_l_H"), falsifier("F4", "S.m6_na_added"), falsifier("F5_B", "B.m6_l_H"),
    falsifier("F5_Bp", "Bp.m6_l_H"))
writeTsv(fals, file.path(out, "falsifiers.tsv"))
withheld <- c(A = fals$fires[fals$falsifier == "F3"], B = fals$fires[fals$falsifier == "F5_B"],
    Bp = fals$fires[fals$falsifier == "F5_Bp"])
withheld_by <- c(A = "F3", B = "F5", Bp = "F5")
f4 <- fals$fires[fals$falsifier == "F4"]

# ---- the primary decision, the positive control, the classes (6.4) -------------------------------------------------------
pc_limits <- limitsOf(pool_vd, "Phi.S", ALPHA_PC)
pc <- isTRUE(pc_limits[["lower"]] > MARGIN)
candidate_rows <- list()
for (x in CANDIDATES) {
    l <- limitsOf(pool_vd, paste0("Phi.", x), ALPHA_PRIMARY)
    ps <- primaryState(l, MARGIN, 1)
    lp <- limitsOf(pool_vd, paste0("dpow.", x), ALPHA_PRIMARY)
    pw <- costState(lp, POWER_MARGIN, "costs power", "no material cost")
    lg <- limitsOf(pool_vd, paste0("dfdr.", x), ALPHA_PRIMARY)
    gd <- costState(lg, GUARD_MARGIN, "raised", "not raised")
    cls <- classOf(ps[["state"]], pw[["state"]], gd[["state"]], withheld[[x]])
    candidate_rows[[x]] <- data.frame(arm = x, R = pool_vd$R, phi = l[["estimate"]], phi_lower = l[["lower"]],
        phi_upper = l[["upper"]], phi_jack_se = l[["jack_se"]], primary = if (withheld[[x]]) paste0("withheld (",
            withheld_by[[x]], ")") else ps[["state"]], primary_jack = ps[["state_jack"]], primary_note = ps[["note"]],
        primary_differs_jack = ps[["state"]] != ps[["state_jack"]], dpow = lp[["estimate"]], dpow_lower = lp[["lower"]],
        dpow_upper = lp[["upper"]], power = if (withheld[[x]]) paste0("withheld (", withheld_by[[x]], ")") else
            pw[["state"]], power_jack = pw[["state_jack"]], dfdr = lg[["estimate"]], dfdr_lower = lg[["lower"]],
        dfdr_upper = lg[["upper"]], guard = if (withheld[[x]]) paste0("withheld (", withheld_by[[x]], ")") else
            gd[["state"]], guard_jack = gd[["state_jack"]], class = cls, withheld = withheld[[x]],
        discriminates_from_S = pc, stringsAsFactors = FALSE)
}
labels <- do.call(rbind, candidate_rows)
labels$pc_phi_S <- pc_limits[["estimate"]]
labels$pc_lower_975 <- pc_limits[["lower"]]
labels$pc_pass <- pc
writeTsv(labels, file.path(out, "labels.tsv"))

setting_labels <- list()
for (key in names(settings)) for (x in CANDIDATES) {
    parts <- strsplit(key, ":", fixed = TRUE)[[1L]]
    l <- limitsOf(settings[[key]], paste0("Phi.", x), ALPHA_SETTING)
    ps <- primaryState(l, MARGIN, 1)
    deciding <- parts[[2L]] %in% VD && !(parts[[1L]] == "D" && parts[[2L]] %in% c("vd5", "vd10"))
    setting_labels[[length(setting_labels) + 1L]] <- data.frame(cell = parts[[1L]], mask = parts[[2L]], arm = x,
        R = settings[[key]]$R, phi = l[["estimate"]], lower = l[["lower"]], upper = l[["upper"]], state = ps[["state"]],
        state_jack = ps[["state_jack"]], deciding = deciding, differs_from_pooled = deciding &&
            ps[["state"]] != labels$primary[labels$arm == x], stringsAsFactors = FALSE)
}
setting_labels <- do.call(rbind, setting_labels)
writeTsv(setting_labels, file.path(out, "labels-settings.tsv"))

stop_reason <- if (f4) "F4" else ""
classes <- stats::setNames(labels$class, labels$arm)
outcome <- outcomeOf(stop_reason, pc, classes)
writeTsv(data.frame(outcome = outcome, pc_pass = pc, f1 = f1 > 0, f2 = f2 > 0, f4 = f4, class_A = classes[["A"]],
    class_B = classes[["B"]], class_Bp = classes[["Bp"]], withheld = paste(names(withheld)[withheld], collapse = ";"),
    bootstrap_B = B_BOOT, stringsAsFactors = FALSE), file.path(out, "outcome.tsv"))
writeTsv(checks, file.path(out, "checks.tsv"))

# ---- the decision table (6.9), numeric rows; the fixed rows are the ticket's text -----------------------------------------
est <- function(x, stat) x$estimate[1L, stat]
pv <- function(stat) est(pool_vd, stat)
lo <- function(stat) q1(pool_vd$boot[, stat], 0.025)
hi <- function(stat) q1(pool_vd$boot[, stat], 0.975)
iv <- function(stat, digits = 3) sprintf("%s [%s, %s]", fmt(pv(stat), digits), fmt(lo(stat), digits), fmt(hi(stat), digits))
arm_labels <- c(S = "S (shipped)", A = "A", B = "B", Bp = "B'", N = "N (no filter)", C0 = "C0 (complete data)")
subgroup <- data.frame(Arm = arm_labels, `Phi (95 % limits)` = vapply(names(arm_labels), \(a) if (paste0("Phi.", a) %in%
    colnames(pool_vd$boot)) iv(paste0("Phi.", a)) else "", character(1L)), Class = vapply(names(arm_labels), \(a)
    if (a %in% CANDIDATES) paste0(labels$primary[labels$arm == a], " / ", classes[[a]]) else if (a == "S")
        paste0("PC ", if (pc) "passes" else "fails") else "", character(1L)),
    `rho_l against O` = vapply(names(arm_labels), \(a) iv(paste0(a, ".rho_l")), character(1L)),
    `rho_l0 (C0)` = vapply(names(arm_labels), \(a) if (a == "C0") iv("C0.rho_l0") else "", character(1L)),
    `share of S's excess removed` = vapply(names(arm_labels), \(a) if (a %in% CANDIDATES) iv(paste0("share.", a)) else "",
        character(1L)), `PhiH` = vapply(names(arm_labels), \(a) if (paste0("PhiH.", a) %in% colnames(pool_vd$boot))
        iv(paste0("PhiH.", a)) else "", character(1L)), check.names = FALSE, stringsAsFactors = FALSE)
nearby <- data.frame(Arm = arm_labels, `alpha, true level within 1 log2 above the threshold` = vapply(names(arm_labels),
    \(a) iv(paste0(a, ".alpha_band"), 4), character(1L)), `ratio to C0` = vapply(names(arm_labels), \(a)
    if (a == "C0") "1" else iv(paste0("bandratio.", a)), character(1L)), `power there` = vapply(names(arm_labels),
    \(a) iv(paste0(a, ".pow_band"), 4), character(1L)), `attenuation there` = vapply(names(arm_labels), \(a)
    iv(paste0(a, ".A_band"), 4), character(1L)), check.names = FALSE, stringsAsFactors = FALSE)
whole <- data.frame(Arm = arm_labels, kept = vapply(names(arm_labels), \(a) fmt(pv(paste0(a, ".kept_per_r")), 1), character(1L)),
    tested = vapply(names(arm_labels), \(a) fmt(pv(paste0(a, ".tested_per_r")), 1), character(1L)),
    called = vapply(names(arm_labels), \(a) fmt(pv(paste0(a, ".calls_per_r")), 2), character(1L)),
    power = vapply(names(arm_labels), \(a) iv(paste0(a, ".pow"), 4), character(1L)),
    `power change against S (S minus X)` = vapply(names(arm_labels), \(a) if (paste0("dpow.", a) %in%
        colnames(pool_vd$boot)) iv(paste0("dpow.", a), 4) else "", character(1L)),
    `false share` = vapply(names(arm_labels), \(a) iv(paste0(a, ".fdr_whole"), 4), character(1L)),
    `false share change (X minus S)` = vapply(names(arm_labels), \(a) if (paste0("dfdr.", a) %in%
        colnames(pool_vd$boot)) iv(paste0("dfdr.", a), 4) else "", character(1L)), check.names = FALSE,
    stringsAsFactors = FALSE)
folds <- data.frame(Arm = arm_labels, `m6_l on the added set` = vapply(names(arm_labels), \(a) iv(paste0(a, ".m6_l_added"), 4),
    character(1L)), `m6_l on the kept set` = vapply(names(arm_labels), \(a) iv(paste0(a, ".m6_l_H"), 4), character(1L)),
    `attenuation, whole list` = vapply(names(arm_labels), \(a) iv(paste0(a, ".A_whole"), 4), character(1L)),
    `attenuation, added set` = vapply(names(arm_labels), \(a) iv(paste0(a, ".A_added"), 4), character(1L)),
    check.names = FALSE, stringsAsFactors = FALSE)
changes <- data.frame(Arm = arm_labels[c("A", "B", "Bp", "N")], `share of S's calls also called` = vapply(c("A", "B", "Bp", "N"),
    \(a) iv(paste0("callshare.", a), 4), character(1L)), `new calls per data set` = vapply(c("A", "B", "Bp", "N"), \(a)
    iv(paste0(a, ".new_calls_per_r"), 2), character(1L)), `proteins gained per data set` = vapply(c("A", "B", "Bp", "N"),
    \(a) iv(paste0(a, ".gained_per_r"), 1), character(1L)), `proteins lost per data set` = vapply(c("A", "B", "Bp", "N"),
    \(a) iv(paste0(a, ".lost_per_r"), 1), character(1L)), `kept but not tested in the contrast` = vapply(c("A", "B", "Bp", "N"),
    \(a) iv(paste0(a, ".untested_per_r"), 1), character(1L)), check.names = FALSE, stringsAsFactors = FALSE)
reported <- function(cellname, stats_wanted) {
    t <- table_settings[table_settings$cell == cellname & table_settings$stat %in% stats_wanted, , drop = FALSE]
    t$value <- sprintf("%s [%s, %s]", fmt(t$estimate, 4), fmt(t$lo95, 4), fmt(t$hi95, 4))
    wide <- stats::reshape(t[c("mask", "stat", "value")], idvar = "stat", timevar = "mask", direction = "wide")
    names(wide) <- sub("^value\\.", "", names(wide))
    rownames(wide) <- NULL
    wide
}
av_stats <- c(paste0(rep(c("S", "A", "B", "Bp", "N"), each = 1), ".rho_l"), paste0("Phi.", c("S", "A", "B", "Bp")),
    paste0(c("S", "A", "B", "Bp", "N", "C0"), ".alpha_band"), paste0(c("S", "A", "B", "Bp"), ".pow"),
    paste0(c("S", "A", "B", "Bp"), ".fdr_whole"))
ag_stats <- c(paste0(rep(ARMS_AG, each = 4), ".", c("kept_per_r", "pow", "powbh", "fdr_whole")),
    paste0(rep(ARMS_AG, each = 4), ".", c("pow_down", "powbh_down", "onoff_removed", "onoff_powbh")),
    paste0("d1.", c("B", "Bp"), ".pow"), paste0("d1.", c("B", "Bp"), ".powbh_down"), paste0("d1.", c("B", "Bp"),
        ".onoff_powbh"), paste0("d1.", c("B", "Bp"), ".fdr_whole"))
fixed_rows <- c(
    "| Routes affected | DIA without limpa, LFQ, TMT; the legacy limma test (with and without TREAT) and the declared designs; not DIA with limpa; strict mode: A changes nothing (strict mode keeps only proteins with no value at or below the threshold); B and B' change it to \"no NA\" |",
    "| What users are told and see | A: the lowest values of the table (at most `percentile / 100` of them) are left out of the test; the test treats them as absent although they are known to be low (left-censored at the threshold, not missing at random), so fold changes of changed proteins near the threshold are biased toward 0; the share left out of the test is the percentile itself (1 % at the default, up to 10 % at the tab's maximum), and its costs near the threshold grow with it; only the default is measured. B: the percentile control no longer does anything at this step. B': the percentile becomes a label-blind low-intensity rule as in lipidomics and metabolomics, with its proportion (0.5); B' is the rule lipidomics and metabolomics ship after KIT-217, with proteomics' minimums. At the default percentile B' removes only proteins with at least half of all samples below the table's 1st percentile, so it keeps nearly B's proteins. The choice between B and B' is whether the percentile control keeps a meaning, not a measured difference. S: today's help text |",
    "| Flag column (SD8; decided with the direction) | whether results tables mark proteins kept with a value at or below the threshold; how many per data set: the kept counts of the first table of 'Whole list' and the added-set counts per arm |",
    "| Code change it needs | A: NA at the DA input on both DA entry points, a record field, a replay marker; B: the presence rule in one helper, the control, record and texts, a replay marker; B': B plus the share rule and its proportion parameter |",
    "| Measured where | the simulated design only: one generator whose variance does not depend on intensity (in real label-free and DIA data the lowest intensities are the noisiest, which changes both what A discards and how large the subgroup's excess is under B); TREAT with the trend and robust settings of the bound commit; equal groups of 3 or 4; no normalisation or RUV between the filter and the test; the default percentile and minimums; strict mode not run. Which direction is better on real data is not measured. One reported cell (Av) with intensity-dependent variance shows the direction of that effect in one design |")
writeLines(c("# KIT-211 decision table (generated by score-k211.R; the owner edits the wording, VAL-41)", "",
    sprintf("Outcome: **%s**. Positive control (S's lower 97.5 %% limit of Phi above %.2f): %s.", outcome, MARGIN,
        if (pc) "passes" else "fails (each candidate's primary state is printed, marked 'does not discriminate from S')"),
    sprintf("Replicates: %s. B = %d.", paste(sprintf("%s %d", cells, vapply(replicates, length, integer(1L))), collapse = ", "),
        B_BOOT), if (any(withheld)) sprintf("Withheld: %s.", paste(sprintf("%s (%s)", names(withheld)[withheld],
            withheld_by[names(withheld)[withheld]]), collapse = ", ")) else "No candidate withheld.", "",
    "## Subgroup (the ten deciding settings, vd; ratios of sums; limits two-sided 95 %)", "", mdTable(subgroup), "",
    "## Near the threshold", "", mdTable(nearby), "", "## Whole list", "", mdTable(whole), "", "## Fold change", "",
    mdTable(folds), "", "## What changes in users' results", "", mdTable(changes), "",
    "## Intensity-dependent variance (Av, reported only, no label)", "", mdTable(reported("Av", av_stats)), "",
    "## Groups required (Ag, reported only, no label): 2 groups against 1 under B and B' (VAL-204 (2), (7); VAL-315)", "",
    mdTable(reported("Ag", ag_stats)), "", "## Fixed rows", "", "| Column | Content |", "|---|---|", fixed_rows),
    file.path(out, "decision-table.md"))

provenance <- list(ticket = "KIT-211", script_sha256 = digest::digest(file = normalizePath(sub("^--file=", "",
    grep("^--file=", commandArgs(FALSE), value = TRUE)[[1L]])), algo = "sha256"),
    inputs_sha256 = as.list(stats::setNames(vapply(list.files(rows_root, "\\.(tsv|gz)$", recursive = TRUE, full.names = TRUE),
        \(f) digest::digest(file = f, algo = "sha256"), character(1L)), sub(paste0(rows_root, "/"), "", list.files(rows_root,
            "\\.(tsv|gz)$", recursive = TRUE, full.names = TRUE), fixed = TRUE))),
    replicates = lapply(replicates, range), replicate_counts = lapply(replicates, length),
    bootstrap = list(B = B_BOOT, seed = SEED, kind = "Mersenne-Twister", normal_kind = "Inversion",
        sample_kind = "Rejection", unit = "cell replicate, all masks and arms", quantile_type = 7L,
        draw_order = names(weights)), margin = MARGIN, power_margin = POWER_MARGIN, guard_margin = GUARD_MARGIN,
    alpha = list(primary = ALPHA_PRIMARY, pc = ALPHA_PC, setting = ALPHA_SETTING, falsifier_pool = ALPHA_F_POOL,
        falsifier_setting = ALPHA_F_SETTING), r_version = R.version.string)
jsonlite::write_json(provenance, file.path(out, "provenance-score.json"), auto_unbox = TRUE, pretty = TRUE, digits = NA)
cat(sprintf("KIT-211 score: outcome %s\n", outcome))
