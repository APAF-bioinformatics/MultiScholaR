# KIT-197 (VAL-111, VAL-121): the scoring of the simulation's rows (ticket sections 5.3 to 5.7). Evidence code under
# docs/qa: no package code, no package loaded.
#
# Usage, from the package root (any R with jsonlite and digest):
#   Rscript docs/qa/workflow-execution/kit/kit-197/score-k197.R descriptive <rows dir> <S2 cells.tsv> <out dir>
#   Rscript docs/qa/workflow-execution/kit/kit-197/score-k197.R deciding <rows dir> <out dir>
# <rows dir> holds sim-k197.R's A.tsv.gz to D.tsv.gz and s2metrics.tsv for one replicate range.
# descriptive (S_a, r = 1 to 100, S2's draws): gate G1 first (S2's figures reproduced to a relative 1e-9; on failure
#   g1-mismatch.tsv is written and the script stops with class k197_g1, computing no estimand), then F1, then every
#   estimand of 5.3 per setting and pooled, with intervals and jackknife SEs and no label.
# deciding (S_b, r = 101 to 200, the new draws): F1, every estimand, the labels of 5.6 per setting (vd and mcar) and
#   pooled, F4, the attribution, the outcome (O1 to O5 or Stop) and the sentence of 5.7.
# Ratios are ratios of sums over the replicates of a setting (or of the ten deciding settings for the pool). Unit of
# resampling: a cell's replicate, with all six masks (one index set per cell, shared by every mask and statistic;
# R1-01). Bootstrap: B = 10,000, set.seed(197000001, Mersenne-Twister, Inversion, Rejection) once, then for cells A,
# B, C, D in that order sample.int(R, B * R, replace = TRUE) filled by row into a B x R index matrix; percentile
# limits by stats::quantile(type = 7). Jackknife: delete one replicate of one cell with all its masks; for a setting
# the units are its cell's R replicates, for the pool the 4 R replicates of the four cells.

CELLS <- c("A", "B", "C", "D")
VD <- c("vd5", "vd10", "vd25")
MCAR <- c("mcar5", "mcar10", "mcar25")
SUBSETS <- c("added", "kept_both", "ref", "whole")
STRATA <- c("one_treated", "one_control", "both", "other_only", "none")
POOL <- list(A = c("vd5", "vd10", "vd25"), B = c("vd5", "vd10", "vd25"), C = c("vd5", "vd10", "vd25"), D = "vd25")
B_BOOT <- 10000L
SEED <- 197000001L
NULL_VALUE <- c(B = 1, P = 0)
MARGIN <- c(B = 1.10, P = 0.005)
COMPONENT_STAT <- c(B = "rho_O", P = "E_comp")
LEVEL_PRESENT <- c(setting = 0.05 / 12, pool = 0.025)
LEVEL_ABSENT <- 0.05
LEVEL_F4 <- 0.05 / 24
G1_TOLERANCE <- 1e-9
R5_STATS <- c("rho_O", "E_cal", "E_comp", "m6_na", "m6_ell", "A_added")
G1_MAP <- data.frame(metric = c("kept_o", "kept_h", "added", "calls", "power", "fdp", "added_share_called", "added_fdp"),
    side = c("base", "head", "head-base", "head", "head", "head", "head", "head"),
    s2_metric = c("kept", "kept", "added", "calls", "power", "fdp", "added_share_called", "added_fdp"),
    stringsAsFactors = FALSE)

k197Abort <- function(class, message) stop(structure(class = c(class, "error", "condition"),
    list(message = message, call = NULL)))

num <- function(x) ifelse(is.na(x), "NA", sprintf("%.17g", x))

writeTsv <- function(df, path) {
    for (column in names(df)) if (is.double(df[[column]])) df[[column]] <- num(df[[column]])
    utils::write.table(df, path, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
}

# ---- input ---------------------------------------------------------------------------------------------------------
readInput <- function(dir) {
    character_columns <- c(cell = "character", mask = "character", subset = "character", class = "character",
        stratum = "character")
    rows <- do.call(rbind, lapply(CELLS, \(cellname) utils::read.delim(file.path(dir, paste0(cellname, ".tsv.gz")),
        colClasses = character_columns, stringsAsFactors = FALSE)))
    metrics <- utils::read.delim(file.path(dir, "s2metrics.tsv"), colClasses = c(cell = "character",
        mask = "character", qmethod = "character"), stringsAsFactors = FALSE)
    list(rows = rows, metrics = metrics)
}

# The sums of one (cell, mask) per replicate: an R x K matrix (rows in the order of `replicates`).
atomMatrix <- function(d, m, replicates) {
    key <- factor(d$r, levels = replicates)
    total <- function(keep, column) {
        v <- tapply(d[[column]][keep], key[keep], sum)
        v[is.na(v)] <- 0
        as.numeric(v)
    }
    is_null <- d$class == "null"
    is_planted <- d$class %in% c("up", "down")
    one_sided <- d$stratum %in% c("one_treated", "one_control")
    atoms <- list(one = rep(1, length(replicates)))
    for (s in SUBSETS) {
        in_s <- d$subset == s
        atoms[[paste0("n0_", s)]] <- total(in_s & is_null, "n_tested")
        atoms[[paste0("p05_", s)]] <- total(in_s & is_null, "n_p05")
        atoms[[paste0("p01_", s)]] <- total(in_s & is_null, "n_p01")
        atoms[[paste0("V_", s)]] <- total(in_s & is_null, "n_call")
        atoms[[paste0("n1_", s)]] <- total(in_s & is_planted, "n_tested")
        atoms[[paste0("T_", s)]] <- total(in_s & is_planted, "n_call")
    }
    for (s in c("added", "kept_both")) {
        in_s <- d$subset == s
        atoms[[paste0("err_", s)]] <- total(in_s & is_planted, "sum_err")
        atoms[[paste0("sgn_", s)]] <- total(in_s & is_planted, "sum_sgn_err")
        for (direction in c("up", "down")) {
            atoms[[paste0("n_", direction, "_", s)]] <- total(in_s & d$class == direction, "n_tested")
            atoms[[paste0("err_", direction, "_", s)]] <- total(in_s & d$class == direction, "sum_err")
        }
        for (st in STRATA) for (l in 0:1) {
            keep <- in_s & is_null & d$stratum == st & d$ell == l
            tag <- paste(s, st, l, sep = "_")
            atoms[[paste0("sn0_", tag)]] <- total(keep, "n_tested")
            atoms[[paste0("sp05_", tag)]] <- total(keep, "n_p05")
            atoms[[paste0("srv_", tag)]] <- total(keep, "sum_resvar")
            atoms[[paste0("snrv_", tag)]] <- total(keep, "n_resvar")
        }
    }
    added_null <- d$subset == "added" & is_null
    for (l in 0:1) {
        keep <- added_null & one_sided & d$ell == l
        atoms[[paste0("m6sum_", l)]] <- total(keep, "sum_s_err")
        atoms[[paste0("m6n_", l)]] <- total(keep, "n_tested")
        atoms[[paste0("dcall_", l)]] <- total(keep, "n_call")
        atoms[[paste0("daway_", l)]] <- total(keep, "n_call_away")
        atoms[[paste0("n0_added_l", l)]] <- total(added_null & d$ell == l, "n_tested")
        atoms[[paste0("p05_added_l", l)]] <- total(added_null & d$ell == l, "n_p05")
    }
    keep <- d$subset == "added" & is_planted & one_sided
    atoms$serr_added <- total(keep, "sum_s_err")
    atoms$nserr_added <- total(keep, "n_tested")
    mm <- m[match(replicates, m$r), , drop = FALSE]
    for (s in SUBSETS) atoms[[paste0("n0t_", s)]] <- atoms[[paste0("n0_", s)]] * mm$t_r
    # Revision 1.0's counterfactual (reported only): Delta's tested nulls called at Ref's own rate of replicate r.
    atoms$n0c_added <- atoms$n0_added * ifelse(atoms$n0_ref > 0, atoms$V_ref / atoms$n0_ref, 0)
    atoms$afdp_sum <- ifelse(is.na(mm$added_fdp), 0, mm$added_fdp)
    atoms$afdp_n <- as.numeric(!is.na(mm$added_fdp))
    do.call(cbind, atoms)
}

# ---- estimands (5.3), on a matrix of summed atoms (one row per resample, jackknife unit or the estimate) --------------
fdrModel <- function(pi0, pow, t) pi0 * t / (pi0 * t + (1 - pi0) * pow)

computeStats <- function(S) {
    g <- function(name) S[, name]
    alpha <- function(s, p = "p05_") g(paste0(p, s)) / g(paste0("n0_", s))
    fdr <- function(s) g(paste0("V_", s)) / (g(paste0("V_", s)) + g(paste0("T_", s)))
    cal <- function(s) g(paste0("n0t_", s)) / (g(paste0("n0t_", s)) + g(paste0("T_", s)))
    pi0 <- function(s) g(paste0("n0_", s)) / (g(paste0("n0_", s)) + g(paste0("n1_", s)))
    pow <- function(s) g(paste0("T_", s)) / g(paste0("n1_", s))
    tbar <- function(s) g(paste0("n0t_", s)) / g(paste0("n0_", s))
    a_d <- alpha("added")
    a_o <- alpha("kept_both")
    a_r <- alpha("ref")
    fdr_d <- fdr("added")
    fdr_o <- fdr("kept_both")
    cal_d <- cal("added")
    cal_o <- cal("kept_both")
    p0d <- pi0("added")
    p0o <- pi0("kept_both")
    pwd <- pow("added")
    pwo <- pow("kept_both")
    td <- tbar("added")
    to <- tbar("kept_both")
    fdr_ref_cf <- g("n0c_added") / (g("n0c_added") + g("T_added"))
    out <- list(
        # R1
        rho_O = a_d / a_o, rho_ref = a_d / a_r, alpha_added = a_d, alpha_kept_both = a_o, alpha_ref = a_r,
        kappa_added = a_d / 0.05, kappa_kept_both = a_o / 0.05, kappa_ref = a_r / 0.05,
        rho01_O = alpha("added", "p01_") / alpha("kept_both", "p01_"),
        rho01_ref = alpha("added", "p01_") / alpha("ref", "p01_"), alpha01_added = alpha("added", "p01_"),
        alpha01_kept_both = alpha("kept_both", "p01_"), alpha01_ref = alpha("ref", "p01_"),
        c_added = g("V_added") / g("n0t_added"), c_kept_both = g("V_kept_both") / g("n0t_kept_both"),
        c_ref = g("V_ref") / g("n0t_ref"),
        # R2
        m6_na = g("m6sum_0") / g("m6n_0"), m6_ell = g("m6sum_1") / g("m6n_1"),
        # R3
        fdr_added = fdr_d, fdr_kept_both = fdr_o, fdr_ref = fdr("ref"), fdr_whole = fdr("whole"),
        fdrcal_added = cal_d, fdrcal_kept_both = cal_o, E_cal = (fdr_d - cal_d) - (fdr_o - cal_o),
        E_comp = cal_d - cal_o, E_total = fdr_d - fdr_o,
        comp_power_first = fdrModel(p0o, pwd, to) - fdrModel(p0o, pwo, to),
        comp_null_second = fdrModel(p0d, pwd, td) - fdrModel(p0o, pwd, to),
        comp_null_first = fdrModel(p0d, pwo, td) - fdrModel(p0o, pwo, to),
        comp_power_second = fdrModel(p0d, pwd, td) - fdrModel(p0d, pwo, td),
        E_null_rev10 = fdr_d - fdr_ref_cf, E_rest_rev10 = fdr_ref_cf - fdr("ref"), gap_whole = fdr_d - fdr("whole"),
        added_fdp_mean_of_ratios = g("afdp_sum") / g("afdp_n"), pi0_added = p0d, pi0_kept_both = p0o,
        pi0_diff = p0d - p0o, power_added = pwd, power_kept_both = pwo, tbar_added = td, tbar_kept_both = to,
        channel_share_b1 = (g("p05_added_l1") - a_o * g("n0_added_l1")) / (g("p05_added") - a_o * g("n0_added")),
        # R4
        M5_added = g("err_added") / g("n1_added"), M5_kept_both = g("err_kept_both") / g("n1_kept_both"),
        M5_diff = g("err_added") / g("n1_added") - g("err_kept_both") / g("n1_kept_both"),
        M5_up_added = g("err_up_added") / g("n_up_added"), M5_up_kept_both = g("err_up_kept_both") / g("n_up_kept_both"),
        M5_down_added = g("err_down_added") / g("n_down_added"),
        M5_down_kept_both = g("err_down_kept_both") / g("n_down_kept_both"),
        A_added = g("sgn_added") / g("n1_added"), A_kept_both = g("sgn_kept_both") / g("n1_kept_both"),
        M5_sided_added = g("serr_added") / g("nserr_added"),
        # R6
        D_all = (g("daway_0") + g("daway_1")) / (g("dcall_0") + g("dcall_1")), D_ell0 = g("daway_0") / g("dcall_0"),
        D_ell1 = g("daway_1") / g("dcall_1"),
        # counts per replicate
        n0_added_per_r = g("n0_added") / g("one"), n1_added_per_r = g("n1_added") / g("one"),
        V_added_per_r = g("V_added") / g("one"), T_added_per_r = g("T_added") / g("one"),
        n0_kept_both_per_r = g("n0_kept_both") / g("one"), n0_added_ell1_per_r = g("n0_added_l1") / g("one"))
    # R7: within Delta by stratum and ell (null proteins), against O's rate; residual variance of Delta's and O's nulls
    for (st in STRATA) {
        sum2 <- function(prefix, s) g(paste0(prefix, s, "_", st, "_0")) + g(paste0(prefix, s, "_", st, "_1"))
        for (l in 0:1) {
            tag <- paste0("added_", st, "_", l)
            out[[paste0("st_rho_O_", st, "_ell", l)]] <- (g(paste0("sp05_", tag)) / g(paste0("sn0_", tag))) / a_o
            out[[paste0("st_n0_per_r_", st, "_ell", l)]] <- g(paste0("sn0_", tag)) / g("one")
            out[[paste0("st_rv_added_", st, "_ell", l)]] <- g(paste0("srv_", tag)) / g(paste0("snrv_", tag))
        }
        out[[paste0("st_rho_O_", st)]] <- (sum2("sp05_", "added") / sum2("sn0_", "added")) / a_o
        out[[paste0("st_n0_per_r_", st)]] <- sum2("sn0_", "added") / g("one")
        out[[paste0("st_rv_added_", st)]] <- sum2("srv_", "added") / sum2("snrv_", "added")
        out[[paste0("st_rv_kept_both_", st)]] <- sum2("srv_", "kept_both") / sum2("snrv_", "kept_both")
        out[[paste0("st_n0_kept_both_per_r_", st)]] <- sum2("sn0_", "kept_both") / g("one")
    }
    do.call(cbind, out)
}

statGroup <- function(stat) {
    ifelse(grepl("^st_", stat), "R7", ifelse(grepl("^(M5_|A_)", stat), "R4", ifelse(grepl("^m6_", stat), "R2",
        ifelse(grepl("^D_", stat), "R6", ifelse(grepl("_per_r$", stat), "counts",
            ifelse(grepl("^(rho|alpha|kappa|c_)", stat), "R1", "R3"))))))
}

# ---- resampling ----------------------------------------------------------------------------------------------------
bootstrapWeights <- function(replicate_counts) {
    set.seed(SEED, kind = "Mersenne-Twister", normal.kind = "Inversion", sample.kind = "Rejection")
    lapply(stats::setNames(CELLS, CELLS), \(cellname) {
        R <- replicate_counts[[cellname]]
        idx <- matrix(sample.int(R, B_BOOT * R, replace = TRUE), nrow = B_BOOT, byrow = TRUE)
        t(apply(idx, 1L, tabulate, nbins = R))
    })
}

jackSe <- function(theta) {
    n <- length(theta)
    if (n < 2L || any(!is.finite(theta))) return(NA_real_)
    sqrt((n - 1) / n * sum((theta - mean(theta))^2))
}

quant <- function(boot, p) apply(boot, 2L, stats::quantile, probs = p, type = 7, na.rm = TRUE, names = FALSE)

# Estimate, bootstrap and jackknife draws of every statistic for sums given per cell: `parts` is a list by cell of
# R x K matrices (each the sum over the masks that the setting or pool takes from that cell).
resampleStats <- function(parts, weights) {
    cells <- names(parts)
    total <- Reduce(`+`, lapply(parts, colSums))
    estimate <- computeStats(matrix(total, nrow = 1L, dimnames = list(NULL, names(total))))
    boot <- computeStats(Reduce(`+`, lapply(cells, \(cellname) weights[[cellname]] %*% parts[[cellname]])))
    jack <- computeStats(do.call(rbind, lapply(cells, \(cellname) matrix(total, nrow(parts[[cellname]]),
        length(total), byrow = TRUE, dimnames = list(NULL, names(total))) - parts[[cellname]])))
    list(estimate = estimate, boot = boot, jack = jack, R = sum(vapply(parts, nrow, integer(1L))))
}

summariseStats <- function(x, scope, cellname, mask) {
    stats_names <- colnames(x$boot)
    data.frame(scope = scope, cell = cellname, mask = mask, stat = stats_names, group = statGroup(stats_names),
        R = x$R, estimate = as.numeric(x$estimate[1L, ]), boot_se = apply(x$boot, 2L, stats::sd, na.rm = TRUE),
        jack_se = apply(x$jack, 2L, jackSe), lo95 = quant(x$boot, 0.025), hi95 = quant(x$boot, 0.975),
        boot_nonfinite = colSums(!is.finite(x$boot)), stringsAsFactors = FALSE, row.names = NULL)
}

differenceStats <- function(a, b) {
    keep <- R5_STATS
    list(estimate = a$estimate[, keep, drop = FALSE] - b$estimate[, keep, drop = FALSE],
        boot = a$boot[, keep, drop = FALSE] - b$boot[, keep, drop = FALSE],
        jack = a$jack[, keep, drop = FALSE] - b$jack[, keep, drop = FALSE], R = a$R)
}

# ---- labels (5.6) --------------------------------------------------------------------------------------------------
componentState <- function(x, component, level_present) {
    stat <- COMPONENT_STAT[[component]]
    boot <- x$boot[, stat]
    estimate <- x$estimate[1L, stat]
    jack_se <- jackSe(x$jack[, stat])
    lower <- stats::quantile(boot, level_present, type = 7, na.rm = TRUE, names = FALSE)
    upper <- stats::quantile(boot, 1 - LEVEL_ABSENT, type = 7, na.rm = TRUE, names = FALSE)
    present <- isTRUE(lower > NULL_VALUE[[component]])
    absent <- isTRUE(upper <= MARGIN[[component]])
    lower_j <- estimate - stats::qnorm(1 - level_present) * jack_se
    upper_j <- estimate + stats::qnorm(1 - LEVEL_ABSENT) * jack_se
    present_j <- isTRUE(lower_j > NULL_VALUE[[component]])
    absent_j <- isTRUE(upper_j <= MARGIN[[component]])
    state <- function(p, a) if (a) "absent" else if (p) "present" else "not determined"
    list(estimate = estimate, lower = lower, upper = upper, jack_se = jack_se, lower_jack = lower_j,
        upper_jack = upper_j, state = state(present, absent), state_jack = state(present_j, absent_j),
        note = if (present && absent) "non-zero, within the margin" else "")
}

settingLabel <- function(b, p) {
    if (b == "present" && p == "absent") "bias" else if (p == "present" && b == "absent") "power" else
        if (b == "present" && p == "present") "both" else if (b == "absent" && p == "absent") "neither" else
            "not determined"
}

attribution <- function(state_vd, state_mcar) {
    if (state_vd != "present") return("")
    switch(state_mcar, absent = "value-dependent", present = "not specific to value dependence", "not determined")
}

labelRow <- function(x, mcar, scope, cellname, mask, level_present) {
    b <- componentState(x, "B", level_present)
    p <- componentState(x, "P", level_present)
    label <- settingLabel(b$state, p$state)
    label_jack <- settingLabel(b$state_jack, p$state_jack)
    mb <- if (is.null(mcar)) NULL else componentState(mcar, "B", level_present)
    mp <- if (is.null(mcar)) NULL else componentState(mcar, "P", level_present)
    share <- x$boot[, "channel_share_b1"]
    data.frame(scope = scope, cell = cellname, mask = mask, R = x$R,
        rho_O = b$estimate, rho_O_lower = b$lower, rho_O_upper95 = b$upper, rho_O_jack_se = b$jack_se,
        rho_O_lower_jack = b$lower_jack, rho_O_upper_jack = b$upper_jack, B_state = b$state,
        B_state_jack = b$state_jack, B_note = b$note,
        E_comp = p$estimate, E_comp_lower = p$lower, E_comp_upper95 = p$upper, E_comp_jack_se = p$jack_se,
        E_comp_lower_jack = p$lower_jack, E_comp_upper_jack = p$upper_jack, P_state = p$state,
        P_state_jack = p$state_jack, P_note = p$note, label = label, label_jack = label_jack,
        label_differs_jack = label != label_jack,
        mcar_B_state = if (is.null(mb)) "" else mb$state, mcar_P_state = if (is.null(mp)) "" else mp$state,
        attribution_B = if (is.null(mb)) "" else attribution(b$state, mb$state),
        attribution_P = if (is.null(mp)) "" else attribution(p$state, mp$state),
        channel_share_b1 = if (label %in% c("bias", "both")) x$estimate[1L, "channel_share_b1"] else NA_real_,
        channel_share_lo95 = if (label %in% c("bias", "both")) stats::quantile(share, 0.025, type = 7, na.rm = TRUE,
            names = FALSE) else NA_real_,
        channel_share_hi95 = if (label %in% c("bias", "both")) stats::quantile(share, 0.975, type = 7, na.rm = TRUE,
            names = FALSE) else NA_real_,
        stringsAsFactors = FALSE)
}

# ---- G1 ------------------------------------------------------------------------------------------------------------
gateG1 <- function(metrics, s2_cells, out) {
    s2 <- utils::read.delim(s2_cells, stringsAsFactors = FALSE)
    rows <- list()
    for (cellname in CELLS) for (mask in VD) for (i in seq_len(nrow(G1_MAP))) {
        m <- metrics[metrics$cell == cellname & metrics$mask == mask, , drop = FALSE]
        m <- m[order(m$r), , drop = FALSE]
        x <- as.numeric(m[[G1_MAP$metric[[i]]]])
        x <- x[!is.na(x)]
        ref <- s2[s2$cell == cellname & s2$mask == mask & s2$side == G1_MAP$side[[i]] &
            s2$metric == G1_MAP$s2_metric[[i]], , drop = FALSE]
        mine <- if (length(x)) mean(x) else NA_real_
        theirs <- if (nrow(ref) == 1L) ref$mean else NA_real_
        rows[[length(rows) + 1L]] <- data.frame(cell = cellname, mask = mask, metric = G1_MAP$s2_metric[[i]],
            side = G1_MAP$side[[i]], reproduced = mine, s2 = theirs, R_reproduced = length(x),
            R_s2 = if (nrow(ref) == 1L) ref$R else NA_integer_,
            pass = isTRUE(mine == theirs || abs(mine - theirs) <= G1_TOLERANCE * abs(theirs)) &&
                isTRUE(length(x) == ref$R[1L]), stringsAsFactors = FALSE)
    }
    g1 <- do.call(rbind, rows)
    g1$replicates_1_to_100 <- identical(sort(unique(metrics$r)), 1:100)
    writeTsv(g1, file.path(out, "g1.tsv"))
    if (!all(g1$pass) || !all(g1$replicates_1_to_100)) {
        writeTsv(g1[!g1$pass | !g1$replicates_1_to_100, , drop = FALSE], file.path(out, "g1-mismatch.tsv"))
        k197Abort("k197_g1", sprintf("gate G1 failed: %d of %d comparisons differ", sum(!g1$pass), nrow(g1)))
    }
    g1
}

# ---- the sentence (5.7) --------------------------------------------------------------------------------------------
fillSentence <- function(outcome, pool, labels, differing, R) {
    v <- function(stat) pool$estimate[pool$stat == stat]
    iv <- function(stat, digits) sprintf("[%s, %s]", formatC(pool$lo95[pool$stat == stat], digits, format = "f"),
        formatC(pool$hi95[pool$stat == stat], digits, format = "f"))
    f <- function(x, digits) formatC(x, digits, format = "f")
    pct <- function(x, digits) formatC(100 * x, digits, format = "f")
    pooled <- labels[labels$scope == "pool" & labels$mask == "vd", , drop = FALSE]
    mcar_b <- pooled$mcar_B_state
    mcar_p <- pooled$mcar_P_state
    prefix <- sprintf(paste("In simulation, the proteins the count rule newly keeps carry a higher share of false calls",
        "(about %s %% against %s %% among the proteins kept under both rules)"), pct(v("fdr_added"), 1),
        pct(v("fdr_kept_both"), 1))
    b_clause <- sprintf("unchanged proteins among them reach p below 0.05 in %s %% of cases against %s %%",
        pct(v("alpha_added"), 2), pct(v("alpha_kept_both"), 2))
    share_clause <- sprintf(paste("; %s %% of this excess is carried by unchanged proteins holding a value at or below",
        "the intensity threshold, which the filter counts as missing but the test keeps; the rest is not separated",
        "between selection on variance and unequal group sizes"), pct(pooled$channel_share_b1, 0))
    m6_clause <- sprintf("The missing value itself does not shift the fold change of an unchanged protein (%s log2 %s)",
        f(v("m6_na"), 4), iv("m6_na", 4))
    mcar_b_clause <- switch(mcar_b, present = "; the excess does appear when the same share of values is missing at random",
        absent = "; the excess does not appear when the same share of values is missing at random",
        "; whether the excess appears when the same share of values is missing at random was not determined")
    p_clause <- sprintf(paste("fewer truly changed proteins are found among them (power %s against %s; %s of the",
        "excess, %s in the other order) and %s of them are unchanged (%s against %s; %s)"),
        f(v("power_added"), 3), f(v("power_kept_both"), 3), f(v("comp_power_first"), 4), f(v("comp_power_second"), 4),
        if (pool$lo95[pool$stat == "pi0_diff"] > 0) "more" else "not more", f(v("pi0_added"), 3),
        f(v("pi0_kept_both"), 3), f(v("comp_null_second"), 4))
    mcar_p_clause <- switch(mcar_p, present = "; the excess also appears when the same share of values is missing at random",
        absent = "; the excess does not appear when the same share of values is missing at random",
        "; whether the excess appears when the same share of values is missing at random was not determined")
    r4_clause <- sprintf(paste(" The estimated fold change of truly changed proteins among them is shrunk toward 0",
        "by %s log2 on average (%s among the proteins kept under both rules)."), f(-v("A_added"), 3),
        f(-v("A_kept_both"), 3))
    differ_clause <- if (nrow(differing)) sprintf(" Deciding settings whose own label differs from the pooled one: %s.",
        paste(sprintf("%s (%s)", differing$setting, differing$label), collapse = ", ")) else ""
    closing <- sprintf(paste(" (simulated data, KIT-185's generator with no unwanted variation; groups of 3 or 4;",
        "%d new draws per setting; the size on real data is not measured)"), R)
    body <- switch(outcome,
        O1 = paste0(prefix, " because ", b_clause, share_clause, ". ", m6_clause, mcar_b_clause, "."),
        O2 = paste0(prefix, " because ", p_clause, "; ", b_clause, ", within the declared 10 % margin", mcar_p_clause, "."),
        O3 = paste0(prefix, ": ", f(v("E_cal"), 4), " of the excess because ", b_clause, share_clause, ", and ",
            f(v("E_comp"), 4), " because ", p_clause, ". ", m6_clause, mcar_b_clause, mcar_p_clause, "."),
        O4 = paste0(prefix, "; the cause was not determined at ", R, " new draws per setting (pooled estimates: ",
            "ratio of null p-values below 0.05 ", f(v("rho_O"), 3), " ", iv("rho_O", 3), ", false-share part from ",
            "power and null share ", f(v("E_comp"), 4), " ", iv("E_comp", 4), ")."),
        NULL)
    if (is.null(body)) return("No sentence: outcome O5 or Stop (no NEWS clause until the owner decides; 5.7).")
    paste0(body, r4_clause, differ_clause, closing)
}

# ---- main ----------------------------------------------------------------------------------------------------------
args <- commandArgs(TRUE)
mode <- if (length(args)) args[[1L]] else ""
if (!(mode == "descriptive" && length(args) == 4L) && !(mode == "deciding" && length(args) == 3L))
    k197Abort("k197_usage", "usage: score-k197.R descriptive <rows dir> <S2 cells.tsv> <out dir> | deciding <rows dir> <out dir>")
rows_dir <- args[[2L]]
out <- args[[length(args)]]
dir.create(out, recursive = TRUE, showWarnings = FALSE)
input <- readInput(rows_dir)
metrics <- input$metrics
replicates <- lapply(stats::setNames(CELLS, CELLS), \(cellname) sort(unique(metrics$r[metrics$cell == cellname])))
expected <- sum(vapply(replicates, length, integer(1L))) * 6L
checks <- data.frame(check = "complete", value = nrow(metrics), pass = nrow(metrics) == expected &&
    !anyDuplicated(metrics[c("cell", "mask", "r")]) && all(c(VD, MCAR) %in% metrics$mask),
    detail = sprintf("%d replicate-masks of %d expected", nrow(metrics), expected), stringsAsFactors = FALSE)
if (mode == "descriptive") gateG1(metrics, args[[3L]], out)
f1 <- sum(metrics$f1_o_not_h)
checks <- rbind(checks,
    data.frame(check = "F1", value = f1, pass = f1 == 0, detail = "proteins in O and not in H, summed", stringsAsFactors = FALSE),
    data.frame(check = "t_r_monotone", value = sum(!metrics$mono_ok), pass = all(metrics$mono_ok),
        detail = "replicate-masks where the count of tested p <= t_r differs from the calls", stringsAsFactors = FALSE),
    data.frame(check = "bh_fallback", value = sum(metrics$qmethod != "qvalue"), pass = TRUE,
        detail = "replicate-masks whose q-values fell back to BH (counted, not a stop)", stringsAsFactors = FALSE))
if (file.exists(file.path(rows_dir, "determinism.tsv"))) {
    det <- utils::read.delim(file.path(rows_dir, "determinism.tsv"), stringsAsFactors = FALSE)
    checks <- rbind(checks, data.frame(check = "determinism", value = nrow(det),
        pass = all(det$rows_identical) && all(det$metrics_identical), detail = "cell A, first five replicates run twice",
        stringsAsFactors = FALSE))
}
if (!checks$pass[checks$check == "complete"]) {
    writeTsv(checks, file.path(out, "checks.tsv"))
    k197Abort("k197_instrument", "the rows are incomplete")
}
if (f1 > 0) {
    writeTsv(checks, file.path(out, "checks.tsv"))
    writeTsv(data.frame(outcome = "Stop", reason = "F1: a protein in O and not in H", stringsAsFactors = FALSE),
        file.path(out, "outcome.tsv"))
    k197Abort("k197_falsifier", "F1: a protein kept by the old rule is not kept by the package's rule")
}

# Atoms per cell and mask, then the bootstrap weights.
atoms <- lapply(stats::setNames(CELLS, CELLS), \(cellname) lapply(stats::setNames(c(VD, MCAR), c(VD, MCAR)), \(mask)
    atomMatrix(input$rows[input$rows$cell == cellname & input$rows$mask == mask, , drop = FALSE],
        metrics[metrics$cell == cellname & metrics$mask == mask, , drop = FALSE], replicates[[cellname]])))
weights <- bootstrapWeights(lapply(replicates, length))

settings <- list()
for (cellname in CELLS) for (mask in c(VD, MCAR)) settings[[paste(cellname, mask, sep = ":")]] <-
    resampleStats(stats::setNames(list(atoms[[cellname]][[mask]]), cellname), weights)
poolParts <- function(masks_of) lapply(stats::setNames(names(POOL), names(POOL)), \(cellname)
    Reduce(`+`, lapply(masks_of(POOL[[cellname]]), \(mask) atoms[[cellname]][[mask]])))
pools <- list(vd = resampleStats(poolParts(identity), weights),
    mcar = resampleStats(poolParts(\(m) sub("^vd", "mcar", m)), weights))

table_settings <- do.call(rbind, lapply(names(settings), \(key) {
    parts <- strsplit(key, ":", fixed = TRUE)[[1L]]
    summariseStats(settings[[key]], "setting", parts[[1L]], parts[[2L]])
}))
table_pools <- rbind(summariseStats(pools$vd, "pool", "ten", "vd"), summariseStats(pools$mcar, "pool", "ten", "mcar"))
r5 <- list()
for (cellname in CELLS) for (mask in VD) r5[[length(r5) + 1L]] <- summariseStats(differenceStats(
    settings[[paste(cellname, mask, sep = ":")]], settings[[paste(cellname, sub("^vd", "mcar", mask), sep = ":")]]),
    "setting", cellname, paste0(mask, "-minus-mcar"))
r5[[length(r5) + 1L]] <- summariseStats(differenceStats(pools$vd, pools$mcar), "pool", "ten", "vd-minus-mcar")
table_r5 <- do.call(rbind, r5)
table_r5$group <- "R5"
all_stats <- rbind(table_settings, table_pools)
main_groups <- c("R1", "R2", "R3", "R6", "counts")
writeTsv(table_r5, file.path(out, "r5-mcar.tsv"))
writeTsv(all_stats[all_stats$group == "R7", , drop = FALSE], file.path(out, "strata.tsv"))
writeTsv(all_stats[all_stats$group == "R4", , drop = FALSE], file.path(out, "m5.tsv"))

provenance <- list(ticket = "KIT-197", mode = mode, script_sha256 = digest::digest(file = normalizePath(sub("^--file=", "",
    grep("^--file=", commandArgs(FALSE), value = TRUE)[[1L]])), algo = "sha256"),
    inputs_sha256 = as.list(stats::setNames(vapply(c(paste0(CELLS, ".tsv.gz"), "s2metrics.tsv"), \(f)
        digest::digest(file = file.path(rows_dir, f), algo = "sha256"), character(1L)),
        c(paste0(CELLS, ".tsv.gz"), "s2metrics.tsv"))),
    replicates = lapply(replicates, range), bootstrap = list(B = B_BOOT, seed = SEED, kind = "Mersenne-Twister",
        normal_kind = "Inversion", sample_kind = "Rejection", unit = "cell replicate, all six masks", quantile_type = 7L),
    margins = as.list(MARGIN), level_present = as.list(LEVEL_PRESENT), level_absent = LEVEL_ABSENT, level_f4 = LEVEL_F4,
    r_version = R.version.string)

if (mode == "descriptive") {
    writeTsv(all_stats[all_stats$group %in% main_groups, , drop = FALSE], file.path(out, "descriptive.tsv"))
    writeTsv(checks, file.path(out, "checks.tsv"))
    jsonlite::write_json(provenance, file.path(out, "provenance-score.json"), auto_unbox = TRUE, pretty = TRUE, digits = NA)
    cat("KIT-197 score (descriptive): G1 passed; tables written\n")
    quit(status = 0L)
}

# ---- deciding: labels, F4, outcome -----------------------------------------------------------------------------------
writeTsv(all_stats[all_stats$scope == "setting" & all_stats$group %in% main_groups, , drop = FALSE],
    file.path(out, "settings.tsv"))
writeTsv(table_pools[table_pools$group %in% main_groups, , drop = FALSE], file.path(out, "pooled.tsv"))
m6 <- table_settings[table_settings$stat == "m6_na", c("cell", "mask", "R", "estimate"), drop = FALSE]
m6$lower_bonferroni <- vapply(names(settings), \(key) stats::quantile(settings[[key]]$boot[, "m6_na"], LEVEL_F4 / 2,
    type = 7, na.rm = TRUE, names = FALSE), numeric(1L))
m6$upper_bonferroni <- vapply(names(settings), \(key) stats::quantile(settings[[key]]$boot[, "m6_na"], 1 - LEVEL_F4 / 2,
    type = 7, na.rm = TRUE, names = FALSE), numeric(1L))
m6$excludes_zero <- m6$lower_bonferroni > 0 | m6$upper_bonferroni < 0
f4 <- any(m6$excludes_zero)
writeTsv(m6, file.path(out, "f4-m6na.tsv"))
checks <- rbind(checks, data.frame(check = "F4", value = sum(m6$excludes_zero), pass = !f4,
    detail = "settings whose Bonferroni (24) interval of m6_NA excludes 0", stringsAsFactors = FALSE))

labels <- list()
for (cellname in CELLS) for (mask in c(VD, MCAR)) {
    mcar <- if (mask %in% VD) settings[[paste(cellname, sub("^vd", "mcar", mask), sep = ":")]] else NULL
    labels[[length(labels) + 1L]] <- labelRow(settings[[paste(cellname, mask, sep = ":")]], mcar, "setting", cellname,
        mask, LEVEL_PRESENT[["setting"]])
}
labels[[length(labels) + 1L]] <- labelRow(pools$vd, pools$mcar, "pool", "ten", "vd", LEVEL_PRESENT[["pool"]])
labels[[length(labels) + 1L]] <- labelRow(pools$mcar, NULL, "pool", "ten", "mcar", LEVEL_PRESENT[["pool"]])
labels <- do.call(rbind, labels)
labels$deciding <- labels$scope == "pool" | (labels$mask %in% VD & !(labels$cell == "D" & labels$mask %in% c("vd5", "vd10")))
if (f4) labels$label <- "not assigned (F4)"
writeTsv(labels, file.path(out, "labels.tsv"))

pooled_label <- labels$label[labels$scope == "pool" & labels$mask == "vd"]
outcome <- if (f4) "Stop" else switch(pooled_label, neither = "O5", `not determined` = "O4", bias = "O1", power = "O2",
    both = "O3")
deciding_rows <- labels[labels$scope == "setting" & labels$deciding, , drop = FALSE]
differing <- if (f4) deciding_rows[0L, ] else deciding_rows[deciding_rows$label != pooled_label, , drop = FALSE]
differing$setting <- if (nrow(differing)) paste(differing$cell, differing$mask) else character()
R_new <- min(vapply(replicates, length, integer(1L)))
sentence <- if (outcome == "Stop") "No sentence: Stop (F4; the owner is told)." else
    fillSentence(outcome, table_pools[table_pools$mask == "vd", , drop = FALSE], labels, differing, R_new)
writeTsv(data.frame(outcome = outcome, pooled_label = pooled_label, f1 = f1 > 0, f4 = f4,
    differing_settings = paste(differing$setting, collapse = "; "), stringsAsFactors = FALSE),
    file.path(out, "outcome.tsv"))
writeTsv(checks, file.path(out, "checks.tsv"))
writeLines(c("# KIT-197 deciding scoring (generated by score-k197.R)", "",
    sprintf("Outcome: **%s** (pooled label: %s).", outcome, pooled_label), "", "Sentence (5.7; the owner edits, VAL-41):", "",
    paste0("> ", sentence)), file.path(out, "sentence.md"))
jsonlite::write_json(provenance, file.path(out, "provenance-score.json"), auto_unbox = TRUE, pretty = TRUE, digits = NA)
cat(sprintf("KIT-197 score (deciding): outcome %s\n", outcome))
