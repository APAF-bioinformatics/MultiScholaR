# KIT-217 measurement M, the decision script (ticket KIT-217 r1.3 section 4.1; methods ledger VAL-184; the declaration
# docs/qa/workflow-execution/kit/kit-217/m/declaration.md). Reads the counted rows of sim-k217.R and applies every
# reading rule fixed before the counted runs. Evidence code under docs/qa: no package code; needs base R only.
#
# Usage (from the package root):
#   Rscript docs/qa/workflow-execution/kit/kit-217/m/score-k217.R score <rows dir> <out dir>
#   Rscript docs/qa/workflow-execution/kit/kit-217/m/score-k217.R selftest <out dir>
# score: <rows dir> holds m-rows-<cell>-<first>-<last>.tsv.gz of every cell and det-rows-L2-1-5.tsv.gz (determinism).
# selftest (run in the lifecycle): the same functions on planted rows (a calibrated arm, a 2-fold arm, a 2-fold positive
#   control); writes selftest.tsv and exits with status 1 when a planted label is not recovered.

args <- commandArgs(TRUE)
mode <- args[[1L]]
stopifnot(mode %in% c("score", "selftest"))
out_dir <- args[[length(args)]]
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ---- fixed constants (declaration section 5) -------------------------------------------------------------------------
B <- 100000L
BOOT_SEED <- c(L2 = 217990001, L1 = 217990002, M2 = 217990003, M1 = 217990004, Lv = 217990005, Le = 217990006, Me = 217990007)
DECIDING <- c("L2", "L1", "M2", "M1")
EXCH <- c("vd5", "vd10", "vd25", "mcar10")
GLINK <- c("vdg10", "vdg25")
FAMILY <- c("D", "Dk", "Atau_g2", "Aobs_any", "Aobs_g2")
F2_ARMS <- c("B", "Dk", "Aobs_any", "Aobs_g2", "Aobs_all")
BAND <- 0.075
EXCESS_EDGE <- 0.025
LEVEL_MEETS <- 0.025 / 5
LEVEL_FAILS <- 0.025 / (5 * 4)
POWER_EDGE <- 0.005
FDR_EDGE <- 0.0025
DNC_SHARE <- 0.05

# ---- resampling: one index set per cell, shared by every arm, mask and statistic --------------------------------------
# Row b of the cell's matrix holds how often each replicate is drawn in resample b:
# set.seed(seed); idx = matrix(sample.int(n, n * B, replace = TRUE), nrow = B).
bootWeights <- function(n, seed) {
    set.seed(seed, kind = "Mersenne-Twister", normal.kind = "Inversion", sample.kind = "Rejection")
    idx <- sample.int(n, n * B, replace = TRUE)
    counts <- tabulate((idx - 1L) * B + rep(seq_len(B), times = n), nbins = B * n)
    matrix(as.numeric(counts), B, n)
}
# A statistic: fn applied to column sums of per-replicate series (one n_c x p matrix per cell). Estimate, the B bootstrap
# values (stratified by cell) and the delete-one-replicate jackknife (stratified: var = sum_h (n_h - 1)/n_h sum_j
# (theta_hj - mean_h)^2, df = sum_h (n_h - 1)).
evalStat <- function(series, fn, W, point_only = FALSE) {
    cells <- names(series)
    total <- Reduce(`+`, lapply(series, colSums))
    est <- fn(matrix(total, nrow = 1L))
    if (point_only) return(list(est = est))
    boot <- fn(Reduce(`+`, lapply(cells, function(c) W[[c]] %*% series[[c]])))
    v <- 0
    df <- 0
    for (c in cells) {
        s <- series[[c]]
        nh <- nrow(s)
        th <- fn(matrix(total, nh, ncol(s), byrow = TRUE) - s)
        v <- v + (nh - 1) / nh * sum((th - mean(th))^2)
        df <- df + nh - 1
    }
    list(est = est, boot = boot, se = sqrt(v), df = df, boot_nan = sum(!is.finite(boot)))
}
upperLim <- function(st, level) c(boot = stats::quantile(st$boot, 1 - level, type = 7, names = FALSE, na.rm = TRUE),
    jack = st$est + stats::qt(1 - level, st$df) * st$se)
lowerLim <- function(st, level) c(boot = stats::quantile(st$boot, level, type = 7, names = FALSE, na.rm = TRUE),
    jack = st$est - stats::qt(1 - level, st$df) * st$se)
# the agreement rule: a classification holds only when bootstrap and jackknife give the same one
agree <- function(labels) if (length(unique(labels)) == 1L) labels[[1L]] else "not determined"
classifyBand <- function(st, edge, good, bad) {
    up <- upperLim(st, LEVEL_MEETS)
    lo <- lowerLim(st, LEVEL_FAILS)
    agree(vapply(c("boot", "jack"), function(m) if (up[[m]] < edge) good else if (lo[[m]] >= edge) bad else "not determined",
        character(1L)))
}
ratio <- function(S) S[, 1L] / S[, 2L]
diffRatio <- function(S) S[, 1L] / S[, 2L] - S[, 3L] / S[, 4L]
fmt <- function(x, d = 4) formatC(x, digits = d, format = "f")

# ---- the rows as arrays [replicate, mask, arm] ------------------------------------------------------------------------
NUMERIC <- c("kept", "kept_novalue", "kept_zerodf", "kept_nogroup", "miss_before", "miss_after", "n_tested", "n_null_tested",
    "n_null_rej", "n_null_low", "n_null_low_rej", "n_null_nolow", "n_null_nolow_rej", "n_planted", "n_planted_tested",
    "n_bh_calls", "n_bh_true", "n_bh_false", "fdp_bh", "n_q_calls", "n_q_true", "n_q_false", "fdp_q", "n_changed_vs_S",
    "sum_fc_null", "n_fc_null", "n_up", "n_down", "rm_up", "rm_down", "n_onoff", "rm_onoff", "nK", "nK_rej_X",
    "nK_NF_tested", "nK_rej_NF")
cellData <- function(rows) {
    reps <- sort(unique(as.integer(rows$rep)))
    masks <- unique(rows$mask)
    arms <- unique(rows$arm)
    idx <- cbind(match(as.integer(rows$rep), reps), match(rows$mask, masks), match(rows$arm, arms))
    ok <- array(FALSE, c(length(reps), length(masks), length(arms)), list(reps, masks, arms))
    ok[idx] <- rows$status == "ok"
    values <- lapply(stats::setNames(intersect(NUMERIC, names(rows)), intersect(NUMERIC, names(rows))), function(col) {
        a <- array(NA_real_, dim(ok), dimnames(ok))
        a[idx] <- suppressWarnings(as.numeric(rows[[col]]))
        a
    })
    list(reps = reps, masks = masks, arms = arms, ok = ok, v = values, rows = rows)
}
# per-replicate series of one column for one arm over masks; `with`: arms that must also have completed the (rep, mask)
seriesOf <- function(cd, arm, masks, col, with = character()) {
    masks <- intersect(masks, cd$masks)
    ok <- cd$ok[, masks, arm, drop = FALSE]
    for (other in with) ok <- ok & cd$ok[, masks, other, drop = FALSE]
    x <- if (identical(col, "one")) array(1, dim(ok)) else cd$v[[col]][, masks, arm, drop = FALSE]
    x[!ok | is.na(x)] <- 0
    rowSums(matrix(x, nrow = dim(x)[1L]))
}
# the same statistic built per cell: cols = list(c(arm, col), ...), every column paired on the arms named
statFor <- function(CD, W, cells, masks, cols, fn, pair = character(), point_only = FALSE) {
    series <- lapply(stats::setNames(cells, cells), function(c) {
        cd <- CD[[c]]
        do.call(cbind, lapply(cols, function(ac) seriesOf(cd, ac[[1L]], masks, ac[[2L]], setdiff(pair, ac[[1L]]))))
    })
    evalStat(series, fn, W[cells], point_only)
}
alphaStat <- function(CD, W, cell, arm, masks) statFor(CD, W, cell, masks, list(c(arm, "n_null_rej"), c(arm, "n_null_tested")), ratio)
excessStat <- function(CD, W, cell, arm, ref, masks) statFor(CD, W, cell, masks, list(c(arm, "n_null_rej"), c(arm, "n_null_tested"),
    c(ref, "n_null_rej"), c(ref, "n_null_tested")), diffRatio, pair = c(arm, ref))
dncShare <- function(cd, arm) {
    if (!arm %in% cd$arms) return(NA_real_)
    1 - mean(cd$ok[, , arm])
}

# ---- Q1, the guard, the class (declaration sections 6 and 7) ----------------------------------------------------------
readQ1 <- function(CD, W, arms_all) {
    out <- list()
    nf_meets <- list()
    for (c in DECIDING) {
        st <- alphaStat(CD, W, c, "NF", EXCH)
        nf_meets[[c]] <- all(upperLim(st, 0.025) < BAND)
        for (arm in intersect(arms_all, CD[[c]]$arms)) {
            a <- alphaStat(CD, W, c, arm, EXCH)
            no_result <- isTRUE(dncShare(CD[[c]], arm) > DNC_SHARE)
            if (identical(arm, "NF")) {
                label <- classifyBand(a, BAND, "meets", "fails")
                basis <- "absolute"
            } else if (nf_meets[[c]]) {
                label <- classifyBand(a, BAND, "meets", "fails")
                basis <- "absolute"
            } else {
                label <- classifyBand(excessStat(CD, W, c, arm, "NF", EXCH), EXCESS_EDGE, "meets", "fails")
                basis <- "excess over NF (NF does not meet)"
            }
            if (no_result) label <- "no result (did not complete above 5 %)"
            up <- upperLim(a, LEVEL_MEETS)
            lo <- lowerLim(a, LEVEL_FAILS)
            out[[length(out) + 1L]] <- data.frame(cell = c, arm = arm, alpha = a$est, upper_boot = up[["boot"]], upper_jack = up[["jack"]],
                lower_boot = lo[["boot"]], lower_jack = lo[["jack"]], dnc_share = dncShare(CD[[c]], arm), nf_meets = nf_meets[[c]],
                basis = basis, label = label, stringsAsFactors = FALSE)
        }
    }
    do.call(rbind, out)
}
verdictOver <- function(labels, cells) {
    if (!length(labels)) return("not determined (no reading cell)")
    if (all(labels == "meets")) return("calibrated")
    if (any(labels == "fails")) return(paste0("miscalibrated (", paste(cells[labels == "fails"], collapse = ", "), ")"))
    paste0("not determined (", paste(cells[labels != "meets"], collapse = ", "), ")")
}
guardFor <- function(CD, W, cells, arm) {
    pw <- statFor(CD, W, cells, EXCH, list(c("S", "n_bh_true"), c("S", "n_planted"), c(arm, "n_bh_true"), c(arm, "n_planted")),
        diffRatio, pair = c("S", arm))
    fd <- statFor(CD, W, cells, EXCH, list(c(arm, "fdp_bh"), c(arm, "one"), c("S", "fdp_bh"), c("S", "one")), diffRatio,
        pair = c("S", arm))
    pl <- c(lo = lowerLim(pw, LEVEL_MEETS), up = upperLim(pw, LEVEL_MEETS))
    fl <- c(lo = lowerLim(fd, LEVEL_MEETS), up = upperLim(fd, LEVEL_MEETS))
    power <- agree(vapply(c("boot", "jack"), function(m) if (pl[[paste0("lo.", m)]] > POWER_EDGE) "costs power" else
        if (pl[[paste0("up.", m)]] <= POWER_EDGE) "no material cost" else "not determined", character(1L)))
    fdr <- agree(vapply(c("boot", "jack"), function(m) if (fl[[paste0("lo.", m)]] > FDR_EDGE) "raised" else
        if (fl[[paste0("up.", m)]] <= FDR_EDGE) "not raised" else "not determined", character(1L)))
    data.frame(arm = arm, power_loss_vs_S = pw$est, power_lo_boot = pl[["lo.boot"]], power_up_boot = pl[["up.boot"]],
        power_lo_jack = pl[["lo.jack"]], power_up_jack = pl[["up.jack"]], power = power, fdr_rise_vs_S = fd$est,
        fdr_lo_boot = fl[["lo.boot"]], fdr_up_boot = fl[["up.boot"]], fdr_lo_jack = fl[["lo.jack"]], fdr_up_jack = fl[["up.jack"]],
        fdr = fdr, stringsAsFactors = FALSE)
}
classOf <- function(verdict, power, fdr) {
    if (identical(verdict, "calibrated")) return(if (identical(power, "no material cost") && identical(fdr, "not raised")) "E" else "C")
    if (startsWith(verdict, "miscalibrated")) return("R")
    "U"
}

# ---- scoring ------------------------------------------------------------------------------------------------------------
readRows <- function(file) utils::read.delim(gzfile(file), stringsAsFactors = FALSE, colClasses = "character", na.strings = "NA")

score <- function(rows_dir) {
    files <- list.files(rows_dir, pattern = "^m-rows-.*\\.tsv\\.gz$", full.names = TRUE)
    all_rows <- do.call(rbind, lapply(files, readRows))
    cells <- intersect(c(DECIDING, "Lv", "Le", "Me"), unique(all_rows$cell))
    stopifnot(all(DECIDING %in% cells))
    mask_rows <- all_rows[all_rows$arm == "mask", , drop = FALSE]
    arm_rows <- all_rows[all_rows$arm != "mask", , drop = FALSE]
    CD <- lapply(stats::setNames(cells, cells), function(c) cellData(arm_rows[arm_rows$cell == c, , drop = FALSE]))
    W <- lapply(stats::setNames(cells, cells), function(c) bootWeights(length(CD[[c]]$reps), BOOT_SEED[[c]]))
    notes <- character()
    # ---- stop rules: E1, F1 (counted rows), determinism ----
    e1_fail <- sum(mask_rows$e1 %in% "FALSE")
    f1_fail <- sum(mask_rows$f1 %in% "FALSE")
    det_file <- list.files(rows_dir, pattern = "^det-rows-L2-1-5\\.tsv\\.gz$", full.names = TRUE)
    det_ok <- FALSE
    if (length(det_file) == 1L) {
        det <- readRows(det_file)
        first <- all_rows[all_rows$cell == "L2" & as.integer(all_rows$rep) %in% 1:5, , drop = FALSE]
        key <- function(x) order(as.integer(x$rep), x$mask, x$arm)
        det_ok <- identical(unname(as.list(det[key(det), ])), unname(as.list(first[key(first), ])))
    }
    stops <- data.frame(rule = c("E1 (oracle equals the package's S)", "F1 (subset identities)", "determinism (L2 1-5)"),
        failures = c(e1_fail, f1_fail, if (det_ok) 0L else 1L), stringsAsFactors = FALSE)
    utils::write.table(stops, file.path(out_dir, "stops.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
    writeOutcome <- function(outcome, recommend, extra = character()) {
        table <- data.frame(item = c("outcome", "recommendation", extra[c(TRUE, FALSE)]), value = c(outcome, recommend,
            extra[c(FALSE, TRUE)]), stringsAsFactors = FALSE)
        utils::write.table(table, file.path(out_dir, "outcome.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
        print(table)
        invisible(table)
    }
    if (any(stops$failures > 0L)) return(writeOutcome("Stop", "none: the owner is told; nothing else is read"))
    # ---- F2 per arm and deciding setting (two-sided 99 % Bonferroni over 80 intervals; floor 0.02 log2) ----
    f2 <- list()
    for (c in DECIDING) for (mask in EXCH) for (arm in intersect(F2_ARMS, CD[[c]]$arms)) {
        st <- statFor(CD, W, c, mask, list(c(arm, "sum_fc_null"), c(arm, "n_fc_null")), ratio)
        lv <- 0.01 / 80 / 2
        lo <- lowerLim(st, lv)
        up <- upperLim(st, lv)
        excludes <- (lo > 0) | (up < 0)
        f2[[length(f2) + 1L]] <- data.frame(cell = c, mask = mask, arm = arm, mean_fc = st$est, lo_boot = lo[["boot"]], up_boot = up[["boot"]],
            lo_jack = lo[["jack"]], up_jack = up[["jack"]], fires = any(excludes) && abs(st$est) > 0.02, stringsAsFactors = FALSE)
    }
    f2 <- do.call(rbind, f2)
    utils::write.table(f2, file.path(out_dir, "f2.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
    f2_arms <- unique(f2$arm[f2$fires])
    # ---- Q1 ----
    q1 <- readQ1(CD, W, c("NF", "S", "PC", FAMILY))
    utils::write.table(q1, file.path(out_dir, "q1-cells.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
    # positive-control gate
    pc_pass <- vapply(DECIDING, function(c) {
        row <- q1[q1$cell == c & q1$arm == "NF", ]
        if (isTRUE(row$nf_meets)) all(lowerLim(alphaStat(CD, W, c, "PC", EXCH), 0.025) > BAND) else
            all(lowerLim(excessStat(CD, W, c, "PC", "NF", EXCH), 0.025) > EXCESS_EDGE)
    }, logical(1L))
    if (!all(pc_pass)) return(writeOutcome("N0", "none: the positive control does not register a 2-fold excess; no candidate is labelled",
        c("positive-control gate", paste(DECIDING, ifelse(pc_pass, "pass", "fail"), collapse = "; "))))
    # O0: the reading cells are those where S meets
    s_labels <- stats::setNames(q1$label[q1$arm == "S"], q1$cell[q1$arm == "S"])
    reading <- DECIDING[s_labels[DECIDING] == "meets"]
    if (!length(reading)) return(writeOutcome("O0 (no reading cell)", "none: S meets the band in no deciding cell; the owner decides on the per-cell figures (q1-cells.tsv)",
        c("positive-control gate", paste(DECIDING, ifelse(pc_pass, "pass", "fail"), collapse = "; "))))
    # verdicts, guard, class
    classes <- do.call(rbind, lapply(FAMILY, function(arm) {
        lab <- stats::setNames(q1$label[q1$arm == arm], q1$cell[q1$arm == arm])[reading]
        verdict <- verdictOver(unname(lab), reading)
        g <- if (length(reading)) guardFor(CD, W, reading, arm) else data.frame(arm = arm, power = "not determined", fdr = "not determined")
        cls <- if (arm %in% f2_arms) "no label (F2)" else classOf(verdict, g$power, g$fdr)
        cbind(data.frame(arm = arm, verdict = verdict, class = cls, stringsAsFactors = FALSE), g[, setdiff(names(g), "arm"), drop = FALSE])
    }))
    utils::write.table(classes, file.path(out_dir, "classes.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
    # ---- Q2 (reported, SD-10 (ii-b)) ----
    q2 <- do.call(rbind, lapply(DECIDING, function(c) do.call(rbind, lapply(c("S", "B", FAMILY), function(arm) {
        e <- excessStat(CD, W, c, arm, "NF", GLINK)
        a <- alphaStat(CD, W, c, arm, GLINK)
        label <- classifyBand(e, EXCESS_EDGE, "robust", "not robust")
        data.frame(cell = c, arm = arm, alphaG = a$est, alphaG_NF = alphaStat(CD, W, c, "NF", GLINK)$est, excess = e$est,
            label = label, stringsAsFactors = FALSE)
    }))))
    utils::write.table(q2, file.path(out_dir, "q2.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
    q2v <- vapply(c("S", "B", FAMILY), function(arm) {
        l <- q2$label[q2$arm == arm]
        if (all(l == "robust")) "robust" else if (any(l == "not robust")) "not robust" else "not determined"
    }, character(1L))
    ad <- do.call(rbind, lapply(DECIDING, function(c) {
        st <- statFor(CD, W, c, GLINK, list(c("Aobs_g2", "n_null_rej"), c("Aobs_g2", "n_null_tested"), c("Dk", "n_null_rej"),
            c("Dk", "n_null_tested")), diffRatio, pair = c("Aobs_g2", "Dk"))
        data.frame(cell = c, eG_Aobs_g2_minus_Dk = st$est, lo95_boot = lowerLim(st, 0.025)[["boot"]], up95_boot = upperLim(st, 0.025)[["boot"]],
            lo95_jack = lowerLim(st, 0.025)[["jack"]], up95_jack = upperLim(st, 0.025)[["jack"]], stringsAsFactors = FALSE)
    }))
    utils::write.table(ad, file.path(out_dir, "q2-aobs-dk.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
    q2_flag <- any(q2$label[q2$arm == "Aobs_g2"] == "not robust") && all(q2$label[q2$arm == "Dk"] == "robust")
    # ---- outcome table (the first row that holds applies) ----
    cls <- stats::setNames(classes$class, classes$arm)
    {
        outcome <- if (length(reading) < length(DECIDING)) "O0" else ""
        if (outcome == "O0") notes <- c(notes, sprintf("S does not meet the band in: %s; O1 to O3 read on: %s",
            paste(setdiff(DECIDING, reading), collapse = ", "), if (length(reading)) paste(reading, collapse = ", ") else "no cell"))
        e_arms <- intersect(c("Aobs_any", "Dk"), names(cls)[cls == "E"])
        if (identical(unname(cls[["Aobs_g2"]]), "E")) {
            code <- "O1"
            recommend <- "Aobs_g2"
        } else if (length(e_arms)) {
            code <- "O2"
            recommend <- if (length(e_arms) == 1L) e_arms else {
                st <- statFor(CD, W, reading, EXCH, list(c("Aobs_any", "n_bh_true"), c("Aobs_any", "n_planted"), c("Dk", "n_bh_true"),
                    c("Dk", "n_planted")), diffRatio, pair = c("Aobs_any", "Dk"))
                lo <- lowerLim(st, 0.025)
                up <- upperLim(st, 0.025)
                if (all(lo > 0)) "Aobs_any" else if (all(up < 0)) "Dk" else {
                    changed <- vapply(e_arms, function(a) {
                        s <- statFor(CD, W, reading, EXCH, list(c(a, "n_changed_vs_S"), c(a, "one")), ratio)
                        s$est
                    }, numeric(1L))
                    if (changed[["Aobs_any"]] < changed[["Dk"]]) "Aobs_any" else if (changed[["Dk"]] < changed[["Aobs_any"]]) "Dk" else "Dk"
                }
            }
            recommend <- paste0(recommend, " (the owner chooses)")
        } else {
            code <- "O3"
            recommend <- "B"
        }
        outcome <- paste(c(if (outcome == "O0") "O0", code), collapse = " then ")
    }
    if (length(f2_arms)) notes <- c(notes, paste("F2 fired; no label for:", paste(f2_arms, collapse = ", ")))
    if (q2_flag) notes <- c(notes, "SD-10 (ii-b): A-obs-g2 is 'not robust' in Q2 in some cell and D-k is 'robust' in every cell: the owner is told before choosing; no class changes")
    # ---- large effects (Le, Me): the recommended candidate against D-k and B, per mask ----
    rec_arm <- base::sub(" .*$", "", recommend)
    large <- list()
    for (c in intersect(c("Le", "Me"), cells)) for (mask in CD[[c]]$masks) for (arm in setdiff(CD[[c]]$arms, "NF")) {
        rm <- statFor(CD, W, c, mask, list(c(arm, "rm_up"), c(arm, "n_planted")), ratio, point_only = TRUE)
        rmd <- statFor(CD, W, c, mask, list(c(arm, "rm_down"), c(arm, "n_planted")), ratio, point_only = TRUE)
        pw <- statFor(CD, W, c, mask, list(c(arm, "n_bh_true"), c(arm, "n_planted")), ratio, point_only = TRUE)
        oo <- statFor(CD, W, c, mask, list(c(arm, "rm_onoff"), c(arm, "n_onoff")), ratio, point_only = TRUE)
        large[[length(large) + 1L]] <- data.frame(cell = c, mask = mask, arm = arm, power_bh = pw$est, removed_up = rm$est,
            removed_down = rmd$est, removed_onoff_share = oo$est, stringsAsFactors = FALSE)
        if (identical(arm, rec_arm)) for (ref in c("Dk", "B")) if (!identical(ref, arm)) {
            both <- statFor(CD, W, c, mask, list(c(arm, "rm_up"), c(arm, "rm_down"), c(arm, "n_planted"), c(ref, "rm_up"),
                c(ref, "rm_down"), c(ref, "n_planted")), function(S) (S[, 1L] + S[, 2L]) / S[, 3L] - (S[, 4L] + S[, 5L]) / S[, 6L],
                pair = c(arm, ref))
            lo <- lowerLim(both, 0.025)
            if (all(lo > 0.05)) notes <- c(notes, sprintf("large effects: %s removes a larger share of planted features than %s in %s %s (difference %.3f, 95 %% lower limits %.3f / %.3f)",
                arm, ref, c, mask, both$est, lo[["boot"]], lo[["jack"]]))
        }
    }
    if (length(large)) utils::write.table(do.call(rbind, large), file.path(out_dir, "large-effects.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
    # ---- printed beside, never deciding ----
    beside <- list()
    for (c in cells) {
        masks_c <- if (c %in% DECIDING) EXCH else CD[[c]]$masks
        for (arm in CD[[c]]$arms) {
            a <- alphaStat(CD, W, c, arm, masks_c)
            vsS <- if (arm %in% c("S", "NF") || !"S" %in% CD[[c]]$arms) NULL else excessStat(CD, W, c, arm, "S", masks_c)
            pw <- statFor(CD, W, c, masks_c, list(c(arm, "n_bh_true"), c(arm, "n_planted")), ratio, point_only = TRUE)
            fd <- statFor(CD, W, c, masks_c, list(c(arm, "fdp_bh"), c(arm, "one")), ratio, point_only = TRUE)
            pq <- statFor(CD, W, c, masks_c, list(c(arm, "n_q_true"), c(arm, "n_planted")), ratio, point_only = TRUE)
            fq <- statFor(CD, W, c, masks_c, list(c(arm, "fdp_q"), c(arm, "one")), ratio, point_only = TRUE)
            low <- statFor(CD, W, c, masks_c, list(c(arm, "n_null_low_rej"), c(arm, "n_null_low")), ratio, point_only = TRUE)
            nolow <- statFor(CD, W, c, masks_c, list(c(arm, "n_null_nolow_rej"), c(arm, "n_null_nolow")), ratio, point_only = TRUE)
            kx <- statFor(CD, W, c, masks_c, list(c(arm, "nK_rej_X"), c(arm, "nK"), c(arm, "nK_rej_NF"), c(arm, "nK_NF_tested")),
                function(S) (S[, 1L] / S[, 2L]) / (S[, 3L] / S[, 4L]), point_only = TRUE)
            mean_of <- function(col) statFor(CD, W, c, masks_c, list(c(arm, col), c(arm, "one")), ratio, point_only = TRUE)$est
            beside[[length(beside) + 1L]] <- data.frame(cell = c, masks = paste(masks_c, collapse = ","), arm = arm, alpha = a$est,
                alpha_lo95 = lowerLim(a, 0.025)[["boot"]], alpha_up95 = upperLim(a, 0.025)[["boot"]],
                alpha_minus_S = if (is.null(vsS)) NA else vsS$est,
                alpha_minus_S_lo95 = if (is.null(vsS)) NA else lowerLim(vsS, 0.025)[["boot"]],
                alpha_minus_S_up95 = if (is.null(vsS)) NA else upperLim(vsS, 0.025)[["boot"]],
                power_bh = pw$est, fdr_bh = fd$est, power_q = pq$est, fdr_q = fq$est, alpha_low = low$est, alpha_nolow = nolow$est,
                ratio_low = low$est / nolow$est, fixed_set_ratio = kx$est, kept = mean_of("kept"), kept_novalue = mean_of("kept_novalue"),
                kept_zerodf = mean_of("kept_zerodf"), kept_nogroup = mean_of("kept_nogroup"), miss_before = mean_of("miss_before"),
                miss_after = mean_of("miss_after"), changed_calls_vs_S = mean_of("n_changed_vs_S"),
                mean_fc_null = statFor(CD, W, c, masks_c, list(c(arm, "sum_fc_null"), c(arm, "n_fc_null")), ratio, point_only = TRUE)$est,
                dnc_share = dncShare(CD[[c]], arm), stringsAsFactors = FALSE)
        }
    }
    beside <- do.call(rbind, beside)
    utils::write.table(beside, file.path(out_dir, "beside.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
    # per-setting null rate (16 per candidate), printed
    setting <- do.call(rbind, lapply(cells, function(c) do.call(rbind, lapply(CD[[c]]$masks, function(mask) do.call(rbind,
        lapply(CD[[c]]$arms, function(arm) {
            a <- alphaStat(CD, W, c, arm, mask)
            data.frame(cell = c, mask = mask, arm = arm, alpha = a$est, lo95 = lowerLim(a, 0.025)[["boot"]],
                up95 = upperLim(a, 0.025)[["boot"]], stringsAsFactors = FALSE)
        }))))))
    utils::write.table(setting, file.path(out_dir, "per-setting.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
    writeOutcome(outcome, recommend, c("positive-control gate", paste(DECIDING, ifelse(pc_pass, "pass", "fail"), collapse = "; "),
        "reading cells", paste(reading, collapse = ", "), "Q2 verdicts", paste(names(q2v), q2v, sep = ": ", collapse = "; "),
        "notes", paste(notes, collapse = " | ")))
}

# ---- selftest: planted rows ----------------------------------------------------------------------------------------------
selftest <- function() {
    set.seed(217990100L, kind = "Mersenne-Twister", normal.kind = "Inversion", sample.kind = "Rejection")
    rate <- c(NF = 0.05, S = 0.05, PC = 0.10, B = 0.05, D = 0.05, Dk = 0.10, Atau_any = 0.05, Atau_g2 = 0.05, Aobs_any = 0.05,
        Aobs_g2 = 0.05, Atau_all = 0.05, Aobs_all = 0.05)
    rows <- do.call(rbind, lapply(DECIDING, function(c) {
        lipid <- startsWith(c, "L")
        grid <- expand.grid(arm = if (lipid) names(rate)[1:10] else names(rate), mask = c("none", EXCH, GLINK), rep = 1:200,
            stringsAsFactors = FALSE)
        n_null <- if (lipid) 360L else 100L
        p <- rate[grid$arm] * exp(stats::rnorm(nrow(grid), 0, 0.15))
        data.frame(cell = c, rep = as.character(grid$rep), mask = grid$mask, arm = grid$arm, status = "ok",
            n_null_tested = as.character(n_null), n_null_rej = as.character(stats::rbinom(nrow(grid), n_null, p)), stringsAsFactors = FALSE)
    }))
    CD <- lapply(stats::setNames(DECIDING, DECIDING), function(c) cellData(rows[rows$cell == c, , drop = FALSE]))
    W <- lapply(stats::setNames(DECIDING, DECIDING), function(c) bootWeights(length(CD[[c]]$reps), BOOT_SEED[[c]]))
    q1 <- readQ1(CD, W, c("NF", "S", "PC", "D", "Dk"))
    pc_pass <- vapply(DECIDING, function(c) all(lowerLim(alphaStat(CD, W, c, "PC", EXCH), 0.025) > BAND), logical(1L))
    checks <- data.frame(check = c("weights: every resample draws n replicates", "planted calibrated arm (D): calibrated",
        "planted 2-fold arm (Dk): miscalibrated in every cell", "NF meets in every cell", "positive control (2-fold): gate passes"),
        pass = c(all(rowSums(W[["M2"]]) == length(CD[["M2"]]$reps)),
            identical(verdictOver(q1$label[q1$arm == "D"], q1$cell[q1$arm == "D"]), "calibrated"),
            all(q1$label[q1$arm == "Dk"] == "fails"), all(q1$label[q1$arm == "NF"] == "meets"), all(pc_pass)),
        stringsAsFactors = FALSE)
    utils::write.table(q1, file.path(out_dir, "selftest-q1.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
    utils::write.table(checks, file.path(out_dir, "selftest.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
    print(checks)
    if (!all(checks$pass)) quit(status = 1L)
}

if (identical(mode, "selftest")) selftest() else score(args[[2L]])
