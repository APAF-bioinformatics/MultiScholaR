# KIT-217 measurement M, the decision script (ticket KIT-217 r1.3 section 4.1; methods ledger VAL-184; the declaration
# docs/qa/workflow-execution/kit/kit-217/m/declaration.md). Reads the counted rows of sim-k217.R and applies every
# reading rule fixed before the counted runs. Evidence code under docs/qa: no package code; needs base R only.
#
# Usage (from the package root):
#   Rscript docs/qa/workflow-execution/kit/kit-217/m/score-k217.R score <rows dir> <out dir>
#   Rscript docs/qa/workflow-execution/kit/kit-217/m/score-k217.R selftest <out dir>
#   Rscript docs/qa/workflow-execution/kit/kit-217/m/score-k217.R e7 <out dir>
# score: <rows dir> holds m-rows-<cell>-<first>-<last>.tsv.gz of every cell and det-rows-L2-<first>-<last>.tsv.gz.
# selftest (lifecycle): the decision rules on planted rows, items (a) to (f) of the declaration audit's D-01, each with
#   its expected label; writes selftest.tsv and exits with status 1 on any mismatch.
# e7 (lifecycle, D-02): scores the lifecycle's rows that sim-k217.R wrote under $RUNNER_TEMP/k217-ml-rows (into
#   $RUNNER_TEMP/k217-ml-score, console discarded) and writes only E7's pass or fail to <out dir>/e7.tsv.

args <- commandArgs(TRUE)
mode <- args[[1L]]
stopifnot(mode %in% c("score", "selftest", "e7"))
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
ND <- "not determined"

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
# the agreement rule: a classification holds only when bootstrap and jackknife give the same one; any non-finite
# bootstrap value or limit makes it "not determined" (D-03 (iv))
agree <- function(labels) if (length(unique(labels)) == 1L) labels[[1L]] else ND
classifyBand <- function(st, edge, good, bad) {
    if (!isTRUE(st$boot_nan == 0)) return(ND)
    up <- upperLim(st, LEVEL_MEETS)
    lo <- lowerLim(st, LEVEL_FAILS)
    agree(vapply(c("boot", "jack"), function(m) if (isTRUE(up[[m]] < edge)) good else if (isTRUE(lo[[m]] >= edge)) bad else ND,
        character(1L)))
}
ratio <- function(S) S[, 1L] / S[, 2L]
diffRatio <- function(S) S[, 1L] / S[, 2L] - S[, 3L] / S[, 4L]

# ---- the rows as arrays [replicate, mask, arm] ------------------------------------------------------------------------
NUMERIC <- c("kept", "kept_novalue", "kept_zerodf", "kept_nogroup", "miss_before", "miss_after", "n_tested", "n_null_tested",
    "n_null_rej", "n_null_low", "n_null_low_rej", "n_null_nolow", "n_null_nolow_rej", "n_planted", "n_planted_tested",
    "n_bh_calls", "n_bh_true", "n_bh_false", "fdp_bh", "n_q_calls", "n_q_true", "n_q_false", "fdp_q", "n_changed_vs_S",
    "sum_fc_null", "n_fc_null", "n_up", "n_down", "rm_up", "rm_down", "n_onoff", "rm_onoff", "nK", "nK_rej_X",
    "nK_NF_tested", "nK_rej_NF", "da_failed")
cellData <- function(rows) {
    reps <- sort(unique(as.integer(rows$rep)))
    masks <- unique(rows$mask)
    arms <- unique(rows$arm)
    idx <- cbind(match(as.integer(rows$rep), reps), match(rows$mask, masks), match(rows$arm, arms))
    ok <- array(FALSE, c(length(reps), length(masks), length(arms)), list(reps, masks, arms))
    # D-03 (i): a unit in which the DA failed for an assay is not complete (KIT-188 r1.9)
    failed <- suppressWarnings(as.numeric(if ("da_failed" %in% names(rows)) rows$da_failed else rep(0, nrow(rows))))
    ok[idx] <- rows$status == "ok" & !(is.finite(failed) & failed > 0)
    values <- lapply(stats::setNames(intersect(NUMERIC, names(rows)), intersect(NUMERIC, names(rows))), function(col) {
        a <- array(NA_real_, dim(ok), dimnames(ok))
        a[idx] <- suppressWarnings(as.numeric(rows[[col]]))
        a
    })
    list(reps = reps, masks = masks, arms = arms, ok = ok, v = values)
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
# the same statistic built per cell: cols = list(c(arm, col), ...), every column paired on the arms named. A paired
# statistic also returns the (replicate, mask) units it used, the units available, and per arm the units it completed
# that were dropped because the other arm did not (D-03 (iii)).
statFor <- function(CD, W, cells, masks, cols, fn, pair = character(), point_only = FALSE) {
    series <- lapply(stats::setNames(cells, cells), function(c) {
        cd <- CD[[c]]
        do.call(cbind, lapply(cols, function(ac) seriesOf(cd, ac[[1L]], masks, ac[[2L]], setdiff(pair, ac[[1L]]))))
    })
    st <- evalStat(series, fn, W[cells], point_only)
    if (length(pair)) {
        used <- 0
        total <- 0
        dropped <- stats::setNames(numeric(length(pair)), pair)
        for (c in cells) {
            m <- intersect(masks, CD[[c]]$masks)
            oks <- lapply(pair, function(a) CD[[c]]$ok[, m, a])
            all_ok <- Reduce(`&`, oks)
            used <- used + sum(all_ok)
            total <- total + length(all_ok)
            for (i in seq_along(pair)) dropped[[i]] <- dropped[[i]] + sum(oks[[i]] & !all_ok)
        }
        st$used <- used
        st$total <- total
        st$dropped <- dropped
    }
    st
}
pairNote <- function(st) sprintf("used %d of %d; dropped %s", st$used, st$total,
    paste(names(st$dropped), st$dropped, sep = " ", collapse = ", "))
pairExcluded <- function(st) (st$total - st$used) / st$total > DNC_SHARE
alphaStat <- function(CD, W, cell, arm, masks) statFor(CD, W, cell, masks, list(c(arm, "n_null_rej"), c(arm, "n_null_tested")), ratio)
excessStat <- function(CD, W, cell, arm, ref, masks) statFor(CD, W, cell, masks, list(c(arm, "n_null_rej"), c(arm, "n_null_tested"),
    c(ref, "n_null_rej"), c(ref, "n_null_tested")), diffRatio, pair = c(arm, ref))
# D-03 (ii): the share not complete over the masks a statistic uses, and the largest share in any one of them
dncShares <- function(cd, arm, masks) {
    masks <- intersect(masks, cd$masks)
    if (!arm %in% cd$arms || !length(masks)) return(c(over = NA_real_, max_mask = NA_real_))
    per <- 1 - colMeans(matrix(cd$ok[, masks, arm], ncol = length(masks)))
    c(over = 1 - mean(cd$ok[, masks, arm]), max_mask = max(per))
}
noResult <- function(cd, arm, masks) {
    s <- dncShares(cd, arm, masks)
    isTRUE(s[["over"]] > DNC_SHARE) || isTRUE(s[["max_mask"]] > DNC_SHARE)
}

# ---- Q1, the guard, the class (declaration sections 6 and 7) ----------------------------------------------------------
readQ1 <- function(CD, W, arms_all) {
    out <- list()
    for (c in DECIDING) {
        st <- alphaStat(CD, W, c, "NF", EXCH)
        nf_meets <- isTRUE(st$boot_nan == 0) && isTRUE(all(upperLim(st, 0.025) < BAND))
        for (arm in intersect(arms_all, CD[[c]]$arms)) {
            a <- alphaStat(CD, W, c, arm, EXCH)
            note <- ""
            if (identical(arm, "NF") || nf_meets) {
                label <- classifyBand(a, BAND, "meets", "fails")
                basis <- "absolute"
            } else {
                e <- excessStat(CD, W, c, arm, "NF", EXCH)
                label <- classifyBand(e, EXCESS_EDGE, "meets", "fails")
                basis <- "excess over NF (NF does not meet)"
                note <- pairNote(e)
            }
            if (noResult(CD[[c]], arm, EXCH)) label <- "no result (did not complete above 5 %)"
            up <- upperLim(a, LEVEL_MEETS)
            lo <- lowerLim(a, LEVEL_FAILS)
            d <- dncShares(CD[[c]], arm, EXCH)
            out[[length(out) + 1L]] <- data.frame(cell = c, arm = arm, alpha = a$est, upper_boot = up[["boot"]], upper_jack = up[["jack"]],
                lower_boot = lo[["boot"]], lower_jack = lo[["jack"]], boot_nan = a$boot_nan, dnc_share = d[["over"]],
                dnc_max_mask = d[["max_mask"]], nf_meets = nf_meets, basis = basis, label = label, paired = note, stringsAsFactors = FALSE)
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
    power <- agree(vapply(c("boot", "jack"), function(m) if (isTRUE(pl[[paste0("lo.", m)]] > POWER_EDGE)) "costs power" else
        if (isTRUE(pl[[paste0("up.", m)]] <= POWER_EDGE)) "no material cost" else ND, character(1L)))
    fdr <- agree(vapply(c("boot", "jack"), function(m) if (isTRUE(fl[[paste0("lo.", m)]] > FDR_EDGE)) "raised" else
        if (isTRUE(fl[[paste0("up.", m)]] <= FDR_EDGE)) "not raised" else ND, character(1L)))
    if (pairExcluded(pw) || !isTRUE(pw$boot_nan == 0)) power <- ND
    if (pairExcluded(fd) || !isTRUE(fd$boot_nan == 0)) fdr <- ND
    data.frame(arm = arm, power_loss_vs_S = pw$est, power_lo_boot = pl[["lo.boot"]], power_up_boot = pl[["up.boot"]],
        power_lo_jack = pl[["lo.jack"]], power_up_jack = pl[["up.jack"]], power = power, fdr_rise_vs_S = fd$est,
        fdr_lo_boot = fl[["lo.boot"]], fdr_up_boot = fl[["up.boot"]], fdr_lo_jack = fl[["lo.jack"]], fdr_up_jack = fl[["up.jack"]],
        fdr = fdr, boot_nan = pw$boot_nan + fd$boot_nan, paired = pairNote(pw), stringsAsFactors = FALSE)
}
classOf <- function(verdict, power, fdr) {
    if (identical(verdict, "calibrated")) return(if (identical(power, "no material cost") && identical(fdr, "not raised")) "E" else "C")
    if (startsWith(verdict, "miscalibrated")) return("R")
    "U"
}

# ---- scoring ------------------------------------------------------------------------------------------------------------
readRows <- function(file) utils::read.delim(gzfile(file), stringsAsFactors = FALSE, colClasses = "character", na.strings = "NA")
# run_all (E7 only): every table is computed even when an earlier outcome row holds; the outcome is still the first row
score <- function(rows_dir, out, run_all = FALSE) {
    dir.create(out, recursive = TRUE, showWarnings = FALSE)
    put <- function(x, name) utils::write.table(x, file.path(out, name), sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
    files <- list.files(rows_dir, pattern = "^m-rows-.*\\.tsv\\.gz$", full.names = TRUE)
    all_rows <- do.call(rbind, lapply(files, readRows))
    cells <- intersect(c(DECIDING, "Lv", "Le", "Me"), unique(all_rows$cell))
    stopifnot(all(DECIDING %in% cells))
    mask_rows <- all_rows[all_rows$arm == "mask", , drop = FALSE]
    arm_rows <- all_rows[all_rows$arm != "mask", , drop = FALSE]
    CD <- lapply(stats::setNames(cells, cells), function(c) cellData(arm_rows[arm_rows$cell == c, , drop = FALSE]))
    W <- lapply(stats::setNames(cells, cells), function(c) bootWeights(length(CD[[c]]$reps), BOOT_SEED[[c]]))
    notes <- character()
    forced <- NULL
    writeOutcome <- function(outcome, recommend, extra = character()) {
        table <- data.frame(item = c("outcome", "recommendation", extra[c(TRUE, FALSE)]), value = c(outcome, recommend,
            extra[c(FALSE, TRUE)]), stringsAsFactors = FALSE)
        put(table, "outcome.tsv")
        print(table)
        invisible(table)
    }
    # ---- stop rules: E1, F1 (counted rows), determinism (range read from the det file's name) ----
    e1_fail <- sum(mask_rows$e1 %in% "FALSE")
    f1_fail <- sum(mask_rows$f1 %in% "FALSE")
    det_file <- list.files(rows_dir, pattern = "^det-rows-L2-[0-9]+-[0-9]+\\.tsv\\.gz$", full.names = TRUE)
    det_ok <- FALSE
    if (length(det_file) == 1L) {
        range <- as.integer(strsplit(sub("^det-rows-L2-([0-9]+)-([0-9]+)\\.tsv\\.gz$", "\\1 \\2", basename(det_file)), " ")[[1L]])
        det <- readRows(det_file)
        first <- all_rows[all_rows$cell == "L2" & as.integer(all_rows$rep) %in% seq(range[[1L]], range[[2L]]), , drop = FALSE]
        key <- function(x) order(as.integer(x$rep), x$mask, x$arm)
        det_ok <- nrow(det) > 0L && identical(unname(as.list(det[key(det), ])), unname(as.list(first[key(first), ])))
    }
    stops <- data.frame(rule = c("E1 (oracle equals the package's S)", "F1 (subset identities)", "determinism (L2)"),
        failures = c(e1_fail, f1_fail, if (det_ok) 0L else 1L), stringsAsFactors = FALSE)
    put(stops, "stops.tsv")
    if (any(stops$failures > 0L)) {
        if (!run_all) return(writeOutcome("Stop", "none: the owner is told; nothing else is read"))
        forced <- c("Stop", "none: the owner is told; nothing else is read")
    }
    # ---- F2 per arm and deciding setting (two-sided 99 % Bonferroni over 80 intervals; floor 0.02 log2) ----
    f2 <- list()
    lv <- 0.01 / 80 / 2
    for (c in DECIDING) for (mask in EXCH) for (arm in intersect(F2_ARMS, CD[[c]]$arms)) {
        st <- statFor(CD, W, c, mask, list(c(arm, "sum_fc_null"), c(arm, "n_fc_null")), ratio)
        lo <- lowerLim(st, lv)
        up <- upperLim(st, lv)
        determined <- isTRUE(st$boot_nan == 0) && all(is.finite(c(lo, up)))
        # a setting whose interval is not determined cannot clear the arm: it is counted as firing (no label)
        fires <- if (!determined) TRUE else any((lo > 0) | (up < 0)) && abs(st$est) > 0.02
        d <- dncShares(CD[[c]], arm, mask)
        f2[[length(f2) + 1L]] <- data.frame(cell = c, mask = mask, arm = arm, mean_fc = st$est, lo_boot = lo[["boot"]], up_boot = up[["boot"]],
            lo_jack = lo[["jack"]], up_jack = up[["jack"]], boot_nan = st$boot_nan, dnc_share = d[["over"]], determined = determined,
            fires = fires, stringsAsFactors = FALSE)
    }
    f2 <- do.call(rbind, f2)
    put(f2, "f2.tsv")
    f2_arms <- unique(f2$arm[f2$fires])
    # ---- Q1 and the positive-control gate ----
    q1 <- readQ1(CD, W, c("NF", "S", "PC", FAMILY))
    put(q1, "q1-cells.tsv")
    pc_pass <- vapply(DECIDING, function(c) {
        row <- q1[q1$cell == c & q1$arm == "NF", ]
        st <- if (isTRUE(row$nf_meets)) alphaStat(CD, W, c, "PC", EXCH) else excessStat(CD, W, c, "PC", "NF", EXCH)
        edge <- if (isTRUE(row$nf_meets)) BAND else EXCESS_EDGE
        isTRUE(st$boot_nan == 0) && isTRUE(all(lowerLim(st, 0.025) > edge)) && !noResult(CD[[c]], "PC", EXCH)
    }, logical(1L))
    gate_text <- paste(DECIDING, ifelse(pc_pass, "pass", "fail"), collapse = "; ")
    if (!all(pc_pass) && is.null(forced)) {
        if (!run_all) return(writeOutcome("N0", "none: the positive control does not register a 2-fold excess; no candidate is labelled",
            c("positive-control gate", gate_text)))
        forced <- c("N0", "none: the positive control does not register a 2-fold excess; no candidate is labelled")
    }
    # O0: the reading cells are those where S meets
    s_labels <- stats::setNames(q1$label[q1$arm == "S"], q1$cell[q1$arm == "S"])
    reading <- DECIDING[s_labels[DECIDING] %in% "meets"]
    if (!length(reading) && is.null(forced)) {
        if (!run_all) return(writeOutcome("O0 (no reading cell)", "none: S meets the band in no deciding cell; the owner decides on the per-cell figures (q1-cells.tsv)",
            c("positive-control gate", gate_text)))
        forced <- c("O0 (no reading cell)", "none")
    }
    # verdicts, guard, class
    classes <- do.call(rbind, lapply(FAMILY, function(arm) {
        lab <- stats::setNames(q1$label[q1$arm == arm], q1$cell[q1$arm == arm])[reading]
        verdict <- verdictOver(unname(lab), reading)
        g <- if (length(reading)) guardFor(CD, W, reading, arm) else data.frame(arm = arm, power = ND, fdr = ND)
        cls <- if (arm %in% f2_arms) "no label (F2)" else classOf(verdict, g$power, g$fdr)
        cbind(data.frame(arm = arm, verdict = verdict, class = cls, stringsAsFactors = FALSE), g[, setdiff(names(g), "arm"), drop = FALSE])
    }))
    put(classes, "classes.tsv")
    # ---- Q2 (reported, SD-10 (ii-b)) ----
    q2 <- do.call(rbind, lapply(DECIDING, function(c) do.call(rbind, lapply(c("S", "B", FAMILY), function(arm) {
        e <- excessStat(CD, W, c, arm, "NF", GLINK)
        a <- alphaStat(CD, W, c, arm, GLINK)
        label <- classifyBand(e, EXCESS_EDGE, "robust", "not robust")
        if (pairExcluded(e)) label <- ND
        if (noResult(CD[[c]], arm, GLINK)) label <- "no result (did not complete above 5 %)"
        data.frame(cell = c, arm = arm, alphaG = a$est, alphaG_NF = alphaStat(CD, W, c, "NF", GLINK)$est, excess = e$est,
            boot_nan = e$boot_nan, label = label, paired = pairNote(e), stringsAsFactors = FALSE)
    }))))
    put(q2, "q2.tsv")
    q2v <- vapply(c("S", "B", FAMILY), function(arm) {
        l <- q2$label[q2$arm == arm]
        if (all(l == "robust")) "robust" else if (any(l == "not robust")) "not robust" else ND
    }, character(1L))
    ad <- do.call(rbind, lapply(DECIDING, function(c) {
        st <- statFor(CD, W, c, GLINK, list(c("Aobs_g2", "n_null_rej"), c("Aobs_g2", "n_null_tested"), c("Dk", "n_null_rej"),
            c("Dk", "n_null_tested")), diffRatio, pair = c("Aobs_g2", "Dk"))
        data.frame(cell = c, eG_Aobs_g2_minus_Dk = st$est, lo95_boot = lowerLim(st, 0.025)[["boot"]], up95_boot = upperLim(st, 0.025)[["boot"]],
            lo95_jack = lowerLim(st, 0.025)[["jack"]], up95_jack = upperLim(st, 0.025)[["jack"]], boot_nan = st$boot_nan,
            paired = pairNote(st), stringsAsFactors = FALSE)
    }))
    put(ad, "q2-aobs-dk.tsv")
    q2_flag <- any(q2$label[q2$arm == "Aobs_g2"] == "not robust") && all(q2$label[q2$arm == "Dk"] == "robust")
    # ---- outcome table (the first row that holds applies) ----
    cls <- stats::setNames(classes$class, classes$arm)
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
            notes <- c(notes, paste("O2 power tie-break:", pairNote(st)))
            lo <- lowerLim(st, 0.025)
            up <- upperLim(st, 0.025)
            if (isTRUE(all(lo > 0))) "Aobs_any" else if (isTRUE(all(up < 0))) "Dk" else {
                changed <- vapply(e_arms, function(a) statFor(CD, W, reading, EXCH, list(c(a, "n_changed_vs_S"), c(a, "one")), ratio,
                    point_only = TRUE)$est, numeric(1L))
                if (isTRUE(changed[["Aobs_any"]] < changed[["Dk"]])) "Aobs_any" else "Dk"
            }
        }
        recommend <- paste0(recommend, " (the owner chooses)")
    } else {
        code <- "O3"
        recommend <- "B"
    }
    outcome <- paste(c(if (outcome == "O0") "O0", code), collapse = " then ")
    if (length(f2_arms)) notes <- c(notes, paste("F2 fired or not determined; no label for:", paste(f2_arms, collapse = ", ")))
    if (q2_flag) notes <- c(notes, "SD-10 (ii-b): A-obs-g2 is 'not robust' in Q2 in some cell and D-k is 'robust' in every cell: the owner is told before choosing; no class changes")
    # ---- large effects (Le, Me): the recommended candidate against D-k and B, per mask ----
    rec_arm <- base::sub(" .*$", "", recommend)
    large <- list()
    for (c in intersect(c("Le", "Me"), cells)) for (mask in CD[[c]]$masks) for (arm in CD[[c]]$arms) {
        point <- function(a, b) statFor(CD, W, c, mask, list(c(arm, a), c(arm, b)), ratio, point_only = TRUE)$est
        large[[length(large) + 1L]] <- data.frame(cell = c, mask = mask, arm = arm, power_bh = point("n_bh_true", "n_planted"),
            removed_up = point("rm_up", "n_planted"), removed_down = point("rm_down", "n_planted"),
            removed_onoff_share = point("rm_onoff", "n_onoff"), dnc_share = dncShares(CD[[c]], arm, mask)[["over"]],
            stringsAsFactors = FALSE)
        if (identical(arm, rec_arm)) for (ref in c("Dk", "B")) if (!identical(ref, arm)) {
            both <- statFor(CD, W, c, mask, list(c(arm, "rm_up"), c(arm, "rm_down"), c(arm, "n_planted"), c(ref, "rm_up"),
                c(ref, "rm_down"), c(ref, "n_planted")), function(S) (S[, 1L] + S[, 2L]) / S[, 3L] - (S[, 4L] + S[, 5L]) / S[, 6L],
                pair = c(arm, ref))
            lo <- lowerLim(both, 0.025)
            if (isTRUE(all(lo > 0.05))) notes <- c(notes, sprintf("large effects: %s removes a larger share of planted features than %s in %s %s (difference %.3f, 95 %% lower limits %.3f / %.3f; %s)",
                arm, ref, c, mask, both$est, lo[["boot"]], lo[["jack"]], pairNote(both)))
        }
    }
    if (length(large)) put(do.call(rbind, large), "large-effects.tsv")
    # ---- printed beside, never deciding ----
    beside <- list()
    for (c in cells) {
        masks_c <- if (c %in% DECIDING) EXCH else CD[[c]]$masks
        for (arm in CD[[c]]$arms) {
            a <- alphaStat(CD, W, c, arm, masks_c)
            vsS <- if (arm %in% c("S", "NF") || !"S" %in% CD[[c]]$arms) NULL else excessStat(CD, W, c, arm, "S", masks_c)
            point <- function(cols, fn = ratio) statFor(CD, W, c, masks_c, cols, fn, point_only = TRUE)$est
            mean_of <- function(col) point(list(c(arm, col), c(arm, "one")))
            low <- point(list(c(arm, "n_null_low_rej"), c(arm, "n_null_low")))
            nolow <- point(list(c(arm, "n_null_nolow_rej"), c(arm, "n_null_nolow")))
            d <- dncShares(CD[[c]], arm, masks_c)
            beside[[length(beside) + 1L]] <- data.frame(cell = c, masks = paste(masks_c, collapse = ","), arm = arm, alpha = a$est,
                alpha_lo95 = lowerLim(a, 0.025)[["boot"]], alpha_up95 = upperLim(a, 0.025)[["boot"]], boot_nan = a$boot_nan,
                alpha_minus_S = if (is.null(vsS)) NA else vsS$est,
                alpha_minus_S_lo95 = if (is.null(vsS)) NA else lowerLim(vsS, 0.025)[["boot"]],
                alpha_minus_S_up95 = if (is.null(vsS)) NA else upperLim(vsS, 0.025)[["boot"]],
                paired_vs_S = if (is.null(vsS)) NA else pairNote(vsS),
                power_bh = point(list(c(arm, "n_bh_true"), c(arm, "n_planted"))), fdr_bh = mean_of("fdp_bh"),
                power_q = point(list(c(arm, "n_q_true"), c(arm, "n_planted"))), fdr_q = mean_of("fdp_q"), alpha_low = low,
                alpha_nolow = nolow, ratio_low = low / nolow,
                fixed_set_ratio = point(list(c(arm, "nK_rej_X"), c(arm, "nK"), c(arm, "nK_rej_NF"), c(arm, "nK_NF_tested")),
                    function(S) (S[, 1L] / S[, 2L]) / (S[, 3L] / S[, 4L])),
                kept = mean_of("kept"), kept_novalue = mean_of("kept_novalue"), kept_zerodf = mean_of("kept_zerodf"),
                kept_nogroup = mean_of("kept_nogroup"), miss_before = mean_of("miss_before"), miss_after = mean_of("miss_after"),
                changed_calls_vs_S = mean_of("n_changed_vs_S"), mean_fc_null = point(list(c(arm, "sum_fc_null"), c(arm, "n_fc_null"))),
                dnc_share = d[["over"]], dnc_max_mask = d[["max_mask"]], stringsAsFactors = FALSE)
        }
    }
    put(do.call(rbind, beside), "beside.tsv")
    # per-setting null rate, printed
    setting <- do.call(rbind, lapply(cells, function(c) do.call(rbind, lapply(CD[[c]]$masks, function(mask) do.call(rbind,
        lapply(CD[[c]]$arms, function(arm) {
            a <- alphaStat(CD, W, c, arm, mask)
            data.frame(cell = c, mask = mask, arm = arm, alpha = a$est, lo95 = lowerLim(a, 0.025)[["boot"]],
                up95 = upperLim(a, 0.025)[["boot"]], boot_nan = a$boot_nan, dnc_share = dncShares(CD[[c]], arm, mask)[["over"]],
                stringsAsFactors = FALSE)
        }))))))
    put(setting, "per-setting.tsv")
    if (!is.null(forced)) {
        outcome <- forced[[1L]]
        recommend <- forced[[2L]]
    }
    writeOutcome(outcome, recommend, c("positive-control gate", gate_text, "reading cells", paste(reading, collapse = ", "),
        "Q2 verdicts", paste(names(q2v), q2v, sep = ": ", collapse = "; "), "notes", paste(notes, collapse = " | ")))
}

# ---- selftest: planted rows (declaration audit D-01 (a) to (f)) --------------------------------------------------------
PLANT_R <- c(L2 = 200L, L1 = 200L, M2 = 200L, M1 = 200L, Le = 100L)
# over: list of list(cell = "<cell>" or "*", arm, and any of rate, power, fdp, changed); defaults: null rate 0.05 (PC
# 0.10), BH power 0.60, FDP 0.05, changed calls 0, mean fold change 0. Null rejections are binomial on a replicate-level
# rate rate x exp(N(0, 0.15^2) - 0.15^2 / 2); power, FDP and changed calls are fixed values.
plantRows <- function(over = list(), seed = 217990100L) {
    set.seed(seed, kind = "Mersenne-Twister", normal.kind = "Inversion", sample.kind = "Rejection")
    do.call(rbind, lapply(names(PLANT_R), function(c) {
        lipid <- startsWith(c, "L")
        arms <- c("NF", "S", "PC", "B", "D", "Dk", "Atau_any", "Atau_g2", "Aobs_any", "Aobs_g2", if (!lipid) c("Atau_all", "Aobs_all"))
        masks <- if (identical(c, "Le")) c("vd10", "vd25", "vdg25") else c("none", EXCH, GLINK)
        g <- expand.grid(arm = arms, mask = masks, rep = seq_len(PLANT_R[[c]]), stringsAsFactors = FALSE)
        field <- function(name, default) {
            v <- default
            for (o in over) if (!is.null(o[[name]]) && (o$cell %in% c("*", c))) v[g$arm == o$arm] <- o[[name]]
            v
        }
        n_null <- if (lipid) 360 else 100
        n_planted <- if (lipid) 100 else 30
        rate <- field("rate", ifelse(g$arm == "PC", 0.10, 0.05))
        s <- 0.15
        p <- pmin(1, rate * exp(stats::rnorm(nrow(g), 0, s) - s^2 / 2))
        arm_rows <- data.frame(cell = c, rep = g$rep, mask = g$mask, arm = g$arm, status = "ok", class = NA_character_,
            stringsAsFactors = FALSE)
        for (col in NUMERIC) arm_rows[[col]] <- 0
        arm_rows$kept <- n_null + n_planted
        arm_rows$n_tested <- n_null + n_planted
        arm_rows$n_null_tested <- n_null
        arm_rows$n_null_rej <- stats::rbinom(nrow(g), n_null, p)
        arm_rows$n_planted <- n_planted
        arm_rows$n_bh_true <- round(field("power", rep(0.60, nrow(g))) * n_planted)
        arm_rows$fdp_bh <- field("fdp", rep(0.05, nrow(g)))
        arm_rows$n_changed_vs_S <- field("changed", rep(0, nrow(g)))
        arm_rows$n_fc_null <- n_null
        arm_rows$e1 <- NA
        arm_rows$f1 <- NA
        mg <- unique(g[, c("mask", "rep")])
        mask_rows <- arm_rows[seq_len(nrow(mg)), , drop = FALSE]
        mask_rows[] <- NA
        mask_rows$cell <- c
        mask_rows$rep <- mg$rep
        mask_rows$mask <- mg$mask
        mask_rows$arm <- "mask"
        mask_rows$status <- "ok"
        mask_rows$e1 <- TRUE
        mask_rows$f1 <- TRUE
        rbind(mask_rows, arm_rows)
    }))
}
writePlanted <- function(rows, dir) {
    dir.create(dir, recursive = TRUE, showWarnings = FALSE)
    w <- function(x, name) {
        connection <- gzfile(file.path(dir, name), "w")
        utils::write.table(x, connection, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
        close(connection)
    }
    for (c in unique(rows$cell)) w(rows[rows$cell == c, , drop = FALSE], sprintf("m-rows-%s-1-%d.tsv.gz", c, PLANT_R[[c]]))
    w(rows[rows$cell == "L2" & rows$rep %in% 1:5, , drop = FALSE], "det-rows-L2-1-5.tsv.gz")
    dir
}
selftest <- function() {
    checks <- list()
    check <- function(name, pass) checks[[length(checks) + 1L]] <<- data.frame(check = name, pass = isTRUE(pass), stringsAsFactors = FALSE)
    prepare <- function(rows) {
        arm_rows <- rows[rows$arm != "mask" & rows$cell %in% DECIDING, , drop = FALSE]
        CD <- lapply(stats::setNames(DECIDING, DECIDING), function(c) cellData(arm_rows[arm_rows$cell == c, , drop = FALSE]))
        W <- lapply(stats::setNames(DECIDING, DECIDING), function(c) bootWeights(length(CD[[c]]$reps), BOOT_SEED[[c]]))
        list(CD = CD, W = W)
    }
    # set A: (a), (b), (c), (e)
    A <- prepare(plantRows(list(list(cell = "*", arm = "Dk", rate = 0.10), list(cell = "M2", arm = "Atau_g2", rate = 0.10),
        list(cell = "*", arm = "Aobs_any", rate = 0.075), list(cell = "*", arm = "Aobs_any", power = 0.57),
        list(cell = "*", arm = "Dk", fdp = 0.06))))
    check("weights: every resample draws n replicates", all(rowSums(A$W[["M2"]]) == length(A$CD[["M2"]]$reps)))
    q1 <- readQ1(A$CD, A$W, c("NF", "S", "PC", FAMILY))
    verdictOf <- function(q, arm) verdictOver(q$label[q$arm == arm], q$cell[q$arm == arm])
    check("(a) planted calibrated arm D reads calibrated", identical(verdictOf(q1, "D"), "calibrated"))
    check("(a) planted 2-fold arm Dk fails in every cell", all(q1$label[q1$arm == "Dk"] == "fails"))
    check("(a) NF meets in every cell", all(q1$label[q1$arm == "NF"] == "meets"))
    check("(a) the positive-control gate passes", all(vapply(DECIDING, function(c) all(lowerLim(alphaStat(A$CD, A$W, c, "PC", EXCH), 0.025) > BAND),
        logical(1L))))
    check("(b) Atau_g2 at 0.10 in M2 only reads miscalibrated (M2)", identical(verdictOf(q1, "Atau_g2"), "miscalibrated (M2)"))
    check("(c) Aobs_any at 0.075 reads not determined in every cell", identical(verdictOf(q1, "Aobs_any"), "not determined (L2, L1, M2, M1)"))
    check("(e) Aobs_g2 at S's power: no material cost", identical(guardFor(A$CD, A$W, DECIDING, "Aobs_g2")$power, "no material cost"))
    check("(e) Aobs_any 0.03 below S's power: costs power", identical(guardFor(A$CD, A$W, DECIDING, "Aobs_any")$power, "costs power"))
    check("(e) Dk FDP 0.06 against S's 0.05: false share raised", identical(guardFor(A$CD, A$W, DECIDING, "Dk")$fdr, "raised"))
    rm(A)
    # set B: (d) NF at 0.12 in M2 only
    Bset <- prepare(plantRows(list(list(cell = "M2", arm = "NF", rate = 0.12), list(cell = "M2", arm = "Aobs_g2", rate = 0.12),
        list(cell = "M2", arm = "Dk", rate = 0.18)), seed = 217990101L))
    qb <- readQ1(Bset$CD, Bset$W, c("NF", "S", "PC", FAMILY))
    m2 <- qb[qb$cell == "M2", ]
    check("(d) M2 is read on the excess over NF", all(m2$basis[m2$arm != "NF"] == "excess over NF (NF does not meet)"))
    check("(d) an arm at NF's rate meets in M2", identical(m2$label[m2$arm == "Aobs_g2"], "meets"))
    check("(d) an arm at NF + 0.06 fails in M2", identical(m2$label[m2$arm == "Dk"], "fails"))
    check("(d) the other cells are read absolutely", all(qb$basis[qb$cell != "M2"] == "absolute"))
    rm(Bset)
    # (f) score() end to end on planted directories (console discarded)
    family_at <- function(rate) lapply(FAMILY, function(a) list(cell = "*", arm = a, rate = rate))
    scenarios <- list(
        stop = list(over = list(), e1 = TRUE, outcome = "Stop", recommend = NULL),
        n0 = list(over = list(list(cell = "*", arm = "PC", rate = 0.05)), outcome = "N0", recommend = NULL),
        o0o1 = list(over = list(list(cell = "M2", arm = "S", rate = 0.10)), outcome = "O0 then O1", recommend = "Aobs_g2"),
        o1 = list(over = list(), outcome = "O1", recommend = "Aobs_g2"),
        o2any = list(over = list(list(cell = "*", arm = "Aobs_g2", rate = 0.10), list(cell = "*", arm = "Aobs_any", power = 0.65)),
            outcome = "O2", recommend = "Aobs_any (the owner chooses)"),
        o2dk = list(over = list(list(cell = "*", arm = "Aobs_g2", rate = 0.10), list(cell = "*", arm = "Aobs_any", changed = 5),
            list(cell = "*", arm = "Dk", changed = 2)), outcome = "O2", recommend = "Dk (the owner chooses)"),
        o3 = list(over = family_at(0.10), outcome = "O3", recommend = "B"))
    for (name in names(scenarios)) {
        sc <- scenarios[[name]]
        rows <- plantRows(sc$over, seed = 217990110L)
        if (isTRUE(sc$e1)) rows$e1[which(rows$arm == "mask")[[1L]]] <- FALSE
        dir <- writePlanted(rows, file.path(tempdir(), paste0("k217-selftest-", name)))
        res_dir <- file.path(tempdir(), paste0("k217-selftest-", name, "-out"))
        result <- NULL
        invisible(utils::capture.output(result <- tryCatch(score(dir, res_dir), error = function(e) NULL)))
        got <- if (is.null(result)) c(NA, NA) else result$value[match(c("outcome", "recommendation"), result$item)]
        check(sprintf("(f) %s: outcome %s", name, sc$outcome), identical(got[[1L]], sc$outcome))
        if (!is.null(sc$recommend)) check(sprintf("(f) %s: recommends %s", name, sc$recommend), identical(got[[2L]], sc$recommend))
    }
    checks <- do.call(rbind, checks)
    utils::write.table(checks, file.path(out_dir, "selftest.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
    print(checks)
    if (!all(checks$pass)) quit(status = 1L)
}

# ---- E7: the lifecycle's rows through score(), blind (declaration audit D-02) ------------------------------------------
e7 <- function() {
    temp <- Sys.getenv("RUNNER_TEMP")
    if (!nzchar(temp)) temp <- tempdir()
    rows_dir <- file.path(temp, "k217-ml-rows")
    score_dir <- file.path(temp, "k217-ml-score")
    ran <- FALSE
    invisible(utils::capture.output(ran <- tryCatch({
        score(rows_dir, score_dir, run_all = TRUE)
        TRUE
    }, error = function(e) FALSE)))
    wanted <- c("stops.tsv", "outcome.tsv", "q1-cells.tsv", "classes.tsv", "q2.tsv", "f2.tsv", "beside.tsv", "per-setting.tsv")
    present <- all(file.exists(file.path(score_dir, wanted)))
    det_ok <- present && isTRUE(utils::read.delim(file.path(score_dir, "stops.tsv"))$failures[[3L]] == 0)
    e7 <- data.frame(check = "E7", pass = isTRUE(ran) && present && det_ok, ran = isTRUE(ran), tables = present,
        determinism = det_ok, stringsAsFactors = FALSE)
    utils::write.table(e7, file.path(out_dir, "e7.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
    if (!e7$pass) quit(status = 1L)
}

if (identical(mode, "selftest")) selftest() else if (identical(mode, "e7")) e7() else score(args[[2L]], out_dir)
