# KIT-211 (VAL-142 (2), VAL-171, VAL-315): the missing-value filter and the limma test on values at or below the intensity
# threshold: arms S (shipped), A, B, B' (candidates), N and C0 (references); the cells A to D, Av (reported) and Ag (reported).
# The simulation (ticket sections 6.2, 7.1, 18.3). Evidence code under docs/qa: no package code.
#
# Extends docs/qa/workflow-execution/kit/kit-197/sim-k197.R (SHA-256 71001de9b7b1c547e1e27fd92fb1b57a067121e336fa09e134f3317f1fa7c599,
# merged at ef63c918): the block between the first "verbatim" and the last "end verbatim" marker below is KIT-197's lines 26 to
# 279 unchanged, with exactly two lines added (checked by e0-k211.py): the guarded variance line of cell Av after the variance
# draw, and the guarded planting line of cell Ag directly after `beta[where] <- signs`.
#
# Usage, from the package root of the bound commit on a run-tests runner (full git history):
#   Rscript docs/qa/workflow-execution/kit/kit-211/sim-k211.R <out dir> <first r> <last r> [<replicate list>]
#   <replicate list>  optional file of "<cell> <r>" lines: only those tasks run (a continuation job); the
#                     part files already in <out dir>/parts are kept.
# Environment: K211_LIFECYCLE=1 (checks E2, E3, E7, E9, E10 on every replicate-mask), K211_DETERMINISM=1 (after the main run,
#   cell A replicates <first r> to <first r> + 4 are run again and compared), K211_WORKERS (forked workers; default 4),
#   K211_CELLS (comma list of cells; default A,B,C,D), K211_ARMS (comma list restricting the arms: E1 runs S only),
#   K211_PARTIAL=1 (assemble the tasks that completed when the range is incomplete after its continuation).
# One tree: pkgload::load_all(".", export_all = TRUE). One task per (cell, r) running every mask and arm; each task writes
# parts/<cell>-<r>.rows.tsv and parts/<cell>-<r>.metrics.tsv through a `.tmp` name, the metrics file last, so a task
# whose metrics file exists is complete. When every task of the range is complete the parts are assembled into
# <cell>.tsv.gz (the sufficient statistics), metrics.tsv (per cell, mask, r and arm) and range.tsv, and the parts are removed.
# Seeds: the generator's set.seed(185000000 + 10000 * idx + r), idx 201 to 204 (A to D), 211 (Av), 212 (Ag).

# ==== verbatim: filter-effect.R:23-234 at 6416dbd6 (line 215 changed: logFC, t, ave kept) ====
# ---- verbatim: sim-k185.R:43-52 at 432f703b ----
warning_count <- new.env()
warning_count$n <- 0L
quiet <- function(expr) {
    value <- NULL
    invisible(utils::capture.output(value <- withCallingHandlers(suppressMessages(expr), warning = function(w) {
        warning_count$n <- warning_count$n + 1L
        invokeRestart("muffleWarning")
    }), type = "output"))
    value
}

# ---- verbatim: sim-k185.R:71-71 at 432f703b ----
GROUPS <- c("control", "treated", "g3", "g4", "g5")

# ---- verbatim: sim-k185.R:72-76 at 432f703b ----
cell <- function(idx, omic, features, groups, per, k = 2L, level = 20, rho = 0, dup = FALSE, batches = 0L,
                 batch_runs = 0L, factors = 2L, assays = 1L) {
    list(idx = idx, omic = omic, features = features, groups = groups, per = per, k = k, level = level, rho = rho,
        dup = dup, batches = batches, batch_runs = batch_runs, factors = factors, assays = assays)
}

# ---- verbatim: sim-k185.R:90-91 at 432f703b ----
MED_VAR_PROT <- 4 * 0.25^2 / stats::qchisq(0.5, 4)
MED_VAR_OMIC <- 0.25

# ---- verbatim: sim-k185.R:101-108 at 432f703b ----
maskMech <- function(mask) if (identical(mask, "none")) "none" else sub("[0-9.]+$", "", mask)
maskRate <- function(mask) if (identical(mask, "none")) 0 else as.numeric(sub("^[a-z]+", "", mask)) / 100

# ---- data (KIT-185 5.2) --------------------------------------------------------------------------------------------
scale1 <- function(v) {
    v <- v - mean(v)
    v / stats::sd(v)
}

# ---- verbatim: sim-k185.R:109-199 at 432f703b ----
generateReplicate <- function(cellname, r, strength = 1.5) {
    p <- CELLS[[cellname]]
    set.seed(185000000L + 10000L * p$idx + r, kind = "Mersenne-Twister", normal.kind = "Inversion", sample.kind = "Rejection")
    n <- p$features
    # the runs: groups, units (P5), batches (D5 cells)
    if (p$batches > 0L) {
        per_batch <- if (identical(cellname, "P4b")) 2L else if (identical(cellname, "L1b")) 16L else 3L
        batch <- rep(seq_len(p$batches), each = p$groups * per_batch)
        group <- rep(rep(GROUPS[seq_len(p$groups)], each = per_batch), times = p$batches)
    } else {
        batch <- NULL
        group <- rep(GROUPS[seq_len(p$groups)], each = p$per * if (p$dup) 2L else 1L)
    }
    runs <- length(group)
    unit <- if (p$dup) rep(seq_len(p$groups * p$per), each = 2L) else seq_len(runs)
    run_ids <- sprintf("run_%03d", seq_len(runs))
    group_factor <- factor(group, levels = GROUPS[seq_len(p$groups)])
    median_variance <- if (identical(p$omic, "prot")) MED_VAR_PROT else MED_VAR_OMIC
    # levels
    mu <- stats::rnorm(n, p$level, 2)
    # variances
    s2 <- if (identical(p$omic, "prot")) 4 * 0.25^2 / stats::rchisq(n, 4) else {
        x <- 10 / stats::rchisq(n, 10)
        0.25 * x / stats::median(x)
    }
    if (identical(cellname, "Av")) s2 <- s2 * 2^((20 - mu) / 2)
    # effects: exactly 20 % of the features, random positions, half up, in group 2
    n_effects <- round(0.2 * n)
    where <- sample.int(n, n_effects)
    signs <- sample(rep(c(1, -1), c(round(n_effects / 2), n_effects - round(n_effects / 2))))
    beta <- numeric(n)
    beta[where] <- signs
    if (identical(cellname, "Ag")) beta <- 3 * beta
    # factors: centred, sample variance 1, exactly orthogonal to each other and to the group design (and batches)
    basis <- stats::model.matrix(~group_factor)
    planted <- NULL
    if (p$batches > 0L) {
        dummies <- stats::model.matrix(~ factor(batch))[, -1L, drop = FALSE]
        dummies <- scale(dummies, center = TRUE, scale = FALSE)
        basis <- cbind(basis, dummies)
        colnames(dummies) <- paste0("B", seq_len(ncol(dummies)))
        planted <- dummies
    }
    g <- scale1(as.numeric(group == "treated"))
    factor_scores <- NULL
    for (f in seq_len(p$factors)) {
        raw <- if (!identical(p$omic, "prot") && f == 1L) {
            order <- sample.int(runs)
            exp(-2 * (order - 1) / (runs - 1))
        } else stats::rnorm(runs)
        e <- scale1(stats::lm.fit(basis, raw)$residuals)
        basis <- cbind(basis, e)
        z <- if (f == 1L && p$rho > 0) p$rho * g + sqrt(1 - p$rho^2) * e else e
        factor_scores <- cbind(factor_scores, z)
    }
    if (!is.null(factor_scores)) colnames(factor_scores) <- paste0("F", seq_len(ncol(factor_scores)))
    planted <- cbind(planted, factor_scores)
    rownames(planted) <- run_ids
    # loadings
    tau <- sqrt(strength * median_variance)
    loadings <- if (p$factors > 0L) matrix(stats::rnorm(n * p$factors), n, p$factors) * tau
    offsets <- NULL
    if (p$batches > 0L) {
        offsets <- matrix(stats::rnorm(n * p$batches), n, p$batches) * sqrt(3 * median_variance)
        offsets <- offsets - rowMeans(offsets)
    }
    # residuals
    noise <- if (p$dup) {
        units <- max(unit)
        biological <- matrix(stats::rnorm(n * units), n, units) * sqrt(0.75 * s2)
        technical <- matrix(stats::rnorm(n * runs), n, runs) * sqrt(0.25 * s2)
        biological[, unit, drop = FALSE] + technical
    } else matrix(stats::rnorm(n * runs), n, runs) * sqrt(s2)
    # the uniforms of the masks, after the values (5.2); the batch-block draw W after U and V (B04)
    U <- matrix(stats::runif(n * runs), n, runs)
    V <- matrix(stats::runif(n * runs), n, runs)
    W <- if (p$batches > 0L) matrix(stats::runif(n * p$batches), n, p$batches)
    unwanted <- matrix(0, n, runs)
    if (p$factors > 0L) unwanted <- unwanted + loadings %*% t(factor_scores)
    if (p$batches > 0L) unwanted <- unwanted + offsets[, batch, drop = FALSE]
    y <- mu + outer(beta, as.numeric(group == "treated")) + unwanted + noise
    ids <- switch(p$omic, prot = sprintf("PG_%05d", seq_len(n)),
        lipid = c(sprintf("A1_%04d", seq_len(n / 2)), sprintf("A2_%04d", seq_len(n / 2))), metab = sprintf("M_%04d", seq_len(n)))
    assay <- switch(p$omic, lipid = rep(c("pos", "neg"), each = n / 2), rep("pos", n))
    dimnames(y) <- dimnames(unwanted) <- dimnames(U) <- dimnames(V) <- list(ids, run_ids)
    design <- if (identical(p$omic, "prot")) {
        data.frame(Run = run_ids, group = group, replicates = seq_len(runs), tech_rep_group = paste0("u", unit),
            stringsAsFactors = FALSE)
    } else data.frame(Sample_ID = run_ids, group = group, extract = paste0("e", seq_len(runs)), stringsAsFactors = FALSE)
    list(cell = cellname, p = p, r = r, y = y, unwanted = unwanted, beta = stats::setNames(beta, ids), U = U, V = V,
        W = W, batch = batch, group = group, unit = unit, runs = run_ids, Z = planted, design = design, assay = stats::setNames(assay, ids),
        strength = strength)
}

# ---- verbatim: sim-k185.R:202-234 at 432f703b ----
solveLod <- function(values, rate, shift = 0) {
    stats::uniroot(function(l) mean(stats::plogis(l + shift - values)) - rate, c(min(values) - 60, max(values) + 60),
        tol = 1e-10)$root
}
maskMatrix <- function(sim, mask) {
    y <- sim$y
    if (identical(mask, "none")) return(matrix(FALSE, nrow(y), ncol(y), dimnames = dimnames(y)))
    rate <- maskRate(mask)
    switch(maskMech(mask),
        vd = sim$U < stats::plogis(solveLod(y, rate) - y),
        mcar = sim$V < rate,
        bb = {
            block_mean <- vapply(seq_len(max(sim$batch)), function(b) rowMeans(y[, sim$batch == b, drop = FALSE]), numeric(nrow(y)))
            lod <- solveLod(block_mean, rate)
            block_missing <- sim$W < stats::plogis(lod - block_mean)
            out <- block_missing[, sim$batch, drop = FALSE]
            dimnames(out) <- dimnames(y)
            out
        },
        ra = {
            z1 <- sim$Z[, "F1"]
            top <- which.max(z1)
            shift_of <- function(c) matrix(c * z1, nrow(y), ncol(y), byrow = TRUE)
            share_top <- function(c) {
                lod <- solveLod(y, rate, shift_of(c))
                mean(stats::plogis(lod + c * z1[[top]] - y[, top]))
            }
            upper <- 1
            while (share_top(upper) < 2 * rate && upper < 1e3) upper <- upper * 2
            c_value <- stats::uniroot(function(c) share_top(c) - 2 * rate, c(0, upper), tol = 1e-10)$root
            sim$U < stats::plogis(solveLod(y, rate, shift_of(c_value)) + shift_of(c_value) - y)
        })
}

# ---- verbatim: sim-k185.R:237-244 at 432f703b ----
protObject <- function(mat, design) {
    quant <- as.data.frame(c(list(Protein.Group = rownames(mat)), as.data.frame(mat, check.names = FALSE)),
        check.names = FALSE, stringsAsFactors = FALSE)
    methods::new("ProteinQuantitativeData", protein_quant_table = quant, protein_id_column = "Protein.Group",
        protein_id_table = data.frame(Protein.Group = rownames(mat), stringsAsFactors = FALSE), design_matrix = design,
        sample_id = "Run", group_id = "group", technical_replicate_id = "tech_rep_group",
        args = list(globalParameters = list(workflow_type = "LFQ")))
}

# ---- verbatim: sim-k185.R:429-465 at 432f703b ----
CONTRAST <- "grouptreated-groupcontrol"
qMethod <- function(p) {
    p <- p[!is.na(p) & is.finite(p)]
    if (!length(p) || all(p == 0) || all(p == 1) || length(p) < 3L) return("BH")
    if (inherits(tryCatch(qvalue::qvalue(p), error = function(e) e), "error")) "BH" else "qvalue"
}
daTable <- function(sim, mat = NULL, object = NULL, ruv_factors = NULL) {
    omic <- sim$p$omic
    before <- warning_count$n
    if (identical(omic, "prot")) {
        object <- object %||% protObject(mat, sim$design)
        value <- quiet(differentialAbundanceAnalysis(object, contrasts_tbl = data.frame(contrasts = paste0("treated_vs_control=", CONTRAST),
            stringsAsFactors = FALSE), formula_string = "~ 0 + group", da_q_val_thresh = 0.05, treat_lfc_cutoff = 0,
            eBayes_trend = TRUE, eBayes_robust = TRUE, args_row_id = "Protein.Group", qvalue_column = "fdr_qvalue",
            raw_pvalue_column = "raw_pvalue"))
        tab <- value$contrasts_results$results[[1L]]
        if (is.null(tab) || !nrow(tab)) stop("differentialAbundanceAnalysis returned no table")
        out <- data.frame(id = rownames(tab), p = tab$raw_pvalue, q = tab$fdr_qvalue, logFC = tab$logFC, t = tab$t, ave = tab$AveExpr, stringsAsFactors = FALSE)
        attr(out, "qmethod") <- qMethod(out$p)
        attr(out, "failed") <- 0L
    } else {
        object <- object %||% omicObject(omic, mat, sim$assay, sim$design)
        contrasts <- data.frame(contrasts = CONTRAST, friendly_names = "treated_vs_control", stringsAsFactors = FALSE)
        value <- quiet(if (identical(omic, "lipid")) {
            runLipidsDA(object, contrasts, "~ 0 + group", 0.05, 0, TRUE, TRUE, ruv_factors = ruv_factors)
        } else runMetabolitesDA(object, contrasts, "~ 0 + group", 0.05, 0, TRUE, TRUE, ruv_factors = ruv_factors))
        long <- as.data.frame(value[[if (identical(omic, "lipid")) "da_lipids_long" else "da_metabolites_long"]])
        id_column <- if (identical(omic, "lipid")) "lipid_id" else "metabolite_id"
        out <- data.frame(id = as.character(long[[id_column]]), p = long$P.Value, q = long$fdr_qvalue, assay = long$assay,
            stringsAsFactors = FALSE)
        methods <- vapply(split(out$p, out$assay), qMethod, character(1L))
        attr(out, "qmethod") <- if (!length(methods)) NA_character_ else if (length(unique(methods)) == 1L) methods[[1L]] else "mixed"
        attr(out, "failed") <- length(value$failed_assays)
    }
    attr(out, "warnings") <- warning_count$n - before
    out
}
# ==== end verbatim ====

# ==== verbatim: filter-effect.R:237-237 at 6416dbd6 ====
`%||%` <- function(x, y) if (is.null(x)) y else x
# ==== end verbatim ====

# ==== verbatim: filter-effect.R:244-246 at 6416dbd6 ====
MIN_REPS <- 2
MIN_GROUPS <- 2
PERCENTILE <- 1
# ==== end verbatim ====

# ==== verbatim: filter-effect.R:248-249 at 6416dbd6 ====
k187Abort <- function(class, message) stop(structure(class = c(class, "error", "condition"),
    list(message = message, call = NULL)))
# ==== end verbatim ====

# ==== verbatim: filter-effect.R:263-270 at 6416dbd6 ====
filterKept <- function(mat, design) {
    object <- protObject(mat, design)
    env <- new.env()
    env$config_list <- list()
    object <- quiet(updateMissingValueParameters(object, min_reps_per_group = MIN_REPS, min_groups = MIN_GROUPS, env = env))
    filtered <- quiet(removeRowsWithMissingValuesPercent(theObject = object, proteins_intensity_cutoff_percentile = PERCENTILE))
    as.character(filtered@protein_quant_table$Protein.Group)
}
# ==== end verbatim ====

# ==== verbatim: filter-effect.R:276-286 at 6416dbd6 ====
perGroupKept <- function(masked, group, m = MIN_REPS, g = MIN_GROUPS) {
    threshold <- stats::quantile(masked[is.finite(masked)], probs = PERCENTILE / 100, type = 7, names = FALSE)
    present <- is.finite(masked) & masked > threshold
    passing <- integer(nrow(masked))
    for (level in unique(group)) {
        columns <- group == level
        n <- sum(columns)
        passing <- passing + (rowSums(!present[, columns, drop = FALSE]) <= n - min(n, m))
    }
    rownames(masked)[passing >= g]
}
# ==== end verbatim ====

# ---- KIT-211 (not copied) ------------------------------------------------------------------------------------------
# Cells A to D are KIT-197's (idx 201 to 204). Av (idx 211, reported only) and Ag (idx 212, reported only, VAL-315) have
# cell A's layout (2 x 3); their marked lines are in generateReplicate (Av: the variance of the replicate halves for every
# 4 log2 of level; Ag: the planted changes are 3 times larger, +-3 log2).
CELLS <- list(A = cell(201L, "prot", 3000L, 2L, 3L), B = cell(202L, "prot", 3000L, 4L, 3L),
    C = cell(203L, "prot", 3000L, 5L, 3L), D = cell(204L, "prot", 3000L, 3L, 4L),
    Av = cell(211L, "prot", 3000L, 2L, 3L), Ag = cell(212L, "prot", 3000L, 2L, 3L))
MASKS <- c("vd5", "vd10", "vd25", "mcar5", "mcar10", "mcar25")
CELL_MASKS <- list(A = MASKS, B = MASKS, C = MASKS, D = MASKS, Av = MASKS[1:3], Ag = MASKS[1:3])
# The old rule (ticket 2.4): minimum present values per group m and minimum passing groups g before KIT-187.
OLD_RULE <- list(A = c(m = 3, g = 2), B = c(m = 3, g = 2), C = c(m = 3, g = 2), D = c(m = 2, g = 3),
    Av = c(m = 3, g = 2), Ag = c(m = 3, g = 2))
# Arms (ticket 4): S shipped; A the test also treats values at or below tau as missing; B the filter counts them as
# present; Bp (B') is B plus the share rule; N is no filter, C0 complete data (both references). Ag: B and B' also at one
# group required (B1, Bp1); no N and no C0.
ARMS <- list(main = c("S", "A", "B", "Bp", "N", "C0"), Ag = c("S", "A", "B", "Bp", "B1", "Bp1"))
SHARE <- 0.5
BAND_WIDTH <- 1
SUBSETS <- c("added", "kept_both", "whole", "whole0", "common", "onoff")
SUM_COLUMNS <- c("n", "n_tested", "n_p05", "n_call", "n_bh", "sum_err", "sum_s_err", "sum_sgn_err")
METRIC_COLUMNS <- c("cell", "mask", "r", "arm", "kept", "added", "calls", "t_r", "mono_ok", "qmethod", "threshold", "viol",
    "e10_dev", "e3_pass", "n_planted", "n_down_all", "n1_band_all", "n_onoff_all", "calls_both_s", "calls_new_s", "gained", "lost",
    "seconds")

k211Abort <- function(class, message) stop(structure(class = c(class, "error", "condition"),
    list(message = message, call = NULL)))

num <- function(x) ifelse(is.na(x), "NA", sprintf("%.17g", x))

writeTsv <- function(df, path) {
    for (column in names(df)) if (is.double(df[[column]])) df[[column]] <- num(df[[column]])
    utils::write.table(df, path, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
}

# B (4.3): a value is present when it is finite; counted per group as perGroupKept counts (m, g).
bFilterKept <- function(masked, group, m = MIN_REPS, g = MIN_GROUPS) {
    present <- is.finite(masked)
    passing <- integer(nrow(masked))
    for (level in unique(group)) {
        columns <- group == level
        n <- sum(columns)
        passing <- passing + (rowSums(!present[, columns, drop = FALSE]) <= n - min(n, m))
    }
    rownames(masked)[passing >= g]
}

# B' (4.3a): a protein of the B-kept set is removed when the share of ALL samples of the table holding a finite value
# strictly below the threshold is at least SHARE (lipidIntensityFilteringHelper's rule: kept when the share is below).
shareKept <- function(masked, kept, threshold) {
    x <- masked[kept, , drop = FALSE]
    kept[rowSums(is.finite(x) & x < threshold) / ncol(masked) < SHARE]
}

# E2: the package's filter run with the threshold switched off (the percentile 0 is the table's minimum, a sentinel
# protein whose every value is -1e6), the sentinel removed; the kept set is B's (m, g).
packageBKept <- function(masked, design, m, g) {
    x <- rbind(masked, K211_SENTINEL = rep(-1e6, ncol(masked)))
    object <- protObject(x, design)
    env <- new.env()
    env$config_list <- list()
    object <- quiet(updateMissingValueParameters(object, min_reps_per_group = m, min_groups = g, env = env))
    filtered <- quiet(removeRowsWithMissingValuesPercent(theObject = object, proteins_intensity_cutoff_percentile = 0))
    setdiff(as.character(filtered@protein_quant_table$Protein.Group), "K211_SENTINEL")
}

# E9: lipidomics' and metabolomics' share helper on the B-kept table (its threshold and proportion given).
helperShareKept <- function(masked, kept, threshold) {
    table <- data.frame(Protein.Group = kept, masked[kept, , drop = FALSE], check.names = FALSE, stringsAsFactors = FALSE)
    filtered <- quiet(lipidIntensityFilteringHelper(table, threshold, SHARE, "Protein.Group",
        sample_columns = colnames(masked)))
    as.character(filtered$Protein.Group)
}

# The largest absolute difference of logFC between two DA tables over the given proteins (an identity, E10); 1e300
# when the two tables disagree on which of the proteins have a value.
lfcDeviation <- function(da_x, da_y, ids) {
    a <- da_x$logFC[match(ids, da_x$id)]
    b <- da_y$logFC[match(ids, da_y$id)]
    if (!identical(is.na(a), is.na(b))) return(1e300)
    if (!length(a) || all(is.na(a))) return(0)
    max(abs(a - b), na.rm = TRUE)
}

# Per protein of H (ticket 6.2 step 4); `threshold` is the filter's: the 1st percentile (type 7) of every finite value
# of the masked matrix (perGroupKept; R/func_prot_qc_filtering_helpers.R:229-233). `rm_all` is the complete-data row
# mean (R-h's band: above the threshold and at most BAND_WIDTH log2 above it); ell0 is the complete-data counterpart of
# ell (a complete value at or below the threshold in control or treated); `onoff` marks Ag's on/off planted proteins.
proteinTable <- function(sim, masked, kept_h, kept_o, kept_s, da, threshold, rm_all, onoff) {
    ids <- rownames(masked)
    h_ids <- ids[ids %in% kept_h]
    x <- masked[h_ids, , drop = FALSE]
    control <- sim$group == "control"
    treated <- sim$group == "treated"
    other <- !(control | treated)
    na <- is.na(x)
    a_c <- rowSums(na[, control, drop = FALSE])
    a_t <- rowSums(na[, treated, drop = FALSE])
    a_o <- if (any(other)) rowSums(na[, other, drop = FALSE]) else numeric(nrow(x))
    low <- is.finite(x) & x <= threshold
    ell <- as.integer(rowSums(low[, control | treated, drop = FALSE]) > 0L)
    low0 <- sim$y[h_ids, , drop = FALSE] <= threshold
    ell0 <- as.integer(rowSums(low0[, control | treated, drop = FALSE]) > 0L)
    band <- as.integer(rm_all[h_ids] > threshold & rm_all[h_ids] <= threshold + BAND_WIDTH)
    s <- ifelse(a_t > 0 & a_c == 0, 1, ifelse(a_c > 0 & a_t == 0, -1, 0))
    stratum <- ifelse(s == 1, "one_treated", ifelse(s == -1, "one_control",
        ifelse(a_c > 0 & a_t > 0, "both", ifelse(a_o > 0, "other_only", "none"))))
    m <- match(h_ids, da$id)
    p <- da$p[m]
    q <- da$q[m]
    lfc <- da$logFC[m]
    beta <- unname(sim$beta[h_ids])
    tested <- is.finite(p) & is.finite(lfc)
    called <- tested & !is.na(q) & q < 0.05
    bh <- rep(NA_real_, length(p))
    bh[tested] <- stats::p.adjust(p[tested], "BH")
    err <- ifelse(tested, lfc - beta, 0)
    data.frame(id = h_ids, in_o = h_ids %in% kept_o, in_s = h_ids %in% kept_s, class = ifelse(beta == 0, "null",
        ifelse(beta > 0, "up", "down")), stratum = stratum, ell = ell, ell0 = ell0, band = band,
        onoff = unname(onoff[h_ids]), n = 1, n_tested = as.numeric(tested), n_p05 = as.numeric(tested & p < 0.05),
        n_call = as.numeric(called), n_bh = as.numeric(tested & !is.na(bh) & bh < 0.05), sum_err = err,
        sum_s_err = s * err, sum_sgn_err = sign(beta) * err, p = p, called = called, stringsAsFactors = FALSE)
}

# Sums of one subset of a protein table by the given key columns; the key columns the subset does not use read "." or -1.
aggregateSubset <- function(tab, members, subset, keys) {
    t <- tab[members, , drop = FALSE]
    if (!nrow(t)) return(NULL)
    key <- do.call(paste, c(lapply(keys, function(k) t[[k]]), sep = "|"))
    sums <- rowsum(as.matrix(t[SUM_COLUMNS]), key, reorder = TRUE)
    parts <- do.call(rbind, strsplit(rownames(sums), "|", fixed = TRUE))
    out <- data.frame(subset = rep(subset, nrow(sums)), class = ".", stratum = ".", ell = -1L, band = -1L, ell0 = -1L,
        stringsAsFactors = FALSE)
    for (i in seq_along(keys)) out[[keys[[i]]]] <- if (keys[[i]] %in% c("ell", "band", "ell0")) as.integer(parts[, i]) else
        parts[, i]
    data.frame(out, sums, row.names = NULL, stringsAsFactors = FALSE)
}

armRows <- function(tab, arm, onoff_rows, reference_only) {
    parts <- if (reference_only) list(aggregateSubset(tab, rep(TRUE, nrow(tab)), "whole", c("class", "stratum", "ell",
        "band")), aggregateSubset(tab, rep(TRUE, nrow(tab)), "whole0", c("class", "ell0"))) else
        list(aggregateSubset(tab, !tab$in_o, "added", c("class", "stratum", "ell")),
            aggregateSubset(tab, tab$in_o, "kept_both", c("class", "stratum", "ell")),
            aggregateSubset(tab, rep(TRUE, nrow(tab)), "whole", c("class", "stratum", "ell", "band")),
            aggregateSubset(tab, rep(TRUE, nrow(tab)), "whole0", c("class", "ell0")),
            aggregateSubset(tab, !tab$in_o & tab$in_s, "common", c("class", "ell")),
            if (onoff_rows) aggregateSubset(tab, tab$onoff, "onoff", "class"))
    out <- do.call(rbind, parts)
    cbind(data.frame(arm = rep(arm, nrow(out)), stringsAsFactors = FALSE), out)
}

runTask <- function(cellname, r, lifecycle) {
    sim <- generateReplicate(cellname, r, strength = 0)
    ids <- rownames(sim$y)
    ag <- identical(cellname, "Ag")
    arms <- if (ag) ARMS$Ag else ARMS$main
    if (length(ARMS_FILTER)) arms <- intersect(arms, ARMS_FILTER)
    rule <- OLD_RULE[[cellname]]
    masks <- CELL_MASKS[[cellname]]
    control <- sim$group == "control"
    treated <- sim$group == "treated"
    rm_all <- rowMeans(sim$y)
    planted <- stats::setNames(sim$beta != 0, ids)
    n_planted <- sum(planted)
    c0_da <- NULL
    c0_seconds <- 0
    if ("C0" %in% arms) {
        started <- proc.time()[[3L]]
        c0_da <- daTable(sim, mat = sim$y)
        c0_seconds <- proc.time()[[3L]] - started
    }
    rows <- list()
    metrics <- list()
    checks <- list()
    for (mask in masks) {
        started <- proc.time()[[3L]]
        masked <- sim$y
        masked[maskMatrix(sim, mask)] <- NA_real_
        threshold <- stats::quantile(masked[is.finite(masked)], probs = PERCENTILE / 100, type = 7, names = FALSE)
        kept <- list(S = filterKept(masked, sim$design))
        kept$A <- kept$S
        kept$B <- bFilterKept(masked, sim$group, MIN_REPS, MIN_GROUPS)
        kept$Bp <- shareKept(masked, kept$B, threshold)
        kept$N <- ids
        if (ag) {
            kept$B1 <- bFilterKept(masked, sim$group, MIN_REPS, 1)
            kept$Bp1 <- shareKept(masked, kept$B1, threshold)
        }
        kept_o <- perGroupKept(masked, sim$group, rule[["m"]], rule[["g"]])
        masked_a <- masked
        masked_a[is.finite(masked) & masked <= threshold] <- NA_real_
        da_arms <- setdiff(arms, "C0")
        input <- lapply(stats::setNames(da_arms, da_arms), function(arm) if (identical(arm, "A"))
            masked_a[kept$A, , drop = FALSE] else masked[kept[[arm]], , drop = FALSE])
        das <- lapply(input, function(mat) daTable(sim, mat = mat))
        onoff <- stats::setNames(rep(FALSE, length(ids)), ids)
        if (ag) {
            lod <- solveLod(sim$y, maskRate(mask))
            onoff <- planted & (rowSums(sim$y[, control, drop = FALSE] < lod) == sum(control) |
                rowSums(sim$y[, treated, drop = FALSE] < lod) == sum(treated))
        }
        n1_band_all <- sum(planted & rm_all > threshold & rm_all <= threshold + BAND_WIDTH)
        tabs <- lapply(stats::setNames(da_arms, da_arms), function(arm) proteinTable(sim, masked, kept[[arm]], kept_o,
            kept$S, das[[arm]], threshold, rm_all, onoff))
        if ("C0" %in% arms) tabs$C0 <- proteinTable(sim, sim$y, ids, character(), kept$S, c0_da, threshold, rm_all, onoff)
        calls_s <- tabs$S$id[tabs$S$called]
        violations <- c(S = length(setdiff(kept_o, kept$S)), A = as.integer(!identical(rownames(input$A), kept$S)),
            B = length(setdiff(kept$S, kept$B)), Bp = length(setdiff(kept$Bp, kept$B)), N = 0L, C0 = 0L,
            B1 = length(setdiff(kept$B, kept$B1)), Bp1 = length(setdiff(kept$Bp, kept$Bp1)))
        deviation <- c(S = 0, A = 0, C0 = 0,
            B = if ("B" %in% arms) lfcDeviation(das$B, das$S, kept$S) else 0,
            Bp = if ("Bp" %in% arms) lfcDeviation(das$Bp, das$S, intersect(kept$S, kept$Bp)) else 0,
            N = if ("N" %in% arms) lfcDeviation(das$N, das$B, kept$B) else 0,
            B1 = if ("B1" %in% arms) lfcDeviation(das$B1, das$B, kept$B) else 0,
            Bp1 = if ("Bp1" %in% arms) lfcDeviation(das$Bp1, das$Bp, kept$Bp) else 0)
        expected_na <- is.na(masked[kept$S, , drop = FALSE]) | (is.finite(masked[kept$S, , drop = FALSE]) &
            masked[kept$S, , drop = FALSE] <= threshold)
        e3 <- if ("A" %in% arms) identical(is.na(input$A), expected_na) && identical(input$A[!is.na(input$A)],
            masked[kept$S, , drop = FALSE][!is.na(input$A)]) else TRUE
        arm_rows <- list()
        arm_metrics <- list()
        for (arm in arms) {
            tab <- tabs[[arm]]
            calls_x <- tab$id[tab$called]
            t_r <- if (any(tab$called)) max(tab$p[tab$called]) else 0
            mono_ok <- sum(tab$n_tested == 1 & tab$p <= t_r) == sum(tab$called)
            arm_rows[[arm]] <- armRows(tab, arm, onoff_rows = ag, reference_only = identical(arm, "C0"))
            arm_metrics[[arm]] <- data.frame(cell = cellname, mask = mask, r = r, arm = arm, kept = nrow(tab),
                added = sum(!tab$in_o), calls = length(calls_x), t_r = t_r, mono_ok = mono_ok,
                qmethod = attr(if (identical(arm, "C0")) c0_da else das[[arm]], "qmethod"), threshold = threshold,
                viol = violations[[arm]], e10_dev = deviation[[arm]], e3_pass = if (identical(arm, "A")) e3 else TRUE,
                n_planted = n_planted, n_down_all = sum(sim$beta < 0), n1_band_all = n1_band_all, n_onoff_all = sum(onoff),
                calls_both_s = length(intersect(calls_x, calls_s)), calls_new_s = length(setdiff(calls_x, calls_s)),
                gained = length(setdiff(tab$id, tabs$S$id)), lost = length(setdiff(tabs$S$id, tab$id)),
                seconds = NA_real_, stringsAsFactors = FALSE)
        }
        if (lifecycle) {
            e9_before <- warning_count$n
            package_b <- packageBKept(masked, sim$design, MIN_REPS, MIN_GROUPS)
            e2 <- setequal(kept$B, package_b) && length(kept$B) == length(package_b)
            package_b1 <- if (ag) packageBKept(masked, sim$design, MIN_REPS, 1)
            e2g1 <- if (ag) setequal(kept$B1, package_b1) && length(kept$B1) == length(package_b1) else TRUE
            helper <- helperShareKept(masked, kept$B, threshold)
            e9 <- setequal(helper, kept$Bp) && length(helper) == length(kept$Bp)
            checks[[length(checks) + 1L]] <- data.frame(cell = cellname, mask = mask, r = r, e2_pass = e2,
                e2g1_pass = e2g1, e3_pass = e3, e7_viol = sum(violations[arms]), e9_pass = e9,
                e9_warnings = warning_count$n - e9_before, e10_dev = max(deviation[arms]), stringsAsFactors = FALSE)
        }
        seconds <- proc.time()[[3L]] - started + c0_seconds / length(masks)
        arm_metrics$S$seconds <- seconds
        rows[[length(rows) + 1L]] <- cbind(data.frame(cell = cellname, mask = mask, r = r, stringsAsFactors = FALSE),
            do.call(rbind, arm_rows))
        metrics[[length(metrics) + 1L]] <- do.call(rbind, arm_metrics)[METRIC_COLUMNS]
    }
    list(rows = do.call(rbind, rows), metrics = do.call(rbind, metrics), checks = if (lifecycle) do.call(rbind, checks))
}

partPaths <- function(out, cellname, r) file.path(out, "parts", sprintf("%s-%04d.%s.tsv", cellname, r, c("rows", "metrics")))

writeTask <- function(out, cellname, r, value) {
    paths <- partPaths(out, cellname, r)
    writeTsv(value$rows, paste0(paths[[1L]], ".tmp"))
    writeTsv(value$metrics, paste0(paths[[2L]], ".tmp"))
    file.rename(paste0(paths[[1L]], ".tmp"), paths[[1L]])
    file.rename(paste0(paths[[2L]], ".tmp"), paths[[2L]])
    if (!is.null(value$checks)) writeTsv(value$checks, file.path(out, "parts", sprintf("%s-%04d.checks.tsv", cellname, r)))
}

readPart <- function(path) utils::read.delim(path, colClasses = "character", check.names = FALSE)

# The parts of the tasks given (a data frame of cell and r) into <cell>.tsv.gz, metrics.tsv and range.tsv.
assemble <- function(out, tasks, expected) {
    tasks <- tasks[order(match(tasks$cell, names(CELLS)), tasks$r), , drop = FALSE]
    cells <- unique(tasks$cell)
    for (cellname in cells) {
        rows <- do.call(rbind, lapply(tasks$r[tasks$cell == cellname], \(r) readPart(partPaths(out, cellname, r)[[1L]])))
        connection <- gzfile(file.path(out, paste0(cellname, ".tsv.gz")), "w")
        utils::write.table(rows, connection, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
        close(connection)
    }
    metrics <- do.call(rbind, lapply(seq_len(nrow(tasks)), \(i) readPart(partPaths(out, tasks$cell[[i]], tasks$r[[i]])[[2L]])))
    utils::write.table(metrics, file.path(out, "metrics.tsv"), sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
    writeTsv(data.frame(cell = cells, first = vapply(cells, \(cellname) min(tasks$r[tasks$cell == cellname]), numeric(1L)),
        last = vapply(cells, \(cellname) max(tasks$r[tasks$cell == cellname]), numeric(1L)),
        expected = expected, complete = vapply(cells, \(cellname) sum(tasks$cell == cellname), numeric(1L)),
        stringsAsFactors = FALSE), file.path(out, "range.tsv"))
    invisible(file.remove(unlist(lapply(seq_len(nrow(tasks)), \(i) partPaths(out, tasks$cell[[i]], tasks$r[[i]])))))
}

sha256 <- function(path) digest::digest(file = path, algo = "sha256")

# ---- the controller ------------------------------------------------------------------------------------------------
args <- commandArgs(TRUE)
if (!length(args) %in% 3:4) k211Abort("k211_usage", "usage: sim-k211.R <out dir> <first r> <last r> [<replicate list>]")
started_at <- Sys.time()
lifecycle <- identical(Sys.getenv("K211_LIFECYCLE"), "1")
determinism <- identical(Sys.getenv("K211_DETERMINISM"), "1")
partial <- identical(Sys.getenv("K211_PARTIAL"), "1")
workers <- as.integer(Sys.getenv("K211_WORKERS", "4"))
ARMS_FILTER <- Filter(nzchar, strsplit(Sys.getenv("K211_ARMS"), ",", fixed = TRUE)[[1L]])
if (length(ARMS_FILTER)) lifecycle <- FALSE
first <- as.integer(args[[2L]])
last <- as.integer(args[[3L]])
if (is.na(first) || is.na(last) || first < 1L || last < first) k211Abort("k211_usage", "bad replicate range")
replicates <- seq.int(first, last)
out <- args[[1L]]
dir.create(file.path(out, "parts"), recursive = TRUE, showWarnings = FALSE)
out <- normalizePath(out)
cells <- strsplit(Sys.getenv("K211_CELLS", "A,B,C,D"), ",", fixed = TRUE)[[1L]]
if (!length(cells) || !all(cells %in% names(CELLS))) k211Abort("k211_usage", "K211_CELLS names an unknown cell")
cells <- names(CELLS)[names(CELLS) %in% cells]
tasks <- expand.grid(r = replicates, cell = cells, stringsAsFactors = FALSE)[c("cell", "r")]
if (length(args) == 4L) {
    listed <- utils::read.table(args[[4L]], col.names = c("cell", "r"), colClasses = c("character", "integer"))
    if (!all(listed$cell %in% cells) || !all(listed$r %in% replicates))
        k211Abort("k211_usage", "the replicate list names a task outside the range")
    tasks <- listed
}
script_file <- normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[[1L]]))
bound_commit <- system2("git", c("rev-parse", "HEAD"), stdout = TRUE)
cat("KIT-211 sim: commit", bound_commit, "range", first, last, "tasks", nrow(tasks), "lifecycle", lifecycle,
    "determinism", determinism, "workers", workers, "\n")
suppressPackageStartupMessages(pkgload::load_all(".", export_all = TRUE, helpers = FALSE, quiet = TRUE))
logger::log_threshold(logger::FATAL)

runTasks <- function(tasks, out) parallel::mclapply(seq_len(nrow(tasks)), \(i) tryCatch({
    started <- proc.time()[[3L]]
    value <- runTask(tasks$cell[[i]], tasks$r[[i]], lifecycle)
    writeTask(out, tasks$cell[[i]], tasks$r[[i]], value)
    list(cell = tasks$cell[[i]], r = tasks$r[[i]], seconds = proc.time()[[3L]] - started)
}, error = function(e) list(cell = tasks$cell[[i]], r = tasks$r[[i]], error = conditionMessage(e))),
mc.cores = workers, mc.preschedule = FALSE)

results <- runTasks(tasks, out)
errors <- Filter(\(x) is.null(x) || inherits(x, "try-error") || !is.null(x$error), results)
if (length(errors)) {
    first_error <- if (is.null(errors[[1L]])) "a worker returned nothing (killed)" else if (inherits(errors[[1L]],
        "try-error")) as.character(errors[[1L]]) else errors[[1L]]$error
    writeLines(sprintf("%d task(s) failed; first: %s", length(errors), first_error), file.path(out, "errors.txt"))
    k211Abort("k211_instrument", sprintf("%d task(s) failed; first: %s", length(errors), first_error))
}
all_tasks <- expand.grid(r = replicates, cell = cells, stringsAsFactors = FALSE)[c("cell", "r")]
done <- vapply(seq_len(nrow(all_tasks)), \(i) file.exists(partPaths(out, all_tasks$cell[[i]], all_tasks$r[[i]])[[2L]]),
    logical(1L))
if (!all(done) && !partial) {
    cat("KIT-211 sim: tasks done, range not yet complete (a continuation is expected)\n")
    quit(status = 0L)
}

# Lifecycle checks (E2, E3, E7, E9, E10) gathered before the parts are assembled.
if (lifecycle) {
    check_files <- list.files(file.path(out, "parts"), "\\.checks\\.tsv$", full.names = TRUE)
    checks <- do.call(rbind, lapply(check_files, readPart))
    units <- nrow(checks)
    summary <- data.frame(check = c("E2", "E3", "E7", "E9", "E10"), units = units,
        passed = c(sum(checks$e2_pass == "TRUE" & checks$e2g1_pass == "TRUE"), sum(checks$e3_pass == "TRUE"),
            sum(as.numeric(checks$e7_viol) == 0), sum(checks$e9_pass == "TRUE"), sum(as.numeric(checks$e10_dev) == 0)),
        max_abs_deviation = c(NA, NA, NA, NA, max(as.numeric(checks$e10_dev))),
        detail = c("B's kept set equals the package filter with the threshold switched off (Ag also at 1 group)",
            "A's DA input holds NA exactly at the masked cells and the finite cells at or below the threshold",
            "F1's identities (O in H_S, H_A = H_S, H_S in H_B, H_Bp in H_B, Ag: H_B in H_B1, H_Bp in H_Bp1)",
            "B' equals lipidIntensityFilteringHelper on the B-kept table",
            "logFC identical across arms on identical inputs"), stringsAsFactors = FALSE)
    summary$pass <- summary$passed == summary$units
    writeTsv(summary, file.path(out, "checks-sim.tsv"))
    writeTsv(data.frame(e9_warnings_counted = sum(as.numeric(checks$e9_warnings))), file.path(out, "e9-warnings.tsv"))
    invisible(file.remove(check_files))
}

assemble(out, all_tasks[done, , drop = FALSE], expected = length(replicates))

# Determinism: cell A, the first five replicates of the range, run again and compared with the assembled rows.
if (determinism) {
    det_out <- file.path(out, "determinism")
    dir.create(file.path(det_out, "parts"), recursive = TRUE, showWarnings = FALSE)
    det_r <- replicates[seq_len(min(5L, length(replicates)))]
    again <- runTasks(data.frame(cell = "A", r = det_r, stringsAsFactors = FALSE), det_out)
    if (length(Filter(\(x) is.null(x) || inherits(x, "try-error") || !is.null(x$error), again))) k211Abort("k211_instrument", "the determinism pass failed")
    main_rows <- utils::read.delim(file.path(out, "A.tsv.gz"), colClasses = "character", check.names = FALSE)
    main_metrics <- readPart(file.path(out, "metrics.tsv"))
    det <- do.call(rbind, lapply(det_r, \(r) {
        paths <- partPaths(det_out, "A", r)
        rows <- readPart(paths[[1L]])
        metrics <- readPart(paths[[2L]])
        mine <- main_rows[main_rows$r == as.character(r), , drop = FALSE]
        theirs <- main_metrics[main_metrics$cell == "A" & main_metrics$r == as.character(r), , drop = FALSE]
        rownames(mine) <- rownames(theirs) <- rownames(rows) <- rownames(metrics) <- NULL
        keep <- setdiff(METRIC_COLUMNS, "seconds")
        data.frame(cell = "A", r = r, rows_identical = identical(rows, mine),
            metrics_identical = identical(metrics[keep], theirs[keep]), stringsAsFactors = FALSE)
    }))
    writeTsv(det, file.path(out, "determinism.tsv"))
    invisible(file.remove(unlist(lapply(det_r, \(r) partPaths(det_out, "A", r)))))
}

tables <- list.files(out, "\\.(tsv|gz)$", full.names = TRUE)
metrics_all <- readPart(file.path(out, "metrics.tsv"))
seconds_s <- as.numeric(metrics_all$seconds[metrics_all$arm == "S"])
provenance <- list(ticket = "KIT-211", script = "docs/qa/workflow-execution/kit/kit-211/sim-k211.R",
    script_sha256 = sha256(script_file), bound_commit = bound_commit, lifecycle = lifecycle, determinism = determinism,
    cells = cells, cell_index = vapply(CELLS[cells], `[[`, integer(1L), "idx"), masks = CELL_MASKS[cells],
    arms = ARMS, arms_filter = ARMS_FILTER, share = SHARE, band_width = BAND_WIDTH,
    old_rule = lapply(OLD_RULE, as.list), replicates = c(first = first, last = last), tasks_this_process = nrow(tasks),
    tasks_complete = sum(done), seed_rule = "set.seed(185000000 + 10000 * idx + r, Mersenne-Twister, Inversion, Rejection)",
    workers = workers, seconds_per_cell = as.list(tapply(seconds_s, metrics_all$cell[metrics_all$arm == "S"], sum)),
    seconds_per_replicate_mask = mean(seconds_s), warnings = warning_count$n,
    r_version = R.version.string, limma = as.character(utils::packageVersion("limma")),
    qvalue = as.character(utils::packageVersion("qvalue")), blas = extSoftVersion()[["BLAS"]], lapack = La_library(),
    runner_image = paste(Sys.getenv("ImageOS"), Sys.getenv("ImageVersion")),
    processor = tryCatch(sub("^model name\\s*:\\s*", "", grep("^model name", readLines("/proc/cpuinfo"), value = TRUE)[[1L]]),
        error = function(e) NA),
    started = format(started_at, "%Y-%m-%dT%H:%M:%S%z"), finished = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
    exit_status = 0L, tables_sha256 = as.list(stats::setNames(vapply(tables, sha256, character(1L)), basename(tables))))
jsonlite::write_json(provenance, file.path(out, "provenance-sim.json"), auto_unbox = TRUE, pretty = TRUE, digits = NA)
cat(sprintf("KIT-211 sim written: %d tasks, wall %.1f s\n", nrow(tasks), as.numeric(difftime(Sys.time(), started_at,
    units = "secs"))))
