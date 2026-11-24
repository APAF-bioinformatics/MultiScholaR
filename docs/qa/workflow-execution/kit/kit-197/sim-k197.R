# KIT-197 (VAL-111, VAL-121): why the proteins only the count rule keeps carry a higher false share among their calls.
# The simulation (ticket sections 5.1, 5.2, 7.1). Evidence code under docs/qa: no package code.
#
# Usage, from the package root of the bound commit on a run-tests runner (full git history):
#   Rscript docs/qa/workflow-execution/kit/kit-197/sim-k197.R <out dir> <first r> <last r> [<replicate list>]
#   <replicate list>  optional file of "<cell> <r>" lines: only those tasks run (a continuation job); the
#                     part files already in <out dir>/parts are kept.
# Environment: K197_LIFECYCLE=1 (checks E1, E2, E3 on every task; E1 needs K197_BASE=<40-digit commit>, S2's base),
#   K197_DETERMINISM=1 (after the main run, cell A replicates <first r> to <first r> + 4 are run again and compared),
#   K197_WORKERS (forked workers; default 4), K197_CELLS (comma list of cells; default A,B,C,D; check E6 uses A).
# One tree: pkgload::load_all(".", export_all = TRUE). One task per (cell, r) running all six masks; each task writes
# parts/<cell>-<r>.rows.tsv and parts/<cell>-<r>.metrics.tsv through a `.tmp` name, the metrics file last, so a task
# whose metrics file exists is complete. When every task of the range is complete the parts are assembled into
# <cell>.tsv.gz (the sufficient statistics of 5.2 (6)) and s2metrics.tsv (per cell, mask and r: S2's metrics, t_r,
# the threshold, |O \ H|), and the part files are removed.
#
# Copied verbatim from docs/qa/workflow-execution/kit/kit-187/filter-effect.R (SHA-256
# 1147a5aa480861581210f1a088c06731e046b20e15d8cdcbe6775a5f58906c9c, committed at 4a58a3be, unchanged at 6416dbd6),
# each block between a "verbatim" marker and an "end verbatim" marker: lines 23-234 (KIT-185's generator, masks,
# object and DA call; one line changed, below), 237 (`%||%`), 244-246 (MIN_REPS, MIN_GROUPS, PERCENTILE), 248-249
# (k187Abort), 263-270 (filterKept), 276-286 (perGroupKept).
# Changed, and only these: filter-effect.R:215 in daTable also keeps `logFC = tab$logFC, t = tab$t,
# ave = tab$AveExpr`; CELLS holds A to D only (S2's indices 201 to 204); MASKS adds mcar5, mcar10, mcar25.
# Seeds: the verbatim generator's set.seed(185000000 + 10000 * idx + r), idx = 200 + i: 187000000 + 10000 i + r.

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
    # effects: exactly 20 % of the features, random positions, half up, in group 2
    n_effects <- round(0.2 * n)
    where <- sample.int(n, n_effects)
    signs <- sample(rep(c(1, -1), c(round(n_effects / 2), n_effects - round(n_effects / 2))))
    beta <- numeric(n)
    beta[where] <- signs
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

# ---- KIT-197 (not copied) ------------------------------------------------------------------------------------------
CELLS <- list(A = cell(201L, "prot", 3000L, 2L, 3L), B = cell(202L, "prot", 3000L, 4L, 3L),
    C = cell(203L, "prot", 3000L, 5L, 3L), D = cell(204L, "prot", 3000L, 3L, 4L))
MASKS <- c("vd5", "vd10", "vd25", "mcar5", "mcar10", "mcar25")
# The old rule (ticket 2.4): minimum present values per group m and minimum passing groups g before KIT-187.
OLD_RULE <- list(A = c(m = 3, g = 2), B = c(m = 3, g = 2), C = c(m = 3, g = 2), D = c(m = 2, g = 3))
S2_FILE <- "docs/qa/workflow-execution/kit/kit-187/filter-effect.R"
S2_SHA256 <- "1147a5aa480861581210f1a088c06731e046b20e15d8cdcbe6775a5f58906c9c"
STRATA <- c("one_treated", "one_control", "both", "other_only", "none")
SUBSETS <- c("added", "kept_both", "ref", "whole")
SUM_COLUMNS <- c("n", "n_tested", "n_p05", "n_p01", "n_call", "n_call_away", "sum_err", "sum_s_err", "sum_sgn_err",
    "sum_sq_err", "sum_resvar", "n_resvar", "n_s")
METRIC_COLUMNS <- c("cell", "mask", "r", "kept_h", "kept_o", "added", "calls", "power", "fdp", "added_share_called",
    "added_fdp", "qmethod", "f1_o_not_h", "t_r", "mono_ok", "threshold", "seconds")

k197Abort <- function(class, message) stop(structure(class = c(class, "error", "condition"),
    list(message = message, call = NULL)))

num <- function(x) ifelse(is.na(x), "NA", sprintf("%.17g", x))

writeTsv <- function(df, path) {
    for (column in names(df)) if (is.double(df[[column]])) df[[column]] <- num(df[[column]])
    utils::write.table(df, path, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
}

# Pooled within-group variance of the observed values (the residual variance of `~ 0 + group` on them).
residualVariance <- function(x, group) {
    ss <- numeric(nrow(x))
    df <- numeric(nrow(x))
    for (level in unique(group)) {
        v <- x[, group == level, drop = FALSE]
        k <- rowSums(!is.na(v))
        centre <- rowMeans(v, na.rm = TRUE)
        ss <- ss + ifelse(k > 0L, rowSums((v - centre)^2, na.rm = TRUE), 0)
        df <- df + pmax(k - 1L, 0L)
    }
    ifelse(df > 0, ss / df, NA_real_)
}

# Per protein of H (ticket 5.2 (4)); `threshold` is the filter's: the 1st percentile (type 7) of every finite value
# of the masked matrix (perGroupKept; R/func_prot_qc_filtering_helpers.R:229-233).
proteinTable <- function(sim, masked, kept_h, kept_o, da, threshold) {
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
    err <- ifelse(tested, lfc - beta, 0)
    resvar <- residualVariance(x, sim$group)
    has_rv <- tested & is.finite(resvar)
    data.frame(id = h_ids, in_o = h_ids %in% kept_o, ref_ok = a_c == 0 & a_t == 0,
        class = ifelse(beta == 0, "null", ifelse(beta > 0, "up", "down")), stratum = stratum, ell = ell,
        n = 1, n_tested = as.numeric(tested), n_p05 = as.numeric(tested & p < 0.05),
        n_p01 = as.numeric(tested & p < 0.01), n_call = as.numeric(called),
        n_call_away = as.numeric(called & s * lfc > 0), sum_err = err, sum_s_err = s * err,
        sum_sgn_err = sign(beta) * err, sum_sq_err = err^2, sum_resvar = ifelse(has_rv, resvar, 0),
        n_resvar = as.numeric(has_rv), n_s = as.numeric(tested & s != 0), p = p, called = called, lfc = lfc,
        stringsAsFactors = FALSE)
}

aggregateSubset <- function(tab, members, subset) {
    t <- tab[members, , drop = FALSE]
    if (!nrow(t)) return(NULL)
    key <- paste(t$class, t$stratum, t$ell, sep = "|")
    sums <- rowsum(as.matrix(t[SUM_COLUMNS]), key, reorder = TRUE)
    parts <- do.call(rbind, strsplit(rownames(sums), "|", fixed = TRUE))
    data.frame(subset = subset, class = parts[, 1L], stratum = parts[, 2L], ell = as.integer(parts[, 3L]),
        sums, row.names = NULL, stringsAsFactors = FALSE)
}

# S2's metrics for one replicate-mask, computed as filter-effect.R:449-463 computes them (O standing for S2's base).
s2Metrics <- function(ids, planted, kept_h, kept_o, da) {
    h <- match(kept_h, ids)
    o <- match(kept_o, ids)
    calls <- match(da$id[!is.na(da$q) & da$q < 0.05], ids)
    added <- setdiff(h, o)
    called_added <- intersect(added, calls)
    list(kept_h = length(h), kept_o = length(o), added = length(added), calls = length(calls),
        power = length(intersect(calls, planted)) / length(planted),
        fdp = if (length(calls)) length(setdiff(calls, planted)) / length(calls) else 0,
        added_share_called = if (length(added)) length(called_added) / length(added) else NA_real_,
        added_fdp = if (length(added)) {
            if (length(called_added)) length(setdiff(called_added, planted)) / length(called_added) else 0
        } else NA_real_,
        f1_o_not_h = length(setdiff(o, h)))
}

# S2's own daTable (filter-effect.R:204-234), read from the committed file for check E2's second half.
loadS2DaTable <- function() {
    if (!identical(digest::digest(file = S2_FILE, algo = "sha256"), S2_SHA256))
        k197Abort("k197_instrument", "filter-effect.R does not have S2's hash")
    env <- new.env(parent = globalenv())
    eval(parse(text = readLines(S2_FILE)[204:234]), envir = env)
    env$daTable
}

runTask <- function(cellname, r, lifecycle, e1_dir) {
    sim <- generateReplicate(cellname, r, strength = 0)
    ids <- rownames(sim$y)
    planted <- which(sim$beta != 0)
    rule <- OLD_RULE[[cellname]]
    rows <- vector("list", length(MASKS))
    metrics <- vector("list", length(MASKS))
    checks <- vector("list", length(MASKS))
    for (k in seq_along(MASKS)) {
        started <- proc.time()[[3L]]
        mask <- MASKS[[k]]
        masked <- sim$y
        masked[maskMatrix(sim, mask)] <- NA_real_
        kept_h <- filterKept(masked, sim$design)
        kept_o <- perGroupKept(masked, sim$group, rule[["m"]], rule[["g"]])
        da <- daTable(sim, mat = masked[kept_h, , drop = FALSE])
        threshold <- stats::quantile(masked[is.finite(masked)], probs = PERCENTILE / 100, type = 7, names = FALSE)
        tab <- proteinTable(sim, masked, kept_h, kept_o, da, threshold)
        rows[[k]] <- do.call(rbind, list(aggregateSubset(tab, !tab$in_o, "added"),
            aggregateSubset(tab, tab$in_o, "kept_both"), aggregateSubset(tab, tab$in_o & tab$ref_ok, "ref"),
            aggregateSubset(tab, rep(TRUE, nrow(tab)), "whole")))
        rows[[k]] <- cbind(data.frame(cell = cellname, mask = mask, r = r, stringsAsFactors = FALSE), rows[[k]])
        t_r <- if (any(tab$called)) max(tab$p[tab$called]) else 0
        mono_ok <- sum(tab$n_tested == 1 & tab$p <= t_r) == sum(tab$called)
        metric <- s2Metrics(ids, planted, kept_h, kept_o, da)
        seconds <- proc.time()[[3L]] - started
        metrics[[k]] <- data.frame(cell = cellname, mask = mask, r = r, metric, qmethod = attr(da, "qmethod"),
            t_r = t_r, mono_ok = mono_ok, threshold = threshold, seconds = seconds,
            stringsAsFactors = FALSE)[METRIC_COLUMNS]
        if (lifecycle) checks[[k]] <- lifecycleChecks(sim, masked, kept_h, kept_o, da, tab, metric, cellname, r, mask,
            e1_dir)
    }
    list(rows = do.call(rbind, rows), metrics = do.call(rbind, metrics), checks = do.call(rbind, checks))
}

# E2 (both halves) and E3 on one replicate-mask; E1's inputs saved for the base tree.
lifecycleChecks <- function(sim, masked, kept_h, kept_o, da, tab, metric, cellname, r, mask, e1_dir) {
    x <- masked[tab$id, , drop = FALSE]
    observed <- rowMeans(x[, sim$group == "treated", drop = FALSE], na.rm = TRUE) -
        rowMeans(x[, sim$group == "control", drop = FALSE], na.rm = TRUE)
    finite <- is.finite(tab$lfc)
    e2a <- if (any(finite)) max(abs(tab$lfc[finite] - observed[finite])) else NA_real_
    s2_da <- S2_DATABLE(sim, mat = masked[kept_h, , drop = FALSE])
    e2b <- identical(s2_da$id, da$id) && identical(s2_da$p, da$p) && identical(s2_da$q, da$q)
    saveRDS(list(masked = masked, design = sim$design, kept_o = kept_o),
        file.path(e1_dir, sprintf("%s-%d-%s.rds", cellname, r, mask)))
    data.frame(cell = cellname, mask = mask, r = r, e2a_max_abs = e2a, e2a_finite = sum(finite),
        e2a_pass = isTRUE(e2a <= 1e-10) && sum(finite) == sum(tab$n_tested), e2b_pass = e2b,
        e3_pass = metric$f1_o_not_h == 0L, stringsAsFactors = FALSE)
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

assemble <- function(out, cells, replicates) {
    paths <- lapply(cells, \(cellname) lapply(replicates, \(r) partPaths(out, cellname, r)))
    for (i in seq_along(cells)) {
        rows <- do.call(rbind, lapply(paths[[i]], \(p) readPart(p[[1L]])))
        connection <- gzfile(file.path(out, paste0(cells[[i]], ".tsv.gz")), "w")
        utils::write.table(rows, connection, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
        close(connection)
    }
    metrics <- do.call(rbind, lapply(unlist(paths, recursive = FALSE), \(p) readPart(p[[2L]])))
    utils::write.table(metrics, file.path(out, "s2metrics.tsv"), sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
    invisible(file.remove(unlist(paths)))
}

sha256 <- function(path) digest::digest(file = path, algo = "sha256")

# ---- E1: the base tree's filter on the lifecycle's masked matrices (a separate Rscript) ---------------------------
args <- commandArgs(TRUE)
if (length(args) && identical(args[[1L]], "--base-filter")) {
    suppressPackageStartupMessages(pkgload::load_all(args[[2L]], export_all = TRUE, helpers = FALSE, quiet = TRUE))
    logger::log_threshold(logger::FATAL)
    files <- sort(list.files(args[[3L]], "\\.rds$", full.names = TRUE))
    e1 <- do.call(rbind, lapply(files, \(f) {
        x <- readRDS(f)
        base <- filterKept(x$masked, x$design)
        data.frame(file = basename(f), n_base = length(base), n_script = length(x$kept_o),
            identical = identical(sort(base), sort(x$kept_o)), stringsAsFactors = FALSE)
    }))
    utils::write.table(e1, file.path(args[[3L]], "e1.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
    quit(status = 0L)
}

# ---- the controller ------------------------------------------------------------------------------------------------
if (!length(args) %in% 3:4) k197Abort("k197_usage", "usage: sim-k197.R <out dir> <first r> <last r> [<replicate list>]")
started_at <- Sys.time()
lifecycle <- identical(Sys.getenv("K197_LIFECYCLE"), "1")
determinism <- identical(Sys.getenv("K197_DETERMINISM"), "1")
workers <- as.integer(Sys.getenv("K197_WORKERS", "4"))
first <- as.integer(args[[2L]])
last <- as.integer(args[[3L]])
if (is.na(first) || is.na(last) || first < 1L || last < first) k197Abort("k197_usage", "bad replicate range")
replicates <- seq.int(first, last)
out <- args[[1L]]
dir.create(file.path(out, "parts"), recursive = TRUE, showWarnings = FALSE)
out <- normalizePath(out)
cells <- strsplit(Sys.getenv("K197_CELLS", paste(names(CELLS), collapse = ",")), ",", fixed = TRUE)[[1L]]
if (!length(cells) || !all(cells %in% names(CELLS))) k197Abort("k197_usage", "K197_CELLS names an unknown cell")
cells <- names(CELLS)[names(CELLS) %in% cells]
tasks <- expand.grid(r = replicates, cell = cells, stringsAsFactors = FALSE)[c("cell", "r")]
if (length(args) == 4L) {
    listed <- utils::read.table(args[[4L]], col.names = c("cell", "r"), colClasses = c("character", "integer"))
    if (!all(listed$cell %in% cells) || !all(listed$r %in% replicates))
        k197Abort("k197_usage", "the replicate list names a task outside the range")
    tasks <- listed
}
script_file <- normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[[1L]]))
bound_commit <- system2("git", c("rev-parse", "HEAD"), stdout = TRUE)
cat("KIT-197 sim: commit", bound_commit, "range", first, last, "tasks", nrow(tasks), "lifecycle", lifecycle,
    "determinism", determinism, "workers", workers, "\n")
suppressPackageStartupMessages(pkgload::load_all(".", export_all = TRUE, helpers = FALSE, quiet = TRUE))
logger::log_threshold(logger::FATAL)
e1_dir <- file.path(tempdir(), "k197-e1")
if (lifecycle) {
    dir.create(e1_dir, showWarnings = FALSE)
    S2_DATABLE <- loadS2DaTable()
}

runTasks <- function(tasks, out) parallel::mclapply(seq_len(nrow(tasks)), \(i) tryCatch({
    started <- proc.time()[[3L]]
    value <- runTask(tasks$cell[[i]], tasks$r[[i]], lifecycle, e1_dir)
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
    k197Abort("k197_instrument", sprintf("%d task(s) failed; first: %s", length(errors), first_error))
}
complete <- all(vapply(cells, \(cellname) all(vapply(replicates, \(r)
    file.exists(partPaths(out, cellname, r)[[2L]]), logical(1L))), logical(1L)))
if (!complete) {
    cat("KIT-197 sim: tasks done, range not yet complete (a continuation is expected)\n")
    quit(status = 0L)
}

# Lifecycle checks (E1 to E3) gathered before the parts are assembled.
if (lifecycle) {
    check_files <- list.files(file.path(out, "parts"), "\\.checks\\.tsv$", full.names = TRUE)
    checks <- do.call(rbind, lapply(check_files, readPart))
    base_commit <- Sys.getenv("K197_BASE")
    if (!grepl("^[0-9a-f]{40}$", base_commit)) k197Abort("k197_usage", "K197_BASE must be a 40-digit commit")
    base_root <- file.path(tempdir(), "k197-base")
    invisible(system2("git", c("worktree", "add", "--detach", shQuote(base_root), base_commit), stdout = TRUE,
        stderr = TRUE))
    rscript <- file.path(R.home("bin"), "Rscript")
    status <- system2(rscript, c(shQuote(script_file), "--base-filter", shQuote(base_root), shQuote(e1_dir)),
        stdout = file.path(out, "e1-base.log"), stderr = file.path(out, "e1-base.log"))
    e1 <- if (status == 0L && file.exists(file.path(e1_dir, "e1.tsv"))) readPart(file.path(e1_dir, "e1.tsv")) else NULL
    base_head <- system2("git", c("-C", shQuote(base_root), "rev-parse", "HEAD"), stdout = TRUE)
    summary <- data.frame(check = c("E1", "E2a", "E2b", "E3"),
        units = c(if (is.null(e1)) 0L else nrow(e1), nrow(checks), nrow(checks), nrow(checks)),
        passed = c(if (is.null(e1)) 0L else sum(e1$identical == "TRUE"), sum(checks$e2a_pass == "TRUE"),
            sum(checks$e2b_pass == "TRUE"), sum(checks$e3_pass == "TRUE")),
        detail = c(sprintf("base tree %s, exit %d", base_head, status),
            sprintf("max |logFC - observed mean difference| %s",
                num(max(suppressWarnings(as.numeric(checks$e2a_max_abs)), na.rm = TRUE))),
            "id, p, q of the modified daTable identical to S2's daTable", "O subset of H"),
        stringsAsFactors = FALSE)
    summary$pass <- summary$units == nrow(checks) & summary$passed == summary$units
    writeTsv(summary, file.path(out, "checks-sim.tsv"))
    invisible(file.remove(check_files))
}

assemble(out, cells, replicates)

# Determinism (6): cell A, the first five replicates of the range, run again and compared with the assembled rows.
if (determinism) {
    det_out <- file.path(out, "determinism")
    dir.create(file.path(det_out, "parts"), recursive = TRUE, showWarnings = FALSE)
    det_r <- replicates[seq_len(min(5L, length(replicates)))]
    again <- runTasks(data.frame(cell = "A", r = det_r, stringsAsFactors = FALSE), det_out)
    if (length(Filter(\(x) is.null(x) || inherits(x, "try-error") || !is.null(x$error), again))) k197Abort("k197_instrument", "the determinism pass failed")
    main_rows <- utils::read.delim(file.path(out, "A.tsv.gz"), colClasses = "character", check.names = FALSE)
    main_metrics <- readPart(file.path(out, "s2metrics.tsv"))
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
metrics_all <- readPart(file.path(out, "s2metrics.tsv"))
provenance <- list(ticket = "KIT-197", script = "docs/qa/workflow-execution/kit/kit-197/sim-k197.R",
    script_sha256 = sha256(script_file), s2_script = S2_FILE, s2_script_sha256 = sha256(S2_FILE),
    bound_commit = bound_commit, s2_head = "e6c42d1d8e18f04d34b3b2a850488ef43ab8d9ab",
    s2_base = "604ec4dd49f6eba3c842df09fdd0a057a7c83658", lifecycle = lifecycle, determinism = determinism,
    cells = cells, cell_index = vapply(CELLS[cells], `[[`, integer(1L), "idx"), masks = MASKS,
    old_rule = lapply(OLD_RULE, as.list), replicates = c(first = first, last = last), tasks_this_process = nrow(tasks),
    seed_rule = "set.seed(185000000 + 10000 * idx + r, Mersenne-Twister, Inversion, Rejection); idx = 200 + i, so 187000000 + 10000 i + r",
    workers = workers, seconds_per_cell = as.list(tapply(as.numeric(metrics_all$seconds), metrics_all$cell, sum)),
    seconds_per_replicate_mask = mean(as.numeric(metrics_all$seconds)), warnings = warning_count$n,
    r_version = R.version.string, limma = as.character(utils::packageVersion("limma")),
    qvalue = as.character(utils::packageVersion("qvalue")), blas = extSoftVersion()[["BLAS"]], lapack = La_library(),
    runner_image = paste(Sys.getenv("ImageOS"), Sys.getenv("ImageVersion")),
    processor = tryCatch(sub("^model name\\s*:\\s*", "", grep("^model name", readLines("/proc/cpuinfo"), value = TRUE)[[1L]]),
        error = function(e) NA),
    started = format(started_at, "%Y-%m-%dT%H:%M:%S%z"), finished = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
    exit_status = 0L, tables_sha256 = as.list(stats::setNames(vapply(tables, sha256, character(1L)), basename(tables))))
jsonlite::write_json(provenance, file.path(out, "provenance-sim.json"), auto_unbox = TRUE, pretty = TRUE, digits = NA)
cat(sprintf("KIT-197 sim written: %d tasks, wall %.1f s\n", nrow(tasks), as.numeric(difftime(Sys.time(), started_at,
    units = "secs"))))
