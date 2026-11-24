# KIT-217 measurement M (ticket KIT-217 r1.3 section 4.1; methods ledger VAL-184): the candidate intensity-filter rules
# for lipidomics and metabolomics, measured by the package's own DA on KIT-185's shapes. Evidence code under docs/qa: no
# package code. Runs on a GitHub-hosted runner with MultiScholaR installed from the bound commit (run-tests `scripts`).
#
# Usage (from the package root):
#   Rscript docs/qa/workflow-execution/kit/kit-217/m/sim-k217.R <cells> <first> <last> <out dir>
#   cells        one cell or a comma list of L2, L1, M2, M1 (deciding), Lv, Le, Me (reported)
#   first, last  replicate range (counted 1 to 200, reported cells 1 to 100; lifecycle 401 to 403)
# Environment: K217_LIFECYCLE=1 (blind lifecycle: no estimand is written; only completion per (cell, mask, arm), seconds
#   per replicate, the SHA-256 of each replicate's rows and pass/fail of E1 to E6), K217_FAULT=1 (E4: arm D raises an
#   error in the first cell's first replicate, mask vd10), K217_WORKERS (default 4), K217_JOB (file-name prefix,
#   default "m"; the determinism job uses "det").
#
# The functions under "verbatim" are copied unchanged from KIT-185's generator,
# docs/qa/workflow-execution/kit/kit-185/sim-k185.R at commit 735c2ea5 (the KIT-185 merge; unchanged on v0.9-staging at
# e6791275), at the line ranges named; the file itself cannot be sourced (it stops without its five arguments). The
# three KIT-217 additions are marked "KIT-217 addition" (cells, mask vdg, the Le/Me effect line) plus the Lv variance
# line. KIT-217's cells take idx = 3200 + i, so the copied seed line gives 217000000 + 10000 i + r. strength = 0: no
# unwanted factor (the factor draws are still made, so the stream is KIT-185's layout).

args <- commandArgs(TRUE)
if (length(args) != 4L) stop("usage: Rscript sim-k217.R <cells> <first> <last> <out dir>")
cells_wanted <- strsplit(args[[1L]], ",", fixed = TRUE)[[1L]]
first_rep <- as.integer(args[[2L]])
last_rep <- as.integer(args[[3L]])
out_dir <- args[[4L]]
stopifnot(!is.na(first_rep), !is.na(last_rep), first_rep <= last_rep)
lifecycle <- identical(Sys.getenv("K217_LIFECYCLE"), "1")
fault <- identical(Sys.getenv("K217_FAULT"), "1")
workers <- as.integer(Sys.getenv("K217_WORKERS", "4"))
job_prefix <- Sys.getenv("K217_JOB", "m")
started_at <- Sys.time()
script_file <- normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[[1L]]))

suppressPackageStartupMessages(library(MultiScholaR))
attach(asNamespace("MultiScholaR"), name = "MultiScholaR:namespace", warn.conflicts = FALSE)
logger::log_threshold(logger::FATAL)

# ---- verbatim: sim-k185.R:43-52 at 735c2ea5 ----
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

# ---- verbatim: sim-k185.R:71-76 at 735c2ea5 ----
GROUPS <- c("control", "treated", "g3", "g4", "g5")
cell <- function(idx, omic, features, groups, per, k = 2L, level = 20, rho = 0, dup = FALSE, batches = 0L,
                 batch_runs = 0L, factors = 2L, assays = 1L) {
    list(idx = idx, omic = omic, features = features, groups = groups, per = per, k = k, level = level, rho = rho,
        dup = dup, batches = batches, batch_runs = batch_runs, factors = factors, assays = assays)
}

# ---- verbatim: sim-k185.R:90-91 at 735c2ea5 ----
MED_VAR_PROT <- 4 * 0.25^2 / stats::qchisq(0.5, 4)
MED_VAR_OMIC <- 0.25

# ---- verbatim: sim-k185.R:101-102 at 735c2ea5 ----
maskMech <- function(mask) if (identical(mask, "none")) "none" else sub("[0-9.]+$", "", mask)
maskRate <- function(mask) if (identical(mask, "none")) 0 else as.numeric(sub("^[a-z]+", "", mask)) / 100

# ---- verbatim: sim-k185.R:105-108 at 735c2ea5 ----
scale1 <- function(v) {
    v <- v - mean(v)
    v / stats::sd(v)
}

# ---- copied: sim-k185.R:109-199 at 735c2ea5, with two marked KIT-217 additions (Lv, Le/Me) ----
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
    if (identical(cellname, "Lv")) s2 <- s2 * 2^((20 - mu) / 2)  # KIT-217 addition: cell Lv (KIT-211's Av line)
    # effects: exactly 20 % of the features, random positions, half up, in group 2
    n_effects <- round(0.2 * n)
    where <- sample.int(n, n_effects)
    signs <- sample(rep(c(1, -1), c(round(n_effects / 2), n_effects - round(n_effects / 2))))
    beta <- numeric(n)
    beta[where] <- signs
    if (cellname %in% c("Le", "Me")) beta <- 3 * beta  # KIT-217 addition (iii): cells Le and Me, effects +-3 log2
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

# ---- verbatim: sim-k185.R:202-205 at 735c2ea5 ----
solveLod <- function(values, rate, shift = 0) {
    stats::uniroot(function(l) mean(stats::plogis(l + shift - values)) - rate, c(min(values) - 60, max(values) + 60),
        tol = 1e-10)$root
}

# ---- copied: sim-k185.R:206-234 at 735c2ea5, with one marked KIT-217 addition (mask vdg) ----
maskMatrix <- function(sim, mask) {
    y <- sim$y
    if (identical(mask, "none")) return(matrix(FALSE, nrow(y), ncol(y), dimnames = dimnames(y)))
    rate <- maskRate(mask)
    switch(maskMech(mask),
        vd = sim$U < stats::plogis(solveLod(y, rate) - y),
        mcar = sim$V < rate,
        # KIT-217 addition (ii): the group-linked mask; lod raised by 0.5 log2 for group 1 (control) before solving
        vdg = {
            shift <- matrix(0.5 * (sim$group == "control"), nrow(y), ncol(y), byrow = TRUE)
            lod <- solveLod(y, rate, shift)
            sim$U < stats::plogis(lod + shift - y)
        },
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

# ---- verbatim: sim-k185.R:257-273 at 735c2ea5 ----
ID_COL <- "feature_id"
assayTable <- function(m) {
    table <- tibble::as_tibble(as.data.frame(m, check.names = FALSE))
    table <- tibble::add_column(table, !!ID_COL := rownames(m), .before = 1L)
    tibble::add_column(table, name = rownames(m), .after = 1L)
}
omicSlot <- function(omic) if (identical(omic, "lipid")) "lipid_data" else "metabolite_data"
omicObject <- function(omic, mat, assay, design) {
    tables <- lapply(split(rownames(mat), factor(assay[rownames(mat)], levels = unique(assay))), function(ids) assayTable(mat[ids, , drop = FALSE]))
    if (identical(omic, "lipid")) {
        methods::new("LipidomicsAssayData", lipid_data = tables, lipid_id_column = ID_COL, annotation_id_column = "name",
            design_matrix = design, sample_id = "Sample_ID", group_id = "group", technical_replicate_id = "extract", args = list())
    } else {
        methods::new("MetaboliteAssayData", metabolite_data = tables, metabolite_id_column = ID_COL, annotation_id_column = "name",
            design_matrix = design, sample_id = "Sample_ID", group_id = "group", technical_replicate_id = "extract", args = list())
    }
}

# ---- verbatim: sim-k185.R:303-328 at 735c2ea5 ----
filterArm <- function(sim, masked, use_presence) {
    p <- sim$p
    if (identical(p$omic, "prot")) {
        object <- protObject(masked, sim$design)
        env <- new.env()
        env$config_list <- list()
        object <- quiet(updateMissingValueParameters(object, min_reps_per_group = 2, min_groups = 2, env = env))
        filtered <- quiet(removeRowsWithMissingValuesPercent(theObject = object))
        ids <- as.character(filtered@protein_quant_table$Protein.Group)
        return(list(ids = ids, object = filtered))
    }
    raw <- omicObject(p$omic, 2^masked, sim$assay, sim$design)
    filtered <- quiet(if (identical(p$omic, "lipid")) {
        lipidIntensityFiltering(raw, lipids_intensity_cutoff_percentile = 10, lipids_proportion_of_samples_below_cutoff = 0.5)
    } else {
        metaboliteIntensityFiltering(raw, metabolites_intensity_cutoff_percentile = 10, metabolites_proportion_of_samples_below_cutoff = 0.5)
    })
    ids <- unlist(lapply(methods::slot(filtered, omicSlot(p$omic)), function(t) as.character(t[[ID_COL]])), use.names = FALSE)
    ids <- rownames(masked)[rownames(masked) %in% ids]
    if (use_presence) {
        seen <- !is.na(masked[ids, , drop = FALSE])
        ok <- Reduce(`&`, lapply(unique(sim$group), function(g) rowSums(seen[, sim$group == g, drop = FALSE]) >= 1L))
        ids <- ids[ok]
    }
    list(ids = ids, object = omicObject(p$omic, masked[ids, , drop = FALSE], sim$assay, sim$design))
}

# ---- verbatim: sim-k185.R:429-465 at 735c2ea5 ----
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
        out <- data.frame(id = rownames(tab), p = tab$raw_pvalue, q = tab$fdr_qvalue, stringsAsFactors = FALSE)
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

# ==== KIT-217 (written for this measurement; not from KIT-185) ======================================================

# ---- KIT-217 addition (i): the cells, their masks and arms (KIT-217 4.1) -------------------------------------------
CELLS <- list(
    L2 = cell(3201L, "lipid", 500L, 2L, 6L, assays = 2L), L1 = cell(3202L, "lipid", 500L, 2L, 32L, assays = 2L),
    M2 = cell(3203L, "metab", 150L, 4L, 3L), M1 = cell(3204L, "metab", 150L, 4L, 45L),
    Lv = cell(3205L, "lipid", 500L, 2L, 6L, assays = 2L), Le = cell(3206L, "lipid", 500L, 2L, 6L, assays = 2L),
    Me = cell(3207L, "metab", 150L, 4L, 3L))
stopifnot(all(cells_wanted %in% names(CELLS)))
DECIDING_CELLS <- c("L2", "L1", "M2", "M1")
DECIDING_MASKS <- c("none", "vd5", "vd10", "vd25", "mcar10", "vdg10", "vdg25")
masksOf217 <- function(cellname) switch(cellname, Lv = c("vd10", "vd25"), Le = , Me = c("vd10", "vd25", "vdg25"), DECIDING_MASKS)
# D-k's k: trxFilterMinSampleSize of the cell's design (edgeR's filterByExpr rule; checked in E6)
K_DECLARED <- c(L2 = 6, L1 = 25.4, M2 = 3, M1 = 34.5, Lv = 6, Le = 6, Me = 3)
PERCENTILE <- 10
PROPORTION <- 0.5
M_PER_GROUP <- 2L
COUNTED_R <- c(L2 = 200L, L1 = 200L, M2 = 200L, M1 = 200L, Lv = 100L, Le = 100L, Me = 100L)
armsOf <- function(cellname) {
    arms <- c("NF", "S", "PC", "B", "D", "Dk", "Atau_any", "Atau_g2", "Aobs_any", "Aobs_g2")
    if (CELLS[[cellname]]$groups > 2L) c(arms, "Atau_all", "Aobs_all") else arms
}

# ---- the in-script rules (KIT-217 section 3), at each assay's own threshold ------------------------------------------
# tau: the assay's PERCENTILE-th type-7 quantile of every finite raw value of the assay (L-1), one value per row.
assayThresholds <- function(raw, assay, percentile) {
    tau <- numeric(nrow(raw))
    for (a in unique(assay)) {
        rows <- assay == a
        tau[rows] <- stats::quantile(as.vector(raw[rows, , drop = FALSE]), probs = percentile / 100, na.rm = TRUE,
            names = FALSE)[[1L]]
    }
    tau
}
# raw: features x samples on the raw scale, NA = missing; group: one label per sample; tau: per-row threshold.
# b = finite values strictly below tau; f = finite values; per group f_i (finite) and a_i (finite, not below tau);
# need_i = min(n_i, m). Returns one logical vector per rule (S is the independent oracle of today's rule, E1).
ruleSets <- function(raw, group, tau, proportion, m, k) {
    n <- ncol(raw)
    fin <- is.finite(raw)
    below <- fin & raw < tau
    below[is.na(below)] <- FALSE
    f <- rowSums(fin)
    b <- rowSums(below)
    share <- b / n
    groups <- unique(group)
    need <- vapply(groups, function(g) min(sum(group == g), m), numeric(1L))
    perGroup <- function(x) matrix(vapply(groups, function(g) rowSums(x[, group == g, drop = FALSE]), numeric(nrow(raw))),
        nrow = nrow(raw))
    f_g <- perGroup(fin)
    a_g <- perGroup(fin & !below)
    n_obs <- rowSums(sweep(f_g, 2L, need, ">="))
    n_tau <- rowSums(sweep(a_g, 2L, need, ">="))
    G <- length(groups)
    low_ok <- share < proportion
    list(S = low_ok, B = f >= 1 & low_ok, D = (b + n - f) / n < proportion, Dk = f >= k & low_ok,
        Atau_any = n_tau >= 1L, Atau_g2 = n_tau >= 2L, Atau_all = n_tau >= G,
        Aobs_any = n_obs >= 1L & low_ok, Aobs_g2 = n_obs >= 2L & low_ok, Aobs_all = n_obs >= G & low_ok)
}
# the effective lod per sample of a value-dependent mask (for the on/off classification of planted features)
lodEffective <- function(sim, mask) {
    mech <- maskMech(mask)
    if (!mech %in% c("vd", "vdg")) return(NULL)
    rate <- maskRate(mask)
    if (identical(mech, "vd")) return(rep(solveLod(sim$y, rate), ncol(sim$y)))
    shift <- matrix(0.5 * (sim$group == "control"), nrow(sim$y), ncol(sim$y), byrow = TRUE)
    solveLod(sim$y, rate, shift) + 0.5 * (sim$group == "control")
}
attemptK <- function(expr) {
    tryCatch(list(value = expr, status = "ok", class = NA_character_), error = function(e) list(value = NULL, status = "dnc",
        class = substr(gsub("[\r\n\t]+", " ", paste(class(e)[[1L]], conditionMessage(e), sep = ": ")), 1L, 200L)))
}
# Benjamini-Hochberg on the arm's raw p per assay at 0.05 (VAL-11's matched rule)
bhCalls <- function(tab, ids) {
    m <- match(ids, tab$id)
    p <- tab$p[m]
    a <- tab$assay[m]
    called <- rep(FALSE, length(ids))
    tested <- is.finite(p)
    for (as in unique(a[tested])) {
        w <- which(tested & a == as)
        called[w] <- stats::p.adjust(p[w], method = "BH") < 0.05
    }
    called
}

ROW_COLUMNS <- c("cell", "rep", "mask", "arm", "status", "class", "kept", "kept_novalue", "kept_zerodf", "kept_nogroup",
    "miss_before", "miss_after", "n_tested", "n_null_tested", "n_null_rej", "n_null_low", "n_null_low_rej", "n_null_nolow",
    "n_null_nolow_rej", "n_planted", "n_planted_tested", "n_bh_calls", "n_bh_true", "n_bh_false", "fdp_bh", "n_q_calls",
    "n_q_true", "n_q_false", "fdp_q", "qmethod", "n_changed_vs_S", "sum_fc_null", "n_fc_null", "n_up", "n_down", "rm_up",
    "rm_down", "n_onoff", "rm_onoff", "nK", "nK_rej_X", "nK_NF_tested", "nK_rej_NF", "da_failed", "da_warn", "tau", "e1",
    "e2", "f1", "f1_note", "e3", "e3_maxdiff", "e6fc")
newRowK <- function(...) {
    row <- stats::setNames(rep(list(NA), length(ROW_COLUMNS)), ROW_COLUMNS)
    values <- lapply(list(...), function(v) if (is.null(v) || !length(v)) NA else if (length(v) > 1L) paste(v, collapse = ";") else v)
    row[names(values)] <- values
    row
}

# ---- one (replicate, mask): every arm --------------------------------------------------------------------------------
maskUnit <- function(sim, cellname, r, mask) {
    ids <- rownames(sim$y)
    group <- sim$group
    ctrl <- group == "control"
    trt <- group == "treated"
    comp <- ctrl | trt
    planted <- sim$beta[ids] != 0
    up <- sim$beta[ids] > 0
    k <- K_DECLARED[[cellname]]
    missing <- maskMatrix(sim, mask)
    masked <- sim$y
    masked[missing] <- NA_real_
    raw <- 2^masked
    tau <- assayThresholds(raw, sim$assay[ids], PERCENTILE)
    rules <- ruleSets(raw, group, tau, PROPORTION, M_PER_GROUP, k)
    fin <- is.finite(masked)
    f <- rowSums(fin)
    groups_with <- rowSums(matrix(vapply(unique(group), function(g) rowSums(fin[, group == g, drop = FALSE]) > 0,
        logical(length(ids))), nrow = length(ids)))
    zero_df <- f >= 1 & f == groups_with
    nogroup <- rowSums(fin[, ctrl, drop = FALSE]) == 0 | rowSums(fin[, trt, drop = FALSE]) == 0
    low_comp <- rowSums(fin[, comp, drop = FALSE] & raw[, comp, drop = FALSE] < tau, na.rm = TRUE) > 0
    miss_comp <- rowSums(!fin[, comp, drop = FALSE]) > 0
    fc <- rowMeans(masked[, trt, drop = FALSE], na.rm = TRUE) - rowMeans(masked[, ctrl, drop = FALSE], na.rm = TRUE)
    lod <- lodEffective(sim, mask)
    onoff <- if (is.null(lod)) rep(NA, length(ids)) else {
        under <- sweep(sim$y, 2L, lod, "<")
        rowSums(under[, ctrl, drop = FALSE]) == sum(ctrl) | rowSums(under[, trt, drop = FALSE]) == sum(trt)
    }
    # the kept sets: S from the package's S4 method (KIT-185's filterArm), the others from the in-script rules
    S_pkg <- attemptK(filterArm(sim, masked, FALSE)$ids)
    sets <- list(NF = ids, S = S_pkg$value)
    for (rule in setdiff(names(rules), "S")) sets[[rule]] <- ids[rules[[rule]]]
    # E1 (oracle), F1 (subset identities), E2 (mask none)
    e1 <- e2 <- f1 <- NA
    f1_fail <- character()
    if (identical(S_pkg$status, "ok")) {
        S <- S_pkg$value
        e1 <- setequal(S, ids[rules$S]) && !anyDuplicated(S)
        family <- c("B", "D", "Dk", "Aobs_any", "Aobs_g2", "Aobs_all")
        for (a in family) if (!all(sets[[a]] %in% S)) f1_fail <- c(f1_fail, paste0(a, " not in S"))
        if (!all(sets$Aobs_g2 %in% sets$Aobs_any)) f1_fail <- c(f1_fail, "Aobs_g2 not in Aobs_any")
        if (identical(mask, "none")) {
            same <- vapply(family, function(a) setequal(sets[[a]], S), logical(1L))
            e2 <- all(same)
            if (!e2) f1_fail <- c(f1_fail, paste0(family[!same], " differs from S on none"))
        }
        f1 <- !length(f1_fail)
    }
    # the DA of every arm (NF first: PC is read from it)
    runDA <- function(arm, kept) {
        if (fault && identical(arm, "D") && identical(mask, "vd10") && r == first_rep && identical(cellname, cells_wanted[[1L]])) {
            stop("deliberate fault (K217_FAULT, E4)")
        }
        daTable(sim, masked[kept, sim$runs, drop = FALSE])
    }
    results <- list(NF = attemptK(runDA("NF", sets$NF)))
    pNF <- rep(NA_real_, length(ids))
    if (identical(results$NF$status, "ok")) {
        pNF <- results$NF$value$p[match(ids, results$NF$value$id)]
        testedNF <- is.finite(pNF)
        threshold <- stats::median(abs(fc[testedNF]))
        sets$PC <- ids[testedNF & abs(fc) >= threshold]
    }
    for (arm in setdiff(armsOf(cellname), "NF")) {
        results[[arm]] <- if (is.null(sets[[arm]])) {
            list(value = NULL, status = "dnc", class = "kept set unavailable (S or NF did not complete)")
        } else attemptK(runDA(arm, sets[[arm]]))
    }
    called <- lapply(results, function(x) if (identical(x$status, "ok")) bhCalls(x$value, ids))
    # E3: B's tested rows and raw p equal S's
    e3 <- e3_max <- NA
    if (identical(results$B$status, "ok") && identical(results$S$status, "ok")) {
        pB <- results$B$value$p[match(ids, results$B$value$id)]
        pS <- results$S$value$p[match(ids, results$S$value$id)]
        same_tested <- identical(is.finite(pB), is.finite(pS))
        e3_max <- if (any(is.finite(pB))) max(abs(pB - pS)[is.finite(pB)]) else 0
        e3 <- same_tested && isTRUE(e3_max <= 1e-12)
    }
    # E6 (fc): the in-script fold change equals the DA's logFC on NF's tested rows (lifecycle, mask vd10)
    e6fc <- NA
    if (lifecycle && identical(mask, "vd10") && identical(results$NF$status, "ok")) {
        e6fc <- tryCatch({
            object <- omicObject(sim$p$omic, masked[ids, sim$runs, drop = FALSE], sim$assay, sim$design)
            contrasts <- data.frame(contrasts = CONTRAST, friendly_names = "treated_vs_control", stringsAsFactors = FALSE)
            value <- quiet(if (identical(sim$p$omic, "lipid")) runLipidsDA(object, contrasts, "~ 0 + group", 0.05, 0, TRUE, TRUE) else
                runMetabolitesDA(object, contrasts, "~ 0 + group", 0.05, 0, TRUE, TRUE))
            long <- as.data.frame(value[[if (identical(sim$p$omic, "lipid")) "da_lipids_long" else "da_metabolites_long"]])
            id_column <- if (identical(sim$p$omic, "lipid")) "lipid_id" else "metabolite_id"
            lfc <- long$logFC[match(ids, as.character(long[[id_column]]))]
            ok <- is.finite(pNF)
            any(ok) && max(abs(lfc[ok] - fc[ok])) <= 1e-10
        }, error = function(e) FALSE)
    }
    common <- list(cell = cellname, rep = r, mask = mask)
    rows <- list(do.call(newRowK, c(common, list(arm = "mask", status = "ok", miss_before = mean(missing),
        tau = signif(unique(tau), 15), e1 = e1, e2 = e2, f1 = f1, f1_note = if (length(f1_fail)) f1_fail else NA,
        e3 = e3, e3_maxdiff = e3_max, e6fc = e6fc))))
    for (arm in armsOf(cellname)) {
        res <- results[[arm]]
        if (!identical(res$status, "ok")) {
            rows[[length(rows) + 1L]] <- do.call(newRowK, c(common, list(arm = arm, status = "dnc", class = res$class)))
            next
        }
        tab <- res$value
        inK <- ids %in% sets[[arm]]
        m <- match(ids, tab$id)
        p <- tab$p[m]
        q <- tab$q[m]
        tested <- inK & is.finite(p)
        nt <- tested & !planted
        rej <- nt & p < 0.05
        rej[is.na(rej)] <- FALSE
        calls <- called[[arm]]
        qcall <- tested & is.finite(q) & q < 0.05
        qcall[is.na(qcall)] <- FALSE
        kset <- nt & (miss_comp | low_comp)
        fcok <- inK & !planted & is.finite(fc)
        rows[[length(rows) + 1L]] <- do.call(newRowK, c(common, list(arm = arm, status = "ok", kept = sum(inK),
            kept_novalue = sum(inK & f == 0), kept_zerodf = sum(inK & zero_df), kept_nogroup = sum(inK & nogroup),
            miss_before = mean(missing), miss_after = if (any(inK)) mean(!fin[inK, , drop = FALSE]) else NA_real_,
            n_tested = sum(tested), n_null_tested = sum(nt), n_null_rej = sum(rej), n_null_low = sum(nt & low_comp),
            n_null_low_rej = sum(rej & low_comp), n_null_nolow = sum(nt & !low_comp), n_null_nolow_rej = sum(rej & !low_comp),
            n_planted = sum(planted), n_planted_tested = sum(tested & planted), n_bh_calls = sum(calls),
            n_bh_true = sum(calls & planted), n_bh_false = sum(calls & !planted),
            fdp_bh = if (any(calls)) sum(calls & !planted) / sum(calls) else 0, n_q_calls = sum(qcall),
            n_q_true = sum(qcall & planted), n_q_false = sum(qcall & !planted),
            fdp_q = if (any(qcall)) sum(qcall & !planted) / sum(qcall) else 0, qmethod = attr(tab, "qmethod"),
            n_changed_vs_S = if (is.null(called$S)) NA else sum(calls != called$S),
            sum_fc_null = sum(fc[fcok]), n_fc_null = sum(fcok), n_up = sum(planted & up), n_down = sum(planted & !up),
            rm_up = sum(planted & up & !inK), rm_down = sum(planted & !up & !inK),
            n_onoff = if (is.null(lod)) NA else sum(planted & onoff), rm_onoff = if (is.null(lod)) NA else sum(planted & onoff & !inK),
            nK = sum(kset), nK_rej_X = sum(kset & rej), nK_NF_tested = sum(kset & is.finite(pNF)),
            nK_rej_NF = sum(kset & is.finite(pNF) & pNF < 0.05, na.rm = TRUE), da_failed = attr(tab, "failed"),
            da_warn = attr(tab, "warnings"))))
    }
    rows
}

# ---- one replicate: every mask ---------------------------------------------------------------------------------------
unitRep <- function(cellname, r) {
    t0 <- proc.time()[[3L]]
    warning_count$n <- 0L
    sim <- generateReplicate(cellname, r, strength = 0)
    rows <- list()
    for (mask in masksOf217(cellname)) {
        rows <- c(rows, tryCatch(maskUnit(sim, cellname, r, mask), error = function(e) {
            cls <- substr(gsub("[\r\n\t]+", " ", paste(class(e)[[1L]], conditionMessage(e), sep = ": ")), 1L, 200L)
            lapply(c("mask", armsOf(cellname)), function(a) newRowK(cell = cellname, rep = r, mask = mask, arm = a,
                status = "dnc", class = cls))
        }))
    }
    table <- do.call(rbind, lapply(rows, as.data.frame, stringsAsFactors = FALSE))[, ROW_COLUMNS]
    list(rows = table, seconds = proc.time()[[3L]] - t0, warnings = warning_count$n)
}

# ---- E6: the in-script rules on 2.4's ten-feature table, T3 and T4; D-k's k from the package -----------------------
checks_rows <- list()
addCheckK <- function(check, cell, pass, n = NA, note = "") {
    checks_rows[[length(checks_rows) + 1L]] <<- data.frame(check = check, cell = cell, pass = pass, n = n, note = note,
        stringsAsFactors = FALSE)
}
patternMatrix <- function(patterns) {
    rows <- lapply(patterns, function(s) {
        ch <- strsplit(gsub("|", "", s, fixed = TRUE), "")[[1L]]
        c(`+` = 1000, `-` = 10, `.` = NA_real_)[ch]
    })
    m <- do.call(rbind, rows)
    dimnames(m) <- list(names(patterns), NULL)
    m
}
checkE6 <- function() {
    table24 <- patternMatrix(c(F1 = "++++++|++++++", F2 = "------|++++++", F3 = "......|......", F4 = "++++++|......",
        F5 = "+++...|+++...", F6 = "+.....|+.....", F7 = "++....|++++++", F8 = "++....|......", F9 = "---...|+++...",
        F10 = "-----+|-----+"))
    rules <- ruleSets(table24, rep(c("A", "B"), each = 6L), rep(100, 10L), 0.5, 2L, 6)
    expected <- list(S = c(1, 3, 4, 5, 6, 7, 8, 9), B = c(1, 4, 5, 6, 7, 8, 9), D = c(1, 7), Dk = c(1, 4, 5, 7, 9),
        Atau_any = c(1, 2, 4, 5, 7, 8, 9), Atau_g2 = c(1, 5, 7), Atau_all = c(1, 5, 7), Aobs_any = c(1, 4, 5, 7, 8, 9),
        Aobs_g2 = c(1, 5, 7, 9), Aobs_all = c(1, 5, 7, 9))
    for (rule in names(expected)) {
        got <- which(rules[[rule]])
        addCheckK("E6_table24", rule, identical(as.numeric(got), as.numeric(expected[[rule]])), length(got))
    }
    t3 <- patternMatrix(c(G1 = "+++|+++|...|...", G2 = "+++|...|...|...", G3 = "++.|++.|++.|++.", G4 = "+..|+..|+..|+.."))
    r3 <- ruleSets(t3, rep(c("a", "b", "c", "d"), each = 3L), rep(100, 4L), 0.5, 2L, 3)
    addCheckK("E6_T3", "g2", identical(rownames(t3)[r3$Aobs_g2], c("G1", "G3")))
    addCheckK("E6_T3", "g1", identical(rownames(t3)[r3$Aobs_any], c("G1", "G2", "G3")))
    addCheckK("E6_T3", "g4", identical(rownames(t3)[r3$Aobs_all], "G3"))
    t4 <- patternMatrix(c(U1 = "+++|++++..", U2 = "++.|++++++"))
    r4 <- ruleSets(t4, rep(c("a", "b"), c(3L, 6L)), rep(100, 2L), 0.5, 4L, 3)
    addCheckK("E6_T4", "g2", identical(rownames(t4)[r4$Aobs_g2], "U1"))
    for (cellname in names(CELLS)) {
        p <- CELLS[[cellname]]
        k <- tryCatch(trxFilterMinSampleSize(rep(GROUPS[seq_len(p$groups)], each = p$per)), error = function(e) NA_real_)
        addCheckK("E6_k", cellname, isTRUE(abs(k - K_DECLARED[[cellname]]) < 1e-12), note = "trxFilterMinSampleSize")
    }
}

# the rows of one job, one gzipped table per cell: <prefix>-rows-<cell>-<first>-<last>.tsv.gz
writeRows <- function(all_rows, dir, prefix, cells, first, last) {
    dir.create(dir, recursive = TRUE, showWarnings = FALSE)
    for (cellname in cells) {
        table <- all_rows[all_rows$cell == cellname, , drop = FALSE]
        table <- table[order(table$rep, match(table$mask, masksOf217(cellname)), match(table$arm, c("mask", armsOf(cellname)))), ]
        connection <- gzfile(file.path(dir, sprintf("%s-rows-%s-%d-%d.tsv.gz", prefix, cellname, first, last)), "w")
        utils::write.table(table, connection, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
        close(connection)
    }
}

# ---- the job ---------------------------------------------------------------------------------------------------------
units <- expand.grid(rep = seq(first_rep, last_rep), cell = cells_wanted, stringsAsFactors = FALSE)
job <- sprintf("%s-%s-%d-%d%s", job_prefix, paste(cells_wanted, collapse = "_"), first_rep, last_rep, if (lifecycle) "-l1" else "")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
cat(sprintf("KIT-217 %s: %d replicate units, %d workers, lifecycle %s, fault %s\n", job, nrow(units), workers, lifecycle, fault))
checkE6()
runOne <- function(i) {
    u <- units[i, ]
    result <- tryCatch(unitRep(u$cell, u$rep), error = function(e) {
        cls <- substr(gsub("[\r\n\t]+", " ", paste(class(e)[[1L]], conditionMessage(e), sep = ": ")), 1L, 200L)
        rows <- do.call(rbind, lapply(masksOf217(u$cell), function(mask) do.call(rbind, lapply(c("mask", armsOf(u$cell)),
            function(a) as.data.frame(newRowK(cell = u$cell, rep = u$rep, mask = mask, arm = a, status = "dnc", class = cls),
                stringsAsFactors = FALSE)))))[, ROW_COLUMNS]
        list(rows = rows, seconds = NA_real_, warnings = NA_integer_)
    })
    result
}
results <- parallel::mclapply(seq_len(nrow(units)), runOne, mc.cores = workers, mc.preschedule = FALSE)
for (i in seq_along(results)) if (!is.list(results[[i]]) || is.null(results[[i]]$rows)) {
    u <- units[i, ]
    cls <- substr(paste("worker failed:", paste(as.character(results[[i]]), collapse = " ")), 1L, 200L)
    rows <- do.call(rbind, lapply(masksOf217(u$cell), function(mask) do.call(rbind, lapply(c("mask", armsOf(u$cell)),
        function(a) as.data.frame(newRowK(cell = u$cell, rep = u$rep, mask = mask, arm = a, status = "dnc", class = cls),
            stringsAsFactors = FALSE)))))[, ROW_COLUMNS]
    results[[i]] <- list(rows = rows, seconds = NA_real_, warnings = NA_integer_)
}
sha <- function(f) strsplit(system2("sha256sum", f, stdout = TRUE), " ")[[1L]][[1L]]
rowsSha <- function(rows) {
    temp <- tempfile(fileext = ".tsv")
    utils::write.table(rows, temp, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
    on.exit(unlink(temp))
    sha(temp)
}
timing <- do.call(rbind, lapply(seq_along(results), function(i) data.frame(cell = units$cell[[i]], rep = units$rep[[i]],
    seconds = results[[i]]$seconds, warnings = results[[i]]$warnings,
    rows_sha256 = if (lifecycle) rowsSha(results[[i]]$rows) else NA_character_, stringsAsFactors = FALSE)))
utils::write.table(timing, file.path(out_dir, sprintf("timing-%s.tsv", job)), sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
all_rows <- do.call(rbind, lapply(results, `[[`, "rows"))
arm_rows <- all_rows[all_rows$arm != "mask", , drop = FALSE]
mask_rows <- all_rows[all_rows$arm == "mask", , drop = FALSE]
# checks on every job (E1 and F1 are stop rules in counted runs; E3 reported)
for (cellname in cells_wanted) {
    mr <- mask_rows[mask_rows$cell == cellname & mask_rows$status == "ok", , drop = FALSE]
    addCheckK("E1", cellname, all(mr$e1 %in% TRUE) && nrow(mr) > 0L, nrow(mr), "package S set equals the oracle")
    addCheckK("F1", cellname, all(mr$f1 %in% TRUE) && nrow(mr) > 0L, nrow(mr), "subset identities (and none: equal to S)")
    if (cellname %in% DECIDING_CELLS) {
        none <- mr[mr$mask == "none", , drop = FALSE]
        addCheckK("E2", cellname, all(none$e2 %in% TRUE) && nrow(none) > 0L, nrow(none), "mask none: B, D, Dk, Aobs equal S")
    }
    addCheckK("E3", cellname, all(mr$e3 %in% TRUE), sum(!is.na(mr$e3)), "B's tested rows and raw p equal S's (1e-12)")
    if (lifecycle) addCheckK("E6_fc", cellname, all(mr$e6fc[mr$mask == "vd10"] %in% TRUE), sum(mr$mask == "vd10"),
        "in-script fold change equals the DA's logFC (1e-10)")
    units_total <- sum(all_rows$cell == cellname)
    addCheckK("units", cellname, TRUE, units_total, sprintf("did not complete: %d of %d rows", sum(all_rows$cell == cellname &
        all_rows$status != "ok"), units_total))
}
if (lifecycle) {
    # E4: the deliberate fault is recorded as did not complete and the job went on
    if (fault) {
        hit <- arm_rows[arm_rows$cell == cells_wanted[[1L]] & arm_rows$rep == first_rep & arm_rows$mask == "vd10" &
            arm_rows$arm == "D", , drop = FALSE]
        others <- arm_rows[!(arm_rows$cell == cells_wanted[[1L]] & arm_rows$rep == first_rep & arm_rows$mask == "vd10" &
            arm_rows$arm == "D"), , drop = FALSE]
        expected_rows <- sum(vapply(cells_wanted, function(cn) length(masksOf217(cn)) * length(armsOf(cn)), numeric(1L))) *
            (last_rep - first_rep + 1L)
        addCheckK("E4", cells_wanted[[1L]], nrow(hit) == 1L && identical(hit$status, "dnc") &&
            grepl("deliberate fault", hit$class, fixed = TRUE) && nrow(arm_rows) == expected_rows && sum(timing$seconds > 0, na.rm = TRUE) ==
            nrow(units), nrow(arm_rows), sprintf("arm rows %d of %d; other rows not completed %d", nrow(arm_rows), expected_rows,
            sum(others$status != "ok")))
    }
    # E5: the counted jobs' planned minutes from these seconds (11 minutes of set-up, 4 workers); the timing rule
    plan <- do.call(rbind, lapply(cells_wanted, function(cn) {
        s <- mean(timing$seconds[timing$cell == cn], na.rm = TRUE)
        data.frame(cell = cn, seconds_per_replicate = round(s, 2), counted_r = COUNTED_R[[cn]],
            compute_minutes = round(s * COUNTED_R[[cn]] / workers / 60, 1), stringsAsFactors = FALSE)
    }))
    jobs <- data.frame(job = c(DECIDING_CELLS, "reported", "determinism"), minutes = c(
        11 + plan$compute_minutes[match(DECIDING_CELLS, plan$cell)],
        11 + sum(plan$compute_minutes[plan$cell %in% c("Lv", "Le", "Me")]),
        11 + plan$seconds_per_replicate[plan$cell == "L2"] * 5 / workers / 60), stringsAsFactors = FALSE)
    jobs$minutes <- round(jobs$minutes, 1)
    utils::write.table(plan, file.path(out_dir, sprintf("timing-plan-%s.tsv", job)), sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
    utils::write.table(jobs, file.path(out_dir, sprintf("job-plan-%s.tsv", job)), sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
    addCheckK("E5", "all", all(jobs$minutes <= 150, na.rm = TRUE) && !anyNA(jobs$minutes), nrow(jobs),
        sprintf("planned job minutes, the largest %.1f (timing rule: a job above 150 is split at the same R)", max(jobs$minutes)))
    # blind: completion only
    completion <- arm_rows[, c("cell", "rep", "mask", "arm", "status", "class")]
    utils::write.table(completion, file.path(out_dir, sprintf("completion-%s.tsv", job)), sep = "\t", quote = FALSE,
        row.names = FALSE, na = "NA")
    # D-02: the counted writing path, run on the lifecycle's rows, to $RUNNER_TEMP only (never under the uploaded out
    # dir); and L2's first replicate run a second time as a determinism file, for E7 (score-k217.R e7)
    temp <- Sys.getenv("RUNNER_TEMP")
    if (!nzchar(temp)) temp <- tempdir()
    ml_rows <- file.path(temp, "k217-ml-rows")
    writeRows(all_rows, ml_rows, "m", cells_wanted, first_rep, last_rep)
    if ("L2" %in% cells_wanted) {
        again <- runOne(which(units$cell == "L2" & units$rep == first_rep))
        writeRows(again$rows, ml_rows, "det", "L2", first_rep, first_rep)
    }
} else {
    writeRows(all_rows, out_dir, job_prefix, cells_wanted, first_rep, last_rep)
}
checks <- do.call(rbind, checks_rows)
utils::write.table(checks, file.path(out_dir, sprintf("checks-%s.tsv", job)), sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
gate <- checks$check %in% c("E1", "E2", "E4", "E6_table24", "E6_T3", "E6_T4", "E6_k", "E6_fc")
cat(sprintf("KIT-217 %s gate G-M (E1, E2, E4, E6): %s\n", job, if (lifecycle) (if (all(checks$pass[gate])) "pass" else "FAIL") else "not a lifecycle job"))
# provenance
library_table <- utils::installed.packages()[, c("Package", "Version")]
key_packages <- c("MultiScholaR", "limma", "qvalue", "edgeR", "dplyr", "tidyr", "tibble")
tables <- list.files(out_dir, full.names = TRUE, pattern = "\\.(tsv|gz)$")
provenance <- list(job = job, cells = cells_wanted, replicates = c(first_rep, last_rep), lifecycle = lifecycle, fault = fault,
    workers = workers, commit = tryCatch(system2("git", c("rev-parse", "HEAD"), stdout = TRUE), error = function(e) NA),
    r_version = R.version.string, blas = extSoftVersion()[["BLAS"]], lapack = La_library(),
    runner_image = paste(Sys.getenv("ImageOS"), Sys.getenv("ImageVersion")),
    processor = tryCatch(sub("^model name\\s*:\\s*", "", grep("^model name", readLines("/proc/cpuinfo"), value = TRUE)[[1L]]), error = function(e) NA),
    seed_rule = "set.seed(185000000 + 10000 * (3200 + i) + r, kind = 'Mersenne-Twister', normal.kind = 'Inversion', sample.kind = 'Rejection') = 217000000 + 10000 i + r",
    cell_index = vapply(CELLS[cells_wanted], `[[`, integer(1L), "idx"),
    key_packages = as.list(stats::setNames(library_table[match(key_packages, library_table[, "Package"]), "Version"], key_packages)),
    script_sha256 = sha(script_file), started = format(started_at, "%Y-%m-%dT%H:%M:%S%z"),
    finished = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"), replicate_units = nrow(units),
    rows_not_completed = sum(all_rows$status != "ok"),
    tables_sha256 = stats::setNames(as.list(vapply(tables, sha, character(1L))), basename(tables)))
jsonlite::write_json(provenance, file.path(out_dir, sprintf("provenance-%s.json", job)), auto_unbox = TRUE, pretty = TRUE, digits = NA)
cat(sprintf("KIT-217 %s finished: %d replicate units, wall %.1f s\n", job, nrow(units), as.numeric(difftime(Sys.time(), started_at, units = "secs"))))
