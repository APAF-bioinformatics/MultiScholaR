# KIT-188 (VAL-72, VAL-81, VAL-84): calibration of the default route with legacy RUV-III Automatic, lipidomics,
# metabolomics (Part L) and proteomics (Part P). Evidence code under docs/qa: no package code. Runs on a GitHub-hosted
# runner with MultiScholaR installed from the bound commit (run-tests `scripts`), from the package root.
#
# Usage:
#   Rscript sim-k188.R <part> <cells> <first> <last> <out dir>       the job (parent)
#   Rscript sim-k188.R --unit <part> <cell> <r> <job> <scenario> <mask> <rds>   one unit (child process)
#   part   L (lipidomics and metabolomics: L12, L64, M12, M168) or P (proteomics: P0, P1, P16, P2, P3, P4b)
#   cells  one cell or a comma list
#   first, last  data sets (S: 1 to 400; lifecycle: 401 to 403)
# Environment:
#   K188_JOBS      comma list of job types (Part L default deciding,reported,munit,ba; Part P deciding,reported,ba):
#                  deciding (LA, GA, NO, OR), reported (LA, OR, LM, FP, SC, LCm in Part L; LM, OR in Part P),
#                  munit (Part M unit check, L12 and M12), ba (mask BA, L64, M168, P4b), lp (LP: the LA request on the
#                  bound commit with KIT-190's markers at their earlier tokens, 5.3; r1.21), las (LAs: the LA request with
#                  KIT-224's marker at its earlier token, 34), l1c (the planted cases of 35.6, L12, one data set, K188_PLANT).
#                  L24 (38.5, reported) runs the deciding job type only. A job type that does not apply to a cell is not run.
#   K188_PLANT     L1c: comma list of the planted cases (itsd, itsdpos, emptysample, ruvskip, ruvdrop, daerror, normplant,
#                  exitplant; default all), one child process each; needs K188_JOBS=l1c and K188_LIFECYCLE=1.
#   K188_PREDICTIONS  L1c and E7 (R-marked predictions of 35.6 and 7.4, read at the bound commit): absolute path of the TSV
#                  (columns key = "<plant or check>/<field>", expected) written into lifecycle.md before the run.
#   K188_LIFECYCLE 1: blind lifecycle (L1, L1p): every estimand column is written NA; the rows' SHA-256 is written; the
#                  checks E1 to E4, E6, E7, E9 (and E8, E2, E5r) are run and written to checks-<job>.tsv.
#   K188_FAULT     1: E7 (crash): one deliberate error in the first unit's first journey.
#   K188_WORKERS   child processes at a time (default 4).
#   K188_P_REPS    Part P lifecycle: comma list cell=last data set (default "P0=403,P2=403"; others run first only).
#   K188_JOB       a name for the job's files (default from part, cells and range).
#
# One data set (Part L: per job type) or one (data set, scenario, mask) (Part P) per Rscript child process, four at a
# time, so that no session owner (duckdb) is forked. Every random draw of a data set is made in generate() before any
# package function is called; seeds set.seed(188000000 + 100000 c + 10000 s + r, kind = "Mersenne-Twister",
# normal.kind = "Inversion", sample.kind = "Rejection"), c the cell (L12 1, L64 2, M12 3, M168 4, P0 5, P1 6, P16 7,
# P2 8, P3 9, P4b 10), s the scenario (N0 1, N1 2, E0 3, E1 4, N2 5, E2 6), r the data set.
# The design, estimands and reading rules are KIT-188 revision 1.4, sections 5, 7, 18 and the declaration.
# Construction choices this file fixes (5.2 leaves them open): each unwanted component is centred with mean square 1
# over the runs (population variance); components are made orthogonal to each other in the order run factor or
# plate/batch basis, drift, sample factor; the drift's order is, of 2,000 random orders (within plate in L64), the one
# whose drift has the smallest share of variance explained by the groups, and is then residualised; labelled controls
# draw their levels and variances as features (the variance scaled by the features' median); every draw is made in
# every scenario (a scenario without unwanted variation multiplies it by 0).

# Partial-match guard (smoke check send 3, the coordinator's ruling; repaired as KIT-275's): `$` on a list or data frame matches a name by its prefix when
# none is exact; that warning becomes an error (`stop`). A `stop()` raised from a GLOBAL calling handler runs with every local handler disabled, so it
# reaches no `tryCatch` and ends the process with nothing written. The guard is therefore a function (`guardWarning`) installed as a LOCAL calling
# handler around every journey, every unit and inside `quiet()`, where its stop reaches the enclosing `tryCatch`; the global handler stays only as a
# last net, with an `options(error=)` writer that leaves a trace file before the non-zero exit (the parent reads it and records an instrument unit row).
# Where each path stands: a partial match inside a journey (inside `quiet()` or the journey's own `tryCatch`) -> local handler -> the journey becomes
# `instrument`; inside a unit outside any journey -> the child's local handler -> the unit row is `instrument`; in the parent or the scorer -> the global
# net only: the process ends with the message and the trace. A stop from the guard is an instrument defect by definition (VAL-54).
options(warnPartialMatchDollar = TRUE)
guardWarning <- function(w) {
    if (grepl("partial match of", tolower(conditionMessage(w)), fixed = TRUE)) {
        options(k188.guard_trace = paste(vapply(sys.calls(), function(cl) substr(paste(deparse(cl)[[1L]], collapse = ""), 1L, 80L), ""), collapse = " <- "))
        stop(paste("partial-match guard:", conditionMessage(w)), call. = FALSE)
    }
}
globalCallingHandlers(warning = guardWarning)
options(error = function() { # an error that escaped every handler: write what is known, then end non-zero
    cl <- commandArgs(TRUE)
    path <- if (identical(cl[1L], "--unit")) sub("\\.rds$", "-trace.txt", cl[8L]) else file.path(cl[5L], "parent-trace.txt")
    try(dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE), silent = TRUE)
    try(writeLines(c(geterrmessage(), getOption("k188.guard_trace", "")), path), silent = TRUE)
    quit(save = "no", status = 1L)
})
args <- commandArgs(TRUE)
child <- length(args) >= 1L && identical(args[[1L]], "--unit")
script_file <- normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[[1L]]))
started_at <- Sys.time()

suppressPackageStartupMessages(library(MultiScholaR))
attach(asNamespace("MultiScholaR"), name = "MultiScholaR:namespace", warn.conflicts = FALSE)
for (helper in c("helper-workflow-omics-owner.R", "helper-workflow-lipid-declarations.R", "helper-workflow-metab-declarations.R",
    "helper-kit172.R")) {
    sys.source(file.path("tests", "testthat", helper), envir = globalenv())
}
logger::log_threshold(logger::FATAL)
invisible(RUVIIIC::omp_set_num_threads(1L))
options(MultiScholaR.ruv_design_panel = FALSE)

# ---- constants -----------------------------------------------------------------------------------------------------
GROUPS <- LETTERS[1:5]
SCENARIOS <- c("N0", "N1", "E0", "E1", "N2", "E2")
MED_VAR_L <- 0.25
MED_VAR_P <- 4 * 0.25^2 / stats::qchisq(0.5, 4)
cellDef <- function(idx, part, omic, assays, per_assay, groups, per, level, sd, unwanted, blocks = 0L, n2 = FALSE) {
    list(idx = idx, part = part, omic = omic, assays = assays, per_assay = per_assay, groups = groups, per = per, level = level,
        sd = sd, unwanted = unwanted, blocks = blocks, n2 = n2)
}
CELLS <- list(
    L12 = cellDef(1L, "L", "lipid", c("LCMS_Neg", "LCMS_Pos"), 250L, 2L, 6L, 20, 2, "drift_factor", n2 = TRUE),
    L64 = cellDef(2L, "L", "lipid", c("LCMS_Neg", "LCMS_Pos"), 277L, 2L, 32L, 20, 2, "plate_drift", blocks = 2L),
    M12 = cellDef(3L, "L", "metab", "LCMS_Pos", 150L, 4L, 3L, 2.2, 3, "drift_factor", n2 = TRUE),
    M168 = cellDef(4L, "L", "metab", "LCMS_Pos", 148L, 4L, 42L, 2.2, 3, "batch_drift", blocks = 14L),
    P0 = cellDef(5L, "P", "prot", "proteins", 3000L, 2L, 6L, 20, 2, "drift_factor", n2 = TRUE),
    P1 = cellDef(6L, "P", "prot", "proteins", 6000L, 2L, 4L, 20, 2, "drift_factor"),
    P16 = cellDef(7L, "P", "prot", "proteins", 6000L, 4L, 4L, 20, 2, "drift_factor"),
    P2 = cellDef(8L, "P", "prot", "proteins", 3000L, 2L, 3L, 20, 2, "drift_factor"),
    P3 = cellDef(9L, "P", "prot", "proteins", 6000L, 5L, 4L, 20, 2, "drift_factor"),
    P4b = cellDef(10L, "P", "prot", "proteins", 5000L, 5L, 10L, 20, 2, "batch_drift", blocks = 5L),
    # 38.5 (r1.19, VAL-243; reported, cap 3 at 24 runs): L12 with 12 runs per group, seed index c = 11 (5.2), no N2 or E2
    L24 = cellDef(11L, "L", "lipid", c("LCMS_Neg", "LCMS_Pos"), 250L, 2L, 12L, 20, 2, "drift_factor"),
    # 5.8b: KIT-172's concentration fixture, both classes, one fixed data set (reported; job k172c only)
    K172c = cellDef(0L, "L", "lipid", "LCMS_Pos", 515L, 2L, 6L, NA, NA, "none"))
scenariosOf <- function(cellname) if (isTRUE(CELLS[[cellname]]$n2)) SCENARIOS else SCENARIOS[1:4]
MASKS <- c("none", "vd5", "vd20")
REFUSALS <- c("multischolar_refusal", "multischolar_ruv_grouping_unusable", "omics_ruv_insufficient_controls", "prot_ruv_insufficient_controls", "prot_ruv_empty_result",
    "workflow_execution_invalid", "ruv_factor_count_unavailable", "prot_ruv_factors_unavailable", "omics_ruv_factors_refused",
    "omics_ruv_no_residual_df", "ruv_k_refused", "ruv_k_above_limit")
# the columns that carry an estimand or a Part M value (blinded in the lifecycle)
BLIND <- c("n_null", "null_p05", "null_p01", "calls", "false_calls", "true_calls", "any_call", "fdp", "power", "n_planted",
    "bh_calls", "bh_false", "bh_any", "bh_fdp", "bh_power", "err_med", "err_p99", "err_rmse", "err_max", "gross1", "slope",
    "n_common", "s5", "or_s5", "gross1_common", "or_gross1_common", "ctl_planted_share", "ctl_labelled", "mdf_null05",
    "mdf_null05_reduced", "mdf_df_le_k", "mdf_tested", "mamp_sigma_min", "mamp_e_max", "mamp_planted_share",
    "mamp_cf_n", "mamp_cf_singular", "mamp_cf_log10kappa", "tie_share", "tie_block", "tie_block_onegroup", "perm_jaccard",
    "mu_itsd_dev", "mu_log2_dev", "mu_loess_dev", "mu_ctl_same", "mu_k_same", "mu_pct_same", "mu_ruv_dev", "mu_logfc_dev",
    "mu_t_dev", "mu_call_changes", "mu_t_dev_notrend", "logfc_max", "mamp_slope")

warning_count <- new.env()
warning_count$n <- 0L
warning_count$log <- character()
warning_count$assays <- character()
# r1.21 (R20-P08, 40.3 (4)): each warning's class, and each classed message's (a notice raised by rlang::inform), per
# journey and assay: "class@assay", the assay being the one the text names ("" when it names none)
quiet <- function(expr) {
    value <- NULL
    noteOf <- function(cond) {
        text <- conditionMessage(cond)
        hit <- warning_count$assays[vapply(warning_count$assays, grepl, logical(1L), x = text, fixed = TRUE)]
        warning_count$log <- c(warning_count$log, paste(class(cond)[[1L]], if (length(hit)) hit[[1L]] else "", sep = "@"))
    }
    invisible(utils::capture.output(value <- withCallingHandlers(expr, warning = function(w) {
        guardWarning(w)
        warning_count$n <- warning_count$n + 1L
        noteOf(w)
        invokeRestart("muffleWarning")
    }, message = function(m) {
        if (!identical(class(m)[[1L]], "simpleMessage")) noteOf(m)
        invokeRestart("muffleMessage")
    }), type = "output"))
    value
}
warningClasses <- function() {
    if (!length(warning_count$log)) return("")
    counts <- table(warning_count$log)
    paste(sprintf("%s=%d", names(counts), as.integer(counts)), collapse = ";")
}
standardise <- function(v) {
    v <- v - mean(v)
    v / sqrt(mean(v^2))
}
residStd <- function(v, basis) standardise(stats::lm.fit(basis, v)$residuals)
sha <- function(f) strsplit(system2("sha256sum", shQuote(f), stdout = TRUE), " ")[[1L]][[1L]]
byName <- function(l, n) if (length(n) == 1L && !is.na(n) && n %in% names(l)) l[[n]] else NULL
oneLine <- function(x, n = 300L) substr(gsub("[\r\n\t]+", " ", paste(x, collapse = " ")), 1L, n)

# ---- replay markers (40.3 (5), R20-P02) ------------------------------------------------------------------------------------
# The marker options are process-wide and the arms of a child run in sequence, so each journey sets its arm's markers
# itself, restores them on exit and, at its start, asserts that no other marker leaks. The earlier tokens an arm sets are
# the ticket's (5.3: LP `ruv_control_grouping` "dropped" and `ruv_k_bound` "none"; 34: LAs `ruv_feature_levels` "shifted"),
# and each is asserted equal to the marker table's `earlier` at the bound commit; a difference stops the unit.
ARM_MARKERS <- list(
    LP = list(MultiScholaR.ruv_control_grouping = "dropped", MultiScholaR.ruv_k_bound = "none"),
    LAs = list(MultiScholaR.ruv_feature_levels = "shifted"),
    # E10u (the coordinator's ruling on SA-P07): a check-only LA journey with the bound off, whose sweep k is the untouched sweep's
    LAu = list(MultiScholaR.ruv_k_bound = "none"))
armMarkers <- function(arm) {
    set <- ARM_MARKERS[[arm]] %||% list()
    for (option in names(set)) {
        marker <- Filter(function(m) identical(m$option, option), workflowReplayRuleMarkers)[[1L]]
        stopifnot(identical(marker$earlier, set[[option]]))
    }
    set
}
markerLeaks <- function(set = list()) {
    leaks <- character()
    for (marker in workflowReplayRuleMarkers) {
        if (marker$option %in% names(set)) next
        value <- getOption(marker$option)
        if (!is.null(value) && !identical(value, marker$current)) leaks <- c(leaks, sprintf("%s=%s", marker$option, value))
    }
    leaks
}
withMarkers <- function(set, expr) {
    old <- options(set)
    on.exit(options(old), add = TRUE)
    expr
}

# ---- the generator (5.2, 18.3) -------------------------------------------------------------------------------------
generate <- function(cellname, scenario, r) {
    p <- CELLS[[cellname]]
    set.seed(188000000 + 100000 * p$idx + 10000 * match(scenario, SCENARIOS) + r, kind = "Mersenne-Twister",
        normal.kind = "Inversion", sample.kind = "Rejection")
    groups <- GROUPS[seq_len(p$groups)]
    if (p$blocks > 0L) {
        per_gb <- p$per %/% p$blocks
        block <- rep(seq_len(p$blocks), each = p$groups * per_gb)
        group <- rep(rep(groups, each = per_gb), times = p$blocks)
    } else {
        block <- rep(1L, p$groups * p$per)
        group <- rep(groups, each = p$per)
    }
    n_runs <- length(group)
    runs <- sprintf("R%03d", seq_len(n_runs))
    omic <- p$omic
    n_assays <- length(p$assays)
    n_feat <- p$per_assay * n_assays
    n_lbl <- if (identical(omic, "prot")) 0L else 5L * n_assays
    n_is <- if (identical(omic, "lipid")) 5L * n_assays else 0L
    medvar <- if (identical(omic, "prot")) MED_VAR_P else MED_VAR_L
    # 1 levels and residual variances (features, then labelled controls)
    mu <- stats::rnorm(n_feat, p$level, p$sd)
    if (identical(omic, "prot")) {
        s2 <- 4 * 0.25^2 / stats::rchisq(n_feat, 4)
    } else {
        x <- 10 / stats::rchisq(n_feat, 10)
        s2 <- 0.25 * x / stats::median(x)
    }
    mu_lbl <- stats::rnorm(n_lbl, p$level, p$sd)
    x_lbl <- 10 / stats::rchisq(n_lbl, 10)
    s2_lbl <- if (n_lbl > 0L) 0.25 * x_lbl / stats::median(x) else numeric()
    # 2 effects: exactly 20 % of the features, half up; the group (2 groups: B; more: one of B.. drawn per feature)
    n_eff <- round(0.2 * n_feat)
    where <- sample.int(n_feat, n_eff)
    signs <- sample(rep(c(1, -1), c(round(n_eff / 2), n_eff - round(n_eff / 2))))
    target <- if (p$groups == 2L) rep(2L, n_eff) else sample(2:p$groups, n_eff, replace = TRUE)
    beta <- matrix(0, n_feat, p$groups, dimnames = list(NULL, groups))
    if (scenario %in% c("E0", "E1", "E2")) beta[cbind(where, target)] <- signs
    # 3 the run factor of N2 and E2 (phi = 1/3 for each contrast with group A)
    level_a <- logical(n_runs)
    for (g in groups) {
        ix <- which(group == g)
        n_a <- if (p$groups == 2L) (if (g == "A") 2L * length(ix) %/% 3L else length(ix) %/% 3L) else
            (if (g == "A") length(ix) %/% 3L else 2L * length(ix) %/% 3L)
        level_a[ix[sample.int(length(ix), n_a)]] <- TRUE
    }
    # 4 drift: of 2,000 random orders (within plate in L64), the one whose drift is least explained by the groups
    G <- stats::model.matrix(~ factor(group, levels = groups))
    within <- identical(p$unwanted, "plate_drift")
    candidates <- replicate(2000L, if (within) {
        o <- integer(n_runs)
        for (b in unique(block)) {
            ix <- which(block == b)
            o[ix] <- sample.int(length(ix))
        }
        o
    } else sample.int(n_runs))
    span <- if (within) as.numeric(table(block)[as.character(block)]) else rep(n_runs, n_runs)
    drifts <- exp(-2 * (candidates - 1) / (span - 1))
    centred <- sweep(drifts, 2L, colMeans(drifts))
    fitted <- G %*% solve(crossprod(G), crossprod(G, drifts))
    explained <- colSums(sweep(fitted, 2L, colMeans(drifts))^2) / colSums(centred^2)
    drift_raw <- drifts[, which.min(explained)]
    # 5 the sample factor
    factor_raw <- stats::rnorm(n_runs)
    # 6 loadings (features, standards, labelled controls)
    C <- switch(p$unwanted, drift_factor = 2L, plate_drift = 2L, batch_drift = p$blocks)
    sd_load <- sqrt(3 * medvar / C)
    A_feat <- matrix(stats::rnorm(n_feat * C), n_feat, C) * sd_load
    A_is <- matrix(stats::rnorm(n_is * C), n_is, C) * sd_load
    A_lbl <- matrix(stats::rnorm(n_lbl * C), n_lbl, C) * sd_load
    # 7 the sample loading
    l <- stats::rnorm(n_runs, 0, 0.5)
    # 8 residuals
    E_feat <- matrix(stats::rnorm(n_feat * n_runs), n_feat, n_runs) * sqrt(s2)
    E_lbl <- matrix(stats::rnorm(n_lbl * n_runs), n_lbl, n_runs) * sqrt(s2_lbl)
    E_is <- matrix(stats::rnorm(n_is * n_runs), n_is, n_runs) * 0.25
    # 9 the uniforms of the masks
    U <- matrix(stats::runif(n_feat * n_runs), n_feat, n_runs)
    # 10 the row order per assay (features and labelled controls; standards last)
    feat_assay <- rep(p$assays, each = p$per_assay)
    lbl_assay <- rep(p$assays, each = 5L)[seq_len(n_lbl)]
    perms <- lapply(stats::setNames(p$assays, p$assays), function(a) sample.int(sum(feat_assay == a) + sum(lbl_assay == a)))
    # 11 mask BA: a third of the features with a masked value at VD 20 %, each losing one whole plate or batch
    ba_u <- stats::runif(n_feat)
    ba_block <- sample.int(max(p$blocks, 1L), n_feat, replace = TRUE)

    # the unwanted components of the scenario
    Z <- NULL
    if (scenario %in% c("N1", "E1", "N2", "E2")) {
        basis <- G
        if (scenario %in% c("N2", "E2")) {
            Z <- cbind(run = standardise(as.numeric(level_a)))
            basis <- cbind(basis, Z)
        } else if (p$blocks > 0L) {
            D <- stats::model.matrix(~ factor(block))[, -1L, drop = FALSE]
            D <- sweep(D, 2L, colMeans(D))
            Q <- qr.Q(qr(D)) * sqrt(n_runs)
            colnames(Q) <- paste0("block", seq_len(ncol(Q)))
            Z <- Q
            basis <- cbind(basis, Q)
        }
        z_drift <- residStd(drift_raw, basis)
        Z <- cbind(Z, drift = z_drift)
        basis <- cbind(basis, z_drift)
        if (identical(p$unwanted, "drift_factor") && !scenario %in% c("N2", "E2")) Z <- cbind(Z, factor = residStd(factor_raw, basis))
        stopifnot(ncol(Z) == C)
    }
    unwanted <- function(A) if (is.null(Z) || !nrow(A)) matrix(0, nrow(A), n_runs) else A %*% t(Z)
    L <- matrix(l, n_feat, n_runs, byrow = TRUE)
    group_index <- match(group, groups)
    unw_feat <- unwanted(A_feat)
    unw_lbl <- unwanted(A_lbl)
    oracle_feat <- mu + beta[, group_index, drop = FALSE] + E_feat
    oracle_lbl <- if (n_lbl > 0L) mu_lbl + E_lbl else matrix(0, 0L, n_runs)
    y_feat <- mu + beta[, group_index, drop = FALSE] + unw_feat + L + E_feat
    y_lbl <- if (n_lbl > 0L) mu_lbl + unw_lbl + matrix(l, n_lbl, n_runs, byrow = TRUE) + E_lbl else matrix(0, 0L, n_runs)
    y_is <- if (n_is > 0L) 22 + unwanted(A_is) + matrix(l, n_is, n_runs, byrow = TRUE) + E_is else matrix(0, 0L, n_runs)

    feat_ids <- if (identical(omic, "prot")) sprintf("P%05d", seq_len(n_feat)) else sprintf("F%04d", seq_len(n_feat))
    lbl_ids <- if (n_lbl > 0L) sprintf(if (identical(omic, "lipid")) "LBL_%03d" else "CTL_%03d", seq_len(n_lbl)) else character()
    is_ids <- if (n_is > 0L) sprintf("IS_%02d", seq_len(n_is)) else character()
    is_assay <- rep(p$assays, each = 5L)[seq_len(n_is)]
    dimnames(y_feat) <- dimnames(oracle_feat) <- dimnames(unw_feat) <- dimnames(U) <- list(feat_ids, runs)
    dimnames(y_lbl) <- dimnames(oracle_lbl) <- dimnames(unw_lbl) <- list(lbl_ids, runs)
    dimnames(y_is) <- list(is_ids, runs)
    # the truth per contrast (group g against A)
    contrast_groups <- groups[-1L]
    theta <- vapply(seq_along(contrast_groups), function(i) beta[, i + 1L] - beta[, 1L], numeric(n_feat))
    theta <- rbind(matrix(theta, n_feat, length(contrast_groups)), matrix(0, n_lbl, length(contrast_groups)))
    dimnames(theta) <- list(c(feat_ids, lbl_ids), contrast_groups)
    # the table order per assay
    order_of <- lapply(stats::setNames(p$assays, p$assays), function(a) {
        rows <- c(feat_ids[feat_assay == a], lbl_ids[lbl_assay == a])
        c(rows[perms[[a]]], is_ids[is_assay == a])
    })
    design <- data.frame(Run = runs, group = group, batch = paste0("b", block),
        replicates = stats::ave(seq_len(n_runs), group, FUN = seq_along), tech_rep_group = runs, stringsAsFactors = FALSE)
    list(cell = cellname, scenario = scenario, r = r, p = p, groups = groups, group = group, block = block, runs = runs,
        design = design, feat_ids = feat_ids, lbl_ids = lbl_ids, is_ids = is_ids,
        assay_of = stats::setNames(c(feat_assay, lbl_assay, is_assay), c(feat_ids, lbl_ids, is_ids)),
        y = rbind(y_feat, y_lbl, y_is), oracle = rbind(oracle_feat, oracle_lbl), unwanted = rbind(unw_feat, unw_lbl),
        l = l, theta = theta, planted = stats::setNames(c(rowSums(beta != 0) > 0, rep(FALSE, n_lbl)), c(feat_ids, lbl_ids)),
        U = U, order_of = order_of, ba_u = ba_u, ba_block = ba_block, Z = Z, C = C)
}

# the mask (TRUE = missing) on the feature rows; standards and labelled controls are never missing
solveLod <- function(values, rate) {
    stats::uniroot(function(lod) mean(stats::plogis(lod - values)) - rate, c(min(values) - 60, max(values) + 60), tol = 1e-10)$root
}
maskOf <- function(sim, mask) {
    y <- sim$y[sim$feat_ids, , drop = FALSE]
    if (identical(mask, "none")) return(matrix(FALSE, nrow(y), ncol(y), dimnames = dimnames(y)))
    rate <- if (identical(mask, "ba")) 0.20 else as.numeric(sub("^vd", "", mask)) / 100
    out <- sim$U < stats::plogis(solveLod(y, rate) - y)
    if (identical(mask, "ba")) {
        with_missing <- which(rowSums(out) > 0)
        chosen <- with_missing[order(sim$ba_u[with_missing])][seq_len(round(length(with_missing) / 3))]
        for (j in chosen) out[j, sim$block == sim$ba_block[[j]]] <- TRUE
    }
    out
}
maskedValues <- function(sim, mask, oracle = FALSE) {
    m <- maskOf(sim, mask)
    y <- if (oracle) sim$oracle else sim$y
    y[sim$feat_ids, ][m] <- NA_real_
    y
}

# ---- package objects -------------------------------------------------------------------------------------------------
contrastsOf <- function(sim) {
    g <- sim$groups[-1L]
    data.frame(contrasts = sprintf("group%s-groupA", g), friendly_names = sprintf("%s_vs_A", g),
        full_format = sprintf("%s_vs_A=group%s-groupA", g, g), stringsAsFactors = FALSE)
}
# 35.2 (1), r1.16: a `Name` column holds every row's name and is the annotation column, so the ITSD step's regex mode (the
# default user's path, `^IS_` in the annotation column) finds the standards; the class column and the slot regex are as
# before. `plant` (35.6, L1c only; the plant changes the generator's input, never the route's state): "itsd" builds the
# table as r1.15 did (the class column is the annotation column, so ITSD finds nothing); "itsdpos" gives the standards of
# the first assay (LCMS_Pos in L cells; here the assay named LCMS_Pos) the names STD_01..STD_05 in the `Name` column only.
omicObject <- function(sim, log2_values, ids = NULL, scale = 1, plant = NULL) {
    omic <- sim$p$omic
    tables <- lapply(sim$p$assays, function(a) {
        rows <- sim$order_of[[a]]
        if (!is.null(ids)) rows <- rows[rows %in% ids]
        values <- 2^log2_values[rows, sim$runs, drop = FALSE] * scale
        annotation <- ifelse(startsWith(rows, "IS_"), "ISTD", "PC")
        names_col <- rows
        if (identical(plant, "itsdpos") && identical(a, "LCMS_Pos")) {
            is_rows <- startsWith(rows, "IS_")
            names_col[is_rows] <- sprintf("STD_%02d", seq_len(sum(is_rows)))
        }
        table <- if (identical(omic, "lipid")) data.frame(LipidName = rows, Name = names_col, LipidClass = annotation,
            stringsAsFactors = FALSE) else data.frame(Metabolite = rows, Name = names_col, Ontology = annotation,
            stringsAsFactors = FALSE)
        cbind(table, stats::setNames(as.data.frame(unname(values)), sim$runs))
    })
    names(tables) <- sim$p$assays
    annotation_column <- if (identical(plant, "itsd")) (if (identical(omic, "lipid")) "LipidClass" else "Ontology") else "Name"
    args <- list(deAnalysisParameters = list(formula_string = "~ 0 + group", de_q_val_thresh = 0.05, treat_lfc_cutoff = 0))
    if (identical(omic, "lipid")) {
        methods::new("LipidomicsAssayData", lipid_data = tables, design_matrix = sim$design, lipid_id_column = "LipidName",
            annotation_id_column = annotation_column, sample_id = "Run", group_id = "group", technical_replicate_id = "tech_rep_group",
            database_identifier_type = "Unknown", internal_standard_regex = "^IS_", args = args)
    } else {
        methods::new("MetaboliteAssayData", metabolite_data = tables, design_matrix = sim$design, metabolite_id_column = "Metabolite",
            annotation_id_column = annotation_column, sample_id = "Run", group_id = "group", technical_replicate_id = "tech_rep_group",
            database_identifier_type = "Unknown", internal_standard_regex = "^IS_", args = args)
    }
}
bindingOf <- function(omic) if (identical(omic, "lipid")) lipidomicsDeclarationBinding else metabolomicsDeclarationBinding
stateMatrices <- function(owner, binding, state) {
    object <- owner$getState(state)
    if (is.null(object)) return(NULL)
    id <- methods::slot(object, binding$id_slot)
    runs <- as.character(object@design_matrix[[object@sample_id]])
    lapply(methods::slot(object, binding$data_slot), function(t) {
        t <- as.data.frame(t, check.names = FALSE)
        m <- as.matrix(t[, intersect(runs, names(t)), drop = FALSE])
        rownames(m) <- as.character(t[[id]])
        m
    })
}

# ---- route defaults, read from the configuration and the tab builders at the bound commit (E3) -----------------------
builderValue <- function(fn, input_id, arg = "value") {
    text <- paste(deparse(fn, width.cutoff = 500L), collapse = " ")
    hit <- regmatches(text, regexec(sprintf('ns\\("%s"\\).*?%s = ([^,)]+)', input_id, arg), text, perl = TRUE))[[1L]]
    if (length(hit) < 2L) return(NA_character_)
    trimws(gsub('"', "", hit[[2L]]))
}
routeDefaults <- function(omic) {
    ruv_section <- omicNormRuvSection
    common <- list(auto_min = as.numeric(builderValue(ruv_section, "auto_percentage_min")),
        auto_max = as.numeric(builderValue(ruv_section, "auto_percentage_max")),
        metric = builderValue(ruv_section, "separation_metric", "selected"),
        penalty = as.numeric(builderValue(ruv_section, "k_penalty_weight")),
        adaptive = as.logical(builderValue(ruv_section, "adaptive_k_penalty")),
        correlation = as.numeric(builderValue(omicNormCorrelationFilterTab, "min_pearson_correlation_threshold")))
    if (identical(omic, "prot")) {
        norm <- omicGuiProteomicsNorm()$panels
        qc <- mod_prot_qc_protein_intensity_ui
        return(c(common, list(method = norm$between_sample$default, ruv_mode = norm$ruv_iii$default,
            manual_pct = norm$ruv_iii$percentage, manual_k = norm$ruv_iii$k,
            min_reps = as.numeric(builderValue(qc, "min_reps_per_group")), min_groups = as.numeric(builderValue(qc, "min_groups")),
            percentile = as.numeric(builderValue(qc, "proteins_intensity_cutoff_percentile")),
            strict = as.logical(builderValue(qc, "use_strict_mode")))))
    }
    norm <- (if (identical(omic, "lipid")) omicGuiLipidNorm() else omicGuiMetabNorm())$panels
    c(common, list(method = norm$between_sample$default, offset = norm$log2$offset, itsd = unname(norm$itsd$aggregation[[1L]]),
        ruv_mode = norm$ruv_iii$default, manual_pct = norm$ruv_iii$percentage, manual_k = norm$ruv_iii$k,
        percentile = as.numeric(builderValue(mod_omic_qc_intensity_ui, "intensity_cutoff_percentile")),
        proportion = as.numeric(builderValue(mod_omic_qc_intensity_ui, "proportion_below_cutoff")),
        # r1.20 (KIT-217; 7.4 E3): the per-group count's two settings, from the QC tab's inputs
        min_values_per_group = as.numeric(builderValue(mod_omic_qc_intensity_ui, "min_values_per_group")),
        min_groups = as.numeric(builderValue(mod_omic_qc_intensity_ui, "min_groups"))))
}
EXPECTED_DEFAULTS <- list(
    lipid = list(auto_min = 1, auto_max = 20, metric = "max_difference", penalty = 0.5, adaptive = TRUE, correlation = 0.5,
        method = "cyclicloess", offset = 0, itsd = "median", ruv_mode = "automatic", manual_pct = 5, manual_k = 2,
        percentile = 10, proportion = 0.5, min_values_per_group = 2, min_groups = 1),
    prot = list(auto_min = 1, auto_max = 20, metric = "max_difference", penalty = 0.5, adaptive = TRUE, correlation = 0.5,
        method = "cyclicloess", ruv_mode = "automatic", manual_pct = 5, manual_k = NULL, min_reps = 2, min_groups = 2,
        percentile = 1, strict = FALSE))
EXPECTED_DEFAULTS$metab <- EXPECTED_DEFAULTS$lipid
# 5.2 (r1.20; VAL-241 OD-1 (c)), the import's default pattern as KIT-260 ships it, read at the bound commit by E3 (a difference
# stops: the route under test changed)
BOUND_MARKERS <- c("duplicate_resolution_rank", "correlation_exclusion", # on staging (KIT-255 at aa119308)
    "normalisation_fallback", "protein_normalisation_fallback", "itsd_outcome", "ruv_fallback") # KIT-242 (itsd_outcome: also KIT-260), all on staging at e4a6d4ae
ITSD_DEFAULT_PATTERN <- "^IS_|_d[0-9]+$|ISTD|^IS-|\\(d[0-9]+\\)"
checkE3 <- function(omic) {
    got <- routeDefaults(omic)
    want <- EXPECTED_DEFAULTS[[omic]]
    differ <- names(want)[!vapply(names(want), function(n) identical(got[[n]], want[[n]]) ||
        (is.numeric(want[[n]]) && is.numeric(got[[n]]) && isTRUE(all.equal(got[[n]], want[[n]]))), logical(1L))]
    patterns <- c("(?i)itsd|istd|internal.?standard", "(?i)^rt_|_rt$|retention.?time", "(?i)^is_|_is$|^isd_|_isd$",
        "(?i)_d[0-9]+$|_13c|_15n|deuterated|labeled")
    names_checked <- c(sprintf("LBL_%03d", 1:10), sprintf("CTL_%03d", 1:5))
    matched <- names_checked[vapply(names_checked, function(x) any(vapply(patterns, grepl, logical(1L), x = x, perl = TRUE)), logical(1L))]
    # r1.21 (R20-P06): the cyclic-loess minimum's state is recorded and either state passes; every marker holds its current
    # token (40.3 (5)); 5.2 (r1.20, VAL-241 OD-1 (c)): the import's default ITSD pattern as KIT-260 ships it
    loess_minimum <- tryCatch(omicsCyclicLoessMinimumActive(), error = function(e) NA)
    leaked <- markerLeaks()
    imported <- tryCatch(if (identical(omic, "metab")) getMetabolomicsColumnDefaults("msdial")$is_pattern else
        getLipidomicsColumnDefaults("msdial")$is_pattern, error = function(e) NA_character_)
    pattern_ok <- identical(omic, "prot") || identical(imported, ITSD_DEFAULT_PATTERN)
    got$cyclic_loess_minimum <- loess_minimum
    got$itsd_import_pattern <- imported
    # BIND-1: the replay markers the bound commit must hold, so that a run on the wrong base stops here. All six are on staging at e4a6d4ae
    # (KIT-255: `correlation_exclusion`; the duplicate-resolution rank marker; KIT-242: `normalisation_fallback`, `protein_normalisation_fallback`,
    # `itsd_outcome` (one marker for KIT-242 and KIT-260), `ruv_fallback`; R/utils_workflow_analysis_record.R).
    missing_markers <- setdiff(BOUND_MARKERS, names(workflowReplayRuleMarkers))
    list(pass = !length(differ) && !length(matched) && !length(leaked) && pattern_ok && !length(missing_markers), got = got,
        note = sprintf("differ: %s; labelled names matching an ITSD pattern: %s; markers off their current token: %s; import pattern as shipped %s; replay markers missing at this commit: %s",
            paste(differ, collapse = ",") %||% "", paste(matched, collapse = ","), paste(leaked, collapse = ","), pattern_ok, paste(missing_markers, collapse = ",")))
}

# ---- Part L journeys (5.3) ------------------------------------------------------------------------------------------
REQUESTS_L <- function(d) list(
    LA = list(gate = FALSE, ruv = list(ruv_mode = "automatic", ruv_grouping_variable = "group", percentage_min = d$auto_min,
        percentage_max = d$auto_max, separation_metric = d$metric, k_penalty_weight = d$penalty, adaptive_k_penalty = d$adaptive,
        manual_k = as.integer(d$manual_k), manual_percentage = d$manual_pct)),
    LM = list(gate = FALSE, ruv = list(ruv_mode = "manual", ruv_grouping_variable = "group", percentage_min = d$auto_min,
        percentage_max = d$auto_max, separation_metric = d$metric, k_penalty_weight = d$penalty, adaptive_k_penalty = d$adaptive,
        manual_k = as.integer(d$manual_k), manual_percentage = d$manual_pct)),
    GA = list(gate = TRUE, ruv = list(ruv_design = "groups_all_features", ruv_grouping_variable = "group", ruv_k = NULL,
        ruv_k_rule = "num_sv_be", ruv_k_seed = RUV_FACTOR_COUNT_SEED, ruv_control_ids = NULL)),
    FP = list(gate = TRUE, ruv = list(ruv_design = "first_pass_nulls", ruv_grouping_variable = "group", ruv_k = NULL,
        ruv_k_rule = "sweep_bounded", ruv_k_seed = NULL, ruv_control_ids = NULL)),
    NO = list(gate = FALSE, ruv = NULL),
    NOnone = list(gate = TRUE, ruv = list(ruv_design = "none", ruv_grouping_variable = "group", ruv_k = NULL, ruv_k_rule = NULL,
        ruv_k_seed = NULL, ruv_control_ids = NULL)))

expectedOperations <- function(binding, itsd, ruv, correlation, filter = TRUE) {
    prefix <- binding$prefix
    c(if (filter) paste0(prefix, "_intensity_filter"), paste0(prefix, "_duplicate_resolution"), paste0(prefix, "_qc_finalize"),
        if (itsd) paste0(prefix, "_itsd_normalisation"), paste0(prefix, "_log2"), paste0(prefix, "_between_sample_normalisation"),
        if (ruv) paste0(prefix, "_ruv"), if (correlation) paste0(prefix, "_pearson_selection"))
}

# NOTE (smoke check send 2, defect 2): `$` on a list matches a name by its prefix when no name is exact, so a journey's `da` (absent when the
# DA did not return) was read as `da_reached`, and `e1` could read `e12` or `e13`. The journey's fields are therefore named so that none is
# a prefix of another (`reached_da`, `itsd_check`, `norm_check`), and the new fields below must follow the same rule.
# ---- the per-journey checks E12 to E14 (35.2 (2), (3); 7.4) ------------------------------------------------------------
#' E12 (journeys that run ITSD): per assay (a) no row starting `IS_` remains after the ITSD step, (b) the step's own record
#' count equals the `IS_` rows the assay held at `<prefix>_qc_complete` (0 when the record holds none), (c) at least one
#' non-standard value differs from its `qc_complete` value. `itsd_failed` lists the failing letters.
checkItsdActed <- function(m_qc, m_itsd, record) {
    counts <- record$itsd_counts_per_assay
    do.call(rbind, lapply(names(m_qc), function(a) {
        before <- m_qc[[a]]
        after <- m_itsd[[a]]
        standards <- startsWith(rownames(before), "IS_")
        expected <- sum(standards)
        left <- if (is.null(after)) NA_integer_ else sum(startsWith(rownames(after), "IS_"))
        used <- as.integer((unlist(counts[[a]]) %||% 0L)[[1L]])
        rows <- setdiff(intersect(rownames(before), rownames(after)), rownames(before)[standards])
        cols <- if (is.null(after)) character() else intersect(colnames(before), colnames(after))
        change <- if (length(rows) && length(cols)) suppressWarnings(max(abs(log2(after[rows, cols, drop = FALSE]) -
            log2(before[rows, cols, drop = FALSE])), na.rm = TRUE)) else NA_real_
        if (!is.finite(change)) change <- NA_real_
        failed <- c(if (!isTRUE(left == 0L)) "a", if (!isTRUE(used == expected)) "b", if (!isTRUE(change > 0)) "c")
        data.frame(assay = a, itsd_is_left = left, itsd_n_used = used, itsd_n_expected = expected, itsd_max_change = change,
            itsd_ok = !length(failed), itsd_failed = paste(failed, collapse = ","), stringsAsFactors = FALSE)
    }))
}
#' E13 (every journey with a between-sample method other than `none`, not OR): per assay at least one finite value of the
#' normalised matrix differs from the log2 matrix, values matched by feature id and run.
checkNormActed <- function(m_norm, m_log2) {
    do.call(rbind, lapply(names(m_log2), function(a) {
        x <- m_norm[[a]]
        y <- m_log2[[a]]
        change <- NA_real_
        if (!is.null(x) && !is.null(y)) {
            rows <- intersect(rownames(x), rownames(y))
            cols <- intersect(colnames(x), colnames(y))
            diff <- abs(x[rows, cols, drop = FALSE] - y[rows, cols, drop = FALSE])
            diff <- diff[is.finite(diff)]
            if (length(diff)) change <- max(diff)
        }
        data.frame(assay = a, norm_max_change = change, norm_ok = isTRUE(change > 0), stringsAsFactors = FALSE)
    }))
}
#' The same matrix, rows and columns matched by name, at tolerance 0 (E14; the product's own status uses `identical()`).
sameByName <- function(a, b) {
    if (is.null(a) || is.null(b)) return(FALSE)
    if (!setequal(rownames(a), rownames(b)) || !setequal(colnames(a), colnames(b))) return(FALSE)
    isTRUE(all.equal(unname(a[rownames(a), colnames(a), drop = FALSE]), unname(b[rownames(a), colnames(a), drop = FALSE]),
        tolerance = 0, check.attributes = FALSE))
}
recordOutcomes <- function(record, state) {
    steps <- Filter(function(s) identical(s$state, state), record$steps %||% list())
    if (length(steps)) steps[[1L]]$outcomes else NULL
}
#' r1.21 (R20-P11, 35.2 (3)): the samples the correlation filter left out of one assay only ("assay:sample"), from the
#' matrices before the filter and after it; a sample missing from every assay is not listed.
oneAssayRemovals <- function(before, after) {
    gone <- lapply(stats::setNames(names(before), names(before)), function(a) {
        x <- before[[a]]
        y <- after[[a]]
        present <- colnames(x)[colSums(is.finite(x)) > 0L]
        kept <- if (is.null(y)) character() else colnames(y)[colSums(is.finite(y)) > 0L]
        setdiff(present, kept)
    })
    hits <- character()
    for (a in names(gone)) for (s in gone[[a]]) {
        if (length(gone) < 2L || !all(vapply(setdiff(names(gone), a), function(o) s %in% gone[[o]], logical(1L)))) hits <- c(hits, paste(a, s, sep = ":"))
    }
    if (length(gone) < 2L) character() else hits
}
eqField <- function(x, y) isTRUE(all.equal(as.character(unlist(x)), as.character(unlist(y))))
pickAssay <- function(l, a) {
    v <- l[[a]]
    if (is.null(v) || !length(v)) NA else v[[1L]]
}
#' The per-assay columns of a journey (35.2 (2), (3)): E12 and E13 results, the RUV outcome (E14) and the DA outcome, each
#' recorded whatever happens later. `ruv_status` is `not_reached`, `failed` (unclassed error at step ruv), `refused` (classed),
#' else the factors' status (declared designs) or, for the legacy step, `dropped`, `uncorrected` (by name, tolerance 0) or
#' `corrected`; the product's own `ruv_status` is written beside it.
journeyAssayColumns <- function(out, request, assays) {
    requested <- !is.null(request$ruv)
    declared <- requested && isTRUE(request$gate)
    lapply(stats::setNames(assays, assays), function(a) {
        e12 <- if (!is.null(out$itsd_check)) out$itsd_check[out$itsd_check$assay == a, , drop = FALSE]
        e13 <- if (!is.null(out$norm_check)) out$norm_check[out$norm_check$assay == a, , drop = FALSE]
        field <- function(x, column) if (!is.null(x) && nrow(x)) x[[column]][[1L]] else NA
        trace <- out$optimisation[[a]]
        res <- trace$optimization_results
        product <- out$ruv_value
        status <- if (!requested) NA_character_ else if (!isTRUE(out$ruv_reached)) "not_reached" else
            if (identical(out$error_step, "ruv")) (if (isTRUE(out$error_refused)) "refused" else "failed") else
            if (declared) as.character(out$ruv_factors[[a]]$status %||% "no_status") else
            if (is.null(out$m_out) || !a %in% names(out$m_out)) "dropped" else
            if (sameByName(out$m_in[[a]], out$m_out[[a]])) "uncorrected" else "corrected"
        product_status <- if (requested && !declared) as.character(pickAssay(product$ruv_status, a)) else NA_character_
        reason <- NA_character_
        if (requested && !declared && !is.null(trace) && !isTRUE(trace$success) && is.data.frame(res) && nrow(res)) {
            counts <- table(res$status)
            reason <- paste(sprintf("%s:%d", names(counts), as.integer(counts)), collapse = ";")
            first <- res$error_reason[nzchar(as.character(res$error_reason))]
            if (length(first)) reason <- paste0(reason, " | ", substr(first[[1L]], 1L, 300L))
            if (any(res$status == "skipped_k_exceeds_complete_controls")) reason <- paste0(reason, " | ", oneLine(trace$error))
            reason <- oneLine(reason, 600L)
        }
        record <- out$record_ruv
        legacy_return <- requested && !declared && status %in% c("corrected", "uncorrected", "dropped")
        da_failed <- identical(out$error_step, "da")
        failed_here <- if (!is.null(out$da_errors)) out$da_errors[out$da_errors$assay == a, , drop = FALSE]
        da_status <- if (da_failed) "failed" else if (!isTRUE(out$reached_da)) "not_reached" else
            if (is.null(out$da) || !nrow(out$da) || !a %in% out$da$assay || (!is.null(failed_here) && nrow(failed_here))) "no_da" else "ok"
        list(itsd_is_left = field(e12, "itsd_is_left"), itsd_n_used = field(e12, "itsd_n_used"), itsd_n_expected = field(e12, "itsd_n_expected"),
            itsd_max_change = field(e12, "itsd_max_change"), itsd_ok = field(e12, "itsd_ok"),
            norm_max_change = field(e13, "norm_max_change"), norm_ok = field(e13, "norm_ok"),
            ruv_status = status, ruv_status_product = product_status,
            ruv_status_agree = if (legacy_return && !is.na(product_status)) identical(status, product_status) else NA,
            ruv_changed = if (legacy_return) identical(status, "corrected") else NA,
            ruv_k_found = if (requested && !declared && !is.null(trace)) (if (isTRUE(trace$success)) as.integer(trace$best_k) else NA_integer_) else NA_integer_,
            ruv_n_controls = if (requested && !declared && !is.null(trace)) sum(trace$control_genes_index, na.rm = TRUE) else NA_integer_,
            ruv_no_k_reason = reason,
            ruv_k_sweep = if (requested && !declared) pickAssay(product$k_sweep, a) else NA,
            ruv_k_used = if (requested && !declared) pickAssay(product$best_k, a) else NA,
            ruv_k_cap = if (requested && !declared) pickAssay(product$k_cap, a) else NA,
            ruv_k_samples = if (requested && !declared) pickAssay(product$k_samples, a) else NA,
            ruv_k_notice = if (requested && !declared) oneLine(pickAssay(product$k_notice, a) %|NA|% "") else NA_character_,
            ruv_record_applied = if (!is.null(record)) as.logical(record$ruv_applied %||% NA) else NA,
            ruv_record_match = if (!is.null(record) && requested && !declared) all(eqField(record$k_sweep_per_assay[[a]], product$k_sweep[[a]]),
                eqField(record$k_cap_per_assay[[a]], product$k_cap[[a]]), eqField(record$k_per_assay[[a]], product$best_k[[a]]),
                eqField(record$ruv_status_per_assay[[a]], product$ruv_status[[a]])) else NA,
            da_status = da_status, da_class = if (da_failed) out$error_class else NA_character_,
            da_message = if (da_failed) out$error_message else if (!is.null(failed_here) && nrow(failed_here)) failed_here$error[[1L]] else NA_character_,
            warning_classes = out$warning_classes %||% NA_character_, corr_removed_one_assay = out$corr_one_assay %||% NA_character_,
            itsd_record_fields = out$itsd_record_fields %||% NA_character_)
    })
}
`%|NA|%` <- function(x, y) if (is.null(x) || (length(x) == 1L && is.na(x))) y else x

#' One journey on its own session owner. Returns the statuses, the RUV choice, the matrices the checks need and the DA.
#' `plant` (35.6, L1c only): "itsd" and "itsdpos" change the table the route is given (omicObject); "normplant" changes the
#' matrix E13 compares; "daerror" stops at step da. A plant never changes the route's own state.
journeyL <- function(sim, log2_values, arm, request, d, root, filter = TRUE, itsd = isTRUE(sim$itsd) || identical(sim$p$omic, "lipid"),
    method = d$method, correlation = TRUE, ids = NULL, scale = 1, keep_states = FALSE, fault = FALSE, plant = NULL) {
    t0 <- proc.time()[[3L]]
    warning_count$n <- 0L
    warning_count$log <- character()
    warning_count$assays <- sim$p$assays
    binding <- bindingOf(sim$p$omic)
    prefix <- binding$prefix
    # 40.3 (5): the arm's markers, asserted before anything is built (a leak stops the unit, not the journey)
    markers <- armMarkers(arm)
    leaks <- markerLeaks(markers)
    if (length(leaks)) stop(sprintf("K188 marker leak at the start of a %s journey: %s", arm, paste(leaks, collapse = "; ")))
    out <- list(arm = arm, status = "ok", class = NA_character_, step = NA_character_)
    old <- options(c(list(MultiScholaR.ruv_design_panel = isTRUE(request$gate)), markers))
    on.exit(options(old), add = TRUE)
    stopifnot(all(vapply(names(markers), function(o) identical(getOption(o), markers[[o]]), logical(1L))))
    fixture <- NULL
    step <- "owner"
    omic_plant <- if (!is.null(plant) && plant %in% c("itsd", "itsdpos")) plant else NULL
    result <- tryCatch(withCallingHandlers({
        object <- omicObject(sim, log2_values, ids = ids, scale = scale, plant = omic_plant)
        fixture <- quiet(if (identical(sim$p$omic, "lipid")) lipidDeclarationsOwner(root, object) else metabDeclarationsOwner(root, object))
        fixture$workflow$contrasts_tbl <- contrastsOf(sim)
        owner <- fixture$owner
        config <- fixture$workflow$config_list
        if (fault) stop("K188 deliberate fault (E7)")
        # the planted guard (journey level): a deliberate `$` prefix match inside the journey; its stop must end the journey as `instrument`
        if (identical(plant, "guardjourney")) {
            probe <- list(partial_probe_value = 1)
            probe$partial_probe
        }
        if (filter) {
            step <- "intensity_filter"
            # 40.3 (3): KIT-217's two settings are passed from E3's read values, not left to the function's defaults
            quiet(runWorkflowOmicsReferenceIntensityFilter(binding, owner, d$percentile, d$proportion, config,
                min_values_per_group = as.integer(d$min_values_per_group), min_groups = as.integer(d$min_groups)))
        }
        step <- "duplicate_resolution"
        quiet(runWorkflowOmicsReferenceDuplicateResolution(binding, owner, config))
        step <- "qc_finalize"
        quiet(commitWorkflowOmicsStoreStep(owner, paste0(prefix, "_qc_finalize"), paste0(prefix, "_qc_complete"), config, "QC complete"))
        m_qc <- stateMatrices(owner, binding, paste0(prefix, "_qc_complete"))
        out$kept <- unlist(lapply(m_qc, rownames), use.names = FALSE)
        out$rows_filter <- vapply(m_qc, nrow, integer(1L))
        if (itsd) {
            step <- "itsd"
            quiet(runWorkflowOmicsReferenceItsdNormalisation(binding, owner, d$itsd, NULL, config))
            out$m_itsd <- stateMatrices(owner, binding, paste0(prefix, "_itsd_norm"))
            out$rows_itsd <- vapply(out$m_itsd, nrow, integer(1L))
            itsd_record <- owner$getState(paste0(prefix, "_itsd_norm"))@args$ITSDNormalization
            out$itsd_record_fields <- paste(names(itsd_record), collapse = ",")
            out$itsd_features_record <- itsd_record[["itsd_features_per_assay"]]
            out$itsd_check <- checkItsdActed(m_qc, out$m_itsd, itsd_record)
            if (!all(out$itsd_check$itsd_ok)) {
                bad <- out$itsd_check[!out$itsd_check$itsd_ok, , drop = FALSE]
                out$status <- "instrument"
                out$class <- paste(sprintf("itsd_not_applied:%s:%s", bad$assay, bad$itsd_failed), collapse = " | ")
            }
        }
        step <- "log2"
        quiet(runWorkflowOmicsReferenceLog2(binding, owner, d$offset, config))
        out$m_log2 <- stateMatrices(owner, binding, paste0(prefix, "_log2"))
        step <- "between_sample"
        quiet(runWorkflowOmicsReferenceBetweenSampleNormalisation(binding, owner, method, config))
        out$m_in <- stateMatrices(owner, binding, paste0(prefix, "_normalized"))
        out$rows_ruv_input <- vapply(out$m_in, nrow, integer(1L))
        if (!identical(method, "none") && !identical(arm, "OR")) {
            compared <- out$m_in
            # `normplant` (35.6): the matrix E13 compares is the log2 matrix for one assay; the route's output is untouched
            if (identical(plant, "normplant")) compared[["LCMS_Pos"]] <- out$m_log2[["LCMS_Pos"]]
            out$norm_check <- checkNormActed(compared, out$m_log2)
            if (!all(out$norm_check$norm_ok)) {
                bad <- out$norm_check[!out$norm_check$norm_ok, , drop = FALSE]
                text <- paste(sprintf("norm_not_applied:%s", bad$assay), collapse = " | ")
                out$class <- if (identical(out$status, "instrument")) paste(out$class, text, sep = " | ") else text
                out$status <- "instrument"
            }
        }
        if (!is.null(request$ruv)) {
            step <- "ruv"
            out$ruv_reached <- TRUE
            run <- quiet(runWorkflowOmicsReferenceRuv(binding, owner, request$ruv, config))
            out$ruv_value <- run$value[setdiff(names(run$value), c("final", "optimisation", "factors", "separation_plots"))]
            out$ruv_factors <- lapply(run$value$factors, function(f) f[intersect(names(f), c("status", "reason", "k", "control_ids",
                "control_count", "ruv_k_estimate"))])
            out$optimisation <- lapply(run$optimisation, function(o) o[intersect(names(o), c("success", "best_k", "best_percentage",
                "control_genes_index", "optimization_results", "error"))])
            out$m_out <- stateMatrices(owner, binding, paste0(prefix, "_ruv_corrected"))
            # 35.2 (3): the analysis record's fields beside the product's (they must equal it)
            record <- tryCatch(quiet(buildWorkflowAnalysisRecord(owner)), error = function(e) NULL)
            out$record_ruv <- recordOutcomes(record, paste0(prefix, "_ruv_corrected"))
            out$record_itsd <- recordOutcomes(record, paste0(prefix, "_itsd_norm"))
            out$record_norm <- recordOutcomes(record, paste0(prefix, "_normalized"))
        } else if (!is.null(plant)) {
            # the planted cases without a RUV request (`emptysample`) read the record's ITSD and normalisation outcomes too (bound commit, KIT-242)
            record <- tryCatch(quiet(buildWorkflowAnalysisRecord(owner)), error = function(e) NULL)
            out$record_itsd <- recordOutcomes(record, paste0(prefix, "_itsd_norm"))
            out$record_norm <- recordOutcomes(record, paste0(prefix, "_normalized"))
        }
        if (correlation) {
            step <- "correlation"
            quiet(runWorkflowOmicsReferenceCorrelation(binding, owner, "group", d$correlation, config))
            out$corr_one_assay <- paste(oneAssayRemovals(out$m_out %||% out$m_in, stateMatrices(owner, binding, owner$current_state)),
                collapse = ";")
        }
        final <- owner$getState(owner$current_state)
        out$runs_final <- as.character(final@design_matrix[[final@sample_id]])
        step <- "da"
        out$reached_da <- TRUE
        if (identical(plant, "daerror")) stop("K188 planted DA error (E15)")
        da <- quiet(runWorkflowOmicsReferenceDaAnalysis(binding, owner, getWorkflowStateReference(owner),
            list(formula_string = "~ 0 + group", da_q_val_thresh = 0.05, treat_lfc_cutoff = 0), contrastsOf(sim)))
        long <- as.data.frame(da$results[[binding$long_table]])
        id <- intersect(c("lipid_id", "metabolite_id"), names(long))[[1L]]
        column <- function(names) {
            hit <- names[names %in% names(long)]
            if (length(hit)) as.numeric(long[[hit[[1L]]]]) else rep(NA_real_, nrow(long))
        }
        out$da <- data.frame(assay = as.character(long$assay), comparison = as.character(long$comparison),
            feature = as.character(long[[id]]), logFC = column("logFC"), t = column("t"), p = column(c("P.Value", "raw_pvalue")),
            q = column("fdr_qvalue"), bh = column("fdr_value_bh"), stringsAsFactors = FALSE)
        # 5.4 (r1.9): a contrast's per-assay DA failure, from the manifest's errors channel (KIT-128), is `no_da` for that family
        errs <- lapply(da$manifest$results, function(e) if (length(e$errors)) data.frame(group = groupOfComparison(e$comparison %||% e$contrast %||% ""),
            assay = names(e$errors), error = vapply(e$errors, function(x) oneLine(paste(as.character(unlist(x)), collapse = " ")), ""),
            stringsAsFactors = FALSE))
        out$da_errors <- do.call(rbind, errs)
        # E7 (no_da): a deliberate per-assay failure in the first journey of the fault unit that reaches the DA with two assays
        if (isTRUE(getOption("k188.fault_noda", FALSE)) && length(unique(out$da$assay)) >= 2L) {
            options(k188.fault_noda = FALSE)
            a2 <- sort(unique(out$da$assay))[[2L]]
            g1 <- groupOfComparison(out$da$comparison)[[1L]]
            out$da <- out$da[!(out$da$assay == a2 & groupOfComparison(out$da$comparison) == g1), , drop = FALSE]
            out$da_errors <- rbind(out$da_errors, data.frame(group = g1, assay = a2, error = "K188 deliberate per-assay DA failure (E7)",
                stringsAsFactors = FALSE))
        }
        # E1: the provenance of every committed state names the declared operation and the binding's default variant
        snapshot <- owner$.context$snapshot
        history <- unlist(snapshot$state_history)
        provenance <- snapshot$step_provenance[history]
        operations <- vapply(provenance, function(x) as.character(x$operation %||% NA), "")
        variants <- vapply(provenance, function(x) as.character(x$variant %||% NA), "")
        after_design <- which(operations == "design_builder")
        steps_run <- if (length(after_design)) seq.int(after_design[[1L]] + 1L, length(operations)) else integer()
        steps_run <- steps_run[steps_run <= length(operations)]
        want <- expectedOperations(binding, itsd, !is.null(request$ruv), correlation, filter)
        want_variant <- vapply(want, function(op) as.character(binding$default_variants[op] %||% NA), "")
        want_variant[is.na(want_variant)] <- "whole"
        got_ops <- operations[steps_run]
        got_variant <- variants[steps_run]
        out$e1 <- list(pass = identical(unname(got_ops), unname(want)) && identical(unname(got_variant), unname(want_variant)),
            note = oneLine(sprintf("ops %s; variants %s", paste(got_ops, collapse = ","), paste(got_variant, collapse = ","))))
        NULL
    }, warning = guardWarning), error = function(e) {
        cls <- intersect(class(e), REFUSALS)
        text <- if (length(cls)) cls[[1L]] else oneLine(paste(class(e)[[1L]], conditionMessage(e), sep = ": "))
        # a partial-match guard stop is an INSTRUMENT defect (VAL-54), never a route failure: it must not enter the 5 % gate
        if (grepl("partial-match guard", conditionMessage(e), fixed = TRUE)) {
            out$status <<- "instrument"
            out$class <<- oneLine(paste("partial_match_guard:", conditionMessage(e)))
        } else if (identical(out$status, "instrument")) {
            # r1.17 (R16-P07): a journey already `instrument` keeps its status; the later error goes to `class`
            out$class <<- paste(out$class, text, sep = " | ")
        } else {
            out$status <<- if (length(cls)) "refused" else "crashed"
            out$class <<- text
        }
        out$step <<- step
        out$error_step <<- step
        out$error_refused <<- length(cls) > 0L
        out$error_class <<- class(e)[[1L]]
        out$error_message <<- oneLine(conditionMessage(e))
        NULL
    }, finally = {
        if (!is.null(fixture)) try(closeDIAStateAdapter(fixture$owner), silent = TRUE)
        unlink(root, recursive = TRUE)
    })
    out$warning_classes <- warningClasses()
    out$assay_cols <- journeyAssayColumns(out, request, sim$p$assays)
    # what the legacy RUV step returned (2.7, E14): the status comparison runs whatever happens later; the journey's own
    # status keeps its precedence (a crash, a refusal or an instrument failure stands)
    if (identical(out$status, "ok") && !is.null(request$ruv) && !isTRUE(request$gate) && !is.null(out$m_in)) {
        statuses <- vapply(out$assay_cols, function(x) x$ruv_status, "")
        if (any(statuses == "dropped")) out$status <- "dropped"
        if (all(statuses == "uncorrected")) out$status <- "uncorrected"
    }
    if (identical(out$status, "ok") && (is.null(out$da) || !nrow(out$da))) out$status <- "no_da"
    out$warnings <- warning_count$n
    out$seconds <- proc.time()[[3L]] - t0
    out
}

#' The legacy choice per assay: the k used (the RUV product's `best_k`, which is what the route removed and what Part M,
#' LCm, E8, M-amp and M-df take), the sweep's k (the optimisation trace's `best_k`, the uncapped one; r1.21, R20-P04), the
#' winning percentage, the control ids (2.2).
legacyChoice <- function(journey) {
    o <- journey$optimisation
    if (is.null(o)) return(NULL)
    lapply(stats::setNames(names(o), names(o)), function(a) {
        x <- o[[a]]
        ctl <- x$control_genes_index
        k_sweep <- if (isTRUE(x$success)) as.integer(x$best_k) else NA_integer_
        k_used <- suppressWarnings(as.integer(unlist(journey[["ruv_value"]][["best_k"]][[a]])[1L]))
        list(k = if (length(k_used) && !is.na(k_used)) k_used else k_sweep, k_sweep = k_sweep, pct = x$best_percentage %||% NA_real_,
            controls = if (is.logical(ctl)) names(ctl)[which(ctl)] else character(), trace = x$optimization_results)
    })
}
declaredChoice <- function(journey) {
    f <- journey$ruv_factors
    if (is.null(f)) return(NULL)
    lapply(f, function(x) list(k = as.integer(x$k %||% NA_integer_), pct = NA_real_, controls = as.character(x$control_ids %||% character()),
        status = x$status %||% NA_character_, estimate = x$ruv_k_estimate %||% NA_integer_))
}

# ---- Part P journeys (18.4; the protein-level LFQ route as kit-177/route.R:82-118 runs it) -------------------------
pConfig <- function() {
    path <- system.file("config", "config.ini", package = "MultiScholaR")
    utils::capture.output(config <- suppressMessages(readConfigFile(file = path, create_workers = FALSE)))
    config
}
REQUESTS_P <- function(d) {
    base <- list(norm_method = d$method, auto_percentage_min = d$auto_min, auto_percentage_max = d$auto_max,
        separation_metric = d$metric, k_penalty_weight = d$penalty, adaptive_k_penalty = d$adaptive,
        ruv_percentage = d$manual_pct, ruv_k = d$manual_k, grouping = "group")
    list(LA = list(gate = FALSE, norm = c(base, list(ruv_mode = "automatic"))),
        LM = list(gate = FALSE, norm = c(base, list(ruv_mode = "manual"))),
        GA = list(gate = TRUE, norm = c(base, list(ruv_mode = "automatic", ruv_design = "groups_all_features", ruv_k_rule = "num_sv_be",
            ruv_k_seed = NULL))),
        NO = list(gate = FALSE, norm = c(base, list(ruv_mode = "skip"))),
        NOnone = list(gate = TRUE, norm = c(base, list(ruv_mode = "skip", ruv_design = "none", ruv_k_rule = NULL, ruv_k_seed = NULL))))
}
protMatrixOf <- function(object) {
    q <- as.data.frame(object@protein_quant_table, check.names = FALSE)
    runs <- as.character(object@design_matrix[[object@sample_id]])
    m <- as.matrix(q[, intersect(runs, names(q)), drop = FALSE])
    rownames(m) <- as.character(q[[object@protein_id_column]])
    m
}
journeyP <- function(sim, log2_values, arm, request, d, root, config, filter = TRUE, correlation = TRUE, ids = NULL, fault = FALSE) {
    t0 <- proc.time()[[3L]]
    warning_count$n <- 0L
    out <- list(arm = arm, status = "ok", class = NA_character_, step = NA_character_)
    old <- options(MultiScholaR.ruv_design_panel = isTRUE(request$gate))
    on.exit(options(old), add = TRUE)
    manager <- NULL
    step <- "owner"
    tryCatch({
        dir.create(root, recursive = TRUE, showWarnings = FALSE)
        manager <- newDIAStateAdapter(newWorkflowEntityStore(root, paste0("k188-", basename(root)), "LFQ"),
            WorkflowState$new(audit_enabled = FALSE))
        manager$.context$reference_states <- TRUE
        manager$setWorkflowType("LFQ")
        paths <- list(source_dir = file.path(root, "src"), protein_qc_dir = file.path(root, "pqc"), da_output_dir = file.path(root, "da"),
            publication_graphs_dir = file.path(root, "pub"), time_dir = file.path(root, "time"), results_dir = file.path(root, "res"))
        for (pth in paths) dir.create(pth, recursive = TRUE, showWarnings = FALSE)
        manager$.context$experiment_paths <- paths
        if (fault) stop("K188 deliberate fault (E7)")
        rows <- sim$order_of[[1L]]
        if (!is.null(ids)) rows <- rows[rows %in% ids]
        values <- 2^log2_values[rows, sim$runs, drop = FALSE]
        # 18.4 (r1.10): the bound commit's own design-step transform (KIT-177: `protDesignPlainLog2`, plain log2), else the
        # base's log2(x + 1)
        raw <- tibble::as_tibble(cbind(data.frame(Protein.Ids = rows, stringsAsFactors = FALSE),
            stats::setNames(as.data.frame(unname(values)), sim$runs)))
        if (exists("protDesignPlainLog2", mode = "function")) {
            made <- protDesignPlainLog2(raw, "Protein.Ids")
            wide <- tibble::as_tibble(made$table)
            fixed <- c(list(pivot = "protein x run", transform = "log2(x); 0 missing"), as.list(made$counts))
        } else {
            wide <- raw
            for (run in sim$runs) wide[[run]] <- log2(wide[[run]] + 1)
            fixed <- list(pivot = "protein x run", transform = "log2(x + 1)")
        }
        out$transform <- fixed$transform
        design <- sim$design[, c("Run", "group", "batch", "replicates")]
        design$tech_rep_group <- paste(design$group, design$replicates, sep = "_")
        s4 <- ProteinQuantitativeData(protein_quant_table = wide, protein_id_column = "Protein.Ids",
            protein_id_table = dplyr::distinct(wide, Protein.Ids), design_matrix = design, sample_id = "Run", group_id = "group",
            technical_replicate_id = "tech_rep_group", args = config)
        setDIAStepProvenance(manager, "design_builder", fixed = fixed)
        manager$saveState("protein_s4_initial", s4, config, "Initial Protein S4 object created in the KIT-188 script")
        workflow <- list2env(list(state_manager = manager, qc_params = list(protein_qc = list()), config_list = list()))
        contrasts <- data.frame(contrasts = sprintf("%s.vs.A", sim$groups[-1L]),
            full_format = sprintf("%s.vs.A=group%s-groupA", sim$groups[-1L], sim$groups[-1L]),
            friendly_names = sprintf("%s.vs.A", sim$groups[-1L]), stringsAsFactors = FALSE)
        withDIAStateScope(manager, function() {
            finish <- function(value) finishDIAReferenceProteinStep(workflow, value, publish = FALSE)
            step <<- "accession_cleanup"
            quiet(finish(computeDIAReferenceAccessionCleanup(manager, getDIAStateReference(manager), list(delimiter = ";",
                aggregation_method = "mean"))))
            if (filter) {
                step <<- "intensity_filter"
                quiet(finish(computeDIAReferenceProteinIntensityFilter(manager, getDIAStateReference(manager), list(strict_mode = d$strict,
                    min_reps_per_group = d$min_reps, min_groups = d$min_groups, proteins_intensity_cutoff_percentile = d$percentile))))
            }
            step <<- "duplicate_removal"
            quiet(finish(computeDIAReferenceProteinDuplicateRemoval(manager, getDIAStateReference(manager), list(aggregation_method = "mean"))))
            step <<- "replicate_filter"
            quiet(completeDIAProteinAction(workflow, computeDIAReferenceProteinReplicateFilter(manager, getDIAStateReference(manager),
                list(grouping_variable = "group", parallel_cores = 1L), variant = "rows"), publish = FALSE))
            before <- manager$getState(manager$current_state)
            out$kept <<- as.character(before@protein_quant_table$Protein.Ids)
            step <<- "normalisation"
            run <- quiet(computeDIAReferenceNormalization(manager, getDIAStateReference(manager), request$norm))
            commit <- diaNormalizationCommit(manager, request$norm, run)
            if (!is.null(commit)) commitDIAReferenceState(manager, commit$state_name, commit$checkpoint, commit$config, commit$description)
            saved <- run$value
            out$ruv_value <<- list(best_k = saved$best_k, percentage = saved$percentage, controls = saved$controls,
                ruv_applied = saved$ruv_result$ruv_applied %||% NA)
            # the RUV input, made in the script as the kernel makes it (E8, E9)
            if (!identical(request$norm$ruv_mode, "skip")) {
                out$m_in <<- list(proteins = protMatrixOf(quiet(normaliseBetweenSamples(before, normalisation_method = request$norm$norm_method))))
                out$m_out <<- list(proteins = protMatrixOf(manager$getState(manager$current_state)))
            }
            step <<- "correlation"
            if (correlation) {
                correlated <- quiet(computeDIAReferenceCorrelation(manager, saved$final, "group", d$correlation))
                commit <- diaCorrelationCommit(manager, correlated$final, d$correlation, FALSE)
            } else commit <- diaCorrelationCommit(manager, saved$final, 0, TRUE)
            if (!is.null(commit)) commitDIAReferenceState(manager, commit$state_name, commit$checkpoint, commit$config, commit$description)
            final <- manager$getState(manager$current_state)
            out$runs_final <<- as.character(final@design_matrix[[final@sample_id]])
            step <<- "da"
            da_data <- new.env(parent = emptyenv())
            da_data$current_s4_object <- getDIAStateReference(manager)
            da_data$contrasts_available <- contrasts$contrasts
            workflow$contrasts_tbl <- contrasts
            manifest <- quiet(runDIAReferenceDaAnalysis(list(formula_string = "~ 0 + group", da_q_val_thresh = 0.05, treat_lfc_cutoff = 0),
                da_data, workflow, paths, contrasts = contrasts, publish = FALSE))
            tables <- lapply(names(manifest$results), function(name) {
                outputs <- readDIAContrastOutput(manifest$results[[name]], c("da_proteins_long", "contrasts_results_table"))
                pick <- function(table, names) {
                    hit <- names[names %in% names(table)]
                    if (length(hit)) table[[hit[[1L]]]] else NULL
                }
                da <- as.data.frame(outputs$da_proteins_long)
                res <- as.data.frame(outputs$contrasts_results_table %||% data.frame())
                feature <- as.character(pick(da, c("Protein.Ids", "uniprot_acc")))
                t <- pick(da, c("t", "t.mod"))
                if (is.null(t) && nrow(res)) {
                    rid <- as.character(pick(res, c("Protein.Ids", "uniprot_acc")))
                    tv <- pick(res, c("t", "t.mod"))
                    if (!is.null(tv) && length(rid)) t <- as.numeric(tv)[match(feature, rid)]
                }
                q <- as.numeric(pick(da, "fdr_qvalue"))
                p <- as.numeric(pick(da, c("raw_pvalue", "P.Value")))
                data.frame(assay = "proteins", comparison = name, feature = feature, logFC = as.numeric(pick(da, c("log2FC", "logFC"))),
                    t = if (is.null(t)) NA_real_ else as.numeric(t), p = p, q = if (length(q)) q else NA_real_,
                    bh = stats::p.adjust(p, "BH"), stringsAsFactors = FALSE)
            })
            out$da <<- do.call(rbind, tables)
            history <- as.character(unlist(manager$.context$snapshot$state_history))
            history <- history[seq_along(history) > match("protein_s4_initial", history, nomatch = 0L)]
            ops <- vapply(manager$.context$snapshot$step_provenance[history], function(x) as.character(x$operation %||% NA), "")
            want <- c("normalization-ruv", if (correlation) "pearson-selection")
            out$e1 <<- list(pass = !any(ops == "legacy", na.rm = TRUE) && !anyNA(ops) && all(want %in% ops),
                note = oneLine(sprintf("ops %s", paste(ops, collapse = ","))))
        })
        NULL
    }, error = function(e) {
        cls <- intersect(class(e), REFUSALS)
        out$status <<- if (length(cls)) "refused" else "crashed"
        out$class <<- if (length(cls)) cls[[1L]] else oneLine(paste(class(e)[[1L]], conditionMessage(e), sep = ": "))
        out$step <<- step
        NULL
    }, finally = {
        if (!is.null(manager)) try(closeDIAStateAdapter(manager), silent = TRUE)
        unlink(root, recursive = TRUE)
    })
    if (identical(out$status, "ok") && (is.null(out$da) || !nrow(out$da))) out$status <- "no_da"
    out$warnings <- warning_count$n
    out$seconds <- proc.time()[[3L]] - t0
    out
}

# ---- estimands (5.5), one family = (assay, contrast) ---------------------------------------------------------------
groupOfComparison <- function(x) {
    x <- as.character(x)
    full <- regmatches(x, regexec("group([A-E])-groupA", x))
    short <- regmatches(x, regexec("^([A-E])(_vs_A|\\.vs\\.A)", x))
    vapply(seq_along(x), function(i) if (length(full[[i]]) >= 2L) full[[i]][[2L]] else if (length(short[[i]]) >= 2L) short[[i]][[2L]] else
        NA_character_, "")
}
#' The estimands of one family. `da` the arm's rows of the family (feature, logFC, p, q); `T_ids` the scored set;
#' `theta` the truth on T_ids; `or_err` the oracle's |error| by feature (NULL for the oracle itself).
scoreFamily <- function(da, T_ids, theta, or_err = NULL) {
    da <- da[da$feature %in% T_ids & is.finite(da$p), , drop = FALSE]
    da <- da[!duplicated(da$feature), , drop = FALSE]
    th <- theta[da$feature]
    called <- is.finite(da$q) & da$q < 0.05
    bh <- stats::p.adjust(da$p, "BH")
    bh_called <- bh < 0.05
    null <- th == 0
    planted_T <- sum(theta[T_ids] != 0)
    err <- abs(da$logFC - th)
    err[!is.finite(err)] <- NA_real_
    share <- function(a, b) if (b > 0) a / b else 0
    out <- list(n_T = length(T_ids), n_tested = nrow(da), lost = length(T_ids) - nrow(da),
        n_null = sum(null), null_p05 = if (any(null)) mean(da$p[null] < 0.05) else NA_real_,
        null_p01 = if (any(null)) mean(da$p[null] < 0.01) else NA_real_,
        calls = sum(called), false_calls = sum(called & null), true_calls = sum(called & !null), any_call = as.integer(any(called)),
        fdp = share(sum(called & null), sum(called)), n_planted = planted_T,
        power = if (planted_T > 0) sum(called & !null) / planted_T else NA_real_,
        bh_calls = sum(bh_called), bh_false = sum(bh_called & null), bh_any = as.integer(any(bh_called)),
        bh_fdp = share(sum(bh_called & null), sum(bh_called)), bh_power = if (planted_T > 0) sum(bh_called & !null) / planted_T else NA_real_,
        err_med = if (any(!is.na(err))) stats::median(err, na.rm = TRUE) else NA_real_,
        err_p99 = if (any(!is.na(err))) unname(stats::quantile(err, 0.99, na.rm = TRUE)) else NA_real_,
        err_rmse = if (any(!is.na(err))) sqrt(mean(err^2, na.rm = TRUE)) else NA_real_,
        err_max = if (any(!is.na(err))) max(err, na.rm = TRUE) else NA_real_,
        # M-warn's fourth signal (5.8): the largest |logFC| the family returned
        logfc_max = if (any(is.finite(da$logFC))) max(abs(da$logFC), na.rm = TRUE) else NA_real_,
        gross1 = if (any(!is.na(err))) mean(err > 1, na.rm = TRUE) else NA_real_,
        slope = if (any(!null)) sum(th[!null] * da$logFC[!null], na.rm = TRUE) / sum(th[!null]^2) else NA_real_,
        q_method = if (!nrow(da) || !any(is.finite(da$q))) NA_character_ else if (isTRUE(all.equal(da$q, bh, check.attributes = FALSE)))
            "BH_or_pi0_1" else "qvalue",
        n_common = NA_integer_, s5 = NA_integer_, or_s5 = NA_integer_, gross1_common = NA_real_, or_gross1_common = NA_real_)
    if (is.null(or_err)) {
        out$or_s5 <- as.integer(any(err > 5, na.rm = TRUE))
    } else {
        common <- intersect(da$feature[!is.na(err)], names(or_err)[!is.na(or_err)])
        e_arm <- err[match(common, da$feature)]
        e_or <- or_err[common]
        out$n_common <- length(common)
        out$s5 <- as.integer(any(e_arm > 5 & e_or <= 5))
        out$gross1_common <- if (length(common)) mean(e_arm > 1) else NA_real_
        out$or_gross1_common <- if (length(common)) mean(e_or > 1) else NA_real_
    }
    out$err_by_feature <- stats::setNames(err, da$feature)
    out
}
ASSAY_COLUMNS <- c("itsd_is_left", "itsd_n_used", "itsd_n_expected", "itsd_max_change", "itsd_ok", "norm_max_change", "norm_ok",
    "ruv_status", "ruv_status_product", "ruv_status_agree", "ruv_changed", "ruv_k_found", "ruv_n_controls", "ruv_no_k_reason",
    "ruv_k_sweep", "ruv_k_used", "ruv_k_cap", "ruv_k_samples", "ruv_k_notice", "ruv_record_applied", "ruv_record_match",
    "da_status", "da_class", "da_message", "warning_classes", "corr_removed_one_assay", "itsd_record_fields")
ROW_COLUMNS <- c("part", "cell", "job", "scenario", "mask", "rep", "arm", "arm_scale", "scoreset", "assay", "contrast", "status", "class", "step",
    "seconds", "warnings", "samples_removed", "k", "pct", "n_controls", "complete_controls", "k_estimate", "e1_pass",
    "n_T", "n_tested", "lost", BLIND, "q_method", "miss_share", "miss_after", "tie_rule_match",
    # r1.16 to r1.21 (35.2): the per-assay E12 to E14 columns, the product's and the record's k fields, the DA outcome
    ASSAY_COLUMNS, "smoke")
newRow <- function(...) {
    row <- stats::setNames(as.list(rep(NA, length(ROW_COLUMNS))), ROW_COLUMNS)
    given <- list(...)
    row[intersect(names(given), ROW_COLUMNS)] <- given[intersect(names(given), ROW_COLUMNS)]
    row
}

#' Rows of one arm: one per (assay, contrast) family, scored on T (and `extra_sets`, named lists of ids to drop).
armRows <- function(sim, journey, T_ids, common, choice, or_errs = NULL, extra_sets = list(), drop_from_T = character()) {
    rows <- list()
    contrasts <- sim$groups[-1L]
    assays <- sim$p$assays
    da <- journey$da
    if (!is.null(da)) da$group <- groupOfComparison(da$comparison)
    sets <- c(list(T = setdiff(T_ids, drop_from_T)), lapply(extra_sets, function(x) setdiff(T_ids, c(drop_from_T, x))))
    errs <- list()
    for (set in names(sets)) for (a in assays) for (g in contrasts) {
        ids <- sets[[set]][sim$assay_of[sets[[set]]] == a]
        ch <- choice[[a]]
        failed <- journey$da_errors
        failed <- if (!is.null(failed)) failed[failed$assay == a & failed$group %in% g, , drop = FALSE]
        family_status <- if (!is.null(failed) && nrow(failed)) "no_da" else journey$status
        family_class <- if (!is.null(failed) && nrow(failed)) failed$error[[1L]] else journey$class
        base <- c(common, list(arm = journey$arm, scoreset = set, assay = a, contrast = g, status = family_status, class = family_class,
            step = journey$step, seconds = journey$seconds, warnings = journey$warnings,
            samples_removed = if (is.null(journey$runs_final)) NA_integer_ else length(sim$runs) - length(journey$runs_final),
            k = ch$k %||% NA_integer_, pct = ch$pct %||% NA_real_, n_controls = length(ch$controls %||% character()),
            complete_controls = ch$complete %||% NA_integer_, k_estimate = ch$estimate %||% NA_integer_,
            e1_pass = if (is.null(journey$e1)) NA else journey$e1$pass,
            miss_after = if (is.null(journey$m_in[[a]])) NA_real_ else mean(!is.finite(journey$m_in[[a]])),
            ctl_planted_share = if (length(ch$controls)) mean(sim$planted[intersect(ch$controls, names(sim$planted))]) else NA_real_,
            ctl_labelled = if (length(ch$controls)) sum(grepl("^(LBL|CTL)_", ch$controls)) else NA_integer_),
            journey$assay_cols[[a]])
        if (is.null(da)) {
            rows[[length(rows) + 1L]] <- do.call(newRow, c(base, list(n_T = length(ids), n_tested = 0L, lost = length(ids))))
            next
        }
        fam <- da[da$assay == a & da$group %in% g, , drop = FALSE]
        or_err <- if (!is.null(or_errs)) or_errs[[paste(set, a, g)]] else NULL
        sc <- scoreFamily(fam, ids, stats::setNames(sim$theta[ids, g], ids), or_err)
        errs[[paste(set, a, g)]] <- sc$err_by_feature
        sc$err_by_feature <- NULL
        rows[[length(rows) + 1L]] <- do.call(newRow, c(base, sc))
    }
    list(rows = rows, errs = errs)
}
choiceOf <- function(journey, legacy) {
    choice <- if (legacy) legacyChoice(journey) else declaredChoice(journey)
    if (is.null(choice)) return(NULL)
    for (a in names(choice)) {
        m <- journey$m_in[[a]]
        ctl <- intersect(choice[[a]]$controls, rownames(m))
        choice[[a]]$complete <- if (is.null(m)) NA_integer_ else sum(rowSums(!is.finite(m[ctl, , drop = FALSE])) == 0L)
    }
    choice
}

# ---- checks ---------------------------------------------------------------------------------------------------------
checks_env <- new.env()
checks_env$rows <- list()
addCheck <- function(check, cell, rep, condition, value, threshold, pass, note = "") {
    checks_env$rows[[length(checks_env$rows) + 1L]] <- data.frame(check = check, cell = cell, rep = rep, condition = condition,
        value = as.character(value), threshold = as.character(threshold), pass = isTRUE(pass), note = oneLine(note),
        stringsAsFactors = FALSE)
    # a unit that is killed still leaves its checks: each is also appended to the unit's progress file (coordinator's ruling, defect 3)
    progress <- getOption("k188.progress")
    if (!is.null(progress)) {
        known <- file.exists(progress)
        try(utils::write.table(checks_env$rows[[length(checks_env$rows)]], progress, sep = "\t", quote = FALSE, row.names = FALSE,
            col.names = !known, append = known), silent = TRUE)
    }
}
#' E9 (Part L, today's base): the legacy controls per assay and percentage are the first round(N_v x pct / 100) rows of
#' the RUV input table among the rows with at least two values, in table order (2.2).
positional <- function(m, pct) {
    valid <- rownames(m)[rowSums(is.finite(m)) > 1L]
    valid[seq_len(round(length(valid) * pct / 100, 0))]
}
checkE9L <- function(sim, journey, condition, arm, d) {
    choice <- legacyChoice(journey)
    if (is.null(choice) || is.null(journey$m_in)) return(addCheck("E9", sim$cell, sim$r, condition, NA, "same sets", FALSE,
        sprintf("%s: no RUV choice (status %s)", arm, journey$status)))
    notes <- character()
    pass <- TRUE
    for (a in names(choice)) {
        m <- journey$m_in[[a]]
        ch <- choice[[a]]
        written <- sim$order_of[[a]][sim$order_of[[a]] %in% rownames(m)]
        notes <- c(notes, sprintf("%s %s: RUV input rows in the written order %s", arm, a, identical(written, rownames(m))))
        if (!is.na(ch$pct)) {
            want <- positional(m, ch$pct)
            same <- setequal(ch$controls, want)
            valid <- rownames(m)[rowSums(is.finite(m)) > 1L]
            pos <- match(ch$controls, valid)
            pass <- pass && same
            notes <- c(notes, sprintf("%s %s pct %s: %d controls, positions among valid rows %s, expected 1..%d: %s", arm, a, ch$pct,
                length(ch$controls), if (length(pos)) paste(range(pos), collapse = "-") else "none", length(want), same))
        }
        trace <- ch$trace
        if (!is.null(trace) && nrow(trace)) {
            n_v <- sum(rowSums(is.finite(m)) > 1L)
            okrow <- !is.na(trace$realized_num_controls)
            counts <- all(trace$realized_num_controls[okrow] == round(n_v * trace$percentage_requested[okrow] / 100, 0))
            pass <- pass && counts
            notes <- c(notes, sprintf("%s %s: counts per percentage = round(N_v pct / 100) for %d of %d percentages: %s", arm, a,
                sum(okrow), nrow(trace), counts))
        }
    }
    # every percentage directly, through the control method the sweep calls, on the RUV input state
    if (identical(arm, "LP") && !is.null(journey$m_in)) {
        object <- omicObjectFromMatrices(sim, journey$m_in)
        for (pct in seq(d$auto_min, d$auto_max)) {
            # LP's own marker (the earlier control rule) is set for the direct call, as the journey sets it
            got <- tryCatch(withMarkers(armMarkers("LP"), quiet(getNegCtrlMetabAnova(object, ruv_grouping_variable = "group",
                percentage_as_neg_ctrl = pct))), error = function(e) e)
            if (inherits(got, "error")) {
                pass <- FALSE
                notes <- c(notes, sprintf("pct %d: error %s", pct, conditionMessage(got)))
                next
            }
            for (a in names(journey$m_in)) {
                same <- setequal(names(got[[a]])[which(got[[a]])], positional(journey$m_in[[a]], pct))
                if (!same) notes <- c(notes, sprintf("pct %d %s: not positional", pct, a))
                pass <- pass && same
            }
        }
    }
    # r1.21 (R20-P02, 40.3 (10)): LP's k is the uncapped sweep's k (the trace's `best_k`); a difference stops before the declaration
    if (identical(arm, "LP")) for (a in names(choice)) {
        same_k <- identical(choice[[a]]$k, choice[[a]]$k_sweep)
        pass <- pass && same_k
        notes <- c(notes, sprintf("%s %s: k %s, uncapped sweep k %s: %s", arm, a, choice[[a]]$k, choice[[a]]$k_sweep, same_k))
    }
    addCheck("E9", sim$cell, sim$r, condition, NA, "same sets", pass, paste(notes, collapse = "; "))
}
#' The RUV input as an object of the omic's class (the state's matrices, in their table order).
omicObjectFromMatrices <- function(sim, mats) {
    omic <- sim$p$omic
    tables <- lapply(mats, function(m) {
        table <- if (identical(omic, "lipid")) data.frame(LipidName = rownames(m), LipidClass = "PC", stringsAsFactors = FALSE) else
            data.frame(Metabolite = rownames(m), Ontology = "PC", stringsAsFactors = FALSE)
        cbind(table, as.data.frame(m, check.names = FALSE))
    })
    if (identical(omic, "lipid")) {
        methods::new("LipidomicsAssayData", lipid_data = tables, design_matrix = sim$design, lipid_id_column = "LipidName",
            annotation_id_column = "LipidClass", sample_id = "Run", group_id = "group", technical_replicate_id = "tech_rep_group",
            database_identifier_type = "Unknown", internal_standard_regex = "^IS_", args = list())
    } else {
        methods::new("MetaboliteAssayData", metabolite_data = tables, design_matrix = sim$design, metabolite_id_column = "Metabolite",
            annotation_id_column = "Ontology", sample_id = "Run", group_id = "group", technical_replicate_id = "tech_rep_group",
            database_identifier_type = "Unknown", internal_standard_regex = "^IS_", args = list())
    }
}
#' E9 (Part P): the LA controls are the top proteins by the helper's rule on `group` with `qvalue` (18.1).
#' The helper's per-protein p (`:555-565`, zero-df proteins dropped by `unlist`), and which were observed in one group.
anovaP <- function(m, group) {
    grps <- group
    ps <- apply(m, 1L, function(x) {
        if (length(unique(grps[!is.na(x)])) > 1L) summary(stats::aov(as.numeric(x) ~ grps))[[1L]][["Pr(>F)"]][1L] else NA_real_
    })
    if (is.list(ps)) ps <- unlist(ps)
    one <- names(ps)[is.na(ps)]
    ps[is.na(ps)] <- 1
    list(p = ps, one_group = one)
}
#' The helper's selection (`:569-582`) on candidates q in table order: the first n of q > 0.05 sorted by q decreasing,
#' ties in table order; and 18.4a (a): the share of the selection decided by table order inside the boundary tie block.
selectByQ <- function(q, n_rows, pct) {
    filtered <- q[q > 0.05]
    n <- min(round(n_rows * pct / 100, 0), length(filtered))
    sorted <- sort(filtered, decreasing = TRUE)
    selected <- names(sorted[seq_len(n)])
    if (n == 0L) return(list(selected = character(), share = NA_real_, block = NA_integer_, block_ids = character()))
    qb <- sorted[[n]]
    block_ids <- names(filtered)[filtered == qb]
    in_selection <- sum(sorted[seq_len(n)] == qb)
    list(selected = selected, share = if (length(block_ids) > in_selection) in_selection / n else 0, block = length(block_ids),
        block_ids = block_ids)
}
jaccard <- function(a, b) if (!length(union(a, b))) NA_real_ else length(intersect(a, b)) / length(union(a, b))
#' 18.4a for one proteomics legacy journey: the shipped rule replicated, (a) the tie share, (b) the Jaccard overlap with the
#' selection on rows in the seeded second order; and whether the replica equals the package's selection.
tieCount <- function(sim, journey, pct, method) {
    m <- journey$m_in$proteins
    if (is.null(m) || !is.finite(pct)) return(NULL)
    group <- sim$design$group[match(colnames(m), sim$design$Run)]
    a <- anovaP(m, group)
    q <- if (identical(method, "qvalue")) qvalue::qvalue(a$p)$qvalues else qvalue::qvalue(a$p, pi0 = 1)$qvalues
    names(q) <- names(a$p)
    shipped <- selectByQ(q, nrow(m), pct)
    set.seed(288000000 + 100000 * sim$p$idx + 10000 * match(sim$scenario, SCENARIOS) + sim$r, kind = "Mersenne-Twister",
        normal.kind = "Inversion", sample.kind = "Rejection")
    order2 <- rownames(m)[sample.int(nrow(m))]
    permuted <- selectByQ(q[order2[order2 %in% names(q)]], nrow(m), pct)
    v <- journey$ruv_value
    got <- if (is.logical(v$controls)) names(v$controls)[which(v$controls)] else character()
    list(tie_share = shipped$share, tie_block = shipped$block, tie_block_onegroup = sum(shipped$block_ids %in% a$one_group),
        perm_jaccard = jaccard(shipped$selected, permuted$selected), tie_rule_match = setequal(got, shipped$selected),
        selected_bh = if (identical(method, "qvalue")) setequal(got, selectByQ(stats::setNames(qvalue::qvalue(a$p, pi0 = 1)$qvalues, names(a$p)),
            nrow(m), pct)$selected) else NA)
}
anovaControls <- function(m, group, pct) {
    grps <- group
    ps <- apply(m, 1L, function(x) {
        if (length(unique(grps[!is.na(x)])) > 1L) summary(stats::aov(as.numeric(x) ~ grps))[[1L]][["Pr(>F)"]][1L] else NA_real_
    })
    # as the helper: a feature with zero residual df has no Pr(>F) (NULL) and `unlist` drops it (KIT-190 audit A11);
    # the count stays round(nrow x pct / 100) over every row, as the helper's
    if (is.list(ps)) ps <- unlist(ps)
    ps[is.na(ps)] <- 1
    pick <- function(q) {
        names(q) <- names(ps)
        filtered <- q[q > 0.05]
        n <- min(round(nrow(m) * pct / 100, 0), length(filtered))
        names(sort(filtered, decreasing = TRUE)[seq_len(n)])
    }
    list(qvalue = pick(qvalue::qvalue(unlist(ps))$qvalues), BH = pick(qvalue::qvalue(unlist(ps), pi0 = 1)$qvalues))
}
checkE9P <- function(sim, journey, condition, tie = NULL) {
    if (!is.null(tie)) return(addCheck("E9", sim$cell, sim$r, condition, NA, "same set as the rule with qvalue", isTRUE(tie$tie_rule_match),
        sprintf("pct %s: equal to the rule with qvalue %s, with BH %s", journey$ruv_value$percentage, tie$tie_rule_match, tie$selected_bh)))
    v <- journey$ruv_value
    if (is.null(v) || is.null(journey$m_in)) return(addCheck("E9", sim$cell, sim$r, condition, NA, "same sets", FALSE,
        sprintf("LA: no RUV choice (status %s)", journey$status)))
    m <- journey$m_in$proteins
    group <- sim$design$group[match(colnames(m), sim$design$Run)]
    got <- if (is.logical(v$controls)) names(v$controls)[which(v$controls)] else character()
    want <- anovaControls(m, group, v$percentage)
    want_q <- want$qvalue
    want_bh <- want$BH
    addCheck("E9", sim$cell, sim$r, condition, length(got), "same set as the rule with qvalue", setequal(got, want_q),
        sprintf("pct %s: %d controls; equal to the rule with qvalue %s, with BH %s", v$percentage, length(got), setequal(got, want_q),
            setequal(got, want_bh)))
}
#' E10 (L1r, on the declaration's base; LA and LM, per assay; 7.4 r1.20). (i) the controls at the sweep's chosen percentage are
#' those of rule (d) computed here on the RUV input: candidates ranked by the group ANOVA's p, highest first, ties by feature
#' identifier in byte order, features observed in one group and (as the helper's `unlist` drops them) features with no residual
#' degrees of freedom left out, the first round(N_v x pct / 100); BH's q (the control method's default `ruv_fdr_method`),
#' the qvalue form recorded beside. (ii) LA: the sweep's k (the optimisation trace's) equals the product's `k_sweep` and the k
#' used is min(k_sweep, `calculateAdaptiveMaxK`(the assay's samples)); LM's typed k is unchanged; K172c's LA has k used at most
#' 2. (iii) the analysis record's fields equal the product's (`ruv_record_match`) and the cap equals `calculateAdaptiveMaxK`.
ruleDControls <- function(m, group, pct, fdr = "BH") {
    m[!is.finite(m)] <- NA
    f <- m[rowSums(!is.na(m)) > 1L, colSums(!is.na(m)) > 1L, drop = FALSE]
    grps <- group[match(colnames(f), names(group))]
    ps <- apply(f, 1L, function(x) {
        if (length(unique(grps[!is.na(x)])) > 1L) summary(stats::aov(as.numeric(x) ~ grps))[[1L]][["Pr(>F)"]][1L] else NA_real_
    })
    if (is.list(ps)) ps <- unlist(ps)
    ps[is.na(ps)] <- 1
    q <- if (identical(fdr, "BH")) qvalue::qvalue(ps, pi0 = 1)$qvalues else qvalue::qvalue(ps)$qvalues
    names(q) <- names(ps)
    seen <- apply(f, 1L, function(x) length(unique(grps[!is.na(x)])))
    candidates <- intersect(names(q)[q > 0.05], names(seen)[seen > 1L])
    ranked <- candidates[order(-ps[candidates], candidates, method = "radix")]
    ranked[seq_len(min(round(nrow(f) * pct / 100, 0), length(ranked)))]
}
checkE10 <- function(sim, journey, arm, condition, d) {
    choice <- legacyChoice(journey)
    if (is.null(choice) || is.null(journey$m_in)) return(addCheck("E10", sim$cell, sim$r, paste(condition, arm), NA, "rule (d) and the bound",
        FALSE, sprintf("%s: no RUV choice (status %s)", arm, journey$status)))
    group <- stats::setNames(sim$design$group, sim$design$Run)
    notes <- character()
    pass <- TRUE
    for (a in names(choice)) {
        ch <- choice[[a]]
        m <- journey$m_in[[a]]
        cols <- journey$assay_cols[[a]]
        cap <- calculateAdaptiveMaxK(ncol(m))
        if (!is.na(ch$pct) && length(ch$controls)) {
            same <- setequal(ch$controls, ruleDControls(m, group, ch$pct, "BH"))
            same_q <- setequal(ch$controls, ruleDControls(m, group, ch$pct, "qvalue"))
            pass <- pass && same
            notes <- c(notes, sprintf("%s %s pct %s: controls equal rule (d) with BH %s (with qvalue %s)", arm, a, ch$pct, same, same_q))
        }
        sweep_ok <- TRUE # the sweep's k is compared with an untouched sweep in E10u (SA-P07, the coordinator's ruling)
        used_ok <- is.na(ch$k) || (if (identical(arm, "LM")) isTRUE(ch$k == as.integer(d$manual_k)) else is.na(ch$k_sweep) ||
            isTRUE(ch$k == min(ch$k_sweep, cap)))
        if (startsWith(sim$cell, "K172c") && identical(arm, "LA") && !is.na(ch$k)) used_ok <- used_ok && ch$k <= 2L
        record_ok <- isTRUE(cols$ruv_record_match) && isTRUE(as.integer(cols$ruv_k_cap) == cap)
        pass <- pass && sweep_ok && used_ok && record_ok
        notes <- c(notes, sprintf("%s %s: k_sweep %s (product %s), k used %s, cap %d: bound %s, record %s", arm, a, ch$k_sweep, cols$ruv_k_sweep,
            ch$k, cap, used_ok, record_ok))
    }
    addCheck("E10", sim$cell, sim$r, paste(condition, arm), NA, "rule (d), the bound, the record", pass, paste(notes, collapse = "; "))
}
#' E10u (SA-P07): LA's sweep (the optimisation trace: percentage and k) equals the untouched sweep's, run on the same data set as a
#' check-only journey with `MultiScholaR.ruv_k_bound = "none"`; that journey's k used equals its sweep k; LA's product `k_sweep`
#' equals the untouched sweep's k.
checkE10u <- function(sim, la, lau, condition) {
    cl <- legacyChoice(la)
    cu <- legacyChoice(lau)
    notes <- character()
    pass <- !is.null(cl) && !is.null(cu)
    for (a in names(cl)) {
        same <- identical(cl[[a]]$k_sweep, cu[[a]]$k_sweep) && identical(cl[[a]]$pct, cu[[a]]$pct) &&
            identical(cu[[a]]$k, cu[[a]]$k_sweep) && identical(as.integer(la$assay_cols[[a]]$ruv_k_sweep), cu[[a]]$k_sweep)
        pass <- pass && isTRUE(same)
        notes <- c(notes, sprintf("%s: LA sweep k %s pct %s, untouched sweep k %s pct %s, product k_sweep %s: %s", a, cl[[a]]$k_sweep, cl[[a]]$pct,
            cu[[a]]$k_sweep, cu[[a]]$pct, la$assay_cols[[a]]$ruv_k_sweep, same))
    }
    addCheck("E10u", sim$cell, sim$r, condition, NA, "the sweep is untouched by the bound", pass, paste(notes, collapse = "; "))
}
#' E11 (L1r; 34): the LAs and LA journeys of a data set share the RUV input, the controls and k; their corrected matrices
#' differ by one constant per feature on observed cells (range at most 1e-9); the fold changes agree to 1e-10 where the same
#' samples are kept (KIT-224 S1 (i), (ii)).
checkE11 <- function(sim, la, las, condition) {
    notes <- character()
    pass <- TRUE
    cl <- legacyChoice(la)
    cs <- legacyChoice(las)
    for (a in sim$p$assays) {
        same_in <- isTRUE(all.equal(la$m_in[[a]], las$m_in[[a]], tolerance = 0, check.attributes = TRUE))
        same_choice <- setequal(cl[[a]]$controls, cs[[a]]$controls) && identical(cl[[a]]$k, cs[[a]]$k)
        x <- la$m_out[[a]]
        y <- las$m_out[[a]]
        spread <- NA_real_
        if (!is.null(x) && !is.null(y)) {
            rows <- intersect(rownames(x), rownames(y))
            cols <- intersect(colnames(x), colnames(y))
            delta <- x[rows, cols, drop = FALSE] - y[rows, cols, drop = FALSE]
            spread <- suppressWarnings(max(apply(delta, 1L, function(v) if (any(is.finite(v))) diff(range(v[is.finite(v)])) else 0)))
        }
        pass <- pass && same_in && same_choice && (is.na(spread) || spread <= 1e-9)
        notes <- c(notes, sprintf("%s: same RUV input %s, same controls and k %s, per-feature range of the difference %s", a, same_in, same_choice,
            signif(spread, 4)))
    }
    dev <- NA_real_
    if (!is.null(la$da) && !is.null(las$da) && identical(la$runs_final, las$runs_final)) {
        at <- match(paste(las$da$assay, las$da$comparison, las$da$feature), paste(la$da$assay, la$da$comparison, la$da$feature))
        dev <- suppressWarnings(max(abs(las[["da"]][["logFC"]] - la[["da"]][["logFC"]][at]), na.rm = TRUE))
        pass <- pass && is.finite(dev) && dev <= 1e-10
    }
    addCheck("E11", sim$cell, sim$r, condition, signif(dev, 4), "constant per feature at most 1e-9; logFC 1e-10", pass, paste(notes, collapse = "; "))
}
#' E8 and M-amp: e_j derived from RUV-III on the complete RUV input with LA's k and controls, against the route.
checkE8 <- function(sim, journey, choice, condition) {
    out <- list()
    for (a in names(choice)) {
        ch <- choice[[a]]
        Yin <- journey$m_in[[a]]
        Yout <- journey$m_out[[a]]
        if (is.na(ch$k) || is.null(Yin) || is.null(Yout) || any(!is.finite(Yin))) next
        Y <- t(Yin)
        group <- sim$design$group[match(rownames(Y), sim$design$Run)]
        M <- stats::model.matrix(~ 0 + factor(group))
        R <- Y - M %*% solve(crossprod(M), crossprod(M, Y))
        U <- svd(R, nu = ch$k, nv = 0L)$u[, seq_len(ch$k), drop = FALSE]
        alpha <- crossprod(U, Y)
        ctl <- colnames(Y) %in% ch$controls
        ac <- alpha[, ctl, drop = FALSE]
        W <- Y[, ctl, drop = FALSE] %*% t(ac) %*% solve(ac %*% t(ac))
        change_ref <- -(W %*% alpha)
        common <- intersect(colnames(Y), rownames(Yout))
        runs <- intersect(rownames(Y), colnames(Yout))
        change_route <- t(Yout[common, runs, drop = FALSE]) - Y[runs, common, drop = FALSE]
        worst <- 0
        e_max <- 0
        for (g in sim$groups[-1L]) {
            lambda <- (group[match(runs, rownames(Y))] == g) / sum(group[match(runs, rownames(Y))] == g) -
                (group[match(runs, rownames(Y))] == "A") / sum(group[match(runs, rownames(Y))] == "A")
            d_ref <- drop(crossprod(lambda, change_ref[runs, common, drop = FALSE]))
            d_route <- drop(crossprod(lambda, change_route))
            worst <- max(worst, abs(d_ref - d_route), na.rm = TRUE)
            e_max <- max(e_max, abs(d_ref), na.rm = TRUE)
        }
        tolerance <- if (length(sim$runs) > 20L) 1e-6 else 1e-8
        addCheck("E8", sim$cell, sim$r, paste(condition, a), signif(worst, 4), tolerance, worst <= tolerance,
            sprintf("k %d, controls %d (reported, never re-run)", ch$k, sum(ctl)))
        # P13-1 (REPORTED ONLY, decides nothing; the coordinator's ruling (a) on KIT-275 P13): the slope of LA's corrected contrast on its uncorrected
        # contrast across features, per contrast from the matrices above (mask none: no missing values), averaged over the contrasts; about 1 means the
        # fold changes are unchanged by the correction
        slopes <- vapply(sim$groups[-1L], function(g) {
            r_ix <- group[match(runs, rownames(Y))]
            lambda <- (r_ix == g) / sum(r_ix == g) - (r_ix == "A") / sum(r_ix == "A")
            u <- drop(crossprod(lambda, Y[runs, common, drop = FALSE]))
            cc <- drop(crossprod(lambda, t(Yout[common, runs, drop = FALSE])))
            if (stats::var(u) > 0) unname(stats::coef(stats::lm(cc ~ u))[[2L]]) else NA_real_
        }, numeric(1L))
        out[[a]] <- list(mamp_slope = mean(slopes, na.rm = TRUE), mamp_sigma_min = min(svd(ac)$d), mamp_e_max = e_max,
            mamp_planted_share = mean(sim$planted[intersect(ch$controls, names(sim$planted))]))
    }
    out
}
#' M-amp, closed-form branch (KIT-193 audit U3; reported): the features whose observed runs leave exactly k residual
#' dimensions, those with fewer complete controls than observed runs (Y_c Y_c' singular), and the largest log10 condition
#' number of Y_c Y_c' (runs by runs, on the observed RUV input values), per assay, with the arm's k and controls.
cfBranch <- function(sim, journey, choice) {
    out <- list()
    for (a in names(choice)) {
        ch <- choice[[a]]
        Y <- journey$m_in[[a]]
        if (is.null(Y) || is.na(ch$k %||% NA)) next
        obs <- is.finite(Y)
        group <- sim$design$group[match(colnames(Y), sim$design$Run)]
        sets <- vapply(unique(group), function(g) rowSums(obs[, group == g, drop = FALSE]) > 0L, logical(nrow(Y)))
        if (is.null(dim(sets))) sets <- matrix(sets, nrow = nrow(Y))
        rd <- rowSums(obs) - rowSums(sets)
        cf <- which(rd == ch$k)
        ctl <- intersect(ch$controls, rownames(Y))
        singular <- 0L
        worst <- NA_real_
        for (j in cf) {
            runs <- which(obs[j, ])
            cj <- ctl[rowSums(!obs[ctl, runs, drop = FALSE]) == 0L]
            if (length(cj) < length(runs)) singular <- singular + 1L
            Yc <- t(Y[cj, runs, drop = FALSE])
            kap <- if (length(cj)) tryCatch(kappa(Yc %*% t(Yc), exact = TRUE), error = function(e) Inf) else Inf
            worst <- max(worst, log10(kap), na.rm = TRUE)
        }
        out[[a]] <- list(mamp_cf_n = length(cf), mamp_cf_singular = singular, mamp_cf_log10kappa = worst)
    }
    out
}
#' M-df: the corrected matrix refitted as the route does, and with the residual variance re-estimated on df - k.
mdf <- function(sim, journey, choice) {
    out <- list()
    for (a in names(choice)) {
        ch <- choice[[a]]
        Y <- journey$m_out[[a]]
        if (is.null(Y) || is.na(ch$k)) next
        group <- factor(sim$design$group[match(colnames(Y), sim$design$Run)], levels = sim$groups)
        design <- stats::model.matrix(~ 0 + group)
        colnames(design) <- levels(group)
        cm <- limma::makeContrasts(contrasts = sprintf("%s-A", sim$groups[-1L]), levels = colnames(design))
        nulls <- function(fit, rows) {
            eb <- limma::eBayes(limma::contrasts.fit(fit, cm), trend = TRUE, robust = TRUE)
            p <- eb$p.value
            p <- p[rownames(p) %in% rownames(sim$theta), , drop = FALSE]
            th <- sim$theta[rownames(p), , drop = FALSE]
            mean(p[th == 0] < 0.05, na.rm = TRUE)
        }
        fit <- quiet(limma::lmFit(Y, design))
        route <- quiet(nulls(fit))
        ok <- is.finite(fit$df.residual) & fit$df.residual > ch$k & is.finite(fit$sigma)
        reduced <- fit[ok, ]
        reduced$sigma <- reduced$sigma * sqrt(reduced$df.residual / (reduced$df.residual - ch$k))
        reduced$df.residual <- reduced$df.residual - ch$k
        out[[a]] <- list(mdf_null05 = route, mdf_null05_reduced = quiet(nulls(reduced)), mdf_df_le_k = sum(!ok), mdf_tested = nrow(Y))
    }
    out
}
#' E4: the estimand functions on a constructed family.
checkE4 <- function() {
    ids <- sprintf("F%02d", 1:10)
    theta <- stats::setNames(c(1, 1, 0, 0, 0, 0, 0, 0, 0, 0), ids)
    err <- c(0.5, 0, 6, 0, 5.5, 0, 0, 0, 1.5, 0)
    da <- data.frame(feature = ids, logFC = theta + err, p = c(0.001, 0.2, 0.01, 0.04, 0.3, 0.5, 0.6, 0.7, 0.8, 0.9),
        q = c(0.01, 0.3, 0.04, 0.2, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95), stringsAsFactors = FALSE)
    or_err <- stats::setNames(c(0, 0, 6, 0, 0, 0, 0, 0, 0, 0), ids)
    a <- scoreFamily(da, ids, theta, or_err)
    want <- list(n_null = 8L, null_p05 = 0.25, null_p01 = 0, calls = 2L, false_calls = 1L, any_call = 1L, fdp = 0.5, power = 0.5,
        gross1 = 0.3, s5 = 1L, n_common = 10L, gross1_common = 0.3, or_gross1_common = 0.1, err_max = 6, slope = 1.25)
    da2 <- da
    da2$logFC[5] <- 0
    b <- scoreFamily(da2, ids, theta, or_err)
    o <- scoreFamily(da, ids, theta, NULL)
    bad <- names(want)[!vapply(names(want), function(n) isTRUE(all.equal(as.numeric(a[[n]]), as.numeric(want[[n]]), tolerance = 1e-12)),
        logical(1L))]
    pass <- !length(bad) && identical(b$s5, 0L) && identical(o$or_s5, 1L)
    none <- scoreFamily(da[0, ], ids, theta, NULL)
    pass <- pass && identical(none$lost, 10L) && identical(none$fdp, 0) && identical(none$any_call, 0L)
    addCheck("E4", "constructed", 0L, "family", NA, "exact", pass, sprintf("differ: %s; s5 without the feature %s; oracle s5 %s",
        paste(bad, collapse = ","), b$s5, o$or_s5))
    # E4b (18.4a): five candidates, three tied at the boundary, two selected: share 1, block 3; a second order c, b, ...
    # selects c, b: Jaccard 1/3; no tie at the boundary: share 0
    q <- c(a = 0.9, b = 0.9, c = 0.9, d = 0.5, e = 0.3)
    one <- selectByQ(q, 40, 5)
    two <- selectByQ(q[c("c", "b", "a", "d", "e")], 40, 5)
    three <- selectByQ(c(a = 0.9, b = 0.8, c = 0.7), 40, 5)
    pass_b <- identical(one$selected, c("a", "b")) && identical(one$share, 1) && identical(one$block, 3L) &&
        identical(jaccard(one$selected, two$selected), 1 / 3) && identical(three$share, 0) && identical(three$selected, c("a", "b"))
    addCheck("E4b", "constructed", 0L, "ties", NA, "exact", pass_b, sprintf("share %s block %s jaccard %s; no-tie share %s", one$share,
        one$block, jaccard(one$selected, two$selected), three$share))
}

# ---- one unit ------------------------------------------------------------------------------------------------------
ARMS_OF <- list(L = list(deciding = c("LA", "GA", "NO", "OR"), reported = c("LA", "OR", "LM", "FP", "SC", "LCm"),
    munit = c("LA", "LM", "GA", "NO"), ba = c("LA", "GA", "NO", "OR"), lp = c("LP", "OR"), las = c("LAs", "OR"),
    k172c = c("LA", "LM", "GA", "NO", "OR")),
    P = list(deciding = c("LA", "GA", "NO", "OR"), reported = c("LM", "OR"), ba = c("LA", "GA", "NO", "OR")))
conditionsOf <- function(part, cellname, job) {
    scen <- scenariosOf(cellname)
    switch(job,
        deciding = expand.grid(mask = MASKS, scenario = scen, stringsAsFactors = FALSE),
        reported = , lp = , las = expand.grid(mask = MASKS, scenario = SCENARIOS[1:4], stringsAsFactors = FALSE),
        munit = expand.grid(mask = c("none", "vd20"), scenario = c("N0", "E1"), stringsAsFactors = FALSE),
        ba = expand.grid(mask = "ba", scenario = SCENARIOS[1:4], stringsAsFactors = FALSE),
        k172c = data.frame(mask = "none", scenario = "N0", stringsAsFactors = FALSE),
        l1c = data.frame(mask = "none", scenario = "N1", stringsAsFactors = FALSE),
        guard = data.frame(mask = "none", scenario = "N0", stringsAsFactors = FALSE))
}
jobApplies <- function(cellname, job) {
    if (identical(cellname, "K172c") || identical(job, "k172c")) return(identical(cellname, "K172c") && identical(job, "k172c"))
    # 35.6: L1c is one job on L12; 38.5: L24 is a reported cell with its deciding-type job only (LA, GA, NO, OR; R = 100); LP
    # (5.3) and LAs (34) run in the four shapes of the base table
    if (identical(job, "guard")) return(identical(cellname, "M12"))
    if (identical(job, "l1c")) return(identical(cellname, "L12"))
    if (identical(cellname, "L24")) return(identical(job, "deciding"))
    switch(job, munit = cellname %in% c("L12", "M12"), ba = cellname %in% c("L64", "M168", "P4b"),
        lp = , las = cellname %in% c("L12", "L64", "M12", "M168"), TRUE)
}

runConditionL <- function(sim, mask, job, d, lifecycle, fault_now, first_condition) {
    rows <- list()
    common <- list(part = "L", cell = sim$cell, job = job, scenario = sim$scenario, mask = mask, rep = sim$r)
    values <- maskedValues(sim, mask)
    miss <- mean(is.na(values[sim$feat_ids, ]))
    common$miss_share <- miss
    requests <- REQUESTS_L(d)
    arms <- ARMS_OF$L[[job]]
    # 34 (E11): in the lifecycle the LAs job also runs LA on the same data set, for the comparison only (its rows are not written)
    if (lifecycle && identical(job, "las")) arms <- c("LA", arms)
    # E10u: one check-only LA journey per L1r cell (N1, mask none) with the bound off; its rows are not written
    e10u <- lifecycle && identical(job, "deciding") && identical(sim$scenario, "N1") && identical(mask, "none")
    if (e10u) arms <- append(arms, "LAu", after = match("LA", arms))
    tmp <- function(arm) tempfile(sprintf("k188-%s-%s-", sim$cell, arm))
    journeys <- list()
    choices <- list()
    T_ids <- NULL
    or_errs <- NULL
    condition <- paste(sim$scenario, mask)
    for (arm in arms) {
        fault <- isTRUE(fault_now) && identical(arm, arms[[1L]])
        if (identical(arm, "OR")) {
            if (is.null(T_ids)) next
            oracle <- maskedValues(sim, mask, oracle = TRUE)
            j <- journeyL(sim, oracle, "OR", list(gate = FALSE, ruv = NULL), d, tmp(arm), filter = FALSE, itsd = FALSE, method = "none",
                correlation = FALSE, ids = T_ids)
            res <- armRows(sim, j, T_ids, common, NULL)
            or_errs <- res$errs
            rows <- c(rows, res$rows)
            journeys$OR <- j
            next
        }
        request <- if (arm %in% names(requests)) requests[[arm]] else if (arm %in% c("LP", "LAs", "LAu")) requests$LA else NULL
        drop <- character()
        extra <- list()
        if (identical(arm, "SC")) {
            request <- list(gate = TRUE, ruv = list(ruv_design = "supplied_controls", ruv_grouping_variable = "group", ruv_k = NULL,
                ruv_k_rule = "num_sv_be", ruv_k_seed = RUV_FACTOR_COUNT_SEED, ruv_control_ids = sim$lbl_ids))
            drop <- sim$lbl_ids
        }
        if (identical(arm, "LCm")) {
            la <- choices$LA
            if (is.null(la) || any(is.na(vapply(la, `[[`, integer(1L), "k"))) || any(lengths(lapply(la, `[[`, "controls")) == 0L)) {
                rows[[length(rows) + 1L]] <- do.call(newRow, c(common, list(arm = "LCm", status = "not_run", class = "LA gave no k or controls")))
                next
            }
            request <- list(gate = TRUE, ruv = list(ruv_design = "supplied_controls", ruv_grouping_variable = "group",
                ruv_k = lapply(la, `[[`, "k"), ruv_k_rule = "user", ruv_k_seed = RUV_FACTOR_COUNT_SEED,
                ruv_control_ids = unlist(lapply(la, `[[`, "controls"), use.names = FALSE)))
            extra <- list(T_noLAc = unlist(lapply(la, `[[`, "controls"), use.names = FALSE))
        }
        j <- journeyL(sim, values, arm, request, d, tmp(arm), fault = fault)
        legacy <- arm %in% c("LA", "LM", "LP", "LAs", "LAu")
        ch <- choiceOf(j, legacy)
        choices[[arm]] <- ch
        if (is.null(T_ids) && !is.null(j$kept)) T_ids <- setdiff(j$kept, sim$is_ids)
        if (!is.null(T_ids) && !is.null(j$kept) && !setequal(setdiff(j$kept, sim$is_ids), T_ids)) {
            addCheck("T_same", sim$cell, sim$r, condition, NA, "the filter keeps the same set in every arm", FALSE, arm)
        }
        journeys[[arm]] <- j
        if (identical(arm, "LA") && identical(job, "reported")) extra <- list(T_noLAc = unlist(lapply(ch, `[[`, "controls"), use.names = FALSE))
        res <- armRows(sim, j, T_ids %||% character(), common, ch, NULL, extra, drop)
        if (!(lifecycle && identical(job, "las") && identical(arm, "LA")) && !identical(arm, "LAu")) rows <- c(rows, res$rows)
        # E9 describes the positional selection, which only LP runs on the bound commit (40.3 (10)); E10 checks the rule that ships
        if (lifecycle && identical(arm, "LP") && j$status %in% c("ok", "dropped", "uncorrected")) checkE9L(sim, j, condition, arm, d)
        if (lifecycle && arm %in% c("LA", "LM") && j$status %in% c("ok", "dropped", "uncorrected")) checkE10(sim, j, arm, condition, d)
        # E14 (7.4): every row of a journey that reached the RUV step has a per-assay status, whatever its step; the script's
        # status equals the product's, but for a strictness of identical() (a product finding for KIT-242, not a stop)
        if (lifecycle && !is.null(request$ruv) && isTRUE(j$ruv_reached)) {
            statuses <- vapply(j$assay_cols, function(x) as.character(x$ruv_status), "")
            agree <- vapply(j$assay_cols, function(x) x$ruv_status_agree, NA)
            strict <- vapply(j$assay_cols, function(x) identical(x$ruv_status, "uncorrected") && identical(x$ruv_status_product, "corrected"), NA)
            ok <- !anyNA(statuses) && all(nzchar(statuses)) && all(is.na(agree) | agree | strict)
            if (!ok || first_condition) addCheck("E14", sim$cell, sim$r, paste(condition, arm), paste(statuses, collapse = "/"),
                "a status per assay; script equals product", ok, sprintf("product %s; identical() strictness only: %s",
                paste(vapply(j$assay_cols, function(x) as.character(x$ruv_status_product), ""), collapse = "/"), any(strict, na.rm = TRUE)))
        }
        if (lifecycle && !is.null(j$e1)) {
            if (!isTRUE(j$e1$pass) || first_condition) addCheck("E1", sim$cell, sim$r, paste(condition, arm), NA, "every step as declared",
                j$e1$pass, j$e1$note)
        } else if (lifecycle && is.null(j$e1) && identical(j$status, "crashed") && !grepl("deliberate fault", j$class, fixed = TRUE)) {
            addCheck("G1_crash", sim$cell, sim$r, paste(condition, arm), NA, "recorded, counted (G1)", FALSE, sprintf("journey %s at %s: %s",
                j$status, j$step, j$class))
        }
    }
    # the oracle's |error| per family, then each arm's attributable and paired figures (S5, F2)
    if (!is.null(or_errs)) {
        for (i in seq_along(rows)) {
            r <- rows[[i]]
            if (identical(r$arm, "OR") || is.na(r$assay) || is.null(journeys[[r$arm]]$da)) next
            fam <- journeys[[r$arm]]$da
            fam$group <- groupOfComparison(fam$comparison)
            fam <- fam[fam$assay == r$assay & fam$group %in% r$contrast, , drop = FALSE]
            ids <- T_ids[sim$assay_of[T_ids] == r$assay]
            if (identical(r$arm, "SC")) ids <- setdiff(ids, sim$lbl_ids)
            if (identical(r$scoreset, "T_noLAc")) ids <- setdiff(ids, unlist(lapply(choices$LA, `[[`, "controls")))
            sc <- scoreFamily(fam, ids, stats::setNames(sim$theta[ids, r$contrast], ids), or_errs[[paste("T", r$assay, r$contrast)]])
            rows[[i]][c("n_common", "s5", "gross1_common", "or_gross1_common")] <- sc[c("n_common", "s5", "gross1_common", "or_gross1_common")]
        }
    }
    # E8 and M-amp (LA, mask none)
    if (identical(mask, "none") && !is.null(journeys$LA) && job %in% c("deciding") && !is.null(choices$LA)) {
        amp <- checkE8(sim, journeys$LA, choices$LA, condition)
        for (i in seq_along(rows)) if (identical(rows[[i]]$arm, "LA") && !is.null(byName(amp, rows[[i]]$assay))) {
            rows[[i]][names(byName(amp, rows[[i]]$assay))] <- byName(amp, rows[[i]]$assay)
        }
    }
    # M-amp, closed-form branch (LA, deciding, every mask; reported)
    if (identical(job, "deciding") && !is.null(journeys$LA) && !is.null(choices$LA)) {
        cf <- cfBranch(sim, journeys$LA, choices$LA)
        for (i in seq_along(rows)) if (identical(rows[[i]]$arm, "LA") && !is.null(byName(cf, rows[[i]]$assay))) rows[[i]][names(byName(cf, rows[[i]]$assay))] <- byName(cf, rows[[i]]$assay)
    }
    # M-df (LA and LM, N0 and N1, mask none; L12 and M12)
    if (identical(job, "reported") && identical(mask, "none") && sim$scenario %in% c("N0", "N1") && sim$cell %in% c("L12", "M12")) {
        for (arm in intersect(c("LA", "LM"), names(choices))) {
            res <- mdf(sim, journeys[[arm]], choices[[arm]])
            for (i in seq_along(rows)) if (identical(rows[[i]]$arm, arm) && !is.null(byName(res, rows[[i]]$assay))) {
                rows[[i]][names(byName(res, rows[[i]]$assay))] <- byName(res, rows[[i]]$assay)
            }
        }
    }
    # R16-P13 (35.6): the metabolomics ITSD identity check, M12 data set 401, N0 and E1, masks none and VD 20 %, NO and LA: ITSD ticked
    # (regex, no standard in the data) against ITSD off gives identical DA tables
    if (lifecycle && identical(job, "deciding") && identical(sim$cell, "M12") && identical(sim$r, 401L) && sim$scenario %in% c("N0", "E1") &&
        mask %in% c("none", "vd20")) {
        for (arm in c("NO", "LA")) if (!is.null(journeys[[arm]]$da)) {
            j2 <- journeyL(sim, values, arm, requests[[arm]], d, tmp(arm), itsd = TRUE)
            if (is.null(j2$da)) {
                addCheck("E3_itsd_identity", sim$cell, sim$r, paste(condition, arm), NA, "identical", FALSE, sprintf("ITSD ticked: %s at %s: %s", j2$status, j2$step, j2$class))
                next
            }
            key <- function(x) x[order(x$assay, x$comparison, x$feature), c("assay", "comparison", "feature", "logFC", "t", "p", "q")]
            a <- key(journeys[[arm]]$da)
            b <- key(j2$da)
            rownames(a) <- rownames(b) <- NULL
            addCheck("E3_itsd_identity", sim$cell, sim$r, paste(condition, arm), NA, "identical DA tables, ITSD ticked against off", identical(a, b),
                sprintf("rows %d and %d; itsd record fields %s", nrow(a), nrow(b), j2$itsd_record_fields %||% ""))
        }
    }
    if (e10u && !is.null(journeys$LA) && !is.null(journeys$LAu)) checkE10u(sim, journeys$LA, journeys$LAu, condition)
    if (lifecycle && identical(job, "las") && !is.null(journeys$LA) && !is.null(journeys$LAs)) checkE11(sim, journeys$LA, journeys$LAs, condition)
    # E2: RUV off as Skip (gate off) and as the design none (gate on) give the same DA (lifecycle; E1 vd20)
    if (lifecycle && identical(job, "deciding") && identical(sim$scenario, "E1") && identical(mask, "vd20") && !is.null(journeys$NO$da)) {
        j2 <- journeyL(sim, values, "NOnone", requests$NOnone, d, tmp("NOnone"))
        if (is.null(j2$da)) {
            addCheck("E2", sim$cell, sim$r, condition, NA, "identical", FALSE, sprintf("design none: %s at %s: %s", j2$status, j2$step, j2$class))
        } else {
            key <- function(x) x[order(x$assay, x$comparison, x$feature), c("assay", "comparison", "feature", "logFC", "t", "p", "q")]
            a <- key(journeys$NO$da)
            b <- key(j2$da)
            rownames(a) <- rownames(b) <- NULL
            same <- identical(a, b)
            diff <- if (nrow(a) == nrow(b)) max(abs(as.matrix(a[, 4:7]) - as.matrix(b[, 4:7])), na.rm = TRUE) else NA
            addCheck("E2", sim$cell, sim$r, condition, signif(diff, 4), "identical", same, sprintf("rows %d and %d", nrow(a), nrow(b)))
        }
    }
    list(rows = rows, journeys = journeys, choices = choices, T_ids = T_ids)
}

#' M-unit: the s = 1 journeys and the journeys at 1,000 and 0.001 (after the filter, on the features kept at s = 1).
runConditionMunit <- function(sim, mask, d) {
    rows <- list()
    common <- list(part = "L", cell = sim$cell, job = "munit", scenario = sim$scenario, mask = mask, rep = sim$r)
    values <- maskedValues(sim, mask)
    requests <- REQUESTS_L(d)
    itsd <- identical(sim$p$omic, "lipid")
    for (arm in ARMS_OF$L$munit) {
        legacy <- arm %in% c("LA", "LM")
        ref <- journeyL(sim, values, arm, requests[[arm]], d, tempfile("k188-mu-"), keep_states = TRUE)
        ref_choice <- choiceOf(ref, legacy)
        T_ids <- setdiff(ref$kept %||% character(), sim$is_ids)
        rows <- c(rows, armRows(sim, ref, T_ids, c(common, list(arm_scale = "1")), ref_choice)$rows)
        for (s in c(1000, 0.001)) {
            if (is.null(ref$kept)) next
            j <- journeyL(sim, values, arm, requests[[arm]], d, tempfile("k188-mu-"), filter = FALSE, ids = ref$kept, scale = s,
                keep_states = TRUE)
            ch <- choiceOf(j, legacy)
            res <- armRows(sim, j, T_ids, common, ch)$rows
            dev <- function(a, b, shift) {
                if (is.null(a) || is.null(b)) return(NA_real_)
                max(vapply(intersect(names(a), names(b)), function(x) {
                    r <- intersect(rownames(a[[x]]), rownames(b[[x]]))
                    cc <- intersect(colnames(a[[x]]), colnames(b[[x]]))
                    max(abs(a[[x]][r, cc] - b[[x]][r, cc] - shift), na.rm = TRUE)
                }, numeric(1L)))
            }
            centre <- function(m) lapply(m, function(x) x - rowMeans(x, na.rm = TRUE))
            log_ratio <- function(m) lapply(m, log2)
            same_ctl <- if (legacy && !is.null(ch) && !is.null(ref_choice)) all(vapply(names(ch), function(a) setequal(ch[[a]]$controls,
                ref_choice[[a]]$controls), logical(1L))) else NA
            same_k <- if (!is.null(ch) && !is.null(ref_choice)) identical(lapply(ch, `[[`, "k"), lapply(ref_choice, `[[`, "k")) else NA
            same_pct <- if (legacy && !is.null(ch) && !is.null(ref_choice)) identical(lapply(ch, `[[`, "pct"), lapply(ref_choice, `[[`, "pct")) else NA
            da_dev <- function(column) {
                if (is.null(j$da) || is.null(ref$da)) return(NA_real_)
                k1 <- paste(ref$da$assay, ref$da$comparison, ref$da$feature)
                k2 <- paste(j$da$assay, j$da$comparison, j$da$feature)
                at <- match(k2, k1)
                max(abs(j$da[[column]] - ref$da[[column]][at]), na.rm = TRUE)
            }
            calls_changed <- if (is.null(j$da) || is.null(ref$da)) NA_integer_ else {
                k1 <- paste(ref$da$assay, ref$da$comparison, ref$da$feature)
                at <- match(paste(j$da$assay, j$da$comparison, j$da$feature), k1)
                c1 <- ifelse(is.finite(ref$da$q[at]) & ref$da$q[at] < 0.05, sign(ref$da$logFC[at]), 0)
                c2 <- ifelse(is.finite(j$da$q) & j$da$q < 0.05, sign(j$da$logFC), 0)
                sum(c1 != c2, na.rm = TRUE)
            }
            notrend <- NA_real_
            if (identical(arm, "LM") && identical(mask, "none") && !is.null(j$m_out) && !is.null(ref$m_out)) {
                tstat <- function(mats) unlist(lapply(mats, function(Y) {
                    group <- factor(sim$design$group[match(colnames(Y), sim$design$Run)], levels = sim$groups)
                    design <- stats::model.matrix(~ 0 + group)
                    colnames(design) <- levels(group)
                    cm <- limma::makeContrasts(contrasts = sprintf("%s-A", sim$groups[-1L]), levels = colnames(design))
                    eb <- limma::eBayes(limma::contrasts.fit(limma::lmFit(Y, design), cm), trend = FALSE, robust = TRUE)
                    stats::setNames(as.vector(eb$t), paste(rep(rownames(eb$t), ncol(eb$t)), rep(colnames(eb$t), each = nrow(eb$t))))
                }))
                t1 <- quiet(tstat(ref$m_out))
                ts <- quiet(tstat(j$m_out))
                notrend <- max(abs(ts - t1[names(ts)]), na.rm = TRUE)
            }
            mu <- list(arm_scale = as.character(s), mu_itsd_dev = if (itsd) dev(log_ratio(j$m_itsd), log_ratio(ref$m_itsd), log2(s)) else NA_real_,
                mu_log2_dev = dev(j$m_log2, ref$m_log2, log2(s)), mu_loess_dev = dev(j$m_in, ref$m_in, log2(s)),
                mu_ctl_same = same_ctl, mu_k_same = same_k, mu_pct_same = same_pct,
                mu_ruv_dev = if (is.null(j$m_out)) NA_real_ else dev(centre(j$m_out), centre(ref$m_out), 0),
                mu_logfc_dev = da_dev("logFC"), mu_t_dev = da_dev("t"), mu_call_changes = calls_changed, mu_t_dev_notrend = notrend)
            for (i in seq_along(res)) res[[i]][intersect(names(mu), ROW_COLUMNS)] <- mu[intersect(names(mu), ROW_COLUMNS)]
            rows <- c(rows, res)
        }
    }
    rows
}

runConditionP <- function(sim, mask, job, d, lifecycle, fault_now, config) {
    rows <- list()
    common <- list(part = "P", cell = sim$cell, job = job, scenario = sim$scenario, mask = mask, rep = sim$r)
    values <- maskedValues(sim, mask)
    common$miss_share <- mean(is.na(values[sim$feat_ids, ]))
    requests <- REQUESTS_P(d)
    journeys <- list()
    T_ids <- NULL
    or_errs <- NULL
    condition <- paste(sim$scenario, mask)
    arms <- ARMS_OF$P[[job]]
    # VAL-132 (2): NOn (the design none) runs in full only when the Skip/none gate failed on the bound commit
    if (identical(job, "deciding") && identical(Sys.getenv("K188_NON"), "1")) arms <- c(setdiff(arms, "OR"), "NOn", "OR")
    for (arm in arms) {
        fault <- isTRUE(fault_now) && identical(arm, arms[[1L]])
        if (identical(arm, "OR")) {
            if (is.null(T_ids)) next
            j <- journeyP(sim, maskedValues(sim, mask, oracle = TRUE), "OR", list(gate = FALSE, norm = c(requests$NO$norm[setdiff(names(requests$NO$norm),
                "norm_method")], list(norm_method = "none"))), d, tempfile("k188-p-"), config, filter = FALSE, correlation = FALSE, ids = T_ids)
            res <- armRows(sim, j, T_ids, common, NULL)
            or_errs <- res$errs
            rows <- c(rows, res$rows)
            journeys$OR <- j
            next
        }
        j <- journeyP(sim, values, arm, requests[[if (identical(arm, "NOn")) "NOnone" else arm]], d, tempfile("k188-p-"), config, fault = fault)
        v <- j$ruv_value
        ch <- list(proteins = list(k = suppressWarnings(as.integer(v[["best_k"]] %||% NA)), pct = suppressWarnings(as.numeric(v[["percentage"]] %||% NA)),
            controls = if (is.logical(v$controls)) names(v$controls)[which(v$controls)] else character()))
        if (!is.null(j$m_in)) {
            m <- j$m_in$proteins
            ch$proteins$complete <- sum(rowSums(!is.finite(m[intersect(ch$proteins$controls, rownames(m)), , drop = FALSE])) == 0L)
        }
        if (is.null(T_ids) && !is.null(j$kept)) T_ids <- j$kept
        journeys[[arm]] <- j
        new_rows <- armRows(sim, j, T_ids %||% character(), common, ch)$rows
        # 18.4a (VAL-126): ties and row order, per legacy journey (LA with qvalue, LM with BH)
        tie <- NULL
        if (arm %in% c("LA", "LP", "LM") && identical(j$status, "ok")) {
            tie <- tryCatch(tieCount(sim, j, suppressWarnings(as.numeric(v[["percentage"]] %||% d$manual_pct)), if (arm == "LM") "BH" else "qvalue"),
                error = function(e) list(error = oneLine(conditionMessage(e))))
            if (!is.null(tie$error)) {
                addCheck("TIE", sim$cell, sim$r, paste(condition, arm), NA, "computed", FALSE, tie$error)
                tie <- NULL
            } else for (i in seq_along(new_rows)) new_rows[[i]][c("tie_share", "tie_block", "tie_block_onegroup", "perm_jaccard", "tie_rule_match")] <-
                tie[c("tie_share", "tie_block", "tie_block_onegroup", "perm_jaccard", "tie_rule_match")]
            if (lifecycle && !is.null(tie)) addCheck("TIE", sim$cell, sim$r, paste(condition, arm), NA, "computed; replica equals the selection",
                isTRUE(tie$tie_rule_match), "values blinded")
        }
        rows <- c(rows, new_rows)
        if (lifecycle && identical(arm, "LA") && identical(j$status, "ok")) checkE9P(sim, j, condition, tie)
        if (lifecycle && !is.null(j$transform) && identical(arm, "LA")) addCheck("TRANSFORM", sim$cell, sim$r, condition, j$transform,
            "the bound commit's design step", TRUE)
        if (lifecycle && !is.null(j$e1) && (!isTRUE(j$e1$pass) || identical(arm, "LA"))) {
            addCheck("E1", sim$cell, sim$r, paste(condition, arm), NA, "every step as declared", j$e1$pass, j$e1$note)
        }
        if (identical(arm, "LA") && identical(mask, "none") && identical(j$status, "ok") && !is.na(ch$proteins$k) && identical(job, "deciding")) {
            amp <- checkE8(sim, j, ch, condition)
            for (i in seq_along(rows)) if (identical(rows[[i]]$arm, "LA") && !is.null(byName(amp, "proteins"))) rows[[i]][names(byName(amp, "proteins"))] <- byName(amp, "proteins")
        }
        if (identical(arm, "LA") && identical(job, "deciding") && !is.na(ch$proteins$k)) {
            cf <- cfBranch(sim, j, ch)
            for (i in seq_along(rows)) if (identical(rows[[i]]$arm, "LA") && !is.null(byName(cf, "proteins"))) rows[[i]][names(byName(cf, "proteins"))] <- byName(cf, "proteins")
        }
    }
    if (!is.null(or_errs)) {
        for (i in seq_along(rows)) {
            r <- rows[[i]]
            if (identical(r$arm, "OR") || is.na(r$assay) || is.null(journeys[[r$arm]]$da)) next
            fam <- journeys[[r$arm]]$da
            fam$group <- groupOfComparison(fam$comparison)
            fam <- fam[fam$group %in% r$contrast, , drop = FALSE]
            sc <- scoreFamily(fam, T_ids, stats::setNames(sim$theta[T_ids, r$contrast], T_ids), or_errs[[paste("T", r$assay, r$contrast)]])
            rows[[i]][c("n_common", "s5", "gross1_common", "or_gross1_common")] <- sc[c("n_common", "s5", "gross1_common", "or_gross1_common")]
        }
    }
    if (lifecycle && identical(job, "deciding") && identical(sim$scenario, "E1") && identical(mask, "vd20") && !is.null(journeys$NO$da)) {
        j2 <- journeyP(sim, values, "NOnone", requests$NOnone, d, tempfile("k188-p-"), config)
        if (is.null(j2$da)) {
            addCheck("GATE_skip_none", sim$cell, sim$r, condition, NA, "identical or at most 1e-10", FALSE, sprintf("design none: %s at %s: %s",
                j2$status, j2$step, j2$class))
        } else {
            key <- function(x) x[order(x$comparison, x$feature), c("comparison", "feature", "logFC", "t", "p", "q")]
            a <- key(journeys$NO$da)
            b <- key(j2$da)
            rownames(a) <- rownames(b) <- NULL
            diff <- if (nrow(a) == nrow(b)) max(abs(as.matrix(a[, 3:6]) - as.matrix(b[, 3:6])), na.rm = TRUE) else NA
            # VAL-132 (2): the pre-declared gate of the declaration (identical, or within 1e-10)
            addCheck("GATE_skip_none", sim$cell, sim$r, condition, signif(diff, 4), "identical or at most 1e-10",
                identical(a, b) || isTRUE(nrow(a) == nrow(b) && identical(a[, 1:2], b[, 1:2]) && diff <= 1e-10),
                sprintf("rows %d and %d; identical %s; if it fails anywhere NOn runs in full (K188_NON=1)", nrow(a), nrow(b), identical(a, b)))
        }
    }
    list(rows = rows, journeys = journeys)
}

#' 5.8b: KIT-172's concentration fixture (`kit172ConcentrationValues()`) as a data set of one class: every feature a true
#' null, no unwanted part, the oracle the same values without ITSD, normalisation or RUV; ITSD on its `IS_` standards in
#' both classes (as T11 and KIT-190 ran it).
k172cSim <- function(class) {
    values <- kit172ConcentrationValues()
    design <- kit172Design()
    ids <- rownames(values)
    is_ids <- ids[startsWith(ids, "IS_")]
    feat <- setdiff(ids, is_ids)
    y <- log2(values)
    p <- list(idx = 0L, part = "L", omic = class, assays = "LCMS_Pos", per_assay = length(feat), groups = 2L, per = 6L, blocks = 0L)
    list(cell = paste0("K172c-", class), scenario = "N0", r = 0L, p = p, itsd = TRUE, groups = c("A", "B"), group = design$group,
        block = rep(1L, nrow(design)), runs = design$Run, design = design, feat_ids = feat, lbl_ids = character(), is_ids = is_ids,
        assay_of = stats::setNames(rep("LCMS_Pos", length(ids)), ids), y = y, oracle = y[feat, , drop = FALSE],
        unwanted = 0 * y[feat, , drop = FALSE], l = rep(0, nrow(design)),
        theta = matrix(0, length(feat), 1L, dimnames = list(feat, "B")), planted = stats::setNames(rep(FALSE, length(feat)), feat),
        U = matrix(1, length(feat), nrow(design), dimnames = list(feat, design$Run)), order_of = list(LCMS_Pos = ids), C = 0L)
}

#' E7 (refusal): the manual legacy journey with k one above its controls observed in every sample.
checkE7Refusal <- function(sim, d, config = NULL) {
    values <- maskedValues(sim, "none")
    if (identical(sim$p$part, "P")) {
        requests <- REQUESTS_P(d)
        lm <- journeyP(sim, values, "LM", requests$LM, d, tempfile("k188-p-"), config)
        k <- (lm$ruv_value %||% list())
        ctl <- if (is.logical(k$controls)) names(k$controls)[which(k$controls)] else character()
        complete <- if (is.null(lm$m_in)) NA_integer_ else sum(rowSums(!is.finite(lm$m_in$proteins[intersect(ctl, rownames(lm$m_in$proteins)), ,
            drop = FALSE])) == 0L)
        request <- requests$LM
        request$norm$ruv_k <- complete + 1L
        j <- journeyP(sim, values, "LM", request, d, tempfile("k188-p-"), config)
        want <- "prot_ruv_insufficient_controls"
    } else {
        requests <- REQUESTS_L(d)
        lm <- journeyL(sim, values, "LM", requests$LM, d, tempfile("k188-e7-"))
        ch <- choiceOf(lm, TRUE)
        complete <- if (is.null(ch)) NA_integer_ else max(vapply(ch, `[[`, integer(1L), "complete"))
        request <- requests$LM
        request$ruv$manual_k <- as.integer(complete + 1L)
        j <- journeyL(sim, values, "LM", request, d, tempfile("k188-e7-"))
        want <- "omics_ruv_insufficient_controls"
        # r1.21 (R20-P07, 7.4 E7): one LM journey with an invalid typed k (-1) stops by rule (VAL-242: an invalid argument stops);
        # its class is recorded as read at the bound commit, and compared with what lifecycle.md wrote (K188_PREDICTIONS)
        bad <- requests$LM
        bad$ruv$manual_k <- -1L
        jb <- journeyL(sim, values, "LM", bad, d, tempfile("k188-e7-"))
        observed <- c(status = jb$status, class = as.character(jb$class %||% NA), step = as.character(jb$step %||% NA))
        written <- predictionsMet("E7_invalid_k", observed)
        addCheck("E7_invalid_k", sim$cell, sim$r, "N0 none", jb$status, "a stop (refused or crashed), class recorded",
            jb$status %in% c("refused", "crashed") && written$n > 0L && !length(written$unmet),
            sprintf("k -1: status %s, class %s, step %s; %s", jb$status, jb$class, jb$step, paste(written$unmet, collapse = "; ")))
        # the r1.20 plant (k one above the complete controls) is kept as an observation of Q6.2 (c): what it meets is recorded
        # and compared with lifecycle.md's text only where that wrote one
        observed <- c(status = j$status, class = as.character(j$class %||% NA), step = as.character(j$step %||% NA))
        written <- predictionsMet("E7_refusal", observed)
        return(addCheck("E7_refusal", sim$cell, sim$r, "N0 none", complete + 1L, "observation (Q6.2 (c)); as written", written$n > 0L && !length(written$unmet),
            sprintf("k %s: status %s, class %s, step %s; %s", complete + 1L, j$status, j$class, j$step, paste(written$unmet, collapse = "; "))))
    }
    addCheck("E7_refusal", sim$cell, sim$r, "N0 none", complete + 1L, want, identical(j$status, "refused") && identical(j$class, want),
        sprintf("k %s: status %s, class %s, step %s", complete + 1L, j$status, j$class, j$step))
}

# ---- E12's format-default pattern assertion (35.2 (2), 5.2; once per job, defaults read from the bound commit) ------------
#' Every format default the import proposes matches exactly the `IS_` rows of the generated names (LipidSearch's searches the
#' class column, 5.2), through the matcher the ITSD step uses (stringr::str_detect).
checkE12Patterns <- function() {
    names_all <- c(sprintf("F%04d", 1:250), sprintf("LBL_%03d", 1:10), sprintf("CTL_%03d", 1:5), sprintf("IS_%02d", 1:5))
    class_all <- ifelse(startsWith(names_all, "IS_"), "ISTD", "PC")
    for (omic in c("lipid", "metab")) for (format in c("msdial", "progenesis", "xcms", "compound_discoverer", "lipidsearch")) {
        defaults <- tryCatch(if (identical(omic, "lipid")) getLipidomicsColumnDefaults(format) else getMetabolomicsColumnDefaults(format),
            error = function(e) NULL)
        pattern <- defaults$is_pattern
        if (is.null(pattern) || is.na(pattern)) next
        hit <- stringr::str_detect(if (identical(format, "lipidsearch")) class_all else names_all, pattern)
        addCheck("E12_patterns", paste(omic, format), 0L, "format default", pattern, "matches exactly the IS_ rows",
            identical(which(hit), which(startsWith(names_all, "IS_"))), sprintf("%d rows matched", sum(hit)))
    }
}

# ---- L1c: the checks fire on planted failures and stay silent otherwise (35.6; 40.3 (6)) -----------------------------------
# K188_PLANT (comma list; default all) with K188_JOBS = l1c, K188_LIFECYCLE = 1, cell L12, one data set (the first). One unit
# per plant, each in its own child process. The plants act on the generator's input or on the script's comparison inputs
# (itsd, itsdpos, emptysample, ruvdrop: the table the route is given; ruvskip: the request's percentages; normplant: the matrix
# E13 reads; daerror: a stop at step da; exitplant: the child's own exit), never on the route's state.
PLANTS <- c("itsd", "itsdpos", "emptysample", "ruvskip", "ruvdrop", "daerror", "normplant", "exitplant")
# The predictions of 35.6 that do not depend on code that was not on staging when this was written (each cited): the key is the
# journey's flattened field (`status`, `class`, `step`, or `<column>.<assay>`); "re:" a Perl regular expression, "range:a:b"
# an inclusive range, "NA" a missing value. The R-marked predictions (35.6 r1.20: KIT-242's and KIT-260's field names,
# classes and notices) are read at the bound commit, written into lifecycle.md and passed in K188_PREDICTIONS.
PLANT_EXPECT <- list(
    # 35.6 (r1.20) `itsd`: E12 fires in both assays (a, b, c), E13 silent
    itsd = c(status = "instrument", class = "re:^itsd_not_applied:LCMS_Neg:a,b,c \\| itsd_not_applied:LCMS_Pos:a,b,c$",
        itsd_ok.LCMS_Neg = "FALSE", itsd_ok.LCMS_Pos = "FALSE", itsd_is_left.LCMS_Neg = "5", itsd_is_left.LCMS_Pos = "5",
        itsd_n_used.LCMS_Neg = "0", itsd_n_used.LCMS_Pos = "0", norm_ok.LCMS_Neg = "TRUE", norm_ok.LCMS_Pos = "TRUE"),
    # `itsdpos`: E12 fires for LCMS_Pos alone
    itsdpos = c(status = "instrument", class = "itsd_not_applied:LCMS_Pos:a,b,c", itsd_ok.LCMS_Pos = "FALSE",
        itsd_ok.LCMS_Neg = "TRUE", itsd_is_left.LCMS_Pos = "5", itsd_n_used.LCMS_Pos = "0", itsd_is_left.LCMS_Neg = "0",
        itsd_n_used.LCMS_Neg = "5"),
    # `emptysample`: E13 silent in both assays, no instrument row (35.6 r1.20; its later steps are recorded, not predicted)
    emptysample = c(status = "re:^(?!instrument$)", itsd_ok.LCMS_Neg = "TRUE", itsd_ok.LCMS_Pos = "TRUE", norm_ok.LCMS_Neg = "TRUE",
        norm_ok.LCMS_Pos = "TRUE"),
    # `ruvskip`: no k in either assay, kept uncorrected (Q6.2 (c)); the reason's string starts as the trace's status counts
    ruvskip = c(status = "uncorrected", ruv_status.LCMS_Neg = "uncorrected", ruv_status.LCMS_Pos = "uncorrected",
        ruv_k_found.LCMS_Neg = "NA", ruv_k_found.LCMS_Pos = "NA", ruv_n_controls.LCMS_Neg = "0", ruv_n_controls.LCMS_Pos = "0",
        ruv_status_product.LCMS_Neg = "uncorrected", ruv_status_product.LCMS_Pos = "uncorrected",
        ruv_changed.LCMS_Neg = "FALSE", ruv_changed.LCMS_Pos = "FALSE", ruv_no_k_reason.LCMS_Neg = "re:^skipped_insufficient_controls:1",
        ruv_no_k_reason.LCMS_Pos = "re:^skipped_insufficient_controls:1"),
    # `ruvdrop`: 11 to 15 rows reach the RUV input in LCMS_Pos (derived by reading, 35.6 r1.20); LCMS_Pos kept uncorrected;
    # LCMS_Neg corrected with a k used of 1 or 2 (an observation: a different result is a defect of the plant)
    ruvdrop = c(rows_ruv_input.LCMS_Pos = "range:11:15", ruv_status.LCMS_Pos = "uncorrected", ruv_k_found.LCMS_Pos = "NA",
        ruv_status.LCMS_Neg = "corrected", ruv_k_used.LCMS_Neg = "range:1:2"),
    # `daerror`: ruvskip's RUV columns, then a crash at step da with the deliberate message
    daerror = c(status = "crashed", step = "da", da_status.LCMS_Neg = "failed", da_status.LCMS_Pos = "failed",
        da_message.LCMS_Neg = "re:K188 planted DA error", ruv_status.LCMS_Neg = "uncorrected", ruv_status.LCMS_Pos = "uncorrected",
        ruv_k_found.LCMS_Neg = "NA", ruv_n_controls.LCMS_Neg = "0", ruv_status_product.LCMS_Neg = "uncorrected",
        ruv_status_product.LCMS_Pos = "uncorrected"),
    # `normplant`: E13 fires for LCMS_Pos alone
    normplant = c(status = "instrument", class = "norm_not_applied:LCMS_Pos", norm_ok.LCMS_Pos = "FALSE", norm_ok.LCMS_Neg = "TRUE",
        itsd_ok.LCMS_Neg = "TRUE", itsd_ok.LCMS_Pos = "TRUE"))
matchExpectation <- function(observed, expected) {
    observed <- unname(observed)
    if (startsWith(expected, "re:")) return(!is.na(observed) && grepl(substring(expected, 4L), observed, perl = TRUE))
    if (startsWith(expected, "range:")) {
        bounds <- as.numeric(strsplit(substring(expected, 7L), ":", fixed = TRUE)[[1L]])
        value <- suppressWarnings(as.numeric(observed))
        return(is.finite(value) && value >= bounds[[1L]] && value <= bounds[[2L]])
    }
    if (identical(expected, "NA")) return(is.na(observed))
    !is.na(observed) && identical(as.character(observed), expected)
}
# SA.8 (tamper-evident): the file is committed and pushed before the request, passed by a path inside the checkout, and its
# SHA-256 is written into the provenance
predictionsFileOk <- function() {
    path <- Sys.getenv("K188_PREDICTIONS")
    nzchar(path) && file.exists(path) && startsWith(normalizePath(path), paste0(normalizePath(getwd()), "/"))
}
loadPredictions <- function() {
    path <- Sys.getenv("K188_PREDICTIONS")
    if (!nzchar(path)) return(NULL)
    utils::read.delim(path, stringsAsFactors = FALSE, colClasses = "character", na.strings = character())
}
#' The R-marked predictions written into lifecycle.md (a TSV with columns key = "<id>/<field>" and expected): each row's
#' field of `observed` against its expectation. Returns the count of rows written for `id` and the unmet ones.
predictionsMet <- function(id, observed) {
    written <- loadPredictions()
    if (is.null(written)) return(list(n = 0L, unmet = character()))
    mine <- written[startsWith(written$key, paste0(id, "/")), , drop = FALSE]
    unmet <- character()
    for (i in seq_len(nrow(mine))) {
        field <- sub("^[^/]*/", "", mine$key[[i]])
        if (!matchExpectation(observed[field], mine$expected[[i]])) {
            unmet <- c(unmet, sprintf("%s: observed %s, written %s", field, observed[field], mine$expected[[i]]))
        }
    }
    list(n = nrow(mine), unmet = unmet)
}
flattenJourney <- function(j, assays) {
    flat <- c(status = j$status, class = as.character(j$class %||% NA), step = as.character(j$step %||% NA),
        warning_classes = as.character(j$warning_classes %||% NA))
    for (a in assays) for (name in names(j$assay_cols[[a]])) flat[paste0(name, ".", a)] <- as.character(j$assay_cols[[a]][[name]])
    for (what in c("rows_filter", "rows_itsd", "rows_ruv_input")) for (a in names(j[[what]])) flat[paste0(what, ".", a)] <- as.character(j[[what]][[a]])
    # the bound commit's records (KIT-242, KIT-260): the analysis record's outcomes of the ITSD, normalisation and RUV steps, per assay (exact `[[`)
    one <- function(x) if (is.null(x) || !length(x)) NA_character_ else paste(as.character(unlist(x)), collapse = ",")
    for (a in assays) {
        itsd <- j[["record_itsd"]]; norm <- j[["record_norm"]]; ruv <- j[["record_ruv"]]
        flat[paste0("itsd_record_applied.", a)] <- one(itsd[["itsd_applied_per_assay"]][[a]])
        flat[paste0("itsd_record_reason.", a)] <- one(itsd[["itsd_not_applied_reason_per_assay"]][[a]])
        flat[paste0("itsd_record_count.", a)] <- one(itsd[["itsd_counts_per_assay"]][[a]])
        flat[paste0("itsd_record_ids.", a)] <- one(j[["itsd_features_record"]][[a]])
        flat[paste0("norm_record_excluded.", a)] <- one(norm[["normalisation_samples_excluded_per_assay"]][[a]])
        flat[paste0("ruv_record_applied.", a)] <- one(ruv[["ruv_applied_per_assay"]][[a]])
        flat[paste0("ruv_record_reason.", a)] <- one(ruv[["ruv_not_applied_reason_per_assay"]][[a]])
    }
    flat[["itsd_record_applied"]] <- one(j[["record_itsd"]][["itsd_applied"]])
    flat[["norm_record_removed"]] <- one(j[["record_norm"]][["normalisation_samples_removed"]])
    flat
}
runPlantUnit <- function(plant, cellname, r, d) {
    sim <- generate(cellname, "N1", r)
    values <- maskedValues(sim, "none")
    requests <- REQUESTS_L(d)
    arm <- if (identical(plant, "emptysample")) "NO" else "LA"
    request <- requests[[arm]]
    ids <- NULL
    in_assay <- function(a) names(sim$assay_of)[sim$assay_of == a]
    if (plant %in% c("ruvskip", "daerror")) {
        request$ruv$percentage_min <- 1
        request$ruv$percentage_max <- 1
    }
    if (identical(plant, "ruvdrop")) {
        pos <- in_assay("LCMS_Pos")
        ids <- c(in_assay("LCMS_Neg"), head(intersect(sim$feat_ids, pos), 10L), intersect(sim$lbl_ids, pos), intersect(sim$is_ids, pos))
    }
    if (identical(plant, "emptysample")) values[setdiff(in_assay("LCMS_Pos"), sim$is_ids), sim$runs[[1L]]] <- NA_real_
    j <- journeyL(sim, values, arm, request, d, tempfile(sprintf("k188-%s-", plant)), ids = ids, plant = plant)
    common <- list(part = "L", cell = sim$cell, job = "l1c", scenario = sim$scenario, mask = "none", rep = sim$r)
    rows <- armRows(sim, j, setdiff(j$kept %||% character(), sim$is_ids), common, choiceOf(j, arm %in% c("LA", "LM")))$rows
    flat <- flattenJourney(j, sim$p$assays)
    expected <- PLANT_EXPECT[[plant]]
    unmet <- names(expected)[!vapply(names(expected), function(f) matchExpectation(flat[f], expected[[f]]), logical(1L))]
    addCheck("E15", sim$cell, r, plant, if (length(unmet)) "unmet" else "met", "as 35.6 predicts", !length(unmet),
        paste(sprintf("%s: observed %s, predicted %s", unmet, flat[unmet], expected[unmet]), collapse = "; "))
    written <- predictionsMet(plant, flat)
    if (written$n) addCheck("E15R", sim$cell, r, plant, written$n, "as written in lifecycle.md (R items)", !length(written$unmet),
        paste(written$unmet, collapse = "; "))
    list(rows = do.call(rbind, lapply(rows, function(x) as.data.frame(x[ROW_COLUMNS], stringsAsFactors = FALSE))),
        checks = do.call(rbind, checks_env$rows), plant = data.frame(plant = plant, key = names(flat), value = unname(flat), stringsAsFactors = FALSE),
        seconds = NA_real_)
}

#' The guard's planted units (the metabolomics cell, one ordinary journey each): `guardjourney` raises a `$` prefix match inside the journey, `guardunit`
#' inside the unit outside any journey, `guardcontrol` is an ordinary journey. The journey and the unit must end as `instrument` and the control as `ok`.
runGuardUnit <- function(variant, cellname, r, d) {
    if (identical(variant, "guardunit")) {
        probe <- list(partial_probe_value = 1)
        probe$partial_probe
    }
    sim <- generate(cellname, "N0", r)
    values <- maskedValues(sim, "none")
    j <- journeyL(sim, values, "NO", REQUESTS_L(d)[["NO"]], d, tempfile("k188-guard-"), plant = variant)
    common <- list(part = "L", cell = sim$cell, job = "guard", scenario = sim$scenario, mask = "none", rep = sim$r)
    rows <- armRows(sim, j, setdiff(j[["kept"]], sim$is_ids), common, NULL)$rows
    guarded <- identical(variant, "guardjourney")
    ok <- if (guarded) identical(j[["status"]], "instrument") && startsWith(as.character(j[["class"]]), "partial_match_guard:") else identical(j[["status"]], "ok")
    addCheck("GUARD_planted", sim$cell, r, variant, j[["status"]], if (guarded) "instrument" else "ok", ok, as.character(j[["class"]]))
    list(rows = do.call(rbind, lapply(rows, function(x) as.data.frame(x[ROW_COLUMNS], stringsAsFactors = FALSE))), checks = do.call(rbind, checks_env$rows), seconds = NA_real_)
}
runUnit <- function(part, cellname, r, job, scenario, mask, lifecycle, fault_now, plant = "") {
    t_unit <- proc.time()[[3L]]
    omic <- CELLS[[cellname]]$omic
    d <- routeDefaults(omic)
    if (identical(job, "guard")) return(runGuardUnit(plant, cellname, r, d))
    if (identical(job, "l1c")) {
        # 35.6 `exitplant` (R20-P10): this unit's child ends its process before writing rows; the parent writes the unit row
        if (identical(plant, "exitplant")) quit(save = "no", status = 3L)
        return(runPlantUnit(plant, cellname, r, d))
    }
    config <- if (identical(part, "P")) pConfig() else NULL
    conditions <- conditionsOf(part, cellname, job)
    if (!identical(scenario, "all")) conditions <- conditions[conditions$scenario == scenario & conditions$mask == mask, , drop = FALSE]
    rows <- list()
    if (identical(cellname, "K172c")) {
        if (fault_now) options(k188.fault_noda = TRUE)
        for (class in c("lipid", "metab")) {
            sim <- k172cSim(class)
            rows <- c(rows, runConditionL(sim, "none", job, routeDefaults(class), lifecycle, FALSE, TRUE)$rows)
        }
        conditions <- conditions[0, , drop = FALSE]
    }
    if (fault_now && identical(part, "L")) options(k188.fault_noda = TRUE)
    for (s in unique(conditions$scenario)) {
        sim <- generate(cellname, s, r)
        if (lifecycle) {
            e6 <- max(abs(sim$oracle + sim$unwanted + matrix(sim$l, nrow(sim$oracle), length(sim$runs), byrow = TRUE) -
                sim$y[rownames(sim$oracle), , drop = FALSE]))
            addCheck("E6", cellname, r, s, signif(e6, 4), 1e-12, e6 <= 1e-12)
        }
        for (m in conditions$mask[conditions$scenario == s]) {
            first <- !length(rows)
            if (identical(job, "munit")) {
                rows <- c(rows, runConditionMunit(sim, m, d))
            } else if (identical(part, "P")) {
                rows <- c(rows, runConditionP(sim, m, job, d, lifecycle, fault_now && first, config)$rows)
            } else {
                rows <- c(rows, runConditionL(sim, m, job, d, lifecycle, fault_now && first, first)$rows)
            }
        }
    }
    if (lifecycle && job %in% c("reported") && (identical(scenario, "all") || (identical(scenario, "N0") && identical(mask, "none"))) &&
        isTRUE(getOption("k188.e7_unit", FALSE))) {
        checkE7Refusal(generate(cellname, "N0", r), d, config)
    }
    table <- do.call(rbind, lapply(rows, function(x) as.data.frame(x[ROW_COLUMNS], stringsAsFactors = FALSE)))
    list(rows = table, checks = if (length(checks_env$rows)) do.call(rbind, checks_env$rows),
        seconds = proc.time()[[3L]] - t_unit)
}

# ---- child process ---------------------------------------------------------------------------------------------------
if (child) {
    stopifnot(length(args) == 8L)
    lifecycle <- identical(Sys.getenv("K188_LIFECYCLE"), "1")
    options(k188.e7_unit = identical(Sys.getenv("K188_E7_UNIT"), "1"))
    fault_now <- identical(Sys.getenv("K188_FAULT_UNIT"), "1")
    options(k188.progress = sub("\\.rds$", "-progress.tsv", args[[8L]]))
    # an error that ends the unit also carries the call stack at the point of the error (deparsed, one line each, truncated), so the
    # failing line is located (coordinator's ruling, defect 2)
    stack <- NULL
    # the planted `crashnet` unit: a non-guard error raised OUTSIDE every handler, so only `options(error=)` sees it (the trace writer; the parent records `crashed`)
    if (identical(Sys.getenv("K188_PLANT_UNIT"), "crashnet")) stop("K188 planted non-guard error outside every handler")
    result <- tryCatch(withCallingHandlers(runUnit(args[[2L]], args[[3L]], as.integer(args[[4L]]), args[[5L]], args[[6L]], args[[7L]], lifecycle, fault_now,
        Sys.getenv("K188_PLANT_UNIT")), warning = guardWarning, error = function(e) {
        stack <<- vapply(sys.calls(), function(cl) substr(oneLine(paste(deparse(cl, width.cutoff = 200L), collapse = " "), 200L), 1L, 200L), "")
    }),
        error = function(e) list(error = oneLine(paste(class(e)[[1L]], conditionMessage(e))), stack = stack,
            checks = if (length(checks_env$rows)) do.call(rbind, checks_env$rows)))
    saveRDS(result, args[[8L]])
    quit(save = "no", status = 0L)
}

# ---- the job (parent) ------------------------------------------------------------------------------------------------
if (length(args) != 5L) stop("usage: Rscript sim-k188.R <part> <cells> <first> <last> <out dir>")
part <- args[[1L]]
cells_wanted <- strsplit(args[[2L]], ",", fixed = TRUE)[[1L]]
first_rep <- as.integer(args[[3L]])
last_rep <- as.integer(args[[4L]])
out_dir <- args[[5L]]
stopifnot(part %in% c("L", "P"), all(cells_wanted %in% names(CELLS)), all(vapply(CELLS[cells_wanted], `[[`, "", "part") == part),
    !is.na(first_rep), !is.na(last_rep), first_rep <= last_rep)
lifecycle <- identical(Sys.getenv("K188_LIFECYCLE"), "1")
# the check SMOKE_counted_refusal: K188_SMOKE=1 marks an instrument check; a counted request (not a lifecycle run) that sets it is refused
# before any unit runs. The flag is in every output: the rows' `smoke` column and the provenance; score-k188.R refuses such rows.
smoke <- identical(Sys.getenv("K188_SMOKE"), "1")
if (smoke && !lifecycle) stop("SMOKE_counted_refusal: K188_SMOKE=1 is set on a counted (non-lifecycle) request; refused")
fault <- identical(Sys.getenv("K188_FAULT"), "1")
workers <- as.integer(Sys.getenv("K188_WORKERS", "4"))
jobs <- strsplit(Sys.getenv("K188_JOBS", if (identical(part, "L")) "deciding,reported,munit,ba" else "deciding,reported,ba"), ",",
    fixed = TRUE)[[1L]]
job_name <- Sys.getenv("K188_JOB", sprintf("%s-%s-%d-%d%s", tolower(part), paste(cells_wanted, collapse = "_"), first_rep, last_rep,
    if (lifecycle) "-l1" else ""))
dir.create(file.path(out_dir, "rows"), recursive = TRUE, showWarnings = FALSE)
out_dir <- normalizePath(out_dir)
p_last <- {
    spec <- strsplit(strsplit(Sys.getenv("K188_P_REPS", "P0=403,P2=403"), ",", fixed = TRUE)[[1L]], "=", fixed = TRUE)
    stats::setNames(as.integer(vapply(spec, `[[`, "", 2L)), vapply(spec, `[[`, "", 1L))
}
# 35.6: K188_PLANT (comma list, default all eight) names the planted cases of the L1c job (K188_JOBS = l1c)
plants_wanted <- strsplit(Sys.getenv("K188_PLANT", paste(PLANTS, collapse = ",")), ",", fixed = TRUE)[[1L]]
# K188_GUARD (comma list, default all three): the planted guard units of the job `guard` (the metabolomics cell M12)
guards_wanted <- strsplit(Sys.getenv("K188_GUARD", "guardjourney,guardunit,guardcontrol,crashnet"), ",", fixed = TRUE)[[1L]]
if ("guard" %in% jobs) stopifnot(lifecycle, identical(part, "L"), all(guards_wanted %in% c("guardjourney", "guardunit", "guardcontrol", "crashnet")))
if ("l1c" %in% jobs) {
    stopifnot(lifecycle, identical(part, "L"), all(plants_wanted %in% PLANTS))
    # 40.3 (8), 35.6 (r1.20): the R-marked predictions are read at the bound commit and written into lifecycle.md before the
    # run; this stops a run sent without them
    if (!predictionsFileOk()) {
        stop("K188_PREDICTIONS (an absolute path to the TSV of the R-marked predictions written into lifecycle.md) is not set or not found")
    }
}
units <- list()
for (cellname in cells_wanted) for (job in jobs) {
    if (!jobApplies(cellname, job)) next
    if (identical(job, "guard")) {
        for (plant in guards_wanted) units[[length(units) + 1L]] <- list(cell = cellname, r = first_rep, job = job, scenario = "N0", mask = "none", plant = plant)
        next
    }
    if (identical(job, "l1c")) {
        for (plant in plants_wanted) units[[length(units) + 1L]] <- list(cell = cellname, r = first_rep, job = job, scenario = "N1",
            mask = "none", plant = plant)
        next
    }
    reps <- seq(first_rep, last_rep)
    if (identical(part, "P") && lifecycle) reps <- reps[reps <= (if (cellname %in% names(p_last)) p_last[[cellname]] else first_rep)]
    for (r in reps) {
        if (identical(part, "P")) {
            cond <- conditionsOf(part, cellname, job)
            for (i in seq_len(nrow(cond))) units[[length(units) + 1L]] <- list(cell = cellname, r = r, job = job, scenario = cond$scenario[[i]],
                mask = cond$mask[[i]], plant = "")
        } else units[[length(units) + 1L]] <- list(cell = cellname, r = r, job = job, scenario = "all", mask = "all", plant = "")
    }
}
units <- do.call(rbind, lapply(units, as.data.frame, stringsAsFactors = FALSE))
# largest first (the four-way packing), the deliberate fault and the refusal check in fixed units
weight <- c(L12 = 1, L24 = 1.5, L64 = 2.5, M12 = 0.5, M168 = 3, P0 = 1, P1 = 1.3, P16 = 2, P2 = 0.5, P3 = 2.6, P4b = 5)[units$cell] *
    c(deciding = 1, reported = 1.8, munit = 1.1, ba = 0.8, lp = 0.7, las = 0.7, l1c = 0.3, guard = 0.1)[units$job]
units <- units[order(-weight), , drop = FALSE]
units$fault <- FALSE
units$e7 <- FALSE
if (lifecycle && fault) units$fault[[1L]] <- TRUE
if (lifecycle) for (cellname in cells_wanted) {
    hit <- which(units$cell == cellname & units$job == "reported" & units$r == first_rep &
        (units$scenario == "all" | (units$scenario == "N0" & units$mask == "none")))
    if (length(hit)) units$e7[[hit[[1L]]]] <- TRUE
}
# SA-P10: E7's R-marked keys are required whenever a unit runs E7 (a lifecycle run of a reported job)
if (lifecycle && identical(part, "L") && any(units$e7) && !predictionsFileOk()) stop("K188_PREDICTIONS (a path inside the checkout) is required: E7's R-marked keys are compared with it")
writeProvenance <- function(final) {
    tables <- if (final) list.files(out_dir, recursive = TRUE, full.names = TRUE, pattern = "\\.(tsv|gz)$")
    library_table <- utils::installed.packages()[, c("Package", "Version")]
    key_packages <- c("MultiScholaR", "RUVIIIC", "ruv", "limma", "qvalue", "sva", "duckdb", "dplyr")
    provenance <- list(job = job_name, part = part, cells = cells_wanted, data_sets = c(first_rep, last_rep), jobs = jobs, lifecycle = lifecycle,
        fault = fault, workers = workers, units = nrow(units),
        commit = tryCatch(system2("git", c("rev-parse", "HEAD"), stdout = TRUE), error = function(e) NA),
        r_version = R.version.string, blas = extSoftVersion()[["BLAS"]], lapack = La_library(), omp_threads = Sys.getenv("OMP_NUM_THREADS"),
        runner_image = paste(Sys.getenv("ImageOS"), Sys.getenv("ImageVersion")),
        processor = tryCatch(sub("^model name\\s*:\\s*", "", grep("^model name", readLines("/proc/cpuinfo"), value = TRUE)[[1L]]), error = function(e) NA),
        cores = parallel::detectCores(),
        predictions_sha256 = if (nzchar(Sys.getenv("K188_PREDICTIONS")) && file.exists(Sys.getenv("K188_PREDICTIONS"))) sha(Sys.getenv("K188_PREDICTIONS")) else NA,
        smoke = smoke,
        prot_design_transform = if (exists("protDesignPlainLog2", mode = "function")) "log2(x); 0 missing" else "log2(x + 1)",
        seed_rule = "set.seed(188000000 + 100000 c + 10000 s + r, kind = 'Mersenne-Twister', normal.kind = 'Inversion', sample.kind = 'Rejection')",
        key_packages = as.list(stats::setNames(library_table[match(key_packages, library_table[, "Package"]), "Version"], key_packages)),
        script_sha256 = sha(script_file), started = format(started_at, "%Y-%m-%dT%H:%M:%S%z"), finished = if (final) format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z") else NA,
        wall_seconds = round(as.numeric(difftime(Sys.time(), started_at, units = "secs")), 1),
        units_with_error = if (final) sum(nzchar(timing$error)) else NA, exit_status = if (final) 0L else NA,
        tables_sha256 = if (final) stats::setNames(as.list(vapply(tables, sha, character(1L))), sub(paste0("^", out_dir, "/?"), "", tables)) else list())
    jsonlite::write_json(provenance, file.path(out_dir, sprintf("provenance-%s.json", job_name)), auto_unbox = TRUE, pretty = TRUE, digits = NA)
}

cat(sprintf("KIT-188 %s: %d units, %d workers, lifecycle %s, jobs %s\n", job_name, nrow(units), workers, lifecycle, paste(jobs, collapse = ",")))
# E12's format-default pattern assertion, once per job (35.2 (2))
if (identical(part, "L")) checkE12Patterns()
if (lifecycle) {
    checkE4()
    for (omic in unique(vapply(CELLS[cells_wanted], `[[`, "", "omic"))) {
        e3 <- checkE3(omic)
        addCheck("E3", omic, 0L, "defaults", NA, "equal to 2.1 and 18.4", e3$pass, paste(e3$note, "|",
            paste(names(e3$got), vapply(e3$got, function(x) paste(format(x), collapse = "/"), ""), sep = "=", collapse = "; ")))
    }
}
# written before any unit runs, so a job that is killed still leaves its provenance and the parent's checks (defect 3)
writeProvenance(FALSE)
if (length(checks_env$rows)) utils::write.table(do.call(rbind, checks_env$rows), file.path(out_dir, sprintf("checks-%s-start.tsv", job_name)),
    sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
runChild <- function(i) {
    u <- units[i, ]
    # SA-P01: the planted case's name is part of every file name, so the eight L1c units (one data set, one cell) stay apart
    suffix <- if (nzchar(u$plant)) paste0("-", u$plant) else ""
    rds <- file.path(out_dir, "rows", sprintf("unit-%s-%s-%d-%s-%s%s.rds", u$cell, u$job, u$r, u$scenario, u$mask, suffix))
    env <- c(sprintf("K188_LIFECYCLE=%s", if (lifecycle) "1" else "0"), sprintf("K188_FAULT_UNIT=%s", if (u$fault) "1" else "0"),
        sprintf("K188_E7_UNIT=%s", if (u$e7) "1" else "0"), sprintf("K188_PLANT_UNIT=%s", u$plant))
    t0 <- Sys.time()
    code <- system2("Rscript", c(shQuote(script_file), "--unit", part, u$cell, u$r, u$job, u$scenario, u$mask, shQuote(rds)),
        env = env, stdout = FALSE, stderr = FALSE)
    seconds <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
    trace_file <- sub("\\.rds$", "-trace.txt", rds)
    res <- if (file.exists(rds)) readRDS(rds) else list(error = paste(sprintf("child process exit %d, no result", code),
        if (file.exists(trace_file)) oneLine(paste(readLines(trace_file, n = 1L), collapse = " "), 200L) else ""))
    unlink(rds)
    res$wall <- seconds
    res$exit <- code
    # written as the unit completes, so that a job stopped at its timeout still uploads what it finished (7.1)
    tag <- sprintf("%s-%s-%d-%s-%s%s", u$cell, u$job, u$r, u$scenario, u$mask, suffix)
    # a unit with NO result file and NO trace (killed by the runner or the job's limit, a signal, out of memory) cannot be told from a failure of
    # the analysis route: it is recorded as an `instrument` unit row, never `crashed` (the coordinator's ruling; the scorer stops on it, and
    # whether it is re-sent or reclassified `crashed` with its cause is decided on the evidence, never by this script). The planted `exitplant`
    # ends so by design and stays `crashed`.
    no_result_no_trace <- !file.exists(rds) && !file.exists(trace_file) && !identical(u$plant, "exitplant")
    if (no_result_no_trace) res$error <- sprintf("no_result_no_trace: child exit %d%s", code, if (code > 128L) sprintf(" (signal %d)", code - 128L) else "")
    rows <- res$rows
    if (!is.null(res$stack)) writeLines(res$stack, file.path(out_dir, "rows", sprintf("unit-trace-%s.txt", tag)))
    if (is.null(rows) || !nrow(rows)) {
        rows <- as.data.frame(do.call(newRow, list(part = part, cell = u$cell, job = u$job, scenario = u$scenario, mask = u$mask, rep = u$r,
            arm = "unit", status = if (no_result_no_trace || grepl("partial-match guard", res$error %||% "", fixed = TRUE)) "instrument" else "crashed",
            class = res$error %||% "no rows")), stringsAsFactors = FALSE)
    }
    rows$smoke <- if (smoke) "1" else NA_character_
    numeric_columns <- vapply(rows, is.numeric, logical(1L))
    rows[numeric_columns] <- lapply(rows[numeric_columns], signif, 7L)
    hash <- NA_character_
    if (lifecycle) {
        temp <- tempfile(fileext = ".tsv")
        utils::write.table(rows, temp, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
        hash <- sha(temp)
        unlink(temp)
        rows[BLIND] <- NA
    }
    utils::write.table(rows, file.path(out_dir, "rows", sprintf("part-%s.tsv", tag)), sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
    timing <- data.frame(cell = u$cell, job = u$job, rep = u$r, scenario = u$scenario, mask = u$mask, wall_seconds = round(res$wall %||% NA, 1),
        unit_seconds = round(res$seconds %||% NA, 1), exit = res$exit %||% NA, error = res$error %||% "", rows = nrow(rows), rows_sha256 = hash,
        fault = u$fault, e7 = u$e7, stringsAsFactors = FALSE)
    utils::write.table(timing, file.path(out_dir, "rows", sprintf("unit-timing-%s.tsv", tag)), sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
    if (!is.null(res$checks)) utils::write.table(res$checks, file.path(out_dir, "rows", sprintf("unit-checks-%s.tsv", tag)), sep = "\t",
        quote = FALSE, row.names = FALSE, na = "NA")
    journeys <- unique(rows[c("cell", "job", "rep", "arm", "arm_scale", "scenario", "mask", "status", "class")])
    journeys[] <- lapply(journeys, function(x) ifelse(is.na(x), "", as.character(x)))
    status <- stats::aggregate(list(n = rep(1L, nrow(journeys))), journeys[c("cell", "job", "arm", "arm_scale", "scenario", "mask",
        "status", "class")], length)
    # 35.6 `exitplant` (R20-P10): the parent wrote one `arm = unit` row, crashed, for a child that exited with status 3
    if (identical(u$plant, "exitplant")) {
        res$checks <- data.frame(check = "E15", cell = u$cell, rep = u$r, condition = "exitplant", value = as.character(code),
            threshold = "exit 3, one unit row, crashed", pass = identical(as.integer(code), 3L) && nrow(rows) == 1L &&
                identical(as.character(rows$arm), "unit") && identical(as.character(rows$status), "crashed"),
            note = sprintf("exit %s; %d row(s)", code, nrow(rows)), stringsAsFactors = FALSE)
    }
    # the planted unit-level guard: the unit ended outside any journey; its row must be `instrument` with the guard's text
    if (identical(u$plant, "guardunit")) {
        res$checks <- data.frame(check = "GUARD_planted", cell = u$cell, rep = u$r, condition = "guardunit", value = as.character(rows$status[[1L]]),
            threshold = "instrument", pass = nrow(rows) == 1L && identical(as.character(rows$arm), "unit") && identical(as.character(rows$status), "instrument") &&
                grepl("partial-match guard", as.character(rows$class), fixed = TRUE), note = oneLine(as.character(rows$class)), stringsAsFactors = FALSE)
    }
    if (identical(u$plant, "crashnet")) {
        res$checks <- data.frame(check = "GUARD_planted", cell = u$cell, rep = u$r, condition = "crashnet", value = as.character(rows$status[[1L]]),
            threshold = "crashed, with a trace", pass = nrow(rows) == 1L && identical(as.character(rows$arm), "unit") && identical(as.character(rows$status), "crashed") &&
                grepl("K188 planted non-guard error", as.character(rows$class), fixed = TRUE) && file.exists(trace_file), note = oneLine(as.character(rows$class)),
            stringsAsFactors = FALSE)
    }
    silence <- if (!identical(u$job, "l1c") && !identical(u$job, "guard") && !is.null(rows$status)) data.frame(n_instrument = sum(rows$status %in% "instrument"),
        n_itsd_false = sum(as.character(rows$itsd_ok) %in% "FALSE"), n_norm_false = sum(as.character(rows$norm_ok) %in% "FALSE"))
    list(timing = timing, checks = res$checks, status = status, plant = res$plant, silence = silence)
}
results <- parallel::mclapply(seq_len(nrow(units)), runChild, mc.cores = workers, mc.preschedule = FALSE)
checks <- if (length(checks_env$rows)) do.call(rbind, checks_env$rows) else NULL
timing <- list()
status_rows <- list()
plant_rows <- list()
for (i in seq_len(nrow(units))) {
    res <- results[[i]]
    if (!is.list(res)) {
        u <- units[i, ]
        res <- list(timing = data.frame(cell = u$cell, job = u$job, rep = u$r, scenario = u$scenario, mask = u$mask, wall_seconds = NA,
            unit_seconds = NA, exit = NA, error = oneLine(paste("worker failed:", paste(as.character(res), collapse = " "))), rows = 0L,
            rows_sha256 = NA, fault = u$fault, e7 = u$e7, stringsAsFactors = FALSE))
    }
    timing[[i]] <- res$timing
    if (!is.null(res$checks)) checks <- rbind(checks, res$checks)
    status_rows[[i]] <- res$status
    if (!is.null(res$plant)) plant_rows[[length(plant_rows) + 1L]] <- res$plant
}
unlink(list.files(file.path(out_dir, "rows"), pattern = "^unit-(timing|checks)-", full.names = TRUE))
timing <- do.call(rbind, timing)
# a unit that ended with an error is a failed check in a lifecycle run (the planted exit and guard units end so by design)
for (i in which(nzchar(timing$error))) checks <- rbind(checks, data.frame(check = "unit_completed", cell = timing$cell[[i]], rep = timing$rep[[i]],
    condition = paste(timing$job[[i]], units$plant[[i]]), value = NA, threshold = "ran to the end", pass = units$plant[[i]] %in% c("exitplant", "guardunit", "crashnet"),
    note = timing$error[[i]], stringsAsFactors = FALSE))
statuses <- do.call(rbind, status_rows)
statuses <- stats::aggregate(statuses["n"], statuses[c("cell", "job", "arm", "arm_scale", "scenario", "mask", "status", "class")], sum)
# E7 (crash): the deliberate fault was recorded crashed and the job went on
if (lifecycle && fault) {
    crashed <- statuses[statuses$status == "crashed" & grepl("deliberate fault", statuses$class, fixed = TRUE), , drop = FALSE]
    checks <- rbind(checks, data.frame(check = "E7_crash", cell = units$cell[[1L]], rep = units$r[[1L]], condition = "first unit",
        value = as.character(sum(crashed$n)), threshold = ">= 1 rows", pass = nrow(crashed) > 0L && nrow(timing) == nrow(units) &&
            all(nzchar(timing$error) == FALSE), note = sprintf("units completed %d of %d", sum(!nzchar(timing$error)), nrow(units)),
        stringsAsFactors = FALSE))
}
# E7 (no_da): the deliberate per-assay DA failure was scored `no_da` for its family and the journey went on
if (lifecycle && fault && identical(part, "L")) {
    noda <- statuses[statuses$status == "no_da" & grepl("deliberate per-assay", statuses$class, fixed = TRUE), , drop = FALSE]
    checks <- rbind(checks, data.frame(check = "E7_noda", cell = units$cell[[1L]], rep = units$r[[1L]], condition = "first unit",
        value = as.character(sum(noda$n)), threshold = ">= 1 journeys", pass = nrow(noda) > 0L, note = "families of the failed assay scored no_da",
        stringsAsFactors = FALSE))
}
# E5 (reported jobs' LA and OR rows equal the deciding jobs', seconds removed): compared by their unblinded hashes per row
if (lifecycle && all(c("deciding", "reported") %in% jobs) && identical(part, "L")) {
    read_unit <- function(cellname, job, r) utils::read.delim(file.path(out_dir, "rows", sprintf("part-%s-%s-%d-all-all.tsv", cellname, job, r)),
        stringsAsFactors = FALSE, colClasses = "character")
    for (cellname in cells_wanted) for (r in seq(first_rep, last_rep)) {
        a <- tryCatch(read_unit(cellname, "deciding", r), error = function(e) NULL)
        b <- tryCatch(read_unit(cellname, "reported", r), error = function(e) NULL)
        if (is.null(a) || is.null(b)) next
        pick <- function(x) {
            x <- x[x$arm %in% c("LA", "OR") & x$scenario %in% SCENARIOS[1:4] & x$scoreset == "T",
                setdiff(names(x), c("job", "seconds", "warnings", "n_common", "s5", "gross1_common", "or_gross1_common", "mamp_sigma_min", "mamp_e_max",
                    "mamp_planted_share", "mdf_null05", "mdf_null05_reduced", "mdf_df_le_k", "mdf_tested")), drop = FALSE]
            x[do.call(order, x[c("arm", "scenario", "mask", "assay", "contrast")]), , drop = FALSE]
        }
        pa <- pick(a)
        pb <- pick(b)
        rownames(pa) <- rownames(pb) <- NULL
        checks <- rbind(checks, data.frame(check = "E5r", cell = cellname, rep = r, condition = "LA and OR rows, deciding against reported",
            value = as.character(nrow(pa)), threshold = "identical (written columns)", pass = identical(pa, pb), note = "lifecycle: blinded columns not compared",
            stringsAsFactors = FALSE))
    }
}
# SA-P12 (7.4 E5, 34): LP's and LAs' oracle rows equal the deciding job's (same commit; `seconds` removed)
if (lifecycle && identical(part, "L") && "deciding" %in% jobs && any(c("lp", "las") %in% jobs)) {
    read_job <- function(cellname, job, r) tryCatch(utils::read.delim(file.path(out_dir, "rows", sprintf("part-%s-%s-%d-all-all.tsv", cellname, job, r)),
        stringsAsFactors = FALSE, colClasses = "character"), error = function(e) NULL)
    for (job2 in intersect(c("lp", "las"), jobs)) for (cellname in cells_wanted) for (r in seq(first_rep, last_rep)) {
        a <- read_job(cellname, "deciding", r)
        b <- read_job(cellname, job2, r)
        if (is.null(a) || is.null(b)) next
        pick <- function(x) {
            x <- x[x$arm %in% "OR" & x$scenario %in% SCENARIOS[1:4] & x$scoreset == "T",
                setdiff(names(x), c("job", "seconds", "warnings", "n_common", "s5", "gross1_common", "or_gross1_common", "mamp_sigma_min", "mamp_e_max",
                    "mamp_planted_share", "mdf_null05", "mdf_null05_reduced", "mdf_df_le_k", "mdf_tested")), drop = FALSE]
            x[do.call(order, x[c("scenario", "mask", "assay", "contrast")]), , drop = FALSE]
        }
        pa <- pick(a)
        pb <- pick(b)
        rownames(pa) <- rownames(pb) <- NULL
        checks <- rbind(checks, data.frame(check = "E5r", cell = cellname, rep = r, condition = paste("OR rows, deciding against", job2),
            value = as.character(nrow(pa)), threshold = "identical (written columns)", pass = identical(pa, pb), note = "lifecycle: blinded columns not compared",
            stringsAsFactors = FALSE))
    }
}
# SA-P16 (35.6): L1r's silence: no instrument row, no itsd_ok or norm_ok FALSE, E14's rows all passing, over every non-l1c unit
silent <- do.call(rbind, lapply(results, function(x) if (is.list(x)) x$silence))
if (lifecycle && !is.null(silent) && !("l1c" %in% jobs)) {
    e14_failed <- if (is.null(checks)) 0L else sum(checks$check == "E14" & !checks$pass)
    checks <- rbind(checks, data.frame(check = "E15_silence", cell = "all", rep = 0L, condition = "every journey of the run",
        value = as.character(sum(silent$n_instrument, silent$n_itsd_false, silent$n_norm_false, e14_failed)), threshold = "0 instrument rows, 0 itsd_ok or norm_ok FALSE, 0 failing E14 rows",
        pass = !sum(silent$n_instrument, silent$n_itsd_false, silent$n_norm_false, e14_failed), note = sprintf("%d units", nrow(silent)),
        stringsAsFactors = FALSE))
}
# assemble one table per cell, the parts removed after
for (cellname in cells_wanted) {
    files <- list.files(file.path(out_dir, "rows"), pattern = sprintf("^part-%s-.*\\.tsv$", cellname), full.names = TRUE)
    if (!length(files)) next
    table <- do.call(rbind, lapply(files, utils::read.delim, stringsAsFactors = FALSE, na.strings = "NA", colClasses = "character"))
    out_file <- file.path(out_dir, "rows", sprintf("%s-%s-%d-%d%s.tsv.gz", tolower(part), cellname, first_rep, last_rep, if (lifecycle) "-l1" else ""))
    connection <- gzfile(out_file, "w")
    utils::write.table(table, connection, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
    close(connection)
    unlink(files)
}
utils::write.table(timing, file.path(out_dir, sprintf("timing-%s.tsv", job_name)), sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
utils::write.table(statuses, file.path(out_dir, sprintf("status-%s.tsv", job_name)), sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
if (length(plant_rows)) utils::write.table(do.call(rbind, plant_rows), file.path(out_dir, sprintf("plants-%s.tsv", job_name)), sep = "\t",
    quote = FALSE, row.names = FALSE, na = "NA")
if (!is.null(checks) && nrow(checks)) utils::write.table(checks, file.path(out_dir, sprintf("checks-%s.tsv", job_name)), sep = "\t",
    quote = FALSE, row.names = FALSE, na = "NA")
writeProvenance(TRUE)
cat(sprintf("KIT-188 %s finished: %d units (%d with an error), wall %.1f min\n", job_name, nrow(units), sum(nzchar(timing$error)),
    as.numeric(difftime(Sys.time(), started_at, units = "mins"))))
# 35.6 (r1.20): a prediction not met (E15, E15R, or the R-marked E7 observations) stops the lifecycle run; everything is written first
# K188_SMOKE=1 (an instrument check with no expectation, declared before it is sent): the checks are written, none is judged
if (!is.null(checks) && !identical(Sys.getenv("K188_SMOKE"), "1") &&
    any(!checks$pass[checks$check %in% c("E15", "E15R", "E7_invalid_k", "E7_refusal", "E15_silence", "GUARD_planted", "unit_completed")])) {
    cat("KIT-188: a prediction of 35.6 or 7.4 was not met; the lifecycle run stops here (reported, 35.6 r1.20)\n")
    quit(save = "no", status = 2L)
}
