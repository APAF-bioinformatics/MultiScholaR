# KIT-275 (VAL-330, VAL-344, VAL-347): the sample-correlation filter after RUV and one anomalous sample (lipidomics, with a
# metabolomics identity check). A measurement; evidence code under docs/qa; no package code. Runs on a GitHub-hosted runner with
# MultiScholaR installed from the bound commit (run-tests `scripts`), from the package root. The generator and the journeys are a
# copy of KIT-188's sim-k188.R at 0866982f with the marked changes of generator-diff.md (marks K275-*); nothing else of it is carried.
#
# Usage:
#   Rscript sim-k275.R <tiers> <design> <first> <last> <out dir>       the job (parent)
#   Rscript sim-k275.R --unit <tier> <design> <r> <rds>                one data set of one tier (child process)
#   tiers   one tier or a comma list of D, K, F, S, V, M, T, E, ID (tier V waits for KIT-242, section 13 of the ticket)
#   design  D3, D4 or D6 (ID: D4); first, last  data sets (counted tiers 1 to 200 or 1 to 100; ID 1 to 20; L1 and the smoke 9101 to 9103)
# Environment:
#   K275_WORKERS   child processes at a time (default 4)
#   K275_LIFECYCLE 1: blind lifecycle (L1): every estimand column of the rows is written NA, the rows' SHA-256 is written, the
#                  checks E1 to E12 are run and written
#   K275_REPEAT    lifecycle: the tiers whose first data set is run a second time for E9 (default all; the smoke send sets one small tier)
#   K275_FAULT     1: E7: one deliberate error in the first journey of the first unit, and the refused journey of E7 (GK3 at D3 on VD)
#   K275_SMOKE     1: the smoke send's flag; written to the provenance and to every row; the scorer refuses any input carrying it for an
#                  outcome line; a counted request that sets it (without K275_LIFECYCLE=1) is refused here (SMOKE_counted_refusal)
#   K275_JOB       a name for the job's files (default from the tiers, design and range)
# One data set per child process, all its versions and arms; no session owner (duckdb) is forked. Every random draw of a data set is
# made in generate() before any package function is called; seeds in section 5.3 of the ticket.

# Partial-match guard (lesson of KIT-188's smoke checks 1 to 3): `$` on a list or data frame matches a name by its prefix when none is
# exact; that warning is an error here, in every process and inside `quiet()` (which would otherwise muffle it). A stop from it is an
# instrument defect by definition (VAL-54). Every read of a field that may be absent below uses `[[` with exact matching.
options(warnPartialMatchDollar = TRUE) # K275-guard
# A `stop()` raised from a GLOBAL calling handler runs with every local handler disabled (R's rule for the handler stack), so it cannot reach a
# `tryCatch`: it ends the process with nothing written (smoke 2's planted guard unit). The guard is therefore also installed as a LOCAL calling
# handler (`guardWarning`) around each journey and around each unit, where its stop reaches the enclosing `tryCatch`; the global handler stays as
# a last net, with an error writer that leaves a trace before the process ends.
guardWarning <- function(w) { # K275-guard
    if (grepl("partial match of", tolower(conditionMessage(w)), fixed = TRUE)) {
        options(k275.guard_trace = paste(vapply(sys.calls(), function(cl) substr(paste(deparse(cl)[[1L]], collapse = ""), 1L, 80L), ""), collapse = " <- "))
        stop(paste("partial-match guard:", conditionMessage(w)), call. = FALSE)
    }
}
globalCallingHandlers(warning = guardWarning) # K275-guard
if (identical(commandArgs(TRUE)[1L], "--unit")) options(error = function() { # K275-guard   # an error that escaped every handler: write what is known, then end non-zero
    dir <- dirname(Sys.getenv("K275_CHECKS_PATH"))
    tag <- sub("\\.tsv$", "", sub("^checks-", "", basename(Sys.getenv("K275_CHECKS_PATH"))))
    try(writeLines(c(geterrmessage(), getOption("k275.guard_trace", "")), file.path(dir, sprintf("unit-trace-%s.txt", tag))), silent = TRUE)
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
cellDef <- function(idx, part, omic, assays, per_assay, groups, per, level, sd, unwanted, blocks = 0L, n2 = FALSE, no_is = FALSE) { # K275-cells
    list(idx = idx, part = part, omic = omic, assays = assays, per_assay = per_assay, groups = groups, per = per, level = level,
        sd = sd, unwanted = unwanted, blocks = blocks, n2 = n2, no_is = no_is) # K275-cells
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
    K172c = cellDef(0L, "L", "lipid", "LCMS_Pos", 515L, 2L, 6L, NA, NA, "none"),
    # K275-cells: the designs of KIT-275 (2 groups of 3, 4, 6; two assays of 250 features; one unwanted component; no internal standards), the
    # metabolomics identity cell D4M (the same seed index as D4, so the draws are the same), and the exactness fixtures of G1 (one assay, 2,000 features)
    D3 = cellDef(1L, "L", "lipid", c("LCMS_Neg", "LCMS_Pos"), 250L, 2L, 3L, 20, 2, "drift_one", no_is = TRUE),
    D4 = cellDef(2L, "L", "lipid", c("LCMS_Neg", "LCMS_Pos"), 250L, 2L, 4L, 20, 2, "drift_one", no_is = TRUE),
    D6 = cellDef(3L, "L", "lipid", c("LCMS_Neg", "LCMS_Pos"), 250L, 2L, 6L, 20, 2, "drift_one", no_is = TRUE),
    D4M = cellDef(2L, "L", "metab", c("LCMS_Neg", "LCMS_Pos"), 250L, 2L, 4L, 20, 2, "drift_one", no_is = TRUE),
    X3 = cellDef(5L, "L", "lipid", "LCMS_Pos", 2000L, 2L, 3L, 20, 2, "drift_one", no_is = TRUE),
    X4 = cellDef(6L, "L", "lipid", "LCMS_Pos", 2000L, 2L, 4L, 20, 2, "drift_one", no_is = TRUE),
    X6 = cellDef(7L, "L", "lipid", "LCMS_Pos", 2000L, 2L, 6L, 20, 2, "drift_one", no_is = TRUE))
scenariosOf <- function(cellname) if (isTRUE(CELLS[[cellname]]$n2)) SCENARIOS else SCENARIOS[1:4]
MASKS <- c("none", "vd5", "vd20")
REFUSALS <- c("multischolar_refusal", "multischolar_ruv_grouping_unusable", "omics_ruv_insufficient_controls", "prot_ruv_insufficient_controls", "prot_ruv_empty_result", # K275-REFUSALS
    "workflow_execution_invalid", "ruv_factor_count_unavailable", "prot_ruv_factors_unavailable", "omics_ruv_factors_refused",
    "omics_ruv_no_residual_df", "ruv_k_refused", "ruv_k_above_limit")

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
        guardWarning(w) # K275-guard
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
    # K275-markers: the control selection before KIT-190 (positional controls; group-blind), for LG1 and LG2
    LG1 = list(MultiScholaR.ruv_control_grouping = "dropped"),
    LG2 = list(MultiScholaR.ruv_control_grouping = "dropped"))
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
generate <- function(cellname, scenario, r, anomaly = NULL) { # K275-sig
    p <- CELLS[[cellname]]
    set.seed(275000000 + 100000 * p$idx + 10000 * match(scenario, SCENARIOS) + r, kind = "Mersenne-Twister", # K275-seed
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
    n_is <- if (identical(omic, "lipid") && !isTRUE(p$no_is)) 5L * n_assays else 0L # K275-noIS
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
    C <- switch(p$unwanted, drift_one = 1L, drift_factor = 2L, plate_drift = 2L, batch_drift = p$blocks) # K275-C
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
    # K275-plant (section 5.3 of the ticket): the planting of one anomalous sample (two for V2) in assay LCMS_Pos, after every draw of
    # generate(); its own stream, the session's random-number state restored, so no draw of any version differs from V0 except the planted entries
    planted_delta <- NULL
    if (!is.null(anomaly)) {   # anomaly = list(form = "g" | "p" | "s" | "t", samples = c(n_runs) or c(n_runs, per))
        old <- if (exists(".Random.seed", envir = globalenv())) get(".Random.seed", envir = globalenv()) else NULL
        set.seed(276000000 + 100000 * p$idx + r, kind = "Mersenne-Twister", normal.kind = "Inversion", sample.kind = "Rejection")
        z_feat <- list(stats::rnorm(n_feat), stats::rnorm(n_feat)); z_lbl <- list(stats::rnorm(n_lbl), stats::rnorm(n_lbl))  # drawn in this order, all, always
        u_feat <- stats::runif(n_feat)   # drawn last, all, always: the sparse form's choice of features
        cm <- if (identical(anomaly$form, "p")) 0.5 else 0
        planted_delta <- list()
        for (s in seq_along(anomaly$samples)) {
            col <- anomaly$samples[[s]]
            rf <- which(feat_assay == "LCMS_Pos"); rl <- which(lbl_assay == "LCMS_Pos")
            d_feat <- (p$level + cm * (mu[rf] - p$level) + sqrt(1 - cm^2) * p$sd * z_feat[[s]][rf]) - mu[rf]
            if (identical(anomaly$form, "t")) d_feat[!(seq_along(rf) %in% order(u_feat[rf])[seq_len(round(0.1 * length(rf)))])] <- 0
            y_feat[rf, col] <- y_feat[rf, col] + d_feat
            d_lbl <- if (!anomaly$form %in% c("s", "t")) (p$level + cm * (mu_lbl[rl] - p$level) + sqrt(1 - cm^2) * p$sd * z_lbl[[s]][rl]) - mu_lbl[rl] else numeric(length(rl))
            y_lbl[rl, col] <- y_lbl[rl, col] + d_lbl
            planted_delta[[s]] <- c(stats::setNames(d_feat, feat_ids[rf]), stats::setNames(d_lbl, lbl_ids[rl]))
        }
        if (!is.null(old)) assign(".Random.seed", old, envir = globalenv()) else rm(".Random.seed", envir = globalenv())
    }
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
        U = U, order_of = order_of, ba_u = ba_u, ba_block = ba_block, Z = Z, C = C, planted_delta = planted_delta) # K275-ret
}

# the mask (TRUE = missing) on the feature rows; standards and labelled controls are never missing
solveLod <- function(values, rate) {
    stats::uniroot(function(lod) mean(stats::plogis(lod - values)) - rate, c(min(values) - 60, max(values) + 60), tol = 1e-10)$root
}
maskOf <- function(sim, mask, lod_values = sim$y[sim$feat_ids, , drop = FALSE]) { # K275-mask: the limit of detection is solved on `lod_values` (V0's values)
    y <- sim$y[sim$feat_ids, , drop = FALSE]
    if (identical(mask, "none")) return(matrix(FALSE, nrow(y), ncol(y), dimnames = dimnames(y)))
    rate <- if (identical(mask, "ba")) 0.20 else as.numeric(sub("^vd", "", mask)) / 100
    out <- sim$U < stats::plogis(solveLod(lod_values, rate) - y) # K275-mask
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
# K275-E3: the ticket's E3 is the route defaults of 5.2 only (ITSD is off in every journey here, so KIT-188's clauses on the ITSD
# import pattern and on labelled names matching an ITSD pattern are not carried: at the bound commit the import pattern is still the one
# before KIT-260, "w/o MS2:"); the replay markers hold their current token (a leak stops the unit).
checkE3 <- function(omic) {
    got <- routeDefaults(omic)
    want <- EXPECTED_DEFAULTS[[omic]]
    differ <- names(want)[!vapply(names(want), function(n) identical(got[[n]], want[[n]]) ||
        (is.numeric(want[[n]]) && is.numeric(got[[n]]) && isTRUE(all.equal(got[[n]], want[[n]]))), logical(1L))]
    leaked <- markerLeaks()
    got$cyclic_loess_minimum <- tryCatch(omicsCyclicLoessMinimumActive(), error = function(e) NA) # recorded, either state passes (as KIT-188)
    need <- c("ruv_control_grouping", "ruv_k_bound", "ruv_feature_levels", "correlation_exclusion", strsplit(Sys.getenv("K275_REQUIRED_MARKERS"), ",", fixed = TRUE)[[1L]]) # K275-E3
    absent <- setdiff(need[nzchar(need)], names(workflowReplayRuleMarkers)) # K275-E3
    ce_ok <- identical(workflowReplayRuleMarkers[["correlation_exclusion"]][["current"]], "within_assay") # K275-E3
    got$correlation_exclusion <- workflowReplayRuleMarkers[["correlation_exclusion"]][["current"]] # K275-E3
    list(pass = !length(differ) && !length(leaked) && !length(absent) && ce_ok, got = got, # K275-E3
        note = sprintf("differ: %s; markers off their current token: %s; markers absent: %s; correlation_exclusion within_assay: %s", paste(differ, collapse = ","),
            paste(leaked, collapse = ","), paste(absent, collapse = ","), ce_ok))
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
    identical(unname(a[rownames(a), colnames(a), drop = FALSE]), unname(b[rownames(a), colnames(a), drop = FALSE]))
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
        e12 <- if (!is.null(out$itsd_check)) out$itsd_check[out$itsd_check$assay == a, , drop = FALSE] # K275-prefix
        e13 <- if (!is.null(out$norm_check)) out$norm_check[out$norm_check$assay == a, , drop = FALSE] # K275-prefix
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
            if (any(res$status %in% c("skipped_k_exceeds_complete_controls", "skipped_no_controls_selected"))) reason <- paste0(reason, " | ", oneLine(trace$error))
        }
        record <- out$record_ruv
        legacy_return <- requested && !declared && status %in% c("corrected", "uncorrected", "dropped")
        da_failed <- identical(out$error_step, "da")
        failed_here <- if (!is.null(out$da_errors)) out$da_errors[out$da_errors$assay == a, , drop = FALSE]
        da_status <- if (da_failed) "failed" else if (!isTRUE(out$reached_da)) "not_reached" else # K275-prefix
            if (is.null(out$da) || !nrow(out$da) || !a %in% out$da$assay || (!is.null(failed_here) && nrow(failed_here))) "no_da" else "ok"
        list(itsd_is_left = field(e12, "itsd_is_left"), itsd_n_used = field(e12, "itsd_n_used"), itsd_n_expected = field(e12, "itsd_n_expected"),
            itsd_max_change = field(e12, "itsd_max_change"), itsd_ok = field(e12, "itsd_ok"),
            norm_max_change = field(e13, "norm_max_change"), norm_ok = field(e13, "norm_ok"),
            ruv_status = status, ruv_status_product = product_status,
            ruv_status_agree = if (legacy_return && !is.na(product_status)) identical(status, product_status) else NA,
            ruv_changed = if (legacy_return) identical(status, "corrected") else NA,
            ruv_k_found = if (requested && !declared && !is.null(trace)) (if (isTRUE(trace$success)) as.integer(trace$best_k) else NA_integer_) else NA_integer_,
            ruv_n_controls = if (requested && !declared && !is.null(trace) && is.logical(trace[["control_genes_index"]])) sum(trace[["control_genes_index"]], na.rm = TRUE) else NA_integer_, # K275-assaycols
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
    result <- tryCatch(withCallingHandlers({ # K275-guard
        object <- omicObject(sim, log2_values, ids = ids, scale = scale, plant = omic_plant)
        fixture <- quiet(if (identical(sim$p$omic, "lipid")) lipidDeclarationsOwner(root, object) else metabDeclarationsOwner(root, object))
        fixture$workflow$contrasts_tbl <- contrastsOf(sim)
        owner <- fixture$owner
        config <- fixture$workflow$config_list
        if (fault) stop("K188 deliberate fault (E7)")
        if (isTRUE(getOption("k275.guard_now"))) { # K275-guard: the planted partial `$` match on a data frame (the self-check of the guard)
            options(k275.guard_now = FALSE)
            out$guard_planted <- TRUE
            invisible(data.frame(abcdef = 1)$abc)
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
            out$itsd_check <- checkItsdActed(m_qc, out$m_itsd, itsd_record) # K275-prefix
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
        step <- "pairs_pre" # K275-pairs
        out$pairs_pre <- quiet(pearsonCorForSamplePairs(theObject = owner$getState(paste0(prefix, "_normalized")), correlation_group = "group")) # K275-pairs
        step <- "between_sample"
        out$rows_ruv_input <- vapply(out$m_in, nrow, integer(1L))
        if (!identical(method, "none") && !identical(arm, "OR")) {
            compared <- out$m_in
            # `normplant` (35.6): the matrix E13 compares is the log2 matrix for one assay; the route's output is untouched
            if (identical(plant, "normplant")) compared[["LCMS_Pos"]] <- out$m_log2[["LCMS_Pos"]]
            out$norm_check <- checkNormActed(compared, out$m_log2) # K275-prefix
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
                "control_count", "ruv_k_estimate", "W"))]) # K275-W
            out$optimisation <- lapply(run$optimisation, function(o) o[intersect(names(o), c("success", "best_k", "best_percentage",
                "control_genes_index", "optimization_results", "error"))])
            out$m_out <- stateMatrices(owner, binding, paste0(prefix, "_ruv_corrected"))
            step <- "pairs_post" # K275-pairs
            out$pairs_post <- quiet(pearsonCorForSamplePairs(theObject = owner$getState(paste0(prefix, "_ruv_corrected")), correlation_group = "group")) # K275-pairs
            step <- "ruv"
            # 35.2 (3): the analysis record's fields beside the product's (they must equal it)
            record <- tryCatch(quiet(buildWorkflowAnalysisRecord(owner)), error = function(e) if (grepl("partial-match guard", conditionMessage(e), fixed = TRUE)) stop(e) else NULL) # K275-guard
            out$record_ruv <- recordOutcomes(record, paste0(prefix, "_ruv_corrected"))
        }
        if (correlation) {
            step <- "correlation"
            quiet(runWorkflowOmicsReferenceCorrelation(binding, owner, "group", d$correlation, config))
            out$corr_one_assay <- paste(oneAssayRemovals(out$m_out %||% out$m_in, stateMatrices(owner, binding, owner$current_state)),
                collapse = ";")
            # K275-pairs: per assay, the samples that held a value before the filter and none after it (the package filter's own removed set, E5)
            after_filter <- stateMatrices(owner, binding, owner$current_state)
            out$corr_removed <- lapply(stats::setNames(names(after_filter), names(after_filter)), function(a) {
                before <- out$m_out[[a]] %||% out$m_in[[a]] # K275-pairs
                had <- colnames(before)[colSums(is.finite(before)) > 0L]
                kept <- colnames(after_filter[[a]])[colSums(is.finite(after_filter[[a]])) > 0L]
                setdiff(had, kept)
            }) # K275-pairs
        }
        final <- owner$getState(owner$current_state)
        out$runs_final <- as.character(final@design_matrix[[final@sample_id]])
        left <- unique(as.character(final@design_matrix[[final@group_id]]))
        if (!all(sim$groups %in% left)) { # K275-nogroup: the filter removed a whole group from the design: the DA has nothing to test (a result, `no_da`)
            out$class <- paste0("filter_removed_group:", paste(setdiff(sim$groups, left), collapse = ","))
        } else {
        step <- "da"
        out$reached_da <- TRUE # K275-prefix
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
        } # K275-nogroup
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
    }, warning = guardWarning), error = function(e) { # K275-guard
        cls <- intersect(class(e), REFUSALS)
        guard <- grepl("partial-match guard", conditionMessage(e), fixed = TRUE) # K275-guard
        text <- if (length(cls)) cls[[1L]] else oneLine(paste(class(e)[[1L]], conditionMessage(e), sep = ": "))
        # r1.17 (R16-P07): a journey already `instrument` keeps its status; the later error goes to `class`
        if (guard) { # K275-guard: a guard stop is an instrument defect, never a crash
            out$status <<- "instrument"
            out$class <<- oneLine(paste("partial_match_guard:", conditionMessage(e)))
        } else if (identical(out$status, "instrument")) {
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
        k_used <- suppressWarnings(as.integer(unlist(journey[["ruv_value"]][["best_k"]][[a]])[1L])) # K275-prefix
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

# ---- estimands (5.5), one family = (assay, contrast) ---------------------------------------------------------------
groupOfComparison <- function(x) {
    x <- as.character(x)
    full <- regmatches(x, regexec("group([A-E])-groupA", x))
    short <- regmatches(x, regexec("^([A-E])(_vs_A|\\.vs\\.A)", x))
    vapply(seq_along(x), function(i) if (length(full[[i]]) >= 2L) full[[i]][[2L]] else if (length(short[[i]]) >= 2L) short[[i]][[2L]] else
        NA_character_, "")
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
}
#' E9 (Part L, today's base): the legacy controls per assay and percentage are the first round(N_v x pct / 100) rows of
#' the RUV input table among the rows with at least two values, in table order (2.2).
positional <- function(m, pct) {
    valid <- rownames(m)[rowSums(is.finite(m)) > 1L]
    valid[seq_len(round(length(valid) * pct / 100, 0))]
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

# ====================================================================================================================
# K275: everything below is added (not copied): versions, requests, the reads of the filter's statistic, the estimand columns of a
# row, the tiers and the job mechanics. Sections 4 to 8 of the ticket.
# ====================================================================================================================

Y_ASSAY <- "LCMS_Pos"
X_ASSAY <- "LCMS_Neg"
THETAS <- c(0.5, 0.7, 0.9, 0.95)
TAGS <- c("05", "07", "09", "095")
VERSIONS <- c("V0", "V1g", "V1p", "V1s", "V1t", "V2", "VD")

# ---- tier plans (section 4): version, arm; LX is added in tier D from GK2 (V1g) and GK1 (VD) ------------------------------
planOf <- function(version, arms) data.frame(version = version, arm = arms, stringsAsFactors = FALSE)
TIER_PLAN <- list(
    D = rbind(planOf("V0", c("NO", "LM2", "GK2", "LA", "GA")), planOf("V1g", c("NO", "LM2", "GK2", "LA", "GA")),
        planOf("VD", c("LM1", "GK1", "LA", "GA"))),
    K = rbind(planOf("V0", c("LM1", "LM3", "GK1", "GK3", "LM2@20", "LG2")), planOf("V1g", c("LM1", "LM3", "GK1", "GK3", "LM2@20", "LG2")),
        planOf("VD", c("NO", "LM2", "GK2", "LM1@20", "LG1"))),
    F = planOf("V1p", c("NO", "LM2", "GK2", "LA", "GA")),
    S = rbind(planOf("V0", "GSC"), planOf("V1g", "GSC"), planOf("V1s", c("GSC", "GK2"))),
    V = rbind(planOf("V0", c("NO", "LM2", "GK2", "LA", "GA")), planOf("V1g", c("NO", "LM2", "GK2", "LA", "GA")),
        planOf("VD", c("LM1", "GK1"))),
    M = planOf("V2", c("NO", "LM1", "LM2", "LM3", "GK1", "GK2", "GK3")),
    T = planOf("V1t", c("NO", "LM2", "LG2", "GK2", "LA", "GA")),
    E = NULL,   # the emulation arms are built in runEmulation (two passes)
    ID = planOf("V1g", c("NO", "LM2", "GK2", "LA", "GA")))
TIERS <- c("D", "K", "F", "S", "V", "M", "T", "E", "ID")
tierDesigns <- function(tier) if (identical(tier, "M")) "D6" else if (identical(tier, "ID")) "D4" else c("D3", "D4", "D6")

# ---- arms (section 4) ------------------------------------------------------------------------------------------------
armSpec <- function(arm) {
    switch(arm,
        NO = list(kind = "NO"), LA = list(kind = "LA"), GA = list(kind = "GA"),
        GSC = list(kind = "GSC", k = 2L),
        PEL2 = list(kind = "LM", k = 2L), PEG2 = list(kind = "GK", k = 2L),
        {
            m <- regmatches(arm, regexec("^(LM|LG|GK)([0-9])(@([0-9]+))?$", arm))[[1L]]
            if (length(m) != 5L) stop(sprintf("K275: unknown arm %s", arm))
            list(kind = m[[2L]], k = as.integer(m[[3L]]), pct = if (nzchar(m[[5L]])) as.numeric(m[[5L]]) else NA_real_)
        })
}
armClass <- function(arm) if (arm %in% c("PEL2", "PEG2")) arm else armSpec(arm)[["kind"]]
isDeclared <- function(arm) armSpec(arm)[["kind"]] %in% c("GA", "GK", "GSC")
requestK275 <- function(arm, d, sim) {
    base <- REQUESTS_L(d)
    spec <- armSpec(arm)
    switch(spec[["kind"]],
        NO = base[["NO"]], LA = base[["LA"]], GA = base[["GA"]],
        LM = , LG = {
            rq <- base[["LM"]]
            rq[["ruv"]][["manual_k"]] <- spec[["k"]]
            rq[["ruv"]][["manual_percentage"]] <- if (is.na(spec[["pct"]] %||% NA_real_)) d$manual_pct else spec[["pct"]]
            rq
        },
        GK = {
            rq <- base[["GA"]]
            rq[["ruv"]][["ruv_k"]] <- spec[["k"]]
            rq[["ruv"]][["ruv_k_rule"]] <- "user"
            rq
        },
        GSC = list(gate = TRUE, ruv = list(ruv_design = "supplied_controls", ruv_grouping_variable = "group", ruv_k = spec[["k"]],
            ruv_k_rule = "user", ruv_k_seed = RUV_FACTOR_COUNT_SEED, ruv_control_ids = sim$lbl_ids)))
}

# ---- versions of a data set (section 4; 5.3) ---------------------------------------------------------------------------
maskedLog2 <- function(sim, mask, lod_values) {
    values <- sim$y
    if (!identical(mask, "none")) {
        block <- values[sim$feat_ids, , drop = FALSE]
        block[maskOf(sim, mask, lod_values)] <- NA_real_
        values[sim$feat_ids, ] <- block
    }
    values
}
#' The version with the named runs removed from the design and from both assays (VD, VE).
dropSamples <- function(sim, drop_runs) {
    keep <- !sim$runs %in% drop_runs
    sim$runs <- sim$runs[keep]
    sim$group <- sim$group[keep]
    sim$block <- sim$block[keep]
    sim$design <- sim$design[keep, , drop = FALSE]
    for (field in c("y", "oracle", "unwanted", "U")) sim[[field]] <- sim[[field]][, keep, drop = FALSE]
    sim$l <- sim$l[keep]
    if (!is.null(sim$Z)) sim$Z <- sim$Z[keep, , drop = FALSE]
    sim
}
makeVersion <- function(cellname, r, version, mask = "none", v0 = NULL) {
    per <- CELLS[[cellname]]$per
    n_runs <- CELLS[[cellname]]$groups * per
    a_run <- sprintf("R%03d", n_runs)
    if (identical(version, "VD")) {
        sim <- dropSamples(v0[["sim"]], a_run)
        return(list(sim = sim, values = v0[["values"]][, sim$runs, drop = FALSE]))
    }
    anomaly <- switch(version, V0 = NULL, V1g = list(form = "g", samples = n_runs), V1p = list(form = "p", samples = n_runs),
        V1s = list(form = "s", samples = n_runs), V1t = list(form = "t", samples = n_runs),
        V2 = list(form = "g", samples = c(n_runs, per)), stop(sprintf("K275: unknown version %s", version)))
    sim <- generate(cellname, "N1", r, anomaly)
    reference <- if (is.null(v0)) sim else v0[["sim"]]
    list(sim = sim, values = maskedLog2(sim, mask, reference$y[reference$feat_ids, , drop = FALSE]))
}
#' The anomaly's share in the contrast B - A (x_j = planted_delta_j / n, + for group B, - for group A), over Y's features, for each form.
plantedShares <- function(cellname, r) {
    per <- CELLS[[cellname]]$per
    n_runs <- CELLS[[cellname]]$groups * per
    one <- function(form, samples) {
        deltas <- generate(cellname, "N1", r, list(form = form, samples = samples))$planted_delta
        keep <- grepl("^F", names(deltas[[1L]]))
        x <- deltas[[1L]][keep] / per
        if (length(deltas) > 1L) x <- (deltas[[1L]][keep] - deltas[[2L]][keep]) / per
        list(x = x, delta = if (length(deltas) > 1L) sqrt(deltas[[1L]]^2 + deltas[[2L]]^2) else deltas[[1L]])
    }
    list(g = one("g", n_runs), p = one("p", n_runs), s = one("s", n_runs), t = one("t", n_runs), two = one("g", c(n_runs, per)))
}

# ---- the filter's statistic (2.1, 6.1 A): the best correlation of a sample with any other sample of its group ------------
bestCorrelation <- function(pairs, samples, id = "Run") {
    out <- stats::setNames(rep(NA_real_, length(samples)), samples)
    if (is.null(pairs) || !nrow(pairs)) return(out)
    x <- as.character(pairs[[paste0(id, ".x")]])
    y <- as.character(pairs[[paste0(id, ".y")]])
    value <- as.numeric(pairs[["pearson_correlation"]])
    for (s in samples) {
        v <- value[(x == s | y == s) & is.finite(value)]
        if (length(v)) out[[s]] <- max(v)
    }
    out
}
failsAt <- function(c_s, theta) is.finite(c_s) & c_s < theta
#' The samples the filter fails at `theta` (a sample with no finite correlation is not tested, KIT-255).
failList <- function(c_s, theta) names(c_s)[failsAt(c_s, theta)]

#' A6: the pseudo-replicate ratio PR = d_a / median of d_s over the other samples, d_s the RMS over the features of a sample's value
#' minus the mean of its group-mates' values.
pseudoReplicateRatio <- function(M, design, a_run) {
    runs <- colnames(M)
    if (!a_run %in% runs) return(NA_real_)
    group <- design$group[match(runs, design$Run)]
    d <- vapply(seq_along(runs), function(i) {
        mates <- runs[group == group[[i]] & runs != runs[[i]]]
        if (!length(mates)) return(NA_real_)
        sqrt(mean((M[, runs[[i]]] - rowMeans(M[, mates, drop = FALSE], na.rm = TRUE))^2, na.rm = TRUE))
    }, numeric(1L))
    names(d) <- runs
    value <- d[[a_run]] / stats::median(d[names(d) != a_run], na.rm = TRUE)
    if (is.finite(value)) value else NA_real_
}
#' A5: the share of the anomaly's energy in the controls the journey used (mean planted_delta^2 over the controls that are affected rows,
#' over the mean over all affected rows).
kappaControls <- function(delta, controls) {
    if (is.null(delta) || !length(controls)) return(NA_real_)
    den <- mean(delta^2)
    use <- intersect(names(delta), controls)
    if (!length(use) || !is.finite(den) || den == 0) return(NA_real_)
    mean(delta[use]^2) / den
}

#' B2: the unmoderated variance of the contrast B - A per feature by ordinary least squares on the matrix and model the DA read
#' (`~ 0 + group`, with the RUV factors W as covariates on the declared route), per feature on its observed values:
#' v_j = RSS_j / df_j * c'(X'X)^-1 c. Returns the sum over the features, their number and the median residual degrees of freedom.
unmoderatedReading <- function(M, design, features, W = NULL) {
    out <- list(sum = NA_real_, n = 0L, df = NA_real_)
    if (is.null(M)) return(out)
    samples <- intersect(colnames(M), design$Run)
    group <- factor(design$group[match(samples, design$Run)], levels = c("A", "B"))
    X <- stats::model.matrix(~ 0 + group)
    dimnames(X) <- list(samples, c("A", "B"))
    cc <- c(-1, 1)
    if (!is.null(W)) {
        X <- cbind(X, W[samples, , drop = FALSE])
        cc <- c(cc, rep(0, ncol(W)))
    }
    features <- intersect(features, rownames(M))
    v <- rep(NA_real_, length(features))
    df <- rep(NA_real_, length(features))
    for (i in seq_along(features)) {
        y <- M[features[[i]], samples]
        ok <- is.finite(y)
        Xo <- X[ok, , drop = FALSE]
        if (sum(ok) <= ncol(Xo)) next
        fit <- stats::lm.fit(Xo, y[ok])
        if (fit$rank < ncol(Xo)) next
        df[[i]] <- sum(ok) - fit$rank
        v[[i]] <- sum(fit$residuals^2) / df[[i]] * drop(crossprod(cc, solve(crossprod(Xo), cc)))
    }
    list(sum = sum(v, na.rm = TRUE), n = sum(is.finite(v)), df = if (any(is.finite(df))) stats::median(df, na.rm = TRUE) else NA_real_)
}

# ---- the row of a journey (section 6.1; the columns are all the scorer reads) --------------------------------------------
DONE <- c("ok", "refused_after_exclusion", "dropped", "uncorrected", "no_da")
ROW_COLUMNS <- c("tier", "design", "cell", "omic", "r", "version", "arm", "smoke", "status", "class", "step", "seconds", "warnings",
    "k_used", "k_found", "k_cap", "pct", "n_controls", "ruv_status_y", "ruv_status_x", "n_runs", "n_runs_final", "ruv_status_product_y",
    "corr_removed_y", "corr_removed_x", "e5_y", "e5_x",
    "cy_pre_a", "cy_post_a", "cx_pre_a", "cx_post_a", "qa_pre", "qa_post",
    paste0("fa_pre_", TAGS), paste0("fa_post_", TAGS), paste0("ff_pre_", TAGS), paste0("ff_post_", TAGS), "ffx_pre_05", "ffx_post_05",
    "pr", "kappa_c",
    "n_tested_y", "tested_sig_y", "any_call_y", "null_p05_y", "se2_mod_y", "n_se_y", "se2_unmod_y", "n_unmod_y", "df_y",
    "n_tested_x", "any_call_x", "null_p05_x",
    "dot_g", "dot_p", "dot_t", "dot_two", "den_g", "den_p", "den_t", "den_two",
    "pe_n_excluded", "pe_a_excluded", "pe_excluded_x")
newRow <- function(...) {
    given <- list(...)
    stopifnot(all(names(given) %in% ROW_COLUMNS))
    row <- stats::setNames(as.list(rep(NA, length(ROW_COLUMNS))), ROW_COLUMNS)
    row[names(given)] <- given
    row
}
rowFrame <- function(row) {
    stopifnot(setequal(names(row), ROW_COLUMNS))   # every row of a file has the same columns in the same order
    row <- row[ROW_COLUMNS]
    row[] <- lapply(row, function(x) if (length(x) != 1L) paste(x, collapse = ";") else x)
    as.data.frame(row, stringsAsFactors = FALSE)
}
daFamily <- function(j, assay) {
    da <- j[["da"]]
    if (is.null(da) || !nrow(da)) return(NULL)
    da[da$assay == assay & startsWith(da$feature, "F"), , drop = FALSE]
}
daNumbers <- function(d, shares) {
    ok <- is.finite(d$p)
    se <- abs(d$logFC / d$t)
    use_se <- is.finite(se) & is.finite(d$t) & d$t != 0
    number <- as.integer(sub("^F", "", d$feature))
    dots <- lapply(shares[c("g", "p", "t", "two")], function(s) sum(s[["x"]][d$feature[ok]] * d$logFC[ok], na.rm = TRUE))
    list(n = sum(ok), sig = sprintf("%d_%.0f_%.0f", sum(ok), sum(number[ok]), sum(number[ok]^2)),
        any = as.integer(any(is.finite(d$q) & d$q < 0.05)), p05 = if (any(ok)) mean(d$p[ok] < 0.05) else NA_real_,
        se2 = sum(se[use_se]^2), nse = sum(use_se), dots = dots)
}
#' The row of one journey. `ctx`: tier, design, cell, omic, r, per, shares, smoke. `ruv`: whether the request asked for RUV.
journeyRow <- function(ctx, j, sim, version, arm, ruv, filter_read = TRUE, extra = list()) {
    a_run <- sprintf("R%03d", 2L * ctx$per)
    has_a <- a_run %in% sim$runs
    others <- setdiff(sim$runs, a_run)
    pre <- j[["pairs_pre"]]
    post <- if (!ruv) pre else j[["pairs_post"]]
    cy_pre <- bestCorrelation(pre[[Y_ASSAY]], sim$runs)
    cy_post <- bestCorrelation(post[[Y_ASSAY]], sim$runs)
    cx_pre <- bestCorrelation(pre[[X_ASSAY]], sim$runs)
    cx_post <- bestCorrelation(post[[X_ASSAY]], sim$runs)
    row <- newRow(tier = ctx$tier, design = ctx$design, cell = ctx$cell, omic = ctx$omic, r = ctx$r, version = version, arm = arm,
        smoke = ctx$smoke, status = j$status, class = j$class, step = j$step, seconds = round(j$seconds, 2), warnings = j$warnings,
        n_runs = length(sim$runs), n_runs_final = if (is.null(j[["runs_final"]])) NA_integer_ else length(j[["runs_final"]]))
    row[["cy_pre_a"]] <- if (has_a) cy_pre[[a_run]] else NA_real_
    row[["cy_post_a"]] <- if (has_a) cy_post[[a_run]] else NA_real_
    row[["cx_pre_a"]] <- if (has_a) cx_pre[[a_run]] else NA_real_
    row[["cx_post_a"]] <- if (has_a) cx_post[[a_run]] else NA_real_
    smallest <- function(c_s) if (has_a && all(is.finite(c_s)) && length(others)) as.integer(c_s[[a_run]] < min(c_s[others])) else NA_integer_
    row[["qa_pre"]] <- smallest(cy_pre)
    row[["qa_post"]] <- smallest(cy_post)
    for (i in seq_along(THETAS)) {
        row[[paste0("fa_pre_", TAGS[[i]])]] <- if (has_a && !is.null(pre)) as.integer(failsAt(cy_pre[[a_run]], THETAS[[i]])) else NA_integer_
        row[[paste0("fa_post_", TAGS[[i]])]] <- if (has_a && !is.null(post)) as.integer(failsAt(cy_post[[a_run]], THETAS[[i]])) else NA_integer_
        row[[paste0("ff_pre_", TAGS[[i]])]] <- if (!is.null(pre)) sum(failsAt(cy_pre[others], THETAS[[i]])) else NA_integer_
        row[[paste0("ff_post_", TAGS[[i]])]] <- if (!is.null(post)) sum(failsAt(cy_post[others], THETAS[[i]])) else NA_integer_
    }
    row[["ffx_pre_05"]] <- if (!is.null(pre)) sum(failsAt(cx_pre[others], 0.5)) else NA_integer_
    row[["ffx_post_05"]] <- if (!is.null(post)) sum(failsAt(cx_post[others], 0.5)) else NA_integer_
    # the RUV choice
    ch <- if (!isTRUE(ruv)) NULL else if (isDeclared(arm)) declaredChoice(j) else legacyChoice(j)
    cols <- j[["assay_cols"]]
    cy <- ch[[Y_ASSAY]]
    row[["k_used"]] <- if (is.null(cy)) NA_integer_ else (cy[["k"]] %||% NA_integer_)
    row[["k_found"]] <- if (is.null(cy) || isDeclared(arm)) NA_integer_ else (cy[["k_sweep"]] %||% NA_integer_)
    row[["k_cap"]] <- if (is.null(cols[[Y_ASSAY]])) NA else (cols[[Y_ASSAY]][["ruv_k_cap"]] %||% NA)
    row[["pct"]] <- if (is.null(cy)) NA_real_ else (cy[["pct"]] %||% NA_real_)
    row[["n_controls"]] <- if (is.null(cy)) NA_integer_ else length(cy[["controls"]])
    row[["ruv_status_y"]] <- if (is.null(cols[[Y_ASSAY]])) NA_character_ else as.character(cols[[Y_ASSAY]][["ruv_status"]])
    row[["ruv_status_product_y"]] <- if (is.null(cols[[Y_ASSAY]])) NA_character_ else as.character(cols[[Y_ASSAY]][["ruv_status_product"]])
    row[["ruv_status_x"]] <- if (is.null(cols[[X_ASSAY]])) NA_character_ else as.character(cols[[X_ASSAY]][["ruv_status"]])
    # E5: the fail list the pairs table the filter read gives, against the package filter's own removed set (per assay)
    removed <- j[["corr_removed"]]
    read_y <- if (isTRUE(ruv)) cy_post else cy_pre
    read_x <- if (isTRUE(ruv)) cx_post else cx_pre
    row[["corr_removed_y"]] <- if (is.null(removed)) NA_character_ else paste(sort(removed[[Y_ASSAY]]), collapse = ";")
    row[["corr_removed_x"]] <- if (is.null(removed)) NA_character_ else paste(sort(removed[[X_ASSAY]]), collapse = ";")
    row[["e5_y"]] <- if (is.null(removed) || is.null(pre)) NA else setequal(failList(read_y, ctx$threshold), removed[[Y_ASSAY]])
    row[["e5_x"]] <- if (is.null(removed) || is.null(pre)) NA else setequal(failList(read_x, ctx$threshold), removed[[X_ASSAY]])
    # PR and kappa_c on the matrix the filter reads
    my <- if (isTRUE(ruv)) j[["m_out"]][[Y_ASSAY]] else j[["m_in"]][[Y_ASSAY]]
    row[["pr"]] <- if (is.null(my)) NA_real_ else pseudoReplicateRatio(my, sim$design, a_run)
    delta <- ctx$shares[[switch(version, V2 = "two", V1p = "p", V1s = "s", V1t = "t", "g")]][["delta"]]
    row[["kappa_c"]] <- if (is.null(cy) || !version %in% c("V1g", "V1p", "V1s", "V1t", "V2")) NA_real_ else kappaControls(delta, cy[["controls"]])
    # the DA of Y and X
    dy <- daFamily(j, Y_ASSAY)
    if (!is.null(dy) && nrow(dy)) {
        num <- daNumbers(dy, ctx$shares)
        row[["n_tested_y"]] <- num$n
        row[["tested_sig_y"]] <- num$sig
        row[["any_call_y"]] <- num$any
        row[["null_p05_y"]] <- num$p05
        row[["se2_mod_y"]] <- num$se2
        row[["n_se_y"]] <- num$nse
        for (k in names(num$dots)) row[[paste0("dot_", k)]] <- num$dots[[k]]
        mat <- if (isTRUE(ruv) && !isDeclared(arm)) j[["m_out"]][[Y_ASSAY]] else j[["m_in"]][[Y_ASSAY]]
        if (!is.null(mat) && !is.null(removed)) mat[, intersect(removed[[Y_ASSAY]], colnames(mat))] <- NA_real_
        w <- if (isDeclared(arm)) j[["ruv_factors"]][[Y_ASSAY]][["W"]] else NULL
        un <- tryCatch(unmoderatedReading(mat, sim$design, dy$feature[is.finite(dy$p)], w), error = function(e) list(sum = NA_real_, n = 0L, df = NA_real_))
        row[["se2_unmod_y"]] <- un$sum
        row[["n_unmod_y"]] <- un$n
        row[["df_y"]] <- un$df
    }
    dx <- daFamily(j, X_ASSAY)
    if (!is.null(dx) && nrow(dx)) {
        row[["n_tested_x"]] <- sum(is.finite(dx$p))
        row[["any_call_x"]] <- as.integer(any(is.finite(dx$q) & dx$q < 0.05))
        row[["null_p05_x"]] <- mean(dx$p[is.finite(dx$p)] < 0.05)
    }
    for (k in c("g", "p", "t", "two")) row[[paste0("den_", k)]] <- sum(ctx$shares[[k]][["x"]]^2)
    for (name in names(extra)) row[[name]] <- extra[[name]]
    row
}
#' LX (4): the legacy test (the DA kernel pair, `~ 0 + group`) on the corrected matrix of a declared journey (GK2 on V1g, GK1 on VD).
lxRow <- function(ctx, j, sim, version) {
    m <- j[["m_out"]][[Y_ASSAY]]
    row <- newRow(tier = ctx$tier, design = ctx$design, cell = ctx$cell, omic = ctx$omic, r = ctx$r, version = version, arm = "LX",
        smoke = ctx$smoke, status = "ok", class = NA_character_, step = NA_character_, n_runs = length(sim$runs))
    result <- tryCatch({
        fit <- quiet(fitLipidDaAssayModel(m, as.data.frame(sim$design), "~ 0 + group", "Run"))
        tab <- quiet(applyLipidDaContrasts(fit, "groupB-groupA", NA, TRUE, TRUE))$results[["groupB-groupA"]]
        data.frame(assay = Y_ASSAY, feature = rownames(tab), logFC = tab$logFC, t = tab$t, p = tab$P.Value, q = tab$fdr_qvalue,
            stringsAsFactors = FALSE)
    }, error = function(e) e)
    if (inherits(result, "error") || is.null(m)) {
        row[["status"]] <- "crashed"
        row[["class"]] <- if (inherits(result, "error")) oneLine(conditionMessage(result)) else "no corrected matrix"
        return(row)
    }
    dy <- result[startsWith(result$feature, "F"), , drop = FALSE]
    num <- daNumbers(dy, ctx$shares)
    row[["n_tested_y"]] <- num$n
    row[["tested_sig_y"]] <- num$sig
    row[["any_call_y"]] <- num$any
    row[["null_p05_y"]] <- num$p05
    row[["se2_mod_y"]] <- num$se2
    row[["n_se_y"]] <- num$nse
    for (k in names(num$dots)) row[[paste0("dot_", k)]] <- num$dots[[k]]
    un <- unmoderatedReading(m, sim$design, dy$feature[is.finite(dy$p)], NULL)
    row[["se2_unmod_y"]] <- un$sum
    row[["n_unmod_y"]] <- un$n
    row[["df_y"]] <- un$df
    for (k in c("g", "p", "t", "two")) row[[paste0("den_", k)]] <- sum(ctx$shares[[k]][["x"]]^2)
    row
}

# ---- progress files: checks and rows are appended as each unit runs, so a killed job still leaves them (lesson 3) -------------
appendTable <- function(path, df) {
    dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
    new <- !file.exists(path)
    utils::write.table(df, path, sep = "\t", quote = FALSE, row.names = FALSE, col.names = new, append = !new, na = "NA")
}
checks_env$written <- 0L
flushChecks <- function(path) {
    n <- length(checks_env$rows)
    if (n > checks_env$written) {
        appendTable(path, do.call(rbind, checks_env$rows[(checks_env$written + 1L):n]))
        checks_env$written <- n
    }
}

#' E12: under LG's marker the control method, called directly on the journey's RUV input, returns the first round(rows x 5 / 100) rows
#' in table order, in every assay (KIT-188's check of its LP arm, sim-k188.R:1136-1148).
checkE12 <- function(ctx, sim, j, arm, version) {
    m_in <- j[["m_in"]]
    condition <- paste(version, arm)
    if (is.null(m_in)) return(addCheck("E12", ctx$cell, ctx$r, condition, NA, "positional controls", FALSE, "no RUV input"))
    pct <- ctx$defaults$manual_pct
    object <- omicObjectFromMatrices(sim, m_in)
    got <- tryCatch(withMarkers(armMarkers(arm), quiet(getNegCtrlMetabAnova(object, ruv_grouping_variable = "group",
        percentage_as_neg_ctrl = pct))), error = function(e) if (grepl("partial-match guard", conditionMessage(e), fixed = TRUE)) stop(e) else e) # C04
    if (inherits(got, "error")) return(addCheck("E12", ctx$cell, ctx$r, condition, NA, "positional controls", FALSE,
        paste(class(got)[[1L]], conditionMessage(got))))
    same <- vapply(names(m_in), function(a) setequal(names(got[[a]])[which(got[[a]])], positional(m_in[[a]], pct)), logical(1L))
    addCheck("E12", ctx$cell, ctx$r, condition, NA, "first round(rows x 5 / 100) rows in table order", all(same), paste(names(same), same, collapse = "; "))
}

#' S02: after an exclusion, a refusal is `refused_after_exclusion` only when too few samples are left (VAL-347) and the refusal's class says so.
tooFewLeft <- function(sim, k = 2L) {
    counts <- table(sim$design[["group"]])
    length(counts) < 2L || any(counts < 2L) || (nrow(sim$design) - length(counts) - 1L) < k
}
FEW_CLASSES <- c("ruv_k_above_limit", "ruv_k_refused", "omics_ruv_no_residual_df", "omics_ruv_insufficient_controls", "multischolar_ruv_grouping_unusable")
#' S04: what L1 accepts as a completed journey: ok; a lost group (`no_da`, R-1); the emulation's results after an exclusion.
l1Accepts <- function(row, sim, arm) {
    status <- row[["status"]]
    class <- as.character(row[["class"]])
    identical(status, "ok") || (identical(status, "no_da") && !is.na(class) && startsWith(class, "filter_removed_group:")) ||
        (arm %in% c("PEL2", "PEG2") && status %in% c("refused_after_exclusion", "dropped", "uncorrected") && tooFewLeft(sim)) ||
        (identical(row[["tier"]], "V") && status %in% c("uncorrected", "dropped"))   # RG6: KIT-242 keeps an assay without RUV (a result)
}
#' One journey: the request of the arm, the journey, its row (written at once), its lifecycle checks.
runArm <- function(ctx, v, version, arm, lifecycle, write, correlation = TRUE, fault = FALSE, extra = list(), note_refusal = NULL) {
    sim <- v[["sim"]]
    request <- requestK275(arm, ctx$defaults, sim)
    j <- journeyL(sim, v[["values"]], arm, request, ctx$defaults, tempfile(sprintf("k275-%s-%s-", version, arm)), filter = TRUE,
        itsd = FALSE, correlation = correlation, fault = fault)
    ruv <- !is.null(request[["ruv"]])
    row <- journeyRow(ctx, j, sim, version, arm, ruv, extra = extra)
    if (identical(row[["status"]], "instrument") && startsWith(as.character(row[["class"]]), "partial_match_guard:")) {   # a guard stop inside a journey: leave its trace
        tag <- sub("\\.tsv$", "", sub("^checks-", "", basename(Sys.getenv("K275_CHECKS_PATH"))))
        try(writeLines(c(as.character(row[["class"]]), getOption("k275.guard_trace", "")), file.path(dirname(Sys.getenv("K275_CHECKS_PATH")), sprintf("unit-trace-%s-journey.txt", tag))), silent = TRUE)
    }
    if (isTRUE(note_refusal) && identical(row[["status"]], "refused") && tooFewLeft(sim) && row[["class"]] %in% FEW_CLASSES) row[["status"]] <- "refused_after_exclusion"
    write(row)
    condition <- paste(version, arm)
    if (lifecycle) {
        if (isTRUE(j[["guard_planted"]])) {
            addCheck("GUARD_planted", ctx$cell, ctx$r, condition, row[["status"]], "instrument", identical(row[["status"]], "instrument") &&
                startsWith(as.character(row[["class"]]), "partial_match_guard:"), "")
        } else if (!isTRUE(fault)) addCheck("L1_status", ctx$cell, ctx$r, condition, row[["status"]], "ok or a stated result", l1Accepts(row, sim, arm),
            paste(row[["step"]], row[["class"]]))
        for (e in c("e5_y", "e5_x")) if (!is.na(row[[e]])) addCheck("E5", ctx$cell, ctx$r, paste(condition, e), NA, "fail list equals the filter's removed set",
            isTRUE(row[[e]]))
    }
    if (startsWith(arm, "LG") && !isTRUE(fault) && (lifecycle || ctx$r <= 20L) && identical(j$status, "ok")) checkE12(ctx, sim, j, arm, version)
    list(j = j, row = row)
}

#' The pre-RUV emulation (VAL-344 OQ-3, VAL-347): pass 1 is the NO journey of the version, its pairs table at 0.5 gives the samples failed
#' in Y; pass 2 is the version VE (those samples removed from the design and both assays, as VD is built), the post-RUV filter off,
#' then LM2 (PEL2) or GK2 (PEG2) and the DA; only Y is read.
runEmulation <- function(ctx, getV, lifecycle, write, fault_now) {
    a_run <- sprintf("R%03d", 2L * ctx$per)
    first <- TRUE
    for (version in c("V0", "V1g")) {
        v <- getV(version)
        one <- runArm(ctx, v, version, "NO", lifecycle, write, fault = fault_now && first)
        first <- FALSE
        pre <- one$j[["pairs_pre"]]
        c_y <- bestCorrelation(pre[[Y_ASSAY]], v$sim$runs)
        c_x <- bestCorrelation(pre[[X_ASSAY]], v$sim$runs)
        failed <- failList(c_y, ctx$threshold)
        ve <- list(sim = dropSamples(v$sim, failed), values = v$values[, !colnames(v$values) %in% failed, drop = FALSE])
        if (lifecycle) {
            removed <- one$j[["corr_removed"]][[Y_ASSAY]]
            addCheck("E11", ctx$cell, ctx$r, paste(version, "pass 1"), NA, "fail list equals the NO journey's filter set",
                !is.null(removed) && setequal(failed, removed))
            addCheck("E11", ctx$cell, ctx$r, paste(version, "VE"), NA, "VE equals the version without those samples (exact)",
                identical(ve$values, v$values[, setdiff(colnames(v$values), failed), drop = FALSE]) &&
                    identical(ve$sim$design$Run, setdiff(v$sim$design$Run, failed)) &&
                    identical(ve$sim$design, v$sim$design[!v$sim$design$Run %in% failed, , drop = FALSE]))
        }
        extra <- list(pe_n_excluded = length(failed), pe_a_excluded = a_run %in% failed, pe_excluded_x = sum(failsAt(c_x, ctx$threshold)))
        for (arm in c("PEL2", "PEG2")) {
            two <- runArm(ctx, ve, version, arm, lifecycle, write, correlation = FALSE, extra = extra, note_refusal = length(failed) > 0L)
            if (lifecycle && identical(two$row[["status"]], "refused")) addCheck("E11", ctx$cell, ctx$r, paste(version, arm), NA,
                "a refusal that is not 'too few samples left' is an instrument defect", FALSE, two$row[["class"]])
        }
    }
}

#' The data checks of tier D's unit (E1, E2): determinism and the planting line, the version VD, and the other assay.
checkVersionsE1 <- function(ctx, getV, journeys) {
    cellname <- ctx$cell
    v0 <- getV("V0")
    v1 <- getV("V1g")
    vd <- getV("VD")
    n_runs <- CELLS[[cellname]]$groups * ctx$per
    a_run <- sprintf("R%03d", n_runs)
    same_twice <- identical(generate(cellname, "N1", ctx$r), generate(cellname, "N1", ctx$r)) &&
        identical(generate(cellname, "N1", ctx$r, list(form = "g", samples = n_runs)), generate(cellname, "N1", ctx$r, list(form = "g", samples = n_runs)))
    addCheck("E1", cellname, ctx$r, "determinism", NA, "byte-identical", same_twice)
    y0 <- v0$sim$y
    y1 <- v1$sim$y
    planted_rows <- which(v1$sim$assay_of[rownames(y1)] == Y_ASSAY)
    planted_cell <- matrix(FALSE, nrow(y1), ncol(y1))
    planted_cell[planted_rows, match(a_run, v1$sim$runs)] <- TRUE
    dev <- if (identical(dim(y0), dim(y1))) max(abs((y1 - y0)[!planted_cell])) else NA_real_
    addCheck("E1", cellname, ctx$r, "V1g equals V0 except the planted entries", dev, "exactly 0", isTRUE(dev == 0) &&
        identical(v0$sim$oracle, v1$sim$oracle) && any(y1 != y0))
    xs <- ctx$shares
    g <- xs[["g"]][["delta"]]
    t_delta <- xs[["t"]][["delta"]]
    s_delta <- xs[["s"]][["delta"]]
    t_nz <- names(t_delta)[t_delta != 0]
    addCheck("E1", cellname, ctx$r, "V1t is V1g on its 25 features and 0 elsewhere", length(t_nz), "25 features; equal to V1g", length(t_nz) == 25L &&
        identical(t_delta[t_nz], g[t_nz]) && all(t_delta[setdiff(names(t_delta), t_nz)] == 0))
    lab <- grepl("^(LBL|CTL)_", names(s_delta))
    addCheck("E1", cellname, ctx$r, "V1s is V1g on the features and 0 on the labelled controls", NA, "exact", identical(s_delta[!lab], g[!lab]) &&
        all(s_delta[lab] == 0))
    addCheck("E1", cellname, ctx$r, "VD is V0 without a (exact)", NA, "exact", identical(vd$values, v0$values[, setdiff(colnames(v0$values), a_run), drop = FALSE]) &&
        identical(vd$sim$design, v0$sim$design[v0$sim$design$Run != a_run, , drop = FALSE]))
    n0 <- journeys[["V0 NO"]]
    n1 <- journeys[["V1g NO"]]
    if (!is.null(n0) && !is.null(n1) && !is.null(n0[["m_in"]]) && !is.null(n1[["m_in"]])) {
        addCheck("E1", cellname, ctx$r, "the planted sample's Y values entering the pipeline differ, X's do not", NA, "Y differs, X identical",
            !identical(n0[["m_in"]][[Y_ASSAY]], n1[["m_in"]][[Y_ASSAY]]) && identical(n0[["m_in"]][[X_ASSAY]], n1[["m_in"]][[X_ASSAY]]))
        key <- function(j) {
            da <- j[["da"]]
            x <- da[da$assay == X_ASSAY, , drop = FALSE]
            x <- x[order(x$feature), , drop = FALSE]
            rownames(x) <- NULL
            x
        }
        addCheck("E2", cellname, ctx$r, "NO: the DA of X is identical across V0 and V1g", NA, "exact", !is.null(n0[["da"]]) && !is.null(n1[["da"]]) &&
            identical(key(n0), key(n1)))
    }
}

#' One data set of one tier (a child process): the versions, the arms, the rows and the checks, written as it goes.
runUnit <- function(tier, design, r, lifecycle, fault_now, e7r, smoke, rows_path, checks_path) {
    cellname <- if (identical(tier, "ID")) "D4M" else design
    p <- CELLS[[cellname]]
    defaults <- routeDefaults(p$omic)
    ctx <- list(tier = tier, design = design, cell = cellname, omic = p$omic, r = r, per = p$per, smoke = smoke, threshold = defaults$correlation,
        defaults = defaults, shares = plantedShares(cellname, r))
    mask <- if (identical(tier, "V")) "vd5" else "none"
    if (isTRUE(getOption("k275.guard_unit"))) invisible(data.frame(abcdef = 1)$abc)   # the planted unit-level partial match (outside any journey)
    v0 <- makeVersion(cellname, r, "V0", mask)
    versions <- list(V0 = v0)
    getV <- function(version) {
        if (is.null(versions[[version]])) versions[[version]] <<- makeVersion(cellname, r, version, mask, v0)
        versions[[version]]
    }
    write <- function(row) {
        appendTable(rows_path, rowFrame(row))
        flushChecks(checks_path)
    }
    journeys <- list()
    if (identical(tier, "E")) {
        runEmulation(ctx, getV, lifecycle, write, fault_now)
    } else {
        plan <- TIER_PLAN[[tier]]
        first <- TRUE
        for (i in seq_len(nrow(plan))) {
            version <- plan$version[[i]]
            arm <- plan$arm[[i]]
            one <- runArm(ctx, getV(version), version, arm, lifecycle, write, fault = fault_now && first)
            first <- FALSE
            journeys[[paste(version, arm)]] <- one$j
        }
        if (identical(tier, "D")) {
            # LX: the legacy test on the declared design's corrected matrix, on V1g (from GK2) and on VD (from GK1)
            for (pair in list(c("V1g", "GK2"), c("VD", "GK1"))) {
                j <- journeys[[paste(pair[[1L]], pair[[2L]])]]
                if (is.null(j)) next
                lx <- lxRow(ctx, j, getV(pair[[1L]])[["sim"]], pair[[1L]])
                write(lx)
                if (lifecycle) addCheck("L1_status", cellname, r, paste(pair[[1L]], "LX"), lx[["status"]], "ok", identical(lx[["status"]], "ok"), lx[["class"]])
            }
            checkVersionsE1(ctx, getV, journeys)   # E1 and E2 are identity checks of the counted journeys too (F4)
        }
    }
    if (identical(tier, "V")) {   # S14: tier V's mask clauses of E1
        v1 <- getV("V1g")
        a_col <- match(sprintf("R%03d", 2L * ctx$per), v0$sim$runs)
        m1 <- maskOf(v1$sim, "vd5", v0$sim$y[v0$sim$feat_ids, , drop = FALSE])
        m0 <- maskOf(v0$sim, "vd5")
        a_cell <- matrix(FALSE, nrow(m1), ncol(m1))
        a_cell[v1$sim$assay_of[rownames(m1)] == Y_ASSAY, a_col] <- TRUE
        addCheck("E1", cellname, r, "V mask: V1g equals V0's except in a's Y column", NA, "exact", identical(m1[!a_cell], m0[!a_cell]))
        vd <- getV("VD")
        addCheck("E1", cellname, r, "V mask: VD is V0's masked values without a", NA, "exact", identical(vd$values, v0$values[, setdiff(colnames(v0$values), colnames(v0$values)[a_col]), drop = FALSE]))
    }
    # E7: the journey the package refuses (GK3 at D3 on VD: k = 3 above the limit 5 - 2 - 1 = 2); not a denominator, written as a refused row
    if (lifecycle && e7r) {
        one <- runArm(ctx, getV("VD"), "VD", "GK3", FALSE, write)
        addCheck("E7_refusal", cellname, r, "VD GK3", one$row[["status"]], "refused, class ruv_k_above_limit",
            identical(one$row[["status"]], "refused") && identical(one$row[["class"]], "ruv_k_above_limit"))
    }
    if (lifecycle && fault_now) {
        file_rows <- utils::read.delim(rows_path, stringsAsFactors = FALSE, colClasses = "character", na.strings = "NA")
        addCheck("E7_crash", cellname, r, "first journey of the first unit", NA, "recorded crashed, with its step", any(file_rows$status == "crashed" &
            grepl("deliberate fault", file_rows$class, fixed = TRUE) & nzchar(file_rows$step)))
    }
    flushChecks(checks_path)
    invisible(NULL)
}

# ---- the child process ---------------------------------------------------------------------------------------------------
runChildMain <- function() {
    stopifnot(length(args) == 5L)
    flag <- function(name) identical(Sys.getenv(name), "1")
    rows_path <- Sys.getenv("K275_ROWS_PATH")
    checks_path <- Sys.getenv("K275_CHECKS_PATH")
    last_calls <- NULL
    options(k275.guard_now = flag("K275_GUARD_UNIT"), k275.guard_unit = flag("K275_GUARD_UNIT_LEVEL"))
    result <- tryCatch(withCallingHandlers({
        runUnit(args[[2L]], args[[3L]], as.integer(args[[4L]]), flag("K275_LIFECYCLE"), flag("K275_FAULT_UNIT"), flag("K275_E7R_UNIT"),
            flag("K275_SMOKE"), rows_path, checks_path)
        list(ok = TRUE)
    }, warning = guardWarning, error = function(e) last_calls <<- paste(vapply(sys.calls(), function(cl) oneLine(deparse(cl)[[1L]], 80L), ""), collapse = " <- ")),
        error = function(e) list(error = oneLine(paste(class(e)[[1L]], conditionMessage(e))), trace = last_calls %||% getOption("k275.guard_trace", "")))
    flushChecks(checks_path)
    saveRDS(result, args[[5L]])
    quit(save = "no", status = 0L)
}

# ---- the job (parent) -------------------------------------------------------------------------------------------------------
readText <- function(path) utils::read.delim(path, stringsAsFactors = FALSE, colClasses = "character", na.strings = "NA")
runParentMain <- function() {
    if (length(args) != 5L) stop("usage: Rscript sim-k275.R <tiers> <designs> <first> <last> <out dir>")
    lifecycle <- identical(Sys.getenv("K275_LIFECYCLE"), "1")
    smoke <- identical(Sys.getenv("K275_SMOKE"), "1")
    fault <- identical(Sys.getenv("K275_FAULT"), "1")
    # the smoke send's flag may never reach a counted request
    if (smoke && !lifecycle) stop("SMOKE_counted_refusal: K275_SMOKE=1 is refused without K275_LIFECYCLE=1 (a counted request never carries the smoke flag)")
    # `all` is every tier but V unless tier V is bound to its own commit (L1-V runs tier V alone, section 10)
    tiers <- if (identical(args[[1L]], "all")) (if (grepl("^[0-9a-f]{40}$", Sys.getenv("K275_TIER_V_BOUND"))) TIERS else setdiff(TIERS, "V")) else
        strsplit(args[[1L]], ",", fixed = TRUE)[[1L]]
    designs <- strsplit(args[[2L]], ",", fixed = TRUE)[[1L]]
    first_rep <- as.integer(args[[3L]])
    last_rep <- as.integer(args[[4L]])
    stopifnot(all(tiers %in% TIERS), all(designs %in% c("D3", "D4", "D6")), !is.na(first_rep), !is.na(last_rep), first_rep <= last_rep)
    # tier V waits for a commit holding KIT-242 Parts A and B (section 13): it runs only when that commit is named
    if ("V" %in% tiers && !grepl("^[0-9a-f]{40}$", Sys.getenv("K275_TIER_V_BOUND"))) {
        stop("tier V waits for KIT-242 Parts A and B (ticket section 13): set K275_TIER_V_BOUND to the 40-digit commit that holds them")
    }
    workers <- as.integer(Sys.getenv("K275_WORKERS", "4"))
    job_name <- Sys.getenv("K275_JOB", sprintf("%s-%s-%d-%d%s", paste(tiers, collapse = "_"), paste(designs, collapse = "_"), first_rep, last_rep,
        if (lifecycle) "-l1" else ""))
    out_dir <- args[[5L]]
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    out_dir <- normalizePath(out_dir)
    # in the blind lifecycle the row files stay under tempdir() (nothing of the measurement leaves the runner)
    rows_dir <- if (lifecycle) file.path(tempdir(), "rows") else file.path(out_dir, "rows")
    progress_dir <- file.path(out_dir, "progress")
    dir.create(rows_dir, recursive = TRUE, showWarnings = FALSE)
    dir.create(progress_dir, recursive = TRUE, showWarnings = FALSE)
    only_r <- Sys.getenv("K275_ONLY_R")   # counted mode: the continuation of a timed-out job runs exactly its missing data sets, once
    if (nzchar(only_r) && lifecycle) stop("K275_ONLY_R is a counted-mode setting (the continuation); it is refused in the lifecycle")
    units <- do.call(rbind, lapply(tiers, function(tier) do.call(rbind, lapply(intersect(designs, tierDesigns(tier)), function(design)
        data.frame(tier = tier, design = design, r = seq(first_rep, last_rep), repeat_run = FALSE, stringsAsFactors = FALSE)))))
    if (is.null(units) || !nrow(units)) stop("K275: no unit: the tiers and designs given have no common cell")
    if (nzchar(only_r)) {
        keep_r <- as.integer(strsplit(only_r, ",", fixed = TRUE)[[1L]])
        units <- units[units$r %in% keep_r, , drop = FALSE]
        if (!nrow(units)) stop("K275: K275_ONLY_R names no data set of the range")
        if (!grepl("-cont$", job_name)) job_name <- paste0(job_name, "-cont")
    }
    if (lifecycle) {
        # E9: the first data set of every (tier, design) run a second time end to end
        again <- units[!duplicated(units[c("tier", "design")]), , drop = FALSE]
        repeat_tiers <- Sys.getenv("K275_REPEAT", "all")   # the smoke send repeats one small tier only; L1 repeats all
        if (!identical(repeat_tiers, "all")) again <- again[again$tier %in% strsplit(repeat_tiers, ",", fixed = TRUE)[[1L]], , drop = FALSE]
        if (nrow(again)) {
            again$repeat_run <- TRUE
            units <- rbind(units, again)
        }
    }
    journeys_per <- vapply(seq_len(nrow(units)), function(i) if (identical(units$tier[[i]], "E")) 6L else nrow(TIER_PLAN[[units$tier[[i]]]]), integer(1L))
    units <- units[order(-journeys_per, units$design, units$r), , drop = FALSE]
    units$tag <- sprintf("%s-%s-%d%s", units$tier, units$design, units$r, ifelse(units$repeat_run, "-rep", ""))
    units$fault <- FALSE
    units$e7r <- FALSE
    units$guard <- FALSE
    units$guard_unit <- FALSE
    if (lifecycle && identical(Sys.getenv("K275_FAULT_GUARD"), "1")) {   # the planted partial-match fault: the guard's own self-check
        hit <- which(units$tier == "D" & !units$repeat_run & units$r != first_rep)   # C01: every plant sits on a data set E9 does not repeat
        if (length(hit)) units$guard[[hit[[1L]]]] <- TRUE
    }
    if (lifecycle && nzchar(Sys.getenv("K275_FAULT_GUARD_UNIT"))) {   # the planted unit-level variant: on the first unit of the tier named (not the journey-level one)
        hit <- which(units$tier == Sys.getenv("K275_FAULT_GUARD_UNIT") & !units$guard & !units$repeat_run & units$r != first_rep)
        if (length(hit)) units$guard_unit[[hit[[1L]]]] <- TRUE
    }
    if (lifecycle && fault) {
        hit <- which(!units$repeat_run & units$r != first_rep)    # E7's fault: the first journey of the first unit that E9 does not repeat (C01)
        if (length(hit)) units$fault[[hit[[1L]]]] <- TRUE
        hit <- which(units$design == "D3" & units$tier == "D" & !units$repeat_run & units$r != first_rep)
        if (length(hit)) units$e7r[[hit[[1L]]]] <- TRUE
    }
    cat(sprintf("K275 %s: %d units, %d workers, lifecycle %s, smoke %s\n", job_name, nrow(units), workers, lifecycle, smoke))
    # E3: the route defaults at the bound commit, once per omic, before any unit (a mismatch stops the job)
    for (omic in unique(vapply(units$tier, function(t) if (identical(t, "ID")) "metab" else "lipid", ""))) {
        e3 <- checkE3(omic)
        addCheck("E3", omic, 0L, "defaults", NA, "equal to 5.2", e3$pass, paste(e3$note, "|",
            paste(names(e3$got), vapply(e3$got, function(x) paste(format(x), collapse = "/"), ""), sep = "=", collapse = "; ")))
    }
    flushChecks(file.path(out_dir, sprintf("checks-%s-start.tsv", job_name)))
    # a counted job stops on an E3 mismatch (the route under test changed); in the lifecycle it is a failed check
    if (!lifecycle && any(vapply(checks_env$rows, function(x) identical(x$check, "E3") && !x$pass, logical(1L)))) {
        stop("K275: E3 failed: the route defaults at the bound commit differ from section 5.2 (see checks-", job_name, "-start.tsv)")
    }
    declaration <- Sys.getenv("K275_DECLARATION")
    scorer_file <- file.path(dirname(script_file), "score-k275.R")
    library_table <- utils::installed.packages()[, c("Package", "Version")]
    key_packages <- c("MultiScholaR", "RUVIIIC", "ruv", "limma", "qvalue", "duckdb", "dplyr")
    provenance <- function(final) {
        tables <- if (final) list.files(out_dir, recursive = TRUE, full.names = TRUE, pattern = "\\.(tsv|gz|json)$") else character()
        list(job = job_name, tiers = tiers, designs = designs, data_sets = c(first_rep, last_rep), lifecycle = lifecycle, smoke = smoke,
            fault = fault, workers = workers, units = nrow(units), final = final,
            commit = tryCatch(system2("git", c("rev-parse", "HEAD"), stdout = TRUE), error = function(e) NA),
            r_version = R.version.string, blas = extSoftVersion()[["BLAS"]], lapack = La_library(), omp_threads = Sys.getenv("OMP_NUM_THREADS"),
            runner_image = paste(Sys.getenv("ImageOS"), Sys.getenv("ImageVersion")),
            processor = tryCatch(sub("^model name\\s*:\\s*", "", grep("^model name", readLines("/proc/cpuinfo"), value = TRUE)[[1L]]), error = function(e) NA),
            cores = parallel::detectCores(),
            seed_rule = "set.seed(275000000 + 100000 c + 10000 s + r, kind = 'Mersenne-Twister', normal.kind = 'Inversion', sample.kind = 'Rejection'); planting 276000000 + 100000 c + r; bootstrap 277000000 + 100 c + t",
            key_packages = as.list(stats::setNames(library_table[match(key_packages, library_table[, "Package"]), "Version"], key_packages)),
            script_sha256 = sha(script_file), scorer_sha256 = if (file.exists(scorer_file)) sha(scorer_file) else NA,
            declaration_file = declaration, declaration_sha256 = if (nzchar(declaration) && file.exists(declaration)) sha(declaration) else NA,
            started = format(started_at, "%Y-%m-%dT%H:%M:%S%z"), finished = if (final) format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z") else NA,
            tables_sha256 = stats::setNames(as.list(vapply(tables, sha, character(1L))), sub(paste0("^", out_dir, "/?"), "", tables)))
    }
    # the provenance is written before any unit runs (lesson 3), and again at the end
    jsonlite::write_json(provenance(FALSE), file.path(out_dir, sprintf("provenance-%s.json", job_name)), auto_unbox = TRUE, pretty = TRUE, digits = NA)
    runChild <- function(i) {
        u <- units[i, ]
        rds <- file.path(progress_dir, sprintf("unit-%s.rds", u$tag))
        rows_path <- file.path(rows_dir, sprintf("part-%s.tsv", u$tag))
        checks_path <- file.path(progress_dir, sprintf("checks-%s.tsv", u$tag))
        env <- c(sprintf("K275_LIFECYCLE=%s", if (lifecycle) "1" else "0"), sprintf("K275_FAULT_UNIT=%s", if (u$fault) "1" else "0"),
            sprintf("K275_E7R_UNIT=%s", if (u$e7r) "1" else "0"), sprintf("K275_SMOKE=%s", if (smoke) "1" else "0"),
            sprintf("K275_GUARD_UNIT=%s", if (u$guard) "1" else "0"), sprintf("K275_GUARD_UNIT_LEVEL=%s", if (u$guard_unit) "1" else "0"), sprintf("K275_ROWS_PATH=%s", rows_path), sprintf("K275_CHECKS_PATH=%s", checks_path))
        t0 <- Sys.time()
        code <- system2("Rscript", c(shQuote(script_file), "--unit", u$tier, u$design, u$r, shQuote(rds)), env = env, stdout = FALSE, stderr = FALSE)
        seconds <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
        trace_file <- file.path(progress_dir, sprintf("unit-trace-%s.txt", u$tag))
        res <- if (file.exists(rds)) readRDS(rds) else list(error = paste(sprintf("child process exit %d, no result", code),
            if (file.exists(trace_file)) oneLine(paste(readLines(trace_file, n = 1L), collapse = " "), 200L) else ""))
        unlink(rds)
        if (!is.null(res[["trace"]])) writeLines(c(res[["error"]] %||% "", res[["trace"]]), trace_file)
        data.frame(tag = u$tag, tier = u$tier, design = u$design, r = u$r, wall_seconds = round(seconds, 1), exit = code,
            error = res[["error"]] %||% "", stringsAsFactors = FALSE)
    }
    results <- parallel::mclapply(seq_len(nrow(units)), runChild, mc.cores = workers, mc.preschedule = FALSE)
    timing <- do.call(rbind, lapply(seq_along(results), function(i) if (is.data.frame(results[[i]])) results[[i]] else
        data.frame(tag = units$tag[[i]], tier = units$tier[[i]], design = units$design[[i]], r = units$r[[i]], wall_seconds = NA, exit = NA,
            error = oneLine(paste("worker failed:", paste(as.character(results[[i]]), collapse = " "))), stringsAsFactors = FALSE)))
    # a unit that left no rows is recorded as a crashed unit row (never a denominator)
    for (i in seq_len(nrow(units))) {
        path <- file.path(rows_dir, sprintf("part-%s.tsv", units$tag[[i]]))
        if (!file.exists(path) || !nzchar(readLines(path, n = 1L))) {
            appendTable(path, rowFrame(newRow(tier = units$tier[[i]], design = units$design[[i]], r = units$r[[i]], arm = "unit", smoke = smoke,
                status = if (grepl("partial-match guard", timing$error[[i]], fixed = TRUE)) "instrument" else "crashed", class = timing$error[[i]])))
        }
    }
    # the checks: the parent's, then every unit's progress file
    unit_checks <- lapply(list.files(progress_dir, pattern = "^checks-.*\\.tsv$", full.names = TRUE), readText)
    parent_checks <- if (length(checks_env$rows)) do.call(rbind, lapply(checks_env$rows, function(x) as.data.frame(lapply(x, as.character), stringsAsFactors = FALSE)))
    checks <- do.call(rbind, c(list(parent_checks), unit_checks))
    # E9: the second run of each first data set equals the first (rows without their seconds)
    if (lifecycle) {
        for (i in which(units$repeat_run)) {
            twin <- which(units$tier == units$tier[[i]] & units$design == units$design[[i]] & units$r == units$r[[i]] & !units$repeat_run)
            a <- readText(file.path(rows_dir, sprintf("part-%s.tsv", units$tag[[twin]])))
            b <- readText(file.path(rows_dir, sprintf("part-%s.tsv", units$tag[[i]])))
            keep <- setdiff(names(a), "seconds")
            canon <- function(x) {
                x <- x[do.call(order, x[c("version", "arm")]), keep, drop = FALSE]
                rownames(x) <- NULL
                x
            }
            same <- identical(canon(a), canon(b))
            checks <- rbind(checks, data.frame(check = "E9", cell = units$tier[[i]], rep = units$r[[i]], condition = units$design[[i]], value = NA,
                threshold = "identical rows", pass = isTRUE(same), note = "the second run equals the first", stringsAsFactors = FALSE))
        }
    }
    # the unit rows' status: a unit that crashed or errored is a failed check in the lifecycle
    for (i in which(nzchar(timing$error))) checks <- rbind(checks, data.frame(check = "unit_completed", cell = timing$tier[[i]], rep = timing$r[[i]],
        condition = timing$design[[i]], value = NA, threshold = "ran to the end", pass = FALSE, note = timing$error[[i]], stringsAsFactors = FALSE))
    utils::write.table(timing, file.path(out_dir, sprintf("timing-units-%s.tsv", job_name)), sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
    if (lifecycle) {
        # L1 writes only checks, timing by arm class and design, hashes and the scorer's parse and coverage; no estimand
        tally <- do.call(rbind, lapply(split(checks, checks$check), function(x) data.frame(check = x$check[[1L]], units_compared = nrow(x),
            units_passed = sum(x$pass == "TRUE"), pass = all(x$pass == "TRUE"), max_abs_deviation = suppressWarnings(max(as.numeric(x[["value"]]),
                na.rm = TRUE)), stringsAsFactors = FALSE)))
        tally$max_abs_deviation[!is.finite(tally$max_abs_deviation)] <- NA
        utils::write.table(tally, file.path(out_dir, sprintf("checks-%s.tsv", job_name)), sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
        files <- list.files(rows_dir, pattern = "^part-.*\\.tsv$", full.names = TRUE)
        all_rows <- do.call(rbind, lapply(files, readText))
        all_rows$class_of_arm <- vapply(all_rows$arm, function(a) if (identical(a, "LX") || identical(a, "unit")) a else armClass(a), "")
        by <- split(all_rows, paste(all_rows$class_of_arm, all_rows$design))
        timing_arm <- do.call(rbind, lapply(by, function(x) data.frame(arm_class = x$class_of_arm[[1L]], design = x$design[[1L]], journeys = nrow(x),
            completed = sum(x$status %in% DONE), mean_seconds = round(mean(as.numeric(x$seconds[x$status %in% DONE]), na.rm = TRUE), 1))))
        timing_arm <- timing_arm[!timing_arm$arm_class %in% c("unit"), , drop = FALSE]
        # G2_coverage: a class with no completed journey at a design fails (L1 must measure every class; tier V only when bound)
        need <- expand.grid(arm_class = c("NO", "LM", "LG", "LA", "GA", "GK", "GSC", "PEL2", "PEG2"), design = unique(units$design[units$tier %in% c("D", "K", "S", "E", "T", "F")]), stringsAsFactors = FALSE)
        need <- need[(need$arm_class %in% c("LG") & need$design %in% unique(units$design[units$tier %in% "K"])) | (need$arm_class %in% c("GSC") & need$design %in% unique(units$design[units$tier %in% "S"])) |
            (need$arm_class %in% c("PEL2", "PEG2") & need$design %in% unique(units$design[units$tier %in% "E"])) | need$arm_class %in% c("NO", "LM", "LA", "GA", "GK"), , drop = FALSE]
        for (i in seq_len(nrow(need))) {
            hit <- timing_arm[timing_arm$arm_class == need$arm_class[[i]] & timing_arm$design == need$design[[i]], , drop = FALSE]
            checks <- rbind(checks, data.frame(check = "G2_coverage", cell = need$arm_class[[i]], rep = 0L, condition = need$design[[i]], value = NA, threshold = ">= 1 completed journey",
                pass = nrow(hit) > 0L && hit$completed[[1L]] > 0L, note = "", stringsAsFactors = FALSE))
        }
        tally <- do.call(rbind, lapply(split(checks, checks$check), function(x) data.frame(check = x$check[[1L]], units_compared = nrow(x),
            units_passed = sum(x$pass == "TRUE"), pass = all(x$pass == "TRUE"), max_abs_deviation = suppressWarnings(max(as.numeric(x[["value"]]), na.rm = TRUE)), stringsAsFactors = FALSE)))
        tally$max_abs_deviation[!is.finite(tally$max_abs_deviation)] <- NA
        utils::write.table(tally, file.path(out_dir, sprintf("checks-%s.tsv", job_name)), sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
        utils::write.table(timing_arm, file.path(out_dir, sprintf("timing-%s.tsv", job_name)), sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
        hashes <- data.frame(item = c(basename(files), "script", "scorer", "bound_commit"), value = c(vapply(files, sha, ""), sha(script_file),
            if (file.exists(scorer_file)) sha(scorer_file) else NA, tryCatch(system2("git", c("rev-parse", "HEAD"), stdout = TRUE), error = function(e) NA)))
        utils::write.table(hashes, file.path(out_dir, sprintf("hashes-%s.tsv", job_name)), sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
        # the scorer, blind, on this job's own rows: only parse success and coverage (and the identity check E8) are written
        blind <- system2("Rscript", c(shQuote(scorer_file), "--blind", shQuote(rows_dir), shQuote(out_dir), shQuote(job_name)), stdout = FALSE, stderr = FALSE)
        cat(sprintf("K275 %s: blind scorer exit %d\n", job_name, blind))
        failed <- tally[!tally$pass, , drop = FALSE]
        for (i in seq_len(nrow(failed))) {
            x <- checks[checks$check == failed$check[[i]] & checks$pass != "TRUE", , drop = FALSE]
            cat(sprintf("K275 %s: FAILED check %s (%d of %d); first: %s\n", job_name, failed$check[[i]], nrow(x), failed$units_compared[[i]],
                oneLine(paste(x$condition[[1L]], x$note[[1L]]), 160L)))
        }
        cat(sprintf("K275 %s: %d checks, %d failed\n", job_name, nrow(tally), nrow(failed)))
    } else {
        utils::write.table(checks, file.path(out_dir, sprintf("checks-%s.tsv", job_name)), sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
        for (tier in unique(units$tier)) for (design in unique(units$design[units$tier == tier])) {
            files <- list.files(rows_dir, pattern = sprintf("^part-%s-%s-.*\\.tsv$", tier, design), full.names = TRUE)
            if (!length(files)) next
            table <- do.call(rbind, lapply(files, readText))
            out_file <- file.path(out_dir, "rows", sprintf("%s-%s-%d-%d.tsv.gz", tier, design, first_rep, last_rep))
            connection <- gzfile(out_file, "w")
            utils::write.table(table, connection, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
            close(connection)
            unlink(files)
        }
    }
    jsonlite::write_json(provenance(TRUE), file.path(out_dir, sprintf("provenance-%s.json", job_name)), auto_unbox = TRUE, pretty = TRUE, digits = NA)
    cat(sprintf("K275 %s finished: %d units (%d with an error), wall %.1f min\n", job_name, nrow(units), sum(nzchar(timing$error)),
        as.numeric(difftime(Sys.time(), started_at, units = "mins"))))
    invisible(NULL)
}

# `options(k275.library = TRUE)` loads the functions only (fixture-k275.R); a job or a child runs below
if (!isTRUE(getOption("k275.library"))) {
    if (child) runChildMain() else runParentMain()
}
