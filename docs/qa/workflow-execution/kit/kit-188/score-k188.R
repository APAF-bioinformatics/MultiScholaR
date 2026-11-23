# KIT-188: the scoring (sections 5.5 to 5.7, 18.5), run once on a runner over the committed rows of S. Evidence code under
# docs/qa: no package code. Hashed in the declaration; a defect found after scoring is fixed in a dated section and all
# rows re-scored once (section 10, Attempts).
#
# Usage: Rscript score-k188.R <rows dir> <out dir> [<part L family> <part P family>]
#   rows dir  the assembled rows (*.tsv.gz) of every S job (deciding, reported, munit, ba, lp)
#   family    the Bonferroni family sizes (default 720 for Part L under D6 (a), 936 for Part P)
# Writes: conditions.tsv (per part, arm, cell, scenario, mask, estimand: R, mean, SD, Monte Carlo SE, limits,
#   classification), verdicts.tsv (deciding arm per omic group), outcome.txt, dataset.tsv (per data set and arm),
#   statuses.tsv, partm-df.tsv, partm-amp.tsv, partm-unit.tsv, partm-repair.tsv.

args <- commandArgs(TRUE)
stopifnot(length(args) >= 2L)
rows_dir <- args[[1L]]
out_dir <- args[[2L]]
family <- c(L = if (length(args) >= 3L) as.numeric(args[[3L]]) else 720, P = if (length(args) >= 4L) as.numeric(args[[4L]]) else 936)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

EDGES <- c(S1 = 0.075, S2 = 0.075, S3 = 0.075, S5 = 0.01)
F1_EDGES <- c(0.05, 0.20)
F2_EDGES <- c(0.005, 0.05)
DECIDING <- c("LA", "GA", "NO")
OMIC_GROUP <- c(L12 = "lipid_metab", L64 = "lipid_metab", M12 = "lipid_metab", M168 = "lipid_metab", P0 = "proteomics",
    P1 = "proteomics", P16 = "proteomics", P2 = "proteomics", P3 = "proteomics", P4b = "proteomics")

files <- list.files(rows_dir, pattern = "\\.tsv\\.gz$", full.names = TRUE, recursive = TRUE)
rows <- do.call(rbind, lapply(files, function(f) utils::read.delim(gzfile(f), stringsAsFactors = FALSE, na.strings = "NA")))
num <- function(x) suppressWarnings(as.numeric(x))
for (column in c("rep", "n_T", "n_tested", "lost", "n_null", "null_p05", "null_p01", "calls", "false_calls", "true_calls", "any_call",
    "fdp", "power", "n_planted", "bh_calls", "bh_false", "bh_any", "bh_fdp", "bh_power", "gross1", "n_common", "s5", "or_s5",
    "gross1_common", "or_gross1_common", "err_max", "k", "samples_removed")) rows[[column]] <- num(rows[[column]])
rows$arm_scale[is.na(rows$arm_scale)] <- "1"
main <- rows[rows$job %in% c("deciding", "ba", "reported", "lp") & rows$arm_scale == "1" & rows$scoreset %in% "T", , drop = FALSE]

# ---- per data set and arm (the unit, VAL-11) -------------------------------------------------------------------------
pooled <- function(num, den) if (sum(den, na.rm = TRUE) > 0) sum(num, na.rm = TRUE) / sum(den, na.rm = TRUE) else NA_real_
key <- c("part", "cell", "job", "scenario", "mask", "rep", "arm")
groups <- split(main, main[key], drop = TRUE)
dataset <- do.call(rbind, lapply(groups, function(x) {
    statuses <- unique(x$status)
    data.frame(x[1L, key], status = paste(sort(statuses), collapse = "+"),
        S1 = pooled(x$null_p05 * x$n_null, x$n_null),
        S1_01 = pooled(x$null_p01 * x$n_null, x$n_null),
        S2 = mean(x$any_call, na.rm = TRUE), S2_pooled = as.numeric(any(x$any_call > 0, na.rm = TRUE)),
        S3 = if (sum(x$calls, na.rm = TRUE) > 0) sum(x$false_calls, na.rm = TRUE) / sum(x$calls, na.rm = TRUE) else 0,
        S5 = if (identical(x$arm[[1L]], "OR")) max(x$or_s5, na.rm = TRUE) else suppressWarnings(max(x$s5, na.rm = TRUE)),
        power = pooled(x$true_calls, x$n_planted),
        gross_common = pooled(x$gross1_common * x$n_common, x$n_common),
        or_gross_common = pooled(x$or_gross1_common * x$n_common, x$n_common),
        S1_bh = NA_real_, S2_bh = mean(x$bh_any, na.rm = TRUE),
        S3_bh = if (sum(x$bh_calls, na.rm = TRUE) > 0) sum(x$bh_false, na.rm = TRUE) / sum(x$bh_calls, na.rm = TRUE) else 0,
        lost = sum(x$lost, na.rm = TRUE), err_max = suppressWarnings(max(x$err_max, na.rm = TRUE)),
        k = paste(unique(x$k), collapse = ";"), samples_removed = max(x$samples_removed, na.rm = TRUE),
        crashed = any(x$status == "crashed"), no_da = any(x$status == "no_da"), dropped = any(x$status == "dropped"),
        stringsAsFactors = FALSE)
}))
dataset$S5[!is.finite(dataset$S5)] <- NA
# F1 and F2, paired with the oracle of the same job and data set
or <- dataset[dataset$arm == "OR", c("part", "cell", "job", "scenario", "mask", "rep", "power")]
names(or)[names(or) == "power"] <- "or_power"
dataset <- merge(dataset, or, by = c("part", "cell", "job", "scenario", "mask", "rep"), all.x = TRUE)
dataset$F1 <- dataset$or_power - dataset$power
dataset$F2 <- dataset$gross_common - dataset$or_gross_common
utils::write.table(dataset, file.path(out_dir, "dataset.tsv"), sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")

# ---- per condition: means, Monte Carlo SEs, limits, classification ------------------------------------------------
estimandsOf <- function(scenario) if (startsWith(scenario, "N")) c("S1", "S2", "S5", "F2") else c("S3", "S5", "F1", "F2")
clopper <- function(x, n, level) c(lower = if (x == 0) 0 else stats::qbeta(1 - level, x, n - x + 1),
    upper = if (x == n) 1 else stats::qbeta(level, x + 1, n - x))
conditions <- list()
cond_key <- c("part", "cell", "job", "scenario", "mask", "arm")
for (x in split(dataset, dataset[cond_key], drop = TRUE)) {
    part <- x$part[[1L]]
    R <- nrow(x)
    for (e in estimandsOf(x$scenario[[1L]])) {
        v <- x[[e]]
        v <- v[is.finite(v)]
        n <- length(v)
        row <- data.frame(x[1L, cond_key], estimand = e, R = R, n = n, mean = if (n) mean(v) else NA, sd = if (n > 1L) stats::sd(v) else NA,
            stringsAsFactors = FALSE)
        if (identical(e, "S5")) {
            events <- sum(v)
            cp <- clopper(events, n, 0.975)
            row$se <- NA
            row$upper <- cp[["upper"]]
            row$lower_bonf <- clopper(events, n, 1 - 0.025 / family[[part]])[["lower"]]
        } else {
            row$se <- row$sd / sqrt(n)
            row$upper <- row$mean + stats::qt(0.975, n - 1L) * row$se
            row$lower <- row$mean - stats::qt(0.975, n - 1L) * row$se
            row$lower_bonf <- row$mean - stats::qt(1 - 0.025 / family[[part]], n - 1L) * row$se
        }
        row$no_result <- mean(x$crashed | x$no_da | x$dropped) > 0.05
        conditions[[length(conditions) + 1L]] <- row
    }
}
conditions <- do.call(rbind, lapply(conditions, function(r) { for (c in c("se", "upper", "lower", "lower_bonf")) if (is.null(r[[c]])) r[[c]] <- NA; r }))
# the oracle check (5.7): where OR does not meet S1, S2 or S3, the arms are classified on the paired excess over OR
or_cond <- conditions[conditions$arm == "OR" & conditions$estimand %in% c("S1", "S2", "S3"), ]
excess_needed <- with(or_cond, paste(part, cell, job, scenario, mask, estimand)[!(upper < EDGES[estimand])])
classify <- function(r) {
    if (isTRUE(r$no_result)) return("no result")
    if (r$estimand %in% names(EDGES)) {
        edge <- EDGES[[r$estimand]]
        if (!is.finite(r$upper)) return("not determined")
        if (r$upper < edge) return("meets")
        if (is.finite(r$lower_bonf) && r$lower_bonf >= edge) return("fails")
        return("not determined")
    }
    edges <- if (identical(r$estimand, "F1")) F1_EDGES else F2_EDGES
    words <- if (identical(r$estimand, "F1")) c("comparable", "lower", "much lower") else c("comparable", "less precise", "much less precise")
    if (!is.finite(r$mean)) return("not determined")
    words[[1L + sum(r$mean >= edges)]]
}
for (i in seq_len(nrow(conditions))) {
    r <- conditions[i, ]
    k <- paste(r$part, r$cell, r$job, r$scenario, r$mask, r$estimand)
    if (r$estimand %in% c("S1", "S2", "S3") && k %in% excess_needed && r$arm != "OR") {
        x <- dataset[dataset$part == r$part & dataset$cell == r$cell & dataset$job == r$job & dataset$scenario == r$scenario &
            dataset$mask == r$mask, ]
        a <- x[x$arm == r$arm, c("rep", r$estimand)]
        o <- x[x$arm == "OR", c("rep", r$estimand)]
        d <- merge(a, o, by = "rep")
        diff <- d[[2L]] - d[[3L]]
        diff <- diff[is.finite(diff)]
        se <- stats::sd(diff) / sqrt(length(diff))
        upper <- mean(diff) + stats::qt(0.975, length(diff) - 1L) * se
        lower <- mean(diff) - stats::qt(1 - 0.025 / family[[r$part]], length(diff) - 1L) * se
        conditions$classification[i] <- if (isTRUE(r$no_result)) "no result" else if (upper < 0.025) "meets (excess over OR)" else
            if (lower >= 0.025) "fails (excess over OR)" else "not determined (excess over OR)"
        conditions$oracle_check[i] <- TRUE
    } else {
        conditions$classification[i] <- classify(r)
        conditions$oracle_check[i] <- FALSE
    }
}
utils::write.table(conditions, file.path(out_dir, "conditions.tsv"), sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")

# ---- verdicts and outcomes (deciding jobs only) ---------------------------------------------------------------------
deciding <- conditions[conditions$job %in% "deciding" & conditions$arm %in% DECIDING, ]
deciding$group <- OMIC_GROUP[deciding$cell]
safety <- deciding[deciding$estimand %in% names(EDGES), ]
verdicts <- do.call(rbind, lapply(split(safety, safety[c("group", "arm")], drop = TRUE), function(x) {
    cls <- x$classification
    verdict <- if (any(startsWith(cls, "fails"))) "unsafe" else if (all(startsWith(cls, "meets"))) "safe" else "not shown safe"
    data.frame(group = x$group[[1L]], arm = x$arm[[1L]], verdict = verdict, meets = sum(startsWith(cls, "meets")),
        fails = sum(startsWith(cls, "fails")), not_determined = sum(startsWith(cls, "not determined")), no_result = sum(cls == "no result"),
        stringsAsFactors = FALSE)
}))
utils::write.table(verdicts, file.path(out_dir, "verdicts.tsv"), sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
outcome <- character()
for (g in unique(verdicts$group)) {
    v <- stats::setNames(verdicts$verdict[verdicts$group == g], verdicts$arm[verdicts$group == g])
    safe <- function(arm) identical(unname(v[arm]), "safe")
    eff <- deciding[deciding$group == g & deciding$estimand %in% c("F1", "F2"), ]
    # K1's efficiency clause: LA neither much lower nor much less precise than GA where GA is safe (paired via OR)
    la_worse <- FALSE
    if (safe("LA") && safe("GA")) {
        la <- eff[eff$arm == "LA", ]
        ga <- eff[eff$arm == "GA", ]
        m <- merge(la, ga, by = c("cell", "scenario", "mask", "estimand"), suffixes = c("_la", "_ga"))
        edge <- ifelse(m$estimand == "F1", F1_EDGES[[2L]], F2_EDGES[[2L]])
        la_worse <- any(m$mean_la - m$mean_ga >= edge, na.rm = TRUE)
    }
    k <- if (safe("LA") && !la_worse) "K1" else if (safe("GA")) "K2" else if (safe("NO")) "K3" else "K4"
    outcome <- c(outcome, sprintf("%s: %s (LA %s, GA %s, NO %s%s)", g, k, v["LA"], v["GA"], v["NO"],
        if (la_worse) "; K1's efficiency clause failed" else ""))
}
writeLines(outcome, file.path(out_dir, "outcome.txt"))

# ---- statuses and Part M --------------------------------------------------------------------------------------------
st <- unique(rows[c("part", "cell", "job", "scenario", "mask", "rep", "arm", "arm_scale", "status", "class")])
st$class[is.na(st$class)] <- ""
utils::write.table(stats::aggregate(list(n = rep(1L, nrow(st))), st[c("part", "cell", "job", "arm", "arm_scale", "status", "class")], length),
    file.path(out_dir, "statuses.tsv"), sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
partm <- function(columns, filter, name) {
    x <- rows[filter, c("cell", "job", "scenario", "mask", "rep", "arm", "arm_scale", "assay", "contrast", columns), drop = FALSE]
    utils::write.table(x, file.path(out_dir, name), sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
}
partm(c("k", "mdf_null05", "mdf_null05_reduced", "mdf_df_le_k", "mdf_tested"), !is.na(num(rows$mdf_null05)), "partm-df.tsv")
partm(c("k", "n_controls", "mamp_sigma_min", "mamp_e_max", "mamp_planted_share", "s5", "null_p05"), !is.na(num(rows$mamp_e_max)), "partm-amp.tsv")
partm(c("k", "n_controls", "mamp_cf_n", "mamp_cf_singular", "mamp_cf_log10kappa", "s5", "err_max"), !is.na(num(rows$mamp_cf_n)),
    "partm-closed-form.tsv")
partm(c("mu_itsd_dev", "mu_log2_dev", "mu_loess_dev", "mu_ctl_same", "mu_k_same", "mu_pct_same", "mu_ruv_dev", "mu_logfc_dev", "mu_t_dev",
    "mu_call_changes", "mu_t_dev_notrend"), rows$job == "munit", "partm-unit.tsv")
repair <- conditions[conditions$arm %in% c("LA", "LP") & conditions$job %in% c("deciding", "lp") & conditions$scenario %in% c("N0", "N1", "E0", "E1"), ]
utils::write.table(repair, file.path(out_dir, "partm-repair.tsv"), sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
cat(sprintf("KIT-188 scoring: %d rows, %d data-set rows, %d conditions\n", nrow(rows), nrow(dataset), nrow(conditions)))
