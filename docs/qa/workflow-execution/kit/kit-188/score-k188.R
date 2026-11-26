# KIT-188: the scoring (sections 5.5 to 5.7, 18.5; r1.19 to r1.21: 38.2, 38.3, 40.3), run once on a runner over the committed rows
# of S. Evidence code under docs/qa: no package code. Hashed in the declaration; a defect found after scoring is fixed in a dated
# section and all rows re-scored once (section 10, Attempts).
#
# Usage: Rscript score-k188.R <rows dir> <out dir> [<part L family> <part P family>]
#   rows dir  the assembled rows (*.tsv.gz) of every S job (deciding, reported, munit, ba, lp, las)
#   family    the Bonferroni family sizes (default 720 for Part L under D6 (a), 936 for Part P)
# Environment (read once, written to settings.txt; both are set ONLY by E16, the scorer on a constructed row set, 7.4; a counted run never
# sets them):
#   K188_SCORE_R      "<deciding>,<other>": the declared data-set ranges (400,100 as 5.6 and 38.3 declare)
#   K188_SCORE_CELLS  "<cell>:<job>,...": the declared (cell, job) set of Part L's S (38.3 (ii), SA-P14)
# phi (38.2 (6), VAL-335, r1.22): events are data sets with an assay not corrected because the method could not and no gate status
# (gate first); phi, the deciding figure, divides by the data sets whose journey reached the RUV step; phi_all (all declared data
# sets) is an informational column, not deciding.
# Writes: conditions.tsv (per part, arm, cell, scenario, mask, estimand: R, mean, SD, Monte Carlo SE, limits, classification,
#   r_used, r_out, n_gate, class_best, class_worst, sensitive, n_paired_dropped, phi columns, lan0_share, label), verdicts.tsv,
#   outcome.txt, dataset.tsv (per data set and arm), statuses.tsv, statuses-assay.tsv, phi.tsv, capbind.tsv, partm-df.tsv,
#   partm-amp.tsv, partm-closed-form.tsv, partm-unit.tsv, partm-repair.tsv, partm-levels.tsv, partm-warn.tsv, settings.txt.

# Partial-match guard (smoke check send 3, the coordinator's ruling): `$` on a list or data frame matches a name by its prefix when none is exact;
# that warning is an error here, in every process of every run (lifecycle, L1r, L1c and counted runs). A stop from it is an instrument defect.
options(warnPartialMatchDollar = TRUE)
globalCallingHandlers(warning = function(w) {
    if (grepl("partial match of", conditionMessage(w), fixed = TRUE)) stop(paste("partial-match guard:", conditionMessage(w)), call. = FALSE)
})
args <- commandArgs(TRUE)
stopifnot(length(args) >= 2L)
rows_dir <- args[[1L]]
out_dir <- args[[2L]]
family <- c(L = if (length(args) >= 3L) as.numeric(args[[3L]]) else 720, P = if (length(args) >= 4L) as.numeric(args[[4L]]) else 936)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

EDGES <- c(S1 = 0.075, S2 = 0.075, S3 = 0.075, S5 = 0.01)
F1_EDGES <- c(0.05, 0.20)
F2_EDGES <- c(0.005, 0.05)
PHI_EDGE <- 0.05 # 38.2 (7): K1's clause on the non-correction share's one-sided 97.5 % Clopper-Pearson upper limit
GATE_EDGE <- 0.05 # 38.2 (2): more than 5 % of a condition's declared data sets out
DECIDING <- c("LA", "GA", "NO")
# L24 (38.5) is reported: it has no entry here, so it enters no family, verdict or outcome
OMIC_GROUP <- c(L12 = "lipid_metab", L64 = "lipid_metab", M12 = "lipid_metab", M168 = "lipid_metab", P0 = "proteomics",
    P1 = "proteomics", P16 = "proteomics", P2 = "proteomics", P3 = "proteomics", P4b = "proteomics")
LEGACY_ARMS <- c("LA", "LP", "LAs", "LM") # 38.2 (6): the arms whose non-correction is "no k found"
DECLARED_ARMS <- c("LA", "GA", "NO") # a data set's arms per job type, as the sim's ARMS_OF declares them (OR added by the unit expansion)
JOB_ARMS <- list(L = list(deciding = c("LA", "GA", "NO", "OR"), reported = c("LA", "LM", "FP", "SC", "LCm", "OR"), ba = c("LA", "GA", "NO", "OR"),
    lp = c("LP", "OR"), las = c("LAs", "OR")), P = list(deciding = c("LA", "GA", "NO", "OR"), reported = c("LM", "OR"), ba = c("LA", "GA", "NO", "OR")))
SCORED_JOBS <- c("deciding", "ba", "reported", "lp", "las")
SCENARIOS <- c("N0", "N1", "E0", "E1", "N2", "E2")
MASKS <- c("none", "vd5", "vd20")
settings <- c(score_r = Sys.getenv("K188_SCORE_R", "400,100"), score_cells = Sys.getenv("K188_SCORE_CELLS", "declared (5.6, 38.3)"))
ranges <- as.integer(strsplit(Sys.getenv("K188_SCORE_R", "400,100"), ",", fixed = TRUE)[[1L]])
declaredRange <- function(job, cell) seq_len(if (identical(job, "deciding") && !identical(cell, "L24")) ranges[[1L]] else ranges[[2L]])
writeLines(sprintf("%s = %s", names(settings), settings), file.path(out_dir, "settings.txt"))

files <- list.files(rows_dir, pattern = "\\.tsv\\.gz$", full.names = TRUE, recursive = TRUE)
# every column is read as text (a column of the letter T, or of TRUE and FALSE, must not become logical: the smoke check's defect 1);
# the numeric columns are converted through num() below
rows <- do.call(rbind, lapply(files, function(f) utils::read.delim(gzfile(f), stringsAsFactors = FALSE, na.strings = "NA", colClasses = "character")))
# the scorer refuses to produce an outcome from any input that carries the smoke flag (an instrument check, K188_SMOKE=1)
if ("smoke" %in% names(rows) && any(rows$smoke %in% "1")) stop("KIT-188 scoring refuses: an input carries the smoke flag (K188_SMOKE=1, an instrument check); no outcome is produced from it")
num <- function(x) suppressWarnings(as.numeric(x))
# rows written before r1.16 carry none of the per-assay columns
for (column in c("itsd_ok", "norm_ok", "ruv_status", "da_status", "corr_removed_one_assay", "ruv_k_found", "ruv_k_sweep", "ruv_k_used", "ruv_k_cap",
    "n_controls", "complete_controls", "k_estimate", "logfc_max")) if (!column %in% names(rows)) rows[[column]] <- NA
for (column in c("rep", "n_T", "n_tested", "lost", "n_null", "null_p05", "null_p01", "calls", "false_calls", "true_calls", "any_call",
    "fdp", "power", "n_planted", "bh_calls", "bh_false", "bh_any", "bh_fdp", "bh_power", "gross1", "n_common", "s5", "or_s5",
    "gross1_common", "or_gross1_common", "err_max", "k", "samples_removed", "n_controls", "complete_controls", "k_estimate", "logfc_max",
    "ruv_k_found", "ruv_k_sweep", "ruv_k_used", "ruv_k_cap")) rows[[column]] <- num(rows[[column]])
rows$arm_scale[is.na(rows$arm_scale)] <- "1"
isFalse <- function(x) as.character(x) %in% "FALSE"

# ---- (i) a counted run holds no instrument failure (38.3 (i); 35.2, 7.4: the columns, not only the status) ----------------
instrument <- isFalse(rows$itsd_ok) | isFalse(rows$norm_ok) | rows$status %in% "instrument"
if (any(instrument)) stop(sprintf("KIT-188 scoring stops: %d row(s) have itsd_ok or norm_ok FALSE or status instrument (an instrument defect, VAL-54, or the bound commit's own behaviour, reported before scoring: 7.4); first: %s",
    sum(instrument), paste(utils::head(unique(paste(rows$cell, rows$job, rows$rep, rows$arm, rows$class)[instrument]), 3L), collapse = " | ")))

# ---- (ii) unit rows expanded to every arm and condition of their job, as `crashed`; the declared range checked (38.3 (ii)) ----
unit_rows <- rows[rows$arm %in% "unit" & rows$job %in% SCORED_JOBS, , drop = FALSE]
rows <- rows[!rows$arm %in% "unit", , drop = FALSE]
condsOf <- function(part, cell, job, u_scenario = NA, u_mask = NA) {
    if (identical(part, "P")) return(data.frame(scenario = u_scenario, mask = u_mask, stringsAsFactors = FALSE))
    scenarios <- if (identical(job, "deciding") && cell %in% c("L12", "M12")) SCENARIOS else SCENARIOS[1:4]
    expand.grid(scenario = scenarios, mask = if (identical(job, "ba")) "ba" else MASKS, stringsAsFactors = FALSE)
}
expanded <- do.call(rbind, c(list(NULL), lapply(seq_len(nrow(unit_rows)), function(i) {
    u <- unit_rows[i, ]
    arms <- JOB_ARMS[[u$part]][[u$job]]
    cond <- condsOf(u$part, u$cell, u$job, u$scenario, u$mask)
    grid <- merge(cond, data.frame(arm = arms, stringsAsFactors = FALSE))
    data.frame(part = u$part, cell = u$cell, job = u$job, scenario = grid$scenario, mask = grid$mask, rep = u$rep, arm = grid$arm,
        arm_scale = "1", scoreset = "T", status = "crashed", class = paste("unit:", u$class), stringsAsFactors = FALSE)
})))
present <- unique(c(paste(rows$cell, rows$job, rows$scenario, rows$mask, rows$arm, rows$rep)[rows$job %in% SCORED_JOBS & rows$arm_scale == "1"],
    if (!is.null(expanded)) paste(expanded$cell, expanded$job, expanded$scenario, expanded$mask, expanded$arm, expanded$rep)))
gaps <- character()
# SA-P14, 38.3 (ii): the declared (cell, job) set of Part L's S, not what is present; a whole absent (cell, job) is a gap
DECLARED_CELL_JOBS <- rbind(data.frame(cell = c("L12", "L64", "M12", "M168", "L24"), job = "deciding"),
    data.frame(cell = rep(c("L12", "L64", "M12", "M168"), 3L), job = rep(c("reported", "lp", "las"), each = 4L)),
    data.frame(cell = c("L64", "M168"), job = "ba"))
if (nzchar(Sys.getenv("K188_SCORE_CELLS"))) {
    spec <- strsplit(strsplit(Sys.getenv("K188_SCORE_CELLS"), ",", fixed = TRUE)[[1L]], ":", fixed = TRUE)
    DECLARED_CELL_JOBS <- data.frame(cell = vapply(spec, `[[`, "", 1L), job = vapply(spec, `[[`, "", 2L), stringsAsFactors = FALSE)
}
cell_jobs <- if (any(rows$part == "L") || any(unit_rows$part == "L")) DECLARED_CELL_JOBS else DECLARED_CELL_JOBS[0, ]
for (i in seq_len(nrow(cell_jobs))) {
    cell <- cell_jobs$cell[[i]]
    job <- cell_jobs$job[[i]]
    cond <- condsOf("L", cell, job)
    # the oracle has no row where no arm reached the QC step; its absence is a paired-set drop (38.2 (5)), not a timeout gap
    arms <- setdiff(JOB_ARMS$L[[job]], "OR")
    declared <- do.call(paste, c(expand.grid(cell, job, paste(cond$scenario, cond$mask), arms, declaredRange(job, cell), stringsAsFactors = FALSE), sep = " "))
    missing <- setdiff(declared, present)
    if (length(missing)) gaps <- c(gaps, sprintf("%s %s: %d (job, cell, scenario, mask, arm, data set) absent without a unit row, first %s", cell, job, length(missing), missing[[1L]]))
}
if (length(gaps)) stop(paste("KIT-188 scoring stops: a data set of the declared range is absent without a unit row (a timeout gap: an instrument defect, VAL-54; its units are re-sent before scoring):",
    paste(gaps, collapse = "; ")))
if (!is.null(expanded)) {
    for (column in setdiff(names(rows), names(expanded))) expanded[[column]] <- NA
    rows <- rbind(rows, expanded[names(rows)])
}
main <- rows[rows$job %in% SCORED_JOBS & rows$arm_scale == "1" & rows$scoreset %in% "T", , drop = FALSE]

# ---- (iii) per data set and arm (the unit, VAL-11): the gate status, whether a DA table exists, the non-correction event -------
pooled <- function(num, den) if (sum(den, na.rm = TRUE) > 0) sum(num, na.rm = TRUE) / sum(den, na.rm = TRUE) else NA_real_
key <- c("part", "cell", "job", "scenario", "mask", "rep", "arm")
groups <- split(main, main[key], drop = TRUE)
dataset <- do.call(rbind, lapply(groups, function(x) {
    statuses <- unique(x$status)
    arm <- x$arm[[1L]]
    # a family has a DA table unless its journey crashed, was refused, had none, did not run (LCm) or is a unit crash
    has_da <- any(!x$status %in% c("crashed", "refused", "no_da", "not_run"))
    gate <- any(x$status %in% c("crashed", "refused", "no_da", "dropped"))
    ruv <- x$ruv_status[!is.na(x$ruv_status)]
    reached <- length(ruv) > 0L && any(ruv != "not_reached")
    event <- if (arm %in% LEGACY_ARMS) any(x$ruv_status %in% c("uncorrected", "dropped") & is.na(x$ruv_k_found)) else
        if (arm %in% c("GA", "FP", "SC", "LCm")) any(x$ruv_status %in% "refused" & !x$status %in% c("crashed", "refused")) else NA
    data.frame(x[1L, key], status = paste(sort(statuses), collapse = "+"), has_da = has_da, gate = gate, ruv_reached = reached, phi_event = event,
        all_uncorrected = arm %in% LEGACY_ARMS && length(ruv) > 0L && all(ruv == "uncorrected"),
        not_corrected = any(x$ruv_status %in% c("uncorrected", "dropped", "none", "refused")), da_failed = any(x$da_status %in% "failed"),
        cap_bound = arm %in% LEGACY_ARMS && any(x$ruv_k_sweep > x$ruv_k_used, na.rm = TRUE),
        S1 = pooled(x$null_p05 * x$n_null, x$n_null),
        S1_01 = pooled(x$null_p01 * x$n_null, x$n_null),
        S2 = mean(x$any_call, na.rm = TRUE), S2_pooled = as.numeric(any(x$any_call > 0, na.rm = TRUE)),
        S3 = if (sum(x$calls, na.rm = TRUE) > 0) sum(x$false_calls, na.rm = TRUE) / sum(x$calls, na.rm = TRUE) else 0,
        S5 = if (identical(arm, "OR")) suppressWarnings(max(x$or_s5, na.rm = TRUE)) else suppressWarnings(max(x$s5, na.rm = TRUE)),
        power = pooled(x$true_calls, x$n_planted), calls = sum(x$calls, na.rm = TRUE), false_calls = sum(x$false_calls, na.rm = TRUE),
        gross_common = pooled(x$gross1_common * x$n_common, x$n_common),
        or_gross_common = pooled(x$or_gross1_common * x$n_common, x$n_common),
        S1_bh = NA_real_, S2_bh = mean(x$bh_any, na.rm = TRUE),
        S3_bh = if (sum(x$bh_calls, na.rm = TRUE) > 0) sum(x$bh_false, na.rm = TRUE) / sum(x$bh_calls, na.rm = TRUE) else 0,
        lost = sum(x$lost, na.rm = TRUE), err_max = suppressWarnings(max(x$err_max, na.rm = TRUE)), logfc_max = suppressWarnings(max(x$logfc_max, na.rm = TRUE)),
        k = paste(unique(x$k), collapse = ";"), k_sweep = paste(unique(x$ruv_k_sweep), collapse = ";"), k_cap = paste(unique(x$ruv_k_cap), collapse = ";"),
        k_num = suppressWarnings(mean(x$k, na.rm = TRUE)), n_controls_mean = suppressWarnings(mean(x$n_controls, na.rm = TRUE)),
        complete_min = suppressWarnings(min(x$complete_controls, na.rm = TRUE)), n_controls_min = suppressWarnings(min(x$n_controls, na.rm = TRUE)),
        k_est_min = suppressWarnings(min(x$k_estimate, na.rm = TRUE)),
        samples_removed = suppressWarnings(max(x$samples_removed, na.rm = TRUE)), corr_one_assay = paste(unique(stats::na.omit(x$corr_removed_one_assay[nzchar(x$corr_removed_one_assay)])), collapse = ";"),
        crashed = any(x$status == "crashed"), no_da = any(x$status == "no_da"), dropped = any(x$status == "dropped"), refused = any(x$status == "refused"),
        stringsAsFactors = FALSE)
}))
for (column in c("S5", "err_max", "logfc_max", "k_num", "n_controls_mean", "complete_min", "n_controls_min", "k_est_min", "samples_removed")) {
    dataset[[column]][!is.finite(dataset[[column]])] <- NA
}
# (vii) F1 and F2, paired with the oracle of the same job and data set, on the data sets where both returned a DA table (38.2 (5))
or <- dataset[dataset$arm == "OR", c("part", "cell", "job", "scenario", "mask", "rep", "power", "has_da")]
names(or)[names(or) == "power"] <- "or_power"
names(or)[names(or) == "has_da"] <- "or_has_da"
dataset <- merge(dataset, or, by = c("part", "cell", "job", "scenario", "mask", "rep"), all.x = TRUE)
paired <- dataset$has_da & dataset$or_has_da %in% TRUE
dataset$F1 <- ifelse(paired, dataset$or_power - dataset$power, NA_real_)
dataset$F2 <- ifelse(paired, dataset$gross_common - dataset$or_gross_common, NA_real_)
utils::write.table(dataset, file.path(out_dir, "dataset.tsv"), sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")

# ---- (iv) per condition: means, Monte Carlo SEs, limits, classification --------------------------------------------------
# The estimands are taken over the data sets that returned a DA table (38.2 (3)); R is the declared count, r_used those with a
# table, r_out the rest, n_gate the data sets that count toward "no result" (38.2 (2)).
estimandsOf <- function(scenario) if (startsWith(scenario, "N")) c("S1", "S2", "S5", "F2") else c("S3", "S5", "F1", "F2")
clopper <- function(x, n, level) c(lower = if (x == 0) 0 else stats::qbeta(1 - level, x, n - x + 1),
    upper = if (x == n) 1 else stats::qbeta(level, x + 1, n - x))
# a classification of safety values (S1, S2, S3: means; S5: events) as 5.7 defines it; used for the best and worst case (38.2 (4))
classifyValues <- function(v, e, part, no_result) {
    if (isTRUE(no_result)) return("no result")
    n <- length(v)
    if (n < 2L) return("not determined")
    if (identical(e, "S5")) {
        upper <- clopper(sum(v), n, 0.975)[["upper"]]
        lower_bonf <- clopper(sum(v), n, 1 - 0.025 / family[[part]])[["lower"]]
    } else {
        se <- stats::sd(v) / sqrt(n)
        upper <- mean(v) + stats::qt(0.975, n - 1L) * se
        lower_bonf <- mean(v) - stats::qt(1 - 0.025 / family[[part]], n - 1L) * se
    }
    if (upper < EDGES[[e]]) "meets" else if (lower_bonf >= EDGES[[e]]) "fails" else "not determined"
}
conditions <- list()
cond_key <- c("part", "cell", "job", "scenario", "mask", "arm")
for (x_all in split(dataset, dataset[cond_key], drop = TRUE)) {
    part <- x_all$part[[1L]]
    x <- x_all[x_all$has_da, , drop = FALSE]
    R <- nrow(x_all)
    n_gate <- sum(x_all$gate)
    no_result <- n_gate / R > GATE_EDGE
    for (e in estimandsOf(x_all$scenario[[1L]])) {
        v <- x[[e]]
        v <- v[is.finite(v)]
        n <- length(v)
        row <- data.frame(x_all[1L, cond_key], estimand = e, R = R, r_used = nrow(x), r_out = R - nrow(x), n_gate = n_gate, n = n,
            mean = if (n) mean(v) else NA, sd = if (n > 1L) stats::sd(v) else NA, se = NA_real_, upper = NA_real_, lower = NA_real_,
            lower_bonf = NA_real_, class_best = NA_character_, class_worst = NA_character_, n_paired_dropped = NA_integer_,
            stringsAsFactors = FALSE)
        if (identical(e, "S5")) {
            events <- sum(v)
            cp <- clopper(events, n, 0.975)
            row$upper <- cp[["upper"]]
            row$lower_bonf <- clopper(events, n, 1 - 0.025 / family[[part]])[["lower"]]
        } else {
            row$se <- row$sd / sqrt(n)
            row$upper <- row$mean + stats::qt(0.975, n - 1L) * row$se
            row$lower <- row$mean - stats::qt(0.975, n - 1L) * row$se
            row$lower_bonf <- row$mean - stats::qt(1 - 0.025 / family[[part]], n - 1L) * row$se
        }
        # 38.2 (4): the classification again with the data sets that returned no table at best case (every safety value 0) and at
        # worst case (every safety value 1); reported, deciding nothing
        if (e %in% c("S1", "S2", "S3", "S5") && row$r_out > 0L) {
            row$class_best <- classifyValues(c(v, rep(0, row$r_out)), e, part, no_result)
            row$class_worst <- classifyValues(c(v, rep(1, row$r_out)), e, part, no_result)
        }
        if (e %in% c("F1", "F2")) row$n_paired_dropped <- R - n
        row$no_result <- no_result
        # 38.3 (r1.21, R20-P13): the smallest upper limit S5 can reach at this R (zero events); S5 cannot meet where it is not below the edge
        row$s5_best_upper <- if (identical(e, "S5") && n > 0L) clopper(0, n, 0.975)[["upper"]] else NA_real_
        conditions[[length(conditions) + 1L]] <- row
    }
}
conditions <- do.call(rbind, conditions)
conditions$sensitive <- !is.na(conditions$class_best) & conditions$class_best != conditions$class_worst
conditions$role <- ifelse(conditions$job %in% "deciding" & conditions$cell %in% names(OMIC_GROUP) & conditions$arm %in% c(DECIDING, "OR"), "deciding", "reported")
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
# the paired excess over the oracle (5.7's oracle check): meets below 0.025 by the one-sided 97.5 % limit, fails at or above it by the
# Bonferroni lower limit; with fewer than two paired data sets it is not determined (SA-P04)
excessClass <- function(diff, part, no_result) {
    if (isTRUE(no_result)) return("no result")
    diff <- diff[is.finite(diff)]
    if (length(diff) < 2L) return("not determined (excess over OR)")
    se <- stats::sd(diff) / sqrt(length(diff))
    upper <- mean(diff) + stats::qt(0.975, length(diff) - 1L) * se
    lower <- mean(diff) - stats::qt(1 - 0.025 / family[[part]], length(diff) - 1L) * se
    if (upper < 0.025) "meets (excess over OR)" else if (lower >= 0.025) "fails (excess over OR)" else "not determined (excess over OR)"
}
for (i in seq_len(nrow(conditions))) {
    r <- conditions[i, ]
    k <- paste(r$part, r$cell, r$job, r$scenario, r$mask, r$estimand)
    if (r$estimand %in% c("S1", "S2", "S3") && k %in% excess_needed && r$arm != "OR") {
        x <- dataset[dataset$part == r$part & dataset$cell == r$cell & dataset$job == r$job & dataset$scenario == r$scenario &
            dataset$mask == r$mask, ]
        # 38.2 (5): the paired excess uses the data sets where the arm and the oracle both returned a DA table
        a <- x[x$arm == r$arm & x$has_da, c("rep", r$estimand)]
        o <- x[x$arm == "OR" & x$has_da, c("rep", r$estimand)]
        d <- merge(a, o, by = "rep")
        diff <- d[[2L]] - d[[3L]]
        diff <- diff[is.finite(diff)]
        conditions$classification[i] <- excessClass(diff, r$part, r$no_result)
        conditions$oracle_check[i] <- TRUE
        conditions$n_paired_dropped[i] <- r$R - length(diff)
        # 38.2 (4), r1.23 (SA-P13): the arm's data sets out under 38.2 (3) put at 0 (best) and 1 (worst), paired with the oracle's value
        # on that data set where the oracle has a DA table
        out_o <- merge(data.frame(rep = x$rep[x$arm == r$arm & !x$has_da]), o, by = "rep")[[2L]]
        out_o <- out_o[is.finite(out_o)]
        if (r$r_out > 0L) {
            conditions$class_best[i] <- excessClass(c(diff, 0 - out_o), r$part, r$no_result)
            conditions$class_worst[i] <- excessClass(c(diff, 1 - out_o), r$part, r$no_result)
        } else conditions$class_best[i] <- conditions$class_worst[i] <- NA_character_
        conditions$sensitive[i] <- !is.na(conditions$class_best[i]) && conditions$class_best[i] != conditions$class_worst[i]
    } else {
        conditions$classification[i] <- classify(r)
        conditions$oracle_check[i] <- FALSE
    }
}
# 38.3 (r1.21, R20-P13): at the R of the reported cells and arms S5 cannot meet; the classification is descriptive
conditions$label <- ifelse(conditions$estimand == "S5" & is.finite(conditions$s5_best_upper) & conditions$s5_best_upper >= EDGES[["S5"]],
    sprintf("cannot meet S5 at this R (zero events give an upper limit of %.4f)", conditions$s5_best_upper),
    ifelse(conditions$role == "reported", "reported (descriptive)", ""))

# ---- (v) phi, the non-correction share, with its one-sided 97.5 % Clopper-Pearson upper limit, both denominators (38.2 (6)) ----
# event: an assay not corrected because the method could not (legacy: uncorrected or dropped with no k found; declared designs:
# the factors refused, not an estimated count of 0); a stop at the RUV step is not an event and counts in the gate. phi_all
# divides by the R declared data sets (r1.19), phi by those whose journey reached the RUV step (the draft rule, 40.1).
phi_rows <- list()
for (x in split(dataset[!dataset$arm %in% c("OR", "NO"), ], dataset[!dataset$arm %in% c("OR", "NO"), cond_key], drop = TRUE)) {
    events <- sum(x$phi_event %in% TRUE & !x$gate) # gate first (38.2 (6), VAL-335): a data set enters at most one of the two
    declared <- nrow(x)
    reached <- sum(x$ruv_reached)
    phi_rows[[length(phi_rows) + 1L]] <- data.frame(x[1L, cond_key], events = events, r_declared = declared, r_reached = reached,
        events_also_gate = sum(x$phi_event %in% TRUE & x$gate), phi_all = if (declared) events / declared else NA_real_,
        phi_all_upper = if (declared) clopper(events, declared, 0.975)[["upper"]] else NA_real_,
        phi = if (reached) events / reached else NA_real_, phi_upper = if (reached) clopper(events, reached, 0.975)[["upper"]] else NA_real_,
        share_not_corrected = mean(x$not_corrected), share_da_failed = mean(x$da_failed), lan0_share = mean(x$all_uncorrected),
        stringsAsFactors = FALSE)
}
phi_tab <- do.call(rbind, phi_rows)
phi_tab$lan0_mark <- phi_tab$arm %in% LEGACY_ARMS & phi_tab$lan0_share > GATE_EDGE # 38.2 (9)
utils::write.table(phi_tab, file.path(out_dir, "phi.tsv"), sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
conditions <- merge(conditions, phi_tab, by = cond_key, all.x = TRUE, sort = FALSE)
utils::write.table(conditions, file.path(out_dir, "conditions.tsv"), sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
# 5.5 (r1.17, r1.21): the share of LA journeys whose sweep k exceeded the cap, per condition (reported, deciding nothing)
capbind <- do.call(rbind, lapply(split(dataset[dataset$arm %in% LEGACY_ARMS & dataset$ruv_reached, ], dataset[dataset$arm %in% LEGACY_ARMS & dataset$ruv_reached, cond_key], drop = TRUE),
    function(x) data.frame(x[1L, cond_key], journeys = nrow(x), cap_bound = sum(x$cap_bound), share = mean(x$cap_bound), stringsAsFactors = FALSE)))
if (!is.null(capbind)) utils::write.table(capbind, file.path(out_dir, "capbind.tsv"), sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")

# ---- (viii) verdicts and outcomes (deciding jobs only; L24, LP, LAs and the reported arms enter none) -----------------------
deciding <- conditions[conditions$role %in% "deciding" & conditions$arm %in% DECIDING, ]
deciding$group <- OMIC_GROUP[deciding$cell]
safety <- deciding[deciding$estimand %in% names(EDGES), ]
verdicts <- do.call(rbind, lapply(split(safety, safety[c("group", "arm")], drop = TRUE), function(x) {
    cls <- x$classification
    verdict <- if (any(startsWith(cls, "fails"))) "unsafe" else if (all(startsWith(cls, "meets"))) "safe" else "not shown safe"
    data.frame(group = x$group[[1L]], arm = x$arm[[1L]], verdict = verdict, meets = sum(startsWith(cls, "meets")),
        fails = sum(startsWith(cls, "fails")), not_determined = sum(startsWith(cls, "not determined")), no_result = sum(cls == "no result"),
        sensitive = sum(x$sensitive, na.rm = TRUE), stringsAsFactors = FALSE)
}))
utils::write.table(verdicts, file.path(out_dir, "verdicts.tsv"), sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
phiPass <- function(g, arm, column) {
    p <- phi_tab[phi_tab$job == "deciding" & phi_tab$arm == arm & OMIC_GROUP[phi_tab$cell] %in% g, ]
    nrow(p) > 0L && all(is.finite(p[[column]]) & p[[column]] < PHI_EDGE)
}
phiMax <- function(g, arm, column) {
    p <- phi_tab[phi_tab$job == "deciding" & phi_tab$arm == arm & OMIC_GROUP[phi_tab$cell] %in% g, ]
    if (nrow(p)) max(p[[column]], na.rm = TRUE) else NA_real_
}
outcome <- character()
denominators <- "phi_upper" # VAL-335: phi over the data sets that reached the RUV step is THE figure
for (g in unique(verdicts$group)) {
    v <- stats::setNames(verdicts$verdict[verdicts$group == g], verdicts$arm[verdicts$group == g])
    safe <- function(arm) identical(unname(v[arm]), "safe")
    eff <- deciding[deciding$group == g & deciding$estimand %in% c("F1", "F2"), ]
    # K1's efficiency clause (R20-P12): in any condition, when GA's verdict is safe, LA is neither much lower nor much less precise
    # than GA; compared per data set on the data sets where LA, GA and the oracle all returned a DA table (38.2 (5))
    la_worse <- FALSE
    if (safe("LA") && safe("GA")) {
        for (cell in unique(eff$cell)) for (scenario in unique(eff$scenario[eff$cell == cell])) for (mask in unique(eff$mask[eff$cell == cell])) {
            d <- dataset[dataset$job == "deciding" & dataset$cell == cell & dataset$scenario == scenario & dataset$mask == mask, ]
            common <- Reduce(intersect, list(d$rep[d$arm == "LA" & d$has_da], d$rep[d$arm == "GA" & d$has_da], d$rep[d$arm == "OR" & d$has_da]))
            for (e in intersect(estimandsOf(scenario), c("F1", "F2"))) {
                la <- d[d$arm == "LA" & d$rep %in% common, e]
                ga <- d[d$arm == "GA" & d$rep %in% common, e]
                edge <- if (e == "F1") F1_EDGES[[2L]] else F2_EDGES[[2L]]
                if (length(common) && mean(la, na.rm = TRUE) - mean(ga, na.rm = TRUE) >= edge) la_worse <- TRUE
            }
        }
    }
    for (column in denominators) {
        phi_pass <- phiPass(g, "LA", column)
        k <- if (safe("LA") && !la_worse && phi_pass) "K1" else
            if (safe("LA") && !la_worse && !phi_pass && !safe("GA")) "no outcome by rule: K1 blocked by non-correction" else # VAL-248, option (1)
            if (safe("GA")) "K2" else if (safe("NO")) "K3" else "K4"
        outcome <- c(outcome, sprintf("%s: %s (LA %s, GA %s, NO %s%s%s)%s", g, k, v["LA"], v["GA"], v["NO"],
            if (la_worse) "; K1's efficiency clause failed" else "", if (!phi_pass) "; K1's non-correction clause failed" else "",
            sprintf(" [phi over the data sets that reached the RUV step%s; LA's largest upper limit %.3f%s]", "", phiMax(g, "LA", column),
                if (identical(k, "K2")) sprintf("; GA's largest upper limit %.3f, per condition in phi.tsv", phiMax(g, "GA", column)) else "")))
    }
}
writeLines(outcome, file.path(out_dir, "outcome.txt"))

# ---- statuses and Part M --------------------------------------------------------------------------------------------
st <- unique(rows[c("part", "cell", "job", "scenario", "mask", "rep", "arm", "arm_scale", "status", "class")])
st$class[is.na(st$class)] <- ""
utils::write.table(stats::aggregate(list(n = rep(1L, nrow(st))), st[c("part", "cell", "job", "arm", "arm_scale", "status", "class")], length),
    file.path(out_dir, "statuses.tsv"), sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
# 35.2 (`score-k188.R`): the per-assay RUV and DA statuses as counts per arm, condition and ruv_status (a second table, so that
# statuses.tsv keeps its meaning of one count per journey)
sa <- unique(rows[rows$job %in% SCORED_JOBS & rows$arm_scale == "1", c("part", "cell", "job", "scenario", "mask", "rep", "arm", "assay", "ruv_status", "da_status")])
sa$ruv_status[is.na(sa$ruv_status)] <- ""
sa$da_status[is.na(sa$da_status)] <- ""
utils::write.table(stats::aggregate(list(n = rep(1L, nrow(sa))), sa[c("part", "cell", "job", "scenario", "mask", "arm", "ruv_status", "da_status")], length),
    file.path(out_dir, "statuses-assay.tsv"), sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
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
# 18.4a (VAL-126): proteomics ties and row order per legacy journey, with their distribution (reported, no threshold)
ties <- rows[rows$part == "P" & rows$arm %in% c("LA", "LP", "LM") & !is.na(num(rows$tie_share)) & rows$scoreset %in% "T",
    c("cell", "job", "scenario", "mask", "rep", "arm", "tie_share", "tie_block", "tie_block_onegroup", "perm_jaccard", "tie_rule_match")]
ties <- unique(ties)
utils::write.table(ties, file.path(out_dir, "partp-ties.tsv"), sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
if (nrow(ties)) {
    tie_summary <- do.call(rbind, lapply(split(ties, ties[c("cell", "arm")], drop = TRUE), function(x) data.frame(cell = x$cell[[1L]],
        arm = x$arm[[1L]], journeys = nrow(x), share_median = stats::median(num(x$tie_share)), share_p90 = stats::quantile(num(x$tie_share), 0.9, names = FALSE),
        share_above_0 = mean(num(x$tie_share) > 0), jaccard_median = stats::median(num(x$perm_jaccard)),
        jaccard_p10 = stats::quantile(num(x$perm_jaccard), 0.1, names = FALSE), jaccard_below_1 = mean(num(x$perm_jaccard) < 1),
        replica_matches = mean(as.logical(x$tie_rule_match)), stringsAsFactors = FALSE)))
    utils::write.table(tie_summary, file.path(out_dir, "partp-ties-summary.tsv"), sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
}
# 5.8b: the named case K172c (reported)
k172c <- rows[rows$job %in% "k172c", c("cell", "arm", "assay", "contrast", "status", "class", "n_T", "n_tested", "lost", "calls",
    "false_calls", "null_p05", "any_call", "k", "ruv_k_sweep", "ruv_k_used", "ruv_k_cap", "pct", "n_controls", "complete_controls", "err_max")]
utils::write.table(k172c, file.path(out_dir, "partl-k172c.tsv"), sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")

# M-repair (5.8) and the LAs-LA table (34): two arms paired per data set on the data sets both returned a DA table, with Monte Carlo
# intervals; the LP-LA difference is KIT-190's change with the rest of the route as on the bound commit (r1.21, R20-P02)
pairedTable <- function(arm_a, arm_b, jobs_a, jobs_b, metrics) {
    a <- dataset[dataset$arm == arm_a & dataset$job %in% jobs_a & dataset$has_da, ]
    b <- dataset[dataset$arm == arm_b & dataset$job %in% jobs_b & dataset$has_da, ]
    d <- merge(a, b, by = c("part", "cell", "scenario", "mask", "rep"), suffixes = c("_a", "_b"))
    d <- d[d$scenario %in% c("N0", "N1", "E0", "E1"), ]
    if (!nrow(d)) return(NULL)
    do.call(rbind, lapply(split(d, d[c("part", "cell", "scenario", "mask")], drop = TRUE), function(x) do.call(rbind, lapply(metrics, function(m) {
        delta <- x[[paste0(m, "_a")]] - x[[paste0(m, "_b")]]
        delta <- delta[is.finite(delta)]
        n <- length(delta)
        se <- if (n > 1L) stats::sd(delta) / sqrt(n) else NA_real_
        data.frame(x[1L, c("part", "cell", "scenario", "mask")], metric = m, arm_a = arm_a, arm_b = arm_b, n = n,
            mean_a = mean(x[[paste0(m, "_a")]], na.rm = TRUE), mean_b = mean(x[[paste0(m, "_b")]], na.rm = TRUE), mean_diff = if (n) mean(delta) else NA_real_,
            se = se, lower = if (n > 1L) mean(delta) - stats::qt(0.975, n - 1L) * se else NA_real_,
            upper = if (n > 1L) mean(delta) + stats::qt(0.975, n - 1L) * se else NA_real_, stringsAsFactors = FALSE)
    }))))
}
METRICS <- c("S1", "S2", "S3", "power", "calls", "samples_removed", "k_num", "n_controls_mean")
repair <- pairedTable("LA", "LP", "deciding", "lp", METRICS)
if (!is.null(repair)) utils::write.table(repair, file.path(out_dir, "partm-repair.tsv"), sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
levels_table <- pairedTable("LAs", "LA", "las", "deciding", METRICS)
if (!is.null(levels_table)) utils::write.table(levels_table, file.path(out_dir, "partm-levels.tsv"), sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")

# M-warn (5.8; informs D-default (d)): signals LA's journeys already record, against "the data set has a false call" (N scenarios) and
# "FDP above 0.10 or an attributable extreme fold-change error" (E scenarios), on L12 and M12, data sets 1 to 100. The first signal
# is the sweep's k above the cap (k_sweep > k_cap; the k used never exceeds the cap, r1.21, R20-P04).
la <- dataset[dataset$arm == "LA" & dataset$job == "deciding" & dataset$cell %in% c("L12", "M12") & dataset$rep <= 100L & dataset$has_da &
    dataset$scenario %in% c("N0", "N1", "E0", "E1") & dataset$ruv_reached, ]
if (nrow(la)) {
    ga <- dataset[dataset$arm == "GA" & dataset$job == "deciding", c("cell", "scenario", "mask", "rep", "k_est_min")]
    names(ga)[5L] <- "ga_estimate_min"
    la <- merge(la, ga, by = c("cell", "scenario", "mask", "rep"), all.x = TRUE)
    la$sig_k_above_cap <- la$cap_bound # k_sweep > k_cap in any assay (the data-set flag, SA.6 (10))
    la$sig_controls_below_10 <- la$n_controls_min < 10
    la$sig_complete_below_5 <- la$complete_min < 5
    la$sig_logfc_above_5 <- la$logfc_max > 5
    la$sig_ga_zero_la_removed <- la$ga_estimate_min == 0 & la$k_num >= 1
    la$target <- ifelse(startsWith(la$scenario, "N"), la$false_calls > 0, la$S3 > 0.10 | la$S5 %in% 1)
    warn <- do.call(rbind, lapply(c("sig_k_above_cap", "sig_controls_below_10", "sig_complete_below_5", "sig_logfc_above_5", "sig_ga_zero_la_removed"), function(sig) {
        do.call(rbind, lapply(c("N", "E"), function(kind) {
            x <- la[startsWith(la$scenario, kind) & !is.na(la[[sig]]), ]
            pos <- x[x$target %in% TRUE, ]
            neg <- x[x$target %in% FALSE, ]
            ci <- function(hits, n) if (n > 0L) clopper(hits, n, 0.975) else c(lower = NA_real_, upper = NA_real_)
            data.frame(signal = sig, scenarios = kind, n_positive = nrow(pos), n_negative = nrow(neg),
                sensitivity = if (nrow(pos)) mean(pos[[sig]]) else NA_real_, sens_lower = ci(sum(pos[[sig]]), nrow(pos))[["lower"]], sens_upper = ci(sum(pos[[sig]]), nrow(pos))[["upper"]],
                specificity = if (nrow(neg)) mean(!neg[[sig]]) else NA_real_, spec_lower = ci(sum(!neg[[sig]]), nrow(neg))[["lower"]], spec_upper = ci(sum(!neg[[sig]]), nrow(neg))[["upper"]],
                stringsAsFactors = FALSE)
        }))
    }))
    utils::write.table(warn, file.path(out_dir, "partm-warn.tsv"), sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
}
cat(sprintf("KIT-188 scoring: %d rows, %d data-set rows, %d conditions\n", nrow(rows), nrow(dataset), nrow(conditions)))
