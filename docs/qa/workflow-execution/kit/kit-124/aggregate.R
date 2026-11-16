# KIT-124 aggregation: the run-level outputs of simulators S1 to S4 against the ticket's decision rule (revision 1.2,
# W3) with the revision 1.3 additions (S1/S3, S4 and S2 pilot decisions; the confirming-study declaration). Evidence
# only: reads simulator outputs, writes tables and a sidecar. Base R only (stats, utils, tools); single-threaded.
#
# Usage (any subset of simulators; a directory is a simulator's out dir, e.g. kit-124/s1-confirm):
#   Rscript docs/qa/workflow-execution/kit/kit-124/aggregate.R <out dir> [S1=<dir>] [S2=<dir>] [S3=<dir>] [S4=<dir>]
# Reads, per simulator, the per-cell restart files when the out dir has them (S2 cells/<cell>.rds, S3
# cells/s3-<cell>.rds, S4 cells/<cell>.tsv with .notes.tsv and .timing.tsv), otherwise the top-level tables
# (s<i>-runs.tsv and the message table: s1-runs.log, s2-errors.tsv, s3-events.tsv, s4-notes.tsv; s2-data.tsv,
# s4-timing.tsv for the realised missing share and factor correlation).
# Writes to <out dir>:
#   k124-cells.tsv        one row per simulator x cell x arm x subset: runs, statuses, measures with Monte Carlo
#                         standard errors and one-sided 97.5% bounds, each rule's status, the row's verdict
#   k124-cell-status.tsv  one row per simulator x cell: class (core / stress / reported), run / not run, void
#   k124-arms.tsv         one row per simulator x arm x subset: the design verdict over the judged core cells, with
#                         the cells behind it; stress-cell results as conditions of use
#   k124-omics.tsv        one row per simulator (omic): eligible designs in the ticket's order of preference, arm 1
#                         on rules 1, 2, 6, 7, and the ticket's recommendation branch
#   k124-aggregate.meta.txt  sidecar (F-16): commit, script hash, R and package versions, command, input file hashes,
#                         and each simulator's own meta lines (commit, command, seed rule)
#
# ---------------------------------------------------------------------------------------------------------------------
# HOW THE RULES ARE READ (W3 and revision 1.3). Where the ticket is ambiguous, the most literal reading is taken and
# the point is listed under AMBIGUITIES.
#   Unit: one (simulator, cell, arm, subset). Runs are the repetitions; the Monte Carlo standard error (MCSE) is over
#   runs. Bounds are one-sided 97.5%: for shares of runs (FDP > 0.10, any call, crash) the exact Clopper-Pearson bound;
#   for per-run means (null rate, mean FDP, power, bias, RMSE, paired differences) mean +/- qnorm(0.975) * MCSE.
#   Statuses: ok; designed refusal (the simulator's `refused` on a package/simulator design arm); method refusal (a
#   comparator method's own stop: `refused` on arms 12, 13, 14, 14E, and rev 1.3's RUVinv stop "fourth-smallest
#   eigenvalue" recorded by S1/S2/S4 as `crashed` on arm 13); crash (every other `crashed`). Refusals of either kind
#   are scored as the run's no-RUV result (R-07, rev 1.3): the simulators already write arm 1's measures on
#   `refused` rows; a reclassified RUVinv row (written with no measures) takes arm 1's row of the same cell, run and
#   subset. Crashed rows have no measures and are left out of every measure; they count in rule 6.
#   Rule 1  null rate at 0.05 and 0.01: pass if both upper bounds <= 0.075 and 0.015; fail if either lower bound is
#           above its limit; else inconclusive.
#   Rule 2  mean FDP: pass if upper <= 0.075; fail if lower > 0.075.
#   Rule 3  share of runs with FDP > 0.10: limit = arm 2's point share in the same cell and subset + 0.05; pass if
#           upper <= limit; fail if lower > limit.
#   Rule 4  cells with a factor (all but C3, S13, S14): paired by run against arm 1 of the same cell and subset; pass
#           if the lower bound of the power gain > 0 or the lower bound of the reduction in the FDP > 0.10 indicator
#           > 0; fail if both upper bounds are <= 0 (a power gain that cannot be computed, as in C4, counts as not
#           above 0); else inconclusive.
#   Rule 5  cells without a factor (C3, S13, S14): power loss against arm 1, paired; pass if upper <= 0.02; fail if
#           lower > 0.02.
#   Rule 6  crash share (crashes / rows): pass if upper <= 0.01; fail if lower > 0.01.
#   Rule 7  C4 only: share of runs with any call; pass if upper <= 0.075; fail if lower > 0.075.
#   Row verdict: fail if any applicable rule fails; pass if all pass; else inconclusive.
#   Void: a (cell, subset) where arm 2 (the true factor) fails rule 1 or 2 is void; it is reported and judges no
#   design. A cell where arm 2 is only inconclusive on rule 1 or 2 is flagged (arm2_rule12), not void.
#   Core cells: C0 to C7, except the cells revision 1.3 makes "reported, not core": S3 C6 and S4 C6.
#   Design verdict (per simulator, arm and subset): over the core cells that are not void and where the arm runs:
#   fail if a judged cell fails; pass if every judged cell passes and no core cell is missing; not run if no core
#   cell has the arm; else inconclusive. Stress cells never fail a design: their failing rules are listed as
#   conditions of use.
#   Recommendation (per omic): RUV designs are arms 3 to 13, 15, 16, S3's R and S4's Qs and Qm; arms 1, 2, 14, 14E,
#   17, I and L are not RUV designs. Eligible = an RUV design whose primary-subset verdict is pass. Order: tier 1 the
#   method paper's usage (15, 16 with technical replicates; 8, 9 where controls are supplied); tier 2 arm 7; tier 3
#   the others by the lower mean rule-3 share over the judged core cells. Branch: eligible designs exist; else, arm 1
#   passes rules 1, 2, 6, 7 in every judged core cell -> "RUV off by default"; else "neither passes".
#   Primary subset: `full` where the arm has it, else `complete` (arms 12 to 14 in S1 and S4).
#   S2 arm 16 repeats arm 7 in C6 (rev 1.3): reported once, as arm 7; arm 16's S2 rows are dropped.
#
# AMBIGUITIES (literal reading taken; listed again in the report):
#   A1  "exact binomial bounds for shares": the null rate is called a share of true nulls, but a pooled binomial over
#       tests assumes independent tests (the unwanted factor makes them dependent; Efron 2007 in Standards) and the
#       number of tested nulls is not recoverable from S4's rows (true effects counted tested or not, contrasts
#       pooled). The measures heading says "each with its Monte Carlo standard error over runs": the null rate and
#       mean FDP are run means with MCSE; exact binomial bounds are used for shares of runs.
#   A2  Rule 3's limit "the true-factor arm's share plus 0.05": arm 2's point share, not its bound.
#   A3  Rule 4 "lower bound ... above 0": paired normal bounds on per-run differences (not binomial). Arm 1 against
#       itself has a gain of exactly 0 and so fails rule 4 in every core cell with a factor; the ticket's arm-1 test for
#       "RUV off" uses rules 1, 2, 6, 7 only, reported separately.
#   A4  "Fail ... with the whole interval outside": the interval is the two one-sided 97.5% bounds; for a "> 0" limit,
#       outside means upper <= 0.
#   A5  "A cell where the true-factor arm fails rules 1 or 2 is void": "fails" read in the ticket's own term (whole
#       interval outside); arm 2 inconclusive is flagged only. Void is decided per subset.
#   A6  "Refused runs carry their data set's no-RUV result": applied to method refusals too (rev 1.3 calls RUVinv's
#       stop a method refusal, and S4 already scores its RUVinv refusals as arm 1). The RUVinv reclassification
#       matches the message on arm 13 in any cell (rev 1.3 says "at 3 per group"; the message is the criterion).
#       In S4 a row is reclassified only if every crash note of that row is the RUVinv message; the whole row (both
#       assays) then takes arm 1's measures.
#   A7  Crashed rows are left out of the measures (no value to score) and counted in rule 6; the ticket does not say
#       how a crash enters the other measures.
#   A8  Which subset judges: R-13 scores arms on the full set and again on the complete subset; both are judged and
#       reported, the omic summary uses the primary subset (full, else complete).
#   A9  Which arms are "RUV designs" for the recommendation, and "the lower share in rule 3" (taken as the mean over
#       the judged core cells): not defined in the ticket; read as above.
#   A10 Rules 4 and 5 "cells with / without a factor": the no-factor cells are C3, S13, S14 in every simulator (S4's
#       "no unwanted variation"); S8 (known batch plus factor) counts as a factor cell.
#   A11 Rule 7 is applied in C4 only (the ticket's "Global null (C4)"), to every arm.
#   A12 Rule 4 in C4: C4 has a factor, so rule 4 applies literally; with no true effects the power gain is undefined
#       and only the reduction in the FDP > 0.10 share can pass it. Where arm 1 already makes (almost) no false calls
#       in C4, no design can show a reduction and rule 4 fails in C4 for every arm. `verdict` keeps the literal
#       reading; `verdict_without_rule4_in_C4` (k124-arms.tsv) and `eligible_without_rule4_in_C4` (k124-omics.tsv)
#       give the reading in which C4 is judged by rules 1-3, 6, 7 only. The choice is the user's.
#   A13 Design verdict "not judged": the arm has rows but no judged core cell (arms 15, 16 where C6 is reported, not
#       core; arm 17, stress only; or every core cell void). "not run": no rows.
# ---------------------------------------------------------------------------------------------------------------------

z975 <- stats::qnorm(0.975)
cells_all <- c(paste0("C", 0:7), paste0("S", 1:14))
no_factor_cells <- c("C3", "S13", "S14")
reported_not_core <- list(S3 = "C6", S4 = "C6")
comparator_arms <- c("12", "13", "14", "14E")
ruv_designs <- c(as.character(c(3:13, 15, 16)), "R", "Qs", "Qm")
omic_label <- c(S1 = "proteomics, classic route", S2 = "proteomics, DIA_LIMPA route",
    S3 = "transcriptomics", S4 = "lipidomics and metabolomics (one class simulated; sharing check)")
omic_notes <- c(
    S1 = "",
    S2 = "engine fastlimpadpc (rev 1.3); factor-group correlation set in expectation (realised per run in s2-data.tsv); S5 crashes are the route's on complete data (parked row 25); arm 16 = arm 7 in C6, reported once",
    S3 = "C6 reported, not core (rev 1.3; parked row 23); no 50%-missing cell (S6 not run; S5, S14 low counts); arms 4, 12, 13 not run (arm 4 = arm 3)",
    S4 = "C6 reported, not core (rev 1.3; parked rows 23, 24); missing share reported as realised (s4-timing.tsv); drift balanced against group by construction")
measure_cols <- c("k", "n_features", "n_true", "n_calls", "null_p05", "null_p01", "fdp", "power", "bias", "rmse", "seconds")

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2L) stop("usage: Rscript aggregate.R <out dir> [S1=<dir>] [S2=<dir>] [S3=<dir>] [S4=<dir>]")
out_dir <- args[[1L]]
sim_args <- args[-1L]
sim_dirs <- stats::setNames(sub("^[^=]*=", "", sim_args), sub("=.*$", "", sim_args))
if (!all(names(sim_dirs) %in% c("S1", "S2", "S3", "S4"))) stop("simulator arguments are S1=<dir> .. S4=<dir>")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
inputs <- character()

# ---- reading ----------------------------------------------------------------------------------------------------
readTsv <- function(path, header = TRUE, col_names = NULL) {
    inputs <<- c(inputs, normalizePath(path))
    if (is.null(col_names)) {
        return(utils::read.delim(path, header = header, colClasses = "character", quote = "", comment.char = "",
            fill = TRUE, na.strings = "NA"))
    }
    utils::read.delim(path, header = header, colClasses = "character", quote = "", comment.char = "",
        fill = TRUE, na.strings = "NA", col.names = col_names)
}
readRds <- function(path) {
    inputs <<- c(inputs, normalizePath(path))
    readRDS(path)
}
bindRows <- function(parts) {
    parts <- Filter(\(p) !is.null(p) && nrow(p) > 0L, parts)
    if (!length(parts)) return(NULL)
    parts <- lapply(parts, \(p) { p[] <- lapply(p, as.character); p })
    do.call(rbind, parts)
}
cellFiles <- function(dir, pattern) {
    files <- list.files(file.path(dir, "cells"), pattern = pattern, full.names = TRUE)
    files[!grepl("\\.tmp$", files)]
}
noMessages <- function() data.frame(cell = character(), run = character(), arm = character(), subset = character(),
    status = character(), message = character())

loadSim <- function(sim, dir) {
    messages <- noMessages(); context <- NULL
    if (sim == "S1") {
        runs <- readTsv(file.path(dir, "s1-runs.tsv"))
        log_path <- file.path(dir, "s1-runs.log")
        if (file.exists(log_path) && file.size(log_path) > 0) {
            lg <- readTsv(log_path, header = FALSE, col_names = c("cell", "run", "arm", "status", "message"))
            messages <- data.frame(cell = lg$cell, run = lg$run, arm = lg$arm, subset = NA_character_,
                status = lg$status, message = lg$message)
        }
    } else if (sim == "S2") {
        files <- cellFiles(dir, "\\.rds$")
        if (length(files)) {
            parts <- lapply(files, readRds)
            runs <- bindRows(lapply(parts, `[[`, "runs"))
            errors <- bindRows(lapply(parts, `[[`, "errors"))
            data <- bindRows(lapply(parts, `[[`, "data"))
        } else {
            runs <- readTsv(file.path(dir, "s2-runs.tsv"))
            errors <- readTsv(file.path(dir, "s2-errors.tsv"))
            data <- readTsv(file.path(dir, "s2-data.tsv"))
        }
        if (!is.null(errors) && nrow(errors)) messages <- data.frame(cell = errors$cell, run = errors$run,
            arm = errors$arm, subset = NA_character_, status = errors$status, message = errors$message)
        if (!is.null(data)) context <- data.frame(cell = data$cell, run = data$run,
            missing_share = as.numeric(data$missing_share), factor_group_cor = as.numeric(data$factor_group_cor))
    } else if (sim == "S3") {
        files <- cellFiles(dir, "^s3-.*\\.rds$")
        if (length(files)) {
            parts <- lapply(files, readRds)
            runs <- bindRows(lapply(parts, `[[`, "rows"))
            events <- bindRows(lapply(parts, `[[`, "events"))
        } else {
            runs <- readTsv(file.path(dir, "s3-runs.tsv"))
            events <- readTsv(file.path(dir, "s3-events.tsv"))
        }
        if (!is.null(events) && nrow(events)) messages <- data.frame(cell = events$cell, run = events$run,
            arm = events$arm, subset = NA_character_, status = events$status, message = events$message)
    } else if (sim == "S4") {
        files <- cellFiles(dir, "^[CS][0-9]+\\.tsv$")
        if (length(files)) {
            runs <- bindRows(lapply(files, readTsv))
            notes <- bindRows(lapply(cellFiles(dir, "\\.notes\\.tsv$"), readTsv))
            timing <- bindRows(lapply(cellFiles(dir, "\\.timing\\.tsv$"), readTsv))
        } else {
            runs <- readTsv(file.path(dir, "s4-runs.tsv"))
            notes <- readTsv(file.path(dir, "s4-notes.tsv"))
            timing <- readTsv(file.path(dir, "s4-timing.tsv"))
        }
        if (!is.null(notes) && nrow(notes)) messages <- data.frame(cell = notes$cell, run = notes$run, arm = notes$arm,
            subset = notes$subset, status = notes$status, message = notes$reason)
        if (!is.null(timing)) context <- data.frame(cell = timing$cell, run = timing$run,
            missing_share = as.numeric(timing$missing_share), factor_group_cor = NA_real_)
    }
    if (is.null(runs) || !nrow(runs)) return(NULL)
    # a table being appended by a live run can end in a partial line: keep complete rows only
    complete_line <- !is.na(runs$status) & runs$status %in% c("ok", "refused", "crashed") & !is.na(runs$seconds)
    complete_line[runs$status == "crashed"] <- TRUE
    runs <- runs[complete_line & !is.na(runs$subset) & nzchar(runs$subset), , drop = FALSE]
    runs$arm <- as.character(runs$arm)
    for (m in measure_cols) runs[[m]] <- suppressWarnings(as.numeric(runs[[m]]))
    list(runs = runs, messages = messages, context = context, meta = metaLines(sim, dir))
}
metaLines <- function(sim, dir) {
    path <- file.path(dir, sprintf("s%s-runs.meta.txt", substring(sim, 2L)))
    if (!file.exists(path)) return(character())
    inputs <<- c(inputs, normalizePath(path))
    lines <- readLines(path, warn = FALSE)
    grep("^(commit|script|command|seed rule|R:|R version|packages|DPC engine|runs|cells|started|finished|== )", lines,
        value = TRUE)
}

# ---- statuses -----------------------------------------------------------------------------------------------------
classify <- function(runs, messages, sim) {
    if (sim == "S2") runs <- runs[runs$arm != "16", , drop = FALSE]
    cls <- ifelse(runs$status == "ok", "ok", ifelse(runs$status == "refused",
        ifelse(runs$arm %in% comparator_arms, "refused_method", "refused_designed"), "crashed"))
    ruvinv <- grepl("fourth-smallest eigenvalue", messages$message, fixed = TRUE) & messages$arm == "13" &
        messages$status == "crashed"
    crash_msgs <- messages[messages$status == "crashed", , drop = FALSE]
    keyOf <- function(x, with_subset) paste(x$cell, x$run, x$arm, if (with_subset) x$subset else rep("", nrow(x)), sep = "|")
    with_subset <- sim == "S4"
    all_ruvinv <- tapply(ruvinv[messages$status == "crashed"], keyOf(crash_msgs, with_subset), all)
    row_keys <- keyOf(runs, with_subset)
    hit <- cls == "crashed" & runs$arm == "13" & row_keys %in% names(all_ruvinv)[all_ruvinv]
    if (any(hit)) {
        ref <- runs[runs$arm == "1", , drop = FALSE]
        idx <- match(paste(runs$cell[hit], runs$run[hit], runs$subset[hit]), paste(ref$cell, ref$run, ref$subset))
        for (m in setdiff(measure_cols, c("k", "seconds"))) runs[[m]][hit] <- ref[[m]][idx]
        cls[hit] <- "refused_method"
    }
    runs$class <- cls
    runs
}

# ---- bounds -------------------------------------------------------------------------------------------------------
cpUpper <- function(x, n) ifelse(n == 0, NA_real_, ifelse(x >= n, 1, stats::qbeta(0.975, x + 1, n - x)))
cpLower <- function(x, n) ifelse(n == 0, NA_real_, ifelse(x <= 0, 0, stats::qbeta(0.025, x, n - x + 1)))
meanSe <- function(v) {
    v <- v[is.finite(v)]
    n <- length(v)
    if (!n) return(c(mean = NA_real_, se = NA_real_, n = 0))
    c(mean = mean(v), se = if (n > 1L) stats::sd(v) / sqrt(n) else NA_real_, n = n)
}
judgeUpper <- function(lower, upper, limit) {
    ifelse(is.na(upper) | is.na(limit), NA_character_, ifelse(upper <= limit, "pass", ifelse(lower > limit, "fail",
        "inconclusive")))
}

# ---- one simulator ------------------------------------------------------------------------------------------------
cellTable <- function(sim, runs, context) {
    runs$key <- paste(runs$cell, runs$subset, runs$run, sep = "|")
    arm1 <- runs[runs$arm == "1", , drop = FALSE]
    groups <- split(seq_len(nrow(runs)), paste(runs$cell, runs$arm, runs$subset, sep = "\t"))
    rows <- lapply(groups, function(i) {
        r <- runs[i, , drop = FALSE]
        scored <- r[r$class != "crashed" & !is.na(r$fdp), , drop = FALSE]
        n_rows <- nrow(r); n_sc <- nrow(scored)
        nul5 <- meanSe(scored$null_p05); nul1 <- meanSe(scored$null_p01); fdp <- meanSe(scored$fdp)
        pw <- meanSe(scored$power); bias <- meanSe(scored$bias); rmse <- meanSe(scored$rmse)
        x10 <- sum(scored$fdp > 0.10); xany <- sum(scored$n_calls > 0); xcr <- sum(r$class == "crashed")
        a1 <- arm1[match(scored$key, arm1$key), , drop = FALSE]
        ok1 <- a1$class != "crashed" & !is.na(a1$fdp)
        gain <- meanSe((scored$power - a1$power)[ok1])
        red <- meanSe(((a1$fdp > 0.10) - (scored$fdp > 0.10))[ok1])
        data.frame(simulator = sim, cell = r$cell[[1L]], arm = r$arm[[1L]], subset = r$subset[[1L]],
            runs = n_rows, ok = sum(r$class == "ok"), refused_designed = sum(r$class == "refused_designed"),
            refused_method = sum(r$class == "refused_method"), crashed = xcr, scored = n_sc,
            mean_k = if (n_sc) mean(scored$k, na.rm = TRUE) else NA_real_,
            mean_n_features = if (n_sc) mean(scored$n_features) else NA_real_,
            null05 = nul5[["mean"]], null05_se = nul5[["se"]], null01 = nul1[["mean"]], null01_se = nul1[["se"]],
            mean_fdp = fdp[["mean"]], mean_fdp_se = fdp[["se"]],
            fdp10_share = if (n_sc) x10 / n_sc else NA_real_, fdp10_lower = cpLower(x10, n_sc),
            fdp10_upper = cpUpper(x10, n_sc),
            anycall_share = if (n_sc) xany / n_sc else NA_real_, anycall_lower = cpLower(xany, n_sc),
            anycall_upper = cpUpper(xany, n_sc),
            power = pw[["mean"]], power_se = pw[["se"]], bias = bias[["mean"]], bias_se = bias[["se"]],
            rmse = rmse[["mean"]], rmse_se = rmse[["se"]],
            power_gain = gain[["mean"]], power_gain_se = gain[["se"]], paired_runs = red[["n"]],
            fdp10_reduction = red[["mean"]], fdp10_reduction_se = red[["se"]],
            crash_share = xcr / n_rows, crash_lower = cpLower(xcr, n_rows), crash_upper = cpUpper(xcr, n_rows),
            stringsAsFactors = FALSE)
    })
    t <- do.call(rbind, rows)
    rownames(t) <- NULL
    lo <- function(m, se) m - z975 * ifelse(is.na(se), 0, se)
    hi <- function(m, se) m + z975 * ifelse(is.na(se), 0, se)
    # one run gives no MCSE: its bound is then not computable
    noSe <- function(se, m) is.na(se) & !is.na(m)
    t$null05_upper <- ifelse(noSe(t$null05_se, t$null05), NA, hi(t$null05, t$null05_se))
    t$null01_upper <- ifelse(noSe(t$null01_se, t$null01), NA, hi(t$null01, t$null01_se))
    t$mean_fdp_upper <- ifelse(noSe(t$mean_fdp_se, t$mean_fdp), NA, hi(t$mean_fdp, t$mean_fdp_se))
    t$power_gain_lower <- lo(t$power_gain, t$power_gain_se); t$power_gain_upper <- hi(t$power_gain, t$power_gain_se)
    t$fdp10_reduction_lower <- lo(t$fdp10_reduction, t$fdp10_reduction_se)
    t$fdp10_reduction_upper <- hi(t$fdp10_reduction, t$fdp10_reduction_se)

    # rules
    r1a <- judgeUpper(lo(t$null05, t$null05_se), t$null05_upper, 0.075)
    r1b <- judgeUpper(lo(t$null01, t$null01_se), t$null01_upper, 0.015)
    t$rule1 <- ifelse(is.na(r1a) | is.na(r1b), NA, ifelse(r1a == "fail" | r1b == "fail", "fail",
        ifelse(r1a == "pass" & r1b == "pass", "pass", "inconclusive")))
    t$rule2 <- judgeUpper(lo(t$mean_fdp, t$mean_fdp_se), t$mean_fdp_upper, 0.075)
    a2 <- t[t$arm == "2", c("cell", "subset", "fdp10_share", "rule1", "rule2")]
    m2 <- match(paste(t$cell, t$subset), paste(a2$cell, a2$subset))
    t$rule3_limit <- a2$fdp10_share[m2] + 0.05
    t$rule3 <- judgeUpper(t$fdp10_lower, t$fdp10_upper, t$rule3_limit)
    factor_cell <- !(t$cell %in% no_factor_cells)
    gain_pass <- !is.na(t$power_gain_lower) & t$power_gain_lower > 0
    red_pass <- !is.na(t$fdp10_reduction_lower) & t$fdp10_reduction_lower > 0
    gain_out <- is.na(t$power_gain_upper) | t$power_gain_upper <= 0
    red_out <- is.na(t$fdp10_reduction_upper) | t$fdp10_reduction_upper <= 0
    t$rule4 <- ifelse(!factor_cell | t$paired_runs == 0, NA, ifelse(gain_pass | red_pass, "pass",
        ifelse(gain_out & red_out, "fail", "inconclusive")))
    t$rule5 <- ifelse(factor_cell | t$paired_runs == 0, NA,
        judgeUpper(-t$power_gain_upper, -t$power_gain_lower, 0.02))
    t$rule6 <- judgeUpper(t$crash_lower, t$crash_upper, 0.01)
    t$rule7 <- ifelse(t$cell == "C4", judgeUpper(t$anycall_lower, t$anycall_upper, 0.075), NA)
    rules <- paste0("rule", 1:7)
    t$verdict <- apply(t[, rules], 1L, function(v) {
        v <- v[!is.na(v)]
        if (!length(v)) "not run" else if (any(v == "fail")) "fail" else if (all(v == "pass")) "pass" else "inconclusive"
    })
    t$arm2_rule1 <- a2$rule1[m2]; t$arm2_rule2 <- a2$rule2[m2]
    t$void <- !is.na(m2) & ((t$arm2_rule1 %in% "fail") | (t$arm2_rule2 %in% "fail"))
    t$arm2_rule12 <- ifelse(is.na(m2), NA, ifelse(t$void, "fail", ifelse(t$arm2_rule1 %in% "pass" &
        t$arm2_rule2 %in% "pass", "pass", "inconclusive")))
    t$cell_role <- ifelse(t$cell %in% reported_not_core[[sim]], "reported", ifelse(startsWith(t$cell, "C"), "core",
        "stress"))
    if (!is.null(context)) {
        ms <- tapply(context$missing_share, context$cell, mean, na.rm = TRUE)
        fc <- tapply(context$factor_group_cor, context$cell, mean, na.rm = TRUE)
        t$realised_missing_share <- as.numeric(ms[t$cell]); t$realised_factor_group_cor <- as.numeric(fc[t$cell])
    } else {
        t$realised_missing_share <- NA_real_; t$realised_factor_group_cor <- NA_real_
    }
    t[order(match(t$cell, cells_all), suppressWarnings(as.numeric(t$arm)), t$arm, t$subset), ]
}

armTable <- function(sim, t) {
    run_cells <- unique(t$cell)
    core_cells <- setdiff(paste0("C", 0:7), reported_not_core[[sim]])
    core_missing <- setdiff(core_cells, run_cells)
    keys <- unique(t[, c("arm", "subset")])
    out <- lapply(seq_len(nrow(keys)), function(j) {
        a <- keys$arm[[j]]; s <- keys$subset[[j]]
        x <- t[t$arm == a & t$subset == s, , drop = FALSE]
        core <- x[x$cell_role == "core", , drop = FALSE]
        judged <- core[!core$void, , drop = FALSE]
        stress <- x[x$cell_role == "stress", , drop = FALSE]
        failing <- function(y) {
            if (!nrow(y)) return("")
            f <- y[y$verdict == "fail", , drop = FALSE]
            if (!nrow(f)) return("")
            paste(vapply(seq_len(nrow(f)), function(i) paste0(f$cell[[i]], "(", paste(sub("rule", "",
                paste0("rule", 1:7)[which(unlist(f[i, paste0("rule", 1:7)]) %in% "fail")]), collapse = ","), ")"), ""),
                collapse = " ")
        }
        # design verdict from the judged cells' rule statuses; `drop_c4_rule4` gives the A12 alternative reading
        over <- function(rules, drop_c4_rule4 = FALSE) {
            if (!nrow(judged)) return(if (nrow(x)) "not judged" else "not run")
            m <- as.matrix(judged[, rules, drop = FALSE])
            if (drop_c4_rule4 && "rule4" %in% rules) m[judged$cell == "C4", "rule4"] <- NA
            v <- m[!is.na(m)]
            if (any(v == "fail")) "fail" else if (all(v == "pass") && !length(core_missing)) "pass" else "inconclusive"
        }
        rules <- paste0("rule", 1:7)
        verdict <- over(rules)
        verdict_alt <- over(rules, drop_c4_rule4 = TRUE)
        arm1_1267 <- if (a != "1") NA_character_ else over(c("rule1", "rule2", "rule6", "rule7"))
        data.frame(simulator = sim, arm = a, subset = s,
            primary = s == "full" || !any(t$arm == a & t$subset == "full"),
            ruv_design = a %in% ruv_designs, verdict = verdict, verdict_without_rule4_in_C4 = verdict_alt,
            judged_core_cells = paste(judged$cell, collapse = ","),
            void_core_cells = paste(core$cell[core$void], collapse = ","),
            core_cells_not_run = paste(core_missing, collapse = ","),
            failing_core_cells = failing(judged),
            inconclusive_core_cells = paste(judged$cell[judged$verdict == "inconclusive"], collapse = ","),
            stress_cells_failing = failing(stress[!stress$void, , drop = FALSE]),
            stress_cells_void = paste(stress$cell[stress$void], collapse = ","),
            reported_cells = paste(x$cell[x$cell_role == "reported"], collapse = ","),
            mean_fdp10_share_core = if (nrow(judged)) mean(judged$fdp10_share, na.rm = TRUE) else NA_real_,
            power_C0 = if (any(x$cell == "C0")) x$power[x$cell == "C0"] else NA_real_,
            refused_designed = sum(x$refused_designed), refused_method = sum(x$refused_method),
            crashed = sum(x$crashed), runs = sum(x$runs), arm1_rules_1_2_6_7 = arm1_1267,
            stringsAsFactors = FALSE)
    })
    a <- do.call(rbind, out)
    a[order(suppressWarnings(as.numeric(a$arm)), a$arm, a$subset), ]
}

omicRow <- function(sim, a) {
    p <- a[a$primary, , drop = FALSE]
    eligible <- function(verdicts) {
        e <- p[p$ruv_design & verdicts == "pass", , drop = FALSE]
        tier <- ifelse(e$arm %in% c("15", "16", "8", "9"), 1L, ifelse(e$arm == "7", 2L, 3L))
        e[order(tier, ifelse(tier == 3L, e$mean_fdp10_share_core, 0)), , drop = FALSE]
    }
    elig <- eligible(p$verdict)
    elig_alt <- eligible(p$verdict_without_rule4_in_C4)
    arm1 <- p[p$arm == "1", , drop = FALSE]
    arm1_1267 <- if (nrow(arm1)) arm1$arm1_rules_1_2_6_7[[1L]] else NA_character_
    branch <- if (nrow(elig)) "eligible designs exist" else if (identical(arm1_1267, "pass"))
        "RUV off by default (no RUV design passes; arm 1 passes rules 1, 2, 6, 7)" else
        if (identical(arm1_1267, "fail")) "neither passes" else
        "no RUV design passes; arm 1 on rules 1, 2, 6, 7 not decided (inconclusive or not run)"
    count <- function(v) sum(p$ruv_design & p$verdict == v)
    data.frame(simulator = sim, omic = omic_label[[sim]],
        eligible_in_order = paste(elig$arm, collapse = " > "),
        ruv_designs_pass = count("pass"), ruv_designs_fail = count("fail"),
        ruv_designs_inconclusive = count("inconclusive"), ruv_designs_not_judged = count("not judged"),
        ruv_designs_not_run = count("not run"),
        arm1_rules_1_2_6_7 = arm1_1267, branch = branch,
        eligible_without_rule4_in_C4 = paste(elig_alt$arm, collapse = " > "),
        void_core_cells = if (nrow(arm1)) arm1$void_core_cells[[1L]] else "",
        core_cells_not_run = if (nrow(p)) p$core_cells_not_run[[1L]] else "",
        notes = omic_notes[[sim]], stringsAsFactors = FALSE)
}

statusTable <- function(sim, t) {
    do.call(rbind, lapply(cells_all, function(cl) {
        x <- t[t$cell == cl, , drop = FALSE]
        role <- if (cl %in% reported_not_core[[sim]]) "reported" else if (startsWith(cl, "C")) "core" else "stress"
        full_or_first <- x[x$arm == "1", , drop = FALSE]
        data.frame(simulator = sim, cell = cl, role = role, status = if (nrow(x)) "run" else "not run",
            runs = if (nrow(full_or_first)) max(full_or_first$runs) else 0L,
            void_subsets = paste(unique(x$subset[x$void]), collapse = ","),
            arm2_rule12 = if (any(x$arm == "2")) paste(unique(paste0(x$subset[x$arm == "2"], ":", x$arm2_rule12[x$arm == "2"])), collapse = " ") else "",
            realised_missing_share = if (nrow(x)) x$realised_missing_share[[1L]] else NA_real_,
            realised_factor_group_cor = if (nrow(x)) x$realised_factor_group_cor[[1L]] else NA_real_,
            stringsAsFactors = FALSE)
    }))
}

# ---- run ----------------------------------------------------------------------------------------------------------
cells_out <- list(); arms_out <- list(); omics_out <- list(); status_out <- list(); metas <- list()
for (sim in intersect(c("S1", "S2", "S3", "S4"), names(sim_dirs))) {
    loaded <- loadSim(sim, sim_dirs[[sim]])
    if (is.null(loaded)) { message(sim, ": no rows"); next }
    runs <- classify(loaded$runs, loaded$messages, sim)
    t <- cellTable(sim, runs, loaded$context)
    a <- armTable(sim, t)
    cells_out[[sim]] <- t; arms_out[[sim]] <- a; omics_out[[sim]] <- omicRow(sim, a)
    status_out[[sim]] <- statusTable(sim, t); metas[[sim]] <- loaded$meta
    message(sprintf("%s: %d rows, cells %s", sim, nrow(runs), paste(unique(runs$cell), collapse = ",")))
}

writeOut <- function(x, name) {
    x[] <- lapply(x, \(v) if (is.double(v)) signif(v, 6) else v)
    utils::write.table(x, file.path(out_dir, name), sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
}
writeOut(do.call(rbind, unname(cells_out)), "k124-cells.tsv")
writeOut(do.call(rbind, unname(status_out)), "k124-cell-status.tsv")
writeOut(do.call(rbind, unname(arms_out)), "k124-arms.tsv")
writeOut(do.call(rbind, unname(omics_out)), "k124-omics.tsv")

# ---- sidecar (F-16) -----------------------------------------------------------------------------------------------
script <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
script <- if (length(script)) normalizePath(script) else NA_character_
git <- function(...) tryCatch(system2("git", c("-C", dirname(script), ...), stdout = TRUE, stderr = FALSE),
    error = \(e) NA_character_, warning = \(w) NA_character_)
inputs <- unique(inputs)
hashes <- if (length(inputs)) tools::sha256sum(inputs) else character()
lines <- c(
    "KIT-124 aggregation (ticket revision 1.2 rules with the revision 1.3 additions)",
    paste("script:", script),
    paste("script sha256:", if (!is.na(script)) unname(tools::sha256sum(script)) else NA),
    paste("commit:", paste(git("rev-parse", "HEAD"), collapse = " ")),
    paste("script differs from commit:", length(git("status", "--porcelain", "--", basename(script))) > 0L),
    paste("command:", paste(c("Rscript", script, commandArgs(TRUE)), collapse = " ")),
    paste("working directory:", getwd()),
    paste("R:", R.version.string, ";", R.version$platform),
    paste("packages:", paste(vapply(c("base", "stats", "utils", "tools"), \(p)
        paste(p, as.character(utils::packageVersion(p))), ""), collapse = "; ")),
    paste("library paths:", paste(.libPaths(), collapse = ":")),
    paste("run at:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
    "inputs (sha256, bytes, path):",
    sprintf("  %s %s %s", hashes, file.size(inputs), inputs))
for (sim in names(metas)) lines <- c(lines, sprintf("simulator %s meta (%s):", sim, sim_dirs[[sim]]),
    paste0("  ", metas[[sim]]))
writeLines(lines, file.path(out_dir, "k124-aggregate.meta.txt"))
message("written: ", out_dir)
