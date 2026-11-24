# KIT-189 scoring (5.5 to 5.7; VAL-86): intervals, the 108 classifications, the outcome, the default and each
# option's trade-off text, from the committed rows. Hashed in the declaration; run once on a runner (SC).
#
# Usage, from the package root:
#   Rscript docs/qa/workflow-execution/kit/kit-189/score-k189.R <rows dir> <out dir> [<blind summary file>]
#   rows dir  every *.tsv.gz and complete part-<cell>-<r>.tsv under it (recursively) is read; a .tmp part never is
#   out dir   cells.tsv, classifications.tsv, paired.tsv, mechanism.tsv, checks-score.tsv, outcome.tsv, options.tsv
#   blind summary file  when given (the lifecycle run), nothing about the results is printed: the tables go to
#             <out dir> (not uploaded) and only their names, sizes and SHA-256 are written to this file

args <- commandArgs(TRUE)
stopifnot(length(args) %in% c(2L, 3L))
rows_dir <- args[[1L]]
out_dir <- args[[2L]]
blind_file <- if (length(args) == 3L) args[[3L]] else NULL
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
if (!is.null(blind_file)) dir.create(dirname(blind_file), showWarnings = FALSE, recursive = TRUE)
say <-function(...) if (is.null(blind_file)) writeLines(sprintf(...))

LINE <- 0.075
N_CLASS <- 108
N_POWER <- 18
MARGIN <- 0.01
DECIDING_ARMS <- c("Q", "BH", "BHf")
NULL_CELLS <- 1:9
DECIDING_PLANTED <- c(10L, 11L, 12L, 15L, 16L, 17L, 20L, 21L, 22L)
CELLS <- data.frame(
    cell = 1:24,
    background = c(rep(c("B1", "B2", "B3"), each = 3L), rep(c("B1", "B2", "B3"), each = 5L)),
    scenario = c(rep("null", 9L), rep(c("small", "mid", "large", "mid", "mid"), 3L)),
    L = c(rep(c(25L, 126L, 400L), 3L), rep(c(126L, 126L, 126L, 25L, 400L), 3L)),
    R_declared = c(rep(2000L, 9L), rep(c(400L, 400L, 400L, 200L, 200L), 3L)),
    stringsAsFactors = FALSE)
BG_LABEL <- c(B1 = "E. coli full proteome", B2 = "E. coli 2,000-protein subset", B3 = "yeast")
PROP_LABEL <- c(none = "no parent terms", is_a_part_of_regulates = "parent terms")
ARM_NAME <- c(Q = "The q-value (as shipped)", BH = "Benjamini-Hochberg over the terms with at least one list protein",
    BHf = "Benjamini-Hochberg over the whole size-filtered family")
ASSUMPTION <- c(
    Q = paste("Storey's estimate of the share of true nulls assumes null p-values uniform on (0, 1); the route adjusts",
        "only the terms with at least one list protein, which breaks that assumption for small sets"),
    BH = paste("Benjamini-Hochberg over the terms with at least one list protein; that family is selected on the",
        "p-value (p below 1), for which Benjamini-Hochberg has no general guarantee"),
    BHf = paste("Benjamini-Hochberg over every term of 10 to 500 proteins; its guarantee needs independent or positively",
        "dependent tests, which overlapping Gene Ontology sets do not ensure"))

# ---- read ----------------------------------------------------------------------------------------------------------
files <- list.files(rows_dir, pattern = "(\\.tsv\\.gz|^part-[0-9]+-[0-9]+\\.tsv)$", recursive = TRUE, full.names = TRUE)
stopifnot(length(files) > 0L)
rows <- do.call(rbind, lapply(files, function(f) utils::read.delim(f, colClasses = "character", na.strings = "NA", quote = "")))
num <- c("cell", "L", "replicate", "m", "m_kept", "calls", "calls_null", "calls_nonnull", "calls_below15",
    "calls_material", "planted_size", "planted_s", "n_p_ge_05", "n_p_ge_095", "n_p_eq_1", "pi0_hat", "pi0_spread",
    "chain_pi0", "median_er_null", "er_dev", "seconds")
lgl <- c("planted_tested", "planted_called", "q_all_na", "chain_equal", "e2_ok")
for (k in num) rows[[k]] <- as.numeric(rows[[k]])
for (k in lgl) rows[[k]] <- as.logical(rows[[k]])

# ---- integrity (each replicate read once; complete; paired) --------------------------------------------------------
key <- paste(rows$cell, rows$replicate, rows$propagation, rows$arm)
per_rep <- table(paste(rows$cell, rows$replicate))
sha_per_rep <- tapply(rows$list_sha, paste(rows$cell, rows$replicate), function(x) length(unique(x)))
q_ok <- rows[rows$arm == "Q" & rows$status == "ok", , drop = FALSE]
checks <- data.frame(check = c("rows read once", "replicates complete (14 rows)", "one list per replicate",
    "files read", "E2 flags on every ok call", "nesting BHf in BH, BH in Q", "chain reproduces the q-values"),
    pass = c(!anyDuplicated(key), all(per_rep == 14L), all(sha_per_rep == 1L), TRUE, all(q_ok$e2_ok),
        !any(grepl("nest", q_ok$e2_fail)), all(q_ok$chain_equal)),
    detail = c(sprintf("%d rows, %d duplicated keys", nrow(rows), sum(duplicated(key))),
        sprintf("%d replicates, %d incomplete", length(per_rep), sum(per_rep != 14L)),
        sprintf("%d replicates with more than one list", sum(sha_per_rep != 1L)), sprintf("%d files", length(files)),
        sprintf("%d of %d ok calls fail a check (%s)", sum(!q_ok$e2_ok), nrow(q_ok),
            paste(unique(q_ok$e2_fail[!q_ok$e2_ok]), collapse = ", ")),
        sprintf("%d calls", sum(grepl("nest", q_ok$e2_fail))), sprintf("%d calls differ", sum(!q_ok$chain_equal))),
    stringsAsFactors = FALSE)
if (anyDuplicated(key)) stop("a replicate's rows are present twice: scoring stops (KIT-189 10, Attempts)")

# ---- interval helpers (5.6) ----------------------------------------------------------------------------------------
cpUpper <- function(x, n, conf) if (n == 0) NA_real_ else if (x >= n) 1 else stats::qbeta(conf, x + 1, n - x)
cpLower <- function(x, n, alpha) if (n == 0) NA_real_ else if (x <= 0) 0 else stats::qbeta(alpha, x, n - x + 1)
meanStats <- function(x, alpha_b) {
    n <- length(x)
    if (n < 2L) return(c(mean = if (n) mean(x) else NA_real_, se = NA_real_, lo95 = NA_real_, hi95 = NA_real_,
        upper1 = NA_real_, lowerB = NA_real_, resolvable = NA_real_))
    m <- mean(x)
    se <- stats::sd(x) / sqrt(n)
    t2 <- stats::qt(0.975, n - 1L)
    c(mean = m, se = se, lo95 = m - t2 * se, hi95 = m + t2 * se, upper1 = m + stats::qt(0.95, n - 1L) * se,
        lowerB = m - stats::qt(1 - alpha_b, n - 1L) * se, resolvable = stats::qt(1 - alpha_b, n - 1L) * se)
}
rateStats <- function(x, n, alpha_b) {
    p <- if (n) x / n else NA_real_
    first_fail <- if (n) which(vapply(0:n, function(k) isTRUE(cpLower(k, n, alpha_b) > LINE), logical(1)))[1L] - 1L else NA
    c(x = x, n = n, estimate = p, mc_error = if (n) sqrt(p * (1 - p) / n) else NA_real_,
        lo95 = if (n) cpLower(x, n, 0.025) else NA_real_, hi95 = if (n) cpUpper(x, n, 0.975) else NA_real_,
        upper1 = cpUpper(x, n, 0.95), lowerB = cpLower(x, n, alpha_b),
        resolvable = if (is.na(first_fail)) NA_real_ else first_fail / n - LINE)
}
classify <- function(upper1, lowerB, crashed_share) {
    if (is.na(upper1) || is.na(lowerB)) return("not determined")
    if (crashed_share > 0.05) return("not determined (crashed calls)")
    if (upper1 < LINE) "controls" else if (lowerB > LINE) "fails" else "not determined"
}

# ---- per cell, propagation value and arm (5.5 (a) to (e)) ----------------------------------------------------------
groups <- split(rows, list(rows$cell, rows$propagation, rows$arm), drop = TRUE)
cell_rows <- lapply(groups, function(d) {
    cell <- d$cell[[1L]]
    info <- CELLS[CELLS$cell == cell, ]
    crashed <- sum(d$status == "crashed")
    v <- d[d$status != "crashed", , drop = FALSE]
    R <- nrow(v)
    out <- data.frame(cell = cell, background = info$background, scenario = info$scenario, L = info$L,
        propagation = d$propagation[[1L]], arm = d$arm[[1L]], R_declared = info$R_declared, R = R, crashed = crashed,
        empty = sum(v$status == "empty"), q_all_na_lists = sum(v$q_all_na, na.rm = TRUE), stringsAsFactors = FALSE)
    if (identical(info$scenario, "null")) {
        a <- rateStats(sum(v$calls > 0), R, 0.05 / N_CLASS)
        b <- meanStats(v$calls, 0.05 / N_CLASS)
        c10 <- rateStats(sum(v$calls >= 10), R, 0.05 / N_CLASS)
        out <- cbind(out, data.frame(estimand = "a_any_call", estimate = a[["estimate"]], mc_error = a[["mc_error"]],
            lo95 = a[["lo95"]], hi95 = a[["hi95"]], upper1 = a[["upper1"]], lowerB = a[["lowerB"]],
            resolvable_excess = a[["resolvable"]], lists_with_call = a[["x"]],
            false_terms_mean = b[["mean"]], false_terms_mc_error = b[["se"]], false_terms_lo95 = b[["lo95"]],
            false_terms_hi95 = b[["hi95"]], false_terms_median = stats::median(v$calls),
            false_terms_p90 = unname(stats::quantile(v$calls, 0.9, type = 7)),
            false_terms_p99 = unname(stats::quantile(v$calls, 0.99, type = 7)), false_terms_max = max(v$calls),
            large_events = c10[["x"]], large_rate = c10[["estimate"]], large_lo95 = c10[["lo95"]],
            large_hi95 = c10[["hi95"]], large_upper1 = c10[["upper1"]], stringsAsFactors = FALSE))
    } else {
        fdp <- ifelse(v$calls > 0, v$calls_null / v$calls, 0)
        fdp15 <- ifelse(v$calls > 0, v$calls_below15 / v$calls, 0)
        d1 <- meanStats(fdp, 0.05 / N_CLASS)
        d15 <- meanStats(fdp15, 0.05 / N_CLASS)
        tested <- v$planted_tested %in% TRUE
        e1 <- rateStats(sum(v$planted_called[tested] %in% TRUE), sum(tested), 0.05 / N_CLASS)
        e2 <- meanStats(v$calls_material, 0.05 / N_CLASS)
        e3 <- rateStats(sum(v$calls_nonnull > 0), R, 0.05 / N_CLASS)
        out <- cbind(out, data.frame(estimand = "d_fdp", estimate = d1[["mean"]], mc_error = d1[["se"]],
            lo95 = d1[["lo95"]], hi95 = d1[["hi95"]], upper1 = d1[["upper1"]], lowerB = d1[["lowerB"]],
            resolvable_excess = d1[["resolvable"]], fdp15 = d15[["mean"]], fdp15_mc_error = d15[["se"]],
            fdp15_upper1 = d15[["upper1"]], planted_tested = sum(tested), power_planted = e1[["estimate"]],
            power_planted_lo95 = e1[["lo95"]], power_planted_hi95 = e1[["hi95"]], power_material_mean = e2[["mean"]],
            power_material_mc_error = e2[["se"]], power_any_nonnull = e3[["estimate"]],
            power_any_nonnull_lo95 = e3[["lo95"]], power_any_nonnull_hi95 = e3[["hi95"]],
            median_er_null_terms = stats::median(v$median_er_null, na.rm = TRUE), stringsAsFactors = FALSE))
    }
    out$classification <- classify(out$upper1, out$lowerB, if (nrow(d)) crashed / nrow(d) else 0)
    out
})
bindFill <- function(xs) {
    cols <- unique(unlist(lapply(xs, names)))
    do.call(rbind, lapply(xs, function(x) {
        for (k in setdiff(cols, names(x))) x[[k]] <- NA
        x[, cols]
    }))
}
cells <- bindFill(cell_rows)
cells <- cells[order(cells$cell, cells$propagation != "none", match(cells$arm, c("Q", "BH", "BHf", "Qf", "BY", "Holm", "raw"))), ]
rownames(cells) <- NULL

deciding <- cells[cells$arm %in% DECIDING_ARMS & (cells$cell %in% NULL_CELLS | cells$cell %in% DECIDING_PLANTED), ]
classifications <- deciding[, c("cell", "background", "scenario", "L", "propagation", "arm", "estimand", "R",
    "estimate", "upper1", "lowerB", "classification")]

# ---- paired (5.5 (e), (f); 5.7 power gain and large-event excess) --------------------------------------------------
wide <- function(d, arm) {
    x <- d[d$arm == arm, ]
    x[order(x$replicate), ]
}
paired_rows <- lapply(split(rows, list(rows$cell, rows$propagation), drop = TRUE), function(d) {
    cell <- d$cell[[1L]]
    info <- CELLS[CELLS$cell == cell, ]
    q <- wide(d, "Q")
    out <- list()
    for (other in c("BH", "BHf")) {
        o <- wide(d, other)
        stopifnot(identical(q$replicate, o$replicate))
        okk <- q$status != "crashed" & o$status != "crashed"
        diff <- q$calls[okk] - o$calls[okk]
        row <- data.frame(cell = cell, background = info$background, scenario = info$scenario, L = info$L,
            propagation = d$propagation[[1L]], comparison = paste0("Q-", other), lists = sum(okk),
            calls_diff_mean = mean(diff), calls_diff_max = max(diff), calls_diff_min = min(diff),
            lists_differ = mean(diff != 0),
            negative_on_non_na = sum(diff < 0 & !q$q_all_na[okk]), stringsAsFactors = FALSE)
        if (identical(info$scenario, "null")) {
            bq <- sum(q$calls[okk] > 0 & o$calls[okk] == 0)
            co <- sum(q$calls[okk] == 0 & o$calls[okk] > 0)
            ind <- as.numeric(q$calls[okk] > 0) - as.numeric(o$calls[okk] > 0)
            ms <- meanStats(ind, 0.05 / N_CLASS)
            non_na <- okk & !q$q_all_na
            ex <- rateStats(sum(q$calls[non_na] >= 10 & o$calls[non_na] < 10), sum(non_na), 0.05)
            row <- cbind(row, data.frame(any_call_diff = ms[["mean"]], any_call_diff_lo95 = ms[["lo95"]],
                any_call_diff_hi95 = ms[["hi95"]], discordant_q_only = bq, discordant_other_only = co,
                mcnemar_share_q_only_lo95 = if (bq + co) cpLower(bq, bq + co, 0.025) else NA_real_,
                mcnemar_share_q_only_hi95 = if (bq + co) cpUpper(bq, bq + co, 0.975) else NA_real_,
                large_excess = ex[["estimate"]], large_excess_x = ex[["x"]], large_excess_n = ex[["n"]],
                large_excess_upper1 = ex[["upper1"]], stringsAsFactors = FALSE))
        } else {
            tested <- okk & (q$planted_tested %in% TRUE)
            dd <- as.numeric(q$planted_called[tested] %in% TRUE) - as.numeric(o$planted_called[tested] %in% TRUE)
            ms <- meanStats(dd, 0.05 / N_POWER)
            mm <- meanStats(q$calls_material[okk] - o$calls_material[okk], 0.05 / N_POWER)
            row <- cbind(row, data.frame(power_diff = ms[["mean"]], power_diff_lo95 = ms[["lo95"]],
                power_diff_hi95 = ms[["hi95"]], power_diff_lowerB18 = ms[["lowerB"]], power_diff_lists = length(dd),
                material_terms_diff = mm[["mean"]], material_terms_diff_lo95 = mm[["lo95"]],
                material_terms_diff_hi95 = mm[["hi95"]], stringsAsFactors = FALSE))
        }
        out[[other]] <- row
    }
    bindFill(out)
})
paired <- bindFill(paired_rows)
paired <- paired[order(paired$cell, paired$propagation != "none", paired$comparison), ]
rownames(paired) <- NULL

# ---- mechanism (5.5 (g)) -------------------------------------------------------------------------------------------
mechanism <- bindFill(lapply(split(rows, list(rows$cell, rows$propagation), drop = TRUE), function(d) {
    q <- wide(d, "Q")
    qf <- wide(d, "Qf")
    q <- q[q$status == "ok", ]
    qf <- qf[qf$status == "ok", ]
    stages <- c("1", "2", "3", "4", "na", "all_one", "all_zero")
    stage_share <- vapply(stages, function(s) mean(q$chain_stage == s), numeric(1))
    stage_any <- vapply(stages, function(s) if (any(q$chain_stage == s)) mean(q$calls[q$chain_stage == s] > 0) else NA_real_, numeric(1))
    stage_large <- vapply(stages, function(s) if (any(q$chain_stage == s)) mean(q$calls[q$chain_stage == s] >= 10) else NA_real_, numeric(1))
    qfs <- vapply(stages, function(s) mean(qf$chain_stage == s), numeric(1))
    low <- !is.na(q$pi0_hat) & q$pi0_hat < 0.5
    matched <- match(q$replicate[low], qf$replicate)
    pq <- stats::quantile(q$pi0_hat, c(0, 0.1, 0.5, 0.9, 1), na.rm = TRUE, names = FALSE)
    data.frame(cell = d$cell[[1L]], propagation = d$propagation[[1L]], calls_ok = nrow(q),
        m_median = stats::median(q$m), m_kept_median = stats::median(q$m_kept),
        p_ge_095_median = stats::median(q$n_p_ge_095), p_ge_05_median = stats::median(q$n_p_ge_05),
        p_eq_1_median = stats::median(q$n_p_eq_1), q_all_na = sum(q$q_all_na),
        pi0_min = pq[[1L]], pi0_p10 = pq[[2L]], pi0_median = pq[[3L]], pi0_p90 = pq[[4L]], pi0_max = pq[[5L]],
        pi0_below_05 = sum(low), qf_pi0_median = stats::median(qf$pi0_hat, na.rm = TRUE),
        below_05_whole_family_ge_09 = sum(qf$pi0_hat[matched] >= 0.9, na.rm = TRUE),
        stage_share = paste(stages, sprintf("%.4f", stage_share), sep = ":", collapse = " "),
        any_call_by_stage = paste(stages, sprintf("%.4f", stage_any), sep = ":", collapse = " "),
        large_by_stage = paste(stages, sprintf("%.4f", stage_large), sep = ":", collapse = " "),
        qf_stage_share = paste(stages, sprintf("%.4f", qfs), sep = ":", collapse = " "),
        stringsAsFactors = FALSE)
}))
mechanism <- mechanism[order(mechanism$cell, mechanism$propagation != "none"), ]

# ---- outcome (5.7) and the readings (rule F, text T; VAL-86) -------------------------------------------------------
fails <- function(arm) any(classifications$classification[classifications$arm == arm] == "fails")
controls <- function(arm) all(classifications$classification[classifications$arm == arm] == "controls")
outcome <- if (fails("BH")) "O3" else if (fails("Q")) "O1" else if (controls("Q") && controls("BH") && controls("BHf")) "O2" else
    if (fails("BHf")) "O5" else "O4"
power_deciding <- paired[paired$comparison == "Q-BH" & paired$cell %in% DECIDING_PLANTED, ]
material <- any(power_deciding$power_diff_lowerB18 >= 0.05, na.rm = TRUE)
excess <- paired[paired$comparison == "Q-BH" & paired$cell %in% NULL_CELLS, ]
large_ok <- nrow(excess) == 18L && all(excess$large_excess_upper1 < MARGIN)
rule_f <- if (controls("BHf")) "BHf" else "BH"
role <- c(Q = "reported", BH = "reported", BHf = "reported")
branch <- ""
o2Reading <- function(bh) {
    if (!material) {
        role[[bh]] <<- "default"
        role[["Q"]] <<- "offered"
        if (identical(bh, "BHf")) role[["BH"]] <<- "offered"
        "no material power gain: BH by default"
    } else if (large_ok) {
        role[["Q"]] <<- "default"
        role[[bh]] <<- "offered"
        "material power gain and the large-event excess below 0.01 in every null setting: the q-value by default"
    } else {
        role[[bh]] <<- "default"
        role[["Q"]] <<- "offered"
        "material power gain but the large-event condition not met: BH by default"
    }
}
o4Reading <- function(bh) {
    if (controls("Q") && material && large_ok) {
        role[["Q"]] <<- "default"
        role[[bh]] <<- "offered"
        "the q-value controls everywhere with a material power gain and the large-event condition: the q-value by default"
    } else {
        role[[bh]] <<- "default"
        role[["Q"]] <<- "offered"
        "BH by default"
    }
}
if (identical(outcome, "O3")) {
    if (controls("BHf")) {
        role[] <- c("offered", "offered", "default")
        branch <- "BH fails; BH over the whole family controls everywhere: it is the default (a package change)"
    } else {
        role[] <- "stated (no default)"
        branch <- "BH fails and BH over the whole family does not control everywhere: no default; the owner decides on a further ticket"
    }
} else if (identical(outcome, "O1")) {
    role[[rule_f]] <- "default"
    role[["Q"]] <- "offered"
    branch <- sprintf("the q-value fails; BH by default over %s (rule F); the q-value offered with its measured failure",
        if (identical(rule_f, "BHf")) "the whole size-filtered family" else "the terms with at least one list protein")
} else if (identical(outcome, "O2")) {
    branch <- o2Reading(rule_f)
} else if (identical(outcome, "O5")) {
    role[["BHf"]] <- "not recommended"
    branch <- paste("BH over the whole family fails:", if (controls("Q") && controls("BH")) o2Reading("BH") else o4Reading("BH"))
} else {
    branch <- o4Reading(rule_f)
}
default_arm <- names(role)[role == "default"]

pct <- function(x) ifelse(is.na(x), "NA", sprintf("%.1f %%", 100 * x))
cellName <- function(cell, prop) {
    info <- CELLS[CELLS$cell == cell, ]
    what <- if (identical(info$scenario, "null")) sprintf("lists of %d", info$L) else
        sprintf("%s planted process, lists of %d", info$scenario, info$L)
    sprintf("%s, %s, %s", BG_LABEL[[info$background]], what, PROP_LABEL[[prop]])
}
textFor <- function(arm) {
    n <- classifications[classifications$arm == arm & classifications$cell %in% NULL_CELLS, ]
    pl <- classifications[classifications$arm == arm & classifications$cell %in% DECIDING_PLANTED, ]
    nc <- cells[cells$arm == arm & cells$cell %in% NULL_CELLS, ]
    pc <- cells[cells$arm == arm & cells$cell %in% DECIDING_PLANTED, ]
    failList <- function(x) if (!any(x$classification == "fails")) "" else sprintf(" (%s)", paste(sprintf("%s: lower limit %s",
        mapply(cellName, x$cell[x$classification == "fails"], x$propagation[x$classification == "fails"]),
        pct(x$lowerB[x$classification == "fails"])), collapse = "; "))
    big <- nc[which.max(nc$large_rate), ]
    txt <- sprintf(paste0("%s. Assumption: %s. On random protein lists (no true signal; 18 settings: E. coli full proteome ",
        "and a 2,000-protein subset, yeast; lists of 25, 126 and 400 proteins; with and without parent terms), %s to %s of ",
        "lists had at least one significant term (highest one-sided 95 %% upper limit %s; line 7.5 %%): it controlled false ",
        "calls in %d of 18 settings, failed in %d%s and was not determined in %d. Lists with 10 or more significant terms: ",
        "at most %s (%s). On lists of 126 with one planted process, mean false discovery proportion %s to %s (highest upper ",
        "limit %s; controlled in %d of 18 settings, failed in %d%s, not determined in %d), the planted term called in %s to ",
        "%s of lists where it was tested"),
        ARM_NAME[[arm]], ASSUMPTION[[arm]], pct(min(n$estimate)), pct(max(n$estimate)), pct(max(n$upper1)),
        sum(n$classification == "controls"), sum(n$classification == "fails"), failList(n),
        sum(!n$classification %in% c("controls", "fails")), pct(big$large_rate), cellName(big$cell, big$propagation),
        pct(min(pl$estimate)), pct(max(pl$estimate)), pct(max(pl$upper1)), sum(pl$classification == "controls"),
        sum(pl$classification == "fails"), failList(pl), sum(!pl$classification %in% c("controls", "fails")),
        pct(min(pc$power_planted, na.rm = TRUE)), pct(max(pc$power_planted, na.rm = TRUE)))
    if (identical(arm, "Q")) {
        txt <- paste0(txt, sprintf(paste0(", and against Benjamini-Hochberg over the same terms it called the planted term ",
            "at most %s more often (paired; highest one-sided lower limit over 18 settings %s; a material gain needs 5.0 %%)"),
            pct(max(power_deciding$power_diff, na.rm = TRUE)), pct(max(power_deciding$power_diff_lowerB18, na.rm = TRUE))))
    }
    paste0(txt, ".")
}
options <- data.frame(arm = DECIDING_ARMS, name = unname(ARM_NAME[DECIDING_ARMS]), role = unname(role[DECIDING_ARMS]),
    text = vapply(DECIDING_ARMS, textFor, character(1)), stringsAsFactors = FALSE)
outcome_table <- data.frame(outcome = outcome, reading = branch, default = if (length(default_arm)) default_arm else "none",
    rule_f = rule_f, material_power_gain = material, large_event_condition = large_ok,
    q_fails = fails("Q"), bh_fails = fails("BH"), bhf_fails = fails("BHf"), q_controls = controls("Q"),
    bh_controls = controls("BH"), bhf_controls = controls("BHf"),
    classifications = nrow(classifications), controls = sum(classifications$classification == "controls"),
    fails_n = sum(classifications$classification == "fails"),
    not_determined = sum(!classifications$classification %in% c("controls", "fails")), stringsAsFactors = FALSE)

# ---- write ---------------------------------------------------------------------------------------------------------
write <- function(x, name) utils::write.table(x, file.path(out_dir, name), sep = "\t", quote = FALSE, row.names = FALSE,
    na = "NA")
write(cells, "cells.tsv")
write(classifications, "classifications.tsv")
write(paired, "paired.tsv")
write(mechanism, "mechanism.tsv")
write(checks, "checks-score.tsv")
write(outcome_table, "outcome.tsv")
write(options, "options.tsv")
written <- file.path(out_dir, c("cells.tsv", "classifications.tsv", "paired.tsv", "mechanism.tsv", "checks-score.tsv",
    "outcome.tsv", "options.tsv"))
if (!is.null(blind_file)) {
    summary <- data.frame(file = basename(written), lines = vapply(written, function(f) length(readLines(f)), integer(1)),
        sha256 = vapply(written, function(f) digest::digest(file = f, algo = "sha256"), character(1)),
        stringsAsFactors = FALSE)
    summary <- rbind(summary, data.frame(file = "(input)", lines = nrow(rows), sha256 = sprintf("%d files", length(files))))
    utils::write.table(summary, blind_file, sep = "\t", quote = FALSE, row.names = FALSE)
    writeLines(sprintf("score: completed blind on %d rows from %d files; %d classifications", nrow(rows), length(files),
        nrow(classifications)))
} else {
    say("outcome %s: %s; default %s", outcome, branch, paste(default_arm, collapse = ","))
    print(outcome_table)
    print(checks)
}
