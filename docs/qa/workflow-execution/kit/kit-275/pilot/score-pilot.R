# KIT-275 revision 2.1 PILOT (ticket R2.3): the report and the rule of pilot-k275.R. NON-COUNTED, unblinded. Base R and stats only.
# Adapted from diag/score-diag.R (the table reading, Clopper-Pearson, hinges), with the rule lines and the outcome added.
# Usage: Rscript score-pilot.R <rows dir> <out dir>
#   <rows dir> holds pilot-rows.tsv (the finished job) or, when the job was killed, part-*.tsv; the assembled file wins.
# Per design and anomaly (V1a, V1b): the data sets with a valid journey (status ok and a finite planted correlation), how many have the planted
# sample's best within-group correlation below 0.5 with the two-sided 95 % Clopper-Pearson interval, the five-number summary (Tukey's hinges) of the
# planted sample's and of the fixed genuine sample's best correlation, the genuine sample's count below 0.5, the not-ok journeys.
# The rule (R2.3, fixed before the pilot): V1a is accepted iff it is failed (below 0.5) in at least 45 of 50 data sets in EVERY design; if not,
# V1b by the same rule; if neither, the matter returns to the owner. A cell with fewer than 50 valid journeys is an instrument defect (audit P-1):
# it is never a rejection, the rule line says "incomplete" and the outcome is "no decision" (re-send once under VAL-54), except that a complete,
# accepted V1a decides whatever V1b holds. Writes pilot-summary.tsv and pilot-decision.tsv and prints them.
options(warnPartialMatchDollar = TRUE)
globalCallingHandlers(warning = function(w) {
    if (grepl("partial match of", tolower(conditionMessage(w)), fixed = TRUE)) stop(paste("partial-match guard:", conditionMessage(w)), call. = FALSE)
})
args <- commandArgs(TRUE)
stopifnot(length(args) == 2L)
rows_dir <- args[[1L]]
out_dir <- args[[2L]]
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
DESIGNS <- c("D3", "D4", "D6")
ANOMALIES <- c("V1a", "V1b")
EXPECTED <- 50L
ACCEPT_AT <- 45L
THETA <- 0.5
readText <- function(path) utils::read.delim(path, stringsAsFactors = FALSE, colClasses = "character", na.strings = "NA")
assembled <- file.path(rows_dir, "pilot-rows.tsv")
files <- if (file.exists(assembled)) assembled else list.files(rows_dir, pattern = "^part-.*\\.tsv$", full.names = TRUE)
if (!length(files)) stop("score-pilot: no rows in ", rows_dir)
rows <- do.call(rbind, lapply(files, readText))
need <- c("design", "t", "anomaly", "status", "cy_planted", "cy_genuine")
if (!all(need %in% names(rows))) stop("score-pilot: columns missing: ", paste(setdiff(need, names(rows)), collapse = ","))
num <- function(x) suppressWarnings(as.numeric(x))
cpInterval <- function(x, n, alpha = 0.05) {
    if (n < 1L) return(c(NA_real_, NA_real_))
    c(if (x == 0L) 0 else stats::qbeta(alpha / 2, x, n - x + 1L), if (x == n) 1 else stats::qbeta(1 - alpha / 2, x + 1L, n - x))
}
f6 <- function(x) if (is.na(x)) "NA" else sprintf("%.6f", x)
summary_rows <- list()
for (anomaly in ANOMALIES) for (design in DESIGNS) {
    x <- rows[rows[["design"]] == design & rows[["anomaly"]] == anomaly, , drop = FALSE]
    planted <- num(x[["cy_planted"]])
    genuine <- num(x[["cy_genuine"]])
    valid <- x[["status"]] == "ok" & is.finite(planted)
    pv <- planted[valid]
    gv <- genuine[valid & is.finite(genuine)]
    n_valid <- length(pv)
    n_below <- sum(pv < THETA)
    ci <- cpInterval(n_below, n_valid)
    fp <- if (length(pv)) stats::fivenum(pv) else rep(NA_real_, 5L)
    fg <- if (length(gv)) stats::fivenum(gv) else rep(NA_real_, 5L)
    summary_rows[[length(summary_rows) + 1L]] <- data.frame(anomaly = anomaly, design = design, n_rows = nrow(x), n_valid = n_valid,
        n_not_ok = sum(x[["status"]] != "ok"), coverage = if (n_valid == EXPECTED) "complete" else sprintf("%d of %d", n_valid, EXPECTED),
        planted_below_05 = n_below, planted_below_share = f6(if (n_valid) n_below / n_valid else NA_real_), cp_lower = f6(ci[[1L]]), cp_upper = f6(ci[[2L]]),
        planted_min = f6(fp[[1L]]), planted_hinge1 = f6(fp[[2L]]), planted_median = f6(fp[[3L]]), planted_hinge3 = f6(fp[[4L]]), planted_max = f6(fp[[5L]]),
        genuine_n = length(gv), genuine_below_05 = sum(gv < THETA),
        genuine_min = f6(fg[[1L]]), genuine_hinge1 = f6(fg[[2L]]), genuine_median = f6(fg[[3L]]), genuine_hinge3 = f6(fg[[4L]]), genuine_max = f6(fg[[5L]]),
        stringsAsFactors = FALSE)
}
summary_table <- do.call(rbind, summary_rows)
utils::write.table(summary_table, file.path(out_dir, "pilot-summary.tsv"), sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
# the rule, per anomaly: accepted / rejected / incomplete (audit P-1)
ruleOf <- function(anomaly) {
    s <- summary_table[summary_table[["anomaly"]] == anomaly, , drop = FALSE]
    s <- s[match(DESIGNS, s[["design"]]), , drop = FALSE]
    complete <- all(s[["n_valid"]] == EXPECTED)
    verdict <- if (!complete) "incomplete" else if (all(s[["planted_below_05"]] >= ACCEPT_AT)) "accepted" else "rejected"
    line <- sprintf("%s: %s -> accepted iff all >= %d: %s%s", anomaly, paste(sprintf("%s %d/%d", s[["design"]], s[["planted_below_05"]], s[["n_valid"]]), collapse = ", "),
        ACCEPT_AT, verdict, if (complete) "" else " (a cell has fewer than 50 valid journeys: an instrument defect, not a rejection)")
    list(verdict = verdict, line = line)
}
rule_a <- ruleOf("V1a")
rule_b <- ruleOf("V1b")
outcome <- if (identical(rule_a[["verdict"]], "accepted")) "V1a accepted" else if (identical(rule_a[["verdict"]], "incomplete")) "no decision: V1a incomplete (instrument defect; re-send once under VAL-54)" else
    if (identical(rule_b[["verdict"]], "accepted")) "V1b accepted (V1a rejected)" else if (identical(rule_b[["verdict"]], "incomplete")) "no decision: V1b incomplete (instrument defect; re-send once under VAL-54)" else
    "neither accepted: back to the owner (R2.3 (3)); no other candidate, no further pilot"
decision <- data.frame(item = c("rule_V1a", "rule_V1b", "outcome"), value = c(rule_a[["line"]], rule_b[["line"]], outcome), stringsAsFactors = FALSE)
utils::write.table(decision, file.path(out_dir, "pilot-decision.tsv"), sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
unit_rows <- sum(rows[["anomaly"]] == "unit")
cat(sprintf("score-pilot: %d rows from %d file(s); %d unit rows; planted sample below 0.5 (valid data sets, two-sided 95%% Clopper-Pearson):\n", nrow(rows), length(files), unit_rows))
for (i in seq_len(nrow(summary_table))) {
    s <- summary_table[i, ]
    cat(sprintf("  %s %s valid %2d of %d  below 0.5: %2d (%s; %s to %s)  planted min/h1/med/h3/max %s %s %s %s %s  genuine %s %s %s %s %s  genuine below 0.5: %d of %d  not ok %d\n",
        s[["anomaly"]], s[["design"]], s[["n_valid"]], EXPECTED, s[["planted_below_05"]], s[["planted_below_share"]], s[["cp_lower"]], s[["cp_upper"]],
        s[["planted_min"]], s[["planted_hinge1"]], s[["planted_median"]], s[["planted_hinge3"]], s[["planted_max"]],
        s[["genuine_min"]], s[["genuine_hinge1"]], s[["genuine_median"]], s[["genuine_hinge3"]], s[["genuine_max"]], s[["genuine_below_05"]], s[["genuine_n"]], s[["n_not_ok"]]))
}
cat(rule_a[["line"]], "\n", sep = "")
cat(rule_b[["line"]], "\n", sep = "")
cat("outcome (R2.3): ", outcome, "\n", sep = "")
invisible(NULL)
