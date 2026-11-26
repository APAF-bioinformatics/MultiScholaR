# KIT-275 diagnostic (VAL-392 (1)): the report of diag-k275.R. NON-COUNTED; descriptive; no label, no stop rule, nothing is decided by it.
# Base R and stats only. Usage: Rscript score-diag.R <rows dir> <out dir>
#   <rows dir> holds diag-rows.tsv (the finished job) or, when the job was killed, part-*.tsv (the unit files); the assembled file wins.
# Per design and normalisation: the data sets with a valid journey (status ok and a finite planted correlation), how many have the planted
# sample's best within-group correlation below 0.5 with the two-sided 95 % Clopper-Pearson interval, the five-number summary (Tukey's hinges,
# stats::fivenum) of the planted sample's best correlation and of the genuine sample's, and the not-ok journeys. Tables are read as TEXT.
# Writes diag-summary.tsv and prints it.
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
NORMS <- c("cyclicloess", "quantile", "none")
EXPECTED <- 50L
THETA <- 0.5
readText <- function(path) utils::read.delim(path, stringsAsFactors = FALSE, colClasses = "character", na.strings = "NA")
assembled <- file.path(rows_dir, "diag-rows.tsv")
files <- if (file.exists(assembled)) assembled else list.files(rows_dir, pattern = "^part-.*\\.tsv$", full.names = TRUE)
if (!length(files)) stop("score-diag: no rows in ", rows_dir)
rows <- do.call(rbind, lapply(files, readText))
need <- c("design", "t", "norm", "status", "cy_planted", "cy_genuine")
if (!all(need %in% names(rows))) stop("score-diag: columns missing: ", paste(setdiff(need, names(rows)), collapse = ","))
num <- function(x) suppressWarnings(as.numeric(x))
cpInterval <- function(x, n, alpha = 0.05) {
    if (n < 1L) return(c(NA_real_, NA_real_))
    c(if (x == 0L) 0 else stats::qbeta(alpha / 2, x, n - x + 1L), if (x == n) 1 else stats::qbeta(1 - alpha / 2, x + 1L, n - x))
}
f6 <- function(x) if (is.na(x)) "NA" else sprintf("%.6f", x)
summary_rows <- list()
for (design in DESIGNS) for (norm in NORMS) {
    x <- rows[rows[["design"]] == design & rows[["norm"]] == norm, , drop = FALSE]
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
    summary_rows[[length(summary_rows) + 1L]] <- data.frame(design = design, norm = norm, n_rows = nrow(x), n_valid = n_valid,
        n_not_ok = sum(x[["status"]] != "ok"), coverage = if (n_valid == EXPECTED) "complete" else sprintf("%d of %d", n_valid, EXPECTED),
        planted_below_05 = n_below, planted_below_share = f6(if (n_valid) n_below / n_valid else NA_real_), cp_lower = f6(ci[[1L]]), cp_upper = f6(ci[[2L]]),
        planted_min = f6(fp[[1L]]), planted_hinge1 = f6(fp[[2L]]), planted_median = f6(fp[[3L]]), planted_hinge3 = f6(fp[[4L]]), planted_max = f6(fp[[5L]]),
        genuine_n = length(gv), genuine_below_05 = sum(gv < THETA),
        genuine_min = f6(fg[[1L]]), genuine_hinge1 = f6(fg[[2L]]), genuine_median = f6(fg[[3L]]), genuine_hinge3 = f6(fg[[4L]]), genuine_max = f6(fg[[5L]]),
        stringsAsFactors = FALSE)
}
summary_table <- do.call(rbind, summary_rows)
utils::write.table(summary_table, file.path(out_dir, "diag-summary.tsv"), sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
unit_rows <- sum(rows[["norm"]] == "unit")
cat(sprintf("score-diag: %d rows from %d file(s); %d unit rows; cells below, planted sample (of valid data sets, two-sided 95%% Clopper-Pearson):\n",
    nrow(rows), length(files), unit_rows))
for (i in seq_len(nrow(summary_table))) {
    s <- summary_table[i, ]
    cat(sprintf("  %s %-12s valid %2d of %d  below 0.5: %2d (%s; %s to %s)  planted min/h1/med/h3/max %s %s %s %s %s  genuine %s %s %s %s %s  not ok %d\n",
        s[["design"]], s[["norm"]], s[["n_valid"]], EXPECTED, s[["planted_below_05"]], s[["planted_below_share"]], s[["cp_lower"]], s[["cp_upper"]],
        s[["planted_min"]], s[["planted_hinge1"]], s[["planted_median"]], s[["planted_hinge3"]], s[["planted_max"]],
        s[["genuine_min"]], s[["genuine_hinge1"]], s[["genuine_median"]], s[["genuine_hinge3"]], s[["genuine_max"]], s[["n_not_ok"]]))
}
invisible(NULL)
