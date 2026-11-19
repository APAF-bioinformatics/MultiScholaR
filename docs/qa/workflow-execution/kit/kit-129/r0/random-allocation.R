# KIT-129 R0 part B, supplement (decision M37): the scoring of the random-allocation run, fixed before the run.
# Descriptive: no pass mark, no threshold set or changed, no verdict. Base R only (stats, utils). Runs on a hosted runner
# (.github/workflows/k129-score.yml), never on the laptop.
#
# Usage: Rscript --vanilla random-allocation.R <purpose: lifecycle | score> <in dir> <label> <cells> <first run> <runs> <out dir>
#   <in dir> holds the run's artifacts as `gh run download` writes them: <in dir>/k124-<label>-s<i>-<cell>/ with
#   s<i>-runs.tsv (the study-mode table) and s<i>-assoc.tsv (the random-allocation cells' statistic), i = 1, 3, 4.
#   <cells> a comma list (score: C0R,C4R,S1R); every cell of every simulator must hold exactly runs <first run> to
#   <first run> + <runs> - 1, or the script stops (a set is complete as declared).
#   purpose `score` is the declared reading: it refuses when the results are already committed in the checkout (read once).
#
# Thresholds: the validated ones (decisions M24, M27), used as they are; fires = the largest factor association strictly
# above the threshold; a row without a W cannot fire and stays in the denominator; crashed rows are left out and counted.
THRESHOLDS <- data.frame(
    simulator = c("S1", "S1", "S1", "S3", "S4", "S4"),
    arm = c("5", "7", "8", "7", "7", "8"),
    threshold = c(0.16376377, 0.058077853, 0.070068989, 0.038965628, 0.088611623, 0.19494193),
    set_at = c("5% point", "5% point", "5% point", "10% point", "10% point", "5% point"),
    source = c("threshold-followup-S1.tsv", "threshold-followup-S1.tsv", "threshold-followup-S1.tsv",
        "threshold-by-design-S3.tsv", "threshold-by-design-S4.tsv", "threshold-followup-S4.tsv"),
    stringsAsFactors = FALSE)
COMPARISON_ARMS <- c("1", "2")     # no correction; the ceiling (the true factor in the model)
POOLS <- list(`C0R+C4R` = c("C0R", "C4R"))   # the calibration's own unrelated designs with a factor, together
RECORD <- "docs/qa/workflow-execution/kit/kit-129/r0"
#
# Quantities (per simulator, cell or pool, and design):
#   association table: runs; rows (S4: one per assay) and crashed rows; the TRUE in-sample association of the factor
#     (one value per run): quartiles, median, 90th percentile, share of runs above the threshold; the ESTIMATED largest
#     association over the rows with a W: the same quantiles; rows with a W; share of the non-crashed rows above the
#     threshold (the warning's firing share, as the calibration counted it; S4 over both assays).
#   errors table: for the design and for arms 1 and 2 on the same runs, on the `full` subset, rows that are not crashed
#     (refused rows carry the run's no-RUV measures, as aggregate.R scores them): runs, mean null_p05 (type I at 0.05),
#     mean null_p01 (type I at 0.01), mean fdp (mean false discovery proportion at the 5% target); overall, split by
#     whether the run exceeded the design's threshold (S4: in either assay; runs without the design's statistic are in
#     `all` only), and by thirds of the true association (cut at its 1/3 and 2/3 quantiles within the cell or pool).
#   Quantiles are R's default (type 7).

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 7L) stop("usage: Rscript --vanilla random-allocation.R <lifecycle|score> <in dir> <label> <cells> <first run> <runs> <out dir>")
purpose <- args[[1L]]; in_dir <- args[[2L]]; label <- args[[3L]]
cells <- strsplit(args[[4L]], ",", fixed = TRUE)[[1L]]
runs_wanted <- as.integer(args[[5L]]) - 1L + seq_len(as.integer(args[[6L]]))
out_dir <- args[[7L]]
stopifnot(purpose %in% c("lifecycle", "score"), length(cells) >= 1L, length(runs_wanted) >= 1L)
if (purpose == "score") {
    done <- file.path(RECORD, c("random-allocation-results.md", "random-allocation-association.tsv", "random-allocation-errors.tsv"))
    if (any(file.exists(done))) stop("the random-allocation results are committed already: the set is read once (one attempt)")
}
if (dir.exists(out_dir) && length(list.files(out_dir))) stop("the out dir is not empty")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

readTsv <- function(path) utils::read.delim(path, colClasses = "character", quote = "", comment.char = "", na.strings = "NA")
num <- function(x) suppressWarnings(as.numeric(x))
problems <- character()
tests <- list(); stats_rows <- list()
for (sim in c("S1", "S3", "S4")) {
    s <- tolower(sim)
    for (cell in cells) {
        dir <- file.path(in_dir, sprintf("k124-%s-%s-%s", label, s, cell))
        runs_path <- file.path(dir, sprintf("%s-runs.tsv", s)); side_path <- file.path(dir, sprintf("%s-assoc.tsv", s))
        if (!file.exists(runs_path) || !file.exists(side_path)) {
            problems <- c(problems, sprintf("%s %s: %s or %s missing", sim, cell, basename(runs_path), basename(side_path)))
            next
        }
        r <- readTsv(runs_path); a <- readTsv(side_path)
        for (part in list(list(name = "runs", x = r), list(name = "assoc", x = a))) {
            x <- part$x
            if (!nrow(x) || !all(x$cell == cell)) problems <- c(problems, sprintf("%s %s %s: not the cell alone", sim, cell, part$name))
            if (!identical(sort(unique(as.integer(x$run))), runs_wanted))
                problems <- c(problems, sprintf("%s %s %s: runs are not %d to %d", sim, cell, part$name, min(runs_wanted), max(runs_wanted)))
        }
        tests[[length(tests) + 1L]] <- r
        stats_rows[[length(stats_rows) + 1L]] <- a
    }
}
if (length(problems)) stop(paste(c("the set is not complete as declared:", problems), collapse = "\n"))
tests <- do.call(rbind, tests); side <- do.call(rbind, stats_rows)
tests <- tests[tests$subset == "full", , drop = FALSE]

tally <- function(x, name) {
    out <- stats::aggregate(list(rows = rep(1L, nrow(x))), x[c("simulator", "cell", "arm", "status")], sum)
    out$table <- name
    out[c("simulator", "table", "cell", "arm", "status", "rows")]
}
status <- rbind(tally(tests, "runs (full subset)"), tally(side, "association"))

q <- function(x, p) if (length(x)) unname(stats::quantile(x, p, type = 7L)) else NA_real_
r4 <- function(x) signif(x, 4L)
groups_of <- function(sim) {
    present <- cells
    g <- stats::setNames(lapply(present, function(c) c), present)
    for (name in names(POOLS)) if (all(POOLS[[name]] %in% present)) g[[name]] <- POOLS[[name]]
    g
}
association <- list(); errors <- list()
for (k in seq_len(nrow(THRESHOLDS))) {
    sim <- THRESHOLDS$simulator[[k]]; design <- THRESHOLDS$arm[[k]]; thr <- THRESHOLDS$threshold[[k]]
    for (gname in names(groups_of(sim))) {
        members <- groups_of(sim)[[gname]]
        a <- side[side$simulator == sim & side$cell %in% members, , drop = FALSE]
        a$key <- paste(a$cell, a$run)
        truth <- a[!is.na(num(a$true_association)) & !duplicated(a$key), c("key", "true_association")]
        true_value <- stats::setNames(num(truth$true_association), truth$key)
        d <- a[a$arm == design, , drop = FALSE]
        alive <- d$status != "crashed"
        has_w <- d$status == "ok"
        est <- num(d$group_association[has_w])
        above_row <- has_w & num(d$group_association) > thr
        above_row[is.na(above_row)] <- FALSE
        exceeded <- tapply(above_row[alive], d$key[alive], any)
        association[[length(association) + 1L]] <- data.frame(simulator = sim, cell = gname, design = design,
            threshold = thr, threshold_set_at = THRESHOLDS$set_at[[k]], runs = length(true_value),
            rows = nrow(d), crashed_rows = sum(!alive), rows_with_W = sum(has_w),
            true_q25 = r4(q(true_value, 0.25)), true_median = r4(q(true_value, 0.5)), true_q75 = r4(q(true_value, 0.75)),
            true_p90 = r4(q(true_value, 0.9)), true_share_above_threshold = r4(mean(true_value > thr)),
            est_q25 = r4(q(est, 0.25)), est_median = r4(q(est, 0.5)), est_q75 = r4(q(est, 0.75)), est_p90 = r4(q(est, 0.9)),
            rows_above = sum(above_row[alive]), share_above = r4(sum(above_row[alive]) / sum(alive)),
            runs_exceeded = sum(exceeded), stringsAsFactors = FALSE)
        cuts <- q(true_value, c(1 / 3, 2 / 3))
        third <- stats::setNames(as.integer(cut(true_value, c(-Inf, cuts, Inf), right = TRUE)), names(true_value))
        t <- tests[tests$simulator == sim & tests$cell %in% members & tests$status != "crashed", , drop = FALSE]
        t$key <- paste(t$cell, t$run)
        t$exceeded <- unname(exceeded[t$key])
        t$third <- unname(third[t$key])
        splits <- list(all = rep(TRUE, nrow(t)), above = t$exceeded %in% TRUE, not_above = t$exceeded %in% FALSE,
            third_1 = t$third %in% 1L, third_2 = t$third %in% 2L, third_3 = t$third %in% 3L)
        for (split in names(splits)) for (arm in c(design, COMPARISON_ARMS)) {
            u <- t[splits[[split]] & t$arm == arm, , drop = FALSE]
            errors[[length(errors) + 1L]] <- data.frame(simulator = sim, cell = gname, design = design, split = split,
                split_true_range = switch(split, third_1 = sprintf("<= %.4g", cuts[[1L]]),
                    third_2 = sprintf("%.4g to %.4g", cuts[[1L]], cuts[[2L]]), third_3 = sprintf("> %.4g", cuts[[2L]]), ""),
                arm = arm, runs = nrow(u), type1_05 = r4(mean(num(u$null_p05), na.rm = TRUE)),
                type1_01 = r4(mean(num(u$null_p01), na.rm = TRUE)), mean_fdp = r4(mean(num(u$fdp), na.rm = TRUE)),
                stringsAsFactors = FALSE)
        }
    }
}
association <- do.call(rbind, association); errors <- do.call(rbind, errors)
errors[] <- lapply(errors, function(x) if (is.numeric(x)) ifelse(is.nan(x), NA, x) else x)
write <- function(x, name) utils::write.table(x, file.path(out_dir, name), sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
write(association, "random-allocation-association.tsv")
write(errors, "random-allocation-errors.tsv")
write(status, "random-allocation-status.tsv")
writeLines(c(paste("purpose:", purpose), paste("label:", label), paste("cells:", paste(cells, collapse = ",")),
    sprintf("runs: %d to %d", min(runs_wanted), max(runs_wanted)), paste("R:", R.version.string),
    paste("finished:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"))), file.path(out_dir, "random-allocation.meta.txt"))
cat("written:", paste(list.files(out_dir), collapse = ", "), "\n")
