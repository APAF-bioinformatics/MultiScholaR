# KIT-184 (ticket 8.6): the cost of deriving the removed factors in the shared writer, on a synthetic pair of 20,000
# features by 500 runs (random log2 values; the corrected matrix the normalised one less two directions over the runs
# times random coefficients). Each case in a fresh process: wall time of the function, the process's resident memory
# before the call (VmRSS, after the two tables are built) and its peak after it (VmHWM). Ceilings: 60 seconds, and 1 GiB
# of peak above the reading before the call. A case passing a ceiling stops with exit status 1 after the table is
# written (the ticket then returns the figures; the QR fallback is not built otherwise).
#
# Usage (from the package root, MultiScholaR installed): Rscript measure-removed-factors.R <out dir>
# Writes <out dir>/removed-factors-cost.tsv (case, features, runs, k, seconds, rss_before_mib, hwm_after_mib,
# increase_mib, status, rank, share).

arguments <- commandArgs(trailingOnly = TRUE)
stopifnot(length(arguments) == 1L)
out_dir <- arguments[[1L]]
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
out_dir <- normalizePath(out_dir)
rscript <- file.path(R.home("bin"), "Rscript")
child <- '
arguments <- commandArgs(trailingOnly = TRUE)
features <- as.integer(arguments[[1L]]); runs <- as.integer(arguments[[2L]]); k <- as.integer(arguments[[3L]])
status <- function(field) {
    line <- grep(paste0("^", field, ":"), readLines("/proc/self/status"), value = TRUE)
    as.numeric(gsub("[^0-9]", "", line)) / 1024
}
set.seed(20261848L)
ids <- sprintf("r%04d", seq_len(runs))
normalised <- matrix(stats::rnorm(features * runs, 20, 2), features, runs, dimnames = list(NULL, ids))
w <- matrix(stats::rnorm(runs * k), runs, k)
corrected <- normalised - matrix(stats::rnorm(features * k), features, k) %*% t(w) - stats::runif(features)
table <- function(m) data.frame(feature_id = sprintf("F%05d", seq_len(features)), m, check.names = FALSE,
    stringsAsFactors = FALSE)
a <- table(normalised); b <- table(corrected)
rm(normalised, corrected); invisible(gc())
before <- status("VmRSS")
seconds <- system.time(out <- MultiScholaR:::integrationInputRemovedFactors(a, b, k, ids))[["elapsed"]]
after <- status("VmHWM")
cat(sprintf("%d\\t%d\\t%d\\t%.3f\\t%.1f\\t%.1f\\t%.1f\\t%s\\t%s\\t%.12g\\n", features, runs, k, seconds, before, after,
    after - before, out$record$status, format(out$record$rank), out$record$share))
'
script <- file.path(tempdir(), "k184-measure-child.R")
writeLines(child, script)
cases <- list(c(case = "small", features = 2000L, runs = 60L, k = 2L), c(case = "ceiling", features = 20000L, runs = 500L, k = 2L),
    c(case = "ceiling_k4", features = 20000L, runs = 500L, k = 4L))
rows <- lapply(cases, function(case) {
    line <- system2(rscript, c("--vanilla", shQuote(script), case[["features"]], case[["runs"]], case[["k"]]), stdout = TRUE)
    fields <- strsplit(line[[length(line)]], "\t", fixed = TRUE)[[1L]]
    stopifnot(length(fields) == 10L)
    data.frame(case = case[["case"]], features = as.integer(fields[[1L]]), runs = as.integer(fields[[2L]]),
        k = as.integer(fields[[3L]]), seconds = as.numeric(fields[[4L]]), rss_before_mib = as.numeric(fields[[5L]]),
        hwm_after_mib = as.numeric(fields[[6L]]), increase_mib = as.numeric(fields[[7L]]), status = fields[[8L]],
        rank = fields[[9L]], share = as.numeric(fields[[10L]]), stringsAsFactors = FALSE)
})
table <- do.call(rbind, rows)
utils::write.table(table, file.path(out_dir, "removed-factors-cost.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
print(table)
ceiling_rows <- table[table$features == 20000L, , drop = FALSE]
within <- all(ceiling_rows$seconds <= 60) && all(ceiling_rows$increase_mib <= 1024)
writeLines(sprintf("within the ceilings of ticket 8.6 (60 s; 1 GiB above the reading before the call): %s", within),
    file.path(out_dir, "removed-factors-cost-verdict.txt"))
if (!within) quit(status = 1L)
