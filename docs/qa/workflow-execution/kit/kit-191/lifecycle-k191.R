# KIT-191 L1 (ticket 5.6, 11): the permutation test's time and memory, each case in a fresh process. A measurement of
# the instrument's cost, not a result about the method. Per case: the inputs built (G features x n biological samples,
# random values, a factor, `sets` sets with sizes in the stated range, seed 20261191), `VmRSS` read immediately before
# the call, the call timed (elapsed), `VmHWM` read after; the increase is VmHWM after minus that VmRSS. Then a declared
# invalid call (permutations = 0), recorded as refused. Stop rule: a case whose increase exceeds 1 GB stops the ticket.
# The coefficients of the predicted time a * n * G * B + c * G * B + s * sum(|set|) * B are each twice the seconds per
# unit of their own term in the case it dominates: a from case (c), c from case (d), s from case (b) (the larger of the
# case's two numbers of reorderings).
#
# Usage (from the installed head's package root): Rscript lifecycle-k191.R <out dir>
# Writes <out dir>/lifecycle.tsv, <out dir>/coefficients.tsv, <out dir>/refused.tsv.

arguments <- commandArgs(trailingOnly = TRUE)
stopifnot(length(arguments) == 1L)
out_dir <- arguments[[1L]]
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
out_dir <- normalizePath(out_dir)
rscript <- file.path(R.home("bin"), "Rscript")

cases <- data.frame(
    case = c("a", "b1", "b2", "c1", "c2", "d"),
    features = c(4000L, 4000L, 4000L, 20000L, 20000L, 20000L),
    samples = c(12L, 12L, 12L, 500L, 500L, 4L),
    sets = c(150L, 5000L, 5000L, 2000L, 2000L, 20L),
    smallest = c(5L, 5L, 5L, 5L, 5L, 5L),
    largest = c(50L, 500L, 500L, 500L, 500L, 5L),
    permutations = c(1000L, 1000L, 10000L, 1000L, 10000L, 10000L),
    stringsAsFactors = FALSE)

child <- file.path(tempdir(), "k191-l1-child.R")
writeLines(c(
    "a <- commandArgs(trailingOnly = TRUE)",
    "g <- as.integer(a[[1L]]); n <- as.integer(a[[2L]]); k <- as.integer(a[[3L]]); lo <- as.integer(a[[4L]])",
    "hi <- as.integer(a[[5L]]); B <- as.integer(a[[6L]]); out <- a[[7L]]",
    "status <- function(field) as.numeric(sub('[^0-9]*([0-9]+).*', '\\\\1', grep(paste0('^', field, ':'),",
    "    readLines('/proc/self/status'), value = TRUE)))",
    "suppressMessages(library(MultiScholaR))",
    "ns <- asNamespace('MultiScholaR')",
    "set.seed(20261191L)",
    "y <- matrix(stats::rnorm(g * n), g, n)",
    "ids <- sprintf('F%06d', seq_len(g)); rownames(y) <- ids",
    "f <- stats::rnorm(n)",
    "sizes <- if (lo == hi) rep(lo, k) else sample(lo:hi, k, replace = TRUE)",
    "members <- lapply(sizes, function(m) sample(ids, m))",
    "term2gene <- data.frame(term = rep(sprintf('S%05d', seq_len(k)), sizes), gene = unlist(members), stringsAsFactors = FALSE)",
    "term2name <- data.frame(term = sprintf('S%05d', seq_len(k)), name = sprintf('S%05d', seq_len(k)), stringsAsFactors = FALSE)",
    "values <- list(y = y, f = f, counts = list(universe = g, samples = n))",
    "settings <- list(min_size = 1L, max_size = 1e6, permutations = B, time_limit = Inf)",
    "invisible(gc())",
    "rss <- status('VmRSS')",
    "seconds <- system.time(result <- ns$integrationEnrichmentPermutation(values, term2gene, term2name, settings, 7L))[['elapsed']]",
    "hwm <- status('VmHWM')",
    "ok <- nrow(result$terms) == k && all(is.finite(result$terms$pvalue))",
    "write.table(data.frame(seconds = seconds, rss_before_kb = rss, hwm_after_kb = hwm, members = sum(sizes), terms_ok = ok),",
    "    out, sep = '\\t', quote = FALSE, row.names = FALSE)"), child)

rows <- list()
for (i in seq_len(nrow(cases))) {
    x <- cases[i, ]
    result <- file.path(tempdir(), paste0("l1-", x$case, ".tsv"))
    log <- file.path(out_dir, paste0("case-", x$case, ".log"))
    started <- Sys.time()
    code <- system2(rscript, c("--vanilla", child, x$features, x$samples, x$sets, x$smallest, x$largest, x$permutations,
        result), stdout = log, stderr = log)
    measured <- if (code == 0L && file.exists(result)) utils::read.delim(result) else
        data.frame(seconds = NA_real_, rss_before_kb = NA_real_, hwm_after_kb = NA_real_, members = NA_real_, terms_ok = FALSE)
    row <- cbind(x, measured, exit = code, wall_process_seconds = as.numeric(difftime(Sys.time(), started, units = "secs")))
    row$increase_mb <- (row$hwm_after_kb - row$rss_before_kb) / 1024
    rows[[i]] <- row
    cat(sprintf("case %s: exit %d, %.1f s, increase %.0f MB\n", x$case, code, row$seconds, row$increase_mb))
    if (isTRUE(row$increase_mb > 1024)) {
        cat("STOP: case", x$case, "increased the process's memory by more than 1 GB (stop rule, ticket 5.6)\n")
        break
    }
}
lifecycle <- do.call(rbind, rows)
utils::write.table(lifecycle, file.path(out_dir, "lifecycle.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

# The coefficients: twice the seconds per unit of the dominating term.
per <- function(cases_of, unit) {
    rows <- lifecycle[lifecycle$case %in% cases_of, , drop = FALSE]
    max(rows$seconds / unit(rows))
}
coefficients <- data.frame(
    coefficient = c("a", "c", "s"),
    from_cases = c("c1, c2", "d", "b1, b2"),
    seconds_per_unit = c(
        per(c("c1", "c2"), function(r) r$samples * r$features * r$permutations),
        per("d", function(r) r$features * r$permutations),
        per(c("b1", "b2"), function(r) r$members * r$permutations)))
coefficients$value <- signif(2 * coefficients$seconds_per_unit, 3)
utils::write.table(coefficients, file.path(out_dir, "coefficients.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

# The declared invalid call.
refused <- tryCatch({
    suppressMessages(library(MultiScholaR))
    y <- matrix(stats::rnorm(400), 40, 10, dimnames = list(sprintf("F%02d", 1:40), NULL))
    asNamespace("MultiScholaR")$integrationEnrichmentPermutation(list(y = y, f = stats::rnorm(10), counts = list()),
        data.frame(term = "S", gene = sprintf("F%02d", 1:5)), data.frame(term = "S", name = "S"),
        list(min_size = 1L, max_size = 100L, permutations = 0L), 1L)
    "not refused"
}, integration_enrichment_parameters_invalid = function(condition) paste("refused:", conditionMessage(condition)))
utils::write.table(data.frame(call = "permutations = 0", outcome = refused), file.path(out_dir, "refused.tsv"), sep = "\t",
    quote = FALSE, row.names = FALSE)
print(lifecycle)
print(coefficients)
cat(refused, "\n")
stopifnot(all(lifecycle$exit == 0L), all(lifecycle$terms_ok), nrow(lifecycle) == nrow(cases), startsWith(refused, "refused"))
