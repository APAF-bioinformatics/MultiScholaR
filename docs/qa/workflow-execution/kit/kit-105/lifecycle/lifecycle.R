# KIT-105 Simulation Setup: one entry lane alone, under strict skips, with counting traces on the new code paths.
# Usage: Rscript lifecycle.R <tree> <lane> <fn1,fn2,...>
args <- commandArgs(TRUE)
tree <- args[1]; lane <- args[2]; fns <- strsplit(args[3], ",", fixed = TRUE)[[1]]
Sys.setenv(MULTISCHOLAR_STRICT_SKIPS = "1")
suppressMessages(pkgload::load_all(tree, quiet = TRUE, helpers = TRUE))
setwd(file.path(tree, "tests", "testthat"))
ns <- asNamespace("MultiScholaR")
calls <- new.env()
for (fn in fns) {
  assign(fn, 0L, envir = calls)
  tracer <- bquote(assign(.(fn), get(.(fn), envir = .(calls)) + 1L, envir = .(calls)))
  suppressMessages(trace(fn, tracer = tracer, where = ns, print = FALSE))
}
started <- Sys.time()
reporter <- testthat::ListReporter$new()
# KIT-110 W2 (instrumentation repair of the kit-105 copy): the ListReporter needs a current file, else
# `add_result` fails with "attempt to use zero-length variable name"; the file is started and ended here.
multi <- testthat::MultiReporter$new(list(testthat::SummaryReporter$new(), reporter))
res <- testthat::with_reporter(multi, {
  multi$start_file(paste0("lifecycle-", lane))
  testthat::test_that(paste("lifecycle", lane), {
    guiEntryExpectLane(lane)
  })
  multi$end_file()
})
r <- as.data.frame(reporter$get_results())
ok <- nrow(r) > 0 && sum(r$failed) == 0 && !any(r$error) && sum(r$skipped) == 0
# (KIT-105 lead: the first version read isTRUE(<reporter object>), always FALSE; see entry-suite.md)
cat(sprintf("LANE %s identical_to_oracle=%s elapsed=%.0fs\n", lane, ok,
  as.numeric(difftime(Sys.time(), started, units = "secs"))))
for (fn in fns) cat(sprintf("CALLS %s %d\n", fn, get(fn, envir = calls)))
