# one.R <pkg> <test file> <desc>: one test_that block, its counts and its failures.
args <- commandArgs(TRUE)
suppressMessages(pkgload::load_all(args[1], quiet = TRUE, export_all = TRUE, helpers = TRUE))
res <- testthat::test_file(args[2], reporter = "silent", stop_on_failure = FALSE, desc = args[3])
df <- as.data.frame(res)
cat(sprintf("COUNTS %s [%s] pass %d fail %d error %d skip %d\n", basename(args[2]), args[3], sum(df$passed), sum(df$failed), sum(df$error), sum(df$skipped)))
for (r in res) for (x in r$results) if (!inherits(x, "expectation_success")) cat("---", class(x)[[1]], substr(conditionMessage(x), 1, 400), "\n")
