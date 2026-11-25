# KIT-222 audit-repair run: the head's test-kit222-replicate-keys.R on the base commit (failing-first evidence for the
# rows added after the tests-only commit). A detached worktree of <base> is loaded with pkgload; the head's test file
# is copied into it and run there. Writes <out>/base-tests.tsv; the rows' "KIT222-OBS" lines go to the console.
# Usage: Rscript base-tests.R <base commit> <out directory>
args <- commandArgs(trailingOnly = TRUE)
base <- args[[1L]]
out <- args[[2L]]
dir.create(out, recursive = TRUE, showWarnings = FALSE)
worktree <- file.path(tempdir(), "kit222-base")
stopifnot(system2("git", c("worktree", "add", "--detach", worktree, base)) == 0L)
file.copy("tests/testthat/test-kit222-replicate-keys.R", file.path(worktree, "tests", "testthat"), overwrite = TRUE)
pkgload::load_all(worktree, export_all = TRUE, helpers = TRUE, quiet = TRUE)
results <- testthat::test_file(file.path(worktree, "tests", "testthat", "test-kit222-replicate-keys.R"),
    reporter = testthat::ListReporter$new(), stop_on_failure = FALSE)
table <- as.data.frame(results)
table <- table[, intersect(c("test", "nb", "failed", "skipped", "error", "warning"), names(table))]
utils::write.table(table, file.path(out, "base-tests.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
print(table)
