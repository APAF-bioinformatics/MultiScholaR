# KIT-069 Simulation Setup: run one module-CI test file from a source tree and hash every file the
# test process wrote under its tempdir() before the process exits (the lanes' projects live there).
# usage: Rscript mci-hash.R <tree> <test file> <out.hashes> [<keep dir>: a copy of the temp tree for the content comparison]
args <- commandArgs(TRUE)
tree <- args[1]; test_file <- args[2]; out <- args[3]; keep <- if (length(args) >= 4) args[4] else NULL
suppressMessages(pkgload::load_all(tree, quiet = TRUE, export_all = TRUE))
res <- testthat::test_file(test_file, reporter = "summary", stop_on_failure = FALSE)
df <- as.data.frame(res)
cat(sprintf("\nCOUNTS %s pass %d fail %d error %d skip %d\n", basename(test_file),
  sum(df$passed), sum(df$failed), sum(df$error), sum(df$skipped)))
root <- tempdir()
files <- list.files(root, recursive = TRUE, all.files = TRUE, full.names = TRUE)
files <- files[!dir.exists(files)]
rel <- substring(files, nchar(root) + 2L)
hashes <- vapply(files, function(f) digest::digest(f, algo = "sha256", file = TRUE), character(1))
writeLines(paste(hashes, rel, sep = "  "), out)
cat(sprintf("HASHED %d files -> %s\n", length(files), out))
if (!is.null(keep)) { dir.create(keep, recursive = TRUE, showWarnings = FALSE); file.copy(list.files(root, full.names = TRUE, all.files = TRUE, no.. = TRUE), keep, recursive = TRUE); cat("KEPT ->", keep, "\n") }
