# KIT-216 R1 (section 4): the base-made archives. Runs the new test file's row "R1 archives" on this commit with
# KIT216_ARCHIVE_OUT set, so that the latest archive the Summary writes after one and two enrichment runs (through the
# application's commit, then `refreshDIANativeSummaryArchive`) is kept in this script's output folder.
out <- file.path("docs", "qa", "workflow-execution", "kit", "kit-216", "r1-archives")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
Sys.setenv(KIT216_ARCHIVE_OUT = normalizePath(out))
suppressMessages(pkgload::load_all(".", quiet = TRUE, export_all = TRUE, helpers = TRUE))
testthat::test_file(file.path("tests", "testthat", "test-kit216-prot-session-after-enrichment.R"),
    desc = "KIT-216 R1 archives (section 4): the latest archive the Summary writes after one and two enrichment runs",
    reporter = "summary", stop_on_failure = FALSE)
files <- list.files(out, full.names = TRUE)
cat(sprintf("KIT-216 R1 archive file %s bytes=%d sha256=%s\n", basename(files), file.size(files),
    vapply(files, function(f) digest::digest(file = f, algo = "sha256"), "")), sep = "")
if (length(files) != 2L) quit(save = "no", status = 1L)
