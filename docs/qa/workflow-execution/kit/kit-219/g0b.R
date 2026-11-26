# KIT-219 G0 (iii) and (iv) (ticket section 10; coordinator's ruling 3 of 2026-11-27). Run from the package root on the
# fix commit, after it is installed:  Rscript docs/qa/workflow-execution/kit/kit-219/g0b.R
#
# (iii) per W1 fixture of tests/testthat/test-omic-summary-module.R (the lipid and metab feature setups for each inputs
#       case, the proteomics in-memory branches and the proteomics owner): whether KIT-212's guard lets Copy run, and which
#       files of `pathway_dir` the package rule would leave out. A case that Copy reaches and that leaves a file out is a
#       listed leaf (the ticket's OD-3). The four GUI-entry lanes are not rebuilt here: their stored captures are compared by
#       R2 (a moved leaf is a new failure on head).
# (iv)  per R2 file, every line that mentions the publication folder's enrichment plots or the enrichment workbook's sheets.
out_dir <- file.path("docs", "qa", "workflow-execution", "kit", "kit-219", "g0b-out")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
pkgload::load_all(".", export_all = TRUE, helpers = TRUE, quiet = TRUE)
Sys.setenv(KIT219_G0B_OUT = normalizePath(out_dir))
Sys.setenv(MULTISCHOLAR_TRX_MANIFEST_NOTE = "1")
result <- testthat::test_file(file.path("docs", "qa", "workflow-execution", "kit", "kit-219", "g0b-survey.R"),
    stop_on_failure = FALSE, reporter = "summary")
print(as.data.frame(result)[, c("test", "failed", "skipped", "error")])
cat("g0b done\n")
