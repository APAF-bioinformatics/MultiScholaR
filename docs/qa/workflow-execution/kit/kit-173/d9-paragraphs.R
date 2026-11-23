# KIT-173 D9: the enrichment paragraph a g:Profiler run prints before and after this ticket, for the owner.
#
# Usage, from the package root on a run-tests runner, after the fix commit is installed:
#   Rscript docs/qa/workflow-execution/kit/kit-173/d9-paragraphs.R <out dir>
#
# The LFQ template's `FormatGOEnrichment` (the other three templates' paragraphs are the same text), at the base
# (tests/testdata/go-propagation/report-base/LFQ_report.R) and on this commit, for a g:Profiler run (recorded method
# "gprofiler2", q 0.05, cutoffs 0) with a recorded correction of gSCS, of fdr, and without a record (T11(d)'s inputs).
# Writes d9-<case>-base.txt and d9-<case>-head.txt.

args <- commandArgs(trailingOnly = TRUE)
out <- if (length(args) >= 1L) args[[1L]] else "k173-d9-out"
dir.create(out, showWarnings = FALSE, recursive = TRUE)
formatter <- function(path) {
    if (grepl("\\.rmd$", path)) {
        lines <- readLines(path, warn = FALSE, encoding = "UTF-8")
        start <- grep("^```\\{r GO Enrichment Section", lines)
        end <- start + which(lines[-seq_len(start)] == "```")[[1L]]
        chunk <- tempfile(fileext = ".Rmd")
        writeLines(lines[start:end], chunk, useBytes = TRUE)
        code <- tempfile(fileext = ".R")
        suppressMessages(knitr::purl(chunk, output = code, quiet = TRUE, documentation = 0))
        path <- code
    }
    exprs <- parse(path, keep.source = FALSE, encoding = "UTF-8")
    env <- new.env(parent = globalenv())
    for (e in exprs) if (is.call(e) && identical(e[[2L]], as.name("FormatGOEnrichment"))) eval(e, env)
    env$FormatGOEnrichment
}
base <- formatter(file.path("tests", "testdata", "go-propagation", "report-base", "LFQ_report.R"))
head <- formatter(file.path("inst", "reports", "proteomics", "LFQ_report.rmd"))
numbers <- list(ui_enrichment_organism_selected = "9606", ui_enrichment_up_lfc = 0, ui_enrichment_down_lfc = 0,
    ui_enrichment_q_value_cutoff = 0.05, ui_enrichment_database_source = "gprofiler2")
cases <- list(gSCS = list(enrichment_correction_method = "gSCS"), fdr = list(enrichment_correction_method = "fdr"),
    no_record = list())
for (name in names(cases)) {
    writeLines(as.character(base(c(numbers, cases[[name]]))), file.path(out, sprintf("d9-%s-base.txt", name)))
    writeLines(as.character(head(c(numbers, cases[[name]]))), file.path(out, sprintf("d9-%s-head.txt", name)))
}
print(list.files(out))
