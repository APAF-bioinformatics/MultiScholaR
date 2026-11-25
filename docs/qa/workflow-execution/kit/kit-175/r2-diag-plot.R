# KIT-175 R2 diagnosis (no result judged): the text behind the `gprofiler_plot` digest that R2 found changed in
# test-omic-enrich-module.R's refusals/pre_da session (steps run_before_da, da, open). Runs that session's drive as the
# test does, with the helpers and the test file's definitions sourced (test_that a no-op), and writes each step's
# masked plot text in full instead of its digest.
# Usage: Rscript docs/qa/workflow-execution/kit/kit-175/r2-diag-plot.R <out dir>
args <- commandArgs(trailingOnly = TRUE)
out <- if (length(args)) args[[1L]] else "k175-diag-out"
dir.create(out, showWarnings = FALSE, recursive = TRUE)
suppressPackageStartupMessages({library(MultiScholaR); library(testthat)})
h <- new.env(parent = asNamespace("MultiScholaR"))
for (f in sort(list.files("tests/testthat", "^(setup|helper).*\\.R$", full.names = TRUE))) {
    tryCatch(sys.source(f, envir = h), error = function(e) message("skipped ", f, ": ", conditionMessage(e)))
}
h$test_that <- function(...) invisible(NULL)
sys.source("tests/testthat/test-omic-enrich-module.R", envir = h)
h$guiEntryOutputValue <- function(value, roots) {
    if (is.character(value)) return(h$guiEntryMask(paste(value, collapse = "\n"), roots))
    paste(utils::capture.output(utils::str(value)), collapse = "\n")
}
seen <- h$omicEnrichW1Case("refusals")
for (step in names(seen)) {
    values <- seen[[step]]$values
    if (is.null(values)) next
    for (name in intersect(c("gprofiler_plot", "clusterprofiler_plot"), names(values))) {
        text <- values[[name]]
        writeLines(c(sprintf("== %s %s sha256=%s", step, name,
            digest::digest(text, algo = "sha256", serialize = FALSE)), text), file.path(out,
            sprintf("%s__%s.txt", gsub("[^A-Za-z0-9_]", "_", step), name)))
    }
}
writeLines(names(seen), file.path(out, "steps.txt"))
