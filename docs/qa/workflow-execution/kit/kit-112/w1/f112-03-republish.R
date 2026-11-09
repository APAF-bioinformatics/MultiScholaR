# KIT-112 W1, draft finding F112-03 (reproduction, start commit's code, no mocks): on one DIA owner, the DA with two
# contrasts (TREAT 0.5, q 0.01; Glimma publishes a page per contrast and saves it as a rendered asset), then the DA
# with the first contrast only on the same head (its product a hit, the manifest new). The in-process publication
# fails in `saveSelfContainedWidget`: `generateProtDAVolcanoPlotGlimma` finds the saved page (KIT-016) and returns a
# `shiny` iframe tag, which is not an htmlwidget, so `htmlwidgets::saveWidget` is called on it and stops.
#   Rscript docs/qa/workflow-execution/kit/kit-112/w1/f112-03-republish.R <package source>
args <- commandArgs(TRUE)
tree <- args[[1]]
setwd(tree)
suppressMessages(pkgload::load_all(tree, quiet = TRUE, export_all = TRUE, helpers = TRUE))
library(testthat)
for (expression in parse(file.path(tree, "tests", "testthat", "test-workflow-saved-da.R"))) {
    if (is.call(expression) && identical(expression[[1L]], as.name("<-")) && grepl("^savedDaW1", deparse(expression[[2L]])))
        eval(expression, globalenv())
}
fixture <- savedDaW1Owner("k112-debug", envir = globalenv())
run <- function(contrasts, q = 0.05, lfc = 0) {
    fixture$workflow$contrasts_tbl <- contrasts
    settings <- list(formula_string = "~ 0 + group", da_q_val_thresh = q, treat_lfc_cutoff = lfc)
    calls <- NULL
    body <- function() utils::capture.output(suppressWarnings(suppressMessages(
        runDIAReferenceDaAnalysis(settings, fixture$da, fixture$workflow, fixture$paths, contrasts, NULL))))
    out <- tryCatch(withCallingHandlers(body(), error = function(e) calls <<- sys.calls()),
        error = function(e) conditionMessage(e))
    cat("RESULT:", if (is.character(out) && length(out) == 1L && !is.null(calls)) out else "ok", "\n")
    if (!is.null(calls)) {
        cs <- vapply(calls, function(x) paste(deparse(x, nlines = 1L), collapse = ""), "")
        cat(paste0("  ", substr(utils::tail(cs, 25L), 1L, 160L)), sep = "\n")
    }
}
cat("== two contrasts TREAT\n"); run(savedDaW1Contrasts(2L), 0.01, 0.5)
cat("== one contrast TREAT hit\n"); run(savedDaW1Contrasts(1L), 0.01, 0.5)
