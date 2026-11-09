# KIT-112 W2 probe (scratch): F112-03 on the managed path. Two contrasts at TREAT 0.5, q 0.01 through the worker and
# the publication step (pages saved as assets by the controller), then the first contrast alone (product a hit, the
# manifest new). Reports whether the publication completes and whether the page's asset lookup stays unambiguous.
#   Rscript probe-f11203-worker.R <tree>
args <- commandArgs(TRUE)
tree <- args[[1]]
suppressMessages(pkgload::load_all(tree, quiet = TRUE, export_all = TRUE, helpers = TRUE))
library(testthat)
for (expression in parse(file.path(tree, "tests", "testthat", "test-workflow-saved-da.R"))) {
    if (is.call(expression) && identical(expression[[1L]], as.name("<-")) && grepl("^savedDaW[12]", deparse(expression[[2L]])))
        eval(expression, globalenv())
}
testthat::with_mocked_bindings({
    fixture <- savedDaW1Owner("k112-f11203-worker", envir = globalenv())
    dispatch <- savedDaW2Dispatcher(fixture, NULL)
    manager <- fixture$workflow$state_manager
    settings <- list(formula_string = "~ 0 + group", da_q_val_thresh = 0.01, treat_lfc_cutoff = 0.5)
    step <- function(n) tryCatch({
        utils::capture.output(m <- suppressWarnings(suppressMessages(savedDaW2Managed(fixture, dispatch, settings, savedDaW1Contrasts(n)))))
        "ok"
    }, error = function(e) conditionMessage(e))
    cat("two contrasts:", step(2L), "\n")
    cat("one contrast (hit, new manifest):", step(1L), "\n")
    snapshot <- withDIAStateScope(manager, function() manager$.context$snapshot)
    manifest <- resolveWorkflowStateValue(manager$.context, snapshot$da_results)
    cat("head outputs:", paste(names(unlist(snapshot$completed_outputs[[getDIAStateReference(manager)$entity_id]])), collapse = ","), "\n")
    lookup <- tryCatch(withDIAStateScope(manager, function() {
        identity <- diaInteractiveVolcanoIdentity(manager, manifest, "treated.vs.control", 0.01)
        findSavedWorkflowEntities(manager$.context$store, list(kind = "rendered_asset", identity = identity))$status
    }), error = function(e) paste("error:", conditionMessage(e)))
    cat("asset lookup for the reused page:", lookup, "\n")
}, workflowExecutionBuildIdentity = function() strrep("f", 64),
    finishWorkflowExecutionJob = function(owner, status) { owner$job <- NULL; list(status = status) },
    .package = "MultiScholaR")
