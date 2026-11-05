# KIT-094 test-first, per site: on the base commit the stub binding cannot even register (the core builds the 13 lipid
# and metab steps from any binding), so every stub test stops at the registry. This diagnostic run lets only
# `omicsWorkflowOperations` honour a binding's `steps`, leaving every other site as the base commit has it, and runs the
# test file again: each remaining failure is one site that does not yet read the binding. Diagnostic only; run on the
# base commit.
#
# Usage: Rscript test-first-sites.R <package root>
args <- commandArgs(TRUE)
suppressMessages(pkgload::load_all(args[[1L]], quiet = TRUE, export_all = TRUE, helpers = TRUE))
original <- get("omicsWorkflowOperations", asNamespace("MultiScholaR"))
testthat::local_mocked_bindings(omicsWorkflowOperations = function(binding) {
    if (!is.function(binding$steps)) return(original(binding))
    declarations <- binding$steps(binding)
    stats::setNames(declarations, vapply(declarations, `[[`, "", "name"))
}, .package = "MultiScholaR", .env = globalenv())
results <- testthat::test_file(file.path(args[[1L]], "tests", "testthat", "test-workflow-omics-binding-steps.R"),
    reporter = "silent", stop_on_failure = FALSE)
frame <- as.data.frame(results)
cat(sprintf("COUNTS test-workflow-omics-binding-steps.R (steps honoured only) pass %d fail %d error %d skip %d\n",
    sum(frame$passed), sum(frame$failed), sum(frame$error), sum(frame$skipped)))
for (result in results) for (expectation in result$results) {
    if (inherits(expectation, c("expectation_success", "expectation_skip"))) next
    cat("---", result$test, "\n", substr(gsub("\\s+", " ", conditionMessage(expectation)), 1L, 400L), "\n")
}
