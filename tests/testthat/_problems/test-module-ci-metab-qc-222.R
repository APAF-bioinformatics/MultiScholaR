# Extracted from test-module-ci-metab-qc.R:222

# test -------------------------------------------------------------------------
current <- module_ci_metab_qc_object(layout = "combined")
workflow <- module_ci_metab_qc_workflow(current)
output <- new.env(parent = emptyenv())
captured <- new.env(parent = emptyenv())
runMetabQcS4FinalizeWorkflow(
    workflowData = workflow,
    omicType = "metabolomics",
    filterPlot = function(value = NULL) {
      if (!missing(value)) captured$filter_plot <- value
      captured$filter_plot
    },
    output = output,
    updateTrackingPlotFn = function(currentS4, omicType, setFilterPlotFn, ...) {
      captured$tracking <- list(currentS4 = currentS4, omicType = omicType)
      setFilterPlotFn("qc-progress-plot")
      "qc-progress-plot"
    },
    reportFinalizeSuccessFn = function(currentS4, history, output, ...) {
      captured$success <- list(currentS4 = currentS4, history = history)
      output$finalize_results <- "finalized"
      invisible("finalized")
    },
    reportFinalizeErrorFn = function(error, ...) {
      captured$error <- conditionMessage(error)
      invisible(captured$error)
    }
  )
expect_identical(workflow$tab_status$quality_control, "complete")
expect_identical(captured$filter_plot, "qc-progress-plot")
expect_identical(captured$tracking$omicType, "metabolomics")
expect_true("metab_qc_complete" %in% workflow$state_manager$getHistory())
expect_identical(output$finalize_results, "finalized")
expect_true(any(captured$success$history == "metab_qc_complete"))
invalid_workflow <- module_ci_metab_qc_workflow(current)
invalid_workflow$state_manager$saveState("metab_not_s4", list(not = "s4"), invalid_workflow$config_list, "not an S4")
runMetabQcS4FinalizeWorkflow(
    workflowData = invalid_workflow,
    omicType = "metabolomics",
    filterPlot = function(...) NULL,
    output = new.env(parent = emptyenv()),
    reportFinalizeSuccessFn = function(...) stop("should not finalize"),
    reportFinalizeErrorFn = function(error, ...) {
      captured$invalid_error <- conditionMessage(error)
      invisible(NULL)
    }
  )
expect_match(captured$invalid_error, "Current state is not a MetaboliteAssayData object", fixed = TRUE)
