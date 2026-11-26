# KIT-268 G0 positive control: the tied_pvalues case of helper-module-ci-lipid-da.R through the real completeOmicDaRun.
# The DA kernel runLipidsDA runs with the helper's stub whose qvalue_warnings is non-empty (tied_pvalues); its result
# (assay -> comparison -> TRUE, as the real kernel returns it) goes to completeOmicDaRun, as in test-module-ci-lipid-da.R.
# No expectation on product behaviour: the observation is the trace's row (g0-observe.R).
test_that("G0 positive control: tied_pvalues lipid DA through completeOmicDaRun", {
  object <- module_ci_lipid_da_object(layout = "triple")
  contrasts <- module_ci_lipid_da_contrasts("multi_group")
  module_ci_lipid_da_bind_orchestration_stubs("tied_pvalues")
  suppressMessages(capture.output(
    results <- runLipidsDA(theObject = object, contrasts_tbl = contrasts, formula_string = "~ 0 + group",
      da_q_val_thresh = 0.05, treat_lfc_cutoff = 0)
  ))
  config <- omicGuiConfig("lipidomics")
  capture <- module_ci_lipid_da_capture_ui()
  workflow <- module_ci_lipid_da_workflow(object)
  da_env <- module_ci_lipid_da_data_env()
  testthat::local_mocked_bindings(updateSelectInput = capture$updateSelectInput, showNotification = capture$showNotification,
    removeNotification = capture$removeNotification, .package = "shiny")
  suppressMessages(completeOmicDaRun(results = results, tab = da_env, workflowData = workflow, session = "g0-control",
    experimentPaths = list(da_output_dir = NULL, publication_graphs_dir = NULL), da = config$da, binding = config$binding$binding))
  testthat::succeed("control ran")
})
