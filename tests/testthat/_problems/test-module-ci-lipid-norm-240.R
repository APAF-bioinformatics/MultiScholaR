# Extracted from test-module-ci-lipid-norm.R:240

# prequel ----------------------------------------------------------------------
library(testthat)

# test -------------------------------------------------------------------------
object <- module_ci_lipid_norm_object(layout = "combined", positive_only = TRUE)
norm_data <- module_ci_lipid_norm_data(object)
norm_data$post_norm_obj <- object
norm_data$normalization_complete <- TRUE
workflow <- module_ci_lipid_norm_workflow(object)
skip_state <- handleLipidNormSkipCorrelationFilter(
    workflowData = workflow,
    normData = norm_data,
    addLog = function(...) invisible(NULL),
    reqFn = function(value) value,
    showNotificationFn = function(...) invisible(NULL)
  )
expect_true(skip_state)
expect_true(norm_data$correlation_filtering_complete)
expect_identical(norm_data$correlation_filtered_obj, object)
expect_true("lipid_norm_complete" %in% workflow$state_manager$getHistory())
expect_identical(workflow$tab_status$normalization, "complete")
run_apply_case <- function(case) {
    # KIT-072: the declared pairs operation calls the namespace kernel, so the case's pairs are bound there.
    testthat::local_mocked_bindings(pearsonCorForSamplePairs = function(...) {
      module_ci_lipid_norm_correlation_pairs(case, names(object@lipid_data))
    }, .package = "MultiScholaR")
    case_workflow <- module_ci_lipid_norm_workflow(object)
    case_norm_data <- module_ci_lipid_norm_data(object)
    case_norm_data$normalization_complete <- TRUE
    case_norm_data$post_norm_obj <- object
    suppressWarnings(suppressMessages(handleLipidNormApplyCorrelationFilter(
      input = module_ci_lipid_norm_input(correlation_threshold = 0.75),
      workflowData = case_workflow,
      normData = case_norm_data,
      addLog = function(...) invisible(NULL),
      reqFn = function(value) value,
      showNotificationFn = function(...) invisible(NULL),
      removeNotificationFn = function(...) invisible(NULL),
      logInfoFn = function(...) invisible(NULL),
      logErrorFn = function(...) invisible(NULL)
    ))) -> result
    list(result = result, workflow = case_workflow, norm_data = case_norm_data)
  }
pass_all <- run_apply_case("pass_all")
expect_true(pass_all$result)
expect_setequal(pass_all$norm_data$correlation_filtered_obj@design_matrix$Run, c("S1", "S2", "S3", "S4"))
