# Extracted from test-module-ci-metab-norm.R:279

# test -------------------------------------------------------------------------
object <- module_ci_metab_norm_object(layout = "combined", positive_only = TRUE)
norm_data <- module_ci_metab_norm_data(object)
norm_data$post_norm_obj <- object
workflow <- module_ci_metab_norm_workflow(object)
skip_state <- runMetabNormSkipCorrelationObserverEntry(
    workflowData = workflow,
    normData = norm_data,
    showNotificationFn = function(...) invisible(NULL)
  )
expect_identical(skip_state$status, "success")
expect_true(norm_data$correlation_filtering_complete)
expect_identical(norm_data$correlation_filtered_obj, object)
expect_true("metab_norm_complete" %in% workflow$state_manager$getHistory())
run_apply_case <- function(case) {
    # KIT-072: the declared pairs operation calls the namespace kernel, so the case's pairs are bound there.
    testthat::local_mocked_bindings(pearsonCorForSamplePairs = function(...) {
      module_ci_metab_norm_correlation_pairs(case, names(object@metabolite_data))
    }, .package = "MultiScholaR")
    case_workflow <- module_ci_metab_norm_workflow(object)
    case_norm_data <- module_ci_metab_norm_data(object)
    observer_state <- list(
      currentS4 = object,
      threshold = 0.75,
      groupingVariable = "tech_rep_group",
      notificationId = paste0("corr_", case)
    )
    suppressWarnings(suppressMessages(dispatchMetabNormApplyCorrelation(
      workflowData = case_workflow,
      normData = case_norm_data,
      observerState = observer_state,
      showNotificationFn = function(...) invisible(NULL),
      removeNotificationFn = function(...) invisible(NULL),
      reqFn = function(value) value
    ))) -> result
    list(result = result, workflow = case_workflow, norm_data = case_norm_data)
  }
pass_all <- run_apply_case("pass_all")
expect_identical(pass_all$result$status, "success")
expect_setequal(pass_all$norm_data$correlation_filtered_obj@design_matrix$Run, c("S1", "S2", "S3", "S4"))
