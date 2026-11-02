# Extracted from test-module-ci-metab-norm.R:207

# test -------------------------------------------------------------------------
object <- module_ci_metab_norm_object(layout = "combined", positive_only = TRUE)
workflow <- module_ci_metab_norm_workflow(object)
norm_data <- module_ci_metab_norm_data(object)
skip_state <- runMetabNormRuvProgressApplyShell(
    currentS4 = object,
    totalSteps = 4,
    ruvMode = "skip",
    autoPercentageMin = 5,
    autoPercentageMax = 10,
    ruvGroupingVariable = "tech_rep_group",
    separationMetric = "max_difference",
    kPenaltyWeight = 0.5,
    adaptiveKPenalty = TRUE,
    manualK = 2L,
    manualPercentage = 40,
    experimentPaths = list(metabolite_qc_dir = tempdir()),
    groupingVariable = "group",
    shapeVariable = "batch",
    workflowData = workflow,
    normData = norm_data,
    incProgressFn = function(...) invisible(NULL)
  )
expect_identical(skip_state$skipLogEntry, "RUV-III skipped")
expect_true(norm_data$ruv_complete)
expect_identical(norm_data$ruv_corrected_obj, object)
optimizer_stub <- function(theObject, ruv_mode, params, experiment_paths) {
    if (identical(ruv_mode, "manual") && params$manual_k < 1) {
      stop("invalid k")
    }
    module_ci_metab_norm_ruv_results(theObject, mode = ruv_mode)
  }
identity_ruv <- function(theObject, ruv_grouping_variable, ruv_number_k, ctrl) {
    theObject@args$ruvIII_C_Varying <- list(ruv_grouping_variable = ruv_grouping_variable,
      ruv_number_k = ruv_number_k, ctrl = ctrl)
    theObject
  }
run_ruv <- function(mode, optimizer, ruv_fn = identity_ruv, manual_k = 2L, manual_percentage = 40,
                      percentage_min = 5, percentage_max = 10, adaptive = FALSE) {
    testthat::local_mocked_bindings(runPerAssayRuvOptimization = optimizer, ruvIII_C_Varying = ruv_fn,
      .package = "MultiScholaR")
    case_workflow <- module_ci_metab_norm_workflow(object)
    case_norm_data <- module_ci_metab_norm_data(object)
    state <- runMetabNormRuvProgressApplyShell(
      currentS4 = object, totalSteps = 4, ruvMode = mode,
      autoPercentageMin = percentage_min, autoPercentageMax = percentage_max,
      ruvGroupingVariable = "tech_rep_group", separationMetric = "max_difference", kPenaltyWeight = 0.5,
      adaptiveKPenalty = adaptive, manualK = manual_k, manualPercentage = manual_percentage,
      experimentPaths = list(metabolite_qc_dir = tempdir()), groupingVariable = "group", shapeVariable = "batch",
      workflowData = case_workflow, normData = case_norm_data, incProgressFn = function(...) invisible(NULL),
      runRuvQcStepFn = function(currentS4, ...) list(currentS4 = currentS4)
    )
    list(state = state, workflow = case_workflow, norm_data = case_norm_data)
  }
manual <- run_ruv("manual", optimizer_stub)
manual_k <- manual$state$optimizationState$bestKPerAssay
expect_named(manual_k, c("LCMS_Pos", "GCMS"))
expect_identical(unname(unlist(manual_k)), c(2L, 2L))
expect_true(all(vapply(extractCtrlPerAssay(manual$state$optimizationState$ruvResults), function(x) {
    is.logical(x) && !is.null(names(x)) && sum(x) == 6L
  }, logical(1))))
auto <- run_ruv("automatic", optimizer_stub, percentage_min = 10, percentage_max = 12, adaptive = TRUE)
expect_identical(unname(unlist(auto$state$optimizationState$bestKPerAssay)), c(1L, 2L))
expect_equal(auto$state$optimizationState$ruvResults$LCMS_Pos$best_percentage, 13)
expect_s3_class(auto$state$optimizationState$ruvResults$GCMS$optimization_results, "data.frame")
expect_s4_class(auto$state$currentS4, "MetaboliteAssayData")
expect_true(auto$norm_data$ruv_complete)
expect_true("metab_ruv_corrected" %in% auto$workflow$state_manager$getHistory())
expect_identical(auto$state$currentS4@args$ruvIII_C_Varying$ruv_number_k$LCMS_Pos, 1L)
