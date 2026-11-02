# Extracted from test-module-ci-metab-norm.R:52

# test -------------------------------------------------------------------------
raw_object <- module_ci_metab_norm_object(layout = "combined", positive_only = FALSE)
logged <- suppressWarnings(suppressMessages(logTransformAssays(raw_object, offset = 1)))
logged_mat <- module_ci_metab_norm_matrix(logged, "LCMS_Pos")
expect_s4_class(logged, "MetaboliteAssayData")
module_ci_metab_norm_assert_alignment(logged, c("LCMS_Pos", "GCMS"))
expect_equal(logged_mat["M_NEG", "S1"], 0)
expect_equal(logged_mat["M_ZERO", "S1"], 0)
expect_true(is.na(logged_mat["M_MISSING", "S1"]))
expect_true(isTRUE(logged@args$log_transformed))
expect_identical(logged@args$log_transform_offset, 1)
expect_error(logTransformAssays(raw_object, offset = 0), "single positive numeric")
for (method_name in c("none", "cyclicloess", "quantile", "scale")) {
    object <- module_ci_metab_norm_object(
      layout = "combined",
      positive_only = TRUE,
      args = module_ci_metab_norm_args(normalisation_method = method_name)
    )
    normalized <- suppressWarnings(suppressMessages(
      normaliseBetweenSamples(object, normalisation_method = method_name)
    ))
    mat <- module_ci_metab_norm_matrix(normalized, "LCMS_Pos")

    expect_s4_class(normalized, "MetaboliteAssayData")
    module_ci_metab_norm_assert_alignment(normalized, c("LCMS_Pos", "GCMS"))
    expect_identical(normalized@args$normalisation_method, method_name)
    expect_false(any(is.nan(mat)))
    expect_false(any(is.infinite(mat)))

    if (identical(method_name, "none")) {
      expect_equal(mat["M_STABLE", "S1"], 100)
      expect_equal(mat["M_ZERO", "S1"], 1)
    }
  }
workflow <- module_ci_metab_norm_workflow(module_ci_metab_norm_object(positive_only = TRUE))
norm_data <- module_ci_metab_norm_data(workflow$state_manager$getState())
shell_state <- suppressWarnings(suppressMessages(runMetabNormBetweenSampleProgressApplyShell(
    currentS4 = workflow$state_manager$getState(),
    totalSteps = 4,
    normMethod = "quantile",
    workflowData = workflow,
    normData = norm_data,
    incProgressFn = function(...) invisible(NULL)
  )))
expect_s4_class(shell_state$currentS4, "MetaboliteAssayData")
expect_true(norm_data$normalization_complete)
expect_true("metab_normalized" %in% workflow$state_manager$getHistory())
module_ci_metab_norm_assert_alignment(shell_state$currentS4, c("LCMS_Pos", "GCMS"))
