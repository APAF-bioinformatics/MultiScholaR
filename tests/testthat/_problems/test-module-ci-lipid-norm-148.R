# Extracted from test-module-ci-lipid-norm.R:148

# prequel ----------------------------------------------------------------------
library(testthat)

# test -------------------------------------------------------------------------
object <- module_ci_lipid_norm_object(layout = "combined", positive_only = TRUE)
skip_pipeline <- suppressWarnings(suppressMessages(module_ci_lipid_norm_run_pipeline(
    input = module_ci_lipid_norm_input(norm_method = "none", ruv_mode = "skip"),
    object = object
  )))
expect_true(skip_pipeline$result)
expect_true(skip_pipeline$norm_data$ruv_complete)
expect_identical(skip_pipeline$norm_data$ruv_corrected_obj, skip_pipeline$norm_data$post_norm_obj)
expect_false("lipid_ruv_corrected" %in% skip_pipeline$workflow$state_manager$getHistory())
manual_pipeline <- suppressWarnings(suppressMessages(module_ci_lipid_norm_run_pipeline(
    input = module_ci_lipid_norm_input(norm_method = "none", ruv_mode = "manual"),
    object = object
  )))
expect_true(manual_pipeline$result)
expect_true(manual_pipeline$norm_data$ruv_complete)
expect_true("lipid_ruv_corrected" %in% manual_pipeline$workflow$state_manager$getHistory())
expect_identical(manual_pipeline$norm_data$ruv_corrected_obj@args$ruv_number_k$LCMS_Pos, 2L)
