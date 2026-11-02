# Extracted from test-module-ci-lipid-norm.R:49

# prequel ----------------------------------------------------------------------
library(testthat)

# test -------------------------------------------------------------------------
raw_object <- module_ci_lipid_norm_object(layout = "combined", positive_only = FALSE)
logged <- suppressWarnings(suppressMessages(logTransformAssays(raw_object, offset = 1)))
logged_mat <- module_ci_lipid_norm_matrix(logged, "LCMS_Pos")
expect_s4_class(logged, "LipidomicsAssayData")
module_ci_lipid_norm_assert_alignment(logged, c("LCMS_Pos", "GCMS"))
expect_equal(logged_mat["L_NEG", "S1"], 0)
expect_equal(logged_mat["L_ZERO", "S1"], 0)
expect_true(is.na(logged_mat["L_MISSING", "S1"]))
expect_true(isTRUE(logged@args$log_transformed))
expect_identical(logged@args$log_transform_offset, 1)
expect_error(logTransformAssays(raw_object, offset = 0), "single positive numeric")
for (method_name in c("none", "cyclicloess", "quantile", "scale")) {
    object <- module_ci_lipid_norm_object(
      layout = "combined",
      positive_only = TRUE,
      args = module_ci_lipid_norm_args(normalisation_method = method_name)
    )
    normalized <- suppressWarnings(suppressMessages(
      normaliseBetweenSamples(object, normalisation_method = method_name)
    ))
    mat <- module_ci_lipid_norm_matrix(normalized, "LCMS_Pos")

    expect_s4_class(normalized, "LipidomicsAssayData")
    module_ci_lipid_norm_assert_alignment(normalized, c("LCMS_Pos", "GCMS"))
    expect_identical(normalized@args$normalisation_method, method_name)
    expect_false(any(is.nan(mat)))
    expect_false(any(is.infinite(mat)))

    if (identical(method_name, "none")) {
      expect_equal(mat["L_STABLE", "S1"], 100)
      expect_equal(mat["L_ZERO", "S1"], 1)
    }
  }
pipeline <- suppressWarnings(suppressMessages(module_ci_lipid_norm_run_pipeline(
    input = module_ci_lipid_norm_input(norm_method = "quantile", ruv_mode = "skip"),
    object = module_ci_lipid_norm_object(positive_only = TRUE)
  )))
expect_true(pipeline$result)
expect_true(pipeline$norm_data$normalization_complete)
expect_true(pipeline$norm_data$ruv_complete)
expect_true("lipid_normalized" %in% pipeline$workflow$state_manager$getHistory())
module_ci_lipid_norm_assert_alignment(pipeline$norm_data$post_norm_obj, c("LCMS_Pos", "GCMS"))
