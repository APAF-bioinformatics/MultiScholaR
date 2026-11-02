# Extracted from test-module-ci-lipid-norm.R:100

# prequel ----------------------------------------------------------------------
library(testthat)

# test -------------------------------------------------------------------------
object <- module_ci_lipid_norm_object(layout = "combined", positive_only = TRUE)
regex_norm <- suppressWarnings(suppressMessages(normaliseUntransformedData(
    object,
    method = "ITSD",
    itsd_aggregation = "median"
  )))
expect_false(any(grepl("^IS_", module_ci_lipid_norm_feature_ids(regex_norm, "LCMS_Pos"))))
expect_false(any(grepl("^IS_", module_ci_lipid_norm_feature_ids(regex_norm, "GCMS"))))
expect_true(isTRUE(regex_norm@args$ITSDNormalization$applied))
expect_identical(regex_norm@args$ITSDNormalization$itsd_aggregation, "median")
expect_identical(regex_norm@args$ITSDNormalization$itsd_counts_per_assay$LCMS_Pos, 2L)
expect_identical(regex_norm@args$ITSDNormalization$itsd_counts_per_assay$GCMS, 2L)
manual_ids <- lapply(object@lipid_data, function(assay) {
    c("IS_A", "IS_B")
  })
manual_norm <- suppressWarnings(suppressMessages(normaliseUntransformedData(
    object,
    method = "ITSD",
    itsd_feature_ids = manual_ids,
    itsd_aggregation = "mean"
  )))
expect_identical(manual_norm@args$ITSDNormalization$itsd_features_per_assay$LCMS_Pos, c("IS_A", "IS_B"))
expect_identical(manual_norm@args$ITSDNormalization$itsd_counts_per_assay$GCMS, 2L)
manual_indices <- lapply(object@lipid_data, function(assay) {
    match(c("IS_A", "IS_B"), assay$lipid_id)
  })
workflow <- module_ci_lipid_norm_workflow(object)
norm_data <- module_ci_lipid_norm_data(object)
norm_data$itsd_selections <- manual_indices
itsd_pipeline <- suppressWarnings(suppressMessages(module_ci_lipid_norm_run_pipeline(
    input = module_ci_lipid_norm_input(
      norm_method = "none",
      ruv_mode = "skip",
      apply_itsd = TRUE,
      itsd_aggregation = "mean"
    ),
    object = object,
    workflow = workflow,
    norm_data = norm_data
  )))
expect_true(itsd_pipeline$result)
expect_s4_class(norm_data$post_itsd_obj, "LipidomicsAssayData")
expect_true("lipid_itsd_norm" %in% workflow$state_manager$getHistory())
expect_identical(norm_data$post_itsd_obj@args$ITSDNormalization$itsd_features_per_assay$LCMS_Pos, c("IS_A", "IS_B"))
