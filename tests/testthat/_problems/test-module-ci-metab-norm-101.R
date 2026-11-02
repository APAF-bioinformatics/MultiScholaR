# Extracted from test-module-ci-metab-norm.R:101

# test -------------------------------------------------------------------------
object <- module_ci_metab_norm_object(layout = "combined", positive_only = TRUE)
regex_norm <- suppressWarnings(suppressMessages(normaliseUntransformedData(
    object,
    method = "ITSD",
    itsd_aggregation = "median"
  )))
expect_false(any(grepl("^IS_", module_ci_metab_norm_feature_ids(regex_norm, "LCMS_Pos"))))
expect_false(any(grepl("^IS_", module_ci_metab_norm_feature_ids(regex_norm, "GCMS"))))
expect_true(isTRUE(regex_norm@args$ITSDNormalization$applied))
expect_identical(regex_norm@args$ITSDNormalization$itsd_aggregation, "median")
expect_identical(regex_norm@args$ITSDNormalization$itsd_counts_per_assay$LCMS_Pos, 2L)
expect_identical(regex_norm@args$ITSDNormalization$itsd_counts_per_assay$GCMS, 2L)
manual_indices <- lapply(object@metabolite_data, function(assay) {
    match(c("IS_A", "IS_B"), assay$database_identifier)
  })
manual_ids <- resolveMetabNormManualItsdFeatureIds(
    currentS4 = object,
    itsdSelections = manual_indices,
    buildSelectionTableFn = function(assay_data, metabolite_id_col, annotation_cols) {
      data.frame(feature_id = assay_data[[metabolite_id_col]], stringsAsFactors = FALSE)
    }
  )
expect_identical(manual_ids$LCMS_Pos, c("IS_A", "IS_B"))
expect_identical(manual_ids$GCMS, c("IS_A", "IS_B"))
workflow <- module_ci_metab_norm_workflow(object)
norm_data <- module_ci_metab_norm_data(object)
itsd_state <- suppressWarnings(suppressMessages(runMetabNormItsdProgressApplyShell(
    currentS4 = object,
    totalSteps = 4,
    applyItsd = TRUE,
    itsdAggregation = "mean",
    itsdSelections = manual_indices,
    workflowData = workflow,
    normData = norm_data,
    incProgressFn = function(...) invisible(NULL),
    resolveManualFeatureIdsFn = function(...) manual_ids
  )))
expect_true(itsd_state$applied)
expect_s4_class(norm_data$post_itsd_obj, "MetaboliteAssayData")
expect_true("metab_itsd_norm" %in% workflow$state_manager$getHistory())
expect_identical(norm_data$post_itsd_obj@args$ITSDNormalization$itsd_features_per_assay$LCMS_Pos, c("IS_A", "IS_B"))
