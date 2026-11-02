# Extracted from test-module-ci-lipid-qc.R:104

# prequel ----------------------------------------------------------------------
library(testthat)

# test -------------------------------------------------------------------------
assay <- module_ci_lipid_qc_assay(include_duplicates = TRUE, include_itsd = FALSE)
detected <- findLipidDuplicateFeatureIDs(module_ci_lipid_qc_object(layout = "gc", include_duplicates = TRUE, include_itsd = FALSE))
expect_true("L_DUP" %in% detected$GCMS$lipid_id)
expect_true("L_TIE" %in% detected$GCMS$lipid_id)
resolved_ids <- resolveLipidDuplicateFeaturesByIntensity(
    assay_tibble = assay,
    id_col = "lipid_id",
    sample_cols = module_ci_lipid_qc_samples()
  )
expect_false(anyDuplicated(resolved_ids$lipid_id) > 0L)
expect_identical(
    resolved_ids$S1[resolved_ids$lipid_id == "L_DUP"],
    90
  )
expect_identical(
    resolved_ids$S1[resolved_ids$lipid_id == "L_TIE"],
    50
  )
resolved_annotations <- resolveLipidDuplicateFeaturesByIntensity(
    assay_tibble = assay,
    id_col = "lipid",
    sample_cols = module_ci_lipid_qc_samples()
  )
same_annotation <- resolved_annotations[resolved_annotations$lipid == "SameAnnotation", , drop = FALSE]
expect_equal(nrow(same_annotation), 1L)
expect_identical(same_annotation$lipid_id, "L_ANN_B")
workflow <- module_ci_lipid_qc_workflow(module_ci_lipid_qc_object(layout = "combined", include_duplicates = TRUE, include_itsd = FALSE))
dispatched <- handleLipidDuplicateResolution(
    workflowData = workflow,
    omicType = "lipidomics",
    reqFn = function(value) value,
    updateFilteringFn = function(...) "qc-plot",
    logWarnFn = function(...) invisible(NULL)
  )
expect_identical(names(dispatched$statsList), c("LCMS_Pos", "GCMS"))
expect_true(all(vapply(dispatched$statsList, `[[`, numeric(1), "removed") >= 2))
expect_true("lipid_duplicates_resolved" %in% workflow$state_manager$getHistory())
expect_identical(names(dispatched$currentS4@lipid_data), c("LCMS_Pos", "GCMS"))
