library(testthat)

test_that("MCI-020.1 intensity filter matrix covers thresholds, zeros, missingness, non-finite values, and per-assay effects", {
  assay <- module_ci_lipid_qc_assay(include_duplicates = FALSE, include_itsd = FALSE)

  strict <- lipidIntensityFilteringHelper(
    assay_table = assay,
    min_lipid_intensity_threshold = 50,
    lipids_proportion_of_samples_below_cutoff = 0.5,
    lipid_id_column = "lipid_id"
  )
  expect_identical(module_ci_lipid_qc_ids(strict), c("L_keep", "L_missing"))
  expect_false("L_boundary" %in% module_ci_lipid_qc_ids(strict))
  expect_false("L_zero" %in% module_ci_lipid_qc_ids(strict))

  relaxed <- lipidIntensityFilteringHelper(
    assay_table = assay,
    min_lipid_intensity_threshold = 50,
    lipids_proportion_of_samples_below_cutoff = 0.51,
    lipid_id_column = "lipid_id"
  )
  expect_true("L_boundary" %in% module_ci_lipid_qc_ids(relaxed))

  all_pass <- lipidIntensityFilteringHelper(
    assay_table = assay,
    min_lipid_intensity_threshold = -1,
    lipids_proportion_of_samples_below_cutoff = 1,
    lipid_id_column = "lipid_id"
  )
  expect_equal(nrow(all_pass), nrow(assay))

  all_fail <- lipidIntensityFilteringHelper(
    assay_table = assay,
    min_lipid_intensity_threshold = 10000,
    lipids_proportion_of_samples_below_cutoff = 0.01,
    lipid_id_column = "lipid_id"
  )
  expect_equal(nrow(all_fail), 0L)

  expect_error(
    lipidIntensityFilteringHelper(
      assay_table = module_ci_lipid_qc_assay(include_duplicates = FALSE, include_itsd = FALSE, nonfinite = TRUE),
      min_lipid_intensity_threshold = 50,
      lipids_proportion_of_samples_below_cutoff = 0.5,
      lipid_id_column = "lipid_id"
    ),
    "non-finite sample intensity",
    fixed = TRUE
  )

  s4 <- module_ci_lipid_qc_object(layout = "combined", include_duplicates = FALSE, include_itsd = FALSE)
  filtered <- lipidIntensityFiltering(
    theObject = s4,
    lipids_intensity_cutoff_percentile = 20,
    lipids_proportion_of_samples_below_cutoff = 0.5
  )
  expect_s4_class(filtered, "LipidomicsAssayData")
  expect_identical(names(filtered@lipid_data), c("LCMS_Pos", "GCMS"))
  expect_true(all(module_ci_lipid_qc_feature_counts(filtered) <= module_ci_lipid_qc_feature_counts(s4)))
})

test_that("MCI-020.2 duplicate matrix covers IDs, annotations/classes, ties, conflicts, and resolution stats", {
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
  expect_true("lipid_duplicates_resolved" %in% names(workflow$state_manager$saved()))
  expect_identical(names(dispatched$currentS4@lipid_data), c("LCMS_Pos", "GCMS"))
})

test_that("MCI-020.3 ITSD matrix covers global regex, assay-specific matches, no matches, multiple matches, and invalid selections", {
  no_itsd <- module_ci_lipid_qc_object(layout = "gc", include_itsd = FALSE)
  empty_metrics <- getLipidInternalStandardMetrics(
    assay_data = no_itsd@lipid_data$GCMS,
    is_pattern = "^IS_",
    lipid_id_col = "lipid_id",
    sample_id_col = "Run",
    sample_columns = module_ci_lipid_qc_samples()
  )
  expect_equal(nrow(empty_metrics), 0L)

  s4 <- module_ci_lipid_qc_object(layout = "combined", include_itsd = TRUE)
  metrics_by_assay <- lapply(s4@lipid_data, function(assay) {
    getLipidInternalStandardMetrics(
      assay_data = assay,
      is_pattern = "^IS_",
      lipid_id_col = "lipid_id",
      sample_id_col = "Run",
      sample_columns = module_ci_lipid_qc_samples()
    )
  })
  expect_identical(names(metrics_by_assay), c("LCMS_Pos", "GCMS"))
  expect_true(all(vapply(metrics_by_assay, nrow, integer(1)) >= 2L))
  expect_true(all(vapply(metrics_by_assay, function(metrics) all(metrics$cv < 5), logical(1))))

  assay_specific <- getLipidInternalStandardMetrics(
    assay_data = s4@lipid_data$LCMS_Pos,
    is_pattern = "^IS_By_Name$",
    lipid_id_col = "lipid",
    sample_id_col = "Run",
    sample_columns = module_ci_lipid_qc_samples()
  )
  expect_identical(assay_specific$is_id, "IS_By_Name")

  no_match_with_fallback <- getLipidInternalStandardMetrics(
    assay_data = s4@lipid_data$LCMS_Pos,
    is_pattern = "^NO_MATCH$",
    lipid_id_col = "lipid_id",
    sample_id_col = "Run",
    sample_columns = module_ci_lipid_qc_samples()
  )
  expect_true(nrow(no_match_with_fallback) >= 2L)

  invalid_pattern <- suppressWarnings(getLipidInternalStandardMetrics(
    assay_data = s4@lipid_data$LCMS_Pos,
    is_pattern = "[",
    lipid_id_col = "lipid_id",
    sample_id_col = "Run",
    sample_columns = module_ci_lipid_qc_samples()
  ))
  expect_true(nrow(invalid_pattern) >= 2L)

  missing_column <- suppressWarnings(getLipidInternalStandardMetrics(
    assay_data = s4@lipid_data$LCMS_Pos,
    is_pattern = "^IS_",
    lipid_id_col = "missing_lipid_id",
    sample_id_col = "Run",
    sample_columns = module_ci_lipid_qc_samples()
  ))
  expect_equal(nrow(missing_column), 0L)
})

test_that("MCI-020.4 S4 plotting matrix covers PCA, Pearson, RLE, density, small-n, constants, and missing group labels", {
  s4 <- module_ci_lipid_qc_object(layout = "gc", include_duplicates = FALSE, include_itsd = FALSE)

  pca <- suppressWarnings(suppressMessages(plotPca(s4, grouping_variable = "group", title = "Module CI PCA")))
  expect_true(is.list(pca))
  expect_named(pca, "GCMS")

  rle <- suppressWarnings(suppressMessages(plotRle(s4, grouping_variable = "group")))
  expect_true(is.list(rle))
  expect_named(rle, "GCMS")

  density <- suppressWarnings(suppressMessages(plotDensity(s4, grouping_variable = "group", title = "Module CI Density")))
  expect_true(is.list(density))

  pearson <- suppressWarnings(suppressMessages(plotPearson(s4, tech_rep_remove_regex = "", correlation_group = "tech_rep_group")))
  expect_true(is.list(pearson))

  small_n <- calculateLipidCVs(
    assay_data = module_ci_lipid_qc_assay(include_duplicates = FALSE, include_itsd = FALSE),
    design_matrix = module_ci_lipid_qc_design(small_n = TRUE),
    group_id_col = "group",
    replicate_id_col = "tech_rep_group",
    sample_id_col = "Run",
    lipid_id_col = "lipid_id",
    sample_columns = module_ci_lipid_qc_samples()
  )
  expect_true(all(is.na(small_n$cv)))

  constant <- module_ci_lipid_qc_object(layout = "gc", include_duplicates = FALSE, include_itsd = FALSE, constant = TRUE)
  constant_pearson <- suppressWarnings(plotPearson(constant, tech_rep_remove_regex = "", correlation_group = "tech_rep_group"))
  expect_true(is.list(constant_pearson))

  missing_group <- module_ci_lipid_qc_object(layout = "gc", include_duplicates = FALSE, include_itsd = FALSE, missing_group = TRUE)
  expect_error(
    validateLipidQcS4FinalizeState(missing_group, reqFn = function(value) value),
    "blank assignments",
    fixed = TRUE
  )
})

test_that("MCI-020.5 finalization validates lipid S4 state and unlocks normalization handoff", {
  current <- module_ci_lipid_qc_object(layout = "combined")
  workflow <- module_ci_lipid_qc_workflow(current)
  output <- new.env(parent = emptyenv())
  captured <- new.env(parent = emptyenv())

  runLipidQcS4FinalizeWorkflow(
    workflowData = workflow,
    omicType = "lipidomics",
    filterPlot = function(value = NULL) {
      if (!missing(value)) captured$filter_plot <- value
      captured$filter_plot
    },
    output = output,
    updateTrackingPlotFn = function(currentS4, omicType, setFilterPlotFn, ...) {
      captured$tracking <- list(currentS4 = currentS4, omicType = omicType)
      setFilterPlotFn("qc-progress-plot")
      "qc-progress-plot"
    },
    reportFinalizeSuccessFn = function(currentS4, history, output, ...) {
      captured$success <- list(currentS4 = currentS4, history = history)
      output$finalize_results <- "finalized"
      invisible("finalized")
    },
    reportFinalizeErrorFn = function(error, ...) {
      captured$error <- conditionMessage(error)
      invisible(captured$error)
    }
  )

  expect_identical(workflow$tab_status$quality_control, "complete")
  expect_identical(workflow$tab_status$normalization, "pending")
  expect_identical(captured$filter_plot, "qc-progress-plot")
  expect_identical(captured$tracking$omicType, "lipidomics")
  expect_true("lipid_qc_complete" %in% names(workflow$state_manager$saved()))
  expect_identical(output$finalize_results, "finalized")
  expect_true(any(captured$success$history == "lipid_qc_complete"))

  invalid_workflow <- module_ci_lipid_qc_workflow(current)
  invalid_workflow$state_manager <- module_ci_lipid_qc_state_manager(list(not = "s4"))
  runLipidQcS4FinalizeWorkflow(
    workflowData = invalid_workflow,
    omicType = "lipidomics",
    filterPlot = function(...) NULL,
    output = new.env(parent = emptyenv()),
    reportFinalizeSuccessFn = function(...) stop("should not finalize"),
    reportFinalizeErrorFn = function(error, ...) {
      captured$invalid_error <- conditionMessage(error)
      invisible(NULL)
    }
  )
  expect_match(captured$invalid_error, "Current state is not a LipidomicsAssayData object", fixed = TRUE)
  expect_identical(invalid_workflow$tab_status$quality_control, "pending")

  empty_design <- current
  empty_design@design_matrix <- empty_design@design_matrix[0, , drop = FALSE]
  expect_error(
    validateLipidQcS4FinalizeState(empty_design, reqFn = function(value) value),
    "no design matrix rows",
    fixed = TRUE
  )

  empty_assay <- current
  empty_assay@lipid_data$LCMS_Pos <- empty_assay@lipid_data$LCMS_Pos[0, , drop = FALSE]
  expect_error(
    validateLipidQcS4FinalizeState(empty_assay, reqFn = function(value) value),
    "no assay data",
    fixed = TRUE
  )
})

test_that("MCI-020.6 browser smoke exposes every lipid QC tab control", {
  ui_fragments <- list(
    orchestrator = htmltools::renderTags(mod_lipid_qc_ui("lipid-qc"))$html,
    intensity = htmltools::renderTags(mod_lipid_qc_intensity_ui("lipid-qc-intensity"))$html,
    duplicates = htmltools::renderTags(mod_lipid_qc_duplicates_ui("lipid-qc-duplicates"))$html,
    itsd = htmltools::renderTags(mod_lipid_qc_itsd_ui("lipid-qc-itsd"))$html,
    finalize = htmltools::renderTags(mod_lipid_qc_s4_ui("lipid-qc-s4"))$html
  )

  expect_match(ui_fragments$orchestrator, "dynamic_qc_tabs", fixed = TRUE)
  expect_match(ui_fragments$intensity, "apply_filter", fixed = TRUE)
  expect_match(ui_fragments$intensity, "revert_filter", fixed = TRUE)
  expect_match(ui_fragments$duplicates, "detect_duplicates", fixed = TRUE)
  expect_match(ui_fragments$duplicates, "resolve_duplicates", fixed = TRUE)
  expect_match(ui_fragments$itsd, "analyze_is", fixed = TRUE)
  expect_match(ui_fragments$itsd, "is_pattern", fixed = TRUE)
  expect_match(ui_fragments$finalize, "finalize_qc", fixed = TRUE)
  expect_match(ui_fragments$finalize, "assay_stats_table", fixed = TRUE)
})
