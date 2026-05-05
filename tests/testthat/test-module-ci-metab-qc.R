test_that("MCI-014.1 intensity filter matrix covers thresholds, zeros, missingness, and per-assay effects", {
  assay <- module_ci_metab_qc_assay(include_duplicates = FALSE, include_itsd = FALSE)

  strict <- metaboliteIntensityFilteringHelper(
    assay_table = assay,
    min_metabolite_intensity_threshold = 50,
    metabolites_proportion_of_samples_below_cutoff = 0.5,
    metabolite_id_column = "database_identifier"
  )
  expect_identical(module_ci_metab_qc_ids(strict), c("M_keep", "M_missing"))
  expect_false("M_boundary" %in% module_ci_metab_qc_ids(strict))
  expect_false("M_zero" %in% module_ci_metab_qc_ids(strict))

  relaxed <- metaboliteIntensityFilteringHelper(
    assay_table = assay,
    min_metabolite_intensity_threshold = 50,
    metabolites_proportion_of_samples_below_cutoff = 0.51,
    metabolite_id_column = "database_identifier"
  )
  expect_true("M_boundary" %in% module_ci_metab_qc_ids(relaxed))

  all_pass <- metaboliteIntensityFilteringHelper(
    assay_table = assay,
    min_metabolite_intensity_threshold = -1,
    metabolites_proportion_of_samples_below_cutoff = 1,
    metabolite_id_column = "database_identifier"
  )
  expect_equal(nrow(all_pass), nrow(assay))

  all_fail <- metaboliteIntensityFilteringHelper(
    assay_table = assay,
    min_metabolite_intensity_threshold = 10000,
    metabolites_proportion_of_samples_below_cutoff = 0.01,
    metabolite_id_column = "database_identifier"
  )
  expect_equal(nrow(all_fail), 0L)

  s4 <- module_ci_metab_qc_object(layout = "combined", include_duplicates = FALSE, include_itsd = FALSE)
  filtered <- metaboliteIntensityFiltering(
    theObject = s4,
    metabolites_intensity_cutoff_percentile = 20,
    metabolites_proportion_of_samples_below_cutoff = 0.5
  )
  expect_s4_class(filtered, "MetaboliteAssayData")
  expect_identical(names(filtered@metabolite_data), c("LCMS_Pos", "GCMS"))
  expect_true(all(module_ci_metab_qc_feature_counts(filtered) <= module_ci_metab_qc_feature_counts(s4)))
})

test_that("MCI-014.2 duplicate matrix covers IDs, annotations, ties, conflicts, and aggregation stats", {
  assay <- module_ci_metab_qc_assay(include_duplicates = TRUE, include_itsd = FALSE)

  detected <- findMetabDuplicateFeatureIDs(module_ci_metab_qc_object(layout = "gc", include_duplicates = TRUE, include_itsd = FALSE))
  expect_true("M_DUP" %in% detected$GCMS$database_identifier)
  expect_true("M_TIE" %in% detected$GCMS$database_identifier)

  resolved_ids <- resolveDuplicateFeaturesByIntensity(
    assay_tibble = assay,
    id_col = "database_identifier",
    sample_cols = module_ci_metab_qc_samples()
  )
  expect_false(anyDuplicated(resolved_ids$database_identifier) > 0L)
  expect_identical(
    resolved_ids$S1[resolved_ids$database_identifier == "M_DUP"],
    90
  )
  expect_identical(
    resolved_ids$S1[resolved_ids$database_identifier == "M_TIE"],
    50
  )

  resolved_annotations <- resolveDuplicateFeaturesByIntensity(
    assay_tibble = assay,
    id_col = "metabolite",
    sample_cols = module_ci_metab_qc_samples()
  )
  same_annotation <- resolved_annotations[resolved_annotations$metabolite == "SameAnnotation", , drop = FALSE]
  expect_equal(nrow(same_annotation), 1L)
  expect_identical(same_annotation$database_identifier, "M_ANN_B")

  dispatched <- resolveMetabDuplicateAssayData(
    assayList = list(LCMS_Pos = assay, GCMS = assay),
    metaboliteIdCol = "database_identifier",
    logWarnFn = function(...) invisible(NULL)
  )
  expect_identical(names(dispatched$resolvedAssayList), c("LCMS_Pos", "GCMS"))
  expect_true(all(vapply(dispatched$statsList, `[[`, numeric(1), "removed") >= 2))
})

test_that("MCI-014.3 ITSD matrix covers absent, valid, multi-assay, no-match, and malformed metadata", {
  no_itsd <- module_ci_metab_qc_object(layout = "gc", include_itsd = FALSE)
  expect_error(
    analyzeMetabQcItsdData(no_itsd, inputPattern = "^IS_"),
    "No internal standards found",
    fixed = TRUE
  )

  valid <- analyzeMetabQcItsdData(
    module_ci_metab_qc_object(layout = "combined", include_itsd = TRUE),
    inputPattern = "^IS_",
    logInfoFn = function(...) invisible(NULL)
  )
  expect_identical(valid$pattern, "^IS_")
  expect_identical(names(valid$metricsByAssay), c("LCMS_Pos", "GCMS"))
  expect_true(valid$nIsTotal >= 4L)
  expect_true(all(valid$metrics$cv < 5))
  expect_true(all(c("LCMS_Pos", "GCMS") %in% unique(valid$longData$assay)))

  no_match <- module_ci_metab_qc_object(layout = "gc", include_itsd = TRUE)
  fallback <- analyzeMetabQcItsdData(
    no_match,
    inputPattern = "^NO_MATCH$",
    logInfoFn = function(...) invisible(NULL)
  )
  expect_identical(fallback$pattern, "^NO_MATCH$")
  expect_true(fallback$nIsTotal >= 2L)

  malformed <- module_ci_metab_qc_object(layout = "gc", include_itsd = TRUE)
  malformed@metabolite_data$GCMS <- malformed@metabolite_data$GCMS[
    ,
    c("database_identifier", module_ci_metab_qc_samples()),
    drop = FALSE
  ]
  malformed@metabolite_id_column <- "missing_id"
  malformed@annotation_id_column <- "missing_annotation"
  expect_error(
    analyzeMetabQcItsdData(malformed, inputPattern = "^IS_"),
    "Columns searched:",
    fixed = TRUE
  )
})

test_that("MCI-014.4 S4 plotting matrix covers PCA, Pearson, RLE, density, small-n, constants, and missing groups", {
  s4 <- module_ci_metab_qc_object(layout = "gc", include_duplicates = FALSE, include_itsd = FALSE)

  pca <- suppressMessages(plotPca(s4, grouping_variable = "group", title = "Module CI PCA"))
  expect_true(is.list(pca))
  expect_named(pca, "GCMS")

  rle <- suppressMessages(plotRle(s4, grouping_variable = "group"))
  expect_true(is.list(rle))
  expect_named(rle, "GCMS")

  density <- suppressWarnings(suppressMessages(plotDensity(s4, grouping_variable = "group", title = "Module CI Density")))
  expect_true(is.list(density))

  pearson <- suppressWarnings(suppressMessages(plotPearson(s4, tech_rep_remove_regex = "", correlation_group = "tech_rep_group")))
  expect_true(is.list(pearson))

  small_n <- calculateMetaboliteCVs(
    assay_data = module_ci_metab_qc_assay(include_duplicates = FALSE, include_itsd = FALSE),
    design_matrix = module_ci_metab_qc_design(small_n = TRUE),
    group_id_col = "group",
    replicate_id_col = "tech_rep_group",
    sample_id_col = "Run",
    metabolite_id_col = "database_identifier",
    sample_columns = module_ci_metab_qc_samples()
  )
  expect_true(all(is.na(small_n$cv)))

  constant <- module_ci_metab_qc_object(layout = "gc", include_duplicates = FALSE, include_itsd = FALSE, constant = TRUE)
  constant_pearson <- suppressWarnings(plotPearson(constant, tech_rep_remove_regex = "", correlation_group = "tech_rep_group"))
  expect_true(is.list(constant_pearson))

  expect_error(
    plotPca(s4, grouping_variable = "missing_group"),
    "not found in design_matrix",
    fixed = TRUE
  )
})

test_that("MCI-014.5 finalization validates current S4 state and unlocks normalization handoff", {
  current <- module_ci_metab_qc_object(layout = "combined")
  workflow <- module_ci_metab_qc_workflow(current)
  output <- new.env(parent = emptyenv())
  captured <- new.env(parent = emptyenv())

  runMetabQcS4FinalizeWorkflow(
    workflowData = workflow,
    omicType = "metabolomics",
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
  expect_identical(captured$filter_plot, "qc-progress-plot")
  expect_identical(captured$tracking$omicType, "metabolomics")
  expect_true("metab_qc_complete" %in% names(workflow$state_manager$saved()))
  expect_identical(output$finalize_results, "finalized")
  expect_true(any(captured$success$history == "metab_qc_complete"))

  invalid_workflow <- module_ci_metab_qc_workflow(current)
  invalid_workflow$state_manager <- module_ci_metab_qc_state_manager(list(not = "s4"))
  runMetabQcS4FinalizeWorkflow(
    workflowData = invalid_workflow,
    omicType = "metabolomics",
    filterPlot = function(...) NULL,
    output = new.env(parent = emptyenv()),
    reportFinalizeSuccessFn = function(...) stop("should not finalize"),
    reportFinalizeErrorFn = function(error, ...) {
      captured$invalid_error <- conditionMessage(error)
      invisible(NULL)
    }
  )
  expect_match(captured$invalid_error, "Current state is not a MetaboliteAssayData object", fixed = TRUE)
  expect_identical(invalid_workflow$tab_status$quality_control, "pending")

  empty_design <- current
  empty_design@design_matrix <- empty_design@design_matrix[0, , drop = FALSE]
  expect_error(
    validateMetabQcS4FinalizeState(empty_design, reqFn = function(value) value),
    "no design matrix rows",
    fixed = TRUE
  )
})

test_that("MCI-014.6 browser smoke exposes every metabolomics QC tab control", {
  ui_fragments <- list(
    orchestrator = htmltools::renderTags(mod_metab_qc_ui("metab-qc"))$html,
    intensity = htmltools::renderTags(mod_metab_qc_intensity_ui("metab-qc-intensity"))$html,
    duplicates = htmltools::renderTags(mod_metab_qc_duplicates_ui("metab-qc-duplicates"))$html,
    itsd = htmltools::renderTags(mod_metab_qc_itsd_ui("metab-qc-itsd"))$html,
    finalize = htmltools::renderTags(mod_metab_qc_s4_ui("metab-qc-s4"))$html
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
