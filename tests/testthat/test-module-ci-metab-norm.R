test_that("MCI-015.1 normalization method matrix preserves assay names, samples, and edge values", {
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
  expect_true("metab_normalized" %in% names(workflow$state_manager$saved()))
  module_ci_metab_norm_assert_alignment(shell_state$currentS4, c("LCMS_Pos", "GCMS"))
})

test_that("MCI-015.2 ITSD matrix serializes regex, manual, absent, and invalid ITSD choices", {
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
  expect_true("metab_itsd_norm" %in% names(workflow$state_manager$saved()))
  expect_identical(norm_data$post_itsd_obj@args$ITSDNormalization$itsd_features_per_assay$LCMS_Pos, c("IS_A", "IS_B"))

  no_itsd <- module_ci_metab_norm_object(layout = "gc", include_itsd = FALSE, positive_only = TRUE)
  no_itsd_warnings <- character()
  no_itsd_norm <- withCallingHandlers(
    suppressMessages(normaliseUntransformedData(no_itsd, method = "ITSD")),
    warning = function(w) {
      no_itsd_warnings <<- c(no_itsd_warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  expect_true(any(grepl("No ITSD features identified", no_itsd_warnings, fixed = TRUE)))
  expect_identical(
    module_ci_metab_norm_feature_ids(no_itsd_norm, "GCMS"),
    module_ci_metab_norm_feature_ids(no_itsd, "GCMS")
  )

  invalid_manual <- suppressWarnings(suppressMessages(normaliseUntransformedData(
    object,
    method = "ITSD",
    itsd_feature_ids = list(LCMS_Pos = "NO_MATCH", GCMS = "NO_MATCH")
  )))
  expect_true("IS_A" %in% module_ci_metab_norm_feature_ids(invalid_manual, "LCMS_Pos"))
  expect_error(
    normaliseUntransformedData(object, method = "ITSD", itsd_aggregation = "mode"),
    "must be one of"
  )
})

test_that("MCI-015.3 RUV matrix covers skip, manual, automatic, invalid k, controls, and per-assay payloads", {
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

  manual_norm_data <- module_ci_metab_norm_data(object)
  manual_opt <- runMetabNormRuvOptimizationStep(
    currentS4 = object,
    ruvMode = "manual",
    autoPercentageMin = 5,
    autoPercentageMax = 10,
    ruvGroupingVariable = "tech_rep_group",
    separationMetric = "max_difference",
    kPenaltyWeight = 0.5,
    adaptiveKPenalty = FALSE,
    manualK = 2L,
    manualPercentage = 40,
    experimentPaths = list(metabolite_qc_dir = tempdir()),
    normData = manual_norm_data,
    runPerAssayRuvOptimizationFn = optimizer_stub
  )
  expect_named(manual_opt$bestKPerAssay, c("LCMS_Pos", "GCMS"))
  expect_identical(unname(unlist(manual_opt$bestKPerAssay)), c(2L, 2L))
  expect_true(all(vapply(manual_opt$ctrlPerAssay, function(x) {
    is.logical(x) && !is.null(names(x)) && sum(x) == 6L
  }, logical(1))))

  auto_opt <- runMetabNormRuvOptimizationStep(
    currentS4 = object,
    ruvMode = "automatic",
    autoPercentageMin = 10,
    autoPercentageMax = 12,
    ruvGroupingVariable = "tech_rep_group",
    separationMetric = "max_difference",
    kPenaltyWeight = 0.5,
    adaptiveKPenalty = TRUE,
    manualK = 2L,
    manualPercentage = 40,
    experimentPaths = list(metabolite_qc_dir = tempdir()),
    normData = module_ci_metab_norm_data(object),
    runPerAssayRuvOptimizationFn = optimizer_stub
  )
  expect_identical(unname(unlist(auto_opt$bestKPerAssay)), c(1L, 2L))
  expect_equal(auto_opt$ruvResults$LCMS_Pos$best_percentage, 13)
  expect_s3_class(auto_opt$ruvResults$GCMS$optimization_results, "data.frame")

  failed_results <- module_ci_metab_norm_ruv_results(object)
  failed_results$LCMS_Pos <- list(
    success = FALSE,
    best_k = NA_integer_,
    best_percentage = NA_real_,
    control_genes_index = NULL,
    error = "too few controls"
  )
  failed_opt <- runMetabNormRuvOptimizationStep(
    currentS4 = object,
    ruvMode = "manual",
    autoPercentageMin = 5,
    autoPercentageMax = 10,
    ruvGroupingVariable = "tech_rep_group",
    separationMetric = "max_difference",
    kPenaltyWeight = 0.5,
    adaptiveKPenalty = FALSE,
    manualK = 2L,
    manualPercentage = 1,
    experimentPaths = list(metabolite_qc_dir = tempdir()),
    normData = module_ci_metab_norm_data(object),
    runPerAssayRuvOptimizationFn = function(...) failed_results
  )
  expect_true(is.na(failed_opt$bestKPerAssay$LCMS_Pos))
  expect_null(failed_opt$ctrlPerAssay$LCMS_Pos)
  expect_identical(failed_opt$bestKPerAssay$GCMS, 2L)

  correction_workflow <- module_ci_metab_norm_workflow(object)
  correction_norm_data <- module_ci_metab_norm_data(object)
  correction <- runMetabNormRuvCorrectionStep(
    currentS4 = object,
    ruvGroupingVariable = "tech_rep_group",
    bestKPerAssay = auto_opt$bestKPerAssay,
    ctrlPerAssay = auto_opt$ctrlPerAssay,
    workflowData = correction_workflow,
    normData = correction_norm_data,
    ruvIII_C_VaryingFn = function(theObject, ruv_grouping_variable, ruv_number_k, ctrl) {
      theObject@args$ruvIII_C_Varying <- list(
        ruv_grouping_variable = ruv_grouping_variable,
        ruv_number_k = ruv_number_k,
        ctrl = ctrl
      )
      theObject
    }
  )
  expect_s4_class(correction$currentS4, "MetaboliteAssayData")
  expect_true(correction_norm_data$ruv_complete)
  expect_true("metab_ruv_corrected" %in% names(correction_workflow$state_manager$saved()))
  expect_identical(correction$currentS4@args$ruvIII_C_Varying$ruv_number_k$LCMS_Pos, 1L)

  invalid_workflow <- module_ci_metab_norm_workflow(object)
  invalid_norm_data <- module_ci_metab_norm_data(object)
  expect_error(
    runMetabNormRuvCorrectionStep(
      currentS4 = object,
      ruvGroupingVariable = "tech_rep_group",
      bestKPerAssay = list(LCMS_Pos = 0L, GCMS = 2L),
      ctrlPerAssay = module_ci_metab_norm_ctrl(object),
      workflowData = invalid_workflow,
      normData = invalid_norm_data,
      ruvIII_C_VaryingFn = function(...) stop("invalid k")
    ),
    "invalid k"
  )
  expect_false(invalid_norm_data$ruv_complete)
  expect_false("metab_ruv_corrected" %in% names(invalid_workflow$state_manager$saved()))
})

test_that("MCI-015.4 correlation filtering matrix covers skip, pass/fail, boundary, small-n, and invalid inputs", {
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
  expect_true("metab_norm_complete" %in% names(workflow$state_manager$saved()))

  run_apply_case <- function(case) {
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
      reqFn = function(value) value,
      calculateCorrelationsFn = function(...) {
        module_ci_metab_norm_correlation_pairs(case, names(object@metabolite_data))
      }
    ))) -> result
    list(result = result, workflow = case_workflow, norm_data = case_norm_data)
  }

  pass_all <- run_apply_case("pass_all")
  expect_identical(pass_all$result$status, "success")
  expect_setequal(pass_all$norm_data$correlation_filtered_obj@design_matrix$Run, c("S1", "S2", "S3", "S4"))
  expect_true("metab_correlation_filtered" %in% names(pass_all$workflow$state_manager$saved()))

  boundary <- run_apply_case("boundary")
  expect_identical(boundary$result$status, "success")
  expect_setequal(boundary$norm_data$correlation_filtered_obj@design_matrix$Run, c("S1", "S2", "S3", "S4"))

  fail_one <- run_apply_case("fail_one")
  expect_identical(fail_one$result$status, "success")
  expect_identical(fail_one$norm_data$correlation_filtered_obj@design_matrix$Run, c("S1", "S2"))
  expect_identical(module_ci_metab_norm_sample_cols(fail_one$norm_data$correlation_filtered_obj, "LCMS_Pos"), c("S1", "S2"))

  fail_many <- suppressMessages(filterSamplesByMetaboliteCorrelationThreshold(
    object,
    pearson_correlation_per_pair = module_ci_metab_norm_correlation_pairs("fail_many", names(object@metabolite_data)),
    min_pearson_correlation_threshold = 0.75
  ))
  expect_identical(module_ci_metab_norm_sample_cols(fail_many, "LCMS_Pos"), character())
  expect_equal(nrow(fail_many@design_matrix), 0L)

  small_n_warnings <- character()
  small_n <- withCallingHandlers(
    suppressMessages(filterSamplesByMetaboliteCorrelationThreshold(
      object,
      pearson_correlation_per_pair = module_ci_metab_norm_correlation_pairs("small_n", names(object@metabolite_data)),
      min_pearson_correlation_threshold = 0.75
    )),
    warning = function(w) {
      small_n_warnings <<- c(small_n_warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  expect_true(any(grepl("No correlation results provided for assay", small_n_warnings, fixed = TRUE)))
  expect_setequal(small_n@design_matrix$Run, c("S1", "S2", "S3", "S4"))

  invalid_workflow <- module_ci_metab_norm_workflow(object)
  invalid_norm_data <- module_ci_metab_norm_data(object)
  invalid_state <- dispatchMetabNormApplyCorrelation(
    workflowData = invalid_workflow,
    normData = invalid_norm_data,
    observerState = list(
      currentS4 = object,
      threshold = 0.75,
      groupingVariable = "tech_rep_group",
      notificationId = "invalid_corr"
    ),
    showNotificationFn = function(...) invisible(NULL),
    removeNotificationFn = function(...) invisible(NULL),
    reqFn = function(value) value,
    calculateCorrelationsFn = function(...) 1,
    logErrorFn = function(...) invisible(NULL)
  )
  expect_identical(invalid_state$status, "error")
  expect_false(invalid_norm_data$correlation_filtering_complete)
  expect_false("metab_correlation_filtered" %in% names(invalid_workflow$state_manager$saved()))
})

test_that("MCI-015.5 session export matrix preserves metabolomics, RUV, ITSD, and feature-count payloads", {
  session_data <- module_ci_metab_norm_export_variant("automatic")

  expect_identical(session_data$omic_type, "metabolomics")
  expect_s4_class(session_data$current_s4_object, "MetaboliteAssayData")
  expect_identical(session_data$assay_names, c("LCMS_Pos", "GCMS"))
  expect_true(session_data$normalization_complete)
  expect_true(session_data$ruv_complete)
  expect_true(session_data$correlation_filtering_complete)
  expect_identical(session_data$normalization_method, "quantile")
  expect_identical(session_data$ruv_mode, "automatic")
  expect_true(isTRUE(session_data$itsd_applied))
  expect_identical(session_data$itsd_aggregation, "median")
  expect_identical(session_data$current_s4_object@args$ITSDNormalization$itsd_counts_per_assay$LCMS_Pos, 2L)
  expect_identical(session_data$ruv_optimization_results$LCMS_Pos$best_k, 1L)
  expect_equal(session_data$feature_counts$LCMS_Pos$features, nrow(session_data$current_s4_object@metabolite_data$LCMS_Pos))

  source_dir <- tempfile("mci-015-metab-norm-export-")
  dir.create(source_dir, recursive = TRUE)
  fixed_time <- as.POSIXct("2026-05-05 12:34:56", tz = "UTC")
  files <- saveMetabNormExportSessionRdsFiles(
    sessionData = session_data,
    sourceDir = source_dir,
    timeFn = function() fixed_time,
    logInfoFn = function(...) invisible(NULL)
  )
  saveMetabNormExportMetadataFiles(
    sessionData = session_data,
    sourceDir = source_dir,
    logInfoFn = function(...) invisible(NULL),
    logWarnFn = function(...) invisible(NULL)
  )
  summary <- saveMetabNormExportSummaryFile(
    sessionData = session_data,
    sourceDir = source_dir,
    sessionFilename = files$sessionFilename,
    timeFn = function() fixed_time,
    logInfoFn = function(...) invisible(NULL)
  )

  expect_true(file.exists(files$sessionFilepath))
  expect_true(file.exists(files$latestFilepath))
  expect_true(file.exists(file.path(source_dir, "metab_ruv_optimization_results.RDS")))
  expect_true(file.exists(file.path(source_dir, "metab_itsd_selections.RDS")))
  expect_true(file.exists(file.path(source_dir, "metab_qc_params.RDS")))
  expect_true(file.exists(summary$summaryFilepath))
  expect_match(summary$summaryContent, "Metabolomics Normalized Session Data Export Summary", fixed = TRUE)
  expect_match(summary$summaryContent, "- Method: quantile", fixed = TRUE)
  expect_match(summary$summaryContent, "- RUV mode: automatic", fixed = TRUE)
  expect_match(summary$summaryContent, "LCMS_Pos: k=1", fixed = TRUE)
  expect_s4_class(readRDS(files$latestFilepath)$current_s4_object, "MetaboliteAssayData")
})

test_that("MCI-015.6 DA reload consumes every accepted metabolomics export variant", {
  for (variant in c("skip", "manual", "automatic")) {
    session_data <- module_ci_metab_norm_export_variant(variant)
    source_dir <- tempfile(paste0("mci-015-da-reload-", variant, "-"))
    dir.create(source_dir, recursive = TRUE)
    session_file <- file.path(source_dir, "metab_filtered_session_data_latest.rds")
    saveRDS(session_data, session_file)

    workflow <- module_ci_metab_norm_workflow(session_data$current_s4_object)
    da_data <- list2env(list(), parent = emptyenv())
    session <- shiny::MockShinySession$new()
    capture <- new.env(parent = emptyenv())
    capture$select <- list()
    capture$text <- list()

    result <- runMetabDaLoadSessionObserverShell(
      sessionFile = session_file,
      workflowData = workflow,
      daData = da_data,
      session = session,
      showNotification = function(...) invisible(NULL),
      removeNotification = function(...) invisible(NULL),
      logInfo = function(...) invisible(NULL),
      logError = function(...) invisible(NULL),
      restoreState = function(sessionData, sessionFile, workflowData, daData, session, debugLog) {
        restoreMetabDaLoadedSessionState(
          sessionData = sessionData,
          sessionFile = sessionFile,
          workflowData = workflowData,
          daData = daData,
          session = session,
          updateSelectInput = function(session, inputId, choices = NULL, selected = NULL, ...) {
            capture$select[[inputId]] <- list(choices = choices, selected = selected)
            invisible(NULL)
          },
          updateTextAreaInput = function(session, inputId, value = NULL, ...) {
            capture$text[[inputId]] <- value
            invisible(NULL)
          },
          logInfo = function(...) invisible(NULL),
          logWarn = function(...) invisible(NULL),
          debugLog = debugLog
        )
      }
    )

    expect_identical(result$status, "success")
    expect_s4_class(da_data$current_s4_object, "MetaboliteAssayData")
    expect_identical(da_data$assays_available, c("LCMS_Pos", "GCMS"))
    expect_identical(da_data$formula_from_s4, "~ 0 + group + batch")
    expect_identical(capture$text$formula_string, "~ 0 + group + batch")
    expect_identical(capture$select$volcano_contrast$choices, "KO_vs_WT")
    expect_identical(capture$select$table_assay$choices, c("All", "LCMS_Pos", "GCMS"))
    expect_true(session_data$r6_current_state_name %in% names(workflow$state_manager$saved()))
  }
})
