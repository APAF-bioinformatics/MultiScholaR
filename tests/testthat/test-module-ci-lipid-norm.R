library(testthat)

test_that("MCI-021.1 normalization method matrix preserves assay names, samples, and edge values", {
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
  expect_true("lipid_normalized" %in% names(pipeline$workflow$state_manager$saved()))
  module_ci_lipid_norm_assert_alignment(pipeline$norm_data$post_norm_obj, c("LCMS_Pos", "GCMS"))
})

test_that("MCI-021.2 ITSD matrix serializes regex, manual, absent, and invalid ITSD choices", {
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
  expect_true("lipid_itsd_norm" %in% names(workflow$state_manager$saved()))
  expect_identical(norm_data$post_itsd_obj@args$ITSDNormalization$itsd_features_per_assay$LCMS_Pos, c("IS_A", "IS_B"))

  no_itsd <- module_ci_lipid_norm_object(layout = "gc", include_itsd = FALSE, positive_only = TRUE)
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
    module_ci_lipid_norm_feature_ids(no_itsd_norm, "GCMS"),
    module_ci_lipid_norm_feature_ids(no_itsd, "GCMS")
  )

  invalid_manual <- suppressWarnings(suppressMessages(normaliseUntransformedData(
    object,
    method = "ITSD",
    itsd_feature_ids = list(LCMS_Pos = "NO_MATCH", GCMS = "NO_MATCH")
  )))
  expect_true("IS_A" %in% module_ci_lipid_norm_feature_ids(invalid_manual, "LCMS_Pos"))
  expect_error(
    normaliseUntransformedData(object, method = "ITSD", itsd_aggregation = "mode"),
    "must be one of"
  )
})

test_that("MCI-021.3 RUV matrix covers skip, manual, automatic, invalid k, controls, and per-assay payloads", {
  object <- module_ci_lipid_norm_object(layout = "combined", positive_only = TRUE)

  skip_pipeline <- suppressWarnings(suppressMessages(module_ci_lipid_norm_run_pipeline(
    input = module_ci_lipid_norm_input(norm_method = "none", ruv_mode = "skip"),
    object = object
  )))
  expect_true(skip_pipeline$result)
  expect_true(skip_pipeline$norm_data$ruv_complete)
  expect_identical(skip_pipeline$norm_data$ruv_corrected_obj, skip_pipeline$norm_data$post_norm_obj)
  expect_false("lipid_ruv_corrected" %in% names(skip_pipeline$workflow$state_manager$saved()))

  manual_pipeline <- suppressWarnings(suppressMessages(module_ci_lipid_norm_run_pipeline(
    input = module_ci_lipid_norm_input(norm_method = "none", ruv_mode = "manual"),
    object = object
  )))
  expect_true(manual_pipeline$result)
  expect_true(manual_pipeline$norm_data$ruv_complete)
  expect_true("lipid_ruv_corrected" %in% names(manual_pipeline$workflow$state_manager$saved()))
  expect_identical(manual_pipeline$norm_data$ruv_corrected_obj@args$ruv_number_k$LCMS_Pos, 2L)
  expect_true(all(vapply(manual_pipeline$norm_data$ruv_corrected_obj@args$ctrl, function(x) {
    is.logical(x) && !is.null(names(x)) && sum(x) == 6L
  }, logical(1))))

  auto_pipeline <- suppressWarnings(suppressMessages(module_ci_lipid_norm_run_pipeline(
    input = module_ci_lipid_norm_input(norm_method = "none", ruv_mode = "automatic"),
    object = object
  )))
  expect_true(auto_pipeline$result)
  expect_identical(auto_pipeline$norm_data$ruv_corrected_obj@args$ruv_number_k$LCMS_Pos, 1L)
  expect_identical(auto_pipeline$norm_data$ruv_corrected_obj@args$ruv_number_k$GCMS, 2L)
  expect_equal(auto_pipeline$norm_data$ruv_optimization_results$LCMS_Pos$best_percentage, 13)
  expect_s3_class(auto_pipeline$norm_data$ruv_optimization_results$GCMS$optimization_results, "data.frame")

  failed_results <- module_ci_lipid_norm_ruv_results(object)
  failed_results$LCMS_Pos <- list(
    success = FALSE,
    best_k = NA_integer_,
    best_percentage = NA_real_,
    control_genes_index = NULL,
    error = "too few controls"
  )
  expect_true(is.na(extractLipidBestKPerAssay(failed_results)$LCMS_Pos))
  expect_null(extractLipidCtrlPerAssay(failed_results)$LCMS_Pos)
  expect_identical(extractLipidBestKPerAssay(failed_results)$GCMS, 2L)

  invalid_workflow <- module_ci_lipid_norm_workflow(object)
  invalid_norm_data <- module_ci_lipid_norm_data(object)
  invalid_pipeline <- suppressWarnings(suppressMessages(module_ci_lipid_norm_run_pipeline(
    input = module_ci_lipid_norm_input(norm_method = "none", ruv_mode = "manual"),
    object = object,
    workflow = invalid_workflow,
    norm_data = invalid_norm_data,
    optimizer_fn = function(...) failed_results,
    ruv_fn = function(theObject, ruv_number_k, ...) {
      if (any(is.na(unlist(ruv_number_k))) || any(unlist(ruv_number_k) < 1)) {
        stop("invalid k")
      }
      theObject
    }
  )))
  expect_false(invalid_pipeline$result)
  expect_false(invalid_norm_data$ruv_complete)
  expect_false("lipid_ruv_corrected" %in% names(invalid_workflow$state_manager$saved()))
})

test_that("MCI-021.4 correlation filtering matrix covers skip, pass/fail, boundary, small-n, and invalid inputs", {
  object <- module_ci_lipid_norm_object(layout = "combined", positive_only = TRUE)
  norm_data <- module_ci_lipid_norm_data(object)
  norm_data$post_norm_obj <- object
  norm_data$normalization_complete <- TRUE
  workflow <- module_ci_lipid_norm_workflow(object)

  skip_state <- handleLipidNormSkipCorrelationFilter(
    workflowData = workflow,
    normData = norm_data,
    addLog = function(...) invisible(NULL),
    reqFn = function(value) value,
    showNotificationFn = function(...) invisible(NULL)
  )
  expect_true(skip_state)
  expect_true(norm_data$correlation_filtering_complete)
  expect_identical(norm_data$correlation_filtered_obj, object)
  expect_true("lipid_norm_complete" %in% names(workflow$state_manager$saved()))
  expect_identical(workflow$tab_status$normalization, "complete")

  run_apply_case <- function(case) {
    case_workflow <- module_ci_lipid_norm_workflow(object)
    case_norm_data <- module_ci_lipid_norm_data(object)
    case_norm_data$normalization_complete <- TRUE
    case_norm_data$post_norm_obj <- object
    suppressWarnings(suppressMessages(handleLipidNormApplyCorrelationFilter(
      input = module_ci_lipid_norm_input(correlation_threshold = 0.75),
      workflowData = case_workflow,
      normData = case_norm_data,
      addLog = function(...) invisible(NULL),
      reqFn = function(value) value,
      showNotificationFn = function(...) invisible(NULL),
      removeNotificationFn = function(...) invisible(NULL),
      pearsonCorForSamplePairsFn = function(...) {
        module_ci_lipid_norm_correlation_pairs(case, names(object@lipid_data))
      },
      logInfoFn = function(...) invisible(NULL),
      logErrorFn = function(...) invisible(NULL)
    ))) -> result
    list(result = result, workflow = case_workflow, norm_data = case_norm_data)
  }

  pass_all <- run_apply_case("pass_all")
  expect_true(pass_all$result)
  expect_setequal(pass_all$norm_data$correlation_filtered_obj@design_matrix$Run, c("S1", "S2", "S3", "S4"))
  expect_true("lipid_correlation_filtered" %in% names(pass_all$workflow$state_manager$saved()))

  boundary <- run_apply_case("boundary")
  expect_true(boundary$result)
  expect_setequal(boundary$norm_data$correlation_filtered_obj@design_matrix$Run, c("S1", "S2", "S3", "S4"))

  fail_one <- run_apply_case("fail_one")
  expect_true(fail_one$result)
  expect_identical(fail_one$norm_data$correlation_filtered_obj@design_matrix$Run, c("S1", "S2"))
  expect_identical(module_ci_lipid_norm_sample_cols(fail_one$norm_data$correlation_filtered_obj, "LCMS_Pos"), c("S1", "S2"))

  fail_many <- suppressMessages(filterSamplesByLipidCorrelationThreshold(
    object,
    pearson_correlation_per_pair = module_ci_lipid_norm_correlation_pairs("fail_many", names(object@lipid_data)),
    min_pearson_correlation_threshold = 0.75
  ))
  expect_identical(module_ci_lipid_norm_sample_cols(fail_many, "LCMS_Pos"), character())
  expect_equal(nrow(fail_many@design_matrix), 0L)

  small_n_warnings <- character()
  small_n <- withCallingHandlers(
    suppressMessages(filterSamplesByLipidCorrelationThreshold(
      object,
      pearson_correlation_per_pair = module_ci_lipid_norm_correlation_pairs("small_n", names(object@lipid_data)),
      min_pearson_correlation_threshold = 0.75
    )),
    warning = function(w) {
      small_n_warnings <<- c(small_n_warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  expect_true(any(grepl("No correlation results provided for assay", small_n_warnings, fixed = TRUE)))
  expect_setequal(small_n@design_matrix$Run, c("S1", "S2", "S3", "S4"))

  invalid_workflow <- module_ci_lipid_norm_workflow(object)
  invalid_norm_data <- module_ci_lipid_norm_data(object)
  invalid_norm_data$normalization_complete <- TRUE
  invalid_norm_data$post_norm_obj <- object
  invalid_state <- handleLipidNormApplyCorrelationFilter(
    input = module_ci_lipid_norm_input(correlation_threshold = 0.75),
    workflowData = invalid_workflow,
    normData = invalid_norm_data,
    addLog = function(...) invisible(NULL),
    reqFn = function(value) value,
    showNotificationFn = function(...) invisible(NULL),
    removeNotificationFn = function(...) invisible(NULL),
    pearsonCorForSamplePairsFn = function(...) 1,
    logInfoFn = function(...) invisible(NULL),
    logErrorFn = function(...) invisible(NULL)
  )
  expect_false(invalid_state)
  expect_false(invalid_norm_data$correlation_filtering_complete)
  expect_false("lipid_correlation_filtered" %in% names(invalid_workflow$state_manager$saved()))
})

test_that("MCI-021.5 session export matrix preserves lipidomics, RUV, ITSD, and feature-count payloads", {
  export <- module_ci_lipid_norm_export_variant("automatic")
  session_data <- export$session_data

  expect_identical(session_data$omic_type, "lipidomics")
  expect_s4_class(session_data$current_s4_object, "LipidomicsAssayData")
  expect_identical(session_data$assay_names, c("LCMS_Pos", "GCMS"))
  expect_true(session_data$normalization_complete)
  expect_true(session_data$ruv_complete)
  expect_true(session_data$correlation_filtering_complete)
  expect_identical(session_data$normalization_method, "quantile")
  expect_identical(session_data$ruv_mode, "automatic")
  expect_true(isTRUE(session_data$itsd_applied))
  expect_identical(session_data$itsd_aggregation, "median")
  expect_identical(session_data$current_s4_object@args$ITSDNormalization$itsd_counts_per_assay$LCMS_Pos, 2L)
  expect_identical(session_data$current_s4_object@args$ruv_number_k$LCMS_Pos, 1L)
  expect_identical(session_data$ruv_optimization_results$LCMS_Pos$best_k, 1L)
  expect_equal(session_data$feature_counts$LCMS_Pos$features, nrow(session_data$current_s4_object@lipid_data$LCMS_Pos))
  expect_equal(session_data$feature_counts$LCMS_Pos$samples, 4L)

  expect_true(file.exists(export$artifacts$session_filepath))
  expect_true(file.exists(export$latest_file))
  expect_true(file.exists(file.path(export$source_dir, "lipid_ruv_optimization_results.RDS")))
  expect_true(file.exists(file.path(export$source_dir, "lipid_itsd_selections.RDS")))
  expect_true(file.exists(file.path(export$source_dir, "lipid_qc_params.RDS")))
  expect_true(file.exists(export$summary_file))
  summary_text <- paste(export$summary, collapse = "\n")
  expect_match(summary_text, "Lipidomics Normalized Session Data Export Summary", fixed = TRUE)
  expect_match(summary_text, "- Method: quantile", fixed = TRUE)
  expect_match(summary_text, "- RUV mode: automatic", fixed = TRUE)
  expect_match(summary_text, "LCMS_Pos: k=1", fixed = TRUE)
  expect_s4_class(readRDS(export$latest_file)$current_s4_object, "LipidomicsAssayData")
})

test_that("MCI-021.6 DA reload consumes every accepted lipid export variant", {
  for (variant in c("skip", "manual", "automatic")) {
    export <- module_ci_lipid_norm_export_variant(variant)
    workflow <- module_ci_lipid_norm_workflow(export$session_data$current_s4_object)
    da_data <- list2env(list(), parent = emptyenv())
    session <- shiny::MockShinySession$new()
    capture <- new.env(parent = emptyenv())
    capture$select <- list()
    capture$text <- list()

    update_select <- function(session, inputId, choices = NULL, selected = NULL, ...) {
      capture$select[[inputId]] <- list(choices = choices, selected = selected)
      invisible(NULL)
    }
    update_text <- function(session, inputId, value = NULL, ...) {
      capture$text[[inputId]] <- value
      invisible(NULL)
    }

    result <- loadLipidDaSessionFromFile(
      sessionFile = export$latest_file,
      workflowData = workflow,
      daData = da_data,
      session = session,
      showLoadingFn = function() list(notificationId = "loading_session"),
      finalizeSuccessFn = function(...) list(successMessage = "Session loaded successfully!"),
      finalizeErrorFn = function(errorMessage, ...) stop(errorMessage),
      restorePostReadFn = function(sessionData, workflowData, daData, session, debugLog) {
        restoreLipidDaPostReadSessionState(
          sessionData = sessionData,
          workflowData = workflowData,
          daData = daData,
          session = session,
          debugLog = debugLog,
          restoreContrastsFn = function(sessionData, workflowData, daData, session, ...) {
            restoreLipidDaContrastsFromSession(
              sessionData = sessionData,
              workflowData = workflowData,
              daData = daData,
              session = session,
              updateSelectInputFn = update_select,
              logInfo = function(...) invisible(NULL)
            )
          },
          restoreAssaysFn = function(sessionData, daData, session, debugLog, ...) {
            restoreLipidDaAssaysFromSession(
              sessionData = sessionData,
              daData = daData,
              session = session,
              updateSelectInputFn = update_select,
              debugLog = debugLog
            )
          },
          restoreFormulaFn = function(sessionData, daData, session, debugLog, ...) {
            restoreLipidDaFormulaFromSession(
              sessionData = sessionData,
              daData = daData,
              session = session,
              updateTextAreaInputFn = update_text,
              debugLog = debugLog,
              hasArgsMethodFn = function() TRUE,
              logWarn = function(...) invisible(NULL)
            )
          }
        )
      }
    )

    expect_identical(result$status, "success")
    expect_s4_class(da_data$current_s4_object, "LipidomicsAssayData")
    expect_identical(da_data$assays_available, c("LCMS_Pos", "GCMS"))
    expect_identical(da_data$formula_from_s4, "~ 0 + group + batch")
    expect_identical(capture$text$formula_string, "~ 0 + group + batch")
    expect_identical(capture$select$volcano_contrast$choices, "KO_vs_WT")
    expect_identical(capture$select$table_assay$choices, c("All", "LCMS_Pos", "GCMS"))
    expect_true(export$session_data$r6_current_state_name %in% names(workflow$state_manager$saved()))
  }
})
