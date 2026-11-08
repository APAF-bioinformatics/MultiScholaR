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

  # KIT-105 W5: the between-sample step through the tab's in-process normalisation, the core's
  # runWorkflowOmicsInProcessNormalisation (runMetabNormBetweenSampleProgressApplyShell was removed with
  # mod_metab_norm_pipeline_helpers.R).
  workflow <- module_ci_metab_norm_workflow(module_ci_metab_norm_object(positive_only = TRUE))
  norm_data <- module_ci_metab_norm_data(workflow$state_manager$getState())
  pipeline <- suppressWarnings(suppressMessages(module_ci_metab_norm_run_pipeline(
    input = module_ci_metab_norm_ruv_input(module_ci_metab_norm_input(norm_method = "quantile", ruv_mode = "skip")),
    workflow = workflow,
    norm_data = norm_data
  )))

  expect_true(pipeline$result)
  expect_true(isSavedWorkflowReference(norm_data$post_norm_obj)) # KIT-075: the step carries the state's reference
  expect_true(norm_data$normalization_complete)
  expect_true("metab_normalized" %in% workflow$state_manager$getHistory())
  module_ci_metab_norm_assert_alignment(module_ci_metab_norm_state(workflow, norm_data$post_norm_obj), c("LCMS_Pos", "GCMS"))
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
  # KIT-105 W5: the manual selection is resolved, and ITSD applied, by the tab's in-process normalisation, the core's
  # runWorkflowOmicsInProcessNormalisation (D8; resolveMetabNormManualItsdFeatureIds and
  # runMetabNormItsdProgressApplyShell were removed with mod_metab_norm_pipeline_helpers.R). The ids it resolves are
  # the ones the ITSD step receives.
  workflow <- module_ci_metab_norm_workflow(object)
  norm_data <- module_ci_metab_norm_data(object)
  norm_data$itsd_selections <- manual_indices
  itsd_state <- suppressWarnings(suppressMessages(module_ci_metab_norm_run_pipeline(
    input = module_ci_metab_norm_ruv_input(module_ci_metab_norm_input(apply_itsd = TRUE, itsd_aggregation = "mean")),
    object = object,
    workflow = workflow,
    norm_data = norm_data,
    itsd_table_fn = function(assay_data, metabolite_id_col, annotation_cols) {
      data.frame(feature_id = assay_data[[metabolite_id_col]], stringsAsFactors = FALSE)
    }
  )))
  manual_ids <- itsd_state$itsd_feature_ids

  expect_identical(manual_ids$LCMS_Pos, c("IS_A", "IS_B"))
  expect_identical(manual_ids$GCMS, c("IS_A", "IS_B"))

  expect_true(itsd_state$result)
  expect_true("ITSD normalization complete" %in% itsd_state$logs)
  expect_true(isSavedWorkflowReference(norm_data$post_itsd_obj)) # KIT-075: the state's reference
  expect_true("metab_itsd_norm" %in% workflow$state_manager$getHistory())
  expect_identical(module_ci_metab_norm_state(workflow, norm_data$post_itsd_obj)@args$ITSDNormalization$itsd_features_per_assay$LCMS_Pos, c("IS_A", "IS_B"))

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

  # KIT-105 W5: RUV-III through the tab's in-process normalisation, the core's runWorkflowOmicsInProcessNormalisation
  # (runMetabNormRuvProgressApplyShell was removed with mod_metab_norm_pipeline_helpers.R). The step's input is the
  # normalised state the run carries, so the skip keeps that state; the k per assay is the correction's argument.
  skip_state <- suppressWarnings(suppressMessages(module_ci_metab_norm_run_pipeline(
    input = module_ci_metab_norm_ruv_input(module_ci_metab_norm_input(ruv_mode = "skip")),
    object = object
  )))
  expect_true("RUV-III skipped" %in% skip_state$logs)
  expect_true(skip_state$norm_data$ruv_complete)
  expect_identical(skip_state$norm_data$ruv_corrected_obj, skip_state$norm_data$post_norm_obj)

  optimizer_stub <- function(theObject, ruv_mode, params, experiment_paths) {
    if (identical(ruv_mode, "manual") && params$manual_k < 1) {
      stop("invalid k")
    }
    module_ci_metab_norm_ruv_results(theObject, mode = ruv_mode)
  }

  # KIT-072: the optimisation and the correction are the one declared operation `metab_ruv` on the owner (the
  # removed optimisation/correction step helpers were the test double's path); the lane's kernels are bound in the
  # namespace the declaration calls.
  identity_ruv <- function(theObject, ruv_grouping_variable, ruv_number_k, ctrl) {
    theObject@args$ruvIII_C_Varying <- list(ruv_grouping_variable = ruv_grouping_variable,
      ruv_number_k = ruv_number_k, ctrl = ctrl)
    theObject
  }
  run_ruv <- function(mode, optimizer, ruv_fn = identity_ruv, manual_k = 2L, manual_percentage = 40,
                      percentage_min = 5, percentage_max = 10, adaptive = FALSE) {
    case_workflow <- module_ci_metab_norm_workflow(object)
    case_norm_data <- module_ci_metab_norm_data(object)
    state <- suppressWarnings(suppressMessages(module_ci_metab_norm_run_pipeline(
      input = module_ci_metab_norm_ruv_input(module_ci_metab_norm_input(ruv_mode = mode),
        percentage_min = percentage_min, percentage_max = percentage_max, adaptive = adaptive,
        manual_k = manual_k, manual_percentage = manual_percentage),
      object = object, workflow = case_workflow, norm_data = case_norm_data, optimizer_fn = optimizer, ruv_fn = ruv_fn
    )))
    corrected <- if (!is.null(case_norm_data$ruv_corrected_obj)) {
      module_ci_metab_norm_state(case_workflow, case_norm_data$ruv_corrected_obj)
    }
    list(state = state, workflow = case_workflow, norm_data = case_norm_data,
      best_k = if (!is.null(corrected)) corrected@args$ruvIII_C_Varying$ruv_number_k)
  }

  manual <- run_ruv("manual", optimizer_stub)
  manual_k <- manual$best_k
  expect_named(manual_k, c("LCMS_Pos", "GCMS"))
  expect_identical(unname(unlist(manual_k)), c(2L, 2L))
  expect_true(all(vapply(extractCtrlPerAssay(manual$norm_data$ruv_optimization_results), function(x) {
    is.logical(x) && !is.null(names(x)) && sum(x) == 6L
  }, logical(1))))

  auto <- run_ruv("automatic", optimizer_stub, percentage_min = 10, percentage_max = 12, adaptive = TRUE)
  expect_identical(unname(unlist(auto$best_k)), c(1L, 2L))
  expect_equal(auto$norm_data$ruv_optimization_results$LCMS_Pos$best_percentage, 13)
  expect_s3_class(auto$norm_data$ruv_optimization_results$GCMS$optimization_results, "data.frame")
  # The correction: committed under `metab_ruv_corrected` with the per-assay k of the optimisation.
  expect_true(isSavedWorkflowReference(auto$norm_data$ruv_corrected_obj)) # KIT-075: the state's reference
  expect_true(auto$norm_data$ruv_complete)
  expect_true("metab_ruv_corrected" %in% auto$workflow$state_manager$getHistory())
  expect_identical(module_ci_metab_norm_state(auto$workflow, auto$norm_data$ruv_corrected_obj)@args$ruvIII_C_Varying$ruv_number_k$LCMS_Pos, 1L)

  failed_results <- module_ci_metab_norm_ruv_results(object)
  failed_results$LCMS_Pos <- list(
    success = FALSE,
    best_k = NA_integer_,
    best_percentage = NA_real_,
    control_genes_index = NULL,
    error = "too few controls"
  )
  failed <- run_ruv("manual", function(...) failed_results, manual_percentage = 1)
  expect_true(is.na(failed$best_k$LCMS_Pos))
  expect_null(extractCtrlPerAssay(failed_results)$LCMS_Pos)
  expect_identical(failed$best_k$GCMS, 2L)

  invalid <- run_ruv("manual", optimizer_stub, ruv_fn = function(...) stop("invalid k"), manual_k = 0L)
  invalid_norm_data <- invalid$norm_data
  # The core reports the failure in the tab's log (and notification) instead of raising it.
  expect_false(invalid$state$result)
  expect_true("ERROR: invalid k" %in% invalid$state$logs)
  expect_false(invalid_norm_data$ruv_complete)
  expect_false("metab_ruv_corrected" %in% invalid$workflow$state_manager$getHistory())
})

test_that("MCI-015.4 correlation filtering matrix covers skip, pass/fail, boundary, small-n, and invalid inputs", {
  object <- module_ci_metab_norm_object(layout = "combined", positive_only = TRUE)
  norm_data <- module_ci_metab_norm_data(object)
  norm_data$post_norm_obj <- object
  norm_data$normalization_complete <- TRUE
  workflow <- module_ci_metab_norm_workflow(object)
  # KIT-105 W5: the shared frame's correlation handlers (mod_omic_norm.R) under the metab configuration
  # (runMetabNormSkipCorrelationObserverEntry and dispatchMetabNormApplyCorrelation were removed with
  # mod_metab_norm_session_helpers.R). They require a finished normalisation and notify through shiny by name, so the
  # lane's `normData` is marked normalised and the notifications are bound to quiet doubles here.
  config <- omicGuiConfig("metabolomics")
  testthat::local_mocked_bindings(showNotification = function(...) invisible(NULL), removeNotification = function(...) invisible(NULL),
    .package = "shiny")

  skip_state <- skipOmicNormCorrelation(
    workflowData = workflow,
    normData = norm_data,
    addLog = function(...) invisible(NULL),
    config = config
  )
  expect_true(skip_state)
  expect_true(norm_data$correlation_filtering_complete)
  expect_identical(norm_data$correlation_filtered_obj, object)
  expect_true("metab_norm_complete" %in% workflow$state_manager$getHistory())

  run_apply_case <- function(case) {
    # KIT-072: the declared pairs operation calls the namespace kernel, so the case's pairs are bound there.
    testthat::local_mocked_bindings(pearsonCorForSamplePairs = function(...) {
      module_ci_metab_norm_correlation_pairs(case, names(object@metabolite_data))
    }, .package = "MultiScholaR")
    case_workflow <- module_ci_metab_norm_workflow(object)
    case_norm_data <- module_ci_metab_norm_data(object)
    case_norm_data$normalization_complete <- TRUE
    case_norm_data$post_norm_obj <- object
    suppressWarnings(suppressMessages(applyOmicNormCorrelation(
      input = module_ci_metab_norm_input(correlation_threshold = 0.75),
      workflowData = case_workflow,
      normData = case_norm_data,
      addLog = function(...) invisible(NULL),
      config = config
    ))) -> result
    list(result = result, workflow = case_workflow, norm_data = case_norm_data)
  }

  pass_all <- run_apply_case("pass_all")
  expect_true(pass_all$result)
  expect_setequal(module_ci_metab_norm_state(pass_all$workflow, pass_all$norm_data$correlation_filtered_obj)@design_matrix$Run, c("S1", "S2", "S3", "S4"))
  expect_true("metab_correlation_filtered" %in% pass_all$workflow$state_manager$getHistory())

  boundary <- run_apply_case("boundary")
  expect_true(boundary$result)
  expect_setequal(module_ci_metab_norm_state(boundary$workflow, boundary$norm_data$correlation_filtered_obj)@design_matrix$Run, c("S1", "S2", "S3", "S4"))

  fail_one <- run_apply_case("fail_one")
  expect_true(fail_one$result)
  expect_identical(module_ci_metab_norm_state(fail_one$workflow, fail_one$norm_data$correlation_filtered_obj)@design_matrix$Run, c("S1", "S2"))
  expect_identical(module_ci_metab_norm_sample_cols(module_ci_metab_norm_state(fail_one$workflow, fail_one$norm_data$correlation_filtered_obj), "LCMS_Pos"), c("S1", "S2"))

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
  invalid_norm_data$normalization_complete <- TRUE
  invalid_norm_data$post_norm_obj <- object
  testthat::local_mocked_bindings(pearsonCorForSamplePairs = function(...) 1, .package = "MultiScholaR")
  invalid_state <- suppressMessages(applyOmicNormCorrelation(
    input = module_ci_metab_norm_input(correlation_threshold = 0.75),
    workflowData = invalid_workflow,
    normData = invalid_norm_data,
    addLog = function(...) invisible(NULL),
    config = config
  ))
  expect_false(invalid_state)
  expect_false(invalid_norm_data$correlation_filtering_complete)
  expect_false("metab_correlation_filtered" %in% invalid_workflow$state_manager$getHistory())
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
  # KIT-071: the session RDS pair is published under the owner's request (the journal and completion record).
  # KIT-105 W5: the tab's saver is the core's (saveWorkflowOmicsNormSessionArtifacts) under metab's export values, with
  # the summary's time fixed (saveMetabNormExportSessionRdsFiles, saveMetabNormExportMetadataFiles and
  # saveMetabNormExportSummaryFile were removed with mod_metab_norm_session_helpers.R).
  files <- saveWorkflowOmicsNormSessionArtifacts(
    binding = metabolomicsDeclarationBinding,
    sessionData = session_data,
    sourceDir = source_dir,
    publicationRequest = list(request_id = "mci-015-5", timestamp = format(fixed_time, "%Y%m%d_%H%M%S")),
    export = omicGuiConfig("metabolomics")$norm$export,
    summaryTime = fixed_time,
    logInfoFn = function(...) invisible(NULL),
    logWarnFn = function(...) invisible(NULL)
  )
  summary <- list(summaryFilepath = files$summaryFilepath,
    summaryContent = paste(readLines(files$summaryFilepath, warn = FALSE), collapse = "\n"))

  expect_true(file.exists(files$sessionFilepath))
  expect_true(file.exists(files$latestFilepath))
  expect_identical(files$publication$status, "complete")
  expect_true(file.exists(file.path(source_dir, "export_complete_mci-015-5.json")))
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

    # KIT-072: a file from another session adopts into the owner, which attaches to the module's session.
    result <- shiny::withReactiveDomain(session, runMetabDaLoadSessionObserverShell(
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
    ))

    expect_identical(result$status, "success")
    # KIT-030: the DA tab holds the state name; the object is read through the state manager.
    expect_s4_class(workflow$state_manager$getState(da_data$current_state_name), "MetaboliteAssayData")
    expect_identical(da_data$assays_available, c("LCMS_Pos", "GCMS"))
    expect_identical(da_data$formula_from_s4, "~ 0 + group + batch")
    expect_identical(capture$text$formula_string, "~ 0 + group + batch")
    expect_identical(capture$select$volcano_contrast$choices, "KO_vs_WT")
    expect_identical(capture$select$table_assay$choices, c("All", "LCMS_Pos", "GCMS"))
    expect_true(session_data$r6_current_state_name %in% workflow$state_manager$getHistory())
  }
})
