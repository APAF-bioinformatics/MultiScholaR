library(testthat)

module_ci_prot_norm_da_load_module <- function(da_data, workflow_data, experiment_paths) {
  force(da_data)
  force(workflow_data)
  force(experiment_paths)

  function(id) {
    shiny::moduleServer(id, function(input, output, session) {
      da_server_load_session_handler(
        input = input,
        output = output,
        session = session,
        da_data = da_data,
        workflow_data = workflow_data,
        experiment_paths = experiment_paths
      )
      da_data
    })
  }
}

test_that("MCI-008.1 normalization method matrix preserves numeric sanity and log contracts", {
  s4_methods <- module_ci_prot_norm_s4_methods()
  offsets <- c(none = 0, cyclicloess = 10, quantile = 20, scale = 30)

  for (method_name in names(offsets)) {
    object <- module_ci_prot_norm_object(args = module_ci_prot_norm_args(normalisation_method = method_name))
    normalized <- suppressMessages(s4_methods$normalise(object, normalisation_method = method_name))
    mat <- module_ci_prot_norm_matrix(normalized)

    expect_s4_class(normalized, "ProteinQuantitativeData")
    expect_identical(normalized@args$module_ci_param_audit$normalisation_method, method_name)
    expect_equal(mat["P_STABLE", "A1"], 10 + offsets[[method_name]])
    expect_equal(mat["P_ZERO", "A1"], 0 + offsets[[method_name]])
    expect_equal(mat["P_NEG", "A1"], -1 + offsets[[method_name]])
    expect_true(is.na(mat["P_INF", "A1"]))
    expect_true(is.na(mat["P_NAN", "A1"]))
    module_ci_prot_norm_assert_sanity(normalized, allow_na = TRUE)
  }

  already_logged <- suppressMessages(s4_methods$normalise(
    module_ci_prot_norm_object(args = module_ci_prot_norm_args(normalisation_method = "none")),
    normalisation_method = "none"
  ))
  expect_equal(module_ci_prot_norm_matrix(already_logged)["P_SHIFT", "B2"], 130)

  raw_matrix <- matrix(
    c(1, 0, 2, NA, -1, Inf),
    nrow = 2,
    byrow = TRUE,
    dimnames = list(c("P1", "P2"), c("S1", "S2", "S3"))
  )
  logged <- suppressWarnings(log2Transformation(raw_matrix))
  expect_identical(dimnames(logged), dimnames(raw_matrix))
  expect_equal(logged["P1", "S1"], log2(1 + 0.01))
  expect_true(is.infinite(logged["P1", "S2"]) && logged["P1", "S2"] < 0)
  expect_true(is.na(logged["P2", "S1"]))
  expect_true(is.infinite(logged["P2", "S3"]) && logged["P2", "S3"] > 0)
})

test_that("MCI-008.2 RUV skip, manual, automatic, controls, k, and optimizer boundaries are covered", {
  s4_methods <- module_ci_prot_norm_s4_methods()
  ctrl <- module_ci_prot_norm_control_index()
  object <- module_ci_prot_norm_object(args = module_ci_prot_norm_args(ctrl = ctrl, ruv_number_k = 2L))

  cancor <- suppressMessages(s4_methods$cancor(object))
  expect_identical(names(cancor$ctl), names(ctrl))
  expect_equal(cancor$Y["A1", "P_NA"], -99)
  expect_identical(cancor$X, object@design_matrix$group)

  expect_error(
    suppressMessages(s4_methods$cancor(module_ci_prot_norm_object(args = module_ci_prot_norm_args(
      ctrl = ctrl,
      ruv_number_k = 2L,
      ruv_grouping_variable = "missing_group"
    )))),
    "is not a column in the design matrix",
    fixed = TRUE
  )
  expect_error(
    suppressMessages(s4_methods$cancor(module_ci_prot_norm_object(args = module_ci_prot_norm_args(
      ctrl = ctrl,
      ruv_number_k = 0L
    )), num_components_to_impute = 0L)),
    "value is invalid",
    fixed = TRUE
  )
  expect_error(
    suppressMessages(s4_methods$cancor(module_ci_prot_norm_object(args = module_ci_prot_norm_args(
      ctrl = c(P1 = TRUE, P2 = TRUE, P3 = TRUE, P4 = TRUE)
    )))),
    "less than 5",
    fixed = TRUE
  )

  finite_object <- module_ci_prot_norm_object(
    protein_quant_table = module_ci_prot_norm_table(module_ci_prot_norm_values()[grep("^P_(STABLE|SHIFT|CONST|CTRL_)", rownames(module_ci_prot_norm_values())), ]),
    args = module_ci_prot_norm_args(rownames(module_ci_prot_norm_values())[grep("^P_(STABLE|SHIFT|CONST|CTRL_)", rownames(module_ci_prot_norm_values()))])
  )
  corrected <- suppressMessages(s4_methods$ruviii(finite_object, ruv_number_k = 2L))
  expect_s4_class(corrected, "ProteinQuantitativeData")
  expect_equal(module_ci_prot_norm_matrix(corrected)["P_STABLE", "A1"], 10.2)
  module_ci_prot_norm_assert_sanity(corrected, allow_na = FALSE)
  expect_error(
    suppressMessages(s4_methods$ruviii(finite_object, ruv_number_k = 0L)),
    "ruv_number_k = 0 value is invalid",
    fixed = TRUE
  )

  workflow <- module_ci_prot_norm_workflow(object)
  norm_data <- module_ci_prot_norm_data(object)
  skip_result <- applyProtNormSkippedRuvState(
    normalizedS4 = object,
    normMethod = "none",
    normData = norm_data,
    workflowData = workflow,
    sourceDir = NULL,
    saveRdsFn = function(...) invisible(NULL),
    messageFn = function(...) invisible(NULL)
  )
  expect_true(skip_result$ruv_skipped)
  expect_true(norm_data$ruv_complete)
  expect_identical(workflow$state_manager$current_state, "normalized")
  expect_false(workflow$state_manager$states$normalized@args$globalParameters$use_limpa)

  captured <- list()
  manual <- resolveProtNormRuvParameters(
    normalizedS4 = object,
    input = module_ci_prot_norm_input(ruv_mode = "manual", ruv_percentage = 40, ruv_k = 3L),
    normData = norm_data,
    workflowData = workflow,
    sourceDir = NULL,
    getRuvGroupingVariableFn = function() "group",
    getNegCtrlProtAnovaFn = function(...) ctrl,
    ruvCancorFn = function(...) module_ci_prot_norm_cancor_plot(),
    updateAuditTrailFn = function(ruvK, controlGenesIndex, percentageAsNegCtrl, modeLabel) {
      captured$manual_audit <<- list(ruvK = ruvK, pct = percentageAsNegCtrl, mode = modeLabel, controls = controlGenesIndex)
    },
    messageFn = function(...) invisible(NULL)
  )
  expect_identical(manual$ruvK, 3L)
  expect_identical(manual$percentageAsNegCtrl, 40)
  expect_identical(names(manual$controlGenesIndex), names(ctrl))
  expect_identical(captured$manual_audit$mode, "manual")
  expect_identical(norm_data$ruv_optimization_result$separation_metric_used, "manual")

  bad_manual_input <- module_ci_prot_norm_input(ruv_mode = "manual", ruv_k = 0L)
  expect_error(
    resolveProtNormRuvParameters(
      normalizedS4 = object,
      input = bad_manual_input,
      normData = norm_data,
      workflowData = workflow,
      getRuvGroupingVariableFn = function() "group",
      getNegCtrlProtAnovaFn = function(...) ctrl,
      messageFn = function(...) invisible(NULL)
    ),
    "ruv_k = 0 value is invalid",
    fixed = TRUE
  )

  auto_result <- list(
    best_percentage = 11,
    best_k = 2L,
    best_control_genes_index = ctrl,
    best_cancor_plot = module_ci_prot_norm_cancor_plot(),
    optimization_results = data.frame(
      percentage_requested = c(10, 11, 12),
      best_k = c(3L, 2L, 2L),
      composite_score = c(0.2, 0.5, 0.4),
      status = "ok"
    ),
    best_separation_score = 0.5,
    best_composite_score = 0.5,
    separation_metric_used = "max_difference",
    k_penalty_weight = 0.5,
    adaptive_k_penalty_used = TRUE
  )
  automatic <- resolveProtNormRuvParameters(
    normalizedS4 = object,
    input = module_ci_prot_norm_input(ruv_mode = "automatic"),
    normData = norm_data,
    workflowData = workflow,
    sourceDir = NULL,
    getRuvGroupingVariableFn = function() "group",
    withProgressFn = function(message, detail = NULL, value = 0, expr) force(expr),
    findBestNegCtrlPercentageFn = function(...) auto_result,
    updateAuditTrailFn = function(ruvK, controlGenesIndex, percentageAsNegCtrl, modeLabel) {
      captured$auto_audit <<- list(ruvK = ruvK, pct = percentageAsNegCtrl, mode = modeLabel, controls = controlGenesIndex)
    },
    messageFn = function(...) invisible(NULL)
  )
  expect_identical(automatic$ruvK, 2L)
  expect_identical(automatic$percentageAsNegCtrl, 11)
  expect_identical(workflow$ruv_optimization_result$best_percentage, 11)
  expect_identical(captured$auto_audit$mode, "automatic")

  bad_input <- module_ci_prot_norm_input(ruv_mode = "automatic")
  bad_input$auto_percentage_min <- 20
  bad_input$auto_percentage_max <- 10
  expect_error(
    resolveProtNormRuvParameters(
      normalizedS4 = object,
      input = bad_input,
      normData = norm_data,
      workflowData = workflow,
      getRuvGroupingVariableFn = function() "group",
      messageFn = function(...) invisible(NULL)
    ),
    "Minimum percentage must be less than maximum percentage",
    fixed = TRUE
  )

  expect_identical(findBestKElbow(module_ci_prot_norm_cancor_plot(c(0.01, 0.02))), 1L)
  expect_identical(findBestKElbow(module_ci_prot_norm_cancor_plot(c(0.10, 0.25, 0.24))), 2L)
  expect_true(is.na(findBestKElbow(list(data = data.frame(K = 1, cc = 1)))))
  expect_identical(calculateCompositeScore(1, 1, 0.5, 1), 1)
  expect_lt(calculateCompositeScore(1, 2, 0.5, 1), 1)
})

test_that("MCI-008.3 correlation filtering covers pass, fail, boundary, small-n, and missing-sample cases", {
  object <- module_ci_prot_norm_object()
  filter_case <- function(case) {
    suppressMessages(filterSamplesByProteinCorrelationThresholdHelper(
      pearson_correlation_per_pair = module_ci_prot_norm_correlation_pairs(case),
      protein_intensity_table = object@protein_quant_table,
      min_pearson_correlation_threshold = 0.75,
      filename_column_x = Run.x,
      filename_column_y = Run.y,
      protein_id_column = Protein.Ids,
      correlation_column = pearson_correlation
    ))
  }
  sample_cols <- function(tbl) setdiff(colnames(tbl), "Protein.Ids")

  expect_setequal(sample_cols(filter_case("pass_all")), c("A1", "A2", "B1", "B2"))
  expect_setequal(sample_cols(filter_case("boundary")), c("A1", "A2", "B1", "B2"))
  expect_setequal(sample_cols(filter_case("fail_one")), c("A1", "A2"))
  expect_identical(sample_cols(filter_case("fail_many")), character(0))
  expect_false("Ghost" %in% colnames(filter_case("missing_sample")))
  expect_setequal(sample_cols(filter_case("missing_sample")), c("A1", "A2", "B1", "B2"))

  one_sample_values <- module_ci_prot_norm_values(samples = c("A1", "A2", "B1", "B2"))[, "A1", drop = FALSE]
  small_object <- module_ci_prot_norm_object(
    protein_quant_table = module_ci_prot_norm_table(one_sample_values),
    design = module_ci_prot_norm_design(samples = "A1", groups = "A", replicates = "R1", batches = "Batch1"),
    args = module_ci_prot_norm_args(rownames(one_sample_values))
  )
  small_result <- suppressMessages(filterSamplesByProteinCorrelationThresholdHelper(
    pearson_correlation_per_pair = module_ci_prot_norm_empty_pairs(),
    protein_intensity_table = small_object@protein_quant_table,
    min_pearson_correlation_threshold = 0.75,
    filename_column_x = Run.x,
    filename_column_y = Run.y,
    protein_id_column = Protein.Ids,
    correlation_column = pearson_correlation
  ))
  expect_identical(sample_cols(small_result), "A1")

  norm_data <- module_ci_prot_norm_data(object)
  filtered <- suppressWarnings(suppressMessages(runProtNormCorrelationFilterStep(
    ruvS4 = object,
    correlationVec = module_ci_prot_norm_correlation_pairs("fail_one"),
    correlationThreshold = 0.75,
    normData = norm_data,
    gcFn = function(...) invisible(NULL),
    messageFn = function(...) invisible(NULL)
  )))
  expect_s4_class(filtered, "ProteinQuantitativeData")
  expect_identical(norm_data$correlation_filtered_obj, filtered)
  expect_setequal(setdiff(colnames(filtered@protein_quant_table), "Protein.Ids"), c("A1", "A2"))

  vector_data <- module_ci_prot_norm_data(object)
  vector_result <- runProtNormCorrelationVectorStep(
    ruvS4 = object,
    correlationThreshold = 0.75,
    normData = vector_data,
    getRuvGroupingVariableFn = function() "replicates",
    pearsonCorForSamplePairsFn = function(object, tech_rep_remove_regex, correlation_group) {
      expect_identical(tech_rep_remove_regex, "pool")
      expect_identical(correlation_group, "replicates")
      module_ci_prot_norm_correlation_pairs("pass_all")
    },
    messageFn = function(...) invisible(NULL)
  )
  expect_identical(vector_data$correlation_threshold, 0.75)
  expect_equal(nrow(vector_result), 2L)
})

test_that("MCI-008.4 QC plot generation covers pre, post, RUV, composite, and degenerate fixtures", {
  object <- module_ci_prot_norm_object()
  qc_dir <- tempfile("mci008-qc-")
  dir.create(qc_dir, recursive = TRUE)
  aesthetics <- list(color_var = "group", shape_var = "batch")
  plot_stub <- function(...) module_ci_prot_norm_ggplot_stub("qc")

  pre_paths <- generateProtNormPreNormalizationQcArtifacts(
    stateManager = module_ci_prot_norm_state_manager(list(protein_replicate_filtered = object), "protein_replicate_filtered"),
    qcDir = qc_dir,
    aesthetics = aesthetics,
    reqFn = function(x) x,
    gcFn = function(...) invisible(NULL),
    plotPcaFn = plot_stub,
    plotRleFn = plot_stub,
    buildDensityFn = plot_stub,
    buildCorrelationFn = plot_stub,
    saveArtifactFn = module_ci_prot_norm_artifact_saver,
    messageFn = function(...) invisible(NULL)
  )
  expect_setequal(names(pre_paths$post_filtering), c("pca", "rle", "density", "correlation"))
  expect_true(all(file.exists(unlist(pre_paths$post_filtering))))

  post_paths <- generateProtNormPostNormalizationQcArtifacts(
    normalizedS4 = object,
    qcDir = qc_dir,
    aesthetics = aesthetics,
    qcPlotPaths = pre_paths,
    gcFn = function(...) invisible(NULL),
    plotPcaFn = plot_stub,
    plotRleFn = plot_stub,
    buildDensityFn = plot_stub,
    buildCorrelationFn = plot_stub,
    saveArtifactFn = module_ci_prot_norm_artifact_saver,
    messageFn = function(...) invisible(NULL)
  )
  expect_setequal(names(post_paths$post_normalization), c("pca", "rle", "density", "correlation"))
  expect_true(all(file.exists(unlist(post_paths$post_normalization))))

  ruv_paths <- generateProtNormRuvCorrectedQcArtifacts(
    ruvCorrectedS4 = object,
    qcDir = qc_dir,
    aesthetics = aesthetics,
    qcPlotPaths = post_paths,
    gcFn = function(...) invisible(NULL),
    plotPcaFn = plot_stub,
    plotRleFn = plot_stub,
    buildDensityFn = plot_stub,
    buildCorrelationFn = plot_stub,
    saveArtifactFn = module_ci_prot_norm_artifact_saver,
    messageFn = function(...) invisible(NULL)
  )
  expect_setequal(names(ruv_paths$ruv_corrected), c("pca", "rle", "density", "correlation"))
  expect_true(all(file.exists(unlist(ruv_paths$ruv_corrected))))

  skip_inputs <- resolveProtNormCompositeFigureInputs("skip", qc_dir)
  ruv_inputs <- resolveProtNormCompositeFigureInputs("manual", qc_dir)
  expect_identical(skip_inputs$ncol, 2)
  expect_identical(ruv_inputs$ncol, 3)
  expect_true(any(grepl("ruv_corrected_cancor.png", ruv_inputs$plotFiles, fixed = TRUE)))

  norm_data <- module_ci_prot_norm_data(object, ruv_result = list(best_cancor_plot = module_ci_prot_norm_cancor_plot()))
  runProtNormStep6RuvQc(
    ruvMode = "manual",
    step6Object = object,
    normData = norm_data,
    qcDir = qc_dir,
    generateRuvCorrectedQcFn = function(ruvCorrectedS4) {
      norm_data$qc_plot_paths <- ruv_paths
    },
    ggsaveFn = function(filename, plot, width, height, dpi) {
      writeLines("cancor", filename)
      invisible(filename)
    },
    messageFn = function(...) invisible(NULL)
  )
  expect_true(file.exists(file.path(qc_dir, "ruv_corrected_cancor.png")))
  expect_identical(norm_data$qc_plot_paths$ruv_corrected$cancor, file.path(qc_dir, "ruv_corrected_cancor.png"))

  constant_values <- module_ci_prot_norm_values()
  constant_values[,] <- 42
  constant_object <- module_ci_prot_norm_object(module_ci_prot_norm_table(constant_values))
  expect_s3_class(suppressMessages(plotPca(constant_object, grouping_variable = "group", label_column = NULL, title = "constant")), "ggplot")
  expect_s3_class(suppressMessages(plotRle(constant_object, grouping_variable = "group")), "ggplot")
  expect_s3_class(suppressMessages(plotDensity(constant_object, grouping_variable = "group", title = "constant")), "ggplot")
})

test_that("MCI-008.5 filtered-session export preserves R6, params, filenames, counts, and DA-ready state", {
  object <- module_ci_prot_norm_object(args = module_ci_prot_norm_args(
    workflow_type = "DIA",
    report_template = "DIANN_limpa_report.rmd",
    use_limpa = TRUE
  ))
  workflow <- module_ci_prot_norm_workflow(
    current_object = object,
    workflow_type = "DIA",
    report_template = "DIANN_limpa_report.rmd",
    use_limpa = TRUE
  )
  norm_data <- module_ci_prot_norm_data(object)

  session_data <- collectProtNormExportSessionData(
    workflowData = workflow,
    normData = norm_data,
    input = module_ci_prot_norm_input(norm_method = "quantile", ruv_mode = "automatic"),
    timeFn = function() as.POSIXct("2026-05-05 12:34:56", tz = "UTC"),
    messageFn = function(...) invisible(NULL)
  )

  expect_identical(session_data$r6_current_state_name, "correlation_filtered")
  expect_setequal(names(session_data$r6_complete_states), c("protein_replicate_filtered", "normalized", "ruv_corrected", "correlation_filtered"))
  expect_identical(session_data$r6_state_history, c("protein_replicate_filtered", "normalized", "ruv_corrected", "correlation_filtered"))
  expect_identical(session_data$workflow_type, "DIA")
  expect_identical(session_data$report_template, "DIANN_limpa_report.rmd")
  expect_true(session_data$limpa_applied)
  expect_identical(session_data$normalization_method, "quantile")
  expect_identical(session_data$ruv_mode, "automatic")
  expect_true(session_data$ruv_applied)
  expect_identical(session_data$ruv_k, 2L)
  expect_identical(session_data$final_protein_count, nrow(object@protein_quant_table))
  expect_identical(session_data$final_sample_count, 4L)

  source_dir <- tempfile("mci008-export-")
  dir.create(source_dir, recursive = TRUE)
  artifacts <- saveProtNormExportArtifacts(
    sessionData = session_data,
    sourceDir = source_dir,
    timeFn = function() as.POSIXct("2026-05-05 12:34:56", tz = "UTC"),
    saveRdsFn = saveRDS,
    writeLinesFn = writeLines,
    messageFn = function(...) invisible(NULL)
  )
  expect_identical(artifacts$sessionFilename, "filtered_session_data_20260505_123456.rds")
  expect_identical(basename(artifacts$latestFilepath), "filtered_session_data_latest.rds")
  expect_true(file.exists(artifacts$sessionFilepath))
  expect_true(file.exists(artifacts$latestFilepath))
  expect_true(file.exists(artifacts$summaryFilepath))
  reloaded <- readRDS(artifacts$latestFilepath)
  expect_identical(reloaded$normalization_method, "quantile")
  expect_identical(reloaded$report_template, "DIANN_limpa_report.rmd")
  summary_text <- readLines(artifacts$summaryFilepath, warn = FALSE)
  expect_true(any(grepl("Workflow Type: DIA", summary_text, fixed = TRUE)))
  expect_true(any(grepl("Report Template: DIANN_limpa_report.rmd", summary_text, fixed = TRUE)))

  output <- new.env(parent = emptyenv())
  metrics <- completeProtNormCorrelationWorkflow(
    finalS4ForDe = object,
    workflowData = workflow,
    output = output,
    normData = norm_data,
    correlationThreshold = 0.75,
    skipped = FALSE,
    successNotification = "ok",
    completionMessage = "complete",
    messagePrefix = "*** MCI-008",
    showNotificationFn = function(...) invisible(NULL),
    renderTextFn = function(x) x,
    messageFn = function(...) invisible(NULL)
  )
  expect_identical(metrics$finalProteinCount, nrow(object@protein_quant_table))
  expect_identical(workflow$state_manager$current_state, "correlation_filtered")
  expect_true(norm_data$correlation_filtering_complete)
  expect_match(output$correlation_filter_summary, "Threshold: 0.75", fixed = TRUE)
  expect_identical(workflow$tab_status$differential_expression, "pending")
})

test_that("MCI-008.6 DA reload smoke consumes every exported normalization session variant", {
  skip_if_not_installed("shiny")

  global_names <- c("contrasts_tbl", "config_list", "uniprot_dat_cln")
  old_globals <- lapply(global_names, function(name) {
    if (exists(name, envir = .GlobalEnv, inherits = FALSE)) {
      list(exists = TRUE, value = get(name, envir = .GlobalEnv, inherits = FALSE))
    } else {
      list(exists = FALSE, value = NULL)
    }
  })
  names(old_globals) <- global_names
  withr::defer({
    for (name in global_names) {
      if (exists(name, envir = .GlobalEnv, inherits = FALSE)) {
        rm(list = name, envir = .GlobalEnv)
      }
      if (isTRUE(old_globals[[name]]$exists)) {
        assign(name, old_globals[[name]]$value, envir = .GlobalEnv)
      }
    }
  })

  notifications <- list()
  text_updates <- list()
  testthat::local_mocked_bindings(
    withProgress = function(message, value, expr) force(expr),
    incProgress = function(amount, detail = NULL) invisible(NULL),
    showNotification = function(message, type = NULL, duration = NULL, ...) {
      notifications[[length(notifications) + 1L]] <<- list(message = message, type = type, duration = duration)
      invisible(NULL)
    },
    updateTextAreaInput = function(session, inputId, value = NULL, ...) {
      text_updates[[inputId]] <<- value
      invisible(NULL)
    },
    .package = "shiny"
  )
  testthat::local_mocked_bindings(
    log_info = function(...) invisible(NULL),
    log_warn = function(...) invisible(NULL),
    log_error = function(...) invisible(NULL),
    .package = "logger"
  )

  variants <- list(
    skip = module_ci_prot_norm_export_variant(ruv_mode = "skip", workflow_type = "DIA", norm_method = "none"),
    manual = module_ci_prot_norm_export_variant(ruv_mode = "manual", workflow_type = "TMT", norm_method = "scale"),
    automatic_limpa = module_ci_prot_norm_export_variant(
      ruv_mode = "automatic",
      workflow_type = "DIA",
      report_template = "DIANN_limpa_report.rmd",
      use_limpa = TRUE,
      norm_method = "quantile"
    )
  )

  for (variant_name in names(variants)) {
    source_dir <- tempfile(paste0("mci008-da-load-", variant_name, "-"))
    dir.create(source_dir, recursive = TRUE)
    saveRDS(variants[[variant_name]], file.path(source_dir, "filtered_session_data_latest.rds"))
    saveRDS(
      data.frame(Entry = c("P_STABLE", "P_SHIFT"), gene_names = c("STABLE", "SHIFT"), stringsAsFactors = FALSE),
      file.path(source_dir, "uniprot_dat_cln.RDS")
    )

    da_data <- shiny::reactiveValues(
      current_s4_object = NULL,
      contrasts_available = NULL,
      formula_from_s4 = NULL
    )
    state_manager <- new.env(parent = emptyenv())
    state_manager$states <- list()
    state_manager$state_history <- character()
    state_manager$current_state <- NULL
    workflow_data <- shiny::reactiveValues(
      state_manager = state_manager,
      tab_status = list(normalization = "pending", differential_expression = "disabled"),
      state_update_trigger = NULL
    )

    shiny::testServer(
      module_ci_prot_norm_da_load_module(
        da_data = da_data,
        workflow_data = workflow_data,
        experiment_paths = list(source_dir = source_dir)
      ),
      {
        session$setInputs(load_filtered_session = 1)
        session$flushReact()

        expect_s4_class(da_data$current_s4_object, "ProteinQuantitativeData")
        expect_identical(workflow_data$state_manager$current_state, "correlation_filtered")
        expect_setequal(names(workflow_data$state_manager$states), names(variants[[variant_name]]$r6_complete_states))
        expect_identical(workflow_data$tab_status$normalization, "complete")
        expect_identical(workflow_data$tab_status$differential_expression, "pending")
        expect_s3_class(workflow_data$state_update_trigger, "POSIXct")
        expect_identical(workflow_data$config_list$globalParameters$workflow_type, variants[[variant_name]]$workflow_type)
        expect_identical(workflow_data$ruv_optimization_result, variants[[variant_name]]$ruv_optimization_result)
        expect_equal(nrow(workflow_data$uniprot_dat_cln), 2L)
      }
    )
  }

  expect_true(any(vapply(notifications, function(notification) identical(notification$type, "message"), logical(1))))
})
