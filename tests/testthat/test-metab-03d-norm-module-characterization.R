# fidelity-coverage-compare: shared
library(testthat)
library(shiny)

source("helpers-scoped-mocked-bindings.R")

if (!methods::isClass("MetaboliteAssayData")) {
  methods::setClass(
    "MetaboliteAssayData",
    slots = c(
      metabolite_data = "list",
      design_matrix = "data.frame",
      sample_id = "character",
      group_id = "character",
      metabolite_id_column = "character",
      annotation_id_column = "character",
      args = "list"
    ),
    prototype = list(
      metabolite_data = list(),
      design_matrix = data.frame(),
      sample_id = "Run",
      group_id = "group",
      metabolite_id_column = "metabolite_id",
      annotation_id_column = "annotation_id",
      args = list()
    )
  )
}

makeMetabNormCharacterizationData <- function(label = "metab_norm_fixture") {
  methods::new(
    "MetaboliteAssayData",
    metabolite_data = list(
      Plasma = data.frame(
        metabolite_id = c(paste0(label, "_m1"), paste0(label, "_m2")),
        annotation_id = c("a1", "a2"),
        S1 = c(10, 20),
        S2 = c(30, 40),
        stringsAsFactors = FALSE
      )
    ),
    design_matrix = data.frame(
      Run = c("S1", "S2"),
      group = c("A", "B"),
      batch = c("B1", "B2"),
      stringsAsFactors = FALSE
    ),
    sample_id = "Run",
    group_id = "group",
    metabolite_id_column = "metabolite_id",
    annotation_id_column = "annotation_id",
    args = list(label = label)
  )
}

multiScholaRNamespace <- function() {
  asNamespace("MultiScholaR")
}

hasMultiScholaRBinding <- function(name) {
  exists(name, envir = multiScholaRNamespace(), inherits = FALSE)
}

getMultiScholaRBinding <- function(name) {
  get(name, envir = multiScholaRNamespace(), inherits = FALSE)
}

register_binding_teardown(
  multiScholaRNamespace(),
  c(
    "generateMetabQcPlots",
    "logTransformAssays",
    "savePlot",
    "plotPca",
    "pearsonCorForSamplePairs",
    "filterSamplesByMetaboliteCorrelationThreshold"
  ),
  .local_envir = environment()
)

skipIfMissingMultiScholaRBindings <- function(...) {
  names <- unlist(list(...), use.names = FALSE)
  missing <- names[!vapply(names, hasMultiScholaRBinding, logical(1))]
  if (length(missing) > 0) {
    skip(paste("requires extracted helper bindings:", paste(missing, collapse = ", ")))
  }
}

makeMetabNormCharacterizationHarness <- function(
  current_state = "metab_post_norm",
  current_s4 = makeMetabNormCharacterizationData()
) {
  capture <- new.env(parent = emptyenv())
  capture$saved_states <- list()

  root_dir <- tempfile("metab-norm-characterization-")
  dir.create(root_dir, recursive = TRUE)
  metabolite_qc_dir <- file.path(root_dir, "qc")
  source_dir <- file.path(root_dir, "source")
  dir.create(metabolite_qc_dir, recursive = TRUE)
  dir.create(source_dir, recursive = TRUE)

  state_manager <- new.env(parent = emptyenv())
  state_manager$current_state <- current_state
  state_manager$getState <- function(state = current_state) current_s4
  state_manager$saveState <- function(...) {
    capture$saved_states[[length(capture$saved_states) + 1L]] <<- list(...)
  }

  workflow_data <- new.env(parent = emptyenv())
  workflow_data$state_manager <- state_manager
  workflow_data$design_matrix <- data.frame(
    sample_id = c("S1", "S2"),
    group = c("A", "B"),
    batch = c("B1", "B2"),
    stringsAsFactors = FALSE
  )
  workflow_data$contrasts_tbl <- data.frame(
    friendly_names = "A vs B",
    stringsAsFactors = FALSE
  )
  workflow_data$config_list <- list(method = "TIC")
  workflow_data$metabolite_counts <- list(Plasma = list(features = 2, samples = 2))
  workflow_data$qc_params <- list(minimum = 1)
  workflow_data$tab_status <- list(
    quality_control = "pending",
    normalization = "pending",
    differential_expression = "locked"
  )

  experiment_paths <- list(
    metabolite_qc_dir = metabolite_qc_dir,
    source_dir = source_dir,
    export_dir = source_dir
  )

  list(
    workflow_data = workflow_data,
    experiment_paths = experiment_paths,
    session = shiny::MockShinySession$new(),
    current_s4 = current_s4,
    capture = capture
  )
}

launchMetabNormModule <- function(harness, selected_tab = NULL) {
  shiny::withReactiveDomain(harness$session, {
    mod_metab_norm_server(
      id = "norm",
      workflow_data = harness$workflow_data,
      experiment_paths = harness$experiment_paths,
      omic_type = "metabolomics",
      experiment_label = "Metabolomics",
      selected_tab = selected_tab
    )
  })

  harness$session$getReturned()
}

setMetabNormInputs <- function(session, ...) {
  values <- list(...)
  names(values) <- paste0("norm-", names(values))
  do.call(session$setInputs, values)
  shiny:::flushReact()
}

test_that("mod_metab_norm_server preserves selected-tab pre-QC hydration behavior", {
  harness <- makeMetabNormCharacterizationHarness(current_state = "metab_import_complete")
  capture <- new.env(parent = emptyenv())
  capture$pre_qc_calls <- 0L
  selected_tab_value <- shiny::reactiveVal("overview")

  scoped_mocked_bindings(
    generateMetabQcPlots = function(...) {
      capture$pre_qc_calls <- capture$pre_qc_calls + 1L
      invisible(NULL)
    },
    .env = asNamespace("MultiScholaR")
  )

  launchMetabNormModule(
    harness,
    selected_tab = function() selected_tab_value()
  )

  expect_identical(capture$pre_qc_calls, 0L)
  selected_tab_value("norm")
  shiny:::flushReact()

  expect_identical(capture$pre_qc_calls, 1L)
})

test_that("mod_metab_norm_server preserves normalization observer success behavior", {
  harness <- makeMetabNormCharacterizationHarness()

  scoped_mocked_bindings(
    generateMetabQcPlots = function(...) invisible(NULL),
    logTransformAssays = function(theObject, ...) theObject,
    savePlot = function(...) invisible(NULL),
    plotPca = function(...) list(Plasma = ggplot2::ggplot()),
    .env = asNamespace("MultiScholaR")
  )

  launchMetabNormModule(harness)

  setMetabNormInputs(
    harness$session,
    color_variable = "group",
    shape_variable = "batch",
    apply_itsd = FALSE,
    norm_method = "none",
    ruv_mode = "skip",
    run_normalization = 1
  )

  expect_true(any(vapply(harness$capture$saved_states, \(x) identical(x$state_name, "metab_normalized"), logical(1))))
})

test_that("mod_metab_norm_server preserves normalization observer error behavior", {
  harness <- makeMetabNormCharacterizationHarness()

  scoped_mocked_bindings(
    generateMetabQcPlots = function(...) invisible(NULL),
    logTransformAssays = function(...) stop("normalization boom"),
    plotPca = function(...) list(Plasma = ggplot2::ggplot()),
    .env = asNamespace("MultiScholaR")
  )

  launchMetabNormModule(harness)

  setMetabNormInputs(
    harness$session,
    apply_itsd = FALSE,
    norm_method = "none",
    ruv_mode = "skip",
    run_normalization = 1
  )

  expect_false(any(vapply(harness$capture$saved_states, \(x) identical(x$state_name, "metab_normalized"), logical(1))))
})

test_that("mod_metab_norm_server preserves apply-correlation public behavior", {
  harness <- makeMetabNormCharacterizationHarness()
  filtered_s4 <- makeMetabNormCharacterizationData(label = "filtered")

  scoped_mocked_bindings(
    generateMetabQcPlots = function(...) invisible(NULL),
    logTransformAssays = function(theObject, ...) theObject,
    savePlot = function(...) invisible(NULL),
    plotPca = function(...) list(Plasma = ggplot2::ggplot()),
    pearsonCorForSamplePairs = function(...) {
      list(Plasma = data.frame(pearson_correlation = 0.96, stringsAsFactors = FALSE))
    },
    filterSamplesByMetaboliteCorrelationThreshold = function(...) filtered_s4,
    .env = asNamespace("MultiScholaR")
  )

  launchMetabNormModule(harness)

  setMetabNormInputs(
    harness$session,
    apply_itsd = FALSE,
    norm_method = "none",
    ruv_mode = "skip",
    run_normalization = 1
  )

  setMetabNormInputs(
    harness$session,
    ruv_grouping_variable = "batch",
    min_pearson_correlation_threshold = 0.85,
    apply_correlation_filter = 1
  )

  expect_true(any(vapply(harness$capture$saved_states, \(x) identical(x$state_name, "metab_correlation_filtered"), logical(1))))
  expect_identical(
    harness$workflow_data$tab_status,
    list(
      quality_control = "complete",
      normalization = "complete",
      differential_expression = "locked"
    )
  )
})

test_that("mod_metab_norm_server preserves skip-correlation public behavior", {
  harness <- makeMetabNormCharacterizationHarness()
  scoped_mocked_bindings(
    generateMetabQcPlots = function(...) invisible(NULL),
    logTransformAssays = function(theObject, ...) theObject,
    savePlot = function(...) invisible(NULL),
    plotPca = function(...) list(Plasma = ggplot2::ggplot()),
    .env = asNamespace("MultiScholaR")
  )

  launchMetabNormModule(harness)

  setMetabNormInputs(
    harness$session,
    apply_itsd = FALSE,
    norm_method = "none",
    ruv_mode = "skip",
    run_normalization = 1
  )

  setMetabNormInputs(harness$session, skip_correlation_filter = 1)

  expect_true(any(vapply(harness$capture$saved_states, \(x) identical(x$state_name, "metab_norm_complete"), logical(1))))
  expect_identical(
    harness$workflow_data$tab_status,
    list(
      quality_control = "complete",
      normalization = "complete",
      differential_expression = "locked"
    )
  )
})

test_that("mod_metab_norm_server preserves export prereq and success behavior", {
  harness <- makeMetabNormCharacterizationHarness(current_state = "metab_norm_complete")
  scoped_mocked_bindings(
    generateMetabQcPlots = function(...) invisible(NULL),
    logTransformAssays = function(theObject, ...) theObject,
    savePlot = function(...) invisible(NULL),
    plotPca = function(...) list(Plasma = ggplot2::ggplot()),
    .env = asNamespace("MultiScholaR")
  )
  launchMetabNormModule(harness)

  setMetabNormInputs(
    harness$session,
    apply_itsd = FALSE,
    norm_method = "none",
    ruv_mode = "skip",
    itsd_aggregation = "median",
    log_offset = 1,
    min_pearson_correlation_threshold = 0.8,
    ruv_grouping_variable = "group",
    export_session = 1
  )

  expect_false(file.exists(file.path(harness$experiment_paths$source_dir, "metab_filtered_session_data_latest.rds")))

  setMetabNormInputs(harness$session, run_normalization = 1)
  setMetabNormInputs(harness$session, skip_correlation_filter = 1)

  setMetabNormInputs(harness$session, export_session = 2)

  expect_true(file.exists(file.path(harness$experiment_paths$source_dir, "metab_filtered_session_data_latest.rds")))
  expect_true(file.exists(file.path(harness$experiment_paths$source_dir, "metab_filtered_session_summary.txt")))
})

test_that("mod_metab_norm_server preserves reset behavior", {
  harness <- makeMetabNormCharacterizationHarness()
  scoped_mocked_bindings(
    generateMetabQcPlots = function(...) invisible(NULL),
    logTransformAssays = function(theObject, ...) theObject,
    savePlot = function(...) invisible(NULL),
    plotPca = function(...) list(Plasma = ggplot2::ggplot()),
    .env = asNamespace("MultiScholaR")
  )
  launchMetabNormModule(harness)

  setMetabNormInputs(
    harness$session,
    apply_itsd = FALSE,
    norm_method = "none",
    ruv_mode = "skip",
    run_normalization = 1
  )

  setMetabNormInputs(harness$session, reset_normalization = 1)

  expect_true(any(vapply(harness$capture$saved_states, \(x) identical(x$state_name, "metab_reset"), logical(1))))
})

test_that("mod_metab_norm_server preserves reset error behavior", {
  harness <- makeMetabNormCharacterizationHarness()
  scoped_mocked_bindings(
    generateMetabQcPlots = function(...) invisible(NULL),
    logTransformAssays = function(theObject, ...) theObject,
    savePlot = function(...) invisible(NULL),
    plotPca = function(...) list(Plasma = ggplot2::ggplot()),
    .env = asNamespace("MultiScholaR")
  )

  launchMetabNormModule(harness)

  setMetabNormInputs(
    harness$session,
    apply_itsd = FALSE,
    norm_method = "none",
    ruv_mode = "skip",
    run_normalization = 1
  )
  harness$workflow_data$state_manager$saveState <- function(...) stop("save failed")
  pre_reset_save_count <- length(harness$capture$saved_states)

  setMetabNormInputs(harness$session, reset_normalization = 1)

  expect_identical(length(harness$capture$saved_states), pre_reset_save_count)
})

test_that("metab norm shared render and plotting helpers preserve current contracts", {
  skipIfMissingMultiScholaRBindings(
    "getPlotAesthetics",
    "appendMetabNormNormalizationLog",
    "renderMetabNormNormalizationLog",
    "buildMetabNormCorrelationFilterSummary",
    "resolveMetabNormFinalQcRenderState",
    "buildMetabNormFinalQcPcaPlot",
    "renderMetabNormCorrelationFilterSummary",
    "renderMetabNormFinalQcPlot",
    "renderMetabNormAssayLabel",
    "renderMetabNormQcImageForAssay",
    "buildMetabNormLabelPlot",
    "buildMetabNormTitlePlot",
    "loadMetabNormImageAsPlot"
  )
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("png")

  get_plot_aesthetics <- getMultiScholaRBinding("getPlotAesthetics")
  append_log <- getMultiScholaRBinding("appendMetabNormNormalizationLog")
  render_log <- getMultiScholaRBinding("renderMetabNormNormalizationLog")
  build_corr_summary <- getMultiScholaRBinding("buildMetabNormCorrelationFilterSummary")
  resolve_render_state <- getMultiScholaRBinding("resolveMetabNormFinalQcRenderState")
  build_final_qc_pca <- getMultiScholaRBinding("buildMetabNormFinalQcPcaPlot")
  render_corr_summary <- getMultiScholaRBinding("renderMetabNormCorrelationFilterSummary")
  render_final_qc <- getMultiScholaRBinding("renderMetabNormFinalQcPlot")
  render_assay_label <- getMultiScholaRBinding("renderMetabNormAssayLabel")
  render_qc_image <- getMultiScholaRBinding("renderMetabNormQcImageForAssay")
  build_label_plot <- getMultiScholaRBinding("buildMetabNormLabelPlot")
  build_title_plot <- getMultiScholaRBinding("buildMetabNormTitlePlot")
  load_image_plot <- getMultiScholaRBinding("loadMetabNormImageAsPlot")

  expect_identical(
    get_plot_aesthetics(colorVariable = NULL, shapeVariable = ""),
    list(color_var = "group", shape_var = "group")
  )
  expect_identical(
    get_plot_aesthetics(colorVariable = "batch", shapeVariable = "factor1"),
    list(color_var = "batch", shape_var = "factor1")
  )

  norm_data <- new.env(parent = emptyenv())
  norm_data$normalization_log <- "[00:00:00] existing"
  expect_identical(
    append_log(normData = norm_data, message = "new entry", timestampFn = function() "07:08:09"),
    c("[00:00:00] existing", "[07:08:09] new entry")
  )
  render_text_capture <- function(expr) eval.parent(substitute(expr))
  expect_identical(
    render_log(normData = norm_data, renderTextFn = render_text_capture),
    "[00:00:00] existing\n[07:08:09] new entry"
  )

  filtered_object <- makeMetabNormCharacterizationData(label = "filtered")
  original_object <- makeMetabNormCharacterizationData(label = "original")
  corr_results <- list(
    Plasma = data.frame(pearson_correlation = c(0.8, 0.9), stringsAsFactors = FALSE)
  )
  summary_text <- build_corr_summary(
    corrResults = corr_results,
    filteredObject = filtered_object,
    originalObject = original_object
  )
  expect_match(summary_text, "=== Correlation Filtering Summary ===", fixed = TRUE)
  expect_match(summary_text, "After filtering: 2 samples", fixed = TRUE)

  render_state <- resolve_render_state(
    correlationFilteredObject = filtered_object,
    ruvCorrectedObject = original_object,
    postNormObject = NULL
  )
  expect_identical(render_state$sourceObject, filtered_object)
  expect_identical(render_state$sourceStage, "correlation_filter")
  empty_state <- resolve_render_state(
    correlationFilteredObject = NULL,
    ruvCorrectedObject = NULL,
    postNormObject = NULL
  )
  expect_true(isTRUE(empty_state$isFallback))
  expect_s3_class(empty_state$plot, "ggplot")

  wrapped_plot <- build_final_qc_pca(
    sourceObject = filtered_object,
    colorVar = "group",
    shapeVar = "batch",
    plotPcaFn = function(...) {
      list(
        ggplot2::ggplot(data.frame(x = 1, y = 1), ggplot2::aes(x, y)) + ggplot2::geom_point(),
        ggplot2::ggplot(data.frame(x = 2, y = 2), ggplot2::aes(x, y)) + ggplot2::geom_point()
      )
    },
    wrapPlotsFn = function(plots, ncol) list(plots = plots, ncol = ncol)
  )
  expect_equal(wrapped_plot$ncol, 1)
  expect_length(wrapped_plot$plots, 2)

  norm_state <- list(
    correlation_filtering_complete = TRUE,
    correlation_results = corr_results,
    correlation_filtered_obj = filtered_object,
    ruv_corrected_obj = original_object,
    post_norm_obj = NULL,
    normalization_log = norm_data$normalization_log
  )
  expect_match(
    render_corr_summary(
      normData = norm_state,
      renderTextFn = render_text_capture,
      buildSummaryFn = build_corr_summary
    ),
    "Sample pairs: 2",
    fixed = TRUE
  )

  final_plot <- render_final_qc(
    normData = norm_state,
    colorVariableFn = function() "group",
    shapeVariableFn = function() "batch",
    renderPlotFn = function(expr) eval.parent(substitute(expr)),
    reqFn = function(value) {
      if (!isTRUE(value)) stop("final qc not ready")
    },
    resolveRenderStateFn = resolve_render_state,
    getPlotAestheticsFn = get_plot_aesthetics,
    buildPcaPlotFn = function(sourceObject, colorVar, shapeVar) {
      expect_identical(sourceObject, filtered_object)
      expect_identical(colorVar, "group")
      expect_identical(shapeVar, "batch")
      "final_qc_plot"
    }
  )
  expect_identical(final_plot, "final_qc_plot")

  expect_identical(
    render_assay_label(
      assaySlot = 1,
      getAssayNamesFn = function() c("Plasma", "Urine"),
      renderTextFn = render_text_capture
    ),
    "Assay: Plasma"
  )

  qc_dir <- tempdir()
  img_path <- file.path(qc_dir, "plasma_post_norm_pca.png")
  png::writePNG(array(1, dim = c(1, 1, 4)), img_path)
  image_payload <- render_qc_image(
    assaySlot = 1,
    plotType = "pca",
    stagePrefix = "post_norm",
    normData = list(plot_refresh_trigger = 1, assay_names = c("Plasma")),
    qcDir = qc_dir,
    renderImageFn = function(expr, deleteFile) eval.parent(substitute(expr)),
    fileExistsFn = function(path) file.exists(path),
    filePathFn = function(...) file.path(...),
    sanitizeAssayNameFn = function(assayName) gsub("[^A-Za-z0-9]", "_", tolower(assayName))
  )
  expect_identical(image_payload$src, img_path)
  expect_identical(image_payload$alt, "pca - Plasma")

  label_plot <- build_label_plot("a)")
  title_plot <- build_title_plot("Pre-Normalisation")
  expect_s3_class(label_plot, "ggplot")
  expect_s3_class(title_plot, "ggplot")

  loaded_plot <- load_image_plot(img_path, logWarnFn = function(...) NULL)
  missing_plot <- load_image_plot(file.path(qc_dir, "missing-plot.png"), logWarnFn = function(...) NULL)
  expect_s3_class(loaded_plot, "ggplot")
  expect_s3_class(missing_plot, "ggplot")
})

test_that("metab norm shared export helpers preserve session and file workflows", {
  skipIfMissingMultiScholaRBindings(
    "checkMetabNormExportSessionReady",
    "resolveMetabNormExportSourceDir",
    "collectMetabNormFeatureCountsPerAssay",
    "buildMetabNormExportSessionData",
    "saveMetabNormExportSessionRdsFiles",
    "saveMetabNormExportMetadataFiles",
    "saveMetabNormExportSummaryFile",
    "runMetabNormExportSessionWorkflow"
  )

  check_export_ready <- getMultiScholaRBinding("checkMetabNormExportSessionReady")
  resolve_source_dir <- getMultiScholaRBinding("resolveMetabNormExportSourceDir")
  collect_feature_counts <- getMultiScholaRBinding("collectMetabNormFeatureCountsPerAssay")
  build_session_data <- getMultiScholaRBinding("buildMetabNormExportSessionData")
  save_session_rds <- getMultiScholaRBinding("saveMetabNormExportSessionRdsFiles")
  save_metadata <- getMultiScholaRBinding("saveMetabNormExportMetadataFiles")
  save_summary <- getMultiScholaRBinding("saveMetabNormExportSummaryFile")
  run_export_workflow <- getMultiScholaRBinding("runMetabNormExportSessionWorkflow")

  notifications <- list()
  expect_false(
    check_export_ready(
      normalizationComplete = FALSE,
      showNotificationFn = function(message, type, duration) {
        notifications[[length(notifications) + 1L]] <<- list(message = message, type = type, duration = duration)
      }
    )
  )
  expect_equal(notifications[[1]]$type, "warning")
  expect_equal(
    resolve_source_dir(
      experimentPaths = list(source_dir = "missing-source", export_dir = tempdir()),
      dirExistsFn = function(path) identical(path, tempdir()),
      dirCreateFn = function(...) TRUE
    ),
    tempdir()
  )

  current_s4 <- makeMetabNormCharacterizationData(label = "export")
  feature_counts <- collect_feature_counts(current_s4)
  expect_identical(feature_counts$Plasma, list(features = 2L, samples = 2L))

  workflow_data <- list(
    state_manager = list(
      current_state = "correlation_filtered",
      getState = function(state_name) current_s4
    ),
    contrasts_tbl = data.frame(friendly_names = "B vs A", stringsAsFactors = FALSE),
    design_matrix = current_s4@design_matrix,
    config_list = list(method = "TIC"),
    metabolite_counts = list(Plasma = list(features = 2L, samples = 2L)),
    qc_params = list(minimum = 1)
  )
  norm_data <- list(
    itsd_selections = list(Plasma = c("IS1", "IS2")),
    ruv_optimization_results = list(Plasma = list(success = TRUE, best_k = 2, best_percentage = 15.5, control_genes_index = c(TRUE, FALSE))),
    correlation_results = list(Plasma = data.frame(pearson_correlation = c(0.8, 0.9), stringsAsFactors = FALSE)),
    assay_names = c("Plasma"),
    normalization_complete = TRUE,
    ruv_complete = TRUE,
    correlation_filtering_complete = TRUE
  )
  input_values <- list(
    norm_method = "TIC",
    ruv_mode = "RUVIII-C",
    apply_itsd = TRUE,
    itsd_aggregation = "mean",
    log_offset = 1,
    min_pearson_correlation_threshold = 0.8,
    ruv_grouping_variable = "Batch"
  )

  session_data <- build_session_data(
    workflowData = workflow_data,
    normData = norm_data,
    inputValues = input_values,
    experimentLabel = "Demo Study",
    exportTimestamp = as.POSIXct("2026-04-17 12:34:56", tz = "UTC")
  )
  expect_identical(session_data$r6_current_state_name, "correlation_filtered")
  expect_identical(session_data$feature_counts$Plasma, list(features = 2L, samples = 2L))
  expect_true(isTRUE(session_data$itsd_applied))

  save_capture <- list()
  rds_files <- save_session_rds(
    sessionData = session_data,
    sourceDir = tempdir(),
    timeFn = function() as.POSIXct("2026-04-17 14:15:16", tz = "UTC"),
    formatTimeFn = function(x, fmt) format(x, fmt, tz = "UTC"),
    saveRdsFn = function(object, path) {
      save_capture[[length(save_capture) + 1L]] <<- path
    },
    logInfoFn = function(...) NULL,
    incProgressFn = function(...) NULL
  )
  expect_true(grepl("metab_filtered_session_data_20260417_141516.rds$", rds_files$sessionFilepath))
  expect_length(save_capture, 2L)

  metadata_capture <- character()
  save_metadata(
    sessionData = session_data,
    sourceDir = tempdir(),
    saveRdsFn = function(object, path) {
      metadata_capture <<- c(metadata_capture, path)
    },
    logInfoFn = function(...) NULL,
    logWarnFn = function(...) NULL
  )
  expect_true(any(grepl("metab_ruv_optimization_results.RDS$", metadata_capture)))
  expect_true(any(grepl("metab_itsd_selections.RDS$", metadata_capture)))

  summary_capture <- list()
  summary_result <- save_summary(
    sessionData = session_data,
    sourceDir = tempdir(),
    sessionFilename = "metab_filtered_session_data_demo.rds",
    writeLinesFn = function(text, path) {
      summary_capture <<- list(text = text, path = path)
    },
    timeFn = function() as.POSIXct("2026-04-17 12:34:56", tz = "UTC"),
    formatTimeFn = function(x, fmt) format(x, fmt, tz = "UTC"),
    logInfoFn = function(...) NULL
  )
  expect_true(grepl("metab_filtered_session_summary.txt$", summary_result$summaryFilepath))
  expect_match(summary_capture$text, "Metabolomics Normalized Session Data Export Summary", fixed = TRUE)
  expect_match(summary_capture$text, "Plasma: 2 features, 2 samples", fixed = TRUE)

  workflow_result <- run_export_workflow(
    workflowData = workflow_data,
    normData = norm_data,
    inputValues = input_values,
    experimentLabel = "Demo Study",
    sourceDir = tempdir(),
    withProgressFn = function(message = NULL, value = NULL, expr) force(expr),
    incProgressFn = function(...) NULL,
    buildSessionDataFn = function(...) session_data,
    saveSessionRdsFilesFn = function(...) list(sessionFilename = "session.rds", sessionFilepath = "session.rds", latestFilename = "latest.rds", latestFilepath = "latest.rds"),
    saveMetadataFilesFn = function(...) invisible(NULL),
    saveSummaryFileFn = function(...) invisible(NULL),
    logInfoFn = function(...) NULL
  )
  expect_identical(workflow_result$sessionFilename, "session.rds")
})

test_that("metab norm shared normalization helpers preserve pipeline building blocks", {
  skipIfMissingMultiScholaRBindings(
    "runMetabNormPreNormalizationQcStep",
    "runMetabNormItsdProgressApplyShell",
    "runMetabNormItsdNormalizationStep",
    "runMetabNormLog2TransformationStep",
    "runMetabNormBetweenSampleNormalizationStep",
    "runMetabNormPostNormalizationQcStep",
    "runMetabNormRuvOptimizationStep",
    "runMetabNormRuvCorrectionStep",
    "runMetabNormRuvQcStep",
    "generateMetabNormCompositeFromFiles",
    "runMetabNormCompositeQcFigureStep",
    "runMetabNormCompositeQcRefreshShell"
  )
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("png")

  run_pre_qc_step <- getMultiScholaRBinding("runMetabNormPreNormalizationQcStep")
  run_itsd_progress <- getMultiScholaRBinding("runMetabNormItsdProgressApplyShell")
  run_itsd_step <- getMultiScholaRBinding("runMetabNormItsdNormalizationStep")
  run_log2_step <- getMultiScholaRBinding("runMetabNormLog2TransformationStep")
  run_between_sample_step <- getMultiScholaRBinding("runMetabNormBetweenSampleNormalizationStep")
  run_post_qc_step <- getMultiScholaRBinding("runMetabNormPostNormalizationQcStep")
  run_ruv_optimization_step <- getMultiScholaRBinding("runMetabNormRuvOptimizationStep")
  run_ruv_correction_step <- getMultiScholaRBinding("runMetabNormRuvCorrectionStep")
  run_ruv_qc_step <- getMultiScholaRBinding("runMetabNormRuvQcStep")
  generate_composite_from_files <- getMultiScholaRBinding("generateMetabNormCompositeFromFiles")
  run_composite_step <- getMultiScholaRBinding("runMetabNormCompositeQcFigureStep")
  run_composite_refresh <- getMultiScholaRBinding("runMetabNormCompositeQcRefreshShell")

  current_s4 <- makeMetabNormCharacterizationData(label = "pipeline")
  workflow_data <- list(
    state_manager = list(
      saveState = function(...) invisible(NULL)
    ),
    config_list = list(method = "TIC")
  )
  norm_data <- new.env(parent = emptyenv())
  norm_data$plot_refresh_trigger <- 0

  pre_state <- run_pre_qc_step(
    currentS4 = current_s4,
    totalSteps = 5,
    experimentPaths = list(metabolite_qc_dir = tempdir()),
    groupingVariable = "group",
    shapeVariable = "batch",
    normData = norm_data,
    addLogFn = function(...) NULL,
    incProgressFn = function(...) NULL,
    generateMetabQcPlotsFn = function(...) invisible(NULL)
  )
  expect_identical(pre_state$currentS4, current_s4)
  expect_identical(norm_data$post_filter_obj, current_s4)

  itsd_state <- run_itsd_step(
    currentS4 = current_s4,
    itsdAggregation = "mean",
    itsdFeatureIds = c("IS1"),
    workflowData = workflow_data,
    normData = norm_data,
    addLogFn = function(...) NULL,
    normaliseUntransformedDataFn = function(theObject, method, itsd_aggregation, itsd_feature_ids) theObject
  )
  expect_identical(itsd_state$currentS4, current_s4)

  progress_state <- run_itsd_progress(
    currentS4 = current_s4,
    totalSteps = 5,
    applyItsd = TRUE,
    itsdAggregation = "mean",
    itsdSelections = list(Plasma = c("IS1")),
    workflowData = workflow_data,
    normData = norm_data,
    addLogFn = function(...) NULL,
    incProgressFn = function(...) NULL,
    resolveManualFeatureIdsFn = function(...) c("IS1"),
    runItsdStepFn = function(...) list(currentS4 = current_s4)
  )
  expect_true(isTRUE(progress_state$applied))
  expect_identical(progress_state$currentS4, current_s4)

  log2_state <- run_log2_step(
    currentS4 = current_s4,
    logOffset = 1,
    workflowData = workflow_data,
    normData = norm_data,
    addLogFn = function(...) NULL,
    logTransformAssaysFn = function(theObject, offset) theObject
  )
  expect_identical(log2_state$currentS4, current_s4)

  between_state <- run_between_sample_step(
    currentS4 = current_s4,
    normMethod = "none",
    workflowData = workflow_data,
    normData = norm_data,
    addLogFn = function(...) NULL,
    normaliseBetweenSamplesFn = function(theObject, normalisation_method) theObject
  )
  expect_identical(between_state$currentS4, current_s4)
  expect_true(isTRUE(norm_data$normalization_complete))

  post_qc_state <- run_post_qc_step(
    currentS4 = current_s4,
    experimentPaths = list(metabolite_qc_dir = tempdir()),
    groupingVariable = "group",
    shapeVariable = "batch",
    addLogFn = function(...) NULL,
    generateMetabQcPlotsFn = function(...) invisible(NULL)
  )
  expect_identical(post_qc_state$currentS4, current_s4)

  norm_data$ruv_optimization_results <- NULL
  optimization_state <- run_ruv_optimization_step(
    currentS4 = current_s4,
    ruvMode = "automatic",
    autoPercentageMin = 1,
    autoPercentageMax = 3,
    ruvGroupingVariable = "batch",
    separationMetric = "max_difference",
    kPenaltyWeight = 0.5,
    adaptiveKPenalty = TRUE,
    manualK = 2,
    manualPercentage = 10,
    experimentPaths = list(source_dir = tempdir()),
    normData = norm_data,
    addLogFn = function(...) NULL,
    runPerAssayRuvOptimizationFn = function(...) list(Plasma = list(success = TRUE, best_k = 2, best_percentage = 15.5, control_genes_index = c(TRUE, FALSE))),
    extractBestKPerAssayFn = function(results) list(Plasma = 2),
    extractCtrlPerAssayFn = function(results) list(Plasma = c(TRUE, FALSE))
  )
  expect_identical(optimization_state$bestKPerAssay$Plasma, 2)

  correction_state <- run_ruv_correction_step(
    currentS4 = current_s4,
    ruvGroupingVariable = "batch",
    bestKPerAssay = list(Plasma = 2),
    ctrlPerAssay = list(Plasma = c(TRUE, FALSE)),
    workflowData = workflow_data,
    normData = norm_data,
    addLogFn = function(...) NULL,
    ruvIII_C_VaryingFn = function(theObject, ruv_grouping_variable, ruv_number_k, ctrl) theObject
  )
  expect_identical(correction_state$currentS4, current_s4)
  expect_true(isTRUE(norm_data$ruv_complete))

  ruv_qc_state <- run_ruv_qc_step(
    currentS4 = current_s4,
    totalSteps = 5,
    experimentPaths = list(metabolite_qc_dir = tempdir()),
    groupingVariable = "group",
    shapeVariable = "batch",
    addLogFn = function(...) NULL,
    incProgressFn = function(...) NULL,
    generateMetabQcPlotsFn = function(...) invisible(NULL)
  )
  expect_identical(ruv_qc_state$currentS4, current_s4)

  img_path <- tempfile(fileext = ".png")
  png::writePNG(array(1, dim = c(1, 1, 4)), img_path)
  composite <- generate_composite_from_files(
    plot_files = c(img_path),
    ncol = 1,
    row_labels = list(Plasma_pca = "a)"),
    column_labels = "Pre"
  )
  expect_type(composite, "list")
  expect_true(all(c("plot", "width", "height") %in% names(composite)))

  composite_state <- run_composite_step(
    experimentPaths = list(metabolite_qc_dir = tempdir()),
    assayNames = c("Plasma"),
    ruvMode = "skip",
    omicType = "metabolomics",
    addLogFn = function(...) NULL,
    generateCompositeFromFilesFn = function(...) list(plot = "plot", width = 10, height = 6),
    savePlotFn = function(...) invisible(NULL),
    dirExistsFn = function(path) TRUE,
    logWarnFn = function(...) NULL
  )
  expect_true(isTRUE(composite_state$compositeSaved))

  refresh_state <- run_composite_refresh(
    currentS4 = current_s4,
    experimentPaths = list(metabolite_qc_dir = tempdir()),
    assayNames = c("Plasma"),
    ruvMode = "skip",
    omicType = "metabolomics",
    normData = norm_data,
    addLogFn = function(...) NULL,
    generateCompositeFromFilesFn = function(...) list(plot = "plot", width = 10, height = 6),
    savePlotFn = function(...) invisible(NULL),
    logWarnFn = function(...) NULL,
    runCompositeQcFigureStepFn = function(...) list(compositeSaved = TRUE)
  )
  expect_identical(refresh_state$currentS4, current_s4)
  expect_identical(refresh_state$plotRefreshTrigger, 1)
})
