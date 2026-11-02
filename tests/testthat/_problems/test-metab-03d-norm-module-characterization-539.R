# Extracted from test-metab-03d-norm-module-characterization.R:539

# prequel ----------------------------------------------------------------------
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
  current_s4 = makeMetabNormCharacterizationData(),
  env = parent.frame()
) {
  capture <- new.env(parent = emptyenv())
  capture$saved_states <- list()

  root_dir <- tempfile("metab-norm-characterization-")
  dir.create(root_dir, recursive = TRUE)
  metabolite_qc_dir <- file.path(root_dir, "qc")
  source_dir <- file.path(root_dir, "source")
  dir.create(metabolite_qc_dir, recursive = TRUE)
  dir.create(source_dir, recursive = TRUE)

  workflow_data <- new.env(parent = emptyenv())
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
  # KIT-072: the bus holds a real metabolomics owner (the module's only state manager) with the data committed under
  # `current_state`; the declared operations commit on it.
  omicsTestOwner(current_s4, current_state, workflow = workflow_data, experiment_paths = list(source_dir = source_dir),
    omic_type = "metabolomics", env = env)

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

# test -------------------------------------------------------------------------
omicsStandInReaders()
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
