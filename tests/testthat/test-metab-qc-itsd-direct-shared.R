# fidelity-coverage-compare: shared
library(testthat)

renderDirectMetabItsdUi <- function(ui) {
  htmltools::renderTags(ui)$html
}

makeMetabItsdObject <- function(pattern = "^IS_") {
  methods::new(
    "MetaboliteAssayData",
    metabolite_data = list(
      LCMS_Pos = data.frame(
        database_identifier = c("DB001", "DB002"),
        Name = c("IS_alpha", "Met One"),
        S1 = c(100, 5),
        S2 = c(110, 6),
        check.names = FALSE,
        stringsAsFactors = FALSE
      ),
      LCMS_Neg = data.frame(
        database_identifier = c("DB101", "DB102"),
        Name = c("IS_beta", "Met Two"),
        S1 = c(90, 7),
        S2 = c(95, 8),
        check.names = FALSE,
        stringsAsFactors = FALSE
      )
    ),
    metabolite_id_column = "database_identifier",
    annotation_id_column = "Name",
    database_identifier_type = "database_identifier",
    internal_standard_regex = pattern,
    design_matrix = data.frame(
      Run = c("S1", "S2"),
      group = c("A", "B"),
      stringsAsFactors = FALSE
    ),
    sample_id = "Run",
    group_id = "group",
    technical_replicate_id = NA_character_,
    args = list()
  )
}

fakeItsdMetrics <- function(assay_data, is_pattern, metabolite_id_col, sample_id_col) {
  ids <- assay_data[[metabolite_id_col]]
  matches <- grepl(is_pattern, ids)

  if (!any(matches)) {
    return(data.frame(
      is_id = character(),
      mean_intensity = numeric(),
      cv = numeric(),
      stringsAsFactors = FALSE
    ))
  }

  sample_columns <- c("S1", "S2")
  intensities <- as.matrix(assay_data[matches, sample_columns, drop = FALSE])

  data.frame(
    is_id = ids[matches],
    mean_intensity = rowMeans(intensities),
    cv = c(8, 12)[seq_len(sum(matches))],
    stringsAsFactors = FALSE
  )
}

fakeMetabQuantData <- function(assay_data) {
  list(sample_names = c("S1", "S2"))
}

makeReactiveVal <- function(value = NULL) {
  force(value)

  function(new_value) {
    if (missing(new_value)) {
      return(value)
    }

    value <<- new_value
    invisible(new_value)
  }
}

evaluateRenderExpression <- function(expr) {
  eval(substitute(expr), parent.frame())
}

test_that("metabolomics ITSD analysis preserves assay aggregation and fallback ID scanning", {
  captured_columns <- character()

  analysis <- analyzeMetabQcItsdData(
    currentS4 = makeMetabItsdObject(),
    inputPattern = "^IS_",
    getInternalStandardMetricsFn = function(assay_data, is_pattern, metabolite_id_col, sample_id_col) {
      captured_columns <<- c(captured_columns, metabolite_id_col)
      fakeItsdMetrics(assay_data, is_pattern, metabolite_id_col, sample_id_col)
    },
    getMetaboliteQuantDataFn = fakeMetabQuantData,
    logInfoFn = function(...) invisible(NULL)
  )

  expect_equal(captured_columns, c("database_identifier", "Name", "database_identifier", "Name"))
  expect_equal(analysis$nIsTotal, 2L)
  expect_equal(sort(unique(analysis$metrics$assay)), c("LCMS_Neg", "LCMS_Pos"))
  expect_equal(sort(unique(analysis$longData$assay)), c("LCMS_Neg", "LCMS_Pos"))
  expect_match(analysis$resultText, "Internal Standard Analysis Complete", fixed = TRUE)
  expect_match(analysis$resultText, "LCMS_Pos: 1 internal standards", fixed = TRUE)
  expect_match(analysis$resultText, "LCMS_Neg: 1 internal standards", fixed = TRUE)
})

test_that("metabolomics ITSD analysis preserves current validation errors", {
  expect_error(
    analyzeMetabQcItsdData(
      currentS4 = list(),
      inputPattern = "^IS_"
    ),
    "Current state is not a MetaboliteAssayData object",
    fixed = TRUE
  )

  expect_error(
    analyzeMetabQcItsdData(
      currentS4 = makeMetabItsdObject(pattern = NA_character_),
      inputPattern = NULL
    ),
    "No internal standard pattern provided",
    fixed = TRUE
  )
})

test_that("metabolomics ITSD summary, tabs, and plots preserve current rendering contracts", {
  metrics <- data.frame(
    is_id = c("IS_alpha", "IS_beta"),
    mean_intensity = c(105, 92.5),
    cv = c(8, 28),
    assay = c("LCMS_Pos", "LCMS_Neg"),
    stringsAsFactors = FALSE
  )
  long_data <- data.frame(
    IS_ID = rep(c("IS_alpha", "IS_beta"), each = 2),
    Sample = rep(c("S1", "S2"), times = 2),
    Intensity = c(100, 110, 90, 95),
    assay = rep(c("LCMS_Pos", "LCMS_Neg"), each = 2),
    stringsAsFactors = FALSE
  )

  empty_summary_html <- renderDirectMetabItsdUi(
    buildMetabQcItsdSummaryUi(metrics = NULL)
  )
  expect_match(empty_summary_html, "Click 'Analyze' to detect internal standards.", fixed = TRUE)

  summary_html <- renderDirectMetabItsdUi(
    buildMetabQcItsdSummaryUi(metrics = metrics)
  )
  expect_match(summary_html, "LCMS_Pos: 1 IS (median CV: 8.0%)", fixed = TRUE)
  expect_match(summary_html, "LCMS_Neg: 1 IS (median CV: 28.0%)", fixed = TRUE)

  expect_null(
    buildMetabQcItsdVizTabsUi(
      metrics = NULL,
      ns = function(id) paste0("itsd-", id)
    )
  )

  tabs_html <- renderDirectMetabItsdUi(
    buildMetabQcItsdVizTabsUi(
      metrics = metrics,
      ns = function(id) paste0("itsd-", id)
    )
  )
  expect_match(tabs_html, "CV Distribution", fixed = TRUE)
  expect_match(tabs_html, "itsd-cv_plot", fixed = TRUE)
  expect_match(tabs_html, "itsd-intensity_plot", fixed = TRUE)

  expect_s3_class(buildMetabQcItsdCvPlot(metrics), "ggplot")
  expect_s3_class(buildMetabQcItsdIntensityPlot(long_data), "ggplot")
})

test_that("metabolomics ITSD server body preserves notification and output wiring", {
  current_s4 <- makeMetabItsdObject()
  notifications <- list()
  output <- new.env(parent = emptyenv())

  runMetabQcItsdServerBody(
    input = list(analyze_is = TRUE, is_pattern = "^IS_"),
    output = output,
    session = list(ns = function(id) paste0("itsd-", id)),
    workflowData = list(
      state_manager = list(
        getState = function() current_s4
      )
    ),
    omicType = "metabolomics",
    experimentLabel = "shared-test",
    reactiveValFn = makeReactiveVal,
    observeEventFn = function(eventExpr, handlerExpr, ...) {
      if (isTRUE(eval(substitute(eventExpr), parent.frame()))) {
        eval(substitute(handlerExpr), parent.frame())
      }
      invisible(NULL)
    },
    reqFn = function(...) {
      args <- list(...)
      if (!length(args)) {
        return(invisible(NULL))
      }
      args[[1L]]
    },
    showNotificationFn = function(message, ...) {
      notifications[[length(notifications) + 1L]] <<- list(message = message, ...)
      invisible(NULL)
    },
    removeNotificationFn = function(id) {
      notifications[[length(notifications) + 1L]] <<- list(removed = id)
      invisible(NULL)
    },
    renderTextFn = evaluateRenderExpression,
    renderUiFn = evaluateRenderExpression,
    renderPlotFn = evaluateRenderExpression,
    analyzeMetabQcItsdDataFn = function(currentS4, inputPattern) {
      analyzeMetabQcItsdData(
        currentS4 = currentS4,
        inputPattern = inputPattern,
        getInternalStandardMetricsFn = fakeItsdMetrics,
        getMetaboliteQuantDataFn = fakeMetabQuantData,
        logInfoFn = function(...) invisible(NULL)
      )
    },
    logInfoFn = function(...) invisible(NULL),
    logErrorFn = function(...) invisible(NULL)
  )

  expect_match(output$is_results, "Internal Standard Analysis Complete", fixed = TRUE)
  expect_match(renderDirectMetabItsdUi(output$is_summary), "LCMS_Pos", fixed = TRUE)
  expect_match(renderDirectMetabItsdUi(output$is_viz_tabs), "itsd-cv_plot", fixed = TRUE)
  expect_s3_class(output$cv_plot, "ggplot")
  expect_s3_class(output$intensity_plot, "ggplot")
  expect_length(notifications, 3L)
  expect_identical(notifications[[2]]$removed, "is_analysis_working")
  expect_match(notifications[[3]]$message, "Found 2 internal standards", fixed = TRUE)
})
