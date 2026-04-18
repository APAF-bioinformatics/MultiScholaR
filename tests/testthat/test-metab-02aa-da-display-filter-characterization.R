library(testthat)

repo_root <- normalizePath(file.path("..", ".."), mustWork = TRUE)

loadSelectedFunctions <- function(paths, symbols, env) {
  for (path in paths[file.exists(paths)]) {
    exprs <- parse(file = path, keep.source = TRUE)

    for (expr in exprs) {
      is_assignment <- is.call(expr) &&
        length(expr) >= 3 &&
        as.character(expr[[1]]) %in% c("<-", "=")

      if (!is_assignment || !is.symbol(expr[[2]])) {
        next
      }

      symbol_name <- as.character(expr[[2]])
      if (symbol_name %in% symbols) {
        eval(expr, envir = env)
      }
    }
  }
}

loadSelectedFunctions(
  paths = file.path(repo_root, "R", c(
    "mod_metab_da_display_helpers.R",
    "mod_metab_da_observer_helpers.R",
    "mod_metab_da_registration_helpers.R",
    "mod_metab_da_server_helpers.R",
    "mod_metab_da_server.R",
    "mod_metab_da.R"
  )),
  symbols = c(
    "filterMetabDaDisplayResults",
    "summarizeMetabDaDisplayResults",
    "buildMetabDaSummaryStatsText",
    "buildMetabDaSummaryStatsRenderOutput",
    "shapeMetabDaTableDisplayResults",
    "buildMetabDaResultsDatatable",
    "buildMetabDaResultsTableRenderOutput",
    "buildMetabDaClusterSummaryText",
    "buildMetabDaClusterSummaryRenderOutput",
    "buildMetabDaHeatmapManualSaveWarning",
    "buildMetabDaHeatmapManualSaveWarningRenderOutput",
    "buildMetabDaContrastsDisplayText",
    "buildMetabDaContrastsRenderOutput",
    "buildMetabDaStatusText",
    "buildMetabDaStatusRenderOutput",
    "registerMetabDaOverviewOutputs",
    "buildMetabDaGlimmaCombinedViewInfoBanner",
    "resolveMetabDaGlimmaWidgetOutput",
    "buildMetabDaGlimmaErrorBanner",
    "buildMetabDaGlimmaRenderOutput",
    "buildMetabDaStaticVolcanoPlot",
    "buildMetabDaStaticVolcanoRenderOutput",
    "buildMetabDaHeatmapPlotOutput",
    "buildMetabDaHeatmapRenderOutput",
    "registerMetabDaVisualizationOutputs",
    "registerMetabDaResultsOutputs",
    "registerMetabDaSaveHeatmapObserver",
    "registerMetabDaLoadFilteredSessionObserver",
    "registerMetabDaRunAnalysisObserver",
    "registerMetabDaMainOutputsAndObservers",
    "createMetabDaServerState",
    "registerMetabDaServerBodyOutputs",
    "startMetabDaServerBody",
    "initializeMetabDaServerBody",
    "createMetabDaServerModuleHandler",
    "runMetabDaServerModuleCallback",
    "runMetabDaServerModuleShell",
    "runMetabDaServerModule",
    "runMetabDaServerEntry",
    "mod_metab_da_server",
    "runMetabDaSaveHeatmapObserverShell",
    "runMetabDaSaveHeatmapObserverEntry",
    "updateMetabDaResultsSelectorInputs",
    "writeMetabDaResultsArtifacts",
    "writeMetabDaResultsDownloadCsv",
    "buildMetabDaResultsDownloadHandler",
    "resolveMetabDaAnalysisInputs",
    "runMetabDaAnalysisObserverEntry",
    "runMetabDaAnalysisObserverShell",
    "resolveMetabDaLoadSessionFile",
    "runMetabDaLoadSessionObserverEntry",
    "runMetabDaLoadSessionObserverShell",
    "restoreMetabDaLoadedSessionState"
  ),
  env = environment()
)

helper_env <- environment()

if (!methods::isClass("MetaboliteAssayData")) {
  methods::setClass("MetaboliteAssayData", slots = c(args = "list"))
}

buildMetabDaDisplayResults <- function() {
  data.frame(
    metabolite_id = c("M1", "M2", "M3", "M4"),
    metabolite_name = c("Met One", "Met Two", "Met Three", "Met Four"),
    comparison = c("groupB-groupA", "groupB-groupA", "groupC-groupA", "groupB-groupA"),
    friendly_name = c("B vs A", "B vs A", "C vs A", "B vs A"),
    assay = c("LCMS_Pos", "LCMS_Neg", "LCMS_Pos", "LCMS_Pos"),
    logFC = c(1.2, -1.1, -0.9, 0.2),
    raw_pvalue = c(0.001, 0.004, 0.02, 0.3),
    fdr_qvalue = c(0.01, 0.02, 0.03, 0.20),
    significant = c("Up", "Down", "Down", "NS"),
    stringsAsFactors = FALSE
  )
}

test_that("metabolomics DA display filter preserves contrast and assay narrowing", {
  results <- buildMetabDaDisplayResults()

  filtered <- filterMetabDaDisplayResults(
    results = results,
    selectedContrast = "B vs A",
    selectedAssay = "LCMS_Pos"
  )

  expect_equal(filtered$metabolite_id, c("M1", "M4"))
  expect_true(all(filtered$friendly_name == "B vs A"))
  expect_true(all(filtered$assay == "LCMS_Pos"))
})

test_that("metabolomics DA display filter preserves significance branches and row limits", {
  results <- buildMetabDaDisplayResults()

  significant <- filterMetabDaDisplayResults(
    results = results,
    significanceFilter = "significant",
    daQValThresh = 0.05
  )
  expect_equal(significant$metabolite_id, c("M1", "M2", "M3"))

  up <- filterMetabDaDisplayResults(
    results = results,
    significanceFilter = "up",
    daQValThresh = 0.05,
    treatLfcCutoff = 0.5
  )
  expect_equal(up$metabolite_id, "M1")

  down <- filterMetabDaDisplayResults(
    results = results,
    significanceFilter = "down",
    daQValThresh = 0.05,
    treatLfcCutoff = 0.5
  )
  expect_equal(down$metabolite_id, c("M2", "M3"))

  limited <- filterMetabDaDisplayResults(
    results = results,
    maxRows = 2
  )
  expect_equal(limited$metabolite_id, c("M1", "M2"))
})

test_that("metabolomics DA display filter preserves null and empty exits", {
  expect_null(filterMetabDaDisplayResults(NULL))

  empty <- buildMetabDaDisplayResults()[0, ]
  filtered_empty <- filterMetabDaDisplayResults(
    results = empty,
    selectedContrast = "B vs A",
    significanceFilter = "significant",
    maxRows = 5
  )

  expect_s3_class(filtered_empty, "data.frame")
  expect_equal(nrow(filtered_empty), 0)
})

test_that("metabolomics DA display summary seam preserves counts and threshold", {
  results <- buildMetabDaDisplayResults()

  summary_stats <- summarizeMetabDaDisplayResults(
    results = results,
    daQValThresh = 0.05
  )

  expect_equal(summary_stats$total, 4)
  expect_equal(summary_stats$significant, 3)
  expect_equal(summary_stats$upRegulated, 1)
  expect_equal(summary_stats$downRegulated, 2)
  expect_equal(summary_stats$significantPct, 75)
  expect_equal(summary_stats$qValueThreshold, 0.05)
})

test_that("metabolomics DA display summary seam preserves null and empty exits", {
  expect_null(summarizeMetabDaDisplayResults(NULL))

  empty <- buildMetabDaDisplayResults()[0, ]
  expect_null(summarizeMetabDaDisplayResults(empty))
})

test_that("metabolomics DA summary-stats seam preserves empty and formatted branches", {
  expect_identical(
    buildMetabDaSummaryStatsText(NULL),
    "No results available."
  )

  empty <- buildMetabDaDisplayResults()[0, ]
  expect_identical(
    buildMetabDaSummaryStatsText(empty),
    "No results available."
  )

  summary_text <- buildMetabDaSummaryStatsText(
    results = buildMetabDaDisplayResults(),
    daQValThresh = 0.05
  )

  expect_identical(
    summary_text,
    paste(
      c(
        "Total metabolites: 4",
        "Significant (Q < 0.050): 3 (75.0%)",
        "  Up-regulated: 1",
        "  Down-regulated: 2"
      ),
      collapse = "\n"
    )
  )
})

test_that("metabolomics DA summary-stats render seam preserves req, filter, and text handoff", {
  captured_req <- NULL
  captured_filter <- NULL
  captured_text <- NULL
  da_results_list <- list(da_metabolites_long = "raw-results")

  output <- buildMetabDaSummaryStatsRenderOutput(
    daResultsList = da_results_list,
    selectedContrast = "B vs A",
    selectedAssay = "LCMS_Pos",
    daQValThresh = 0.01,
    requireInputs = function(...) {
      captured_req <<- list(...)
      invisible(TRUE)
    },
    filterResults = function(...) {
      captured_filter <<- list(...)
      "filtered-results"
    },
    buildSummaryText = function(results, daQValThresh = 0.05) {
      captured_text <<- list(results = results, daQValThresh = daQValThresh)
      "summary-stats-output"
    },
    writeOutput = function(text, ...) {
      paste0("written:", text)
    }
  )

  expect_identical(captured_req, list(da_results_list))
  expect_identical(captured_filter$results, "raw-results")
  expect_identical(captured_filter$selectedContrast, "B vs A")
  expect_identical(captured_filter$selectedAssay, "LCMS_Pos")
  expect_identical(captured_text$results, "filtered-results")
  expect_identical(captured_text$daQValThresh, 0.01)
  expect_identical(output, "written:summary-stats-output")
})

test_that("metabolomics DA table display seam preserves display columns and row order", {
  results <- buildMetabDaDisplayResults()
  results$extra_debug <- c("x1", "x2", "x3", "x4")

  display_results <- shapeMetabDaTableDisplayResults(results)

  expect_equal(
    colnames(display_results),
    c(
      "metabolite_id", "metabolite_name", "assay", "logFC",
      "raw_pvalue", "fdr_qvalue", "significant"
    )
  )
  expect_equal(display_results$metabolite_id, results$metabolite_id)
  expect_false("extra_debug" %in% colnames(display_results))
})

test_that("metabolomics DA table display seam preserves null and empty exits", {
  expect_null(shapeMetabDaTableDisplayResults(NULL))

  empty <- buildMetabDaDisplayResults()[0, ]
  display_empty <- shapeMetabDaTableDisplayResults(empty)

  expect_s3_class(display_empty, "data.frame")
  expect_equal(nrow(display_empty), 0)
  expect_equal(
    colnames(display_empty),
    c(
      "metabolite_id", "metabolite_name", "assay", "logFC",
      "raw_pvalue", "fdr_qvalue", "significant"
    )
  )
})

test_that("metabolomics DA datatable seam preserves display widget options and formatting", {
  results <- buildMetabDaDisplayResults()
  display_results <- shapeMetabDaTableDisplayResults(results)

  table_widget <- buildMetabDaResultsDatatable(display_results)
  column_defs <- table_widget$x$options$columnDefs
  named_defs <- Filter(function(def) !is.null(def$name), column_defs)

  expect_s3_class(table_widget, "datatables")
  expect_equal(table_widget$x$filter, "top")
  expect_equal(unlist(table_widget$x$extensions), "Buttons")
  expect_equal(table_widget$x$options$pageLength, 25)
  expect_true(isTRUE(table_widget$x$options$scrollX))
  expect_equal(table_widget$x$options$dom, "Bfrtip")
  expect_equal(unlist(table_widget$x$options$buttons), c("copy", "csv", "excel"))
  expect_setequal(
    vapply(named_defs, `[[`, character(1), "name"),
    c(
      "metabolite_id", "metabolite_name", "assay", "logFC",
      "raw_pvalue", "fdr_qvalue", "significant"
    )
  )
  expect_match(table_widget$x$options$rowCallback, "#ffcccc", fixed = TRUE)
  expect_match(table_widget$x$options$rowCallback, "#cce5ff", fixed = TRUE)
})

test_that("metabolomics DA datatable seam preserves null and empty exits", {
  expect_null(buildMetabDaResultsDatatable(NULL))

  empty <- buildMetabDaDisplayResults()[0, ]
  display_empty <- shapeMetabDaTableDisplayResults(empty)

  expect_null(buildMetabDaResultsDatatable(display_empty))
})

test_that("metabolomics DA results-table render seam preserves req, filter, shape, and datatable handoff", {
  captured_req <- NULL
  captured_filter <- NULL
  captured_shape <- NULL
  captured_datatable <- NULL
  da_results_list <- list(da_metabolites_long = "raw-results")

  output <- buildMetabDaResultsTableRenderOutput(
    daResultsList = da_results_list,
    selectedContrast = "B vs A",
    selectedAssay = "LCMS_Pos",
    significanceFilter = "up",
    daQValThresh = 0.01,
    treatLfcCutoff = 1.5,
    maxRows = 7,
    requireInputs = function(...) {
      captured_req <<- list(...)
      invisible(TRUE)
    },
    filterResults = function(...) {
      captured_filter <<- list(...)
      data.frame(filtered = "yes", stringsAsFactors = FALSE)
    },
    shapeResults = function(results) {
      captured_shape <<- results
      data.frame(shaped = "yes", stringsAsFactors = FALSE)
    },
    buildDatatable = function(results) {
      captured_datatable <<- results
      "datatable-output"
    }
  )

  expect_identical(captured_req, list(da_results_list))
  expect_identical(captured_filter$results, "raw-results")
  expect_identical(captured_filter$selectedContrast, "B vs A")
  expect_identical(captured_filter$selectedAssay, "LCMS_Pos")
  expect_identical(captured_filter$significanceFilter, "up")
  expect_identical(captured_filter$daQValThresh, 0.01)
  expect_identical(captured_filter$treatLfcCutoff, 1.5)
  expect_identical(captured_filter$maxRows, 7)
  expect_identical(captured_shape, data.frame(filtered = "yes", stringsAsFactors = FALSE))
  expect_identical(captured_datatable, data.frame(shaped = "yes", stringsAsFactors = FALSE))
  expect_identical(output, "datatable-output")
})

test_that("metabolomics DA cluster summary seam preserves grouped text output", {
  clusters <- c(MetC = 2, MetA = 1, MetB = 2)

  summary_text <- buildMetabDaClusterSummaryText(clusters)

  expect_identical(
    summary_text,
    paste0(
      "Total Clusters: 2\n",
      "----------------------------------\n",
      "\nCluster 1 (1 metabolites):\n",
      "MetA\n",
      "\nCluster 2 (2 metabolites):\n",
      "MetC, MetB\n"
    )
  )
})

test_that("metabolomics DA cluster summary render seam preserves req gate and text handoff", {
  captured_req <- NULL
  captured_text <- NULL

  output <- buildMetabDaClusterSummaryRenderOutput(
    treeCutMethod = "dynamic",
    clusters = c(MetA = 1, MetB = 2),
    requireInputs = function(...) {
      captured_req <<- list(...)
      invisible(TRUE)
    },
    buildSummaryText = function(clusters, maxMembers = 20) {
      captured_text <<- list(clusters = clusters, maxMembers = maxMembers)
      "cluster-summary-output"
    },
    writeOutput = function(text, ...) {
      paste0("written:", text)
    }
  )

  expect_identical(captured_req, list(TRUE))
  expect_identical(captured_text$clusters, c(MetA = 1, MetB = 2))
  expect_identical(captured_text$maxMembers, 20)
  expect_identical(output, "written:cluster-summary-output")
})

test_that("metabolomics DA heatmap warning seam preserves null and banner exits", {
  expect_null(buildMetabDaHeatmapManualSaveWarning(FALSE))
  expect_null(buildMetabDaHeatmapManualSaveWarning(NULL))

  warning_banner <- buildMetabDaHeatmapManualSaveWarning(TRUE)

  expect_s3_class(warning_banner, "shiny.tag")
  expect_identical(warning_banner$name, "div")
  expect_true(grepl("alert alert-info", warning_banner$attribs$class, fixed = TRUE))
  expect_true(grepl("margin-bottom: 15px", warning_banner$attribs$style, fixed = TRUE))

  banner_text <- paste(capture.output(print(warning_banner)), collapse = "\n")
  expect_match(banner_text, "Heatmaps are NOT saved automatically during analysis.", fixed = TRUE)
  expect_match(banner_text, "Heatmap", fixed = TRUE)
  expect_match(banner_text, "Save Heatmap", fixed = TRUE)
})

test_that("metabolomics DA heatmap warning render seam preserves warning helper handoff", {
  captured_analysis_complete <- NULL

  output <- buildMetabDaHeatmapManualSaveWarningRenderOutput(
    analysisComplete = TRUE,
    buildWarning = function(analysisComplete) {
      captured_analysis_complete <<- analysisComplete
      "warning-banner-output"
    }
  )

  expect_identical(captured_analysis_complete, TRUE)
  expect_identical(output, "warning-banner-output")
})

test_that("metabolomics DA contrasts seam preserves empty, named, and fallback branches", {
  expect_identical(
    buildMetabDaContrastsDisplayText(NULL),
    "No contrasts defined.\nLoad a filtered session or define contrasts in the design tab."
  )

  empty <- data.frame(friendly_names = character(), stringsAsFactors = FALSE)
  expect_identical(
    buildMetabDaContrastsDisplayText(empty),
    "No contrasts defined.\nLoad a filtered session or define contrasts in the design tab."
  )

  friendly_tbl <- data.frame(
    friendly_names = c("B vs A", "C vs A"),
    contrasts = c("groupB-groupA", "groupC-groupA"),
    stringsAsFactors = FALSE
  )
  expect_identical(
    buildMetabDaContrastsDisplayText(friendly_tbl),
    "B vs A\nC vs A"
  )

  contrasts_tbl <- data.frame(
    contrasts = c("groupB-groupA", "groupC-groupA"),
    stringsAsFactors = FALSE
  )
  expect_identical(
    buildMetabDaContrastsDisplayText(contrasts_tbl),
    "groupB-groupA\ngroupC-groupA"
  )

  fallback_tbl <- data.frame(
    label = c("First", "Second"),
    value = c(1, 2),
    stringsAsFactors = FALSE
  )
  expect_identical(
    buildMetabDaContrastsDisplayText(fallback_tbl),
    paste(capture.output(print(fallback_tbl)), collapse = "\n")
  )
})

test_that("metabolomics DA contrasts render seam preserves display-text handoff", {
  captured_tbl <- NULL

  output <- buildMetabDaContrastsRenderOutput(
    contrastsTbl = "raw-contrasts",
    buildDisplayText = function(contrastsTbl) {
      captured_tbl <<- contrastsTbl
      "formatted-contrasts"
    },
    writeOutput = function(text, ...) {
      paste0("written:", text)
    }
  )

  expect_identical(captured_tbl, "raw-contrasts")
  expect_identical(output, "written:formatted-contrasts")
})

test_that("metabolomics DA status seam preserves analysis, ready, and waiting branches", {
  expect_identical(
    buildMetabDaStatusText(
      analysisComplete = TRUE,
      significantCounts = list(
        LCMS_Pos = list(up = 3, down = 1, ns = 12),
        LCMS_Neg = list(up = 0, down = 2, ns = 8)
      )
    ),
    paste(
      c(
        "Analysis Complete",
        "",
        "LCMS_Pos:\n  Up: 3 | Down: 1 | NS: 12",
        "LCMS_Neg:\n  Up: 0 | Down: 2 | NS: 8"
      ),
      collapse = "\n"
    )
  )

  expect_identical(
    buildMetabDaStatusText(analysisComplete = TRUE),
    "Analysis complete."
  )

  expect_identical(
    buildMetabDaStatusText(currentS4Object = structure(list(), class = "MetaboliteAssayData")),
    "Data loaded. Ready to run analysis."
  )

  expect_identical(
    buildMetabDaStatusText(),
    "Waiting for data.\nClick 'Load Filtered Session' to begin."
  )
})

test_that("metabolomics DA status render seam preserves status-count extraction and text handoff", {
  captured_text <- NULL

  output <- buildMetabDaStatusRenderOutput(
    daResultsList = list(
      significant_counts = list(
        LCMS_Pos = list(up = 2, down = 1, ns = 5)
      )
    ),
    analysisComplete = TRUE,
    currentS4Object = "s4-state",
    buildStatusText = function(analysisComplete, currentS4Object, significantCounts) {
      captured_text <<- list(
        analysisComplete = analysisComplete,
        currentS4Object = currentS4Object,
        significantCounts = significantCounts
      )
      "status-output"
    },
    writeOutput = function(text, ...) {
      paste0("written:", text)
    }
  )

  expect_identical(captured_text$analysisComplete, TRUE)
  expect_identical(captured_text$currentS4Object, "s4-state")
  expect_identical(
    captured_text$significantCounts,
    list(LCMS_Pos = list(up = 2, down = 1, ns = 5))
  )
  expect_identical(output, "written:status-output")
})

test_that("metabolomics DA overview-output seam preserves contrasts and status registration fan-out", {
  output <- new.env(parent = emptyenv())
  workflow_state <- new.env(parent = emptyenv())
  workflow_state$contrasts_tbl <- "contrast-state"
  da_state <- new.env(parent = emptyenv())
  da_state$da_results_list <- "results-state"
  da_state$analysis_complete <- TRUE
  da_state$current_s4_object <- "current-s4"

  captured_contrasts <- NULL
  captured_status <- NULL

  registration <- registerMetabDaOverviewOutputs(
    output = output,
    workflowData = workflow_state,
    daData = da_state,
    renderPrint = function(expr) list(kind = "print", value = expr),
    buildContrastsOutput = function(...) {
      captured_contrasts <<- list(...)
      "contrasts-output"
    },
    buildStatusOutput = function(...) {
      captured_status <<- list(...)
      "status-output"
    }
  )

  expect_identical(registration, output)
  expect_identical(captured_contrasts$contrastsTbl, "contrast-state")
  expect_identical(output$contrasts_display, list(kind = "print", value = "contrasts-output"))
  expect_identical(captured_status$daResultsList, "results-state")
  expect_identical(captured_status$analysisComplete, TRUE)
  expect_identical(captured_status$currentS4Object, "current-s4")
  expect_identical(output$da_status, list(kind = "print", value = "status-output"))
})

test_that("metabolomics DA Glimma combined-view seam preserves assay gate and banner text", {
  expect_null(buildMetabDaGlimmaCombinedViewInfoBanner("LCMS_Pos"))

  combined_banner <- buildMetabDaGlimmaCombinedViewInfoBanner("Combined")
  null_banner <- buildMetabDaGlimmaCombinedViewInfoBanner(NULL)

  expect_s3_class(combined_banner, "shiny.tag")
  expect_identical(combined_banner$name, "div")
  expect_true(grepl("alert alert-info", combined_banner$attribs$class, fixed = TRUE))
  expect_true(grepl("margin: 20px", combined_banner$attribs$style, fixed = TRUE))
  expect_identical(combined_banner, null_banner)

  banner_text <- paste(capture.output(print(combined_banner)), collapse = "\n")
  expect_match(banner_text, "Combined View", fixed = TRUE)
  expect_match(
    banner_text,
    "Interactive Glimma plots require a single assay selection.",
    fixed = TRUE
  )
  expect_match(banner_text, "static volcano plot", fixed = TRUE)
})

test_that("metabolomics DA Glimma widget seam preserves warning fallback and widget passthrough", {
  warning_banner <- resolveMetabDaGlimmaWidgetOutput(NULL)
  widget <- shiny::div(class = "glimma-holder", "ready")

  expect_s3_class(warning_banner, "shiny.tag")
  expect_identical(warning_banner$name, "div")
  expect_true(grepl("alert alert-warning", warning_banner$attribs$class, fixed = TRUE))

  banner_text <- paste(capture.output(print(warning_banner)), collapse = "\n")
  expect_match(banner_text, "Could not generate Glimma plot.", fixed = TRUE)
  expect_match(banner_text, "static plot option", fixed = TRUE)

  expect_identical(resolveMetabDaGlimmaWidgetOutput(widget), widget)
})

test_that("metabolomics DA Glimma error seam preserves danger banner text", {
  error_banner <- buildMetabDaGlimmaErrorBanner("plot backend failed")

  expect_s3_class(error_banner, "shiny.tag")
  expect_identical(error_banner$name, "div")
  expect_true(grepl("alert alert-danger", error_banner$attribs$class, fixed = TRUE))

  banner_text <- paste(capture.output(print(error_banner)), collapse = "\n")
  expect_match(banner_text, "Error generating plot: plot backend failed", fixed = TRUE)
})

test_that("metabolomics DA Glimma render seam preserves combined-view short circuit", {
  captured_req <- NULL
  generated <- FALSE
  warning_banner <- shiny::div(class = "combined-banner", "combined-only")

  output <- buildMetabDaGlimmaRenderOutput(
    daResultsList = "results-state",
    selectedContrast = "B vs A",
    selectedAssay = "Combined",
    requireInputs = function(...) {
      captured_req <<- list(...)
      invisible(TRUE)
    },
    buildCombinedViewInfo = function(selectedAssay) {
      expect_identical(selectedAssay, "Combined")
      warning_banner
    },
    generatePlot = function(...) {
      generated <<- TRUE
      stop("should not generate in combined view")
    }
  )

  expect_identical(captured_req, list("results-state", "B vs A"))
  expect_false(generated)
  expect_identical(output, warning_banner)
})

test_that("metabolomics DA Glimma render seam preserves generator handoff and widget resolution", {
  captured_req <- NULL
  captured_generate <- NULL
  captured_widget <- NULL

  output <- buildMetabDaGlimmaRenderOutput(
    daResultsList = "results-state",
    selectedContrast = "B vs A",
    selectedAssay = "LCMS_Pos",
    daQValThresh = 0.01,
    requireInputs = function(...) {
      captured_req <<- list(...)
      invisible(TRUE)
    },
    buildCombinedViewInfo = function(...) NULL,
    generatePlot = function(...) {
      captured_generate <<- list(...)
      "glimma-widget"
    },
    resolveWidgetOutput = function(widget) {
      captured_widget <<- widget
      list(resolved = widget)
    }
  )

  expect_identical(captured_req, list("results-state", "B vs A"))
  expect_identical(captured_generate$da_results_list, "results-state")
  expect_identical(captured_generate$selected_contrast, "B vs A")
  expect_identical(captured_generate$selected_assay, "LCMS_Pos")
  expect_identical(captured_generate$da_q_val_thresh, 0.01)
  expect_identical(captured_widget, "glimma-widget")
  expect_identical(output, list(resolved = "glimma-widget"))
})

test_that("metabolomics DA Glimma render seam preserves error logging and fallback banner", {
  captured_log <- NULL

  output <- buildMetabDaGlimmaRenderOutput(
    daResultsList = "results-state",
    selectedContrast = "B vs A",
    selectedAssay = "LCMS_Pos",
    requireInputs = function(...) invisible(TRUE),
    buildCombinedViewInfo = function(...) NULL,
    generatePlot = function(...) stop("backend failed"),
    buildErrorBanner = function(message) list(type = "banner", message = message),
    logError = function(message) {
      captured_log <<- message
      invisible(NULL)
    }
  )

  expect_identical(captured_log, "Glimma error: backend failed")
  expect_identical(
    output,
    list(type = "banner", message = "backend failed")
  )
})

test_that("metabolomics DA static volcano seam preserves generator inputs and defaults", {
  captured <- NULL
  stub_results <- list(da_metabolites_long = data.frame(dummy = 1))

  old_generator <- get0(
    "generateMetabDAVolcanoStatic",
    envir = helper_env,
    inherits = FALSE
  )
  on.exit({
    if (is.null(old_generator)) {
      rm("generateMetabDAVolcanoStatic", envir = helper_env)
    } else {
      assign("generateMetabDAVolcanoStatic", old_generator, envir = helper_env)
    }
  }, add = TRUE)

  assign(
    "generateMetabDAVolcanoStatic",
    function(...) {
      captured <<- list(...)
      "static-volcano-plot"
    },
    envir = helper_env
  )

  plot <- buildMetabDaStaticVolcanoPlot(
    daResultsList = stub_results,
    selectedContrast = "B vs A",
    selectedAssay = "LCMS_Pos",
    daQValThresh = 0.01,
    treatLfcCutoff = 1.5
  )

  expect_identical(plot, "static-volcano-plot")
  expect_identical(captured$da_results_list, stub_results)
  expect_identical(captured$selected_contrast, "B vs A")
  expect_identical(captured$selected_assay, "LCMS_Pos")
  expect_identical(captured$da_q_val_thresh, 0.01)
  expect_identical(captured$lfc_threshold, 1.5)
  expect_identical(captured$show_labels, TRUE)
  expect_identical(captured$n_labels, 15)
})

test_that("metabolomics DA static volcano render seam preserves req gate and plot handoff", {
  captured_req <- NULL
  captured_plot <- NULL

  output <- buildMetabDaStaticVolcanoRenderOutput(
    daResultsList = "results-state",
    selectedContrast = "B vs A",
    selectedAssay = "LCMS_Pos",
    daQValThresh = 0.01,
    treatLfcCutoff = 1.5,
    requireInputs = function(...) {
      captured_req <<- list(...)
      invisible(TRUE)
    },
    buildPlot = function(...) {
      captured_plot <<- list(...)
      "static-render-output"
    }
  )

  expect_identical(captured_req, list("results-state", "B vs A"))
  expect_identical(captured_plot$daResultsList, "results-state")
  expect_identical(captured_plot$selectedContrast, "B vs A")
  expect_identical(captured_plot$selectedAssay, "LCMS_Pos")
  expect_identical(captured_plot$daQValThresh, 0.01)
  expect_identical(captured_plot$treatLfcCutoff, 1.5)
  expect_identical(output, "static-render-output")
})

test_that("metabolomics DA heatmap seam preserves generator inputs and stored plot state", {
  captured <- NULL
  stub_results <- list(da_metabolites_long = data.frame(dummy = 1))
  da_state <- new.env(parent = emptyenv())

  old_generator <- get0(
    "generateMetabDAHeatmap",
    envir = helper_env,
    inherits = FALSE
  )
  on.exit({
    if (is.null(old_generator)) {
      rm("generateMetabDAHeatmap", envir = helper_env)
    } else {
      assign("generateMetabDAHeatmap", old_generator, envir = helper_env)
    }
  }, add = TRUE)

  assign(
    "generateMetabDAHeatmap",
    function(...) {
      captured <<- list(...)
      list(
        plot = "heatmap-plot",
        row_clusters = c(Met1 = 1, Met2 = 2),
        col_clusters = c(S1 = 1, S2 = 1)
      )
    },
    envir = helper_env
  )

  plot <- buildMetabDaHeatmapPlotOutput(
    daResultsList = stub_results,
    selectedContrast = "B vs A",
    selectedAssay = "LCMS_Pos",
    topN = 25,
    clusteringMethod = "average",
    distanceMethod = "manhattan",
    clusterRows = FALSE,
    clusterCols = TRUE,
    scaleData = "column",
    colorScheme = "Viridis",
    showMetaboliteNames = FALSE,
    daQValThresh = 0.01,
    treeCutMethod = "dynamic",
    nClusters = 4,
    cutHeight = 0.8,
    minClusterSize = 3,
    daData = da_state
  )

  expect_identical(plot, "heatmap-plot")
  expect_identical(captured$da_results_list, stub_results)
  expect_identical(captured$selected_contrast, "B vs A")
  expect_identical(captured$selected_assay, "LCMS_Pos")
  expect_identical(captured$top_n, 25)
  expect_identical(captured$clustering_method, "average")
  expect_identical(captured$distance_method, "manhattan")
  expect_identical(captured$cluster_rows, FALSE)
  expect_identical(captured$cluster_cols, TRUE)
  expect_identical(captured$scale_data, "column")
  expect_identical(captured$color_scheme, "Viridis")
  expect_identical(captured$show_metabolite_names, FALSE)
  expect_identical(captured$da_q_val_thresh, 0.01)
  expect_identical(captured$tree_cut_method, "dynamic")
  expect_identical(captured$n_clusters, 4)
  expect_identical(captured$cut_height, 0.8)
  expect_identical(captured$min_cluster_size, 3)
  expect_identical(da_state$current_row_clusters, c(Met1 = 1, Met2 = 2))
  expect_identical(da_state$current_col_clusters, c(S1 = 1, S2 = 1))
  expect_identical(da_state$current_heatmap_plot, "heatmap-plot")
})

test_that("metabolomics DA heatmap seam preserves passthrough and null exits", {
  old_generator <- get0(
    "generateMetabDAHeatmap",
    envir = helper_env,
    inherits = FALSE
  )
  on.exit({
    if (is.null(old_generator)) {
      rm("generateMetabDAHeatmap", envir = helper_env)
    } else {
      assign("generateMetabDAHeatmap", old_generator, envir = helper_env)
    }
  }, add = TRUE)

  assign(
    "generateMetabDAHeatmap",
    function(...) NULL,
    envir = helper_env
  )
  expect_null(
    buildMetabDaHeatmapPlotOutput(
      daResultsList = list(da_metabolites_long = data.frame(dummy = 1)),
      selectedContrast = "B vs A"
    )
  )

  assign(
    "generateMetabDAHeatmap",
    function(...) "prebuilt-heatmap-plot",
    envir = helper_env
  )
  expect_identical(
    buildMetabDaHeatmapPlotOutput(
      daResultsList = list(da_metabolites_long = data.frame(dummy = 1)),
      selectedContrast = "B vs A"
    ),
    "prebuilt-heatmap-plot"
  )
})

test_that("metabolomics DA heatmap render seam preserves req gate and plot handoff", {
  captured_req <- NULL
  captured_plot <- NULL

  output <- buildMetabDaHeatmapRenderOutput(
    daResultsList = "results-state",
    selectedContrast = "B vs A",
    selectedAssay = "LCMS_Pos",
    topN = 25,
    clusteringMethod = "average",
    distanceMethod = "manhattan",
    heatmapClustering = "column",
    scaleData = "column",
    colorScheme = "Viridis",
    showMetaboliteNames = FALSE,
    daQValThresh = 0.01,
    treeCutMethod = "dynamic",
    nClusters = 4,
    cutHeight = 0.8,
    minClusterSize = 3,
    daData = "da-state",
    requireInputs = function(...) {
      captured_req <<- list(...)
      invisible(TRUE)
    },
    buildPlot = function(...) {
      captured_plot <<- list(...)
      "heatmap-render-output"
    }
  )

  expect_identical(captured_req, list("results-state", "B vs A"))
  expect_identical(captured_plot$daResultsList, "results-state")
  expect_identical(captured_plot$selectedContrast, "B vs A")
  expect_identical(captured_plot$selectedAssay, "LCMS_Pos")
  expect_identical(captured_plot$topN, 25)
  expect_identical(captured_plot$clusteringMethod, "average")
  expect_identical(captured_plot$distanceMethod, "manhattan")
  expect_identical(captured_plot$clusterRows, FALSE)
  expect_identical(captured_plot$clusterCols, TRUE)
  expect_identical(captured_plot$scaleData, "column")
  expect_identical(captured_plot$colorScheme, "Viridis")
  expect_identical(captured_plot$showMetaboliteNames, FALSE)
  expect_identical(captured_plot$daQValThresh, 0.01)
  expect_identical(captured_plot$treeCutMethod, "dynamic")
  expect_identical(captured_plot$nClusters, 4)
  expect_identical(captured_plot$cutHeight, 0.8)
  expect_identical(captured_plot$minClusterSize, 3)
  expect_identical(captured_plot$daData, "da-state")
  expect_identical(output, "heatmap-render-output")
})

test_that("metabolomics DA visualization-output seam preserves registration fan-out", {
  output <- new.env(parent = emptyenv())
  input <- list(
    volcano_contrast = "B vs A",
    volcano_assay = "LCMS_Pos",
    da_q_val_thresh = 0.01,
    treat_lfc_cutoff = 1.5,
    heatmap_contrast = "C vs A",
    heatmap_assay = "LCMS_Neg",
    heatmap_top_n = 25,
    heatmap_cluster_method = "average",
    heatmap_distance_method = "manhattan",
    heatmap_clustering = "column",
    heatmap_scaling = "row",
    heatmap_color_scheme = "Viridis",
    heatmap_show_labels = FALSE,
    heatmap_tree_cut_method = "dynamic",
    heatmap_n_clusters = 4,
    heatmap_cut_height = 0.8,
    heatmap_min_cluster_size = 3
  )
  da_state <- new.env(parent = emptyenv())
  da_state$analysis_complete <- TRUE
  da_state$da_results_list <- "results-state"
  da_state$current_row_clusters <- c(MetA = 1, MetB = 2)

  captured_warning <- NULL
  captured_glimma <- NULL
  captured_static <- NULL
  captured_heatmap <- NULL
  captured_cluster <- NULL

  registration <- registerMetabDaVisualizationOutputs(
    input = input,
    output = output,
    daData = da_state,
    renderUi = function(expr) list(kind = "ui", value = expr),
    renderPlot = function(expr) list(kind = "plot", value = expr),
    renderPrint = function(expr) list(kind = "print", value = expr),
    buildHeatmapWarningOutput = function(analysisComplete) {
      captured_warning <<- analysisComplete
      "warning-output"
    },
    buildGlimmaOutput = function(...) {
      captured_glimma <<- list(...)
      "glimma-output"
    },
    buildStaticVolcanoOutput = function(...) {
      captured_static <<- list(...)
      "static-output"
    },
    buildHeatmapOutput = function(...) {
      captured_heatmap <<- list(...)
      "heatmap-output"
    },
    buildClusterSummaryOutput = function(...) {
      captured_cluster <<- list(...)
      "cluster-output"
    }
  )

  expect_identical(registration, output)
  expect_identical(captured_warning, TRUE)
  expect_identical(output$heatmap_manual_save_warning, list(kind = "ui", value = "warning-output"))
  expect_identical(captured_glimma$daResultsList, "results-state")
  expect_identical(captured_glimma$selectedContrast, "B vs A")
  expect_identical(captured_glimma$selectedAssay, "LCMS_Pos")
  expect_identical(captured_glimma$daQValThresh, 0.01)
  expect_identical(output$volcano_glimma, list(kind = "ui", value = "glimma-output"))
  expect_identical(captured_static$daResultsList, "results-state")
  expect_identical(captured_static$selectedContrast, "B vs A")
  expect_identical(captured_static$selectedAssay, "LCMS_Pos")
  expect_identical(captured_static$daQValThresh, 0.01)
  expect_identical(captured_static$treatLfcCutoff, 1.5)
  expect_identical(output$volcano_static, list(kind = "plot", value = "static-output"))
  expect_identical(captured_heatmap$daResultsList, "results-state")
  expect_identical(captured_heatmap$selectedContrast, "C vs A")
  expect_identical(captured_heatmap$selectedAssay, "LCMS_Neg")
  expect_identical(captured_heatmap$topN, 25)
  expect_identical(captured_heatmap$clusteringMethod, "average")
  expect_identical(captured_heatmap$distanceMethod, "manhattan")
  expect_identical(captured_heatmap$heatmapClustering, "column")
  expect_identical(captured_heatmap$scaleData, "row")
  expect_identical(captured_heatmap$colorScheme, "Viridis")
  expect_identical(captured_heatmap$showMetaboliteNames, FALSE)
  expect_identical(captured_heatmap$daQValThresh, 0.01)
  expect_identical(captured_heatmap$treeCutMethod, "dynamic")
  expect_identical(captured_heatmap$nClusters, 4)
  expect_identical(captured_heatmap$cutHeight, 0.8)
  expect_identical(captured_heatmap$minClusterSize, 3)
  expect_identical(captured_heatmap$daData, da_state)
  expect_identical(output$heatmap_plot, list(kind = "plot", value = "heatmap-output"))
  expect_identical(captured_cluster$treeCutMethod, "dynamic")
  expect_identical(captured_cluster$clusters, c(MetA = 1, MetB = 2))
  expect_identical(output$cluster_summary, list(kind = "print", value = "cluster-output"))
})

test_that("metabolomics DA results-output seam preserves registration fan-out", {
  output <- new.env(parent = emptyenv())
  input <- list(
    table_contrast = "B vs A",
    table_assay = "LCMS_Pos",
    table_significance = "up",
    da_q_val_thresh = 0.01,
    treat_lfc_cutoff = 1.5,
    table_max_rows = 12
  )
  da_state <- new.env(parent = emptyenv())
  da_state$da_results_list <- "results-state"

  captured_summary <- NULL
  captured_table <- NULL
  captured_download <- NULL

  registration <- registerMetabDaResultsOutputs(
    input = input,
    output = output,
    daData = da_state,
    renderPrint = function(expr) list(kind = "print", value = expr),
    renderDt = function(expr) list(kind = "dt", value = expr),
    buildSummaryStatsOutput = function(...) {
      captured_summary <<- list(...)
      "summary-output"
    },
    buildResultsTableOutput = function(...) {
      captured_table <<- list(...)
      "table-output"
    },
    buildResultsDownloadHandler = function(...) {
      captured_download <<- list(...)
      "download-output"
    }
  )

  expect_identical(registration, output)
  expect_identical(captured_summary$daResultsList, "results-state")
  expect_identical(captured_summary$selectedContrast, "B vs A")
  expect_identical(captured_summary$selectedAssay, "LCMS_Pos")
  expect_identical(captured_summary$daQValThresh, 0.01)
  expect_identical(output$da_summary_stats, list(kind = "print", value = "summary-output"))
  expect_identical(captured_table$daResultsList, "results-state")
  expect_identical(captured_table$selectedContrast, "B vs A")
  expect_identical(captured_table$selectedAssay, "LCMS_Pos")
  expect_identical(captured_table$significanceFilter, "up")
  expect_identical(captured_table$daQValThresh, 0.01)
  expect_identical(captured_table$treatLfcCutoff, 1.5)
  expect_identical(captured_table$maxRows, 12)
  expect_identical(output$da_results_table, list(kind = "dt", value = "table-output"))
  expect_identical(captured_download$daData, da_state)
  expect_identical(output$download_da_results, "download-output")
})

test_that("metabolomics DA save-heatmap observer registration seam preserves event and entry fan-out", {
  input <- list(
    save_heatmap = 1,
    heatmap_contrast = "B vs A",
    heatmap_top_n = 25,
    heatmap_cluster_method = "average",
    heatmap_distance_method = "manhattan",
    heatmap_clustering = "column",
    heatmap_scaling = "row",
    heatmap_color_scheme = "Viridis",
    heatmap_tree_cut_method = "dynamic",
    heatmap_n_clusters = 4,
    heatmap_cut_height = 0.8,
    heatmap_min_cluster_size = 3
  )
  da_state <- new.env(parent = emptyenv())
  da_state$current_heatmap_plot <- "heatmap-object"
  da_state$current_row_clusters <- c(MetA = 1, MetB = 2)
  experiment_paths <- list(publication_graphs_dir = "/tmp/metab-publication-graphs")

  captured_event <- NULL
  captured_entry <- NULL

  registration <- registerMetabDaSaveHeatmapObserver(
    input = input,
    daData = da_state,
    experimentPaths = experiment_paths,
    observeEvent = function(eventExpr, handlerExpr) {
      captured_event <<- eventExpr
      eval(substitute(handlerExpr), envir = parent.frame())
      "observer-registered"
    },
    runSaveHeatmapObserverEntry = function(...) {
      captured_entry <<- list(...)
      "save-entry-output"
    }
  )

  expect_identical(registration, "observer-registered")
  expect_identical(captured_event, 1)
  expect_identical(captured_entry$currentHeatmapPlot, "heatmap-object")
  expect_identical(captured_entry$publicationGraphsDir, "/tmp/metab-publication-graphs")
  expect_identical(captured_entry$currentRowClusters, c(MetA = 1, MetB = 2))
  expect_identical(captured_entry$selectedContrast, "B vs A")
  expect_identical(captured_entry$topN, 25)
  expect_identical(captured_entry$clusteringMethod, "average")
  expect_identical(captured_entry$distanceMethod, "manhattan")
  expect_identical(captured_entry$heatmapClustering, "column")
  expect_identical(captured_entry$scaling, "row")
  expect_identical(captured_entry$colorScheme, "Viridis")
  expect_identical(captured_entry$treeCutMethod, "dynamic")
  expect_identical(captured_entry$nClusters, 4)
  expect_identical(captured_entry$cutHeight, 0.8)
  expect_identical(captured_entry$minClusterSize, 3)
})

test_that("metabolomics DA load-filtered-session observer registration seam preserves event, debug log, and entry handoff", {
  captured_event <- NULL
  captured_entry <- NULL
  captured_log <- NULL
  captured_messages <- character()
  start_time <- as.POSIXct("2026-04-16 15:30:00", tz = "UTC")

  workflow_state <- new.env(parent = emptyenv())
  da_state <- new.env(parent = emptyenv())
  experiment_paths <- list(source_dir = "/tmp/metab-source")

  registration <- registerMetabDaLoadFilteredSessionObserver(
    input = list(load_filtered_session = 3),
    workflowData = workflow_state,
    daData = da_state,
    session = "metab-session",
    experimentPaths = experiment_paths,
    observeEvent = function(eventExpr, handlerExpr) {
      captured_event <<- eventExpr
      eval(substitute(handlerExpr), envir = parent.frame())
      "observer-registered"
    },
    runLoadSessionObserverEntry = function(...) {
      captured_entry <<- list(...)
      invisible(list(status = "success"))
    },
    logInfo = function(message) {
      captured_log <<- message
      invisible(NULL)
    },
    messageFn = function(message) {
      captured_messages <<- c(captured_messages, message)
      invisible(NULL)
    },
    sysTime = function() start_time
  )

  expect_identical(registration, "observer-registered")
  expect_identical(captured_event, 3)
  expect_identical(captured_log, "=== LOAD FILTERED SESSION BUTTON CLICKED ===")
  expect_identical(captured_messages, "[D66] --- ENTER load_filtered_session observer ---")
  expect_identical(captured_entry$experimentPaths, experiment_paths)
  expect_identical(captured_entry$workflowData, workflow_state)
  expect_identical(captured_entry$daData, da_state)
  expect_identical(captured_entry$session, "metab-session")
  expect_true(is.function(captured_entry$debugLog))
  expect_identical(captured_entry$startTime, start_time)

  captured_entry$debugLog("hello ", "world")
  expect_identical(tail(captured_messages, 1), "[D66] hello world")
})

test_that("metabolomics DA run-analysis observer registration seam preserves event, log, and entry handoff", {
  captured_event <- NULL
  captured_entry <- NULL
  captured_log <- NULL

  workflow_state <- new.env(parent = emptyenv())
  workflow_state$contrasts_tbl <- data.frame(
    friendly_names = c("B vs A"),
    contrasts = c("groupB-groupA"),
    stringsAsFactors = FALSE
  )

  da_state <- new.env(parent = emptyenv())
  da_state$current_s4_object <- "current-s4"

  registration <- registerMetabDaRunAnalysisObserver(
    input = list(
      run_da_analysis = 2,
      formula_string = "~ 0 + group",
      da_q_val_thresh = 0.01,
      treat_lfc_cutoff = 1.5
    ),
    workflowData = workflow_state,
    daData = da_state,
    session = "metab-session",
    experimentPaths = list(
      da_output_dir = "/tmp/metab-da-output",
      publication_graphs_dir = "/tmp/metab-publication-graphs"
    ),
    observeEvent = function(eventExpr, handlerExpr) {
      captured_event <<- eventExpr
      eval(substitute(handlerExpr), envir = parent.frame())
      "observer-registered"
    },
    runAnalysisObserverEntry = function(...) {
      captured_entry <<- list(...)
      invisible(list(status = "success"))
    },
    logInfo = function(message) {
      captured_log <<- message
      invisible(NULL)
    }
  )

  expect_identical(registration, "observer-registered")
  expect_identical(captured_event, 2)
  expect_identical(captured_log, "=== RUN DA ANALYSIS BUTTON CLICKED ===")
  expect_identical(captured_entry$currentS4Object, "current-s4")
  expect_identical(captured_entry$workflowData, workflow_state)
  expect_identical(captured_entry$formulaString, "~ 0 + group")
  expect_identical(captured_entry$daQValThresh, 0.01)
  expect_identical(captured_entry$treatLfcCutoff, 1.5)
  expect_identical(captured_entry$daData, da_state)
  expect_identical(captured_entry$session, "metab-session")
  expect_identical(
    captured_entry$experimentPaths,
    list(
      da_output_dir = "/tmp/metab-da-output",
      publication_graphs_dir = "/tmp/metab-publication-graphs"
    )
  )
})

test_that("metabolomics DA main outputs-and-observers seam preserves server fan-out", {
  input <- list(run_da_analysis = 1, load_filtered_session = 2)
  output <- new.env(parent = emptyenv())
  workflow_state <- new.env(parent = emptyenv())
  da_state <- new.env(parent = emptyenv())
  experiment_paths <- list(source_dir = "/tmp/metab-source")

  captured_load <- NULL
  captured_run <- NULL
  captured_visualization <- NULL
  captured_save <- NULL
  captured_results <- NULL

  registration <- registerMetabDaMainOutputsAndObservers(
    input = input,
    output = output,
    workflowData = workflow_state,
    daData = da_state,
    session = "metab-session",
    experimentPaths = experiment_paths,
    registerLoadFilteredSessionObserver = function(...) {
      captured_load <<- list(...)
      "load-registration"
    },
    registerRunAnalysisObserver = function(...) {
      captured_run <<- list(...)
      "run-registration"
    },
    registerVisualizationOutputs = function(...) {
      captured_visualization <<- list(...)
      "visualization-registration"
    },
    registerSaveHeatmapObserver = function(...) {
      captured_save <<- list(...)
      "save-registration"
    },
    registerResultsOutputs = function(...) {
      captured_results <<- list(...)
      "results-registration"
    }
  )

  expect_identical(
    registration,
    list(
      loadFilteredSessionObserver = "load-registration",
      runAnalysisObserver = "run-registration",
      visualizationOutputs = "visualization-registration",
      saveHeatmapObserver = "save-registration",
      resultsOutputs = "results-registration"
    )
  )
  expect_identical(captured_load$input, input)
  expect_identical(captured_load$workflowData, workflow_state)
  expect_identical(captured_load$daData, da_state)
  expect_identical(captured_load$session, "metab-session")
  expect_identical(captured_load$experimentPaths, experiment_paths)
  expect_identical(captured_run$input, input)
  expect_identical(captured_run$workflowData, workflow_state)
  expect_identical(captured_run$daData, da_state)
  expect_identical(captured_run$session, "metab-session")
  expect_identical(captured_run$experimentPaths, experiment_paths)
  expect_identical(captured_visualization$input, input)
  expect_identical(captured_visualization$output, output)
  expect_identical(captured_visualization$daData, da_state)
  expect_identical(captured_save$input, input)
  expect_identical(captured_save$daData, da_state)
  expect_identical(captured_save$experimentPaths, experiment_paths)
  expect_identical(captured_results$input, input)
  expect_identical(captured_results$output, output)
  expect_identical(captured_results$daData, da_state)
})

test_that("metabolomics DA server-body seam preserves state defaults and registration fan-out", {
  captured_state <- NULL

  server_state <- createMetabDaServerState(
    createReactiveValues = function(...) {
      captured_state <<- list(...)
      captured_state
    }
  )

  expect_identical(server_state, captured_state)
  expect_identical(
    names(captured_state),
    c(
      "da_results_list",
      "contrasts_available",
      "assays_available",
      "analysis_complete",
      "current_s4_object",
      "formula_from_s4",
      "current_row_clusters",
      "current_col_clusters",
      "current_heatmap_plot"
    )
  )
  expect_true(all(vapply(captured_state[c(
    "da_results_list",
    "contrasts_available",
    "assays_available",
    "current_s4_object",
    "formula_from_s4",
    "current_row_clusters",
    "current_col_clusters",
    "current_heatmap_plot"
  )], is.null, logical(1))))
  expect_false(captured_state$analysis_complete)
})

test_that("metabolomics DA server registration seam preserves overview and main fan-out", {
  captured_overview <- NULL
  captured_main <- NULL

  registration <- registerMetabDaServerBodyOutputs(
    input = "metab-input",
    output = "metab-output",
    session = "metab-session",
    workflowData = "workflow-state",
    experimentPaths = "experiment-paths",
    daData = "server-state",
    registerOverviewOutputs = function(...) {
      captured_overview <<- list(...)
      invisible("overview-registration")
    },
    registerMainOutputsAndObservers = function(...) {
      captured_main <<- list(...)
      invisible("main-registration")
    }
  )

  expect_identical(registration, "main-registration")
  expect_identical(captured_overview$output, "metab-output")
  expect_identical(captured_overview$workflowData, "workflow-state")
  expect_identical(captured_overview$daData, "server-state")
  expect_identical(captured_main$input, "metab-input")
  expect_identical(captured_main$output, "metab-output")
  expect_identical(captured_main$workflowData, "workflow-state")
  expect_identical(captured_main$daData, "server-state")
  expect_identical(captured_main$session, "metab-session")
  expect_identical(captured_main$experimentPaths, "experiment-paths")
})

test_that("metabolomics DA server startup seam preserves state-builder and registration handoff", {
  captured_state <- NULL
  captured_registration <- NULL

  server_state <- startMetabDaServerBody(
    input = "metab-input",
    output = "metab-output",
    session = "metab-session",
    workflowData = "workflow-state",
    experimentPaths = "experiment-paths",
    createServerState = function() {
      captured_state <<- list(server = "state")
      captured_state
    },
    registerServerBodyOutputs = function(...) {
      captured_registration <<- list(...)
      "server-registration"
    }
  )

  expect_identical(server_state, "server-registration")
  expect_identical(captured_state, list(server = "state"))
  expect_identical(captured_registration$input, "metab-input")
  expect_identical(captured_registration$output, "metab-output")
  expect_identical(captured_registration$session, "metab-session")
  expect_identical(captured_registration$workflowData, "workflow-state")
  expect_identical(captured_registration$experimentPaths, "experiment-paths")
  expect_identical(captured_registration$daData, captured_state)
})

test_that("metabolomics DA server-body initializer seam preserves startup handoff", {
  captured_startup <- NULL

  server_state <- initializeMetabDaServerBody(
    input = "metab-input",
    output = "metab-output",
    session = "metab-session",
    workflowData = "workflow-state",
    experimentPaths = "experiment-paths",
    createServerState = function() "unused-state-builder",
    registerServerBodyOutputs = function(...) "unused-registration",
    startServerBody = function(...) {
      captured_startup <<- list(...)
      "startup-state"
    }
  )

  expect_identical(server_state, "startup-state")
  expect_identical(captured_startup$input, "metab-input")
  expect_identical(captured_startup$output, "metab-output")
  expect_identical(captured_startup$session, "metab-session")
  expect_identical(captured_startup$workflowData, "workflow-state")
  expect_identical(captured_startup$experimentPaths, "experiment-paths")
  expect_true(is.function(captured_startup$createServerState))
  expect_true(is.function(captured_startup$registerServerBodyOutputs))
})

test_that("metabolomics DA server-module callback seam preserves initializer handoff", {
  captured_body <- NULL

  server_state <- runMetabDaServerModuleCallback(
    input = "metab-input",
    output = "metab-output",
    session = "metab-session",
    workflowData = "workflow-state",
    experimentPaths = "experiment-paths",
    initializeServerBody = function(...) {
      captured_body <<- list(...)
      "server-body-state"
    }
  )

  expect_identical(server_state, "server-body-state")
  expect_identical(captured_body$input, "metab-input")
  expect_identical(captured_body$output, "metab-output")
  expect_identical(captured_body$session, "metab-session")
  expect_identical(captured_body$workflowData, "workflow-state")
  expect_identical(captured_body$experimentPaths, "experiment-paths")
})

test_that("metabolomics DA server-module handler seam preserves callback handoff", {
  captured_callback <- NULL

  module_handler <- createMetabDaServerModuleHandler(
    workflowData = "workflow-state",
    experimentPaths = "experiment-paths",
    runServerCallback = function(...) {
      captured_callback <<- list(...)
      "server-callback-state"
    },
    initializeServerBody = function(...) "server-body-initializer"
  )

  handler_state <- module_handler(
    input = "metab-input",
    output = "metab-output",
    session = "metab-session"
  )

  expect_identical(handler_state, "server-callback-state")
  expect_identical(captured_callback$input, "metab-input")
  expect_identical(captured_callback$output, "metab-output")
  expect_identical(captured_callback$session, "metab-session")
  expect_identical(captured_callback$workflowData, "workflow-state")
  expect_identical(captured_callback$experimentPaths, "experiment-paths")
  expect_true(is.function(captured_callback$initializeServerBody))
})

test_that("metabolomics DA server-module shell seam preserves handler-builder and moduleServer handoff", {
  captured_id <- NULL
  captured_builder <- NULL
  captured_module <- NULL
  module_handler <- function(...) "server-module-handler"

  server_state <- runMetabDaServerModuleShell(
    id = "metab-da",
    workflowData = "workflow-state",
    experimentPaths = "experiment-paths",
    moduleServer = function(id, module) {
      captured_id <<- id
      captured_module <<- module
      "module-server-state"
    },
    runServerCallback = "server-callback-sentinel",
    initializeServerBody = function(...) "server-body-initializer",
    createServerModuleHandler = function(...) {
      captured_builder <<- list(...)
      module_handler
    }
  )

  expect_identical(server_state, "module-server-state")
  expect_identical(captured_id, "metab-da")
  expect_identical(captured_builder$workflowData, "workflow-state")
  expect_identical(captured_builder$experimentPaths, "experiment-paths")
  expect_identical(captured_builder$runServerCallback, "server-callback-sentinel")
  expect_true(is.function(captured_builder$initializeServerBody))
  expect_identical(captured_module, module_handler)
})

test_that("metabolomics DA server-module seam preserves moduleServer shell handoff", {
  captured_shell <- NULL

  server_state <- runMetabDaServerModule(
    id = "metab-da",
    workflowData = "workflow-state",
    experimentPaths = "experiment-paths",
    moduleServer = "module-server-sentinel",
    runServerCallback = "server-callback-sentinel",
    initializeServerBody = function(...) "unused-server-body",
    runServerModuleShell = function(...) {
      captured_shell <<- list(...)
      "server-module-shell-state"
    }
  )

  expect_identical(server_state, "server-module-shell-state")
  expect_identical(captured_shell$id, "metab-da")
  expect_identical(captured_shell$workflowData, "workflow-state")
  expect_identical(captured_shell$experimentPaths, "experiment-paths")
  expect_identical(captured_shell$moduleServer, "module-server-sentinel")
  expect_identical(captured_shell$runServerCallback, "server-callback-sentinel")
  expect_true(is.function(captured_shell$initializeServerBody))
})

test_that("metabolomics DA server-entry seam preserves public wrapper handoff", {
  captured_module <- NULL

  server_state <- runMetabDaServerEntry(
    id = "metab-da",
    workflowData = "workflow-state",
    experimentPaths = "experiment-paths",
    omicType = "metabolomics",
    experimentLabel = "Metabolomics",
    runServerModule = function(...) {
      captured_module <<- list(...)
      "server-entry-state"
    }
  )

  expect_identical(server_state, "server-entry-state")
  expect_identical(captured_module$id, "metab-da")
  expect_identical(captured_module$workflowData, "workflow-state")
  expect_identical(captured_module$experimentPaths, "experiment-paths")
  expect_false("omicType" %in% names(captured_module))
  expect_false("experimentLabel" %in% names(captured_module))
})

test_that("metabolomics DA public server wrapper preserves server-entry handoff", {
  captured_entry <- NULL
  wrapper_env <- environment(mod_metab_da_server)
  original_run_server_entry <- get("runMetabDaServerEntry", envir = wrapper_env)

  on.exit(
    assign(
      "runMetabDaServerEntry",
      original_run_server_entry,
      envir = wrapper_env
    ),
    add = TRUE
  )

  assign(
    "runMetabDaServerEntry",
    function(...) {
      captured_entry <<- list(...)
      "public-server-wrapper-state"
    },
    envir = wrapper_env
  )

  server_state <- mod_metab_da_server(
    id = "metab-da",
    workflow_data = "workflow-state",
    experiment_paths = "experiment-paths",
    omic_type = "metabolomics",
    experiment_label = "Metabolomics"
  )

  expect_identical(server_state, "public-server-wrapper-state")
  expect_identical(captured_entry$id, "metab-da")
  expect_identical(captured_entry$workflowData, "workflow-state")
  expect_identical(captured_entry$experimentPaths, "experiment-paths")
  expect_identical(captured_entry$omicType, "metabolomics")
  expect_identical(captured_entry$experimentLabel, "Metabolomics")
})

test_that("metabolomics DA heatmap save seam preserves payload, prefix, and notification", {
  captured_save <- NULL
  captured_notification <- NULL
  captured_log <- NULL

  file_prefix <- runMetabDaSaveHeatmapObserverShell(
    currentHeatmapPlot = "heatmap-object",
    currentRowClusters = c(MetA = 1, MetB = 2),
    publicationGraphsDir = "/tmp/metab-publication-graphs",
    selectedContrast = "B vs A / LCMS+",
    topN = 30,
    clusteringMethod = "average",
    distanceMethod = "manhattan",
    heatmapClustering = "column",
    scaling = "column",
    colorScheme = "Viridis",
    treeCutMethod = "dynamic",
    nClusters = 5,
    cutHeight = 0.7,
    minClusterSize = 4,
    saveHeatmapProducts = function(...) {
      captured_save <<- list(...)
      invisible("saved")
    },
    showNotification = function(...) {
      captured_notification <<- list(...)
      invisible(NULL)
    },
    logInfo = function(message) {
      captured_log <<- message
      invisible(NULL)
    }
  )

  expect_identical(file_prefix, "metab_B_vs_A___LCMS_")
  expect_identical(captured_log, "Save Heatmap button clicked")
  expect_identical(captured_save$heatmap_obj, "heatmap-object")
  expect_identical(captured_save$row_clusters, c(MetA = 1, MetB = 2))
  expect_identical(captured_save$output_dir, "/tmp/metab-publication-graphs")
  expect_identical(captured_save$file_prefix, "metab_B_vs_A___LCMS_")
  expect_identical(
    captured_save$params_list,
    list(
      contrast = "B vs A / LCMS+",
      top_n = 30,
      clustering_method = "average",
      distance_method = "manhattan",
      cluster_rows = FALSE,
      cluster_cols = TRUE,
      scaling = "column",
      color_scheme = "Viridis",
      tree_cut_method = "dynamic",
      n_clusters = 5,
      cut_height = 0.7,
      min_cluster_size = 4
    )
  )
  expect_identical(
    captured_notification,
    list(
      "Heatmap and cluster info saved to publication_graphs/Heatmap",
      type = "message",
      duration = 5
    )
  )
})

test_that("metabolomics DA heatmap save observer entry preserves req gate", {
  captured_req <- NULL

  expect_error(
    runMetabDaSaveHeatmapObserverEntry(
      currentHeatmapPlot = "heatmap-object",
      publicationGraphsDir = "/tmp/metab-publication-graphs",
      selectedContrast = "B vs A",
      requireInputs = function(...) {
        captured_req <<- list(...)
        stop("req gate tripped")
      },
      saveHeatmapShell = function(...) stop("shell should not run when req fails")
    ),
    "req gate tripped"
  )

  expect_identical(
    captured_req,
    list("heatmap-object", "/tmp/metab-publication-graphs")
  )
})

test_that("metabolomics DA heatmap save observer entry preserves shell handoff", {
  captured_shell <- NULL

  save_state <- runMetabDaSaveHeatmapObserverEntry(
    currentHeatmapPlot = "heatmap-object",
    publicationGraphsDir = "/tmp/metab-publication-graphs",
    currentRowClusters = c(MetA = 1),
    selectedContrast = "B vs A / LCMS+",
    topN = 30,
    clusteringMethod = "average",
    distanceMethod = "manhattan",
    heatmapClustering = "column",
    scaling = "column",
    colorScheme = "Viridis",
    treeCutMethod = "dynamic",
    nClusters = 5,
    cutHeight = 0.7,
    minClusterSize = 4,
    requireInputs = function(...) invisible(list(...)),
    saveHeatmapShell = function(...) {
      captured_shell <<- list(...)
      invisible("saved-prefix")
    }
  )

  expect_identical(save_state$status, "success")
  expect_identical(save_state$saveState, "saved-prefix")
  expect_identical(captured_shell$currentHeatmapPlot, "heatmap-object")
  expect_identical(captured_shell$publicationGraphsDir, "/tmp/metab-publication-graphs")
  expect_identical(captured_shell$currentRowClusters, c(MetA = 1))
  expect_identical(captured_shell$selectedContrast, "B vs A / LCMS+")
  expect_identical(captured_shell$topN, 30)
  expect_identical(captured_shell$clusteringMethod, "average")
  expect_identical(captured_shell$distanceMethod, "manhattan")
  expect_identical(captured_shell$heatmapClustering, "column")
  expect_identical(captured_shell$scaling, "column")
  expect_identical(captured_shell$colorScheme, "Viridis")
  expect_identical(captured_shell$treeCutMethod, "dynamic")
  expect_identical(captured_shell$nClusters, 5)
  expect_identical(captured_shell$cutHeight, 0.7)
  expect_identical(captured_shell$minClusterSize, 4)
})

test_that("metabolomics DA analysis-input seam preserves reactive and state-manager resolution", {
  get_state_called <- FALSE
  current_s4 <- methods::new("MetaboliteAssayData", args = list())
  contrasts_tbl <- data.frame(
    friendly_names = c("B vs A", "C vs A"),
    contrasts = c("groupB-groupA", "groupC-groupA"),
    stringsAsFactors = FALSE
  )

  analysis_inputs <- resolveMetabDaAnalysisInputs(
    currentS4Object = current_s4,
    workflowData = list(contrasts_tbl = contrasts_tbl),
    getState = function() {
      get_state_called <<- TRUE
      stop("state manager should not be used")
    }
  )

  expect_true(isTRUE(analysis_inputs$ok))
  expect_identical(analysis_inputs$currentS4, current_s4)
  expect_identical(analysis_inputs$contrastsTbl, contrasts_tbl)
  expect_null(analysis_inputs$errorMessage)
  expect_false(get_state_called)

  fallback_s4 <- methods::new("MetaboliteAssayData", args = list(source = "state-manager"))
  get_state_called <- FALSE
  analysis_inputs <- resolveMetabDaAnalysisInputs(
    currentS4Object = NULL,
    workflowData = list(contrasts_tbl = contrasts_tbl),
    getState = function() {
      get_state_called <<- TRUE
      fallback_s4
    }
  )

  expect_true(isTRUE(analysis_inputs$ok))
  expect_true(get_state_called)
  expect_identical(analysis_inputs$currentS4, fallback_s4)
  expect_identical(analysis_inputs$contrastsTbl, contrasts_tbl)
})

test_that("metabolomics DA analysis-input seam preserves missing-data and missing-contrast errors", {
  no_data <- resolveMetabDaAnalysisInputs(
    currentS4Object = NULL,
    workflowData = list(contrasts_tbl = data.frame(contrast = "groupB-groupA")),
    getState = function() NULL
  )

  expect_false(isTRUE(no_data$ok))
  expect_null(no_data$contrastsTbl)
  expect_identical(
    no_data$errorMessage,
    "No metabolomics data loaded. Please load a filtered session first."
  )

  wrong_class <- resolveMetabDaAnalysisInputs(
    currentS4Object = structure(list(), class = "NotMetaboliteAssayData"),
    workflowData = list(contrasts_tbl = data.frame(contrast = "groupB-groupA")),
    getState = NULL
  )

  expect_false(isTRUE(wrong_class$ok))
  expect_identical(
    wrong_class$errorMessage,
    "No metabolomics data loaded. Please load a filtered session first."
  )

  missing_contrasts <- resolveMetabDaAnalysisInputs(
    currentS4Object = methods::new("MetaboliteAssayData", args = list()),
    workflowData = list(contrasts_tbl = data.frame(
      friendly_names = character(),
      contrasts = character(),
      stringsAsFactors = FALSE
    )),
    getState = NULL
  )

  expect_false(isTRUE(missing_contrasts$ok))
  expect_s3_class(missing_contrasts$contrastsTbl, "data.frame")
  expect_equal(nrow(missing_contrasts$contrastsTbl), 0)
  expect_identical(
    missing_contrasts$errorMessage,
    "No contrasts defined. Please define contrasts in the design tab."
  )
})

test_that("metabolomics DA run-analysis observer entry preserves early validation notification", {
  captured_notifications <- list()

  analysis_state <- runMetabDaAnalysisObserverEntry(
    currentS4Object = NULL,
    workflowData = list(contrasts_tbl = data.frame()),
    formulaString = "~ 0 + group",
    daQValThresh = 0.05,
    treatLfcCutoff = 0,
    daData = "da-state",
    session = "metab-session",
    experimentPaths = list(),
    resolveAnalysisInputs = function(currentS4Object, workflowData, getState) {
      expect_null(currentS4Object)
      expect_true(is.list(workflowData))
      expect_null(getState)
      list(
        ok = FALSE,
        currentS4 = NULL,
        contrastsTbl = NULL,
        errorMessage = "No metabolomics data loaded. Please load a filtered session first."
      )
    },
    runAnalysisShell = function(...) stop("shell should not run when validation fails"),
    showNotification = function(...) {
      captured_notifications[[length(captured_notifications) + 1]] <<- list(...)
      invisible(NULL)
    }
  )

  expect_identical(analysis_state$status, "error")
  expect_identical(analysis_state$stage, "resolve_analysis_inputs")
  expect_identical(
    analysis_state$errorMessage,
    "No metabolomics data loaded. Please load a filtered session first."
  )
  expect_identical(
    captured_notifications,
    list(
      list(
        "No metabolomics data loaded. Please load a filtered session first.",
        type = "error",
        duration = 5
      )
    )
  )
})

test_that("metabolomics DA run-analysis observer entry preserves validation-to-shell handoff", {
  captured_shell <- NULL

  workflow_state <- new.env(parent = emptyenv())
  workflow_state$contrasts_tbl <- data.frame(
    friendly_names = c("B vs A"),
    contrasts = c("groupB-groupA"),
    stringsAsFactors = FALSE
  )

  analysis_state <- runMetabDaAnalysisObserverEntry(
    currentS4Object = "current-s4",
    workflowData = workflow_state,
    formulaString = "~ 0 + group",
    daQValThresh = 0.01,
    treatLfcCutoff = 1.5,
    daData = "da-state",
    session = "metab-session",
    experimentPaths = list(da_output_dir = "/tmp/metab-da-output"),
    resolveAnalysisInputs = function(currentS4Object, workflowData, getState) {
      expect_identical(currentS4Object, "current-s4")
      expect_identical(workflowData, workflow_state)
      expect_null(getState)
      list(
        ok = TRUE,
        currentS4 = "resolved-s4",
        contrastsTbl = workflow_state$contrasts_tbl,
        errorMessage = NULL
      )
    },
    runAnalysisShell = function(...) {
      captured_shell <<- list(...)
      invisible(list(status = "success", results = "analysis-results"))
    },
    showNotification = function(...) stop("should not notify on validation success")
  )

  expect_identical(analysis_state$status, "success")
  expect_identical(analysis_state$analysisInputs$currentS4, "resolved-s4")
  expect_identical(
    analysis_state$analysisInputs$contrastsTbl,
    workflow_state$contrasts_tbl
  )
  expect_identical(
    analysis_state$shellState,
    list(status = "success", results = "analysis-results")
  )
  expect_identical(captured_shell$currentS4, "resolved-s4")
  expect_identical(captured_shell$contrastsTbl, workflow_state$contrasts_tbl)
  expect_identical(captured_shell$formulaString, "~ 0 + group")
  expect_identical(captured_shell$daQValThresh, 0.01)
  expect_identical(captured_shell$treatLfcCutoff, 1.5)
  expect_identical(captured_shell$workflowData, workflow_state)
  expect_identical(captured_shell$daData, "da-state")
  expect_identical(captured_shell$session, "metab-session")
  expect_identical(
    captured_shell$experimentPaths$da_output_dir,
    "/tmp/metab-da-output"
  )
})

test_that("metabolomics DA run-analysis shell preserves success-path state, notifications, and follow-up helpers", {
  captured_notifications <- list()
  captured_removed <- character()
  captured_run <- NULL
  captured_selectors <- NULL
  captured_artifacts <- NULL
  captured_info_logs <- character()
  captured_error_logs <- character()

  workflow_state <- new.env(parent = emptyenv())
  workflow_state$tab_status <- list(differential_analysis = "pending", other = "keep")

  da_state <- new.env(parent = emptyenv())
  current_s4 <- methods::new("MetaboliteAssayData", args = list())
  contrasts_tbl <- data.frame(
    friendly_names = c("B vs A"),
    contrasts = c("groupB-groupA"),
    stringsAsFactors = FALSE
  )
  results <- list(
    da_metabolites_long = buildMetabDaDisplayResults(),
    significant_counts = list(LCMS_Pos = list(up = 1, down = 2, ns = 3))
  )

  analysis_state <- runMetabDaAnalysisObserverShell(
    currentS4 = current_s4,
    contrastsTbl = contrasts_tbl,
    formulaString = "~ 0 + group",
    daQValThresh = 0.01,
    treatLfcCutoff = 1.5,
    workflowData = workflow_state,
    daData = da_state,
    session = "metab-session",
    experimentPaths = list(
      da_output_dir = "/tmp/metab-da-output",
      publication_graphs_dir = "/tmp/metab-publication-graphs"
    ),
    runAnalysis = function(...) {
      captured_run <<- list(...)
      results
    },
    updateSelectors = function(...) {
      captured_selectors <<- list(...)
      invisible(list(updated = TRUE))
    },
    writeArtifacts = function(...) {
      captured_artifacts <<- list(...)
      invisible(list(status = "written"))
    },
    showNotification = function(...) {
      captured_notifications[[length(captured_notifications) + 1]] <<- list(...)
      invisible(NULL)
    },
    removeNotification = function(id) {
      captured_removed <<- c(captured_removed, id)
      invisible(NULL)
    },
    logInfo = function(message) {
      captured_info_logs <<- c(captured_info_logs, message)
      invisible(NULL)
    },
    logError = function(message) {
      captured_error_logs <<- c(captured_error_logs, message)
      invisible(NULL)
    }
  )

  expect_identical(analysis_state$status, "success")
  expect_identical(analysis_state$results, results)
  expect_identical(analysis_state$selectorState, list(updated = TRUE))
  expect_identical(analysis_state$diskState, list(status = "written"))
  expect_identical(captured_run$theObject, current_s4)
  expect_identical(captured_run$contrasts_tbl, contrasts_tbl)
  expect_identical(captured_run$formula_string, "~ 0 + group")
  expect_identical(captured_run$da_q_val_thresh, 0.01)
  expect_identical(captured_run$treat_lfc_cutoff, 1.5)
  expect_identical(captured_run$eBayes_trend, TRUE)
  expect_identical(captured_run$eBayes_robust, TRUE)
  expect_identical(da_state$da_results_list, results)
  expect_true(isTRUE(da_state$analysis_complete))
  expect_identical(
    workflow_state$tab_status,
    list(differential_analysis = "complete", other = "keep")
  )
  expect_identical(captured_removed, "da_running")
  expect_identical(captured_selectors$daResultsLong, results$da_metabolites_long)
  expect_identical(captured_selectors$session, "metab-session")
  expect_identical(captured_artifacts$results, results)
  expect_identical(captured_artifacts$experimentPaths$da_output_dir, "/tmp/metab-da-output")
  expect_identical(captured_artifacts$daQValThresh, 0.01)
  expect_identical(captured_artifacts$treatLfcCutoff, 1.5)
  expect_identical(
    captured_notifications,
    list(
      list(
        "Running differential expression analysis...",
        id = "da_running",
        duration = NULL
      ),
      list(
        "Differential expression analysis complete!",
        type = "message",
        duration = 5
      )
    )
  )
  expect_identical(captured_info_logs, "   DA analysis completed successfully")
  expect_length(captured_error_logs, 0)
})

test_that("metabolomics DA run-analysis shell preserves error-path notification and leaves state untouched", {
  captured_notifications <- list()
  captured_removed <- character()
  captured_error_logs <- character()

  workflow_state <- new.env(parent = emptyenv())
  workflow_state$tab_status <- list(differential_analysis = "pending")

  da_state <- new.env(parent = emptyenv())
  da_state$analysis_complete <- FALSE
  current_s4 <- methods::new("MetaboliteAssayData", args = list())
  contrasts_tbl <- data.frame(
    friendly_names = c("B vs A"),
    contrasts = c("groupB-groupA"),
    stringsAsFactors = FALSE
  )

  analysis_state <- runMetabDaAnalysisObserverShell(
    currentS4 = current_s4,
    contrastsTbl = contrasts_tbl,
    formulaString = "~ 0 + group",
    daQValThresh = 0.05,
    treatLfcCutoff = 0,
    workflowData = workflow_state,
    daData = da_state,
    session = "metab-session",
    experimentPaths = list(),
    runAnalysis = function(...) stop("model matrix singular"),
    updateSelectors = function(...) stop("should not update selectors"),
    writeArtifacts = function(...) stop("should not write artifacts"),
    showNotification = function(...) {
      captured_notifications[[length(captured_notifications) + 1]] <<- list(...)
      invisible(NULL)
    },
    removeNotification = function(id) {
      captured_removed <<- c(captured_removed, id)
      invisible(NULL)
    },
    logInfo = function(...) invisible(NULL),
    logError = function(message) {
      captured_error_logs <<- c(captured_error_logs, message)
      invisible(NULL)
    }
  )

  expect_identical(analysis_state$status, "error")
  expect_identical(analysis_state$errorMessage, "model matrix singular")
  expect_identical(workflow_state$tab_status, list(differential_analysis = "pending"))
  expect_identical(da_state$analysis_complete, FALSE)
  expect_false(exists("da_results_list", envir = da_state, inherits = FALSE))
  expect_identical(captured_removed, "da_running")
  expect_identical(
    captured_notifications,
    list(
      list(
        "Running differential expression analysis...",
        id = "da_running",
        duration = NULL
      ),
      list(
        "Analysis error: model matrix singular",
        type = "error",
        duration = 10
      )
    )
  )
  expect_identical(captured_error_logs, "   DA analysis error: model matrix singular")
})

test_that("metabolomics DA load-session path seam preserves fallback resolution and error branches", {
  captured_logs <- character()

  resolved_session <- resolveMetabDaLoadSessionFile(
    experimentPaths = list(
      source_dir = "/tmp/missing",
      export_dir = "/tmp/export"
    ),
    dirExists = function(path) identical(path, "/tmp/export"),
    fileExists = function(path) identical(
      path,
      "/tmp/export/metab_filtered_session_data_latest.rds"
    ),
    filePath = function(...) paste(..., sep = "/"),
    debugLog = function(...) {
      captured_logs <<- c(captured_logs, paste0(...))
      invisible(NULL)
    }
  )

  expect_true(resolved_session$ok)
  expect_identical(resolved_session$sourceDir, "/tmp/export")
  expect_identical(
    resolved_session$sessionFile,
    "/tmp/export/metab_filtered_session_data_latest.rds"
  )
  expect_null(resolved_session$errorMessage)
  expect_true(any(grepl("trying export_dir", captured_logs, fixed = TRUE)))

  missing_dir <- resolveMetabDaLoadSessionFile(
    experimentPaths = list(
      source_dir = NULL,
      export_dir = "/tmp/missing"
    ),
    dirExists = function(...) FALSE,
    fileExists = function(...) stop("should not check files when directories are missing"),
    filePath = function(...) paste(..., sep = "/")
  )

  expect_false(missing_dir$ok)
  expect_null(missing_dir$sessionFile)
  expect_identical(
    missing_dir$errorMessage,
    "Could not find source directory for session data."
  )

  missing_file <- resolveMetabDaLoadSessionFile(
    experimentPaths = list(
      source_dir = "/tmp/source",
      export_dir = "/tmp/export"
    ),
    dirExists = function(path) identical(path, "/tmp/source"),
    fileExists = function(...) FALSE,
    filePath = function(...) paste(..., sep = "/")
  )

  expect_false(missing_file$ok)
  expect_identical(missing_file$sourceDir, "/tmp/source")
  expect_identical(
    missing_file$sessionFile,
    "/tmp/source/metab_filtered_session_data_latest.rds"
  )
  expect_identical(
    missing_file$errorMessage,
    "Session file not found: /tmp/source/metab_filtered_session_data_latest.rds"
  )
})

test_that("metabolomics DA load-session observer entry preserves early error notifications", {
  captured_notifications <- list()

  load_state <- runMetabDaLoadSessionObserverEntry(
    experimentPaths = "experiment-paths",
    workflowData = "workflow-state",
    daData = "da-state",
    session = "metab-session",
    resolveSessionFile = function(experimentPaths, debugLog) {
      expect_identical(experimentPaths, "experiment-paths")
      expect_true(is.function(debugLog))
      list(
        ok = FALSE,
        sessionFile = NULL,
        errorMessage = "Could not find source directory for session data."
      )
    },
    loadSessionShell = function(...) stop("should not load when resolution fails"),
    showNotification = function(...) {
      captured_notifications[[length(captured_notifications) + 1]] <<- list(...)
      invisible(NULL)
    }
  )

  expect_identical(load_state$status, "error")
  expect_identical(load_state$stage, "resolve_session_file")
  expect_null(load_state$sessionResolution$sessionFile)
  expect_identical(
    load_state$errorMessage,
    "Could not find source directory for session data."
  )
  expect_identical(
    captured_notifications,
    list(
      list(
        "Could not find source directory for session data.",
        type = "error",
        duration = 5
      )
    )
  )
})

test_that("metabolomics DA load-session observer entry preserves resolution-to-shell handoff", {
  captured_shell <- NULL

  load_state <- runMetabDaLoadSessionObserverEntry(
    experimentPaths = "experiment-paths",
    workflowData = "workflow-state",
    daData = "da-state",
    session = "metab-session",
    resolveSessionFile = function(experimentPaths, debugLog) {
      expect_identical(experimentPaths, "experiment-paths")
      expect_true(is.function(debugLog))
      list(
        ok = TRUE,
        sessionFile = "/tmp/metab_filtered_session_data_latest.rds",
        errorMessage = NULL
      )
    },
    loadSessionShell = function(sessionFile, workflowData, daData, session, debugLog, startTime) {
      captured_shell <<- list(
        sessionFile = sessionFile,
        workflowData = workflowData,
        daData = daData,
        session = session,
        debugLog = debugLog,
        startTime = startTime
      )
      invisible(list(status = "success", restored = TRUE))
    },
    showNotification = function(...) stop("should not notify on success"),
    startTime = "start-marker"
  )

  expect_identical(load_state$status, "success")
  expect_identical(load_state$sessionResolution$sessionFile, "/tmp/metab_filtered_session_data_latest.rds")
  expect_identical(load_state$loadState, list(status = "success", restored = TRUE))
  expect_identical(captured_shell$sessionFile, "/tmp/metab_filtered_session_data_latest.rds")
  expect_identical(captured_shell$workflowData, "workflow-state")
  expect_identical(captured_shell$daData, "da-state")
  expect_identical(captured_shell$session, "metab-session")
  expect_true(is.function(captured_shell$debugLog))
  expect_identical(captured_shell$startTime, "start-marker")
})

test_that("metabolomics DA load-session observer shell preserves success notifications and restore handoff", {
  captured_notifications <- list()
  captured_removed <- character()
  captured_info_logs <- character()
  captured_error_logs <- character()
  captured_logs <- character()
  captured_restore <- NULL
  session_data <- list(
    current_s4_object = methods::new("MetaboliteAssayData", args = list()),
    contrasts_tbl = data.frame(),
    assay_names = character()
  )

  load_state <- runMetabDaLoadSessionObserverShell(
    sessionFile = "/tmp/metab_filtered_session_data_latest.rds",
    workflowData = "workflow-state",
    daData = "da-state",
    session = "metab-session",
    readSessionData = function(path) {
      expect_identical(path, "/tmp/metab_filtered_session_data_latest.rds")
      session_data
    },
    restoreState = function(sessionData, sessionFile, workflowData, daData, session, debugLog) {
      captured_restore <<- list(
        sessionData = sessionData,
        sessionFile = sessionFile,
        workflowData = workflowData,
        daData = daData,
        session = session,
        debugLog = debugLog
      )
      invisible(list(restored = TRUE))
    },
    showNotification = function(...) {
      captured_notifications[[length(captured_notifications) + 1]] <<- list(...)
      invisible(NULL)
    },
    removeNotification = function(id) {
      captured_removed <<- c(captured_removed, id)
      invisible(NULL)
    },
    logInfo = function(message) {
      captured_info_logs <<- c(captured_info_logs, message)
      invisible(NULL)
    },
    logError = function(message) {
      captured_error_logs <<- c(captured_error_logs, message)
      invisible(NULL)
    },
    debugLog = function(...) {
      captured_logs <<- c(captured_logs, paste0(...))
      invisible(NULL)
    },
    startTime = as.POSIXct("2026-04-16 12:00:00", tz = "UTC"),
    timeNow = function() as.POSIXct("2026-04-16 12:00:02", tz = "UTC")
  )

  expect_identical(load_state$status, "success")
  expect_identical(load_state$sessionData, session_data)
  expect_identical(load_state$restoreState, list(restored = TRUE))
  expect_identical(captured_restore$sessionData, session_data)
  expect_identical(captured_restore$sessionFile, "/tmp/metab_filtered_session_data_latest.rds")
  expect_identical(captured_restore$workflowData, "workflow-state")
  expect_identical(captured_restore$daData, "da-state")
  expect_identical(captured_restore$session, "metab-session")
  expect_true(is.function(captured_restore$debugLog))
  expect_identical(captured_removed, "loading_session")
  expect_identical(
    captured_notifications,
    list(
      list(
        "Loading filtered session...",
        id = "loading_session",
        duration = NULL
      ),
      list(
        "Session loaded successfully!",
        type = "message",
        duration = 3
      )
    )
  )
  expect_identical(
    captured_info_logs,
    "   Loaded session from: /tmp/metab_filtered_session_data_latest.rds"
  )
  expect_length(captured_error_logs, 0)
  expect_true(any(grepl("STEP 5: Reading RDS file...", captured_logs, fixed = TRUE)))
  expect_true(any(grepl("names(session_data): current_s4_object, contrasts_tbl, assay_names", captured_logs, fixed = TRUE)))
  expect_true(any(grepl("STEP 10: Removing notification and showing success...", captured_logs, fixed = TRUE)))
  expect_true(any(grepl("EXIT load_filtered_session (2s) --- SUCCESS", captured_logs, fixed = TRUE)))
})

test_that("metabolomics DA load-session observer shell preserves error notifications and leaves restore untouched", {
  captured_notifications <- list()
  captured_removed <- character()
  captured_error_logs <- character()
  captured_logs <- character()

  load_state <- runMetabDaLoadSessionObserverShell(
    sessionFile = "/tmp/metab_filtered_session_data_latest.rds",
    workflowData = "workflow-state",
    daData = "da-state",
    session = "metab-session",
    readSessionData = function(...) stop("corrupt session"),
    restoreState = function(...) stop("should not restore when read fails"),
    showNotification = function(...) {
      captured_notifications[[length(captured_notifications) + 1]] <<- list(...)
      invisible(NULL)
    },
    removeNotification = function(id) {
      captured_removed <<- c(captured_removed, id)
      invisible(NULL)
    },
    logInfo = function(...) invisible(NULL),
    logError = function(message) {
      captured_error_logs <<- c(captured_error_logs, message)
      invisible(NULL)
    },
    debugLog = function(...) {
      captured_logs <<- c(captured_logs, paste0(...))
      invisible(NULL)
    }
  )

  expect_identical(load_state$status, "error")
  expect_identical(load_state$errorMessage, "corrupt session")
  expect_identical(captured_removed, "loading_session")
  expect_identical(
    captured_notifications,
    list(
      list(
        "Loading filtered session...",
        id = "loading_session",
        duration = NULL
      ),
      list(
        "Error loading session: corrupt session",
        type = "error",
        duration = 10
      )
    )
  )
  expect_identical(
    captured_error_logs,
    "   Error loading session: corrupt session"
  )
  expect_true(any(grepl("FATAL ERROR: corrupt session", captured_logs, fixed = TRUE)))
  expect_true(any(grepl("EXIT load_filtered_session --- FAILED", captured_logs, fixed = TRUE)))
})

test_that("metabolomics DA load-session restore seam preserves state, selectors, and formula sync", {
  captured_save_state <- NULL
  captured_updates <- list()
  captured_formula_update <- NULL
  captured_info_logs <- character()
  captured_warn_logs <- character()

  workflow_state <- new.env(parent = emptyenv())
  workflow_state$state_manager <- list(
    saveState = function(...) {
      captured_save_state <<- list(...)
      invisible(NULL)
    }
  )

  da_state <- new.env(parent = emptyenv())
  s4_object <- methods::new(
    "MetaboliteAssayData",
    args = list(
      daAnalysisParameters = list(formula_string = "~ 0 + group")
    )
  )
  session_data <- list(
    current_s4_object = s4_object,
    r6_current_state_name = "filtered_ready",
    contrasts_tbl = data.frame(
      friendly_names = c("B vs A", "C vs A"),
      contrasts = c("groupB-groupA", "groupC-groupA"),
      stringsAsFactors = FALSE
    ),
    assay_names = c("LCMS_Pos", "LCMS_Neg")
  )

  restore_state <- restoreMetabDaLoadedSessionState(
    sessionData = session_data,
    sessionFile = "/tmp/metab_filtered_session_data_latest.rds",
    workflowData = workflow_state,
    daData = da_state,
    session = "metab-session",
    updateSelectInput = function(session, inputId, choices, selected = NULL) {
      captured_updates[[length(captured_updates) + 1]] <<- list(
        session = session,
        inputId = inputId,
        choices = choices,
        selected = selected
      )
      invisible(NULL)
    },
    updateTextAreaInput = function(session, inputId, value) {
      captured_formula_update <<- list(
        session = session,
        inputId = inputId,
        value = value
      )
      invisible(NULL)
    },
    logInfo = function(message) {
      captured_info_logs <<- c(captured_info_logs, message)
      invisible(NULL)
    },
    logWarn = function(message) {
      captured_warn_logs <<- c(captured_warn_logs, message)
      invisible(NULL)
    }
  )

  expect_identical(restore_state$stateName, "filtered_ready")
  expect_identical(restore_state$contrastChoices, c("B vs A", "C vs A"))
  expect_identical(restore_state$assayChoices, c("Combined", "LCMS_Pos", "LCMS_Neg"))
  expect_identical(restore_state$formulaValue, "~ 0 + group")
  expect_identical(da_state$current_s4_object, s4_object)
  expect_identical(da_state$contrasts_available, session_data$contrasts_tbl)
  expect_identical(da_state$assays_available, c("LCMS_Pos", "LCMS_Neg"))
  expect_identical(da_state$formula_from_s4, "~ 0 + group")
  expect_identical(workflow_state$contrasts_tbl, session_data$contrasts_tbl)
  expect_identical(captured_save_state$state_name, "filtered_ready")
  expect_identical(captured_save_state$s4_data_object, s4_object)
  expect_identical(
    captured_save_state$config_object,
    list(loaded_from = "/tmp/metab_filtered_session_data_latest.rds")
  )
  expect_identical(
    captured_save_state$description,
    "Loaded from filtered session for DA analysis"
  )
  expect_length(captured_updates, 6)
  expect_identical(
    vapply(captured_updates, `[[`, character(1), "inputId"),
    c(
      "volcano_contrast",
      "heatmap_contrast",
      "table_contrast",
      "volcano_assay",
      "heatmap_assay",
      "table_assay"
    )
  )
  expect_identical(captured_updates[[1]]$selected, "B vs A")
  expect_null(captured_updates[[4]]$selected)
  expect_null(captured_updates[[6]]$selected)
  expect_identical(
    captured_formula_update,
    list(
      session = "metab-session",
      inputId = "formula_string",
      value = "~ 0 + group"
    )
  )
  expect_identical(
    captured_info_logs,
    c(
      "   Restored S4 object to state manager",
      "   Restored 2 contrasts"
    )
  )
  expect_length(captured_warn_logs, 0)
})

test_that("metabolomics DA load-session restore seam preserves fallback state and formula-warning branch", {
  captured_save_state <- NULL
  captured_warn_logs <- character()

  workflow_state <- new.env(parent = emptyenv())
  workflow_state$state_manager <- list(
    saveState = function(...) {
      captured_save_state <<- list(...)
      invisible(NULL)
    }
  )

  da_state <- new.env(parent = emptyenv())
  session_data <- list(
    current_s4_object = list(kind = "non-s4"),
    r6_current_state_name = NULL,
    contrasts_tbl = NULL,
    assay_names = NULL
  )

  restore_state <- restoreMetabDaLoadedSessionState(
    sessionData = session_data,
    sessionFile = "/tmp/metab_filtered_session_data_latest.rds",
    workflowData = workflow_state,
    daData = da_state,
    session = "metab-session",
    updateSelectInput = function(...) stop("should not update selectors"),
    updateTextAreaInput = function(...) stop("should not update formula"),
    logInfo = function(...) invisible(NULL),
    logWarn = function(message) {
      captured_warn_logs <<- c(captured_warn_logs, message)
      invisible(NULL)
    }
  )

  expect_identical(restore_state$stateName, "loaded_for_de")
  expect_identical(restore_state$contrastChoices, character())
  expect_identical(restore_state$assayChoices, character())
  expect_null(restore_state$formulaValue)
  expect_identical(da_state$current_s4_object, session_data$current_s4_object)
  expect_identical(captured_save_state$state_name, "loaded_for_de")
  expect_identical(captured_save_state$s4_data_object, session_data$current_s4_object)
  expect_length(captured_warn_logs, 1)
  expect_match(captured_warn_logs[[1]], "^Could not extract formula from S4:")
})

test_that("metabolomics DA selector-update seam preserves choices, selections, and fallback", {
  captured_updates <- list()
  captured_log <- NULL

  selector_state <- updateMetabDaResultsSelectorInputs(
    daResultsLong = buildMetabDaDisplayResults(),
    session = "metab-session",
    updateSelectInput = function(session, inputId, choices, selected = NULL) {
      captured_updates[[length(captured_updates) + 1]] <<- list(
        session = session,
        inputId = inputId,
        choices = choices,
        selected = selected
      )
      invisible(NULL)
    },
    logInfo = function(message) {
      captured_log <<- message
      invisible(NULL)
    }
  )

  expect_length(captured_updates, 6)
  expect_identical(
    vapply(captured_updates, `[[`, character(1), "inputId"),
    c(
      "volcano_contrast",
      "heatmap_contrast",
      "table_contrast",
      "volcano_assay",
      "heatmap_assay",
      "table_assay"
    )
  )
  expect_identical(selector_state$contrastChoices, c("B vs A", "C vs A"))
  expect_identical(selector_state$assayChoices, c("Combined", "LCMS_Pos", "LCMS_Neg"))
  expect_identical(selector_state$tableAssayChoices, c("All", "LCMS_Pos", "LCMS_Neg"))
  expect_identical(captured_updates[[1]]$selected, "B vs A")
  expect_identical(captured_updates[[4]]$selected, "Combined")
  expect_identical(captured_updates[[6]]$selected, "All")
  expect_identical(captured_log, "   Updated dropdowns: 2 contrasts, 2 assays")

  captured_updates <- list()
  fallback_results <- buildMetabDaDisplayResults()[setdiff(
    names(buildMetabDaDisplayResults()),
    "friendly_name"
  )]
  selector_state <- updateMetabDaResultsSelectorInputs(
    daResultsLong = fallback_results,
    session = "metab-session",
    updateSelectInput = function(session, inputId, choices, selected = NULL) {
      captured_updates[[length(captured_updates) + 1]] <<- list(
        session = session,
        inputId = inputId,
        choices = choices,
        selected = selected
      )
      invisible(NULL)
    },
    logInfo = function(...) invisible(NULL)
  )

  expect_identical(
    selector_state$contrastChoices,
    c("groupB-groupA", "groupC-groupA")
  )
  expect_identical(captured_updates[[1]]$choices, c("groupB-groupA", "groupC-groupA"))
  expect_null(updateMetabDaResultsSelectorInputs(NULL, session = "metab-session"))
})

test_that("metabolomics DA results-disk seam preserves success, skip, and warning branches", {
  captured_output <- NULL
  captured_notifications <- list()
  captured_logs <- character()
  results <- list(da_metabolites_long = buildMetabDaDisplayResults())

  disk_state <- writeMetabDaResultsArtifacts(
    results = results,
    experimentPaths = list(
      da_output_dir = "/tmp/metab-da-output",
      publication_graphs_dir = "/tmp/metab-publication-graphs"
    ),
    daQValThresh = 0.01,
    treatLfcCutoff = 1.5,
    outputResults = function(...) {
      captured_output <<- list(...)
      TRUE
    },
    showNotification = function(...) {
      captured_notifications[[length(captured_notifications) + 1]] <<- list(...)
      invisible(NULL)
    },
    logInfo = function(message) {
      captured_logs <<- c(captured_logs, message)
      invisible(NULL)
    },
    logWarn = function(message) {
      captured_logs <<- c(captured_logs, message)
      invisible(NULL)
    }
  )

  expect_identical(disk_state$status, "written")
  expect_true(isTRUE(disk_state$success))
  expect_identical(captured_output$da_results_list, results)
  expect_identical(captured_output$da_output_dir, "/tmp/metab-da-output")
  expect_identical(captured_output$publication_graphs_dir, "/tmp/metab-publication-graphs")
  expect_identical(captured_output$da_q_val_thresh, 0.01)
  expect_identical(captured_output$lfc_threshold, 1.5)
  expect_identical(captured_output$heatmap_top_n, 50)
  expect_identical(captured_output$heatmap_clustering, "both")
  expect_identical(captured_output$heatmap_color_scheme, "RdBu")
  expect_identical(
    captured_notifications[[1]],
    list(
      "DA results saved to disk (tables, volcano plots, heatmaps)",
      type = "message",
      duration = 5
    )
  )
  expect_identical(
    captured_logs,
    c(
      "   Writing DA results to disk...",
      "   da_output_dir = /tmp/metab-da-output",
      "   publication_graphs_dir = /tmp/metab-publication-graphs",
      "   All DA results written to disk successfully"
    )
  )

  captured_notifications <- list()
  captured_logs <- character()
  skip_state <- writeMetabDaResultsArtifacts(
    results = results,
    experimentPaths = list(
      da_output_dir = NULL,
      publication_graphs_dir = "/tmp/metab-publication-graphs"
    ),
    outputResults = function(...) stop("should not be called"),
    showNotification = function(...) {
      captured_notifications[[length(captured_notifications) + 1]] <<- list(...)
      invisible(NULL)
    },
    logInfo = function(message) {
      captured_logs <<- c(captured_logs, message)
      invisible(NULL)
    },
    logWarn = function(message) {
      captured_logs <<- c(captured_logs, message)
      invisible(NULL)
    }
  )

  expect_identical(skip_state$status, "skipped")
  expect_identical(skip_state$daOutputDir, NULL)
  expect_identical(skip_state$publicationGraphsDir, "/tmp/metab-publication-graphs")
  expect_length(captured_notifications, 0)
  expect_identical(
    captured_logs,
    c(
      "   Writing DA results to disk...",
      "   da_output_dir = NULL",
      "   publication_graphs_dir = /tmp/metab-publication-graphs",
      "   Output directories not configured, skipping file output"
    )
  )

  captured_notifications <- list()
  captured_logs <- character()
  error_state <- writeMetabDaResultsArtifacts(
    results = results,
    experimentPaths = list(
      da_output_dir = "/tmp/metab-da-output",
      publication_graphs_dir = "/tmp/metab-publication-graphs"
    ),
    outputResults = function(...) stop("disk full"),
    showNotification = function(...) {
      captured_notifications[[length(captured_notifications) + 1]] <<- list(...)
      invisible(NULL)
    },
    logInfo = function(message) {
      captured_logs <<- c(captured_logs, message)
      invisible(NULL)
    },
    logWarn = function(message) {
      captured_logs <<- c(captured_logs, message)
      invisible(NULL)
    }
  )

  expect_identical(error_state$status, "error")
  expect_identical(error_state$errorMessage, "disk full")
  expect_identical(
    captured_notifications[[1]],
    list(
      "Warning: Could not save results to disk: disk full",
      type = "warning",
      duration = 8
    )
  )
  expect_identical(
    tail(captured_logs, 1),
    "   Could not write DA results to disk: disk full"
  )
})

test_that("metabolomics DA download seam preserves req gate and CSV payload", {
  captured_req <- NULL
  captured_write <- NULL
  results <- buildMetabDaDisplayResults()

  output_file <- writeMetabDaResultsDownloadCsv(
    results = results,
    file = "metab-da-results.csv",
    req = function(value) {
      captured_req <<- value
      value
    },
    writeCsv = function(x, file, row.names) {
      captured_write <<- list(x = x, file = file, row.names = row.names)
      invisible(NULL)
    }
  )

  expect_identical(output_file, "metab-da-results.csv")
  expect_identical(captured_req, results)
  expect_identical(captured_write$x, results)
  expect_identical(captured_write$file, "metab-da-results.csv")
  expect_identical(captured_write$row.names, FALSE)
})

test_that("metabolomics DA download-handler seam preserves filename and CSV handoff", {
  captured_download_handler <- NULL
  captured_write <- NULL
  results <- buildMetabDaDisplayResults()

  handler <- buildMetabDaResultsDownloadHandler(
    daData = list(
      da_results_list = list(
        da_metabolites_long = results
      )
    ),
    downloadHandler = function(filename, content) {
      captured_download_handler <<- list(
        filename = filename,
        content = content
      )
      list(filename = filename, content = content)
    },
    buildFilename = function() {
      "metabolomics-da-results.csv"
    },
    writeDownloadCsv = function(results, file) {
      captured_write <<- list(results = results, file = file)
      invisible(file)
    }
  )

  expect_true(is.list(handler))
  expect_identical(handler$filename(), "metabolomics-da-results.csv")
  expect_identical(
    captured_download_handler$filename(),
    "metabolomics-da-results.csv"
  )

  handler$content("download-target.csv")

  expect_identical(captured_write$results, results)
  expect_identical(captured_write$file, "download-target.csv")
})

test_that("metabolomics DA cluster summary seam preserves empty and truncation exits", {
  expect_identical(
    buildMetabDaClusterSummaryText(NULL),
    "No clusters defined. Enable clustering and tree cutting on the heatmap."
  )

  cluster_names <- paste0("Met", seq_len(22))
  clusters <- rep(3, 22)
  names(clusters) <- cluster_names
  summary_text <- buildMetabDaClusterSummaryText(clusters)

  expect_match(summary_text, "Cluster 3 \\(22 metabolites\\):")
  expect_match(summary_text, "Met20")
  expect_match(summary_text, ", ... and 2 more", fixed = TRUE)
  expect_false(grepl("Met21|Met22", summary_text))
})
