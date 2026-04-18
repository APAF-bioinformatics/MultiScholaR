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
  paths = c(
    file.path(repo_root, "R", "mod_metab_norm_server_helpers.R"),
    file.path(repo_root, "R", "mod_metab_norm.R")
  ),
  symbols = c(
    "updateMetabNormDesignDrivenChoices",
    "initializeMetabNormAssayNames",
    "runMetabNormAutoPreNormalizationQcObserverShell",
    "renderMetabNormItsdSelectionUi",
    "renderMetabNormItsdSelectionTable",
    "registerMetabNormItsdSelectionTracking",
    "runMetabNormItsdSelectionTrackingObserverShell",
    "runMetabNormItsdSelectionTableObserverShell",
    "renderMetabNormRuvQcUi",
    "renderMetabNormRuvCancorPlot",
    "renderMetabNormRuvOptimizationSummary",
    "renderMetabNormRuvResultsTable",
    "runMetabNormRuvBindingObserverShell",
    "runMetabNormQcImageBindingShell",
    "runMetabNormAssayLabelBindingShell",
    "runMetabNormNormalizationObserverWrapper",
    "runMetabNormResetNormalizationObserverWrapper",
    "runMetabNormApplyCorrelationObserverWrapper",
    "runMetabNormSkipCorrelationObserverWrapper",
    "runMetabNormExportSessionObserverWrapper",
    "buildMetabNormCorrelationFilterSummary",
    "resolveMetabNormFinalQcRenderState",
    "buildMetabNormFinalQcPcaPlot",
    "getPlotAesthetics",
    "appendMetabNormNormalizationLog",
    "renderMetabNormNormalizationLog",
    "renderMetabNormCorrelationFilterSummary",
    "renderMetabNormFinalQcPlot",
    "renderMetabNormAssayLabel",
    "renderMetabNormQcImageForAssay",
    "resolveMetabNormExportSourceDir",
    "collectMetabNormFeatureCountsPerAssay",
    "buildMetabNormExportSessionData",
    "saveMetabNormExportSessionRdsFiles",
    "saveMetabNormExportMetadataFiles",
    "saveMetabNormExportSummaryFile",
    "runMetabNormExportSessionWorkflow",
    "checkMetabNormExportSessionReady",
    "dispatchMetabNormExportSession",
    "handleMetabNormExportSessionOutcome",
    "resolveMetabNormManualItsdFeatureIds",
    "runMetabNormPreNormalizationQcStep",
    "runMetabNormItsdProgressApplyShell",
    "runMetabNormItsdNormalizationStep",
    "runMetabNormLog2ProgressApplyShell",
    "runMetabNormLog2TransformationStep",
    "runMetabNormBetweenSampleProgressApplyShell",
    "runMetabNormBetweenSampleNormalizationStep",
    "runMetabNormPostNormalizationQcStep",
    "runMetabNormRuvProgressApplyShell",
    "runMetabNormRuvOptimizationStep",
    "runMetabNormRuvCorrectionStep",
    "runMetabNormRuvQcStep",
    "buildMetabNormLabelPlot",
    "buildMetabNormTitlePlot",
    "loadMetabNormImageAsPlot",
    "generateMetabNormCompositeFromFiles",
    "generateMetabNormPreNormalizationQc",
    "runMetabNormCompositeQcFigureStep",
    "runMetabNormCompositeQcRefreshShell",
    "runMetabNormNormalizationPipelineShell",
    "runMetabNormNormalizationObserverShell",
    "runMetabNormExportSessionObserverShell",
    "runMetabNormResetNormalizationObserverShell",
    "resolveMetabNormSkipCorrelationInputObject",
    "completeMetabNormSkipCorrelationState",
    "handleMetabNormSkipCorrelationOutcome",
    "dispatchMetabNormSkipCorrelation",
    "runMetabNormSkipCorrelationObserverEntry",
    "runMetabNormSkipCorrelationObserverShell",
    "runMetabNormApplyCorrelationObserverShell",
    "runMetabNormApplyCorrelationWorkflow",
    "handleMetabNormApplyCorrelationOutcome",
    "dispatchMetabNormApplyCorrelation",
    "runMetabNormApplyCorrelationObserverEntry"
  ),
  env = environment()
)

if (!methods::isClass("MetabNormSummaryMock")) {
  methods::setClass("MetabNormSummaryMock", slots = c(design_matrix = "data.frame"))
}

if (!methods::isClass("MetaboliteAssayData")) {
  methods::setClass(
    "MetaboliteAssayData",
    slots = c(
      metabolite_data = "list",
      metabolite_id_column = "character",
      annotation_id_column = "character"
    ),
    prototype = list(
      metabolite_data = list(),
      metabolite_id_column = "metabolite_id",
      annotation_id_column = "annotation_id"
    )
  )
}

test_that("metabolomics normalization module plot aesthetics helper preserves defaults and explicit values", {
  expect_identical(
    getPlotAesthetics(colorVariable = NULL, shapeVariable = ""),
    list(color_var = "group", shape_var = "group")
  )

  expect_identical(
    getPlotAesthetics(colorVariable = "batch", shapeVariable = "factor1"),
    list(color_var = "batch", shape_var = "factor1")
  )
})

test_that("metabolomics normalization module design-driven choices helper preserves available input updates", {
  capture <- new.env(parent = emptyenv())
  capture$updates <- list()

  visible <- withVisible(
    updateMetabNormDesignDrivenChoices(
      session = "session-token",
      designMatrix = data.frame(
        sample_id = c("s1", "s2"),
        batch = c("b1", "b2"),
        group = c("g1", "g2"),
        stringsAsFactors = FALSE
      ),
      updateSelectInputFn = function(session, inputId, choices, selected) {
        capture$updates[[length(capture$updates) + 1L]] <- list(
          session = session,
          input_id = inputId,
          choices = choices,
          selected = selected
        )
      }
    )
  )

  expect_false(visible$visible)
  expect_null(visible$value)
  expect_length(capture$updates, 3L)
  expect_identical(
    capture$updates[[1]],
    list(
      session = "session-token",
      input_id = "color_variable",
      choices = c("sample_id", "batch", "group"),
      selected = "group"
    )
  )
  expect_identical(
    capture$updates[[2]],
    list(
      session = "session-token",
      input_id = "shape_variable",
      choices = c("sample_id", "batch", "group"),
      selected = "group"
    )
  )
  expect_identical(
    capture$updates[[3]],
    list(
      session = "session-token",
      input_id = "ruv_grouping_variable",
      choices = c("batch", "group"),
      selected = "group"
    )
  )
})

test_that("metabolomics normalization module assay-name initialization helper preserves detected assay updates", {
  norm_data <- new.env(parent = emptyenv())
  norm_data$assay_names <- NULL
  norm_data$itsd_selections <- list(Plasma = c("IS1", "IS2"))

  capture <- new.env(parent = emptyenv())
  capture$req <- NULL
  capture$manager <- NULL
  capture$info <- NULL
  capture$warn <- NULL

  current_s4 <- methods::new(
    "MetaboliteAssayData",
    metabolite_data = list(
      Plasma = data.frame(
        metabolite_id = "m1",
        annotation_id = "a1",
        sample_1 = 1,
        stringsAsFactors = FALSE
      ),
      Urine = data.frame(
        metabolite_id = "m2",
        annotation_id = "a2",
        sample_1 = 2,
        stringsAsFactors = FALSE
      )
    )
  )

  visible <- withVisible(
    initializeMetabNormAssayNames(
      stateManager = "state-manager",
      normData = norm_data,
      reqFn = function(value) {
        capture$req <- value
        value
      },
      getStateFn = function(manager) {
        capture$manager <- manager
        current_s4
      },
      logInfoFn = function(message) capture$info <- message,
      logWarnFn = function(message) capture$warn <- message
    )
  )

  expect_false(visible$visible)
  expect_identical(visible$value, c("Plasma", "Urine"))
  expect_identical(capture$req, "state-manager")
  expect_identical(capture$manager, "state-manager")
  expect_identical(norm_data$assay_names, c("Plasma", "Urine"))
  expect_identical(norm_data$itsd_selections, list(Plasma = c("IS1", "IS2")))
  expect_identical(capture$info, "Detected assays: Plasma, Urine")
  expect_null(capture$warn)
})

test_that("metabolomics normalization module assay-name initialization helper preserves warning fallback", {
  norm_data <- new.env(parent = emptyenv())
  norm_data$assay_names <- "Existing"
  norm_data$itsd_selections <- list(Existing = "IS1")

  capture <- new.env(parent = emptyenv())
  capture$info <- NULL
  capture$warn <- NULL

  visible <- withVisible(
    initializeMetabNormAssayNames(
      stateManager = "state-manager",
      normData = norm_data,
      reqFn = identity,
      getStateFn = function(manager) {
        expect_identical(manager, "state-manager")
        stop("boom")
      },
      logInfoFn = function(message) capture$info <- message,
      logWarnFn = function(message) capture$warn <- message
    )
  )

  expect_false(visible$visible)
  expect_null(visible$value)
  expect_identical(norm_data$assay_names, "Existing")
  expect_identical(norm_data$itsd_selections, list(Existing = "IS1"))
  expect_null(capture$info)
  expect_identical(capture$warn, "Could not detect assay names: boom")
})

test_that("metabolomics normalization module log append helper preserves timestamped log mutation", {
  norm_data <- new.env(parent = emptyenv())
  norm_data$normalization_log <- "[00:00:00] existing"

  visible <- withVisible(
    appendMetabNormNormalizationLog(
      normData = norm_data,
      message = "new entry",
      timestampFn = function() "07:08:09"
    )
  )

  expect_false(visible$visible)
  expect_identical(
    norm_data$normalization_log,
    c("[00:00:00] existing", "[07:08:09] new entry")
  )
  expect_identical(visible$value, norm_data$normalization_log)
})

test_that("metabolomics normalization module log render helper preserves empty and populated output text", {
  render_text_capture <- function(expr) {
    eval.parent(substitute(expr))
  }

  empty_norm_data <- list(normalization_log = character(0))
  expect_identical(
    renderMetabNormNormalizationLog(
      normData = empty_norm_data,
      renderTextFn = render_text_capture
    ),
    "Normalization log will appear here as you apply steps..."
  )

  populated_norm_data <- list(
    normalization_log = c("[07:08:09] first entry", "[07:08:10] second entry")
  )
  expect_identical(
    renderMetabNormNormalizationLog(
      normData = populated_norm_data,
      renderTextFn = render_text_capture
    ),
    "[07:08:09] first entry\n[07:08:10] second entry"
  )
})

test_that("metabolomics normalization module correlation summary render helper preserves fallback and source selection", {
  render_text_capture <- function(expr) {
    eval.parent(substitute(expr))
  }

  expect_identical(
    renderMetabNormCorrelationFilterSummary(
      normData = list(correlation_filtering_complete = FALSE),
      renderTextFn = render_text_capture,
      buildSummaryFn = function(...) {
        stop("buildSummaryFn should not run before filtering completes")
      }
    ),
    "Apply correlation filter to see results..."
  )

  filtered_object <- list(stage = "filtered")
  ruv_corrected_object <- list(stage = "ruv")
  post_norm_object <- list(stage = "post_norm")
  corr_results <- list(Plasma = data.frame(pearson_correlation = c(0.8, 0.9)))
  capture <- new.env(parent = emptyenv())

  expect_identical(
    renderMetabNormCorrelationFilterSummary(
      normData = list(
        correlation_filtering_complete = TRUE,
        correlation_results = corr_results,
        correlation_filtered_obj = filtered_object,
        ruv_corrected_obj = ruv_corrected_object,
        post_norm_obj = post_norm_object
      ),
      renderTextFn = render_text_capture,
      buildSummaryFn = function(corrResults, filteredObject, originalObject) {
        capture$corr_results <- corrResults
        capture$filtered_object <- filteredObject
        capture$original_object <- originalObject
        "summary from builder"
      }
    ),
    "summary from builder"
  )
  expect_identical(capture$corr_results, corr_results)
  expect_identical(capture$filtered_object, filtered_object)
  expect_identical(capture$original_object, ruv_corrected_object)

  expect_identical(
    renderMetabNormCorrelationFilterSummary(
      normData = list(
        correlation_filtering_complete = TRUE,
        correlation_results = corr_results,
        correlation_filtered_obj = filtered_object,
        ruv_corrected_obj = NULL,
        post_norm_obj = post_norm_object
      ),
      renderTextFn = render_text_capture,
      buildSummaryFn = function(corrResults, filteredObject, originalObject) {
        capture$fallback_original_object <- originalObject
        "summary without ruv"
      }
    ),
    "summary without ruv"
  )
  expect_identical(capture$fallback_original_object, post_norm_object)
})

test_that("metabolomics normalization module final QC render helper preserves readiness gate", {
  render_plot_capture <- function(expr) {
    eval.parent(substitute(expr))
  }
  capture <- new.env(parent = emptyenv())

  expect_error(
    renderMetabNormFinalQcPlot(
      normData = list(
        correlation_filtering_complete = FALSE,
        ruv_complete = FALSE
      ),
      renderPlotFn = render_plot_capture,
      reqFn = function(value) {
        capture$req <- value
        if (!isTRUE(value)) {
          stop("final qc not ready")
        }
      }
    ),
    "final qc not ready"
  )

  expect_false(isTRUE(capture$req))
})

test_that("metabolomics normalization module final QC render helper preserves fallback return", {
  render_plot_capture <- function(expr) {
    eval.parent(substitute(expr))
  }

  fallback_plot <- ggplot2::ggplot(data.frame(x = 1, y = 1), ggplot2::aes(x, y)) +
    ggplot2::geom_point()
  capture <- new.env(parent = emptyenv())

  expect_identical(
    renderMetabNormFinalQcPlot(
      normData = list(
        correlation_filtering_complete = TRUE,
        ruv_complete = FALSE,
        correlation_filtered_obj = "filtered-object",
        ruv_corrected_obj = "ruv-object",
        post_norm_obj = "post-object"
      ),
      colorVariableFn = function() {
        stop("colorVariableFn should not run for fallback renders")
      },
      shapeVariableFn = function() {
        stop("shapeVariableFn should not run for fallback renders")
      },
      renderPlotFn = render_plot_capture,
      resolveRenderStateFn = function(correlationFilteredObject, ruvCorrectedObject, postNormObject) {
        capture$state_args <- list(
          correlation_filtered_object = correlationFilteredObject,
          ruv_corrected_object = ruvCorrectedObject,
          post_norm_object = postNormObject
        )
        list(
          sourceObject = NULL,
          sourceStage = "empty",
          isFallback = TRUE,
          plot = fallback_plot
        )
      },
      getPlotAestheticsFn = function(...) {
        stop("getPlotAestheticsFn should not run for fallback renders")
      },
      buildPcaPlotFn = function(...) {
        stop("buildPcaPlotFn should not run for fallback renders")
      }
    ),
    fallback_plot
  )

  expect_identical(
    capture$state_args,
    list(
      correlation_filtered_object = "filtered-object",
      ruv_corrected_object = "ruv-object",
      post_norm_object = "post-object"
    )
  )
})

test_that("metabolomics normalization module final QC render helper preserves aesthetics and builder delegation", {
  render_plot_capture <- function(expr) {
    eval.parent(substitute(expr))
  }

  capture <- new.env(parent = emptyenv())
  source_object <- structure(list(stage = "correlation_filter"), class = "MetabNormRenderMock")

  expect_identical(
    renderMetabNormFinalQcPlot(
      normData = list(
        correlation_filtering_complete = TRUE,
        ruv_complete = FALSE,
        correlation_filtered_obj = source_object,
        ruv_corrected_obj = "ruv-object",
        post_norm_obj = "post-object"
      ),
      colorVariableFn = function() {
        capture$color_requested <- TRUE
        "Condition"
      },
      shapeVariableFn = function() {
        capture$shape_requested <- TRUE
        "Batch"
      },
      renderPlotFn = render_plot_capture,
      reqFn = function(value) {
        capture$req <- value
        invisible(value)
      },
      resolveRenderStateFn = function(correlationFilteredObject, ruvCorrectedObject, postNormObject) {
        capture$state_args <- list(
          correlation_filtered_object = correlationFilteredObject,
          ruv_corrected_object = ruvCorrectedObject,
          post_norm_object = postNormObject
        )
        list(
          sourceObject = correlationFilteredObject,
          sourceStage = "correlation_filter",
          isFallback = FALSE,
          plot = NULL
        )
      },
      getPlotAestheticsFn = function(colorVariable, shapeVariable) {
        capture$aesthetics_args <- list(
          color_variable = colorVariable,
          shape_variable = shapeVariable
        )
        list(color_var = "resolved-color", shape_var = "resolved-shape")
      },
      buildPcaPlotFn = function(sourceObject, colorVar, shapeVar) {
        capture$builder_args <- list(
          source_object = sourceObject,
          color_var = colorVar,
          shape_var = shapeVar
        )
        "built-final-qc-plot"
      }
    ),
    "built-final-qc-plot"
  )

  expect_true(isTRUE(capture$req))
  expect_true(isTRUE(capture$color_requested))
  expect_true(isTRUE(capture$shape_requested))
  expect_identical(
    capture$state_args,
    list(
      correlation_filtered_object = source_object,
      ruv_corrected_object = "ruv-object",
      post_norm_object = "post-object"
    )
  )
  expect_identical(
    capture$aesthetics_args,
    list(color_variable = "Condition", shape_variable = "Batch")
  )
  expect_identical(
    capture$builder_args,
    list(
      source_object = source_object,
      color_var = "resolved-color",
      shape_var = "resolved-shape"
    )
  )
})

test_that("metabolomics normalization module assay label helper preserves detected and fallback labels", {
  render_text_capture <- function(expr) {
    eval.parent(substitute(expr))
  }

  expect_identical(
    renderMetabNormAssayLabel(
      assaySlot = 1,
      getAssayNamesFn = function() c("Plasma", "Urine"),
      renderTextFn = render_text_capture
    ),
    "Assay: Plasma"
  )

  expect_identical(
    renderMetabNormAssayLabel(
      assaySlot = 2,
      getAssayNamesFn = function() "Plasma",
      renderTextFn = render_text_capture
    ),
    "Assay 2: (detecting...)"
  )

  expect_identical(
    renderMetabNormAssayLabel(
      assaySlot = 3,
      getAssayNamesFn = function() NULL,
      renderTextFn = render_text_capture
    ),
    "Assay 3: (detecting...)"
  )
})

test_that("metabolomics normalization module QC image helper preserves assay resolution and fallbacks", {
  capture <- new.env(parent = emptyenv())
  capture$refresh_reads <- 0L

  render_image_capture <- function(expr, deleteFile = FALSE) {
    capture$delete_file <- deleteFile
    eval.parent(substitute(expr))
  }

  norm_data <- new.env(parent = emptyenv())
  norm_data$assay_names <- c("Plasma Panel", "Urine")
  makeActiveBinding("plot_refresh_trigger", function() {
    capture$refresh_reads <- capture$refresh_reads + 1L
    5L
  }, norm_data)

  expect_identical(
    renderMetabNormQcImageForAssay(
      assaySlot = 1,
      plotType = "pca",
      stagePrefix = "pre_norm",
      normData = norm_data,
      qcDir = "/tmp/metab-qc",
      renderImageFn = render_image_capture,
      fileExistsFn = function(path) {
        capture$existing_path <- path
        TRUE
      }
    ),
    list(
      src = "/tmp/metab-qc/plasma_panel_pre_norm_pca.png",
      contentType = "image/png",
      width = "100%",
      height = "auto",
      alt = "pca - Plasma Panel"
    )
  )
  expect_identical(capture$existing_path, "/tmp/metab-qc/plasma_panel_pre_norm_pca.png")
  expect_identical(capture$refresh_reads, 1L)
  expect_identical(capture$delete_file, FALSE)

  unresolved_norm_data <- list(
    plot_refresh_trigger = 0L,
    assay_names = "Plasma"
  )
  expect_identical(
    renderMetabNormQcImageForAssay(
      assaySlot = 2,
      plotType = "density",
      stagePrefix = "post_norm",
      normData = unresolved_norm_data,
      qcDir = "/tmp/metab-qc",
      renderImageFn = render_image_capture,
      fileExistsFn = function(...) {
        stop("fileExistsFn should not run when the assay slot is unresolved")
      }
    ),
    list(src = "", alt = "Assay not detected yet")
  )

  expect_identical(
    renderMetabNormQcImageForAssay(
      assaySlot = 1,
      plotType = "density",
      stagePrefix = "post_norm",
      normData = list(plot_refresh_trigger = 0L, assay_names = "Plasma"),
      qcDir = NULL,
      renderImageFn = render_image_capture,
      fileExistsFn = function(...) {
        stop("fileExistsFn should not run when the QC directory is missing")
      }
    ),
    list(src = "", alt = "QC directory not configured")
  )

  expect_identical(
    renderMetabNormQcImageForAssay(
      assaySlot = 1,
      plotType = "rle",
      stagePrefix = "ruv_corrected",
      normData = list(plot_refresh_trigger = 0L, assay_names = "Plasma"),
      qcDir = "/tmp/metab-qc",
      renderImageFn = render_image_capture,
      fileExistsFn = function(path) {
        capture$missing_path <- path
        FALSE
      }
    ),
    list(
      src = "",
      alt = "Plot not generated yet: plasma_ruv_corrected_rle.png"
    )
  )
  expect_identical(capture$missing_path, "/tmp/metab-qc/plasma_ruv_corrected_rle.png")
})

test_that("metabolomics normalization module manual ITSD helper returns NULL when no selections are present", {
  capture <- new.env(parent = emptyenv())
  capture$build_called <- FALSE

  current_s4 <- methods::new(
    "MetaboliteAssayData",
    metabolite_data = list(
      Plasma = data.frame(
        metabolite_id = "M1",
        annotation_id = "A1",
        stringsAsFactors = FALSE
      )
    )
  )

  result <- resolveMetabNormManualItsdFeatureIds(
    currentS4 = current_s4,
    itsdSelections = list(Plasma = integer(0), Urine = NULL),
    addLogFn = function(message) {
      stop(sprintf("addLogFn should not be called: %s", message))
    },
    buildSelectionTableFn = function(...) {
      capture$build_called <- TRUE
      data.frame(feature_id = character(), stringsAsFactors = FALSE)
    },
    mapSelectionsFn = function(...) {
      stop("mapSelectionsFn should not be called when selections are empty")
    },
    compactFn = function(...) {
      stop("compactFn should not be called when selections are empty")
    }
  )

  expect_null(result)
  expect_false(capture$build_called)
})

test_that("metabolomics normalization module manual ITSD helper preserves feature-ID lookup and logging", {
  capture <- new.env(parent = emptyenv())
  capture$build <- list()
  capture$log <- character()

  current_s4 <- methods::new(
    "MetaboliteAssayData",
    metabolite_data = list(
      Plasma = data.frame(
        metabolite_id = c("M1", "M2", "M3"),
        annotation_id = c("A1", "A2", "A3"),
        stringsAsFactors = FALSE
      )
    )
  )

  result <- resolveMetabNormManualItsdFeatureIds(
    currentS4 = current_s4,
    itsdSelections = list(Plasma = c(1L, 3L), Missing = 1L, Urine = integer(0)),
    addLogFn = function(message) {
      capture$log <- c(capture$log, message)
    },
    buildSelectionTableFn = function(assay_data, metabolite_id_col, annotation_cols) {
      capture$build <- c(
        capture$build,
        list(list(
          assay_data = assay_data,
          metabolite_id_col = metabolite_id_col,
          annotation_cols = annotation_cols
        ))
      )
      data.frame(feature_id = assay_data[[metabolite_id_col]], stringsAsFactors = FALSE)
    },
    mapSelectionsFn = function(selectionList, fn) {
      output <- lapply(names(selectionList), function(assayName) {
        fn(selectionList[[assayName]], assayName)
      })
      stats::setNames(output, names(selectionList))
    },
    compactFn = function(values) {
      values[!vapply(values, is.null, logical(1))]
    }
  )

  expect_identical(result, list(Plasma = c("M1", "M3")))
  expect_length(capture$build, 1)
  expect_identical(capture$build[[1]]$metabolite_id_col, "metabolite_id")
  expect_identical(capture$build[[1]]$annotation_cols, "annotation_id")
  expect_identical(capture$log, "Assay Plasma : 2 ITSD features selected")
})

test_that("metabolomics normalization module pre-normalization QC helper preserves capture and generation shell", {
  capture <- new.env(parent = emptyenv())
  capture$log <- character()

  current_s4 <- structure(list(stage = "post_filter"), class = "MetabNormPreQcStepMock")
  experiment_paths <- list(metabolite_qc_dir = "qc-dir")
  norm_data <- new.env(parent = emptyenv())

  visible <- withVisible(
    runMetabNormPreNormalizationQcStep(
      currentS4 = current_s4,
      totalSteps = 6,
      experimentPaths = experiment_paths,
      groupingVariable = "batch",
      shapeVariable = "group",
      normData = norm_data,
      addLogFn = function(message) {
        capture$log <- c(capture$log, message)
      },
      incProgressFn = function(amount, detail = NULL) {
        capture$progress <- list(amount = amount, detail = detail)
      },
      generateMetabQcPlotsFn = function(
        theObject,
        experiment_paths,
        stage,
        grouping_variable,
        shape_variable
      ) {
        capture$plot_call <- list(
          theObject = theObject,
          experiment_paths = experiment_paths,
          stage = stage,
          grouping_variable = grouping_variable,
          shape_variable = shape_variable
        )
        invisible(NULL)
      }
    )
  )

  expect_false(visible$visible)
  expect_identical(
    capture$progress,
    list(amount = 1 / 6, detail = "Capturing pre-normalization state...")
  )
  expect_identical(norm_data$post_filter_obj, current_s4)
  expect_identical(
    capture$plot_call,
    list(
      theObject = current_s4,
      experiment_paths = experiment_paths,
      stage = "post_filter",
      grouping_variable = "batch",
      shape_variable = "group"
    )
  )
  expect_identical(
    capture$log,
    c("Post-filtering state captured", "Pre-normalization QC plots generated")
  )
  expect_identical(
    visible$value,
    list(
      currentS4 = current_s4,
      stage = "post_filter",
      progressDetail = "Capturing pre-normalization state...",
      captureLogEntry = "Post-filtering state captured",
      preQcLogEntry = "Pre-normalization QC plots generated"
    )
  )
})

test_that("metabolomics normalization module auto pre-normalization QC seam preserves assay detection and plot refresh shell", {
  capture <- new.env(parent = emptyenv())
  capture$add_log <- character()
  capture$info <- character()
  capture$warn <- character()
  capture$error <- character()
  capture$req <- list()

  current_s4 <- methods::new(
    "MetaboliteAssayData",
    metabolite_data = list(
      Plasma = data.frame(
        metabolite_id = "M1",
        annotation_id = "A1",
        stringsAsFactors = FALSE
      ),
      Urine = data.frame(
        metabolite_id = "M2",
        annotation_id = "A2",
        stringsAsFactors = FALSE
      )
    )
  )
  workflow_data <- list(
    state_manager = list(
      getState = function() current_s4
    )
  )
  norm_data <- new.env(parent = emptyenv())
  norm_data$assay_names <- NULL
  norm_data$plot_refresh_trigger <- 4
  norm_data$pre_norm_qc_generated <- FALSE

  generateMetabNormPreNormalizationQc(
    workflowData = workflow_data,
    experimentPaths = list(metabolite_qc_dir = "qc-dir"),
    normData = norm_data,
    getPlotAestheticsFn = function() {
      capture$aesthetics_called <- TRUE
      list(color_var = "batch", shape_var = "group")
    },
    addLogFn = function(message) {
      capture$add_log <- c(capture$add_log, message)
    },
    reqFn = function(value) {
      capture$req <- c(capture$req, list(value))
      invisible(value)
    },
    generateMetabQcPlotsFn = function(
      theObject,
      experiment_paths,
      stage,
      grouping_variable,
      shape_variable
    ) {
      capture$plot_call <- list(
        theObject = theObject,
        experiment_paths = experiment_paths,
        stage = stage,
        grouping_variable = grouping_variable,
        shape_variable = shape_variable
      )
      invisible(NULL)
    },
    logInfoFn = function(message) {
      capture$info <- c(capture$info, message)
    },
    logWarnFn = function(message) {
      capture$warn <- c(capture$warn, message)
    },
    logErrorFn = function(message) {
      capture$error <- c(capture$error, message)
    }
  )

  expect_true(isTRUE(capture$aesthetics_called))
  expect_identical(capture$req, list(workflow_data$state_manager))
  expect_identical(norm_data$assay_names, c("Plasma", "Urine"))
  expect_identical(norm_data$plot_refresh_trigger, 5)
  expect_true(norm_data$pre_norm_qc_generated)
  expect_identical(
    capture$plot_call,
    list(
      theObject = current_s4,
      experiment_paths = list(metabolite_qc_dir = "qc-dir"),
      stage = "post_filter",
      grouping_variable = "batch",
      shape_variable = "group"
    )
  )
  expect_identical(
    capture$info,
    c(
      "=== GENERATING PRE-NORMALIZATION QC PLOTS ===",
      "Set assay names: Plasma, Urine",
      "Pre-normalization QC plots generated successfully"
    )
  )
  expect_identical(capture$add_log, character())
  expect_identical(capture$warn, character())
  expect_identical(capture$error, character())
})

test_that("metabolomics normalization module auto pre-normalization QC seam preserves error logging tail", {
  capture <- new.env(parent = emptyenv())
  capture$add_log <- character()
  capture$info <- character()
  capture$warn <- character()
  capture$error <- character()

  current_s4 <- methods::new(
    "MetaboliteAssayData",
    metabolite_data = list(
      Plasma = data.frame(
        metabolite_id = "M1",
        annotation_id = "A1",
        stringsAsFactors = FALSE
      )
    )
  )
  workflow_data <- list(
    state_manager = list(
      getState = function() current_s4
    )
  )
  norm_data <- new.env(parent = emptyenv())
  norm_data$assay_names <- NULL
  norm_data$plot_refresh_trigger <- 0
  norm_data$pre_norm_qc_generated <- FALSE

  generateMetabNormPreNormalizationQc(
    workflowData = workflow_data,
    experimentPaths = list(metabolite_qc_dir = "qc-dir"),
    normData = norm_data,
    getPlotAestheticsFn = function() list(color_var = "batch", shape_var = "group"),
    addLogFn = function(message) {
      capture$add_log <- c(capture$add_log, message)
    },
    reqFn = function(value) invisible(value),
    generateMetabQcPlotsFn = function(...) {
      stop("qc failed", call. = FALSE)
    },
    logInfoFn = function(message) {
      capture$info <- c(capture$info, message)
    },
    logWarnFn = function(message) {
      capture$warn <- c(capture$warn, message)
    },
    logErrorFn = function(message) {
      capture$error <- c(capture$error, message)
    }
  )

  expect_identical(norm_data$assay_names, "Plasma")
  expect_identical(norm_data$plot_refresh_trigger, 0)
  expect_false(norm_data$pre_norm_qc_generated)
  expect_identical(
    capture$info,
    c(
      "=== GENERATING PRE-NORMALIZATION QC PLOTS ===",
      "Set assay names: Plasma"
    )
  )
  expect_identical(capture$warn, character())
  expect_identical(
    capture$error,
    "Error generating pre-normalization QC: qc failed"
  )
  expect_identical(capture$add_log, "Error generating Pre-QC: qc failed")
})

test_that("metabolomics normalization module auto pre-normalization QC observer shell preserves tab-gated handoff", {
  capture <- new.env(parent = emptyenv())
  capture$info <- character()
  capture$req <- list()

  current_s4 <- methods::new(
    "MetaboliteAssayData",
    metabolite_data = list(
      Plasma = data.frame(
        metabolite_id = "M1",
        annotation_id = "A1",
        stringsAsFactors = FALSE
      )
    )
  )
  norm_data <- new.env(parent = emptyenv())
  norm_data$pre_norm_qc_generated <- FALSE
  plot_generator <- function(...) invisible(NULL)

  visible <- withVisible(
    runMetabNormAutoPreNormalizationQcObserverShell(
      selectedTab = "norm",
      workflowData = list(state_manager = "state-manager"),
      experimentPaths = list(metabolite_qc_dir = "qc-dir"),
      normData = norm_data,
      colorVariable = "batch",
      shapeVariable = "group",
      addLogFn = function(message) {
        capture$add_log <- c(capture$add_log, message)
      },
      reqFn = function(value) {
        capture$req <- c(capture$req, list(value))
        invisible(value)
      },
      getStateFn = function(manager) {
        capture$manager <- manager
        current_s4
      },
      withProgressFn = function(message, value, expr, detail = NULL) {
        capture$with_progress <- list(
          message = message,
          value = value,
          detail = detail
        )
        expr
      },
      generatePreNormalizationQcFn = function(
        workflowData,
        experimentPaths,
        normData,
        getPlotAestheticsFn,
        addLogFn,
        reqFn,
        generateMetabQcPlotsFn,
        logInfoFn,
        logWarnFn,
        logErrorFn
      ) {
        capture$generate <- list(
          workflowData = workflowData,
          experimentPaths = experimentPaths,
          normData = normData,
          aesthetics = getPlotAestheticsFn(),
          generateMetabQcPlotsFn = generateMetabQcPlotsFn
        )
        invisible(NULL)
      },
      generateMetabQcPlotsFn = plot_generator,
      logInfoFn = function(message) {
        capture$info <- c(capture$info, message)
      },
      logWarnFn = function(message) {
        stop(sprintf("logWarnFn should not be called: %s", message))
      },
      logErrorFn = function(message) {
        stop(sprintf("logErrorFn should not be called: %s", message))
      }
    )
  )

  expect_false(visible$visible)
  expect_identical(capture$req, list("state-manager"))
  expect_identical(capture$manager, "state-manager")
  expect_identical(
    capture$with_progress,
    list(
      message = "Generating Pre-Normalization QC...",
      value = 0.5,
      detail = NULL
    )
  )
  expect_identical(capture$generate$workflowData, list(state_manager = "state-manager"))
  expect_identical(capture$generate$experimentPaths, list(metabolite_qc_dir = "qc-dir"))
  expect_identical(capture$generate$normData, norm_data)
  expect_identical(
    capture$generate$aesthetics,
    list(color_var = "batch", shape_var = "group")
  )
  expect_identical(capture$generate$generateMetabQcPlotsFn, plot_generator)
  expect_identical(
    capture$info,
    c(
      "Normalization tab selected - checking if pre-QC needed",
      "Auto-triggering pre-normalization QC plots"
    )
  )
  expect_true(norm_data$pre_norm_qc_generated)
  expect_identical(
    visible$value,
    list(
      selectedTab = "norm",
      currentS4 = current_s4,
      preNormQcGenerated = TRUE
    )
  )
})

test_that("metabolomics normalization module auto pre-normalization QC observer shell preserves generated short-circuit", {
  capture <- new.env(parent = emptyenv())
  capture$info <- character()

  norm_data <- new.env(parent = emptyenv())
  norm_data$pre_norm_qc_generated <- TRUE

  visible <- withVisible(
    runMetabNormAutoPreNormalizationQcObserverShell(
      selectedTab = "norm",
      workflowData = list(state_manager = "state-manager"),
      experimentPaths = list(metabolite_qc_dir = "qc-dir"),
      normData = norm_data,
      reqFn = function(value) {
        capture$req <- c(capture$req, list(value))
        invisible(value)
      },
      getStateFn = function(manager) {
        capture$manager <- manager
        stop("getStateFn should not run when pre-QC is already generated")
      },
      withProgressFn = function(...) {
        capture$with_progress <- TRUE
        invisible(NULL)
      },
      generatePreNormalizationQcFn = function(...) {
        capture$generate <- TRUE
        invisible(NULL)
      },
      logInfoFn = function(message) {
        capture$info <- c(capture$info, message)
      }
    )
  )

  expect_false(visible$visible)
  expect_null(visible$value)
  expect_identical(
    capture$info,
    "Normalization tab selected - checking if pre-QC needed"
  )
  expect_null(capture$req)
  expect_null(capture$manager)
  expect_null(capture$with_progress)
  expect_null(capture$generate)
  expect_true(norm_data$pre_norm_qc_generated)
})

test_that("metabolomics normalization module ITSD selection render seam preserves assay panel rendering", {
  norm_data <- new.env(parent = emptyenv())
  norm_data$assay_names <- c("LCMS Pos", "GC-MS/Neg")

  capture <- new.env(parent = emptyenv())
  capture$req <- list()
  capture$output_ids <- character()

  render_fn <- renderMetabNormItsdSelectionUi(
    normData = norm_data,
    ns = function(id) paste0("module-", id),
    renderUIFn = function(expr) {
      render_env <- parent.frame()
      render_expr <- substitute(expr)

      function() {
        eval(render_expr, envir = render_env)
      }
    },
    reqFn = function(value) {
      capture$req <- c(capture$req, list(value))
      invisible(value)
    },
    mapFn = function(.x, .f) lapply(.x, .f),
    wellPanelFn = function(...) list(tag = "wellPanel", children = list(...)),
    h5Fn = function(text) list(tag = "h5", text = text),
    dataTableOutputFn = function(outputId) {
      capture$output_ids <- c(capture$output_ids, outputId)
      list(tag = "dataTableOutput", outputId = outputId)
    },
    brFn = function() list(tag = "br"),
    tagListFn = function(children) list(tag = "tagList", children = children)
  )

  expect_true(is.function(render_fn))

  rendered_ui <- render_fn()

  expect_identical(capture$req, list(c("LCMS Pos", "GC-MS/Neg")))
  expect_identical(
    capture$output_ids,
    c("module-itsd_table_lcms_pos", "module-itsd_table_gc_ms_neg")
  )
  expect_identical(rendered_ui$tag, "tagList")
  expect_length(rendered_ui$children, 2L)
  expect_identical(
    rendered_ui$children[[1]]$children[[1]],
    list(tag = "h5", text = "Assay: LCMS Pos")
  )
  expect_identical(
    rendered_ui$children[[1]]$children[[2]],
    list(tag = "dataTableOutput", outputId = "module-itsd_table_lcms_pos")
  )
  expect_identical(rendered_ui$children[[1]]$children[[3]], list(tag = "br"))
  expect_identical(
    rendered_ui$children[[2]]$children[[1]],
    list(tag = "h5", text = "Assay: GC-MS/Neg")
  )
  expect_identical(
    rendered_ui$children[[2]]$children[[2]],
    list(tag = "dataTableOutput", outputId = "module-itsd_table_gc_ms_neg")
  )
})

test_that("metabolomics normalization module ITSD selection table render seam preserves success and null exits", {
  render_data_table_capture <- function(expr) {
    render_env <- parent.frame()
    render_expr <- substitute(expr)

    function() {
      eval(render_expr, envir = render_env)
    }
  }

  success_capture <- new.env(parent = emptyenv())
  selection_table <- data.frame(
    feature_id = c("F1", "F2", "F3"),
    annotation = c("IS A", "Metab B", "IS C"),
    mean_intensity = c(10, 20, 30),
    cv_percent = c(5, 15, 25),
    is_candidate = c(TRUE, FALSE, TRUE),
    stringsAsFactors = FALSE
  )
  assay_data <- data.frame(feature = c("raw-1", "raw-2"))
  current_s4 <- methods::new(
    "MetaboliteAssayData",
    metabolite_data = list("LCMS Pos" = assay_data),
    metabolite_id_column = "feature_id",
    annotation_id_column = "annotation_id"
  )

  success_render <- renderMetabNormItsdSelectionTable(
    assayName = "LCMS Pos",
    currentS4 = current_s4,
    renderDataTableFn = render_data_table_capture,
    buildItsdSelectionTableFn = function(assay_data, metabolite_id_col, annotation_cols) {
      success_capture$build_args <- list(
        assay_data = assay_data,
        metabolite_id_col = metabolite_id_col,
        annotation_cols = annotation_cols
      )
      selection_table
    },
    datatableFn = function(data, selection = NULL, filter = NULL, options = NULL, rownames = NULL) {
      success_capture$datatable_args <- list(
        data = data,
        selection = selection,
        filter = filter,
        options = options,
        rownames = rownames
      )
      list(
        tag = "datatable",
        data = data,
        selection = selection,
        filter = filter,
        options = options,
        rownames = rownames
      )
    },
    formatStyleFn = function(table, columns, backgroundColor = NULL) {
      success_capture$format_style_args <- list(
        columns = columns,
        backgroundColor = backgroundColor
      )
      table$style <- list(columns = columns, backgroundColor = backgroundColor)
      table
    },
    styleEqualFn = function(levels, values) {
      success_capture$style_equal_args <- list(levels = levels, values = values)
      list(levels = levels, values = values)
    },
    formatRoundFn = function(table, columns, digits) {
      success_capture$format_round_args <- list(columns = columns, digits = digits)
      table$round <- list(columns = columns, digits = digits)
      table
    }
  )

  expect_true(is.function(success_render))
  success_value <- success_render()
  expect_identical(success_capture$build_args$assay_data, assay_data)
  expect_identical(success_capture$build_args$metabolite_id_col, "feature_id")
  expect_identical(success_capture$build_args$annotation_cols, "annotation_id")
  expect_identical(success_capture$datatable_args$data, selection_table)
  expect_identical(
    success_capture$datatable_args$selection,
    list(mode = "multiple", selected = c(1L, 3L))
  )
  expect_identical(success_capture$datatable_args$filter, "top")
  expect_identical(
    success_capture$datatable_args$options,
    list(
      pageLength = 10,
      scrollX = TRUE,
      order = list(list(4, "desc"), list(3, "asc"))
    )
  )
  expect_false(success_capture$datatable_args$rownames)
  expect_identical(success_capture$format_style_args$columns, "is_candidate")
  expect_identical(
    success_capture$format_style_args$backgroundColor,
    list(levels = TRUE, values = "#d4edda")
  )
  expect_identical(success_capture$style_equal_args, list(levels = TRUE, values = "#d4edda"))
  expect_identical(
    success_capture$format_round_args,
    list(columns = c("mean_intensity", "cv_percent"), digits = 2)
  )
  expect_identical(
    success_value,
    list(
      tag = "datatable",
      data = selection_table,
      selection = list(mode = "multiple", selected = c(1L, 3L)),
      filter = "top",
      options = list(
        pageLength = 10,
        scrollX = TRUE,
        order = list(list(4, "desc"), list(3, "asc"))
      ),
      rownames = FALSE,
      style = list(
        columns = "is_candidate",
        backgroundColor = list(levels = TRUE, values = "#d4edda")
      ),
      round = list(columns = c("mean_intensity", "cv_percent"), digits = 2)
    )
  )

  failure_capture <- new.env(parent = emptyenv())
  missing_render <- renderMetabNormItsdSelectionTable(
    assayName = "Missing",
    currentS4 = current_s4,
    renderDataTableFn = render_data_table_capture,
    buildItsdSelectionTableFn = function(...) {
      failure_capture$build_called <- TRUE
      stop("buildItsdSelectionTableFn should not run for missing assay data")
    }
  )

  expect_true(is.function(missing_render))
  expect_null(missing_render())
  expect_false(isTRUE(failure_capture$build_called))
})

test_that("metabolomics normalization module ITSD selection tracking seam preserves per-assay observer registration", {
  norm_data <- new.env(parent = emptyenv())
  norm_data$assay_names <- c("LCMS Pos", "GC-MS/Neg")
  norm_data$itsd_selections <- list()

  input_values <- list(
    itsd_table_lcms_pos_rows_selected = c(1L, 3L),
    itsd_table_gc_ms_neg_rows_selected = integer(0)
  )

  capture <- new.env(parent = emptyenv())
  capture$events <- list()
  capture$messages <- character()

  visible <- withVisible(
    registerMetabNormItsdSelectionTracking(
      normData = norm_data,
      input = input_values,
      walkFn = function(.x, .f) lapply(.x, .f),
      observeEventFn = function(eventExpr, handlerExpr, ignoreNULL = TRUE) {
        event_env <- parent.frame()
        event_expr <- substitute(eventExpr)
        handler_env <- parent.frame()
        handler_expr <- substitute(handlerExpr)

        capture$events[[length(capture$events) + 1L]] <- list(
          value = eval(event_expr, envir = event_env),
          ignoreNULL = ignoreNULL
        )

        eval(handler_expr, envir = handler_env)
        invisible(NULL)
      },
      logInfoFn = function(message) {
        capture$messages <- c(capture$messages, message)
      }
    )
  )

  expect_false(visible$visible)
  expect_null(visible$value)
  expect_length(capture$events, 2L)
  expect_identical(
    capture$events[[1]],
    list(value = c(1L, 3L), ignoreNULL = FALSE)
  )
  expect_identical(
    capture$events[[2]],
    list(value = integer(0), ignoreNULL = FALSE)
  )
  expect_identical(
    norm_data$itsd_selections,
    list("LCMS Pos" = c(1L, 3L), "GC-MS/Neg" = integer(0))
  )
  expect_identical(
    capture$messages,
    c(
      "ITSD selection updated for LCMS Pos : 2 features selected",
      "ITSD selection updated for GC-MS/Neg : 0 features selected"
    )
  )
})

test_that("metabolomics normalization module ITSD selection tracking seam preserves null-selection clearing", {
  norm_data <- new.env(parent = emptyenv())
  norm_data$assay_names <- "LCMS Pos"
  norm_data$itsd_selections <- list("LCMS Pos" = c(2L, 4L))

  capture <- new.env(parent = emptyenv())
  capture$input_ids <- character()
  capture$messages <- character()

  visible <- withVisible(
    registerMetabNormItsdSelectionTracking(
      normData = norm_data,
      input = list(),
      walkFn = function(.x, .f) lapply(.x, .f),
      selectionGetter = function(input, inputId) {
        capture$input_ids <- c(capture$input_ids, inputId)
        NULL
      },
      observeEventFn = function(eventExpr, handlerExpr, ignoreNULL = TRUE) {
        event_env <- parent.frame()
        event_expr <- substitute(eventExpr)
        handler_env <- parent.frame()
        handler_expr <- substitute(handlerExpr)

        capture$event <- list(
          value = eval(event_expr, envir = event_env),
          ignoreNULL = ignoreNULL
        )

        eval(handler_expr, envir = handler_env)
        invisible(NULL)
      },
      logInfoFn = function(message) {
        capture$messages <- c(capture$messages, message)
      }
    )
  )

  expect_false(visible$visible)
  expect_null(visible$value)
  expect_identical(capture$input_ids, rep("itsd_table_lcms_pos_rows_selected", 2L))
  expect_identical(capture$event, list(value = NULL, ignoreNULL = FALSE))
  expect_length(norm_data$itsd_selections, 0L)
  expect_false("LCMS Pos" %in% names(norm_data$itsd_selections))
  expect_identical(
    capture$messages,
    "ITSD selection updated for LCMS Pos : 0 features selected"
  )
})

test_that("metabolomics normalization module ITSD selection tracking observer shell preserves req and registration handoff", {
  norm_data <- new.env(parent = emptyenv())
  norm_data$assay_names <- c("LCMS Pos", "GC-MS/Neg")

  capture <- new.env(parent = emptyenv())
  capture$req <- NULL
  capture$registration <- NULL

  visible <- withVisible(
    runMetabNormItsdSelectionTrackingObserverShell(
      normData = norm_data,
      input = "input-bag",
      reqFn = function(value) {
        capture$req <- value
        invisible(value)
      },
      registerItsdSelectionTrackingFn = function(normData, input) {
        capture$registration <- list(normData = normData, input = input)
        invisible("registered")
      }
    )
  )

  expect_false(visible$visible)
  expect_null(visible$value)
  expect_identical(capture$req, c("LCMS Pos", "GC-MS/Neg"))
  expect_identical(
    capture$registration,
    list(normData = norm_data, input = "input-bag")
  )
})

test_that("metabolomics normalization module ITSD selection tracking observer shell preserves req gate exit", {
  norm_data <- new.env(parent = emptyenv())
  norm_data$assay_names <- NULL

  capture <- new.env(parent = emptyenv())
  capture$req <- NULL
  capture$registered <- FALSE

  expect_error(
    runMetabNormItsdSelectionTrackingObserverShell(
      normData = norm_data,
      input = "input-bag",
      reqFn = function(value) {
        capture$req <- value
        stop("req gate tripped")
      },
      registerItsdSelectionTrackingFn = function(...) {
        capture$registered <- TRUE
        invisible(NULL)
      }
    ),
    "req gate tripped"
  )

  expect_null(capture$req)
  expect_false(capture$registered)
})

test_that("metabolomics normalization module ITSD selection table observer shell preserves per-assay bindings", {
  norm_data <- new.env(parent = emptyenv())
  norm_data$assay_names <- c("LCMS Pos", "GC-MS/Neg")

  workflow_data <- new.env(parent = emptyenv())
  workflow_data$state_manager <- "state-manager"

  output <- new.env(parent = emptyenv())

  capture <- new.env(parent = emptyenv())
  capture$req <- list()
  capture$render_calls <- list()

  current_s4 <- methods::new(
    "MetaboliteAssayData",
    metabolite_data = list(
      "LCMS Pos" = data.frame(feature = "F1"),
      "GC-MS/Neg" = data.frame(feature = "F2")
    ),
    metabolite_id_column = "feature_id",
    annotation_id_column = "annotation_id"
  )

  visible <- withVisible(
    runMetabNormItsdSelectionTableObserverShell(
      normData = norm_data,
      workflowData = workflow_data,
      output = output,
      reqFn = function(value) {
        capture$req <- c(capture$req, list(value))
        invisible(value)
      },
      getStateFn = function(stateManager) {
        capture$state_manager <- stateManager
        current_s4
      },
      walkFn = function(.x, .f) lapply(.x, .f),
      renderItsdSelectionTableFn = function(assayName, currentS4) {
        capture$render_calls[[length(capture$render_calls) + 1L]] <- list(
          assayName = assayName,
          currentS4 = currentS4
        )
        paste("rendered", assayName)
      }
    )
  )

  expect_false(visible$visible)
  expect_null(visible$value)
  expect_identical(capture$req, list(c("LCMS Pos", "GC-MS/Neg"), "state-manager"))
  expect_identical(capture$state_manager, "state-manager")
  expect_identical(
    capture$render_calls,
    list(
      list(assayName = "LCMS Pos", currentS4 = current_s4),
      list(assayName = "GC-MS/Neg", currentS4 = current_s4)
    )
  )
  expect_identical(output[["itsd_table_lcms_pos"]], "rendered LCMS Pos")
  expect_identical(output[["itsd_table_gc_ms_neg"]], "rendered GC-MS/Neg")
})

test_that("metabolomics normalization module ITSD selection table observer shell preserves getter-error and invalid-state exits", {
  norm_data <- new.env(parent = emptyenv())
  norm_data$assay_names <- "LCMS Pos"

  workflow_data <- new.env(parent = emptyenv())
  workflow_data$state_manager <- "state-manager"

  error_output <- new.env(parent = emptyenv())
  error_capture <- new.env(parent = emptyenv())

  error_visible <- withVisible(
    runMetabNormItsdSelectionTableObserverShell(
      normData = norm_data,
      workflowData = workflow_data,
      output = error_output,
      walkFn = function(.x, .f) lapply(.x, .f),
      renderItsdSelectionTableFn = function(...) {
        error_capture$render_called <- TRUE
        stop("renderItsdSelectionTableFn should not run after getState error")
      },
      getStateFn = function(...) {
        stop("state manager unavailable")
      }
    )
  )

  expect_false(error_visible$visible)
  expect_null(error_visible$value)
  expect_false(isTRUE(error_capture$render_called))
  expect_identical(ls(error_output), character(0))

  invalid_output <- new.env(parent = emptyenv())
  invalid_capture <- new.env(parent = emptyenv())

  invalid_visible <- withVisible(
    runMetabNormItsdSelectionTableObserverShell(
      normData = norm_data,
      workflowData = workflow_data,
      output = invalid_output,
      walkFn = function(.x, .f) lapply(.x, .f),
      getStateFn = function(...) list(not = "an-s4"),
      renderItsdSelectionTableFn = function(...) {
        invalid_capture$render_called <- TRUE
        stop("renderItsdSelectionTableFn should not run for invalid state")
      }
    )
  )

  expect_false(invalid_visible$visible)
  expect_null(invalid_visible$value)
  expect_false(isTRUE(invalid_capture$render_called))
  expect_identical(ls(invalid_output), character(0))
})

test_that("metabolomics normalization module RUV QC render seam preserves assay plot, summary, and table UI", {
  norm_data <- new.env(parent = emptyenv())
  norm_data$assay_names <- c("LCMS Pos", "GC-MS/Neg")

  capture <- new.env(parent = emptyenv())
  capture$req <- list()
  capture$plot_ids <- character()
  capture$summary_ids <- character()
  capture$table_ids <- character()

  render_fn <- renderMetabNormRuvQcUi(
    normData = norm_data,
    ns = function(id) paste0("module-", id),
    renderUIFn = function(expr) {
      render_env <- parent.frame()
      render_expr <- substitute(expr)

      function() {
        eval(render_expr, envir = render_env)
      }
    },
    reqFn = function(value) {
      capture$req <- c(capture$req, list(value))
      invisible(value)
    },
    mapFn = function(.x, .f) lapply(.x, .f),
    tagListFn = function(...) list(tag = "tagList", children = list(...)),
    fluidRowFn = function(...) list(tag = "fluidRow", children = list(...)),
    columnFn = function(width, ...) list(tag = "column", width = width, children = list(...)),
    h5Fn = function(text, style = NULL) list(tag = "h5", text = text, style = style),
    h6Fn = function(text) list(tag = "h6", text = text),
    jquiResizableFn = function(widget) list(tag = "jqui_resizable", child = widget),
    plotOutputFn = function(outputId, height = NULL) {
      capture$plot_ids <- c(capture$plot_ids, outputId)
      list(tag = "plotOutput", outputId = outputId, height = height)
    },
    wellPanelFn = function(...) list(tag = "wellPanel", children = list(...)),
    verbatimTextOutputFn = function(outputId) {
      capture$summary_ids <- c(capture$summary_ids, outputId)
      list(tag = "verbatimTextOutput", outputId = outputId)
    },
    brFn = function() list(tag = "br"),
    dataTableOutputFn = function(outputId) {
      capture$table_ids <- c(capture$table_ids, outputId)
      list(tag = "dataTableOutput", outputId = outputId)
    },
    hrFn = function() list(tag = "hr")
  )

  expect_true(is.function(render_fn))

  rendered_ui <- render_fn()

  expect_identical(capture$req, list(c("LCMS Pos", "GC-MS/Neg")))
  expect_identical(
    capture$plot_ids,
    c("module-cancor_plot_lcms_pos", "module-cancor_plot_gc_ms_neg")
  )
  expect_identical(
    capture$summary_ids,
    c("module-ruv_summary_lcms_pos", "module-ruv_summary_gc_ms_neg")
  )
  expect_identical(
    capture$table_ids,
    c("module-ruv_table_lcms_pos", "module-ruv_table_gc_ms_neg")
  )
  expect_identical(rendered_ui$tag, "tagList")
  expect_length(rendered_ui$children, 2L)
  expect_identical(
    rendered_ui$children[[1]]$children[[1]]$children[[1]]$children[[1]],
    list(
      tag = "h5",
      text = "Assay: LCMS Pos",
      style = "border-bottom: 1px solid #ddd; padding-bottom: 5px;"
    )
  )
  expect_identical(
    rendered_ui$children[[1]]$children[[2]]$children[[1]]$children[[1]]$child,
    list(
      tag = "plotOutput",
      outputId = "module-cancor_plot_lcms_pos",
      height = "400px"
    )
  )
  expect_identical(
    rendered_ui$children[[1]]$children[[2]]$children[[2]]$children[[1]]$children[[1]],
    list(tag = "h6", text = "Optimization Summary")
  )
  expect_identical(
    rendered_ui$children[[1]]$children[[2]]$children[[2]]$children[[1]]$children[[2]],
    list(tag = "verbatimTextOutput", outputId = "module-ruv_summary_lcms_pos")
  )
  expect_identical(
    rendered_ui$children[[1]]$children[[2]]$children[[2]]$children[[1]]$children[[4]],
    list(tag = "h6", text = "Results Table")
  )
  expect_identical(
    rendered_ui$children[[2]]$children[[2]]$children[[2]]$children[[1]]$children[[5]],
    list(tag = "dataTableOutput", outputId = "module-ruv_table_gc_ms_neg")
  )
})

test_that("metabolomics normalization module RUV cancor plot render helper preserves success and fallback plot selection", {
  render_plot_capture <- function(expr) {
    render_env <- parent.frame()
    render_expr <- substitute(expr)

    function() {
      eval(render_expr, envir = render_env)
    }
  }

  success_plot <- ggplot2::ggplot(data.frame(x = 1, y = 1), ggplot2::aes(x, y)) +
    ggplot2::geom_point()
  success_capture <- new.env(parent = emptyenv())

  success_render <- renderMetabNormRuvCancorPlot(
    assayName = "LCMS Pos",
    normData = list(
      ruv_optimization_results = list(
        "LCMS Pos" = list(success = TRUE, cancor_plot = success_plot)
      )
    ),
    renderPlotFn = render_plot_capture,
    buildFallbackPlotFn = function() {
      success_capture$fallback_called <- TRUE
      "fallback"
    }
  )

  expect_true(is.function(success_render))
  expect_identical(success_render(), success_plot)
  expect_false(isTRUE(success_capture$fallback_called))

  fallback_capture <- new.env(parent = emptyenv())
  fallback_capture$fallback_called <- FALSE

  fallback_render <- renderMetabNormRuvCancorPlot(
    assayName = "GC-MS/Neg",
    normData = list(
      ruv_optimization_results = list(
        "GC-MS/Neg" = list(success = FALSE, error = "optimization failed")
      )
    ),
    renderPlotFn = render_plot_capture,
    buildFallbackPlotFn = function() {
      fallback_capture$fallback_called <- TRUE
      list(kind = "fallback-plot")
    }
  )

  expect_true(is.function(fallback_render))
  expect_identical(fallback_render(), list(kind = "fallback-plot"))
  expect_true(fallback_capture$fallback_called)
})

test_that("metabolomics normalization module RUV optimization summary render helper preserves success, failure, and pending text", {
  render_text_capture <- function(expr) {
    eval.parent(substitute(expr))
  }

  expect_identical(
    renderMetabNormRuvOptimizationSummary(
      assayName = "LCMS Pos",
      normData = list(
        ruv_optimization_results = list(
          "LCMS Pos" = list(
            success = TRUE,
            best_k = 2,
            best_percentage = 0.35,
            separation_score = 0.98765,
            control_genes_index = c(TRUE, FALSE, TRUE, NA)
          )
        )
      ),
      renderTextFn = render_text_capture
    ),
    paste0(
      "Best k: 2\n",
      "Best %: 0.35\n",
      "Separation: 0.9876\n",
      "Controls: 2"
    )
  )

  expect_identical(
    renderMetabNormRuvOptimizationSummary(
      assayName = "GC-MS/Neg",
      normData = list(
        ruv_optimization_results = list(
          "GC-MS/Neg" = list(success = FALSE, error = "optimization failed")
        )
      ),
      renderTextFn = render_text_capture
    ),
    "Failed: optimization failed"
  )

  expect_identical(
    renderMetabNormRuvOptimizationSummary(
      assayName = "Missing",
      normData = list(ruv_optimization_results = list()),
      renderTextFn = render_text_capture
    ),
    "Not yet computed"
  )
})

test_that("metabolomics normalization module RUV results table render helper preserves success and null exits", {
  render_data_table_capture <- function(expr) {
    render_env <- parent.frame()
    render_expr <- substitute(expr)

    function() {
      eval(render_expr, envir = render_env)
    }
  }

  success_capture <- new.env(parent = emptyenv())
  optimization_results <- data.frame(
    k = c(1, 2),
    score = c(0.12, 0.34),
    stringsAsFactors = FALSE
  )

  success_render <- renderMetabNormRuvResultsTable(
    assayName = "LCMS Pos",
    normData = list(
      ruv_optimization_results = list(
        "LCMS Pos" = list(
          success = TRUE,
          optimization_results = optimization_results
        )
      )
    ),
    renderDataTableFn = render_data_table_capture,
    datatableFn = function(data, options = NULL, rownames = NULL) {
      success_capture$args <- list(
        data = data,
        options = options,
        rownames = rownames
      )
      list(tag = "datatable", data = data, options = options, rownames = rownames)
    }
  )

  expect_true(is.function(success_render))
  expect_identical(
    success_render(),
    list(
      tag = "datatable",
      data = optimization_results,
      options = list(pageLength = 5, dom = "t"),
      rownames = FALSE
    )
  )
  expect_identical(success_capture$args$data, optimization_results)
  expect_identical(success_capture$args$options, list(pageLength = 5, dom = "t"))
  expect_false(success_capture$args$rownames)

  failure_capture <- new.env(parent = emptyenv())
  failure_render <- renderMetabNormRuvResultsTable(
    assayName = "GC-MS/Neg",
    normData = list(
      ruv_optimization_results = list(
        "GC-MS/Neg" = list(success = FALSE, error = "optimization failed")
      )
    ),
    renderDataTableFn = render_data_table_capture,
    datatableFn = function(...) {
      failure_capture$called <- TRUE
      "unexpected"
    }
  )

  expect_true(is.function(failure_render))
  expect_null(failure_render())
  expect_false(isTRUE(failure_capture$called))

  pending_render <- renderMetabNormRuvResultsTable(
    assayName = "Missing",
    normData = list(ruv_optimization_results = list()),
    renderDataTableFn = render_data_table_capture,
    datatableFn = function(...) {
      stop("datatableFn should not run for missing optimization results")
    }
  )

  expect_true(is.function(pending_render))
  expect_null(pending_render())
})

test_that("metabolomics normalization module RUV binding observer shell preserves per-assay plot, summary, and table bindings", {
  norm_data <- new.env(parent = emptyenv())
  norm_data$assay_names <- c("LCMS Pos", "GC-MS/Neg")
  norm_data$ruv_optimization_results <- list(
    "LCMS Pos" = list(success = TRUE),
    "GC-MS/Neg" = list(success = TRUE)
  )

  output <- new.env(parent = emptyenv())
  capture <- new.env(parent = emptyenv())
  capture$req <- list()
  capture$plot_calls <- list()
  capture$summary_calls <- list()
  capture$table_calls <- list()

  visible <- withVisible(
    runMetabNormRuvBindingObserverShell(
      normData = norm_data,
      output = output,
      reqFn = function(value) {
        capture$req <- c(capture$req, list(value))
        invisible(value)
      },
      walkFn = function(.x, .f) lapply(.x, .f),
      renderRuvCancorPlotFn = function(assayName, normData) {
        capture$plot_calls[[length(capture$plot_calls) + 1L]] <- list(
          assayName = assayName,
          normData = normData
        )
        paste("plot", assayName)
      },
      renderRuvOptimizationSummaryFn = function(assayName, normData) {
        capture$summary_calls[[length(capture$summary_calls) + 1L]] <- list(
          assayName = assayName,
          normData = normData
        )
        paste("summary", assayName)
      },
      renderRuvResultsTableFn = function(assayName, normData) {
        capture$table_calls[[length(capture$table_calls) + 1L]] <- list(
          assayName = assayName,
          normData = normData
        )
        paste("table", assayName)
      }
    )
  )

  expect_false(visible$visible)
  expect_null(visible$value)
  expect_identical(capture$req, list(c("LCMS Pos", "GC-MS/Neg"), TRUE))
  expect_identical(
    capture$plot_calls,
    list(
      list(assayName = "LCMS Pos", normData = norm_data),
      list(assayName = "GC-MS/Neg", normData = norm_data)
    )
  )
  expect_identical(capture$summary_calls, capture$plot_calls)
  expect_identical(capture$table_calls, capture$plot_calls)
  expect_identical(output[["cancor_plot_lcms_pos"]], "plot LCMS Pos")
  expect_identical(output[["ruv_summary_lcms_pos"]], "summary LCMS Pos")
  expect_identical(output[["ruv_table_lcms_pos"]], "table LCMS Pos")
  expect_identical(output[["cancor_plot_gc_ms_neg"]], "plot GC-MS/Neg")
  expect_identical(output[["ruv_summary_gc_ms_neg"]], "summary GC-MS/Neg")
  expect_identical(output[["ruv_table_gc_ms_neg"]], "table GC-MS/Neg")
})

test_that("metabolomics normalization module RUV binding observer shell preserves empty-results req gate", {
  norm_data <- new.env(parent = emptyenv())
  norm_data$assay_names <- "LCMS Pos"
  norm_data$ruv_optimization_results <- list()

  output <- new.env(parent = emptyenv())
  capture <- new.env(parent = emptyenv())
  capture$req <- list()

  expect_error(
    runMetabNormRuvBindingObserverShell(
      normData = norm_data,
      output = output,
      reqFn = function(value) {
        capture$req <- c(capture$req, list(value))
        if (isFALSE(value) || is.null(value) || length(value) == 0L) {
          stop("req gate tripped")
        }
        invisible(value)
      },
      walkFn = function(.x, .f) lapply(.x, .f),
      renderRuvCancorPlotFn = function(...) {
        stop("renderRuvCancorPlotFn should not run when req gate fails")
      },
      renderRuvOptimizationSummaryFn = function(...) {
        stop("renderRuvOptimizationSummaryFn should not run when req gate fails")
      },
      renderRuvResultsTableFn = function(...) {
        stop("renderRuvResultsTableFn should not run when req gate fails")
      }
    ),
    "req gate tripped"
  )

  expect_identical(capture$req, list("LCMS Pos", FALSE))
  expect_identical(ls(output), character(0))
})

test_that("metabolomics normalization module assay label binding shell preserves static output fan-out", {
  output <- new.env(parent = emptyenv())
  capture <- new.env(parent = emptyenv())
  capture$calls <- list()
  get_assay_names <- function() c("Plasma", "Urine")

  visible <- withVisible(
    runMetabNormAssayLabelBindingShell(
      output = output,
      getAssayNamesFn = get_assay_names,
      renderAssayLabelFn = function(assaySlot, getAssayNamesFn) {
        capture$calls[[length(capture$calls) + 1L]] <- list(
          assaySlot = assaySlot,
          getAssayNamesFn = getAssayNamesFn
        )
        paste("label", assaySlot)
      }
    )
  )

  expect_false(visible$visible)
  expect_null(visible$value)
  expect_identical(
    lapply(capture$calls, `[[`, "assaySlot"),
    list(1, 2, 1, 2, 1, 2, 1, 2)
  )
  expect_true(all(vapply(
    capture$calls,
    function(call) identical(call$getAssayNamesFn, get_assay_names),
    logical(1)
  )))
  expect_identical(output$assay1_label_pca, "label 1")
  expect_identical(output$assay2_label_pca, "label 2")
  expect_identical(output$assay1_label_density, "label 1")
  expect_identical(output$assay2_label_density, "label 2")
  expect_identical(output$assay1_label_rle, "label 1")
  expect_identical(output$assay2_label_rle, "label 2")
  expect_identical(output$assay1_label_correlation, "label 1")
  expect_identical(output$assay2_label_correlation, "label 2")
})

test_that("metabolomics normalization module QC image binding shell preserves static output fan-out", {
  output <- new.env(parent = emptyenv())
  capture <- new.env(parent = emptyenv())
  capture$calls <- list()

  expected_bindings <- list(
    list(outputId = "pca_post_filter_assay1", assaySlot = 1, plotType = "pca", stagePrefix = "pre_norm"),
    list(outputId = "pca_post_norm_assay1", assaySlot = 1, plotType = "pca", stagePrefix = "post_norm"),
    list(outputId = "pca_ruv_corrected_assay1", assaySlot = 1, plotType = "pca", stagePrefix = "ruv_corrected"),
    list(outputId = "pca_post_filter_assay2", assaySlot = 2, plotType = "pca", stagePrefix = "pre_norm"),
    list(outputId = "pca_post_norm_assay2", assaySlot = 2, plotType = "pca", stagePrefix = "post_norm"),
    list(outputId = "pca_ruv_corrected_assay2", assaySlot = 2, plotType = "pca", stagePrefix = "ruv_corrected"),
    list(outputId = "density_post_filter_assay1", assaySlot = 1, plotType = "density", stagePrefix = "pre_norm"),
    list(outputId = "density_post_norm_assay1", assaySlot = 1, plotType = "density", stagePrefix = "post_norm"),
    list(outputId = "density_ruv_corrected_assay1", assaySlot = 1, plotType = "density", stagePrefix = "ruv_corrected"),
    list(outputId = "density_post_filter_assay2", assaySlot = 2, plotType = "density", stagePrefix = "pre_norm"),
    list(outputId = "density_post_norm_assay2", assaySlot = 2, plotType = "density", stagePrefix = "post_norm"),
    list(outputId = "density_ruv_corrected_assay2", assaySlot = 2, plotType = "density", stagePrefix = "ruv_corrected"),
    list(outputId = "rle_post_filter_assay1", assaySlot = 1, plotType = "rle", stagePrefix = "pre_norm"),
    list(outputId = "rle_post_norm_assay1", assaySlot = 1, plotType = "rle", stagePrefix = "post_norm"),
    list(outputId = "rle_ruv_corrected_assay1", assaySlot = 1, plotType = "rle", stagePrefix = "ruv_corrected"),
    list(outputId = "rle_post_filter_assay2", assaySlot = 2, plotType = "rle", stagePrefix = "pre_norm"),
    list(outputId = "rle_post_norm_assay2", assaySlot = 2, plotType = "rle", stagePrefix = "post_norm"),
    list(outputId = "rle_ruv_corrected_assay2", assaySlot = 2, plotType = "rle", stagePrefix = "ruv_corrected"),
    list(outputId = "correlation_post_filter_assay1", assaySlot = 1, plotType = "correlation", stagePrefix = "pre_norm"),
    list(outputId = "correlation_post_norm_assay1", assaySlot = 1, plotType = "correlation", stagePrefix = "post_norm"),
    list(outputId = "correlation_ruv_corrected_assay1", assaySlot = 1, plotType = "correlation", stagePrefix = "ruv_corrected"),
    list(outputId = "correlation_post_filter_assay2", assaySlot = 2, plotType = "correlation", stagePrefix = "pre_norm"),
    list(outputId = "correlation_post_norm_assay2", assaySlot = 2, plotType = "correlation", stagePrefix = "post_norm"),
    list(outputId = "correlation_ruv_corrected_assay2", assaySlot = 2, plotType = "correlation", stagePrefix = "ruv_corrected")
  )

  norm_data <- structure(list(assay_names = c("Plasma", "Urine")), class = "MetabNormQcShellMock")
  qc_dir <- "/tmp/metab-qc"

  visible <- withVisible(
    runMetabNormQcImageBindingShell(
      output = output,
      normData = norm_data,
      qcDir = qc_dir,
      renderQcImageFn = function(assaySlot, plotType, stagePrefix, normData, qcDir) {
        capture$calls[[length(capture$calls) + 1L]] <- list(
          assaySlot = assaySlot,
          plotType = plotType,
          stagePrefix = stagePrefix,
          normData = normData,
          qcDir = qcDir
        )
        paste(plotType, stagePrefix, assaySlot, sep = "::")
      }
    )
  )

  expect_false(visible$visible)
  expect_null(visible$value)
  expect_length(capture$calls, length(expected_bindings))

  for (index in seq_along(expected_bindings)) {
    expected <- expected_bindings[[index]]
    actual <- capture$calls[[index]]

    expect_identical(actual$assaySlot, expected$assaySlot)
    expect_identical(actual$plotType, expected$plotType)
    expect_identical(actual$stagePrefix, expected$stagePrefix)
    expect_identical(actual$normData, norm_data)
    expect_identical(actual$qcDir, qc_dir)
    expect_identical(
      output[[expected$outputId]],
      paste(expected$plotType, expected$stagePrefix, expected$assaySlot, sep = "::")
    )
  }
})

test_that("metabolomics normalization module ITSD progress/apply shell preserves progress, manual-ID resolution, and step handoff", {
  capture <- new.env(parent = emptyenv())
  capture$log <- character()

  current_s4 <- structure(list(stage = "pre_itsd"), class = "MetabNormItsdShellMock")
  updated_s4 <- structure(list(stage = "post_itsd"), class = "MetabNormItsdShellMock")
  workflow_data <- list(state_manager = "state-manager", config_list = list(method = "median"))
  norm_data <- new.env(parent = emptyenv())
  itsd_feature_ids <- list(Plasma = c("M1", "M3"))

  visible <- withVisible(
    runMetabNormItsdProgressApplyShell(
      currentS4 = current_s4,
      totalSteps = 6,
      applyItsd = TRUE,
      itsdAggregation = "median",
      itsdSelections = list(Plasma = c(1L, 3L)),
      workflowData = workflow_data,
      normData = norm_data,
      addLogFn = function(message) {
        capture$log <- c(capture$log, message)
      },
      incProgressFn = function(amount, detail = NULL) {
        capture$progress <- list(amount = amount, detail = detail)
      },
      resolveManualFeatureIdsFn = function(currentS4, itsdSelections, addLogFn) {
        capture$resolve <- list(
          currentS4 = currentS4,
          itsdSelections = itsdSelections,
          addLogFn = addLogFn
        )
        itsd_feature_ids
      },
      runItsdStepFn = function(currentS4, itsdAggregation, itsdFeatureIds, workflowData, normData, addLogFn) {
        capture$run_step <- list(
          currentS4 = currentS4,
          itsdAggregation = itsdAggregation,
          itsdFeatureIds = itsdFeatureIds,
          workflowData = workflowData,
          normData = normData,
          addLogFn = addLogFn
        )
        list(currentS4 = updated_s4, stateName = "metab_itsd_norm")
      }
    )
  )

  expect_false(visible$visible)
  expect_identical(
    capture$progress,
    list(amount = 1 / 6, detail = "Applying ITSD normalization...")
  )
  expect_identical(
    capture$log,
    "Applying ITSD normalization (aggregation: median )"
  )
  expect_identical(
    capture$resolve$currentS4,
    current_s4
  )
  expect_identical(
    capture$resolve$itsdSelections,
    list(Plasma = c(1L, 3L))
  )
  expect_true(is.function(capture$resolve$addLogFn))
  expect_identical(
    capture$run_step,
    list(
      currentS4 = current_s4,
      itsdAggregation = "median",
      itsdFeatureIds = itsd_feature_ids,
      workflowData = workflow_data,
      normData = norm_data,
      addLogFn = capture$resolve$addLogFn
    )
  )
  expect_identical(visible$value$currentS4, updated_s4)
  expect_true(visible$value$applied)
  expect_identical(visible$value$progressDetail, "Applying ITSD normalization...")
  expect_identical(
    visible$value$applyLogEntry,
    "Applying ITSD normalization (aggregation: median )"
  )
  expect_null(visible$value$skippedLogEntry)
  expect_identical(visible$value$itsdFeatureIds, itsd_feature_ids)
  expect_identical(
    visible$value$itsdState,
    list(currentS4 = updated_s4, stateName = "metab_itsd_norm")
  )
})

test_that("metabolomics normalization module ITSD progress/apply shell preserves skipped branch", {
  capture <- new.env(parent = emptyenv())
  capture$log <- character()

  current_s4 <- structure(list(stage = "pre_itsd"), class = "MetabNormItsdShellSkipMock")
  workflow_data <- list(state_manager = "state-manager")
  norm_data <- new.env(parent = emptyenv())

  visible <- withVisible(
    runMetabNormItsdProgressApplyShell(
      currentS4 = current_s4,
      totalSteps = 6,
      applyItsd = FALSE,
      itsdAggregation = "median",
      itsdSelections = list(Plasma = 1L),
      workflowData = workflow_data,
      normData = norm_data,
      addLogFn = function(message) {
        capture$log <- c(capture$log, message)
      },
      incProgressFn = function(amount, detail = NULL) {
        capture$progress <- list(amount = amount, detail = detail)
      },
      resolveManualFeatureIdsFn = function(...) {
        stop("resolveManualFeatureIdsFn should not be called")
      },
      runItsdStepFn = function(...) {
        stop("runItsdStepFn should not be called")
      }
    )
  )

  expect_false(visible$visible)
  expect_identical(
    capture$progress,
    list(amount = 1 / 6, detail = "Applying ITSD normalization...")
  )
  expect_identical(capture$log, "ITSD normalization skipped")
  expect_identical(visible$value$currentS4, current_s4)
  expect_false(visible$value$applied)
  expect_identical(visible$value$progressDetail, "Applying ITSD normalization...")
  expect_null(visible$value$applyLogEntry)
  expect_identical(visible$value$skippedLogEntry, "ITSD normalization skipped")
  expect_null(visible$value$itsdFeatureIds)
  expect_null(visible$value$itsdState)
})

test_that("metabolomics normalization module ITSD step helper preserves apply/saveState shell", {
  capture <- new.env(parent = emptyenv())
  capture$log <- character()

  current_s4 <- structure(list(stage = "pre_itsd"), class = "MetabNormItsdStepMock")
  updated_s4 <- structure(list(stage = "post_itsd"), class = "MetabNormItsdStepMock")
  workflow_data <- new.env(parent = emptyenv())
  norm_data <- new.env(parent = emptyenv())

  workflow_data$config_list <- list(method = "median")
  workflow_data$state_manager <- list(
    saveState = function(state_name, s4_data_object, config_object, description) {
      capture$save_state <- list(
        state_name = state_name,
        s4_data_object = s4_data_object,
        config_object = config_object,
        description = description
      )
    }
  )

  visible <- withVisible(
    runMetabNormItsdNormalizationStep(
      currentS4 = current_s4,
      itsdAggregation = "median",
      itsdFeatureIds = list(Plasma = c("M1", "M3")),
      workflowData = workflow_data,
      normData = norm_data,
      addLogFn = function(message) {
        capture$log <- c(capture$log, message)
      },
      normaliseUntransformedDataFn = function(theObject, method, itsd_aggregation, itsd_feature_ids) {
        capture$normalise <- list(
          theObject = theObject,
          method = method,
          itsd_aggregation = itsd_aggregation,
          itsd_feature_ids = itsd_feature_ids
        )
        updated_s4
      }
    )
  )

  expect_false(visible$visible)
  expect_identical(
    capture$normalise,
    list(
      theObject = current_s4,
      method = "ITSD",
      itsd_aggregation = "median",
      itsd_feature_ids = list(Plasma = c("M1", "M3"))
    )
  )
  expect_identical(norm_data$post_itsd_obj, updated_s4)
  expect_identical(
    capture$save_state,
    list(
      state_name = "metab_itsd_norm",
      s4_data_object = updated_s4,
      config_object = workflow_data$config_list,
      description = "ITSD normalization (aggregation: median )"
    )
  )
  expect_identical(capture$log, "ITSD normalization complete")
  expect_identical(
    visible$value,
    list(
      currentS4 = updated_s4,
      stateName = "metab_itsd_norm",
      description = "ITSD normalization (aggregation: median )",
      logEntry = "ITSD normalization complete"
    )
  )
})

test_that("metabolomics normalization module log2 step helper preserves apply/saveState shell", {
  capture <- new.env(parent = emptyenv())
  capture$log <- character()

  current_s4 <- structure(list(stage = "pre_log2"), class = "MetabNormLog2StepMock")
  updated_s4 <- structure(list(stage = "post_log2"), class = "MetabNormLog2StepMock")
  workflow_data <- new.env(parent = emptyenv())
  norm_data <- new.env(parent = emptyenv())

  workflow_data$config_list <- list(offset = 0.5)
  workflow_data$state_manager <- list(
    saveState = function(state_name, s4_data_object, config_object, description) {
      capture$save_state <- list(
        state_name = state_name,
        s4_data_object = s4_data_object,
        config_object = config_object,
        description = description
      )
    }
  )

  visible <- withVisible(
    runMetabNormLog2TransformationStep(
      currentS4 = current_s4,
      logOffset = 0.5,
      workflowData = workflow_data,
      normData = norm_data,
      addLogFn = function(message) {
        capture$log <- c(capture$log, message)
      },
      logTransformAssaysFn = function(theObject, offset) {
        capture$transform <- list(
          theObject = theObject,
          offset = offset
        )
        updated_s4
      }
    )
  )

  expect_false(visible$visible)
  expect_identical(
    capture$transform,
    list(
      theObject = current_s4,
      offset = 0.5
    )
  )
  expect_identical(norm_data$post_log2_obj, updated_s4)
  expect_identical(
    capture$save_state,
    list(
      state_name = "metab_log2",
      s4_data_object = updated_s4,
      config_object = workflow_data$config_list,
      description = "Log2 transformation (offset: 0.5 )"
    )
  )
  expect_identical(capture$log, "Log2 transformation complete")
  expect_identical(
    visible$value,
    list(
      currentS4 = updated_s4,
      stateName = "metab_log2",
      description = "Log2 transformation (offset: 0.5 )",
      logEntry = "Log2 transformation complete"
    )
  )
})

test_that("metabolomics normalization module log2 progress/apply shell preserves progress and step handoff", {
  capture <- new.env(parent = emptyenv())
  capture$log <- character()

  current_s4 <- structure(list(stage = "pre_log2"), class = "MetabNormLog2ShellMock")
  updated_s4 <- structure(list(stage = "post_log2"), class = "MetabNormLog2ShellMock")
  workflow_data <- list(state_manager = "state-manager", config_list = list(offset = 0.5))
  norm_data <- new.env(parent = emptyenv())

  visible <- withVisible(
    runMetabNormLog2ProgressApplyShell(
      currentS4 = current_s4,
      totalSteps = 6,
      logOffset = 0.5,
      workflowData = workflow_data,
      normData = norm_data,
      addLogFn = function(message) {
        capture$log <- c(capture$log, message)
      },
      incProgressFn = function(amount, detail = NULL) {
        capture$progress <- list(amount = amount, detail = detail)
      },
      runLog2StepFn = function(currentS4, logOffset, workflowData, normData, addLogFn) {
        capture$run_step <- list(
          currentS4 = currentS4,
          logOffset = logOffset,
          workflowData = workflowData,
          normData = normData,
          addLogFn = addLogFn
        )
        list(currentS4 = updated_s4, stateName = "metab_log2")
      }
    )
  )

  expect_false(visible$visible)
  expect_identical(
    capture$progress,
    list(amount = 1 / 6, detail = "Applying log2 transformation...")
  )
  expect_identical(
    capture$log,
    "Applying log2 transformation (offset: 0.5 )"
  )
  expect_identical(
    capture$run_step,
    list(
      currentS4 = current_s4,
      logOffset = 0.5,
      workflowData = workflow_data,
      normData = norm_data,
      addLogFn = capture$run_step$addLogFn
    )
  )
  expect_true(is.function(capture$run_step$addLogFn))
  expect_identical(visible$value$currentS4, updated_s4)
  expect_identical(visible$value$progressDetail, "Applying log2 transformation...")
  expect_identical(
    visible$value$applyLogEntry,
    "Applying log2 transformation (offset: 0.5 )"
  )
  expect_identical(
    visible$value$log2State,
    list(currentS4 = updated_s4, stateName = "metab_log2")
  )
})

test_that("metabolomics normalization module between-sample progress/apply shell preserves progress and step handoff", {
  capture <- new.env(parent = emptyenv())
  capture$log <- character()

  current_s4 <- structure(list(stage = "post_log2"), class = "MetabNormBetweenSampleShellMock")
  updated_s4 <- structure(list(stage = "post_norm"), class = "MetabNormBetweenSampleShellMock")
  workflow_data <- list(state_manager = "state-manager", config_list = list(method = "median"))
  norm_data <- new.env(parent = emptyenv())

  visible <- withVisible(
    runMetabNormBetweenSampleProgressApplyShell(
      currentS4 = current_s4,
      totalSteps = 6,
      normMethod = "median",
      workflowData = workflow_data,
      normData = norm_data,
      addLogFn = function(message) {
        capture$log <- c(capture$log, message)
      },
      incProgressFn = function(amount, detail = NULL) {
        capture$progress <- list(amount = amount, detail = detail)
      },
      runBetweenSampleStepFn = function(currentS4, normMethod, workflowData, normData, addLogFn) {
        capture$run_step <- list(
          currentS4 = currentS4,
          normMethod = normMethod,
          workflowData = workflowData,
          normData = normData,
          addLogFn = addLogFn
        )
        list(currentS4 = updated_s4, stateName = "metab_normalized")
      }
    )
  )

  expect_false(visible$visible)
  expect_identical(
    capture$progress,
    list(amount = 1 / 6, detail = "Applying between-sample normalization...")
  )
  expect_identical(
    capture$log,
    "Applying between-sample normalization (method: median )"
  )
  expect_identical(
    capture$run_step,
    list(
      currentS4 = current_s4,
      normMethod = "median",
      workflowData = workflow_data,
      normData = norm_data,
      addLogFn = capture$run_step$addLogFn
    )
  )
  expect_true(is.function(capture$run_step$addLogFn))
  expect_identical(visible$value$currentS4, updated_s4)
  expect_identical(
    visible$value$progressDetail,
    "Applying between-sample normalization..."
  )
  expect_identical(
    visible$value$applyLogEntry,
    "Applying between-sample normalization (method: median )"
  )
  expect_identical(
    visible$value$betweenSampleState,
    list(currentS4 = updated_s4, stateName = "metab_normalized")
  )
})

test_that("metabolomics normalization module between-sample progress/apply shell preserves no-op logging branch", {
  capture <- new.env(parent = emptyenv())
  capture$log <- character()

  current_s4 <- structure(list(stage = "post_log2"), class = "MetabNormBetweenSampleShellMock")
  workflow_data <- list(state_manager = "state-manager", config_list = list(method = "none"))
  norm_data <- new.env(parent = emptyenv())

  visible <- withVisible(
    runMetabNormBetweenSampleProgressApplyShell(
      currentS4 = current_s4,
      totalSteps = 6,
      normMethod = "none",
      workflowData = workflow_data,
      normData = norm_data,
      addLogFn = function(message) {
        capture$log <- c(capture$log, message)
      },
      incProgressFn = function(amount, detail = NULL) {
        capture$progress <- list(amount = amount, detail = detail)
      },
      runBetweenSampleStepFn = function(currentS4, normMethod, workflowData, normData, addLogFn) {
        capture$run_step <- list(
          currentS4 = currentS4,
          normMethod = normMethod,
          workflowData = workflowData,
          normData = normData,
          addLogFn = addLogFn
        )
        list(currentS4 = currentS4, stateName = "metab_normalized")
      }
    )
  )

  expect_false(visible$visible)
  expect_identical(
    capture$progress,
    list(amount = 1 / 6, detail = "Applying between-sample normalization...")
  )
  expect_identical(capture$log, character())
  expect_identical(
    capture$run_step,
    list(
      currentS4 = current_s4,
      normMethod = "none",
      workflowData = workflow_data,
      normData = norm_data,
      addLogFn = capture$run_step$addLogFn
    )
  )
  expect_true(is.function(capture$run_step$addLogFn))
  expect_identical(visible$value$currentS4, current_s4)
  expect_identical(
    visible$value$progressDetail,
    "Applying between-sample normalization..."
  )
  expect_null(visible$value$applyLogEntry)
  expect_identical(
    visible$value$betweenSampleState,
    list(currentS4 = current_s4, stateName = "metab_normalized")
  )
})

test_that("metabolomics normalization module between-sample step helper preserves apply/saveState shell", {
  capture <- new.env(parent = emptyenv())
  capture$log <- character()

  current_s4 <- structure(list(stage = "post_log2"), class = "MetabNormBetweenSampleStepMock")
  updated_s4 <- structure(list(stage = "post_norm"), class = "MetabNormBetweenSampleStepMock")
  workflow_data <- new.env(parent = emptyenv())
  norm_data <- new.env(parent = emptyenv())

  workflow_data$config_list <- list(method = "median")
  workflow_data$state_manager <- list(
    saveState = function(state_name, s4_data_object, config_object, description) {
      capture$save_state <- list(
        state_name = state_name,
        s4_data_object = s4_data_object,
        config_object = config_object,
        description = description
      )
    }
  )

  visible <- withVisible(
    runMetabNormBetweenSampleNormalizationStep(
      currentS4 = current_s4,
      normMethod = "median",
      workflowData = workflow_data,
      normData = norm_data,
      addLogFn = function(message) {
        capture$log <- c(capture$log, message)
      },
      normaliseBetweenSamplesFn = function(theObject, normalisation_method) {
        capture$normalise <- list(
          theObject = theObject,
          normalisation_method = normalisation_method
        )
        updated_s4
      }
    )
  )

  expect_false(visible$visible)
  expect_identical(
    capture$normalise,
    list(
      theObject = current_s4,
      normalisation_method = "median"
    )
  )
  expect_identical(norm_data$post_norm_obj, updated_s4)
  expect_identical(
    capture$save_state,
    list(
      state_name = "metab_normalized",
      s4_data_object = updated_s4,
      config_object = workflow_data$config_list,
      description = "Between-sample normalization (method: median )"
    )
  )
  expect_identical(capture$log, "Between-sample normalization complete")
  expect_true(isTRUE(norm_data$normalization_complete))
  expect_identical(
    visible$value,
    list(
      currentS4 = updated_s4,
      stateName = "metab_normalized",
      description = "Between-sample normalization (method: median )",
      logEntry = "Between-sample normalization complete"
    )
  )
})

test_that("metabolomics normalization module between-sample step helper preserves no-op passthrough", {
  capture <- new.env(parent = emptyenv())
  capture$log <- character()
  capture$normalise_called <- FALSE

  current_s4 <- structure(list(stage = "post_log2"), class = "MetabNormBetweenSampleStepMock")
  workflow_data <- new.env(parent = emptyenv())
  norm_data <- new.env(parent = emptyenv())

  workflow_data$config_list <- list(method = "none")
  workflow_data$state_manager <- list(
    saveState = function(state_name, s4_data_object, config_object, description) {
      capture$save_state <- list(
        state_name = state_name,
        s4_data_object = s4_data_object,
        config_object = config_object,
        description = description
      )
    }
  )

  visible <- withVisible(
    runMetabNormBetweenSampleNormalizationStep(
      currentS4 = current_s4,
      normMethod = "none",
      workflowData = workflow_data,
      normData = norm_data,
      addLogFn = function(message) {
        capture$log <- c(capture$log, message)
      },
      normaliseBetweenSamplesFn = function(theObject, normalisation_method) {
        capture$normalise_called <- TRUE
        stop(sprintf(
          "normaliseBetweenSamplesFn should not run for method 'none' (got %s)",
          normalisation_method
        ))
      }
    )
  )

  expect_false(visible$visible)
  expect_false(capture$normalise_called)
  expect_identical(norm_data$post_norm_obj, current_s4)
  expect_identical(
    capture$save_state,
    list(
      state_name = "metab_normalized",
      s4_data_object = current_s4,
      config_object = workflow_data$config_list,
      description = "Between-sample normalization (method: none )"
    )
  )
  expect_identical(capture$log, "Between-sample normalization complete")
  expect_true(isTRUE(norm_data$normalization_complete))
  expect_identical(
    visible$value,
    list(
      currentS4 = current_s4,
      stateName = "metab_normalized",
      description = "Between-sample normalization (method: none )",
      logEntry = "Between-sample normalization complete"
    )
  )
})

test_that("metabolomics normalization module post-normalization QC helper preserves generation shell", {
  capture <- new.env(parent = emptyenv())
  capture$log <- character()

  current_s4 <- structure(list(stage = "post_norm_input"), class = "MetabNormPostQcStepMock")
  experiment_paths <- list(metabolite_qc_dir = "qc-dir")

  visible <- withVisible(
    runMetabNormPostNormalizationQcStep(
      currentS4 = current_s4,
      experimentPaths = experiment_paths,
      groupingVariable = "batch",
      shapeVariable = "group",
      addLogFn = function(message) {
        capture$log <- c(capture$log, message)
      },
      generateMetabQcPlotsFn = function(
        theObject,
        experiment_paths,
        stage,
        grouping_variable,
        shape_variable
      ) {
        capture$plot_call <- list(
          theObject = theObject,
          experiment_paths = experiment_paths,
          stage = stage,
          grouping_variable = grouping_variable,
          shape_variable = shape_variable
        )
        invisible(NULL)
      }
    )
  )

  expect_false(visible$visible)
  expect_identical(
    capture$plot_call,
    list(
      theObject = current_s4,
      experiment_paths = experiment_paths,
      stage = "post_norm",
      grouping_variable = "batch",
      shape_variable = "group"
    )
  )
  expect_identical(capture$log, "Post-normalization QC plots generated")
  expect_identical(
    visible$value,
    list(
      currentS4 = current_s4,
      stage = "post_norm",
      logEntry = "Post-normalization QC plots generated"
    )
  )
})

test_that("metabolomics normalization module RUV progress/apply shell preserves progress and branch handoff", {
  capture <- new.env(parent = emptyenv())
  capture$log <- character()

  current_s4 <- structure(list(stage = "post_norm"), class = "MetabNormRuvShellMock")
  updated_s4 <- structure(list(stage = "ruv_corrected"), class = "MetabNormRuvShellMock")
  experiment_paths <- list(metabolite_qc_dir = "qc-dir")
  workflow_data <- list(state_manager = "state-manager", config_list = list(ruv_mode = "RUVIII-C"))
  norm_data <- new.env(parent = emptyenv())
  optimization_state <- list(
    currentS4 = current_s4,
    bestKPerAssay = list(Plasma = 2),
    ctrlPerAssay = list(Plasma = c("M1", "M2"))
  )
  correction_state <- list(currentS4 = updated_s4, stateName = "metab_ruv_corrected")
  ruv_qc_state <- list(currentS4 = updated_s4, stage = "ruv_corrected")

  visible <- withVisible(
    runMetabNormRuvProgressApplyShell(
      currentS4 = current_s4,
      totalSteps = 6,
      ruvMode = "RUVIII-C",
      autoPercentageMin = 0.1,
      autoPercentageMax = 0.4,
      ruvGroupingVariable = "batch",
      separationMetric = "variance",
      kPenaltyWeight = 2,
      adaptiveKPenalty = TRUE,
      manualK = 3,
      manualPercentage = 0.25,
      experimentPaths = experiment_paths,
      groupingVariable = "batch",
      shapeVariable = "group",
      workflowData = workflow_data,
      normData = norm_data,
      addLogFn = function(message) {
        capture$log <- c(capture$log, message)
      },
      incProgressFn = function(amount, detail = NULL) {
        capture$progress <- list(amount = amount, detail = detail)
      },
      runRuvOptimizationStepFn = function(
        currentS4,
        ruvMode,
        autoPercentageMin,
        autoPercentageMax,
        ruvGroupingVariable,
        separationMetric,
        kPenaltyWeight,
        adaptiveKPenalty,
        manualK,
        manualPercentage,
        experimentPaths,
        normData,
        addLogFn
      ) {
        capture$optimization_call <- list(
          currentS4 = currentS4,
          ruvMode = ruvMode,
          autoPercentageMin = autoPercentageMin,
          autoPercentageMax = autoPercentageMax,
          ruvGroupingVariable = ruvGroupingVariable,
          separationMetric = separationMetric,
          kPenaltyWeight = kPenaltyWeight,
          adaptiveKPenalty = adaptiveKPenalty,
          manualK = manualK,
          manualPercentage = manualPercentage,
          experimentPaths = experimentPaths,
          normData = normData,
          addLogFn = addLogFn
        )
        optimization_state
      },
      runRuvCorrectionStepFn = function(
        currentS4,
        ruvGroupingVariable,
        bestKPerAssay,
        ctrlPerAssay,
        workflowData,
        normData,
        addLogFn
      ) {
        capture$correction_call <- list(
          currentS4 = currentS4,
          ruvGroupingVariable = ruvGroupingVariable,
          bestKPerAssay = bestKPerAssay,
          ctrlPerAssay = ctrlPerAssay,
          workflowData = workflowData,
          normData = normData,
          addLogFn = addLogFn
        )
        correction_state
      },
      runRuvQcStepFn = function(
        currentS4,
        totalSteps,
        experimentPaths,
        groupingVariable,
        shapeVariable,
        addLogFn
      ) {
        capture$qc_call <- list(
          currentS4 = currentS4,
          totalSteps = totalSteps,
          experimentPaths = experimentPaths,
          groupingVariable = groupingVariable,
          shapeVariable = shapeVariable,
          addLogFn = addLogFn
        )
        ruv_qc_state
      }
    )
  )

  expect_false(visible$visible)
  expect_identical(
    capture$progress,
    list(amount = 1 / 6, detail = "Running RUV-III batch correction...")
  )
  expect_identical(capture$log, "Running RUV-III (mode: RUVIII-C )")
  expect_identical(
    capture$optimization_call,
    list(
      currentS4 = current_s4,
      ruvMode = "RUVIII-C",
      autoPercentageMin = 0.1,
      autoPercentageMax = 0.4,
      ruvGroupingVariable = "batch",
      separationMetric = "variance",
      kPenaltyWeight = 2,
      adaptiveKPenalty = TRUE,
      manualK = 3,
      manualPercentage = 0.25,
      experimentPaths = experiment_paths,
      normData = norm_data,
      addLogFn = capture$optimization_call$addLogFn
    )
  )
  expect_true(is.function(capture$optimization_call$addLogFn))
  expect_identical(
    capture$correction_call,
    list(
      currentS4 = current_s4,
      ruvGroupingVariable = "batch",
      bestKPerAssay = list(Plasma = 2),
      ctrlPerAssay = list(Plasma = c("M1", "M2")),
      workflowData = workflow_data,
      normData = norm_data,
      addLogFn = capture$correction_call$addLogFn
    )
  )
  expect_true(is.function(capture$correction_call$addLogFn))
  expect_identical(
    capture$qc_call,
    list(
      currentS4 = updated_s4,
      totalSteps = 6,
      experimentPaths = experiment_paths,
      groupingVariable = "batch",
      shapeVariable = "group",
      addLogFn = capture$qc_call$addLogFn
    )
  )
  expect_true(is.function(capture$qc_call$addLogFn))
  expect_identical(
    visible$value,
    list(
      currentS4 = updated_s4,
      progressDetail = "Running RUV-III batch correction...",
      applyLogEntry = "Running RUV-III (mode: RUVIII-C )",
      skipLogEntry = NULL,
      optimizationState = optimization_state,
      correctionState = correction_state,
      ruvQcState = ruv_qc_state
    )
  )
})

test_that("metabolomics normalization module RUV progress/apply shell preserves skip branch state capture", {
  capture <- new.env(parent = emptyenv())
  capture$log <- character()
  capture$optimization_called <- FALSE
  capture$correction_called <- FALSE
  capture$qc_called <- FALSE

  current_s4 <- structure(list(stage = "post_norm"), class = "MetabNormRuvShellMock")
  experiment_paths <- list(metabolite_qc_dir = "qc-dir")
  workflow_data <- list(state_manager = "state-manager", config_list = list(ruv_mode = "skip"))
  norm_data <- new.env(parent = emptyenv())

  visible <- withVisible(
    runMetabNormRuvProgressApplyShell(
      currentS4 = current_s4,
      totalSteps = 6,
      ruvMode = "skip",
      autoPercentageMin = 0.1,
      autoPercentageMax = 0.4,
      ruvGroupingVariable = "batch",
      separationMetric = "variance",
      kPenaltyWeight = 2,
      adaptiveKPenalty = TRUE,
      manualK = 3,
      manualPercentage = 0.25,
      experimentPaths = experiment_paths,
      groupingVariable = "batch",
      shapeVariable = "group",
      workflowData = workflow_data,
      normData = norm_data,
      addLogFn = function(message) {
        capture$log <- c(capture$log, message)
      },
      incProgressFn = function(amount, detail = NULL) {
        capture$progress <- list(amount = amount, detail = detail)
      },
      runRuvOptimizationStepFn = function(...) {
        capture$optimization_called <- TRUE
        stop("runRuvOptimizationStepFn should not run when RUV is skipped")
      },
      runRuvCorrectionStepFn = function(...) {
        capture$correction_called <- TRUE
        stop("runRuvCorrectionStepFn should not run when RUV is skipped")
      },
      runRuvQcStepFn = function(...) {
        capture$qc_called <- TRUE
        stop("runRuvQcStepFn should not run when RUV is skipped")
      }
    )
  )

  expect_false(visible$visible)
  expect_identical(
    capture$progress,
    list(amount = 1 / 6, detail = "Running RUV-III batch correction...")
  )
  expect_identical(capture$log, "RUV-III skipped")
  expect_false(capture$optimization_called)
  expect_false(capture$correction_called)
  expect_false(capture$qc_called)
  expect_identical(norm_data$ruv_corrected_obj, current_s4)
  expect_true(isTRUE(norm_data$ruv_complete))
  expect_identical(
    visible$value,
    list(
      currentS4 = current_s4,
      progressDetail = "Running RUV-III batch correction...",
      applyLogEntry = NULL,
      skipLogEntry = "RUV-III skipped",
      optimizationState = NULL,
      correctionState = NULL,
      ruvQcState = NULL
    )
  )
})

test_that("metabolomics normalization module RUV optimization helper preserves opening shell", {
  capture <- new.env(parent = emptyenv())
  capture$log <- character()

  current_s4 <- structure(list(stage = "post_norm"), class = "MetabNormRuvOptimizationStepMock")
  experiment_paths <- list(metabolite_qc_dir = "qc-dir")
  norm_data <- new.env(parent = emptyenv())
  ruv_results <- list(Plasma = list(best_k = 2), Urine = list(best_k = 4))
  best_k_list <- list(Plasma = 2, Urine = 4)
  ctrl_list <- list(Plasma = c("M1", "M2"), Urine = "M3")

  visible <- withVisible(
    runMetabNormRuvOptimizationStep(
      currentS4 = current_s4,
      ruvMode = "RUVIII-C",
      autoPercentageMin = 0.1,
      autoPercentageMax = 0.4,
      ruvGroupingVariable = "batch",
      separationMetric = "variance",
      kPenaltyWeight = 2,
      adaptiveKPenalty = TRUE,
      manualK = 3,
      manualPercentage = 0.25,
      experimentPaths = experiment_paths,
      normData = norm_data,
      addLogFn = function(message) {
        capture$log <- c(capture$log, message)
      },
      runPerAssayRuvOptimizationFn = function(theObject, ruv_mode, params, experiment_paths) {
        capture$run_optimization <- list(
          theObject = theObject,
          ruv_mode = ruv_mode,
          params = params,
          experiment_paths = experiment_paths
        )
        ruv_results
      },
      extractBestKPerAssayFn = function(results) {
        capture$extract_best_k <- results
        best_k_list
      },
      extractCtrlPerAssayFn = function(results) {
        capture$extract_ctrl <- results
        ctrl_list
      }
    )
  )

  expect_false(visible$visible)
  expect_identical(
    capture$run_optimization,
    list(
      theObject = current_s4,
      ruv_mode = "RUVIII-C",
      params = list(
        percentage_min = 0.1,
        percentage_max = 0.4,
        ruv_grouping_variable = "batch",
        separation_metric = "variance",
        k_penalty_weight = 2,
        adaptive_k_penalty = TRUE,
        manual_k = 3,
        manual_percentage = 0.25
      ),
      experiment_paths = experiment_paths
    )
  )
  expect_identical(norm_data$ruv_optimization_results, ruv_results)
  expect_identical(capture$extract_best_k, ruv_results)
  expect_identical(capture$extract_ctrl, ruv_results)
  expect_identical(
    capture$log,
    "RUV optimization complete. Best k per assay: Plasma = 2, Urine = 4"
  )
  expect_identical(
    visible$value,
    list(
      currentS4 = current_s4,
      ruvParams = list(
        percentage_min = 0.1,
        percentage_max = 0.4,
        ruv_grouping_variable = "batch",
        separation_metric = "variance",
        k_penalty_weight = 2,
        adaptive_k_penalty = TRUE,
        manual_k = 3,
        manual_percentage = 0.25
      ),
      ruvResults = ruv_results,
      bestKPerAssay = best_k_list,
      ctrlPerAssay = ctrl_list,
      logEntry = "RUV optimization complete. Best k per assay: Plasma = 2, Urine = 4"
    )
  )
})

test_that("metabolomics normalization module RUV correction helper preserves apply/saveState shell", {
  capture <- new.env(parent = emptyenv())
  capture$log <- character()

  current_s4 <- structure(list(stage = "post_norm"), class = "MetabNormRuvCorrectionStepMock")
  updated_s4 <- structure(list(stage = "ruv_corrected"), class = "MetabNormRuvCorrectionStepMock")
  workflow_data <- new.env(parent = emptyenv())
  norm_data <- new.env(parent = emptyenv())
  best_k_list <- list(Plasma = 2, Urine = 4)
  ctrl_list <- list(Plasma = c("M1", "M2"), Urine = "M3")

  workflow_data$config_list <- list(ruv_mode = "RUVIII-C")
  workflow_data$state_manager <- list(
    saveState = function(state_name, s4_data_object, config_object, description) {
      capture$save_state <- list(
        state_name = state_name,
        s4_data_object = s4_data_object,
        config_object = config_object,
        description = description
      )
    }
  )

  visible <- withVisible(
    runMetabNormRuvCorrectionStep(
      currentS4 = current_s4,
      ruvGroupingVariable = "batch",
      bestKPerAssay = best_k_list,
      ctrlPerAssay = ctrl_list,
      workflowData = workflow_data,
      normData = norm_data,
      addLogFn = function(message) {
        capture$log <- c(capture$log, message)
      },
      ruvIII_C_VaryingFn = function(theObject, ruv_grouping_variable, ruv_number_k, ctrl) {
        capture$apply_ruv <- list(
          theObject = theObject,
          ruv_grouping_variable = ruv_grouping_variable,
          ruv_number_k = ruv_number_k,
          ctrl = ctrl
        )
        updated_s4
      }
    )
  )

  expect_false(visible$visible)
  expect_identical(
    capture$apply_ruv,
    list(
      theObject = current_s4,
      ruv_grouping_variable = "batch",
      ruv_number_k = best_k_list,
      ctrl = ctrl_list
    )
  )
  expect_identical(norm_data$ruv_corrected_obj, updated_s4)
  expect_true(isTRUE(norm_data$ruv_complete))
  expect_identical(
    capture$save_state,
    list(
      state_name = "metab_ruv_corrected",
      s4_data_object = updated_s4,
      config_object = workflow_data$config_list,
      description = "RUV-III batch correction complete"
    )
  )
  expect_identical(capture$log, "RUV-III correction applied")
  expect_identical(
    visible$value,
    list(
      currentS4 = updated_s4,
      stateName = "metab_ruv_corrected",
      description = "RUV-III batch correction complete",
      logEntry = "RUV-III correction applied"
    )
  )
})

test_that("metabolomics normalization module RUV QC helper preserves progress and generation shell", {
  capture <- new.env(parent = emptyenv())
  capture$log <- character()

  current_s4 <- structure(list(stage = "ruv_corrected"), class = "MetabNormRuvQcStepMock")
  experiment_paths <- list(metabolite_qc_dir = "qc-dir")

  visible <- withVisible(
    runMetabNormRuvQcStep(
      currentS4 = current_s4,
      totalSteps = 6,
      experimentPaths = experiment_paths,
      groupingVariable = "batch",
      shapeVariable = "group",
      addLogFn = function(message) {
        capture$log <- c(capture$log, message)
      },
      incProgressFn = function(amount, detail = NULL) {
        capture$progress <- list(amount = amount, detail = detail)
      },
      generateMetabQcPlotsFn = function(
        theObject,
        experiment_paths,
        stage,
        grouping_variable,
        shape_variable
      ) {
        capture$plot_call <- list(
          theObject = theObject,
          experiment_paths = experiment_paths,
          stage = stage,
          grouping_variable = grouping_variable,
          shape_variable = shape_variable
        )
        invisible(NULL)
      }
    )
  )

  expect_false(visible$visible)
  expect_identical(
    capture$progress,
    list(amount = 1 / 6, detail = "Generating RUV QC plots...")
  )
  expect_identical(
    capture$plot_call,
    list(
      theObject = current_s4,
      experiment_paths = experiment_paths,
      stage = "ruv_corrected",
      grouping_variable = "batch",
      shape_variable = "group"
    )
  )
  expect_identical(capture$log, "RUV QC plots generated")
  expect_identical(
    visible$value,
    list(
      currentS4 = current_s4,
      stage = "ruv_corrected",
      progressDetail = "Generating RUV QC plots...",
      logEntry = "RUV QC plots generated"
    )
  )
})

test_that("metabolomics normalization module composite image assembly seam preserves layout orchestration", {
  capture <- new.env(parent = emptyenv())
  capture$info <- character()
  capture$wraps <- list()

  plot_files <- c(
    "/tmp/plasma-pre.png",
    "/tmp/plasma-post.png",
    "/tmp/plasma-ruv.png",
    "/tmp/urine-pre.png",
    "/tmp/urine-post.png",
    "/tmp/urine-ruv.png"
  )

  visible <- withVisible(
    generateMetabNormCompositeFromFiles(
      plot_files = plot_files,
      ncol = 3,
      row_labels = list(
        plasma_pca = c("a)", "b)", "c)"),
        urine_pca = c("d)", "e)", "f)")
      ),
      column_labels = c("Pre", "Post", "RUV"),
      logInfoFn = function(message) {
        capture$info <- c(capture$info, message)
      },
      logErrorFn = function(message) {
        capture$error <- c(capture$error, message)
      },
      warningFn = function(message) {
        capture$warn <- c(capture$warn, message)
      },
      requireNamespaceFn = function(pkg, quietly = TRUE) {
        capture$packages <- c(capture$packages, pkg)
        TRUE
      },
      buildLabelPlotFn = function(title) paste0("label:", title),
      buildTitlePlotFn = function(title) paste0("title:", title),
      loadImageAsPlotFn = function(file_path) paste0("image:", basename(file_path)),
      wrapPlotsFn = function(plots, ncol) {
        wrapped <- list(plots = plots, ncol = ncol)
        capture$wraps <- c(capture$wraps, list(wrapped))
        wrapped
      },
      plotLayoutFn = function(heights) {
        capture$heights <- heights
        list(heights = heights)
      },
      combineLayoutFn = function(plot, layout) {
        capture$combined <- list(plot = plot, layout = layout)
        list(plot = plot, layout = layout)
      },
      fileExistsFn = function(path) {
        capture$exists <- c(capture$exists, path)
        TRUE
      },
      gcFn = function() {
        capture$gc_called <- TRUE
        invisible(NULL)
      }
    )
  )

  expect_true(visible$visible)
  expect_identical(capture$packages, c("patchwork", "ggplot2", "png"))
  expect_null(capture$error)
  expect_null(capture$warn)
  expect_identical(
    capture$info,
    c(
      "[generateCompositeFromFiles] Generating composite from 6 files...",
      "[generateCompositeFromFiles] Added column titles",
      "[generateCompositeFromFiles] Added row: plasma_pca",
      "[generateCompositeFromFiles] Added row: urine_pca",
      "[generateCompositeFromFiles] Combining plot sections..."
    )
  )
  expect_length(capture$wraps, 6)
  expect_identical(
    capture$wraps[[1]],
    list(plots = list("title:Pre", "title:Post", "title:RUV"), ncol = 3)
  )
  expect_identical(
    capture$wraps[[2]],
    list(plots = list("label:a)", "label:b)", "label:c)"), ncol = 3)
  )
  expect_identical(
    capture$wraps[[3]],
    list(
      plots = list("image:plasma-pre.png", "image:plasma-post.png", "image:plasma-ruv.png"),
      ncol = 3
    )
  )
  expect_identical(capture$heights, c(0.2, 0.1, 1, 0.1, 1))
  expect_true(capture$gc_called)
  expect_identical(
    visible$value,
    list(
      plot = list(
        plot = capture$wraps[[6]],
        layout = list(heights = c(0.2, 0.1, 1, 0.1, 1))
      ),
      width = 13,
      height = 14
    )
  )
})

test_that("metabolomics normalization module composite QC helper preserves generation and save shell", {
  capture <- new.env(parent = emptyenv())
  capture$log <- character()
  capture$warn <- NULL

  experiment_paths <- list(metabolite_qc_dir = "/tmp/metab-qc")
  assay_names <- c("Plasma Panel", "Urine Panel")

  visible <- withVisible(
    runMetabNormCompositeQcFigureStep(
      experimentPaths = experiment_paths,
      assayNames = assay_names,
      ruvMode = "per_assay",
      omicType = "metabolite",
      addLogFn = function(message) {
        capture$log <- c(capture$log, message)
      },
      generateCompositeFromFilesFn = function(plot_files, ncol, row_labels, column_labels) {
        capture$composite <- list(
          plot_files = plot_files,
          ncol = ncol,
          row_labels = row_labels,
          column_labels = column_labels
        )
        list(plot = "combined-plot", width = 13, height = 19)
      },
      savePlotFn = function(plot, base_path, plot_name, width, height, dpi, limitsize) {
        capture$save <- list(
          plot = plot,
          base_path = base_path,
          plot_name = plot_name,
          width = width,
          height = height,
          dpi = dpi,
          limitsize = limitsize
        )
        invisible(NULL)
      },
      dirExistsFn = function(path) {
        capture$dir_path <- path
        TRUE
      },
      logWarnFn = function(message) {
        capture$warn <- c(capture$warn, message)
      }
    )
  )

  expect_false(visible$visible)
  expect_identical(capture$dir_path, "/tmp/metab-qc")
  expect_identical(
    capture$composite$column_labels,
    c("Pre-Normalisation", "Post-Normalisation", "RUV-Corrected")
  )
  expect_identical(capture$composite$ncol, 3)
  expect_identical(length(capture$composite$plot_files), 24L)
  expect_identical(
    capture$composite$plot_files[1:3],
    c(
      "/tmp/metab-qc/plasma_panel_pre_norm_pca.png",
      "/tmp/metab-qc/plasma_panel_post_norm_pca.png",
      "/tmp/metab-qc/plasma_panel_ruv_corrected_pca.png"
    )
  )
  expect_identical(
    capture$composite$row_labels[1:2],
    list(
      plasma_panel_pca = c("a)", "b)", "c)"),
      plasma_panel_density = c("d)", "e)", "f)")
    )
  )
  expect_identical(
    capture$save,
    list(
      plot = "combined-plot",
      base_path = "/tmp/metab-qc",
      plot_name = "metabolite_composite_QC_figure",
      width = 13,
      height = 19,
      dpi = 150,
      limitsize = FALSE
    )
  )
  expect_null(capture$warn)
  expect_identical(
    capture$log,
    c(
      "Generating composite QC figure...",
      "Composite QC figure saved to: /tmp/metab-qc/composite_QC_figure"
    )
  )
  expect_identical(
    visible$value,
    list(
      qcDir = "/tmp/metab-qc",
      ncolComposite = 3,
      columnLabels = c("Pre-Normalisation", "Post-Normalisation", "RUV-Corrected"),
      allPlotFiles = capture$composite$plot_files,
      allRowLabels = capture$composite$row_labels,
      logEntry = "Generating composite QC figure...",
      compositeSaved = TRUE
    )
  )
})

test_that("metabolomics normalization module composite QC refresh shell preserves final tail handoff", {
  capture <- new.env(parent = emptyenv())
  capture$log <- character()

  current_s4 <- structure(list(stage = "ruv_corrected"), class = "MetabNormCompositeRefreshMock")
  experiment_paths <- list(metabolite_qc_dir = "/tmp/metab-qc")
  assay_names <- c("Plasma", "Urine")
  norm_data <- new.env(parent = emptyenv())
  norm_data$plot_refresh_trigger <- 4

  add_log <- function(message) {
    capture$log <- c(capture$log, message)
  }
  composite_generator <- function(...) {
    "unused-generator"
  }
  save_plot <- function(...) {
    "unused-save"
  }
  warn_log <- function(...) {
    "unused-warn"
  }
  composite_state <- list(qcDir = "/tmp/metab-qc", compositeSaved = TRUE)

  visible <- withVisible(
    runMetabNormCompositeQcRefreshShell(
      currentS4 = current_s4,
      experimentPaths = experiment_paths,
      assayNames = assay_names,
      ruvMode = "per_assay",
      omicType = "metabolite",
      normData = norm_data,
      addLogFn = add_log,
      generateCompositeFromFilesFn = composite_generator,
      savePlotFn = save_plot,
      logWarnFn = warn_log,
      runCompositeQcFigureStepFn = function(
        experimentPaths,
        assayNames,
        ruvMode,
        omicType,
        addLogFn,
        generateCompositeFromFilesFn,
        savePlotFn,
        logWarnFn
      ) {
        capture$call <- list(
          experimentPaths = experimentPaths,
          assayNames = assayNames,
          ruvMode = ruvMode,
          omicType = omicType,
          addLogFn = addLogFn,
          generateCompositeFromFilesFn = generateCompositeFromFilesFn,
          savePlotFn = savePlotFn,
          logWarnFn = logWarnFn
        )
        addLogFn("Composite refresh shell delegated")
        composite_state
      }
    )
  )

  expect_false(visible$visible)
  expect_identical(
    capture$call,
    list(
      experimentPaths = experiment_paths,
      assayNames = assay_names,
      ruvMode = "per_assay",
      omicType = "metabolite",
      addLogFn = add_log,
      generateCompositeFromFilesFn = composite_generator,
      savePlotFn = save_plot,
      logWarnFn = warn_log
    )
  )
  expect_identical(capture$log, "Composite refresh shell delegated")
  expect_identical(norm_data$plot_refresh_trigger, 5)
  expect_identical(
    visible$value,
    list(
      currentS4 = current_s4,
      compositeQcState = composite_state,
      plotRefreshTrigger = 5
    )
  )
})

test_that("metabolomics normalization module summary helper preserves correlation summary formatting", {
  corr_results <- list(
    Plasma = data.frame(
      pearson_correlation = c(0.5, 0.7, 0.9)
    ),
    Urine = data.frame(
      pearson_correlation = c(0.8, 0.85)
    )
  )

  original_object <- methods::new(
    "MetabNormSummaryMock",
    design_matrix = data.frame(sample_id = paste0("S", 1:5))
  )
  filtered_object <- methods::new(
    "MetabNormSummaryMock",
    design_matrix = data.frame(sample_id = paste0("S", 1:3))
  )

  expect_identical(
    buildMetabNormCorrelationFilterSummary(
      corrResults = corr_results,
      filteredObject = filtered_object,
      originalObject = original_object
    ),
    paste0(
      "=== Correlation Filtering Summary ===\n",
      "\n[Plasma]\n  Sample pairs: 3\n  Correlation: mean=0.700, min=0.500, max=0.900",
      "\n[Urine]\n  Sample pairs: 2\n  Correlation: mean=0.825, min=0.800, max=0.850",
      "\n\n[Sample Filtering]\n  Original: 5 samples\n  After filtering: 3 samples\n  Removed: 2 samples"
    )
  )
})

test_that("metabolomics normalization module summary helper preserves empty-result fallback", {
  expect_identical(
    buildMetabNormCorrelationFilterSummary(corrResults = list()),
    "No correlation results available."
  )

  expect_identical(
    buildMetabNormCorrelationFilterSummary(corrResults = NULL),
    "No correlation results available."
  )
})

test_that("metabolomics normalization module final QC render-state helper preserves source priority", {
  post_norm_object <- structure(list(stage = "post_norm"), class = "MetabNormRenderMock")
  ruv_corrected_object <- structure(list(stage = "ruv_corrected"), class = "MetabNormRenderMock")
  correlation_filtered_object <- structure(list(stage = "correlation_filter"), class = "MetabNormRenderMock")

  expect_identical(
    resolveMetabNormFinalQcRenderState(
      correlationFilteredObject = correlation_filtered_object,
      ruvCorrectedObject = ruv_corrected_object,
      postNormObject = post_norm_object
    ),
    list(
      sourceObject = correlation_filtered_object,
      sourceStage = "correlation_filter",
      isFallback = FALSE,
      plot = NULL
    )
  )

  expect_identical(
    resolveMetabNormFinalQcRenderState(
      correlationFilteredObject = NULL,
      ruvCorrectedObject = ruv_corrected_object,
      postNormObject = post_norm_object
    ),
    list(
      sourceObject = ruv_corrected_object,
      sourceStage = "ruv_corrected",
      isFallback = FALSE,
      plot = NULL
    )
  )

  expect_identical(
    resolveMetabNormFinalQcRenderState(
      correlationFilteredObject = NULL,
      ruvCorrectedObject = NULL,
      postNormObject = post_norm_object
    ),
    list(
      sourceObject = post_norm_object,
      sourceStage = "post_norm",
      isFallback = FALSE,
      plot = NULL
    )
  )
})

test_that("metabolomics normalization module final QC render-state helper preserves empty fallback", {
  render_state <- resolveMetabNormFinalQcRenderState(
    correlationFilteredObject = NULL,
    ruvCorrectedObject = NULL,
    postNormObject = NULL
  )

  expect_identical(render_state$sourceObject, NULL)
  expect_identical(render_state$sourceStage, "empty")
  expect_true(isTRUE(render_state$isFallback))
  expect_s3_class(render_state$plot, "ggplot")
})

test_that("metabolomics normalization module export-source helper preserves source-dir priority", {
  capture <- new.env(parent = emptyenv())
  capture$dir_exists <- character()
  capture$dir_create <- list()

  resolved_dir <- resolveMetabNormExportSourceDir(
    experimentPaths = list(
      source_dir = "existing-source",
      export_dir = "existing-export"
    ),
    dirExistsFn = function(path) {
      capture$dir_exists <- c(capture$dir_exists, path)
      identical(path, "existing-source")
    },
    dirCreateFn = function(...) {
      capture$dir_create <- c(capture$dir_create, list(list(...)))
      TRUE
    }
  )

  expect_identical(resolved_dir, "existing-source")
  expect_identical(capture$dir_exists, "existing-source")
  expect_length(capture$dir_create, 0)
})

test_that("metabolomics normalization module export-source helper preserves export-dir fallback creation", {
  capture <- new.env(parent = emptyenv())
  capture$dir_exists <- character()
  capture$dir_create <- list()

  resolved_dir <- resolveMetabNormExportSourceDir(
    experimentPaths = list(
      source_dir = "missing-source",
      export_dir = "new-export-dir"
    ),
    dirExistsFn = function(path) {
      capture$dir_exists <- c(capture$dir_exists, path)
      FALSE
    },
    dirCreateFn = function(path, recursive = FALSE) {
      capture$dir_create <- list(path = path, recursive = recursive)
      TRUE
    }
  )

  expect_identical(resolved_dir, "new-export-dir")
  expect_identical(capture$dir_exists, c("missing-source", "new-export-dir"))
  expect_identical(
    capture$dir_create,
    list(path = "new-export-dir", recursive = TRUE)
  )
})

test_that("metabolomics normalization module export-source helper preserves missing-directory error", {
  expect_error(
    resolveMetabNormExportSourceDir(
      experimentPaths = list(
        source_dir = "missing-source",
        export_dir = NULL
      ),
      dirExistsFn = function(path) {
        FALSE
      },
      dirCreateFn = function(...) {
        stop("dirCreateFn should not be called when export_dir is NULL")
      }
    ),
    "Could not find a valid directory to save session data."
  )
})

test_that("metabolomics normalization module feature-count helper preserves invalid-object fallback", {
  expect_null(collectMetabNormFeatureCountsPerAssay(NULL))
  expect_null(
    collectMetabNormFeatureCountsPerAssay(
      structure(list(), class = "NotMetaboliteAssayData")
    )
  )
})

test_that("metabolomics normalization module feature-count helper preserves per-assay counts", {
  current_s4 <- methods::new(
    "MetaboliteAssayData",
    metabolite_data = list(
      Plasma = data.frame(
        metabolite_id = c("M1", "M1", "M2"),
        annotation_id = c("A1", "A1", "A2"),
        Sample_1 = c(10, 11, 12),
        Sample_2 = c(20, 21, 22),
        stringsAsFactors = FALSE
      ),
      Urine = NULL
    ),
    metabolite_id_column = "metabolite_id",
    annotation_id_column = "annotation_id"
  )

  feature_counts <- collectMetabNormFeatureCountsPerAssay(current_s4)

  expect_named(feature_counts, c("Plasma", "Urine"))
  expect_equal(feature_counts$Plasma, list(features = 2, samples = 2))
  expect_equal(feature_counts$Urine, list(features = 0, samples = 0))
})

test_that("metabolomics normalization module export-session payload helper preserves state payload assembly", {
  capture <- new.env(parent = emptyenv())
  current_s4 <- methods::new(
    "MetaboliteAssayData",
    metabolite_data = list(
      Plasma = data.frame(
        metabolite_id = c("M1", "M2"),
        annotation_id = c("A1", "A2"),
        Sample_1 = c(10, 12),
        Sample_2 = c(20, 22),
        stringsAsFactors = FALSE
      )
    ),
    metabolite_id_column = "metabolite_id",
    annotation_id_column = "annotation_id"
  )
  workflow_data <- list(
    state_manager = list(
      current_state = "correlation_filtered",
      getState = function(state_name) {
        capture$state_name <- state_name
        current_s4
      }
    ),
    contrasts_tbl = data.frame(friendly_names = "Treatment-Control", stringsAsFactors = FALSE),
    design_matrix = data.frame(sample = c("Sample_1", "Sample_2"), stringsAsFactors = FALSE),
    config_list = list(method = "TIC"),
    metabolite_counts = c(Plasma = 2L),
    qc_params = list(min_present = 0.8)
  )
  norm_data <- list(
    itsd_selections = list(Plasma = "M1"),
    ruv_optimization_results = list(Plasma = list(success = TRUE, best_k = 2)),
    correlation_results = list(Plasma = data.frame(pearson_correlation = c(0.9, 0.95))),
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
    min_pearson_correlation_threshold = 0.85,
    ruv_grouping_variable = "Batch"
  )
  export_timestamp <- as.POSIXct("2026-04-17 12:34:56", tz = "UTC")

  session_data <- buildMetabNormExportSessionData(
    workflowData = workflow_data,
    normData = norm_data,
    inputValues = input_values,
    experimentLabel = "Demo Study",
    exportTimestamp = export_timestamp
  )

  expect_identical(capture$state_name, "correlation_filtered")
  expect_identical(session_data$r6_current_state_name, "correlation_filtered")
  expect_identical(session_data$current_s4_object, current_s4)
  expect_identical(session_data$feature_counts$Plasma, list(features = 2L, samples = 2L))
  expect_identical(session_data$itsd_applied, TRUE)
  expect_identical(session_data$itsd_aggregation, "mean")
  expect_identical(session_data$export_timestamp, export_timestamp)
  expect_identical(session_data$omic_type, "metabolomics")
  expect_identical(session_data$experiment_label, "Demo Study")
  expect_identical(session_data$qc_params, workflow_data$qc_params)
  expect_identical(session_data$normalization_complete, TRUE)
  expect_identical(session_data$ruv_complete, TRUE)
  expect_identical(session_data$correlation_filtering_complete, TRUE)
})

test_that("metabolomics normalization module export-session payload helper preserves ITSD fallback", {
  workflow_data <- list(
    state_manager = list(
      current_state = "post_norm",
      getState = function(state_name) {
        structure(list(state = state_name), class = "NotMetaboliteAssayData")
      }
    ),
    contrasts_tbl = NULL,
    design_matrix = NULL,
    config_list = NULL,
    metabolite_counts = NULL,
    qc_params = NULL
  )
  norm_data <- list(
    itsd_selections = list(),
    ruv_optimization_results = list(),
    correlation_results = list(),
    assay_names = c("Plasma"),
    normalization_complete = TRUE,
    ruv_complete = FALSE,
    correlation_filtering_complete = FALSE
  )
  input_values <- list(
    norm_method = "Median",
    ruv_mode = "skip",
    apply_itsd = FALSE,
    itsd_aggregation = "median",
    log_offset = 0.5,
    min_pearson_correlation_threshold = NULL,
    ruv_grouping_variable = NULL
  )

  session_data <- buildMetabNormExportSessionData(
    workflowData = workflow_data,
    normData = norm_data,
    inputValues = input_values,
    experimentLabel = "Fallback Study",
    exportTimestamp = as.POSIXct("2026-04-17 13:00:00", tz = "UTC")
  )

  expect_null(session_data$feature_counts)
  expect_true(is.na(session_data$itsd_aggregation))
  expect_identical(session_data$itsd_applied, FALSE)
  expect_identical(session_data$ruv_complete, FALSE)
  expect_identical(session_data$correlation_filtering_complete, FALSE)
})

test_that("metabolomics normalization module export-session RDS helper preserves main/latest writes", {
  capture <- new.env(parent = emptyenv())
  capture$saves <- list()
  capture$info <- character()
  capture$progress <- list()
  session_data <- list(export = "payload")

  visible <- withVisible(
    saveMetabNormExportSessionRdsFiles(
      sessionData = session_data,
      sourceDir = "/tmp/metab-export",
      timeFn = function() {
        as.POSIXct("2026-04-17 14:15:16", tz = "UTC")
      },
      formatTimeFn = function(x, fmt) {
        format(x, fmt, tz = "UTC")
      },
      saveRdsFn = function(object, path) {
        capture$saves <- c(capture$saves, list(list(object = object, path = path)))
      },
      logInfoFn = function(message) {
        capture$info <- c(capture$info, message)
      },
      incProgressFn = function(value, detail = NULL) {
        capture$progress <- c(
          capture$progress,
          list(list(value = value, detail = detail))
        )
      }
    )
  )

  expect_false(visible$visible)
  expect_identical(
    visible$value,
    list(
      sessionFilename = "metab_filtered_session_data_20260417_141516.rds",
      sessionFilepath = "/tmp/metab-export/metab_filtered_session_data_20260417_141516.rds",
      latestFilename = "metab_filtered_session_data_latest.rds",
      latestFilepath = "/tmp/metab-export/metab_filtered_session_data_latest.rds"
    )
  )
  expect_length(capture$saves, 2)
  expect_true(all(vapply(capture$saves, function(entry) identical(entry$object, session_data), logical(1))))
  expect_identical(
    vapply(capture$saves, `[[`, character(1), "path"),
    c(
      "/tmp/metab-export/metab_filtered_session_data_20260417_141516.rds",
      "/tmp/metab-export/metab_filtered_session_data_latest.rds"
    )
  )
  expect_identical(
    capture$info,
    c(
      "*** EXPORT: Session data saved to: /tmp/metab-export/metab_filtered_session_data_20260417_141516.rds ***",
      "*** EXPORT: Latest version saved to: /tmp/metab-export/metab_filtered_session_data_latest.rds ***"
    )
  )
  expect_identical(
    capture$progress,
    list(
      list(value = 0.3, detail = "Saving to file..."),
      list(value = 0.1, detail = "Creating latest version...")
    )
  )
})

test_that("metabolomics normalization module export metadata helper preserves redundancy writes", {
  capture <- new.env(parent = emptyenv())
  capture$saves <- list()
  capture$info <- character()
  capture$warn <- character()

  visible <- withVisible(
    saveMetabNormExportMetadataFiles(
      sessionData = list(
        ruv_optimization_results = list(Plasma = list(success = TRUE, best_k = 2)),
        itsd_selections = list(Plasma = "M1"),
        qc_params = list(min_present = 0.8)
      ),
      sourceDir = "/tmp/metab-export",
      saveRdsFn = function(object, path) {
        capture$saves <- c(capture$saves, list(list(object = object, path = path)))
      },
      logInfoFn = function(message) {
        capture$info <- c(capture$info, message)
      },
      logWarnFn = function(message) {
        capture$warn <- c(capture$warn, message)
      }
    )
  )

  expect_null(visible$value)
  expect_false(visible$visible)
  expect_length(capture$saves, 3)
  expect_identical(
    vapply(capture$saves, `[[`, character(1), "path"),
    c(
      "/tmp/metab-export/metab_ruv_optimization_results.RDS",
      "/tmp/metab-export/metab_itsd_selections.RDS",
      "/tmp/metab-export/metab_qc_params.RDS"
    )
  )
  expect_identical(
    capture$info,
    c(
      "*** EXPORT: Saved metab_ruv_optimization_results.RDS ***",
      "*** EXPORT: Saved metab_itsd_selections.RDS ***",
      "*** EXPORT: Saved metab_qc_params.RDS ***"
    )
  )
  expect_identical(capture$warn, character())
})

test_that("metabolomics normalization module export metadata helper preserves warning fallback", {
  capture <- new.env(parent = emptyenv())
  capture$paths <- character()
  capture$info <- character()
  capture$warn <- character()

  visible <- withVisible(
    saveMetabNormExportMetadataFiles(
      sessionData = list(
        ruv_optimization_results = list(Plasma = list(success = TRUE, best_k = 2)),
        itsd_selections = list(Plasma = "M1"),
        qc_params = list(min_present = 0.8)
      ),
      sourceDir = "/tmp/metab-export",
      saveRdsFn = function(object, path) {
        capture$paths <- c(capture$paths, path)
        stop("disk full")
      },
      logInfoFn = function(message) {
        capture$info <- c(capture$info, message)
      },
      logWarnFn = function(message) {
        capture$warn <- c(capture$warn, message)
      }
    )
  )

  expect_null(visible$value)
  expect_false(visible$visible)
  expect_identical(
    capture$paths,
    "/tmp/metab-export/metab_ruv_optimization_results.RDS"
  )
  expect_identical(capture$info, character())
  expect_identical(
    capture$warn,
    "*** WARNING: Some metadata files could not be saved: disk full ***"
  )
})

test_that("metabolomics normalization module export summary helper preserves summary writes", {
  capture <- new.env(parent = emptyenv())
  capture$content <- NULL
  capture$path <- NULL
  capture$info <- character()

  visible <- withVisible(
    saveMetabNormExportSummaryFile(
      sessionData = list(
        ruv_optimization_results = list(
          Plasma = list(
            success = TRUE,
            best_k = 2,
            best_percentage = 15.5,
            control_genes_index = c(TRUE, FALSE, TRUE)
          ),
          Urine = list(
            success = FALSE,
            best_k = 1,
            best_percentage = 10,
            control_genes_index = c(TRUE)
          )
        ),
        feature_counts = list(
          Plasma = list(features = 2, samples = 3),
          Urine = list(features = 1, samples = 2)
        ),
        normalization_method = "TIC",
        itsd_applied = TRUE,
        itsd_aggregation = "mean",
        log_offset = 1,
        ruv_mode = "RUVIII-C",
        ruv_grouping_variable = "Batch",
        correlation_threshold = 0.8,
        contrasts_tbl = data.frame(
          friendly_names = c("B vs A", "C vs A"),
          stringsAsFactors = FALSE
        )
      ),
      sourceDir = "/tmp/metab-export",
      sessionFilename = "metab_filtered_session_data_20260417_123456.rds",
      writeLinesFn = function(text, path) {
        capture$content <- text
        capture$path <- path
      },
      timeFn = function() {
        as.POSIXct("2026-04-17 12:34:56", tz = "UTC")
      },
      logInfoFn = function(message) {
        capture$info <- c(capture$info, message)
      }
    )
  )

  expect_false(visible$visible)
  expect_identical(
    visible$value,
    list(
      summaryContent = capture$content,
      summaryFilepath = "/tmp/metab-export/metab_filtered_session_summary.txt"
    )
  )
  expect_identical(capture$path, "/tmp/metab-export/metab_filtered_session_summary.txt")
  expect_identical(
    capture$content,
    paste0(
      "Metabolomics Normalized Session Data Export Summary\n",
      "===================================================\n\n",
      "Export Timestamp: 2026-04-17 12:34:56\n",
      "Session File: metab_filtered_session_data_20260417_123456.rds\n\n",
      "Data Summary:\n",
      "  Plasma: 2 features, 3 samples\n",
      "  Urine: 1 features, 2 samples\n\n",
      "Normalization Parameters:\n",
      "- Method: TIC\n",
      "- ITSD applied: Yes\n",
      "- ITSD aggregation: mean\n",
      "- Log2 offset: 1\n",
      "- RUV mode: RUVIII-C\n",
      "- RUV grouping variable: Batch\n",
      "- Correlation threshold: 0.8\n\n",
      "RUV Optimization Results (per-assay):\n",
      "  Plasma: k=2, %=15.5, controls=2\n\n",
      "Contrasts:\n",
      "B vs A\n",
      "C vs A\n\n",
      "This data is ready for differential expression analysis.\n",
      "Use 'Load Filtered Session' in the DE tab to import.\n"
    )
  )
  expect_identical(
    capture$info,
    "*** EXPORT: Summary saved to: /tmp/metab-export/metab_filtered_session_summary.txt ***"
  )
})

test_that("metabolomics normalization module export summary helper preserves empty fallbacks", {
  capture <- new.env(parent = emptyenv())
  capture$content <- NULL
  capture$path <- NULL
  capture$info <- character()

  visible <- withVisible(
    saveMetabNormExportSummaryFile(
      sessionData = list(
        ruv_optimization_results = list(),
        feature_counts = list(),
        normalization_method = "Median",
        itsd_applied = FALSE,
        itsd_aggregation = NA_character_,
        log_offset = 0.5,
        ruv_mode = "skip",
        ruv_grouping_variable = "None",
        correlation_threshold = NULL,
        contrasts_tbl = NULL
      ),
      sourceDir = "/tmp/metab-export",
      sessionFilename = "metab_filtered_session_data_latest.rds",
      writeLinesFn = function(text, path) {
        capture$content <- text
        capture$path <- path
      },
      timeFn = function() {
        as.POSIXct("2026-04-17 13:00:00", tz = "UTC")
      },
      logInfoFn = function(message) {
        capture$info <- c(capture$info, message)
      }
    )
  )

  expect_false(visible$visible)
  expect_identical(
    visible$value,
    list(
      summaryContent = capture$content,
      summaryFilepath = "/tmp/metab-export/metab_filtered_session_summary.txt"
    )
  )
  expect_identical(capture$path, "/tmp/metab-export/metab_filtered_session_summary.txt")
  expect_identical(
    capture$content,
    paste0(
      "Metabolomics Normalized Session Data Export Summary\n",
      "===================================================\n\n",
      "Export Timestamp: 2026-04-17 13:00:00\n",
      "Session File: metab_filtered_session_data_latest.rds\n\n",
      "Data Summary:\n\n",
      "Normalization Parameters:\n",
      "- Method: Median\n",
      "- ITSD applied: No\n",
      "- ITSD aggregation: N/A\n",
      "- Log2 offset: 0.5\n",
      "- RUV mode: skip\n",
      "- RUV grouping variable: None\n",
      "- Correlation threshold: N/A\n\n",
      "RUV Optimization Results (per-assay):\n",
      "  (RUV skipped or not applied)\n\n",
      "Contrasts:\n",
      "None defined\n\n",
      "This data is ready for differential expression analysis.\n",
      "Use 'Load Filtered Session' in the DE tab to import.\n"
    )
  )
  expect_identical(
    capture$info,
    "*** EXPORT: Summary saved to: /tmp/metab-export/metab_filtered_session_summary.txt ***"
  )
})

test_that("metabolomics normalization module export workflow helper preserves withProgress orchestration", {
  capture <- new.env(parent = emptyenv())
  capture$info <- character()
  capture$progress <- list()

  workflow_data <- list(
    state_manager = list(current_state = "correlation_filtered")
  )
  norm_data <- list(
    assay_names = c("Plasma", "Urine")
  )
  input_values <- list(
    norm_method = "TIC"
  )
  session_data <- list(
    contrasts_tbl = data.frame(friendly_names = "B vs A", stringsAsFactors = FALSE),
    export = "payload"
  )

  visible <- withVisible(
    runMetabNormExportSessionWorkflow(
      workflowData = workflow_data,
      normData = norm_data,
      inputValues = input_values,
      experimentLabel = "Demo Study",
      sourceDir = "/tmp/metab-export",
      withProgressFn = function(message = NULL, detail = NULL, value = NULL, expr) {
        capture$with_progress <- list(message = message, detail = detail, value = value)
        result <- withVisible(force(expr))
        if (isTRUE(result$visible)) {
          result$value
        } else {
          invisible(result$value)
        }
      },
      incProgressFn = function(value, detail = NULL) {
        capture$progress <- c(
          capture$progress,
          list(list(value = value, detail = detail))
        )
      },
      buildSessionDataFn = function(workflowData, normData, inputValues, experimentLabel) {
        capture$build <- list(
          workflowData = workflowData,
          normData = normData,
          inputValues = inputValues,
          experimentLabel = experimentLabel
        )
        session_data
      },
      saveSessionRdsFilesFn = function(sessionData, sourceDir, incProgressFn) {
        capture$save_rds <- list(sessionData = sessionData, sourceDir = sourceDir)
        incProgressFn(0.3, detail = "Saving to file...")
        incProgressFn(0.1, detail = "Creating latest version...")
        list(
          sessionFilename = "metab_filtered_session_data_demo.rds",
          sessionFilepath = "/tmp/metab-export/metab_filtered_session_data_demo.rds",
          latestFilename = "metab_filtered_session_data_latest.rds",
          latestFilepath = "/tmp/metab-export/metab_filtered_session_data_latest.rds"
        )
      },
      saveMetadataFilesFn = function(sessionData, sourceDir) {
        capture$save_metadata <- list(sessionData = sessionData, sourceDir = sourceDir)
        NULL
      },
      saveSummaryFileFn = function(sessionData, sourceDir, sessionFilename) {
        capture$save_summary <- list(
          sessionData = sessionData,
          sourceDir = sourceDir,
          sessionFilename = sessionFilename
        )
        NULL
      },
      logInfoFn = function(message) {
        capture$info <- c(capture$info, message)
      }
    )
  )

  expect_false(visible$visible)
  expect_identical(
    visible$value,
    list(
      sessionFilename = "metab_filtered_session_data_demo.rds",
      sessionFilepath = "/tmp/metab-export/metab_filtered_session_data_demo.rds",
      latestFilename = "metab_filtered_session_data_latest.rds",
      latestFilepath = "/tmp/metab-export/metab_filtered_session_data_latest.rds"
    )
  )
  expect_identical(
    capture$with_progress,
    list(
      message = "Exporting normalized session data...",
      detail = NULL,
      value = 0
    )
  )
  expect_identical(capture$build$workflowData, workflow_data)
  expect_identical(capture$build$normData, norm_data)
  expect_identical(capture$build$inputValues, input_values)
  expect_identical(capture$build$experimentLabel, "Demo Study")
  expect_identical(capture$save_rds, list(sessionData = session_data, sourceDir = "/tmp/metab-export"))
  expect_identical(capture$save_metadata, list(sessionData = session_data, sourceDir = "/tmp/metab-export"))
  expect_identical(
    capture$save_summary,
    list(
      sessionData = session_data,
      sourceDir = "/tmp/metab-export",
      sessionFilename = "metab_filtered_session_data_demo.rds"
    )
  )
  expect_identical(
    capture$info,
    c(
      "*** EXPORT: Gathered session data successfully ***",
      "*** EXPORT: Assays: Plasma, Urine ***",
      "*** EXPORT: Contrasts available: 1 ***"
    )
  )
  expect_identical(
    capture$progress,
    list(
      list(value = 0.2, detail = "Gathering data..."),
      list(value = 0.3, detail = "Saving to file..."),
      list(value = 0.1, detail = "Creating latest version..."),
      list(value = 0.1, detail = "Saving metadata files..."),
      list(value = 0.2, detail = "Creating summary...")
    )
  )
})

test_that("metabolomics normalization module export ready helper preserves warning gate", {
  capture <- new.env(parent = emptyenv())
  capture$notifications <- list()

  expect_false(
    checkMetabNormExportSessionReady(
      normalizationComplete = FALSE,
      showNotificationFn = function(message, type = NULL, duration = NULL) {
        capture$notifications <- c(
          capture$notifications,
          list(list(message = message, type = type, duration = duration))
        )
        invisible(NULL)
      }
    )
  )

  expect_identical(
    capture$notifications,
    list(list(
      message = "Please complete normalization before exporting session data.",
      type = "warning",
      duration = 5
    ))
  )
})

test_that("metabolomics normalization module export ready helper preserves completed state passthrough", {
  expect_true(
    checkMetabNormExportSessionReady(
      normalizationComplete = TRUE,
      showNotificationFn = function(...) {
        stop("showNotificationFn should not be called")
      }
    )
  )
})

test_that("metabolomics normalization module export dispatch helper preserves success shell", {
  capture <- new.env(parent = emptyenv())
  capture$log <- character()

  visible <- withVisible(
    dispatchMetabNormExportSession(
      workflowData = list(state_manager = "workflow-state"),
      normData = list(assay_names = "Plasma"),
      inputValues = list(norm_method = "TIC"),
      experimentPaths = list(source_dir = "/tmp/source"),
      experimentLabel = "Demo Study",
      addLogFn = function(message) {
        capture$log <- c(capture$log, message)
      },
      resolveSourceDirFn = function(experimentPaths) {
        capture$resolve_paths <- experimentPaths
        "/tmp/metab-export"
      },
      runWorkflowFn = function(workflowData, normData, inputValues, experimentLabel, sourceDir) {
        capture$workflow <- list(
          workflowData = workflowData,
          normData = normData,
          inputValues = inputValues,
          experimentLabel = experimentLabel,
          sourceDir = sourceDir
        )
        list(
          sessionFilename = "metab_filtered_session_data_demo.rds",
          sessionFilepath = "/tmp/metab-export/metab_filtered_session_data_demo.rds"
        )
      },
      handleOutcomeFn = function(sessionFilename = NULL, sessionFilepath = NULL, error = NULL, addLogFn) {
        capture$outcome <- list(
          sessionFilename = sessionFilename,
          sessionFilepath = sessionFilepath,
          error = error
        )
        addLogFn("dispatch-success")
        invisible(list(
          status = "success",
          sessionFilename = sessionFilename,
          sessionFilepath = sessionFilepath
        ))
      }
    )
  )

  expect_false(visible$visible)
  expect_identical(
    visible$value,
    list(
      status = "success",
      sessionFilename = "metab_filtered_session_data_demo.rds",
      sessionFilepath = "/tmp/metab-export/metab_filtered_session_data_demo.rds"
    )
  )
  expect_identical(capture$resolve_paths, list(source_dir = "/tmp/source"))
  expect_identical(
    capture$workflow,
    list(
      workflowData = list(state_manager = "workflow-state"),
      normData = list(assay_names = "Plasma"),
      inputValues = list(norm_method = "TIC"),
      experimentLabel = "Demo Study",
      sourceDir = "/tmp/metab-export"
    )
  )
  expect_identical(
    capture$outcome,
    list(
      sessionFilename = "metab_filtered_session_data_demo.rds",
      sessionFilepath = "/tmp/metab-export/metab_filtered_session_data_demo.rds",
      error = NULL
    )
  )
  expect_identical(capture$log, "dispatch-success")
})

test_that("metabolomics normalization module export dispatch helper preserves error shell", {
  capture <- new.env(parent = emptyenv())
  capture$log <- character()

  visible <- withVisible(
    dispatchMetabNormExportSession(
      workflowData = list(state_manager = "workflow-state"),
      normData = list(assay_names = "Plasma"),
      inputValues = list(norm_method = "TIC"),
      experimentPaths = list(export_dir = "/tmp/metab-export"),
      experimentLabel = "Demo Study",
      addLogFn = function(message) {
        capture$log <- c(capture$log, message)
      },
      resolveSourceDirFn = function(experimentPaths) {
        capture$resolve_paths <- experimentPaths
        "/tmp/metab-export"
      },
      runWorkflowFn = function(...) {
        stop("disk full")
      },
      handleOutcomeFn = function(sessionFilename = NULL, sessionFilepath = NULL, error = NULL, addLogFn) {
        capture$outcome <- list(
          sessionFilename = sessionFilename,
          sessionFilepath = sessionFilepath,
          error = conditionMessage(error)
        )
        addLogFn("dispatch-error")
        invisible(list(
          status = "error",
          errorMessage = conditionMessage(error)
        ))
      }
    )
  )

  expect_false(visible$visible)
  expect_identical(
    visible$value,
    list(
      status = "error",
      errorMessage = "disk full"
    )
  )
  expect_identical(capture$resolve_paths, list(export_dir = "/tmp/metab-export"))
  expect_identical(
    capture$outcome,
    list(
      sessionFilename = NULL,
      sessionFilepath = NULL,
      error = "disk full"
    )
  )
  expect_identical(capture$log, "dispatch-error")
})

test_that("metabolomics normalization module export outcome helper preserves success notification tail", {
  capture <- new.env(parent = emptyenv())
  capture$log <- character()
  capture$notifications <- list()
  capture$info <- character()

  visible <- withVisible(
    handleMetabNormExportSessionOutcome(
      sessionFilename = "metab_filtered_session_data_demo.rds",
      sessionFilepath = "/tmp/metab-export/metab_filtered_session_data_demo.rds",
      addLogFn = function(message) {
        capture$log <- c(capture$log, message)
      },
      showNotificationFn = function(message, type = NULL, duration = NULL) {
        capture$notifications <- c(
          capture$notifications,
          list(list(message = message, type = type, duration = duration))
        )
        invisible(NULL)
      },
      logInfoFn = function(message) {
        capture$info <- c(capture$info, message)
      },
      logErrorFn = function(message) {
        stop(sprintf("logErrorFn should not be called: %s", message))
      }
    )
  )

  expect_false(visible$visible)
  expect_identical(
    visible$value,
    list(
      status = "success",
      sessionFilename = "metab_filtered_session_data_demo.rds",
      sessionFilepath = "/tmp/metab-export/metab_filtered_session_data_demo.rds"
    )
  )
  expect_identical(
    capture$log,
    "Exported comprehensive session data to: /tmp/metab-export/metab_filtered_session_data_demo.rds"
  )
  expect_identical(
    capture$notifications,
    list(list(
      message = paste0(
        "Session data exported successfully!\n",
        "Saved as: metab_filtered_session_data_demo.rds\n",
        "See summary file for details."
      ),
      type = "message",
      duration = 10
    ))
  )
  expect_identical(capture$info, "=== EXPORT NORMALIZED SESSION COMPLETED SUCCESSFULLY ===")
})

test_that("metabolomics normalization module export outcome helper preserves error notification tail", {
  capture <- new.env(parent = emptyenv())
  capture$log <- character()
  capture$notifications <- list()
  capture$error <- character()

  visible <- withVisible(
    handleMetabNormExportSessionOutcome(
      error = simpleError("disk full"),
      addLogFn = function(message) {
        capture$log <- c(capture$log, message)
      },
      showNotificationFn = function(message, type = NULL, duration = NULL) {
        capture$notifications <- c(
          capture$notifications,
          list(list(message = message, type = type, duration = duration))
        )
        invisible(NULL)
      },
      logInfoFn = function(message) {
        stop(sprintf("logInfoFn should not be called: %s", message))
      },
      logErrorFn = function(message) {
        capture$error <- c(capture$error, message)
      }
    )
  )

  expect_false(visible$visible)
  expect_identical(
    visible$value,
    list(
      status = "error",
      errorMessage = "disk full"
    )
  )
  expect_identical(capture$error, "*** ERROR in session export: disk full ***")
  expect_identical(capture$log, "Export error: disk full")
  expect_identical(
    capture$notifications,
    list(list(
      message = "Export error: disk full",
      type = "error",
      duration = 10
    ))
  )
})

test_that("metabolomics normalization module export observer shell preserves readiness gate and dispatch handoff", {
  capture <- new.env(parent = emptyenv())
  capture$info <- character()
  capture$req <- list()
  capture$log <- character()

  workflow_data <- list(state_manager = "state-manager")
  norm_data <- list(normalization_complete = TRUE)
  input_values <- list(export_session = 1)
  experiment_paths <- list(source_dir = "/tmp/export-source")

  visible <- withVisible(
    runMetabNormExportSessionObserverShell(
      workflowData = workflow_data,
      normData = norm_data,
      inputValues = input_values,
      experimentPaths = experiment_paths,
      experimentLabel = "Metab Experiment",
      addLogFn = function(message) {
        capture$log <- c(capture$log, message)
      },
      logInfoFn = function(message) {
        capture$info <- c(capture$info, message)
      },
      reqFn = function(value) {
        capture$req <- c(capture$req, list(value))
        invisible(value)
      },
      checkReadyFn = function(normalizationComplete) {
        capture$ready <- normalizationComplete
        TRUE
      },
      dispatchExportSessionFn = function(
        workflowData,
        normData,
        inputValues,
        experimentPaths,
        experimentLabel,
        addLogFn
      ) {
        capture$dispatch <- list(
          workflowData = workflowData,
          normData = normData,
          inputValues = inputValues,
          experimentPaths = experimentPaths,
          experimentLabel = experimentLabel
        )
        addLogFn("export-shell-log")
        list(status = "delegated", source = "export-observer-shell")
      }
    )
  )

  expect_false(visible$visible)
  expect_identical(capture$info, "=== EXPORT NORMALIZED SESSION BUTTON CLICKED ===")
  expect_identical(capture$req, list("state-manager"))
  expect_true(capture$ready)
  expect_identical(
    capture$dispatch,
    list(
      workflowData = workflow_data,
      normData = norm_data,
      inputValues = input_values,
      experimentPaths = experiment_paths,
      experimentLabel = "Metab Experiment"
    )
  )
  expect_identical(capture$log, "export-shell-log")
  expect_identical(
    visible$value,
    list(status = "delegated", source = "export-observer-shell")
  )
})

test_that("metabolomics normalization module export observer shell preserves not-ready short-circuit", {
  capture <- new.env(parent = emptyenv())
  capture$req <- list()
  capture$dispatchCalled <- FALSE

  visible <- withVisible(
    runMetabNormExportSessionObserverShell(
      workflowData = list(state_manager = "state-manager"),
      normData = list(normalization_complete = FALSE),
      inputValues = list(export_session = 1),
      experimentPaths = list(export_dir = "/tmp/export-fallback"),
      experimentLabel = "Metab Experiment",
      logInfoFn = function(message) {
        capture$info <- c(capture$info, message)
      },
      reqFn = function(value) {
        capture$req <- c(capture$req, list(value))
        invisible(value)
      },
      checkReadyFn = function(normalizationComplete) {
        capture$ready <- normalizationComplete
        FALSE
      },
      dispatchExportSessionFn = function(...) {
        capture$dispatchCalled <- TRUE
        invisible(NULL)
      }
    )
  )

  expect_false(visible$visible)
  expect_null(visible$value)
  expect_identical(capture$info, "=== EXPORT NORMALIZED SESSION BUTTON CLICKED ===")
  expect_identical(capture$req, list("state-manager"))
  expect_false(capture$ready)
  expect_false(capture$dispatchCalled)
})

test_that("metabolomics normalization module normalization pipeline shell preserves orchestration handoff", {
  capture <- new.env(parent = emptyenv())
  capture$order <- character()
  capture$req <- list()

  initial_s4 <- structure(list(stage = "initial"), class = "MetabNormPipelineShellMock")
  pre_norm_s4 <- structure(list(stage = "pre"), class = "MetabNormPipelineShellMock")
  itsd_s4 <- structure(list(stage = "itsd"), class = "MetabNormPipelineShellMock")
  log2_s4 <- structure(list(stage = "log2"), class = "MetabNormPipelineShellMock")
  between_s4 <- structure(list(stage = "between"), class = "MetabNormPipelineShellMock")
  post_norm_s4 <- structure(list(stage = "post"), class = "MetabNormPipelineShellMock")
  ruv_s4 <- structure(list(stage = "ruv"), class = "MetabNormPipelineShellMock")
  final_s4 <- structure(list(stage = "final"), class = "MetabNormPipelineShellMock")

  workflow_data <- list(
    state_manager = list(
      getState = function() {
        capture$order <- c(capture$order, "getState")
        initial_s4
      }
    )
  )
  input_values <- list(
    apply_itsd = TRUE,
    itsd_aggregation = "median",
    log_offset = 1,
    norm_method = "median",
    ruv_mode = "per_assay",
    auto_percentage_min = 0.1,
    auto_percentage_max = 0.3,
    ruv_grouping_variable = "group",
    separation_metric = "pearson",
    k_penalty_weight = 0.25,
    adaptive_k_penalty = FALSE,
    ruv_k = 2,
    ruv_percentage = 0.2
  )
  experiment_paths <- list(metabolite_qc_dir = "/tmp/metab-qc")
  norm_data <- new.env(parent = emptyenv())
  norm_data$assay_names <- c("Plasma", "Urine")
  norm_data$itsd_selections <- list(Plasma = 1:2)

  composite_generator <- function(...) {
    "unused-generator"
  }
  save_plot <- function(...) {
    "unused-save"
  }
  warn_log <- function(message) {
    capture$warn <- c(capture$warn, message)
  }
  composite_qc_state <- list(qcDir = "/tmp/metab-qc", compositeSaved = TRUE)

  visible <- withVisible(
    runMetabNormNormalizationPipelineShell(
      workflowData = workflow_data,
      inputValues = input_values,
      experimentPaths = experiment_paths,
      omicType = "metabolite",
      normData = norm_data,
      getPlotAestheticsFn = function() {
        capture$order <- c(capture$order, "aesthetics")
        list(color_var = "group", shape_var = "batch")
      },
      addLogFn = function(message) {
        capture$log <- c(capture$log, message)
      },
      reqFn = function(value) {
        capture$order <- c(capture$order, "req")
        capture$req <- c(capture$req, list(value))
        invisible(value)
      },
      generateCompositeFromFilesFn = composite_generator,
      savePlotFn = save_plot,
      logWarnFn = warn_log,
      runPreNormalizationQcStepFn = function(currentS4, totalSteps, experimentPaths, groupingVariable, shapeVariable, normData, addLogFn) {
        capture$order <- c(capture$order, "pre")
        capture$pre <- list(
          currentS4 = currentS4,
          totalSteps = totalSteps,
          experimentPaths = experimentPaths,
          groupingVariable = groupingVariable,
          shapeVariable = shapeVariable,
          normData = normData,
          addLogFn = addLogFn
        )
        list(currentS4 = pre_norm_s4)
      },
      runItsdProgressApplyShellFn = function(currentS4, totalSteps, applyItsd, itsdAggregation, itsdSelections, workflowData, normData, addLogFn) {
        capture$order <- c(capture$order, "itsd")
        capture$itsd <- list(
          currentS4 = currentS4,
          totalSteps = totalSteps,
          applyItsd = applyItsd,
          itsdAggregation = itsdAggregation,
          itsdSelections = itsdSelections,
          workflowData = workflowData,
          normData = normData,
          addLogFn = addLogFn
        )
        list(currentS4 = itsd_s4)
      },
      runLog2ProgressApplyShellFn = function(currentS4, totalSteps, logOffset, workflowData, normData, addLogFn) {
        capture$order <- c(capture$order, "log2")
        capture$log2 <- list(
          currentS4 = currentS4,
          totalSteps = totalSteps,
          logOffset = logOffset,
          workflowData = workflowData,
          normData = normData,
          addLogFn = addLogFn
        )
        list(currentS4 = log2_s4)
      },
      runBetweenSampleProgressApplyShellFn = function(currentS4, totalSteps, normMethod, workflowData, normData, addLogFn) {
        capture$order <- c(capture$order, "between")
        capture$between <- list(
          currentS4 = currentS4,
          totalSteps = totalSteps,
          normMethod = normMethod,
          workflowData = workflowData,
          normData = normData,
          addLogFn = addLogFn
        )
        list(currentS4 = between_s4)
      },
      runPostNormalizationQcStepFn = function(currentS4, experimentPaths, groupingVariable, shapeVariable, addLogFn) {
        capture$order <- c(capture$order, "post")
        capture$post <- list(
          currentS4 = currentS4,
          experimentPaths = experimentPaths,
          groupingVariable = groupingVariable,
          shapeVariable = shapeVariable,
          addLogFn = addLogFn
        )
        list(currentS4 = post_norm_s4)
      },
      runRuvProgressApplyShellFn = function(
        currentS4,
        totalSteps,
        ruvMode,
        autoPercentageMin,
        autoPercentageMax,
        ruvGroupingVariable,
        separationMetric,
        kPenaltyWeight,
        adaptiveKPenalty,
        manualK,
        manualPercentage,
        experimentPaths,
        groupingVariable,
        shapeVariable,
        workflowData,
        normData,
        addLogFn
      ) {
        capture$order <- c(capture$order, "ruv")
        capture$ruv <- list(
          currentS4 = currentS4,
          totalSteps = totalSteps,
          ruvMode = ruvMode,
          autoPercentageMin = autoPercentageMin,
          autoPercentageMax = autoPercentageMax,
          ruvGroupingVariable = ruvGroupingVariable,
          separationMetric = separationMetric,
          kPenaltyWeight = kPenaltyWeight,
          adaptiveKPenalty = adaptiveKPenalty,
          manualK = manualK,
          manualPercentage = manualPercentage,
          experimentPaths = experimentPaths,
          groupingVariable = groupingVariable,
          shapeVariable = shapeVariable,
          workflowData = workflowData,
          normData = normData,
          addLogFn = addLogFn
        )
        list(currentS4 = ruv_s4)
      },
      runCompositeQcRefreshShellFn = function(
        currentS4,
        experimentPaths,
        assayNames,
        ruvMode,
        omicType,
        normData,
        addLogFn,
        generateCompositeFromFilesFn,
        savePlotFn,
        logWarnFn
      ) {
        capture$order <- c(capture$order, "composite")
        capture$composite <- list(
          currentS4 = currentS4,
          experimentPaths = experimentPaths,
          assayNames = assayNames,
          ruvMode = ruvMode,
          omicType = omicType,
          normData = normData,
          addLogFn = addLogFn,
          generateCompositeFromFilesFn = generateCompositeFromFilesFn,
          savePlotFn = savePlotFn,
          logWarnFn = logWarnFn
        )
        list(
          currentS4 = final_s4,
          compositeQcState = composite_qc_state,
          plotRefreshTrigger = 12
        )
      }
    )
  )

  expect_false(visible$visible)
  expect_identical(
    capture$order,
    c("getState", "req", "aesthetics", "pre", "itsd", "log2", "between", "post", "ruv", "composite")
  )
  expect_identical(capture$req, list(initial_s4))
  expect_identical(capture$pre$currentS4, initial_s4)
  expect_identical(capture$pre$totalSteps, 6)
  expect_identical(capture$pre$experimentPaths, experiment_paths)
  expect_identical(capture$pre$groupingVariable, "group")
  expect_identical(capture$pre$shapeVariable, "batch")
  expect_identical(capture$pre$normData, norm_data)
  expect_identical(capture$itsd$currentS4, pre_norm_s4)
  expect_identical(capture$itsd$applyItsd, TRUE)
  expect_identical(capture$itsd$itsdAggregation, "median")
  expect_identical(capture$itsd$itsdSelections, norm_data$itsd_selections)
  expect_identical(capture$log2$currentS4, itsd_s4)
  expect_identical(capture$log2$logOffset, 1)
  expect_identical(capture$between$currentS4, log2_s4)
  expect_identical(capture$between$normMethod, "median")
  expect_identical(capture$post$currentS4, between_s4)
  expect_identical(capture$post$groupingVariable, "group")
  expect_identical(capture$post$shapeVariable, "batch")
  expect_identical(capture$ruv$currentS4, post_norm_s4)
  expect_identical(capture$ruv$ruvMode, "per_assay")
  expect_identical(capture$ruv$ruvGroupingVariable, "group")
  expect_identical(capture$ruv$separationMetric, "pearson")
  expect_identical(capture$ruv$kPenaltyWeight, 0.25)
  expect_identical(capture$ruv$manualK, 2)
  expect_identical(capture$ruv$manualPercentage, 0.2)
  expect_identical(capture$composite$currentS4, ruv_s4)
  expect_identical(capture$composite$experimentPaths, experiment_paths)
  expect_identical(capture$composite$assayNames, c("Plasma", "Urine"))
  expect_identical(capture$composite$ruvMode, "per_assay")
  expect_identical(capture$composite$omicType, "metabolite")
  expect_identical(capture$composite$normData, norm_data)
  expect_identical(capture$composite$generateCompositeFromFilesFn, composite_generator)
  expect_identical(capture$composite$savePlotFn, save_plot)
  expect_identical(capture$composite$logWarnFn, warn_log)
  expect_identical(
    visible$value,
    list(
      currentS4 = final_s4,
      totalSteps = 6,
      compositeQcState = composite_qc_state,
      plotRefreshTrigger = 12
    )
  )
})

test_that("metabolomics normalization module normalization pipeline shell preserves current-state req gate", {
  capture <- new.env(parent = emptyenv())
  capture$req <- list()
  capture$pre_called <- FALSE

  expect_error(
    runMetabNormNormalizationPipelineShell(
      workflowData = list(
        state_manager = list(
          getState = function() NULL
        )
      ),
      inputValues = list(),
      experimentPaths = list(),
      omicType = "metabolite",
      normData = list(),
      getPlotAestheticsFn = function() {
        capture$aesthetics_called <- TRUE
        list(color_var = "group", shape_var = "batch")
      },
      reqFn = function(value) {
        capture$req <- c(capture$req, list(value))
        if (is.null(value)) {
          stop("state required", call. = FALSE)
        }
        invisible(value)
      },
      runPreNormalizationQcStepFn = function(...) {
        capture$pre_called <- TRUE
        invisible(NULL)
      }
    ),
    "state required"
  )

  expect_identical(capture$req, list(NULL))
  expect_false(isTRUE(capture$aesthetics_called))
  expect_false(capture$pre_called)
})

test_that("metabolomics normalization module run-normalization observer wrapper preserves pipeline input assembly", {
  capture <- new.env(parent = emptyenv())

  add_log <- function(message) message
  show_notification <- function(...) NULL
  req_fn <- function(value) value
  with_progress <- function(...) NULL
  composite_generator <- function(...) NULL
  save_plot <- function(...) NULL
  warn_log <- function(...) NULL

  input <- list(
    apply_itsd = TRUE,
    itsd_aggregation = "median",
    log_offset = 1.5,
    norm_method = "quantile",
    ruv_mode = "manual",
    auto_percentage_min = 2,
    auto_percentage_max = 18,
    ruv_grouping_variable = "batch",
    separation_metric = "auc",
    k_penalty_weight = 0.7,
    adaptive_k_penalty = FALSE,
    ruv_k = 4,
    ruv_percentage = 9,
    color_variable = "factor1",
    shape_variable = "batch"
  )

  workflow_data <- list(state_manager = "state-manager")
  experiment_paths <- list(metabolite_qc_dir = "qc-dir")
  norm_data <- list(stage = "post-filter")

  visible <- withVisible(
    runMetabNormNormalizationObserverWrapper(
      workflowData = workflow_data,
      input = input,
      experimentPaths = experiment_paths,
      omicType = "metabolite",
      normData = norm_data,
      addLogFn = add_log,
      showNotificationFn = show_notification,
      reqFn = req_fn,
      withProgressFn = with_progress,
      getPlotAestheticsFn = function(colorVariable, shapeVariable) {
        capture$plot_args <- list(
          colorVariable = colorVariable,
          shapeVariable = shapeVariable
        )
        list(color_var = toupper(colorVariable), shape_var = toupper(shapeVariable))
      },
      runObserverShellFn = function(workflowData, addLogFn, showNotificationFn, reqFn, withProgressFn, runPipelineFn) {
        capture$observer <- list(
          workflowData = workflowData,
          addLogFn = addLogFn,
          showNotificationFn = showNotificationFn,
          reqFn = reqFn,
          withProgressFn = withProgressFn
        )
        pipeline_state <- runPipelineFn()
        invisible(list(observer = "shell", pipelineState = pipeline_state))
      },
      runPipelineShellFn = function(workflowData, inputValues, experimentPaths, omicType, normData, getPlotAestheticsFn, addLogFn, reqFn, generateCompositeFromFilesFn, savePlotFn, logWarnFn) {
        capture$pipeline <- list(
          workflowData = workflowData,
          inputValues = inputValues,
          experimentPaths = experimentPaths,
          omicType = omicType,
          normData = normData,
          addLogFn = addLogFn,
          reqFn = reqFn,
          generateCompositeFromFilesFn = generateCompositeFromFilesFn,
          savePlotFn = savePlotFn,
          logWarnFn = logWarnFn,
          aesthetics = getPlotAestheticsFn()
        )
        invisible(list(
          inputValues = inputValues,
          aesthetics = capture$pipeline$aesthetics
        ))
      },
      generateCompositeFromFilesFn = composite_generator,
      savePlotFn = save_plot,
      logWarnFn = warn_log
    )
  )

  expect_false(visible$visible)
  expect_identical(capture$observer$workflowData, workflow_data)
  expect_identical(capture$observer$addLogFn, add_log)
  expect_identical(capture$observer$showNotificationFn, show_notification)
  expect_identical(capture$observer$reqFn, req_fn)
  expect_identical(capture$observer$withProgressFn, with_progress)
  expect_identical(capture$pipeline$workflowData, workflow_data)
  expect_identical(
    capture$pipeline$inputValues,
    list(
      apply_itsd = TRUE,
      itsd_aggregation = "median",
      log_offset = 1.5,
      norm_method = "quantile",
      ruv_mode = "manual",
      auto_percentage_min = 2,
      auto_percentage_max = 18,
      ruv_grouping_variable = "batch",
      separation_metric = "auc",
      k_penalty_weight = 0.7,
      adaptive_k_penalty = FALSE,
      ruv_k = 4,
      ruv_percentage = 9
    )
  )
  expect_identical(capture$pipeline$experimentPaths, experiment_paths)
  expect_identical(capture$pipeline$omicType, "metabolite")
  expect_identical(capture$pipeline$normData, norm_data)
  expect_identical(capture$pipeline$addLogFn, add_log)
  expect_identical(capture$pipeline$reqFn, req_fn)
  expect_identical(capture$pipeline$generateCompositeFromFilesFn, composite_generator)
  expect_identical(capture$pipeline$savePlotFn, save_plot)
  expect_identical(capture$pipeline$logWarnFn, warn_log)
  expect_identical(
    capture$plot_args,
    list(colorVariable = "factor1", shapeVariable = "batch")
  )
  expect_identical(
    capture$pipeline$aesthetics,
    list(color_var = "FACTOR1", shape_var = "BATCH")
  )
  expect_identical(
    visible$value,
    list(
      observer = "shell",
      pipelineState = list(
        inputValues = capture$pipeline$inputValues,
        aesthetics = list(color_var = "FACTOR1", shape_var = "BATCH")
      )
    )
  )
})

test_that("metabolomics normalization module reset observer wrapper preserves dependency forwarding", {
  capture <- new.env(parent = emptyenv())

  add_log <- function(message) message
  show_notification <- function(...) NULL
  req_fn <- function(value) value

  workflow_data <- list(state_manager = "state-manager")
  norm_data <- list(post_filter_obj = "post-filter")

  visible <- withVisible(
    runMetabNormResetNormalizationObserverWrapper(
      workflowData = workflow_data,
      normData = norm_data,
      addLogFn = add_log,
      showNotificationFn = show_notification,
      reqFn = req_fn,
      runObserverShellFn = function(workflowData, normData, addLogFn, showNotificationFn, reqFn) {
        capture$args <- list(
          workflowData = workflowData,
          normData = normData,
          addLogFn = addLogFn,
          showNotificationFn = showNotificationFn,
          reqFn = reqFn
        )
        invisible(list(status = "success", stateSaved = TRUE))
      }
    )
  )

  expect_false(visible$visible)
  expect_identical(
    capture$args,
    list(
      workflowData = workflow_data,
      normData = norm_data,
      addLogFn = add_log,
      showNotificationFn = show_notification,
      reqFn = req_fn
    )
  )
  expect_identical(
    visible$value,
    list(status = "success", stateSaved = TRUE)
  )
})

test_that("metabolomics normalization module apply-correlation observer wrapper preserves dependency forwarding", {
  capture <- new.env(parent = emptyenv())

  add_log <- function(message) message
  show_notification <- function(...) NULL
  remove_notification <- function(id) id
  req_fn <- function(value) value
  run_observer_entry <- function(...) NULL

  workflow_data <- list(state_manager = "state-manager")
  input <- list(
    min_pearson_correlation_threshold = 0.85,
    ruv_grouping_variable = "Batch"
  )
  norm_data <- list(ruv_corrected_obj = "ruv-object")

  visible <- withVisible(
    runMetabNormApplyCorrelationObserverWrapper(
      workflowData = workflow_data,
      input = input,
      normData = norm_data,
      addLogFn = add_log,
      showNotificationFn = show_notification,
      removeNotificationFn = remove_notification,
      reqFn = req_fn,
      runObserverShellFn = function(
        workflowData,
        normData,
        threshold,
        groupingVariable,
        addLogFn,
        showNotificationFn,
        removeNotificationFn,
        reqFn,
        runObserverEntryFn
      ) {
        capture$args <- list(
          workflowData = workflowData,
          normData = normData,
          threshold = threshold,
          groupingVariable = groupingVariable,
          addLogFn = addLogFn,
          showNotificationFn = showNotificationFn,
          removeNotificationFn = removeNotificationFn,
          reqFn = reqFn,
          runObserverEntryFn = runObserverEntryFn
        )
        invisible(list(status = "success", dispatched = TRUE))
      },
      runObserverEntryFn = run_observer_entry
    )
  )

  expect_false(visible$visible)
  expect_identical(
    capture$args,
    list(
      workflowData = workflow_data,
      normData = norm_data,
      threshold = 0.85,
      groupingVariable = "Batch",
      addLogFn = add_log,
      showNotificationFn = show_notification,
      removeNotificationFn = remove_notification,
      reqFn = req_fn,
      runObserverEntryFn = run_observer_entry
    )
  )
  expect_identical(
    visible$value,
    list(status = "success", dispatched = TRUE)
  )
})

test_that("metabolomics normalization module skip-correlation observer wrapper preserves dependency forwarding", {
  capture <- new.env(parent = emptyenv())

  add_log <- function(message) message
  show_notification <- function(...) NULL
  req_fn <- function(value) value
  run_observer_entry <- function(...) NULL

  workflow_data <- list(state_manager = "state-manager")
  norm_data <- list(ruv_corrected_obj = "ruv-object")

  visible <- withVisible(
    runMetabNormSkipCorrelationObserverWrapper(
      workflowData = workflow_data,
      normData = norm_data,
      addLogFn = add_log,
      showNotificationFn = show_notification,
      reqFn = req_fn,
      runObserverShellFn = function(
        workflowData,
        normData,
        addLogFn,
        showNotificationFn,
        reqFn,
        runObserverEntryFn
      ) {
        capture$args <- list(
          workflowData = workflowData,
          normData = normData,
          addLogFn = addLogFn,
          showNotificationFn = showNotificationFn,
          reqFn = reqFn,
          runObserverEntryFn = runObserverEntryFn
        )
        invisible(list(status = "success", dispatched = TRUE))
      },
      runObserverEntryFn = run_observer_entry
    )
  )

  expect_false(visible$visible)
  expect_identical(
    capture$args,
    list(
      workflowData = workflow_data,
      normData = norm_data,
      addLogFn = add_log,
      showNotificationFn = show_notification,
      reqFn = req_fn,
      runObserverEntryFn = run_observer_entry
    )
  )
  expect_identical(
    visible$value,
    list(status = "success", dispatched = TRUE)
  )
})

test_that("metabolomics normalization module export observer wrapper preserves input assembly and dependency forwarding", {
  capture <- new.env(parent = emptyenv())

  add_log <- function(message) message
  log_info <- function(message) message
  req_fn <- function(value) value
  check_ready <- function(normalizationComplete) normalizationComplete
  dispatch_export <- function(...) NULL

  workflow_data <- list(state_manager = "state-manager")
  input <- list(
    export_session = 1,
    norm_method = "log2",
    ruv_mode = "adaptive",
    apply_itsd = TRUE,
    itsd_aggregation = "mean",
    log_offset = 0.5,
    min_pearson_correlation_threshold = 0.9,
    ruv_grouping_variable = "Batch"
  )
  norm_data <- list(normalization_complete = TRUE)
  experiment_paths <- list(source_dir = "/tmp/export-source")

  visible <- withVisible(
    runMetabNormExportSessionObserverWrapper(
      workflowData = workflow_data,
      input = input,
      normData = norm_data,
      experimentPaths = experiment_paths,
      experimentLabel = "Metab Experiment",
      addLogFn = add_log,
      logInfoFn = log_info,
      reqFn = req_fn,
      runObserverShellFn = function(
        workflowData,
        normData,
        inputValues,
        experimentPaths,
        experimentLabel,
        addLogFn,
        logInfoFn,
        reqFn,
        checkReadyFn,
        dispatchExportSessionFn
      ) {
        capture$args <- list(
          workflowData = workflowData,
          normData = normData,
          inputValues = inputValues,
          experimentPaths = experimentPaths,
          experimentLabel = experimentLabel,
          addLogFn = addLogFn,
          logInfoFn = logInfoFn,
          reqFn = reqFn,
          checkReadyFn = checkReadyFn,
          dispatchExportSessionFn = dispatchExportSessionFn
        )
        invisible(list(status = "success", dispatched = TRUE))
      },
      checkReadyFn = check_ready,
      dispatchExportSessionFn = dispatch_export
    )
  )

  expect_false(visible$visible)
  expect_identical(
    capture$args,
    list(
      workflowData = workflow_data,
      normData = norm_data,
      inputValues = list(
        export_session = 1,
        norm_method = "log2",
        ruv_mode = "adaptive",
        apply_itsd = TRUE,
        itsd_aggregation = "mean",
        log_offset = 0.5,
        min_pearson_correlation_threshold = 0.9,
        ruv_grouping_variable = "Batch"
      ),
      experimentPaths = experiment_paths,
      experimentLabel = "Metab Experiment",
      addLogFn = add_log,
      logInfoFn = log_info,
      reqFn = req_fn,
      checkReadyFn = check_ready,
      dispatchExportSessionFn = dispatch_export
    )
  )
  expect_identical(
    visible$value,
    list(status = "success", dispatched = TRUE)
  )
})

test_that("metabolomics normalization module run-normalization observer shell preserves progress shell and completion tail", {
  capture <- new.env(parent = emptyenv())
  capture$req <- list()
  capture$log <- character()
  capture$notifications <- list()

  visible <- withVisible(
    runMetabNormNormalizationObserverShell(
      workflowData = list(state_manager = "state-manager"),
      addLogFn = function(message) {
        capture$log <- c(capture$log, message)
      },
      showNotificationFn = function(message, type = NULL, duration = NULL) {
        capture$notifications <- c(
          capture$notifications,
          list(list(message = message, type = type, duration = duration))
        )
        invisible(NULL)
      },
      reqFn = function(value) {
        capture$req <- c(capture$req, list(value))
        invisible(value)
      },
      withProgressFn = function(message, value, expr, detail = NULL) {
        capture$with_progress <- list(
          message = message,
          value = value,
          detail = detail
        )
        expr
      },
      runPipelineFn = function() {
        capture$log <- c(capture$log, "pipeline-body")
        list(status = "pipeline", steps = 6)
      },
      logErrorFn = function(message) {
        stop(sprintf("logErrorFn should not be called: %s", message))
      }
    )
  )

  expect_false(visible$visible)
  expect_identical(capture$req, list("state-manager"))
  expect_identical(
    capture$with_progress,
    list(
      message = "Running normalization pipeline...",
      value = 0,
      detail = NULL
    )
  )
  expect_identical(
    capture$log,
    c(
      "Starting normalization pipeline...",
      "pipeline-body",
      "Normalization pipeline complete!"
    )
  )
  expect_identical(
    capture$notifications,
    list(list(
      message = "Normalization pipeline complete!",
      type = "message",
      duration = NULL
    ))
  )
  expect_identical(
    visible$value,
    list(
      status = "success",
      pipelineState = list(status = "pipeline", steps = 6)
    )
  )
})

test_that("metabolomics normalization module run-normalization observer shell preserves error notification tail", {
  capture <- new.env(parent = emptyenv())
  capture$req <- list()
  capture$log <- character()
  capture$notifications <- list()
  capture$error <- character()

  visible <- withVisible(
    runMetabNormNormalizationObserverShell(
      workflowData = list(state_manager = "state-manager"),
      addLogFn = function(message) {
        capture$log <- c(capture$log, message)
      },
      showNotificationFn = function(message, type = NULL, duration = NULL) {
        capture$notifications <- c(
          capture$notifications,
          list(list(message = message, type = type, duration = duration))
        )
        invisible(NULL)
      },
      reqFn = function(value) {
        capture$req <- c(capture$req, list(value))
        invisible(value)
      },
      withProgressFn = function(message, value, expr, detail = NULL) {
        capture$with_progress <- list(
          message = message,
          value = value,
          detail = detail
        )
        expr
      },
      runPipelineFn = function() {
        stop("pipeline failed", call. = FALSE)
      },
      logErrorFn = function(message) {
        capture$error <- c(capture$error, message)
        invisible(NULL)
      }
    )
  )

  expect_false(visible$visible)
  expect_identical(capture$req, list("state-manager"))
  expect_identical(
    capture$with_progress,
    list(
      message = "Running normalization pipeline...",
      value = 0,
      detail = NULL
    )
  )
  expect_identical(
    capture$log,
    c(
      "Starting normalization pipeline...",
      "ERROR: pipeline failed"
    )
  )
  expect_identical(capture$error, "Normalization pipeline error: pipeline failed")
  expect_identical(
    capture$notifications,
    list(list(
      message = "Error: pipeline failed",
      type = "error",
      duration = NULL
    ))
  )
  expect_identical(
    visible$value,
    list(
      status = "error",
      errorMessage = "pipeline failed"
    )
  )
})

test_that("metabolomics normalization module reset observer shell preserves state-save and reset tail", {
  capture <- new.env(parent = emptyenv())
  capture$req <- list()
  capture$log <- character()
  capture$notifications <- list()
  capture$save_state <- NULL

  state_manager <- new.env(parent = emptyenv())
  state_manager$saveState <- function(state_name, s4_data_object, config_object, description) {
    capture$save_state <- list(
      state_name = state_name,
      s4_data_object = s4_data_object,
      config_object = config_object,
      description = description
    )
  }

  workflow_data <- new.env(parent = emptyenv())
  workflow_data$state_manager <- state_manager
  workflow_data$config_list <- list(method = "TIC")

  post_filter_object <- structure(list(stage = "post_filter"), class = "MetabNormResetMock")
  norm_data <- new.env(parent = emptyenv())
  norm_data$post_filter_obj <- post_filter_object
  norm_data$normalization_complete <- TRUE
  norm_data$ruv_complete <- TRUE
  norm_data$correlation_filtering_complete <- TRUE
  norm_data$post_norm_obj <- "post-norm"
  norm_data$ruv_corrected_obj <- "ruv-corrected"
  norm_data$correlation_filtered_obj <- "corr-filtered"
  norm_data$ruv_optimization_results <- list(Plasma = data.frame(best_k = 2))

  visible <- withVisible(
    runMetabNormResetNormalizationObserverShell(
      workflowData = workflow_data,
      normData = norm_data,
      addLogFn = function(message) {
        capture$log <- c(capture$log, message)
      },
      showNotificationFn = function(message, type = NULL, duration = NULL) {
        capture$notifications <- c(
          capture$notifications,
          list(list(message = message, type = type, duration = duration))
        )
        invisible(NULL)
      },
      reqFn = function(value) {
        capture$req <- c(capture$req, list(value))
        invisible(value)
      }
    )
  )

  expect_false(visible$visible)
  expect_identical(capture$req, list(state_manager))
  expect_identical(
    capture$save_state,
    list(
      state_name = "metab_reset",
      s4_data_object = post_filter_object,
      config_object = workflow_data$config_list,
      description = "Reset to pre-normalization state"
    )
  )
  expect_false(norm_data$normalization_complete)
  expect_false(norm_data$ruv_complete)
  expect_false(norm_data$correlation_filtering_complete)
  expect_null(norm_data$post_norm_obj)
  expect_null(norm_data$ruv_corrected_obj)
  expect_null(norm_data$correlation_filtered_obj)
  expect_identical(norm_data$ruv_optimization_results, list())
  expect_identical(capture$log, "Reset to pre-normalization state")
  expect_identical(
    capture$notifications,
    list(list(
      message = "Reset to pre-normalization state",
      type = "message",
      duration = NULL
    ))
  )
  expect_identical(
    visible$value,
    list(
      status = "success",
      stateSaved = TRUE,
      logEntry = "Reset to pre-normalization state"
    )
  )
})

test_that("metabolomics normalization module reset observer shell preserves error notification tail", {
  capture <- new.env(parent = emptyenv())
  capture$req <- list()
  capture$log <- character()
  capture$notifications <- list()

  state_manager <- new.env(parent = emptyenv())
  state_manager$saveState <- function(...) {
    stop("save failed", call. = FALSE)
  }

  workflow_data <- new.env(parent = emptyenv())
  workflow_data$state_manager <- state_manager
  workflow_data$config_list <- list(method = "TIC")

  post_filter_object <- structure(list(stage = "post_filter"), class = "MetabNormResetMock")
  norm_data <- new.env(parent = emptyenv())
  norm_data$post_filter_obj <- post_filter_object
  norm_data$normalization_complete <- TRUE
  norm_data$ruv_complete <- TRUE
  norm_data$correlation_filtering_complete <- TRUE
  norm_data$post_norm_obj <- "post-norm"
  norm_data$ruv_corrected_obj <- "ruv-corrected"
  norm_data$correlation_filtered_obj <- "corr-filtered"
  norm_data$ruv_optimization_results <- list(Plasma = data.frame(best_k = 2))

  visible <- withVisible(
    runMetabNormResetNormalizationObserverShell(
      workflowData = workflow_data,
      normData = norm_data,
      addLogFn = function(message) {
        capture$log <- c(capture$log, message)
      },
      showNotificationFn = function(message, type = NULL, duration = NULL) {
        capture$notifications <- c(
          capture$notifications,
          list(list(message = message, type = type, duration = duration))
        )
        invisible(NULL)
      },
      reqFn = function(value) {
        capture$req <- c(capture$req, list(value))
        invisible(value)
      }
    )
  )

  expect_false(visible$visible)
  expect_identical(capture$req, list(state_manager))
  expect_identical(capture$log, "ERROR during reset: save failed")
  expect_identical(
    capture$notifications,
    list(list(
      message = "Error: save failed",
      type = "error",
      duration = NULL
    ))
  )
  expect_true(norm_data$normalization_complete)
  expect_true(norm_data$ruv_complete)
  expect_true(norm_data$correlation_filtering_complete)
  expect_identical(norm_data$post_norm_obj, "post-norm")
  expect_identical(norm_data$ruv_corrected_obj, "ruv-corrected")
  expect_identical(norm_data$correlation_filtered_obj, "corr-filtered")
  expect_identical(
    norm_data$ruv_optimization_results,
    list(Plasma = data.frame(best_k = 2))
  )
  expect_identical(
    visible$value,
    list(
      status = "error",
      errorMessage = "save failed"
    )
  )
})

test_that("metabolomics normalization module skip-correlation input helper preserves source priority", {
  post_norm_object <- structure(list(stage = "post_norm"), class = "MetabNormSkipInputMock")
  ruv_corrected_object <- structure(list(stage = "ruv_corrected"), class = "MetabNormSkipInputMock")

  expect_identical(
    resolveMetabNormSkipCorrelationInputObject(
      ruvCorrectedObject = ruv_corrected_object,
      postNormObject = post_norm_object
    ),
    ruv_corrected_object
  )

  expect_identical(
    resolveMetabNormSkipCorrelationInputObject(
      ruvCorrectedObject = NULL,
      postNormObject = post_norm_object
    ),
    post_norm_object
  )
})

test_that("metabolomics normalization module skip-correlation input helper preserves empty fallback", {
  expect_null(
    resolveMetabNormSkipCorrelationInputObject(
      ruvCorrectedObject = NULL,
      postNormObject = NULL
    )
  )
})

test_that("metabolomics normalization module skip-correlation state helper preserves state-save/status shell", {
  capture <- new.env(parent = emptyenv())
  workflow_data <- new.env(parent = emptyenv())
  current_s4 <- structure(list(stage = "post_norm"), class = "MetabNormSkipStateMock")

  workflow_data$state_manager <- list(
    saveState = function(state_name, s4_data_object, config_object, description) {
      capture$save_state <- list(
        state_name = state_name,
        s4_data_object = s4_data_object,
        config_object = config_object,
        description = description
      )
    }
  )
  workflow_data$config_list <- list(method = "TIC")
  workflow_data$tab_status <- list(
    setup_import = "complete",
    quality_control = "pending",
    normalization = "pending",
    differential_expression = "disabled"
  )

  visible <- withVisible(
    completeMetabNormSkipCorrelationState(
      workflowData = workflow_data,
      currentS4 = current_s4
    )
  )

  expect_false(visible$visible)
  expect_identical(
    capture$save_state,
    list(
      state_name = "metab_norm_complete",
      s4_data_object = current_s4,
      config_object = workflow_data$config_list,
      description = "Normalization complete (correlation filtering skipped)"
    )
  )
  expect_identical(
    visible$value,
    list(
      setup_import = "complete",
      quality_control = "complete",
      normalization = "complete",
      differential_expression = "disabled"
    )
  )
  expect_identical(workflow_data$tab_status, visible$value)
})

test_that("metabolomics normalization module skip-correlation outcome helper preserves success feedback tail", {
  capture <- new.env(parent = emptyenv())
  capture$log <- character()
  capture$notifications <- list()

  visible <- withVisible(
    handleMetabNormSkipCorrelationOutcome(
      addLogFn = function(message) {
        capture$log <- c(capture$log, message)
      },
      showNotificationFn = function(message, type = NULL, duration = NULL) {
        capture$notifications <- c(
          capture$notifications,
          list(list(message = message, type = type, duration = duration))
        )
        invisible(NULL)
      }
    )
  )

  expect_false(visible$visible)
  expect_identical(
    visible$value,
    list(
      status = "success",
      logEntry = "Correlation filtering skipped - ready for DE analysis",
      notificationMessage = "Normalization complete! Proceeding to DE analysis."
    )
  )
  expect_identical(capture$log, "Correlation filtering skipped - ready for DE analysis")
  expect_identical(
    capture$notifications,
    list(list(
      message = "Normalization complete! Proceeding to DE analysis.",
      type = "message",
      duration = NULL
    ))
  )
})

test_that("metabolomics normalization module skip-correlation dispatch helper preserves empty fallback", {
  visible <- withVisible(
    dispatchMetabNormSkipCorrelation(
      workflowData = list(state_manager = "workflow-state"),
      currentS4 = NULL,
      completeStateFn = function(...) {
        stop("completeStateFn should not be called when currentS4 is NULL")
      },
      handleOutcomeFn = function(...) {
        stop("handleOutcomeFn should not be called when currentS4 is NULL")
      }
    )
  )

  expect_false(visible$visible)
  expect_null(visible$value)
})

test_that("metabolomics normalization module skip-correlation dispatch helper preserves delegation shell", {
  capture <- new.env(parent = emptyenv())
  capture$log <- character()
  capture$notifications <- list()
  workflow_data <- list(state_manager = "workflow-state")
  current_s4 <- structure(list(stage = "ruv_corrected"), class = "MetabNormSkipDispatchMock")

  visible <- withVisible(
    dispatchMetabNormSkipCorrelation(
      workflowData = workflow_data,
      currentS4 = current_s4,
      addLogFn = function(message) {
        capture$log <- c(capture$log, message)
      },
      showNotificationFn = function(message, type = NULL, duration = NULL) {
        capture$notifications <- c(
          capture$notifications,
          list(list(message = message, type = type, duration = duration))
        )
        invisible(NULL)
      },
      completeStateFn = function(workflowData, currentS4) {
        capture$state <- list(
          workflowData = workflowData,
          currentS4 = currentS4
        )
        list(
          quality_control = "complete",
          normalization = "complete"
        )
      },
      handleOutcomeFn = function(addLogFn, showNotificationFn) {
        addLogFn("dispatch-skip-success")
        showNotificationFn("skip dispatch ready", type = "message", duration = 5)
        list(
          status = "success",
          logEntry = "dispatch-skip-success",
          notificationMessage = "skip dispatch ready"
        )
      }
    )
  )

  expect_false(visible$visible)
  expect_identical(
    capture$state,
    list(
      workflowData = workflow_data,
      currentS4 = current_s4
    )
  )
  expect_identical(capture$log, "dispatch-skip-success")
  expect_identical(
    capture$notifications,
    list(list(
      message = "skip dispatch ready",
      type = "message",
      duration = 5
    ))
  )
  expect_identical(
    visible$value,
    list(
      status = "success",
      updatedStatus = list(
        quality_control = "complete",
        normalization = "complete"
      ),
      outcome = list(
        status = "success",
        logEntry = "dispatch-skip-success",
        notificationMessage = "skip dispatch ready"
      )
    )
  )
})

test_that("metabolomics normalization module skip-correlation observer entry preserves empty-resolution handoff", {
  capture <- new.env(parent = emptyenv())

  visible <- withVisible(
    runMetabNormSkipCorrelationObserverEntry(
      workflowData = "workflow-state",
      normData = list(
        ruv_corrected_obj = "ruv-object",
        post_norm_obj = "post-object"
      ),
      resolveInputObjectFn = function(ruvCorrectedObject, postNormObject) {
        capture$resolved <- list(
          ruvCorrectedObject = ruvCorrectedObject,
          postNormObject = postNormObject
        )
        NULL
      },
      dispatchSkipCorrelationFn = function(workflowData, currentS4, addLogFn, showNotificationFn) {
        capture$dispatch <- list(
          workflowData = workflowData,
          currentS4 = currentS4,
          addLogFn = addLogFn,
          showNotificationFn = showNotificationFn
        )
        NULL
      }
    )
  )

  expect_false(visible$visible)
  expect_null(visible$value)
  expect_identical(
    capture$resolved,
    list(
      ruvCorrectedObject = "ruv-object",
      postNormObject = "post-object"
    )
  )
  expect_identical(capture$dispatch$workflowData, "workflow-state")
  expect_null(capture$dispatch$currentS4)
  expect_true(is.function(capture$dispatch$addLogFn))
  expect_true(is.function(capture$dispatch$showNotificationFn))
})

test_that("metabolomics normalization module skip-correlation observer entry preserves resolution-to-dispatch handoff", {
  capture <- new.env(parent = emptyenv())
  capture$log <- character()
  capture$notifications <- list()
  current_s4 <- structure(list(stage = "ruv_corrected"), class = "MetabNormSkipObserverMock")

  visible <- withVisible(
    runMetabNormSkipCorrelationObserverEntry(
      workflowData = "workflow-state",
      normData = list(
        ruv_corrected_obj = "ruv-object",
        post_norm_obj = "post-object"
      ),
      addLogFn = function(message) {
        capture$log <- c(capture$log, message)
      },
      showNotificationFn = function(message, type = NULL, duration = NULL) {
        capture$notifications <- c(
          capture$notifications,
          list(list(message = message, type = type, duration = duration))
        )
        invisible(NULL)
      },
      resolveInputObjectFn = function(ruvCorrectedObject, postNormObject) {
        capture$resolved <- list(
          ruvCorrectedObject = ruvCorrectedObject,
          postNormObject = postNormObject
        )
        current_s4
      },
      dispatchSkipCorrelationFn = function(workflowData, currentS4, addLogFn, showNotificationFn) {
        capture$dispatch <- list(
          workflowData = workflowData,
          currentS4 = currentS4
        )
        addLogFn("observer-entry-log")
        showNotificationFn("observer-entry-notification", type = "message", duration = 5)
        list(status = "success", source = "observer-entry")
      }
    )
  )

  expect_false(visible$visible)
  expect_identical(
    capture$resolved,
    list(
      ruvCorrectedObject = "ruv-object",
      postNormObject = "post-object"
    )
  )
  expect_identical(
    capture$dispatch,
    list(
      workflowData = "workflow-state",
      currentS4 = current_s4
    )
  )
  expect_identical(capture$log, "observer-entry-log")
  expect_identical(
    capture$notifications,
    list(list(
      message = "observer-entry-notification",
      type = "message",
      duration = 5
    ))
  )
  expect_identical(
    visible$value,
    list(status = "success", source = "observer-entry")
  )
})

test_that("metabolomics normalization module skip-correlation observer shell preserves readiness gate and entry handoff", {
  capture <- new.env(parent = emptyenv())
  capture$req <- list()
  capture$log <- character()
  capture$notifications <- list()

  workflow_data <- list(state_manager = "state-manager")
  norm_data <- list(
    ruv_complete = FALSE,
    normalization_complete = TRUE
  )

  visible <- withVisible(
    runMetabNormSkipCorrelationObserverShell(
      workflowData = workflow_data,
      normData = norm_data,
      addLogFn = function(message) {
        capture$log <- c(capture$log, message)
      },
      showNotificationFn = function(message, type = NULL, duration = NULL) {
        capture$notifications <- c(
          capture$notifications,
          list(list(message = message, type = type, duration = duration))
        )
        invisible(NULL)
      },
      reqFn = function(value) {
        capture$req <- c(capture$req, list(value))
        invisible(value)
      },
      runObserverEntryFn = function(workflowData, normData, addLogFn, showNotificationFn) {
        capture$entry <- list(
          workflowData = workflowData,
          normData = normData
        )
        addLogFn("skip-shell-log")
        showNotificationFn("skip-shell-notification", type = "message", duration = 5)

        list(status = "delegated", source = "skip-observer-shell")
      }
    )
  )

  expect_false(visible$visible)
  expect_identical(capture$req, list("state-manager", TRUE))
  expect_identical(capture$entry$workflowData, workflow_data)
  expect_identical(capture$entry$normData, norm_data)
  expect_identical(capture$log, "skip-shell-log")
  expect_identical(
    capture$notifications,
    list(list(
      message = "skip-shell-notification",
      type = "message",
      duration = 5
    ))
  )
  expect_identical(
    visible$value,
    list(status = "delegated", source = "skip-observer-shell")
  )
})

test_that("metabolomics normalization module skip-correlation observer shell preserves incomplete-state req gate", {
  capture <- new.env(parent = emptyenv())
  capture$req <- list()
  capture$entryCalled <- FALSE

  expect_error(
    runMetabNormSkipCorrelationObserverShell(
      workflowData = list(state_manager = "state-manager"),
      normData = list(
        ruv_complete = FALSE,
        normalization_complete = FALSE
      ),
      reqFn = function(value) {
        capture$req <- c(capture$req, list(value))
        if (is.null(value) || identical(value, FALSE)) {
          stop("req failed", call. = FALSE)
        }
        invisible(value)
      },
      runObserverEntryFn = function(...) {
        capture$entryCalled <- TRUE
        invisible(NULL)
      }
    ),
    "req failed"
  )

  expect_identical(capture$req, list("state-manager", FALSE))
  expect_false(capture$entryCalled)
})

test_that("metabolomics normalization module apply-correlation observer entry preserves empty-resolution setup-to-dispatch handoff", {
  capture <- new.env(parent = emptyenv())
  capture$log <- character()
  capture$notifications <- list()

  visible <- withVisible(
    runMetabNormApplyCorrelationObserverEntry(
      workflowData = "workflow-state",
      normData = list(
        ruv_corrected_obj = "ruv-object",
        post_norm_obj = "post-object"
      ),
      threshold = 0.8,
      groupingVariable = "Batch",
      addLogFn = function(message) {
        capture$log <- c(capture$log, message)
      },
      showNotificationFn = function(message, ...) {
        capture$notifications <- c(
          capture$notifications,
          list(c(list(message = message), list(...)))
        )
        invisible(NULL)
      },
      removeNotificationFn = function(id) {
        capture$removed <- c(capture$removed, id)
        invisible(NULL)
      },
      resolveInputObjectFn = function(ruvCorrectedObject, postNormObject) {
        capture$resolved <- list(
          ruvCorrectedObject = ruvCorrectedObject,
          postNormObject = postNormObject
        )
        NULL
      },
      dispatchApplyCorrelationFn = function(
        workflowData,
        normData,
        observerState,
        addLogFn,
        showNotificationFn,
        removeNotificationFn
      ) {
        capture$dispatch <- list(
          workflowData = workflowData,
          normData = normData,
          observerState = observerState,
          addLogFn = addLogFn,
          showNotificationFn = showNotificationFn,
          removeNotificationFn = removeNotificationFn
        )

        list(status = "delegated", source = "observer-entry")
      }
    )
  )

  expect_false(visible$visible)
  expect_identical(
    capture$resolved,
    list(
      ruvCorrectedObject = "ruv-object",
      postNormObject = "post-object"
    )
  )
  expect_identical(capture$log, "Applying correlation filter (threshold: 0.8 )")
  expect_identical(
    capture$notifications,
    list(list(
      message = "Applying correlation filter...",
      id = "corr_working",
      duration = NULL
    ))
  )
  expect_identical(capture$dispatch$workflowData, "workflow-state")
  expect_identical(
    capture$dispatch$normData,
    list(
      ruv_corrected_obj = "ruv-object",
      post_norm_obj = "post-object"
    )
  )
  expect_identical(
    capture$dispatch$observerState,
    list(
      currentS4 = NULL,
      threshold = 0.8,
      groupingVariable = "Batch",
      logEntry = "Applying correlation filter (threshold: 0.8 )",
      notificationId = "corr_working"
    )
  )
  expect_true(is.function(capture$dispatch$addLogFn))
  expect_true(is.function(capture$dispatch$showNotificationFn))
  expect_true(is.function(capture$dispatch$removeNotificationFn))
  expect_identical(
    visible$value,
    list(status = "delegated", source = "observer-entry")
  )
})

test_that("metabolomics normalization module apply-correlation observer entry preserves resolved setup-to-dispatch workflow", {
  capture <- new.env(parent = emptyenv())
  capture$log <- character()
  capture$notifications <- list()
  capture$removed <- character()
  current_s4 <- structure(list(stage = "ruv_corrected"), class = "MetabNormApplyObserverMock")
  norm_data <- list(
    ruv_corrected_obj = "ruv-object",
    post_norm_obj = "post-object"
  )

  visible <- withVisible(
    runMetabNormApplyCorrelationObserverEntry(
      workflowData = "workflow-state",
      normData = norm_data,
      threshold = 0.9,
      groupingVariable = "Condition",
      addLogFn = function(message) {
        capture$log <- c(capture$log, message)
      },
      showNotificationFn = function(message, ...) {
        capture$notifications <- c(
          capture$notifications,
          list(c(list(message = message), list(...)))
        )
        invisible(NULL)
      },
      removeNotificationFn = function(id) {
        capture$removed <- c(capture$removed, id)
        invisible(NULL)
      },
      resolveInputObjectFn = function(ruvCorrectedObject, postNormObject) {
        capture$resolved <- list(
          ruvCorrectedObject = ruvCorrectedObject,
          postNormObject = postNormObject
        )
        current_s4
      },
      dispatchApplyCorrelationFn = function(
        workflowData,
        normData,
        observerState,
        addLogFn,
        showNotificationFn,
        removeNotificationFn
      ) {
        capture$dispatch <- list(
          workflowData = workflowData,
          normData = normData,
          observerState = observerState
        )
        addLogFn("observer-entry-log")
        showNotificationFn("observer-entry-notification", type = "message", duration = 5)
        removeNotificationFn(observerState$notificationId)

        list(status = "success", source = "observer-entry")
      }
    )
  )

  expect_false(visible$visible)
  expect_identical(
    capture$resolved,
    list(
      ruvCorrectedObject = "ruv-object",
      postNormObject = "post-object"
    )
  )
  expect_identical(capture$dispatch$workflowData, "workflow-state")
  expect_identical(capture$dispatch$normData, norm_data)
  expect_identical(
    capture$dispatch$observerState,
    list(
      currentS4 = current_s4,
      threshold = 0.9,
      groupingVariable = "Condition",
      logEntry = "Applying correlation filter (threshold: 0.9 )",
      notificationId = "corr_working"
    )
  )
  expect_identical(
    capture$log,
    c("Applying correlation filter (threshold: 0.9 )", "observer-entry-log")
  )
  expect_identical(
    capture$notifications,
    list(
      list(
        message = "Applying correlation filter...",
        id = "corr_working",
        duration = NULL
      ),
      list(
        message = "observer-entry-notification",
        type = "message",
        duration = 5
      )
    )
  )
  expect_identical(capture$removed, "corr_working")
  expect_identical(visible$value, list(status = "success", source = "observer-entry"))
})

test_that("metabolomics normalization module apply-correlation observer shell preserves readiness gate and entry handoff", {
  capture <- new.env(parent = emptyenv())
  capture$req <- list()
  capture$log <- character()
  capture$notifications <- list()
  capture$removed <- character()

  workflow_data <- list(state_manager = "state-manager")
  norm_data <- list(
    ruv_complete = TRUE,
    normalization_complete = FALSE
  )

  visible <- withVisible(
    runMetabNormApplyCorrelationObserverShell(
      workflowData = workflow_data,
      normData = norm_data,
      threshold = 0.9,
      groupingVariable = "Condition",
      addLogFn = function(message) {
        capture$log <- c(capture$log, message)
      },
      showNotificationFn = function(message, ...) {
        capture$notifications <- c(
          capture$notifications,
          list(c(list(message = message), list(...)))
        )
        invisible(NULL)
      },
      removeNotificationFn = function(id) {
        capture$removed <- c(capture$removed, id)
        invisible(NULL)
      },
      reqFn = function(value) {
        capture$req <- c(capture$req, list(value))
        invisible(value)
      },
      runObserverEntryFn = function(
        workflowData,
        normData,
        threshold,
        groupingVariable,
        addLogFn,
        showNotificationFn,
        removeNotificationFn
      ) {
        capture$entry <- list(
          workflowData = workflowData,
          normData = normData,
          threshold = threshold,
          groupingVariable = groupingVariable
        )
        addLogFn("observer-shell-log")
        showNotificationFn("observer-shell-notification", type = "message", duration = 5)
        removeNotificationFn("corr_working")

        list(status = "delegated", source = "observer-shell")
      }
    )
  )

  expect_false(visible$visible)
  expect_identical(capture$req, list("state-manager", TRUE))
  expect_identical(capture$entry$workflowData, workflow_data)
  expect_identical(capture$entry$normData, norm_data)
  expect_identical(capture$entry$threshold, 0.9)
  expect_identical(capture$entry$groupingVariable, "Condition")
  expect_identical(capture$log, "observer-shell-log")
  expect_identical(
    capture$notifications,
    list(list(
      message = "observer-shell-notification",
      type = "message",
      duration = 5
    ))
  )
  expect_identical(capture$removed, "corr_working")
  expect_identical(visible$value, list(status = "delegated", source = "observer-shell"))
})

test_that("metabolomics normalization module apply-correlation observer shell preserves incomplete-state req gate", {
  capture <- new.env(parent = emptyenv())
  capture$req <- list()
  capture$entryCalled <- FALSE

  expect_error(
    runMetabNormApplyCorrelationObserverShell(
      workflowData = list(state_manager = "state-manager"),
      normData = list(
        ruv_complete = FALSE,
        normalization_complete = FALSE
      ),
      threshold = 0.75,
      groupingVariable = "Batch",
      reqFn = function(value) {
        capture$req <- c(capture$req, list(value))
        if (is.null(value) || identical(value, FALSE)) {
          stop("req failed", call. = FALSE)
        }
        invisible(value)
      },
      runObserverEntryFn = function(...) {
        capture$entryCalled <- TRUE
        invisible(NULL)
      }
    ),
    "req failed"
  )

  expect_identical(capture$req, list("state-manager", FALSE))
  expect_false(capture$entryCalled)
})

test_that("metabolomics normalization module apply-correlation outcome helper preserves success completion tail", {
  capture <- new.env(parent = emptyenv())
  capture$saved <- NULL
  capture$log <- character()
  capture$notifications <- list()
  capture$removed <- character()
  capture$error <- character()

  filtered_s4 <- structure(list(stage = "filtered"), class = "MetabNormApplyObserverMock")
  corr_results <- list(
    Plasma = data.frame(pearson_correlation = c(0.82, 0.91))
  )

  state_manager <- new.env(parent = emptyenv())
  state_manager$saveState <- function(...) {
    capture$saved <- list(...)
    invisible(NULL)
  }

  workflow_data <- new.env(parent = emptyenv())
  workflow_data$state_manager <- state_manager
  workflow_data$config_list <- list(config = "value")
  workflow_data$tab_status <- list(
    quality_control = "pending",
    normalization = "pending",
    differential_expression = "locked"
  )

  norm_data <- new.env(parent = emptyenv())
  norm_data$correlation_results <- NULL
  norm_data$correlation_filtered_obj <- NULL
  norm_data$correlation_filtering_complete <- FALSE

  observer_state <- list(
    threshold = 0.85,
    notificationId = "corr_working"
  )

  visible <- withVisible(
    handleMetabNormApplyCorrelationOutcome(
      workflowData = workflow_data,
      normData = norm_data,
      observerState = observer_state,
      corrResults = corr_results,
      filteredS4 = filtered_s4,
      addLogFn = function(message) {
        capture$log <- c(capture$log, message)
      },
      showNotificationFn = function(message, type = NULL, duration = NULL) {
        capture$notifications <- c(
          capture$notifications,
          list(list(message = message, type = type, duration = duration))
        )
        invisible(NULL)
      },
      removeNotificationFn = function(id) {
        capture$removed <- c(capture$removed, id)
        invisible(NULL)
      },
      logErrorFn = function(message) {
        capture$error <- c(capture$error, message)
        invisible(NULL)
      }
    )
  )

  expect_false(visible$visible)
  expect_identical(norm_data$correlation_results, corr_results)
  expect_identical(norm_data$correlation_filtered_obj, filtered_s4)
  expect_true(isTRUE(norm_data$correlation_filtering_complete))
  expect_identical(
    capture$saved,
    list(
      state_name = "metab_correlation_filtered",
      s4_data_object = filtered_s4,
      config_object = workflow_data$config_list,
      description = "Correlation filtering (threshold: 0.85 )"
    )
  )
  expect_identical(
    workflow_data$tab_status,
    list(
      quality_control = "complete",
      normalization = "complete",
      differential_expression = "locked"
    )
  )
  expect_identical(capture$log, "Correlation filtering complete")
  expect_identical(
    capture$notifications,
    list(list(
      message = "Correlation filtering complete! Ready for DE analysis.",
      type = "message",
      duration = NULL
    ))
  )
  expect_identical(capture$removed, "corr_working")
  expect_identical(capture$error, character())
  expect_identical(
    visible$value,
    list(
      status = "success",
      corrResults = corr_results,
      filteredS4 = filtered_s4,
      updatedStatus = workflow_data$tab_status
    )
  )
})

test_that("metabolomics normalization module apply-correlation outcome helper preserves error cleanup", {
  capture <- new.env(parent = emptyenv())
  capture$saved <- NULL
  capture$log <- character()
  capture$notifications <- list()
  capture$removed <- character()
  capture$error <- character()

  state_manager <- new.env(parent = emptyenv())
  state_manager$saveState <- function(...) {
    capture$saved <- list(...)
    invisible(NULL)
  }

  workflow_data <- new.env(parent = emptyenv())
  workflow_data$state_manager <- state_manager
  workflow_data$config_list <- list(config = "value")
  workflow_data$tab_status <- list(
    quality_control = "pending",
    normalization = "pending"
  )

  norm_data <- new.env(parent = emptyenv())
  norm_data$correlation_results <- NULL
  norm_data$correlation_filtered_obj <- NULL
  norm_data$correlation_filtering_complete <- FALSE

  observer_state <- list(
    threshold = 0.9,
    notificationId = "corr_working"
  )

  visible <- withVisible(
    handleMetabNormApplyCorrelationOutcome(
      workflowData = workflow_data,
      normData = norm_data,
      observerState = observer_state,
      error = simpleError("correlation exploded"),
      addLogFn = function(message) {
        capture$log <- c(capture$log, message)
      },
      showNotificationFn = function(message, type = NULL, duration = NULL) {
        capture$notifications <- c(
          capture$notifications,
          list(list(message = message, type = type, duration = duration))
        )
        invisible(NULL)
      },
      removeNotificationFn = function(id) {
        capture$removed <- c(capture$removed, id)
        invisible(NULL)
      },
      logErrorFn = function(message) {
        capture$error <- c(capture$error, message)
        invisible(NULL)
      }
    )
  )

  expect_false(visible$visible)
  expect_null(norm_data$correlation_results)
  expect_null(norm_data$correlation_filtered_obj)
  expect_false(isTRUE(norm_data$correlation_filtering_complete))
  expect_null(capture$saved)
  expect_identical(
    workflow_data$tab_status,
    list(
      quality_control = "pending",
      normalization = "pending"
    )
  )
  expect_identical(
    capture$log,
    "ERROR in correlation filtering: correlation exploded"
  )
  expect_identical(
    capture$notifications,
    list(list(
      message = "Error: correlation exploded",
      type = "error",
      duration = NULL
    ))
  )
  expect_identical(capture$removed, "corr_working")
  expect_identical(
    capture$error,
    "Correlation filtering error: correlation exploded"
  )
  expect_identical(
    visible$value,
    list(
      status = "error",
      errorMessage = "correlation exploded"
    )
  )
})

test_that("metabolomics normalization module apply-correlation workflow helper preserves req-calc-filter handoff", {
  capture <- new.env(parent = emptyenv())
  capture$req <- list()
  capture$calc <- NULL
  capture$filter <- NULL
  capture$info <- character()

  current_s4 <- structure(list(stage = "ruv_corrected"), class = "MetabNormApplyObserverMock")
  filtered_s4 <- structure(list(stage = "filtered"), class = "MetabNormApplyObserverMock")
  corr_results <- list(
    Plasma = data.frame(pearson_correlation = c(0.82, 0.91))
  )

  observer_state <- list(
    currentS4 = current_s4,
    threshold = 0.85,
    groupingVariable = "Batch"
  )

  visible <- withVisible(
    runMetabNormApplyCorrelationWorkflow(
      observerState = observer_state,
      reqFn = function(value) {
        capture$req <- c(capture$req, list(value))
        value
      },
      calculateCorrelationsFn = function(theObject, correlation_group) {
        capture$calc <- list(
          theObject = theObject,
          correlation_group = correlation_group
        )
        corr_results
      },
      filterSamplesFn = function(theObject, pearson_correlation_per_pair, min_pearson_correlation_threshold) {
        capture$filter <- list(
          theObject = theObject,
          pearson_correlation_per_pair = pearson_correlation_per_pair,
          min_pearson_correlation_threshold = min_pearson_correlation_threshold
        )
        filtered_s4
      },
      logInfoFn = function(message) {
        capture$info <- c(capture$info, message)
        invisible(NULL)
      }
    )
  )

  expect_false(visible$visible)
  expect_identical(capture$req, list(current_s4))
  expect_identical(
    capture$calc,
    list(
      theObject = current_s4,
      correlation_group = "Batch"
    )
  )
  expect_identical(
    capture$filter,
    list(
      theObject = current_s4,
      pearson_correlation_per_pair = corr_results,
      min_pearson_correlation_threshold = 0.85
    )
  )
  expect_identical(capture$info, "Calculating Pearson correlations per sample pair...")
  expect_identical(
    visible$value,
    list(
      corrResults = corr_results,
      filteredS4 = filtered_s4
    )
  )
})

test_that("metabolomics normalization module apply-correlation dispatch preserves workflow-to-outcome handoff", {
  capture <- new.env(parent = emptyenv())
  capture$workflow <- NULL
  capture$outcome <- NULL
  capture$log <- character()
  capture$notifications <- list()
  capture$removed <- character()
  capture$error <- character()

  current_s4 <- structure(list(stage = "ruv_corrected"), class = "MetabNormApplyObserverMock")
  filtered_s4 <- structure(list(stage = "filtered"), class = "MetabNormApplyObserverMock")
  corr_results <- list(
    Plasma = data.frame(pearson_correlation = c(0.82, 0.91))
  )

  workflow_data <- list(state_manager = "workflow-state")
  norm_data <- list(
    correlation_results = NULL,
    correlation_filtered_obj = NULL,
    correlation_filtering_complete = FALSE
  )
  observer_state <- list(
    currentS4 = current_s4,
    threshold = 0.85,
    groupingVariable = "Batch",
    notificationId = "corr_working"
  )

  req_fn <- function(value) value
  calculate_fn <- function(theObject, correlation_group) corr_results
  filter_fn <- function(theObject, pearson_correlation_per_pair, min_pearson_correlation_threshold) filtered_s4
  log_info_fn <- function(message) invisible(NULL)
  log_error_fn <- function(message) {
    capture$error <- c(capture$error, message)
    invisible(NULL)
  }

  visible <- withVisible(
    dispatchMetabNormApplyCorrelation(
      workflowData = workflow_data,
      normData = norm_data,
      observerState = observer_state,
      addLogFn = function(message) {
        capture$log <- c(capture$log, message)
      },
      showNotificationFn = function(message, type = NULL, duration = NULL) {
        capture$notifications <- c(
          capture$notifications,
          list(list(message = message, type = type, duration = duration))
        )
        invisible(NULL)
      },
      removeNotificationFn = function(id) {
        capture$removed <- c(capture$removed, id)
        invisible(NULL)
      },
      reqFn = req_fn,
      calculateCorrelationsFn = calculate_fn,
      filterSamplesFn = filter_fn,
      runWorkflowFn = function(observerState, reqFn, calculateCorrelationsFn, filterSamplesFn, logInfoFn) {
        capture$workflow <- list(
          observerState = observerState,
          reqFn = reqFn,
          calculateCorrelationsFn = calculateCorrelationsFn,
          filterSamplesFn = filterSamplesFn,
          logInfoFn = logInfoFn
        )
        invisible(list(
          corrResults = corr_results,
          filteredS4 = filtered_s4
        ))
      },
      handleOutcomeFn = function(
        workflowData,
        normData,
        observerState,
        corrResults = NULL,
        filteredS4 = NULL,
        error = NULL,
        addLogFn,
        showNotificationFn,
        removeNotificationFn,
        logErrorFn
      ) {
        capture$outcome <- list(
          workflowData = workflowData,
          normData = normData,
          observerState = observerState,
          corrResults = corrResults,
          filteredS4 = filteredS4,
          error = error,
          logErrorFn = logErrorFn
        )
        addLogFn("dispatch-outcome-log")
        showNotificationFn("dispatch-outcome-notification", type = "message", duration = 5)
        removeNotificationFn(observerState$notificationId)
        invisible(list(status = "success", source = "dispatch-outcome"))
      },
      logInfoFn = log_info_fn,
      logErrorFn = log_error_fn
    )
  )

  expect_false(visible$visible)
  expect_identical(
    capture$workflow,
    list(
      observerState = observer_state,
      reqFn = req_fn,
      calculateCorrelationsFn = calculate_fn,
      filterSamplesFn = filter_fn,
      logInfoFn = log_info_fn
    )
  )
  expect_identical(capture$outcome$workflowData, workflow_data)
  expect_identical(capture$outcome$normData, norm_data)
  expect_identical(capture$outcome$observerState, observer_state)
  expect_identical(capture$outcome$corrResults, corr_results)
  expect_identical(capture$outcome$filteredS4, filtered_s4)
  expect_null(capture$outcome$error)
  expect_identical(capture$outcome$logErrorFn, log_error_fn)
  expect_identical(capture$log, "dispatch-outcome-log")
  expect_identical(
    capture$notifications,
    list(list(
      message = "dispatch-outcome-notification",
      type = "message",
      duration = 5
    ))
  )
  expect_identical(capture$removed, "corr_working")
  expect_identical(capture$error, character())
  expect_identical(
    visible$value,
    list(status = "success", source = "dispatch-outcome")
  )
})

test_that("metabolomics normalization module apply-correlation dispatch preserves workflow-error-to-outcome handoff", {
  capture <- new.env(parent = emptyenv())
  capture$workflow <- NULL
  capture$outcome <- NULL
  capture$log <- character()
  capture$notifications <- list()
  capture$removed <- character()
  capture$error <- character()

  current_s4 <- structure(list(stage = "ruv_corrected"), class = "MetabNormApplyObserverMock")

  workflow_data <- list(state_manager = "workflow-state")
  norm_data <- list(
    correlation_results = NULL,
    correlation_filtered_obj = NULL,
    correlation_filtering_complete = FALSE
  )
  observer_state <- list(
    currentS4 = current_s4,
    threshold = 0.9,
    groupingVariable = "Condition",
    notificationId = "corr_working"
  )

  req_fn <- function(value) value
  calculate_fn <- function(theObject, correlation_group) NULL
  filter_fn <- function(theObject, pearson_correlation_per_pair, min_pearson_correlation_threshold) NULL
  log_info_fn <- function(message) invisible(NULL)
  log_error_fn <- function(message) {
    capture$error <- c(capture$error, message)
    invisible(NULL)
  }

  visible <- withVisible(
    dispatchMetabNormApplyCorrelation(
      workflowData = workflow_data,
      normData = norm_data,
      observerState = observer_state,
      addLogFn = function(message) {
        capture$log <- c(capture$log, message)
      },
      showNotificationFn = function(message, type = NULL, duration = NULL) {
        capture$notifications <- c(
          capture$notifications,
          list(list(message = message, type = type, duration = duration))
        )
        invisible(NULL)
      },
      removeNotificationFn = function(id) {
        capture$removed <- c(capture$removed, id)
        invisible(NULL)
      },
      reqFn = req_fn,
      calculateCorrelationsFn = calculate_fn,
      filterSamplesFn = filter_fn,
      runWorkflowFn = function(observerState, reqFn, calculateCorrelationsFn, filterSamplesFn, logInfoFn) {
        capture$workflow <- list(
          observerState = observerState,
          reqFn = reqFn,
          calculateCorrelationsFn = calculateCorrelationsFn,
          filterSamplesFn = filterSamplesFn,
          logInfoFn = logInfoFn
        )
        stop("correlation exploded")
      },
      handleOutcomeFn = function(
        workflowData,
        normData,
        observerState,
        corrResults = NULL,
        filteredS4 = NULL,
        error = NULL,
        addLogFn,
        showNotificationFn,
        removeNotificationFn,
        logErrorFn
      ) {
        capture$outcome <- list(
          workflowData = workflowData,
          normData = normData,
          observerState = observerState,
          corrResults = corrResults,
          filteredS4 = filteredS4,
          errorMessage = if (inherits(error, "condition")) conditionMessage(error) else as.character(error),
          logErrorFn = logErrorFn
        )
        addLogFn("dispatch-error-log")
        showNotificationFn("dispatch-error-notification", type = "error", duration = 5)
        removeNotificationFn(observerState$notificationId)
        logErrorFn("dispatch-error-forwarded")
        invisible(list(status = "error", errorMessage = conditionMessage(error)))
      },
      logInfoFn = log_info_fn,
      logErrorFn = log_error_fn
    )
  )

  expect_false(visible$visible)
  expect_identical(
    capture$workflow,
    list(
      observerState = observer_state,
      reqFn = req_fn,
      calculateCorrelationsFn = calculate_fn,
      filterSamplesFn = filter_fn,
      logInfoFn = log_info_fn
    )
  )
  expect_identical(capture$outcome$workflowData, workflow_data)
  expect_identical(capture$outcome$normData, norm_data)
  expect_identical(capture$outcome$observerState, observer_state)
  expect_null(capture$outcome$corrResults)
  expect_null(capture$outcome$filteredS4)
  expect_identical(capture$outcome$errorMessage, "correlation exploded")
  expect_identical(capture$outcome$logErrorFn, log_error_fn)
  expect_identical(capture$log, "dispatch-error-log")
  expect_identical(
    capture$notifications,
    list(list(
      message = "dispatch-error-notification",
      type = "error",
      duration = 5
    ))
  )
  expect_identical(capture$removed, "corr_working")
  expect_identical(capture$error, "dispatch-error-forwarded")
  expect_identical(
    visible$value,
    list(
      status = "error",
      errorMessage = "correlation exploded"
    )
  )
})

test_that("metabolomics normalization module final QC PCA helper preserves argument forwarding", {
  expected_plot <- ggplot2::ggplot(data.frame(x = 1, y = 1), ggplot2::aes(x, y)) +
    ggplot2::geom_point()
  capture <- new.env(parent = emptyenv())
  source_object <- structure(list(stage = "post_norm"), class = "MetabNormRenderMock")

  rendered_plot <- buildMetabNormFinalQcPcaPlot(
    sourceObject = source_object,
    colorVar = "Condition",
    shapeVar = "Batch",
    plotPcaFn = function(object, grouping_variable = NULL, shape_variable = NULL, title = NULL) {
      capture$object <- object
      capture$grouping_variable <- grouping_variable
      capture$shape_variable <- shape_variable
      capture$title <- title
      expected_plot
    }
  )

  expect_identical(capture$object, source_object)
  expect_identical(capture$grouping_variable, "Condition")
  expect_identical(capture$shape_variable, "Batch")
  expect_identical(capture$title, "Final QC - PCA")
  expect_identical(rendered_plot, expected_plot)
})

test_that("metabolomics normalization module final QC PCA helper preserves multi-plot wrapping", {
  plot_a <- ggplot2::ggplot(data.frame(x = 1, y = 1), ggplot2::aes(x, y)) +
    ggplot2::geom_point()
  plot_b <- ggplot2::ggplot(data.frame(x = 2, y = 2), ggplot2::aes(x, y)) +
    ggplot2::geom_point()

  wrapped_plot <- buildMetabNormFinalQcPcaPlot(
    sourceObject = structure(list(stage = "post_norm"), class = "MetabNormRenderMock"),
    plotPcaFn = function(...) {
      list(plot_a, plot_b)
    },
    wrapPlotsFn = function(plots, ncol) {
      list(plots = plots, ncol = ncol)
    }
  )

  expect_identical(wrapped_plot$ncol, 1)
  expect_identical(wrapped_plot$plots, list(plot_a, plot_b))
})

test_that("metabolomics normalization module final QC PCA helper preserves single-plot list normalization", {
  expected_plot <- ggplot2::ggplot(data.frame(x = 1, y = 1), ggplot2::aes(x, y)) +
    ggplot2::geom_point()

  rendered_plot <- buildMetabNormFinalQcPcaPlot(
    sourceObject = structure(list(stage = "post_norm"), class = "MetabNormRenderMock"),
    plotPcaFn = function(...) {
      list(expected_plot)
    },
    wrapPlotsFn = function(...) {
      stop("wrapPlotsFn should not be used for a single plot")
    }
  )

  expect_identical(rendered_plot, expected_plot)
})

test_that("metabolomics normalization module final QC PCA helper preserves empty-shape fallback", {
  fallback_plot <- buildMetabNormFinalQcPcaPlot(
    sourceObject = structure(list(stage = "post_norm"), class = "MetabNormRenderMock"),
    plotPcaFn = function(...) {
      "unexpected return shape"
    }
  )

  expect_s3_class(fallback_plot, "ggplot")
  expect_length(fallback_plot$layers, 0)
})

test_that("metabolomics normalization module final QC PCA helper preserves error fallback", {
  error_plot <- buildMetabNormFinalQcPcaPlot(
    sourceObject = structure(list(stage = "post_norm"), class = "MetabNormRenderMock"),
    plotPcaFn = function(...) {
      stop("boom")
    }
  )

  expect_s3_class(error_plot, "ggplot")
  expect_identical(ggplot2::ggplot_build(error_plot)$data[[1]]$label, "Error: boom")
})
