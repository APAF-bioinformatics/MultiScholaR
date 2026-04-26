# fidelity-coverage-compare: shared
# testthat for Proteomics Normalisation module helper contracts

skipIfMissingMultiScholaRBindings <- function(...) {
  missing <- setdiff(c(...), ls(envir = asNamespace("MultiScholaR")))
  if (length(missing) > 0) {
    testthat::skip(sprintf("Target-only extracted helper(s) not present: %s", paste(missing, collapse = ", ")))
  }
}

skipIfMissingMultiScholaRBindings(
  "buildProtNormQcImagePayload",
  "renderProtNormQcImage",
  "runProtNormNormalizationWorkflow"
)

test_that("buildProtNormQcImagePayload returns image payload when file exists", {
  img_path <- tempfile(fileext = ".png")
  writeLines("placeholder", img_path)

  payload <- buildProtNormQcImagePayload(dirname(img_path), basename(img_path), "QC plot")

  expect_equal(payload$src, img_path)
  expect_equal(payload$contentType, "image/png")
  expect_equal(payload$width, "100%")
  expect_equal(payload$height, "auto")
  expect_equal(payload$alt, "QC plot")
})

test_that("buildProtNormQcImagePayload returns placeholder when file is missing", {
  payload <- buildProtNormQcImagePayload(tempdir(), "missing.png", "QC plot")

  expect_equal(payload$src, "")
  expect_equal(payload$alt, "Plot not generated yet")
})

test_that("renderProtNormQcImage delegates through renderImage with QC payload builder", {
  norm_data <- new.env(parent = emptyenv())
  norm_data$plot_refresh_trigger <- Sys.time()
  payload_call <- NULL

  rendered <- renderProtNormQcImage(
    filename = "pre_norm_pca.png",
    altText = "PCA Post-Filtering",
    normData = norm_data,
    proteinQcDir = "/tmp/qc",
    renderImageFn = function(expr, deleteFile) {
      list(result = force(expr), deleteFile = deleteFile)
    },
    buildQcImagePayloadFn = function(proteinQcDir, filename, altText) {
      payload_call <<- list(
        proteinQcDir = proteinQcDir,
        filename = filename,
        altText = altText
      )
      list(src = file.path(proteinQcDir, filename), alt = altText)
    }
  )

  expect_equal(payload_call$proteinQcDir, "/tmp/qc")
  expect_equal(payload_call$filename, "pre_norm_pca.png")
  expect_equal(payload_call$altText, "PCA Post-Filtering")
  expect_equal(rendered$result$src, "/tmp/qc/pre_norm_pca.png")
  expect_equal(rendered$result$alt, "PCA Post-Filtering")
  expect_false(rendered$deleteFile)
})

test_that("registerProtNormQcImageOutputs assigns all QC image outputs", {
  output <- new.env(parent = emptyenv())
  norm_data <- new.env(parent = emptyenv())
  calls <- list()

  registerProtNormQcImageOutputs(
    output = output,
    normData = norm_data,
    proteinQcDir = "/tmp/qc",
    renderQcImageFn = function(filename, altText, normData, proteinQcDir) {
      calls[[length(calls) + 1]] <<- list(
        filename = filename,
        altText = altText,
        normData = normData,
        proteinQcDir = proteinQcDir
      )
      list(filename = filename, altText = altText)
    }
  )

  expect_length(calls, 12)
  expect_identical(calls[[1]]$normData, norm_data)
  expect_equal(calls[[1]]$filename, "pre_norm_pca.png")
  expect_equal(calls[[12]]$filename, "ruv_corrected_correlation.png")
  expect_equal(calls[[12]]$altText, "Correlation RUV Corrected")
  expect_equal(calls[[1]]$proteinQcDir, "/tmp/qc")
  expect_equal(output$pca_post_filtering$filename, "pre_norm_pca.png")
  expect_equal(output$correlation_ruv_corrected$filename, "ruv_corrected_correlation.png")
})

test_that("renderProtNormFilteringSummaryText delegates through renderText and summary helper", {
  norm_data <- new.env(parent = emptyenv())
  norm_data$filtering_summary_text <- "ready"
  captured_text <- NULL

  rendered <- renderProtNormFilteringSummaryText(
    normData = norm_data,
    renderTextFn = function(expr) force(expr),
    getSummaryTextFn = function(text) {
      captured_text <<- text
      paste("summary:", text)
    }
  )

  expect_equal(captured_text, "ready")
  expect_equal(rendered, "summary: ready")
})

test_that("renderProtNormFinalQcPlot falls back to placeholder text when no plot is available", {
  norm_data <- new.env(parent = emptyenv())
  norm_data$final_qc_plot <- NULL
  norm_data$final_filtering_plot <- NULL
  messages <- character()
  calls <- list()

  renderProtNormFinalQcPlot(
    normData = norm_data,
    renderPlotFn = function(expr) force(expr),
    resolveFinalQcRenderStateFn = function(finalQcPlot, finalFilteringPlot) {
      list(mode = "empty", finalQcPlot = finalQcPlot, finalFilteringPlot = finalFilteringPlot)
    },
    plotNewFn = function() {
      calls$plot_new <<- TRUE
    },
    textFn = function(x, y, labels, cex) {
      calls$text <<- list(x = x, y = y, labels = labels, cex = cex)
    },
    messageFn = function(text) {
      messages <<- c(messages, text)
    }
  )

  expect_true(isTRUE(calls$plot_new))
  expect_equal(calls$text$labels, "Apply correlation filter to generate final QC plot")
  expect_true(any(grepl("Rendering", messages, fixed = TRUE)))
})

test_that("renderProtNormRuvOptimizationTable formats highlighted rows when results exist", {
  format_call <- NULL

  rendered <- renderProtNormRuvOptimizationTable(
    normData = new.env(parent = emptyenv()),
    renderDataTableFn = function(expr) force(expr),
    prepareOptimizationResultsTableFn = function(x) {
      list(
        hasResults = TRUE,
        bestPercentage = 2,
        results = data.frame(percentage = c(1, 2, 3))
      )
    },
    datatableFn = function(data, options, rownames) {
      list(data = data, options = options, rownames = rownames)
    },
    formatStyleFn = function(table, column, target, backgroundColor) {
      format_call <<- list(
        table = table,
        column = column,
        target = target,
        backgroundColor = backgroundColor
      )
      list(formatted = TRUE, table = table)
    },
    styleEqualFn = function(levels, values) {
      list(levels = levels, values = values)
    }
  )

  expect_equal(format_call$column, "percentage")
  expect_equal(format_call$target, "row")
  expect_equal(format_call$backgroundColor$levels, 2)
  expect_equal(format_call$backgroundColor$values, "#e6f3ff")
  expect_true(isTRUE(rendered$formatted))
})

test_that("registerProtNormRenderOutputs assigns all remaining render outputs", {
  output <- new.env(parent = emptyenv())
  norm_data <- new.env(parent = emptyenv())
  calls <- list()

  registerProtNormRenderOutputs(
    output = output,
    normData = norm_data,
    ruvMode = "automatic",
    groupingVariable = "group",
    renderPostNormFilteringSummaryFn = function(normData) {
      calls$post_norm <<- normData
      "post_norm_plot"
    },
    renderFilteringSummaryTextFn = function(normData) {
      calls$filtering_text <<- normData
      "filtering_text"
    },
    renderFinalQcPlotFn = function(normData) {
      calls$final_qc <<- normData
      "final_qc"
    },
    renderRuvCanonicalCorrelationPlotFn = function(normData) {
      calls$ruv_plot <<- normData
      "ruv_plot"
    },
    renderRuvOptimizationSummaryFn = function(normData, ruvMode, groupingVariable) {
      calls$ruv_summary <<- list(
        normData = normData,
        ruvMode = ruvMode,
        groupingVariable = groupingVariable
      )
      "ruv_summary"
    },
    renderRuvOptimizationTableFn = function(normData) {
      calls$ruv_table <<- normData
      "ruv_table"
    }
  )

  expect_identical(calls$post_norm, norm_data)
  expect_identical(calls$final_qc, norm_data)
  expect_identical(calls$ruv_plot, norm_data)
  expect_identical(calls$ruv_table, norm_data)
  expect_identical(calls$ruv_summary$normData, norm_data)
  expect_equal(calls$ruv_summary$ruvMode, "automatic")
  expect_equal(calls$ruv_summary$groupingVariable, "group")
  expect_equal(output$post_norm_filtering_summary, "post_norm_plot")
  expect_equal(output$filtering_summary_text, "filtering_text")
  expect_equal(output$final_qc_plot, "final_qc")
  expect_equal(output$ruv_canonical_correlation_plot, "ruv_plot")
  expect_equal(output$ruv_optimization_summary, "ruv_summary")
  expect_equal(output$ruv_optimization_table, "ruv_table")
})

test_that("getProtNormFilteringSummaryText returns summary or default text", {
  expect_equal(
    getProtNormFilteringSummaryText("ready"),
    "ready"
  )
  expect_equal(
    getProtNormFilteringSummaryText(NULL),
    "Filtering summary will be available after normalization and RUV correction."
  )
})

test_that("buildProtNormRuvOptimizationSummary formats automatic and manual summaries", {
  optimization_result <- list(
    best_percentage = 5,
    best_k = 2,
    best_control_genes_index = c(TRUE, FALSE, TRUE),
    best_separation_score = 0.1234,
    best_composite_score = 0.5678,
    separation_metric_used = "max_difference",
    k_penalty_weight = 0.5,
    adaptive_k_penalty_used = TRUE,
    sample_size = 4
  )

  automatic_summary <- buildProtNormRuvOptimizationSummary(
    ruvOptimizationResult = optimization_result,
    ruvMode = "automatic",
    groupingVariable = "group"
  )
  manual_summary <- buildProtNormRuvOptimizationSummary(
    ruvOptimizationResult = optimization_result,
    ruvMode = "manual",
    groupingVariable = "batch"
  )

  expect_match(automatic_summary, "Automatic RUV Optimization Results", fixed = TRUE)
  expect_match(automatic_summary, "Best percentage: 5.0%", fixed = TRUE)
  expect_match(automatic_summary, "RUV grouping variable: group", fixed = TRUE)
  expect_match(manual_summary, "Manual RUV Parameters", fixed = TRUE)
  expect_match(manual_summary, "RUV grouping variable: batch", fixed = TRUE)
})

test_that("prepareProtNormOptimizationResultsTable rounds numeric columns and preserves best percentage", {
  optimization_result <- list(
    best_percentage = 10,
    optimization_results = data.frame(
      percentage = c(5, 10),
      separation_score = c(0.123456, 0.987654),
      composite_score = c(1.23456, 9.87654)
    )
  )

  table_data <- prepareProtNormOptimizationResultsTable(optimization_result)

  expect_true(table_data$hasResults)
  expect_equal(table_data$bestPercentage, 10)
  expect_equal(table_data$results$separation_score, c(0.1235, 0.9877))
  expect_equal(table_data$results$composite_score, c(1.2346, 9.8765))
})

test_that("prepareProtNormOptimizationResultsTable returns placeholder table without results", {
  table_data <- prepareProtNormOptimizationResultsTable(NULL)

  expect_false(table_data$hasResults)
  expect_equal(table_data$results$Message, "Run normalization to see optimization results")
  expect_null(table_data$bestPercentage)
})

test_that("resolveProtNormFinalQcRenderState selects the expected render mode", {
  combined_state <- resolveProtNormFinalQcRenderState(
    finalQcPlot = "pca",
    finalFilteringPlot = "filter"
  )
  filtering_only_state <- resolveProtNormFinalQcRenderState(
    finalQcPlot = NULL,
    finalFilteringPlot = "filter"
  )
  pca_only_state <- resolveProtNormFinalQcRenderState(
    finalQcPlot = "pca",
    finalFilteringPlot = NULL
  )
  placeholder_state <- resolveProtNormFinalQcRenderState(
    finalQcPlot = NULL,
    finalFilteringPlot = NULL
  )

  expect_equal(combined_state$mode, "combined")
  expect_equal(combined_state$finalQcPlot, "pca")
  expect_equal(combined_state$finalFilteringPlot, "filter")
  expect_equal(filtering_only_state$mode, "filtering_only")
  expect_equal(filtering_only_state$finalFilteringPlot, "filter")
  expect_equal(pca_only_state$mode, "pca_only")
  expect_equal(pca_only_state$finalQcPlot, "pca")
  expect_equal(placeholder_state$mode, "placeholder")
})

test_that("getProtNormRuvCanonicalCorrelationPlot returns plot only when available", {
  expect_null(getProtNormRuvCanonicalCorrelationPlot(NULL))
  expect_null(getProtNormRuvCanonicalCorrelationPlot(list(best_cancor_plot = NULL)))
  expect_equal(
    getProtNormRuvCanonicalCorrelationPlot(list(best_cancor_plot = "plot")),
    "plot"
  )
})

test_that("getProtNormDefaultCorrelationFilterSummaryText returns the reset default", {
  expect_equal(
    getProtNormDefaultCorrelationFilterSummaryText(),
    "No correlation filtering applied yet"
  )
})

test_that("updateProtNormDesignDrivenChoices updates plot and RUV choices from design matrix", {
  calls <- list()
  capture_update <- function(session, inputId, choices, selected) {
    calls[[length(calls) + 1]] <<- list(
      inputId = inputId,
      choices = choices,
      selected = selected
    )
  }
  messages <- character()
  capture_message <- function(text) {
    messages <<- c(messages, text)
  }
  design_matrix <- data.frame(
    group = c("A", "B"),
    batch = c("b1", "b2"),
    sample_id = c("s1", "s2"),
    stringsAsFactors = FALSE
  )

  updateProtNormDesignDrivenChoices(
    session = "session",
    designMatrix = design_matrix,
    updateSelectInputFn = capture_update,
    messageFn = capture_message
  )

  expect_length(calls, 3)
  expect_equal(calls[[1]]$inputId, "color_variable")
  expect_equal(calls[[1]]$selected, "group")
  expect_true(all(c("group", "batch", "sample_id") %in% calls[[1]]$choices))
  expect_equal(calls[[2]]$inputId, "shape_variable")
  expect_equal(calls[[2]]$selected, "group")
  expect_equal(calls[[3]]$inputId, "ruv_grouping_variable")
  expect_equal(calls[[3]]$selected, "group")
  expect_true(any(grepl("Batch column detected", messages, fixed = TRUE)))
  expect_true(any(grepl("Updated grouping variable choices", messages, fixed = TRUE)))
})

test_that("updateProtNormDesignDrivenChoices is a no-op without design matrix", {
  calls <- list()
  updateProtNormDesignDrivenChoices(
    session = "session",
    designMatrix = NULL,
    updateSelectInputFn = function(...) {
      calls[[length(calls) + 1]] <<- list(...)
    }
  )

  expect_length(calls, 0)
})

test_that("getProtNormPlotAesthetics and getProtNormRuvGroupingVariable apply defaults", {
  expect_equal(
    getProtNormPlotAesthetics(NULL, "") ,
    list(color_var = "group", shape_var = "group")
  )
  expect_equal(
    getProtNormPlotAesthetics("batch", "factor1"),
    list(color_var = "batch", shape_var = "factor1")
  )
  expect_equal(getProtNormRuvGroupingVariable(NULL), "group")
  expect_equal(getProtNormRuvGroupingVariable(""), "group")
  expect_equal(getProtNormRuvGroupingVariable("batch"), "batch")
})

test_that("buildProtNormLabelPlot and buildProtNormTitlePlot return ggplot objects", {
  label_plot <- buildProtNormLabelPlot("a)")
  title_plot <- buildProtNormTitlePlot("Pre-Normalisation")

  expect_s3_class(label_plot, "ggplot")
  expect_s3_class(title_plot, "ggplot")
})

test_that("loadProtNormImageAsPlot returns ggplot placeholder for missing files", {
  plot_obj <- loadProtNormImageAsPlot(file.path(tempdir(), "missing.png"))

  expect_s3_class(plot_obj, "ggplot")
})

test_that("generateProtNormCompositeFromFiles returns a composite payload for simple input", {
  skip_if_not_installed("patchwork")
  skip_if_not_installed("png")
  skip_if_not_installed("ggplot2")

  img_path <- tempfile(fileext = ".png")
  png::writePNG(array(1, dim = c(1, 1, 4)), img_path)

  composite <- generateProtNormCompositeFromFiles(
    plotFiles = c(img_path),
    ncol = 1,
    rowLabels = list(pca = "a)"),
    columnLabels = "Pre"
  )

  expect_type(composite, "list")
  expect_true(all(c("plot", "width", "height") %in% names(composite)))
  expect_s3_class(composite$plot, "patchwork")
  expect_gt(composite$width, 0)
  expect_gt(composite$height, 0)
})

test_that("shouldProtNormAutoGeneratePreQc matches the expected trigger conditions", {
  expect_true(shouldProtNormAutoGeneratePreQc("normalization", "protein_replicate_filtered", FALSE))
  expect_true(shouldProtNormAutoGeneratePreQc("normalization", "normalized", FALSE))
  expect_false(shouldProtNormAutoGeneratePreQc("qc", "protein_replicate_filtered", FALSE))
  expect_false(shouldProtNormAutoGeneratePreQc("normalization", "intensity_filtered", FALSE))
  expect_false(shouldProtNormAutoGeneratePreQc("normalization", "protein_replicate_filtered", TRUE))
})

test_that("notifyProtNormNormalizationPrereqWarning logs and notifies consistently", {
  messages <- character()
  notification <- NULL

  result <- notifyProtNormNormalizationPrereqWarning(
    currentState = "intensity_filtered",
    showNotificationFn = function(message, type, duration) {
      notification <<- list(message = message, type = type, duration = duration)
    },
    messageFn = function(text) {
      messages <<- c(messages, text)
    }
  )

  expect_false(isTRUE(result))
  expect_true(any(grepl("intensity_filtered", messages, fixed = TRUE)))
  expect_equal(
    notification$message,
    "Please complete the Quality Control filtering steps before accessing normalization."
  )
  expect_equal(notification$type, "warning")
  expect_equal(notification$duration, 5)
})

test_that("handleProtNormPreQcGenerationError logs and notifies consistently", {
  messages <- character()
  notification <- NULL

  result <- handleProtNormPreQcGenerationError(
    error = structure(list(message = "boom"), class = c("simpleError", "error", "condition")),
    showNotificationFn = function(message, type, duration) {
      notification <<- list(message = message, type = type, duration = duration)
    },
    messageFn = function(text) {
      messages <<- c(messages, text)
    }
  )

  expect_false(isTRUE(result))
  expect_true(any(grepl("boom", messages, fixed = TRUE)))
  expect_equal(notification$message, "Error generating pre-normalization QC: boom")
  expect_equal(notification$type, "error")
  expect_equal(notification$duration, 10)
})

test_that("runProtNormTabEntryWorkflow auto-generates pre-normalization QC when eligible", {
  workflow_data <- new.env(parent = emptyenv())
  workflow_data$state_manager <- new.env(parent = emptyenv())
  workflow_data$state_manager$current_state <- "protein_replicate_filtered"
  norm_data <- new.env(parent = emptyenv())
  norm_data$pre_norm_qc_generated <- FALSE
  generated <- FALSE
  messages <- character()

  result <- runProtNormTabEntryWorkflow(
    selectedTab = "normalization",
    workflowData = workflow_data,
    normData = norm_data,
    generatePreNormalizationQcFn = function() {
      generated <<- TRUE
    },
    messageFn = function(text) {
      messages <<- c(messages, text)
    }
  )

  expect_true(isTRUE(result))
  expect_true(generated)
  expect_true(isTRUE(norm_data$pre_norm_qc_generated))
  expect_true(any(grepl("AUTO-TRIGGERING PRE-NORMALIZATION QC", messages, fixed = TRUE)))
})

test_that("regenerateProtNormQcForAestheticChange runs only the available regeneration steps", {
  calls <- character()
  capture_message <- function(text) {
    calls <<- c(calls, text)
  }
  pre_fn <- function() {
    calls <<- c(calls, "pre")
  }
  post_fn <- function(x) {
    calls <<- c(calls, paste("post", x))
  }
  ruv_fn <- function(x) {
    calls <<- c(calls, paste("ruv", x))
  }
  norm_data <- list(
    pre_norm_qc_generated = TRUE,
    normalization_complete = TRUE,
    normalized_protein_obj = "norm_obj",
    ruv_complete = TRUE,
    ruv_normalized_obj = "ruv_obj"
  )

  regenerateProtNormQcForAestheticChange(
    normData = norm_data,
    generatePreNormalizationQcFn = pre_fn,
    generatePostNormalizationQcFn = post_fn,
    generateRuvCorrectedQcFn = ruv_fn,
    messageFn = capture_message
  )

  expect_true("pre" %in% calls)
  expect_true("post norm_obj" %in% calls)
  expect_true("ruv ruv_obj" %in% calls)
})

test_that("initializeProtNormQcPlotPaths preserves existing paths or creates defaults", {
  existing_paths <- list(existing = TRUE)

  expect_identical(initializeProtNormQcPlotPaths(existing_paths), existing_paths)

  default_paths <- initializeProtNormQcPlotPaths(NULL)
  expect_true(all(c("post_filtering", "post_normalization", "ruv_corrected") %in% names(default_paths)))
})

test_that("saveProtNormQcPlotArtifact returns NULL without qc dir and saves when available", {
  saved <- list()
  save_stub <- function(filename, plot, width, height, dpi) {
    saved <<- list(filename = filename, plot = plot, width = width, height = height, dpi = dpi)
  }

  expect_null(saveProtNormQcPlotArtifact(NULL, "plot.png", "plot", 8, 6, savePlotFn = save_stub))

  artifact_path <- saveProtNormQcPlotArtifact(tempdir(), "plot.png", "plot", 8, 6, savePlotFn = save_stub)
  expect_equal(artifact_path, file.path(tempdir(), "plot.png"))
  expect_equal(saved$filename, artifact_path)
  expect_equal(saved$plot, "plot")
  expect_equal(saved$width, 8)
  expect_equal(saved$height, 6)
  expect_equal(saved$dpi, 150)
})

test_that("recordProtNormQcPlotPath initializes and stores stage paths", {
  paths <- recordProtNormQcPlotPath(NULL, "post_filtering", "pca", "file.png")

  expect_equal(paths$post_filtering$pca, "file.png")
  expect_true("post_normalization" %in% names(paths))
})

test_that("buildProtNormDensityPlot uses PCA and box helpers with the expected arguments", {
  pca_calls <- list()
  box_calls <- list()
  pca_stub <- function(...) {
    pca_calls <<- list(...)
    "pca_plot"
  }
  box_stub <- function(...) {
    box_calls <<- list(...)
    "density_plot"
  }
  aesthetics <- list(color_var = "group", shape_var = "batch")

  density_plot <- buildProtNormDensityPlot(
    s4Object = "s4",
    aesthetics = aesthetics,
    plotPcaFn = pca_stub,
    plotPcaBoxFn = box_stub
  )

  expect_equal(density_plot, "density_plot")
  expect_equal(pca_calls[[1]], "s4")
  expect_equal(pca_calls$grouping_variable, "group")
  expect_equal(pca_calls$shape_variable, "batch")
  expect_equal(box_calls[[1]], "pca_plot")
  expect_equal(box_calls$grouping_variable, "group")
  expect_true(box_calls$show_legend)
})

test_that("buildProtNormCorrelationPlot runs with and without progress wrappers", {
  pearson_calls <- list()
  pearson_stub <- function(...) {
    pearson_calls <<- list(...)
    "corr_plot"
  }

  plain_plot <- buildProtNormCorrelationPlot(
    s4Object = "s4",
    colorVar = "group",
    isRunningFn = function() FALSE,
    plotPearsonFn = pearson_stub
  )
  expect_equal(plain_plot, "corr_plot")
  expect_equal(pearson_calls[[1]], "s4")
  expect_equal(pearson_calls$correlation_group, "group")
  expect_true(pearson_calls$exclude_pool_samples)

  progress_calls <- character()
  pearson_calls <- list()
  progress_plot <- buildProtNormCorrelationPlot(
    s4Object = "s4",
    colorVar = "batch",
    isRunningFn = function() TRUE,
    withProgressFn = function(message, detail, value, expr) {
      progress_calls <<- c(progress_calls, message, detail, as.character(value))
      force(expr)
    },
    incProgressFn = function(value) {
      progress_calls <<- c(progress_calls, paste("inc", value))
    },
    plotPearsonFn = pearson_stub
  )
  expect_equal(progress_plot, "corr_plot")
  expect_equal(pearson_calls$correlation_group, "batch")
  expect_true(any(grepl("Generating Pearson correlation plot", progress_calls, fixed = TRUE)))
  expect_true(any(grepl("inc 0.5", progress_calls, fixed = TRUE)))
})

test_that("resolveProtNormQcStateObject returns state info and errors on missing object", {
  state_manager <- new.env(parent = emptyenv())
  state_manager$current_state <- "normalized"
  state_manager$getState <- function(state) {
    expect_equal(state, "normalized")
    "s4_obj"
  }

  resolved <- resolveProtNormQcStateObject(
    stateManager = state_manager,
    reqFn = function(x) x,
    messageFn = function(...) NULL
  )
  expect_equal(resolved$currentState, "normalized")
  expect_equal(resolved$currentS4, "s4_obj")

  state_manager$getState <- function(state) NULL
  expect_error(
    resolveProtNormQcStateObject(
      stateManager = state_manager,
      reqFn = function(x) x,
      messageFn = function(...) NULL
    ),
    "No S4 object available for QC plot generation"
  )
})

test_that("generateProtNormPostNormalizationQcArtifacts records post-normalization outputs", {
  saved_files <- character()
  save_stub <- function(qcDir, filename, plotObject, width, height, dpi = 150, savePlotFn = NULL) {
    saved_files <<- c(saved_files, filename)
    file.path(qcDir, filename)
  }
  record_stub <- function(qcPlotPaths, stage, plotType, path) {
    qcPlotPaths <- initializeProtNormQcPlotPaths(qcPlotPaths)
    qcPlotPaths[[stage]][[plotType]] <- path
    qcPlotPaths
  }

  result <- generateProtNormPostNormalizationQcArtifacts(
    normalizedS4 = "norm_obj",
    qcDir = tempdir(),
    aesthetics = list(color_var = "group", shape_var = "batch"),
    qcPlotPaths = NULL,
    messageFn = function(...) NULL,
    gcFn = function(...) NULL,
    plotPcaFn = function(...) "pca",
    plotRleFn = function(...) "rle",
    buildDensityFn = function(...) "density",
    buildCorrelationFn = function(...) "correlation",
    saveArtifactFn = save_stub,
    recordPathFn = record_stub
  )

  expect_equal(
    saved_files,
    c("post_norm_pca.png", "post_norm_rle.png", "post_norm_density.png", "post_norm_correlation.png")
  )
  expect_equal(result$post_normalization$pca, file.path(tempdir(), "post_norm_pca.png"))
  expect_equal(result$post_normalization$correlation, file.path(tempdir(), "post_norm_correlation.png"))
})

test_that("generateProtNormRuvCorrectedQcArtifacts records ruv-corrected outputs", {
  saved_files <- character()
  save_stub <- function(qcDir, filename, plotObject, width, height, dpi = 150, savePlotFn = NULL) {
    saved_files <<- c(saved_files, filename)
    file.path(qcDir, filename)
  }
  record_stub <- function(qcPlotPaths, stage, plotType, path) {
    qcPlotPaths <- initializeProtNormQcPlotPaths(qcPlotPaths)
    qcPlotPaths[[stage]][[plotType]] <- path
    qcPlotPaths
  }

  result <- generateProtNormRuvCorrectedQcArtifacts(
    ruvCorrectedS4 = "ruv_obj",
    qcDir = tempdir(),
    aesthetics = list(color_var = "group", shape_var = "batch"),
    qcPlotPaths = NULL,
    messageFn = function(...) NULL,
    gcFn = function(...) NULL,
    plotPcaFn = function(...) "pca",
    plotRleFn = function(...) "rle",
    buildDensityFn = function(...) "density",
    buildCorrelationFn = function(...) "correlation",
    saveArtifactFn = save_stub,
    recordPathFn = record_stub
  )

  expect_equal(
    saved_files,
    c("ruv_corrected_pca.png", "ruv_corrected_rle.png", "ruv_corrected_density.png", "ruv_corrected_correlation.png")
  )
  expect_equal(result$ruv_corrected$density, file.path(tempdir(), "ruv_corrected_density.png"))
})

test_that("generateProtNormPreNormalizationQcArtifacts records pre-normalization outputs", {
  if (!methods::isClass("MockProteinQcData")) {
    methods::setClass(
      "MockProteinQcData",
      slots = c(
        protein_quant_table = "data.frame",
        protein_id_column = "character"
      )
    )
  }

  mock_s4 <- methods::new(
    "MockProteinQcData",
    protein_quant_table = data.frame(Protein.Ids = c("p1", "p2"), S1 = c(1, 2), S2 = c(3, 4)),
    protein_id_column = "Protein.Ids"
  )
  state_manager <- new.env(parent = emptyenv())
  state_manager$current_state <- "protein_replicate_filtered"
  state_manager$getState <- function(state) mock_s4

  saved_files <- character()
  save_stub <- function(qcDir, filename, plotObject, width, height, dpi = 150, savePlotFn = NULL) {
    saved_files <<- c(saved_files, filename)
    file.path(qcDir, filename)
  }
  record_stub <- function(qcPlotPaths, stage, plotType, path) {
    qcPlotPaths <- initializeProtNormQcPlotPaths(qcPlotPaths)
    qcPlotPaths[[stage]][[plotType]] <- path
    qcPlotPaths
  }

  result <- generateProtNormPreNormalizationQcArtifacts(
    stateManager = state_manager,
    qcDir = tempdir(),
    aesthetics = list(color_var = "group", shape_var = "batch"),
    qcPlotPaths = NULL,
    reqFn = function(x) x,
    messageFn = function(...) NULL,
    gcFn = function(...) NULL,
    plotPcaFn = function(...) "pca",
    plotRleFn = function(...) "rle",
    buildDensityFn = function(...) "density",
    buildCorrelationFn = function(...) "correlation",
    saveArtifactFn = save_stub,
    recordPathFn = record_stub
  )

  expect_equal(
    saved_files,
    c("pre_norm_pca.png", "pre_norm_rle.png", "pre_norm_density.png", "pre_norm_correlation.png")
  )
  expect_equal(result$post_filtering$rle, file.path(tempdir(), "pre_norm_rle.png"))
})

test_that("prepareProtNormNormalizationRun returns current state and seeds the QC composite", {
  state_manager <- new.env(parent = emptyenv())
  state_manager$current_state <- "protein_replicate_filtered"
  state_manager$getState <- function(state) {
    expect_equal(state, "protein_replicate_filtered")
    "s4_obj"
  }

  norm_data <- new.env(parent = emptyenv())

  run_context <- prepareProtNormNormalizationRun(
    stateManager = state_manager,
    normData = norm_data,
    reqFn = function(x) x,
    initialiseGridFn = function() "grid_obj",
    messageFn = function(...) NULL
  )

  expect_equal(run_context$currentState, "protein_replicate_filtered")
  expect_equal(run_context$currentS4, "s4_obj")
  expect_equal(norm_data$QC_composite_figure, "grid_obj")
})

test_that("runProtNormBetweenSamplesStep normalizes data and saves the pre-RUV matrix", {
  if (!methods::isClass("MockNormStepData")) {
    methods::setClass(
      "MockNormStepData",
      slots = c(protein_quant_table = "data.frame")
    )
  }

  current_s4 <- methods::new(
    "MockNormStepData",
    protein_quant_table = data.frame(Protein.Ids = "p1", S1 = 1)
  )
  normalized_s4 <- methods::new(
    "MockNormStepData",
    protein_quant_table = data.frame(Protein.Ids = "p1", S1 = 2)
  )
  norm_data <- new.env(parent = emptyenv())
  assigned_config <- NULL
  checkpoints <- list()
  saved_matrix <- list()

  result <- runProtNormBetweenSamplesStep(
    currentS4 = current_s4,
    normMethod = "quantile",
    normData = norm_data,
    proteinQcDir = tempdir(),
    normaliseBetweenSamplesFn = function(object, normalisation_method) {
      expect_identical(object, current_s4)
      expect_equal(normalisation_method, "quantile")
      normalized_s4
    },
    captureCheckpointFn = function(object, checkpointId, label) {
      checkpoints <<- list(object = object, checkpointId = checkpointId, label = label)
    },
    existsFn = function(name, envir) {
      expect_equal(name, "config_list")
      TRUE
    },
    getFn = function(name, envir) {
      list(normaliseBetweenSamples = list(method = "cyclicloess"))
    },
    assignFn = function(name, value, envir) {
      assigned_config <<- value
    },
    saveMatrixFn = function(data, path) {
      saved_matrix <<- list(data = data, path = path)
    },
    messageFn = function(...) NULL
  )

  expect_identical(result, normalized_s4)
  expect_identical(norm_data$normalized_protein_obj, normalized_s4)
  expect_identical(checkpoints$object, normalized_s4)
  expect_equal(checkpoints$checkpointId, "cp05")
  expect_equal(checkpoints$label, "normalised")
  expect_equal(assigned_config$normaliseBetweenSamples$method, "quantile")
  expect_equal(saved_matrix$data, normalized_s4@protein_quant_table)
  expect_equal(
    saved_matrix$path,
    file.path(tempdir(), "normalized_protein_matrix_pre_ruv.tsv")
  )
})

test_that("runProtNormPostNormalizationQcStep marks normalization complete after QC generation", {
  norm_data <- new.env(parent = emptyenv())
  calls <- character()

  runProtNormPostNormalizationQcStep(
    normalizedS4 = "norm_obj",
    normData = norm_data,
    generatePostNormalizationQcFn = function(object) {
      calls <<- c(calls, object)
    },
    messageFn = function(...) NULL
  )

  expect_equal(calls, "norm_obj")
  expect_true(isTRUE(norm_data$normalization_complete))
})

test_that("buildProtNormSkippedRuvResult returns the expected skip payload", {
  result <- buildProtNormSkippedRuvResult()

  expect_true(isTRUE(result$ruv_skipped))
  expect_true(is.na(result$best_percentage))
  expect_true(is.na(result$best_k))
  expect_equal(result$skip_reason, "User selected skip due to dataset constraints")
})

test_that("applyProtNormSkippedRuvState stores skip results and saves normalized state", {
  norm_data <- new.env(parent = emptyenv())
  workflow_data <- new.env(parent = emptyenv())
  saved_state <- NULL
  workflow_data$state_manager <- new.env(parent = emptyenv())
  workflow_data$state_manager$saveState <- function(...) {
    saved_state <<- list(...)
  }
  saved_rds <- list()
  skip_result <- list(
    best_percentage = NA,
    best_k = NA,
    best_control_genes_index = NA,
    best_cancor_plot = NULL,
    ruv_skipped = TRUE,
    skip_reason = "skip"
  )

  applyProtNormSkippedRuvState(
    normalizedS4 = "norm_obj",
    normMethod = "cyclicloess",
    normData = norm_data,
    workflowData = workflow_data,
    sourceDir = tempdir(),
    skipResult = skip_result,
    saveRdsFn = function(object, path) {
      saved_rds <<- list(object = object, path = path)
    },
    messageFn = function(...) NULL
  )

  expect_equal(norm_data$ruv_normalized_obj, "norm_obj")
  expect_true(isTRUE(norm_data$ruv_complete))
  expect_true(is.na(norm_data$best_k))
  expect_true(is.na(norm_data$control_genes_index))
  expect_identical(norm_data$ruv_optimization_result, skip_result)
  expect_identical(workflow_data$ruv_optimization_result, skip_result)
  expect_identical(saved_rds$object, skip_result)
  expect_equal(saved_rds$path, file.path(tempdir(), "ruv_optimization_results.RDS"))
  expect_equal(saved_state$state_name, "normalized")
  expect_equal(saved_state$s4_data_object, "norm_obj")
  expect_equal(saved_state$config_object$norm_method, "cyclicloess")
  expect_equal(saved_state$config_object$ruv_mode, "skip")
  expect_false(saved_state$config_object$ruv_applied)
})

test_that("buildProtNormManualRuvResult returns the manual-mode result payload", {
  result <- buildProtNormManualRuvResult(
    percentageAsNegCtrl = 10,
    ruvK = 3,
    controlGenesIndex = c(TRUE, FALSE, TRUE),
    cancorPlot = "plot_obj"
  )

  expect_equal(result$best_percentage, 10)
  expect_equal(result$best_k, 3)
  expect_identical(result$best_control_genes_index, c(TRUE, FALSE, TRUE))
  expect_equal(result$best_cancor_plot, "plot_obj")
  expect_equal(result$optimization_results$num_controls, 2)
  expect_equal(result$separation_metric_used, "manual")
})

test_that("updateProtNormRuvAuditTrail updates config_list when available", {
  assigned_config <- NULL

  updateProtNormRuvAuditTrail(
    ruvK = 4,
    controlGenesIndex = c(TRUE, FALSE),
    percentageAsNegCtrl = 12,
    modeLabel = "automatic",
    existsFn = function(name, envir) TRUE,
    getFn = function(name, envir) list(existing = TRUE),
    assignFn = function(name, value, envir) {
      assigned_config <<- value
    },
    updateRuvParametersFn = function(config, ruv_k, control_genes_index, percentage_as_neg_ctrl) {
      config$ruv_k <- ruv_k
      config$controls <- control_genes_index
      config$percentage <- percentage_as_neg_ctrl
      config
    },
    messageFn = function(...) NULL
  )

  expect_equal(assigned_config$ruv_k, 4)
  expect_identical(assigned_config$controls, c(TRUE, FALSE))
  expect_equal(assigned_config$percentage, 12)
})

test_that("resolveProtNormRuvParameters resolves automatic mode and persists results", {
  if (!methods::isClass("MockRuvParamData")) {
    methods::setClass(
      "MockRuvParamData",
      slots = c(protein_quant_table = "data.frame")
    )
  }

  normalized_s4 <- methods::new(
    "MockRuvParamData",
    protein_quant_table = data.frame(Protein.Ids = c("p1", "p2"), S1 = c(1, 2), S2 = c(3, 4))
  )
  norm_data <- new.env(parent = emptyenv())
  workflow_data <- new.env(parent = emptyenv())
  persisted <- NULL
  audit_calls <- list()

  result <- resolveProtNormRuvParameters(
    normalizedS4 = normalized_s4,
    input = list(
      ruv_mode = "automatic",
      auto_percentage_min = 1,
      auto_percentage_max = 3,
      separation_metric = "max_difference",
      k_penalty_weight = 0.5,
      adaptive_k_penalty = TRUE
    ),
    normData = norm_data,
    workflowData = workflow_data,
    sourceDir = tempdir(),
    getRuvGroupingVariableFn = function() "group",
    withProgressFn = function(message, detail, value, expr) force(expr),
    findBestNegCtrlPercentageFn = function(...) {
      list(
        best_percentage = 3,
        best_k = 2,
        best_control_genes_index = c(TRUE, FALSE),
        best_cancor_plot = NULL,
        optimization_results = data.frame(percentage = 1:3)
      )
    },
    updateAuditTrailFn = function(ruvK, controlGenesIndex, percentageAsNegCtrl, modeLabel, ...) {
      audit_calls <<- list(
        ruvK = ruvK,
        controlGenesIndex = controlGenesIndex,
        percentageAsNegCtrl = percentageAsNegCtrl,
        modeLabel = modeLabel
      )
    },
    persistRuvResultFn = function(ruvResult, workflowData, sourceDir, resultLabel, ...) {
      persisted <<- list(ruvResult = ruvResult, sourceDir = sourceDir, resultLabel = resultLabel)
      workflowData$ruv_optimization_result <- ruvResult
    },
    messageFn = function(...) NULL
  )

  expect_equal(result$percentageAsNegCtrl, 3)
  expect_equal(result$ruvK, 2)
  expect_identical(result$controlGenesIndex, c(TRUE, FALSE))
  expect_equal(norm_data$best_k, 2)
  expect_identical(norm_data$control_genes_index, c(TRUE, FALSE))
  expect_equal(audit_calls$modeLabel, "automatic")
  expect_equal(persisted$resultLabel, "RUV optimization results")
  expect_equal(persisted$sourceDir, tempdir())
  expect_equal(workflow_data$ruv_optimization_result$best_percentage, 3)
})

test_that("resolveProtNormRuvParameters resolves manual mode and stores manual results", {
  normalized_s4 <- "norm_obj"
  norm_data <- new.env(parent = emptyenv())
  workflow_data <- new.env(parent = emptyenv())
  persisted <- NULL

  result <- resolveProtNormRuvParameters(
    normalizedS4 = normalized_s4,
    input = list(
      ruv_mode = "manual",
      ruv_percentage = 7,
      ruv_k = 5
    ),
    normData = norm_data,
    workflowData = workflow_data,
    sourceDir = tempdir(),
    getRuvGroupingVariableFn = function() "batch",
    getNegCtrlProtAnovaFn = function(...) c(TRUE, TRUE, FALSE),
    ruvCancorFn = function(...) "cancor_plot",
    updateAuditTrailFn = function(...) NULL,
    persistRuvResultFn = function(ruvResult, workflowData, sourceDir, resultLabel, ...) {
      persisted <<- list(ruvResult = ruvResult, resultLabel = resultLabel)
      workflowData$ruv_optimization_result <- ruvResult
    },
    messageFn = function(...) NULL
  )

  expect_equal(result$percentageAsNegCtrl, 7)
  expect_equal(result$ruvK, 5)
  expect_identical(result$controlGenesIndex, c(TRUE, TRUE, FALSE))
  expect_equal(norm_data$best_k, 5)
  expect_identical(norm_data$control_genes_index, c(TRUE, TRUE, FALSE))
  expect_equal(norm_data$ruv_optimization_result$best_cancor_plot, "cancor_plot")
  expect_equal(persisted$resultLabel, "manual RUV results")
})

test_that("applyProtNormRuvCorrectionStep applies correction and captures cp06", {
  norm_data <- new.env(parent = emptyenv())
  norm_data$best_k <- 2
  norm_data$control_genes_index <- c(TRUE, FALSE)
  checkpoint <- NULL

  result <- applyProtNormRuvCorrectionStep(
    normalizedS4 = "norm_obj",
    normData = norm_data,
    getRuvGroupingVariableFn = function() "group",
    ruvIII_C_VaryingFn = function(object, ruv_grouping_variable, ruv_number_k, ctrl) {
      expect_equal(object, "norm_obj")
      expect_equal(ruv_grouping_variable, "group")
      expect_equal(ruv_number_k, 2)
      expect_identical(ctrl, c(TRUE, FALSE))
      "ruv_obj"
    },
    captureCheckpointFn = function(object, checkpointId, label) {
      checkpoint <<- list(object = object, checkpointId = checkpointId, label = label)
    },
    messageFn = function(...) NULL
  )

  expect_equal(result, "ruv_obj")
  expect_equal(norm_data$ruv_normalized_obj, "ruv_obj")
  expect_equal(checkpoint$checkpointId, "cp06")
  expect_equal(checkpoint$label, "ruv_corrected")
})

test_that("finalizeProtNormRuvCleanupStep tracks protein counts and saves normalized state", {
  norm_data <- new.env(parent = emptyenv())
  norm_data$best_k <- 6
  norm_data$ruv_optimization_result <- list(best_percentage = 8)
  workflow_data <- new.env(parent = emptyenv())
  saved_state <- NULL
  workflow_data$state_manager <- new.env(parent = emptyenv())
  workflow_data$state_manager$saveState <- function(...) {
    saved_state <<- list(...)
  }

  result <- finalizeProtNormRuvCleanupStep(
    ruvCorrectedS4 = "ruv_obj",
    input = list(
      norm_method = "cyclicloess",
      ruv_mode = "automatic",
      ruv_percentage = 5
    ),
    normData = norm_data,
    workflowData = workflow_data,
    omicType = "proteomics",
    experimentLabel = "exp1",
    removeRowsWithMissingValuesPercentFn = function(theObject) {
      expect_equal(theObject, "ruv_obj")
      structure(list(clean = TRUE), class = "clean_obj")
    },
    countDistinctProteinsFn = function(s4Object) {
      expect_true(inherits(s4Object, "clean_obj"))
      42
    },
    updateProteinFilteringFn = function(...) "filter_plot",
    messageFn = function(...) NULL
  )

  expect_true(inherits(result, "clean_obj"))
  expect_equal(workflow_data$protein_counts$after_ruv_filtering, 42)
  expect_equal(norm_data$post_norm_filtering_plot, "filter_plot")
  expect_match(norm_data$filtering_summary_text, "Post-RUV filtering: 42 proteins", fixed = TRUE)
  expect_identical(norm_data$ruv_normalized_obj, result)
  expect_equal(saved_state$state_name, "ruv_corrected")
  expect_equal(saved_state$config_object$norm_method, "cyclicloess")
  expect_equal(saved_state$config_object$ruv_k, 6)
  expect_equal(saved_state$config_object$percentage_as_neg_ctrl, 8)
})

test_that("resolveProtNormStep6QcObject prefers step5 output and falls back to norm_data", {
  norm_data <- new.env(parent = emptyenv())
  norm_data$ruv_normalized_obj <- "fallback_obj"

  expect_equal(
    resolveProtNormStep6QcObject(
      step5Object = "step5_obj",
      normData = norm_data,
      messageFn = function(...) NULL
    ),
    "step5_obj"
  )

  expect_equal(
    resolveProtNormStep6QcObject(
      step5Object = NULL,
      normData = norm_data,
      messageFn = function(...) NULL
    ),
    "fallback_obj"
  )
})

test_that("runProtNormStep6RuvQc runs RUV QC generation and saves cancor plot", {
  norm_data <- new.env(parent = emptyenv())
  norm_data$ruv_optimization_result <- list(best_cancor_plot = "plot_obj")
  norm_data$qc_plot_paths <- NULL
  calls <- list()

  runProtNormStep6RuvQc(
    ruvMode = "automatic",
    step6Object = "ruv_obj",
    normData = norm_data,
    qcDir = tempdir(),
    generateRuvCorrectedQcFn = function(object) {
      calls$generated <<- object
    },
    ggsaveFn = function(filename, plot, width, height, dpi) {
      calls$saved <<- list(filename = filename, plot = plot, width = width, height = height, dpi = dpi)
    },
    initPathsFn = function(paths) initializeProtNormQcPlotPaths(paths),
    messageFn = function(...) NULL
  )

  expect_equal(calls$generated, "ruv_obj")
  expect_equal(calls$saved$filename, file.path(tempdir(), "ruv_corrected_cancor.png"))
  expect_equal(calls$saved$plot, "plot_obj")
  expect_equal(norm_data$qc_plot_paths$ruv_corrected$cancor, file.path(tempdir(), "ruv_corrected_cancor.png"))
})

test_that("runProtNormStep6RuvQc skips generation when RUV was skipped", {
  norm_data <- new.env(parent = emptyenv())
  norm_data$ruv_optimization_result <- list(best_cancor_plot = "plot_obj")
  generated <- FALSE

  runProtNormStep6RuvQc(
    ruvMode = "skip",
    step6Object = "ruv_obj",
    normData = norm_data,
    qcDir = tempdir(),
    generateRuvCorrectedQcFn = function(object) {
      generated <<- TRUE
    },
    messageFn = function(...) NULL
  )

  expect_false(generated)
})

test_that("resolveProtNormCompositeFigureInputs returns expected skip and full layouts", {
  skip_inputs <- resolveProtNormCompositeFigureInputs("skip", tempdir())
  full_inputs <- resolveProtNormCompositeFigureInputs("automatic", tempdir())

  expect_equal(skip_inputs$ncol, 2)
  expect_length(skip_inputs$plotFiles, 8)
  expect_equal(skip_inputs$columnLabels, c("Pre-Normalisation", "Post-Normalisation"))

  expect_equal(full_inputs$ncol, 3)
  expect_length(full_inputs$plotFiles, 15)
  expect_true(any(is.na(full_inputs$plotFiles)))
  expect_equal(full_inputs$columnLabels, c("Pre-Normalisation", "Post-Normalisation", "RUV-Corrected"))
})

test_that("generateProtNormCompositeQcFigure saves the composite when generated", {
  saved <- list()

  generateProtNormCompositeQcFigure(
    ruvMode = "automatic",
    qcDir = tempdir(),
    omicType = "proteomics",
    resolveInputsFn = function(ruvMode, qcDir) {
      list(
        ncol = 3,
        plotFiles = c(file.path(qcDir, "a.png")),
        rowLabels = list(pca = c("a)")),
        columnLabels = c("Pre")
      )
    },
    generateCompositeFn = function(plotFiles, ncol, rowLabels, columnLabels) {
      list(plot = "composite_plot", width = 12, height = 8)
    },
    savePlotFn = function(plot, base_path, plot_name, formats = c("pdf", "png"), width, height, ...) {
      saved <<- list(plot = plot, base_path = base_path, plot_name = plot_name, width = width, height = height)
    },
    messageFn = function(...) NULL
  )

  expect_equal(saved$plot, "composite_plot")
  expect_equal(saved$base_path, tempdir())
  expect_equal(saved$plot_name, "proteomics_composite_QC_figure")
  expect_equal(saved$width, 12)
  expect_equal(saved$height, 8)
})

test_that("finalizeProtNormWorkflowState clears composite state and marks RUV complete", {
  norm_data <- new.env(parent = emptyenv())
  norm_data$QC_composite_figure <- "grid_obj"
  norm_data$ruv_complete <- FALSE
  gc_called <- FALSE

  finalizeProtNormWorkflowState(
    normData = norm_data,
    gcFn = function() {
      gc_called <<- TRUE
    },
    messageFn = function(...) NULL
  )

  expect_null(norm_data$QC_composite_figure)
  expect_true(norm_data$ruv_complete)
  expect_true(gc_called)
})

test_that("buildProtNormCompletionNotification returns mode-specific text", {
  expect_match(
    buildProtNormCompletionNotification("skip"),
    "RUV skipped",
    fixed = TRUE
  )
  expect_match(
    buildProtNormCompletionNotification("automatic"),
    "Normalization and RUV correction completed",
    fixed = TRUE
  )
})

test_that("handleProtNormNormalizationError logs and notifies consistently", {
  messages <- character()
  notification <- NULL

  result <- handleProtNormNormalizationError(
    error = structure(list(message = "boom"), class = c("simpleError", "error", "condition")),
    showNotificationFn = function(message, type, duration) {
      notification <<- list(message = message, type = type, duration = duration)
    },
    messageFn = function(text) {
      messages <<- c(messages, text)
    }
  )

  expect_false(isTRUE(result))
  expect_true(any(grepl("boom", messages, fixed = TRUE)))
  expect_equal(notification$message, "Error in normalization: boom")
  expect_equal(notification$type, "error")
  expect_equal(notification$duration, 10)
})

test_that("runProtNormNormalizationWorkflow sequences the automatic normalization shell", {
  workflow_data <- new.env(parent = emptyenv())
  workflow_data$state_manager <- "state_manager"
  norm_data <- new.env(parent = emptyenv())
  calls <- list()
  progress <- character()
  notification <- NULL

  result <- runProtNormNormalizationWorkflow(
    input = list(norm_method = "cyclicloess", ruv_mode = "automatic"),
    workflowData = workflow_data,
    normData = norm_data,
    experimentPaths = list(protein_qc_dir = tempdir(), source_dir = tempdir()),
    omicType = "proteomics",
    experimentLabel = "exp1",
    checkMemoryUsageFn = function(threshold_gb, context) {
      calls$memory <<- list(threshold_gb = threshold_gb, context = context)
    },
    generatePostNormalizationQcFn = function(object) {
      calls$post_qc_generated <<- object
    },
    generateRuvCorrectedQcFn = function(object) {
      calls$ruv_qc_generated <<- object
    },
    getRuvGroupingVariableFn = function() "group",
    prepareNormalizationRunFn = function(stateManager, normData) {
      calls$prepare <<- list(stateManager = stateManager, normData = normData)
      list(currentS4 = "current_obj")
    },
    runBetweenSamplesStepFn = function(currentS4, normMethod, normData, proteinQcDir) {
      calls$between <<- list(
        currentS4 = currentS4,
        normMethod = normMethod,
        normData = normData,
        proteinQcDir = proteinQcDir
      )
      "normalized_obj"
    },
    runPostNormalizationQcStepFn = function(normalizedS4, normData, generatePostNormalizationQcFn) {
      calls$post_step <<- list(normalizedS4 = normalizedS4, normData = normData)
      generatePostNormalizationQcFn(normalizedS4)
    },
    resolveRuvParametersFn = function(normalizedS4, input, normData, workflowData, sourceDir, getRuvGroupingVariableFn) {
      calls$resolve_ruv <<- list(
        normalizedS4 = normalizedS4,
        sourceDir = sourceDir,
        grouping = getRuvGroupingVariableFn()
      )
    },
    applyRuvCorrectionStepFn = function(normalizedS4, normData, getRuvGroupingVariableFn) {
      calls$apply_ruv <<- list(
        normalizedS4 = normalizedS4,
        grouping = getRuvGroupingVariableFn()
      )
      "ruv_obj"
    },
    finalizeRuvCleanupStepFn = function(ruvCorrectedS4, input, normData, workflowData, omicType, experimentLabel) {
      calls$cleanup <<- list(
        ruvCorrectedS4 = ruvCorrectedS4,
        omicType = omicType,
        experimentLabel = experimentLabel
      )
      "ruv_clean"
    },
    resolveStep6QcObjectFn = function(step5Object, normData) {
      calls$resolve_step6 <<- list(step5Object = step5Object, normData = normData)
      "step6_obj"
    },
    runStep6RuvQcFn = function(ruvMode, step6Object, normData, qcDir, generateRuvCorrectedQcFn) {
      calls$step6_qc <<- list(ruvMode = ruvMode, step6Object = step6Object, qcDir = qcDir)
      generateRuvCorrectedQcFn(step6Object)
    },
    generateCompositeQcFigureFn = function(ruvMode, qcDir, omicType) {
      calls$composite <<- list(ruvMode = ruvMode, qcDir = qcDir, omicType = omicType)
    },
    finalizeWorkflowStateFn = function(normData) {
      calls$finalize <<- normData
    },
    buildCompletionNotificationFn = function(ruvMode) {
      calls$completion_message <<- ruvMode
      "done"
    },
    withProgressFn = function(message, value, expr) force(expr),
    incProgressFn = function(value, detail) {
      progress <<- c(progress, detail)
    },
    showNotificationFn = function(message, type, duration) {
      notification <<- list(message = message, type = type, duration = duration)
    },
    messageFn = function(...) NULL
  )

  expect_true(isTRUE(result))
  expect_equal(calls$memory$context, "Normalization Start")
  expect_identical(calls$prepare$stateManager, "state_manager")
  expect_equal(calls$between$currentS4, "current_obj")
  expect_equal(calls$between$normMethod, "cyclicloess")
  expect_equal(calls$post_qc_generated, "normalized_obj")
  expect_equal(calls$resolve_ruv$grouping, "group")
  expect_equal(calls$apply_ruv$grouping, "group")
  expect_equal(calls$cleanup$ruvCorrectedS4, "ruv_obj")
  expect_equal(calls$resolve_step6$step5Object, "ruv_clean")
  expect_equal(calls$step6_qc$ruvMode, "automatic")
  expect_equal(calls$ruv_qc_generated, "step6_obj")
  expect_equal(calls$composite$omicType, "proteomics")
  expect_identical(calls$finalize, norm_data)
  expect_equal(calls$completion_message, "automatic")
  expect_equal(
    progress,
    c(
      "Normalizing between samples...",
      "Generating post-normalization QC plots...",
      "Determining RUV parameters...",
      "Applying RUV-III batch correction...",
      "Generating RUV-corrected QC plots..."
    )
  )
  expect_equal(notification$message, "done")
  expect_equal(notification$type, "message")
  expect_equal(notification$duration, 10)
})

test_that("resolveProtNormCorrelationResultFilenames selects the expected save names", {
  skipped_files <- resolveProtNormCorrelationResultFilenames(
    normRuvOptimizationResult = list(ruv_skipped = TRUE),
    workflowRuvOptimizationResult = NULL,
    messageFn = function(...) NULL
  )
  applied_files <- resolveProtNormCorrelationResultFilenames(
    normRuvOptimizationResult = list(ruv_skipped = FALSE),
    workflowRuvOptimizationResult = NULL,
    messageFn = function(...) NULL
  )

  expect_true(skipped_files$ruvWasSkipped)
  expect_equal(skipped_files$tsvFilename, "normalised_results_cln_with_replicates.tsv")
  expect_equal(applied_files$rdsFilename, "ruv_normalised_results_cln_with_replicates.RDS")
})

test_that("saveProtNormCorrelationResults persists the final matrix with conditional filenames", {
  if (!methods::isClass("MockCorrelationResultData")) {
    methods::setClass(
      "MockCorrelationResultData",
      slots = c(
        protein_quant_table = "data.frame",
        protein_id_column = "character"
      )
    )
  }

  final_s4 <- methods::new(
    "MockCorrelationResultData",
    protein_quant_table = data.frame(Protein.Ids = "p1", S1 = 1),
    protein_id_column = "Protein.Ids"
  )
  saved_tsv <- NULL
  saved_rds <- NULL

  file_names <- saveProtNormCorrelationResults(
    finalS4ForDe = final_s4,
    experimentPaths = list(protein_qc_dir = tempdir()),
    normRuvOptimizationResult = list(ruv_skipped = TRUE),
    workflowRuvOptimizationResult = NULL,
    writeTsvFn = function(data, path) {
      saved_tsv <<- list(data = data, path = path)
    },
    saveRdsFn = function(object, path) {
      saved_rds <<- list(object = object, path = path)
    },
    messageFn = function(...) NULL
  )

  expect_equal(saved_tsv$data, final_s4@protein_quant_table)
  expect_equal(saved_tsv$path, file.path(tempdir(), "normalised_results_cln_with_replicates.tsv"))
  expect_identical(saved_rds$object, final_s4)
  expect_equal(saved_rds$path, file.path(tempdir(), "normalised_results_cln_with_replicates.RDS"))
  expect_true(file_names$ruvWasSkipped)
})

test_that("updateProtNormFinalFilteringPlot and updateProtNormFinalQcPlot store final QC artifacts", {
  if (!methods::isClass("MockCorrelationArtifactData")) {
    methods::setClass(
      "MockCorrelationArtifactData",
      slots = c(
        protein_quant_table = "data.frame",
        protein_id_column = "character"
      )
    )
  }

  final_s4 <- methods::new(
    "MockCorrelationArtifactData",
    protein_quant_table = data.frame(Protein.Ids = c("p1", "p2"), S1 = c(1, 2)),
    protein_id_column = "Protein.Ids"
  )
  norm_data <- new.env(parent = emptyenv())

  filtering_plot <- updateProtNormFinalFilteringPlot(
    finalS4ForDe = final_s4,
    normData = norm_data,
    omicType = "proteomics",
    experimentLabel = "exp1",
    updateProteinFilteringFn = function(data, step_name, omic_type, experiment_label, return_grid, overwrite) {
      expect_equal(data, final_s4@protein_quant_table)
      expect_equal(step_name, "12_correlation_filtered")
      expect_equal(omic_type, "proteomics")
      expect_equal(experiment_label, "exp1")
      expect_true(return_grid)
      expect_true(overwrite)
      "filter_plot"
    },
    messageFn = function(...) NULL
  )

  qc_plot <- updateProtNormFinalQcPlot(
    finalS4ForDe = final_s4,
    normData = norm_data,
    title = "Final Data",
    getPlotAestheticsFn = function() list(color_var = "group", shape_var = "batch"),
    plotPcaFn = function(object, grouping_variable, label_column, shape_variable, title, font_size) {
      expect_identical(object, final_s4)
      expect_equal(grouping_variable, "group")
      expect_equal(shape_variable, "batch")
      expect_equal(title, "Final Data")
      expect_equal(font_size, 8)
      "qc_plot"
    },
    messageFn = function(...) NULL
  )

  expect_equal(filtering_plot, "filter_plot")
  expect_equal(norm_data$final_filtering_plot, "filter_plot")
  expect_equal(qc_plot, "qc_plot")
  expect_equal(norm_data$final_qc_plot, "qc_plot")
})

test_that("finalizeProtNormCorrelationWorkflowState stores DE-ready state and counts", {
  if (!methods::isClass("MockCorrelationWorkflowData")) {
    methods::setClass(
      "MockCorrelationWorkflowData",
      slots = c(
        protein_quant_table = "data.frame",
        protein_id_column = "character"
      )
    )
  }

  final_s4 <- methods::new(
    "MockCorrelationWorkflowData",
    protein_quant_table = data.frame(
      Protein.Ids = c("p1", "p1", "p2"),
      S1 = c(1, 2, 3),
      S2 = c(4, 5, 6)
    ),
    protein_id_column = "Protein.Ids"
  )
  workflow_data <- new.env(parent = emptyenv())
  workflow_data$state_manager <- new.env(parent = emptyenv())
  saved_state <- NULL
  workflow_data$state_manager$saveState <- function(...) {
    saved_state <<- list(...)
  }
  workflow_data$tab_status <- list(normalization = "pending", differential_expression = "disabled")
  workflow_data$state_update_trigger <- as.POSIXct("2026-04-11 10:00:00", tz = "UTC")

  metrics <- finalizeProtNormCorrelationWorkflowState(
    finalS4ForDe = final_s4,
    workflowData = workflow_data,
    correlationThreshold = 0.65,
    skipped = FALSE,
    timeFn = function() as.POSIXct("2026-04-11 10:05:00", tz = "UTC"),
    messageFn = function(...) NULL,
    catFn = function(...) NULL
  )

  expect_identical(workflow_data$ruv_normalised_for_da_analysis_obj, final_s4)
  expect_equal(saved_state$state_name, "correlation_filtered")
  expect_equal(saved_state$config_object$min_pearson_correlation_threshold, 0.65)
  expect_equal(saved_state$description, "Applied final sample correlation filter (chunk 28)")
  expect_equal(workflow_data$tab_status$normalization, "complete")
  expect_equal(workflow_data$tab_status$differential_expression, "pending")
  expect_equal(workflow_data$protein_counts$final_for_de, 2)
  expect_equal(metrics$finalProteinCount, 2)
  expect_equal(metrics$finalSampleCount, 2)
})

test_that("buildProtNormCorrelationSummaryText returns apply and skip summaries", {
  completed_summary <- buildProtNormCorrelationSummaryText(
    finalProteinCount = 12,
    finalSampleCount = 4,
    correlationThreshold = 0.75,
    skipped = FALSE
  )
  skipped_summary <- buildProtNormCorrelationSummaryText(
    finalProteinCount = 12,
    finalSampleCount = 4,
    skipped = TRUE
  )

  expect_match(completed_summary, "Threshold: 0.75", fixed = TRUE)
  expect_match(completed_summary, "Proteins remaining: 12", fixed = TRUE)
  expect_match(skipped_summary, "Correlation filtering SKIPPED.", fixed = TRUE)
  expect_match(skipped_summary, "All samples retained.", fixed = TRUE)
})

test_that("resolveProtNormCorrelationInputObject validates presence and logs dimensions", {
  if (!methods::isClass("MockCorrelationInputData")) {
    methods::setClass(
      "MockCorrelationInputData",
      slots = c(
        protein_quant_table = "data.frame",
        protein_id_column = "character"
      )
    )
  }

  mock_s4 <- methods::new(
    "MockCorrelationInputData",
    protein_quant_table = data.frame(Protein.Ids = c("p1", "p2"), S1 = c(1, 2), S2 = c(3, 4)),
    protein_id_column = "Protein.Ids"
  )
  messages <- character()

  result <- resolveProtNormCorrelationInputObject(
    ruvNormalizedObj = mock_s4,
    startMessage = "Starting Correlation Filter Flow",
    missingObjectMessage = "missing",
    messageFn = function(text) {
      messages <<- c(messages, text)
    }
  )

  expect_identical(result, mock_s4)
  expect_true(any(grepl("Starting Correlation Filter Flow", messages, fixed = TRUE)))
  expect_true(any(grepl("2 proteins x 2 samples", messages, fixed = TRUE)))
  expect_error(
    resolveProtNormCorrelationInputObject(
      ruvNormalizedObj = NULL,
      startMessage = "start",
      missingObjectMessage = "missing",
      messageFn = function(...) NULL
    ),
    "missing"
  )
})

test_that("runProtNormCorrelationVectorStep stores the correlation vector and threshold", {
  if (!methods::isClass("MockCorrelationVectorData")) {
    methods::setClass(
      "MockCorrelationVectorData",
      slots = c(
        protein_quant_table = "data.frame",
        protein_id_column = "character"
      )
    )
  }

  mock_s4 <- methods::new(
    "MockCorrelationVectorData",
    protein_quant_table = data.frame(Protein.Ids = c("p1", "p2"), S1 = c(1, 2), S2 = c(3, 4)),
    protein_id_column = "Protein.Ids"
  )
  norm_data <- new.env(parent = emptyenv())
  correlation_vec <- data.frame(sample = "S1", correlation = 0.8)
  time_values <- as.POSIXct(c("2026-04-11 10:00:00", "2026-04-11 10:00:03"), tz = "UTC")
  idx <- 0

  result <- runProtNormCorrelationVectorStep(
    ruvS4 = mock_s4,
    correlationThreshold = 0.7,
    normData = norm_data,
    getRuvGroupingVariableFn = function() "group",
    pearsonCorForSamplePairsFn = function(object, tech_rep_remove_regex, correlation_group) {
      expect_identical(object, mock_s4)
      expect_equal(tech_rep_remove_regex, "pool")
      expect_equal(correlation_group, "group")
      correlation_vec
    },
    timeFn = function() {
      idx <<- idx + 1
      time_values[[idx]]
    },
    messageFn = function(...) NULL
  )

  expect_equal(result, correlation_vec)
  expect_equal(norm_data$correlation_vector, correlation_vec)
  expect_equal(norm_data$correlation_threshold, 0.7)
})

test_that("runProtNormCorrelationFilterStep stores the filtered object and cleans up", {
  if (!methods::isClass("MockCorrelationFilterData")) {
    methods::setClass(
      "MockCorrelationFilterData",
      slots = c(
        protein_quant_table = "data.frame",
        protein_id_column = "character"
      )
    )
  }

  mock_s4 <- methods::new(
    "MockCorrelationFilterData",
    protein_quant_table = data.frame(Protein.Ids = c("p1", "p2"), S1 = c(1, 2), S2 = c(3, 4)),
    protein_id_column = "Protein.Ids"
  )
  filtered_s4 <- methods::new(
    "MockCorrelationFilterData",
    protein_quant_table = data.frame(Protein.Ids = "p1", S1 = 1, S2 = 3),
    protein_id_column = "Protein.Ids"
  )
  norm_data <- new.env(parent = emptyenv())
  gc_called <- FALSE
  time_values <- as.POSIXct(c("2026-04-11 10:00:00", "2026-04-11 10:00:02"), tz = "UTC")
  idx <- 0

  result <- runProtNormCorrelationFilterStep(
    ruvS4 = mock_s4,
    correlationVec = data.frame(sample = "S1", correlation = 0.8),
    correlationThreshold = 0.75,
    normData = norm_data,
    filterSamplesFn = function(ruvS4, pearson_correlation_per_pair, min_pearson_correlation_threshold) {
      expect_identical(ruvS4, mock_s4)
      expect_equal(min_pearson_correlation_threshold, 0.75)
      expect_equal(pearson_correlation_per_pair$correlation, 0.8)
      filtered_s4
    },
    gcFn = function() {
      gc_called <<- TRUE
    },
    timeFn = function() {
      idx <<- idx + 1
      time_values[[idx]]
    },
    messageFn = function(...) NULL
  )

  expect_identical(result, filtered_s4)
  expect_identical(norm_data$correlation_filtered_obj, filtered_s4)
  expect_true(gc_called)
})

test_that("prepareProtNormSkippedCorrelationState clears correlation metadata", {
  if (!methods::isClass("MockSkippedCorrelationData")) {
    methods::setClass(
      "MockSkippedCorrelationData",
      slots = c(
        protein_quant_table = "data.frame",
        protein_id_column = "character"
      )
    )
  }

  mock_s4 <- methods::new(
    "MockSkippedCorrelationData",
    protein_quant_table = data.frame(Protein.Ids = "p1", S1 = 1),
    protein_id_column = "Protein.Ids"
  )
  norm_data <- new.env(parent = emptyenv())
  norm_data$correlation_vector <- "old"
  norm_data$correlation_threshold <- 0.7
  norm_data$correlation_filtered_obj <- "old_obj"
  gc_called <- FALSE

  result <- prepareProtNormSkippedCorrelationState(
    ruvS4 = mock_s4,
    normData = norm_data,
    gcFn = function() {
      gc_called <<- TRUE
    },
    messageFn = function(...) NULL
  )

  expect_identical(result, mock_s4)
  expect_null(norm_data$correlation_vector)
  expect_null(norm_data$correlation_threshold)
  expect_identical(norm_data$correlation_filtered_obj, mock_s4)
  expect_true(gc_called)
})

test_that("runProtNormApplyCorrelationWorkflow sequences the apply correlation shell", {
  if (!methods::isClass("MockApplyCorrelationWorkflowData")) {
    methods::setClass(
      "MockApplyCorrelationWorkflowData",
      slots = c(
        protein_quant_table = "data.frame",
        protein_id_column = "character"
      )
    )
  }

  ruv_s4 <- methods::new(
    "MockApplyCorrelationWorkflowData",
    protein_quant_table = data.frame(Protein.Ids = c("p1", "p2"), S1 = c(1, 2)),
    protein_id_column = "Protein.Ids"
  )
  final_s4 <- methods::new(
    "MockApplyCorrelationWorkflowData",
    protein_quant_table = data.frame(Protein.Ids = "p1", S1 = 1),
    protein_id_column = "Protein.Ids"
  )
  norm_data <- new.env(parent = emptyenv())
  norm_data$ruv_optimization_result <- list(ruv_skipped = FALSE)
  workflow_data <- new.env(parent = emptyenv())
  workflow_data$ruv_optimization_result <- list(ruv_skipped = FALSE)
  progress_calls <- list()
  calls <- list()

  result <- runProtNormApplyCorrelationWorkflow(
    ruvS4 = ruv_s4,
    correlationThreshold = 0.7,
    normData = norm_data,
    experimentPaths = list(protein_qc_dir = tempdir()),
    workflowData = workflow_data,
    omicType = "proteomics",
    experimentLabel = "exp1",
    getRuvGroupingVariableFn = function() "group",
    getPlotAestheticsFn = function() list(color_var = "group", shape_var = "batch"),
    runVectorStepFn = function(...) {
      calls$vector <<- list(...)
      data.frame(sample = "S1", correlation = 0.9)
    },
    runFilterStepFn = function(...) {
      calls$filter <<- list(...)
      final_s4
    },
    updateFinalFilteringPlotFn = function(...) {
      calls$updateFiltering <<- list(...)
      "filter_plot"
    },
    updateFinalQcPlotFn = function(...) {
      calls$updateQc <<- list(...)
      "qc_plot"
    },
    saveCorrelationResultsFn = function(...) {
      calls$save <<- list(...)
      invisible(NULL)
    },
    withProgressFn = function(message, value, expr) force(expr),
    incProgressFn = function(value, detail) {
      progress_calls[[length(progress_calls) + 1]] <<- list(value = value, detail = detail)
    },
    gcFn = function() {
      calls$gc <<- TRUE
    },
    messageFn = function(...) NULL
  )

  expect_identical(result, final_s4)
  expect_equal(progress_calls[[1]]$detail, "Calculating sample correlations...")
  expect_equal(progress_calls[[2]]$detail, "Filtering low-correlation samples...")
  expect_equal(progress_calls[[3]]$detail, "Updating tracking...")
  expect_equal(progress_calls[[4]]$detail, "Saving results...")
  expect_equal(calls$vector$correlationThreshold, 0.7)
  expect_identical(calls$filter$correlationVec$correlation, 0.9)
  expect_identical(calls$updateFiltering$finalS4ForDe, final_s4)
  expect_identical(calls$updateQc$finalS4ForDe, final_s4)
  expect_identical(calls$save$finalS4ForDe, final_s4)
  expect_true(isTRUE(calls$gc))
})

test_that("runProtNormSkipCorrelationWorkflow sequences the skip correlation shell", {
  if (!methods::isClass("MockSkipCorrelationWorkflowData")) {
    methods::setClass(
      "MockSkipCorrelationWorkflowData",
      slots = c(
        protein_quant_table = "data.frame",
        protein_id_column = "character"
      )
    )
  }

  ruv_s4 <- methods::new(
    "MockSkipCorrelationWorkflowData",
    protein_quant_table = data.frame(Protein.Ids = c("p1", "p2"), S1 = c(1, 2)),
    protein_id_column = "Protein.Ids"
  )
  norm_data <- new.env(parent = emptyenv())
  norm_data$ruv_optimization_result <- list(ruv_skipped = TRUE)
  workflow_data <- new.env(parent = emptyenv())
  workflow_data$ruv_optimization_result <- list(ruv_skipped = TRUE)
  progress_calls <- list()
  calls <- list()

  result <- runProtNormSkipCorrelationWorkflow(
    ruvS4 = ruv_s4,
    normData = norm_data,
    experimentPaths = list(protein_qc_dir = tempdir()),
    workflowData = workflow_data,
    omicType = "proteomics",
    experimentLabel = "exp1",
    getPlotAestheticsFn = function() list(color_var = "group", shape_var = "batch"),
    prepareSkippedStateFn = function(...) {
      calls$prepare <<- list(...)
      ruv_s4
    },
    updateFinalFilteringPlotFn = function(...) {
      calls$updateFiltering <<- list(...)
      "filter_plot"
    },
    updateFinalQcPlotFn = function(...) {
      calls$updateQc <<- list(...)
      "qc_plot"
    },
    saveCorrelationResultsFn = function(...) {
      calls$save <<- list(...)
      invisible(NULL)
    },
    withProgressFn = function(message, value, expr) force(expr),
    incProgressFn = function(value, detail) {
      progress_calls[[length(progress_calls) + 1]] <<- list(value = value, detail = detail)
    },
    gcFn = function() {
      calls$gc <<- TRUE
    },
    messageFn = function(...) NULL
  )

  expect_identical(result, ruv_s4)
  expect_equal(progress_calls[[1]]$detail, "Bypassing filter...")
  expect_equal(progress_calls[[2]]$detail, "Updating tracking...")
  expect_equal(progress_calls[[3]]$detail, "Saving results...")
  expect_identical(calls$prepare$ruvS4, ruv_s4)
  expect_identical(calls$updateFiltering$finalS4ForDe, ruv_s4)
  expect_identical(calls$updateQc$finalS4ForDe, ruv_s4)
  expect_identical(calls$save$finalS4ForDe, ruv_s4)
  expect_true(isTRUE(calls$gc))
})

test_that("completeProtNormCorrelationWorkflow finalizes summary state and notification", {
  if (!methods::isClass("MockCompleteCorrelationWorkflowData")) {
    methods::setClass(
      "MockCompleteCorrelationWorkflowData",
      slots = c(
        protein_quant_table = "data.frame",
        protein_id_column = "character"
      )
    )
  }

  final_s4 <- methods::new(
    "MockCompleteCorrelationWorkflowData",
    protein_quant_table = data.frame(Protein.Ids = c("p1", "p2"), S1 = c(1, 2), S2 = c(3, 4)),
    protein_id_column = "Protein.Ids"
  )
  workflow_data <- new.env(parent = emptyenv())
  workflow_data$state_manager <- new.env(parent = emptyenv())
  workflow_data$state_manager$saveState <- function(...) NULL
  workflow_data$tab_status <- list(normalization = "pending", differential_expression = "disabled")
  workflow_data$state_update_trigger <- as.POSIXct("2026-04-11 10:00:00", tz = "UTC")
  output <- new.env(parent = emptyenv())
  norm_data <- new.env(parent = emptyenv())
  notification <- NULL

  completeProtNormCorrelationWorkflow(
    finalS4ForDe = final_s4,
    workflowData = workflow_data,
    output = output,
    normData = norm_data,
    correlationThreshold = 0.8,
    skipped = FALSE,
    successNotification = "done",
    completionMessage = "completed",
    messagePrefix = "*** CORRELATION",
    showNotificationFn = function(message, type, duration) {
      notification <<- list(message = message, type = type, duration = duration)
    },
    renderTextFn = function(text) text,
    messageFn = function(...) NULL
  )

  expect_true(isTRUE(norm_data$correlation_filtering_complete))
  expect_match(output$correlation_filter_summary, "Threshold: 0.80", fixed = TRUE)
  expect_equal(notification$message, "done")
  expect_equal(notification$type, "message")
  expect_equal(workflow_data$tab_status$differential_expression, "pending")
})

test_that("handleProtNormCorrelationError logs and notifies consistently", {
  messages <- character()
  notification <- NULL

  handleProtNormCorrelationError(
    error = structure(list(message = "boom"), class = c("simpleError", "error", "condition")),
    logPrefix = "Error in correlation filtering:",
    notificationPrefix = "Error in correlation filtering:",
    showNotificationFn = function(message, type, duration) {
      notification <<- list(message = message, type = type, duration = duration)
    },
    messageFn = function(text) {
      messages <<- c(messages, text)
    }
  )

  expect_true(any(grepl("boom", messages, fixed = TRUE)))
  expect_equal(notification$message, "Error in correlation filtering: boom")
  expect_equal(notification$type, "error")
  expect_equal(notification$duration, 10)
})

test_that("runProtNormApplyCorrelationObserver sequences the apply observer shell", {
  norm_data <- new.env(parent = emptyenv())
  norm_data$ruv_normalized_obj <- "ruv_obj"
  calls <- list()

  result <- runProtNormApplyCorrelationObserver(
    input = list(min_pearson_correlation_threshold = 0.8),
    output = "output",
    workflowData = "workflow",
    normData = norm_data,
    experimentPaths = "paths",
    omicType = "proteomics",
    experimentLabel = "exp1",
    getRuvGroupingVariableFn = function() "group",
    getPlotAestheticsFn = function() list(color_var = "group"),
    checkMemoryUsageFn = function(threshold_gb, context) {
      calls$memory <<- list(threshold_gb = threshold_gb, context = context)
    },
    resolveCorrelationInputObjectFn = function(ruvNormalizedObj, startMessage, missingObjectMessage) {
      calls$resolve <<- list(
        ruvNormalizedObj = ruvNormalizedObj,
        startMessage = startMessage,
        missingObjectMessage = missingObjectMessage
      )
      "resolved_ruv"
    },
    runApplyCorrelationWorkflowFn = function(ruvS4, correlationThreshold, normData, experimentPaths, workflowData, omicType, experimentLabel, getRuvGroupingVariableFn, getPlotAestheticsFn) {
      calls$apply <<- list(
        ruvS4 = ruvS4,
        correlationThreshold = correlationThreshold,
        normData = normData,
        experimentPaths = experimentPaths,
        workflowData = workflowData,
        omicType = omicType,
        experimentLabel = experimentLabel,
        grouping = getRuvGroupingVariableFn(),
        aesthetics = getPlotAestheticsFn()
      )
      "final_obj"
    },
    completeCorrelationWorkflowFn = function(finalS4ForDe, workflowData, output, normData, correlationThreshold, skipped, successNotification, completionMessage, messagePrefix) {
      calls$complete <<- list(
        finalS4ForDe = finalS4ForDe,
        workflowData = workflowData,
        output = output,
        normData = normData,
        correlationThreshold = correlationThreshold,
        skipped = skipped,
        successNotification = successNotification,
        completionMessage = completionMessage,
        messagePrefix = messagePrefix
      )
      invisible(NULL)
    },
    messageFn = function(...) NULL
  )

  expect_equal(calls$memory$context, "Correlation Filtering Start")
  expect_equal(calls$resolve$ruvNormalizedObj, "ruv_obj")
  expect_equal(calls$apply$ruvS4, "resolved_ruv")
  expect_equal(calls$apply$correlationThreshold, 0.8)
  expect_equal(calls$apply$grouping, "group")
  expect_equal(calls$apply$experimentLabel, "exp1")
  expect_equal(calls$complete$finalS4ForDe, "final_obj")
  expect_false(calls$complete$skipped)
  expect_equal(calls$complete$correlationThreshold, 0.8)
  expect_equal(calls$complete$messagePrefix, "*** CORRELATION")
  expect_equal(result, "final_obj")
})

test_that("runProtNormSkipCorrelationObserver sequences the skip observer shell", {
  norm_data <- new.env(parent = emptyenv())
  norm_data$ruv_normalized_obj <- "ruv_obj"
  calls <- list()

  result <- runProtNormSkipCorrelationObserver(
    output = "output",
    workflowData = "workflow",
    normData = norm_data,
    experimentPaths = "paths",
    omicType = "proteomics",
    experimentLabel = "exp1",
    getPlotAestheticsFn = function() list(color_var = "group"),
    checkMemoryUsageFn = function(threshold_gb, context) {
      calls$memory <<- list(threshold_gb = threshold_gb, context = context)
    },
    resolveCorrelationInputObjectFn = function(ruvNormalizedObj, startMessage, missingObjectMessage) {
      calls$resolve <<- list(
        ruvNormalizedObj = ruvNormalizedObj,
        startMessage = startMessage,
        missingObjectMessage = missingObjectMessage
      )
      "resolved_ruv"
    },
    runSkipCorrelationWorkflowFn = function(ruvS4, normData, experimentPaths, workflowData, omicType, experimentLabel, getPlotAestheticsFn) {
      calls$skip <<- list(
        ruvS4 = ruvS4,
        normData = normData,
        experimentPaths = experimentPaths,
        workflowData = workflowData,
        omicType = omicType,
        experimentLabel = experimentLabel,
        aesthetics = getPlotAestheticsFn()
      )
      "final_obj"
    },
    completeCorrelationWorkflowFn = function(finalS4ForDe, workflowData, output, normData, correlationThreshold, skipped, successNotification, completionMessage, messagePrefix) {
      calls$complete <<- list(
        finalS4ForDe = finalS4ForDe,
        workflowData = workflowData,
        output = output,
        normData = normData,
        correlationThreshold = correlationThreshold,
        skipped = skipped,
        successNotification = successNotification,
        completionMessage = completionMessage,
        messagePrefix = messagePrefix
      )
      invisible(NULL)
    },
    messageFn = function(...) NULL
  )

  expect_equal(calls$memory$context, "Skip Correlation Filtering Start")
  expect_equal(calls$resolve$ruvNormalizedObj, "ruv_obj")
  expect_equal(calls$skip$ruvS4, "resolved_ruv")
  expect_equal(calls$skip$experimentLabel, "exp1")
  expect_equal(calls$complete$finalS4ForDe, "final_obj")
  expect_true(calls$complete$skipped)
  expect_equal(calls$complete$correlationThreshold, 0)
  expect_equal(calls$complete$messagePrefix, "*** SKIP CORRELATION")
  expect_equal(result, "final_obj")
})

test_that("canProtNormExportFilteredSession and resolveProtNormExportSourceDir enforce export readiness", {
  expect_true(canProtNormExportFilteredSession(TRUE, "obj"))
  expect_false(canProtNormExportFilteredSession(FALSE, "obj"))
  expect_false(canProtNormExportFilteredSession(TRUE, NULL))

  expect_equal(
    resolveProtNormExportSourceDir(
      experimentPaths = list(source_dir = tempdir()),
      dirExistsFn = function(path) TRUE
    ),
    tempdir()
  )
  expect_error(
    resolveProtNormExportSourceDir(
      experimentPaths = list(source_dir = tempdir()),
      dirExistsFn = function(path) FALSE
    ),
    "Could not find the source directory"
  )
})

test_that("collectProtNormExportSessionData builds the export payload with counts and workflow type", {
  if (!methods::isClass("MockProtNormExportData")) {
    methods::setClass(
      "MockProtNormExportData",
      slots = c(
        protein_quant_table = "data.frame",
        protein_id_column = "character",
        args = "list"
      )
    )
  }

  current_s4 <- methods::new(
    "MockProtNormExportData",
    protein_quant_table = data.frame(Protein.Ids = c("p1", "p1", "p2"), S1 = c(1, 2, 3), S2 = c(4, 5, 6)),
    protein_id_column = "Protein.Ids",
    args = list(globalParameters = list(workflow_type = "DIA"))
  )
  workflow_data <- new.env(parent = emptyenv())
  workflow_data$state_manager <- new.env(parent = emptyenv())
  workflow_data$state_manager$current_state <- "correlation_filtered"
  workflow_data$state_manager$getState <- function(state) current_s4
  workflow_data$contrasts_tbl <- data.frame(friendly_names = c("A vs B"), stringsAsFactors = FALSE)
  workflow_data$design_matrix <- data.frame(group = c("A", "B"), stringsAsFactors = FALSE)
  workflow_data$config_list <- list(globalParameters = list(workflow_type = "TMT"))
  workflow_data$fasta_metadata <- list(db = "ref")
  workflow_data$accession_cleanup_results <- list(clean = TRUE)
  workflow_data$ruv_optimization_result <- list(best_k = 2)
  workflow_data$qc_params <- list(min = 1)
  workflow_data$protein_counts <- list(final = 2)
  workflow_data$mixed_species_analysis <- list(enabled = TRUE)
  norm_data <- new.env(parent = emptyenv())
  norm_data$correlation_filtered_obj <- "filtered_obj"
  norm_data$best_k <- 3
  norm_data$correlation_threshold <- 0.7

  session_data <- collectProtNormExportSessionData(
    workflowData = workflow_data,
    normData = norm_data,
    input = list(norm_method = "cyclicloess", ruv_mode = "automatic"),
    timeFn = function() as.POSIXct("2026-04-11 12:00:00", tz = "UTC"),
    messageFn = function(...) NULL
  )

  expect_equal(session_data$r6_current_state_name, "correlation_filtered")
  expect_identical(session_data$current_s4_object, current_s4)
  expect_equal(session_data$workflow_type, "TMT")
  expect_equal(session_data$normalization_method, "cyclicloess")
  expect_true(session_data$ruv_applied)
  expect_equal(session_data$ruv_k, 3)
  expect_equal(session_data$correlation_threshold, 0.7)
  expect_equal(session_data$final_protein_count, 2)
  expect_equal(session_data$final_sample_count, 2)
})

test_that("buildProtNormExportSummaryContent and saveProtNormExportArtifacts create expected export outputs", {
  session_data <- list(
    final_protein_count = 12,
    final_sample_count = 4,
    contrasts_tbl = data.frame(friendly_names = c("A vs B", "A vs C"), stringsAsFactors = FALSE),
    normalization_method = "cyclicloess",
    ruv_mode = "automatic",
    ruv_k = 2,
    correlation_threshold = 0.75,
    accession_cleanup_results = list(clean = TRUE),
    qc_params = list(min = 1),
    fasta_metadata = list(db = "ref"),
    ruv_optimization_result = list(best_k = 2),
    mixed_species_analysis = list(enabled = TRUE)
  )
  summary_text <- buildProtNormExportSummaryContent(
    sessionData = session_data,
    sessionFilename = "filtered_session_data_20260411_120000.rds",
    timeFn = function() as.POSIXct("2026-04-11 12:00:00", tz = "UTC"),
    formatTimeFn = function(x, format) {
      if (identical(format, "%Y-%m-%d %H:%M:%S")) {
        "2026-04-11 12:00:00"
      } else {
        "20260411_120000"
      }
    }
  )

  expect_match(summary_text, "Proteins: 12", fixed = TRUE)
  expect_match(summary_text, "A vs B\nA vs C", fixed = TRUE)

  saved_rds <- list()
  written_summary <- NULL
  artifacts <- saveProtNormExportArtifacts(
    sessionData = session_data,
    sourceDir = tempdir(),
    timeFn = function() as.POSIXct("2026-04-11 12:00:00", tz = "UTC"),
    formatTimeFn = function(x, format) {
      if (identical(format, "%Y%m%d_%H%M%S")) {
        "20260411_120000"
      } else {
        "2026-04-11 12:00:00"
      }
    },
    saveRdsFn = function(object, path) {
      saved_rds[[length(saved_rds) + 1]] <<- path
    },
    writeLinesFn = function(text, con) {
      written_summary <<- list(text = text, con = con)
    },
    fileExistsFn = function(path) FALSE,
    messageFn = function(...) NULL
  )

  expect_equal(artifacts$sessionFilename, "filtered_session_data_20260411_120000.rds")
  expect_equal(artifacts$sessionFilepath, file.path(tempdir(), "filtered_session_data_20260411_120000.rds"))
  expect_equal(artifacts$latestFilepath, file.path(tempdir(), "filtered_session_data_latest.rds"))
  expect_equal(artifacts$summaryFilepath, file.path(tempdir(), "filtered_session_summary.txt"))
  expect_true(any(grepl("filtered_session_data_20260411_120000.rds", saved_rds, fixed = TRUE)))
  expect_true(any(grepl("accession_cleanup_results.RDS", saved_rds, fixed = TRUE)))
  expect_equal(written_summary$con, file.path(tempdir(), "filtered_session_summary.txt"))
})

test_that("runProtNormExportSessionWorkflow sequences collection and artifact saving", {
  progress_calls <- list()
  result <- runProtNormExportSessionWorkflow(
    workflowData = "workflow",
    normData = "norm",
    input = list(norm_method = "cyclicloess", ruv_mode = "automatic"),
    sourceDir = tempdir(),
    withProgressFn = function(message, value, expr) force(expr),
    incProgressFn = function(value, detail) {
      progress_calls[[length(progress_calls) + 1]] <<- list(value = value, detail = detail)
    },
    collectSessionDataFn = function(workflowData, normData, input, messageFn) {
      expect_equal(workflowData, "workflow")
      expect_equal(normData, "norm")
      list(id = 1)
    },
    saveExportArtifactsFn = function(sessionData, sourceDir, messageFn) {
      expect_equal(sessionData$id, 1)
      expect_equal(sourceDir, tempdir())
      list(sessionFilename = "filtered_session_data_test.rds")
    },
    messageFn = function(...) NULL
  )

  expect_equal(progress_calls[[1]]$detail, "Gathering data...")
  expect_equal(progress_calls[[2]]$detail, "Saving to file...")
  expect_equal(progress_calls[[3]]$detail, "Creating latest version...")
  expect_equal(progress_calls[[4]]$detail, "Saving metadata files...")
  expect_equal(progress_calls[[5]]$detail, "Creating summary...")
  expect_equal(result$sessionData$id, 1)
  expect_equal(result$exportArtifacts$sessionFilename, "filtered_session_data_test.rds")
})

test_that("notifyProtNormExportSessionPrereqWarning logs and notifies consistently", {
  messages <- character()
  notification <- NULL

  result <- notifyProtNormExportSessionPrereqWarning(
    showNotificationFn = function(message, type, duration) {
      notification <<- list(message = message, type = type, duration = duration)
    },
    messageFn = function(text) {
      messages <<- c(messages, text)
    }
  )

  expect_false(isTRUE(result))
  expect_true(any(grepl("correlation filtering is not complete", messages, fixed = TRUE)))
  expect_equal(notification$message, "Please complete correlation filtering before exporting session data.")
  expect_equal(notification$type, "warning")
  expect_equal(notification$duration, 5)
})

test_that("runProtNormExportObserver sequences prereq check and export success path", {
  norm_data <- new.env(parent = emptyenv())
  norm_data$correlation_filtering_complete <- TRUE
  norm_data$correlation_filtered_obj <- "filtered_obj"
  notification <- NULL
  calls <- list()
  messages <- character()

  result <- runProtNormExportObserver(
    input = list(norm_method = "cyclicloess"),
    workflowData = "workflow",
    normData = norm_data,
    experimentPaths = "paths",
    canExportFilteredSessionFn = function(correlationFilteringComplete, correlationFilteredObj) {
      calls$can_export <<- list(
        correlationFilteringComplete = correlationFilteringComplete,
        correlationFilteredObj = correlationFilteredObj
      )
      TRUE
    },
    notifyExportPrereqFn = function() {
      calls$notify <<- TRUE
      invisible(FALSE)
    },
    resolveExportSourceDirFn = function(experimentPaths) {
      calls$resolve_dir <<- experimentPaths
      "/tmp/source"
    },
    runExportSessionWorkflowFn = function(workflowData, normData, input, sourceDir) {
      calls$export <<- list(
        workflowData = workflowData,
        normData = normData,
        input = input,
        sourceDir = sourceDir
      )
      list(exportArtifacts = list(sessionFilename = "filtered_session_data_test.rds"))
    },
    showNotificationFn = function(message, type, duration) {
      notification <<- list(message = message, type = type, duration = duration)
    },
    messageFn = function(text) {
      messages <<- c(messages, text)
    }
  )

  expect_true(calls$can_export$correlationFilteringComplete)
  expect_equal(calls$can_export$correlationFilteredObj, "filtered_obj")
  expect_equal(calls$resolve_dir, "paths")
  expect_equal(calls$export$sourceDir, "/tmp/source")
  expect_null(calls$notify)
  expect_match(notification$message, "filtered_session_data_test.rds", fixed = TRUE)
  expect_equal(notification$type, "message")
  expect_equal(notification$duration, 10)
  expect_true(any(grepl("COMPLETED SUCCESSFULLY", messages, fixed = TRUE)))
  expect_equal(result$exportArtifacts$sessionFilename, "filtered_session_data_test.rds")
})

test_that("handleProtNormExportError logs and notifies consistently", {
  messages <- character()
  notification <- NULL

  handleProtNormExportError(
    error = structure(list(message = "boom"), class = c("simpleError", "error", "condition")),
    showNotificationFn = function(message, type, duration) {
      notification <<- list(message = message, type = type, duration = duration)
    },
    messageFn = function(text) {
      messages <<- c(messages, text)
    }
  )

  expect_true(any(grepl("boom", messages, fixed = TRUE)))
  expect_equal(notification$message, "Error exporting session: boom")
  expect_equal(notification$type, "error")
  expect_equal(notification$duration, 10)
})

test_that("resolveProtNormPreNormalizationState selects the most recent eligible state", {
  history <- c(
    "raw_data_s4",
    "qvalue_filtered",
    "precursor_rollup",
    "sample_filtered",
    "protein_replicate_filtered",
    "normalised"
  )

  expect_equal(
    resolveProtNormPreNormalizationState(history),
    "protein_replicate_filtered"
  )
  expect_null(resolveProtNormPreNormalizationState(c("normalised", "ruv_corrected")))
})

test_that("revertProtNormStateManagerToPreNormalization reverts when a valid state exists", {
  workflow_data <- new.env(parent = emptyenv())
  workflow_data$state_manager <- new.env(parent = emptyenv())
  workflow_data$state_manager$getHistory <- function() {
    c("raw_data_s4", "sample_filtered", "protein_replicate_filtered", "normalised")
  }

  reverted_to <- NULL
  messages <- character()
  workflow_data$state_manager$revertToState <- function(state) {
    reverted_to <<- state
    "reverted_object"
  }

  result <- revertProtNormStateManagerToPreNormalization(
    workflowData = workflow_data,
    messageFn = function(text) {
      messages <<- c(messages, text)
    }
  )

  expect_equal(result$previousState, "protein_replicate_filtered")
  expect_equal(result$revertedS4, "reverted_object")
  expect_equal(reverted_to, "protein_replicate_filtered")
  expect_true(any(grepl("Reverted R6 state manager", messages, fixed = TRUE)))
})

test_that("revertProtNormStateManagerToPreNormalization handles missing state manager gracefully", {
  workflow_data <- new.env(parent = emptyenv())
  workflow_data$state_manager <- NULL
  messages <- character()

  result <- revertProtNormStateManagerToPreNormalization(
    workflowData = workflow_data,
    messageFn = function(text) {
      messages <<- c(messages, text)
    }
  )

  expect_null(result$previousState)
  expect_null(result$revertedS4)
  expect_true(any(grepl("state_manager is NULL", messages, fixed = TRUE)))
})

test_that("resetProtNormReactiveState clears normalization state and tab status", {
  norm_data <- new.env(parent = emptyenv())
  norm_data$normalization_complete <- TRUE
  norm_data$ruv_complete <- TRUE
  norm_data$correlation_filtering_complete <- TRUE
  norm_data$normalized_protein_obj <- "normalized"
  norm_data$ruv_normalized_obj <- "ruv"
  norm_data$correlation_filtered_obj <- "filtered"
  norm_data$best_k <- 2
  norm_data$control_genes_index <- c(TRUE, FALSE)
  norm_data$correlation_vector <- data.frame(x = 1)
  norm_data$correlation_threshold <- 0.8
  norm_data$final_qc_plot <- "pca"
  norm_data$final_filtering_plot <- "filtering"
  norm_data$post_norm_filtering_plot <- "post"
  norm_data$filtering_summary_text <- "summary"
  norm_data$ruv_optimization_result <- list(best_k = 2)
  norm_data$qc_plots <- list(
    post_normalization = list(pca = "plot"),
    ruv_corrected = list(pca = "plot")
  )

  workflow_data <- new.env(parent = emptyenv())
  workflow_data$ruv_normalised_for_da_analysis_obj <- "da_obj"
  workflow_data$tab_status <- list(normalization = "complete", differential_expression = "pending")

  resetProtNormReactiveState(
    normData = norm_data,
    workflowData = workflow_data
  )

  expect_false(norm_data$normalization_complete)
  expect_false(norm_data$ruv_complete)
  expect_false(norm_data$correlation_filtering_complete)
  expect_null(norm_data$normalized_protein_obj)
  expect_null(norm_data$ruv_normalized_obj)
  expect_null(norm_data$correlation_filtered_obj)
  expect_null(norm_data$best_k)
  expect_null(norm_data$control_genes_index)
  expect_null(norm_data$correlation_vector)
  expect_null(norm_data$correlation_threshold)
  expect_null(norm_data$final_qc_plot)
  expect_null(norm_data$final_filtering_plot)
  expect_null(norm_data$post_norm_filtering_plot)
  expect_null(norm_data$filtering_summary_text)
  expect_null(norm_data$ruv_optimization_result)
  expect_equal(
    norm_data$qc_plots$post_normalization,
    list(pca = NULL, density = NULL, rle = NULL, correlation = NULL)
  )
  expect_equal(
    norm_data$qc_plots$ruv_corrected,
    list(pca = NULL, density = NULL, rle = NULL, correlation = NULL)
  )
  expect_null(workflow_data$ruv_normalised_for_da_analysis_obj)
  expect_equal(workflow_data$tab_status$normalization, "pending")
  expect_equal(workflow_data$tab_status$differential_expression, "disabled")
})

test_that("resetProtNormOutputs restores the default reset text", {
  output <- new.env(parent = emptyenv())

  resetProtNormOutputs(
    output = output,
    ruvMode = "automatic",
    groupingVariable = "group",
    renderTextFn = function(text) text
  )

  expect_equal(
    output$correlation_filter_summary,
    "No correlation filtering applied yet"
  )
  expect_equal(
    output$filtering_summary_text,
    "Filtering summary will be available after normalization and RUV correction."
  )
  expect_match(
    output$ruv_optimization_summary,
    "Run normalization to see RUV optimization results",
    fixed = TRUE
  )
})

test_that("runProtNormResetWorkflow sequences reset helpers and notification", {
  workflow_data <- new.env(parent = emptyenv())
  norm_data <- new.env(parent = emptyenv())
  output <- new.env(parent = emptyenv())
  calls <- list()
  notification <- NULL
  messages <- character()

  result <- runProtNormResetWorkflow(
    workflowData = workflow_data,
    normData = norm_data,
    output = output,
    ruvMode = "manual",
    groupingVariable = "condition",
    revertStateManagerFn = function(workflowData, messageFn) {
      calls$revert <<- list(workflowData = workflowData)
      list(previousState = "sample_filtered", revertedS4 = "state")
    },
    resetReactiveStateFn = function(normData, workflowData) {
      calls$resetReactive <<- list(normData = normData, workflowData = workflowData)
      invisible(NULL)
    },
    resetOutputsFn = function(output, ruvMode, groupingVariable) {
      calls$resetOutputs <<- list(
        output = output,
        ruvMode = ruvMode,
        groupingVariable = groupingVariable
      )
      invisible(NULL)
    },
    buildNotificationMessageFn = function(previousState) {
      calls$buildNotification <<- previousState
      sprintf("reset to %s", previousState)
    },
    showNotificationFn = function(message, type, duration) {
      notification <<- list(message = message, type = type, duration = duration)
    },
    messageFn = function(text) {
      messages <<- c(messages, text)
    }
  )

  expect_identical(calls$revert$workflowData, workflow_data)
  expect_identical(calls$resetReactive$normData, norm_data)
  expect_identical(calls$resetOutputs$output, output)
  expect_equal(calls$resetOutputs$ruvMode, "manual")
  expect_equal(calls$resetOutputs$groupingVariable, "condition")
  expect_equal(calls$buildNotification, "sample_filtered")
  expect_equal(notification$message, "reset to sample_filtered")
  expect_equal(notification$type, "warning")
  expect_equal(notification$duration, 5)
  expect_true(any(grepl("reset completed successfully", messages, fixed = TRUE)))
  expect_equal(result$previousState, "sample_filtered")
})

test_that("handleProtNormResetError logs and notifies consistently", {
  messages <- character()
  notification <- NULL

  handleProtNormResetError(
    error = structure(list(message = "boom"), class = c("simpleError", "error", "condition")),
    showNotificationFn = function(message, type, duration) {
      notification <<- list(message = message, type = type, duration = duration)
    },
    messageFn = function(text) {
      messages <<- c(messages, text)
    }
  )

  expect_true(any(grepl("boom", messages, fixed = TRUE)))
  expect_equal(notification$message, "Error resetting normalization: boom")
  expect_equal(notification$type, "error")
  expect_equal(notification$duration, 10)
})

test_that("runProtNormResetObserver delegates to reset workflow and error handler", {
  calls <- list()

  result <- runProtNormResetObserver(
    workflowData = "workflow",
    normData = "norm",
    output = "output",
    ruvMode = "automatic",
    groupingVariable = "group",
    runResetWorkflowFn = function(workflowData, normData, output, ruvMode, groupingVariable) {
      calls$run <<- list(
        workflowData = workflowData,
        normData = normData,
        output = output,
        ruvMode = ruvMode,
        groupingVariable = groupingVariable
      )
      "reset_done"
    },
    handleResetErrorFn = function(error) {
      calls$error <<- error
    },
    messageFn = function(...) NULL
  )

  expect_equal(calls$run$workflowData, "workflow")
  expect_equal(calls$run$normData, "norm")
  expect_equal(calls$run$output, "output")
  expect_equal(calls$run$ruvMode, "automatic")
  expect_equal(calls$run$groupingVariable, "group")
  expect_null(calls$error)
  expect_equal(result, "reset_done")
})

test_that("checkProtNormMemoryUsage reports memory usage and warns only above threshold", {
  messages <- character()
  warnings <- character()

  low_result <- checkProtNormMemoryUsage(
    threshold_gb = 8,
    context = "low",
    gcFn = function() matrix(c(0, 1024, 0, 0), ncol = 2, byrow = TRUE),
    warningFn = function(text) {
      warnings <<- c(warnings, text)
    },
    messageFn = function(text) {
      messages <<- c(messages, text)
    }
  )
  high_result <- checkProtNormMemoryUsage(
    threshold_gb = 1,
    context = "high",
    gcFn = function() matrix(c(0, 4096, 0, 0), ncol = 2, byrow = TRUE),
    warningFn = function(text) {
      warnings <<- c(warnings, text)
    },
    messageFn = function(text) {
      messages <<- c(messages, text)
    }
  )

  expect_equal(low_result, 1)
  expect_equal(high_result, 4)
  expect_true(any(grepl("MEMORY CHECK \\[low\\]", messages)))
  expect_true(any(grepl("HIGH MEMORY WARNING \\[high\\]", messages)))
  expect_true(any(grepl("HIGH MEMORY WARNING \\[high\\]", warnings)))
})

test_that("createProtNormReactiveState initializes the normalization reactive defaults", {
  state <- createProtNormReactiveState(
    reactiveValuesFn = function(...) list(...)
  )

  expect_false(state$pre_norm_qc_generated)
  expect_false(state$normalization_complete)
  expect_false(state$ruv_complete)
  expect_false(state$correlation_filtering_complete)
  expect_null(state$QC_composite_figure)
  expect_equal(state$plot_refresh_trigger, 0)
  expect_null(state$normalized_protein_obj)
  expect_null(state$ruv_normalized_obj)
  expect_null(state$correlation_filtered_obj)
  expect_null(state$ruv_optimization_result)
  expect_named(
    state$qc_plots,
    c("post_filtering", "post_normalization", "ruv_corrected")
  )
  expect_named(
    state$qc_plots$post_filtering,
    c("pca", "density", "rle", "correlation")
  )
})

test_that("registerProtNormServerObservers wires the observer shell through helper delegates", {
  calls <- list(
    observe = 0,
    events = character()
  )
  input <- list(
    color_variable = "group",
    shape_variable = "batch",
    run_normalization = 1,
    apply_correlation_filter = 1,
    skip_correlation_filter = 1,
    reset_normalization = 1,
    export_filtered_session = 1,
    ruv_mode = "manual"
  )
  output <- new.env(parent = emptyenv())
  session <- "session"
  selected_tab <- function() "normalization"

  registerProtNormServerObservers(
    input = input,
    output = output,
    session = session,
    selectedTab = selected_tab,
    workflowData = list(design_matrix = "design", state_manager = "state"),
    normData = "norm_data",
    experimentPaths = list(protein_qc_dir = "/tmp/qc", source_dir = "/tmp/source"),
    omicType = "proteomics",
    experimentLabel = "label",
    generatePreNormalizationQcFn = function() NULL,
    generatePostNormalizationQcFn = function(x) x,
    generateRuvCorrectedQcFn = function(x) x,
    getPlotAestheticsFn = function() list(color_var = "group", shape_var = "batch"),
    getRuvGroupingVariableFn = function() "group",
    checkMemoryUsageFn = "check_memory",
    updateDesignDrivenChoicesFn = function(session, designMatrix) {
      calls$observe <<- calls$observe + 1
      calls$design <<- list(session = session, designMatrix = designMatrix)
    },
    runTabEntryWorkflowFn = function(selectedTab, workflowData, normData, generatePreNormalizationQcFn) {
      calls$tab_entry <<- list(
        selectedTab = selectedTab,
        workflowData = workflowData,
        normData = normData,
        generatePreNormalizationQcFn = generatePreNormalizationQcFn
      )
    },
    regenerateQcForAestheticChangeFn = function(normData, generatePreNormalizationQcFn, generatePostNormalizationQcFn, generateRuvCorrectedQcFn) {
      calls$regenerate <<- list(
        normData = normData,
        generatePreNormalizationQcFn = generatePreNormalizationQcFn,
        generatePostNormalizationQcFn = generatePostNormalizationQcFn,
        generateRuvCorrectedQcFn = generateRuvCorrectedQcFn
      )
    },
    runNormalizationWorkflowFn = function(input, workflowData, normData, experimentPaths, omicType, experimentLabel, checkMemoryUsageFn, generatePostNormalizationQcFn, generateRuvCorrectedQcFn, getRuvGroupingVariableFn) {
      calls$normalize <<- list(
        input = input,
        workflowData = workflowData,
        normData = normData,
        experimentPaths = experimentPaths,
        omicType = omicType,
        experimentLabel = experimentLabel,
        checkMemoryUsageFn = checkMemoryUsageFn,
        generatePostNormalizationQcFn = generatePostNormalizationQcFn,
        generateRuvCorrectedQcFn = generateRuvCorrectedQcFn,
        getRuvGroupingVariableFn = getRuvGroupingVariableFn
      )
    },
    runApplyCorrelationObserverFn = function(input, output, workflowData, normData, experimentPaths, omicType, experimentLabel, getRuvGroupingVariableFn, getPlotAestheticsFn, checkMemoryUsageFn) {
      calls$apply <<- list(
        input = input,
        output = output,
        workflowData = workflowData,
        normData = normData,
        experimentPaths = experimentPaths,
        omicType = omicType,
        experimentLabel = experimentLabel,
        getRuvGroupingVariableFn = getRuvGroupingVariableFn,
        getPlotAestheticsFn = getPlotAestheticsFn,
        checkMemoryUsageFn = checkMemoryUsageFn
      )
    },
    runSkipCorrelationObserverFn = function(output, workflowData, normData, experimentPaths, omicType, experimentLabel, getPlotAestheticsFn, checkMemoryUsageFn) {
      calls$skip <<- list(
        output = output,
        workflowData = workflowData,
        normData = normData,
        experimentPaths = experimentPaths,
        omicType = omicType,
        experimentLabel = experimentLabel,
        getPlotAestheticsFn = getPlotAestheticsFn,
        checkMemoryUsageFn = checkMemoryUsageFn
      )
    },
    runResetObserverFn = function(workflowData, normData, output, ruvMode, groupingVariable) {
      calls$reset <<- list(
        workflowData = workflowData,
        normData = normData,
        output = output,
        ruvMode = ruvMode,
        groupingVariable = groupingVariable
      )
    },
    runExportObserverFn = function(input, workflowData, normData, experimentPaths) {
      calls$export <<- list(
        input = input,
        workflowData = workflowData,
        normData = normData,
        experimentPaths = experimentPaths
      )
    },
    observeFn = function(expr) {
      eval(substitute(expr), parent.frame())
      invisible("observe")
    },
    observeEventFn = function(eventExpr, handlerExpr, ignoreInit = FALSE) {
      calls$events <<- c(calls$events, paste(deparse(substitute(eventExpr)), collapse = ""))
      force(ignoreInit)
      eval(substitute(handlerExpr), parent.frame())
      invisible("observeEvent")
    }
  )

  expect_equal(calls$observe, 1)
  expect_equal(calls$design$session, session)
  expect_equal(calls$design$designMatrix, "design")
  expect_equal(calls$tab_entry$selectedTab, "normalization")
  expect_equal(calls$regenerate$normData, "norm_data")
  expect_equal(calls$normalize$checkMemoryUsageFn, "check_memory")
  expect_equal(calls$normalize$omicType, "proteomics")
  expect_equal(calls$apply$experimentLabel, "label")
  expect_equal(calls$skip$checkMemoryUsageFn, "check_memory")
  expect_equal(calls$reset$ruvMode, "manual")
  expect_equal(calls$reset$groupingVariable, "group")
  expect_equal(calls$export$normData, "norm_data")
  expect_length(calls$events, 7)
})
