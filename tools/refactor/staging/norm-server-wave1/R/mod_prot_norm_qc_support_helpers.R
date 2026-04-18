buildProtNormQcImagePayload <- function(proteinQcDir, filename, altText) {
  img_path <- if (!is.null(proteinQcDir)) {
    file.path(proteinQcDir, filename)
  } else {
    ""
  }

  if (img_path != "" && file.exists(img_path)) {
    return(list(
      src = img_path,
      contentType = "image/png",
      width = "100%",
      height = "auto",
      alt = altText
    ))
  }

  list(src = "", alt = "Plot not generated yet")
}

getProtNormFilteringSummaryText <- function(filteringSummaryText) {
  if (!is.null(filteringSummaryText)) {
    return(filteringSummaryText)
  }

  "Filtering summary will be available after normalization and RUV correction."
}

buildProtNormRuvOptimizationSummary <- function(ruvOptimizationResult, ruvMode, groupingVariable) {
  if (is.null(ruvOptimizationResult)) {
    return("Run normalization to see RUV optimization results")
  }

  if (ruvMode == "automatic") {
    return(sprintf(
      paste(
        "Automatic RUV Optimization Results:\n\n* Best percentage: %.1f%%\n* Best k value: %d\n* Separation score: %.4f\n* Composite score: %.4f\n* Control genes: %d\n* RUV grouping variable: %s\n* Separation metric: %s\n* K penalty weight: %.1f\n* Adaptive penalty: %s\n* Sample size: %s",
        collapse = ""
      ),
      ruvOptimizationResult$best_percentage,
      ruvOptimizationResult$best_k,
      ifelse(is.null(ruvOptimizationResult$best_separation_score), 0, ruvOptimizationResult$best_separation_score),
      ifelse(is.null(ruvOptimizationResult$best_composite_score), 0, ruvOptimizationResult$best_composite_score),
      sum(ruvOptimizationResult$best_control_genes_index, na.rm = TRUE),
      groupingVariable,
      ruvOptimizationResult$separation_metric_used,
      ifelse(is.null(ruvOptimizationResult$k_penalty_weight), 0, ruvOptimizationResult$k_penalty_weight),
      ifelse(is.null(ruvOptimizationResult$adaptive_k_penalty_used), FALSE, ruvOptimizationResult$adaptive_k_penalty_used),
      ifelse(is.null(ruvOptimizationResult$sample_size), "N/A", ruvOptimizationResult$sample_size)
    ))
  }

  sprintf(
    paste(
      "Manual RUV Parameters:\n\n* Selected percentage: %.1f%%\n* Selected k value: %d\n* Control genes: %d\n* RUV grouping variable: %s\n* Mode: Manual selection",
      collapse = ""
    ),
    ruvOptimizationResult$best_percentage,
    ruvOptimizationResult$best_k,
    sum(ruvOptimizationResult$best_control_genes_index, na.rm = TRUE),
    groupingVariable
  )
}

prepareProtNormOptimizationResultsTable <- function(ruvOptimizationResult) {
  if (is.null(ruvOptimizationResult) || is.null(ruvOptimizationResult$optimization_results)) {
    return(list(
      hasResults = FALSE,
      results = data.frame(Message = "Run normalization to see optimization results"),
      bestPercentage = NULL
    ))
  }

  results_df <- ruvOptimizationResult$optimization_results

  if ("separation_score" %in% colnames(results_df)) {
    results_df$separation_score <- round(results_df$separation_score, 4)
  }
  if ("composite_score" %in% colnames(results_df)) {
    results_df$composite_score <- round(results_df$composite_score, 4)
  }

  list(
    hasResults = TRUE,
    results = results_df,
    bestPercentage = ruvOptimizationResult$best_percentage
  )
}

resolveProtNormFinalQcRenderState <- function(finalQcPlot, finalFilteringPlot) {
  if (!is.null(finalQcPlot) && !is.null(finalFilteringPlot)) {
    return(list(
      mode = "combined",
      finalQcPlot = finalQcPlot,
      finalFilteringPlot = finalFilteringPlot
    ))
  }

  if (!is.null(finalFilteringPlot)) {
    return(list(
      mode = "filtering_only",
      finalFilteringPlot = finalFilteringPlot
    ))
  }

  if (!is.null(finalQcPlot)) {
    return(list(
      mode = "pca_only",
      finalQcPlot = finalQcPlot
    ))
  }

  list(mode = "placeholder")
}

getProtNormRuvCanonicalCorrelationPlot <- function(ruvOptimizationResult) {
  if (!is.null(ruvOptimizationResult) && !is.null(ruvOptimizationResult$best_cancor_plot)) {
    return(ruvOptimizationResult$best_cancor_plot)
  }

  NULL
}

getProtNormDefaultCorrelationFilterSummaryText <- function() {
  "No correlation filtering applied yet"
}

updateProtNormDesignDrivenChoices <- function(session, designMatrix, updateSelectInputFn = shiny::updateSelectInput, messageFn = message) {
  if (is.null(designMatrix)) {
    return(invisible(NULL))
  }

  design_cols <- colnames(designMatrix)

  plot_available_vars <- intersect(
    design_cols,
    c("group", "factor1", "factor2", "batch", "technical_replicate_id", "sample_id")
  )

  if (length(plot_available_vars) > 0) {
    default_plot_var <- if ("group" %in% plot_available_vars) "group" else plot_available_vars[1]

    updateSelectInputFn(
      session,
      "color_variable",
      choices = plot_available_vars,
      selected = default_plot_var
    )
    updateSelectInputFn(
      session,
      "shape_variable",
      choices = plot_available_vars,
      selected = default_plot_var
    )
  }

  ruv_available_vars <- intersect(design_cols, c("group", "factor1", "factor2", "batch"))

  if (length(ruv_available_vars) > 0) {
    default_ruv <- if ("batch" %in% ruv_available_vars && any(!is.na(designMatrix$batch))) {
      messageFn("*** RUV: Batch column detected with values - suggesting 'batch' as RUV variable to set to shapes ***")
      "group"
    } else if ("group" %in% ruv_available_vars) {
      "group"
    } else {
      ruv_available_vars[1]
    }

    updateSelectInputFn(
      session,
      "ruv_grouping_variable",
      choices = ruv_available_vars,
      selected = default_ruv
    )

    messageFn(sprintf("*** RUV: Updated grouping variable choices to: %s ***", paste(ruv_available_vars, collapse = ", ")))
  }

  invisible(NULL)
}

getProtNormPlotAesthetics <- function(colorVariable, shapeVariable) {
  list(
    color_var = if (is.null(colorVariable) || colorVariable == "") "group" else colorVariable,
    shape_var = if (is.null(shapeVariable) || shapeVariable == "") "group" else shapeVariable
  )
}

getProtNormRuvGroupingVariable <- function(ruvGroupingVariable) {
  if (is.null(ruvGroupingVariable) || ruvGroupingVariable == "") {
    return("group")
  }

  ruvGroupingVariable
}

buildProtNormLabelPlot <- function(title) {
  ggplot2::ggplot() +
    ggplot2::annotate("text", x = 0, y = 0.5, label = title, size = 5, hjust = 0) +
    ggplot2::xlim(0, 1) +
    ggplot2::theme_void() +
    ggplot2::theme(
      plot.margin = ggplot2::margin(5, 5, 5, 5),
      panel.background = ggplot2::element_blank()
    )
}

buildProtNormTitlePlot <- function(title) {
  ggplot2::ggplot() +
    ggplot2::annotate("text", x = 0.5, y = 0.5, label = title, size = 6, fontface = "bold", hjust = 0.5) +
    ggplot2::xlim(0, 1) +
    ggplot2::theme_void() +
    ggplot2::theme(
      plot.margin = ggplot2::margin(5, 5, 10, 5),
      panel.background = ggplot2::element_blank()
    )
}

loadProtNormImageAsPlot <- function(filePath, messageFn = message) {
  if (is.na(filePath) || !file.exists(filePath)) {
    return(ggplot2::ggplot() + ggplot2::theme_void())
  }

  tryCatch({
    img <- png::readPNG(filePath)
    g <- grid::rasterGrob(img, interpolate = TRUE)
    ggplot2::ggplot() +
      ggplot2::annotation_custom(g, xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf) +
      ggplot2::theme_void()
  }, error = function(e) {
    messageFn(sprintf("   [generateCompositeFromFiles] Could not load image: %s", filePath))
    ggplot2::ggplot() + ggplot2::theme_void()
  })
}

generateProtNormCompositeFromFiles <- function(
  plotFiles,
  ncol = 3,
  rowLabels = NULL,
  columnLabels = NULL,
  messageFn = message,
  warningFn = warning
) {
  messageFn(sprintf("   [generateCompositeFromFiles] Generating composite from %d files...", length(plotFiles)))

  if (!requireNamespace("patchwork", quietly = TRUE)) {
    warningFn("patchwork package required for composite generation")
    return(NULL)
  }
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    warningFn("ggplot2 package required for composite generation")
    return(NULL)
  }
  if (!requireNamespace("png", quietly = TRUE)) {
    warningFn("png package required for composite generation")
    return(NULL)
  }

  tryCatch({
    n_files <- length(plotFiles)
    n_plot_types <- n_files / ncol

    if (is.null(rowLabels)) {
      all_labels <- letters[1:n_files]
      rowLabels <- split(paste0(all_labels, ")"), rep(1:n_plot_types, each = ncol))
      names(rowLabels) <- paste0("row", seq_len(n_plot_types))
    }

    plot_sections <- list()
    height_values <- c()

    if (!is.null(columnLabels)) {
      if (length(columnLabels) != ncol) {
        warningFn("Number of column labels does not match ncol")
      } else {
        title_plots <- lapply(columnLabels, buildProtNormTitlePlot)
        plot_sections <- append(plot_sections, list(
          patchwork::wrap_plots(title_plots, ncol = ncol)
        ))
        height_values <- c(height_values, 0.2)
        messageFn("   [generateCompositeFromFiles] Added column titles")
      }
    }

    row_names <- names(rowLabels)

    for (i in seq_along(row_names)) {
      row_name <- row_names[i]
      labels <- rowLabels[[row_name]]

      start_idx <- (i - 1) * ncol + 1
      end_idx <- min(i * ncol, n_files)
      row_files <- plotFiles[start_idx:end_idx]

      has_files <- any(!is.na(row_files) & sapply(row_files, function(f) !is.na(f) && file.exists(f)))

      if (has_files || row_name == "cancor") {
        label_plots <- lapply(labels, buildProtNormLabelPlot)
        image_plots <- lapply(row_files, loadProtNormImageAsPlot, messageFn = messageFn)

        plot_sections <- append(plot_sections, list(
          patchwork::wrap_plots(label_plots, ncol = ncol),
          patchwork::wrap_plots(image_plots, ncol = ncol)
        ))
        height_values <- c(height_values, 0.1, 1)

        messageFn(sprintf("   [generateCompositeFromFiles] Added row: %s", row_name))
      } else {
        messageFn(sprintf("   [generateCompositeFromFiles] Skipping empty row: %s", row_name))
      }
    }

    if (length(plot_sections) == 0) {
      warningFn("No valid plot sections to combine")
      return(NULL)
    }

    messageFn("   [generateCompositeFromFiles] Combining plot sections...")
    combined_plot <- patchwork::wrap_plots(plot_sections, ncol = 1) +
      patchwork::plot_layout(heights = height_values)

    plot_width <- 4 + (ncol * 3)
    plot_height <- 4 + (length(height_values) * 2)

    rm(plot_sections)
    gc()

    list(plot = combined_plot, width = plot_width, height = plot_height)
  }, error = function(e) {
    messageFn(paste("   [generateCompositeFromFiles] Error:", e$message))
    NULL
  })
}

shouldProtNormAutoGeneratePreQc <- function(selectedTabValue, currentState, preNormQcGenerated) {
  valid_states_for_norm_tab <- c("protein_replicate_filtered", "normalized", "ruv_corrected", "correlation_filtered")

  !is.null(selectedTabValue) &&
    selectedTabValue == "normalization" &&
    currentState %in% valid_states_for_norm_tab &&
    !preNormQcGenerated
}

regenerateProtNormQcForAestheticChange <- function(
  normData,
  generatePreNormalizationQcFn,
  generatePostNormalizationQcFn,
  generateRuvCorrectedQcFn,
  messageFn = message
) {
  if (normData$pre_norm_qc_generated) {
    messageFn("Regenerating pre-normalization QC plots with new aesthetics...")
    tryCatch({
      generatePreNormalizationQcFn()
    }, error = function(e) {
      messageFn(paste("Error regenerating pre-normalization QC:", e$message))
    })
  }

  if (normData$normalization_complete && !is.null(normData$normalized_protein_obj)) {
    messageFn("Regenerating post-normalization QC plots with new aesthetics...")
    tryCatch({
      generatePostNormalizationQcFn(normData$normalized_protein_obj)
    }, error = function(e) {
      messageFn(paste("Error regenerating post-normalization QC:", e$message))
    })
  }

  if (normData$ruv_complete && !is.null(normData$ruv_normalized_obj)) {
    messageFn("Regenerating RUV-corrected QC plots with new aesthetics...")
    tryCatch({
      generateRuvCorrectedQcFn(normData$ruv_normalized_obj)
    }, error = function(e) {
      messageFn(paste("Error regenerating RUV-corrected QC:", e$message))
    })
  }

  invisible(NULL)
}

initializeProtNormQcPlotPaths <- function(qcPlotPaths = NULL) {
  if (!is.null(qcPlotPaths)) {
    return(qcPlotPaths)
  }

  list(
    post_filtering = list(),
    post_normalization = list(),
    ruv_corrected = list()
  )
}

saveProtNormQcPlotArtifact <- function(
  qcDir,
  filename,
  plotObject,
  width,
  height,
  dpi = 150,
  savePlotFn = ggplot2::ggsave
) {
  if (is.null(qcDir)) {
    return(NULL)
  }

  artifactPath <- file.path(qcDir, filename)
  savePlotFn(artifactPath, plotObject, width = width, height = height, dpi = dpi)
  artifactPath
}

recordProtNormQcPlotPath <- function(qcPlotPaths, stage, plotType, path) {
  qcPlotPaths <- initializeProtNormQcPlotPaths(qcPlotPaths)

  if (!is.null(path)) {
    qcPlotPaths[[stage]][[plotType]] <- path
  }

  qcPlotPaths
}

buildProtNormDensityPlot <- function(
  s4Object,
  aesthetics,
  plotPcaFn = plotPca,
  plotPcaBoxFn = plotPcaBox
) {
  pcaPlot <- plotPcaFn(
    s4Object,
    grouping_variable = aesthetics$color_var,
    label_column = "",
    shape_variable = aesthetics$shape_var,
    title = "",
    font_size = 8
  )

  plotPcaBoxFn(
    pcaPlot,
    grouping_variable = aesthetics$color_var,
    show_legend = TRUE
  )
}

buildProtNormCorrelationPlot <- function(
  s4Object,
  colorVar,
  isRunningFn = shiny::isRunning,
  withProgressFn = shiny::withProgress,
  incProgressFn = shiny::incProgress,
  plotPearsonFn = plotPearson
) {
  if (isRunningFn()) {
    return(withProgressFn(
      message = "Generating Pearson correlation plot...",
      detail = "Calculating pairwise correlations...",
      value = 0.5,
      {
        result <- plotPearsonFn(
          s4Object,
          correlation_group = colorVar,
          exclude_pool_samples = TRUE
        )
        incProgressFn(0.5)
        result
      }
    ))
  }

  plotPearsonFn(
    s4Object,
    correlation_group = colorVar,
    exclude_pool_samples = TRUE
  )
}

