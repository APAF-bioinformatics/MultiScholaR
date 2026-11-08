# KIT-105 W3: proteomics' normalisation core. Back-references #10 to #22 (the functions the DIA reference path, the
# DIA sessions and the lipid/metab RUV declaration call), their module-file callees and the DIA_LIMPA pass-through,
# moved unchanged and with their names kept from R/mod_prot_norm_{workflow,ruv,qc_generation,qc_support,correlation,
# session}_helpers.R (design 2.2, 4.7). The tab's adapters, renderers and observers stay in its module files; KIT-111
# moves the run out of the GUI process.

# --- RUV-III parameters and results (#11, #12, #10; the manual result #11 calls) ---

buildProtNormManualRuvResult <- function(
  percentageAsNegCtrl,
  ruvK,
  controlGenesIndex,
  cancorPlot
) {
  list(
    best_percentage = percentageAsNegCtrl,
    best_k = ruvK,
    best_control_genes_index = controlGenesIndex,
    best_cancor_plot = cancorPlot,
    optimization_results = data.frame(
      percentage = percentageAsNegCtrl,
      separation_score = NA,
      best_k = ruvK,
      composite_score = NA,
      num_controls = sum(controlGenesIndex, na.rm = TRUE),
      valid_plot = TRUE
    ),
    separation_metric_used = "manual",
    k_penalty_weight = NA,
    adaptive_k_penalty_used = FALSE
  )
}

validateProtNormRuvK <- function(ruvK, paramLabel = "ruv_number_k") {
  if (is.null(ruvK) ||
      length(ruvK) != 1L ||
      is.na(ruvK) ||
      !is.numeric(ruvK) ||
      ruvK < 1) {
    stop(sprintf(
      "The %s = %s value is invalid.",
      paramLabel,
      paste(ruvK, collapse = ",")
    ), call. = FALSE)
  }

  as.integer(ruvK)
}

# Kernel-only RUV parameter determination (KIT-004): the same control selection, k validation and
# canonical-correlation plot as resolveProtNormRuvParameters, without progress, messages, holders or
# the config_list audit trail. The module keeps the adapter below; the operation runner calls this.
resolveProtRuvParameters <- function(
  normalizedS4,
  settings,
  ruvGroupingVariable,
  findBestNegCtrlPercentageFn = findBestNegCtrlPercentage,
  getNegCtrlProtAnovaFn = getNegCtrlProtAnova,
  ruvCancorFn = ruvCancor,
  buildManualRuvResultFn = buildProtNormManualRuvResult
) {
  if (identical(settings$ruv_mode, "automatic")) {
    percentage_min <- if (is.null(settings$auto_percentage_min) || is.na(settings$auto_percentage_min)) 1 else settings$auto_percentage_min
    percentage_max <- if (is.null(settings$auto_percentage_max) || is.na(settings$auto_percentage_max)) 20 else settings$auto_percentage_max
    if (percentage_min >= percentage_max) {
      stop("Minimum percentage must be less than maximum percentage")
    }
    result <- findBestNegCtrlPercentageFn(
      normalised_protein_matrix_obj = normalizedS4,
      percentage_range = seq(percentage_min, percentage_max, by = 1),
      separation_metric = settings$separation_metric,
      k_penalty_weight = settings$k_penalty_weight,
      adaptive_k_penalty = settings$adaptive_k_penalty,
      ruv_grouping_variable = ruvGroupingVariable,
      verbose = TRUE
    )
    ruv_k <- validateProtNormRuvK(result$best_k, "ruv_k")
    return(list(result = result, k = ruv_k, controls = result$best_control_genes_index,
      percentage = result$best_percentage))
  }
  percentage <- settings$ruv_percentage
  ruv_k <- if (is.null(settings$ruv_k) || is.na(settings$ruv_k)) 3 else settings$ruv_k
  ruv_k <- validateProtNormRuvK(ruv_k, "ruv_k")
  controls <- getNegCtrlProtAnovaFn(
    normalizedS4,
    ruv_grouping_variable = ruvGroupingVariable,
    percentage_as_neg_ctrl = percentage,
    ruv_qval_cutoff = 0.05,
    ruv_fdr_method = "BH"
  )
  # The adapter tolerates a failed canonical-correlation plot and leaves the result unset.
  result <- tryCatch({
    cancor_plot <- ruvCancorFn(normalizedS4, ctrl = controls, num_components_to_impute = 2,
      ruv_grouping_variable = ruvGroupingVariable)
    buildManualRuvResultFn(percentageAsNegCtrl = percentage, ruvK = ruv_k,
      controlGenesIndex = controls, cancorPlot = cancor_plot)
  }, error = function(e) {
    message(paste("Warning: Could not generate canonical correlation plot for manual mode:", e$message))
    NULL
  })
  list(result = result, k = ruv_k, controls = controls, percentage = result$best_percentage)
}

buildProtNormSkippedRuvResult <- function() {
  list(
    best_percentage = NA,
    best_k = NA,
    best_control_genes_index = NA,
    best_cancor_plot = NULL,
    ruv_skipped = TRUE,
    skip_reason = "User selected skip due to dataset constraints"
  )
}

# --- QC plot products (#13 to #16) ---

#' Remove captured plotting environments while preserving plot data and layers.
#' @param plotObject Existing ggplot or patchwork object.
#' @return Plot without ownership of producer-local scientific inputs.
#' @noRd
stripProtNormQcPlotOwners <- function(plotObject) {
  if (inherits(plotObject, c("gg", "ggplot"))) {
    plotObject$plot_env <- rlang::base_env()
  }
  if (!is.null(plotObject$patches) && is.list(plotObject$patches$plots)) {
    plotObject$patches$plots <- lapply(
      plotObject$patches$plots,
      stripProtNormQcPlotOwners
    )
  }
  plotObject
}

#' Produce one existing normalization QC plot as a bounded reusable product.
#' @param s4Object Existing ProteinQuantitativeData.
#' @param aesthetics Existing color/shape selection.
#' @param plotType One of pca, rle, density or correlation.
#' @param plotPcaFn,plotRleFn,buildDensityFn,buildCorrelationFn Existing producers.
#' @return Product containing ordered population metadata, plot data and baseline plot.
#' @noRd
buildProtNormQcPlotProduct <- function(
  s4Object,
  aesthetics,
  plotType,
  plotPcaFn = plotPca,
  plotRleFn = plotRle,
  buildDensityFn = buildProtNormDensityPlot,
  buildCorrelationFn = buildProtNormCorrelationPlot
) {
  if (!methods::is(s4Object, "ProteinQuantitativeData") ||
      !is.list(aesthetics) ||
      !identical(sort(names(aesthetics)), c("color_var", "shape_var")) ||
      !is.character(plotType) || length(plotType) != 1L || is.na(plotType) ||
      !plotType %in% c("pca", "rle", "density", "correlation")) {
    cli::cli_abort("invalid protein normalization QC product request",
      class = "saved_entity_invalid")
  }
  previous_plot <- ggplot2::get_last_plot()
  on.exit(ggplot2::set_last_plot(previous_plot), add = TRUE)
  custom_producer <- switch(plotType,
    pca = !identical(plotPcaFn, plotPca),
    rle = !identical(plotRleFn, plotRle),
    density = !identical(buildDensityFn, buildProtNormDensityPlot),
    correlation = !identical(buildCorrelationFn, buildProtNormCorrelationPlot))
  quant <- methods::slot(s4Object, "protein_quant_table")
  key <- methods::slot(s4Object, "protein_id_column")
  design <- methods::slot(s4Object, "design_matrix")
  sample_id <- methods::slot(s4Object, "sample_id")
  samples <- as.character(design[[sample_id]])
  annotations <- design[match(samples, design[[sample_id]]), , drop = FALSE]
  plot <- switch(plotType,
    pca = plotPcaFn(s4Object, grouping_variable = aesthetics$color_var,
      label_column = "", shape_variable = aesthetics$shape_var,
      title = "", font_size = 8),
    rle = plotRleFn(s4Object, group = aesthetics$color_var,
      yaxis_limit = c(-6, 6)),
    density = buildDensityFn(s4Object, aesthetics),
    correlation = buildCorrelationFn(s4Object, aesthetics$color_var)
  )
  data <- if (!is.null(plot$data) && is.data.frame(plot$data)) {
    as.data.frame(plot$data)
  } else {
    data.frame()
  }
  population_role <- if (identical(plotType, "density")) {
    "density_independent_pca_call"
  } else {
    paste0(plotType, "_existing_producer")
  }
  coordinate_limits <- if (!is.null(plot$coordinates$limits)) {
    plot$coordinates$limits
  } else {
    list(x = NULL, y = NULL)
  }
  product <- list(
    schema = "prot_norm_qc_plot_product_v1",
    plot_type = plotType,
    population_role = population_role,
    feature_ids = as.character(quant[[key]]),
    sample_ids = samples,
    annotations = annotations,
    scientific_parameters = list(
      color_variable = aesthetics$color_var,
      shape_variable = aesthetics$shape_var,
      pca_cv_percentile = 0.90,
      rle_yaxis_limit = c(-6, 6),
      correlation_exclude_pool_samples = TRUE
    ),
    render_metadata = list(
      x_label = if (is.null(plot$labels$x)) "" else as.character(plot$labels$x),
      y_label = if (is.null(plot$labels$y)) "" else as.character(plot$labels$y),
      x_limits = coordinate_limits$x,
      y_limits = coordinate_limits$y
    ),
    data = data,
    plot = NULL
  )
  # The reconstruction frame receives calculated display data, never this S4 or
  # the native plot's captured matrices/fits. Preserve explicitly injected seams.
  if (custom_producer) {
    product$plot <- stripProtNormQcPlotOwners(plot)
  } else {
    plot <- NULL
    product$plot <- stripProtNormQcPlotOwners(rebuildWorkflowDiaQcPlot(product))
  }
  product
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

# --- The composite QC figure and the workflow's completion (#21, #22; the composite's callees) ---

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
    stripProtNormQcPlotOwners(ggplot2::ggplot() +
      ggplot2::annotation_custom(g, xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf) +
      ggplot2::theme_void())
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
      all_labels <- letters[seq_len(n_files)] # KIT-064 (roar BUG-16): as the lipid copy already does
      rowLabels <- split(paste0(all_labels, ")"), rep(seq_len(n_plot_types), each = ncol))
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
        # The single RUV control plot takes a shorter row than the three-panel rows.
        height_values <- c(height_values, 0.1, if (row_name == "cancor") 0.75 else 1)

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

resolveProtNormCompositeFigureInputs <- function(ruvMode, qcDir) {
  if (ruvMode == "skip") {
    ncol_composite <- 2
    file_names <- c(
      "pre_norm_pca.png", "post_norm_pca.png",
      "pre_norm_density.png", "post_norm_density.png",
      "pre_norm_rle.png", "post_norm_rle.png",
      "pre_norm_correlation.png", "post_norm_correlation.png"
    )
    row_labels <- list(
      pca = c("a)", "b)"),
      density = c("c)", "d)"),
      rle = c("e)", "f)"),
      correlation = c("g)", "h)")
    )
    column_labels <- c("Pre-Normalisation", "Post-Normalisation")
  } else {
    ncol_composite <- 3
    file_names <- c(
      "pre_norm_pca.png", "post_norm_pca.png", "ruv_corrected_pca.png",
      "pre_norm_density.png", "post_norm_density.png", "ruv_corrected_density.png",
      "pre_norm_rle.png", "post_norm_rle.png", "ruv_corrected_rle.png",
      "pre_norm_correlation.png", "post_norm_correlation.png", "ruv_corrected_correlation.png",
      # The RUV control plot stands alone on its row: centred (the middle column), not under the RUV column.
      NA, "ruv_corrected_cancor.png", NA
    )
    row_labels <- list(
      pca = c("a)", "b)", "c)"),
      density = c("d)", "e)", "f)"),
      rle = c("g)", "h)", "i)"),
      correlation = c("j)", "k)", "l)"),
      cancor = c("", "m)", "")
    )
    column_labels <- c("Pre-Normalisation", "Post-Normalisation", "RUV-Corrected")
  }

  plot_files <- if (!is.null(qcDir)) {
    vapply(
      file_names,
      function(fn) {
        if (is.na(fn)) NA_character_ else file.path(qcDir, fn)
      },
      character(1)
    )
  } else {
    NULL
  }

  list(
    ncol = ncol_composite,
    plotFiles = plot_files,
    rowLabels = row_labels,
    columnLabels = column_labels
  )
}

generateProtNormCompositeQcFigure <- function(
  ruvMode,
  qcDir,
  omicType,
  resolveInputsFn = resolveProtNormCompositeFigureInputs,
  generateCompositeFn = generateProtNormCompositeFromFiles,
  savePlotFn = savePlot,
  messageFn = message
) {
  messageFn("*** STEP 6B: Generating composite QC figure from saved images ***")

  if (is.null(qcDir)) {
    return(invisible(NULL))
  }

  # ggplot2 retains the last constructed plot independently of our GUI holders.
  # Keep the caller's plot while releasing this export's raster composite on exit.
  previous_plot <- ggplot2::get_last_plot()
  on.exit(ggplot2::set_last_plot(previous_plot), add = TRUE)
  composite_inputs <- resolveInputsFn(ruvMode, qcDir)

  tryCatch({
    composite_res <- generateCompositeFn(
      plotFiles = composite_inputs$plotFiles,
      ncol = composite_inputs$ncol,
      rowLabels = composite_inputs$rowLabels,
      columnLabels = composite_inputs$columnLabels
    )

    if (!is.null(composite_res)) {
      savePlotFn(
        composite_res$plot,
        base_path = qcDir,
        plot_name = paste0(omicType, "_composite_QC_figure"),
        width = composite_res$width,
        height = composite_res$height,
        dpi = 150,
        limitsize = FALSE
      )
      messageFn("*** STEP 6B: Composite QC figure saved to protein_qc_dir ***")
    }
  }, error = function(e) {
    messageFn(paste("Warning: Could not generate composite QC figure:", e$message))
  })

  invisible(NULL)
}

finalizeProtNormWorkflowState <- function(
  normData,
  gcFn = gc,
  messageFn = message
) {
  normData$QC_composite_figure <- NULL
  messageFn("*** STEP 6: Clearing redundant plot objects from memory ***")
  gcFn()

  normData$ruv_complete <- TRUE
  messageFn("*** STEP 6: RUV-corrected workflow completed ***")
  messageFn("*** STEP 7: Enabling correlation filtering step ***")
  messageFn("*** STEP 7: Normalization and RUV workflow completed - ready for correlation filtering ***")

  invisible(NULL)
}

# --- The DIA_LIMPA route's normalisation pass-through (beside runDIAReferenceNormalization's callees) ---

#' The DIA_LIMPA route's normalisation action (KIT-041): the between-sample normalisation and the peptide-count
#' filter ran in the protein QC stage, so the tab's normalisation button designates the current state as the
#' correlation step's input, draws the post-normalisation QC from it and completes; the RUV covariates run after
#' the correlation filter.
#' @noRd
runDIAReferenceLimpaNormalisationPassThrough <- function(settings, workflow, norm_data, post_qc) {
  manager <- workflow$state_manager
  withDIAStateScope(manager, function() {
    assertDIAExportReady(manager, allow_action = "normalization")
    source <- getDIAStateReference(manager)
    norm_data$normalization_input_reference <- source
    norm_data$normalized_protein_obj <- source
    norm_data$ruv_normalized_obj <- source
    norm_data$ruv_optimization_result <- NULL
    norm_data$best_k <- NA
    norm_data$control_genes_index <- NA
    post_qc(source)
    norm_data$normalization_complete <- TRUE
    norm_data$filtering_summary_text <- paste0(
      "DIA with limpa: between-sample normalisation (", diaStateHistory(manager) |> (\(h) if ("protein_normalised" %in% h) "applied in the protein QC stage" else "not run")(),
      ") and the peptide count filter are protein QC steps; RUV-III factors are computed as covariates after the correlation filter.")
    TRUE
  })
}

# --- Correlation result files (#17, #18) ---

resolveProtNormCorrelationResultFilenames <- function(
  normRuvOptimizationResult = NULL,
  workflowRuvOptimizationResult = NULL,
  messagePrefix = "*** CORRELATION",
  messageFn = message
) {
  ruv_was_skipped <- isTRUE(normRuvOptimizationResult$ruv_skipped) ||
    isTRUE(workflowRuvOptimizationResult$ruv_skipped)

  if (ruv_was_skipped) {
    messageFn(sprintf("%s: RUV was skipped - using 'normalised' file names ***", messagePrefix))
    return(list(
      ruvWasSkipped = TRUE,
      tsvFilename = "normalised_results_cln_with_replicates.tsv",
      rdsFilename = "normalised_results_cln_with_replicates.RDS"
    ))
  }

  messageFn(sprintf("%s: RUV was applied - using 'ruv_normalised' file names ***", messagePrefix))

  list(
    ruvWasSkipped = FALSE,
    tsvFilename = "ruv_normalised_results_cln_with_replicates.tsv",
    rdsFilename = "ruv_normalised_results_cln_with_replicates.RDS"
  )
}

saveProtNormCorrelationResults <- function(
  finalS4ForDe,
  experimentPaths,
  normRuvOptimizationResult = NULL,
  workflowRuvOptimizationResult = NULL,
  messagePrefix = "*** CORRELATION STEP 4",
  writeTsvFn = readr::write_tsv,
  saveRdsFn = saveRDS,
  messageFn = message,
  compatibilityExport = !isDIANativeSession(getDIAActiveManager())
) {
  file_names <- resolveProtNormCorrelationResultFilenames(
    normRuvOptimizationResult = normRuvOptimizationResult,
    workflowRuvOptimizationResult = workflowRuvOptimizationResult,
    messagePrefix = messagePrefix,
    messageFn = messageFn
  )

  if (!is.null(experimentPaths) && "protein_qc_dir" %in% names(experimentPaths)) {
    writeTsvFn(
      finalS4ForDe@protein_quant_table,
      file.path(experimentPaths$protein_qc_dir, file_names$tsvFilename)
    )

    if (isTRUE(compatibilityExport)) {
      saveRdsFn(
        finalS4ForDe,
        file.path(experimentPaths$protein_qc_dir, file_names$rdsFilename)
      )
    } else {
      file_names$rdsFilename <- NULL
    }

    messageFn(sprintf(
      "%s: Saved files: %s ***",
      messagePrefix,
      paste(c(file_names$tsvFilename, file_names$rdsFilename), collapse = ", ")
    ))
  }

  invisible(file_names)
}

# --- Session export metadata and summary (#19, #20) ---

buildProtNormExportSummaryContent <- function(
  sessionData,
  sessionFilename,
  timeFn = Sys.time,
  formatTimeFn = format
) {
  sprintf(
    paste(
      "Filtered Session Data Export Summary",
      "=====================================",
      "",
      "Export Timestamp: %s",
      "Session File: %s",
      "",
      "Data Summary:",
      "- Proteins: %d",
      "- Samples: %d",
      "- Contrasts: %d",
      "- Normalization: %s",
      "- RUV Mode: %s",
      "- RUV K: %s",
      "- Correlation Threshold: %s",
      "- Workflow Type: %s",
      "- Report Template: %s",
      "- limpa DPC-Quant: %s",
      "- limpa Quantification Method: %s",
      "",
      "Contrasts:",
      "%s",
      "",
      "This data is ready for differential expression analysis.",
      "Use 'Load Filtered Session' in the DE tab to import.",
      sep = "\n"
    ),
    formatTimeFn(timeFn(), "%Y-%m-%d %H:%M:%S"),
    sessionFilename,
    sessionData$final_protein_count,
    sessionData$final_sample_count,
    ifelse(is.null(sessionData$contrasts_tbl), 0, nrow(sessionData$contrasts_tbl)),
    sessionData$normalization_method,
    sessionData$ruv_mode,
    ifelse(is.null(sessionData$ruv_k), "NA", sessionData$ruv_k),
    ifelse(is.null(sessionData$correlation_threshold), "NA", sessionData$correlation_threshold),
    ifelse(is.null(sessionData$workflow_type), "NA", sessionData$workflow_type),
    ifelse(is.null(sessionData$report_template), "NA", sessionData$report_template),
    ifelse(isTRUE(sessionData$limpa_applied), "yes", "no"),
    ifelse(
      is.null(sessionData$limpa_dpc_quant_results$quantification_method),
      "NA",
      sessionData$limpa_dpc_quant_results$quantification_method
    ),
    if (!is.null(sessionData$contrasts_tbl)) paste(sessionData$contrasts_tbl$friendly_names, collapse = "\n") else "None"
  )
}

saveProtNormExportMetadataFiles <- function(
  sessionData,
  sourceDir,
  saveRdsFn = saveRDS,
  fileExistsFn = file.exists,
  messageFn = message,
  strict = FALSE
) {
  tryCatch({
    if (!is.null(sessionData$accession_cleanup_results)) {
      saveRdsFn(sessionData$accession_cleanup_results, file.path(sourceDir, "accession_cleanup_results.RDS"))
      messageFn("*** EXPORT: Saved accession_cleanup_results.RDS ***")
    }

    if (!is.null(sessionData$qc_params)) {
      saveRdsFn(sessionData$qc_params, file.path(sourceDir, "qc_params.RDS"))
      messageFn("*** EXPORT: Saved qc_params.RDS ***")
    }

    if (!is.null(sessionData$fasta_metadata)) {
      fasta_metadata_file <- file.path(sourceDir, "fasta_metadata.RDS")
      if (!fileExistsFn(fasta_metadata_file)) {
        saveRdsFn(sessionData$fasta_metadata, fasta_metadata_file)
        messageFn("*** EXPORT: Saved fasta_metadata.RDS ***")
      }
    }

    if (!is.null(sessionData$ruv_optimization_result)) {
      ruv_file <- file.path(sourceDir, "ruv_optimization_results.RDS")
      if (!fileExistsFn(ruv_file)) {
        saveRdsFn(sessionData$ruv_optimization_result, ruv_file)
        messageFn("*** EXPORT: Saved ruv_optimization_results.RDS ***")
      }
    }

    if (!is.null(sessionData$mixed_species_analysis)) {
      mixed_species_file <- file.path(sourceDir, "mixed_species_analysis.RDS")
      saveRdsFn(sessionData$mixed_species_analysis, mixed_species_file)
      messageFn(sprintf(
        "*** EXPORT: Saved mixed_species_analysis.RDS (enabled: %s) ***",
        isTRUE(sessionData$mixed_species_analysis$enabled)
      ))
    }

    if (!is.null(sessionData$limpa_parameters)) {
      saveRdsFn(sessionData$limpa_parameters, file.path(sourceDir, "limpa_parameters.RDS"))
      messageFn("*** EXPORT: Saved limpa_parameters.RDS ***")
    }

    if (!is.null(sessionData$limpa_dpc_quant_results)) {
      saveRdsFn(sessionData$limpa_dpc_quant_results, file.path(sourceDir, "limpa_dpc_quant_results.RDS"))
      messageFn("*** EXPORT: Saved limpa_dpc_quant_results.RDS ***")
    }
  }, error = function(e) {
    messageFn(sprintf("*** WARNING: Some metadata files could not be saved: %s ***", e$message))
    if (strict) stop(e)
  })

  invisible(NULL)
}
