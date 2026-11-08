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

  if (ruvMode %in% c("automatic", "covariates", "correction")) { # KIT-041/056: the DIA_LIMPA sweep shares the summary
    return(sprintf(
      paste(
        paste0(if (ruvMode == "covariates") "RUV-III Covariates (negative-control sweep):" else if (ruvMode == "correction") "RUV-III Correction (negative-control sweep):" else "Automatic RUV Optimization Results:", "\n\n* Best percentage: %.1f%%\n* Best k value: %d\n* Separation score: %.4f\n* Composite score: %.4f\n* Control genes: %d\n* RUV grouping variable: %s\n* Separation metric: %s\n* K penalty weight: %.1f\n* Adaptive penalty: %s\n* Sample size: %s"),
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
  withDIAStateScope(getDIAActiveManager(), function() {
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
      generatePostNormalizationQcFn(getDIAOwnedInput(normData, "normalized_protein_obj"))
    }, error = function(e) {
      messageFn(paste("Error regenerating post-normalization QC:", e$message))
    })
  }

  if (normData$ruv_complete && !is.null(normData$ruv_normalized_obj)) {
    messageFn("Regenerating RUV-corrected QC plots with new aesthetics...")
    tryCatch({
      generateRuvCorrectedQcFn(getDIAOwnedInput(normData, "ruv_normalized_obj"))
    }, error = function(e) {
      messageFn(paste("Error regenerating RUV-corrected QC:", e$message))
    })
  }

  invisible(NULL)
  })
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
  manager <- getDIAActiveManager()
  product_reference <- attr(plotObject, "dia_product_reference", exact = TRUE)
  if (isSavedWorkflowReference(product_reference)) {
    dimensions <- as.integer(c(width, height) * dpi)
    aesthetics <- attr(plotObject, "dia_render_aesthetics", exact = TRUE)
    aesthetics <- aesthetics[!vapply(aesthetics, is.null, logical(1))]
    identity <- newSavedWorkflowRenderIdentity(product_reference,
      c(aesthetics, list(dpi = dpi)),
      "png", dimensions, manager$.context$software)
    hit <- resolveSavedWorkflowAsset(manager$.context$store, product_reference, identity)
    if (identical(hit$status, "hit")) {
      writeBin(hit$value$bytes, artifactPath)
      return(artifactPath)
    }
    savePlotFn(artifactPath, plotObject, width = width, height = height, dpi = dpi)
    asset <- newSavedWorkflowRenderedAsset(readBin(artifactPath, "raw", n = file.size(artifactPath)),
      "image/png", dimensions)
    reference <- saveWorkflowDiaQcAsset(manager$.context$store, asset,
      basename(tempfile("dia-qc-image-")), identity)
    rootDIAProduct(manager, reference, list(product_reference))
    return(artifactPath)
  }
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
