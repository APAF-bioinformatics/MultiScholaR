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
  messageFn = message,
  dispatchFiguresFn = NULL
) {
  # KIT-122: on a managed session the worker saves each existing stage's figures under the new aesthetics first, and
  # the regeneration below runs in its delivery (this function again, without the dispatcher).
  if (is.function(dispatchFiguresFn)) {
    post <- isTRUE(normData$normalization_complete) && !is.null(normData$normalized_protein_obj)
    ruv <- isTRUE(normData$ruv_complete) && !is.null(normData$ruv_normalized_obj)
    if (!isTRUE(normData$pre_norm_qc_generated) && !post && !ruv) return(invisible(NULL))
    dispatchFiguresFn(
      sources = function(manager) c(
        if (isTRUE(normData$pre_norm_qc_generated)) list(resolveProtNormQcStateObject(manager, messageFn = messageFn)$currentS4),
        if (isTRUE(normData$normalization_complete)) list(normData$normalized_protein_obj),
        if (isTRUE(normData$ruv_complete)) list(normData$ruv_normalized_obj)),
      then = function() regenerateProtNormQcForAestheticChange(normData, generatePreNormalizationQcFn,
        generatePostNormalizationQcFn, generateRuvCorrectedQcFn, messageFn),
      failed = function(e) messageFn(paste("Error regenerating the QC figures:", conditionMessage(e))))
    return(invisible(NULL))
  }
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

recordProtNormQcPlotPath <- function(qcPlotPaths, stage, plotType, path) {
  qcPlotPaths <- initializeProtNormQcPlotPaths(qcPlotPaths)

  if (!is.null(path)) {
    qcPlotPaths[[stage]][[plotType]] <- path
  }

  qcPlotPaths
}
