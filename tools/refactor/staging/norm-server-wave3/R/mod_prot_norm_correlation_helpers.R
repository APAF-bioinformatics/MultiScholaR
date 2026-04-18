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
  messageFn = message
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

    saveRdsFn(
      finalS4ForDe,
      file.path(experimentPaths$protein_qc_dir, file_names$rdsFilename)
    )

    messageFn(sprintf(
      "%s: Saved files: %s, %s ***",
      messagePrefix,
      file_names$tsvFilename,
      file_names$rdsFilename
    ))
  }

  invisible(file_names)
}

updateProtNormFinalFilteringPlot <- function(
  finalS4ForDe,
  normData,
  omicType,
  experimentLabel,
  updateProteinFilteringFn = updateProteinFiltering,
  messagePrefix = "*** CORRELATION STEP 3",
  messageFn = message
) {
  final_filtering_plot <- tryCatch({
    messageFn("   [mod_prot_norm] Calling updateProteinFiltering...")
    result <- updateProteinFilteringFn(
      data = finalS4ForDe@protein_quant_table,
      step_name = "12_correlation_filtered",
      omic_type = omicType,
      experiment_label = experimentLabel,
      return_grid = TRUE,
      overwrite = TRUE
    )
    messageFn("   [mod_prot_norm] updateProteinFiltering returned.")
    result
  }, error = function(e) {
    messageFn(paste("*** WARNING: updateProteinFiltering failed:", e$message, "***"))
    messageFn(sprintf("%s: Continuing without filtering plot update ***", messagePrefix))
    NULL
  })

  normData$final_filtering_plot <- final_filtering_plot

  invisible(final_filtering_plot)
}

updateProtNormFinalQcPlot <- function(
  finalS4ForDe,
  normData,
  title,
  getPlotAestheticsFn,
  plotPcaFn = plotPca,
  messageFn = message
) {
  tryCatch({
    aesthetics <- getPlotAestheticsFn()
    normData$final_qc_plot <- plotPcaFn(
      finalS4ForDe,
      grouping_variable = aesthetics$color_var,
      label_column = "",
      shape_variable = aesthetics$shape_var,
      title = title,
      font_size = 8
    )
  }, error = function(e) {
    messageFn(paste("Error generating final QC plot:", e$message))
  })

  invisible(normData$final_qc_plot)
}

finalizeProtNormCorrelationWorkflowState <- function(
  finalS4ForDe,
  workflowData,
  correlationThreshold = NULL,
  skipped = FALSE,
  timeFn = Sys.time,
  messagePrefix = "*** CORRELATION",
  messageFn = message,
  catFn = cat
) {
  workflowData$ruv_normalised_for_da_analysis_obj <- finalS4ForDe

  config_object <- if (skipped) {
    list(min_pearson_correlation_threshold = 0, skipped = TRUE)
  } else {
    list(min_pearson_correlation_threshold = correlationThreshold)
  }

  description <- if (skipped) {
    "Skipped final sample correlation filter"
  } else {
    "Applied final sample correlation filter (chunk 28)"
  }

  workflowData$state_manager$saveState(
    state_name = "correlation_filtered",
    s4_data_object = finalS4ForDe,
    config_object = config_object,
    description = description
  )

  trigger_label <- if (skipped) {
    "--- Entering STATE UPDATE TRIGGER setting (SKIP MODE) ---\n"
  } else {
    "--- Entering STATE UPDATE TRIGGER setting ---\n"
  }
  exit_label <- if (skipped) {
    "--- Exiting STATE UPDATE TRIGGER setting (SKIP MODE) ---\n"
  } else {
    "--- Exiting STATE UPDATE TRIGGER setting ---\n"
  }

  catFn(trigger_label)
  catFn("   STATE UPDATE Step: Setting workflow_data$state_update_trigger...\n")
  old_trigger_value <- workflowData$state_update_trigger
  new_trigger_value <- timeFn()
  workflowData$state_update_trigger <- new_trigger_value
  if (!skipped) {
    catFn(sprintf("   STATE UPDATE Step: Old trigger value = %s\n", old_trigger_value))
    catFn(sprintf("   STATE UPDATE Step: New trigger value = %s\n", new_trigger_value))
  }
  catFn("   STATE UPDATE Step: state_update_trigger SET SUCCESSFULLY\n")

  catFn("   STATE UPDATE Step: Setting tab status...\n")
  updated_status <- workflowData$tab_status
  updated_status$normalization <- "complete"
  updated_status$differential_expression <- "pending"
  workflowData$tab_status <- updated_status
  catFn(sprintf("   STATE UPDATE Step: normalization status = %s\n", workflowData$tab_status$normalization))
  if (!skipped) {
    catFn(sprintf("   STATE UPDATE Step: differential_expression status = %s\n", workflowData$tab_status$differential_expression))
  }
  catFn(exit_label)

  final_protein_count <- length(unique(finalS4ForDe@protein_quant_table$Protein.Ids))
  final_sample_count <- length(setdiff(
    colnames(finalS4ForDe@protein_quant_table),
    finalS4ForDe@protein_id_column
  ))

  if (is.null(workflowData$protein_counts)) {
    workflowData$protein_counts <- list()
  }
  workflowData$protein_counts$final_for_de <- final_protein_count

  messageFn(sprintf("%s: Tracked final protein count for DE: %d ***", messagePrefix, final_protein_count))

  list(
    finalProteinCount = final_protein_count,
    finalSampleCount = final_sample_count
  )
}

buildProtNormCorrelationSummaryText <- function(
  finalProteinCount,
  finalSampleCount,
  correlationThreshold = NULL,
  skipped = FALSE
) {
  if (skipped) {
    return(sprintf(
      paste(
        "Correlation filtering SKIPPED.",
        "",
        "All samples retained.",
        "Proteins remaining: %d",
        "Samples remaining: %d",
        "",
        "Ready for differential expression analysis.",
        sep = "\n"
      ),
      finalProteinCount,
      finalSampleCount
    ))
  }

  sprintf(
    paste(
      "Correlation filtering completed successfully!",
      "",
      "Threshold: %.2f",
      "Proteins remaining: %d",
      "Samples remaining: %d",
      "",
      "Ready for differential expression analysis.",
      sep = "\n"
    ),
    correlationThreshold,
    finalProteinCount,
    finalSampleCount
  )
}

resolveProtNormCorrelationInputObject <- function(
  ruvNormalizedObj,
  startMessage,
  missingObjectMessage,
  messageFn = message
) {
  if (is.null(ruvNormalizedObj)) {
    stop(missingObjectMessage)
  }

  messageFn(sprintf("--- DEBUG66 [mod_prot_norm]: %s ---", startMessage))
  messageFn(sprintf(
    "   [mod_prot_norm] Input Object (ruv_s4) Dimensions: %d proteins x %d samples",
    nrow(ruvNormalizedObj@protein_quant_table),
    ncol(ruvNormalizedObj@protein_quant_table) - 1
  ))

  ruvNormalizedObj
}

runProtNormCorrelationVectorStep <- function(
  ruvS4,
  correlationThreshold,
  normData,
  getRuvGroupingVariableFn,
  pearsonCorForSamplePairsFn = pearsonCorForSamplePairs,
  timeFn = Sys.time,
  messageFn = message
) {
  messageFn("*** CORRELATION STEP 1: Calculating sample correlations ***")
  messageFn("   [mod_prot_norm] Calling pearsonCorForSamplePairs...")

  start_time_step1 <- timeFn()
  correlation_vec <- pearsonCorForSamplePairsFn(
    ruvS4,
    tech_rep_remove_regex = "pool",
    correlation_group = getRuvGroupingVariableFn()
  )
  end_time_step1 <- timeFn()

  messageFn(sprintf(
    "   [mod_prot_norm] pearsonCorForSamplePairs returned. Duration: %.2f secs",
    as.numeric(difftime(end_time_step1, start_time_step1, units = "secs"))
  ))
  messageFn(sprintf("   [mod_prot_norm] Correlation Vector Size: %d rows", nrow(correlation_vec)))

  normData$correlation_vector <- correlation_vec
  normData$correlation_threshold <- correlationThreshold
  messageFn("*** CORRELATION STEP 1: Sample correlations calculated ***")

  correlation_vec
}

runProtNormCorrelationFilterStep <- function(
  ruvS4,
  correlationVec,
  correlationThreshold,
  normData,
  filterSamplesFn = filterSamplesByProteinCorrelationThreshold,
  gcFn = gc,
  timeFn = Sys.time,
  messageFn = message
) {
  messageFn("*** CORRELATION STEP 2: Applying correlation threshold filter ***")
  messageFn("   [mod_prot_norm] Calling filterSamplesByProteinCorrelationThreshold...")

  start_time_step2 <- timeFn()
  final_s4_for_de <- filterSamplesFn(
    ruvS4,
    pearson_correlation_per_pair = correlationVec,
    min_pearson_correlation_threshold = correlationThreshold
  )
  end_time_step2 <- timeFn()

  messageFn(sprintf(
    "   [mod_prot_norm] filterSamplesByProteinCorrelationThreshold returned. Duration: %.2f secs",
    as.numeric(difftime(end_time_step2, start_time_step2, units = "secs"))
  ))

  normData$correlation_filtered_obj <- final_s4_for_de
  messageFn("*** CORRELATION STEP 2: Correlation filtering applied ***")
  messageFn("*** CORRELATION STEP 2: Running garbage collection ***")
  gcFn()

  final_s4_for_de
}

prepareProtNormSkippedCorrelationState <- function(
  ruvS4,
  normData,
  gcFn = gc,
  messageFn = message
) {
  messageFn("*** SKIP CORRELATION: Bypassing sample filtering ***")

  normData$correlation_vector <- NULL
  normData$correlation_threshold <- NULL
  normData$correlation_filtered_obj <- ruvS4

  messageFn("*** SKIP CORRELATION: Object passed through without filtering ***")
  gcFn()

  ruvS4
}

runProtNormApplyCorrelationWorkflow <- function(
  ruvS4,
  correlationThreshold,
  normData,
  experimentPaths,
  workflowData,
  omicType,
  experimentLabel,
  getRuvGroupingVariableFn,
  getPlotAestheticsFn,
  runVectorStepFn = runProtNormCorrelationVectorStep,
  runFilterStepFn = runProtNormCorrelationFilterStep,
  updateFinalFilteringPlotFn = updateProtNormFinalFilteringPlot,
  updateFinalQcPlotFn = updateProtNormFinalQcPlot,
  saveCorrelationResultsFn = saveProtNormCorrelationResults,
  withProgressFn = shiny::withProgress,
  incProgressFn = shiny::incProgress,
  gcFn = gc,
  messageFn = message
) {
  withProgressFn(message = "Applying correlation filter...", value = 0, {
    incProgressFn(0.3, detail = "Calculating sample correlations...")
    correlation_vec <- runVectorStepFn(
      ruvS4 = ruvS4,
      correlationThreshold = correlationThreshold,
      normData = normData,
      getRuvGroupingVariableFn = getRuvGroupingVariableFn,
      messageFn = messageFn
    )

    incProgressFn(0.4, detail = "Filtering low-correlation samples...")
    final_s4_for_de <- runFilterStepFn(
      ruvS4 = ruvS4,
      correlationVec = correlation_vec,
      correlationThreshold = correlationThreshold,
      normData = normData,
      messageFn = messageFn
    )

    incProgressFn(0.2, detail = "Updating tracking...")
    messageFn("*** CORRELATION STEP 3: Updating protein filtering tracking ***")

    updateFinalFilteringPlotFn(
      finalS4ForDe = final_s4_for_de,
      normData = normData,
      omicType = omicType,
      experimentLabel = experimentLabel,
      messagePrefix = "*** CORRELATION STEP 3",
      messageFn = messageFn
    )

    updateFinalQcPlotFn(
      finalS4ForDe = final_s4_for_de,
      normData = normData,
      title = "Final Correlation-Filtered Data",
      getPlotAestheticsFn = getPlotAestheticsFn,
      messageFn = messageFn
    )

    incProgressFn(0.1, detail = "Saving results...")
    messageFn("*** CORRELATION STEP 4: Saving final results ***")
    messageFn("*** CORRELATION STEP 4: Running garbage collection before saves ***")
    gcFn()

    tryCatch({
      saveCorrelationResultsFn(
        finalS4ForDe = final_s4_for_de,
        experimentPaths = experimentPaths,
        normRuvOptimizationResult = normData$ruv_optimization_result,
        workflowRuvOptimizationResult = workflowData$ruv_optimization_result,
        messagePrefix = "*** CORRELATION STEP 4",
        messageFn = messageFn
      )
    }, error = function(e) {
      messageFn(paste("Warning: Could not save files:", e$message))
    })

    final_s4_for_de
  })
}

runProtNormSkipCorrelationWorkflow <- function(
  ruvS4,
  normData,
  experimentPaths,
  workflowData,
  omicType,
  experimentLabel,
  getPlotAestheticsFn,
  prepareSkippedStateFn = prepareProtNormSkippedCorrelationState,
  updateFinalFilteringPlotFn = updateProtNormFinalFilteringPlot,
  updateFinalQcPlotFn = updateProtNormFinalQcPlot,
  saveCorrelationResultsFn = saveProtNormCorrelationResults,
  withProgressFn = shiny::withProgress,
  incProgressFn = shiny::incProgress,
  gcFn = gc,
  messageFn = message
) {
  withProgressFn(message = "Skipping filter and saving...", value = 0, {
    incProgressFn(0.3, detail = "Bypassing filter...")
    final_s4_for_de <- prepareSkippedStateFn(
      ruvS4 = ruvS4,
      normData = normData,
      gcFn = gcFn,
      messageFn = messageFn
    )

    incProgressFn(0.2, detail = "Updating tracking...")
    messageFn("*** SKIP CORRELATION: Updating protein filtering tracking ***")

    updateFinalFilteringPlotFn(
      finalS4ForDe = final_s4_for_de,
      normData = normData,
      omicType = omicType,
      experimentLabel = experimentLabel,
      messagePrefix = "*** SKIP CORRELATION",
      messageFn = messageFn
    )

    updateFinalQcPlotFn(
      finalS4ForDe = final_s4_for_de,
      normData = normData,
      title = "Final Data (Correlation Filter Skipped)",
      getPlotAestheticsFn = getPlotAestheticsFn,
      messageFn = messageFn
    )

    incProgressFn(0.2, detail = "Saving results...")
    messageFn("*** SKIP CORRELATION: Saving final results ***")
    gcFn()

    tryCatch({
      saveCorrelationResultsFn(
        finalS4ForDe = final_s4_for_de,
        experimentPaths = experimentPaths,
        normRuvOptimizationResult = normData$ruv_optimization_result,
        workflowRuvOptimizationResult = workflowData$ruv_optimization_result,
        messagePrefix = "*** SKIP CORRELATION",
        messageFn = messageFn
      )
    }, error = function(e) {
      messageFn(paste("Warning: Could not save files:", e$message))
    })

    final_s4_for_de
  })
}

completeProtNormCorrelationWorkflow <- function(
  finalS4ForDe,
  workflowData,
  output,
  normData,
  correlationThreshold = NULL,
  skipped = FALSE,
  successNotification,
  completionMessage,
  messagePrefix,
  showNotificationFn = shiny::showNotification,
  renderTextFn = shiny::renderText,
  messageFn = message
) {
  correlation_metrics <- finalizeProtNormCorrelationWorkflowState(
    finalS4ForDe = finalS4ForDe,
    workflowData = workflowData,
    correlationThreshold = correlationThreshold,
    skipped = skipped,
    messagePrefix = messagePrefix,
    messageFn = messageFn
  )

  correlation_summary <- buildProtNormCorrelationSummaryText(
    finalProteinCount = correlation_metrics$finalProteinCount,
    finalSampleCount = correlation_metrics$finalSampleCount,
    correlationThreshold = correlationThreshold,
    skipped = skipped
  )
  output$correlation_filter_summary <- renderTextFn(correlation_summary)

  normData$correlation_filtering_complete <- TRUE

  showNotificationFn(
    successNotification,
    type = "message",
    duration = 5
  )

  messageFn(completionMessage)

  invisible(correlation_metrics)
}

handleProtNormCorrelationError <- function(
  error,
  logPrefix,
  notificationPrefix,
  showNotificationFn = shiny::showNotification,
  messageFn = message
) {
  messageFn(paste(logPrefix, error$message))
  showNotificationFn(
    paste(notificationPrefix, error$message),
    type = "error",
    duration = 10
  )

  invisible(NULL)
}

