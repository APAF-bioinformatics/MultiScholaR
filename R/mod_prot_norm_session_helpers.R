canProtNormExportFilteredSession <- function(
  correlationFilteringComplete,
  correlationFilteredObj
) {
  isTRUE(correlationFilteringComplete) && !is.null(correlationFilteredObj)
}

resolveProtNormExportSourceDir <- function(
  experimentPaths,
  dirExistsFn = dir.exists
) {
  source_dir <- experimentPaths$source_dir
  if (is.null(source_dir) || !dirExistsFn(source_dir)) {
    stop("Could not find the source directory to save session data.")
  }

  source_dir
}

collectProtNormExportSessionData <- function(
  workflowData,
  normData,
  input,
  timeFn = Sys.time,
  messageFn = message
) {
  stateSnapshot <- workflowStateLegacySnapshot(workflowData$state_manager)
  current_state_name <- stateSnapshot$r6_current_state_name
  current_s4_object <- workflowData$state_manager$getState(current_state_name)
  r6_complete_states <- stateSnapshot$r6_complete_states
  if (is.null(r6_complete_states)) {
    r6_complete_states <- list()
  }
  r6_state_history <- stateSnapshot$r6_state_history
  if (is.null(r6_state_history)) {
    r6_state_history <- current_state_name
  }
  stateManifest <- workflowStateManifest(workflowData$state_manager)

  workflow_type <- if (!is.null(workflowData$config_list) &&
    !is.null(workflowData$config_list$globalParameters) &&
    !is.null(workflowData$config_list$globalParameters$workflow_type)) {
    workflowData$config_list$globalParameters$workflow_type
  } else if (!is.null(current_s4_object@args$globalParameters$workflow_type)) {
    current_s4_object@args$globalParameters$workflow_type
  } else {
    "DIA"
  }

  current_args <- tryCatch(current_s4_object@args, error = function(e) list())
  current_global_params <- if (!is.null(current_args$globalParameters)) {
    current_args$globalParameters
  } else {
    list()
  }
  config_global_params <- if (!is.null(workflowData$config_list) &&
    !is.null(workflowData$config_list$globalParameters)) {
    workflowData$config_list$globalParameters
  } else {
    list()
  }
  report_template <- config_global_params$report_template
  if (is.null(report_template) || !nzchar(report_template)) {
    report_template <- current_global_params$report_template
  }
  limpa_dpc_quant_results <- current_args$limpa_dpc_quant_results
  limpa_parameters <- current_args$proteinMissingValueImputationLimpa
  limpa_applied <- isTRUE(config_global_params$use_limpa) ||
    isTRUE(current_global_params$use_limpa) ||
    !is.null(limpa_dpc_quant_results) ||
    !is.null(limpa_parameters)

  session_data <- list(
    r6_current_state_name = current_state_name,
    r6_complete_states = r6_complete_states,
    r6_state_history = r6_state_history,
    current_s4_object = current_s4_object,
    correlation_filtered_s4 = normData$correlation_filtered_obj,
    contrasts_tbl = workflowData$contrasts_tbl,
    design_matrix = workflowData$design_matrix,
    config_list = workflowData$config_list,
    export_timestamp = timeFn(),
    normalization_method = input$norm_method,
    ruv_mode = input$ruv_mode,
    ruv_applied = (input$ruv_mode != "skip"),
    ruv_k = normData$best_k,
    correlation_threshold = normData$correlation_threshold,
    workflow_type = workflow_type,
    report_template = report_template,
    use_limpa = limpa_applied,
    limpa_applied = limpa_applied,
    limpa_parameters = limpa_parameters,
    limpa_dpc_quant_results = limpa_dpc_quant_results,
    fasta_metadata = workflowData$fasta_metadata,
    accession_cleanup_results = workflowData$accession_cleanup_results,
    ruv_optimization_result = workflowData$ruv_optimization_result,
    qc_params = workflowData$qc_params,
    protein_counts = workflowData$protein_counts,
    mixed_species_analysis = workflowData$mixed_species_analysis,
    final_protein_count = countDistinctProteinQuantIdentities(current_s4_object),
    final_sample_count = length(setdiff(
      colnames(current_s4_object@protein_quant_table),
      current_s4_object@protein_id_column
    ))
  )
  if (!is.null(stateManifest)) {
    session_data$workflow_state_manifest <- stateManifest
  }

  messageFn("*** EXPORT: Gathered session data successfully ***")
  messageFn(sprintf("*** EXPORT: Final protein count: %d ***", session_data$final_protein_count))
  messageFn(sprintf("*** EXPORT: Final sample count: %d ***", session_data$final_sample_count))
  messageFn(sprintf(
    "*** EXPORT: Contrasts available: %d ***",
    ifelse(is.null(session_data$contrasts_tbl), 0, nrow(session_data$contrasts_tbl))
  ))

  session_data
}

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
  messageFn = message
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
  })

  invisible(NULL)
}

saveProtNormExportArtifacts <- function(
  sessionData,
  sourceDir,
  timeFn = Sys.time,
  formatTimeFn = format,
  saveRdsFn = saveRDS,
  writeLinesFn = writeLines,
  fileExistsFn = file.exists,
  messageFn = message
) {
  timestamp_str <- formatTimeFn(timeFn(), "%Y%m%d_%H%M%S")
  session_filename <- sprintf("filtered_session_data_%s.rds", timestamp_str)
  session_filepath <- file.path(sourceDir, session_filename)

  saveRdsFn(sessionData, session_filepath)
  messageFn(sprintf("*** EXPORT: Session data saved to: %s ***", session_filepath))

  latest_filename <- "filtered_session_data_latest.rds"
  latest_filepath <- file.path(sourceDir, latest_filename)
  saveRdsFn(sessionData, latest_filepath)
  messageFn(sprintf("*** EXPORT: Latest version saved to: %s ***", latest_filepath))

  saveProtNormExportMetadataFiles(
    sessionData = sessionData,
    sourceDir = sourceDir,
    saveRdsFn = saveRdsFn,
    fileExistsFn = fileExistsFn,
    messageFn = messageFn
  )

  summary_content <- buildProtNormExportSummaryContent(
    sessionData = sessionData,
    sessionFilename = session_filename,
    timeFn = timeFn,
    formatTimeFn = formatTimeFn
  )
  summary_filepath <- file.path(sourceDir, "filtered_session_summary.txt")
  writeLinesFn(summary_content, summary_filepath)

  list(
    sessionFilename = session_filename,
    sessionFilepath = session_filepath,
    latestFilepath = latest_filepath,
    summaryFilepath = summary_filepath
  )
}

runProtNormExportSessionWorkflow <- function(
  workflowData,
  normData,
  input,
  sourceDir,
  withProgressFn = shiny::withProgress,
  incProgressFn = shiny::incProgress,
  collectSessionDataFn = collectProtNormExportSessionData,
  saveExportArtifactsFn = saveProtNormExportArtifacts,
  messageFn = message
) {
  withProgressFn(message = "Exporting filtered session data...", value = 0, {
    incProgressFn(0.2, detail = "Gathering data...")
    session_data <- collectSessionDataFn(
      workflowData = workflowData,
      normData = normData,
      input = input,
      messageFn = messageFn
    )

    incProgressFn(0.4, detail = "Saving to file...")
    export_artifacts <- saveExportArtifactsFn(
      sessionData = session_data,
      sourceDir = sourceDir,
      messageFn = messageFn
    )

    incProgressFn(0.2, detail = "Creating latest version...")
    incProgressFn(0.1, detail = "Saving metadata files...")
    incProgressFn(0.1, detail = "Creating summary...")

    list(sessionData = session_data, exportArtifacts = export_artifacts)
  })
}

handleProtNormExportError <- function(
  error,
  showNotificationFn = shiny::showNotification,
  messageFn = message
) {
  messageFn(paste("*** ERROR in session export:", error$message, "***"))
  showNotificationFn(
    paste("Error exporting session:", error$message),
    type = "error",
    duration = 10
  )

  invisible(NULL)
}

resolveProtNormPreNormalizationState <- function(
  history,
  preNormalizationStates = c(
    "protein_replicate_filtered",
    "imputed",
    "replicate_filtered",
    "sample_filtered",
    "protein_peptide_filtered",
    "intensity_filtered",
    "precursor_rollup",
    "qvalue_filtered",
    "raw_data_s4",
    "protein_s4_initial"
  )
) {
  matching_states <- intersect(rev(history), preNormalizationStates)

  if (length(matching_states) == 0) {
    return(NULL)
  }

  matching_states[[1]]
}

revertProtNormStateManagerToPreNormalization <- function(
  workflowData,
  preNormalizationStates = c(
    "protein_replicate_filtered",
    "imputed",
    "replicate_filtered",
    "sample_filtered",
    "protein_peptide_filtered",
    "intensity_filtered",
    "precursor_rollup",
    "qvalue_filtered",
    "raw_data_s4",
    "protein_s4_initial"
  ),
  messageFn = message
) {
  previous_state <- NULL
  reverted_s4 <- NULL

  if (!is.null(workflowData$state_manager)) {
    history <- workflowData$state_manager$getHistory()
    previous_state <- resolveProtNormPreNormalizationState(
      history = history,
      preNormalizationStates = preNormalizationStates
    )

    if (!is.null(previous_state)) {
      reverted_s4 <- workflowData$state_manager$revertToState(previous_state)
      messageFn(sprintf(
        "*** RESET: Reverted R6 state manager to '%s' (actual previous state) ***",
        previous_state
      ))
    } else {
      messageFn("*** WARNING: No valid pre-normalization state found in history ***")
    }
  } else {
    messageFn("*** WARNING: workflow_data$state_manager is NULL - cannot revert state ***")
  }

  list(previousState = previous_state, revertedS4 = reverted_s4)
}

resetProtNormReactiveState <- function(
  normData,
  workflowData
) {
  normData$normalization_complete <- FALSE
  normData$ruv_complete <- FALSE
  normData$correlation_filtering_complete <- FALSE
  normData$normalized_protein_obj <- NULL
  normData$ruv_normalized_obj <- NULL
  normData$correlation_filtered_obj <- NULL
  normData$best_k <- NULL
  normData$control_genes_index <- NULL
  normData$correlation_vector <- NULL
  normData$correlation_threshold <- NULL
  normData$final_qc_plot <- NULL
  normData$final_filtering_plot <- NULL
  normData$post_norm_filtering_plot <- NULL
  normData$filtering_summary_text <- NULL
  normData$ruv_optimization_result <- NULL

  if (is.null(normData$qc_plots)) {
    normData$qc_plots <- list()
  }

  normData$qc_plots$post_normalization <- list(
    pca = NULL,
    density = NULL,
    rle = NULL,
    correlation = NULL
  )
  normData$qc_plots$ruv_corrected <- list(
    pca = NULL,
    density = NULL,
    rle = NULL,
    correlation = NULL
  )

  workflowData$ruv_normalised_for_da_analysis_obj <- NULL
  updated_status <- workflowData$tab_status
  updated_status$normalization <- "pending"
  updated_status$differential_expression <- "disabled"
  workflowData$tab_status <- updated_status

  invisible(NULL)
}

resetProtNormOutputs <- function(
  output,
  ruvMode,
  groupingVariable,
  renderTextFn = shiny::renderText
) {
  output$correlation_filter_summary <- renderTextFn(
    getProtNormDefaultCorrelationFilterSummaryText()
  )
  output$filtering_summary_text <- renderTextFn(
    getProtNormFilteringSummaryText(NULL)
  )
  output$ruv_optimization_summary <- renderTextFn(
    buildProtNormRuvOptimizationSummary(
      ruvOptimizationResult = NULL,
      ruvMode = ruvMode,
      groupingVariable = groupingVariable
    )
  )

  invisible(NULL)
}

buildProtNormResetNotificationMessage <- function(
  previousState
) {
  if (!is.null(previousState)) {
    return(sprintf(
      "Normalization has been reset to pre-normalization state (%s)",
      previousState
    ))
  }

  "Normalization has been reset (no valid previous state found)"
}

runProtNormResetWorkflow <- function(
  workflowData,
  normData,
  output,
  ruvMode,
  groupingVariable,
  revertStateManagerFn = revertProtNormStateManagerToPreNormalization,
  resetReactiveStateFn = resetProtNormReactiveState,
  resetOutputsFn = resetProtNormOutputs,
  buildNotificationMessageFn = buildProtNormResetNotificationMessage,
  showNotificationFn = shiny::showNotification,
  messageFn = message
) {
  reset_result <- revertStateManagerFn(
    workflowData = workflowData,
    messageFn = messageFn
  )

  resetReactiveStateFn(
    normData = normData,
    workflowData = workflowData
  )

  resetOutputsFn(
    output = output,
    ruvMode = ruvMode,
    groupingVariable = groupingVariable
  )

  showNotificationFn(
    buildNotificationMessageFn(reset_result$previousState),
    type = "warning",
    duration = 5
  )

  messageFn("*** RESET: Normalization reset completed successfully ***")

  invisible(reset_result)
}

handleProtNormResetError <- function(
  error,
  showNotificationFn = shiny::showNotification,
  messageFn = message
) {
  messageFn(paste("*** ERROR in normalization reset:", error$message, "***"))
  showNotificationFn(
    paste("Error resetting normalization:", error$message),
    type = "error",
    duration = 10
  )

  invisible(NULL)
}
