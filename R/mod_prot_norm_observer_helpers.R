registerProtNormServerObservers <- function(
  input,
  output,
  session,
  selectedTab = NULL,
  workflowData,
  normData,
  experimentPaths,
  omicType,
  experimentLabel,
  generatePreNormalizationQcFn,
  generatePostNormalizationQcFn,
  generateRuvCorrectedQcFn,
  getPlotAestheticsFn,
  getRuvGroupingVariableFn,
  checkMemoryUsageFn = checkProtNormMemoryUsage,
  updateDesignDrivenChoicesFn = updateProtNormDesignDrivenChoices,
  runTabEntryWorkflowFn = runProtNormTabEntryWorkflow,
  regenerateQcForAestheticChangeFn = regenerateProtNormQcForAestheticChange,
  runNormalizationWorkflowFn = runProtNormNormalizationWorkflow,
  runApplyCorrelationObserverFn = runProtNormApplyCorrelationObserver,
  runSkipCorrelationObserverFn = runProtNormSkipCorrelationObserver,
  runResetObserverFn = runProtNormResetObserver,
  runExportObserverFn = runProtNormExportObserver,
  observeFn = shiny::observe,
  observeEventFn = shiny::observeEvent,
  updateRadioButtonsFn = shiny::updateRadioButtons
) {
  observeFn({
    updateDesignDrivenChoicesFn(session, workflowData$design_matrix)
  })

  if (!is.null(selectedTab)) {
    observeEventFn(selectedTab(), {
      runTabEntryWorkflowFn(
        selectedTab = selectedTab(),
        workflowData = workflowData,
        normData = normData,
        generatePreNormalizationQcFn = generatePreNormalizationQcFn
      )
    }, ignoreInit = TRUE)
  }

  observeEventFn(c(input$color_variable, input$shape_variable), {
    regenerateQcForAestheticChangeFn(
      normData = normData,
      generatePreNormalizationQcFn = generatePreNormalizationQcFn,
      generatePostNormalizationQcFn = generatePostNormalizationQcFn,
      generateRuvCorrectedQcFn = generateRuvCorrectedQcFn
    )
  })

  observeEventFn(input$run_normalization, {
    runNormalizationWorkflowFn(
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
  })

  # KIT-041: on the DIA_LIMPA route the RUV-III factors are covariates computed after the correlation filter (the
  # limpa script's order); the tab's RUV mode offers `covariates` or `skip`, set here once the type is known.
  observeFn({
    workflowData$state_update_trigger
    manager <- workflowData$state_manager
    if (identical(tryCatch(manager$workflow_type, error = function(e) NULL), "DIA_LIMPA") &&
        !identical(shiny::isolate(input$ruv_mode), "covariates") && !identical(shiny::isolate(input$ruv_mode), "skip")) {
      updateRadioButtonsFn(session, "ruv_mode", choices = protNormRuvModeChoices("DIA_LIMPA"), selected = "covariates")
    }
  })

  observeEventFn(input$apply_correlation_filter, {
    runApplyCorrelationObserverFn(
      input = input,
      output = output,
      workflowData = workflowData,
      normData = normData,
      experimentPaths = experimentPaths,
      omicType = omicType,
      experimentLabel = experimentLabel,
      getRuvGroupingVariableFn = getRuvGroupingVariableFn,
      getPlotAestheticsFn = getPlotAestheticsFn,
      checkMemoryUsageFn = checkMemoryUsageFn,
      generateRuvCorrectedQcFn = generateRuvCorrectedQcFn
    )
  })

  observeEventFn(input$skip_correlation_filter, {
    runSkipCorrelationObserverFn(
      output = output,
      workflowData = workflowData,
      normData = normData,
      experimentPaths = experimentPaths,
      omicType = omicType,
      experimentLabel = experimentLabel,
      getPlotAestheticsFn = getPlotAestheticsFn,
      checkMemoryUsageFn = checkMemoryUsageFn,
      input = input,
      getRuvGroupingVariableFn = getRuvGroupingVariableFn,
      generateRuvCorrectedQcFn = generateRuvCorrectedQcFn
    )
  })

  observeEventFn(input$reset_normalization, {
    runResetObserverFn(
      workflowData = workflowData,
      normData = normData,
      output = output,
      ruvMode = input$ruv_mode,
      groupingVariable = getRuvGroupingVariableFn()
    )
  })

  observeEventFn(input$export_filtered_session, {
    runExportObserverFn(
      input = input,
      workflowData = workflowData,
      normData = normData,
      experimentPaths = experimentPaths
    )
  })

  invisible(TRUE)
}

runProtNormModuleServerShell <- function(
  input,
  output,
  session,
  id,
  workflowData,
  experimentPaths,
  omicType,
  experimentLabel,
  selectedTab = NULL,
  messageFn = message,
  createReactiveStateFn = createProtNormReactiveState,
  getPlotAestheticsFn = getProtNormPlotAesthetics,
  getRuvGroupingVariableFn = getProtNormRuvGroupingVariable,
  generatePreNormalizationQcArtifactsFn = generateProtNormPreNormalizationQcArtifacts,
  generatePostNormalizationQcArtifactsFn = generateProtNormPostNormalizationQcArtifacts,
  generateRuvCorrectedQcArtifactsFn = generateProtNormRuvCorrectedQcArtifacts,
  registerServerObserversFn = registerProtNormServerObservers,
  registerQcImageOutputsFn = registerProtNormQcImageOutputs,
  registerRenderOutputsFn = registerProtNormRenderOutputs,
  checkMemoryUsageFn = checkProtNormMemoryUsage
) {
  messageFn("=== NORMALIZATION MODULE SERVER STARTED ===")
  messageFn(sprintf("Module ID: %s", id))
  messageFn(sprintf("workflow_data is NULL: %s", is.null(workflowData)))
  if (!is.null(workflowData$state_manager)) {
    messageFn(sprintf("Current state at module start: %s", workflowData$state_manager$current_state))
  }

  normData <- createReactiveStateFn()
  output$export_s4_compatibility <- shiny::downloadHandler(
    filename = function() paste0(omicType, "_", experimentLabel, "_final_s4.RDS"),
    content = function(file) {
      names <- resolveProtNormCorrelationResultFilenames(normData$ruv_optimization_result,
        workflowData$ruv_optimization_result)
      canonical <- file.path(experimentPaths$protein_qc_dir, names$rdsFilename)
      exportDIAS4Compatibility(workflowData$state_manager, canonical)
      if (!file.copy(canonical, file, overwrite = TRUE)) stop("Could not deliver compatibility export")
      if (!is.null(experimentPaths$integration_dir)) {
        destination <- file.path(experimentPaths$integration_dir,
          paste0(omicType, "_", experimentLabel, "_final_s4.RDS"))
        if (!file.copy(canonical, destination, overwrite = TRUE)) stop("Could not publish compatibility Integration file")
      }
    }
  )


  shiny::observeEvent(workflowData$native_session_restored, {
    manager <- workflowData$state_manager
    if (!usesDIAReferenceStates(manager)) return(invisible(NULL))
    withDIAStateScope(manager, function() {
      context <- manager$.context
      extra <- context$native_extra_references %||% list()
      metadata <- if (!is.null(context$native_metadata_reference))
        resolveDIAStateValue(context, context$native_metadata_reference) else list()
      normData$normalization_input_reference <- extra$normalization_input_reference %||%
        getDIAStateReference(manager, "protein_replicate_filtered")
      normData$normalized_protein_obj <- extra$normalized_protein_obj %||% getDIAStateReference(manager, "normalized")
      normData$ruv_normalized_obj <- extra$ruv_normalized_obj %||%
        getDIAStateReference(manager, "ruv_corrected") %||% getDIAStateReference(manager, "normalized")
      normData$correlation_filtered_obj <- getDIAStateReference(manager, "correlation_filtered")
      active_history <- manager$getHistory()
      normalized_active <- any(c("normalized", "ruv_corrected", "correlation_filtered") %in% active_history)
      normData$normalization_complete <- normalized_active && !is.null(normData$normalized_protein_obj)
      normData$ruv_complete <- normalized_active && !is.null(normData$ruv_normalized_obj)
      normData$correlation_filtering_complete <- "correlation_filtered" %in% active_history &&
        !is.null(normData$correlation_filtered_obj)
      if (!normalized_active) normData$normalized_protein_obj <- normData$ruv_normalized_obj <- NULL
      if (!"correlation_filtered" %in% active_history) normData$correlation_filtered_obj <- NULL
      ruv <- metadata$ruv_optimization_result
      if (!is.null(ruv$best_cancor_plot) && !inherits(ruv$best_cancor_plot, "grob"))
        ruv$best_cancor_plot <- ggplot2::ggplotGrob(stripProtNormQcPlotOwners(ruv$best_cancor_plot))
      normData$ruv_optimization_result <- ruv
      normData$best_k <- metadata$ruv_k
      normData$control_genes_index <- ruv$best_control_genes_index
      normData$correlation_threshold <- metadata$correlation_threshold
      normData$pre_norm_qc_generated <- FALSE
      normData$plot_refresh_trigger <- normData$plot_refresh_trigger + 1L
      if (!is.null(metadata$normalization_method)) shiny::updateSelectInput(session, "norm_method",
        selected = metadata$normalization_method)
      invisible(NULL)
    })
  }, ignoreInit = TRUE)

  getPlotAesthetics <- function() {
    getPlotAestheticsFn(input$color_variable, input$shape_variable)
  }

  getRuvGroupingVariable <- function() {
    getRuvGroupingVariableFn(input$ruv_grouping_variable)
  }

  generatePreNormalizationQc <- function() {
    normData$qc_plot_paths <- generatePreNormalizationQcArtifactsFn(
      stateManager = workflowData$state_manager,
      qcDir = experimentPaths$protein_qc_dir,
      aesthetics = getPlotAesthetics(),
      qcPlotPaths = normData$qc_plot_paths
    )
  }

  generatePostNormalizationQc <- function(normalizedS4) {
    normData$qc_plot_paths <- generatePostNormalizationQcArtifactsFn(
      normalizedS4 = normalizedS4,
      qcDir = normData$publication_paths$protein_qc_dir %||% experimentPaths$protein_qc_dir,
      aesthetics = getPlotAesthetics(),
      qcPlotPaths = normData$qc_plot_paths
    )

    normData$plot_refresh_trigger <- normData$plot_refresh_trigger + 1
  }

  generateRuvCorrectedQc <- function(ruvCorrectedS4) {
    normData$qc_plot_paths <- generateRuvCorrectedQcArtifactsFn(
      ruvCorrectedS4 = ruvCorrectedS4,
      qcDir = normData$publication_paths$protein_qc_dir %||% experimentPaths$protein_qc_dir,
      aesthetics = getPlotAesthetics(),
      qcPlotPaths = normData$qc_plot_paths
    )

    normData$plot_refresh_trigger <- normData$plot_refresh_trigger + 1
  }

  registerServerObserversFn(
    input = input,
    output = output,
    session = session,
    selectedTab = selectedTab,
    workflowData = workflowData,
    normData = normData,
    experimentPaths = experimentPaths,
    omicType = omicType,
    experimentLabel = experimentLabel,
    generatePreNormalizationQcFn = generatePreNormalizationQc,
    generatePostNormalizationQcFn = generatePostNormalizationQc,
    generateRuvCorrectedQcFn = generateRuvCorrectedQc,
    getPlotAestheticsFn = getPlotAesthetics,
    getRuvGroupingVariableFn = getRuvGroupingVariable,
    checkMemoryUsageFn = checkMemoryUsageFn
  )

  registerQcImageOutputsFn(
    output = output,
    normData = normData,
    proteinQcDir = experimentPaths$protein_qc_dir
  )

  registerRenderOutputsFn(
    output = output,
    normData = normData,
    ruvMode = input$ruv_mode,
    groupingVariable = getRuvGroupingVariable()
  )

  invisible(normData)
}

runProtNormModuleServerEntryShell <- function(
  id,
  workflowData,
  experimentPaths,
  omicType,
  experimentLabel,
  selectedTab = NULL,
  moduleServerFn = shiny::moduleServer,
  runModuleServerShellFn = runProtNormModuleServerShell
) {
  moduleServerFn(id, function(input, output, session) {
    runModuleServerShellFn(
      input = input,
      output = output,
      session = session,
      id = id,
      workflowData = workflowData,
      experimentPaths = experimentPaths,
      omicType = omicType,
      experimentLabel = experimentLabel,
      selectedTab = selectedTab
    )
  })
}

runProtNormModuleServerPublicWrapper <- function(
  id,
  workflow_data,
  experiment_paths,
  omic_type,
  experiment_label,
  selected_tab = NULL,
  runModuleServerEntryShellFn = runProtNormModuleServerEntryShell
) {
  runModuleServerEntryShellFn(
    id = id,
    workflowData = workflow_data,
    experimentPaths = experiment_paths,
    omicType = omic_type,
    experimentLabel = experiment_label,
    selectedTab = selected_tab
  )
}

