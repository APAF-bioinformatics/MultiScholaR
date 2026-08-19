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
  observeEventFn = shiny::observeEvent
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
      checkMemoryUsageFn = checkMemoryUsageFn
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
      checkMemoryUsageFn = checkMemoryUsageFn
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
    messageFn(sprintf(
      "Current state at module start: %s",
      workflowStateCurrentName(workflowData$state_manager)
    ))
  }

  normData <- createReactiveStateFn()

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
      qcDir = experimentPaths$protein_qc_dir,
      aesthetics = getPlotAesthetics(),
      qcPlotPaths = normData$qc_plot_paths
    )

    normData$plot_refresh_trigger <- normData$plot_refresh_trigger + 1
  }

  generateRuvCorrectedQc <- function(ruvCorrectedS4) {
    normData$qc_plot_paths <- generateRuvCorrectedQcArtifactsFn(
      ruvCorrectedS4 = ruvCorrectedS4,
      qcDir = experimentPaths$protein_qc_dir,
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
