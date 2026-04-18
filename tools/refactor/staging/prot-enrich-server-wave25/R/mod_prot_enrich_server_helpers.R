runProtEnrichAnalysisBody <- function(input,
                                      enrichmentData,
                                      workflowData,
                                      experimentPaths,
                                      currentAnalysisMethodFn,
                                      prepareAnalysisSetupFn = prepareProtEnrichAnalysisBodySetup,
                                      captureCheckpointFn = .capture_checkpoint,
                                      resolveSelectedDaResultsFn = resolveProtEnrichSelectedDaResults,
                                      resolveRunDependenciesFn = resolveProtEnrichRunDependencies,
                                      resolveOutputDirectoriesFn = resolveProtEnrichOutputDirectories,
                                      createDaResultsForEnrichmentFn = createDAResultsForEnrichment,
                                      resolveUniprotAnnotationsFn = resolveProtEnrichUniprotAnnotations,
                                      resolveAnnotationMatchingFn = resolveProtEnrichAnnotationMatching,
                                      resolveOrganismMappingFn = resolveProtEnrichOrganismMapping,
                                      applyOrganismFilterFn = applyProtEnrichOrganismFilter,
                                      persistOrganismFilterMetadataFn = persistProtEnrichOrganismFilterMetadata,
                                      resolveAnalysisInputColumnsFn = resolveProtEnrichAnalysisInputColumns,
                                      buildProcessEnrichmentsArgsFn = buildProtEnrichProcessEnrichmentsArgs,
                                      prepareProcessExecutionFn = prepareProtEnrichProcessExecution,
                                      processEnrichmentsFn = processEnrichments,
                                      executeProcessEnrichmentsFn = executeProtEnrichProcessEnrichments,
                                      buildAllContrastResultsFn = buildProtEnrichAllContrastResults,
                                      resolveSelectedContrastResultsFn = resolveProtEnrichSelectedContrastResults,
                                      capturePostProcessResultsFn = captureProtEnrichPostProcessResults,
                                      persistAnalysisResultsFn = persistProtEnrichAnalysisResults,
                                      buildAnalysisResultsPayloadFn = buildProtEnrichAnalysisResultsPayload,
                                      propagateResultsArgsFn = propagateProtEnrichResultsArgs,
                                      propagateUiParamsFn = propagateProtEnrichUiParams,
                                      updateStateManagerUiParamsFn = updateProtEnrichStateManagerUiParams,
                                      saveCompletedStateFn = saveProtEnrichCompletedState,
                                      completeTabStatusFn = completeProtEnrichTabStatus,
                                      completeProgressFn = completeProtEnrichProgress,
                                      finalizeAnalysisResultsFn = finalizeProtEnrichAnalysisBodyResults,
                                      showNotificationFn = shiny::showNotification,
                                      incProgressFn = shiny::incProgress,
                                      globalEnv = .GlobalEnv,
                                      catFn = cat) {
  selectedContrast <- input$selected_contrast

  incProgressFn(0.2, detail = "Transforming DE data...")

  incProgressFn(0.3, detail = "Creating DE results S4 object...")

  catFn("   ENRICHMENT Step: Creating DA results S4 object using createDAResultsForEnrichment\n")

  setupConfig <- prepareAnalysisSetupFn(
    selectedContrast = selectedContrast,
    input = input,
    enrichmentData = enrichmentData,
    workflowData = workflowData,
    experimentPaths = experimentPaths,
    resolveSelectedDaResultsFn = resolveSelectedDaResultsFn,
    resolveRunDependenciesFn = resolveRunDependenciesFn,
    resolveOutputDirectoriesFn = resolveOutputDirectoriesFn,
    createDaResultsForEnrichmentFn = createDaResultsForEnrichmentFn,
    resolveUniprotAnnotationsFn = resolveUniprotAnnotationsFn,
    resolveAnnotationMatchingFn = resolveAnnotationMatchingFn,
    resolveOrganismMappingFn = resolveOrganismMappingFn,
    applyOrganismFilterFn = applyOrganismFilterFn,
    persistOrganismFilterMetadataFn = persistOrganismFilterMetadataFn,
    showNotificationFn = showNotificationFn,
    globalEnv = globalEnv,
    catFn = catFn
  )
  rawContrastName <- setupConfig$rawContrastName
  contrastsTbl <- setupConfig$contrastsTbl
  pathwayDir <- setupConfig$pathwayDir
  daResultsForEnrichment <- setupConfig$daResultsForEnrichment
  organismFilterApplied <- setupConfig$organismFilterApplied
  filterStats <- setupConfig$filterStats

  incProgressFn(0.5, detail = "Running enrichment analysis...")

  catFn("   ENRICHMENT Step: Running processEnrichments\n")

  processExecutionConfig <- prepareProcessExecutionFn(
    input = input,
    enrichmentData = enrichmentData,
    daResultsForEnrichment = daResultsForEnrichment,
    pathwayDir = pathwayDir,
    goAnnotations = setupConfig$goAnnotations,
    currentAnalysisMethodFn = currentAnalysisMethodFn,
    resolveAnalysisInputColumnsFn = resolveAnalysisInputColumnsFn,
    buildProcessEnrichmentsArgsFn = buildProcessEnrichmentsArgsFn,
    catFn = catFn
  )
  methodInfo <- processExecutionConfig$methodInfo
  enrichmentArgs <- processExecutionConfig$enrichmentArgs

  enrichmentResults <- executeProcessEnrichmentsFn(
    enrichmentArgs = enrichmentArgs,
    upCutoff = input$up_cutoff,
    downCutoff = input$down_cutoff,
    qCutoff = input$q_cutoff,
    captureCheckpointFn = captureCheckpointFn,
    processEnrichmentsFn = processEnrichmentsFn,
    catFn = catFn
  )

  finalizeAnalysisResultsFn(
    selectedContrast = selectedContrast,
    rawContrastName = rawContrastName,
    organismFilterApplied = organismFilterApplied,
    filterStats = filterStats,
    enrichmentResults = enrichmentResults,
    enrichmentData = enrichmentData,
    workflowData = workflowData,
    input = input,
    methodInfo = methodInfo,
    contrastsTbl = contrastsTbl,
    pathwayDir = pathwayDir,
    buildAllContrastResultsFn = buildAllContrastResultsFn,
    resolveSelectedContrastResultsFn = resolveSelectedContrastResultsFn,
    capturePostProcessResultsFn = capturePostProcessResultsFn,
    persistAnalysisResultsFn = persistAnalysisResultsFn,
    buildAnalysisResultsPayloadFn = buildAnalysisResultsPayloadFn,
    propagateResultsArgsFn = propagateResultsArgsFn,
    propagateUiParamsFn = propagateUiParamsFn,
    updateStateManagerUiParamsFn = updateStateManagerUiParamsFn,
    saveCompletedStateFn = saveCompletedStateFn,
    completeTabStatusFn = completeTabStatusFn,
    completeProgressFn = completeProgressFn,
    incProgressFn = incProgressFn,
    catFn = catFn
  )
}

