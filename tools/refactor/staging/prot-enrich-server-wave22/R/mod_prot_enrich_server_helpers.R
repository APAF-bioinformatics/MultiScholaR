captureProtEnrichPostProcessResults <- function(selectedContrast,
                                                enrichmentResults,
                                                enrichmentData,
                                                contrastsTbl,
                                                methodInfo,
                                                buildAllContrastResultsFn = buildProtEnrichAllContrastResults,
                                                resolveSelectedContrastResultsFn = resolveProtEnrichSelectedContrastResults,
                                                catFn = cat) {
  if (!is.null(enrichmentResults) && !is.null(enrichmentResults@enrichment_data)) {
    allContrastResults <- buildAllContrastResultsFn(
      enrichmentResults = enrichmentResults,
      methodInfo = methodInfo
    )

    enrichmentData$all_enrichment_results <- allContrastResults

    initialContrastState <- resolveSelectedContrastResultsFn(
      selectedContrast,
      allContrastResults,
      contrastsTbl
    )

    if (initialContrastState$found) {
      enrichmentData$gprofiler_results <- initialContrastState$gprofilerResults
      enrichmentData$clusterprofiler_results <- initialContrastState$clusterprofilerResults
      enrichmentData$stringdb_results <- initialContrastState$stringdbResults
      catFn(sprintf(
        "   ENRICHMENT Step: Set initial display to contrast '%s'\n",
        initialContrastState$rawContrastName
      ))
    }

    return(list(
      allContrastResults = allContrastResults,
      initialContrastState = initialContrastState,
      hasResults = TRUE
    ))
  }

  catFn("   ENRICHMENT Step: processEnrichments returned NULL or empty results\n")

  list(
    allContrastResults = NULL,
    initialContrastState = NULL,
    hasResults = FALSE
  )
}

persistProtEnrichAnalysisResults <- function(input,
                                             enrichmentData,
                                             workflowData,
                                             selectedContrast,
                                             enrichmentResults,
                                             methodInfo,
                                             pathwayDir,
                                             buildAnalysisResultsPayloadFn = buildProtEnrichAnalysisResultsPayload,
                                             propagateResultsArgsFn = propagateProtEnrichResultsArgs,
                                             propagateUiParamsFn = propagateProtEnrichUiParams,
                                             updateStateManagerUiParamsFn = updateProtEnrichStateManagerUiParams,
                                             saveCompletedStateFn = saveProtEnrichCompletedState,
                                             completeTabStatusFn = completeProtEnrichTabStatus,
                                             completeProgressFn = completeProtEnrichProgress,
                                             incProgressFn = shiny::incProgress) {
  incProgressFn(0.8, detail = "Storing results...")

  enrichmentData$enrichment_results_full <- enrichmentResults
  enrichmentData$analysis_complete <- TRUE

  workflowData$enrichment_analysis_results <- buildAnalysisResultsPayloadFn(
    gprofilerResults = enrichmentData$gprofiler_results,
    clusterprofilerResults = enrichmentData$clusterprofiler_results,
    stringdbResults = enrichmentData$stringdb_results,
    fullEnrichmentResults = enrichmentData$enrichment_results_full,
    selectedContrast = selectedContrast,
    analysisMethod = methodInfo$method,
    upCutoff = input$up_cutoff,
    downCutoff = input$down_cutoff,
    qCutoff = input$q_cutoff,
    organismTaxid = input$organism_taxid
  )

  argsPropagation <- propagateResultsArgsFn(
    enrichmentResults = enrichmentResults,
    currentS4Object = enrichmentData$current_s4_object,
    selectedContrast = selectedContrast,
    methodInfo = methodInfo,
    upCutoff = input$up_cutoff,
    downCutoff = input$down_cutoff,
    qCutoff = input$q_cutoff,
    organismTaxid = input$organism_taxid,
    pathwayDir = pathwayDir
  )
  enrichmentResults <- argsPropagation$enrichmentResults

  uiParamsPropagation <- propagateUiParamsFn(
    currentS4Object = enrichmentData$current_s4_object,
    workflowData = workflowData,
    selectedContrast = selectedContrast,
    methodInfo = methodInfo,
    upCutoff = input$up_cutoff,
    downCutoff = input$down_cutoff,
    qCutoff = input$q_cutoff,
    organismTaxid = input$organism_taxid
  )
  enrichmentData$current_s4_object <- uiParamsPropagation$currentS4Object

  updateStateManagerUiParamsFn(
    workflowData = workflowData,
    storedUiParams = uiParamsPropagation$storedUiParams
  )

  saveCompletedStateFn(
    workflowData = workflowData,
    enrichmentResults = enrichmentResults,
    selectedContrast = selectedContrast,
    methodInfo = methodInfo,
    upCutoff = input$up_cutoff,
    downCutoff = input$down_cutoff,
    qCutoff = input$q_cutoff,
    organismTaxid = input$organism_taxid,
    pathwayDir = pathwayDir
  )

  completeTabStatusFn(workflowData = workflowData)
  completeProgressFn()

  list(
    enrichmentResults = enrichmentResults,
    analysisComplete = enrichmentData$analysis_complete,
    currentS4Object = enrichmentData$current_s4_object
  )
}

finalizeProtEnrichAnalysisBodyResults <- function(selectedContrast,
                                                  rawContrastName,
                                                  organismFilterApplied,
                                                  filterStats,
                                                  enrichmentResults,
                                                  enrichmentData,
                                                  workflowData,
                                                  input,
                                                  methodInfo,
                                                  contrastsTbl,
                                                  pathwayDir,
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
                                                  incProgressFn = shiny::incProgress,
                                                  catFn = cat) {
  capturePostProcessResultsFn(
    selectedContrast = selectedContrast,
    enrichmentResults = enrichmentResults,
    enrichmentData = enrichmentData,
    contrastsTbl = contrastsTbl,
    methodInfo = methodInfo,
    buildAllContrastResultsFn = buildAllContrastResultsFn,
    resolveSelectedContrastResultsFn = resolveSelectedContrastResultsFn,
    catFn = catFn
  )

  persistenceResult <- persistAnalysisResultsFn(
    input = input,
    enrichmentData = enrichmentData,
    workflowData = workflowData,
    selectedContrast = selectedContrast,
    enrichmentResults = enrichmentResults,
    methodInfo = methodInfo,
    pathwayDir = pathwayDir,
    buildAnalysisResultsPayloadFn = buildAnalysisResultsPayloadFn,
    propagateResultsArgsFn = propagateResultsArgsFn,
    propagateUiParamsFn = propagateUiParamsFn,
    updateStateManagerUiParamsFn = updateStateManagerUiParamsFn,
    saveCompletedStateFn = saveCompletedStateFn,
    completeTabStatusFn = completeTabStatusFn,
    completeProgressFn = completeProgressFn,
    incProgressFn = incProgressFn
  )
  enrichmentResults <- persistenceResult$enrichmentResults

  list(
    selectedContrast = selectedContrast,
    rawContrastName = rawContrastName,
    analysisMethod = methodInfo$method,
    analysisComplete = persistenceResult$analysisComplete,
    organismFilterApplied = organismFilterApplied,
    filterStats = filterStats,
    enrichmentResults = enrichmentResults
  )
}

