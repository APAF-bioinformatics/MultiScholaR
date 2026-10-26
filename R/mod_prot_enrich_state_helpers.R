buildProtEnrichAnalysisResultsPayload <- function(gprofilerResults,
                                                  clusterprofilerResults,
                                                  stringdbResults,
                                                  fullEnrichmentResults,
                                                  selectedContrast,
                                                  analysisMethod,
                                                  upCutoff,
                                                  downCutoff,
                                                  qCutoff,
                                                  organismTaxid) {
  list(
    gprofiler_results = gprofilerResults,
    clusterprofiler_results = clusterprofilerResults,
    stringdb_results = stringdbResults,
    full_enrichment_results = fullEnrichmentResults,
    selected_contrast = selectedContrast,
    analysis_method = analysisMethod,
    parameters = list(
      up_cutoff = upCutoff,
      down_cutoff = downCutoff,
      q_cutoff = qCutoff,
      organism_taxid = organismTaxid
    )
  )
}

propagateProtEnrichResultsArgs <- function(enrichmentResults,
                                           currentS4Object,
                                           selectedContrast,
                                           methodInfo,
                                           upCutoff,
                                           downCutoff,
                                           qCutoff,
                                           organismTaxid,
                                           pathwayDir,
                                           timeFn = Sys.time,
                                           catFn = cat) {
  catFn("   ENRICHMENT Step: Copying @args from original data object...\n")

  dataHasArgs <- tryCatch({
    !is.null(currentS4Object) && !is.null(currentS4Object@args)
  }, error = function(e) {
    FALSE
  })

  resultsHasArgs <- tryCatch({
    !is.null(enrichmentResults@args)
  }, error = function(e) {
    FALSE
  })

  if (dataHasArgs) {
    if (resultsHasArgs) {
      enrichmentResults@args <- currentS4Object@args

      if (is.null(enrichmentResults@args$enrichmentAnalysis)) {
        enrichmentResults@args$enrichmentAnalysis <- list()
      }

      enrichmentResults@args$enrichmentAnalysis <- list(
        selected_contrast = selectedContrast,
        analysis_method = methodInfo$method,
        organism_supported = methodInfo$supported,
        up_cutoff = upCutoff,
        down_cutoff = downCutoff,
        q_cutoff = qCutoff,
        organism_taxid = organismTaxid,
        pathway_dir = pathwayDir
      )

      catFn("   ENRICHMENT Step: Storing UI parameters in @args\n")
      if (is.null(enrichmentResults@args$enrichmentAnalysisUI)) {
        enrichmentResults@args$enrichmentAnalysisUI <- list()
      }

      enrichmentResults@args$enrichmentAnalysisUI <- list(
        up_log2fc_cutoff = upCutoff,
        down_log2fc_cutoff = downCutoff,
        q_value_cutoff = qCutoff,
        organism_taxon_id = organismTaxid,
        analysis_method = methodInfo$method,
        organism_name = methodInfo$species_name,
        organism_supported = methodInfo$supported,
        selected_contrast = selectedContrast,
        timestamp = timeFn()
      )

      catFn("   ENRICHMENT Step: Successfully copied and updated @args\n")
    } else {
      catFn("   ENRICHMENT Step: EnrichmentResults doesn't have @args slot\n")
    }
  } else {
    catFn("   ENRICHMENT Step: Original data object doesn't have @args to copy\n")
  }

  list(
    enrichmentResults = enrichmentResults,
    dataHasArgs = dataHasArgs,
    resultsHasArgs = resultsHasArgs,
    copiedArgs = isTRUE(dataHasArgs) && isTRUE(resultsHasArgs)
  )
}

propagateProtEnrichUiParams <- function(currentS4Object,
                                        workflowData,
                                        selectedContrast,
                                        methodInfo,
                                        upCutoff,
                                        downCutoff,
                                        qCutoff,
                                        organismTaxid,
                                        timeFn = Sys.time,
                                        catFn = cat) {
  dataHasArgs <- tryCatch({
    !is.null(currentS4Object) && !is.null(currentS4Object@args)
  }, error = function(e) {
    FALSE
  })

  if (dataHasArgs) {
    catFn("   ENRICHMENT Step: Storing UI parameters in original data object @args\n")
    if (is.null(currentS4Object@args$enrichmentAnalysisUI)) {
      currentS4Object@args$enrichmentAnalysisUI <- list()
    }

    currentS4Object@args$enrichmentAnalysisUI <- list(
      up_log2fc_cutoff = upCutoff,
      down_log2fc_cutoff = downCutoff,
      q_value_cutoff = qCutoff,
      organism_taxon_id = organismTaxid,
      analysis_method = methodInfo$method,
      organism_name = methodInfo$species_name,
      organism_supported = methodInfo$supported,
      selected_contrast = selectedContrast,
      timestamp = timeFn()
    )

    workflowData$enrichment_ui_params <- list(
      up_log2fc_cutoff = upCutoff,
      down_log2fc_cutoff = downCutoff,
      q_value_cutoff = qCutoff,
      organism_selected = organismTaxid,
      database_source = methodInfo$method,
      organism_name = methodInfo$species_name,
      organism_supported = methodInfo$supported,
      selected_contrast = selectedContrast,
      timestamp = timeFn()
    )
    catFn("   ENRICHMENT Step: Stored UI parameters in workflow_data for sessionSummary\n")
  }

  list(
    currentS4Object = currentS4Object,
    dataHasArgs = dataHasArgs,
    storedUiParams = isTRUE(dataHasArgs),
    workflowUiParams = if (isTRUE(dataHasArgs)) workflowData$enrichment_ui_params else NULL
  )
}

updateProtEnrichStateManagerUiParams <- function(workflowData,
                                                 storedUiParams,
                                                 currentDataStates = c(
                                                   "ruv_covariates",
                                                   "correlation_filtered",
                                                   "ruv_corrected",
                                                   "protein_replicate_filtered"
                                                 ),
                                                 detectFn = purrr::detect,
                                                 catFn = cat) {
  if (!isTRUE(storedUiParams)) {
    return(list(
      attempted = FALSE,
      currentDataState = NULL,
      availableStates = NULL,
      warning = NULL,
      updated = FALSE
    ))
  }

  catFn("   ENRICHMENT Step: Updating R6 state with enrichment UI parameters\n")

  currentDataState <- NULL
  availableStates <- NULL
  warningMessage <- NULL

  tryCatch({
    availableStates <- diaStateHistory(workflowData$state_manager)
    currentDataState <- detectFn(currentDataStates, ~ .x %in% availableStates)

    if (!is.null(currentDataState)) {
      catFn("   ENRICHMENT Step: Skipping state update (updateState method verification needed)\n")
    }
  }, error = function(e) {
    warningMessage <<- e$message
    catFn(sprintf(
      "   ENRICHMENT Step: Warning - could not update state with UI parameters: %s\n",
      e$message
    ))
  })

  list(
    attempted = TRUE,
    currentDataState = currentDataState,
    availableStates = availableStates,
    warning = warningMessage,
    updated = FALSE
  )
}

saveProtEnrichCompletedState <- function(workflowData,
                                         enrichmentResults,
                                         selectedContrast,
                                         methodInfo,
                                         upCutoff,
                                         downCutoff,
                                         qCutoff,
                                         organismTaxid,
                                         pathwayDir,
                                         catFn = cat) {
  catFn("   ENRICHMENT Step: Saving results to R6 state manager...\n")

  warningMessage <- NULL
  saved <- FALSE

  tryCatch({
    saveDIAOwnerState(workflowData$state_manager, 
      state_name = "enrichment_completed",
      s4_data_object = enrichmentResults,
      config_object = list(
        selected_contrast = selectedContrast,
        analysis_method = methodInfo$method,
        organism_supported = methodInfo$supported,
        up_cutoff = upCutoff,
        down_cutoff = downCutoff,
        q_cutoff = qCutoff,
        organism_taxid = organismTaxid,
        pathway_dir = pathwayDir
      ),
      description = paste(
        "Enrichment analysis completed using",
        methodInfo$method,
        "for contrast:",
        selectedContrast
      )
    )

    saved <- TRUE
    catFn("   ENRICHMENT Step: Successfully saved state to R6 state manager\n")
  }, error = function(e) {
    warningMessage <<- e$message
    catFn(sprintf(
      "   ENRICHMENT Step: Warning saving to R6 state manager: %s\n",
      e$message
    ))
  })

  list(
    attempted = TRUE,
    saved = saved,
    warning = warningMessage
  )
}

completeProtEnrichTabStatus <- function(workflowData,
                                        tabName = "enrichment_analysis",
                                        status = "complete") {
  updatedStatus <- workflowData$tab_status
  updatedStatus[[tabName]] <- status
  workflowData$tab_status <- updatedStatus

  list(
    tabName = tabName,
    status = status,
    tabStatus = updatedStatus
  )
}

completeProtEnrichProgress <- function(value = 1.0,
                                       detail = "Complete!",
                                       incProgressFn = shiny::incProgress) {
  incProgressFn(value, detail = detail)

  list(
    value = value,
    detail = detail
  )
}

notifyProtEnrichCompletion <- function(selectedContrast,
                                       type = "message",
                                       duration = 5,
                                       showNotificationFn = shiny::showNotification) {
  notificationMessage <- paste(
    "Enrichment analysis completed successfully for contrast:",
    selectedContrast
  )

  showNotificationFn(
    notificationMessage,
    type = type,
    duration = duration
  )

  list(
    message = notificationMessage,
    type = type,
    duration = duration,
    selectedContrast = selectedContrast
  )
}

notifyProtEnrichAnalysisError <- function(errorMessage,
                                          type = "error",
                                          duration = 10,
                                          showNotificationFn = shiny::showNotification) {
  notificationMessage <- sprintf("Error in enrichment analysis: %s", errorMessage)

  showNotificationFn(
    notificationMessage,
    type = type,
    duration = duration
  )

  list(
    message = notificationMessage,
    type = type,
    duration = duration,
    errorMessage = errorMessage
  )
}

logProtEnrichAnalysisError <- function(errorMessage,
                                       template = "*** ERROR in enrichment analysis: %s ***\n",
                                       catFn = cat) {
  message <- sprintf(template, errorMessage)
  catFn(message)

  list(
    message = message,
    errorMessage = errorMessage
  )
}

messageProtEnrichAnalysisError <- function(errorMessage,
                                           template = "*** ERROR in enrichment analysis: %s",
                                           messageFn = message) {
  formattedMessage <- sprintf(template, errorMessage)
  messageFn(formattedMessage)

  list(
    message = formattedMessage,
    errorMessage = errorMessage
  )
}

reportProtEnrichAnalysisError <- function(errorMessage,
                                          messageErrorFn = messageProtEnrichAnalysisError,
                                          logErrorFn = logProtEnrichAnalysisError,
                                          notifyErrorFn = notifyProtEnrichAnalysisError) {
  messageResult <- messageErrorFn(errorMessage = errorMessage)
  logResult <- logErrorFn(errorMessage = errorMessage)
  notificationResult <- notifyErrorFn(errorMessage = errorMessage)

  list(
    errorMessage = errorMessage,
    messageResult = messageResult,
    logResult = logResult,
    notificationResult = notificationResult
  )
}
captureProtEnrichPostProcessResults <- function(selectedContrast,
                                                enrichmentResults,
                                                enrichmentData,
                                                contrastsTbl,
                                                methodInfo,
                                                buildAllContrastResultsFn = buildProtEnrichAllContrastResults,
                                                resolveSelectedContrastResultsFn = resolveProtEnrichSelectedContrastResults,
                                                catFn = cat,
                                                workflowData = NULL) {
  if (!is.null(enrichmentResults) && !is.null(enrichmentResults@enrichment_data)) {
    allContrastResults <- buildAllContrastResultsFn(
      enrichmentResults = enrichmentResults,
      methodInfo = methodInfo
    )
    # KIT-018: on the on-disk default each contrast's result table is also a saved display table the
    # results tables page from the store (content-addressed; the in-memory frames stay for plots and downloads).
    if (!is.null(workflowData)) {
      for (contrastName in names(allContrastResults)) {
        allContrastResults[[contrastName]]$gprofiler_results_table <- saveDIADisplayTable(
          workflowData$state_manager, allContrastResults[[contrastName]]$gprofiler_results)
        allContrastResults[[contrastName]]$clusterprofiler_results_table <- saveDIADisplayTable(
          workflowData$state_manager, allContrastResults[[contrastName]]$clusterprofiler_results)
      }
    }

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
      enrichmentData$gprofiler_results_table <- initialContrastState$gprofilerResultsTable
      enrichmentData$clusterprofiler_results_table <- initialContrastState$clusterprofilerResultsTable
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
    catFn = catFn,
    workflowData = workflowData
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

#' Restore the enrichment results from the `enrichment_completed` state (KIT-024): after an Open the module holds
#' no results although the state carries the KIT-003 `EnrichmentResults`; on first use they are rebuilt into the
#' per-contrast tables (and saved display tables on the on-disk default) exactly as the analysis builds them.
#' @return TRUE when restored, FALSE when nothing was needed or nothing was available.
#' @noRd
restoreProtEnrichResultsFromState <- function(enrichmentData, workflowData,
                                              buildAllContrastResultsFn = buildProtEnrichAllContrastResults) {
  if (!is.null(shiny::isolate(enrichmentData$all_enrichment_results))) return(FALSE)
  manager <- if (is.null(workflowData)) NULL else shiny::isolate(workflowData$state_manager)
  if (is.null(manager) || !"enrichment_completed" %in% unlist(diaStateHistory(manager))) return(FALSE)
  results <- readDIAStateObject(manager, "enrichment_completed")
  if (!methods::is(results, "EnrichmentResults")) return(FALSE)
  methodInfo <- list(method = diaEnrichmentMethodOf(results))
  allContrastResults <- buildAllContrastResultsFn(enrichmentResults = results, methodInfo = methodInfo)
  for (contrastName in names(allContrastResults)) {
    allContrastResults[[contrastName]]$gprofiler_results_table <- saveDIADisplayTable(manager,
      allContrastResults[[contrastName]]$gprofiler_results)
    allContrastResults[[contrastName]]$clusterprofiler_results_table <- saveDIADisplayTable(manager,
      allContrastResults[[contrastName]]$clusterprofiler_results)
  }
  enrichmentData$all_enrichment_results <- allContrastResults
  enrichmentData$enrichment_results_full <- results
  enrichmentData$analysis_complete <- TRUE
  TRUE
}
