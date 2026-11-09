runProtEnrichObserverShell <- function(selectedContrast,
                                       runAnalysisBody,
                                       withProgressFn = shiny::withProgress,
                                       finalizeFailureFn = finalizeProtEnrichObserverFailure,
                                       finalizeObserverRunFn = finalizeProtEnrichObserverRun,
                                       progressMessage = "Running enrichment analysis...",
                                       progressValue = 0) {
  analysisCompleted <- FALSE
  failureResult <- NULL
  successResult <- NULL

  tryCatch({
    withProgressFn(message = progressMessage, value = progressValue, {
      runAnalysisBody()
    })
    analysisCompleted <- TRUE
  }, error = function(e) {
    failureResult <<- finalizeFailureFn(errorMessage = e$message)
  })

  if (isTRUE(analysisCompleted)) {
    successResult <- finalizeObserverRunFn(
      completed = TRUE,
      selectedContrast = selectedContrast
    )
  }

  list(
    completed = analysisCompleted,
    selectedContrast = selectedContrast,
    failureResult = failureResult,
    successResult = successResult
  )
}

finalizeProtEnrichObserverFailure <- function(errorMessage,
                                              reportAnalysisErrorFn = reportProtEnrichAnalysisError,
                                              removeWorkingNotificationFn = removeProtEnrichWorkingNotification) {
  reportResult <- reportAnalysisErrorFn(errorMessage = errorMessage)
  cleanupResult <- removeWorkingNotificationFn()

  list(
    errorMessage = errorMessage,
    reportResult = reportResult,
    cleanupResult = cleanupResult
  )
}

reportProtEnrichCompletion <- function(selectedContrast,
                                       notifyCompletionFn = notifyProtEnrichCompletion,
                                       logCompletionFn = logProtEnrichCompletion) {
  notificationResult <- notifyCompletionFn(selectedContrast = selectedContrast)
  logResult <- logCompletionFn()

  list(
    selectedContrast = selectedContrast,
    notificationResult = notificationResult,
    logResult = logResult
  )
}

finalizeProtEnrichObserverRun <- function(completed,
                                          selectedContrast = NULL,
                                          reportCompletionFn = reportProtEnrichCompletion,
                                          removeWorkingNotificationFn = removeProtEnrichWorkingNotification) {
  reportResult <- NULL

  if (isTRUE(completed)) {
    reportResult <- reportCompletionFn(selectedContrast = selectedContrast)
  }

  cleanupResult <- removeWorkingNotificationFn()

  list(
    completed = isTRUE(completed),
    selectedContrast = selectedContrast,
    reportResult = reportResult,
    cleanupResult = cleanupResult
  )
}

logProtEnrichCompletion <- function(message = "=== ENRICHMENT ANALYSIS COMPLETED ===\n",
                                    catFn = cat) {
  catFn(message)

  list(message = message)
}

removeProtEnrichWorkingNotification <- function(notificationId = "enrichment_working",
                                                removeNotificationFn = shiny::removeNotification) {
  removeNotificationFn(notificationId)

  list(notificationId = notificationId)
}

handoffProtEnrichObserverRun <- function(input,
                                         enrichmentData,
                                         workflowData,
                                         experimentPaths,
                                         currentAnalysisMethodFn,
                                         showNotificationFn = shiny::showNotification,
                                         runObserverShellFn = runProtEnrichObserverShell,
                                         runAnalysisBodyFn = runProtEnrichAnalysisBody) {
  selectedContrast <- input$selected_contrast
  notificationMessage <- "Running enrichment analysis..."
  notificationId <- "enrichment_working"
  notificationDuration <- NULL

  showNotificationFn(
    notificationMessage,
    id = notificationId,
    duration = notificationDuration
  )

  shellResult <- runObserverShellFn(
    selectedContrast = selectedContrast,
    runAnalysisBody = function() {
      runAnalysisBodyFn(
        input = input,
        enrichmentData = enrichmentData,
        workflowData = workflowData,
        experimentPaths = experimentPaths,
        currentAnalysisMethodFn = currentAnalysisMethodFn
      )
    }
  )

  list(
    selectedContrast = selectedContrast,
    notificationMessage = notificationMessage,
    notificationId = notificationId,
    notificationDuration = notificationDuration,
    shellResult = shellResult
  )
}

runProtEnrichObserverPreflight <- function(input,
                                           enrichmentData,
                                           workflowData,
                                           experimentPaths,
                                           currentAnalysisMethodFn,
                                           reqFn = shiny::req,
                                           handoffObserverRunFn = handoffProtEnrichObserverRun) {
  selectedContrast <- input$selected_contrast
  daResultsData <- enrichmentData$da_results_data

  reqFn(selectedContrast, daResultsData)

  handoffResult <- handoffObserverRunFn(
    input = input,
    enrichmentData = enrichmentData,
    workflowData = workflowData,
    experimentPaths = experimentPaths,
    currentAnalysisMethodFn = currentAnalysisMethodFn
  )

  list(
    selectedContrast = selectedContrast,
    daResultsData = daResultsData,
    handoffResult = handoffResult
  )
}

createProtEnrichSupportedOrganismsReactive <- function(reactiveFn = shiny::reactive,
                                                       buildSupportedOrganismsFn = buildProtEnrichSupportedOrganisms) {
  reactiveFn({
    buildSupportedOrganismsFn()
  })
}

createProtEnrichReactiveValues <- function(reactiveValuesFn = shiny::reactiveValues) {
  reactiveValuesFn(
    enrichment_results = NULL,
    contrasts_available = NULL,
    analysis_complete = FALSE,
    current_s4_object = NULL,
    da_results_data = NULL,
    gprofiler_results = NULL,
    clusterprofiler_results = NULL,
    stringdb_results = NULL,
    analysis_method = NULL,
    organism_supported = NULL,
    all_enrichment_results = list(),
    current_contrast_results = list(),
    enrichment_plots = list()
  )
}
