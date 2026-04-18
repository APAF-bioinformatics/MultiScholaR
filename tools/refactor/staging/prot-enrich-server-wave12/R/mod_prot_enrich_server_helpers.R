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

