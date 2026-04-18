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

