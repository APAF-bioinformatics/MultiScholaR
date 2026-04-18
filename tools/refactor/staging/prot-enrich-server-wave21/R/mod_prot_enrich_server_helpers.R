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

