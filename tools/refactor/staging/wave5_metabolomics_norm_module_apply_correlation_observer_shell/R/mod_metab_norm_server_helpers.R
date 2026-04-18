runMetabNormApplyCorrelationObserverShell <- function(
    workflowData,
    normData,
    threshold,
    groupingVariable,
    addLogFn = function(entry) invisible(entry),
    showNotificationFn = shiny::showNotification,
    removeNotificationFn = shiny::removeNotification,
    reqFn = shiny::req,
    runObserverEntryFn = runMetabNormApplyCorrelationObserverEntry
) {
    reqFn(workflowData$state_manager)
    reqFn(normData$ruv_complete || normData$normalization_complete)

    dispatchState <- runObserverEntryFn(
        workflowData = workflowData,
        normData = normData,
        threshold = threshold,
        groupingVariable = groupingVariable,
        addLogFn = addLogFn,
        showNotificationFn = showNotificationFn,
        removeNotificationFn = removeNotificationFn
    )

    invisible(dispatchState)
}

