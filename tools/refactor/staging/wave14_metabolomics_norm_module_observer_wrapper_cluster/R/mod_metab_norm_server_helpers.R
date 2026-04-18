# Keep the reset-normalization observer wrapper top-level so a later
# extraction wave can move it without reopening the module server wrapper.
runMetabNormResetNormalizationObserverWrapper <- function(
    workflowData,
    normData,
    addLogFn = function(entry) invisible(entry),
    showNotificationFn = shiny::showNotification,
    reqFn = shiny::req,
    runObserverShellFn = runMetabNormResetNormalizationObserverShell
) {
    runObserverShellFn(
        workflowData = workflowData
        , normData = normData
        , addLogFn = addLogFn
        , showNotificationFn = showNotificationFn
        , reqFn = reqFn
    )
}

# Keep the apply-correlation observer wrapper top-level so a later extraction
# wave can move it without reopening the module server wrapper.
runMetabNormApplyCorrelationObserverWrapper <- function(
    workflowData,
    input,
    normData,
    addLogFn = function(entry) invisible(entry),
    showNotificationFn = shiny::showNotification,
    removeNotificationFn = shiny::removeNotification,
    reqFn = shiny::req,
    runObserverShellFn = runMetabNormApplyCorrelationObserverShell,
    runObserverEntryFn = runMetabNormApplyCorrelationObserverEntry
) {
    runObserverShellFn(
        workflowData = workflowData
        , normData = normData
        , threshold = input$min_pearson_correlation_threshold
        , groupingVariable = input$ruv_grouping_variable
        , addLogFn = addLogFn
        , showNotificationFn = showNotificationFn
        , removeNotificationFn = removeNotificationFn
        , reqFn = reqFn
        , runObserverEntryFn = runObserverEntryFn
    )
}

# Keep the skip-correlation observer wrapper top-level so a later extraction
# wave can move it without reopening the module server wrapper.
runMetabNormSkipCorrelationObserverWrapper <- function(
    workflowData,
    normData,
    addLogFn = function(entry) invisible(entry),
    showNotificationFn = shiny::showNotification,
    reqFn = shiny::req,
    runObserverShellFn = runMetabNormSkipCorrelationObserverShell,
    runObserverEntryFn = runMetabNormSkipCorrelationObserverEntry
) {
    runObserverShellFn(
        workflowData = workflowData
        , normData = normData
        , addLogFn = addLogFn
        , showNotificationFn = showNotificationFn
        , reqFn = reqFn
        , runObserverEntryFn = runObserverEntryFn
    )
}

# Keep the export-session observer wrapper top-level so a later extraction
# wave can move it without reopening the module server wrapper.
runMetabNormExportSessionObserverWrapper <- function(
    workflowData,
    input,
    normData,
    experimentPaths,
    experimentLabel,
    addLogFn = function(entry) invisible(entry),
    logInfoFn = logger::log_info,
    reqFn = shiny::req,
    runObserverShellFn = runMetabNormExportSessionObserverShell,
    checkReadyFn = checkMetabNormExportSessionReady,
    dispatchExportSessionFn = dispatchMetabNormExportSession
) {
    inputValues <- list(
        export_session = input$export_session
        , norm_method = input$norm_method
        , ruv_mode = input$ruv_mode
        , apply_itsd = input$apply_itsd
        , itsd_aggregation = input$itsd_aggregation
        , log_offset = input$log_offset
        , min_pearson_correlation_threshold = input$min_pearson_correlation_threshold
        , ruv_grouping_variable = input$ruv_grouping_variable
    )

    runObserverShellFn(
        workflowData = workflowData
        , normData = normData
        , inputValues = inputValues
        , experimentPaths = experimentPaths
        , experimentLabel = experimentLabel
        , addLogFn = addLogFn
        , logInfoFn = logInfoFn
        , reqFn = reqFn
        , checkReadyFn = checkReadyFn
        , dispatchExportSessionFn = dispatchExportSessionFn
    )
}

