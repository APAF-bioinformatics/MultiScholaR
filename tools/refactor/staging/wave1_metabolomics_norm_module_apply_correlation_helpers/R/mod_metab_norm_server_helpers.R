runMetabNormApplyCorrelationWorkflow <- function(
    observerState,
    reqFn = shiny::req,
    calculateCorrelationsFn = pearsonCorForSamplePairs,
    filterSamplesFn = filterSamplesByMetaboliteCorrelationThreshold,
    logInfoFn = logger::log_info
) {
    currentS4 <- observerState$currentS4
    threshold <- observerState$threshold

    reqFn(currentS4)

    logInfoFn("Calculating Pearson correlations per sample pair...")
    corrResults <- calculateCorrelationsFn(
        theObject = currentS4
        , correlation_group = observerState$groupingVariable
    )

    filteredS4 <- filterSamplesFn(
        theObject = currentS4
        , pearson_correlation_per_pair = corrResults
        , min_pearson_correlation_threshold = threshold
    )

    invisible(list(
        corrResults = corrResults
        , filteredS4 = filteredS4
    ))
}

dispatchMetabNormApplyCorrelation <- function(
    workflowData,
    normData,
    observerState,
    addLogFn = function(entry) invisible(entry),
    showNotificationFn = shiny::showNotification,
    removeNotificationFn = shiny::removeNotification,
    reqFn = shiny::req,
    calculateCorrelationsFn = pearsonCorForSamplePairs,
    filterSamplesFn = filterSamplesByMetaboliteCorrelationThreshold,
    runWorkflowFn = runMetabNormApplyCorrelationWorkflow,
    handleOutcomeFn = handleMetabNormApplyCorrelationOutcome,
    logInfoFn = logger::log_info,
    logErrorFn = logger::log_error
) {
    tryCatch({
        workflowState <- runWorkflowFn(
            observerState = observerState
            , reqFn = reqFn
            , calculateCorrelationsFn = calculateCorrelationsFn
            , filterSamplesFn = filterSamplesFn
            , logInfoFn = logInfoFn
        )

        handleOutcomeFn(
            workflowData = workflowData,
            normData = normData,
            observerState = observerState,
            corrResults = workflowState$corrResults,
            filteredS4 = workflowState$filteredS4,
            addLogFn = addLogFn,
            showNotificationFn = showNotificationFn,
            removeNotificationFn = removeNotificationFn,
            logErrorFn = logErrorFn
        )
    }, error = function(e) {
        handleOutcomeFn(
            workflowData = workflowData,
            normData = normData,
            observerState = observerState,
            error = e,
            addLogFn = addLogFn,
            showNotificationFn = showNotificationFn,
            removeNotificationFn = removeNotificationFn,
            logErrorFn = logErrorFn
        )
    })
}

handleMetabNormApplyCorrelationOutcome <- function(
    workflowData,
    normData,
    observerState,
    corrResults = NULL,
    filteredS4 = NULL,
    error = NULL,
    addLogFn = function(entry) invisible(entry),
    showNotificationFn = shiny::showNotification,
    removeNotificationFn = shiny::removeNotification,
    logErrorFn = logger::log_error
) {
    if (is.null(error)) {
        threshold <- observerState$threshold

        normData$correlation_results <- corrResults
        normData$correlation_filtered_obj <- filteredS4
        normData$correlation_filtering_complete <- TRUE

        workflowData$state_manager$saveState(
            state_name = "metab_correlation_filtered"
            , s4_data_object = filteredS4
            , config_object = workflowData$config_list
            , description = paste("Correlation filtering (threshold:", threshold, ")")
        )

        updatedStatus <- workflowData$tab_status
        updatedStatus$quality_control <- "complete"
        updatedStatus$normalization <- "complete"
        workflowData$tab_status <- updatedStatus

        addLogFn("Correlation filtering complete")
        removeNotificationFn(observerState$notificationId)
        showNotificationFn("Correlation filtering complete! Ready for DE analysis.", type = "message")

        return(invisible(list(
            status = "success"
            , corrResults = corrResults
            , filteredS4 = filteredS4
            , updatedStatus = updatedStatus
        )))
    }

    errorMessage <- if (inherits(error, "condition")) conditionMessage(error) else as.character(error)

    addLogFn(paste("ERROR in correlation filtering:", errorMessage))
    logErrorFn(paste("Correlation filtering error:", errorMessage))
    removeNotificationFn(observerState$notificationId)
    showNotificationFn(paste("Error:", errorMessage), type = "error")

    invisible(list(
        status = "error"
        , errorMessage = errorMessage
    ))
}

runMetabNormApplyCorrelationObserverEntry <- function(
    workflowData,
    normData,
    threshold,
    groupingVariable,
    addLogFn = function(entry) invisible(entry),
    showNotificationFn = shiny::showNotification,
    removeNotificationFn = shiny::removeNotification,
    resolveInputObjectFn = resolveMetabNormSkipCorrelationInputObject,
    dispatchApplyCorrelationFn = dispatchMetabNormApplyCorrelation
) {
    currentS4 <- resolveInputObjectFn(
        ruvCorrectedObject = normData$ruv_corrected_obj,
        postNormObject = normData$post_norm_obj
    )
    logEntry <- paste("Applying correlation filter (threshold:", threshold, ")")
    notificationId <- "corr_working"

    addLogFn(logEntry)
    showNotificationFn("Applying correlation filter...", id = notificationId, duration = NULL)

    observerState <- list(
        currentS4 = currentS4,
        threshold = threshold,
        groupingVariable = groupingVariable,
        logEntry = logEntry,
        notificationId = notificationId
    )

    dispatchState <- dispatchApplyCorrelationFn(
        workflowData = workflowData,
        normData = normData,
        observerState = observerState,
        addLogFn = addLogFn,
        showNotificationFn = showNotificationFn,
        removeNotificationFn = removeNotificationFn
    )

    invisible(dispatchState)
}

