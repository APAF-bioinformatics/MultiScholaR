reportMetabDuplicateDetection <- function(
    totalDuplicates
    , logInfoFn = logger::log_info
    , showNotificationFn = shiny::showNotification
    , sprintfFn = sprintf
) {
    logMessage <- paste(
        "Detected"
        , totalDuplicates
        , "duplicate feature IDs across assays"
    )
    notificationMessage <- sprintfFn(
        "Detection complete: %d duplicate IDs found"
        , totalDuplicates
    )
    notificationType <- if (totalDuplicates > 0) "warning" else "message"

    logInfoFn(logMessage)
    showNotificationFn(notificationMessage, type = notificationType)

    invisible(list(
        logMessage = logMessage
        , notificationMessage = notificationMessage
        , notificationType = notificationType
    ))
}

reportMetabDuplicateDetectionError <- function(
    errorCondition
    , prefix = "Error detecting duplicates:"
    , conditionMessageFn = conditionMessage
    , pasteFn = paste
    , logErrorFn = logger::log_error
    , showNotificationFn = shiny::showNotification
) {
    errorMessage <- pasteFn(prefix, conditionMessageFn(errorCondition))

    logErrorFn(errorMessage)
    showNotificationFn(errorMessage, type = "error")

    invisible(list(
        errorMessage = errorMessage
        , notificationType = "error"
    ))
}

runMetabDuplicateDetectionObserverShell <- function(
    stateManager
    , setDuplicateInfoFn
    , detectDuplicatesFn = detectMetabDuplicateFeatures
    , reportDetectionFn = reportMetabDuplicateDetection
    , reportDetectionErrorFn = reportMetabDuplicateDetectionError
) {
    tryCatch({
        detectionResult <- detectDuplicatesFn(
            stateManager = stateManager
        )
        setDuplicateInfoFn(detectionResult$duplicatesList)

        reportDetectionFn(
            totalDuplicates = detectionResult$totalDuplicates
        )

        invisible(list(
            status = "success"
            , duplicatesList = detectionResult$duplicatesList
            , totalDuplicates = detectionResult$totalDuplicates
        ))
    }, error = function(e) {
        reportDetectionErrorFn(
            errorCondition = e
        )

        invisible(list(
            status = "error"
            , errorCondition = e
        ))
    })
}

runMetabDuplicateResolutionObserverShell <- function(
    runResolutionFn
    , output
    , setDuplicateInfoFn
    , renderTextFn = shiny::renderText
    , logInfoFn = logger::log_info
    , logErrorFn = logger::log_error
    , showNotificationFn = shiny::showNotification
    , removeNotificationFn = shiny::removeNotification
) {
    tryCatch({
        resolutionDispatch <- runResolutionFn()

        output$resolution_results <- renderTextFn(
            resolutionDispatch$resultText
        )
        setDuplicateInfoFn(NULL)

        logInfoFn(paste(
            "Resolved duplicates: removed"
            , resolutionDispatch$totalRemoved
            , "rows"
        ))

        removeNotificationFn("metab_dup_resolve_working")
        showNotificationFn(
            sprintf(
                "Duplicates resolved: %d rows removed"
                , resolutionDispatch$totalRemoved
            )
            , type = "message"
        )

        invisible(list(
            status = "success"
            , totalRemoved = resolutionDispatch$totalRemoved
            , resultText = resolutionDispatch$resultText
        ))
    }, error = function(e) {
        msg <- paste("Error resolving duplicates:", e$message)

        logErrorFn(msg)
        showNotificationFn(msg, type = "error", duration = 15)
        removeNotificationFn("metab_dup_resolve_working")

        invisible(list(
            status = "error"
            , errorCondition = e
            , errorMessage = msg
        ))
    })
}

runMetabDuplicateRevertObserverShell <- function(
    runRevertFn
    , output
    , setResolutionStatsFn
    , setDuplicateInfoFn
    , setFilterPlotFn
    , renderTextFn = shiny::renderText
    , logInfoFn = logger::log_info
    , logErrorFn = logger::log_error
    , showNotificationFn = shiny::showNotification
) {
    tryCatch({
        revertDispatch <- runRevertFn()

        output$resolution_results <- renderTextFn(
            revertDispatch$resultText
        )
        setResolutionStatsFn(NULL)
        setDuplicateInfoFn(NULL)
        setFilterPlotFn(NULL)

        logInfoFn(paste(
            "Reverted duplicate resolution to"
            , revertDispatch$previousStateName
        ))
        showNotificationFn("Reverted successfully", type = "message")

        invisible(list(
            status = "success"
            , previousStateName = revertDispatch$previousStateName
            , resultText = revertDispatch$resultText
        ))
    }, error = function(e) {
        msg <- paste("Error reverting:", e$message)

        logErrorFn(msg)
        showNotificationFn(msg, type = "error")

        invisible(list(
            status = "error"
            , errorCondition = e
            , errorMessage = msg
        ))
    })
}

