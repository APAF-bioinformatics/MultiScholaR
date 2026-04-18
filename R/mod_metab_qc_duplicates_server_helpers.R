resolveMetabDuplicateAssayData <- function(
    assayList
    , metaboliteIdCol
    , assayNames = names(assayList)
    , resolveDuplicateFeaturesByIntensityFn = resolveDuplicateFeaturesByIntensity
    , logWarnFn = logger::log_warn
) {
    statsList <- list()

    resolvedAssayList <- lapply(seq_along(assayList), function(i) {
        assayData <- assayList[[i]]
        assayName <- if (!is.null(assayNames)) assayNames[i] else paste0("Assay_", i)

        originalRows <- nrow(assayData)

        sampleCols <- names(assayData)[sapply(assayData, is.numeric)]

        if (length(sampleCols) == 0) {
            logWarnFn(paste("No numeric columns found in assay:", assayName))
            statsList[[assayName]] <<- list(
                original = originalRows
                , resolved = originalRows
                , removed = 0
            )
            return(assayData)
        }

        resolvedAssay <- resolveDuplicateFeaturesByIntensityFn(
            assay_tibble = assayData
            , id_col = metaboliteIdCol
            , sample_cols = sampleCols
        )

        resolvedRows <- nrow(resolvedAssay)

        statsList[[assayName]] <<- list(
            original = originalRows
            , resolved = resolvedRows
            , removed = originalRows - resolvedRows
        )

        resolvedAssay
    })

    names(resolvedAssayList) <- assayNames

    list(
        resolvedAssayList = resolvedAssayList
        , statsList = statsList
    )
}

buildMetabDuplicateResolutionSummary <- function(
    statsList
    , stateName = "metab_duplicates_resolved"
) {
    perAssayText <- if (length(statsList) == 0) {
        "  No assay statistics available"
    } else {
        paste(vapply(names(statsList), function(assayName) {
            stats <- statsList[[assayName]]
            sprintf(
                "  %s: %d -> %d rows (removed %d duplicates)"
                , assayName
                , stats$original
                , stats$resolved
                , stats$removed
            )
        }, character(1)), collapse = "\n")
    }

    totalRemoved <- if (length(statsList) == 0) {
        0
    } else {
        sum(vapply(statsList, `[[`, numeric(1), "removed"))
    }

    list(
        totalRemoved = totalRemoved
        , resultText = paste(
            "Duplicate Resolution Complete"
            , "=============================="
            , "Strategy: Keep feature with highest average intensity"
            , ""
            , "Per-Assay Results:"
            , perAssayText
            , ""
            , sprintf("Total duplicate rows removed: %d", totalRemoved)
            , ""
            , sprintf("State saved as: '%s'", stateName)
            , sep = "\n"
        )
    )
}

detectMetabDuplicateFeatures <- function(
    stateManager
    , duplicateFinderFn = findMetabDuplicateFeatureIDs
    , reqFn = shiny::req
    , inheritsFn = inherits
    , expectedClass = "MetaboliteAssayData"
) {
    reqFn(stateManager)

    currentS4 <- stateManager$getState()
    reqFn(currentS4)

    if (!inheritsFn(currentS4, expectedClass)) {
        stop(sprintf("Current state is not a %s object", expectedClass))
    }

    duplicatesList <- duplicateFinderFn(currentS4)
    totalDuplicates <- sum(vapply(duplicatesList, function(dupDf) {
        if (is.null(dupDf)) {
            return(0L)
        }

        nrow(dupDf)
    }, integer(1)))

    list(
        duplicatesList = duplicatesList
        , totalDuplicates = totalDuplicates
    )
}

revertMetabDuplicateResolution <- function(
    stateManager
    , reqFn = shiny::req
    , historyGetterFn = function(manager) manager$getHistory()
    , revertStateFn = function(manager, stateName) manager$revertToState(stateName)
) {
    reqFn(stateManager)

    history <- historyGetterFn(stateManager)
    if (length(history) <= 1) {
        stop("No previous state to revert to.")
    }

    previousStateName <- history[[length(history) - 1L]]
    revertStateFn(stateManager, previousStateName)

    list(
        previousStateName = previousStateName
        , resultText = paste("Reverted to previous state:", previousStateName)
    )
}

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

prepareMetabDuplicateResolutionState <- function(
    stateManager
    , resolveDuplicateAssayDataFn = resolveMetabDuplicateAssayData
    , reqFn = shiny::req
    , inheritsFn = inherits
    , expectedClass = "MetaboliteAssayData"
) {
    currentS4 <- stateManager$getState()
    reqFn(currentS4)

    if (!inheritsFn(currentS4, expectedClass)) {
        stop(sprintf("Current state is not a %s object", expectedClass))
    }

    resolutionResult <- resolveDuplicateAssayDataFn(
        assayList = currentS4@metabolite_data
        , metaboliteIdCol = currentS4@metabolite_id_column
    )

    currentS4@metabolite_data <- resolutionResult$resolvedAssayList

    list(
        currentS4 = currentS4
        , statsList = resolutionResult$statsList
    )
}

applyMetabDuplicateResolutionState <- function(
    currentS4
    , statsList
    , workflowData
    , omicType
    , setResolutionStatsFn
    , setFilterPlotFn
    , stateName = "metab_duplicates_resolved"
    , description = "Resolved duplicate metabolite features by keeping highest intensity"
    , stepName = "3_Duplicates_Resolved"
    , updateMetaboliteFilteringFn = updateMetaboliteFiltering
    , logWarnFn = logger::log_warn
) {
    setResolutionStatsFn(statsList)

    workflowData$state_manager$saveState(
        state_name = stateName
        , s4_data_object = currentS4
        , config_object = workflowData$config_list
        , description = description
    )

    qcPlot <- tryCatch({
        updateMetaboliteFilteringFn(
            theObject = currentS4
            , step_name = stepName
            , omics_type = omicType
            , return_grid = TRUE
            , overwrite = TRUE
        )
    }, error = function(e) {
        logWarnFn(paste("Could not generate QC plot:", e$message))
        NULL
    })

    setFilterPlotFn(qcPlot)

    list(
        stateName = stateName
        , qcPlot = qcPlot
    )
}

# Keep the remaining resolve-duplicates dispatch top-level so later waves can
# move it without reopening the observer body.
runMetabDuplicateResolutionWorkflow <- function(
    workflowData
    , omicType
    , setResolutionStatsFn
    , setFilterPlotFn
    , prepareResolutionStateFn = prepareMetabDuplicateResolutionState
    , applyResolutionStateFn = applyMetabDuplicateResolutionState
    , buildResolutionSummaryFn = buildMetabDuplicateResolutionSummary
) {
    resolutionPreflight <- prepareResolutionStateFn(
        stateManager = workflowData$state_manager
    )

    resolutionApply <- applyResolutionStateFn(
        currentS4 = resolutionPreflight$currentS4
        , statsList = resolutionPreflight$statsList
        , workflowData = workflowData
        , omicType = omicType
        , setResolutionStatsFn = setResolutionStatsFn
        , setFilterPlotFn = setFilterPlotFn
    )

    resultSummary <- buildResolutionSummaryFn(
        statsList = resolutionPreflight$statsList
        , stateName = resolutionApply$stateName
    )

    list(
        resultText = resultSummary$resultText
        , totalRemoved = resultSummary$totalRemoved
    )
}

runMetabDuplicateResolutionObserver <- function(
    workflowData
    , omicType
    , output
    , setDuplicateInfoFn
    , setResolutionStatsFn
    , setFilterPlotFn
    , reqFn = shiny::req
    , showNotificationFn = shiny::showNotification
    , runResolutionObserverShellFn = runMetabDuplicateResolutionObserverShell
    , runResolutionWorkflowFn = runMetabDuplicateResolutionWorkflow
) {
    reqFn(workflowData$state_manager)

    showNotificationFn(
        "Resolving duplicate features..."
        , id = "metab_dup_resolve_working"
        , duration = NULL
    )

    runResolutionObserverShellFn(
        runResolutionFn = function() {
            runResolutionWorkflowFn(
                workflowData = workflowData
                , omicType = omicType
                , setResolutionStatsFn = setResolutionStatsFn
                , setFilterPlotFn = setFilterPlotFn
            )
        }
        , output = output
        , setDuplicateInfoFn = setDuplicateInfoFn
    )
}

