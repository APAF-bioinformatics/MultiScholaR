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

    currentS4 <- readMetabWorkflowIds(stateManager) # KIT-082: the finder reads the id column and nothing else
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
    , revertStateFn = function(manager, stateName) selectWorkflowStateReference(manager, stateName) # KIT-075
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
    , metrics = NULL
) {
    setResolutionStatsFn(statsList)
    # KIT-032: the declared operation committed the state on the owner; this step publishes (KIT-072: one path).

    qcPlot <- tryCatch({
        # KIT-078: a worker's metrics stand for the object.
        if (is.null(currentS4) && is.null(metrics)) stop("the worker computed no filtering metrics")
        do.call(updateMetaboliteFilteringFn, c(list(
            theObject = currentS4
            , step_name = stepName
            , omics_type = omicType
            , return_grid = TRUE
            , overwrite = TRUE
        ), if (!is.null(metrics)) list(metrics = metrics)))
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
    , applyResolutionStateFn = applyMetabDuplicateResolutionState
    , buildResolutionSummaryFn = buildMetabDuplicateResolutionSummary
    , managedResult = NULL
) {
    # KIT-032: the resolution is the declared operation `metab_duplicate_resolution`, committed under
    # `metab_duplicates_resolved` with its provenance; the apply step then publishes (KIT-072: the owner's path is
    # the only path).
    # KIT-078: under a managed session the worker ran it (`managedResult`, which also carries what the figures read).
    resolved <- if (is.null(managedResult)) runMetabReferenceDuplicateResolution(workflowData$state_manager, workflowData$config_list) else
        list(value = managedResult$value, reference = managedResult$stateReference)
    # KIT-075: the resolved state is read for the QC plot only.
    resolutionPreflight <- list(currentS4 = if (is.null(managedResult)) readMetabWorkflowState(workflowData$state_manager, resolved$reference),
        statsList = resolved$value$stats)
    apply_arguments <- list(
        currentS4 = resolutionPreflight$currentS4
        , statsList = resolutionPreflight$statsList
        , workflowData = workflowData
        , omicType = omicType
        , setResolutionStatsFn = setResolutionStatsFn
        , setFilterPlotFn = setFilterPlotFn
    )
    if (!is.null(managedResult)) apply_arguments$metrics <- readWorkflowOmicsFilteringMetrics(workflowData$state_manager, managedResult)
    resolutionApply <- do.call(applyResolutionStateFn, apply_arguments)
    resultSummary <- buildResolutionSummaryFn(statsList = resolutionPreflight$statsList, stateName = resolutionApply$stateName)
    list(resultText = resultSummary$resultText, totalRemoved = resultSummary$totalRemoved)
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

    # KIT-078: on the worker under a managed session.
    runner <- workflowOmicsManagedActions("metabolomics")
    if (!is.null(runner)) {
        runner$run(list(omicsManagedCommit(omicsDuplicateResolutionStep(metabolomicsDeclarationBinding), workflowData$config_list, function(result) {
            runResolutionObserverShellFn(
                runResolutionFn = function() {
                    runResolutionWorkflowFn(
                        workflowData = workflowData
                        , omicType = omicType
                        , setResolutionStatsFn = setResolutionStatsFn
                        , setFilterPlotFn = setFilterPlotFn
                        , managedResult = result
                    )
                }
                , output = output
                , setDuplicateInfoFn = setDuplicateInfoFn
            )
        })), failed = function(error) shiny::removeNotification("metab_dup_resolve_working"))
        return(invisible(NULL))
    }
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

