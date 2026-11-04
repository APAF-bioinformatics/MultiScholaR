initializeLipidDuplicateServerState <- function(
    reactiveValFn = shiny::reactiveVal
) {
    list(
        duplicateInfo = reactiveValFn(NULL),
        resolutionStats = reactiveValFn(NULL),
        filterPlot = reactiveValFn(NULL)
    )
}

buildLipidDuplicateTablesUi <- function(dupList, ns) {
    if (is.null(dupList)) {
        return(NULL)
    }

    assaysWithDuplicates <- names(dupList)[vapply(
        dupList,
        function(dupDf) !is.null(dupDf) && nrow(dupDf) > 0,
        logical(1)
    )]

    if (length(assaysWithDuplicates) == 0) {
        return(shiny::wellPanel(
            shiny::icon("check-circle", style = "color: green; font-size: 24px;")
            , shiny::h5("No duplicates found in any assay!")
            , shiny::p("All lipid IDs are unique. No resolution needed.")
        ))
    }

    tabList <- lapply(assaysWithDuplicates, function(assayName) {
        shiny::tabPanel(
            assayName
            , shiny::br()
            , DT::DTOutput(ns(sanitizeLipidDuplicateTableOutputId(assayName)))
        )
    })

    do.call(shiny::tabsetPanel, c(list(id = ns("dup_tables_tabs")), tabList))
}

sanitizeLipidDuplicateTableOutputId <- function(assayName) {
    paste0("dup_table_", gsub("[^a-zA-Z0-9]", "_", assayName))
}

buildLipidDuplicateDatatable <- function(dupDf) {
    DT::datatable(
        dupDf
        , options = list(
            pageLength = 10
            , scrollX = TRUE
            , dom = "frtip"
        )
        , rownames = FALSE
        , class = "compact stripe"
    )
}

registerLipidDuplicateTableOutputs <- function(output, dupList, renderTable = NULL) {
    if (is.null(dupList)) {
        return(character())
    }

    if (is.null(renderTable)) {
        renderTable <- function(dupDf) {
            DT::renderDT({
                buildLipidDuplicateDatatable(dupDf)
            })
        }
    }

    registeredIds <- character()

    lapply(names(dupList), function(assayName) {
        dupDf <- dupList[[assayName]]

        if (is.null(dupDf) || nrow(dupDf) == 0) {
            return(NULL)
        }

        outputId <- sanitizeLipidDuplicateTableOutputId(assayName)
        output[[outputId]] <- renderTable(dupDf)
        registeredIds <<- c(registeredIds, outputId)

        NULL
    })

    registeredIds
}

registerLipidDuplicateTableObserver <- function(
    output,
    duplicateInfoVal,
    observeFn = shiny::observe,
    reqFn = shiny::req,
    registerTableOutputsFn = registerLipidDuplicateTableOutputs
) {
    observeFn({
        dupList <- duplicateInfoVal()
        reqFn(dupList)

        registerTableOutputsFn(
            output = output
            , dupList = dupList
        )
    })

    invisible(output)
}

buildLipidDuplicateSummaryUi <- function(dupList) {
    if (is.null(dupList)) {
        return(shiny::p(
            shiny::icon("info-circle")
            , " Click 'Detect Duplicates' to scan for duplicate features."
            , style = "color: #666;"
        ))
    }

    summaryItems <- lapply(names(dupList), function(assayName) {
        dupDf <- dupList[[assayName]]
        duplicateCount <- if (is.null(dupDf)) 0 else nrow(dupDf)

        iconClass <- if (duplicateCount == 0) "check-circle" else "exclamation-triangle"
        iconColor <- if (duplicateCount == 0) "green" else "orange"

        shiny::tags$li(
            shiny::icon(iconClass, style = paste("color:", iconColor))
            , sprintf(" %s: %d duplicate IDs", assayName, duplicateCount)
        )
    })

    shiny::tags$ul(summaryItems, style = "list-style: none; padding-left: 0;")
}

registerLipidDuplicateSummaryOutput <- function(
    output,
    duplicateInfoVal,
    renderUiFn = shiny::renderUI,
    buildSummaryUiFn = buildLipidDuplicateSummaryUi
) {
    output$duplicate_summary <- renderUiFn({
        buildSummaryUiFn(duplicateInfoVal())
    })

    invisible(output)
}

registerLipidDuplicateTablesOutput <- function(
    output,
    duplicateInfoVal,
    ns,
    renderUiFn = shiny::renderUI,
    buildTablesUiFn = buildLipidDuplicateTablesUi
) {
    output$duplicate_tables <- renderUiFn({
        buildTablesUiFn(
            dupList = duplicateInfoVal()
            , ns = ns
        )
    })

    invisible(output)
}

registerLipidDuplicateFilterPlotOutput <- function(
    output,
    filterPlotVal,
    renderPlotFn = shiny::renderPlot,
    reqFn = shiny::req,
    gridDrawFn = grid::grid.draw,
    printFn = print
) {
    output$filter_plot <- renderPlotFn({
        plotObject <- filterPlotVal()
        reqFn(plotObject)

        if (inherits(plotObject, "grob") || inherits(plotObject, "gtable")) {
            gridDrawFn(plotObject)
        } else if (inherits(plotObject, "ggplot")) {
            printFn(plotObject)
        }
    })

    invisible(output)
}

handleLipidDuplicateDetection <- function(
    workflowData,
    reqFn = shiny::req,
    findDuplicatesFn = findLipidDuplicateFeatureIDs,
    logInfoFn = logger::log_info
) {
    reqFn(workflowData$state_manager)

    # KIT-082: the finder reads the id column and nothing else.
    currentS4 <- readLipidWorkflowIds(workflowData$state_manager)
    reqFn(currentS4)

    if (!inherits(currentS4, "LipidomicsAssayData")) {
        stop("Current state is not a LipidomicsAssayData object")
    }

    duplicatesList <- findDuplicatesFn(currentS4)
    totalDuplicates <- sum(vapply(duplicatesList, function(dupDf) {
        if (is.null(dupDf)) {
            return(0)
        }

        nrow(dupDf)
    }, numeric(1)))

    logInfoFn(paste(
        "Detected"
        , totalDuplicates
        , "duplicate feature IDs across assays"
    ))

    list(
        duplicatesList = duplicatesList
        , totalDuplicates = totalDuplicates
        , notificationMessage = sprintf(
            "Detection complete: %d duplicate IDs found"
            , totalDuplicates
        )
        , notificationType = if (totalDuplicates > 0) "warning" else "message"
    )
}

handleLipidDuplicateResolution <- function(
    workflowData,
    omicType,
    reqFn = shiny::req,
    resolveDuplicatesFn = resolveLipidDuplicateFeaturesByIntensity,
    updateFilteringFn = updateLipidFiltering,
    logWarnFn = logger::log_warn,
    managedResult = NULL
) {
    reqFn(workflowData$state_manager)

    # KIT-075: the assay names from the state's metadata, the QC plot from the resolved state read for this action;
    # the result carries the state's reference, not the object.
    current <- readLipidWorkflowMetadata(workflowData$state_manager)
    reqFn(current)

    if (!identical(current$class, "LipidomicsAssayData")) {
        stop("Current state is not a LipidomicsAssayData object")
    }

    assayNames <- current$assay_names

    # KIT-031: the resolution is the declared operation `lipid_duplicate_resolution`, committed under
    # `lipid_duplicates_resolved` with its provenance (KIT-072: the owner's path is the only path).
    # KIT-078: under a managed session the worker ran it (`managedResult`, which also carries what the figures read).
    resolved <- if (is.null(managedResult)) runLipidReferenceDuplicateResolution(workflowData$state_manager, workflowData$config_list) else
        list(value = managedResult$value, reference = managedResult$stateReference)
    counts <- resolved$value$counts
    statsList <- lapply(stats::setNames(assayNames, assayNames), function(assayName) list(
        original = unname(counts$before[[assayName]]), resolved = unname(counts$after[[assayName]]),
        removed = unname(counts$before[[assayName]]) - unname(counts$after[[assayName]])))

    qcPlot <- tryCatch({
        arguments <- if (is.null(managedResult)) list(theObject = readLipidWorkflowState(workflowData$state_manager, resolved$reference)) else {
            metrics <- readWorkflowOmicsFilteringMetrics(workflowData$state_manager, managedResult)
            if (is.null(metrics)) stop("the worker computed no filtering metrics")
            list(theObject = NULL, metrics = metrics)
        }
        do.call(updateFilteringFn, c(arguments, list(
            step_name = "3_Duplicates_Resolved"
            , omics_type = omicType
            , return_grid = TRUE
            , overwrite = TRUE
        )))
    }, error = function(e) {
        logWarnFn(paste("Could not generate QC plot:", e$message))
        NULL
    })

    totalRemoved <- sum(vapply(statsList, function(stat) stat$removed, numeric(1)))
    resultText <- paste(
        "Duplicate Resolution Complete"
        , "=============================="
        , "Strategy: Keep feature with highest average intensity"
        , ""
        , "Per-Assay Results:"
        , paste(vapply(names(statsList), function(assayName) {
            stat <- statsList[[assayName]]
            sprintf(
                "  %s: %d -> %d rows (removed %d duplicates)"
                , assayName
                , stat$original
                , stat$resolved
                , stat$removed
            )
        }, character(1)), collapse = "\n")
        , ""
        , sprintf("Total duplicate rows removed: %d", totalRemoved)
        , ""
        , "State saved as: 'lipid_duplicates_resolved'"
        , sep = "\n"
    )

    list(
        stateReference = resolved$reference
        , statsList = statsList
        , qcPlot = qcPlot
        , totalRemoved = totalRemoved
        , resultText = resultText
    )
}

applyLipidDuplicateDetectionResult <- function(
    detectionResult,
    duplicateInfoVal,
    showNotificationFn = shiny::showNotification
) {
    duplicateInfoVal(detectionResult$duplicatesList)
    showNotificationFn(
        detectionResult$notificationMessage,
        type = detectionResult$notificationType
    )

    list(
        totalDuplicates = detectionResult$totalDuplicates
        , notificationMessage = detectionResult$notificationMessage
        , notificationType = detectionResult$notificationType
    )
}

registerLipidDuplicateDetectObserver <- function(
    input,
    workflowData,
    duplicateInfoVal,
    observeEventFn = shiny::observeEvent,
    reqFn = shiny::req,
    handleDetectionFn = handleLipidDuplicateDetection,
    applyDetectionResultFn = applyLipidDuplicateDetectionResult,
    logErrorFn = logger::log_error,
    showNotificationFn = shiny::showNotification
) {
    observeEventFn(input$detect_duplicates, {
        reqFn(workflowData$state_manager)

        tryCatch({
            detectionResult <- handleDetectionFn(
                workflowData = workflowData
            )

            applyDetectionResultFn(
                detectionResult = detectionResult
                , duplicateInfoVal = duplicateInfoVal
                , showNotificationFn = showNotificationFn
            )
        }, error = function(e) {
            msg <- paste("Error detecting duplicates:", e$message)
            logErrorFn(msg)
            showNotificationFn(msg, type = "error")
        })
    })

    invisible(input)
}

applyLipidDuplicateResolutionResult <- function(
    resolutionResult,
    resolutionStatsVal,
    duplicateInfoVal,
    filterPlotVal,
    output,
    renderTextFn = shiny::renderText,
    logInfoFn = logger::log_info,
    removeNotificationFn = shiny::removeNotification,
    showNotificationFn = shiny::showNotification,
    workingNotificationId = "lipid_dup_resolve_working"
) {
    resolutionStatsVal(resolutionResult$statsList)
    filterPlotVal(resolutionResult$qcPlot)
    output$resolution_results <- renderTextFn(resolutionResult$resultText)
    duplicateInfoVal(NULL)

    successMessage <- sprintf(
        "Duplicates resolved: %d rows removed",
        resolutionResult$totalRemoved
    )

    logInfoFn(paste(
        "Resolved duplicates: removed"
        , resolutionResult$totalRemoved
        , "rows"
    ))
    removeNotificationFn(workingNotificationId)
    showNotificationFn(successMessage, type = "message")

    list(
        successMessage = successMessage
        , notificationType = "message"
        , workingNotificationId = workingNotificationId
    )
}

registerLipidDuplicateResolveObserver <- function(
    input,
    workflowData,
    omicType,
    resolutionStatsVal,
    duplicateInfoVal,
    filterPlotVal,
    output,
    observeEventFn = shiny::observeEvent,
    reqFn = shiny::req,
    showNotificationFn = shiny::showNotification,
    handleResolutionFn = handleLipidDuplicateResolution,
    applyResolutionResultFn = applyLipidDuplicateResolutionResult,
    logErrorFn = logger::log_error,
    removeNotificationFn = shiny::removeNotification,
    workingNotificationId = "lipid_dup_resolve_working"
) {
    observeEventFn(input$resolve_duplicates, {
        reqFn(workflowData$state_manager)

        showNotificationFn(
            "Resolving duplicate features..."
            , id = workingNotificationId
            , duration = NULL
        )

        tryCatch({
            apply_result <- function(resolutionResult) applyResolutionResultFn(
                resolutionResult = resolutionResult
                , resolutionStatsVal = resolutionStatsVal
                , duplicateInfoVal = duplicateInfoVal
                , filterPlotVal = filterPlotVal
                , output = output
                , removeNotificationFn = removeNotificationFn
                , showNotificationFn = showNotificationFn
                , workingNotificationId = workingNotificationId
            )
            # KIT-078: on the worker under a managed session.
            runner <- workflowOmicsManagedActions("lipidomics")
            if (!is.null(runner)) {
                runner$run(list(omicsManagedCommit(omicsDuplicateResolutionStep(lipidomicsDeclarationBinding), workflowData$config_list, function(result) {
                    resolutionResult <- handleResolutionFn(workflowData = workflowData, omicType = omicType, managedResult = result)
                    apply_result(resolutionResult)
                })), failed = function(error) removeNotificationFn(workingNotificationId))
                return(invisible(NULL))
            }
            resolutionResult <- handleResolutionFn(
                workflowData = workflowData
                , omicType = omicType
            )
            apply_result(resolutionResult)
        }, error = function(e) {
            msg <- paste("Error resolving duplicates:", e$message)
            logErrorFn(msg)
            showNotificationFn(msg, type = "error", duration = 15)
            removeNotificationFn(workingNotificationId)
        })
    })

    invisible(input)
}

applyLipidDuplicateRevertResult <- function(
    revertResult,
    resolutionStatsVal,
    duplicateInfoVal,
    filterPlotVal,
    output,
    renderTextFn = shiny::renderText,
    showNotificationFn = shiny::showNotification
) {
    output$resolution_results <- renderTextFn(revertResult$resultText)
    resolutionStatsVal(NULL)
    duplicateInfoVal(NULL)
    filterPlotVal(NULL)
    showNotificationFn(
        revertResult$notificationMessage,
        type = revertResult$notificationType
    )

    list(
        notificationMessage = revertResult$notificationMessage
        , notificationType = revertResult$notificationType
    )
}

registerLipidDuplicateRevertObserver <- function(
    input,
    workflowData,
    resolutionStatsVal,
    duplicateInfoVal,
    filterPlotVal,
    output,
    observeEventFn = shiny::observeEvent,
    handleRevertFn = handleLipidDuplicateRevert,
    applyRevertResultFn = applyLipidDuplicateRevertResult,
    logErrorFn = logger::log_error,
    showNotificationFn = shiny::showNotification
) {
    observeEventFn(input$revert_duplicates, {
        tryCatch({
            revertResult <- handleRevertFn(
                workflowData = workflowData
            )

            applyRevertResultFn(
                revertResult = revertResult
                , resolutionStatsVal = resolutionStatsVal
                , duplicateInfoVal = duplicateInfoVal
                , filterPlotVal = filterPlotVal
                , output = output
                , showNotificationFn = showNotificationFn
            )
        }, error = function(e) {
            msg <- paste("Error reverting:", e$message)
            logErrorFn(msg)
            showNotificationFn(msg, type = "error")
        })
    })

    invisible(input)
}

handleLipidDuplicateRevert <- function(
    workflowData,
    reqFn = shiny::req,
    logInfoFn = logger::log_info
) {
    reqFn(workflowData$state_manager)

    history <- workflowData$state_manager$getHistory()
    if (length(history) <= 1) {
        stop("No previous state to revert to.")
    }

    previousStateName <- history[[length(history) - 1]]
    selectWorkflowStateReference(workflowData$state_manager, previousStateName) # KIT-075: no materialised object

    logInfoFn(paste("Reverted duplicate resolution to", previousStateName))

    list(
        previousStateName = previousStateName
        , resultText = paste("Reverted to previous state:", previousStateName)
        , notificationMessage = "Reverted successfully"
        , notificationType = "message"
    )
}

registerLipidDuplicateServerBindings <- function(
    input,
    output,
    workflowData,
    omicType,
    ns,
    duplicateInfoVal,
    resolutionStatsVal,
    filterPlotVal,
    registerDetectObserverFn = registerLipidDuplicateDetectObserver,
    registerSummaryOutputFn = registerLipidDuplicateSummaryOutput,
    registerTablesOutputFn = registerLipidDuplicateTablesOutput,
    registerTableObserverFn = registerLipidDuplicateTableObserver,
    registerResolveObserverFn = registerLipidDuplicateResolveObserver,
    registerRevertObserverFn = registerLipidDuplicateRevertObserver,
    registerFilterPlotOutputFn = registerLipidDuplicateFilterPlotOutput
) {
    registerDetectObserverFn(
        input = input
        , workflowData = workflowData
        , duplicateInfoVal = duplicateInfoVal
    )

    registerSummaryOutputFn(
        output = output
        , duplicateInfoVal = duplicateInfoVal
    )

    registerTablesOutputFn(
        output = output
        , duplicateInfoVal = duplicateInfoVal
        , ns = ns
    )

    registerTableObserverFn(
        output = output
        , duplicateInfoVal = duplicateInfoVal
    )

    registerResolveObserverFn(
        input = input
        , workflowData = workflowData
        , omicType = omicType
        , resolutionStatsVal = resolutionStatsVal
        , duplicateInfoVal = duplicateInfoVal
        , filterPlotVal = filterPlotVal
        , output = output
    )

    registerRevertObserverFn(
        input = input
        , workflowData = workflowData
        , resolutionStatsVal = resolutionStatsVal
        , duplicateInfoVal = duplicateInfoVal
        , filterPlotVal = filterPlotVal
        , output = output
    )

    registerFilterPlotOutputFn(
        output = output
        , filterPlotVal = filterPlotVal
    )

    invisible(output)
}

runLipidDuplicateModuleServerShell <- function(
    input,
    output,
    session,
    workflowData,
    omicType,
    initializeServerStateFn = initializeLipidDuplicateServerState,
    registerServerBindingsFn = registerLipidDuplicateServerBindings
) {
    ns <- session$ns
    duplicateState <- initializeServerStateFn()

    registerServerBindingsFn(
        input = input
        , output = output
        , workflowData = workflowData
        , omicType = omicType
        , ns = ns
        , duplicateInfoVal = duplicateState$duplicateInfo
        , resolutionStatsVal = duplicateState$resolutionStats
        , filterPlotVal = duplicateState$filterPlot
    )

    invisible(duplicateState)
}

