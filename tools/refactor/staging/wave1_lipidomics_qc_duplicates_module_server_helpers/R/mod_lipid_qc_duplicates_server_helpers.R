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

    currentS4 <- workflowData$state_manager$getState()
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
    logWarnFn = logger::log_warn
) {
    reqFn(workflowData$state_manager)

    currentS4 <- workflowData$state_manager$getState()
    reqFn(currentS4)

    if (!inherits(currentS4, "LipidomicsAssayData")) {
        stop("Current state is not a LipidomicsAssayData object")
    }

    lipidIdCol <- currentS4@lipid_id_column
    assayList <- currentS4@lipid_data
    assayNames <- names(assayList)
    if (is.null(assayNames)) {
        assayNames <- paste0("Assay_", seq_along(assayList))
    }

    resolvedAssayList <- vector("list", length(assayList))
    statsList <- vector("list", length(assayList))
    names(resolvedAssayList) <- assayNames
    names(statsList) <- assayNames

    for (i in seq_along(assayList)) {
        assayData <- assayList[[i]]
        assayName <- assayNames[[i]]
        originalRows <- nrow(assayData)
        sampleCols <- names(assayData)[vapply(assayData, is.numeric, logical(1))]

        if (length(sampleCols) == 0) {
            logWarnFn(paste("No numeric columns found in assay:", assayName))
            resolvedAssay <- assayData
        } else {
            resolvedAssay <- resolveDuplicatesFn(
                assay_tibble = assayData
                , id_col = lipidIdCol
                , sample_cols = sampleCols
            )
        }

        resolvedRows <- nrow(resolvedAssay)
        resolvedAssayList[[assayName]] <- resolvedAssay
        statsList[[assayName]] <- list(
            original = originalRows
            , resolved = resolvedRows
            , removed = originalRows - resolvedRows
        )
    }

    currentS4@lipid_data <- resolvedAssayList

    workflowData$state_manager$saveState(
        state_name = "lipid_duplicates_resolved"
        , s4_data_object = currentS4
        , config_object = workflowData$config_list
        , description = "Resolved duplicate lipid features by keeping highest intensity"
    )

    qcPlot <- tryCatch({
        updateFilteringFn(
            theObject = currentS4
            , step_name = "3_Duplicates_Resolved"
            , omics_type = omicType
            , return_grid = TRUE
            , overwrite = TRUE
        )
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
        currentS4 = currentS4
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
            resolutionResult <- handleResolutionFn(
                workflowData = workflowData
                , omicType = omicType
            )

            applyResolutionResultFn(
                resolutionResult = resolutionResult
                , resolutionStatsVal = resolutionStatsVal
                , duplicateInfoVal = duplicateInfoVal
                , filterPlotVal = filterPlotVal
                , output = output
                , removeNotificationFn = removeNotificationFn
                , showNotificationFn = showNotificationFn
                , workingNotificationId = workingNotificationId
            )
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
    workflowData$state_manager$revertToState(previousStateName)

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

