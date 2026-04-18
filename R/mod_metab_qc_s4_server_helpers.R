# Keep the summary table assembly isolated so the finalization wrapper can
# shed render logic without changing the public server entry point.
buildMetabQcS4DataSummaryUi <- function(
    currentS4,
    paragraphFn = shiny::p,
    iconFn = shiny::icon,
    tagListFn = shiny::tagList,
    tableTagFn = shiny::tags$table,
    tbodyTagFn = shiny::tags$tbody,
    trTagFn = shiny::tags$tr,
    tdTagFn = shiny::tags$td,
    strongFn = shiny::strong,
    formatFn = format,
    uniqueFn = unique,
    isNumericFn = is.numeric
) {
    if (is.null(currentS4) || !inherits(currentS4, "MetaboliteAssayData")) {
        return(paragraphFn(
            iconFn("exclamation-triangle", style = "color: orange;")
            , " No MetaboliteAssayData object available."
        ))
    }

    assay_list <- currentS4@metabolite_data
    n_assays <- length(assay_list)

    total_metabolites <- sum(vapply(assay_list, function(assay) {
        if (currentS4@metabolite_id_column %in% names(assay)) {
            length(uniqueFn(assay[[currentS4@metabolite_id_column]]))
        } else {
            nrow(assay)
        }
    }, numeric(1)))

    sample_cols <- if (n_assays > 0) {
        names(assay_list[[1]])[vapply(assay_list[[1]], isNumericFn, logical(1))]
    } else {
        character(0)
    }
    n_samples <- length(sample_cols)

    dm <- currentS4@design_matrix
    n_groups <- if (nrow(dm) > 0 && currentS4@group_id %in% names(dm)) {
        length(uniqueFn(dm[[currentS4@group_id]]))
    } else {
        NA_integer_
    }

    tagListFn(
        tableTagFn(
            class = "table table-condensed"
            , style = "margin-bottom: 0;"
            , tbodyTagFn(
                trTagFn(
                    tdTagFn(strongFn("Number of Assays:"))
                    , tdTagFn(n_assays)
                )
                , trTagFn(
                    tdTagFn(strongFn("Total Metabolites:"))
                    , tdTagFn(formatFn(total_metabolites, big.mark = ","))
                )
                , trTagFn(
                    tdTagFn(strongFn("Number of Samples:"))
                    , tdTagFn(n_samples)
                )
                , trTagFn(
                    tdTagFn(strongFn("Experimental Groups:"))
                    , tdTagFn(if (!is.na(n_groups)) n_groups else "N/A")
                )
                , trTagFn(
                    tdTagFn(strongFn("Metabolite ID Column:"))
                    , tdTagFn(currentS4@metabolite_id_column)
                )
                , trTagFn(
                    tdTagFn(strongFn("Sample ID Column:"))
                    , tdTagFn(currentS4@sample_id)
                )
            )
        )
    )
}

# Isolate the state-history markup so the server wrapper can hand off another
# render branch without changing the public module entry point.
buildMetabQcS4StateHistoryUi <- function(
    history,
    paragraphFn = shiny::p,
    iconFn = shiny::icon,
    orderedListTagFn = shiny::tags$ol,
    listItemTagFn = shiny::tags$li,
    spanTagFn = shiny::tags$span,
    sprintfFn = sprintf
) {
    if (length(history) == 0) {
        return(paragraphFn(
            iconFn("info-circle")
            , " No processing history available."
            , style = "color: #666;"
        ))
    }

    history_items <- lapply(seq_along(history), function(i) {
        state_name <- history[i]
        is_current <- i == length(history)

        icon_name <- if (is_current) "arrow-right" else "check"
        icon_color <- if (is_current) "blue" else "green"
        text_style <- if (is_current) "font-weight: bold;" else ""

        item_children <- list(
            iconFn(icon_name, style = paste("color:", icon_color))
            , spanTagFn(
                sprintfFn(" %d. %s", i, state_name)
                , style = text_style
            )
        )

        if (is_current) {
            item_children <- c(
                item_children,
                list(spanTagFn(" (current)", style = "color: blue;"))
            )
        }

        do.call(listItemTagFn, item_children)
    })

    do.call(
        orderedListTagFn,
        c(history_items, list(style = "list-style: none; padding-left: 0;"))
    )
}

# Keep the state-history render-path history-fetch shell isolated so the
# finalization wrapper can hand off history retrieval without changing the
# public module entry point.
getMetabQcS4StateHistory <- function(
    stateManager,
    historyGetterFn = function(manager) manager$getHistory(),
    errorHistory = character(0)
) {
    tryCatch(
        historyGetterFn(stateManager),
        error = function(e) errorHistory
    )
}

# Keep the state-history render-path sequencing isolated so the finalization
# wrapper can hand off the remaining UI render branch without changing the
# public module entry point.
buildMetabQcS4StateHistoryRenderOutput <- function(
    stateManager,
    getStateHistoryFn = getMetabQcS4StateHistory,
    buildStateHistoryUiFn = buildMetabQcS4StateHistoryUi
) {
    history <- getStateHistoryFn(stateManager = stateManager)

    buildStateHistoryUiFn(history)
}

# Keep the per-assay statistics datatable assembly isolated so the finalization
# wrapper can hand off another render branch without changing the public server
# entry point.
buildMetabQcS4AssayStatsDatatable <- function(
    currentS4,
    datatableFn = DT::datatable,
    datatableOptions = list(
        dom = "t",
        paging = FALSE,
        ordering = FALSE
    ),
    uniqueFn = unique,
    isNumericFn = is.numeric,
    sumFn = sum,
    isNaFn = is.na,
    roundFn = round
) {
    if (is.null(currentS4) || !inherits(currentS4, "MetaboliteAssayData")) {
        return(NULL)
    }

    assay_list <- currentS4@metabolite_data
    assay_names <- names(assay_list)

    if (is.null(assay_names)) {
        assay_names <- character(0)
    }

    metabolite_id_col <- currentS4@metabolite_id_column

    stats_df <- data.frame(
        Assay = assay_names
        , Metabolites = vapply(assay_list, function(assay) {
            if (metabolite_id_col %in% names(assay)) {
                length(uniqueFn(assay[[metabolite_id_col]]))
            } else {
                nrow(assay)
            }
        }, numeric(1))
        , Samples = vapply(assay_list, function(assay) {
            sumFn(vapply(assay, isNumericFn, logical(1)))
        }, numeric(1))
        , stringsAsFactors = FALSE
    )

    stats_df$Missingness <- vapply(assay_list, function(assay) {
        quant_cols <- names(assay)[vapply(assay, isNumericFn, logical(1))]
        if (length(quant_cols) == 0) {
            return(NA_real_)
        }

        total_cells <- nrow(assay) * length(quant_cols)
        missing_cells <- sumFn(vapply(quant_cols, function(col) {
            sumFn(isNaFn(assay[[col]]) | assay[[col]] == 0)
        }, numeric(1)))

        roundFn((missing_cells / total_cells) * 100, 1)
    }, numeric(1))

    datatableFn(
        stats_df
        , options = datatableOptions
        , rownames = FALSE
        , class = "compact stripe"
    )
}

# Keep the QC-progress plot dispatch isolated so the finalization wrapper can
# hand off the last remaining simple render branch without changing the public
# module entry point.
renderMetabQcS4FilterPlot <- function(
    filterPlot,
    reqFn = shiny::req,
    inheritsFn = inherits,
    gridDrawFn = grid::grid.draw,
    printFn = print
) {
    plotObject <- filterPlot()
    reqFn(plotObject)

    if (inheritsFn(plotObject, "grob") || inheritsFn(plotObject, "gtable")) {
        gridDrawFn(plotObject)
    } else if (inheritsFn(plotObject, "ggplot")) {
        printFn(plotObject)
    }

    invisible(NULL)
}

# Keep the filter-plot render-path sequencing isolated so the finalization
# wrapper can hand off the remaining plot render branch without changing the
# public module entry point.
buildMetabQcS4FilterPlotRenderOutput <- function(
    filterPlot,
    renderFilterPlotFn = renderMetabQcS4FilterPlot
) {
    renderFilterPlotFn(filterPlot = filterPlot)
}

# Keep the finalize-results text assembly isolated so the finalization wrapper
# can hand off observer reporting without changing the public module entry
# point.
buildMetabQcS4FinalizeResultsText <- function(
    currentS4,
    history,
    uniqueFn = unique,
    sumFn = sum,
    pasteFn = paste,
    sprintfFn = sprintf
) {
    assayList <- currentS4@metabolite_data
    metaboliteIdColumn <- currentS4@metabolite_id_column

    totalMetabolites <- sumFn(vapply(assayList, function(assay) {
        if (metaboliteIdColumn %in% names(assay)) {
            length(uniqueFn(assay[[metaboliteIdColumn]]))
        } else {
            nrow(assay)
        }
    }, numeric(1)))

    historyText <- pasteFn(
        sprintfFn("  %d. %s", seq_along(history), history)
        , collapse = "\n"
    )

    pasteFn(
        "QC Finalization Complete"
        , "========================"
        , ""
        , sprintfFn("Total metabolites retained: %d", totalMetabolites)
        , sprintfFn("Number of assays: %d", length(assayList))
        , sprintfFn("Processing steps completed: %d", length(history))
        , ""
        , "Processing History:"
        , historyText
        , ""
        , "State saved as: 'metab_qc_complete'"
        , ""
        , "You can now proceed to the Normalization tab."
        , sep = "\n"
    )
}

# Keep the finalize observer current-state fetch shell isolated so the
# finalization wrapper can hand off state retrieval without changing the public
# module entry point.
getMetabQcS4FinalizeState <- function(
    stateManager,
    stateGetterFn = function(manager) manager$getState()
) {
    stateGetterFn(stateManager)
}

# Keep the finalize observer current-state validation shell isolated so the
# finalization wrapper can hand off `shiny::req()` plus class validation
# without changing the public module entry point.
validateMetabQcS4FinalizeState <- function(
    currentS4,
    reqFn = shiny::req,
    inheritsFn = inherits,
    expectedClass = "MetaboliteAssayData",
    errorMessage = paste("Current state is not a", expectedClass, "object")
) {
    reqFn(currentS4)

    if (!inheritsFn(currentS4, expectedClass)) {
        stop(errorMessage)
    }

    currentS4
}

# Keep the finalization state-save shell isolated so the finalization wrapper
# can hand off persistence without changing the public module entry point.
saveMetabQcS4CompletedState <- function(
    stateManager,
    currentS4,
    configObject,
    stateName = "metab_qc_complete",
    description = "QC processing complete - ready for normalization"
) {
    stateManager$saveState(
        state_name = stateName
        , s4_data_object = currentS4
        , config_object = configObject
        , description = description
    )

    invisible(stateName)
}

# Keep the finalize observer tab-status completion shell isolated so the
# finalization wrapper can hand off reactive status replacement without
# changing the public module entry point.
completeMetabQcS4TabStatus <- function(
    workflowData,
    tabName = "quality_control",
    status = "complete"
) {
    updatedStatus <- workflowData$tab_status
    updatedStatus[[tabName]] <- status
    workflowData$tab_status <- updatedStatus

    invisible(updatedStatus)
}

# Keep the finalize observer history-fetch shell isolated so the finalization
# wrapper can hand off state-history retrieval without changing the public
# module entry point.
getMetabQcS4FinalizeHistory <- function(
    stateManager,
    historyGetterFn = function(manager) manager$getHistory()
) {
    historyGetterFn(stateManager)
}

# Keep the finalize observer QC-tracking plot-refresh shell isolated so the
# finalization wrapper can hand off plot regeneration without changing the
# public module entry point.
updateMetabQcS4TrackingPlot <- function(
    currentS4,
    omicType,
    setFilterPlotFn,
    stepName = "4_QC_Complete",
    updateMetaboliteFilteringFn = updateMetaboliteFiltering
) {
    qcPlot <- tryCatch({
        updateMetaboliteFilteringFn(
            theObject = currentS4
            , step_name = stepName
            , omics_type = omicType
            , return_grid = TRUE
            , overwrite = TRUE
        )
    }, error = function(e) {
        NULL
    })

    setFilterPlotFn(qcPlot)

    invisible(qcPlot)
}

# Keep the finalize observer success-reporting shell isolated so the
# finalization wrapper can hand off render/log/notification work without
# changing the public module entry point.
reportMetabQcS4FinalizeSuccess <- function(
    currentS4,
    history,
    output,
    finalizeResultsOutputName = "finalize_results",
    buildResultsTextFn = buildMetabQcS4FinalizeResultsText,
    renderTextFn = shiny::renderText,
    logInfoFn = logger::log_info,
    showNotificationFn = shiny::showNotification
) {
    resultText <- buildResultsTextFn(
        currentS4 = currentS4
        , history = history
    )

    output[[finalizeResultsOutputName]] <- renderTextFn(resultText)

    logInfoFn("Metabolomics QC finalized successfully")
    showNotificationFn(
        "QC complete! Proceed to Normalization."
        , type = "message"
        , duration = 5
    )

    invisible(resultText)
}

# Keep the finalize observer error-reporting shell isolated so the
# finalization wrapper can hand off failure reporting without changing the
# public module entry point.
reportMetabQcS4FinalizeError <- function(
    error,
    messagePrefix = "Error finalizing QC:",
    pasteFn = paste,
    logErrorFn = logger::log_error,
    showNotificationFn = shiny::showNotification,
    notificationType = "error"
) {
    messageText <- pasteFn(messagePrefix, error$message)

    logErrorFn(messageText)
    showNotificationFn(messageText, type = notificationType)

    invisible(messageText)
}

# Keep the remaining finalize observer orchestration shell isolated so the
# finalization wrapper can hand off the `tryCatch()` envelope and existing
# helper chain without changing the public module entry point.
runMetabQcS4FinalizeWorkflow <- function(
    workflowData,
    omicType,
    filterPlot,
    output,
    getFinalizeStateFn = getMetabQcS4FinalizeState,
    validateFinalizeStateFn = validateMetabQcS4FinalizeState,
    saveCompletedStateFn = saveMetabQcS4CompletedState,
    updateTrackingPlotFn = updateMetabQcS4TrackingPlot,
    completeTabStatusFn = completeMetabQcS4TabStatus,
    getFinalizeHistoryFn = getMetabQcS4FinalizeHistory,
    reportFinalizeSuccessFn = reportMetabQcS4FinalizeSuccess,
    reportFinalizeErrorFn = reportMetabQcS4FinalizeError
) {
    tryCatch({
        currentS4 <- getFinalizeStateFn(
            stateManager = workflowData$state_manager
        )
        currentS4 <- validateFinalizeStateFn(currentS4 = currentS4)

        saveCompletedStateFn(
            stateManager = workflowData$state_manager
            , currentS4 = currentS4
            , configObject = workflowData$config_list
        )

        updateTrackingPlotFn(
            currentS4 = currentS4
            , omicType = omicType
            , setFilterPlotFn = filterPlot
        )

        completeTabStatusFn(workflowData = workflowData)

        history <- getFinalizeHistoryFn(
            stateManager = workflowData$state_manager
        )

        reportFinalizeSuccessFn(
            currentS4 = currentS4
            , history = history
            , output = output
        )
    }, error = function(e) {
        reportFinalizeErrorFn(error = e)
    })
}

# Keep the data-summary render-path current-state fetch shell isolated so the
# finalization wrapper can hand off current-state retrieval without changing
# the public module entry point.
getMetabQcS4DataSummaryState <- function(
    stateManager,
    stateGetterFn = function(manager) manager$getState(),
    errorState = NULL
) {
    tryCatch(
        stateGetterFn(stateManager),
        error = function(e) errorState
    )
}

# Keep the data-summary render-path sequencing isolated so the finalization
# wrapper can hand off the remaining UI render branch without changing the
# public module entry point.
buildMetabQcS4DataSummaryRenderOutput <- function(
    stateManager,
    getDataSummaryStateFn = getMetabQcS4DataSummaryState,
    buildDataSummaryUiFn = buildMetabQcS4DataSummaryUi
) {
    currentS4 <- getDataSummaryStateFn(stateManager = stateManager)

    buildDataSummaryUiFn(currentS4)
}

# Keep the assay-stats render-path current-state fetch shell isolated so the
# finalization wrapper can hand off current-state retrieval without changing
# the public module entry point.
getMetabQcS4AssayStatsState <- function(
    stateManager,
    stateGetterFn = function(manager) manager$getState(),
    errorState = NULL
) {
    tryCatch(
        stateGetterFn(stateManager),
        error = function(e) errorState
    )
}

# Keep the assay-stats render-path sequencing isolated so the finalization
# wrapper can hand off the remaining DT render branch without changing the
# public module entry point.
buildMetabQcS4AssayStatsRenderOutput <- function(
    stateManager,
    getAssayStatsStateFn = getMetabQcS4AssayStatsState,
    buildAssayStatsDatatableFn = buildMetabQcS4AssayStatsDatatable
) {
    currentS4 <- getAssayStatsStateFn(stateManager = stateManager)

    buildAssayStatsDatatableFn(currentS4)
}

# Keep the module server body isolated so the public wrapper can shed the
# observer and render wiring without changing the module entrypoint.
runMetabQcS4ServerBody <- function(
    input,
    output,
    session,
    workflowData,
    omicType = NULL,
    experimentLabel = NULL,
    reactiveValFn = shiny::reactiveVal,
    renderUiFn = shiny::renderUI,
    renderDtFn = DT::renderDT,
    observeEventFn = shiny::observeEvent,
    renderPlotFn = shiny::renderPlot,
    reqFn = shiny::req,
    buildStateHistoryRenderOutputFn = buildMetabQcS4StateHistoryRenderOutput,
    buildDataSummaryRenderOutputFn = buildMetabQcS4DataSummaryRenderOutput,
    buildAssayStatsRenderOutputFn = buildMetabQcS4AssayStatsRenderOutput,
    runFinalizeWorkflowFn = runMetabQcS4FinalizeWorkflow,
    buildFilterPlotRenderOutputFn = buildMetabQcS4FilterPlotRenderOutput
) {
    filterPlot <- reactiveValFn(NULL)

    # Render state history
    output$state_history <- renderUiFn({
        reqFn(workflowData$state_manager)

        buildStateHistoryRenderOutputFn(
            stateManager = workflowData$state_manager
        )
    })

    # Render data summary
    output$data_summary <- renderUiFn({
        reqFn(workflowData$state_manager)

        buildDataSummaryRenderOutputFn(
            stateManager = workflowData$state_manager
        )
    })

    # Render per-assay statistics table
    output$assay_stats_table <- renderDtFn({
        reqFn(workflowData$state_manager)

        buildAssayStatsRenderOutputFn(
            stateManager = workflowData$state_manager
        )
    })

    # Finalize QC
    observeEventFn(input$finalize_qc, {
        reqFn(workflowData$state_manager)

        runFinalizeWorkflowFn(
            workflowData = workflowData
            , omicType = omicType
            , filterPlot = filterPlot
            , output = output
        )
    })

    # Render QC progress plot
    output$filter_plot <- renderPlotFn({
        buildFilterPlotRenderOutputFn(
            filterPlot = filterPlot
        )
    })

    invisible(NULL)
}

