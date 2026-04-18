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

