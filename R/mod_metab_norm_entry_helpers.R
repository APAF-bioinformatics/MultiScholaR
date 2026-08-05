createMetabNormReactiveState <- function(
    createReactiveValuesFn = shiny::reactiveValues
) {
    createReactiveValuesFn(
        assay_names = NULL,
        itsd_selections = list(),
        pre_norm_qc_generated = FALSE,
        normalization_complete = FALSE,
        ruv_complete = FALSE,
        correlation_filtering_complete = FALSE,
        qc_plot_paths = list(),
        post_filter_obj = NULL,
        post_itsd_obj = NULL,
        post_log2_obj = NULL,
        post_norm_obj = NULL,
        ruv_corrected_obj = NULL,
        correlation_filtered_obj = NULL,
        ruv_optimization_results = list(),
        correlation_results = list(),
        plot_refresh_trigger = 0,
        normalization_log = character(0)
    )
}

runMetabNormModuleServerShell <- function(
    input,
    output,
    session,
    id,
    workflowData,
    experimentPaths,
    omicType,
    experimentLabel,
    selectedTab = NULL,
    logInfoFn = logger::log_info,
    createReactiveStateFn = createMetabNormReactiveState,
    appendNormalizationLogFn = appendMetabNormNormalizationLog,
    observeFn = shiny::observe,
    observeEventFn = shiny::observeEvent,
    initializeAssayNamesFn = initializeMetabNormAssayNames,
    runAutoPreNormalizationQcObserverShellFn = runMetabNormAutoPreNormalizationQcObserverShell,
    updateDesignDrivenChoicesFn = updateMetabNormDesignDrivenChoices,
    renderNormalizationLogFn = renderMetabNormNormalizationLog,
    renderItsdSelectionUiFn = renderMetabNormItsdSelectionUi,
    renderRuvQcUiFn = renderMetabNormRuvQcUi,
    runQcImageBindingShellFn = runMetabNormQcImageBindingShell,
    runAssayLabelBindingShellFn = runMetabNormAssayLabelBindingShell,
    runItsdSelectionTableObserverShellFn = runMetabNormItsdSelectionTableObserverShell,
    runItsdSelectionTrackingObserverShellFn = runMetabNormItsdSelectionTrackingObserverShell,
    runNormalizationObserverWrapperFn = runMetabNormNormalizationObserverWrapper,
    runResetNormalizationObserverWrapperFn = runMetabNormResetNormalizationObserverWrapper,
    runRuvBindingObserverShellFn = runMetabNormRuvBindingObserverShell,
    runApplyCorrelationObserverWrapperFn = runMetabNormApplyCorrelationObserverWrapper,
    runSkipCorrelationObserverWrapperFn = runMetabNormSkipCorrelationObserverWrapper,
    renderCorrelationFilterSummaryFn = renderMetabNormCorrelationFilterSummary,
    renderFinalQcPlotFn = renderMetabNormFinalQcPlot,
    runExportSessionObserverWrapperFn = runMetabNormExportSessionObserverWrapper,
    showNotificationFn = shiny::showNotification,
    removeNotificationFn = shiny::removeNotification,
    reqFn = shiny::req,
    withProgressFn = shiny::withProgress
) {
    ns <- session$ns

    logInfoFn("=== METABOLOMICS NORMALIZATION MODULE STARTED ===")
    logInfoFn(paste("Module ID:", id))

    normData <- createReactiveStateFn()

    addLog <- function(message) {
        appendNormalizationLogFn(
            normData = normData,
            message = message
        )
    }

    observeFn({
        initializeAssayNamesFn(
            stateManager = workflowData$state_manager,
            normData = normData
        )
    })

    if (!is.null(selectedTab)) {
        observeEventFn(selectedTab(), {
            runAutoPreNormalizationQcObserverShellFn(
                selectedTab = selectedTab(),
                workflowData = workflowData,
                experimentPaths = experimentPaths,
                normData = normData,
                colorVariable = input$color_variable,
                shapeVariable = input$shape_variable,
                addLogFn = addLog
            )
        }, ignoreInit = FALSE)
    }

    observeFn({
        updateDesignDrivenChoicesFn(
            session = session,
            designMatrix = workflowData$design_matrix
        )
    })

    output$norm_log <- renderNormalizationLogFn(
        normData = normData
    )

    output$itsd_selection_ui <- renderItsdSelectionUiFn(
        normData = normData,
        ns = ns
    )

    output$ruv_qc_ui <- renderRuvQcUiFn(
        normData = normData,
        ns = ns
    )

    runQcImageBindingShellFn(
        output = output,
        normData = normData,
        qcDir = experimentPaths$metabolite_qc_dir
    )

    runAssayLabelBindingShellFn(
        output = output,
        getAssayNamesFn = function() normData$assay_names
    )

    observeFn({
        runItsdSelectionTableObserverShellFn(
            normData = normData,
            workflowData = workflowData,
            output = output
        )
    })

    observeFn({
        runItsdSelectionTrackingObserverShellFn(
            normData = normData,
            input = input
        )
    })

    observeEventFn(input$run_normalization, {
        runNormalizationObserverWrapperFn(
            workflowData = workflowData,
            input = input,
            experimentPaths = experimentPaths,
            omicType = omicType,
            normData = normData,
            addLogFn = addLog,
            showNotificationFn = showNotificationFn,
            reqFn = reqFn,
            withProgressFn = withProgressFn
        )
    })

    observeEventFn(input$reset_normalization, {
        runResetNormalizationObserverWrapperFn(
            workflowData = workflowData,
            normData = normData,
            addLogFn = addLog,
            showNotificationFn = showNotificationFn,
            reqFn = reqFn
        )
    })

    observeFn({
        runRuvBindingObserverShellFn(
            normData = normData,
            output = output
        )
    })

    observeEventFn(input$apply_correlation_filter, {
        runApplyCorrelationObserverWrapperFn(
            workflowData = workflowData,
            input = input,
            normData = normData,
            addLogFn = addLog,
            showNotificationFn = showNotificationFn,
            removeNotificationFn = removeNotificationFn,
            reqFn = reqFn
        )
    })

    observeEventFn(input$skip_correlation_filter, {
        runSkipCorrelationObserverWrapperFn(
            workflowData = workflowData,
            normData = normData,
            addLogFn = addLog,
            showNotificationFn = showNotificationFn,
            reqFn = reqFn
        )
    })

    output$correlation_filter_summary <- renderCorrelationFilterSummaryFn(
        normData = normData
    )

    output$final_qc_plot <- renderFinalQcPlotFn(
        normData = normData,
        colorVariableFn = function() input$color_variable,
        shapeVariableFn = function() input$shape_variable
    )

    observeEventFn(input$export_session, {
        runExportSessionObserverWrapperFn(
            workflowData = workflowData,
            input = input,
            normData = normData,
            experimentPaths = experimentPaths,
            experimentLabel = experimentLabel,
            addLogFn = addLog,
            logInfoFn = logInfoFn,
            reqFn = reqFn
        )
    })

    invisible(normData)
}

runMetabNormModuleServerEntryShell <- function(
    id,
    workflowData,
    experimentPaths,
    omicType,
    experimentLabel,
    selectedTab = NULL,
    moduleServerFn = shiny::moduleServer,
    runModuleServerShellFn = runMetabNormModuleServerShell
) {
    moduleServerFn(id, function(input, output, session) {
        runModuleServerShellFn(
            input = input,
            output = output,
            session = session,
            id = id,
            workflowData = workflowData,
            experimentPaths = experimentPaths,
            omicType = omicType,
            experimentLabel = experimentLabel,
            selectedTab = selectedTab
        )
    })
}

runMetabNormModuleServerPublicWrapper <- function(
    id,
    workflow_data,
    experiment_paths,
    omic_type,
    experiment_label,
    selected_tab = NULL,
    runModuleServerEntryShellFn = runMetabNormModuleServerEntryShell
) {
    runModuleServerEntryShellFn(
        id = id,
        workflowData = workflow_data,
        experimentPaths = experiment_paths,
        omicType = omic_type,
        experimentLabel = experiment_label,
        selectedTab = selected_tab
    )
}

