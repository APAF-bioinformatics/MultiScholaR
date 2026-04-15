createLipidNormReactiveState <- function(
    reactiveValuesFn = shiny::reactiveValues
) {
    reactiveValuesFn(
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

createLipidNormStartupRuntime <- function(
    input,
    experimentPaths,
    normData,
    ns,
    buildAddLogFn = buildLipidNormAddLog,
    buildPlotAestheticsGetterFn = buildLipidNormPlotAestheticsGetter,
    buildCompositeFromFilesGeneratorFn = buildLipidNormCompositeFromFilesGenerator,
    buildQcImageRendererFn = buildLipidNormQcImageRenderer,
    buildLogRendererFn = buildLipidNormLogRenderer,
    buildItsdSelectionUiRendererFn = buildLipidNormItsdSelectionUiRenderer,
    buildRuvQcUiRendererFn = buildLipidNormRuvQcUiRenderer,
    buildAssayLabelRendererFn = buildLipidNormAssayLabelRenderer,
    buildCorrelationFilterSummaryRendererFn = buildLipidNormCorrelationFilterSummaryRenderer,
    buildFinalQcPlotRendererFn = buildLipidNormFinalQcPlotRenderer
) {
    addLog <- buildAddLogFn(
        normData = normData
    )

    getPlotAesthetics <- buildPlotAestheticsGetterFn(
        input = input
    )

    generateCompositeFromFiles <- buildCompositeFromFilesGeneratorFn()

    renderQcImageForAssay <- buildQcImageRendererFn(
        normData = normData,
        experimentPaths = experimentPaths
    )

    renderNormLog <- buildLogRendererFn(
        normData = normData
    )

    renderItsdSelectionUi <- buildItsdSelectionUiRendererFn(
        normData = normData,
        ns = ns
    )

    renderRuvQcUi <- buildRuvQcUiRendererFn(
        normData = normData,
        ns = ns
    )

    renderAssayLabel <- buildAssayLabelRendererFn(
        normData = normData
    )

    renderCorrelationFilterSummary <- buildCorrelationFilterSummaryRendererFn(
        normData = normData
    )

    renderFinalQcPlot <- buildFinalQcPlotRendererFn(
        normData = normData,
        getPlotAestheticsFn = getPlotAesthetics
    )

    list(
        addLog = addLog,
        getPlotAesthetics = getPlotAesthetics,
        generateCompositeFromFiles = generateCompositeFromFiles,
        renderQcImageForAssay = renderQcImageForAssay,
        renderNormLog = renderNormLog,
        renderItsdSelectionUi = renderItsdSelectionUi,
        renderRuvQcUi = renderRuvQcUi,
        renderAssayLabel = renderAssayLabel,
        renderCorrelationFilterSummary = renderCorrelationFilterSummary,
        renderFinalQcPlot = renderFinalQcPlot
    )
}

registerLipidNormPrimaryStartupOutputs <- function(
    output,
    startupRuntime,
    registerLogOutputFn = registerLipidNormLogOutput,
    registerItsdSelectionOutputFn = registerLipidNormItsdSelectionOutput,
    registerRuvQcOutputFn = registerLipidNormRuvQcOutput,
    registerStaticQcImageOutputsFn = registerLipidNormStaticQcImageOutputs,
    registerAssayLabelOutputsFn = registerLipidNormAssayLabelOutputs
) {
    registerLogOutputFn(output, startupRuntime$renderNormLog)
    registerItsdSelectionOutputFn(output, startupRuntime$renderItsdSelectionUi)
    registerRuvQcOutputFn(output, startupRuntime$renderRuvQcUi)
    registerStaticQcImageOutputsFn(output, startupRuntime$renderQcImageForAssay)
    registerAssayLabelOutputsFn(output, startupRuntime$renderAssayLabel)

    invisible(output)
}

registerLipidNormItsdSelectionRuntime <- function(
    input,
    output,
    workflowData,
    normData,
    registerItsdTableOutputsFn = registerLipidNormItsdTableOutputs,
    registerItsdSelectionTrackingFn = registerLipidNormItsdSelectionTracking
) {
    registerItsdTableOutputsFn(
        output = output,
        workflowData = workflowData,
        normData = normData
    )

    registerItsdSelectionTrackingFn(
        input = input,
        normData = normData
    )

    invisible(normData)
}

registerLipidNormPostNormalizationOutputs <- function(
    output,
    normData,
    startupRuntime,
    registerRuvCancorOutputsFn = registerLipidNormRuvCancorOutputs,
    registerCorrelationFilterSummaryOutputFn = registerLipidNormCorrelationFilterSummaryOutput,
    registerFinalQcPlotOutputFn = registerLipidNormFinalQcPlotOutput
) {
    registerRuvCancorOutputsFn(
        output = output,
        normData = normData
    )

    registerCorrelationFilterSummaryOutputFn(
        output = output,
        renderCorrelationFilterSummary = startupRuntime$renderCorrelationFilterSummary
    )

    registerFinalQcPlotOutputFn(
        output = output,
        renderFinalQcPlot = startupRuntime$renderFinalQcPlot
    )

    invisible(output)
}

registerLipidNormStartupObserverRuntime <- function(
    session,
    workflowData,
    experimentPaths,
    normData,
    addLog,
    getPlotAestheticsFn,
    selectedTab = NULL,
    registerAssayNameInitializationObserverFn = registerLipidNormAssayNameInitializationObserver,
    registerSelectedTabPreNormalizationObserverFn = registerLipidNormSelectedTabPreNormalizationObserver,
    registerDesignDrivenChoiceObserverFn = registerLipidNormDesignDrivenChoiceObserver
) {
    registerAssayNameInitializationObserverFn(
        workflowData = workflowData,
        normData = normData
    )

    if (!is.null(selectedTab)) {
        registerSelectedTabPreNormalizationObserverFn(
            selectedTab = selectedTab
            , workflowData = workflowData
            , experimentPaths = experimentPaths
            , normData = normData
            , addLog = addLog
            , getPlotAestheticsFn = getPlotAestheticsFn
        )
    }

    registerDesignDrivenChoiceObserverFn(
        session = session,
        workflowData = workflowData
    )

    invisible(session)
}

registerLipidNormServerRuntime <- function(
    input,
    output,
    session,
    workflowData,
    experimentPaths,
    omicType,
    experimentLabel,
    normData,
    startupRuntime,
    addLog,
    getPlotAestheticsFn,
    generateCompositeFromFilesFn,
    selectedTab = NULL,
    registerStartupObserverRuntimeFn = registerLipidNormStartupObserverRuntime,
    registerPrimaryStartupOutputsFn = registerLipidNormPrimaryStartupOutputs,
    registerItsdSelectionRuntimeFn = registerLipidNormItsdSelectionRuntime,
    registerRunNormalizationObserverFn = registerLipidNormRunNormalizationObserver,
    registerResetNormalizationObserverFn = registerLipidNormResetNormalizationObserver,
    registerPostNormalizationOutputsFn = registerLipidNormPostNormalizationOutputs,
    registerApplyCorrelationFilterObserverFn = registerLipidNormApplyCorrelationFilterObserver,
    registerSkipCorrelationFilterObserverFn = registerLipidNormSkipCorrelationFilterObserver,
    registerExportSessionObserverFn = registerLipidNormExportSessionObserver
) {
    registerStartupObserverRuntimeFn(
        session = session
        , workflowData = workflowData
        , experimentPaths = experimentPaths
        , normData = normData
        , addLog = addLog
        , getPlotAestheticsFn = getPlotAestheticsFn
        , selectedTab = selectedTab
    )

    registerPrimaryStartupOutputsFn(
        output = output
        , startupRuntime = startupRuntime
    )

    registerItsdSelectionRuntimeFn(
        input = input
        , output = output
        , workflowData = workflowData
        , normData = normData
    )

    registerRunNormalizationObserverFn(
        input = input
        , workflowData = workflowData
        , experimentPaths = experimentPaths
        , omicType = omicType
        , normData = normData
        , addLog = addLog
        , getPlotAestheticsFn = getPlotAestheticsFn
        , generateCompositeFromFilesFn = generateCompositeFromFilesFn
    )

    registerResetNormalizationObserverFn(
        input = input
        , workflowData = workflowData
        , normData = normData
        , addLog = addLog
    )

    registerPostNormalizationOutputsFn(
        output = output
        , normData = normData
        , startupRuntime = startupRuntime
    )

    registerApplyCorrelationFilterObserverFn(
        input = input
        , workflowData = workflowData
        , normData = normData
        , addLog = addLog
    )

    registerSkipCorrelationFilterObserverFn(
        input = input
        , workflowData = workflowData
        , normData = normData
        , addLog = addLog
    )

    registerExportSessionObserverFn(
        input = input
        , workflowData = workflowData
        , experimentPaths = experimentPaths
        , experimentLabel = experimentLabel
        , normData = normData
        , addLog = addLog
    )

    invisible(input)
}

