# Keep the main run-normalization observer wrapper top-level so a later
# extraction wave can move it without reopening the module server wrapper.
runMetabNormNormalizationObserverWrapper <- function(
    workflowData,
    input,
    experimentPaths,
    omicType,
    normData,
    addLogFn = function(entry) invisible(entry),
    showNotificationFn = shiny::showNotification,
    reqFn = shiny::req,
    withProgressFn = shiny::withProgress,
    getPlotAestheticsFn = getPlotAesthetics,
    runObserverShellFn = runMetabNormNormalizationObserverShell,
    runPipelineShellFn = runMetabNormNormalizationPipelineShell,
    generateCompositeFromFilesFn = generateMetabNormCompositeFromFiles,
    savePlotFn = savePlot,
    logWarnFn = logger::log_warn
) {
    runObserverShellFn(
        workflowData = workflowData
        , addLogFn = addLogFn
        , showNotificationFn = showNotificationFn
        , reqFn = reqFn
        , withProgressFn = withProgressFn
        , runPipelineFn = function() {
            runPipelineShellFn(
                workflowData = workflowData
                , inputValues = list(
                    apply_itsd = input$apply_itsd
                    , itsd_aggregation = input$itsd_aggregation
                    , log_offset = input$log_offset
                    , norm_method = input$norm_method
                    , ruv_mode = input$ruv_mode
                    , auto_percentage_min = input$auto_percentage_min
                    , auto_percentage_max = input$auto_percentage_max
                    , ruv_grouping_variable = input$ruv_grouping_variable
                    , separation_metric = input$separation_metric
                    , k_penalty_weight = input$k_penalty_weight
                    , adaptive_k_penalty = input$adaptive_k_penalty
                    , ruv_k = input$ruv_k
                    , ruv_percentage = input$ruv_percentage
                )
                , experimentPaths = experimentPaths
                , omicType = omicType
                , normData = normData
                , getPlotAestheticsFn = function() {
                    getPlotAestheticsFn(
                        colorVariable = input$color_variable
                        , shapeVariable = input$shape_variable
                    )
                }
                , addLogFn = addLogFn
                , reqFn = reqFn
                , generateCompositeFromFilesFn = generateCompositeFromFilesFn
                , savePlotFn = savePlotFn
                , logWarnFn = logWarnFn
            )
        }
    )
}

