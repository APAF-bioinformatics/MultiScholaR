# Keep the remaining run-normalization pipeline body top-level so later waves
# can move this orchestration shell without reopening the observer wrapper.
runMetabNormNormalizationPipelineShell <- function(
    workflowData,
    inputValues,
    experimentPaths,
    omicType,
    normData,
    getPlotAestheticsFn,
    addLogFn = function(entry) invisible(entry),
    reqFn = shiny::req,
    generateCompositeFromFilesFn = NULL,
    savePlotFn = savePlot,
    logWarnFn = logger::log_warn,
    runPreNormalizationQcStepFn = runMetabNormPreNormalizationQcStep,
    runItsdProgressApplyShellFn = runMetabNormItsdProgressApplyShell,
    runLog2ProgressApplyShellFn = runMetabNormLog2ProgressApplyShell,
    runBetweenSampleProgressApplyShellFn = runMetabNormBetweenSampleProgressApplyShell,
    runPostNormalizationQcStepFn = runMetabNormPostNormalizationQcStep,
    runRuvProgressApplyShellFn = runMetabNormRuvProgressApplyShell,
    runCompositeQcRefreshShellFn = runMetabNormCompositeQcRefreshShell
) {
    currentS4 <- workflowData$state_manager$getState()
    reqFn(currentS4)

    aesthetics <- getPlotAestheticsFn()
    totalSteps <- 6  # Pre-QC, ITSD, Log2, Norm, RUV, Post-QC

    preNormQcState <- runPreNormalizationQcStepFn(
        currentS4 = currentS4
        , totalSteps = totalSteps
        , experimentPaths = experimentPaths
        , groupingVariable = aesthetics$color_var
        , shapeVariable = aesthetics$shape_var
        , normData = normData
        , addLogFn = addLogFn
    )
    currentS4 <- preNormQcState$currentS4

    itsdShellState <- runItsdProgressApplyShellFn(
        currentS4 = currentS4
        , totalSteps = totalSteps
        , applyItsd = inputValues$apply_itsd
        , itsdAggregation = inputValues$itsd_aggregation
        , itsdSelections = normData$itsd_selections
        , workflowData = workflowData
        , normData = normData
        , addLogFn = addLogFn
    )
    currentS4 <- itsdShellState$currentS4

    log2ShellState <- runLog2ProgressApplyShellFn(
        currentS4 = currentS4
        , totalSteps = totalSteps
        , logOffset = inputValues$log_offset
        , workflowData = workflowData
        , normData = normData
        , addLogFn = addLogFn
    )
    currentS4 <- log2ShellState$currentS4

    betweenSampleShellState <- runBetweenSampleProgressApplyShellFn(
        currentS4 = currentS4
        , totalSteps = totalSteps
        , normMethod = inputValues$norm_method
        , workflowData = workflowData
        , normData = normData
        , addLogFn = addLogFn
    )
    currentS4 <- betweenSampleShellState$currentS4

    postNormQcState <- runPostNormalizationQcStepFn(
        currentS4 = currentS4
        , experimentPaths = experimentPaths
        , groupingVariable = aesthetics$color_var
        , shapeVariable = aesthetics$shape_var
        , addLogFn = addLogFn
    )
    currentS4 <- postNormQcState$currentS4

    ruvShellState <- runRuvProgressApplyShellFn(
        currentS4 = currentS4
        , totalSteps = totalSteps
        , ruvMode = inputValues$ruv_mode
        , autoPercentageMin = inputValues$auto_percentage_min
        , autoPercentageMax = inputValues$auto_percentage_max
        , ruvGroupingVariable = inputValues$ruv_grouping_variable
        , separationMetric = inputValues$separation_metric
        , kPenaltyWeight = inputValues$k_penalty_weight
        , adaptiveKPenalty = inputValues$adaptive_k_penalty
        , manualK = inputValues$ruv_k
        , manualPercentage = inputValues$ruv_percentage
        , experimentPaths = experimentPaths
        , groupingVariable = aesthetics$color_var
        , shapeVariable = aesthetics$shape_var
        , workflowData = workflowData
        , normData = normData
        , addLogFn = addLogFn
    )
    currentS4 <- ruvShellState$currentS4

    compositeQcRefreshState <- runCompositeQcRefreshShellFn(
        currentS4 = currentS4,
        experimentPaths = experimentPaths,
        assayNames = normData$assay_names,
        ruvMode = inputValues$ruv_mode,
        omicType = omicType,
        normData = normData,
        addLogFn = addLogFn,
        generateCompositeFromFilesFn = generateCompositeFromFilesFn,
        savePlotFn = savePlotFn,
        logWarnFn = logWarnFn
    )

    invisible(list(
        currentS4 = compositeQcRefreshState$currentS4,
        totalSteps = totalSteps,
        compositeQcState = compositeQcRefreshState$compositeQcState,
        plotRefreshTrigger = compositeQcRefreshState$plotRefreshTrigger
    ))
}

