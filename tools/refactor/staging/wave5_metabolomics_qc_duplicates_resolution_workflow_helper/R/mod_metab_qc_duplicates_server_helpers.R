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

