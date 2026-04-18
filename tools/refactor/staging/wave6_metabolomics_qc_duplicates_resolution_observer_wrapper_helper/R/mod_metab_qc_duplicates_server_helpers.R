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

