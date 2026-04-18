prepareMetabDuplicateResolutionState <- function(
    stateManager
    , resolveDuplicateAssayDataFn = resolveMetabDuplicateAssayData
    , reqFn = shiny::req
    , inheritsFn = inherits
    , expectedClass = "MetaboliteAssayData"
) {
    currentS4 <- stateManager$getState()
    reqFn(currentS4)

    if (!inheritsFn(currentS4, expectedClass)) {
        stop(sprintf("Current state is not a %s object", expectedClass))
    }

    resolutionResult <- resolveDuplicateAssayDataFn(
        assayList = currentS4@metabolite_data
        , metaboliteIdCol = currentS4@metabolite_id_column
    )

    currentS4@metabolite_data <- resolutionResult$resolvedAssayList

    list(
        currentS4 = currentS4
        , statsList = resolutionResult$statsList
    )
}

applyMetabDuplicateResolutionState <- function(
    currentS4
    , statsList
    , workflowData
    , omicType
    , setResolutionStatsFn
    , setFilterPlotFn
    , stateName = "metab_duplicates_resolved"
    , description = "Resolved duplicate metabolite features by keeping highest intensity"
    , stepName = "3_Duplicates_Resolved"
    , updateMetaboliteFilteringFn = updateMetaboliteFiltering
    , logWarnFn = logger::log_warn
) {
    setResolutionStatsFn(statsList)

    workflowData$state_manager$saveState(
        state_name = stateName
        , s4_data_object = currentS4
        , config_object = workflowData$config_list
        , description = description
    )

    qcPlot <- tryCatch({
        updateMetaboliteFilteringFn(
            theObject = currentS4
            , step_name = stepName
            , omics_type = omicType
            , return_grid = TRUE
            , overwrite = TRUE
        )
    }, error = function(e) {
        logWarnFn(paste("Could not generate QC plot:", e$message))
        NULL
    })

    setFilterPlotFn(qcPlot)

    list(
        stateName = stateName
        , qcPlot = qcPlot
    )
}

