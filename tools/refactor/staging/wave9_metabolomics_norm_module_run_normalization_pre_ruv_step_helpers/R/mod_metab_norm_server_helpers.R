# Keep the log2 progress/apply shell top-level so later waves can move this
# step without reopening the rest of the normalization pipeline body.
runMetabNormLog2ProgressApplyShell <- function(
    currentS4,
    totalSteps,
    logOffset,
    workflowData,
    normData,
    addLogFn = function(entry) invisible(entry),
    incProgressFn = shiny::incProgress,
    runLog2StepFn = runMetabNormLog2TransformationStep
) {
    progressDetail <- "Applying log2 transformation..."
    incProgressFn(1 / totalSteps, detail = progressDetail)

    applyLogEntry <- paste("Applying log2 transformation (offset:", logOffset, ")")
    addLogFn(applyLogEntry)

    log2State <- runLog2StepFn(
        currentS4 = currentS4,
        logOffset = logOffset,
        workflowData = workflowData,
        normData = normData,
        addLogFn = addLogFn
    )

    invisible(list(
        currentS4 = log2State$currentS4,
        progressDetail = progressDetail,
        applyLogEntry = applyLogEntry,
        log2State = log2State
    ))
}

# Keep the log2 apply/saveState block top-level so later waves can move this
# step without reopening the rest of the normalization pipeline body.
runMetabNormLog2TransformationStep <- function(
    currentS4,
    logOffset,
    workflowData,
    normData,
    addLogFn = function(entry) invisible(entry),
    logTransformAssaysFn = logTransformAssays
) {
    stateDescription <- paste("Log2 transformation (offset:", logOffset, ")")

    currentS4 <- logTransformAssaysFn(
        theObject = currentS4,
        offset = logOffset
    )
    normData$post_log2_obj <- currentS4

    workflowData$state_manager$saveState(
        state_name = "metab_log2",
        s4_data_object = currentS4,
        config_object = workflowData$config_list,
        description = stateDescription
    )

    logEntry <- "Log2 transformation complete"
    addLogFn(logEntry)

    invisible(list(
        currentS4 = currentS4,
        stateName = "metab_log2",
        description = stateDescription,
        logEntry = logEntry
    ))
}

# Keep the between-sample normalization progress/apply shell top-level so
# later waves can move this step without reopening the rest of the
# normalization pipeline body.
runMetabNormBetweenSampleProgressApplyShell <- function(
    currentS4,
    totalSteps,
    normMethod,
    workflowData,
    normData,
    addLogFn = function(entry) invisible(entry),
    incProgressFn = shiny::incProgress,
    runBetweenSampleStepFn = runMetabNormBetweenSampleNormalizationStep
) {
    progressDetail <- "Applying between-sample normalization..."
    incProgressFn(1 / totalSteps, detail = progressDetail)

    applyLogEntry <- NULL
    if (normMethod != "none") {
        applyLogEntry <- paste(
            "Applying between-sample normalization (method:",
            normMethod,
            ")"
        )
        addLogFn(applyLogEntry)
    }

    betweenSampleState <- runBetweenSampleStepFn(
        currentS4 = currentS4,
        normMethod = normMethod,
        workflowData = workflowData,
        normData = normData,
        addLogFn = addLogFn
    )

    invisible(list(
        currentS4 = betweenSampleState$currentS4,
        progressDetail = progressDetail,
        applyLogEntry = applyLogEntry,
        betweenSampleState = betweenSampleState
    ))
}

# Keep the between-sample normalization apply/saveState block top-level so
# later waves can move this step without reopening the rest of the
# normalization pipeline body.
runMetabNormBetweenSampleNormalizationStep <- function(
    currentS4,
    normMethod,
    workflowData,
    normData,
    addLogFn = function(entry) invisible(entry),
    normaliseBetweenSamplesFn = normaliseBetweenSamples
) {
    stateDescription <- paste("Between-sample normalization (method:", normMethod, ")")

    if (normMethod != "none") {
        currentS4 <- normaliseBetweenSamplesFn(
            theObject = currentS4,
            normalisation_method = normMethod
        )
    }
    normData$post_norm_obj <- currentS4

    workflowData$state_manager$saveState(
        state_name = "metab_normalized",
        s4_data_object = currentS4,
        config_object = workflowData$config_list,
        description = stateDescription
    )

    logEntry <- "Between-sample normalization complete"
    addLogFn(logEntry)
    normData$normalization_complete <- TRUE

    invisible(list(
        currentS4 = currentS4,
        stateName = "metab_normalized",
        description = stateDescription,
        logEntry = logEntry
    ))
}

# Keep the post-normalization QC generation block top-level so later waves can
# move this step without reopening the rest of the normalization pipeline body.
runMetabNormPostNormalizationQcStep <- function(
    currentS4,
    experimentPaths,
    groupingVariable,
    shapeVariable,
    addLogFn = function(entry) invisible(entry),
    generateMetabQcPlotsFn = generateMetabQcPlots
) {
    generateMetabQcPlotsFn(
        theObject = currentS4,
        experiment_paths = experimentPaths,
        stage = "post_norm",
        grouping_variable = groupingVariable,
        shape_variable = shapeVariable
    )

    logEntry <- "Post-normalization QC plots generated"
    addLogFn(logEntry)

    invisible(list(
        currentS4 = currentS4,
        stage = "post_norm",
        logEntry = logEntry
    ))
}

