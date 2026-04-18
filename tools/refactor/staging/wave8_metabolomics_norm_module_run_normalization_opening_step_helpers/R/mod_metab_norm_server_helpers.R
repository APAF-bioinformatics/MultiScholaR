# Keep manual ITSD feature-ID resolution top-level so later waves can move this
# step without reopening the full normalization pipeline body.
resolveMetabNormManualItsdFeatureIds <- function(
    currentS4,
    itsdSelections,
    addLogFn = function(entry) invisible(entry),
    buildSelectionTableFn = buildItsdSelectionTable,
    mapSelectionsFn = purrr::imap,
    compactFn = purrr::compact
) {
    itsdFeatureIds <- NULL
    hasManualSelections <- any(vapply(itsdSelections, \(selection) length(selection) > 0, logical(1)))

    if (!hasManualSelections) {
        return(NULL)
    }

    metaboliteIdCol <- currentS4@metabolite_id_column
    annotationCol <- currentS4@annotation_id_column

    itsdFeatureIds <- mapSelectionsFn(itsdSelections, \(rowIndices, assayName) {
        if (is.null(rowIndices) || length(rowIndices) == 0) {
            return(NULL)
        }

        assayData <- currentS4@metabolite_data[[assayName]]
        if (is.null(assayData)) {
            return(NULL)
        }

        selectionTable <- buildSelectionTableFn(
            assay_data = assayData,
            metabolite_id_col = metaboliteIdCol,
            annotation_cols = annotationCol
        )

        selectedIds <- selectionTable$feature_id[rowIndices]
        addLogFn(paste("Assay", assayName, ":", length(selectedIds), "ITSD features selected"))
        selectedIds
    })

    itsdFeatureIds <- compactFn(itsdFeatureIds)
    if (length(itsdFeatureIds) == 0) {
        return(NULL)
    }

    itsdFeatureIds
}

# Keep the pre-normalization capture/pre-QC opening block top-level so later
# waves can move this step without reopening the rest of the normalization
# pipeline body.
runMetabNormPreNormalizationQcStep <- function(
    currentS4,
    totalSteps,
    experimentPaths,
    groupingVariable,
    shapeVariable,
    normData,
    addLogFn = function(entry) invisible(entry),
    incProgressFn = shiny::incProgress,
    generateMetabQcPlotsFn = generateMetabQcPlots
) {
    progressDetail <- "Capturing pre-normalization state..."
    incProgressFn(1 / totalSteps, detail = progressDetail)

    normData$post_filter_obj <- currentS4
    captureLogEntry <- "Post-filtering state captured"
    addLogFn(captureLogEntry)

    generateMetabQcPlotsFn(
        theObject = currentS4,
        experiment_paths = experimentPaths,
        stage = "post_filter",
        grouping_variable = groupingVariable,
        shape_variable = shapeVariable
    )

    preQcLogEntry <- "Pre-normalization QC plots generated"
    addLogFn(preQcLogEntry)

    invisible(list(
        currentS4 = currentS4,
        stage = "post_filter",
        progressDetail = progressDetail,
        captureLogEntry = captureLogEntry,
        preQcLogEntry = preQcLogEntry
    ))
}

# Keep the ITSD progress/apply shell top-level so later waves can move this
# step without reopening the rest of the normalization pipeline body.
runMetabNormItsdProgressApplyShell <- function(
    currentS4,
    totalSteps,
    applyItsd,
    itsdAggregation,
    itsdSelections,
    workflowData,
    normData,
    addLogFn = function(entry) invisible(entry),
    incProgressFn = shiny::incProgress,
    resolveManualFeatureIdsFn = resolveMetabNormManualItsdFeatureIds,
    runItsdStepFn = runMetabNormItsdNormalizationStep
) {
    progressDetail <- "Applying ITSD normalization..."
    incProgressFn(1 / totalSteps, detail = progressDetail)

    if (!isTRUE(applyItsd)) {
        skippedLogEntry <- "ITSD normalization skipped"
        addLogFn(skippedLogEntry)

        return(invisible(list(
            currentS4 = currentS4,
            applied = FALSE,
            progressDetail = progressDetail,
            applyLogEntry = NULL,
            skippedLogEntry = skippedLogEntry,
            itsdFeatureIds = NULL,
            itsdState = NULL
        )))
    }

    applyLogEntry <- paste("Applying ITSD normalization (aggregation:", itsdAggregation, ")")
    addLogFn(applyLogEntry)

    itsdFeatureIds <- resolveManualFeatureIdsFn(
        currentS4 = currentS4,
        itsdSelections = itsdSelections,
        addLogFn = addLogFn
    )

    itsdState <- runItsdStepFn(
        currentS4 = currentS4,
        itsdAggregation = itsdAggregation,
        itsdFeatureIds = itsdFeatureIds,
        workflowData = workflowData,
        normData = normData,
        addLogFn = addLogFn
    )

    invisible(list(
        currentS4 = itsdState$currentS4,
        applied = TRUE,
        progressDetail = progressDetail,
        applyLogEntry = applyLogEntry,
        skippedLogEntry = NULL,
        itsdFeatureIds = itsdFeatureIds,
        itsdState = itsdState
    ))
}

# Keep the ITSD apply/saveState block top-level so later waves can move this
# step without reopening the rest of the normalization pipeline body.
runMetabNormItsdNormalizationStep <- function(
    currentS4,
    itsdAggregation,
    itsdFeatureIds = NULL,
    workflowData,
    normData,
    addLogFn = function(entry) invisible(entry),
    normaliseUntransformedDataFn = normaliseUntransformedData
) {
    stateDescription <- paste("ITSD normalization (aggregation:", itsdAggregation, ")")

    currentS4 <- normaliseUntransformedDataFn(
        theObject = currentS4,
        method = "ITSD",
        itsd_aggregation = itsdAggregation,
        itsd_feature_ids = itsdFeatureIds
    )
    normData$post_itsd_obj <- currentS4

    workflowData$state_manager$saveState(
        state_name = "metab_itsd_norm",
        s4_data_object = currentS4,
        config_object = workflowData$config_list,
        description = stateDescription
    )

    logEntry <- "ITSD normalization complete"
    addLogFn(logEntry)

    invisible(list(
        currentS4 = currentS4,
        stateName = "metab_itsd_norm",
        description = stateDescription,
        logEntry = logEntry
    ))
}

