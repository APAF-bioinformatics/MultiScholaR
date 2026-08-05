# Keep the observer shell top-level so later waves can move it without
# reopening the run-normalization body itself.
runMetabNormNormalizationObserverShell <- function(
    workflowData,
    addLogFn = function(entry) invisible(entry),
    showNotificationFn = shiny::showNotification,
    reqFn = shiny::req,
    withProgressFn = shiny::withProgress,
    runPipelineFn,
    logErrorFn = logger::log_error
) {
    reqFn(workflowData$state_manager)
    addLogFn("Starting normalization pipeline...")

    shellState <- withProgressFn(
        message = "Running normalization pipeline..."
        , value = 0
        , {
            tryCatch({
                pipelineState <- runPipelineFn()

                addLogFn("Normalization pipeline complete!")
                showNotificationFn(
                    "Normalization pipeline complete!"
                    , type = "message"
                )

                invisible(list(
                    status = "success"
                    , pipelineState = pipelineState
                ))
            }, error = function(e) {
                errorMessage <- if (inherits(e, "condition")) conditionMessage(e) else as.character(e)

                addLogFn(paste("ERROR:", errorMessage))
                logErrorFn(paste("Normalization pipeline error:", errorMessage))
                showNotificationFn(paste("Error:", errorMessage), type = "error")

                invisible(list(
                    status = "error"
                    , errorMessage = errorMessage
                ))
            })
        }
    )

    invisible(shellState)
}

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

# Keep the RUV-III progress/apply shell top-level so later waves can move
# this step without reopening the rest of the normalization pipeline body.
runMetabNormRuvProgressApplyShell <- function(
    currentS4,
    totalSteps,
    ruvMode,
    autoPercentageMin,
    autoPercentageMax,
    ruvGroupingVariable,
    separationMetric,
    kPenaltyWeight,
    adaptiveKPenalty,
    manualK,
    manualPercentage,
    experimentPaths,
    groupingVariable,
    shapeVariable,
    workflowData,
    normData,
    addLogFn = function(entry) invisible(entry),
    incProgressFn = shiny::incProgress,
    runRuvOptimizationStepFn = runMetabNormRuvOptimizationStep,
    runRuvCorrectionStepFn = runMetabNormRuvCorrectionStep,
    runRuvQcStepFn = runMetabNormRuvQcStep
) {
    progressDetail <- "Running RUV-III batch correction..."
    incProgressFn(1 / totalSteps, detail = progressDetail)

    applyLogEntry <- NULL
    skipLogEntry <- NULL
    optimizationState <- NULL
    correctionState <- NULL
    ruvQcState <- NULL

    if (ruvMode != "skip") {
        applyLogEntry <- paste("Running RUV-III (mode:", ruvMode, ")")
        addLogFn(applyLogEntry)

        optimizationState <- runRuvOptimizationStepFn(
            currentS4 = currentS4,
            ruvMode = ruvMode,
            autoPercentageMin = autoPercentageMin,
            autoPercentageMax = autoPercentageMax,
            ruvGroupingVariable = ruvGroupingVariable,
            separationMetric = separationMetric,
            kPenaltyWeight = kPenaltyWeight,
            adaptiveKPenalty = adaptiveKPenalty,
            manualK = manualK,
            manualPercentage = manualPercentage,
            experimentPaths = experimentPaths,
            normData = normData,
            addLogFn = addLogFn
        )

        correctionState <- runRuvCorrectionStepFn(
            currentS4 = currentS4,
            ruvGroupingVariable = ruvGroupingVariable,
            bestKPerAssay = optimizationState$bestKPerAssay,
            ctrlPerAssay = optimizationState$ctrlPerAssay,
            workflowData = workflowData,
            normData = normData,
            addLogFn = addLogFn
        )
        currentS4 <- correctionState$currentS4

        ruvQcState <- runRuvQcStepFn(
            currentS4 = currentS4,
            totalSteps = totalSteps,
            experimentPaths = experimentPaths,
            groupingVariable = groupingVariable,
            shapeVariable = shapeVariable,
            addLogFn = addLogFn
        )
        currentS4 <- ruvQcState$currentS4
    } else {
        skipLogEntry <- "RUV-III skipped"
        addLogFn(skipLogEntry)
        normData$ruv_corrected_obj <- currentS4
        normData$ruv_complete <- TRUE
    }

    invisible(list(
        currentS4 = currentS4,
        progressDetail = progressDetail,
        applyLogEntry = applyLogEntry,
        skipLogEntry = skipLogEntry,
        optimizationState = optimizationState,
        correctionState = correctionState,
        ruvQcState = ruvQcState
    ))
}

# Keep the opening RUV optimization block top-level so later waves can move
# this step without reopening the rest of the normalization pipeline body.
runMetabNormRuvOptimizationStep <- function(
    currentS4,
    ruvMode,
    autoPercentageMin,
    autoPercentageMax,
    ruvGroupingVariable,
    separationMetric,
    kPenaltyWeight,
    adaptiveKPenalty,
    manualK,
    manualPercentage,
    experimentPaths,
    normData,
    addLogFn = function(entry) invisible(entry),
    runPerAssayRuvOptimizationFn = runPerAssayRuvOptimization,
    extractBestKPerAssayFn = extractBestKPerAssay,
    extractCtrlPerAssayFn = extractCtrlPerAssay
) {
    ruvParams <- list(
        percentage_min = autoPercentageMin,
        percentage_max = autoPercentageMax,
        ruv_grouping_variable = ruvGroupingVariable,
        separation_metric = separationMetric,
        k_penalty_weight = kPenaltyWeight,
        adaptive_k_penalty = adaptiveKPenalty,
        manual_k = manualK,
        manual_percentage = manualPercentage
    )

    ruvResults <- runPerAssayRuvOptimizationFn(
        theObject = currentS4,
        ruv_mode = ruvMode,
        params = ruvParams,
        experiment_paths = experimentPaths
    )
    normData$ruv_optimization_results <- ruvResults

    bestKPerAssay <- extractBestKPerAssayFn(ruvResults)
    ctrlPerAssay <- extractCtrlPerAssayFn(ruvResults)

    logEntry <- paste(
        "RUV optimization complete. Best k per assay:",
        paste(names(bestKPerAssay), "=", unlist(bestKPerAssay), collapse = ", ")
    )
    addLogFn(logEntry)

    invisible(list(
        currentS4 = currentS4,
        ruvParams = ruvParams,
        ruvResults = ruvResults,
        bestKPerAssay = bestKPerAssay,
        ctrlPerAssay = ctrlPerAssay,
        logEntry = logEntry
    ))
}

# Keep the RUV-III apply/saveState block top-level so later waves can move
# this step without reopening the rest of the normalization pipeline body.
runMetabNormRuvCorrectionStep <- function(
    currentS4,
    ruvGroupingVariable,
    bestKPerAssay,
    ctrlPerAssay,
    workflowData,
    normData,
    addLogFn = function(entry) invisible(entry),
    ruvIII_C_VaryingFn = ruvIII_C_Varying
) {
    stateDescription <- "RUV-III batch correction complete"

    currentS4 <- ruvIII_C_VaryingFn(
        theObject = currentS4,
        ruv_grouping_variable = ruvGroupingVariable,
        ruv_number_k = bestKPerAssay,
        ctrl = ctrlPerAssay
    )
    normData$ruv_corrected_obj <- currentS4
    normData$ruv_complete <- TRUE

    workflowData$state_manager$saveState(
        state_name = "metab_ruv_corrected",
        s4_data_object = currentS4,
        config_object = workflowData$config_list,
        description = stateDescription
    )

    logEntry <- "RUV-III correction applied"
    addLogFn(logEntry)

    invisible(list(
        currentS4 = currentS4,
        stateName = "metab_ruv_corrected",
        description = stateDescription,
        logEntry = logEntry
    ))
}

# Keep the RUV QC generation block top-level so later waves can move this step
# without reopening the rest of the normalization pipeline body.
runMetabNormRuvQcStep <- function(
    currentS4,
    totalSteps,
    experimentPaths,
    groupingVariable,
    shapeVariable,
    addLogFn = function(entry) invisible(entry),
    incProgressFn = shiny::incProgress,
    generateMetabQcPlotsFn = generateMetabQcPlots
) {
    progressDetail <- "Generating RUV QC plots..."
    incProgressFn(1 / totalSteps, detail = progressDetail)

    generateMetabQcPlotsFn(
        theObject = currentS4,
        experiment_paths = experimentPaths,
        stage = "ruv_corrected",
        grouping_variable = groupingVariable,
        shape_variable = shapeVariable
    )

    logEntry <- "RUV QC plots generated"
    addLogFn(logEntry)

    invisible(list(
        currentS4 = currentS4,
        stage = "ruv_corrected",
        progressDetail = progressDetail,
        logEntry = logEntry
    ))
}

# Keep the composite QC figure generation block top-level so later waves can
# move this step without reopening the rest of the normalization pipeline body.
runMetabNormCompositeQcFigureStep <- function(
    experimentPaths,
    assayNames,
    ruvMode,
    omicType,
    addLogFn = function(entry) invisible(entry),
    generateCompositeFromFilesFn = NULL,
    savePlotFn = savePlot,
    dirExistsFn = dir.exists,
    logWarnFn = logger::log_warn
) {
    if (is.null(generateCompositeFromFilesFn)) {
        stop("generateCompositeFromFilesFn must be provided")
    }

    logEntry <- "Generating composite QC figure..."
    addLogFn(logEntry)

    tryCatch({
        qcDir <- experimentPaths$metabolite_qc_dir
        if (is.null(qcDir) || !dirExistsFn(qcDir) || is.null(assayNames)) {
            return(invisible(list(
                qcDir = qcDir,
                ncolComposite = NULL,
                columnLabels = NULL,
                allPlotFiles = character(0),
                allRowLabels = list(),
                logEntry = logEntry,
                compositeSaved = FALSE
            )))
        }

        ncolComposite <- if (identical(ruvMode, "skip")) 2 else 3
        columnLabels <- if (identical(ruvMode, "skip")) {
            c("Pre-Normalisation", "Post-Normalisation")
        } else {
            c("Pre-Normalisation", "Post-Normalisation", "RUV-Corrected")
        }

        allPlotFiles <- c()
        allRowLabels <- list()
        labelCounter <- 1L
        plotTypes <- c("pca", "density", "rle", "correlation")

        for (assayName in assayNames) {
            safeName <- gsub("[^A-Za-z0-9]", "_", tolower(assayName))

            for (plotType in plotTypes) {
                if (identical(ruvMode, "skip")) {
                    files <- c(
                        file.path(qcDir, sprintf("%s_pre_norm_%s.png", safeName, plotType)),
                        file.path(qcDir, sprintf("%s_post_norm_%s.png", safeName, plotType))
                    )
                    labels <- c(
                        sprintf("%s)", letters[labelCounter]),
                        sprintf("%s)", letters[labelCounter + 1L])
                    )
                    labelCounter <- labelCounter + 2L
                } else {
                    files <- c(
                        file.path(qcDir, sprintf("%s_pre_norm_%s.png", safeName, plotType)),
                        file.path(qcDir, sprintf("%s_post_norm_%s.png", safeName, plotType)),
                        file.path(qcDir, sprintf("%s_ruv_corrected_%s.png", safeName, plotType))
                    )
                    labels <- c(
                        sprintf("%s)", letters[labelCounter]),
                        sprintf("%s)", letters[labelCounter + 1L]),
                        sprintf("%s)", letters[labelCounter + 2L])
                    )
                    labelCounter <- labelCounter + 3L
                }

                allPlotFiles <- c(allPlotFiles, files)
                rowKey <- sprintf("%s_%s", safeName, plotType)
                allRowLabels[[rowKey]] <- labels
            }
        }

        compositeRes <- generateCompositeFromFilesFn(
            plot_files = allPlotFiles,
            ncol = ncolComposite,
            row_labels = allRowLabels,
            column_labels = columnLabels
        )

        if (!is.null(compositeRes)) {
            savePlotFn(
                compositeRes$plot,
                qcDir,
                paste0(omicType, "_composite_QC_figure"),
                width = compositeRes$width,
                height = compositeRes$height,
                dpi = 150,
                limitsize = FALSE
            )
            addLogFn(sprintf("Composite QC figure saved to: %s", file.path(qcDir, "composite_QC_figure")))
        }

        invisible(list(
            qcDir = qcDir,
            ncolComposite = ncolComposite,
            columnLabels = columnLabels,
            allPlotFiles = allPlotFiles,
            allRowLabels = allRowLabels,
            logEntry = logEntry,
            compositeSaved = !is.null(compositeRes)
        ))
    }, error = function(e) {
        errorMessage <- conditionMessage(e)
        addLogFn(paste("Warning: Could not generate composite QC figure:", errorMessage))
        logWarnFn(paste("Composite QC generation failed:", errorMessage))

        invisible(list(
            qcDir = experimentPaths$metabolite_qc_dir,
            ncolComposite = NULL,
            columnLabels = NULL,
            allPlotFiles = character(0),
            allRowLabels = list(),
            logEntry = logEntry,
            errorMessage = errorMessage,
            compositeSaved = FALSE
        ))
    })
}

# Keep the final composite-QC / plot-refresh tail top-level so later waves can
# move this closing shell without reopening the rest of the normalization
# pipeline body.
runMetabNormCompositeQcRefreshShell <- function(
    currentS4,
    experimentPaths,
    assayNames,
    ruvMode,
    omicType,
    normData,
    addLogFn = function(entry) invisible(entry),
    generateCompositeFromFilesFn = NULL,
    savePlotFn = savePlot,
    logWarnFn = logger::log_warn,
    runCompositeQcFigureStepFn = runMetabNormCompositeQcFigureStep
) {
    compositeQcState <- runCompositeQcFigureStepFn(
        experimentPaths = experimentPaths,
        assayNames = assayNames,
        ruvMode = ruvMode,
        omicType = omicType,
        addLogFn = addLogFn,
        generateCompositeFromFilesFn = generateCompositeFromFilesFn,
        savePlotFn = savePlotFn,
        logWarnFn = logWarnFn
    )

    normData$plot_refresh_trigger <- normData$plot_refresh_trigger + 1

    invisible(list(
        currentS4 = currentS4,
        compositeQcState = compositeQcState,
        plotRefreshTrigger = normData$plot_refresh_trigger
    ))
}

