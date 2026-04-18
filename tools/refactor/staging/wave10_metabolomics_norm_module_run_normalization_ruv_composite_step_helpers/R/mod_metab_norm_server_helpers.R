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

