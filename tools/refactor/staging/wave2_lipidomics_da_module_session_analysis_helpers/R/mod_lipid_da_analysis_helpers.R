bootstrapLipidDaRunAnalysis <- function(
    daData,
    workflowData,
    formulaString,
    daQValThresh,
    treatLfcCutoff,
    session,
    experimentPaths,
    prepareAnalysisContextFn = prepareLipidDaRunAnalysisContext,
    handlePreflightFn = handleLipidDaRunAnalysisPreflight
) {
    analysisContext <- prepareAnalysisContextFn(
        daData = daData,
        workflowData = workflowData
    )

    handlePreflightFn(
        analysisContext = analysisContext,
        formulaString = formulaString,
        daQValThresh = daQValThresh,
        treatLfcCutoff = treatLfcCutoff,
        daData = daData,
        workflowData = workflowData,
        session = session,
        experimentPaths = experimentPaths
    )
}

prepareLipidDaRunAnalysisContext <- function(
    daData,
    workflowData,
    notify = shiny::showNotification,
    logInfo = logger::log_info,
    resolveCurrentStateFn = function(workflowData) workflowData$state_manager$getState()
) {
    logInfo("=== RUN DA ANALYSIS BUTTON CLICKED ===")

    currentS4 <- daData$current_s4_object
    if (is.null(currentS4)) {
        currentS4 <- resolveCurrentStateFn(workflowData)
    }

    if (is.null(currentS4) || !inherits(currentS4, "LipidomicsAssayData")) {
        notify(
            "No lipidomics data loaded. Please load a filtered session first.",
            type = "error",
            duration = 5
        )
        return(NULL)
    }

    contrastsTbl <- workflowData$contrasts_tbl
    if (is.null(contrastsTbl) || nrow(contrastsTbl) == 0) {
        notify(
            "No contrasts defined. Please define contrasts in the design tab.",
            type = "error",
            duration = 5
        )
        return(NULL)
    }

    runningNotification <- list(
        message = "Running differential expression analysis...",
        id = "da_running",
        duration = NULL
    )

    notify(
        runningNotification$message,
        id = runningNotification$id,
        duration = runningNotification$duration
    )

    list(
        currentS4 = currentS4,
        contrastsTbl = contrastsTbl,
        runningNotification = runningNotification
    )
}

handleLipidDaRunAnalysisPreflight <- function(
    analysisContext,
    formulaString,
    daQValThresh,
    treatLfcCutoff,
    daData,
    workflowData,
    session,
    experimentPaths,
    executeRunAnalysisFn = executeLipidDaRunAnalysis
) {
    if (is.null(analysisContext)) {
        return(NULL)
    }

    currentS4 <- analysisContext$currentS4
    contrastsTbl <- analysisContext$contrastsTbl

    executeRunAnalysisFn(
        currentS4 = currentS4,
        contrastsTbl = contrastsTbl,
        formulaString = formulaString,
        daQValThresh = daQValThresh,
        treatLfcCutoff = treatLfcCutoff,
        daData = daData,
        workflowData = workflowData,
        session = session,
        experimentPaths = experimentPaths
    )
}

executeLipidDaRunAnalysis <- function(
    currentS4,
    contrastsTbl,
    formulaString,
    daQValThresh,
    treatLfcCutoff,
    daData,
    workflowData,
    session,
    experimentPaths,
    runLipidsDaFn = runLipidsDA,
    successFinalizerFn = finalizeLipidDaRunAnalysisSuccess,
    errorFinalizerFn = finalizeLipidDaRunAnalysisError
) {
    tryCatch(
        {
            results <- runLipidsDaFn(
                theObject = currentS4,
                contrasts_tbl = contrastsTbl,
                formula_string = formulaString,
                da_q_val_thresh = daQValThresh,
                treat_lfc_cutoff = treatLfcCutoff,
                eBayes_trend = TRUE,
                eBayes_robust = TRUE
            )

            list(
                status = "success",
                results = results,
                finalizerResult = successFinalizerFn(
                    results = results,
                    daData = daData,
                    workflowData = workflowData,
                    session = session,
                    experimentPaths = experimentPaths,
                    daQValThresh = daQValThresh,
                    treatLfcCutoff = treatLfcCutoff
                )
            )
        },
        error = function(e) {
            list(
                status = "error",
                errorMessage = e$message,
                finalizerResult = errorFinalizerFn(
                    errorMessage = e$message
                )
            )
        }
    )
}

finalizeLipidDaRunAnalysisSuccess <- function(
    results,
    daData,
    workflowData,
    session,
    experimentPaths,
    daQValThresh,
    treatLfcCutoff,
    removeNotificationFn = shiny::removeNotification,
    notify = shiny::showNotification,
    updateSelectInputFn = shiny::updateSelectInput,
    logInfo = logger::log_info,
    logWarn = logger::log_warn,
    writeResultsFn = outputLipidDaResultsAllContrasts
) {
    daData$da_results_list <- results
    daData$analysis_complete <- TRUE

    updatedStatus <- workflowData$tab_status
    updatedStatus$differential_analysis <- "complete"
    workflowData$tab_status <- updatedStatus

    removeNotificationFn("da_running")
    notify(
        "Differential expression analysis complete!",
        type = "message",
        duration = 5
    )

    logInfo("   DA analysis completed successfully")

    contrastChoices <- NULL
    assayChoices <- NULL
    tableAssayChoices <- NULL

    if (!is.null(results$da_lipids_long)) {
        contrastChoices <- unique(results$da_lipids_long$friendly_name)
        if (length(contrastChoices) == 0) {
            contrastChoices <- unique(results$da_lipids_long$comparison)
        }

        updateSelectInputFn(session, "volcano_contrast",
            choices = contrastChoices, selected = contrastChoices[1]
        )
        updateSelectInputFn(session, "heatmap_contrast",
            choices = contrastChoices, selected = contrastChoices[1]
        )
        updateSelectInputFn(session, "table_contrast",
            choices = contrastChoices, selected = contrastChoices[1]
        )

        assayChoices <- c("Combined", unique(results$da_lipids_long$assay))
        updateSelectInputFn(session, "volcano_assay",
            choices = assayChoices, selected = "Combined"
        )
        updateSelectInputFn(session, "heatmap_assay",
            choices = assayChoices, selected = "Combined"
        )

        tableAssayChoices <- c("All", unique(results$da_lipids_long$assay))
        updateSelectInputFn(session, "table_assay",
            choices = tableAssayChoices, selected = "All"
        )

        logInfo(sprintf(
            "   Updated dropdowns: %d contrasts, %d assays",
            length(contrastChoices), length(assayChoices) - 1
        ))
    }

    logInfo("   Writing DA results to disk...")

    diskWrite <- tryCatch(
        {
            daOutputDir <- experimentPaths$da_output_dir
            publicationGraphsDir <- experimentPaths$publication_graphs_dir

            logInfo(sprintf("   da_output_dir = %s", daOutputDir))
            logInfo(sprintf("   publication_graphs_dir = %s", publicationGraphsDir))

            if (is.null(daOutputDir) || is.null(publicationGraphsDir)) {
                logWarn("   Output directories not configured, skipping file output")
                return(list(status = "skipped", success = FALSE))
            }

            success <- writeResultsFn(
                da_results_list = results,
                da_output_dir = daOutputDir,
                publication_graphs_dir = publicationGraphsDir,
                da_q_val_thresh = daQValThresh,
                lfc_threshold = treatLfcCutoff,
                heatmap_top_n = 50,
                heatmap_clustering = "both",
                heatmap_color_scheme = "RdBu"
            )

            if (success) {
                logInfo("   All DA results written to disk successfully")
                notify(
                    "DA results saved to disk (tables, volcano plots, heatmaps)",
                    type = "message",
                    duration = 5
                )
            }

            list(status = "completed", success = success)
        },
        error = function(e) {
            logWarn(paste("   Could not write DA results to disk:", e$message))
            notify(
                paste("Warning: Could not save results to disk:", e$message),
                type = "warning",
                duration = 8
            )

            list(status = "error", success = FALSE, message = e$message)
        }
    )

    list(
        updatedStatus = workflowData$tab_status,
        contrastChoices = contrastChoices,
        assayChoices = assayChoices,
        tableAssayChoices = tableAssayChoices,
        diskWrite = diskWrite
    )
}

finalizeLipidDaRunAnalysisError <- function(
    errorMessage,
    removeNotificationFn = shiny::removeNotification,
    notify = shiny::showNotification,
    logError = logger::log_error
) {
    contract <- list(
        errorMessage = errorMessage,
        loggedMessage = sprintf("   DA analysis error: %s", errorMessage),
        removedNotificationId = "da_running",
        notificationMessage = sprintf("Analysis error: %s", errorMessage),
        notificationType = "error",
        notificationDuration = 10
    )

    logError(contract$loggedMessage)
    removeNotificationFn(contract$removedNotificationId)
    notify(
        contract$notificationMessage,
        type = contract$notificationType,
        duration = contract$notificationDuration
    )

    contract
}

