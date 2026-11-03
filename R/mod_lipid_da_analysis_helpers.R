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
    resolveCurrentStateFn = function(workflowData) readLipidWorkflowMetadata(workflowData$state_manager)
) {
    logInfo("=== RUN DA ANALYSIS BUTTON CLICKED ===")

    # KIT-075: the state's metadata; the DA reads its components inside the declaration.
    currentS4 <- lipidDaCurrentS4(daData, workflowData)
    if (is.null(currentS4)) {
        currentS4 <- resolveCurrentStateFn(workflowData)
    }

    if (!is.list(currentS4) || !identical(currentS4$class, "LipidomicsAssayData")) {
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

validateLipidDaAnalysisPreflight <- function(
    currentS4,
    contrastsTbl,
    formulaString = "~ 0 + group",
    validateAlignmentFn = validateLipidDesignAlignment,
    validateContrastsFn = validateLipidDesignContrasts
) {
    alignment <- validateAlignmentFn(
        designMatrix = currentS4@design_matrix,
        assayList = currentS4@lipid_data,
        sampleIdCol = currentS4@sample_id,
        groupCol = currentS4@group_id,
        requireReplicates = FALSE
    )
    contrastValidation <- validateContrastsFn(
        designMatrix = currentS4@design_matrix,
        contrastsTbl = contrastsTbl,
        formulaString = formulaString,
        allowEmpty = FALSE
    )
    errors <- c(alignment$errors, contrastValidation$errors)
    warnings <- c(alignment$warnings, contrastValidation$warnings)

    list(
        valid = length(errors) == 0L,
        errors = unique(errors),
        warnings = unique(warnings),
        alignment = alignment,
        contrast_validation = contrastValidation,
        model_terms = contrastValidation$model_terms
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
    validatePreflightFn = validateLipidDaAnalysisPreflight,
    removeNotificationFn = shiny::removeNotification,
    notify = shiny::showNotification,
    executeRunAnalysisFn = executeLipidDaRunAnalysis
) {
    if (is.null(analysisContext)) {
        return(NULL)
    }

    currentS4 <- analysisContext$currentS4
    contrastsTbl <- analysisContext$contrastsTbl

    # KIT-078: under a managed session the declaration's own preflight runs in the worker, on each contrast.
    if (is.null(workflowOmicsManagedActions("lipidomics")) && is.list(currentS4) && identical(currentS4$class, "LipidomicsAssayData") && is.function(validatePreflightFn)) {
        # KIT-075: the alignment check reads the state for this action only.
        preflight <- validatePreflightFn(
            currentS4 = readLipidWorkflowState(workflowData$state_manager, currentS4$reference),
            contrastsTbl = contrastsTbl,
            formulaString = formulaString
        )

        if (!isTRUE(preflight$valid)) {
            errorMessage <- paste(preflight$errors, collapse = "; ")
            removeNotificationFn("da_running")
            notify(errorMessage, type = "error", duration = 5)
            return(list(
                status = "error",
                stage = "preflight",
                preflight = preflight,
                errorMessage = errorMessage
            ))
        }
    }

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
            # KIT-078: on the worker under a managed session, with the tables and figures.
            runner <- workflowOmicsManagedActions("lipidomics")
            if (!is.null(runner)) {
                manager <- workflowData$state_manager
                runWorkflowOmicsManagedDa(lipidomicsDeclarationBinding, runner, manager, getWorkflowStateReference(manager, daData$current_state_name),
                    list(formula_string = formulaString, da_q_val_thresh = daQValThresh, treat_lfc_cutoff = treatLfcCutoff), contrastsTbl,
                    deliver = function(results, export) successFinalizerFn(results = results, daData = daData, workflowData = workflowData, session = session,
                        experimentPaths = experimentPaths, daQValThresh = daQValThresh, treatLfcCutoff = treatLfcCutoff, workerExport = export),
                    failed = function(error) errorFinalizerFn(errorMessage = conditionMessage(error)))
                return(list(status = "dispatched"))
            }
            # KIT-031: the DA is the declared operation on the loaded state (KIT-076: `lipid_da_contrast`, one per contrast),
            # its manifest `lipid.da_results.v2` accepted on the head; the kernel's list drives the display and the files
            # as before (KIT-072: the owner's path is the only path).
            manager <- if (is.list(workflowData) || is.environment(workflowData)) workflowData$state_manager else NULL
            source <- getWorkflowStateReference(manager, if (is.list(daData) || is.environment(daData)) daData$current_state_name else NULL)
            results <- runLipidReferenceDaAnalysis(manager, source, list(formula_string = formulaString,
                da_q_val_thresh = daQValThresh, treat_lfc_cutoff = treatLfcCutoff), contrastsTbl)$results

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
    writeResultsFn = outputLipidDaResultsAllContrasts,
    workerExport = NULL
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

    # KIT-076: a failed contrast fails alone; the run says so instead of a table that is silently missing.
    if (length(results$failed_contrasts) > 0) {
        notify(
            sprintf("DA failed for: %s", paste(names(results$failed_contrasts), collapse = ", ")),
            type = "warning",
            duration = NULL
        )
    }

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

    # KIT-078: a managed session's worker has written them (`workerExport` is its outcome).
    diskWrite <- if (!is.null(workerExport)) reportWorkflowOmicsDaExport(workerExport, notify, logInfo, logWarn) else tryCatch(
        {
            daOutputDir <- experimentPaths$da_output_dir
            publicationGraphsDir <- experimentPaths$publication_graphs_dir

            logInfo(sprintf("   da_output_dir = %s", daOutputDir))
            logInfo(sprintf("   publication_graphs_dir = %s", publicationGraphsDir))

            if (is.null(daOutputDir) || is.null(publicationGraphsDir)) {
                logWarn("   Output directories not configured, skipping file output")
                list(status = "skipped", success = FALSE)
            } else {
                success <- writeResultsFn(
                    da_results_list = lipidDaResultsWithSource(results, daData, workflowData), # KIT-075: the heatmaps' source, read for the export
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
            }
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
