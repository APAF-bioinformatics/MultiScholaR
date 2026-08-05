runMetabNormApplyCorrelationWorkflow <- function(
    observerState,
    reqFn = shiny::req,
    calculateCorrelationsFn = pearsonCorForSamplePairs,
    filterSamplesFn = filterSamplesByMetaboliteCorrelationThreshold,
    logInfoFn = logger::log_info
) {
    currentS4 <- observerState$currentS4
    threshold <- observerState$threshold

    reqFn(currentS4)

    logInfoFn("Calculating Pearson correlations per sample pair...")
    corrResults <- calculateCorrelationsFn(
        theObject = currentS4
        , correlation_group = observerState$groupingVariable
    )

    filteredS4 <- filterSamplesFn(
        theObject = currentS4
        , pearson_correlation_per_pair = corrResults
        , min_pearson_correlation_threshold = threshold
    )

    invisible(list(
        corrResults = corrResults
        , filteredS4 = filteredS4
    ))
}

dispatchMetabNormApplyCorrelation <- function(
    workflowData,
    normData,
    observerState,
    addLogFn = function(entry) invisible(entry),
    showNotificationFn = shiny::showNotification,
    removeNotificationFn = shiny::removeNotification,
    reqFn = shiny::req,
    calculateCorrelationsFn = pearsonCorForSamplePairs,
    filterSamplesFn = filterSamplesByMetaboliteCorrelationThreshold,
    runWorkflowFn = runMetabNormApplyCorrelationWorkflow,
    handleOutcomeFn = handleMetabNormApplyCorrelationOutcome,
    logInfoFn = logger::log_info,
    logErrorFn = logger::log_error
) {
    tryCatch({
        workflowState <- runWorkflowFn(
            observerState = observerState
            , reqFn = reqFn
            , calculateCorrelationsFn = calculateCorrelationsFn
            , filterSamplesFn = filterSamplesFn
            , logInfoFn = logInfoFn
        )

        handleOutcomeFn(
            workflowData = workflowData,
            normData = normData,
            observerState = observerState,
            corrResults = workflowState$corrResults,
            filteredS4 = workflowState$filteredS4,
            addLogFn = addLogFn,
            showNotificationFn = showNotificationFn,
            removeNotificationFn = removeNotificationFn,
            logErrorFn = logErrorFn
        )
    }, error = function(e) {
        handleOutcomeFn(
            workflowData = workflowData,
            normData = normData,
            observerState = observerState,
            error = e,
            addLogFn = addLogFn,
            showNotificationFn = showNotificationFn,
            removeNotificationFn = removeNotificationFn,
            logErrorFn = logErrorFn
        )
    })
}

handleMetabNormApplyCorrelationOutcome <- function(
    workflowData,
    normData,
    observerState,
    corrResults = NULL,
    filteredS4 = NULL,
    error = NULL,
    addLogFn = function(entry) invisible(entry),
    showNotificationFn = shiny::showNotification,
    removeNotificationFn = shiny::removeNotification,
    logErrorFn = logger::log_error
) {
    if (is.null(error)) {
        threshold <- observerState$threshold

        normData$correlation_results <- corrResults
        normData$correlation_filtered_obj <- filteredS4
        normData$correlation_filtering_complete <- TRUE

        workflowData$state_manager$saveState(
            state_name = "metab_correlation_filtered"
            , s4_data_object = filteredS4
            , config_object = workflowData$config_list
            , description = paste("Correlation filtering (threshold:", threshold, ")")
        )

        updatedStatus <- workflowData$tab_status
        updatedStatus$quality_control <- "complete"
        updatedStatus$normalization <- "complete"
        workflowData$tab_status <- updatedStatus

        addLogFn("Correlation filtering complete")
        removeNotificationFn(observerState$notificationId)
        showNotificationFn("Correlation filtering complete! Ready for DE analysis.", type = "message")

        return(invisible(list(
            status = "success"
            , corrResults = corrResults
            , filteredS4 = filteredS4
            , updatedStatus = updatedStatus
        )))
    }

    errorMessage <- if (inherits(error, "condition")) conditionMessage(error) else as.character(error)

    addLogFn(paste("ERROR in correlation filtering:", errorMessage))
    logErrorFn(paste("Correlation filtering error:", errorMessage))
    removeNotificationFn(observerState$notificationId)
    showNotificationFn(paste("Error:", errorMessage), type = "error")

    invisible(list(
        status = "error"
        , errorMessage = errorMessage
    ))
}

runMetabNormApplyCorrelationObserverEntry <- function(
    workflowData,
    normData,
    threshold,
    groupingVariable,
    addLogFn = function(entry) invisible(entry),
    showNotificationFn = shiny::showNotification,
    removeNotificationFn = shiny::removeNotification,
    resolveInputObjectFn = resolveMetabNormSkipCorrelationInputObject,
    dispatchApplyCorrelationFn = dispatchMetabNormApplyCorrelation
) {
    currentS4 <- resolveInputObjectFn(
        ruvCorrectedObject = normData$ruv_corrected_obj,
        postNormObject = normData$post_norm_obj
    )
    logEntry <- paste("Applying correlation filter (threshold:", threshold, ")")
    notificationId <- "corr_working"

    addLogFn(logEntry)
    showNotificationFn("Applying correlation filter...", id = notificationId, duration = NULL)

    observerState <- list(
        currentS4 = currentS4,
        threshold = threshold,
        groupingVariable = groupingVariable,
        logEntry = logEntry,
        notificationId = notificationId
    )

    dispatchState <- dispatchApplyCorrelationFn(
        workflowData = workflowData,
        normData = normData,
        observerState = observerState,
        addLogFn = addLogFn,
        showNotificationFn = showNotificationFn,
        removeNotificationFn = removeNotificationFn
    )

    invisible(dispatchState)
}

checkMetabNormExportSessionReady <- function(
    normalizationComplete,
    showNotificationFn = shiny::showNotification
) {
    if (isTRUE(normalizationComplete)) {
        return(TRUE)
    }

    showNotificationFn(
        "Please complete normalization before exporting session data."
        , type = "warning"
        , duration = 5
    )

    FALSE
}

dispatchMetabNormExportSession <- function(
    workflowData,
    normData,
    inputValues,
    experimentPaths,
    experimentLabel,
    addLogFn = function(entry) invisible(entry),
    resolveSourceDirFn = resolveMetabNormExportSourceDir,
    runWorkflowFn = runMetabNormExportSessionWorkflow,
    handleOutcomeFn = handleMetabNormExportSessionOutcome
) {
    tryCatch({
        sourceDir <- resolveSourceDirFn(experimentPaths)

        exportFiles <- runWorkflowFn(
            workflowData = workflowData,
            normData = normData,
            inputValues = inputValues,
            experimentLabel = experimentLabel,
            sourceDir = sourceDir
        )

        handleOutcomeFn(
            sessionFilename = exportFiles$sessionFilename,
            sessionFilepath = exportFiles$sessionFilepath,
            addLogFn = addLogFn
        )
    }, error = function(e) {
        handleOutcomeFn(
            error = e,
            addLogFn = addLogFn
        )
    })
}

handleMetabNormExportSessionOutcome <- function(
    sessionFilename = NULL,
    sessionFilepath = NULL,
    error = NULL,
    addLogFn = function(entry) invisible(entry),
    showNotificationFn = shiny::showNotification,
    logInfoFn = logger::log_info,
    logErrorFn = logger::log_error
) {
    if (is.null(error)) {
        addLogFn(paste("Exported comprehensive session data to:", sessionFilepath))
        showNotificationFn(
            sprintf("Session data exported successfully!\nSaved as: %s\nSee summary file for details.", sessionFilename)
            , type = "message"
            , duration = 10
        )
        logInfoFn("=== EXPORT NORMALIZED SESSION COMPLETED SUCCESSFULLY ===")

        return(invisible(list(
            status = "success"
            , sessionFilename = sessionFilename
            , sessionFilepath = sessionFilepath
        )))
    }

    errorMessage <- if (inherits(error, "condition")) conditionMessage(error) else as.character(error)
    logErrorFn(paste("*** ERROR in session export:", errorMessage, "***"))
    addLogFn(paste("Export error:", errorMessage))
    showNotificationFn(paste("Export error:", errorMessage), type = "error", duration = 10)

    invisible(list(
        status = "error"
        , errorMessage = errorMessage
    ))
}

runMetabNormExportSessionObserverShell <- function(
    workflowData,
    normData,
    inputValues,
    experimentPaths,
    experimentLabel,
    addLogFn = function(entry) invisible(entry),
    logInfoFn = logger::log_info,
    reqFn = shiny::req,
    checkReadyFn = checkMetabNormExportSessionReady,
    dispatchExportSessionFn = dispatchMetabNormExportSession
) {
    logInfoFn("=== EXPORT NORMALIZED SESSION BUTTON CLICKED ===")
    reqFn(workflowData$state_manager)

    if (!isTRUE(checkReadyFn(normData$normalization_complete))) {
        return(invisible(NULL))
    }

    dispatchState <- dispatchExportSessionFn(
        workflowData = workflowData,
        normData = normData,
        inputValues = inputValues,
        experimentPaths = experimentPaths,
        experimentLabel = experimentLabel,
        addLogFn = addLogFn
    )

    invisible(dispatchState)
}

runMetabNormResetNormalizationObserverShell <- function(
    workflowData,
    normData,
    addLogFn = function(entry) invisible(entry),
    showNotificationFn = shiny::showNotification,
    reqFn = shiny::req
) {
    reqFn(workflowData$state_manager)
    postFilterObject <- normData$post_filter_obj

    tryCatch({
        if (!is.null(postFilterObject)) {
            workflowData$state_manager$saveState(
                state_name = "metab_reset"
                , s4_data_object = postFilterObject
                , config_object = workflowData$config_list
                , description = "Reset to pre-normalization state"
            )
        }

        normData$normalization_complete <- FALSE
        normData$ruv_complete <- FALSE
        normData$correlation_filtering_complete <- FALSE
        normData$post_norm_obj <- NULL
        normData$ruv_corrected_obj <- NULL
        normData$correlation_filtered_obj <- NULL
        normData$ruv_optimization_results <- list()

        addLogFn("Reset to pre-normalization state")
        showNotificationFn("Reset to pre-normalization state", type = "message")

        invisible(list(
            status = "success"
            , stateSaved = !is.null(postFilterObject)
            , logEntry = "Reset to pre-normalization state"
        ))
    }, error = function(e) {
        errorMessage <- if (inherits(e, "condition")) conditionMessage(e) else as.character(e)

        addLogFn(paste("ERROR during reset:", errorMessage))
        showNotificationFn(paste("Error:", errorMessage), type = "error")

        invisible(list(
            status = "error"
            , errorMessage = errorMessage
        ))
    })
}

resolveMetabNormSkipCorrelationInputObject <- function(ruvCorrectedObject, postNormObject) {
    if (!is.null(ruvCorrectedObject)) {
        return(ruvCorrectedObject)
    }

    if (!is.null(postNormObject)) {
        return(postNormObject)
    }

    NULL
}

completeMetabNormSkipCorrelationState <- function(workflowData, currentS4) {
    workflowData$state_manager$saveState(
        state_name = "metab_norm_complete"
        , s4_data_object = currentS4
        , config_object = workflowData$config_list
        , description = "Normalization complete (correlation filtering skipped)"
    )

    updated_status <- workflowData$tab_status
    updated_status$quality_control <- "complete"
    updated_status$normalization <- "complete"
    workflowData$tab_status <- updated_status

    invisible(updated_status)
}

handleMetabNormSkipCorrelationOutcome <- function(
    addLogFn = function(entry) invisible(entry),
    showNotificationFn = shiny::showNotification
) {
    logEntry <- "Correlation filtering skipped - ready for DE analysis"
    notificationMessage <- "Normalization complete! Proceeding to DE analysis."

    addLogFn(logEntry)
    showNotificationFn(notificationMessage, type = "message")

    invisible(list(
        status = "success"
        , logEntry = logEntry
        , notificationMessage = notificationMessage
    ))
}

dispatchMetabNormSkipCorrelation <- function(
    workflowData,
    currentS4,
    addLogFn = function(entry) invisible(entry),
    showNotificationFn = shiny::showNotification,
    completeStateFn = completeMetabNormSkipCorrelationState,
    handleOutcomeFn = handleMetabNormSkipCorrelationOutcome
) {
    if (is.null(currentS4)) {
        return(invisible(NULL))
    }

    updatedStatus <- completeStateFn(
        workflowData = workflowData,
        currentS4 = currentS4
    )
    outcome <- handleOutcomeFn(
        addLogFn = addLogFn,
        showNotificationFn = showNotificationFn
    )

    status <- NULL
    if (is.list(outcome) && !is.null(outcome$status)) {
        status <- outcome$status
    }

    invisible(list(
        status = status,
        updatedStatus = updatedStatus,
        outcome = outcome
    ))
}

runMetabNormSkipCorrelationObserverEntry <- function(
    workflowData,
    normData,
    addLogFn = function(entry) invisible(entry),
    showNotificationFn = shiny::showNotification,
    resolveInputObjectFn = resolveMetabNormSkipCorrelationInputObject,
    dispatchSkipCorrelationFn = dispatchMetabNormSkipCorrelation
) {
    currentS4 <- resolveInputObjectFn(
        ruvCorrectedObject = normData$ruv_corrected_obj,
        postNormObject = normData$post_norm_obj
    )

    dispatchState <- dispatchSkipCorrelationFn(
        workflowData = workflowData,
        currentS4 = currentS4,
        addLogFn = addLogFn,
        showNotificationFn = showNotificationFn
    )

    if (!is.null(dispatchState)) {
        normData$correlation_filtering_complete <- TRUE
        normData$correlation_filtered_obj <- currentS4
    }

    invisible(dispatchState)
}

runMetabNormSkipCorrelationObserverShell <- function(
    workflowData,
    normData,
    addLogFn = function(entry) invisible(entry),
    showNotificationFn = shiny::showNotification,
    reqFn = shiny::req,
    runObserverEntryFn = runMetabNormSkipCorrelationObserverEntry
) {
    reqFn(workflowData$state_manager)
    reqFn(normData$ruv_complete || normData$normalization_complete)

    dispatchState <- runObserverEntryFn(
        workflowData = workflowData,
        normData = normData,
        addLogFn = addLogFn,
        showNotificationFn = showNotificationFn
    )

    invisible(dispatchState)
}

resolveMetabNormExportSourceDir <- function(experimentPaths, dirExistsFn = dir.exists, dirCreateFn = dir.create) {
    sourceDir <- experimentPaths$source_dir
    if (!is.null(sourceDir) && dirExistsFn(sourceDir)) {
        return(sourceDir)
    }

    sourceDir <- experimentPaths$export_dir
    if (is.null(sourceDir)) {
        stop("Could not find a valid directory to save session data.")
    }

    if (!dirExistsFn(sourceDir)) {
        dirCreateFn(sourceDir, recursive = TRUE)
    }

    sourceDir
}

collectMetabNormFeatureCountsPerAssay <- function(currentS4) {
    if (is.null(currentS4) || !inherits(currentS4, "MetaboliteAssayData")) {
        return(NULL)
    }

    assayNames <- names(currentS4@metabolite_data)

    purrr::map(assayNames, \(assay_name) {
        assay_data <- currentS4@metabolite_data[[assay_name]]
        if (!is.null(assay_data)) {
            n_features <- length(unique(assay_data[[currentS4@metabolite_id_column]]))
            n_samples <- ncol(assay_data) - length(c(currentS4@metabolite_id_column, currentS4@annotation_id_column))
            list(features = n_features, samples = n_samples)
        } else {
            list(features = 0, samples = 0)
        }
    }) |> purrr::set_names(assayNames)
}

buildMetabNormExportSessionData <- function(workflowData, normData, inputValues, experimentLabel, exportTimestamp = Sys.time()) {
    currentStateName <- workflowData$state_manager$current_state
    currentS4 <- workflowData$state_manager$getState(currentStateName)
    featureCountsPerAssay <- collectMetabNormFeatureCountsPerAssay(currentS4)

    list(
        # --- R6 State Info ---
        r6_current_state_name = currentStateName
        , current_s4_object = currentS4

        # --- Workflow artifacts ---
        , contrasts_tbl = workflowData$contrasts_tbl
        , design_matrix = workflowData$design_matrix
        , config_list = workflowData$config_list

        # --- Metabolomics-specific ---
        , itsd_selections = normData$itsd_selections
        , ruv_optimization_results = normData$ruv_optimization_results
        , correlation_results = normData$correlation_results
        , assay_names = normData$assay_names

        # --- Export metadata ---
        , export_timestamp = exportTimestamp
        , omic_type = "metabolomics"
        , experiment_label = experimentLabel

        # --- Normalization parameters ---
        , normalization_method = inputValues$norm_method
        , ruv_mode = inputValues$ruv_mode
        , itsd_applied = inputValues$apply_itsd
        , itsd_aggregation = if (isTRUE(inputValues$apply_itsd)) inputValues$itsd_aggregation else NA
        , log_offset = inputValues$log_offset
        , correlation_threshold = inputValues$min_pearson_correlation_threshold
        , ruv_grouping_variable = inputValues$ruv_grouping_variable

        # --- Feature counts per assay ---
        , feature_counts = featureCountsPerAssay
        , metabolite_counts = workflowData$metabolite_counts

        # --- QC parameters ---
        , qc_params = workflowData$qc_params

        # --- Processing flags ---
        , normalization_complete = normData$normalization_complete
        , ruv_complete = normData$ruv_complete
        , correlation_filtering_complete = normData$correlation_filtering_complete
    )
}

saveMetabNormExportSessionRdsFiles <- function(
    sessionData,
    sourceDir,
    timeFn = Sys.time,
    formatTimeFn = format,
    saveRdsFn = saveRDS,
    logInfoFn = logger::log_info,
    incProgressFn = NULL
) {
    if (!is.null(incProgressFn)) {
        incProgressFn(0.3, detail = "Saving to file...")
    }

    timestampStr <- formatTimeFn(timeFn(), "%Y%m%d_%H%M%S")
    sessionFilename <- sprintf("metab_filtered_session_data_%s.rds", timestampStr)
    sessionFilepath <- file.path(sourceDir, sessionFilename)

    saveRdsFn(sessionData, sessionFilepath)
    logInfoFn(sprintf("*** EXPORT: Session data saved to: %s ***", sessionFilepath))

    if (!is.null(incProgressFn)) {
        incProgressFn(0.1, detail = "Creating latest version...")
    }

    latestFilename <- "metab_filtered_session_data_latest.rds"
    latestFilepath <- file.path(sourceDir, latestFilename)

    saveRdsFn(sessionData, latestFilepath)
    logInfoFn(sprintf("*** EXPORT: Latest version saved to: %s ***", latestFilepath))

    invisible(list(
        sessionFilename = sessionFilename,
        sessionFilepath = sessionFilepath,
        latestFilename = latestFilename,
        latestFilepath = latestFilepath
    ))
}

saveMetabNormExportMetadataFiles <- function(
    sessionData,
    sourceDir,
    saveRdsFn = saveRDS,
    logInfoFn = logger::log_info,
    logWarnFn = logger::log_warn
) {
    tryCatch({
        if (!is.null(sessionData$ruv_optimization_results) && length(sessionData$ruv_optimization_results) > 0) {
            saveRdsFn(
                sessionData$ruv_optimization_results,
                file.path(sourceDir, "metab_ruv_optimization_results.RDS")
            )
            logInfoFn("*** EXPORT: Saved metab_ruv_optimization_results.RDS ***")
        }

        if (!is.null(sessionData$itsd_selections) && length(sessionData$itsd_selections) > 0) {
            saveRdsFn(
                sessionData$itsd_selections,
                file.path(sourceDir, "metab_itsd_selections.RDS")
            )
            logInfoFn("*** EXPORT: Saved metab_itsd_selections.RDS ***")
        }

        if (!is.null(sessionData$qc_params)) {
            saveRdsFn(
                sessionData$qc_params,
                file.path(sourceDir, "metab_qc_params.RDS")
            )
            logInfoFn("*** EXPORT: Saved metab_qc_params.RDS ***")
        }
    }, error = function(e) {
        logWarnFn(sprintf("*** WARNING: Some metadata files could not be saved: %s ***", e$message))
    })

    invisible(NULL)
}

saveMetabNormExportSummaryFile <- function(
    sessionData,
    sourceDir,
    sessionFilename,
    writeLinesFn = writeLines,
    timeFn = Sys.time,
    formatTimeFn = format,
    logInfoFn = logger::log_info
) {
    ruvSummaryLines <- ""
    if (!is.null(sessionData$ruv_optimization_results)) {
        for (assayName in names(sessionData$ruv_optimization_results)) {
            result <- sessionData$ruv_optimization_results[[assayName]]
            if (!is.null(result) && isTRUE(result$success)) {
                ruvSummaryLines <- paste0(ruvSummaryLines, sprintf(
                    "\n  %s: k=%d, %%=%.1f, controls=%d"
                    , assayName
                    , result$best_k
                    , result$best_percentage
                    , sum(result$control_genes_index, na.rm = TRUE)
                ))
            }
        }
    }
    if (ruvSummaryLines == "") ruvSummaryLines <- "\n  (RUV skipped or not applied)"

    featureSummaryLines <- ""
    if (!is.null(sessionData$feature_counts)) {
        for (assayName in names(sessionData$feature_counts)) {
            counts <- sessionData$feature_counts[[assayName]]
            featureSummaryLines <- paste0(featureSummaryLines, sprintf(
                "\n  %s: %d features, %d samples"
                , assayName
                , counts$features
                , counts$samples
            ))
        }
    }

    summaryContent <- sprintf(
        "Metabolomics Normalized Session Data Export Summary\n===================================================\n\nExport Timestamp: %s\nSession File: %s\n\nData Summary:%s\n\nNormalization Parameters:\n- Method: %s\n- ITSD applied: %s\n- ITSD aggregation: %s\n- Log2 offset: %s\n- RUV mode: %s\n- RUV grouping variable: %s\n- Correlation threshold: %s\n\nRUV Optimization Results (per-assay):%s\n\nContrasts:\n%s\n\nThis data is ready for differential expression analysis.\nUse 'Load Filtered Session' in the DE tab to import.\n"
        , formatTimeFn(timeFn(), "%Y-%m-%d %H:%M:%S")
        , sessionFilename
        , featureSummaryLines
        , sessionData$normalization_method
        , ifelse(isTRUE(sessionData$itsd_applied), "Yes", "No")
        , ifelse(is.na(sessionData$itsd_aggregation), "N/A", sessionData$itsd_aggregation)
        , sessionData$log_offset
        , sessionData$ruv_mode
        , sessionData$ruv_grouping_variable
        , ifelse(is.null(sessionData$correlation_threshold), "N/A", sessionData$correlation_threshold)
        , ruvSummaryLines
        , if (!is.null(sessionData$contrasts_tbl)) paste(sessionData$contrasts_tbl$friendly_names, collapse = "\n") else "None defined"
    )

    summaryFilepath <- file.path(sourceDir, "metab_filtered_session_summary.txt")
    writeLinesFn(summaryContent, summaryFilepath)
    logInfoFn(sprintf("*** EXPORT: Summary saved to: %s ***", summaryFilepath))

    invisible(list(
        summaryContent = summaryContent,
        summaryFilepath = summaryFilepath
    ))
}

runMetabNormExportSessionWorkflow <- function(
    workflowData,
    normData,
    inputValues,
    experimentLabel,
    sourceDir,
    withProgressFn = shiny::withProgress,
    incProgressFn = shiny::incProgress,
    buildSessionDataFn = buildMetabNormExportSessionData,
    saveSessionRdsFilesFn = saveMetabNormExportSessionRdsFiles,
    saveMetadataFilesFn = saveMetabNormExportMetadataFiles,
    saveSummaryFileFn = saveMetabNormExportSummaryFile,
    logInfoFn = logger::log_info
) {
    withProgressFn(message = "Exporting normalized session data...", value = 0, expr = {
        incProgressFn(0.2, detail = "Gathering data...")

        sessionData <- buildSessionDataFn(
            workflowData = workflowData,
            normData = normData,
            inputValues = inputValues,
            experimentLabel = experimentLabel
        )

        logInfoFn("*** EXPORT: Gathered session data successfully ***")
        logInfoFn(sprintf("*** EXPORT: Assays: %s ***", paste(normData$assay_names, collapse = ", ")))
        logInfoFn(sprintf(
            "*** EXPORT: Contrasts available: %d ***",
            ifelse(is.null(sessionData$contrasts_tbl), 0, nrow(sessionData$contrasts_tbl))
        ))

        exportFiles <- saveSessionRdsFilesFn(
            sessionData = sessionData,
            sourceDir = sourceDir,
            incProgressFn = incProgressFn
        )

        incProgressFn(0.1, detail = "Saving metadata files...")
        saveMetadataFilesFn(
            sessionData = sessionData,
            sourceDir = sourceDir
        )

        incProgressFn(0.2, detail = "Creating summary...")
        saveSummaryFileFn(
            sessionData = sessionData,
            sourceDir = sourceDir,
            sessionFilename = exportFiles$sessionFilename
        )

        invisible(exportFiles)
    })
}

buildMetabNormCorrelationFilterSummary <- function(corrResults, filteredObject = NULL, originalObject = NULL) {
    if (is.null(corrResults) || length(corrResults) == 0) {
        return("No correlation results available.")
    }

    summaryLines <- c("=== Correlation Filtering Summary ===\n")

    for (assay_name in names(corrResults)) {
        assay_corr <- corrResults[[assay_name]]
        if (!is.null(assay_corr) && nrow(assay_corr) > 0) {
            n_pairs <- nrow(assay_corr)
            mean_corr <- round(mean(assay_corr$pearson_correlation, na.rm = TRUE), 3)
            min_corr <- round(min(assay_corr$pearson_correlation, na.rm = TRUE), 3)
            max_corr <- round(max(assay_corr$pearson_correlation, na.rm = TRUE), 3)

            summaryLines <- c(summaryLines, sprintf(
                "\n[%s]\n  Sample pairs: %d\n  Correlation: mean=%.3f, min=%.3f, max=%.3f",
                assay_name, n_pairs, mean_corr, min_corr, max_corr
            ))
        }
    }

    if (!is.null(originalObject) && !is.null(filteredObject)) {
        original_samples <- nrow(originalObject@design_matrix)
        filtered_samples <- nrow(filteredObject@design_matrix)
        removed <- original_samples - filtered_samples

        summaryLines <- c(summaryLines, sprintf(
            "\n\n[Sample Filtering]\n  Original: %d samples\n  After filtering: %d samples\n  Removed: %d samples",
            original_samples, filtered_samples, removed
        ))
    }

    paste(summaryLines, collapse = "")
}

resolveMetabNormFinalQcRenderState <- function(correlationFilteredObject = NULL, ruvCorrectedObject = NULL, postNormObject = NULL) {
    sourceObject <- if (!is.null(correlationFilteredObject)) {
        correlationFilteredObject
    } else if (!is.null(ruvCorrectedObject)) {
        ruvCorrectedObject
    } else {
        postNormObject
    }

    if (is.null(sourceObject)) {
        return(list(
            sourceObject = NULL
            , sourceStage = "empty"
            , isFallback = TRUE
            , plot = ggplot2::ggplot() +
                ggplot2::annotate("text", x = 0.5, y = 0.5, label = "No data", size = 6) +
                ggplot2::theme_void()
        ))
    }

    list(
        sourceObject = sourceObject
        , sourceStage = if (!is.null(correlationFilteredObject)) {
            "correlation_filter"
        } else if (!is.null(ruvCorrectedObject)) {
            "ruv_corrected"
        } else {
            "post_norm"
        }
        , isFallback = FALSE
        , plot = NULL
    )
}

buildMetabNormFinalQcPcaPlot <- function(sourceObject, colorVar = NULL, shapeVar = NULL, plotPcaFn = plotPca, wrapPlotsFn = patchwork::wrap_plots) {
    tryCatch({
        pcaPlots <- plotPcaFn(
            sourceObject
            , grouping_variable = colorVar
            , shape_variable = shapeVar
            , title = "Final QC - PCA"
        )

        if (is.list(pcaPlots) && length(pcaPlots) > 1) {
            wrapPlotsFn(pcaPlots, ncol = 1)
        } else if (is.list(pcaPlots) && length(pcaPlots) == 1) {
            pcaPlots[[1]]
        } else if (inherits(pcaPlots, "ggplot")) {
            pcaPlots
        } else {
            ggplot2::ggplot() + ggplot2::theme_void()
        }
    }, error = function(e) {
        ggplot2::ggplot() +
            ggplot2::annotate("text", x = 0.5, y = 0.5, label = paste("Error:", e$message), size = 4) +
            ggplot2::theme_void()
    })
}

runMetabNormApplyCorrelationObserverShell <- function(
    workflowData,
    normData,
    threshold,
    groupingVariable,
    addLogFn = function(entry) invisible(entry),
    showNotificationFn = shiny::showNotification,
    removeNotificationFn = shiny::removeNotification,
    reqFn = shiny::req,
    runObserverEntryFn = runMetabNormApplyCorrelationObserverEntry
) {
    reqFn(workflowData$state_manager)
    reqFn(normData$ruv_complete || normData$normalization_complete)

    dispatchState <- runObserverEntryFn(
        workflowData = workflowData,
        normData = normData,
        threshold = threshold,
        groupingVariable = groupingVariable,
        addLogFn = addLogFn,
        showNotificationFn = showNotificationFn,
        removeNotificationFn = removeNotificationFn
    )

    invisible(dispatchState)
}

