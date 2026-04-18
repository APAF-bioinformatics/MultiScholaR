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

# Keep the composite image-assembly helper cluster top-level so later waves can
# move this seam without reopening the module wrapper body.
buildMetabNormLabelPlot <- function(title) {
    ggplot2::ggplot() +
        ggplot2::annotate("text", x = 0, y = 0.5, label = title, size = 5, hjust = 0) +
        ggplot2::xlim(0, 1) +
        ggplot2::theme_void() +
        ggplot2::theme(
            plot.margin = ggplot2::margin(5, 5, 5, 5),
            panel.background = ggplot2::element_blank()
        )
}

buildMetabNormTitlePlot <- function(title) {
    ggplot2::ggplot() +
        ggplot2::annotate("text", x = 0.5, y = 0.5, label = title, size = 6, fontface = "bold", hjust = 0.5) +
        ggplot2::xlim(0, 1) +
        ggplot2::theme_void() +
        ggplot2::theme(
            plot.margin = ggplot2::margin(5, 5, 10, 5),
            panel.background = ggplot2::element_blank()
        )
}

loadMetabNormImageAsPlot <- function(
    file_path,
    fileExistsFn = file.exists,
    readPngFn = png::readPNG,
    rasterGrobFn = grid::rasterGrob,
    logWarnFn = logger::log_warn
) {
    if (is.na(file_path) || !fileExistsFn(file_path)) {
        return(ggplot2::ggplot() + ggplot2::theme_void())
    }

    tryCatch({
        img <- readPngFn(file_path)
        grob <- rasterGrobFn(img, interpolate = TRUE)
        ggplot2::ggplot() +
            ggplot2::annotation_custom(grob, xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf) +
            ggplot2::theme_void()
    }, error = function(e) {
        logWarnFn(sprintf("[generateCompositeFromFiles] Could not load image: %s", file_path))
        ggplot2::ggplot() + ggplot2::theme_void()
    })
}

generateMetabNormCompositeFromFiles <- function(
    plot_files,
    ncol = 3,
    row_labels = NULL,
    column_labels = NULL,
    logInfoFn = logger::log_info,
    logErrorFn = logger::log_error,
    warningFn = warning,
    requireNamespaceFn = requireNamespace,
    buildLabelPlotFn = buildMetabNormLabelPlot,
    buildTitlePlotFn = buildMetabNormTitlePlot,
    loadImageAsPlotFn = loadMetabNormImageAsPlot,
    wrapPlotsFn = patchwork::wrap_plots,
    plotLayoutFn = patchwork::plot_layout,
    combineLayoutFn = function(plot, layout) plot + layout,
    fileExistsFn = file.exists,
    gcFn = gc
) {
    logInfoFn(sprintf("[generateCompositeFromFiles] Generating composite from %d files...", length(plot_files)))

    if (!requireNamespaceFn("patchwork", quietly = TRUE)) {
        warningFn("patchwork package required for composite generation")
        return(NULL)
    }
    if (!requireNamespaceFn("ggplot2", quietly = TRUE)) {
        warningFn("ggplot2 package required for composite generation")
        return(NULL)
    }
    if (!requireNamespaceFn("png", quietly = TRUE)) {
        warningFn("png package required for composite generation")
        return(NULL)
    }

    tryCatch({
        n_files <- length(plot_files)
        n_plot_types <- n_files / ncol

        if (is.null(row_labels)) {
            all_labels <- letters[1:n_files]
            row_labels <- split(paste0(all_labels, ")"), rep(1:n_plot_types, each = ncol))
            names(row_labels) <- paste0("row", seq_len(n_plot_types))
        }

        plot_sections <- list()
        height_values <- c()

        if (!is.null(column_labels) && length(column_labels) == ncol) {
            title_plots <- lapply(column_labels, buildTitlePlotFn)
            plot_sections <- append(plot_sections, list(
                wrapPlotsFn(title_plots, ncol = ncol)
            ))
            height_values <- c(height_values, 0.2)
            logInfoFn("[generateCompositeFromFiles] Added column titles")
        }

        row_names <- names(row_labels)

        for (i in seq_along(row_names)) {
            row_name <- row_names[i]
            labels <- row_labels[[row_name]]

            start_idx <- (i - 1) * ncol + 1
            end_idx <- min(i * ncol, n_files)
            row_files <- plot_files[start_idx:end_idx]

            has_files <- any(!is.na(row_files) & vapply(
                row_files,
                function(path) !is.na(path) && fileExistsFn(path),
                logical(1)
            ))

            if (has_files) {
                label_plots <- lapply(labels, buildLabelPlotFn)
                image_plots <- lapply(row_files, loadImageAsPlotFn)

                plot_sections <- append(plot_sections, list(
                    wrapPlotsFn(label_plots, ncol = ncol),
                    wrapPlotsFn(image_plots, ncol = ncol)
                ))
                height_values <- c(height_values, 0.1, 1)

                logInfoFn(sprintf("[generateCompositeFromFiles] Added row: %s", row_name))
            } else {
                logInfoFn(sprintf("[generateCompositeFromFiles] Skipping empty row: %s", row_name))
            }
        }

        if (length(plot_sections) == 0) {
            warningFn("No valid plot sections to combine")
            return(NULL)
        }

        logInfoFn("[generateCompositeFromFiles] Combining plot sections...")
        combined_plot <- combineLayoutFn(
            wrapPlotsFn(plot_sections, ncol = 1),
            plotLayoutFn(heights = height_values)
        )

        plot_width <- 4 + (ncol * 3)
        plot_height <- 4 + (length(height_values) * 2)

        rm(plot_sections)
        gcFn()

        list(plot = combined_plot, width = plot_width, height = plot_height)
    }, error = function(e) {
        logErrorFn(paste("[generateCompositeFromFiles] Error:", e$message))
        NULL
    })
}

# Keep the auto pre-normalization QC helper top-level so a later exact-source
# wave can move it without reopening the selected-tab observer body.
generateMetabNormPreNormalizationQc <- function(
    workflowData,
    experimentPaths,
    normData,
    getPlotAestheticsFn,
    addLogFn = function(entry) invisible(entry),
    reqFn = shiny::req,
    generateMetabQcPlotsFn = generateMetabQcPlots,
    logInfoFn = logger::log_info,
    logWarnFn = logger::log_warn,
    logErrorFn = logger::log_error
) {
    logInfoFn("=== GENERATING PRE-NORMALIZATION QC PLOTS ===")

    reqFn(workflowData$state_manager)
    current_s4 <- workflowData$state_manager$getState()

    if (is.null(current_s4)) {
        logWarnFn("No S4 object available for QC plot generation")
        return()
    }

    if (inherits(current_s4, "MetaboliteAssayData")) {
        detected_assays <- names(current_s4@metabolite_data)
        if (length(detected_assays) > 0 && is.null(normData$assay_names)) {
            normData$assay_names <- detected_assays
            logInfoFn(paste("Set assay names:", paste(detected_assays, collapse = ", ")))
        }
    }

    aesthetics <- getPlotAestheticsFn()

    tryCatch({
        generateMetabQcPlotsFn(
            theObject = current_s4
            , experiment_paths = experimentPaths
            , stage = "post_filter"
            , grouping_variable = aesthetics$color_var
            , shape_variable = aesthetics$shape_var
        )

        normData$plot_refresh_trigger <- normData$plot_refresh_trigger + 1
        normData$pre_norm_qc_generated <- TRUE
        logInfoFn("Pre-normalization QC plots generated successfully")

    }, error = function(e) {
        logErrorFn(paste("Error generating pre-normalization QC:", e$message))
        addLogFn(paste("Error generating Pre-QC:", e$message))
    })
}

# Keep plot-aesthetics resolution top-level so a later extraction wave can move
# it without reopening the module server wrapper.
getPlotAesthetics <- function(
    colorVariable = NULL,
    shapeVariable = NULL
) {
    list(
        color_var = if (is.null(colorVariable) || colorVariable == "") {
            "group"
        } else {
            colorVariable
        }
        , shape_var = if (is.null(shapeVariable) || shapeVariable == "") {
            "group"
        } else {
            shapeVariable
        }
    )
}

# Keep assay-label rendering top-level so a later extraction wave can move it
# without reopening the module server wrapper's static label bindings.
renderMetabNormAssayLabel <- function(
    assaySlot,
    getAssayNamesFn,
    renderTextFn = shiny::renderText
) {
    renderTextFn({
        assayNames <- getAssayNamesFn()

        if (!is.null(assayNames) && length(assayNames) >= assaySlot) {
            paste("Assay:", assayNames[[assaySlot]])
        } else {
            paste0("Assay ", assaySlot, ": (detecting...)")
        }
    })
}

# Keep QC image rendering top-level so a later extraction wave can move it
# without reopening the module server wrapper's static QC output bindings.
renderMetabNormQcImageForAssay <- function(
    assaySlot,
    plotType,
    stagePrefix,
    normData,
    qcDir,
    renderImageFn = shiny::renderImage,
    fileExistsFn = file.exists,
    filePathFn = file.path,
    sanitizeAssayNameFn = function(assayName) {
        gsub("[^A-Za-z0-9]", "_", tolower(assayName))
    }
) {
    renderImageFn({
        normData$plot_refresh_trigger
        assayNames <- normData$assay_names

        if (is.null(assayNames) || length(assayNames) < assaySlot) {
            return(list(src = "", alt = "Assay not detected yet"))
        }

        assayName <- assayNames[[assaySlot]]
        safeName <- sanitizeAssayNameFn(assayName)
        filename <- paste0(safeName, "_", stagePrefix, "_", plotType, ".png")

        if (is.null(qcDir)) {
            return(list(src = "", alt = "QC directory not configured"))
        }

        imgPath <- filePathFn(qcDir, filename)

        if (fileExistsFn(imgPath)) {
            list(
                src = imgPath
                , contentType = "image/png"
                , width = "100%"
                , height = "auto"
                , alt = paste(plotType, "-", assayName)
            )
        } else {
            list(src = "", alt = paste("Plot not generated yet:", filename))
        }
    }, deleteFile = FALSE)
}

# Keep normalization-log mutation top-level so a later extraction wave can move
# it without reopening the module server wrapper's log state shell.
appendMetabNormNormalizationLog <- function(
    normData,
    message,
    timestampFn = function() {
        format(Sys.time(), "%H:%M:%S")
    }
) {
    timestamp <- timestampFn()
    updatedLog <- c(
        normData$normalization_log
        , sprintf("[%s] %s", timestamp, message)
    )
    normData$normalization_log <- updatedLog

    invisible(updatedLog)
}

# Keep normalization-log rendering top-level so a later extraction wave can move
# it without reopening the module server wrapper's output binding shell.
renderMetabNormNormalizationLog <- function(
    normData,
    renderTextFn = shiny::renderText,
    emptyMessage = "Normalization log will appear here as you apply steps..."
) {
    renderTextFn({
        normalizationLog <- normData$normalization_log

        if (length(normalizationLog) == 0) {
            return(emptyMessage)
        }

        paste(normalizationLog, collapse = "\n")
    })
}

# Keep correlation-filter summary rendering top-level so a later extraction
# wave can move it without reopening the module server wrapper's output shell.
renderMetabNormCorrelationFilterSummary <- function(
    normData,
    renderTextFn = shiny::renderText,
    buildSummaryFn = buildMetabNormCorrelationFilterSummary,
    incompleteMessage = "Apply correlation filter to see results..."
) {
    renderTextFn({
        if (!normData$correlation_filtering_complete) {
            return(incompleteMessage)
        }

        corrResults <- normData$correlation_results
        filteredObject <- normData$correlation_filtered_obj
        originalObject <- if (!is.null(normData$ruv_corrected_obj)) {
            normData$ruv_corrected_obj
        } else {
            normData$post_norm_obj
        }

        buildSummaryFn(
            corrResults = corrResults
            , filteredObject = filteredObject
            , originalObject = originalObject
        )
    })
}

# Keep final-QC rendering top-level so a later extraction wave can move it
# without reopening the module server wrapper's output shell.
renderMetabNormFinalQcPlot <- function(
    normData,
    colorVariableFn = function() NULL,
    shapeVariableFn = function() NULL,
    renderPlotFn = shiny::renderPlot,
    reqFn = shiny::req,
    resolveRenderStateFn = resolveMetabNormFinalQcRenderState,
    getPlotAestheticsFn = getPlotAesthetics,
    buildPcaPlotFn = buildMetabNormFinalQcPcaPlot
) {
    renderPlotFn({
        reqFn(normData$correlation_filtering_complete || normData$ruv_complete)

        finalQcState <- resolveRenderStateFn(
            correlationFilteredObject = normData$correlation_filtered_obj
            , ruvCorrectedObject = normData$ruv_corrected_obj
            , postNormObject = normData$post_norm_obj
        )

        if (isTRUE(finalQcState$isFallback)) {
            return(finalQcState$plot)
        }

        aesthetics <- getPlotAestheticsFn(
            colorVariable = colorVariableFn()
            , shapeVariable = shapeVariableFn()
        )

        buildPcaPlotFn(
            sourceObject = finalQcState$sourceObject
            , colorVar = aesthetics$color_var
            , shapeVar = aesthetics$shape_var
        )
    })
}

# Keep the reset-normalization observer wrapper top-level so a later
# extraction wave can move it without reopening the module server wrapper.
runMetabNormResetNormalizationObserverWrapper <- function(
    workflowData,
    normData,
    addLogFn = function(entry) invisible(entry),
    showNotificationFn = shiny::showNotification,
    reqFn = shiny::req,
    runObserverShellFn = runMetabNormResetNormalizationObserverShell
) {
    runObserverShellFn(
        workflowData = workflowData
        , normData = normData
        , addLogFn = addLogFn
        , showNotificationFn = showNotificationFn
        , reqFn = reqFn
    )
}

# Keep the apply-correlation observer wrapper top-level so a later extraction
# wave can move it without reopening the module server wrapper.
runMetabNormApplyCorrelationObserverWrapper <- function(
    workflowData,
    input,
    normData,
    addLogFn = function(entry) invisible(entry),
    showNotificationFn = shiny::showNotification,
    removeNotificationFn = shiny::removeNotification,
    reqFn = shiny::req,
    runObserverShellFn = runMetabNormApplyCorrelationObserverShell,
    runObserverEntryFn = runMetabNormApplyCorrelationObserverEntry
) {
    runObserverShellFn(
        workflowData = workflowData
        , normData = normData
        , threshold = input$min_pearson_correlation_threshold
        , groupingVariable = input$ruv_grouping_variable
        , addLogFn = addLogFn
        , showNotificationFn = showNotificationFn
        , removeNotificationFn = removeNotificationFn
        , reqFn = reqFn
        , runObserverEntryFn = runObserverEntryFn
    )
}

# Keep the skip-correlation observer wrapper top-level so a later extraction
# wave can move it without reopening the module server wrapper.
runMetabNormSkipCorrelationObserverWrapper <- function(
    workflowData,
    normData,
    addLogFn = function(entry) invisible(entry),
    showNotificationFn = shiny::showNotification,
    reqFn = shiny::req,
    runObserverShellFn = runMetabNormSkipCorrelationObserverShell,
    runObserverEntryFn = runMetabNormSkipCorrelationObserverEntry
) {
    runObserverShellFn(
        workflowData = workflowData
        , normData = normData
        , addLogFn = addLogFn
        , showNotificationFn = showNotificationFn
        , reqFn = reqFn
        , runObserverEntryFn = runObserverEntryFn
    )
}

# Keep the export-session observer wrapper top-level so a later extraction
# wave can move it without reopening the module server wrapper.
runMetabNormExportSessionObserverWrapper <- function(
    workflowData,
    input,
    normData,
    experimentPaths,
    experimentLabel,
    addLogFn = function(entry) invisible(entry),
    logInfoFn = logger::log_info,
    reqFn = shiny::req,
    runObserverShellFn = runMetabNormExportSessionObserverShell,
    checkReadyFn = checkMetabNormExportSessionReady,
    dispatchExportSessionFn = dispatchMetabNormExportSession
) {
    inputValues <- list(
        export_session = input$export_session
        , norm_method = input$norm_method
        , ruv_mode = input$ruv_mode
        , apply_itsd = input$apply_itsd
        , itsd_aggregation = input$itsd_aggregation
        , log_offset = input$log_offset
        , min_pearson_correlation_threshold = input$min_pearson_correlation_threshold
        , ruv_grouping_variable = input$ruv_grouping_variable
    )

    runObserverShellFn(
        workflowData = workflowData
        , normData = normData
        , inputValues = inputValues
        , experimentPaths = experimentPaths
        , experimentLabel = experimentLabel
        , addLogFn = addLogFn
        , logInfoFn = logInfoFn
        , reqFn = reqFn
        , checkReadyFn = checkReadyFn
        , dispatchExportSessionFn = dispatchExportSessionFn
    )
}

# Keep design-driven input-choice updates top-level so a later extraction wave
# can move them without reopening the module server wrapper.
updateMetabNormDesignDrivenChoices <- function(
    session,
    designMatrix,
    updateSelectInputFn = shiny::updateSelectInput
) {
    if (is.null(designMatrix)) {
        return(invisible(NULL))
    }

    designCols <- colnames(designMatrix)

    plotAvailableVars <- intersect(
        designCols,
        c("group", "factor1", "factor2", "batch", "technical_replicate_id", "sample_id")
    )

    if (length(plotAvailableVars) > 0) {
        defaultPlotVar <- if ("group" %in% plotAvailableVars) {
            "group"
        } else {
            plotAvailableVars[[1]]
        }

        updateSelectInputFn(
            session, "color_variable"
            , choices = plotAvailableVars
            , selected = defaultPlotVar
        )
        updateSelectInputFn(
            session, "shape_variable"
            , choices = plotAvailableVars
            , selected = defaultPlotVar
        )
    }

    ruvAvailableVars <- intersect(designCols, c("group", "factor1", "factor2", "batch"))

    if (length(ruvAvailableVars) > 0) {
        defaultRuvVar <- if ("group" %in% ruvAvailableVars) {
            "group"
        } else {
            ruvAvailableVars[[1]]
        }

        updateSelectInputFn(
            session, "ruv_grouping_variable"
            , choices = ruvAvailableVars
            , selected = defaultRuvVar
        )
    }

    invisible(NULL)
}

# Keep assay-name detection top-level so a later extraction wave can move the
# observer shell without reopening the module server wrapper.
initializeMetabNormAssayNames <- function(
    stateManager,
    normData,
    reqFn = shiny::req,
    getStateFn = function(manager) manager$getState(),
    logInfoFn = logger::log_info,
    logWarnFn = logger::log_warn
) {
    reqFn(stateManager)

    tryCatch({
        currentS4 <- getStateFn(stateManager)

        if (!inherits(currentS4, "MetaboliteAssayData")) {
            return(invisible(NULL))
        }

        detectedAssays <- names(currentS4@metabolite_data)
        normData$assay_names <- detectedAssays
        logInfoFn(paste("Detected assays:", paste(detectedAssays, collapse = ", ")))

        for (assayName in normData$assay_names) {
            if (!assayName %in% names(normData$itsd_selections)) {
                normData$itsd_selections[[assayName]] <- NULL
            }
        }

        invisible(detectedAssays)
    }, error = function(e) {
        logWarnFn(paste("Could not detect assay names:", e$message))
        invisible(NULL)
    })
}

# Keep the selected-tab auto pre-normalization QC shell top-level so a later
# extraction wave can move it without reopening the module server wrapper.
runMetabNormAutoPreNormalizationQcObserverShell <- function(
    selectedTab,
    workflowData,
    experimentPaths,
    normData,
    colorVariable = NULL,
    shapeVariable = NULL,
    addLogFn = function(entry) invisible(entry),
    reqFn = shiny::req,
    getStateFn = function(manager) manager$getState(),
    withProgressFn = shiny::withProgress,
    generatePreNormalizationQcFn = generateMetabNormPreNormalizationQc,
    getPlotAestheticsFn = getPlotAesthetics,
    generateMetabQcPlotsFn = generateMetabQcPlots,
    logInfoFn = logger::log_info,
    logWarnFn = logger::log_warn,
    logErrorFn = logger::log_error
) {
    if (is.null(selectedTab) || selectedTab != "norm") {
        return(invisible(NULL))
    }

    logInfoFn("Normalization tab selected - checking if pre-QC needed")

    if (isTRUE(normData$pre_norm_qc_generated)) {
        return(invisible(NULL))
    }

    reqFn(workflowData$state_manager)
    currentS4 <- tryCatch(getStateFn(workflowData$state_manager), error = function(e) NULL)

    if (is.null(currentS4) || !inherits(currentS4, "MetaboliteAssayData")) {
        return(invisible(NULL))
    }

    logInfoFn("Auto-triggering pre-normalization QC plots")

    withProgressFn(
        message = "Generating Pre-Normalization QC..."
        , value = 0.5
        , {
            generatePreNormalizationQcFn(
                workflowData = workflowData
                , experimentPaths = experimentPaths
                , normData = normData
                , getPlotAestheticsFn = function() {
                    getPlotAestheticsFn(
                        colorVariable = colorVariable
                        , shapeVariable = shapeVariable
                    )
                }
                , addLogFn = addLogFn
                , reqFn = reqFn
                , generateMetabQcPlotsFn = generateMetabQcPlotsFn
                , logInfoFn = logInfoFn
                , logWarnFn = logWarnFn
                , logErrorFn = logErrorFn
            )
        }
    )

    normData$pre_norm_qc_generated <- TRUE

    invisible(list(
        selectedTab = selectedTab
        , currentS4 = currentS4
        , preNormQcGenerated = normData$pre_norm_qc_generated
    ))
}

# Keep the dynamic ITSD selection tables render shell top-level so a later
# extraction wave can move it without reopening the module server wrapper.
renderMetabNormItsdSelectionUi <- function(
    normData,
    ns,
    renderUIFn = shiny::renderUI,
    reqFn = shiny::req,
    mapFn = purrr::map,
    wellPanelFn = shiny::wellPanel,
    h5Fn = shiny::h5,
    dataTableOutputFn = DT::dataTableOutput,
    brFn = shiny::br,
    tagListFn = shiny::tagList
) {
    renderUIFn({
        reqFn(normData$assay_names)

        assayUis <- mapFn(normData$assay_names, \(assayName) {
            safeName <- gsub("[^A-Za-z0-9]", "_", tolower(assayName))
            wellPanelFn(
                h5Fn(paste("Assay:", assayName))
                , dataTableOutputFn(ns(paste0("itsd_table_", safeName)))
                , brFn()
            )
        })

        tagListFn(assayUis)
    })
}

# Keep the per-assay ITSD selection table render shell top-level so a later
# extraction wave can move it without reopening the module server wrapper.
renderMetabNormItsdSelectionTable <- function(
    assayName,
    currentS4,
    renderDataTableFn = DT::renderDataTable,
    buildItsdSelectionTableFn = buildItsdSelectionTable,
    datatableFn = DT::datatable,
    formatStyleFn = DT::formatStyle,
    styleEqualFn = DT::styleEqual,
    formatRoundFn = DT::formatRound,
    datatableOptions = list(
        pageLength = 10,
        scrollX = TRUE,
        order = list(list(4, "desc"), list(3, "asc"))
    ),
    candidateHighlight = "#d4edda",
    roundedColumns = c("mean_intensity", "cv_percent"),
    digits = 2,
    filter = "top"
) {
    renderDataTableFn({
        assayData <- currentS4@metabolite_data[[assayName]]
        if (is.null(assayData)) {
            return(NULL)
        }

        metaboliteIdCol <- currentS4@metabolite_id_column
        annotationCol <- currentS4@annotation_id_column

        selectionTable <- buildItsdSelectionTableFn(
            assay_data = assayData,
            metabolite_id_col = metaboliteIdCol,
            annotation_cols = annotationCol
        )

        preselected <- which(selectionTable$is_candidate)

        datatableFn(
            selectionTable,
            selection = list(mode = "multiple", selected = preselected),
            filter = filter,
            options = datatableOptions,
            rownames = FALSE
        ) |>
            formatStyleFn(
                "is_candidate",
                backgroundColor = styleEqualFn(TRUE, candidateHighlight)
            ) |>
            formatRoundFn(columns = roundedColumns, digits = digits)
    })
}

# Keep the per-assay ITSD selection tracking registration top-level so a later
# extraction wave can move it without reopening the module server wrapper.
registerMetabNormItsdSelectionTracking <- function(
    normData,
    input,
    walkFn = purrr::walk,
    observeEventFn = shiny::observeEvent,
    selectionGetter = function(input, inputId) input[[inputId]],
    safeNameFn = function(assayName) gsub("[^A-Za-z0-9]", "_", tolower(assayName)),
    logInfoFn = logger::log_info
) {
    walkFn(normData$assay_names, \(assayName) {
        safeName <- safeNameFn(assayName)
        inputId <- paste0("itsd_table_", safeName, "_rows_selected")

        observeEventFn(selectionGetter(input, inputId), {
            selectedRows <- selectionGetter(input, inputId)
            normData$itsd_selections[[assayName]] <- selectedRows
            logInfoFn(paste(
                "ITSD selection updated for", assayName, ":"
                , length(selectedRows), "features selected"
            ))
        }, ignoreNULL = FALSE)
    })

    invisible(NULL)
}

# Keep the ITSD selection tracking observer shell top-level so a later
# extraction wave can move it without reopening the module server wrapper.
runMetabNormItsdSelectionTrackingObserverShell <- function(
    normData,
    input,
    reqFn = shiny::req,
    registerItsdSelectionTrackingFn = registerMetabNormItsdSelectionTracking
) {
    reqFn(normData$assay_names)

    registerItsdSelectionTrackingFn(
        normData = normData,
        input = input
    )

    invisible(NULL)
}

# Keep the per-assay ITSD selection table binding observer shell top-level so a
# later extraction wave can move it without reopening the module server wrapper.
runMetabNormItsdSelectionTableObserverShell <- function(
    normData,
    workflowData,
    output,
    reqFn = shiny::req,
    getStateFn = function(stateManager) stateManager$getState(),
    walkFn = purrr::walk,
    safeNameFn = function(assayName) gsub("[^A-Za-z0-9]", "_", tolower(assayName)),
    renderItsdSelectionTableFn = renderMetabNormItsdSelectionTable
) {
    reqFn(normData$assay_names)
    reqFn(workflowData$state_manager)

    currentS4 <- tryCatch({
        getStateFn(workflowData$state_manager)
    }, error = function(e) NULL)

    if (is.null(currentS4) || !inherits(currentS4, "MetaboliteAssayData")) {
        return(invisible(NULL))
    }

    walkFn(normData$assay_names, \(assayName) {
        safeName <- safeNameFn(assayName)
        outputId <- paste0("itsd_table_", safeName)

        output[[outputId]] <- renderItsdSelectionTableFn(
            assayName = assayName,
            currentS4 = currentS4
        )
    })

    invisible(NULL)
}

# Keep the dynamic RUV QC plots render shell top-level so a later extraction
# wave can move it without reopening the module server wrapper.
renderMetabNormRuvQcUi <- function(
    normData,
    ns,
    renderUIFn = shiny::renderUI,
    reqFn = shiny::req,
    mapFn = purrr::map,
    tagListFn = shiny::tagList,
    fluidRowFn = shiny::fluidRow,
    columnFn = shiny::column,
    h5Fn = shiny::h5,
    h6Fn = shiny::h6,
    jquiResizableFn = shinyjqui::jqui_resizable,
    plotOutputFn = shiny::plotOutput,
    wellPanelFn = shiny::wellPanel,
    verbatimTextOutputFn = shiny::verbatimTextOutput,
    brFn = shiny::br,
    dataTableOutputFn = DT::dataTableOutput,
    hrFn = shiny::hr
) {
    renderUIFn({
        reqFn(normData$assay_names)

        assayUis <- mapFn(normData$assay_names, \(assayName) {
            safeName <- gsub("[^A-Za-z0-9]", "_", tolower(assayName))

            tagListFn(
                fluidRowFn(
                    columnFn(
                        12,
                        h5Fn(
                            paste("Assay:", assayName),
                            style = "border-bottom: 1px solid #ddd; padding-bottom: 5px;"
                        )
                    )
                ),
                fluidRowFn(
                    columnFn(
                        8,
                        jquiResizableFn(
                            plotOutputFn(
                                ns(paste0("cancor_plot_", safeName)),
                                height = "400px"
                            )
                        )
                    ),
                    columnFn(
                        4,
                        wellPanelFn(
                            h6Fn("Optimization Summary"),
                            verbatimTextOutputFn(ns(paste0("ruv_summary_", safeName))),
                            brFn(),
                            h6Fn("Results Table"),
                            dataTableOutputFn(ns(paste0("ruv_table_", safeName)))
                        )
                    )
                ),
                hrFn()
            )
        })

        do.call(tagListFn, assayUis)
    })
}

# Keep the per-assay RUV cancor plot render shell top-level so a later
# extraction wave can move it without reopening the module server wrapper.
renderMetabNormRuvCancorPlot <- function(
    assayName,
    normData,
    renderPlotFn = shiny::renderPlot,
    buildFallbackPlotFn = function() {
        ggplot2::ggplot() +
            ggplot2::annotate("text", x = 0.5, y = 0.5, label = "No data", size = 6) +
            ggplot2::theme_void()
    }
) {
    renderPlotFn({
        result <- normData$ruv_optimization_results[[assayName]]
        if (!is.null(result) && isTRUE(result$success) && !is.null(result$cancor_plot)) {
            result$cancor_plot
        } else {
            buildFallbackPlotFn()
        }
    })
}

# Keep the per-assay RUV optimization summary render shell top-level so a later
# extraction wave can move it without reopening the module server wrapper.
renderMetabNormRuvOptimizationSummary <- function(
    assayName,
    normData,
    renderTextFn = shiny::renderText,
    pendingMessage = "Not yet computed"
) {
    renderTextFn({
        result <- normData$ruv_optimization_results[[assayName]]

        if (!is.null(result) && isTRUE(result$success)) {
            return(paste0(
                "Best k: ", result$best_k, "\n"
                , "Best %: ", result$best_percentage, "\n"
                , "Separation: ", round(result$separation_score, 4), "\n"
                , "Controls: ", sum(result$control_genes_index, na.rm = TRUE)
            ))
        }

        if (!is.null(result)) {
            return(paste0("Failed: ", result$error))
        }

        pendingMessage
    })
}

# Keep the per-assay RUV results table render shell top-level so a later
# extraction wave can move it without reopening the module server wrapper.
renderMetabNormRuvResultsTable <- function(
    assayName,
    normData,
    renderDataTableFn = DT::renderDataTable,
    datatableFn = DT::datatable,
    datatableOptions = list(pageLength = 5, dom = "t")
) {
    renderDataTableFn({
        result <- normData$ruv_optimization_results[[assayName]]

        if (!is.null(result) && isTRUE(result$success) && !is.null(result$optimization_results)) {
            return(datatableFn(
                result$optimization_results
                , options = datatableOptions
                , rownames = FALSE
            ))
        }

        NULL
    })
}

# Keep the per-assay RUV plot/summary/table binding observer shell top-level
# so a later extraction wave can move it without reopening the module server
# wrapper.
runMetabNormRuvBindingObserverShell <- function(
    normData,
    output,
    reqFn = shiny::req,
    walkFn = purrr::walk,
    safeNameFn = function(assayName) gsub("[^A-Za-z0-9]", "_", tolower(assayName)),
    renderRuvCancorPlotFn = renderMetabNormRuvCancorPlot,
    renderRuvOptimizationSummaryFn = renderMetabNormRuvOptimizationSummary,
    renderRuvResultsTableFn = renderMetabNormRuvResultsTable
) {
    reqFn(normData$assay_names)
    reqFn(length(normData$ruv_optimization_results) > 0)

    walkFn(normData$assay_names, \(assayName) {
        safeName <- safeNameFn(assayName)

        output[[paste0("cancor_plot_", safeName)]] <- renderRuvCancorPlotFn(
            assayName = assayName,
            normData = normData
        )

        output[[paste0("ruv_summary_", safeName)]] <- renderRuvOptimizationSummaryFn(
            assayName = assayName,
            normData = normData
        )

        output[[paste0("ruv_table_", safeName)]] <- renderRuvResultsTableFn(
            assayName = assayName,
            normData = normData
        )
    })

    invisible(NULL)
}

# Keep the static QC image output binding shell top-level so a later extraction
# wave can move it without reopening the module server wrapper.
runMetabNormQcImageBindingShell <- function(
    output,
    normData,
    qcDir,
    renderQcImageFn = renderMetabNormQcImageForAssay
) {
    output$pca_post_filter_assay1 <- renderQcImageFn(1, "pca", "pre_norm", normData, qcDir)
    output$pca_post_norm_assay1 <- renderQcImageFn(1, "pca", "post_norm", normData, qcDir)
    output$pca_ruv_corrected_assay1 <- renderQcImageFn(1, "pca", "ruv_corrected", normData, qcDir)
    output$pca_post_filter_assay2 <- renderQcImageFn(2, "pca", "pre_norm", normData, qcDir)
    output$pca_post_norm_assay2 <- renderQcImageFn(2, "pca", "post_norm", normData, qcDir)
    output$pca_ruv_corrected_assay2 <- renderQcImageFn(2, "pca", "ruv_corrected", normData, qcDir)

    output$density_post_filter_assay1 <- renderQcImageFn(1, "density", "pre_norm", normData, qcDir)
    output$density_post_norm_assay1 <- renderQcImageFn(1, "density", "post_norm", normData, qcDir)
    output$density_ruv_corrected_assay1 <- renderQcImageFn(1, "density", "ruv_corrected", normData, qcDir)
    output$density_post_filter_assay2 <- renderQcImageFn(2, "density", "pre_norm", normData, qcDir)
    output$density_post_norm_assay2 <- renderQcImageFn(2, "density", "post_norm", normData, qcDir)
    output$density_ruv_corrected_assay2 <- renderQcImageFn(2, "density", "ruv_corrected", normData, qcDir)

    output$rle_post_filter_assay1 <- renderQcImageFn(1, "rle", "pre_norm", normData, qcDir)
    output$rle_post_norm_assay1 <- renderQcImageFn(1, "rle", "post_norm", normData, qcDir)
    output$rle_ruv_corrected_assay1 <- renderQcImageFn(1, "rle", "ruv_corrected", normData, qcDir)
    output$rle_post_filter_assay2 <- renderQcImageFn(2, "rle", "pre_norm", normData, qcDir)
    output$rle_post_norm_assay2 <- renderQcImageFn(2, "rle", "post_norm", normData, qcDir)
    output$rle_ruv_corrected_assay2 <- renderQcImageFn(2, "rle", "ruv_corrected", normData, qcDir)

    output$correlation_post_filter_assay1 <- renderQcImageFn(1, "correlation", "pre_norm", normData, qcDir)
    output$correlation_post_norm_assay1 <- renderQcImageFn(1, "correlation", "post_norm", normData, qcDir)
    output$correlation_ruv_corrected_assay1 <- renderQcImageFn(1, "correlation", "ruv_corrected", normData, qcDir)
    output$correlation_post_filter_assay2 <- renderQcImageFn(2, "correlation", "pre_norm", normData, qcDir)
    output$correlation_post_norm_assay2 <- renderQcImageFn(2, "correlation", "post_norm", normData, qcDir)
    output$correlation_ruv_corrected_assay2 <- renderQcImageFn(2, "correlation", "ruv_corrected", normData, qcDir)

    invisible(NULL)
}

# Keep the static assay-label output binding shell top-level so a later
# extraction wave can move it without reopening the module server wrapper.
runMetabNormAssayLabelBindingShell <- function(
    output,
    getAssayNamesFn,
    renderAssayLabelFn = renderMetabNormAssayLabel
) {
    output$assay1_label_pca <- renderAssayLabelFn(1, getAssayNamesFn)
    output$assay2_label_pca <- renderAssayLabelFn(2, getAssayNamesFn)

    output$assay1_label_density <- renderAssayLabelFn(1, getAssayNamesFn)
    output$assay2_label_density <- renderAssayLabelFn(2, getAssayNamesFn)

    output$assay1_label_rle <- renderAssayLabelFn(1, getAssayNamesFn)
    output$assay2_label_rle <- renderAssayLabelFn(2, getAssayNamesFn)

    output$assay1_label_correlation <- renderAssayLabelFn(1, getAssayNamesFn)
    output$assay2_label_correlation <- renderAssayLabelFn(2, getAssayNamesFn)

    invisible(NULL)
}

# Keep the main run-normalization observer wrapper top-level so a later
# extraction wave can move it without reopening the module server wrapper.
runMetabNormNormalizationObserverWrapper <- function(
    workflowData,
    input,
    experimentPaths,
    omicType,
    normData,
    addLogFn = function(entry) invisible(entry),
    showNotificationFn = shiny::showNotification,
    reqFn = shiny::req,
    withProgressFn = shiny::withProgress,
    getPlotAestheticsFn = getPlotAesthetics,
    runObserverShellFn = runMetabNormNormalizationObserverShell,
    runPipelineShellFn = runMetabNormNormalizationPipelineShell,
    generateCompositeFromFilesFn = generateMetabNormCompositeFromFiles,
    savePlotFn = savePlot,
    logWarnFn = logger::log_warn
) {
    runObserverShellFn(
        workflowData = workflowData
        , addLogFn = addLogFn
        , showNotificationFn = showNotificationFn
        , reqFn = reqFn
        , withProgressFn = withProgressFn
        , runPipelineFn = function() {
            runPipelineShellFn(
                workflowData = workflowData
                , inputValues = list(
                    apply_itsd = input$apply_itsd
                    , itsd_aggregation = input$itsd_aggregation
                    , log_offset = input$log_offset
                    , norm_method = input$norm_method
                    , ruv_mode = input$ruv_mode
                    , auto_percentage_min = input$auto_percentage_min
                    , auto_percentage_max = input$auto_percentage_max
                    , ruv_grouping_variable = input$ruv_grouping_variable
                    , separation_metric = input$separation_metric
                    , k_penalty_weight = input$k_penalty_weight
                    , adaptive_k_penalty = input$adaptive_k_penalty
                    , ruv_k = input$ruv_k
                    , ruv_percentage = input$ruv_percentage
                )
                , experimentPaths = experimentPaths
                , omicType = omicType
                , normData = normData
                , getPlotAestheticsFn = function() {
                    getPlotAestheticsFn(
                        colorVariable = input$color_variable
                        , shapeVariable = input$shape_variable
                    )
                }
                , addLogFn = addLogFn
                , reqFn = reqFn
                , generateCompositeFromFilesFn = generateCompositeFromFilesFn
                , savePlotFn = savePlotFn
                , logWarnFn = logWarnFn
            )
        }
    )
}

