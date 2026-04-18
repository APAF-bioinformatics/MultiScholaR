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

