handleLipidNormRunNormalization <- function(
    input,
    workflowData,
    experimentPaths,
    omicType,
    normData,
    addLog,
    getPlotAestheticsFn = getPlotAesthetics,
    reqFn = shiny::req,
    withProgressFn = shiny::withProgress,
    incProgressFn = shiny::incProgress,
    generateLipidQcPlotsFn = generateLipidQcPlots,
    buildLipidItsdSelectionTableFn = buildLipidItsdSelectionTable,
    compactFn = purrr::compact,
    imapFn = purrr::imap,
    normaliseUntransformedDataFn = normaliseUntransformedData,
    logTransformAssaysFn = logTransformAssays,
    normaliseBetweenSamplesFn = normaliseBetweenSamples,
    runLipidPerAssayRuvOptimizationFn = runLipidPerAssayRuvOptimization,
    extractLipidBestKPerAssayFn = extractLipidBestKPerAssay,
    extractLipidCtrlPerAssayFn = extractLipidCtrlPerAssay,
    ruvIII_C_VaryingFn = ruvIII_C_Varying,
    generateCompositeFromFilesFn = buildLipidNormCompositeFromFilesGenerator(),
    savePlotFn = savePlot,
    dirExistsFn = dir.exists,
    filePathFn = file.path,
    showNotificationFn = shiny::showNotification,
    logWarnFn = logger::log_warn,
    logErrorFn = logger::log_error
) {
    reqFn(workflowData$state_manager)

    addLog("Starting normalization pipeline...")

    withProgressFn(
        message = "Running normalization pipeline..."
        , value = 0
        , {
            tryCatch({
                current_s4 <- workflowData$state_manager$getState()
                reqFn(current_s4)

                aesthetics <- getPlotAestheticsFn()
                total_steps <- 6

                incProgressFn(1 / total_steps, detail = "Capturing pre-normalization state...")
                normData$post_filter_obj <- current_s4
                addLog("Post-filtering state captured")

                generateLipidQcPlotsFn(
                    theObject = current_s4
                    , experiment_paths = experimentPaths
                    , stage = "post_filter"
                    , grouping_variable = aesthetics$color_var
                    , shape_variable = aesthetics$shape_var
                )
                addLog("Pre-normalization QC plots generated")

                incProgressFn(1 / total_steps, detail = "Applying ITSD normalization...")
                if (isTRUE(input$apply_itsd)) {
                    addLog(paste("Applying ITSD normalization (aggregation:", input$itsd_aggregation, ")"))

                    itsd_feature_ids <- NULL
                    has_manual_selections <- any(sapply(normData$itsd_selections, \(x) length(x) > 0))

                    if (has_manual_selections) {
                        lipid_id_col <- current_s4@lipid_id_column
                        annotation_col <- current_s4@annotation_id_column
                        itsd_feature_ids <- imapFn(normData$itsd_selections, \(row_indices, assay_name) {
                            if (is.null(row_indices) || length(row_indices) == 0) {
                                return(NULL)
                            }

                            assay_data <- current_s4@lipid_data[[assay_name]]
                            if (is.null(assay_data)) {
                                return(NULL)
                            }

                            selection_table <- buildLipidItsdSelectionTableFn(
                                assay_data = assay_data
                                , lipid_id_col = lipid_id_col
                                , annotation_cols = annotation_col
                            )

                            selected_ids <- selection_table$feature_id[row_indices]
                            addLog(paste("Assay", assay_name, ":", length(selected_ids), "ITSD features selected"))
                            selected_ids
                        })
                        itsd_feature_ids <- compactFn(itsd_feature_ids)
                        if (length(itsd_feature_ids) == 0) {
                            itsd_feature_ids <- NULL
                        }
                    }

                    current_s4 <- normaliseUntransformedDataFn(
                        theObject = current_s4
                        , method = "ITSD"
                        , itsd_aggregation = input$itsd_aggregation
                        , itsd_feature_ids = itsd_feature_ids
                    )
                    normData$post_itsd_obj <- current_s4

                    workflowData$state_manager$saveState(
                        state_name = "lipid_itsd_norm"
                        , s4_data_object = current_s4
                        , config_object = workflowData$config_list
                        , description = paste("ITSD normalization (aggregation:", input$itsd_aggregation, ")")
                    )
                    addLog("ITSD normalization complete")
                } else {
                    addLog("ITSD normalization skipped")
                }

                incProgressFn(1 / total_steps, detail = "Applying log2 transformation...")
                addLog(paste("Applying log2 transformation (offset:", input$log_offset, ")"))

                current_s4 <- logTransformAssaysFn(
                    theObject = current_s4
                    , offset = input$log_offset
                )
                normData$post_log2_obj <- current_s4

                workflowData$state_manager$saveState(
                    state_name = "lipid_log2"
                    , s4_data_object = current_s4
                    , config_object = workflowData$config_list
                    , description = paste("Log2 transformation (offset:", input$log_offset, ")")
                )
                addLog("Log2 transformation complete")

                incProgressFn(1 / total_steps, detail = "Applying between-sample normalization...")
                if (input$norm_method != "none") {
                    addLog(paste("Applying between-sample normalization (method:", input$norm_method, ")"))

                    current_s4 <- normaliseBetweenSamplesFn(
                        theObject = current_s4
                        , normalisation_method = input$norm_method
                    )
                }
                normData$post_norm_obj <- current_s4

                workflowData$state_manager$saveState(
                    state_name = "lipid_normalized"
                    , s4_data_object = current_s4
                    , config_object = workflowData$config_list
                    , description = paste("Between-sample normalization (method:", input$norm_method, ")")
                )
                addLog("Between-sample normalization complete")
                normData$normalization_complete <- TRUE

                generateLipidQcPlotsFn(
                    theObject = current_s4
                    , experiment_paths = experimentPaths
                    , stage = "post_norm"
                    , grouping_variable = aesthetics$color_var
                    , shape_variable = aesthetics$shape_var
                )
                addLog("Post-normalization QC plots generated")

                incProgressFn(1 / total_steps, detail = "Running RUV-III batch correction...")
                if (input$ruv_mode != "skip") {
                    addLog(paste("Running RUV-III (mode:", input$ruv_mode, ")"))

                    ruv_params <- list(
                        percentage_min = input$auto_percentage_min
                        , percentage_max = input$auto_percentage_max
                        , ruv_grouping_variable = input$ruv_grouping_variable
                        , separation_metric = input$separation_metric
                        , k_penalty_weight = input$k_penalty_weight
                        , adaptive_k_penalty = input$adaptive_k_penalty
                        , manual_k = input$ruv_k
                        , manual_percentage = input$ruv_percentage
                    )

                    ruv_results <- runLipidPerAssayRuvOptimizationFn(
                        theObject = current_s4
                        , ruv_mode = input$ruv_mode
                        , params = ruv_params
                        , experiment_paths = experimentPaths
                    )
                    normData$ruv_optimization_results <- ruv_results

                    best_k_list <- extractLipidBestKPerAssayFn(ruv_results)
                    ctrl_list <- extractLipidCtrlPerAssayFn(ruv_results)

                    addLog(paste(
                        "RUV optimization complete. Best k per assay:",
                        paste(names(best_k_list), "=", unlist(best_k_list), collapse = ", ")
                    ))

                    current_s4 <- ruvIII_C_VaryingFn(
                        theObject = current_s4
                        , ruv_grouping_variable = input$ruv_grouping_variable
                        , ruv_number_k = best_k_list
                        , ctrl = ctrl_list
                    )
                    normData$ruv_corrected_obj <- current_s4
                    normData$ruv_complete <- TRUE

                    workflowData$state_manager$saveState(
                        state_name = "lipid_ruv_corrected"
                        , s4_data_object = current_s4
                        , config_object = workflowData$config_list
                        , description = "RUV-III batch correction complete"
                    )
                    addLog("RUV-III correction applied")

                    incProgressFn(1 / total_steps, detail = "Generating RUV QC plots...")
                    generateLipidQcPlotsFn(
                        theObject = current_s4
                        , experiment_paths = experimentPaths
                        , stage = "ruv_corrected"
                        , grouping_variable = aesthetics$color_var
                        , shape_variable = aesthetics$shape_var
                    )
                    addLog("RUV QC plots generated")
                } else {
                    addLog("RUV-III skipped")
                    normData$ruv_corrected_obj <- current_s4
                    normData$ruv_complete <- TRUE
                }

                addLog("Generating composite QC figure...")
                tryCatch({
                    qc_dir <- experimentPaths$lipid_qc_dir
                    if (!is.null(qc_dir) && dirExistsFn(qc_dir) && !is.null(normData$assay_names)) {
                        ncol_composite <- if (input$ruv_mode == "skip") 2 else 3
                        column_labels <- if (input$ruv_mode == "skip") {
                            c("Pre-Normalisation", "Post-Normalisation")
                        } else {
                            c("Pre-Normalisation", "Post-Normalisation", "RUV-Corrected")
                        }

                        all_plot_files <- c()
                        all_row_labels <- list()
                        label_counter <- 1

                        for (assay_name in normData$assay_names) {
                            safe_name <- gsub("[^A-Za-z0-9]", "_", tolower(assay_name))
                            plot_types <- c("pca", "density", "rle", "correlation")

                            for (plot_type in plot_types) {
                                if (input$ruv_mode == "skip") {
                                    files <- c(
                                        filePathFn(qc_dir, sprintf("%s_pre_norm_%s.png", safe_name, plot_type))
                                        , filePathFn(qc_dir, sprintf("%s_post_norm_%s.png", safe_name, plot_type))
                                    )
                                    labels <- c(
                                        sprintf("%s)", letters[label_counter])
                                        , sprintf("%s)", letters[label_counter + 1])
                                    )
                                    label_counter <- label_counter + 2
                                } else {
                                    files <- c(
                                        filePathFn(qc_dir, sprintf("%s_pre_norm_%s.png", safe_name, plot_type))
                                        , filePathFn(qc_dir, sprintf("%s_post_norm_%s.png", safe_name, plot_type))
                                        , filePathFn(qc_dir, sprintf("%s_ruv_corrected_%s.png", safe_name, plot_type))
                                    )
                                    labels <- c(
                                        sprintf("%s)", letters[label_counter])
                                        , sprintf("%s)", letters[label_counter + 1])
                                        , sprintf("%s)", letters[label_counter + 2])
                                    )
                                    label_counter <- label_counter + 3
                                }

                                all_plot_files <- c(all_plot_files, files)
                                row_key <- sprintf("%s_%s", safe_name, plot_type)
                                all_row_labels[[row_key]] <- labels
                            }
                        }

                        composite_res <- generateCompositeFromFilesFn(
                            plot_files = all_plot_files
                            , ncol = ncol_composite
                            , row_labels = all_row_labels
                            , column_labels = column_labels
                        )

                        if (!is.null(composite_res)) {
                            savePlotFn(
                                composite_res$plot
                                , qc_dir
                                , paste0(omicType, "_composite_QC_figure")
                                , width = composite_res$width
                                , height = composite_res$height
                                , dpi = 150
                                , limitsize = FALSE
                            )
                            addLog(sprintf(
                                "Composite QC figure saved to: %s",
                                filePathFn(qc_dir, "composite_QC_figure")
                            ))
                        }
                    }
                }, error = function(e) {
                    addLog(paste("Warning: Could not generate composite QC figure:", e$message))
                    logWarnFn(paste("Composite QC generation failed:", e$message))
                })

                normData$plot_refresh_trigger <- normData$plot_refresh_trigger + 1

                addLog("Normalization pipeline complete!")
                showNotificationFn("Normalization pipeline complete!", type = "message")

                invisible(TRUE)
            }, error = function(e) {
                addLog(paste("ERROR:", e$message))
                logErrorFn(paste("Normalization pipeline error:", e$message))
                showNotificationFn(paste("Error:", e$message), type = "error")

                invisible(FALSE)
            })
        }
    )
}

handleLipidNormExportSession <- function(
    input,
    workflowData,
    experimentPaths,
    experimentLabel,
    normData,
    addLog,
    reqFn = shiny::req,
    showNotificationFn = shiny::showNotification,
    withProgressFn = shiny::withProgress,
    incProgressFn = shiny::incProgress,
    dirExistsFn = dir.exists,
    dirCreateFn = dir.create,
    saveRdsFn = saveRDS,
    writeLinesFn = writeLines,
    getTimeFn = Sys.time,
    logInfoFn = logger::log_info,
    logWarnFn = logger::log_warn,
    logErrorFn = logger::log_error
) {
    logInfoFn("=== EXPORT NORMALIZED SESSION BUTTON CLICKED ===")
    reqFn(workflowData$state_manager)

    if (!normData$normalization_complete) {
        showNotificationFn(
            "Please complete normalization before exporting session data."
            , type = "warning"
            , duration = 5
        )
        return(invisible(FALSE))
    }

    tryCatch({
        source_dir <- experimentPaths$source_dir
        if (is.null(source_dir) || !dirExistsFn(source_dir)) {
            source_dir <- experimentPaths$export_dir
            if (is.null(source_dir)) {
                stop("Could not find a valid directory to save session data.")
            }
            if (!dirExistsFn(source_dir)) {
                dirCreateFn(source_dir, recursive = TRUE)
            }
        }

        export_artifacts <- list(
            session_filepath = NULL
            , session_filename = NULL
        )

        withProgressFn(message = "Exporting normalized session data...", value = 0, {

            # Step 1: Gather all essential data for DE analysis
            incProgressFn(0.2, detail = "Gathering data...")

            # Get current state from R6 manager
            current_state_name <- workflowData$state_manager$current_state
            current_s4 <- workflowData$state_manager$getState(current_state_name)

            # Calculate feature counts per assay
            feature_counts_per_assay <- NULL
            if (!is.null(current_s4) && inherits(current_s4, "LipidomicsAssayData")) {
                feature_counts_per_assay <- purrr::map(names(current_s4@lipid_data), \(assay_name) {
                    assay_data <- current_s4@lipid_data[[assay_name]]
                    if (!is.null(assay_data)) {
                        n_features <- length(unique(assay_data[[current_s4@lipid_id_column]]))
                        n_samples <- ncol(assay_data) - length(c(current_s4@lipid_id_column, current_s4@annotation_id_column))
                        list(features = n_features, samples = n_samples)
                    } else {
                        list(features = 0, samples = 0)
                    }
                }) |> purrr::set_names(names(current_s4@lipid_data))
            }

            export_timestamp <- getTimeFn()

            # Build comprehensive session data
            session_data <- list(
                # --- R6 State Info ---
                r6_current_state_name = current_state_name
                , current_s4_object = current_s4

                # --- Workflow artifacts ---
                , contrasts_tbl = workflowData$contrasts_tbl
                , design_matrix = workflowData$design_matrix
                , config_list = workflowData$config_list

                # --- Lipidomics-specific ---
                , itsd_selections = normData$itsd_selections
                , ruv_optimization_results = normData$ruv_optimization_results
                , correlation_results = normData$correlation_results
                , assay_names = normData$assay_names

                # --- Export metadata ---
                , export_timestamp = export_timestamp
                , omic_type = "lipidomics"
                , experiment_label = experimentLabel

                # --- Normalization parameters ---
                , normalization_method = input$norm_method
                , ruv_mode = input$ruv_mode
                , itsd_applied = input$apply_itsd
                , itsd_aggregation = if (isTRUE(input$apply_itsd)) input$itsd_aggregation else NA
                , log_offset = input$log_offset
                , correlation_threshold = input$min_pearson_correlation_threshold
                , ruv_grouping_variable = input$ruv_grouping_variable

                # --- Feature counts per assay ---
                , feature_counts = feature_counts_per_assay
                , lipid_counts = workflowData$lipid_counts

                # --- QC parameters ---
                , qc_params = workflowData$qc_params

                # --- Processing flags ---
                , normalization_complete = normData$normalization_complete
                , ruv_complete = normData$ruv_complete
                , correlation_filtering_complete = normData$correlation_filtering_complete
            )

            logInfoFn("*** EXPORT: Gathered session data successfully ***")
            logInfoFn(sprintf("*** EXPORT: Assays: %s ***", paste(normData$assay_names, collapse = ", ")))
            logInfoFn(sprintf("*** EXPORT: Contrasts available: %d ***",
                ifelse(is.null(session_data$contrasts_tbl), 0, nrow(session_data$contrasts_tbl))))

            # Step 2: Save to RDS file with timestamp
            incProgressFn(0.3, detail = "Saving to file...")

            timestamp_str <- format(export_timestamp, "%Y%m%d_%H%M%S")
            export_artifacts$session_filename <- sprintf("lipid_filtered_session_data_%s.rds", timestamp_str)
            export_artifacts$session_filepath <- file.path(source_dir, export_artifacts$session_filename)

            saveRdsFn(session_data, export_artifacts$session_filepath)
            logInfoFn(sprintf("*** EXPORT: Session data saved to: %s ***", export_artifacts$session_filepath))

            # Step 3: Save "latest" version for easy access
            incProgressFn(0.1, detail = "Creating latest version...")

            latest_filename <- "lipid_filtered_session_data_latest.rds"
            latest_filepath <- file.path(source_dir, latest_filename)

            saveRdsFn(session_data, latest_filepath)
            logInfoFn(sprintf("*** EXPORT: Latest version saved to: %s ***", latest_filepath))

            # Step 4: Save individual metadata files for redundancy
            incProgressFn(0.1, detail = "Saving metadata files...")

            tryCatch({
                # Save RUV optimization results (per-assay)
                if (!is.null(session_data$ruv_optimization_results) && length(session_data$ruv_optimization_results) > 0) {
                    ruv_file <- file.path(source_dir, "lipid_ruv_optimization_results.RDS")
                    saveRdsFn(session_data$ruv_optimization_results, ruv_file)
                    logInfoFn("*** EXPORT: Saved lipid_ruv_optimization_results.RDS ***")
                }

                # Save ITSD selections
                if (!is.null(session_data$itsd_selections) && length(session_data$itsd_selections) > 0) {
                    itsd_file <- file.path(source_dir, "lipid_itsd_selections.RDS")
                    saveRdsFn(session_data$itsd_selections, itsd_file)
                    logInfoFn("*** EXPORT: Saved lipid_itsd_selections.RDS ***")
                }

                # Save QC parameters
                if (!is.null(session_data$qc_params)) {
                    qc_params_file <- file.path(source_dir, "lipid_qc_params.RDS")
                    saveRdsFn(session_data$qc_params, qc_params_file)
                    logInfoFn("*** EXPORT: Saved lipid_qc_params.RDS ***")
                }

            }, error = function(e) {
                logWarnFn(sprintf("*** WARNING: Some metadata files could not be saved: %s ***", e$message))
            })

            # Step 5: Create human-readable summary file
            incProgressFn(0.2, detail = "Creating summary...")

            # Build RUV summary per assay
            ruv_summary_lines <- ""
            if (!is.null(session_data$ruv_optimization_results)) {
                for (assay_name in names(session_data$ruv_optimization_results)) {
                    result <- session_data$ruv_optimization_results[[assay_name]]
                    if (!is.null(result) && isTRUE(result$success)) {
                        ruv_summary_lines <- paste0(ruv_summary_lines, sprintf(
                            "\n  %s: k=%d, %%=%.1f, controls=%d"
                            , assay_name
                            , result$best_k
                            , result$best_percentage
                            , sum(result$control_genes_index, na.rm = TRUE)
                        ))
                    }
                }
            }
            if (ruv_summary_lines == "") ruv_summary_lines <- "\n  (RUV skipped or not applied)"

            # Build feature counts summary
            feature_summary_lines <- ""
            if (!is.null(feature_counts_per_assay)) {
                for (assay_name in names(feature_counts_per_assay)) {
                    counts <- feature_counts_per_assay[[assay_name]]
                    feature_summary_lines <- paste0(feature_summary_lines, sprintf(
                        "\n  %s: %d features, %d samples"
                        , assay_name
                        , counts$features
                        , counts$samples
                    ))
                }
            }

            summary_content <- sprintf(
                "Lipidomics Normalized Session Data Export Summary\n===================================================\n\nExport Timestamp: %s\nSession File: %s\n\nData Summary:%s\n\nNormalization Parameters:\n- Method: %s\n- ITSD applied: %s\n- ITSD aggregation: %s\n- Log2 offset: %s\n- RUV mode: %s\n- RUV grouping variable: %s\n- Correlation threshold: %s\n\nRUV Optimization Results (per-assay):%s\n\nContrasts:\n%s\n\nThis data is ready for differential expression analysis.\nUse 'Load Filtered Session' in the DE tab to import.\n"
                , format(export_timestamp, "%Y-%m-%d %H:%M:%S")
                , export_artifacts$session_filename
                , feature_summary_lines
                , session_data$normalization_method
                , ifelse(isTRUE(session_data$itsd_applied), "Yes", "No")
                , ifelse(is.na(session_data$itsd_aggregation), "N/A", session_data$itsd_aggregation)
                , session_data$log_offset
                , session_data$ruv_mode
                , session_data$ruv_grouping_variable
                , ifelse(is.null(session_data$correlation_threshold), "N/A", session_data$correlation_threshold)
                , ruv_summary_lines
                , if (!is.null(session_data$contrasts_tbl)) paste(session_data$contrasts_tbl$friendly_names, collapse = "\n") else "None defined"
            )

            summary_filepath <- file.path(source_dir, "lipid_filtered_session_summary.txt")
            writeLinesFn(summary_content, summary_filepath)
            logInfoFn(sprintf("*** EXPORT: Summary saved to: %s ***", summary_filepath))
        })

        addLog(paste("Exported comprehensive session data to:", export_artifacts$session_filepath))
        showNotificationFn(
            sprintf("Session data exported successfully!\nSaved as: %s\nSee summary file for details.", export_artifacts$session_filename)
            , type = "message"
            , duration = 10
        )

        logInfoFn("=== EXPORT NORMALIZED SESSION COMPLETED SUCCESSFULLY ===")
        invisible(export_artifacts)

    }, error = function(e) {
        logErrorFn(paste("*** ERROR in session export:", e$message, "***"))
        addLog(paste("Export error:", e$message))
        showNotificationFn(paste("Export error:", e$message), type = "error", duration = 10)
        invisible(FALSE)
    })
}

handleLipidNormSkipCorrelationFilter <- function(
    workflowData,
    normData,
    addLog,
    reqFn = shiny::req,
    showNotificationFn = shiny::showNotification
) {
    reqFn(workflowData$state_manager)
    reqFn(normData$ruv_complete || normData$normalization_complete)

    current_s4 <- if (!is.null(normData$ruv_corrected_obj)) {
        normData$ruv_corrected_obj
    } else {
        normData$post_norm_obj
    }

    if (is.null(current_s4)) {
        return(invisible(FALSE))
    }

    workflowData$state_manager$saveState(
        state_name = "lipid_norm_complete"
        , s4_data_object = current_s4
        , config_object = workflowData$config_list
        , description = "Normalization complete (correlation filtering skipped)"
    )

    updated_status <- workflowData$tab_status
    updated_status$quality_control <- "complete"
    updated_status$normalization <- "complete"
    workflowData$tab_status <- updated_status

    addLog("Correlation filtering skipped - ready for DE analysis")
    showNotificationFn("Normalization complete! Proceeding to DE analysis.", type = "message")

    invisible(TRUE)
}

handleLipidNormResetNormalization <- function(
    workflowData,
    normData,
    addLog,
    reqFn = shiny::req,
    showNotificationFn = shiny::showNotification
) {
    reqFn(workflowData$state_manager)

    tryCatch({
        if (!is.null(normData$post_filter_obj)) {
            workflowData$state_manager$saveState(
                state_name = "lipid_reset"
                , s4_data_object = normData$post_filter_obj
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

        addLog("Reset to pre-normalization state")
        showNotificationFn("Reset to pre-normalization state", type = "message")

        invisible(TRUE)
    }, error = function(e) {
        addLog(paste("ERROR during reset:", e$message))
        showNotificationFn(paste("Error:", e$message), type = "error")

        invisible(FALSE)
    })
}

handleLipidNormApplyCorrelationFilter <- function(
    input,
    workflowData,
    normData,
    addLog,
    reqFn = shiny::req,
    showNotificationFn = shiny::showNotification,
    removeNotificationFn = shiny::removeNotification,
    pearsonCorForSamplePairsFn = pearsonCorForSamplePairs,
    filterSamplesByLipidCorrelationThresholdFn = filterSamplesByLipidCorrelationThreshold,
    logInfoFn = logger::log_info,
    logErrorFn = logger::log_error
) {
    reqFn(workflowData$state_manager)
    reqFn(normData$ruv_complete || normData$normalization_complete)

    threshold <- input$min_pearson_correlation_threshold
    addLog(paste("Applying correlation filter (threshold:", threshold, ")"))
    showNotificationFn("Applying correlation filter...", id = "corr_working", duration = NULL)

    tryCatch({
        current_s4 <- if (!is.null(normData$ruv_corrected_obj)) {
            normData$ruv_corrected_obj
        } else {
            normData$post_norm_obj
        }
        reqFn(current_s4)

        logInfoFn("Calculating Pearson correlations per sample pair...")
        corr_results <- pearsonCorForSamplePairsFn(
            theObject = current_s4
            , correlation_group = input$ruv_grouping_variable
        )

        normData$correlation_results <- corr_results

        filtered_s4 <- filterSamplesByLipidCorrelationThresholdFn(
            theObject = current_s4
            , pearson_correlation_per_pair = corr_results
            , min_pearson_correlation_threshold = threshold
        )

        normData$correlation_filtered_obj <- filtered_s4
        normData$correlation_filtering_complete <- TRUE

        workflowData$state_manager$saveState(
            state_name = "lipid_correlation_filtered"
            , s4_data_object = filtered_s4
            , config_object = workflowData$config_list
            , description = paste("Correlation filtering (threshold:", threshold, ")")
        )

        updated_status <- workflowData$tab_status
        updated_status$quality_control <- "complete"
        updated_status$normalization <- "complete"
        workflowData$tab_status <- updated_status

        addLog("Correlation filtering complete")
        removeNotificationFn("corr_working")
        showNotificationFn("Correlation filtering complete! Ready for DE analysis.", type = "message")

        invisible(TRUE)
    }, error = function(e) {
        addLog(paste("ERROR in correlation filtering:", e$message))
        logErrorFn(paste("Correlation filtering error:", e$message))
        removeNotificationFn("corr_working")
        showNotificationFn(paste("Error:", e$message), type = "error")

        invisible(FALSE)
    })
}

registerLipidNormRunNormalizationObserver <- function(
    input,
    workflowData,
    experimentPaths,
    omicType,
    normData,
    addLog,
    getPlotAestheticsFn,
    generateCompositeFromFilesFn,
    observeEventFn = shiny::observeEvent,
    handleRunNormalizationFn = handleLipidNormRunNormalization
) {
    observeEventFn(input$run_normalization, {
        handleRunNormalizationFn(
            input = input
            , workflowData = workflowData
            , experimentPaths = experimentPaths
            , omicType = omicType
            , normData = normData
            , addLog = addLog
            , getPlotAestheticsFn = getPlotAestheticsFn
            , generateCompositeFromFilesFn = generateCompositeFromFilesFn
        )
    })

    invisible(input)
}

registerLipidNormResetNormalizationObserver <- function(
    input,
    workflowData,
    normData,
    addLog,
    observeEventFn = shiny::observeEvent,
    handleResetNormalizationFn = handleLipidNormResetNormalization
) {
    observeEventFn(input$reset_normalization, {
        handleResetNormalizationFn(
            workflowData = workflowData
            , normData = normData
            , addLog = addLog
        )
    })

    invisible(input)
}

registerLipidNormApplyCorrelationFilterObserver <- function(
    input,
    workflowData,
    normData,
    addLog,
    observeEventFn = shiny::observeEvent,
    handleApplyCorrelationFilterFn = handleLipidNormApplyCorrelationFilter
) {
    observeEventFn(input$apply_correlation_filter, {
        handleApplyCorrelationFilterFn(
            input = input
            , workflowData = workflowData
            , normData = normData
            , addLog = addLog
        )
    })

    invisible(input)
}

registerLipidNormSkipCorrelationFilterObserver <- function(
    input,
    workflowData,
    normData,
    addLog,
    observeEventFn = shiny::observeEvent,
    handleSkipCorrelationFilterFn = handleLipidNormSkipCorrelationFilter
) {
    observeEventFn(input$skip_correlation_filter, {
        handleSkipCorrelationFilterFn(
            workflowData = workflowData
            , normData = normData
            , addLog = addLog
        )
    })

    invisible(input)
}

registerLipidNormExportSessionObserver <- function(
    input,
    workflowData,
    experimentPaths,
    experimentLabel,
    normData,
    addLog,
    observeEventFn = shiny::observeEvent,
    handleExportSessionFn = handleLipidNormExportSession
) {
    observeEventFn(input$export_session, {
        handleExportSessionFn(
            input = input
            , workflowData = workflowData
            , experimentPaths = experimentPaths
            , experimentLabel = experimentLabel
            , normData = normData
            , addLog = addLog
        )
    })

    invisible(input)
}

registerLipidNormSelectedTabPreNormalizationObserver <- function(
    selectedTab,
    workflowData,
    experimentPaths,
    normData,
    addLog,
    getPlotAestheticsFn,
    observeEventFn = shiny::observeEvent,
    handleSelectedTabPreNormalizationTriggerFn = handleLipidNormSelectedTabPreNormalizationTrigger
) {
    observeEventFn(selectedTab(), {
        handleSelectedTabPreNormalizationTriggerFn(
            selectedTabValue = selectedTab()
            , workflowData = workflowData
            , experimentPaths = experimentPaths
            , normData = normData
            , addLog = addLog
            , getPlotAestheticsFn = getPlotAestheticsFn
        )
    }, ignoreInit = FALSE)

    invisible(selectedTab)
}

