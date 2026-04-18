registerLipidNormDesignDrivenChoiceObserver <- function(
    session,
    workflowData,
    observeFn = shiny::observe,
    updateSelectInputFn = shiny::updateSelectInput
) {
    observeFn({
        if (is.null(workflowData$design_matrix)) {
            return()
        }

        design_cols <- colnames(workflowData$design_matrix)

        plot_available_vars <- intersect(
            design_cols,
            c("group", "factor1", "factor2", "batch", "technical_replicate_id", "sample_id")
        )

        if (length(plot_available_vars) > 0) {
            default_plot_var <- if ("group" %in% plot_available_vars) {
                "group"
            } else {
                plot_available_vars[1]
            }

            updateSelectInputFn(
                session, "color_variable"
                , choices = plot_available_vars
                , selected = default_plot_var
            )
            updateSelectInputFn(
                session, "shape_variable"
                , choices = plot_available_vars
                , selected = default_plot_var
            )
        }

        ruv_available_vars <- intersect(design_cols, c("group", "factor1", "factor2", "batch"))
        if (length(ruv_available_vars) > 0) {
            updateSelectInputFn(
                session, "ruv_grouping_variable"
                , choices = ruv_available_vars
                , selected = if ("group" %in% ruv_available_vars) "group" else ruv_available_vars[1]
            )
        }
    })

    invisible(session)
}

registerLipidNormAssayNameInitializationObserver <- function(
    workflowData,
    normData,
    observeFn = shiny::observe,
    reqFn = shiny::req,
    logInfoFn = logger::log_info,
    logWarnFn = logger::log_warn
) {
    observeFn({
        reqFn(workflowData$state_manager)

        tryCatch({
            current_s4 <- workflowData$state_manager$getState()

            if (!inherits(current_s4, "LipidomicsAssayData")) {
                return()
            }

            detected_assays <- names(current_s4@lipid_data)
            normData$assay_names <- detected_assays
            logInfoFn(paste("Detected assays:", paste(detected_assays, collapse = ", ")))

            for (assay_name in normData$assay_names) {
                if (!assay_name %in% names(normData$itsd_selections)) {
                    normData$itsd_selections[[assay_name]] <- NULL
                }
            }
        }, error = function(e) {
            logWarnFn(paste("Could not detect assay names:", e$message))
        })
    })

    invisible(normData)
}

registerLipidNormRuvCancorOutputs <- function(
    output,
    normData,
    observeFn = shiny::observe,
    reqFn = shiny::req,
    walkFn = purrr::walk,
    renderPlotFn = shiny::renderPlot,
    renderTextFn = shiny::renderText,
    renderDataTableFn = DT::renderDataTable,
    datatableFn = DT::datatable
) {
    observeFn({
        reqFn(normData$assay_names)
        reqFn(length(normData$ruv_optimization_results) > 0)

        walkFn(normData$assay_names, \(assay_name) {
            safe_name <- gsub("[^A-Za-z0-9]", "_", tolower(assay_name))

            output[[paste0("cancor_plot_", safe_name)]] <- renderPlotFn({
                result <- normData$ruv_optimization_results[[assay_name]]
                if (!is.null(result) && isTRUE(result$success) && !is.null(result$cancor_plot)) {
                    result$cancor_plot
                } else {
                    ggplot2::ggplot() +
                        ggplot2::annotate("text", x = 0.5, y = 0.5, label = "No data", size = 6) +
                        ggplot2::theme_void()
                }
            })

            output[[paste0("ruv_summary_", safe_name)]] <- renderTextFn({
                result <- normData$ruv_optimization_results[[assay_name]]
                if (!is.null(result) && isTRUE(result$success)) {
                    paste0(
                        "Best k: ", result$best_k, "\n"
                        , "Best %: ", result$best_percentage, "\n"
                        , "Separation: ", round(result$separation_score, 4), "\n"
                        , "Controls: ", sum(result$control_genes_index, na.rm = TRUE)
                    )
                } else if (!is.null(result)) {
                    paste0("Failed: ", result$error)
                } else {
                    "Not yet computed"
                }
            })

            output[[paste0("ruv_table_", safe_name)]] <- renderDataTableFn({
                result <- normData$ruv_optimization_results[[assay_name]]
                if (!is.null(result) && isTRUE(result$success) && !is.null(result$optimization_results)) {
                    datatableFn(
                        result$optimization_results
                        , options = list(pageLength = 5, dom = "t")
                        , rownames = FALSE
                    )
                } else {
                    NULL
                }
            })
        })
    })

    invisible(output)
}

registerLipidNormItsdTableOutputs <- function(
    output,
    workflowData,
    normData,
    observeFn = shiny::observe,
    reqFn = shiny::req,
    walkFn = purrr::walk,
    renderDataTableFn = DT::renderDataTable,
    datatableFn = DT::datatable,
    formatStyleFn = DT::formatStyle,
    styleEqualFn = DT::styleEqual,
    formatRoundFn = DT::formatRound,
    buildLipidItsdSelectionTableFn = buildLipidItsdSelectionTable
) {
    observeFn({
        reqFn(normData$assay_names)
        reqFn(workflowData$state_manager)

        current_s4 <- tryCatch({
            workflowData$state_manager$getState()
        }, error = function(e) NULL)

        if (is.null(current_s4) || !inherits(current_s4, "LipidomicsAssayData")) {
            return()
        }

        walkFn(normData$assay_names, \(assay_name) {
            safe_name <- gsub("[^A-Za-z0-9]", "_", tolower(assay_name))
            output_id <- paste0("itsd_table_", safe_name)

            output[[output_id]] <- renderDataTableFn({
                assay_data <- current_s4@lipid_data[[assay_name]]
                if (is.null(assay_data)) {
                    return(NULL)
                }

                selection_table <- buildLipidItsdSelectionTableFn(
                    assay_data = assay_data
                    , lipid_id_col = current_s4@lipid_id_column
                    , annotation_cols = current_s4@annotation_id_column
                )

                preselected <- which(selection_table$is_candidate)

                datatableFn(
                    selection_table
                    , selection = list(mode = "multiple", selected = preselected)
                    , filter = "top"
                    , options = list(
                        pageLength = 10
                        , scrollX = TRUE
                        , order = list(list(4, "desc"), list(3, "asc"))
                    )
                    , rownames = FALSE
                ) |>
                    formatStyleFn(
                        "is_candidate"
                        , backgroundColor = styleEqualFn(TRUE, "#d4edda")
                    ) |>
                    formatRoundFn(columns = c("mean_intensity", "cv_percent"), digits = 2)
            })
        })
    })

    invisible(output)
}

registerLipidNormItsdSelectionTracking <- function(
    input,
    normData,
    observeFn = shiny::observe,
    reqFn = shiny::req,
    walkFn = purrr::walk,
    observeEventFn = shiny::observeEvent,
    logInfoFn = logger::log_info
) {
    observeFn({
        reqFn(normData$assay_names)

        walkFn(normData$assay_names, \(assay_name) {
            safe_name <- gsub("[^A-Za-z0-9]", "_", tolower(assay_name))
            input_id <- paste0("itsd_table_", safe_name, "_rows_selected")

            observeEventFn(input[[input_id]], {
                normData$itsd_selections[[assay_name]] <- input[[input_id]]
                logInfoFn(paste(
                    "ITSD selection updated for", assay_name, ":"
                    , length(input[[input_id]]), "features selected"
                ))
            }, ignoreNULL = FALSE)
        })
    })

    invisible(normData)
}

handleLipidNormPreNormalizationQc <- function(
    workflowData,
    experimentPaths,
    normData,
    addLog,
    getPlotAestheticsFn = getPlotAesthetics,
    reqFn = shiny::req,
    generateLipidQcPlotsFn = generateLipidQcPlots,
    logInfoFn = logger::log_info,
    logWarnFn = logger::log_warn,
    logErrorFn = logger::log_error
) {
    logInfoFn("=== GENERATING PRE-NORMALIZATION QC PLOTS ===")

    reqFn(workflowData$state_manager)
    current_s4 <- workflowData$state_manager$getState()

    if (is.null(current_s4)) {
        logWarnFn("No S4 object available for QC plot generation")
        return(invisible(FALSE))
    }

    if (inherits(current_s4, "LipidomicsAssayData")) {
        detected_assays <- names(current_s4@lipid_data)
        if (length(detected_assays) > 0 && is.null(normData$assay_names)) {
            normData$assay_names <- detected_assays
            logInfoFn(paste("Set assay names:", paste(detected_assays, collapse = ", ")))
        }
    }

    aesthetics <- getPlotAestheticsFn()

    tryCatch({
        generateLipidQcPlotsFn(
            theObject = current_s4
            , experiment_paths = experimentPaths
            , stage = "post_filter"
            , grouping_variable = aesthetics$color_var
            , shape_variable = aesthetics$shape_var
        )

        normData$plot_refresh_trigger <- normData$plot_refresh_trigger + 1
        normData$pre_norm_qc_generated <- TRUE
        logInfoFn("Pre-normalization QC plots generated successfully")

        invisible(TRUE)
    }, error = function(e) {
        logErrorFn(paste("Error generating pre-normalization QC:", e$message))
        addLog(paste("Error generating Pre-QC:", e$message))

        invisible(FALSE)
    })
}

handleLipidNormSelectedTabPreNormalizationTrigger <- function(
    selectedTabValue,
    workflowData,
    experimentPaths,
    normData,
    addLog,
    getPlotAestheticsFn = getPlotAesthetics,
    reqFn = shiny::req,
    withProgressFn = shiny::withProgress,
    handlePreNormalizationQcFn = handleLipidNormPreNormalizationQc,
    logInfoFn = logger::log_info
) {
    if (is.null(selectedTabValue) || !identical(selectedTabValue, "norm")) {
        return(invisible(FALSE))
    }

    logInfoFn("Normalization tab selected - checking if pre-QC needed")

    if (normData$pre_norm_qc_generated) {
        return(invisible(FALSE))
    }

    reqFn(workflowData$state_manager)
    current_s4 <- tryCatch(
        workflowData$state_manager$getState(),
        error = function(e) NULL
    )

    if (is.null(current_s4) || !inherits(current_s4, "LipidomicsAssayData")) {
        return(invisible(FALSE))
    }

    logInfoFn("Auto-triggering pre-normalization QC plots")

    withProgressFn(
        message = "Generating Pre-Normalization QC..."
        , value = 0.5
        , {
            handlePreNormalizationQcFn(
                workflowData = workflowData
                , experimentPaths = experimentPaths
                , normData = normData
                , addLog = addLog
                , getPlotAestheticsFn = getPlotAestheticsFn
            )
        }
    )

    normData$pre_norm_qc_generated <- TRUE

    invisible(TRUE)
}

