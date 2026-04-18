registerLipidIntensityAssayResultsOutput <- function(
    output,
    filterStatsVal,
    ns,
    renderUiFn = shiny::renderUI
) {
    output$assay_results_tabs <- renderUiFn({
        stats <- filterStatsVal()
        if (is.null(stats) || nrow(stats) == 0) {
            return(NULL)
        }

        tab_list <- lapply(seq_len(nrow(stats)), function(i) {
            shiny::tabPanel(
                stats$Assay[i]
                , shiny::br()
                , shiny::tags$table(
                    class = "table table-striped table-condensed"
                    , style = "width: auto;"
                    , shiny::tags$tbody(
                        shiny::tags$tr(
                            shiny::tags$td(shiny::strong("Original lipids:"))
                            , shiny::tags$td(stats$Original[i])
                        )
                        , shiny::tags$tr(
                            shiny::tags$td(shiny::strong("After filtering:"))
                            , shiny::tags$td(stats$Filtered[i])
                        )
                        , shiny::tags$tr(
                            shiny::tags$td(shiny::strong("Removed:"))
                            , shiny::tags$td(stats$Removed[i])
                        )
                        , shiny::tags$tr(
                            shiny::tags$td(shiny::strong("Percent retained:"))
                            , shiny::tags$td(paste0(stats$Percent_Retained[i], "%"))
                        )
                    )
                )
            )
        })

        do.call(shiny::tabsetPanel, c(list(id = ns("assay_stats_tabs")), tab_list))
    })

    output
}

registerLipidIntensityFilterPlotOutput <- function(
    output,
    filterPlotVal,
    renderPlotFn = shiny::renderPlot,
    reqFn = shiny::req,
    gridDrawFn = grid::grid.draw,
    printFn = print
) {
    output$filter_plot <- renderPlotFn({
        plot_obj <- filterPlotVal()
        reqFn(plot_obj)
        if (inherits(plot_obj, "grob") || inherits(plot_obj, "gtable")) {
            gridDrawFn(plot_obj)
        } else if (inherits(plot_obj, "ggplot")) {
            printFn(plot_obj)
        } else {
        }
    })

    output
}

registerLipidIntensityRevertObserver <- function(
    input,
    workflowData,
    output,
    filterStatsVal,
    filterPlotVal,
    observeEventFn = shiny::observeEvent,
    renderTextFn = shiny::renderText,
    logInfoFn = logger::log_info,
    logErrorFn = logger::log_error,
    showNotificationFn = shiny::showNotification
) {
    observeEventFn(input$revert_filter, {
        tryCatch({
            history <- workflowData$state_manager$getHistory()
            if (length(history) > 1) {
                prev_state_name <- history[[length(history) - 1]]
                workflowData$state_manager$revertToState(prev_state_name)
                output$filter_results <- renderTextFn(
                    paste("Reverted to previous state:", prev_state_name)
                )
                filterStatsVal(NULL)
                filterPlotVal(NULL)
                logInfoFn(paste("Reverted lipid intensity filter to", prev_state_name))
                showNotificationFn("Reverted successfully", type = "message")
            } else {
                stop("No previous state to revert to.")
            }
        }, error = function(e) {
            msg <- paste("Error reverting:", e$message)
            logErrorFn(msg)
            showNotificationFn(msg, type = "error")
        })
    })

    invisible(input)
}

registerLipidIntensityApplyFilterObserver <- function(
    input,
    workflowData,
    output,
    filterStatsVal,
    filterPlotVal,
    omicType,
    observeEventFn = shiny::observeEvent,
    reqFn = shiny::req,
    showNotificationFn = shiny::showNotification,
    removeNotificationFn = shiny::removeNotification,
    renderTextFn = shiny::renderText,
    logInfoFn = logger::log_info,
    logWarnFn = logger::log_warn,
    logErrorFn = logger::log_error,
    lipidIntensityFilteringFn = lipidIntensityFiltering,
    updateLipidFilteringFn = updateLipidFiltering
) {
    observeEventFn(input$apply_filter, {
        reqFn(workflowData$state_manager)

        showNotificationFn(
            "Applying lipid intensity filter..."
            , id = "lipid_intensity_filter_working"
            , duration = NULL
        )

        tryCatch({
            current_s4 <- workflowData$state_manager$getState()
            reqFn(current_s4)

            if (!inherits(current_s4, "LipidomicsAssayData")) {
                stop("Current state is not a LipidomicsAssayData object")
            }

            original_counts <- sapply(current_s4@lipid_data, function(assay) {
                if (current_s4@lipid_id_column %in% names(assay)) {
                    length(unique(assay[[current_s4@lipid_id_column]]))
                } else {
                    nrow(assay)
                }
            })

            logInfoFn(paste(
                "Applying lipid intensity filter: percentile ="
                , input$intensity_cutoff_percentile
                , ", proportion ="
                , input$proportion_below_cutoff
            ))

            filtered_s4 <- lipidIntensityFilteringFn(
                current_s4
                , lipids_intensity_cutoff_percentile = input$intensity_cutoff_percentile
                , lipids_proportion_of_samples_below_cutoff = input$proportion_below_cutoff
            )

            filtered_counts <- sapply(filtered_s4@lipid_data, function(assay) {
                if (filtered_s4@lipid_id_column %in% names(assay)) {
                    length(unique(assay[[filtered_s4@lipid_id_column]]))
                } else {
                    nrow(assay)
                }
            })

            stats_df <- data.frame(
                Assay = names(original_counts)
                , Original = as.numeric(original_counts)
                , Filtered = as.numeric(filtered_counts)
                , Removed = as.numeric(original_counts) - as.numeric(filtered_counts)
                , stringsAsFactors = FALSE
            )
            stats_df$Percent_Retained <- round(
                (stats_df$Filtered / stats_df$Original) * 100
                , 1
            )
            filterStatsVal(stats_df)

            qc_plot <- tryCatch({
                result <- updateLipidFilteringFn(
                    theObject = filtered_s4
                    , step_name = "2_Intensity_Filtered"
                    , omics_type = omicType
                    , return_grid = TRUE
                    , overwrite = TRUE
                )
                result
            }, error = function(e) {
                logWarnFn(paste("Could not generate QC plot:", e$message))
                NULL
            })

            filterPlotVal(qc_plot)

            workflowData$state_manager$saveState(
                state_name = "lipid_intensity_filtered"
                , s4_data_object = filtered_s4
                , config_object = workflowData$config_list
                , description = sprintf(
                    "Applied lipid intensity filter (percentile: %d%%, proportion: %.2f)"
                    , input$intensity_cutoff_percentile
                    , input$proportion_below_cutoff
                )
            )

            total_original <- sum(stats_df$Original)
            total_filtered <- sum(stats_df$Filtered)
            total_removed <- sum(stats_df$Removed)

            result_text <- paste(
                "Lipid Intensity Filter Applied Successfully"
                , "============================================"
                , sprintf("Intensity cutoff percentile: %d%%", input$intensity_cutoff_percentile)
                , sprintf("Max proportion below cutoff: %.2f", input$proportion_below_cutoff)
                , ""
                , "Per-Assay Results:"
                , paste(sapply(seq_len(nrow(stats_df)), function(i) {
                    sprintf(
                        "  %s: %d -> %d (removed %d, %.1f%% retained)"
                        , stats_df$Assay[i]
                        , stats_df$Original[i]
                        , stats_df$Filtered[i]
                        , stats_df$Removed[i]
                        , stats_df$Percent_Retained[i]
                    )
                }), collapse = "\n")
                , ""
                , sprintf("Total: %d -> %d lipids (removed %d)", total_original, total_filtered, total_removed)
                , ""
                , "State saved as: 'lipid_intensity_filtered'"
                , sep = "\n"
            )

            output$filter_results <- renderTextFn(result_text)

            logInfoFn(paste(
                "Lipid intensity filter applied: removed"
                , total_removed
                , "lipids"
            ))

            removeNotificationFn("lipid_intensity_filter_working")
            showNotificationFn(
                sprintf("Intensity filter applied: %d lipids retained", total_filtered)
                , type = "message"
            )
        }, error = function(e) {
            msg <- paste("Error applying lipid intensity filter:", e$message)
            logErrorFn(msg)
            showNotificationFn(msg, type = "error", duration = 15)
            removeNotificationFn("lipid_intensity_filter_working")
        })
    })

    invisible(input)
}

