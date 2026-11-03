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
                # KIT-075: select the saved state without materialising it (the proteomics GUI's call).
                selectWorkflowStateReference(workflowData$state_manager, prev_state_name)
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
            # KIT-075: the counts come from the declaration's result and the QC plot from the committed state read for
            # this action; no object outlives the observer (its outputs' closures keep this frame).
            current <- readLipidWorkflowMetadata(workflowData$state_manager)
            reqFn(current)

            if (!identical(current$class, "LipidomicsAssayData")) {
                stop("Current state is not a LipidomicsAssayData object")
            }

            logInfoFn(paste(
                "Applying lipid intensity filter: percentile ="
                , input$intensity_cutoff_percentile
                , ", proportion ="
                , input$proportion_below_cutoff
            ))

            # KIT-078: the tab's completion, from the declaration's counts in process or on the worker.
            finish <- function(counts, draw_qc_plot) {
                original_counts <- counts$before
                filtered_counts <- counts$after

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

                qc_plot <- tryCatch(draw_qc_plot(), error = function(e) {
                    logWarnFn(paste("Could not generate QC plot:", e$message))
                    NULL
                })

                filterPlotVal(qc_plot)

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
            }
            draw <- function(arguments) do.call(updateLipidFilteringFn, c(arguments, list(step_name = "2_Intensity_Filtered",
                omics_type = omicType, return_grid = TRUE, overwrite = TRUE)))

            # KIT-031: the filter is the declared operation `lipid_intensity_filter`, committed under
            # `lipid_intensity_filtered` with its provenance (KIT-072: the owner's path is the only path). KIT-078: on the
            # worker under a managed session, which also computes what the figures read from the matrix.
            runner <- workflowOmicsManagedActions("lipidomics")
            if (!is.null(runner)) {
                step <- omicsIntensityFilterStep(lipidomicsDeclarationBinding, input$intensity_cutoff_percentile, input$proportion_below_cutoff)
                runner$run(list(omicsManagedCommit(step, workflowData$config_list, function(result) {
                    finish(result$value$counts, function() {
                        metrics <- readWorkflowOmicsFilteringMetrics(workflowData$state_manager, result)
                        if (is.null(metrics)) stop("the worker computed no filtering metrics")
                        draw(list(theObject = NULL, metrics = metrics))
                    })
                })), failed = function(error) removeNotificationFn("lipid_intensity_filter_working"))
                return(invisible(NULL))
            }
            filter_run <- runLipidReferenceIntensityFilter(workflowData$state_manager, input$intensity_cutoff_percentile,
                input$proportion_below_cutoff, workflowData$config_list)
            finish(filter_run$value$counts, function() {
                draw(list(theObject = readLipidWorkflowState(workflowData$state_manager, filter_run$reference)))
            })
        }, error = function(e) {
            msg <- paste("Error applying lipid intensity filter:", e$message)
            logErrorFn(msg)
            showNotificationFn(msg, type = "error", duration = 15)
            removeNotificationFn("lipid_intensity_filter_working")
        })
    })

    invisible(input)
}

