# ============================================================================
# mod_metab_qc_intensity.R
# ============================================================================
# Purpose: Metabolite intensity filtering Shiny module
#
# This module wraps the metaboliteIntensityFiltering() S4 method to provide
# an interactive UI for filtering metabolites based on intensity thresholds.
# ============================================================================

#' @title Metabolite Intensity Filter Module
#' @description A Shiny module for applying metabolite intensity and missing value filters.
#'              Wraps the metaboliteIntensityFiltering() S4 method.
#' @name mod_metab_qc_intensity
NULL

#' @rdname mod_metab_qc_intensity
#' @export
#' @importFrom shiny NS tagList tabPanel br fluidRow column wellPanel h4 p hr numericInput sliderInput div actionButton verbatimTextOutput uiOutput
#' @importFrom shinyjqui jqui_resizable
mod_metab_qc_intensity_ui <- function(id) {
    ns <- shiny::NS(id)
    
    shiny::tabPanel(
        "Intensity Filter"
        , shiny::br()
        , shiny::fluidRow(
            shiny::column(4
                , shiny::wellPanel(
                    shiny::h4("Metabolite Intensity Filter Parameters")
                    , shiny::p("Filter metabolites based on intensity thresholds. Metabolites with low intensities across many samples will be removed.")
                    , shiny::hr()
                    
                    # Intensity cutoff percentile
                    , shiny::sliderInput(
                        ns("intensity_cutoff_percentile")
                        , "Intensity Cutoff Percentile (%)"
                        , min = 1
                        , max = 50
                        , value = 10
                        , step = 1
                    )
                    , shiny::helpText("Intensities below this percentile are considered 'low' (default: 10%)")
                    
                    # Proportion threshold
                    , shiny::sliderInput(
                        ns("proportion_below_cutoff")
                        , "Max Proportion of Samples Below Threshold"
                        , min = 0.1
                        , max = 1.0
                        , value = 0.5
                        , step = 0.05
                    )
                    , shiny::helpText("Remove metabolites where this proportion of samples are below the intensity threshold (default: 0.5)")
                    
                    , shiny::hr()
                    , shiny::div(
                        shiny::actionButton(
                            ns("apply_filter")
                            , "Apply Filter"
                            , class = "btn-primary"
                            , width = "48%"
                        )
                        , shiny::actionButton(
                            ns("revert_filter")
                            , "Revert"
                            , class = "btn-warning"
                            , width = "48%"
                            , style = "margin-left: 4%"
                        )
                    )
                )
            )
            , shiny::column(8
                , shiny::verbatimTextOutput(ns("filter_results"))
                , shiny::br()
                # Per-assay results tabs
                , shiny::uiOutput(ns("assay_results_tabs"))
                , shiny::br()
                , shinyjqui::jqui_resizable(
                    shiny::plotOutput(ns("filter_plot"), height = "600px", width = "100%")
                )
            )
        )
    )
}

#' @rdname mod_metab_qc_intensity
#' @export
#' @importFrom shiny moduleServer reactiveVal observeEvent req showNotification removeNotification renderText renderPlot renderUI tabsetPanel tabPanel
#' @importFrom logger log_info log_error
#' @importFrom grid grid.draw
mod_metab_qc_intensity_server <- function(id, workflow_data, omic_type, experiment_label) {
    shiny::moduleServer(id, function(input, output, session) {
        ns <- session$ns
        
        filter_plot <- shiny::reactiveVal(NULL)
        filter_stats <- shiny::reactiveVal(NULL)
        
        # Apply intensity filter
        shiny::observeEvent(input$apply_filter, {
            shiny::req(workflow_data$state_manager)
            
            shiny::showNotification(
                "Applying metabolite intensity filter..."
                , id = "metab_intensity_filter_working"
                , duration = NULL
            )
            
            tryCatch({
                # Get current MetaboliteAssayData S4 object from state
                current_s4 <- workflow_data$state_manager$getState()
                shiny::req(current_s4)
                
                # Verify it's a MetaboliteAssayData object
                if (!inherits(current_s4, "MetaboliteAssayData")) {
                    stop("Current state is not a MetaboliteAssayData object")
                }
                
                # Get original metabolite counts per assay
                original_counts <- sapply(current_s4@metabolite_data, function(assay) {
                    if (current_s4@metabolite_id_column %in% names(assay)) {
                        length(unique(assay[[current_s4@metabolite_id_column]]))
                    } else {
                        nrow(assay)
                    }
                })
                
                logger::log_info(paste(
                    "Applying metabolite intensity filter: percentile ="
                    , input$intensity_cutoff_percentile
                    , ", proportion ="
                    , input$proportion_below_cutoff
                ))
                
                # Apply S4 method for filtering
                filtered_s4 <- metaboliteIntensityFiltering(
                    current_s4
                    , metabolites_intensity_cutoff_percentile = input$intensity_cutoff_percentile
                    , metabolites_proportion_of_samples_below_cutoff = input$proportion_below_cutoff
                )
                
                # Get filtered metabolite counts per assay
                filtered_counts <- sapply(filtered_s4@metabolite_data, function(assay) {
                    if (filtered_s4@metabolite_id_column %in% names(assay)) {
                        length(unique(assay[[filtered_s4@metabolite_id_column]]))
                    } else {
                        nrow(assay)
                    }
                })
                
                # Store stats for display
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
                filter_stats(stats_df)
                
                # Update QC tracking visualization
                
                qc_plot <- tryCatch({
                    result <- updateMetaboliteFiltering(
                        theObject = filtered_s4
                        , step_name = "2_Intensity_Filtered"
                        , omics_type = omic_type
                        , return_grid = TRUE
                        , overwrite = TRUE
                    )
                    result
                }, error = function(e) {
                    logger::log_warn(paste("Could not generate QC plot:", e$message))
                    NULL
                })
                
                filter_plot(qc_plot)
                
                # Save state
                workflow_data$state_manager$saveState(
                    state_name = "metab_intensity_filtered"
                    , s4_data_object = filtered_s4
                    , config_object = workflow_data$config_list
                    , description = sprintf(
                        "Applied metabolite intensity filter (percentile: %d%%, proportion: %.2f)"
                        , input$intensity_cutoff_percentile
                        , input$proportion_below_cutoff
                    )
                )
                
                # REMOVED duplicate filter_plot(qc_plot) - was causing potential issues
                
                # Generate summary text
                total_original <- sum(stats_df$Original)
                total_filtered <- sum(stats_df$Filtered)
                total_removed <- sum(stats_df$Removed)
                
                result_text <- paste(
                    "Metabolite Intensity Filter Applied Successfully"
                    , "============================================"
                    , sprintf("Intensity cutoff percentile: %d%%", input$intensity_cutoff_percentile)
                    , sprintf("Max proportion below cutoff: %.2f", input$proportion_below_cutoff)
                    , ""
                    , "Per-Assay Results:"
                    , paste(sapply(seq_len(nrow(stats_df)), function(i) {
                        sprintf(
                            "  %s: %d → %d (removed %d, %.1f%% retained)"
                            , stats_df$Assay[i]
                            , stats_df$Original[i]
                            , stats_df$Filtered[i]
                            , stats_df$Removed[i]
                            , stats_df$Percent_Retained[i]
                        )
                    }), collapse = "\n")
                    , ""
                    , sprintf("Total: %d → %d metabolites (removed %d)", total_original, total_filtered, total_removed)
                    , ""
                    , "State saved as: 'metab_intensity_filtered'"
                    , sep = "\n"
                )
                
                output$filter_results <- shiny::renderText(result_text)
                
                logger::log_info(paste(
                    "Metabolite intensity filter applied: removed"
                    , total_removed
                    , "metabolites"
                ))
                
                shiny::removeNotification("metab_intensity_filter_working")
                shiny::showNotification(
                    sprintf("Intensity filter applied: %d metabolites retained", total_filtered)
                    , type = "message"
                )
                
            }, error = function(e) {
                msg <- paste("Error applying metabolite intensity filter:", e$message)
                logger::log_error(msg)
                shiny::showNotification(msg, type = "error", duration = 15)
                shiny::removeNotification("metab_intensity_filter_working")
            })
        })
        
        # Revert filter
        shiny::observeEvent(input$revert_filter, {
            tryCatch({
                history <- workflow_data$state_manager$getHistory()
                if (length(history) > 1) {
                    prev_state_name <- history[length(history) - 1]
                    workflow_data$state_manager$revertToState(prev_state_name)
                    output$filter_results <- shiny::renderText(
                        paste("Reverted to previous state:", prev_state_name)
                    )
                    filter_stats(NULL)
                    filter_plot(NULL)
                    logger::log_info(paste("Reverted metabolite intensity filter to", prev_state_name))
                    shiny::showNotification("Reverted successfully", type = "message")
                } else {
                    stop("No previous state to revert to.")
                }
            }, error = function(e) {
                msg <- paste("Error reverting:", e$message)
                logger::log_error(msg)
                shiny::showNotification(msg, type = "error")
            })
        })
        
        # Render per-assay results tabs
        output$assay_results_tabs <- shiny::renderUI({
            stats <- filter_stats()
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
                                shiny::tags$td(shiny::strong("Original metabolites:"))
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
        
        # Render filter plot
        output$filter_plot <- shiny::renderPlot({
            plot_obj <- filter_plot()
            shiny::req(plot_obj)
            if (inherits(plot_obj, "grob") || inherits(plot_obj, "gtable")) {
                grid::grid.draw(plot_obj)
            } else if (inherits(plot_obj, "ggplot")) {
                print(plot_obj)
            } else {
            }
        })
    })
}

