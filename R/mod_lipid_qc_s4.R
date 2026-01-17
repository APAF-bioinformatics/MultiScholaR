# ============================================================================
# mod_lipid_qc_s4.R
# ============================================================================
# Purpose: Final S4 state confirmation module for lipidomics QC
#
# This module provides a summary of QC processing and allows the user to
# finalize the QC step before proceeding to normalization.
# ============================================================================

#' @title Lipidomics QC Finalization Module
#' @description A Shiny module for finalizing the QC step and confirming the S4 state.
#'              Displays QC summary statistics and state history.
#' @name mod_lipid_qc_s4
NULL

#' @rdname mod_lipid_qc_s4
#' @export
#' @importFrom shiny NS tagList tabPanel br fluidRow column wellPanel h4 h5 p hr div actionButton verbatimTextOutput uiOutput icon tags plotOutput
#' @importFrom DT DTOutput
#' @importFrom shinyjqui jqui_resizable
mod_lipid_qc_s4_ui <- function(id) {
    ns <- shiny::NS(id)
    
    shiny::tabPanel(
        "Finalize QC"
        , shiny::br()
        , shiny::fluidRow(
            shiny::column(6
                , shiny::wellPanel(
                    shiny::h4("QC Processing Summary")
                    , shiny::p("Review the QC processing steps before finalizing. You can revert to any previous state if needed.")
                    , shiny::hr()
                    
                    # State history
                    , shiny::h5("Processing History")
                    , shiny::uiOutput(ns("state_history"))
                    
                    , shiny::hr()
                    , shiny::div(
                        shiny::actionButton(
                            ns("finalize_qc")
                            , "Finalize QC"
                            , class = "btn-success"
                            , width = "100%"
                            , icon = shiny::icon("check")
                        )
                    )
                )
            )
            , shiny::column(6
                , shiny::wellPanel(
                    shiny::h4("Current Data Summary")
                    , shiny::hr()
                    , shiny::uiOutput(ns("data_summary"))
                )
                , shiny::wellPanel(
                    shiny::h4("Per-Assay Statistics")
                    , shiny::hr()
                    , DT::DTOutput(ns("assay_stats_table"))
                )
            )
        )
        , shiny::fluidRow(
            shiny::column(12
                , shiny::verbatimTextOutput(ns("finalize_results"))
                , shiny::br()
                # QC Progress Grid
                , shinyjqui::jqui_resizable(
                    shiny::plotOutput(ns("filter_plot"), height = "600px", width = "100%")
                )
            )
        )
    )
}

#' @rdname mod_lipid_qc_s4
#' @export
#' @importFrom shiny moduleServer reactiveVal observeEvent req showNotification renderText renderUI observe tags renderPlot
#' @importFrom DT renderDT datatable
#' @importFrom logger log_info log_error
#' @importFrom grid grid.draw
mod_lipid_qc_s4_server <- function(id, workflow_data, omic_type, experiment_label) {
    shiny::moduleServer(id, function(input, output, session) {
        ns <- session$ns

        filter_plot <- shiny::reactiveVal(NULL)

        # Render state history
        output$state_history <- shiny::renderUI({
            shiny::req(workflow_data$state_manager)
            
            history <- tryCatch({
                workflow_data$state_manager$getHistory()
            }, error = function(e) {
                character(0)
            })
            
            if (length(history) == 0) {
                return(shiny::p(
                    shiny::icon("info-circle")
                    , " No processing history available."
                    , style = "color: #666;"
                ))
            }
            
            # Create history items with icons
            history_items <- lapply(seq_along(history), function(i) {
                state_name <- history[i]
                is_current <- i == length(history)
                
                icon_name <- if (is_current) "arrow-right" else "check"
                icon_color <- if (is_current) "blue" else "green"
                text_style <- if (is_current) "font-weight: bold;" else ""
                
                shiny::tags$li(
                    shiny::icon(icon_name, style = paste("color:", icon_color))
                    , shiny::tags$span(
                        sprintf(" %d. %s", i, state_name)
                        , style = text_style
                    )
                    , if (is_current) shiny::tags$span(" (current)", style = "color: blue;") else NULL
                )
            })
            
            shiny::tags$ol(history_items, style = "list-style: none; padding-left: 0;")
        })
        
        # Render data summary
        output$data_summary <- shiny::renderUI({
            shiny::req(workflow_data$state_manager)
            
            current_s4 <- tryCatch({
                workflow_data$state_manager$getState()
            }, error = function(e) {
                NULL
            })
            
            if (is.null(current_s4) || !inherits(current_s4, "LipidomicsAssayData")) {
                return(shiny::p(
                    shiny::icon("exclamation-triangle", style = "color: orange;")
                    , " No LipidomicsAssayData object available."
                ))
            }
            
            # Calculate summary stats
            assay_list <- current_s4@lipid_data
            n_assays <- length(assay_list)
            
            total_lipids <- sum(sapply(assay_list, function(a) {
                if (current_s4@lipid_id_column %in% names(a)) {
                    length(unique(a[[current_s4@lipid_id_column]]))
                } else {
                    nrow(a)
                }
            }))
            
            # Get sample count from first assay
            sample_cols <- if (n_assays > 0) {
                names(assay_list[[1]])[sapply(assay_list[[1]], is.numeric)]
            } else {
                character(0)
            }
            n_samples <- length(sample_cols)
            
            # Design matrix info
            dm <- current_s4@design_matrix
            n_groups <- if (nrow(dm) > 0 && current_s4@group_id %in% names(dm)) {
                length(unique(dm[[current_s4@group_id]]))
            } else {
                NA
            }
            
            shiny::tagList(
                shiny::tags$table(
                    class = "table table-condensed"
                    , style = "margin-bottom: 0;"
                    , shiny::tags$tbody(
                        shiny::tags$tr(
                            shiny::tags$td(shiny::strong("Number of Assays:"))
                            , shiny::tags$td(n_assays)
                        )
                        , shiny::tags$tr(
                            shiny::tags$td(shiny::strong("Total Lipids:"))
                            , shiny::tags$td(format(total_lipids, big.mark = ","))
                        )
                        , shiny::tags$tr(
                            shiny::tags$td(shiny::strong("Number of Samples:"))
                            , shiny::tags$td(n_samples)
                        )
                        , shiny::tags$tr(
                            shiny::tags$td(shiny::strong("Experimental Groups:"))
                            , shiny::tags$td(if (!is.na(n_groups)) n_groups else "N/A")
                        )
                        , shiny::tags$tr(
                            shiny::tags$td(shiny::strong("Lipid ID Column:"))
                            , shiny::tags$td(current_s4@lipid_id_column)
                        )
                        , shiny::tags$tr(
                            shiny::tags$td(shiny::strong("Sample ID Column:"))
                            , shiny::tags$td(current_s4@sample_id)
                        )
                    )
                )
            )
        })
        
        # Render per-assay statistics table
        output$assay_stats_table <- DT::renderDT({
            shiny::req(workflow_data$state_manager)
            
            current_s4 <- tryCatch({
                workflow_data$state_manager$getState()
            }, error = function(e) {
                NULL
            })
            
            if (is.null(current_s4) || !inherits(current_s4, "LipidomicsAssayData")) {
                return(NULL)
            }
            
            assay_list <- current_s4@lipid_data
            lipid_id_col <- current_s4@lipid_id_column
            
            # Build stats data frame
            stats_df <- data.frame(
                Assay = names(assay_list)
                , Lipids = sapply(assay_list, function(a) {
                    if (lipid_id_col %in% names(a)) {
                        length(unique(a[[lipid_id_col]]))
                    } else {
                        nrow(a)
                    }
                })
                , Samples = sapply(assay_list, function(a) {
                    sum(sapply(a, is.numeric))
                })
                , stringsAsFactors = FALSE
            )
            
            # Calculate missingness per assay
            stats_df$Missingness <- sapply(assay_list, function(a) {
                quant_cols <- names(a)[sapply(a, is.numeric)]
                if (length(quant_cols) == 0) return(NA)
                
                total_cells <- nrow(a) * length(quant_cols)
                missing_cells <- sum(sapply(quant_cols, function(col) {
                    sum(is.na(a[[col]]) | a[[col]] == 0)
                }))
                
                round((missing_cells / total_cells) * 100, 1)
            })
            
            DT::datatable(
                stats_df
                , options = list(
                    dom = 't'
                    , paging = FALSE
                    , ordering = FALSE
                )
                , rownames = FALSE
                , class = "compact stripe"
            )
        })
        
        # Finalize QC
        shiny::observeEvent(input$finalize_qc, {
            shiny::req(workflow_data$state_manager)
            
            tryCatch({
                current_s4 <- workflow_data$state_manager$getState()
                shiny::req(current_s4)
                
                if (!inherits(current_s4, "LipidomicsAssayData")) {
                    stop("Current state is not a LipidomicsAssayData object")
                }
                
                # Save final QC state
                workflow_data$state_manager$saveState(
                    state_name = "lipid_qc_complete"
                    , s4_data_object = current_s4
                    , config_object = workflow_data$config_list
                    , description = "QC processing complete - ready for normalization"
                )

                # Update QC tracking visualization

                qc_plot <- tryCatch({
                    result <- updateLipidFiltering(
                        theObject = current_s4
                        , step_name = "4_QC_Complete"
                        , omics_type = omic_type
                        , return_grid = TRUE
                        , overwrite = TRUE
                    )
                    if (!is.null(result)) {
                    }
                    result
                }, error = function(e) {
                    NULL
                })

                filter_plot(qc_plot)

                # Update tab status - must replace entire list to trigger reactivity
                updated_status <- workflow_data$tab_status
                updated_status$quality_control <- "complete"
                workflow_data$tab_status <- updated_status
                
                # Generate summary
                assay_list <- current_s4@lipid_data
                total_lipids <- sum(sapply(assay_list, function(a) {
                    if (current_s4@lipid_id_column %in% names(a)) {
                        length(unique(a[[current_s4@lipid_id_column]]))
                    } else {
                        nrow(a)
                    }
                }))
                
                history <- workflow_data$state_manager$getHistory()
                
                result_text <- paste(
                    "QC Finalization Complete"
                    , "========================"
                    , ""
                    , sprintf("Total lipids retained: %d", total_lipids)
                    , sprintf("Number of assays: %d", length(assay_list))
                    , sprintf("Processing steps completed: %d", length(history))
                    , ""
                    , "Processing History:"
                    , paste(sprintf("  %d. %s", seq_along(history), history), collapse = "\n")
                    , ""
                    , "State saved as: 'lipid_qc_complete'"
                    , ""
                    , "You can now proceed to the Normalization tab."
                    , sep = "\n"
                )
                
                output$finalize_results <- shiny::renderText(result_text)
                
                logger::log_info("Lipidomics QC finalized successfully")
                
                shiny::showNotification(
                    "QC complete! Proceed to Normalization."
                    , type = "message"
                    , duration = 5
                )
                
            }, error = function(e) {
                msg <- paste("Error finalizing QC:", e$message)
                logger::log_error(msg)
                shiny::showNotification(msg, type = "error")
            })
        })

        # Render QC progress plot
        output$filter_plot <- shiny::renderPlot({

            shiny::req(filter_plot())


            if (inherits(filter_plot(), "grob") || inherits(filter_plot(), "gtable")) {
                grid::grid.draw(filter_plot())
            } else if (inherits(filter_plot(), "ggplot")) {
                print(filter_plot())
            } else {
            }
        })
    })
}

