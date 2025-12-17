# ============================================================================
# mod_metab_summary.R
# ============================================================================
# Purpose: Metabolomics session summary and integration export Shiny module
#
# This module provides a summary of the processing session and enables
# export of results for multi-omics integration.
# ============================================================================

#' @title Metabolomics Session Summary Module
#' @description A Shiny module for displaying processing summary and enabling
#'              data export for integration workflows.
#' @name mod_metab_summary
NULL

#' @rdname mod_metab_summary
#' @export
#' @importFrom shiny NS tagList fluidRow column wellPanel h3 h4 h5 p hr br uiOutput verbatimTextOutput actionButton downloadButton icon tags
mod_metab_summary_ui <- function(id) {
    ns <- shiny::NS(id)
    
    shiny::tagList(
        shiny::fluidRow(
            shiny::column(12
                , shiny::wellPanel(
                    shiny::h3("Session Summary & Export")
                    
                    , shiny::fluidRow(
                        # Left column: Processing summary
                        shiny::column(6
                            , shiny::h4("Processing Summary")
                            , shiny::uiOutput(ns("processing_summary"))
                            
                            , shiny::hr()
                            
                            , shiny::h4("State History")
                            , shiny::uiOutput(ns("state_history"))
                            
                            , shiny::hr()
                            
                            , shiny::h4("Data Summary")
                            , shiny::verbatimTextOutput(ns("data_summary"))
                        )
                        
                        # Right column: Export options
                        , shiny::column(6
                            , shiny::h4("Export Options")
                            , shiny::p("Export your processed data for downstream analysis or multi-omics integration.")
                            
                            , shiny::wellPanel(
                                style = "background-color: #f8f9fa;"
                                , shiny::h5(shiny::icon("table"), " Tabular Exports")
                                , shiny::br()
                                
                                , shiny::downloadButton(
                                    ns("download_design")
                                    , "Design Matrix (CSV)"
                                    , class = "btn-secondary"
                                    , style = "width: 100%; margin-bottom: 10px;"
                                )
                                
                                , shiny::downloadButton(
                                    ns("download_normalized")
                                    , "Normalized Data (CSV)"
                                    , class = "btn-secondary"
                                    , style = "width: 100%; margin-bottom: 10px;"
                                )
                                
                                , shiny::downloadButton(
                                    ns("download_de_results")
                                    , "DE Results (CSV)"
                                    , class = "btn-secondary"
                                    , style = "width: 100%;"
                                )
                            )
                            
                            , shiny::wellPanel(
                                style = "background-color: #e8f4fd;"
                                , shiny::h5(shiny::icon("database"), " S4 Object Export")
                                , shiny::br()
                                
                                , shiny::downloadButton(
                                    ns("download_s4")
                                    , "MetaboliteAssayData Object (RDS)"
                                    , class = "btn-primary"
                                    , style = "width: 100%; margin-bottom: 10px;"
                                )
                                
                                , shiny::p(
                                    shiny::tags$small(
                                        "The RDS file contains the complete S4 object for loading in future R sessions."
                                    )
                                )
                            )
                            
                            , shiny::wellPanel(
                                style = "background-color: #e8fde8;"
                                , shiny::h5(shiny::icon("link"), " Multi-Omics Integration")
                                , shiny::br()
                                
                                , shiny::downloadButton(
                                    ns("download_integration")
                                    , "Integration Export (RDS)"
                                    , class = "btn-success"
                                    , style = "width: 100%; margin-bottom: 10px;"
                                )
                                
                                , shiny::p(
                                    shiny::tags$small(
                                        "Formatted for use with MOFA or other multi-omics integration tools."
                                    )
                                )
                            )
                            
                            , shiny::hr()
                            
                            , shiny::h4("Processing Log")
                            , shiny::verbatimTextOutput(ns("processing_log_display"))
                        )
                    )
                )
            )
        )
    )
}

#' @rdname mod_metab_summary
#' @export
#' @importFrom shiny moduleServer reactiveValues reactive observeEvent req renderUI renderText downloadHandler tags
#' @importFrom logger log_info log_error
mod_metab_summary_server <- function(id, workflow_data, experiment_paths, omic_type, experiment_label) {
    shiny::moduleServer(id, function(input, output, session) {
        ns <- session$ns
        
        # Processing summary
        output$processing_summary <- shiny::renderUI({
            tab_status <- workflow_data$tab_status
            
            steps <- list(
                list(name = "Setup & Import", key = "setup_import")
                , list(name = "Design Matrix", key = "design_matrix")
                , list(name = "Quality Control", key = "quality_control")
                , list(name = "Normalization", key = "normalization")
                , list(name = "Differential Analysis", key = "differential_analysis")
            )
            
            step_items <- lapply(steps, function(step) {
                status <- tab_status[[step$key]]
                is_complete <- !is.null(status) && status == "complete"
                
                icon_name <- if (is_complete) "check-circle" else "circle"
                icon_color <- if (is_complete) "green" else "gray"
                
                shiny::tags$li(
                    shiny::icon(icon_name, style = paste("color:", icon_color))
                    , " "
                    , step$name
                    , if (is_complete) shiny::tags$span(" (Complete)", style = "color: green;") else NULL
                )
            })
            
            shiny::tags$ul(step_items, style = "list-style: none; padding-left: 0;")
        })
        
        # State history
        output$state_history <- shiny::renderUI({
            shiny::req(workflow_data$state_manager)
            
            history <- tryCatch({
                workflow_data$state_manager$getHistory()
            }, error = function(e) {
                character(0)
            })
            
            if (length(history) == 0) {
                return(shiny::p("No state history available.", style = "color: #666;"))
            }
            
            history_items <- lapply(seq_along(history), function(i) {
                shiny::tags$li(
                    shiny::tags$code(sprintf("%d. %s", i, history[i]))
                )
            })
            
            shiny::tags$ol(history_items, style = "font-size: 0.9em;")
        })
        
        # Data summary
        output$data_summary <- shiny::renderText({
            shiny::req(workflow_data$state_manager)
            
            current_s4 <- tryCatch({
                workflow_data$state_manager$getState()
            }, error = function(e) NULL)
            
            if (is.null(current_s4) || !inherits(current_s4, "MetaboliteAssayData")) {
                return("No MetaboliteAssayData object available.")
            }
            
            assay_list <- current_s4@metabolite_data
            dm <- current_s4@design_matrix
            
            # Assay info
            assay_info <- sapply(names(assay_list), function(name) {
                assay <- assay_list[[name]]
                n_features <- if (current_s4@metabolite_id_column %in% names(assay)) {
                    length(unique(assay[[current_s4@metabolite_id_column]]))
                } else {
                    nrow(assay)
                }
                sprintf("  %s: %d metabolites", name, n_features)
            })
            
            paste(
                "MetaboliteAssayData Object"
                , "=========================="
                , ""
                , sprintf("Number of assays: %d", length(assay_list))
                , paste(assay_info, collapse = "\n")
                , ""
                , sprintf("Design matrix samples: %d", nrow(dm))
                , sprintf("Groups: %s", paste(unique(dm[[current_s4@group_id]]), collapse = ", "))
                , ""
                , sprintf("Metabolite ID column: %s", current_s4@metabolite_id_column)
                , sprintf("Sample ID column: %s", current_s4@sample_id)
                , sep = "\n"
            )
        })
        
        # Processing log display
        output$processing_log_display <- shiny::renderText({
            log_entries <- workflow_data$processing_log
            
            if (length(log_entries) == 0) {
                return("No processing log available.")
            }
            
            log_text <- lapply(names(log_entries), function(step_name) {
                entry <- log_entries[[step_name]]
                timestamp <- if (!is.null(entry$timestamp)) {
                    format(entry$timestamp, "%Y-%m-%d %H:%M:%S")
                } else {
                    "Unknown"
                }
                
                paste0(
                    "[", step_name, "] ", timestamp
                    , if (!is.null(entry$description)) paste0("\n  ", entry$description) else ""
                )
            })
            
            paste(log_text, collapse = "\n\n")
        })
        
        # Download handlers
        
        # Design matrix
        output$download_design <- shiny::downloadHandler(
            filename = function() {
                paste0("metabolomics_design_matrix_", Sys.Date(), ".csv")
            }
            , content = function(file) {
                current_s4 <- workflow_data$state_manager$getState()
                shiny::req(current_s4)
                write.csv(current_s4@design_matrix, file, row.names = FALSE)
            }
        )
        
        # Normalized data
        output$download_normalized <- shiny::downloadHandler(
            filename = function() {
                paste0("metabolomics_normalized_", Sys.Date(), ".csv")
            }
            , content = function(file) {
                current_s4 <- workflow_data$state_manager$getState()
                shiny::req(current_s4)
                
                # Combine all assays with assay label
                all_data <- lapply(names(current_s4@metabolite_data), function(name) {
                    assay_data <- current_s4@metabolite_data[[name]]
                    assay_data$assay_name <- name
                    return(assay_data)
                })
                
                combined <- do.call(rbind, all_data)
                write.csv(combined, file, row.names = FALSE)
            }
        )
        
        # DE results
        output$download_de_results <- shiny::downloadHandler(
            filename = function() {
                paste0("metabolomics_de_results_", Sys.Date(), ".csv")
            }
            , content = function(file) {
                de_log <- workflow_data$processing_log$differential_analysis
                
                if (!is.null(de_log) && !is.null(de_log$results)) {
                    write.csv(de_log$results, file, row.names = FALSE)
                } else {
                    # Empty file with headers
                    write.csv(
                        data.frame(
                            message = "No DE results available. Run differential analysis first."
                        )
                        , file
                        , row.names = FALSE
                    )
                }
            }
        )
        
        # S4 object
        output$download_s4 <- shiny::downloadHandler(
            filename = function() {
                paste0("metabolomics_s4_object_", Sys.Date(), ".rds")
            }
            , content = function(file) {
                current_s4 <- workflow_data$state_manager$getState()
                shiny::req(current_s4)
                saveRDS(current_s4, file)
            }
        )
        
        # Integration export
        output$download_integration <- shiny::downloadHandler(
            filename = function() {
                paste0("metabolomics_integration_export_", Sys.Date(), ".rds")
            }
            , content = function(file) {
                current_s4 <- workflow_data$state_manager$getState()
                shiny::req(current_s4)
                
                # Create integration-ready format
                integration_data <- list(
                    omic_type = "metabolomics"
                    , data = current_s4@metabolite_data
                    , design_matrix = current_s4@design_matrix
                    , sample_id_col = current_s4@sample_id
                    , group_col = current_s4@group_id
                    , feature_id_col = current_s4@metabolite_id_column
                    , assay_names = names(current_s4@metabolite_data)
                    , processing_log = workflow_data$processing_log
                    , export_timestamp = Sys.time()
                )
                
                saveRDS(integration_data, file)
                
                logger::log_info("Exported metabolomics data for integration")
            }
        )
    })
}

