# ============================================================================
# mod_lipid_qc_duplicates.R
# ============================================================================
# Purpose: Lipid duplicate resolution Shiny module
#
# This module wraps the resolveDuplicateFeaturesByIntensity() function to
# provide an interactive UI for identifying and resolving duplicate lipids.
# ============================================================================

#' @title Lipid Duplicate Resolution Module
#' @description A Shiny module for identifying and resolving duplicate lipid features.
#'              Wraps the resolveDuplicateFeaturesByIntensity() function.
#' @name mod_lipid_qc_duplicates
NULL

#' @rdname mod_lipid_qc_duplicates
#' @export
#' @importFrom shiny NS tagList tabPanel br fluidRow column wellPanel h4 p hr div actionButton verbatimTextOutput uiOutput plotOutput
#' @importFrom DT DTOutput
#' @importFrom shinyjqui jqui_resizable
mod_lipid_qc_duplicates_ui <- function(id) {
    ns <- shiny::NS(id)

    shiny::tabPanel(
        "Duplicate Resolution"
        , shiny::br()
        , shiny::fluidRow(
            shiny::column(4
                , shiny::wellPanel(
                    shiny::h4("Duplicate Feature Resolution")
                    , shiny::p("Identify and resolve duplicate lipid IDs within each assay. Duplicates are resolved by keeping the feature with the highest average intensity across samples.")
                    , shiny::hr()

                    # Detection summary
                    , shiny::h5("Detection Summary")
                    , shiny::uiOutput(ns("duplicate_summary"))

                    , shiny::hr()
                    , shiny::div(
                        shiny::actionButton(
                            ns("detect_duplicates")
                            , "Detect Duplicates"
                            , class = "btn-info"
                            , width = "100%"
                        )
                    )
                    , shiny::br()
                    , shiny::div(
                        shiny::actionButton(
                            ns("resolve_duplicates")
                            , "Resolve Duplicates"
                            , class = "btn-primary"
                            , width = "48%"
                        )
                        , shiny::actionButton(
                            ns("revert_duplicates")
                            , "Revert"
                            , class = "btn-warning"
                            , width = "48%"
                            , style = "margin-left: 4%"
                        )
                    )
                )
            )
            , shiny::column(8
                , shiny::verbatimTextOutput(ns("resolution_results"))
                , shiny::br()
                # Per-assay duplicate tables
                , shiny::uiOutput(ns("duplicate_tables"))
                , shiny::br()
                # QC Progress Grid
                , shinyjqui::jqui_resizable(
                    shiny::plotOutput(ns("filter_plot"), height = "600px", width = "100%")
                )
            )
        )
    )
}

#' @rdname mod_lipid_qc_duplicates
#' @export
#' @importFrom shiny moduleServer reactiveVal observeEvent req showNotification removeNotification renderText renderUI tabsetPanel tabPanel tags renderPlot
#' @importFrom DT renderDT datatable
#' @importFrom logger log_info log_error log_warn
#' @importFrom grid grid.draw
mod_lipid_qc_duplicates_server <- function(id, workflow_data, omic_type, experiment_label) {
    shiny::moduleServer(id, function(input, output, session) {
        ns <- session$ns

        duplicate_info <- shiny::reactiveVal(NULL)
        resolution_stats <- shiny::reactiveVal(NULL)
        filter_plot <- shiny::reactiveVal(NULL)
        
        # Detect duplicates
        shiny::observeEvent(input$detect_duplicates, {
            shiny::req(workflow_data$state_manager)
            
            tryCatch({
                current_s4 <- workflow_data$state_manager$getState()
                shiny::req(current_s4)
                
                if (!inherits(current_s4, "LipidomicsAssayData")) {
                    stop("Current state is not a LipidomicsAssayData object")
                }
                
                # Use the findDuplicateFeatureIDs function
                duplicates_list <- findDuplicateFeatureIDs(current_s4)
                duplicate_info(duplicates_list)
                
                # Count total duplicates
                total_duplicates <- sum(sapply(duplicates_list, function(x) {
                    if (is.null(x)) 0 else nrow(x)
                }))
                
                logger::log_info(paste("Detected", total_duplicates, "duplicate feature IDs across assays"))
                
                shiny::showNotification(
                    sprintf("Detection complete: %d duplicate IDs found", total_duplicates)
                    , type = if (total_duplicates > 0) "warning" else "message"
                )
                
            }, error = function(e) {
                msg <- paste("Error detecting duplicates:", e$message)
                logger::log_error(msg)
                shiny::showNotification(msg, type = "error")
            })
        })
        
        # Render duplicate summary
        output$duplicate_summary <- shiny::renderUI({
            dup_list <- duplicate_info()
            
            if (is.null(dup_list)) {
                return(shiny::p(
                    shiny::icon("info-circle")
                    , " Click 'Detect Duplicates' to scan for duplicate features."
                    , style = "color: #666;"
                ))
            }
            
            summary_items <- lapply(names(dup_list), function(assay_name) {
                dup_df <- dup_list[[assay_name]]
                n_dups <- if (is.null(dup_df)) 0 else nrow(dup_df)
                
                icon_class <- if (n_dups == 0) "check-circle" else "exclamation-triangle"
                icon_color <- if (n_dups == 0) "green" else "orange"
                
                shiny::tags$li(
                    shiny::icon(icon_class, style = paste("color:", icon_color))
                    , sprintf(" %s: %d duplicate IDs", assay_name, n_dups)
                )
            })
            
            shiny::tags$ul(summary_items, style = "list-style: none; padding-left: 0;")
        })
        
        # Render duplicate tables per assay
        output$duplicate_tables <- shiny::renderUI({
            dup_list <- duplicate_info()
            
            if (is.null(dup_list)) {
                return(NULL)
            }
            
            # Filter to assays with duplicates
            assays_with_dups <- names(dup_list)[sapply(dup_list, function(x) !is.null(x) && nrow(x) > 0)]
            
            if (length(assays_with_dups) == 0) {
                return(shiny::wellPanel(
                    shiny::icon("check-circle", style = "color: green; font-size: 24px;")
                    , shiny::h5("No duplicates found in any assay!")
                    , shiny::p("All lipid IDs are unique. No resolution needed.")
                ))
            }
            
            tab_list <- lapply(assays_with_dups, function(assay_name) {
                shiny::tabPanel(
                    assay_name
                    , shiny::br()
                    , DT::DTOutput(ns(paste0("dup_table_", gsub("[^a-zA-Z0-9]", "_", assay_name))))
                )
            })
            
            do.call(shiny::tabsetPanel, c(list(id = ns("dup_tables_tabs")), tab_list))
        })
        
        # Render individual duplicate tables
        shiny::observe({
            dup_list <- duplicate_info()
            shiny::req(dup_list)
            
            lapply(names(dup_list), function(assay_name) {
                dup_df <- dup_list[[assay_name]]
                if (!is.null(dup_df) && nrow(dup_df) > 0) {
                    output_id <- paste0("dup_table_", gsub("[^a-zA-Z0-9]", "_", assay_name))
                    output[[output_id]] <- DT::renderDT({
                        DT::datatable(
                            dup_df
                            , options = list(
                                pageLength = 10
                                , scrollX = TRUE
                                , dom = 'frtip'
                            )
                            , rownames = FALSE
                            , class = "compact stripe"
                        )
                    })
                }
            })
        })
        
        # Resolve duplicates
        shiny::observeEvent(input$resolve_duplicates, {
            shiny::req(workflow_data$state_manager)
            
            shiny::showNotification(
                "Resolving duplicate features..."
                , id = "lipid_dup_resolve_working"
                , duration = NULL
            )
            
            tryCatch({
                current_s4 <- workflow_data$state_manager$getState()
                shiny::req(current_s4)
                
                if (!inherits(current_s4, "LipidomicsAssayData")) {
                    stop("Current state is not a LipidomicsAssayData object")
                }
                
                lipid_id_col <- current_s4@lipid_id_column
                assay_list <- current_s4@lipid_data
                assay_names <- names(assay_list)
                
                # Track stats
                stats_list <- list()
                
                # Process each assay
                resolved_assay_list <- lapply(seq_along(assay_list), function(i) {
                    assay_data <- assay_list[[i]]
                    assay_name <- if (!is.null(assay_names)) assay_names[i] else paste0("Assay_", i)
                    
                    original_rows <- nrow(assay_data)
                    
                    # Identify sample (numeric) columns
                    sample_cols <- names(assay_data)[sapply(assay_data, is.numeric)]
                    
                    if (length(sample_cols) == 0) {
                        logger::log_warn(paste("No numeric columns found in assay:", assay_name))
                        stats_list[[assay_name]] <<- list(
                            original = original_rows
                            , resolved = original_rows
                            , removed = 0
                        )
                        return(assay_data)
                    }
                    
                    # Resolve duplicates
                    resolved_assay <- resolveDuplicateFeaturesByIntensity(
                        assay_tibble = assay_data
                        , id_col = lipid_id_col
                        , sample_cols = sample_cols
                    )
                    
                    resolved_rows <- nrow(resolved_assay)
                    
                    stats_list[[assay_name]] <<- list(
                        original = original_rows
                        , resolved = resolved_rows
                        , removed = original_rows - resolved_rows
                    )
                    
                    return(resolved_assay)
                })
                
                # Restore names
                names(resolved_assay_list) <- assay_names
                
                # Update the S4 object
                current_s4@lipid_data <- resolved_assay_list
                
                resolution_stats(stats_list)
                
                # Save state
                workflow_data$state_manager$saveState(
                    state_name = "lipid_duplicates_resolved"
                    , s4_data_object = current_s4
                    , config_object = workflow_data$config_list
                    , description = "Resolved duplicate lipid features by keeping highest intensity"
                )

                # Update QC tracking visualization
                qc_plot <- tryCatch({
                    updateLipidFiltering(
                        theObject = current_s4
                        , step_name = "3_Duplicates_Resolved"
                        , omics_type = omic_type
                        , return_grid = TRUE
                        , overwrite = TRUE
                    )
                }, error = function(e) {
                    logger::log_warn(paste("Could not generate QC plot:", e$message))
                    NULL
                })
                filter_plot(qc_plot)

                # Generate summary text
                total_removed <- sum(sapply(stats_list, function(x) x$removed))
                
                result_text <- paste(
                    "Duplicate Resolution Complete"
                    , "=============================="
                    , "Strategy: Keep feature with highest average intensity"
                    , ""
                    , "Per-Assay Results:"
                    , paste(sapply(names(stats_list), function(name) {
                        s <- stats_list[[name]]
                        sprintf("  %s: %d → %d rows (removed %d duplicates)"
                            , name, s$original, s$resolved, s$removed)
                    }), collapse = "\n")
                    , ""
                    , sprintf("Total duplicate rows removed: %d", total_removed)
                    , ""
                    , "State saved as: 'lipid_duplicates_resolved'"
                    , sep = "\n"
                )
                
                output$resolution_results <- shiny::renderText(result_text)
                
                # Clear duplicate info after resolution
                duplicate_info(NULL)
                
                logger::log_info(paste("Resolved duplicates: removed", total_removed, "rows"))
                
                shiny::removeNotification("lipid_dup_resolve_working")
                shiny::showNotification(
                    sprintf("Duplicates resolved: %d rows removed", total_removed)
                    , type = "message"
                )
                
            }, error = function(e) {
                msg <- paste("Error resolving duplicates:", e$message)
                logger::log_error(msg)
                shiny::showNotification(msg, type = "error", duration = 15)
                shiny::removeNotification("lipid_dup_resolve_working")
            })
        })
        
        # Revert
        shiny::observeEvent(input$revert_duplicates, {
            tryCatch({
                history <- workflow_data$state_manager$getHistory()
                if (length(history) > 1) {
                    prev_state_name <- history[length(history) - 1]
                    workflow_data$state_manager$revertToState(prev_state_name)
                    output$resolution_results <- shiny::renderText(
                        paste("Reverted to previous state:", prev_state_name)
                    )
                    resolution_stats(NULL)
                    duplicate_info(NULL)
                    filter_plot(NULL)
                    logger::log_info(paste("Reverted duplicate resolution to", prev_state_name))
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

        # Render QC progress plot
        output$filter_plot <- shiny::renderPlot({
            shiny::req(filter_plot())
            if (inherits(filter_plot(), "grob") || inherits(filter_plot(), "gtable")) {
                grid::grid.draw(filter_plot())
            } else if (inherits(filter_plot(), "ggplot")) {
                print(filter_plot())
            }
        })
    })
}

