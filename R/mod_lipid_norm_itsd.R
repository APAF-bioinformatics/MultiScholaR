# MultiScholaR: Interactive Multi-Omics Analysis
# Copyright (C) 2024-2026 Ignatius Pang, William Klare, and APAF-bioinformatics
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

# ============================================================================
# mod_lipid_norm_itsd.R
# ============================================================================
# Purpose: Lipidomics ITSD selection and normalization sub-module
# ============================================================================

#' @title Lipidomics ITSD Module
#' @description A Shiny module for selecting internal standards and applying ITSD normalization.
#' @name mod_lipid_norm_itsd
NULL

#' @rdname mod_lipid_norm_itsd
#' @export
#' @importFrom shiny NS tagList tabPanel br fluidRow column wellPanel h5 hr uiOutput icon helpText
mod_lipid_norm_itsd_ui <- function(id) {
    ns <- shiny::NS(id)
    
    shiny::tabPanel(
        "ITSD Selection"
        , icon = shiny::icon("flask")
        , shiny::br()
        , shiny::fluidRow(
            shiny::column(12
                , shiny::h5("Internal Standard Selection per Assay", style = "text-align: center;")
                , shiny::helpText("Select internal standards for each assay. Pre-selected candidates are auto-detected from naming patterns (ITSD, IS, labeled compounds).")
                , shiny::hr()
                , shiny::uiOutput(ns("itsd_selection_ui"))
            )
        )
    )
}

#' @rdname mod_lipid_norm_itsd
#' @export
#' @importFrom shiny moduleServer reactiveValues observe req renderUI wellPanel h5 br observeEvent
#' @importFrom DT renderDataTable datatable formatStyle styleEqual formatRound
#' @importFrom purrr map walk imap compact
#' @importFrom logger log_info
mod_lipid_norm_itsd_server <- function(id, workflow_data, omic_type, experiment_label) {
    shiny::moduleServer(id, function(input, output, session) {
        ns <- session$ns
        
        # Reactive values to store selections
        itsd_data <- shiny::reactiveValues(
            assay_names = NULL,
            selections = list()
        )
        
        # Detect assay names and initialize selections
        shiny::observe({
            shiny::req(workflow_data$state_manager)
            current_s4 <- workflow_data$state_manager$getState()
            if (!is.null(current_s4) && inherits(current_s4, "LipidomicsAssayData")) {
                itsd_data$assay_names <- names(current_s4@lipid_data)
                
                # Initialize selections for new assays
                for (assay_name in itsd_data$assay_names) {
                    if (!assay_name %in% names(itsd_data$selections)) {
                        itsd_data$selections[[assay_name]] <- NULL
                    }
                }
            }
        })
        
        # Render Selection UIs
        output$itsd_selection_ui <- shiny::renderUI({
            shiny::req(itsd_data$assay_names)
            
            assay_uis <- purrr::map(itsd_data$assay_names, \(assay_name) {
                safe_name <- gsub("[^A-Za-z0-9]", "_", tolower(assay_name))
                shiny::wellPanel(
                    shiny::h5(paste("Assay:", assay_name))
                    , DT::dataTableOutput(ns(paste0("itsd_table_", safe_name)))
                    , shiny::br()
                )
            })
            shiny::tagList(assay_uis)
        })
        
        # Render DT Tables and Handle Selections
        shiny::observe({
            shiny::req(itsd_data$assay_names)
            shiny::req(workflow_data$state_manager)
            
            current_s4 <- workflow_data$state_manager$getState()
            if (is.null(current_s4)) return()
            
            purrr::walk(itsd_data$assay_names, \(assay_name) {
                safe_name <- gsub("[^A-Za-z0-9]", "_", tolower(assay_name))
                output_id <- paste0("itsd_table_", safe_name)
                input_id <- paste0("itsd_table_", safe_name, "_rows_selected")
                
                # Table Rendering
                output[[output_id]] <- DT::renderDataTable({
                    assay_data <- current_s4@lipid_data[[assay_name]]
                    if (is.null(assay_data)) return(NULL)
                    
                    selection_table <- buildItsdSelectionTable(
                        assay_data = assay_data,
                        metabolite_id_col = current_s4@lipid_id_column,
                        annotation_cols = current_s4@annotation_id_column
                    )
                    
                    preselected <- which(selection_table$is_candidate)
                    
                    DT::datatable(
                        selection_table
                        , selection = list(mode = "multiple", selected = preselected)
                        , filter = "top"
                        , options = list(pageLength = 10, scrollX = TRUE)
                        , rownames = FALSE
                    ) |>
                        DT::formatStyle("is_candidate", backgroundColor = DT::styleEqual(TRUE, "#d4edda")) |>
                        DT::formatRound(columns = c("mean_intensity", "cv_percent"), digits = 2)
                })
                
                # Selection Tracking
                shiny::observeEvent(input[[input_id]], {
                    itsd_data$selections[[assay_name]] <- input[[input_id]]
                }, ignoreNULL = FALSE)
            })
        })
        
        # Return selections and helper to parent
        return(list(
            get_selections = shiny::reactive({ itsd_data$selections }),
            get_feature_ids = function(current_s4) {
                if (is.null(itsd_data$selections) || length(itsd_data$selections) == 0) return(NULL)
                
                ids <- purrr::imap(itsd_data$selections, \(row_indices, assay_name) {
                    if (is.null(row_indices) || length(row_indices) == 0) return(NULL)
                    assay_data <- current_s4@lipid_data[[assay_name]]
                    if (is.null(assay_data)) return(NULL)
                    
                    selection_table <- buildItsdSelectionTable(
                        assay_data = assay_data,
                        metabolite_id_col = current_s4@lipid_id_column,
                        annotation_cols = current_s4@annotation_id_column
                    )
                    selection_table$feature_id[row_indices]
                })
                purrr::compact(ids)
            }
        ))
    })
}

# <!-- APAF Bioinformatics | mod_lipid_norm_itsd.R | Approved | 2026-03-17 -->
