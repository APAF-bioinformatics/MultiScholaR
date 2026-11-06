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
# mod_lipid_qc.R
# ============================================================================
# Purpose: Lipidomics QC orchestrator module
#
# This module coordinates all QC sub-modules (intensity filtering, duplicate
# resolution, internal standards, and finalization) in a tabbed interface.
# ============================================================================

#' @title Lipidomics Quality Control Orchestrator Module
#' @description A Shiny module coordinator for the Quality Control step of the
#'              lipidomics workflow. Integrates sub-modules for intensity filtering,
#'              duplicate resolution, internal standard QC, and finalization.
#' @name mod_lipid_qc
#' @param id,workflow_data,omic_type,experiment_label Runtime inputs used by this function; see the usage section for accepted values.
#' @param profile An omic profile (KIT-096): its title, sub-tabs, their servers and class. `NULL`, the default, is lipidomics.
NULL

#' @rdname mod_lipid_qc
#' @export
#' @importFrom shiny NS tagList fluidRow column wellPanel h3 uiOutput
mod_lipid_qc_ui <- function(id, profile = NULL) {
    ns <- shiny::NS(id)
    
    shiny::tagList(
        shiny::fluidRow(
            shiny::column(12
                , shiny::wellPanel(
                    shiny::h3(profile$qc_title %||% "Lipid Quality Control & Filtering")
                    , shiny::uiOutput(ns("dynamic_qc_tabs"))
                )
            )
        )
    )
}

#' @rdname mod_lipid_qc
#' @export
#' @importFrom shiny moduleServer reactive observeEvent req renderUI tabsetPanel tabPanel
#' @importFrom logger log_info log_error
initializeLipidQcSubmodules <- function(workflow_data, omic_type, experiment_label) {
    mod_lipid_qc_intensity_server(
        "intensity"
        , workflow_data
        , omic_type
        , experiment_label
    )
    
    mod_lipid_qc_duplicates_server(
        "duplicates"
        , workflow_data
        , omic_type
        , experiment_label
    )
    
    mod_lipid_qc_itsd_server(
        "itsd"
        , workflow_data
        , omic_type
        , experiment_label
    )
    
    mod_lipid_qc_s4_server(
        "s4_finalize"
        , workflow_data
        , omic_type
        , experiment_label
    )
}

mod_lipid_qc_server <- function(id, workflow_data, experiment_paths, omic_type, experiment_label, qc_trigger = NULL, profile = NULL) {
    shiny::moduleServer(id, function(input, output, session) {
        ns <- session$ns
        # KIT-096: a profile's metadata reader, sub-tabs, sub-module servers and class (lipidomics' otherwise).
        readMetadata <- profile$read_metadata %||% readLipidWorkflowMetadata
        initializeSubmodules <- profile$qc_init %||% initializeLipidQcSubmodules
        
        # Track whether sub-modules have been initialized
        modules_initialized <- shiny::reactiveVal(FALSE)
        
        # Render dynamic QC tabs
        output$dynamic_qc_tabs <- shiny::renderUI({
            # Check if we have data available
            # KIT-075: the state's metadata answers whether there is data; no object is read.
            has_data <- tryCatch({
                !is.null(workflow_data$state_manager) && 
                !is.null(readMetadata(workflow_data$state_manager))
            }, error = function(e) {
                FALSE
            })
            
            if (!has_data) {
                return(shiny::div(
                    class = "alert alert-info"
                    , shiny::icon("info-circle")
                    , " Please complete the 'Design Matrix' step first. QC tabs will appear here once data is available."
                ))
            }
            
            if (!is.null(profile)) return(profile$qc_tabs(ns))

            # Build the QC tabset
            shiny::tabsetPanel(
                id = ns("lipid_qc_tabs")
                , mod_lipid_qc_intensity_ui(ns("intensity"))
                , mod_lipid_qc_duplicates_ui(ns("duplicates"))
                , mod_lipid_qc_itsd_ui(ns("itsd"))
                , mod_lipid_qc_s4_ui(ns("s4_finalize"))
            )
        })
        
        # Initialize sub-module servers when triggered
        shiny::observeEvent(qc_trigger(), {
            shiny::req(qc_trigger() == TRUE)
            
            if (!modules_initialized()) {
                logger::log_info("Lipidomics QC: Initializing sub-module servers")
                initializeSubmodules(workflow_data, omic_type, experiment_label)
                modules_initialized(TRUE)
                logger::log_info("Lipidomics QC: Sub-modules initialized successfully")
            }
        }, ignoreNULL = TRUE, once = TRUE)
        
        # Alternative initialization when state manager becomes available
        shiny::observe({
            shiny::req(workflow_data$state_manager)
            
            # Check if we have a lipidomics S4 object
            current_state <- tryCatch({
                readMetadata(workflow_data$state_manager)
            }, error = function(e) {
                NULL
            })
            
            if (!is.null(current_state) && 
                identical(current_state$class, profile$class %||% "LipidomicsAssayData") && 
                !shiny::isolate(modules_initialized())) {
                
                logger::log_info("Lipidomics QC: Auto-initializing sub-modules (state detected)")
                initializeSubmodules(workflow_data, omic_type, experiment_label)
                modules_initialized(TRUE)
            }
        })
    })
}
