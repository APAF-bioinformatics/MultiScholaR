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
NULL

#' @rdname mod_lipid_qc
#' @export
#' @importFrom shiny NS tagList fluidRow column wellPanel h3 uiOutput
mod_lipid_qc_ui <- function(id) {
    ns <- shiny::NS(id)
    
    shiny::tagList(
        shiny::fluidRow(
            shiny::column(12
                , shiny::wellPanel(
                    shiny::h3("Lipid Quality Control & Filtering")
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
mod_lipid_qc_server <- function(id, workflow_data, experiment_paths, omic_type, experiment_label, qc_trigger = NULL) {
    shiny::moduleServer(id, function(input, output, session) {
        ns <- session$ns
        
        # Track whether sub-modules have been initialized
        modules_initialized <- shiny::reactiveVal(FALSE)
        
        # Render dynamic QC tabs
        output$dynamic_qc_tabs <- shiny::renderUI({
            # Check if we have data available
            has_data <- tryCatch({
                !is.null(workflow_data$state_manager) && 
                !is.null(workflow_data$state_manager$getState())
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
                
                # Initialize all QC sub-modules
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
                
                modules_initialized(TRUE)
                logger::log_info("Lipidomics QC: Sub-modules initialized successfully")
            }
        }, ignoreNULL = TRUE, once = TRUE)
        
        # Alternative initialization when state manager becomes available
        shiny::observe({
            shiny::req(workflow_data$state_manager)
            
            # Check if we have a lipidomics S4 object
            current_state <- tryCatch({
                workflow_data$state_manager$getState()
            }, error = function(e) {
                NULL
            })
            
            if (!is.null(current_state) && 
                inherits(current_state, "LipidomicsAssayData") && 
                !modules_initialized()) {
                
                logger::log_info("Lipidomics QC: Auto-initializing sub-modules (state detected)")
                
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
                
                modules_initialized(TRUE)
            }
        })
    })
}

