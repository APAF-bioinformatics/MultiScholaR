# Import extracted components at top of file
source("ui/proteomics/qc/peptideQCApplet_ui.R")
source("server/proteomics/qc/peptideQCApplet_server.R")
source("ui/proteomics/qc/proteinQCApplet_ui.R")
source("server/proteomics/qc/proteinQCApplet_server.R")
source("server/proteomics/qc/protein/create_s4_from_protein_data_server.R")
source("ui/proteomics/qc/create_s4_from_protein_data_ui.R")

#' @title qualityControlAppletModule
#'
#' @description A Shiny module coordinator for the Quality Control step of the proteomics
#' workflow. This is now a lightweight coordinator that integrates extracted components
#' while maintaining backward compatibility.
#'
#' @name qualityControlAppletModule
NULL

#' @rdname qualityControlAppletModule
#' @export
#' @import shiny
#' @import shinydashboard
qualityControlAppletUI <- function(id) {
  ns <- NS(id)
  shiny::fluidRow(
    shiny::column(12,
      shiny::wellPanel(
        shiny::h3("Quality Control & Filtering"),
        # Replace the conditional logic with a single UI Output
        shiny::uiOutput(ns("dynamic_qc_tabs"))
      )
    )
  )
}

#' @rdname qualityControlAppletModule
#' @export
#' @import shiny
#' @importFrom logger log_info log_error
#' @importFrom grid grid.draw
qualityControlAppletServer <- function(id, workflow_data, experiment_paths, omic_type, experiment_label, qc_trigger = NULL) {
  shiny::moduleServer(id, function(input, output, session) {
    ns <- session$ns
    
    output$dynamic_qc_tabs <- shiny::renderUI({
      # This logic will run once workflow_type is available
      shiny::req(workflow_data$state_manager$workflow_type)
      workflow_type <- workflow_data$state_manager$workflow_type
      
      if (workflow_type == "TMT") {
        # Define and return UI for TMT workflow
        shiny::tabsetPanel(
          id = ns("qc_tabs_tmt"),
          shiny::tabPanel("Initial Protein Processing", create_s4_from_protein_data_ui(ns("create_s4_tmt"))),
          shiny::tabPanel("Protein QC", proteinQCAppletUI(ns("protein_qc")))
        )
      } else { # LFQ and DIA
        # Define and return UI for peptide-based workflows
        shiny::tabsetPanel(
          id = ns("qc_tabs_lfq"),
          shiny::tabPanel("Peptide QC", peptideQCAppletUI(ns("peptide_qc"))),
          shiny::tabPanel("Protein QC", proteinQCAppletUI(ns("protein_qc")))
        )
      }
    })
    
    # == Component Server Calls ==========================================
    
    # Call the appropriate sub-applet orchestrators based on workflow type
    # This now waits for the qc_trigger from the design matrix step.
    observeEvent(qc_trigger(), {
      req(qc_trigger() == TRUE) # Only run when trigger is explicitly set to TRUE
      
      workflow_type <- shiny::isolate(workflow_data$state_manager$workflow_type)
      logger::log_info(paste("Main QC Applet: Routing for workflow type:", workflow_type))
      
      if (workflow_type %in% c("LFQ", "DIA")) {
        # The peptide applet orchestrator will handle all peptide-level steps
        peptideQCAppletServer("peptide_qc", workflow_data, experiment_paths, omic_type, experiment_label)
        
        # The protein applet orchestrator will handle all protein-level steps
        proteinQCAppletServer("protein_qc", workflow_data, experiment_paths, omic_type, experiment_label)
        
      } else if (workflow_type == "TMT") {
        # For TMT, we first create the S4 object from the protein data
        create_s4_from_protein_data_server("create_s4_tmt", workflow_data, omic_type, experiment_label)
        
        # Then, we call the protein QC orchestrator, which will run the common steps
        proteinQCAppletServer("protein_qc", workflow_data, experiment_paths, omic_type, experiment_label)
      }
    }, ignoreNULL = TRUE, once = TRUE) # Run only once when triggered
    
    # Final step: update tab status when the final state is reached
    observeEvent(workflow_data$state_manager$states$protein_replicate_filtered, {
      # This observer is triggered when the final protein state is reached.
      # It can be used to update the UI or perform final actions.
      # For now, we'll just log it.
      logger::log_info("Final protein state reached. Quality Control module is ready.")
    })
    
  })
} 