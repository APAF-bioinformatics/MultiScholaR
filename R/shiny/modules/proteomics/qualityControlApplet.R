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
      message("=== DEBUG66: qualityControlApplet renderUI fired ===")
      
      # This logic will run once workflow_type is available
      shiny::req(workflow_data$state_manager$workflow_type)
      workflow_type <- workflow_data$state_manager$workflow_type
      message(sprintf("   DEBUG66: workflow_type = %s", workflow_type))
      
      if (workflow_type == "TMT") {
        # Check if S4 object already exists from design matrix step
        s4_exists <- tryCatch({
          s4_obj <- workflow_data$state_manager$getState("protein_s4_initial")
          exists <- !is.null(s4_obj)
          message(sprintf("   DEBUG66: TMT - Checking for protein_s4_initial: exists = %s", exists))
          if (exists) {
            message(sprintf("   DEBUG66: TMT - S4 object class = %s", class(s4_obj)))
          }
          exists
        }, error = function(e) {
          message(sprintf("   DEBUG66: TMT - Error checking S4: %s", e$message))
          FALSE
        })
        
        message(sprintf("   DEBUG66: TMT - s4_exists final = %s", s4_exists))
        
        if (s4_exists) {
          # S4 already created in design matrix step, skip to protein QC
          logger::log_info("QC UI: S4 object exists from design matrix. Skipping S4 creation UI.")
          message("   DEBUG66: TMT - Rendering single Protein QC tab")
          shiny::tabsetPanel(
            id = ns("qc_tabs_tmt"),
            shiny::tabPanel("Protein QC", proteinQCAppletUI(ns("protein_qc"), workflow_type = workflow_type))
          )
        } else {
          # Fallback: show S4 creation UI if somehow not created
          logger::log_info("QC UI: S4 object not found. Showing S4 creation UI.")
          message("   DEBUG66: TMT - Rendering S4 creation + Protein QC tabs")
          shiny::tabsetPanel(
            id = ns("qc_tabs_tmt"),
            shiny::tabPanel("Initial Protein Processing", create_s4_from_protein_data_ui(ns("create_s4_tmt"))),
            shiny::tabPanel("Protein QC", proteinQCAppletUI(ns("protein_qc"), workflow_type = workflow_type))
          )
        }
      } else { # LFQ and DIA
        message(sprintf("   DEBUG66: %s - Rendering Peptide + Protein QC tabs", workflow_type))
        # Define and return UI for peptide-based workflows
        shiny::tabsetPanel(
          id = ns("qc_tabs_lfq"),
          shiny::tabPanel("Peptide QC", peptideQCAppletUI(ns("peptide_qc"))),
          shiny::tabPanel("Protein QC", proteinQCAppletUI(ns("protein_qc"), workflow_type = workflow_type))
        )
      }
    })
    
    # == Component Server Calls ==========================================
    
    # Call the appropriate sub-applet orchestrators based on workflow type
    # This now waits for the qc_trigger from the design matrix step.
    observeEvent(qc_trigger(), {
      message("=== DEBUG66: qualityControlApplet observeEvent(qc_trigger) FIRED ===")
      message(sprintf("   DEBUG66: qc_trigger() value = %s", qc_trigger()))
      
      req(qc_trigger() == TRUE) # Only run when trigger is explicitly set to TRUE
      message("   DEBUG66: qc_trigger is TRUE, proceeding with module activation")
      
      workflow_type <- shiny::isolate(workflow_data$state_manager$workflow_type)
      logger::log_info(paste("Main QC Applet: Routing for workflow type:", workflow_type))
      message(sprintf("   DEBUG66: Isolated workflow_type = %s", workflow_type))
      
      if (workflow_type %in% c("LFQ", "DIA")) {
        message(sprintf("   DEBUG66: %s workflow - calling peptide and protein servers", workflow_type))
        # The peptide applet orchestrator will handle all peptide-level steps
        peptideQCAppletServer("peptide_qc", workflow_data, experiment_paths, omic_type, experiment_label)
        message("   DEBUG66: peptideQCAppletServer called")
        
        # The protein applet orchestrator will handle all protein-level steps
        proteinQCAppletServer("protein_qc", workflow_data, experiment_paths, omic_type, experiment_label)
        message("   DEBUG66: proteinQCAppletServer called")
        
      } else if (workflow_type == "TMT") {
        message("   DEBUG66: TMT workflow - checking S4 existence")
        # Check if S4 object already exists from design matrix step
        s4_exists <- tryCatch({
          s4_obj <- workflow_data$state_manager$getState("protein_s4_initial")
          exists <- !is.null(s4_obj)
          message(sprintf("   DEBUG66: TMT observeEvent - S4 exists = %s", exists))
          exists
        }, error = function(e) {
          message(sprintf("   DEBUG66: TMT observeEvent - Error checking S4: %s", e$message))
          FALSE
        })
        
        if (s4_exists) {
          # S4 already created in design matrix step, go straight to protein QC
          logger::log_info("TMT: S4 object exists from design matrix. Skipping S4 creation, starting protein QC.")
          message("   DEBUG66: TMT - S4 exists, calling proteinQCAppletServer only")
          proteinQCAppletServer("protein_qc", workflow_data, experiment_paths, omic_type, experiment_label)
          message("   DEBUG66: TMT - proteinQCAppletServer called")
        } else {
          # S4 doesn't exist (shouldn't happen in normal flow), use the creation module
          logger::log_warn("TMT: S4 object not found. Using S4 creation module.")
          message("   DEBUG66: TMT - S4 MISSING, calling create_s4 + protein servers")
          create_s4_from_protein_data_server("create_s4_tmt", workflow_data, omic_type, experiment_label)
          message("   DEBUG66: TMT - create_s4_from_protein_data_server called")
          proteinQCAppletServer("protein_qc", workflow_data, experiment_paths, omic_type, experiment_label)
          message("   DEBUG66: TMT - proteinQCAppletServer called")
        }
      }
      message("=== DEBUG66: qualityControlApplet observeEvent complete ===")
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