# Import extracted components at top of file
source("ui/proteomics/qc/peptideQCApplet_ui.R")
source("server/proteomics/qc/peptideQCApplet_server.R")
source("ui/proteomics/qc/proteinQCApplet_ui.R")
source("server/proteomics/qc/proteinQCApplet_server.R")
source("R/shiny/server/proteomics/qc/protein/create_s4_from_protein_data_server.R")
source("R/shiny/ui/proteomics/qc/create_s4_from_protein_data_ui.R")

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
qualityControlAppletUI <- function(id, workflow_data) {
  ns <- NS(id)
  
  # Get workflow type to conditionally render UI
  workflow_type <- shiny::isolate(workflow_data$state_manager$workflow_type)
  
  # Define UI for LFQ/DIA workflows (peptide + protein tabs)
  lfq_dia_ui <- shiny::tabsetPanel(
    id = ns("qc_tabs_lfq"),
    shiny::tabPanel("Peptide QC", peptideQCAppletUI(ns("peptide_qc"))),
    shiny::tabPanel("Protein QC", proteinQCAppletUI(ns("protein_qc")))
  )
  
  # Define UI for TMT workflow (S4 creation + protein tabs)
  tmt_ui <- shiny::tabsetPanel(
    id = ns("qc_tabs_tmt"),
    shiny::tabPanel("Initial Protein Processing", create_s4_from_protein_data_ui(ns("create_s4_tmt"))),
    shiny::tabPanel("Protein QC", proteinQCAppletUI(ns("protein_qc")))
  )
  
  shiny::fluidRow(
    shiny::column(12,
      shiny::wellPanel(
        shiny::h3("Quality Control & Filtering"),
        
        # Conditionally render the correct UI based on workflow type
        if (workflow_type == "TMT") {
          tmt_ui
        } else {
          lfq_dia_ui
        }
      )
    )
  )
}

#' @rdname qualityControlAppletModule
#' @export
#' @import shiny
#' @importFrom logger log_info log_error
#' @importFrom grid grid.draw
qualityControlAppletServer <- function(id, workflow_data, experiment_paths, omic_type, experiment_label) {
  shiny::moduleServer(id, function(input, output, session) {
    
    # Create global project_dirs object that updateProteinFiltering expects
    # This maps our experiment_paths to the expected global structure
    if (!exists("project_dirs", envir = .GlobalEnv)) {
      project_dirs <- list()
      project_dirs[[paste0(omic_type, "_", experiment_label)]] <- experiment_paths
      assign("project_dirs", project_dirs, envir = .GlobalEnv)
      logger::log_info("Created global project_dirs object for updateProteinFiltering compatibility")
    }
    
    # == Coordinator-owned logic: Raw Data QC ===========================
    
    # Reactive value to store raw data QC plot
    raw_qc_plot <- reactiveVal(NULL)
    
    observeEvent(input$run_raw_data_qc, {
      shiny::req(workflow_data$state_manager)
      
      # Show a notification while running
      shiny::showNotification("Generating raw data QC plot...", id = "raw_qc_plot_working", duration = NULL)
      
      tryCatch({
        # 1. Retrieve the initial S4 object from the state manager
        logger::log_info("QC Step: Retrieving 'raw_data_s4' object from state manager.")
        s4_object <- workflow_data$state_manager$getState("raw_data_s4")
        
        shiny::req(s4_object)
        
        # 2. Run the updateProteinFiltering function
        # This function generates a plot and saves it, but we also want to capture it.
        logger::log_info("QC Step: Running updateProteinFiltering for '1_Raw Data'.")
        
        plot_grid <- updateProteinFiltering(
          data = s4_object, # Pass the S4 object directly
          step_name = "1_Raw Data",
          omic_type = omic_type,
          experiment_label = experiment_label,
          return_grid = TRUE,
          overwrite = TRUE
        )
        
        raw_qc_plot(plot_grid)
        
        logger::log_info("QC Step: Raw data QC plot generated successfully.")
        shiny::removeNotification("raw_qc_plot_working")
        shiny::showNotification("Plot generated.", type = "message")
        
      }, error = function(e) {
        msg <- paste("Error generating raw data QC plot:", e$message)
        logger::log_error(msg)
        shiny::showNotification(msg, type = "error", duration = 15)
        shiny::removeNotification("raw_qc_plot_working")
      })
    })
    
    # Render the raw data plot
    output$raw_data_plot <- shiny::renderPlot({
      shiny::req(raw_qc_plot())
      # The return object is a grob, so it needs to be drawn
      grid::grid.draw(raw_qc_plot())
    })
    
    # Revert to Raw Data functionality
    observeEvent(input$revert_to_raw, {
      tryCatch({
        # Revert to raw data state and clear all subsequent states
        reverted_s4 <- workflow_data$state_manager$revertToState("raw_data_s4")
        
        logger::log_info("Session reset to raw data state - all filtering steps cleared")
        shiny::showNotification("Session reset to raw data state. All filtering steps have been cleared.", 
          type = "message", duration = 5)
        
      }, error = function(e) {
        msg <- paste("Error resetting session:", e$message)
        logger::log_error(msg)
        shiny::showNotification(msg, type = "error")
      })
    })
    
    # == Component Server Calls ==========================================
    
    # Call the appropriate sub-applet orchestrators based on workflow type
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
    
    # Final step: update tab status when the final state is reached
    observeEvent(workflow_data$state_manager$states$protein_replicate_filtered, {
      # This observer is triggered when the final protein state is reached.
      # It can be used to update the UI or perform final actions.
      # For now, we'll just log it.
      logger::log_info("Final protein state reached. Quality Control module is ready.")
    })
    
  })
} 