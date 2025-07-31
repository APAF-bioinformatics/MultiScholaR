# Import extracted components at top of file
source("ui/proteomics/qc/peptideQCApplet_ui.R")
source("server/proteomics/qc/peptideQCApplet_server.R")
source("ui/proteomics/qc/proteinQCApplet_ui.R")
source("server/proteomics/qc/proteinQCApplet_server.R")

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
        shiny::h3("Quality Control"),
        shiny::p("This section provides tools for assessing and filtering the data based on quality metrics."),
        
        shiny::tabsetPanel(
          id = ns("qc_tabs"),
          
          # == Raw Data QC Sub-Tab (Coordinator-owned) ==================
          shiny::tabPanel(
            "Raw Data Overview",
            icon = shiny::icon("images"),
            shiny::br(),
            shiny::fluidRow(
              shiny::column(3,
                shiny::wellPanel(
                  shiny::h4("Raw Data Visualization"),
                  shiny::p("This generates a plot summarizing the number of proteins identified at various stages of the raw data processing."),
                  shiny::actionButton(ns("run_raw_data_qc"), "Generate Plot", class = "btn-primary", width = "100%"),
                  shiny::hr(),
                  shiny::h5("Session Reset"),
                  shiny::p("Reset the analysis back to the raw data state. This will clear all filtering steps."),
                  shiny::actionButton(ns("revert_to_raw"), "Revert to Raw Data", 
                    class = "btn-warning", width = "100%")
                )
              ),
              shiny::column(9,
                shinyjqui::jqui_resizable(
                  shiny::plotOutput(ns("raw_data_plot"), height = "800px", width = "100%")
                )
              )
            )
          ),
          
          # == Peptide Filtering Component ===============================
          shiny::tabPanel(
            "Peptide Filtering",
            icon = shiny::icon("filter"),
            shiny::br(),
            peptideQCAppletUI(ns("peptide"))
          ),
          
          # == Protein Filtering Component ===============================
          shiny::tabPanel(
            "Protein Filtering",
            icon = shiny::icon("filter"),
            shiny::br(),
            proteinQCAppletUI(ns("protein"))
          )
        )
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
    
    # Peptide filtering component server
    peptideQCAppletServer("peptide", workflow_data, experiment_paths, omic_type, experiment_label)
    
    # Protein filtering component server
    proteinQCAppletServer("protein", workflow_data, experiment_paths, omic_type, experiment_label)
    
  })
} 