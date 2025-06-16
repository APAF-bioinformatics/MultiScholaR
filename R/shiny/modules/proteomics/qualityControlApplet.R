#' @title qualityControlAppletModule
#'
#' @description A Shiny module for the Quality Control step of the proteomics
#' workflow. It contains sub-tabs for different QC operations.
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
          
          # == Raw Data QC Sub-Tab ========================================
          shiny::tabPanel(
            "Raw Data Overview",
            icon = shiny::icon("images"),
            shiny::br(),
            shiny::fluidRow(
              shiny::column(3,
                shiny::wellPanel(
                  shiny::h4("Raw Data Visualization"),
                  shiny::p("This generates a plot summarizing the number of proteins identified at various stages of the raw data processing."),
                  shiny::actionButton(ns("run_raw_data_qc"), "Generate Plot", class = "btn-primary", width = "100%")
                )
              ),
              shiny::column(9,
                shiny::plotOutput(ns("raw_data_plot"), height = "600px")
              )
            )
          ),
          
          # == Peptide Filtering Sub-Tab =================================
          shiny::tabPanel(
            "Peptide Filtering",
            icon = shiny::icon("filter")
            # Placeholder for future content
          ),
          
          # == Protein Filtering Sub-Tab =================================
          shiny::tabPanel(
            "Protein Filtering",
            icon = shiny::icon("filter")
            # Placeholder for future content
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
    
    # Reactive value to store the plot from the raw data QC step
    raw_qc_plot <- reactiveVal(NULL)
    
    # == Raw Data QC Logic ==================================================
    
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
        
        # The function needs 'data_cln' from the global env, which is not ideal.
        # For now, we get it from workflow_data.
        # A better long-term solution would be to refactor updateProteinFiltering to accept the S4 object.
        plot_grid <- updateProteinFiltering(
          data = s4_object, # Pass the S4 object directly
          step_name = "1_Raw Data",
          publication_graphs_dir = NULL, # Disable saving for now
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
    
    # Render the plot
    output$raw_data_plot <- shiny::renderPlot({
      shiny::req(raw_qc_plot())
      # The return object is a grob, so it needs to be drawn
      grid::grid.draw(raw_qc_plot())
    })
    
  })
} 