#' Design Matrix Applet UI
#' 
#' UI for integrating the existing design matrix applet
#' 
#' @param id Module ID
#' 
#' @importFrom shiny NS actionButton uiOutput wellPanel h4 DTOutput
#' @importFrom shiny fluidRow column conditionalPanel div h3 p icon br
#' @export
designMatrixUi <- function(id) {
  ns <- shiny::NS(id)
  
  shiny::fluidRow(
    shiny::column(12,
      shiny::wellPanel(
        shiny::h3("Design Matrix Builder"),
        shiny::p("The design matrix defines your experimental groups and contrasts for differential expression analysis."),
        
        shiny::conditionalPanel(
          condition = paste0("output['", ns("data_available"), "']"),
          shiny::fluidRow(
            shiny::column(12,
              shiny::actionButton(
                ns("launch_design_matrix"),
                "Launch Design Matrix Builder",
                class = "btn-primary btn-lg",
                icon = shiny::icon("edit"),
                width = "100%"
              ),
              shiny::br(),
              shiny::br(),
              shiny::conditionalPanel(
                condition = paste0("output['", ns("design_matrix_exists"), "']"),
                shiny::wellPanel(
                  shiny::h4("Current Design Matrix"),
                  DT::DTOutput(ns("design_matrix_preview")),
                  shiny::br(),
                  shiny::h4("Defined Contrasts"),
                  DT::DTOutput(ns("contrasts_preview"))
                )
              )
            )
          )
        ),
        
        shiny::conditionalPanel(
          condition = paste0("!output['", ns("data_available"), "']"),
          shiny::div(
            class = "alert alert-info",
            shiny::icon("info-circle"),
            " Please complete the Setup & Import step first."
          )
        )
      )
    )
  )
}

#' Design Matrix Applet Server
#' 
#' Server logic for integrating the existing design matrix applet
#' 
#' @param id Module ID
#' @param workflow_data Reactive values object to store workflow data
#' @param experiment_paths List of paths for this experiment
#' @param omic_type The omics type
#' @param experiment_label The experiment label
#' 
#' @importFrom shiny moduleServer reactive observeEvent req renderDT
#' @importFrom shiny showNotification outputOptions
#' @importFrom logger log_info log_error
#' @export
designMatrixServer <- function(id, workflow_data, experiment_paths, omic_type, experiment_label) {
  shiny::moduleServer(id, function(input, output, session) {
    
    # Check if data is available
    output$data_available <- shiny::reactive({
      !is.null(workflow_data$data_tbl) && !is.null(workflow_data$config_list)
    })
    shiny::outputOptions(output, "data_available", suspendWhenHidden = FALSE)
    
    # Check if design matrix exists
    output$design_matrix_exists <- shiny::reactive({
      !is.null(workflow_data$design_matrix)
    })
    shiny::outputOptions(output, "design_matrix_exists", suspendWhenHidden = FALSE)
    
    # Launch design matrix applet
    shiny::observeEvent(input$launch_design_matrix, {
      shiny::req(workflow_data$data_tbl)
      shiny::req(workflow_data$config_list)
      
      log_info("Launching design matrix applet for {omic_type} experiment: {experiment_label}")
      
      # Ensure required objects exist in the parent environment
      # The RunApplet function expects these in the parent frame
      assign("data_tbl", workflow_data$data_tbl, envir = parent.frame(2))
      assign("config_list", workflow_data$config_list, envir = parent.frame(2))
      assign("project_dirs", experiment_paths, envir = parent.frame(2))
      
      tryCatch({
        # Check if RunApplet exists
        if (!exists("RunApplet")) {
          stop("RunApplet function not found. Make sure shiny_applets.R is loaded.")
        }
        
        # Call the existing RunApplet function
        result <- RunApplet(
          applet_type = "designMatrix",
          omic_type = omic_type,
          experiment_label = experiment_label,
          project_dirs_object_name = "project_dirs",
          force = FALSE
        )
        
        # Update workflow data with results
        if (!is.null(result)) {
          workflow_data$design_matrix <- result$design_matrix
          workflow_data$data_cln <- result$data_cln
          workflow_data$contrasts_tbl <- result$contrasts_tbl
          workflow_data$config_list <- result$config_list
          
          # Update processing log
          workflow_data$processing_log$design_matrix <- list(
            timestamp = Sys.time(),
            n_samples = nrow(result$design_matrix),
            n_groups = length(unique(result$design_matrix$group)),
            n_contrasts = if (!is.null(result$contrasts_tbl)) nrow(result$contrasts_tbl) else 0
          )
          
          # Mark this step as complete
          workflow_data$tab_status$design_matrix <- "complete"
          
          shiny::showNotification("Design matrix saved successfully!", type = "success")
        }
        
      }, error = function(e) {
        log_error("Error running design matrix applet: {e$message}")
        shiny::showNotification(
          paste("Error launching design matrix builder:", e$message),
          type = "error",
          duration = NULL
        )
      })
    })
    
    # Display design matrix preview
    output$design_matrix_preview <- DT::renderDT({
      shiny::req(workflow_data$design_matrix)
      workflow_data$design_matrix
    }, options = list(pageLength = 10, scrollX = TRUE))
    
    # Display contrasts preview
    output$contrasts_preview <- DT::renderDT({
      shiny::req(workflow_data$contrasts_tbl)
      workflow_data$contrasts_tbl
    }, options = list(pageLength = 10, scrollX = TRUE))
    
  })
} 