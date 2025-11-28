#' @title Peptide Precursor Rollup Module
#'
#' @description A Shiny module for applying precursor to peptide rollup.
#'
#' @name mod_prot_qc_peptide_rollup
NULL

#' @rdname mod_prot_qc_peptide_rollup
#' @export
#' @importFrom shiny NS tagList tabPanel br fluidRow column wellPanel h4 p hr div actionButton verbatimTextOutput plotOutput
mod_prot_qc_peptide_rollup_ui <- function(id) {
  ns <- shiny::NS(id)
  
  shiny::tabPanel(
    "Precursor Rollup",
    shiny::br(),
    shiny::fluidRow(
      shiny::column(4,
        shiny::wellPanel(
          shiny::h4("Precursor to Peptide Rollup"),
          shiny::p("Aggregate intensity measurements from multiple precursor ions to peptide level."),
          shiny::hr(),
          shiny::p("This step has no configurable parameters - it uses statistical methods to combine precursor measurements."),
          shiny::hr(),
          shiny::div(
            shiny::actionButton(ns("apply_rollup"), "Apply Rollup", 
              class = "btn-primary", width = "48%"),
            shiny::actionButton(ns("revert_rollup"), "Revert", 
              class = "btn-warning", width = "48%", style = "margin-left: 4%")
          )
        )
      ),
      shiny::column(8,
        shiny::verbatimTextOutput(ns("rollup_results")),
        shiny::br(),
        shinyjqui::jqui_resizable(
          shiny::plotOutput(ns("rollup_plot"), height = "800px", width = "100%")
        )
      )
    )
  )
}

#' @rdname mod_prot_qc_peptide_rollup
#' @export
#' @importFrom shiny moduleServer reactiveVal observeEvent req showNotification removeNotification renderText renderPlot
#' @importFrom logger log_info log_error
#' @importFrom grid grid.draw
mod_prot_qc_peptide_rollup_server <- function(id, workflow_data, omic_type, experiment_label) {
  shiny::moduleServer(id, function(input, output, session) {
    
    rollup_plot <- shiny::reactiveVal(NULL)
    
    # Step 2: Precursor Rollup (chunk 11)
    shiny::observeEvent(input$apply_rollup, {
      shiny::req(workflow_data$state_manager)
      
      shiny::showNotification("Applying precursor rollup...", id = "rollup_working", duration = NULL)
      
      tryCatch({
        # Get current S4 object from the active state
        current_s4 <- workflow_data$state_manager$getState()
        shiny::req(current_s4)
        
        logger::log_info("QC Step: Applying precursor rollup")
        
        # Apply S4 transformation (EXISTING S4 CODE - UNCHANGED)
        rolled_up_s4 <- rollUpPrecursorToPeptide(current_s4)
        
        # Track QC parameters in workflow_data
        if (is.null(workflow_data$qc_params)) {
          workflow_data$qc_params <- list()
        }
        if (is.null(workflow_data$qc_params$peptide_qc)) {
          workflow_data$qc_params$peptide_qc <- list()
        }
        
        workflow_data$qc_params$peptide_qc$precursor_rollup <- list(
          applied = TRUE,
          timestamp = Sys.time()
        )
        
        # Save new state in R6 manager
        workflow_data$state_manager$saveState(
          state_name = "precursor_rollup",
          s4_data_object = rolled_up_s4,
          config_object = list(),
          description = "Applied precursor to peptide rollup"
        )
        
        # Generate summary
        protein_count <- rolled_up_s4@peptide_data |>
          dplyr::distinct(Protein.Ids) |>
          nrow()
        
        result_text <- paste(
          "Precursor Rollup Applied Successfully\n",
          "====================================\n",
          sprintf("Proteins remaining: %d\n", protein_count),
          "State saved as: 'precursor_rollup'\n"
        )
        
        output$rollup_results <- shiny::renderText(result_text)
        
        # Update filtering visualization and capture plot
        plot_grid <- updateProteinFiltering(
          data = rolled_up_s4@peptide_data,
          step_name = "3_precursor_rollup",
          omic_type = omic_type,
          experiment_label = experiment_label,
          return_grid = TRUE,
          overwrite = TRUE
        )
        rollup_plot(plot_grid)
        
        logger::log_info("Precursor rollup applied successfully")
        shiny::removeNotification("rollup_working")
        shiny::showNotification("Precursor rollup applied successfully", type = "message")
        
      }, error = function(e) {
        msg <- paste("Error applying precursor rollup:", e$message)
        logger::log_error(msg)
        shiny::showNotification(msg, type = "error", duration = 15)
        shiny::removeNotification("rollup_working")
      })
    })
    
    # Revert Precursor Rollup
    shiny::observeEvent(input$revert_rollup, {
      tryCatch({
        # Revert to the previous state in the history
        history <- workflow_data$state_manager$getHistory()
        # The previous state is the one before the current one.
        if (length(history) > 1) {
          prev_state_name <- history[length(history) - 1]
          reverted_s4 <- workflow_data$state_manager$revertToState(prev_state_name)
          output$rollup_results <- shiny::renderText(paste("Reverted to previous state:", prev_state_name))
          logger::log_info(paste("Reverted precursor rollup to", prev_state_name))
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
    
    # Render rollup plot
    output$rollup_plot <- shiny::renderPlot({
      shiny::req(rollup_plot())
      grid::grid.draw(rollup_plot())
    })
  })
}

