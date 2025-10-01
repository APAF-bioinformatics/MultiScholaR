#' @title precursor_rollup_server
#'
#' @description Server module for applying precursor to peptide rollup.
#' Extracted from peptideQCApplet_server.R for modularity.
#'
#' @param id Module ID
#' @param workflow_data A reactive values object to store workflow data.
#' @param omic_type The omics type (e.g., "proteomics").
#' @param experiment_label The experiment label.
#'
#' @export
#' @import shiny
#' @importFrom logger log_info log_error
precursor_rollup_server <- function(id, workflow_data, omic_type, experiment_label) {
  shiny::moduleServer(id, function(input, output, session) {
    
    rollup_plot <- reactiveVal(NULL)
    
    # Step 2: Precursor Rollup (chunk 11)
    observeEvent(input$apply_rollup, {
      shiny::req(workflow_data$state_manager)
      
      shiny::showNotification("Applying precursor rollup...", id = "rollup_working", duration = NULL)
      
      tryCatch({
        # Get current S4 object from the active state
        current_s4 <- workflow_data$state_manager$getState()
        shiny::req(current_s4)
        
        logger::log_info("QC Step: Applying precursor rollup")
        
        # Apply S4 transformation (EXISTING S4 CODE - UNCHANGED)
        rolled_up_s4 <- rollUpPrecursorToPeptide(current_s4)
        
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
    observeEvent(input$revert_rollup, {
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
