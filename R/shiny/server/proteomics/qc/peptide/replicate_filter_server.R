#' @title replicate_filter_server
#'
#' @description Server module for applying the replicate filter.
#' Extracted from peptideQCApplet_server.R for modularity.
#'
#' @param input Shiny input object from parent module
#' @param output Shiny output object from parent module
#' @param session Shiny session object from parent module
#' @param workflow_data A reactive values object to store workflow data.
#' @param omic_type The omics type (e.g., "proteomics").
#' @param experiment_label The experiment label.
#'
#' @export
#' @import shiny
#' @importFrom logger log_info log_error
replicate_filter_server <- function(input, output, session, workflow_data, omic_type, experiment_label) {
  
  replicate_plot <- reactiveVal(NULL)
    
    # Step 6: Replicate Filter (chunk 15)
    observeEvent(input$apply_replicate_filter, {
      shiny::req(workflow_data$state_manager)
      
      shiny::showNotification("Applying replicate filter...", id = "replicate_working", duration = NULL)
      
      tryCatch({
        # Get current S4 object from the active state
        current_s4 <- workflow_data$state_manager$getState()
        shiny::req(current_s4)
        
        logger::log_info(sprintf("QC Step: Applying replicate filter (column: %s)", input$replicate_group_column))
        
        # Apply S4 transformation (EXISTING S4 CODE - UNCHANGED)
        # Note: This function takes the replicate column as a parameter
        filtered_s4 <- removePeptidesWithOnlyOneReplicate(
          current_s4,
          replicate_group_column = input$replicate_group_column
        )
        
        # Save new state in R6 manager
        workflow_data$state_manager$saveState(
          state_name = "replicate_filtered",
          s4_data_object = filtered_s4,
          config_object = list(
            replicate_group_column = input$replicate_group_column
          ),
          description = "Applied replicate filter (removed single-replicate peptides)"
        )
        
        # Generate summary
        protein_count <- filtered_s4@peptide_data |>
          dplyr::distinct(Protein.Ids) |>
          nrow()
        
        result_text <- paste(
          "Replicate Filter Applied Successfully\n",
          "====================================\n",
          sprintf("Proteins remaining: %d\n", protein_count),
          sprintf("Replicate group column: %s\n", input$replicate_group_column),
          "State saved as: 'replicate_filtered'\n"
        )
        
        output$replicate_results <- shiny::renderText(result_text)
        
        # Update filtering visualization and capture plot
        plot_grid <- updateProteinFiltering(
          data = filtered_s4@peptide_data,
          step_name = "7_replicate_filtered",
          omic_type = omic_type,
          experiment_label = experiment_label,
          return_grid = TRUE,
          overwrite = TRUE
        )
        replicate_plot(plot_grid)
        
        logger::log_info("Replicate filter applied successfully")
        shiny::removeNotification("replicate_working")
        shiny::showNotification("Replicate filter applied successfully", type = "message")
        
      }, error = function(e) {
        msg <- paste("Error applying replicate filter:", e$message)
        logger::log_error(msg)
        shiny::showNotification(msg, type = "error", duration = 15)
        shiny::removeNotification("replicate_working")
      })
    })
    
    # Revert Replicate Filter
    observeEvent(input$revert_replicate, {
      tryCatch({
        # Revert to the previous state in the history
        history <- workflow_data$state_manager$getHistory()
        if (length(history) > 1) {
          prev_state_name <- history[length(history) - 1]
          reverted_s4 <- workflow_data$state_manager$revertToState(prev_state_name)
          output$replicate_results <- shiny::renderText(paste("Reverted to previous state:", prev_state_name))
          logger::log_info(paste("Reverted replicate filter to", prev_state_name))
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

  # Render replicate filter plot
  output$replicate_plot <- shiny::renderPlot({
    shiny::req(replicate_plot())
    grid::grid.draw(replicate_plot())
  })
}
