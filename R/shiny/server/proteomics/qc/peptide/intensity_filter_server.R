#' @title intensity_filter_server
#'
#' @description Server module for applying peptide intensity filtering.
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
intensity_filter_server <- function(input, output, session, workflow_data, omic_type, experiment_label) {
  
  intensity_plot <- reactiveVal(NULL)
    
    # Step 3: Intensity Filter (chunk 12)
    observeEvent(input$apply_intensity_filter, {
      shiny::req(workflow_data$state_manager)
      
      shiny::showNotification("Applying intensity filter...", id = "intensity_working", duration = NULL)
      
      tryCatch({
        # Get current S4 object from the active state
        current_s4 <- workflow_data$state_manager$getState()
        shiny::req(current_s4)
        
        logger::log_info(sprintf("QC Step: Applying intensity filter with cutoff %s%%", input$intensity_cutoff_percentile))
        
        # ✅ FIXED: Use updateConfigParameter to sync S4 object AND global config_list
        current_s4 <- updateConfigParameter(
          theObject = current_s4,
          function_name = "peptideIntensityFiltering",
          parameter_name = "peptides_intensity_cutoff_percentile",
          new_value = input$intensity_cutoff_percentile
        )
        
        current_s4 <- updateConfigParameter(
          theObject = current_s4,
          function_name = "peptideIntensityFiltering",
          parameter_name = "peptides_proportion_of_samples_below_cutoff",
          new_value = input$proportion_samples_below_cutoff
        )
        
        # Apply S4 transformation (EXISTING S4 CODE - UNCHANGED)
        filtered_s4 <- peptideIntensityFiltering(theObject = current_s4)
        
        # Save new state in R6 manager
        workflow_data$state_manager$saveState(
          state_name = "intensity_filtered",
          s4_data_object = filtered_s4,
          config_object = list(
            intensity_cutoff_percentile = input$intensity_cutoff_percentile,
            proportion_samples_below_cutoff = input$proportion_samples_below_cutoff
          ),
          description = "Applied peptide intensity filter"
        )
        
        # Generate summary
        protein_count <- filtered_s4@peptide_data |>
          dplyr::distinct(Protein.Ids) |>
          nrow()
        
        result_text <- paste(
          "Intensity Filter Applied Successfully\n",
          "====================================\n",
          sprintf("Proteins remaining: %d\n", protein_count),
          sprintf("Intensity cutoff percentile: %.1f%%\n", input$intensity_cutoff_percentile),
          sprintf("Max proportion below cutoff: %.1f\n", input$proportion_samples_below_cutoff),
          "State saved as: 'intensity_filtered'\n"
        )
        
        output$intensity_results <- shiny::renderText(result_text)
        
        # Update filtering visualization and capture plot
        plot_grid <- updateProteinFiltering(
          data = filtered_s4@peptide_data,
          step_name = "4_intensity_filtered",
          omic_type = omic_type,
          experiment_label = experiment_label,
          return_grid = TRUE,
          overwrite = TRUE
        )
        intensity_plot(plot_grid)
        
        logger::log_info("Intensity filter applied successfully")
        shiny::removeNotification("intensity_working")
        shiny::showNotification("Intensity filter applied successfully", type = "message")
        
      }, error = function(e) {
        msg <- paste("Error applying intensity filter:", e$message)
        logger::log_error(msg)
        shiny::showNotification(msg, type = "error", duration = 15)
        shiny::removeNotification("intensity_working")
      })
    })
    
    # Revert Intensity Filter
    observeEvent(input$revert_intensity, {
      tryCatch({
        # Revert to the previous state in the history
        history <- workflow_data$state_manager$getHistory()
        if (length(history) > 1) {
          prev_state_name <- history[length(history) - 1]
          reverted_s4 <- workflow_data$state_manager$revertToState(prev_state_name)
          output$intensity_results <- shiny::renderText(paste("Reverted to previous state:", prev_state_name))
          logger::log_info(paste("Reverted intensity filter to", prev_state_name))
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

  # Render intensity filter plot
  output$intensity_plot <- shiny::renderPlot({
    shiny::req(intensity_plot())
    grid::grid.draw(intensity_plot())
  })
}
