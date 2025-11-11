#' @title protein_intensity_filter_server
#'
#' @description Server module for applying protein intensity and missing value filters.
#' Extracted from proteinQCApplet_server.R for modularity.
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
protein_intensity_filter_server <- function(input, output, session, workflow_data, omic_type, experiment_label) {
  
  protein_intensity_filter_plot <- reactiveVal(NULL)
    
    # Step 3: Protein Intensity Filter (chunks 20+21 combined)
    observeEvent(input$apply_protein_intensity_filter, {
      shiny::req(workflow_data$state_manager)
      
      shiny::showNotification("Applying protein intensity filter...", id = "protein_intensity_filter_working", duration = NULL)
      
      tryCatch({
        # Get current ProteinQuantitativeData S4 object from the active state
        current_s4 <- workflow_data$state_manager$getState()
        shiny::req(current_s4)
        
        logger::log_info("Protein Processing: Applying missing value parameters and intensity filter")
        
        # First: Update missing value parameters (chunk 20) using NEW function signature
        current_s4 <- updateMissingValueParameters(
          theObject = current_s4,
          min_reps_per_group = input$min_reps_per_group,
          min_groups = input$min_groups
        )
        
        # Only update the intensity cutoff percentile since it's not calculated by updateMissingValueParameters
        current_s4 <- updateConfigParameter(
          theObject = current_s4,
          function_name = "removeRowsWithMissingValuesPercent",
          parameter_name = "proteins_intensity_cutoff_percentile",
          new_value = input$proteins_intensity_cutoff_percentile
        )
        
        # Apply S4 transformation (EXISTING S4 CODE - UNCHANGED)
        filtered_s4 <- removeRowsWithMissingValuesPercent(current_s4)
        
        # Track QC parameters in workflow_data
        if (is.null(workflow_data$qc_params)) {
          workflow_data$qc_params <- list()
        }
        if (is.null(workflow_data$qc_params$protein_qc)) {
          workflow_data$qc_params$protein_qc <- list()
        }
        
        workflow_data$qc_params$protein_qc$intensity_filter <- list(
          min_reps_per_group = input$min_reps_per_group,
          min_groups = input$min_groups,
          groupwise_percentage_cutoff = current_s4@args$removeRowsWithMissingValuesPercent$groupwise_percentage_cutoff,
          max_groups_percentage_cutoff = current_s4@args$removeRowsWithMissingValuesPercent$max_groups_percentage_cutoff,
          proteins_intensity_cutoff_percentile = input$proteins_intensity_cutoff_percentile,
          timestamp = Sys.time()
        )
        
        # Save new state
        workflow_data$state_manager$saveState(
          state_name = "protein_intensity_filtered",
          s4_data_object = filtered_s4,
          config_object = list(
            min_reps_per_group = input$min_reps_per_group,
            min_groups = input$min_groups,
            groupwise_percentage_cutoff = current_s4@args$removeRowsWithMissingValuesPercent$groupwise_percentage_cutoff,
            max_groups_percentage_cutoff = current_s4@args$removeRowsWithMissingValuesPercent$max_groups_percentage_cutoff,
            proteins_intensity_cutoff_percentile = input$proteins_intensity_cutoff_percentile
          ),
          description = "Applied missing value parameters and protein intensity filter"
        )
        
        # Generate summary
        protein_count <- filtered_s4@protein_quant_table |>
          dplyr::distinct(Protein.Ids) |>
          nrow()
        
        result_text <- paste(
          "Protein Intensity Filter Applied Successfully\n",
          "============================================\n",
          sprintf("Proteins remaining: %d\n", protein_count),
          sprintf("Min replicates per group: %d\n", input$min_reps_per_group),
          sprintf("Min groups required: %d\n", input$min_groups),
          sprintf("Groupwise %% cutoff: %.3f%% (calculated)\n", current_s4@args$removeRowsWithMissingValuesPercent$groupwise_percentage_cutoff),
          sprintf("Max groups %% cutoff: %.3f%% (calculated)\n", current_s4@args$removeRowsWithMissingValuesPercent$max_groups_percentage_cutoff),
          sprintf("Intensity cutoff percentile: %.1f%%\n", input$proteins_intensity_cutoff_percentile),
          "State saved as: 'protein_intensity_filtered'\n"
        )
        
        output$protein_intensity_filter_results <- shiny::renderText(result_text)
        
        # Update filtering visualization and capture plot
        plot_grid <- updateProteinFiltering(
          data = filtered_s4@protein_quant_table,
          step_name = "11_protein_intensity_filtered",
          omic_type = omic_type,
          experiment_label = experiment_label,
          return_grid = TRUE,
          overwrite = TRUE
        )
        protein_intensity_filter_plot(plot_grid)
        
        logger::log_info("Protein intensity filter applied successfully")
        shiny::removeNotification("protein_intensity_filter_working")
        shiny::showNotification("Protein intensity filter applied successfully", type = "message")
        
      }, error = function(e) {
        msg <- paste("Error applying protein intensity filter:", e$message)
        logger::log_error(msg)
        shiny::showNotification(msg, type = "error", duration = 15)
        shiny::removeNotification("protein_intensity_filter_working")
      })
    })
    
    # Revert Protein Intensity Filter
    observeEvent(input$revert_protein_intensity_filter, {
      tryCatch({
        # Revert to the previous state in the history
        history <- workflow_data$state_manager$getHistory()
        if (length(history) > 1) {
          prev_state_name <- history[length(history) - 1]
          reverted_s4 <- workflow_data$state_manager$revertToState(prev_state_name)
          output$protein_intensity_filter_results <- shiny::renderText(paste("Reverted to previous state:", prev_state_name))
          logger::log_info(paste("Reverted protein intensity filter to", prev_state_name))
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
    
  # Render protein intensity filter plot
  output$protein_intensity_filter_plot <- shiny::renderPlot({
    shiny::req(protein_intensity_filter_plot())
    grid::grid.draw(protein_intensity_filter_plot())
  })
}
