#' @title protein_replicate_filter_server
#'
#' @description Server module for applying the protein replicate filter.
#' Extracted from proteinQCApplet_server.R for modularity.
#'
#' @param input Shiny input object from parent module
#' @param output Shiny output object from parent module
#' @param session Shiny session object from parent module
#' @param workflow_data A reactive values object to store workflow data.
#' @param experiment_paths A list of paths for the current experiment.
#' @param omic_type The omics type (e.g., "proteomics").
#' @param experiment_label The experiment label.
#'
#' @export
#' @import shiny
#' @importFrom logger log_info log_error
protein_replicate_filter_server <- function(input, output, session, workflow_data, experiment_paths, omic_type, experiment_label) {
  
  protein_replicate_filter_plot <- reactiveVal(NULL)
    
    # Step 5: Protein Replicate Filter (chunk 23)
    observeEvent(input$apply_protein_replicate_filter, {
      shiny::req(workflow_data$state_manager)
      
      shiny::showNotification("Applying protein replicate filter...", id = "protein_replicate_filter_working", duration = NULL)
      
      tryCatch({
        # Get current ProteinQuantitativeData S4 object from the active state
        current_s4 <- workflow_data$state_manager$getState()
        shiny::req(current_s4)
        
        logger::log_info("Protein Processing: Applying protein replicate filter with {input$parallel_cores} cores")
        
        # Set up parallel processing
        core_utilisation <- new_cluster(input$parallel_cores)
        
        # Apply S4 transformation (EXISTING S4 CODE - UNCHANGED)
        filtered_s4 <- removeProteinsWithOnlyOneReplicate(
          current_s4,
          core_utilisation,
          grouping_variable = input$protein_grouping_variable
        )
        
        # Save filtered data to file (as in original workflow)
        output_file <- file.path(experiment_paths$protein_qc_dir, "remove_proteins_with_only_one_rep.tsv")
        vroom::vroom_write(filtered_s4@protein_quant_table, output_file)
        
        # Save new state
        workflow_data$state_manager$saveState(
          state_name = "protein_replicate_filtered",
          s4_data_object = filtered_s4,
          config_object = list(
            grouping_variable = input$protein_grouping_variable,
            parallel_cores = input$parallel_cores,
            output_file = output_file
          ),
          description = "Applied protein replicate filter (removed single-replicate proteins)"
        )
        
        # Generate summary
        protein_count <- filtered_s4@protein_quant_table |>
          dplyr::distinct(Protein.Ids) |>
          nrow()
        
        result_text <- paste(
          "Protein Replicate Filter Applied Successfully\n",
          "============================================\n",
          sprintf("Proteins remaining: %d\n", protein_count),
          sprintf("Grouping variable: %s\n", input$protein_grouping_variable),
          sprintf("Parallel cores used: %d\n", input$parallel_cores),
          sprintf("Output file: %s\n", basename(output_file)),
          "State saved as: 'protein_replicate_filtered'\n",
          "\nProtein filtering pipeline complete!"
        )
        
        output$protein_replicate_filter_results <- shiny::renderText(result_text)
        
        # Update filtering visualization and capture plot
        plot_grid <- updateProteinFiltering(
          data = filtered_s4@protein_quant_table,
          step_name = "13_protein_replicate_filtered",
          omic_type = omic_type,
          experiment_label = experiment_label,
          return_grid = TRUE,
          overwrite = TRUE
        )
        protein_replicate_filter_plot(plot_grid)
        
        logger::log_info("Protein replicate filter applied successfully")
        shiny::removeNotification("protein_replicate_filter_working")
        shiny::showNotification("Protein replicate filter applied successfully. Pipeline complete!", type = "message")
        
      }, error = function(e) {
        msg <- paste("Error applying protein replicate filter:", e$message)
        logger::log_error(msg)
        shiny::showNotification(msg, type = "error", duration = 15)
        shiny::removeNotification("protein_replicate_filter_working")
      })
    })
    
    # Revert Protein Replicate Filter
    observeEvent(input$revert_protein_replicate_filter, {
      tryCatch({
        # Revert to the previous state in the history
        history <- workflow_data$state_manager$getHistory()
        if (length(history) > 1) {
          prev_state_name <- history[length(history) - 1]
          reverted_s4 <- workflow_data$state_manager$revertToState(prev_state_name)
          output$protein_replicate_filter_results <- shiny::renderText(paste("Reverted to previous state:", prev_state_name))
          logger::log_info(paste("Reverted protein replicate filter to", prev_state_name))
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
    
  # Render protein replicate filter plot
  output$protein_replicate_filter_plot <- shiny::renderPlot({
    shiny::req(protein_replicate_filter_plot())
    grid::grid.draw(protein_replicate_filter_plot())
  })
}
