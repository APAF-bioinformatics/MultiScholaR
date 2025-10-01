#' @title sample_filter_server
#'
#' @description Server module for applying the sample quality filter.
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
sample_filter_server <- function(id, workflow_data, omic_type, experiment_label) {
  shiny::moduleServer(id, function(input, output, session) {
    
    sample_plot <- reactiveVal(NULL)
    
    # Step 5: Sample Quality Filter (chunk 14)
    observeEvent(input$apply_sample_filter, {
      shiny::req(workflow_data$state_manager)
      
      shiny::showNotification("Applying sample quality filter...", id = "sample_working", duration = NULL)
      
      tryCatch({
        # Get current S4 object from the active state
        current_s4 <- workflow_data$state_manager$getState()
        shiny::req(current_s4)
        
        logger::log_info(sprintf("QC Step: Applying sample quality filter (min: %s)", input$min_peptides_per_sample))
        
        # ✅ FIXED: Use updateConfigParameter to sync S4 object AND global config_list
        current_s4 <- updateConfigParameter(
          theObject = current_s4,
          function_name = "filterMinNumPeptidesPerSample",
          parameter_name = "peptides_per_sample_cutoff",
          new_value = input$min_peptides_per_sample
        )
        
        # Apply S4 transformation (EXISTING S4 CODE - UNCHANGED)
        filtered_s4 <- filterMinNumPeptidesPerSample(theObject = current_s4)
        
        # Save new state in R6 manager
        workflow_data$state_manager$saveState(
          state_name = "sample_filtered",
          s4_data_object = filtered_s4,
          config_object = list(
            min_peptides_per_sample = input$min_peptides_per_sample
          ),
          description = "Applied minimum peptides per sample filter"
        )
        
        # Generate summary
        protein_count <- filtered_s4@peptide_data |>
          dplyr::distinct(Protein.Ids) |>
          nrow()
        run_count <- filtered_s4@peptide_data |>
          dplyr::distinct(Run) |>
          nrow()
        
        result_text <- paste(
          "Sample Quality Filter Applied Successfully\n",
          "=========================================\n",
          sprintf("Proteins remaining: %d\n", protein_count),
          sprintf("Samples remaining: %d\n", run_count),
          sprintf("Min peptides per sample: %d\n", input$min_peptides_per_sample),
          "State saved as: 'sample_filtered'\n"
        )
        
        output$sample_results <- shiny::renderText(result_text)
        
        # Update filtering visualization and capture plot
        plot_grid <- updateProteinFiltering(
          data = filtered_s4@peptide_data,
          step_name = "6_sample_filtered",
          omic_type = omic_type,
          experiment_label = experiment_label,
          return_grid = TRUE,
          overwrite = TRUE
        )
        sample_plot(plot_grid)
        
        logger::log_info("Sample quality filter applied successfully")
        shiny::removeNotification("sample_working")
        shiny::showNotification("Sample quality filter applied successfully", type = "message")
        
      }, error = function(e) {
        msg <- paste("Error applying sample quality filter:", e$message)
        logger::log_error(msg)
        shiny::showNotification(msg, type = "error", duration = 15)
        shiny::removeNotification("sample_working")
      })
    })
    
    # Revert Sample Filter
    observeEvent(input$revert_sample, {
      tryCatch({
        # Revert to the previous state in the history
        history <- workflow_data$state_manager$getHistory()
        if (length(history) > 1) {
          prev_state_name <- history[length(history) - 1]
          reverted_s4 <- workflow_data$state_manager$revertToState(prev_state_name)
          output$sample_results <- shiny::renderText(paste("Reverted to previous state:", prev_state_name))
          logger::log_info(paste("Reverted sample filter to", prev_state_name))
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

    # Render sample filter plot
    output$sample_plot <- shiny::renderPlot({
      shiny::req(sample_plot())
      grid::grid.draw(sample_plot())
    })
    
  })
}
