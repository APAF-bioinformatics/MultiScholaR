#' @title protein_peptide_filter_server
#'
#' @description Server module for applying the protein peptide count filter.
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
protein_peptide_filter_server <- function(id, workflow_data, omic_type, experiment_label) {
  shiny::moduleServer(id, function(input, output, session) {
    
    protein_peptide_plot <- reactiveVal(NULL)
    
    # Step 4: Protein Peptide Count Filter (chunk 13)
    observeEvent(input$apply_protein_peptide_filter, {
      shiny::req(workflow_data$state_manager)
      
      shiny::showNotification("Applying protein peptide count filter...", id = "protein_peptide_working", duration = NULL)
      
      tryCatch({
        # Get current S4 object from the active state
        current_s4 <- workflow_data$state_manager$getState()
        shiny::req(current_s4)
        
        logger::log_info(sprintf("QC Step: Applying protein peptide count filter (min: %s)", input$min_peptides_per_protein))
        
        # ✅ FIXED: Use updateConfigParameter to sync S4 object AND global config_list
        # Use config.ini parameter names, not function parameter names
        current_s4 <- updateConfigParameter(
          theObject = current_s4,
          function_name = "filterMinNumPeptidesPerProtein",
          parameter_name = "peptides_per_protein_cutoff",
          new_value = input$min_peptides_per_protein
        )
        
        current_s4 <- updateConfigParameter(
          theObject = current_s4,
          function_name = "filterMinNumPeptidesPerProtein",
          parameter_name = "peptidoforms_per_protein_cutoff",
          new_value = input$min_peptidoforms_per_protein
        )
        
        # Apply S4 transformation (EXISTING S4 CODE - UNCHANGED)
        filtered_s4 <- filterMinNumPeptidesPerProtein(theObject = current_s4)
        
        # Save new state in R6 manager
        workflow_data$state_manager$saveState(
          state_name = "protein_peptide_filtered",
          s4_data_object = filtered_s4,
          config_object = list(
            min_peptides_per_protein = input$min_peptides_per_protein,
            min_peptidoforms_per_protein = input$min_peptidoforms_per_protein
          ),
          description = "Applied minimum peptides per protein filter"
        )
        
        # Generate summary
        protein_count <- filtered_s4@peptide_data |>
          dplyr::distinct(Protein.Ids) |>
          nrow()
        
        result_text <- paste(
          "Protein Peptide Count Filter Applied Successfully\n",
          "===============================================\n",
          sprintf("Proteins remaining: %d\n", protein_count),
          sprintf("Min peptides per protein: %d\n", input$min_peptides_per_protein),
          sprintf("Min peptidoforms per protein: %d\n", input$min_peptidoforms_per_protein),
          "State saved as: 'protein_peptide_filtered'\n"
        )
        
        output$protein_peptide_results <- shiny::renderText(result_text)
        
        # Update filtering visualization and capture plot
        plot_grid <- updateProteinFiltering(
          data = filtered_s4@peptide_data,
          step_name = "5_protein_peptide_filtered",
          omic_type = omic_type,
          experiment_label = experiment_label,
          return_grid = TRUE,
          overwrite = TRUE
        )
        protein_peptide_plot(plot_grid)
        
        logger::log_info("Protein peptide count filter applied successfully")
        shiny::removeNotification("protein_peptide_working")
        shiny::showNotification("Protein peptide count filter applied successfully", type = "message")
        
      }, error = function(e) {
        msg <- paste("Error applying protein peptide count filter:", e$message)
        logger::log_error(msg)
        shiny::showNotification(msg, type = "error", duration = 15)
        shiny::removeNotification("protein_peptide_working")
      })
    })
    
    # Revert Protein Peptide Filter
    observeEvent(input$revert_protein_peptide, {
      tryCatch({
        # Revert to the previous state in the history
        history <- workflow_data$state_manager$getHistory()
        if (length(history) > 1) {
          prev_state_name <- history[length(history) - 1]
          reverted_s4 <- workflow_data$state_manager$revertToState(prev_state_name)
          output$protein_peptide_results <- shiny::renderText(paste("Reverted to previous state:", prev_state_name))
          logger::log_info(paste("Reverted protein peptide filter to", prev_state_name))
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

    # Render protein peptide filter plot
    output$protein_peptide_plot <- shiny::renderPlot({
      shiny::req(protein_peptide_plot())
      grid::grid.draw(protein_peptide_plot())
    })
    
  })
}
