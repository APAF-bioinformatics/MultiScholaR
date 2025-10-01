#' @title imputation_server
#'
#' @description Server module for applying missing value imputation.
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
imputation_server <- function(id, workflow_data, omic_type, experiment_label) {
  shiny::moduleServer(id, function(input, output, session) {
    
    imputation_plot <- reactiveVal(NULL)
    
    # Step 7: Missing Value Imputation (chunk 16)
    observeEvent(input$apply_imputation, {
      shiny::req(workflow_data$state_manager)
      
      shiny::showNotification("Applying missing value imputation...", id = "imputation_working", duration = NULL)
      
      tryCatch({
        # Get current S4 object from the active state
        current_s4 <- workflow_data$state_manager$getState()
        shiny::req(current_s4)
        
        logger::log_info(sprintf("QC Step: Applying missing value imputation (proportion: %s)", input$proportion_missing_values))
        
        # ✅ FIXED: Use updateConfigParameter to sync S4 object AND global config_list
        current_s4 <- updateConfigParameter(
          theObject = current_s4,
          function_name = "peptideMissingValueImputation",
          parameter_name = "proportion_missing_values",
          new_value = input$proportion_missing_values
        )
        
        # Apply S4 transformation (EXISTING S4 CODE - UNCHANGED)
        imputed_s4 <- peptideMissingValueImputation(theObject = current_s4)
        
        # Save new state in R6 manager
        workflow_data$state_manager$saveState(
          state_name = "imputed",
          s4_data_object = imputed_s4,
          config_object = list(
            proportion_missing_values = input$proportion_missing_values
          ),
          description = "Applied missing value imputation using technical replicates"
        )
        
        # Generate summary
        protein_count <- imputed_s4@peptide_data |>
          dplyr::distinct(Protein.Ids) |>
          nrow()
        
        result_text <- paste(
          "Missing Value Imputation Applied Successfully\n",
          "============================================\n",
          sprintf("Proteins remaining: %d\n", protein_count),
          sprintf("Max proportion missing: %.1f\n", input$proportion_missing_values),
          "State saved as: 'imputed'\n",
          "\nReady for peptide-to-protein rollup step."
        )
        
        output$imputation_results <- shiny::renderText(result_text)
        
        # Update filtering visualization and capture plot
        plot_grid <- updateProteinFiltering(
          data = imputed_s4@peptide_data,
          step_name = "8_imputed",
          omic_type = omic_type,
          experiment_label = experiment_label,
          return_grid = TRUE,
          overwrite = TRUE
        )
        imputation_plot(plot_grid)
        
        logger::log_info("Missing value imputation applied successfully")
        shiny::removeNotification("imputation_working")
        shiny::showNotification("Missing value imputation applied successfully", type = "message")
        
      }, error = function(e) {
        msg <- paste("Error applying missing value imputation:", e$message)
        logger::log_error(msg)
        shiny::showNotification(msg, type = "error", duration = 15)
        shiny::removeNotification("imputation_working")
      })
    })
    
    # Revert Imputation
    observeEvent(input$revert_imputation, {
      tryCatch({
        # Revert to the previous state in the history
        history <- workflow_data$state_manager$getHistory()
        if (length(history) > 1) {
          prev_state_name <- history[length(history) - 1]
          reverted_s4 <- workflow_data$state_manager$revertToState(prev_state_name)
          output$imputation_results <- shiny::renderText(paste("Reverted to previous state:", prev_state_name))
          logger::log_info(paste("Reverted imputation to", prev_state_name))
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

    # Render imputation plot
    output$imputation_plot <- shiny::renderPlot({
      shiny::req(imputation_plot())
      grid::grid.draw(imputation_plot())
    })
    
  })
}
