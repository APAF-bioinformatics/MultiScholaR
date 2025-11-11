#' @title qvalue_filter_server
#'
#' @description Server module for applying Q-value and proteotypic peptide filters.
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
qvalue_filter_server <- function(input, output, session, workflow_data, omic_type, experiment_label) {
  
  qvalue_plot <- reactiveVal(NULL)
    
    # Step 1: Q-Value Filter (chunk 10)
    observeEvent(input$apply_qvalue_filter, {
      shiny::req(workflow_data$state_manager)
      
      shiny::showNotification("Applying Q-value filter...", id = "qvalue_working", duration = NULL)
      
      tryCatch({
        # Get current S4 object from state manager - NOW USES CURRENT STATE
        current_s4 <- workflow_data$state_manager$getState()
        shiny::req(current_s4)
        
        # ✅ DIAGNOSTIC: Check config and columns
        logger::log_info(sprintf("Q-value filter: S4 class = %s", class(current_s4)[1]))
        logger::log_info(sprintf("Q-value filter: S4 @peptide_data has %d rows, %d columns", 
                                 nrow(current_s4@peptide_data), ncol(current_s4@peptide_data)))
        
        # ✅ DIAGNOSTIC: Check input_matrix_column_ids from config
        if (!is.null(current_s4@args$srlQvalueProteotypicPeptideClean$input_matrix_column_ids)) {
          ids_vector <- current_s4@args$srlQvalueProteotypicPeptideClean$input_matrix_column_ids
          logger::log_info(sprintf("Q-value filter: input_matrix_column_ids length = %d", length(ids_vector)))
          logger::log_info(sprintf("Q-value filter: input_matrix_column_ids values: [%s]", 
                                   paste(sapply(ids_vector, function(x) paste0("'", x, "'")), collapse = ", ")))
          
          # Check for issues
          has_whitespace <- any(grepl("^\\s|\\s$", ids_vector))
          has_empty <- any(ids_vector == "")
          if (has_whitespace) logger::log_warn("Q-value filter: input_matrix_column_ids contains values with leading/trailing whitespace!")
          if (has_empty) logger::log_warn("Q-value filter: input_matrix_column_ids contains empty strings!")
        }
        
        logger::log_info(sprintf("QC Step: Applying Q-value filter with thresholds %s, %s", input$qvalue_threshold, input$global_qvalue_threshold))
        
        # ✅ FIXED: Use updateConfigParameter to sync S4 object AND global config_list
        current_s4 <- updateConfigParameter(
          theObject = current_s4,
          function_name = "srlQvalueProteotypicPeptideClean",
          parameter_name = "qvalue_threshold",
          new_value = input$qvalue_threshold
        )
        
        current_s4 <- updateConfigParameter(
          theObject = current_s4,
          function_name = "srlQvalueProteotypicPeptideClean",
          parameter_name = "global_qvalue_threshold",
          new_value = input$global_qvalue_threshold
        )
        
        current_s4 <- updateConfigParameter(
          theObject = current_s4,
          function_name = "srlQvalueProteotypicPeptideClean",
          parameter_name = "choose_only_proteotypic_peptide",
          new_value = as.numeric(input$proteotypic_only)
        )
        
        # Apply S4 transformation (EXISTING S4 CODE - UNCHANGED)
        filtered_s4 <- srlQvalueProteotypicPeptideClean(theObject = current_s4)
        
        # Track QC parameters in workflow_data
        if (is.null(workflow_data$qc_params)) {
          workflow_data$qc_params <- list()
        }
        if (is.null(workflow_data$qc_params$peptide_qc)) {
          workflow_data$qc_params$peptide_qc <- list()
        }
        
        workflow_data$qc_params$peptide_qc$qvalue_filter <- list(
          qvalue_threshold = input$qvalue_threshold,
          global_qvalue_threshold = input$global_qvalue_threshold,
          proteotypic_only = input$proteotypic_only,
          timestamp = Sys.time()
        )
        
        # Save new state in R6 manager
        workflow_data$state_manager$saveState(
          state_name = "qvalue_filtered",
          s4_data_object = filtered_s4,
          config_object = list(
            qvalue_threshold = input$qvalue_threshold,
            global_qvalue_threshold = input$global_qvalue_threshold,
            proteotypic_only = input$proteotypic_only
          ),
          description = "Applied Q-value and proteotypic peptide filter"
        )
        
        # Generate filtering summary
        protein_count <- filtered_s4@peptide_data |>
          dplyr::distinct(Protein.Ids) |>
          nrow()
        
        result_text <- paste(
          "Q-Value Filter Applied Successfully\n",
          "================================\n",
          sprintf("Proteins remaining: %d\n", protein_count),
          sprintf("Q-value threshold: %g\n", input$qvalue_threshold),
          sprintf("Global Q-value threshold: %g\n", input$global_qvalue_threshold),
          sprintf("Proteotypic only: %s\n", input$proteotypic_only),
          sprintf("State saved as: 'qvalue_filtered'\n")
        )
        
        output$qvalue_results <- shiny::renderText(result_text)
        
        # Update filtering visualization and capture plot
        plot_grid <- updateProteinFiltering(
          data = filtered_s4@peptide_data,
          step_name = "2_qval_filtered",
          omic_type = omic_type,
          experiment_label = experiment_label,
          return_grid = TRUE,
          overwrite = TRUE
        )
        qvalue_plot(plot_grid)
        
        logger::log_info("Q-value filter applied successfully")
        shiny::removeNotification("qvalue_working")
        shiny::showNotification("Q-value filter applied successfully", type = "message")
        
      }, error = function(e) {
        msg <- paste("Error applying Q-value filter:", e$message)
        logger::log_error(msg)
        shiny::showNotification(msg, type = "error", duration = 15)
        shiny::removeNotification("qvalue_working")
      })
    })
    
    # Revert Q-Value Filter
    observeEvent(input$revert_qvalue, {
      tryCatch({
        # Revert to raw data state - THIS LOGIC IS NOW TOO SIMPLE AND NEEDS REVIEW
        # For now, we revert to the state named 'raw_data_s4' if it exists.
        # A more robust undo/redo system could be implemented later.
        if ("raw_data_s4" %in% workflow_data$state_manager$getHistory()) {
          reverted_s4 <- workflow_data$state_manager$revertToState("raw_data_s4")
          output$qvalue_results <- shiny::renderText("Reverted to raw data state")
          logger::log_info("Reverted to raw data state")
          shiny::showNotification("Reverted to raw data state", type = "message")
        } else {
          stop("Cannot revert: 'raw_data_s4' state not found in history.")
        }
        
      }, error = function(e) {
        msg <- paste("Error reverting:", e$message)
        logger::log_error(msg)
        shiny::showNotification(msg, type = "error")
      })
    })
    
  # Render Q-value filter plot
  output$qvalue_plot <- shiny::renderPlot({
    shiny::req(qvalue_plot())
    grid::grid.draw(qvalue_plot())
  })
}
