#' @title peptideQCAppletServer
#'
#' @description Server component for peptide-level quality control and filtering.
#' Extracted from qualityControlApplet.R to enable format-specific workflows.
#' This module handles all peptide filtering server logic for DIA workflows.
#'
#' @param id Module ID
#' @param workflow_data A reactive values object to store workflow data
#' @param experiment_paths A list of paths for the current experiment
#' @param omic_type The omics type (e.g., "proteomics")
#' @param experiment_label The experiment label
#' @export
#' @import shiny
#' @importFrom logger log_info log_error
#' @importFrom grid grid.draw
peptideQCAppletServer <- function(id, workflow_data, experiment_paths, omic_type, experiment_label) {
  shiny::moduleServer(id, function(input, output, session) {
    
    # Reactive values to store plots from each peptide filtering step
    qvalue_plot <- reactiveVal(NULL)
    rollup_plot <- reactiveVal(NULL)
    intensity_plot <- reactiveVal(NULL)
    protein_peptide_plot <- reactiveVal(NULL)
    sample_plot <- reactiveVal(NULL)
    replicate_plot <- reactiveVal(NULL)
    imputation_plot <- reactiveVal(NULL)
    
    # == Peptide Filtering Server Logic =====================================
    # Extracted from qualityControlApplet.R lines ~686-1540
    
    # Step 1: Q-Value Filter (chunk 10)
    observeEvent(input$apply_qvalue_filter, {
      shiny::req(workflow_data$state_manager)
      
      shiny::showNotification("Applying Q-value filter...", id = "qvalue_working", duration = NULL)
      
      tryCatch({
        # Get current S4 object from state manager
        current_s4 <- workflow_data$state_manager$getState("raw_data_s4")
        shiny::req(current_s4)
        
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
        # Revert to raw data state
        reverted_s4 <- workflow_data$state_manager$revertToState("raw_data_s4")
        
        output$qvalue_results <- shiny::renderText("Reverted to raw data state")
        
        logger::log_info("Reverted to raw data state")
        shiny::showNotification("Reverted to raw data state", type = "message")
        
      }, error = function(e) {
        msg <- paste("Error reverting:", e$message)
        logger::log_error(msg)
        shiny::showNotification(msg, type = "error")
      })
    })
    
    # Step 2: Precursor Rollup (chunk 11)
    observeEvent(input$apply_rollup, {
      shiny::req(workflow_data$state_manager)
      
      shiny::showNotification("Applying precursor rollup...", id = "rollup_working", duration = NULL)
      
      tryCatch({
        # Get current S4 object (should be qvalue_filtered or raw_data_s4)
        current_state <- workflow_data$state_manager$current_state
        current_s4 <- workflow_data$state_manager$getState(current_state)
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
        # Revert to previous state (qvalue_filtered or raw_data_s4)
        history <- workflow_data$state_manager$getHistory()
        if ("qvalue_filtered" %in% history) {
          reverted_s4 <- workflow_data$state_manager$revertToState("qvalue_filtered")
          output$rollup_results <- shiny::renderText("Reverted to Q-value filtered state")
        } else {
          reverted_s4 <- workflow_data$state_manager$revertToState("raw_data_s4")
          output$rollup_results <- shiny::renderText("Reverted to raw data state")
        }
        
        logger::log_info("Reverted precursor rollup")
        shiny::showNotification("Reverted successfully", type = "message")
        
      }, error = function(e) {
        msg <- paste("Error reverting:", e$message)
        logger::log_error(msg)
        shiny::showNotification(msg, type = "error")
      })
    })
    
    # Step 3: Intensity Filter (chunk 12)
    observeEvent(input$apply_intensity_filter, {
      shiny::req(workflow_data$state_manager)
      
      shiny::showNotification("Applying intensity filter...", id = "intensity_working", duration = NULL)
      
      tryCatch({
        # Get current S4 object
        current_state <- workflow_data$state_manager$current_state
        current_s4 <- workflow_data$state_manager$getState(current_state)
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
        # Revert to previous state
        history <- workflow_data$state_manager$getHistory()
        prev_states <- c("precursor_rollup", "qvalue_filtered", "raw_data_s4")
        prev_state <- intersect(rev(history), prev_states)[1]
        
        reverted_s4 <- workflow_data$state_manager$revertToState(prev_state)
        output$intensity_results <- shiny::renderText(paste("Reverted to", prev_state, "state"))
        
        logger::log_info(sprintf("Reverted intensity filter to %s", prev_state))
        shiny::showNotification("Reverted successfully", type = "message")
        
      }, error = function(e) {
        msg <- paste("Error reverting:", e$message)
        logger::log_error(msg)
        shiny::showNotification(msg, type = "error")
      })
    })
    
    # Step 4: Protein Peptide Count Filter (chunk 13)
    observeEvent(input$apply_protein_peptide_filter, {
      shiny::req(workflow_data$state_manager)
      
      shiny::showNotification("Applying protein peptide count filter...", id = "protein_peptide_working", duration = NULL)
      
      tryCatch({
        # Get current S4 object
        current_state <- workflow_data$state_manager$current_state
        current_s4 <- workflow_data$state_manager$getState(current_state)
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
        # Revert to previous state
        history <- workflow_data$state_manager$getHistory()
        prev_states <- c("intensity_filtered", "precursor_rollup", "qvalue_filtered", "raw_data_s4")
        prev_state <- intersect(rev(history), prev_states)[1]
        
        reverted_s4 <- workflow_data$state_manager$revertToState(prev_state)
        output$protein_peptide_results <- shiny::renderText(paste("Reverted to", prev_state, "state"))
        
        logger::log_info(sprintf("Reverted protein peptide filter to %s", prev_state))
        shiny::showNotification("Reverted successfully", type = "message")
        
      }, error = function(e) {
        msg <- paste("Error reverting:", e$message)
        logger::log_error(msg)
        shiny::showNotification(msg, type = "error")
      })
    })
    
    # Step 5: Sample Quality Filter (chunk 14)
    observeEvent(input$apply_sample_filter, {
      shiny::req(workflow_data$state_manager)
      
      shiny::showNotification("Applying sample quality filter...", id = "sample_working", duration = NULL)
      
      tryCatch({
        # Get current S4 object
        current_state <- workflow_data$state_manager$current_state
        current_s4 <- workflow_data$state_manager$getState(current_state)
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
        # Revert to previous state
        history <- workflow_data$state_manager$getHistory()
        prev_states <- c("protein_peptide_filtered", "intensity_filtered", "precursor_rollup", "qvalue_filtered", "raw_data_s4")
        prev_state <- intersect(rev(history), prev_states)[1]
        
        reverted_s4 <- workflow_data$state_manager$revertToState(prev_state)
        output$sample_results <- shiny::renderText(paste("Reverted to", prev_state, "state"))
        
        logger::log_info(sprintf("Reverted sample filter to %s", prev_state))
        shiny::showNotification("Reverted successfully", type = "message")
        
      }, error = function(e) {
        msg <- paste("Error reverting:", e$message)
        logger::log_error(msg)
        shiny::showNotification(msg, type = "error")
      })
    })
    
    # Step 6: Replicate Filter (chunk 15)
    observeEvent(input$apply_replicate_filter, {
      shiny::req(workflow_data$state_manager)
      
      shiny::showNotification("Applying replicate filter...", id = "replicate_working", duration = NULL)
      
      tryCatch({
        # Get current S4 object
        current_state <- workflow_data$state_manager$current_state
        current_s4 <- workflow_data$state_manager$getState(current_state)
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
        # Revert to previous state
        history <- workflow_data$state_manager$getHistory()
        prev_states <- c("sample_filtered", "protein_peptide_filtered", "intensity_filtered", "precursor_rollup", "qvalue_filtered", "raw_data_s4")
        prev_state <- intersect(rev(history), prev_states)[1]
        
        reverted_s4 <- workflow_data$state_manager$revertToState(prev_state)
        output$replicate_results <- shiny::renderText(paste("Reverted to", prev_state, "state"))
        
        logger::log_info(sprintf("Reverted replicate filter to %s", prev_state))
        shiny::showNotification("Reverted successfully", type = "message")
        
      }, error = function(e) {
        msg <- paste("Error reverting:", e$message)
        logger::log_error(msg)
        shiny::showNotification(msg, type = "error")
      })
    })
    
    # Step 7: Missing Value Imputation (chunk 16)
    observeEvent(input$apply_imputation, {
      shiny::req(workflow_data$state_manager)
      
      shiny::showNotification("Applying missing value imputation...", id = "imputation_working", duration = NULL)
      
      tryCatch({
        # Get current S4 object
        current_state <- workflow_data$state_manager$current_state
        current_s4 <- workflow_data$state_manager$getState(current_state)
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
        # Revert to previous state
        history <- workflow_data$state_manager$getHistory()
        prev_states <- c("replicate_filtered", "sample_filtered", "protein_peptide_filtered", "intensity_filtered", "precursor_rollup", "qvalue_filtered", "raw_data_s4")
        prev_state <- intersect(rev(history), prev_states)[1]
        
        reverted_s4 <- workflow_data$state_manager$revertToState(prev_state)
        output$imputation_results <- shiny::renderText(paste("Reverted to", prev_state, "state"))
        
        logger::log_info(sprintf("Reverted imputation to %s", prev_state))
        shiny::showNotification("Reverted successfully", type = "message")
        
      }, error = function(e) {
        msg <- paste("Error reverting:", e$message)
        logger::log_error(msg)
        shiny::showNotification(msg, type = "error")
      })
    })
    
    # == Real-time Calculated Values Display ================================
    
    # Show calculated percentages based on replicate settings
    output$calculated_groupwise_percent <- shiny::renderText({
      shiny::req(input$min_reps_per_group, input$min_groups)
      
      # Get design matrix if available
      current_state <- workflow_data$state_manager$current_state
      if (!is.null(current_state)) {
        current_s4 <- workflow_data$state_manager$getState(current_state)
        if (!is.null(current_s4) && !is.null(current_s4@design_matrix)) {
          design_matrix <- current_s4@design_matrix
          
          # Calculate the same way as updateMissingValueParameters (duplicated for UI display)
          reps_per_group_tbl <- design_matrix |>
            dplyr::group_by(group) |>
            dplyr::summarise(n_reps = n()) |>
            dplyr::ungroup()
          
          group_thresholds <- reps_per_group_tbl |>
            dplyr::mutate(
              adjusted_min_reps = pmin(n_reps, input$min_reps_per_group),
              max_missing = n_reps - adjusted_min_reps,
              missing_percent = round((max_missing / n_reps) * 100, 3)
            )
          
          groupwise_cutoff <- max(group_thresholds$missing_percent)
          return(sprintf("Groupwise Percentage Cutoff: %.3f%% (calculated)", groupwise_cutoff))
        }
      }
      
      return("Groupwise Percentage Cutoff: [Design matrix needed]")
    })
    
    output$calculated_max_groups_percent <- shiny::renderText({
      shiny::req(input$min_reps_per_group, input$min_groups)
      
      # Get design matrix if available
      current_state <- workflow_data$state_manager$current_state
      if (!is.null(current_state)) {
        current_s4 <- workflow_data$state_manager$getState(current_state)
        if (!is.null(current_s4) && !is.null(current_s4@design_matrix)) {
          design_matrix <- current_s4@design_matrix
          
          # Calculate the same way as updateMissingValueParameters (duplicated for UI display)
          total_groups <- design_matrix |> dplyr::distinct(group) |> nrow()
          max_failing_groups <- total_groups - input$min_groups
          max_groups_cutoff <- round((max_failing_groups / total_groups) * 100, 3)
          
          return(sprintf("Max Groups Percentage Cutoff: %.3f%% (calculated)", max_groups_cutoff))
        }
      }
      
      return("Max Groups Percentage Cutoff: [Design matrix needed]")
    })

    # == Plot Rendering Functions ==============================================
    
    # Render Q-value filter plot
    output$qvalue_plot <- shiny::renderPlot({
      shiny::req(qvalue_plot())
      grid::grid.draw(qvalue_plot())
    })
    
    # Render rollup plot
    output$rollup_plot <- shiny::renderPlot({
      shiny::req(rollup_plot())
      grid::grid.draw(rollup_plot())
    })
    
    # Render intensity filter plot
    output$intensity_plot <- shiny::renderPlot({
      shiny::req(intensity_plot())
      grid::grid.draw(intensity_plot())
    })
    
    # Render protein peptide filter plot
    output$protein_peptide_plot <- shiny::renderPlot({
      shiny::req(protein_peptide_plot())
      grid::grid.draw(protein_peptide_plot())
    })
    
    # Render sample filter plot
    output$sample_plot <- shiny::renderPlot({
      shiny::req(sample_plot())
      grid::grid.draw(sample_plot())
    })
    
    # Render replicate filter plot
    output$replicate_plot <- shiny::renderPlot({
      shiny::req(replicate_plot())
      grid::grid.draw(replicate_plot())
    })
    
    # Render imputation plot
    output$imputation_plot <- shiny::renderPlot({
      shiny::req(imputation_plot())
      grid::grid.draw(imputation_plot())
    })
    
  })
} 