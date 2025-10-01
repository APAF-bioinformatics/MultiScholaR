#' @title protein_rollup_server
#'
#' @description Server module for performing peptide-to-protein rollup using the IQ tool.
#' This is specific to LFQ/DIA workflows. Extracted from proteinQCApplet_server.R for modularity.
#'
#' @param id Module ID
#' @param workflow_data A reactive values object to store workflow data.
#' @param experiment_paths A list of paths for the current experiment.
#' @param omic_type The omics type (e.g., "proteomics").
#' @param experiment_label The experiment label.
#'
#' @export
#' @import shiny
#' @importFrom logger log_info log_error
protein_rollup_server <- function(id, workflow_data, experiment_paths, omic_type, experiment_label) {
  shiny::moduleServer(id, function(input, output, session) {
    
    iq_rollup_plot <- reactiveVal(NULL)
    
    # Step 1: IQ Protein Rollup (chunk 17)
    observeEvent(input$apply_iq_rollup, {
      shiny::req(workflow_data$state_manager)
      
      shiny::showNotification("Running IQ protein rollup & creating S4 object...", id = "iq_rollup_working", duration = NULL)
      
      tryCatch({
        # Get the final peptide S4 object (should be 'imputed' state)
        current_state <- workflow_data$state_manager$current_state
        peptide_s4 <- workflow_data$state_manager$getState(current_state)
        shiny::req(peptide_s4)
        
        logger::log_info("Protein Processing: Starting IQ rollup from peptide state")
        
        # Save peptide data to TSV file for IQ input
        peptide_values_imputed_file <- file.path(
          experiment_paths$peptide_qc_dir,
          "peptide_values_imputed.tsv"
        )
        
        # Prepare data for IQ (ensure correct column names and format)
        peptide_data_for_iq <- peptide_s4@peptide_data |>
          dplyr::mutate(
            Q.Value = 0.0009,
            PG.Q.Value = 0.009,
            Peptide.Imputed = ifelse(is.na(Peptide.Imputed), 0, Peptide.Imputed)
          )
        
        vroom::vroom_write(peptide_data_for_iq, peptide_values_imputed_file)
        
        # Run IQ processing
        iq_output_file <- file.path(experiment_paths$protein_qc_dir, "iq_output_file.txt")
        
        iq::process_long_format(
          peptide_values_imputed_file,
          output_filename = iq_output_file,
          sample_id = "Run",
          primary_id = "Protein.Ids",
          secondary_id = "Stripped.Sequence",
          intensity_col = "Peptide.Imputed",
          filter_double_less = c("Q.Value" = "0.01", "PG.Q.Value" = "0.01"),
          normalization = "none"  # Critical - no normalization at this stage
        )
        
        # Wait for IQ output file to be available
        max_wait <- 30  # Maximum 30 seconds
        wait_count <- 0
        while (!file.exists(iq_output_file) && wait_count < max_wait) {
          Sys.sleep(1)
          wait_count <- wait_count + 1
        }
        
        if (!file.exists(iq_output_file)) {
          stop("IQ output file not created within timeout period")
        }
        
        # Read IQ output
        protein_log2_quant <- vroom::vroom(iq_output_file)
        
        logger::log_info("Protein Processing: Creating ProteinQuantitativeData S4 object")
        
        # Create ProteinQuantitativeData S4 object
        protein_obj <- ProteinQuantitativeData(
          protein_quant_table = protein_log2_quant,
          protein_id_column = "Protein.Ids",
          protein_id_table = protein_log2_quant |> dplyr::distinct(Protein.Ids),
          design_matrix = peptide_s4@design_matrix,
          sample_id = "Run",
          group_id = "group",
          technical_replicate_id = "replicates",
          args = peptide_s4@args
        )
        
        # Save S4 object as new state (combined rollup + S4 creation)
        workflow_data$state_manager$saveState(
          state_name = "protein_s4_created",
          s4_data_object = protein_obj,
          config_object = list(
            iq_output_file = iq_output_file,
            peptide_input_file = peptide_values_imputed_file,
            s4_class = "ProteinQuantitativeData",
            protein_id_column = "Protein.Ids"
          ),
          description = "IQ protein rollup completed and ProteinQuantitativeData S4 object created"
        )
        
        # Generate summary
        protein_count <- protein_obj@protein_quant_table |>
          dplyr::distinct(Protein.Ids) |>
          nrow()
        
        result_text <- paste(
          "IQ Protein Rollup & S4 Object Creation Completed Successfully\n",
          "============================================================\n",
          sprintf("Proteins quantified: %d\n", protein_count),
          sprintf("Samples: %d\n", ncol(protein_obj@protein_quant_table) - 1),
          sprintf("Algorithm: MaxLFQ (via IQ tool)\n"),
          sprintf("S4 Class: %s\n", class(protein_obj)[1]),
          sprintf("Design matrix: %s\n", paste(colnames(protein_obj@design_matrix), collapse = ", ")),
          sprintf("Output file: %s\n", basename(iq_output_file)),
          "State saved as: 'protein_s4_created'\n",
          "\nReady for protein accession cleanup."
        )
        
        output$iq_rollup_results <- shiny::renderText(result_text)
        
        # Update filtering visualization and capture plot
        plot_grid <- updateProteinFiltering(
          data = protein_obj@protein_quant_table,
          step_name = "9_protein_s4_created",
          omic_type = omic_type,
          experiment_label = experiment_label,
          return_grid = TRUE,
          overwrite = TRUE
        )
        iq_rollup_plot(plot_grid)
        
        logger::log_info("IQ protein rollup and S4 object creation completed successfully")
        shiny::removeNotification("iq_rollup_working")
        shiny::showNotification("IQ protein rollup & S4 object creation completed successfully", type = "message")
        
      }, error = function(e) {
        msg <- paste("Error in IQ protein rollup & S4 creation:", e$message)
        logger::log_error(msg)
        shiny::showNotification(msg, type = "error", duration = 15)
        shiny::removeNotification("iq_rollup_working")
      })
    })
    
    # Revert IQ Rollup
    observeEvent(input$revert_iq_rollup, {
      tryCatch({
        # Revert to final peptide state (should be 'imputed')
        history <- workflow_data$state_manager$getHistory()
        peptide_states <- c("imputed", "replicate_filtered", "sample_filtered", "protein_peptide_filtered", 
                           "intensity_filtered", "precursor_rollup", "qvalue_filtered", "raw_data_s4")
        prev_state <- intersect(rev(history), peptide_states)[1]
        
        reverted_s4 <- workflow_data$state_manager$revertToState(prev_state)
        output$iq_rollup_results <- shiny::renderText(paste("Reverted to", prev_state, "state"))
        
        logger::log_info("Reverted IQ rollup to {prev_state}")
        shiny::showNotification("Reverted successfully", type = "message")
        
      }, error = function(e) {
        msg <- paste("Error reverting:", e$message)
        logger::log_error(msg)
        shiny::showNotification(msg, type = "error")
      })
    })

    # Render IQ rollup plot
    output$iq_rollup_plot <- shiny::renderPlot({
      shiny::req(iq_rollup_plot())
      grid::grid.draw(iq_rollup_plot())
    })
    
  })
}
