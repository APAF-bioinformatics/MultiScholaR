#' @title proteinQCAppletServer
#'
#' @description Server component for protein-level quality control and filtering.
#' Extracted from qualityControlApplet.R to enable format-specific workflows.
#' This module handles all protein filtering server logic after peptide-to-protein rollup.
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
proteinQCAppletServer <- function(id, workflow_data, experiment_paths, omic_type, experiment_label) {
  shiny::moduleServer(id, function(input, output, session) {
    
    # Reactive values to store plots from each protein filtering step
    iq_rollup_plot <- reactiveVal(NULL)
    accession_cleanup_plot <- reactiveVal(NULL)
    protein_intensity_filter_plot <- reactiveVal(NULL)
    duplicate_removal_plot <- reactiveVal(NULL)
    protein_replicate_filter_plot <- reactiveVal(NULL)
    
    # == Protein Filtering Server Logic ====================================
    # Extracted from qualityControlApplet.R lines ~1550-2120
    
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
        
        # --- START DEBUGGING ---
        message("\n\n--- DEBUG: Inspecting S4 object AFTER Protein Creation ---")
        message(sprintf("S4 object class: %s", class(protein_obj)))

        if (methods::is(protein_obj, "ProteinQuantitativeData")) {
          message(sprintf("Dimensions of protein_obj@protein_quant_table: %d rows, %d cols",
                          nrow(protein_obj@protein_quant_table), ncol(protein_obj@protein_quant_table)))

          head_output <- capture.output(head(protein_obj@protein_quant_table))
          message("Head of protein_obj@protein_quant_table:")
          for (line in head_output) {
            message(line)
          }

          message("Head of protein_obj@design_matrix:")
          head_dm_output <- capture.output(head(protein_obj@design_matrix))
           for (line in head_dm_output) {
            message(line)
          }

        } else {
          message("DEBUG: The created S4 object is NOT of class ProteinQuantitativeData.")
        }
        message("--- END DEBUGGING ---\n\n")
        # --- END DEBUGGING ---
        
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
        
        # CRITICAL FIX: Update state trigger to notify reactive outputs (delimiter detection, etc.)
        cat("\n========== IQ ROLLUP: UPDATING STATE TRIGGER ==========\n")
        old_trigger <- workflow_data$state_update_trigger
        workflow_data$state_update_trigger <- Sys.time()
        cat(sprintf("   Old trigger: %s\n", old_trigger))
        cat(sprintf("   New trigger: %s\n", workflow_data$state_update_trigger))
        cat("========== STATE TRIGGER UPDATED ==========\n\n")
        
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
    
    # Step 2: Protein Accession Cleanup (chunk 19)
    
    # Helper function to detect delimiters in protein IDs
    detectProteinIdDelimiters <- function(protein_ids) {
      cat("=== detectProteinIdDelimiters called ===\n")
      cat(paste("Number of IDs to check:", length(protein_ids), "\n"))
      
      # Common delimiters to check
      common_delimiters <- c(";", ":", ",", "|", "/", "_")
      
      # Count occurrences of each delimiter
      delimiter_counts <- sapply(common_delimiters, function(delim) {
        sum(grepl(delim, protein_ids, fixed = TRUE))
      })
      
      cat("Delimiter counts:\n")
      print(delimiter_counts)
      
      # Return delimiters that appear in at least one protein ID
      found_delimiters <- common_delimiters[delimiter_counts > 0]
      cat(paste("Found delimiters:", paste(found_delimiters, collapse = ", "), "\n"))
      
      # Also count how many IDs contain each delimiter
      delimiter_info <- sapply(found_delimiters, function(delim) {
        count <- sum(grepl(delim, protein_ids, fixed = TRUE))
        paste0(delim, " (", count, " IDs)")
      })
      
      result <- list(
        delimiters = found_delimiters,
        info = delimiter_info
      )
      
      cat("Returning result:\n")
      print(result)
      cat("=== detectProteinIdDelimiters done ===\n")
      
      return(result)
    }
    
    # Reactive to detect delimiters in current data
    detected_delimiter_info <- shiny::reactive({
      cat("\n========== DELIMITER DETECTION REACTIVE TRIGGERED ==========\n")
      
      # Trigger on state changes
      workflow_data$state_update_trigger
      
      cat("State update trigger checked\n")
      cat(paste("workflow_data$state_manager is NULL:", is.null(workflow_data$state_manager), "\n"))
      
      # Only run if we have the protein S4 object
      if (!is.null(workflow_data$state_manager)) {
        cat("State manager exists, checking for states...\n")
        tryCatch({
          # Try to get the most recent protein state
          # Check multiple possible states in order of most recent
          current_s4 <- NULL
          possible_states <- c("protein_replicate_filtered", "protein_duplicate_removed", 
                              "protein_intensity_filtered", "protein_accession_cleaned", 
                              "protein_s4_created")
          
          for (state_name in possible_states) {
            if (state_name %in% workflow_data$state_manager$getHistory()) {
              current_s4 <- workflow_data$state_manager$getState(state_name)
              if (!is.null(current_s4)) {
                message(paste("DEBUG: Found protein S4 object in state:", state_name))
                break
              }
            }
          }
          
          if (!is.null(current_s4) && !is.null(current_s4@protein_quant_table)) {
            protein_id_col <- current_s4@protein_id_column
            message(paste("DEBUG: protein_id_col =", protein_id_col))
            
            if (protein_id_col %in% names(current_s4@protein_quant_table)) {
              protein_ids <- current_s4@protein_quant_table[[protein_id_col]]
              message(paste("DEBUG: Found", length(protein_ids), "protein IDs"))
              message(paste("DEBUG: First 3 IDs:", paste(head(protein_ids, 3), collapse = ", ")))
              
              delimiter_detection <- detectProteinIdDelimiters(protein_ids)
              message(paste("DEBUG: Detected delimiters:", paste(delimiter_detection$delimiters, collapse = ", ")))
              
              if (length(delimiter_detection$delimiters) > 0) {
                # Return both the display string and the most common delimiter
                return(list(
                  display = paste(delimiter_detection$info, collapse = ", "),
                  most_common = delimiter_detection$delimiters[1]  # First one is most common
                ))
              } else {
                return(list(
                  display = "No common delimiters detected (single accessions only)",
                  most_common = NULL
                ))
              }
            } else {
              message(paste("DEBUG: protein_id_col not found in columns:", paste(names(current_s4@protein_quant_table), collapse = ", ")))
            }
          } else {
            message("DEBUG: current_s4 is NULL or protein_quant_table is NULL")
          }
        }, error = function(e) {
          message(paste("DEBUG ERROR in delimiter detection:", e$message))
          return(list(
            display = paste("Unable to detect delimiters:", e$message),
            most_common = NULL
          ))
        })
      } else {
        message("DEBUG: workflow_data$state_manager is NULL")
      }
      
      return(list(
        display = "Data not loaded yet - complete IQ Rollup first",
        most_common = NULL
      ))
    })
    
    # Output detected delimiters
    output$detected_delimiters <- shiny::renderText({
      info <- detected_delimiter_info()
      if (is.list(info)) {
        info$display
      } else {
        info
      }
    })
    
    # Auto-update delimiter input when data is detected
    observe({
      info <- detected_delimiter_info()
      if (is.list(info) && !is.null(info$most_common)) {
        # Only update if the current value is the default ";"
        # This prevents overwriting user's manual changes
        if (input$delimiter == ";") {
          shiny::updateTextInput(session, "delimiter", value = info$most_common)
          logger::log_info(paste("Auto-detected and set delimiter to:", info$most_common))
        }
      }
    })
    
    # Display sample protein IDs
    output$sample_protein_ids <- shiny::renderText({
      # Trigger on state changes
      workflow_data$state_update_trigger
      
      if (!is.null(workflow_data$state_manager)) {
        tryCatch({
          # Try to get the most recent protein state (same logic as delimiter detection)
          current_s4 <- NULL
          possible_states <- c("protein_replicate_filtered", "protein_duplicate_removed", 
                              "protein_intensity_filtered", "protein_accession_cleaned", 
                              "protein_s4_created")
          
          for (state_name in possible_states) {
            if (state_name %in% workflow_data$state_manager$getHistory()) {
              current_s4 <- workflow_data$state_manager$getState(state_name)
              if (!is.null(current_s4)) {
                break
              }
            }
          }
          
          if (!is.null(current_s4) && !is.null(current_s4@protein_quant_table)) {
            protein_id_col <- current_s4@protein_id_column
            
            if (protein_id_col %in% names(current_s4@protein_quant_table)) {
              protein_ids <- current_s4@protein_quant_table[[protein_id_col]]
              
              # Show up to 10 example protein IDs
              # Prioritize showing ones with delimiters first
              has_delim <- grepl("[;:,|/_]", protein_ids)
              
              if (any(has_delim)) {
                # Show IDs with delimiters first
                with_delims <- protein_ids[has_delim]
                examples <- head(with_delims, 7)
                
                # Add a few without delimiters if available
                without_delims <- protein_ids[!has_delim]
                if (length(without_delims) > 0) {
                  examples <- c(examples, head(without_delims, 3))
                }
              } else {
                # No delimiters found, just show first 10
                examples <- head(protein_ids, 10)
              }
              
              return(paste(examples, collapse = "\n"))
            }
          }
        }, error = function(e) {
          return(paste("Unable to display sample IDs:", e$message))
        })
      }
      
      return("Data not loaded yet - complete IQ Rollup first")
    })
    
    observeEvent(input$apply_accession_cleanup, {
      shiny::req(workflow_data$state_manager)
      
      shiny::showNotification("Applying protein accession cleanup...", id = "accession_cleanup_working", duration = NULL)
      
      tryCatch({
        # Get current ProteinQuantitativeData S4 object
        current_s4 <- workflow_data$state_manager$getState("protein_s4_created")
        shiny::req(current_s4)
        
        logger::log_info("Protein Processing: Applying accession cleanup with delimiter: {input$delimiter}")
        
        # Check if aa_seq_tbl_final exists in workflow_data (from FASTA processing)
        if (exists("aa_seq_tbl_final", envir = .GlobalEnv)) {
          aa_seq_tbl_final <- get("aa_seq_tbl_final", envir = .GlobalEnv)
          aa_seq_tbl_final <- aa_seq_tbl_final |>
            dplyr::rename(uniprot_acc = database_id)
          
          # Apply S4 transformation (EXISTING S4 CODE - UNCHANGED)
          cleaned_s4 <- chooseBestProteinAccession(
            theObject = current_s4,
            delim = input$delimiter,
            seqinr_obj = aa_seq_tbl_final,
            seqinr_accession_column = "uniprot_acc",
            replace_zero_with_na = TRUE,
            aggregation_method = input$aggregation_method
          )
          
          cleanup_applied <- TRUE
        } else {
          # No FASTA data available, skip cleanup
          cleaned_s4 <- current_s4
          cleanup_applied <- FALSE
        }
        
        # Save new state
        workflow_data$state_manager$saveState(
          state_name = "protein_accession_cleaned",
          s4_data_object = cleaned_s4,
          config_object = list(
            delimiter = input$delimiter,
            aggregation_method = input$aggregation_method,
            cleanup_applied = cleanup_applied
          ),
          description = if (cleanup_applied) "Applied protein accession cleanup" else "Skipped accession cleanup (no FASTA data)"
        )
        
        # Generate summary
        protein_count <- cleaned_s4@protein_quant_table |>
          dplyr::distinct(Protein.Ids) |>
          nrow()
        
        result_text <- if (cleanup_applied) {
          paste(
            "Protein Accession Cleanup Applied Successfully\n",
            "==============================================\n",
            sprintf("Proteins remaining: %d\n", protein_count),
            sprintf("Delimiter: %s\n", input$delimiter),
            sprintf("Aggregation method: %s\n", input$aggregation_method),
            "State saved as: 'protein_accession_cleaned'\n"
          )
        } else {
          paste(
            "Protein Accession Cleanup Skipped\n",
            "=================================\n",
            sprintf("Proteins remaining: %d\n", protein_count),
            "Reason: No FASTA data available for cleanup\n",
            "State saved as: 'protein_accession_cleaned'\n"
          )
        }
        
        output$accession_cleanup_results <- shiny::renderText(result_text)
        
        # Update filtering visualization and capture plot
        plot_grid <- updateProteinFiltering(
          data = cleaned_s4@protein_quant_table,
          step_name = "10_protein_accession_cleaned",
          omic_type = omic_type,
          experiment_label = experiment_label,
          return_grid = TRUE,
          overwrite = TRUE
        )
        accession_cleanup_plot(plot_grid)
        
        logger::log_info("Protein accession cleanup completed")
        shiny::removeNotification("accession_cleanup_working")
        shiny::showNotification("Protein accession cleanup completed", type = "message")
        
      }, error = function(e) {
        msg <- paste("Error in protein accession cleanup:", e$message)
        logger::log_error(msg)
        shiny::showNotification(msg, type = "error", duration = 15)
        shiny::removeNotification("accession_cleanup_working")
      })
    })
    
    # Revert Accession Cleanup
    observeEvent(input$revert_accession_cleanup, {
      tryCatch({
        # Revert to protein S4 created state
        reverted_s4 <- workflow_data$state_manager$revertToState("protein_s4_created")
        output$accession_cleanup_results <- shiny::renderText("Reverted to protein S4 created state")
        
        logger::log_info("Reverted accession cleanup to protein_s4_created")
        shiny::showNotification("Reverted successfully", type = "message")
        
      }, error = function(e) {
        msg <- paste("Error reverting:", e$message)
        logger::log_error(msg)
        shiny::showNotification(msg, type = "error")
      })
    })
    
    # Step 3: Protein Intensity Filter (chunks 20+21 combined)
    observeEvent(input$apply_protein_intensity_filter, {
      shiny::req(workflow_data$state_manager)
      
      shiny::showNotification("Applying protein intensity filter...", id = "protein_intensity_filter_working", duration = NULL)
      
      tryCatch({
        # Get current ProteinQuantitativeData S4 object
        current_state <- workflow_data$state_manager$current_state
        current_s4 <- workflow_data$state_manager$getState(current_state)
        shiny::req(current_s4)
        
        logger::log_info("Protein Processing: Applying missing value parameters and intensity filter")
        
        # First: Update missing value parameters (chunk 20) using NEW function signature
        current_s4 <- updateMissingValueParameters(
          theObject = current_s4,
          min_reps_per_group = input$min_reps_per_group,
          min_groups = input$min_groups
        )
        
        # ✅ FIXED: Use ONLY the calculated percentages from updateMissingValueParameters
        # The replicate-based inputs (min_reps_per_group, min_groups) automatically calculate the percentages
        # No need to override with UI percentage inputs - they would conflict!
        
        # Only update the intensity cutoff percentile since it's not calculated by updateMissingValueParameters
        current_s4 <- updateConfigParameter(
          theObject = current_s4,
          function_name = "removeRowsWithMissingValuesPercent",
          parameter_name = "proteins_intensity_cutoff_percentile",
          new_value = input$proteins_intensity_cutoff_percentile
        )
        
        # Apply S4 transformation (EXISTING S4 CODE - UNCHANGED)
        filtered_s4 <- removeRowsWithMissingValuesPercent(current_s4)
        
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
        # Revert to previous state
        history <- workflow_data$state_manager$getHistory()
        prev_states <- c("protein_accession_cleaned", "protein_s4_created")
        prev_state <- intersect(rev(history), prev_states)[1]
        
        reverted_s4 <- workflow_data$state_manager$revertToState(prev_state)
        output$protein_intensity_filter_results <- shiny::renderText(paste("Reverted to", prev_state, "state"))
        
        logger::log_info("Reverted protein intensity filter to {prev_state}")
        shiny::showNotification("Reverted successfully", type = "message")
        
      }, error = function(e) {
        msg <- paste("Error reverting:", e$message)
        logger::log_error(msg)
        shiny::showNotification(msg, type = "error")
      })
    })
    
    # Step 4: Duplicate Protein Removal (chunk 22)
    observeEvent(input$apply_duplicate_removal, {
      shiny::req(workflow_data$state_manager)
      
      shiny::showNotification("Removing duplicate proteins...", id = "duplicate_removal_working", duration = NULL)
      
      tryCatch({
        # Get current ProteinQuantitativeData S4 object
        current_state <- workflow_data$state_manager$current_state
        current_s4 <- workflow_data$state_manager$getState(current_state)
        shiny::req(current_s4)
        
        logger::log_info("Protein Processing: Removing duplicate proteins using {input$duplicate_aggregation_method}")
        
        # Identify duplicates first
        duplicates <- current_s4@protein_quant_table |>
          dplyr::group_by(Protein.Ids) |>
          dplyr::filter(dplyr::n() > 1) |>
          dplyr::select(Protein.Ids) |>
          dplyr::distinct() |>
          dplyr::pull(Protein.Ids)
        
        # Apply duplicate removal
        current_s4@protein_quant_table <- current_s4@protein_quant_table |>
          dplyr::group_by(Protein.Ids) |>
          dplyr::summarise(
            dplyr::across(dplyr::matches("\\d+"), ~ get(input$duplicate_aggregation_method)(.x, na.rm = TRUE))
          ) |>
          dplyr::ungroup()
        
        # Save new state
        workflow_data$state_manager$saveState(
          state_name = "duplicates_removed",
          s4_data_object = current_s4,
          config_object = list(
            aggregation_method = input$duplicate_aggregation_method,
            duplicates_found = duplicates,
            num_duplicates = length(duplicates)
          ),
          description = "Removed duplicate proteins by aggregation"
        )
        
        # Generate summary
        protein_count <- current_s4@protein_quant_table |>
          dplyr::distinct(Protein.Ids) |>
          nrow()
        
        result_text <- paste(
          "Duplicate Protein Removal Completed Successfully\n",
          "===============================================\n",
          sprintf("Proteins remaining: %d\n", protein_count),
          sprintf("Duplicates found: %d\n", length(duplicates)),
          sprintf("Aggregation method: %s\n", input$duplicate_aggregation_method),
          "State saved as: 'duplicates_removed'\n"
        )
        
        output$duplicate_removal_results <- shiny::renderText(result_text)
        
        # Update filtering visualization and capture plot
        plot_grid <- updateProteinFiltering(
          data = current_s4@protein_quant_table,
          step_name = "12_duplicates_removed",
          omic_type = omic_type,
          experiment_label = experiment_label,
          return_grid = TRUE,
          overwrite = TRUE
        )
        duplicate_removal_plot(plot_grid)
        
        logger::log_info("Duplicate protein removal completed successfully")
        shiny::removeNotification("duplicate_removal_working")
        shiny::showNotification("Duplicate protein removal completed successfully", type = "message")
        
      }, error = function(e) {
        msg <- paste("Error removing duplicate proteins:", e$message)
        logger::log_error(msg)
        shiny::showNotification(msg, type = "error", duration = 15)
        shiny::removeNotification("duplicate_removal_working")
      })
    })
    
    # Revert Duplicate Removal
    observeEvent(input$revert_duplicate_removal, {
      tryCatch({
        # Revert to previous state
        history <- workflow_data$state_manager$getHistory()
        prev_states <- c("protein_intensity_filtered", "protein_accession_cleaned", "protein_s4_created")
        prev_state <- intersect(rev(history), prev_states)[1]
        
        reverted_s4 <- workflow_data$state_manager$revertToState(prev_state)
        output$duplicate_removal_results <- shiny::renderText(paste("Reverted to", prev_state, "state"))
        
        logger::log_info("Reverted duplicate removal to {prev_state}")
        shiny::showNotification("Reverted successfully", type = "message")
        
      }, error = function(e) {
        msg <- paste("Error reverting:", e$message)
        logger::log_error(msg)
        shiny::showNotification(msg, type = "error")
      })
    })
    
    # Step 5: Protein Replicate Filter (chunk 23)
    observeEvent(input$apply_protein_replicate_filter, {
      shiny::req(workflow_data$state_manager)
      
      shiny::showNotification("Applying protein replicate filter...", id = "protein_replicate_filter_working", duration = NULL)
      
      tryCatch({
        # Get current ProteinQuantitativeData S4 object
        current_state <- workflow_data$state_manager$current_state
        current_s4 <- workflow_data$state_manager$getState(current_state)
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
        # Revert to previous state
        history <- workflow_data$state_manager$getHistory()
        prev_states <- c("duplicates_removed", "protein_intensity_filtered", "protein_accession_cleaned", "protein_s4_created")
        prev_state <- intersect(rev(history), prev_states)[1]
        
        reverted_s4 <- workflow_data$state_manager$revertToState(prev_state)
        output$protein_replicate_filter_results <- shiny::renderText(paste("Reverted to", prev_state, "state"))
        
        logger::log_info("Reverted protein replicate filter to {prev_state}")
        shiny::showNotification("Reverted successfully", type = "message")
        
      }, error = function(e) {
        msg <- paste("Error reverting:", e$message)
        logger::log_error(msg)
        shiny::showNotification(msg, type = "error")
      })
    })
    
    # == Real-time Calculated Values Display ================================
    
    # Show calculated percentages based on replicate settings (for protein intensity filter)
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
    
    # Render IQ rollup plot
    output$iq_rollup_plot <- shiny::renderPlot({
      shiny::req(iq_rollup_plot())
      grid::grid.draw(iq_rollup_plot())
    })
    
    # Render accession cleanup plot
    output$accession_cleanup_plot <- shiny::renderPlot({
      shiny::req(accession_cleanup_plot())
      grid::grid.draw(accession_cleanup_plot())
    })
    
    # Render protein intensity filter plot
    output$protein_intensity_filter_plot <- shiny::renderPlot({
      shiny::req(protein_intensity_filter_plot())
      grid::grid.draw(protein_intensity_filter_plot())
    })
    
    # Render duplicate removal plot
    output$duplicate_removal_plot <- shiny::renderPlot({
      shiny::req(duplicate_removal_plot())
      grid::grid.draw(duplicate_removal_plot())
    })
    
    # Render protein replicate filter plot
    output$protein_replicate_filter_plot <- shiny::renderPlot({
      shiny::req(protein_replicate_filter_plot())
      grid::grid.draw(protein_replicate_filter_plot())
    })
    
  })
} 