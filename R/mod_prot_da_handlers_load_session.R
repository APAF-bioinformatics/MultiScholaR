#' Internal session loading handler for Proteomics DA module
#' @keywords internal
da_server_load_session_handler <- function(input, output, session, da_data, workflow_data, experiment_paths) {
  shiny::observeEvent(input$load_filtered_session, {
    cat("=== LOAD FILTERED SESSION BUTTON CLICKED (DA) ===\n")

    tryCatch(
      {
        # Check for latest session file
        source_dir <- experiment_paths$source_dir
        if (is.null(source_dir) || !dir.exists(source_dir)) {
          stop("Could not find the source directory to load session data.")
        }

        latest_session_file <- file.path(source_dir, "filtered_session_data_latest.rds")

        if (!file.exists(latest_session_file)) {
          stop("No exported session data found. Please export session data from the Normalization tab first.")
        }

        shiny::withProgress(message = "Loading filtered session data...", value = 0, {
          # Step 1: Load session data
          shiny::incProgress(0.3, detail = "Reading session file...")
          session_data <- readRDS(latest_session_file)

          cat("*** LOAD: Session data loaded successfully ***\n")
          cat(sprintf("*** LOAD: Export timestamp: %s ***\n", session_data$export_timestamp))
          cat(sprintf("*** LOAD: R6 current state: %s ***\n", session_data$r6_current_state_name))
          cat(sprintf("*** LOAD: Protein count: %d ***\n", session_data$final_protein_count))
          cat(sprintf("*** LOAD: Sample count: %d ***\n", session_data$final_sample_count))

          # Step 2: Restore R6 state manager completely
          shiny::incProgress(0.3, detail = "Restoring R6 state...")

          # Restore the complete R6 state structure
          workflow_data$state_manager$states <- session_data$r6_complete_states
          workflow_data$state_manager$state_history <- session_data$r6_state_history
          workflow_data$state_manager$current_state <- session_data$r6_current_state_name

          cat("*** LOAD DA: R6 state manager completely restored ***\n")
          cat(sprintf("*** LOAD DA: Restored %d states ***\n", length(session_data$r6_complete_states)))
          cat(sprintf("*** LOAD DA: State history: %s ***\n", paste(unlist(session_data$r6_state_history), collapse = " -> ")))
          cat(sprintf("*** LOAD DA: Current state set to: %s ***\n", session_data$r6_current_state_name))

          # Step 3: Restore workflow data
          shiny::incProgress(0.2, detail = "Restoring workflow data...")

          workflow_data$contrasts_tbl <- session_data$contrasts_tbl
          workflow_data$design_matrix <- session_data$design_matrix
          workflow_data$config_list <- session_data$config_list

          # Set global contrasts_tbl for DA analysis
          if (!is.null(session_data$contrasts_tbl)) {
            assign("contrasts_tbl", session_data$contrasts_tbl, envir = .GlobalEnv)
            cat("*** LOAD DA: Restored contrasts_tbl to global environment ***\n")
          }

          # Set global config_list
          if (!is.null(session_data$config_list)) {
            assign("config_list", session_data$config_list, envir = .GlobalEnv)
            cat("*** LOAD: Restored config_list to global environment ***\n")
          }

          # Step 4: Update DA module state
          shiny::incProgress(0.15, detail = "Updating DA module...")

          da_data$current_s4_object <- session_data$current_s4_object

          # Set contrasts available from loaded session
          if (!is.null(session_data$contrasts_tbl)) {
            if ("contrasts" %in% names(session_data$contrasts_tbl)) {
              da_data$contrasts_available <- session_data$contrasts_tbl$contrasts
            } else {
              da_data$contrasts_available <- session_data$contrasts_tbl[, 1]
            }
          }

          # Extract formula from S4 object
          if (!is.null(session_data$current_s4_object) && "deAnalysisParameters" %in% names(session_data$current_s4_object@args)) {
            if ("formula_string" %in% names(session_data$current_s4_object@args$deAnalysisParameters)) {
              formula_from_s4 <- session_data$current_s4_object@args$deAnalysisParameters$formula_string
              da_data$formula_from_s4 <- formula_from_s4

              # Update UI with formula from S4
              shiny::updateTextAreaInput(
                session,
                "formula_string",
                value = formula_from_s4
              )
            }
          }

          # [OK] NEW: Load UniProt annotations if available
          shiny::incProgress(0.05, detail = "Loading UniProt annotations...")

          cat("*** LOAD: Checking for UniProt annotations ***\n")
          tryCatch(
            {
              # Try to load uniprot_dat_cln.RDS from scripts directory
              scripts_uniprot_path <- file.path(source_dir, "uniprot_dat_cln.RDS")

              if (file.exists(scripts_uniprot_path)) {
                cat(sprintf("*** LOAD: Found uniprot_dat_cln.RDS at %s ***\n", scripts_uniprot_path))

                uniprot_dat_cln <- readRDS(scripts_uniprot_path)

                # Store in workflow data and global environment
                workflow_data$uniprot_dat_cln <- uniprot_dat_cln
                assign("uniprot_dat_cln", uniprot_dat_cln, envir = .GlobalEnv)

                cat(sprintf("*** LOAD: Successfully loaded %d UniProt annotations ***\n", nrow(uniprot_dat_cln)))
                logger::log_info(paste("Loaded UniProt annotations from", scripts_uniprot_path))

                # Show notification to user
                shiny::showNotification(
                  sprintf("UniProt annotations loaded: %d protein annotations available for enrichment analysis", nrow(uniprot_dat_cln)),
                  type = "message",
                  duration = 5
                )
              } else {
                cat(sprintf("*** LOAD: No uniprot_dat_cln.RDS found at %s ***\n", scripts_uniprot_path))
                logger::log_info("No UniProt annotations file found - enrichment analysis may have limited functionality")

                shiny::showNotification(
                  "No UniProt annotations found in session. Enrichment analysis may be limited.",
                  type = "warning",
                  duration = 5
                )
              }
            },
            error = function(e) {
              cat(sprintf("*** LOAD: Error loading UniProt annotations: %s ***\n", e$message))
              logger::log_warn(paste("Error loading UniProt annotations:", e$message))

              shiny::showNotification(
                paste("Warning: Could not load UniProt annotations:", e$message),
                type = "warning",
                duration = 8
              )
            }
          )

          # [OK] NEW: Restore FASTA metadata for complete audit trail
          cat("*** LOAD: Restoring FASTA metadata ***\n")
          if (!is.null(session_data$fasta_metadata)) {
            workflow_data$fasta_metadata <- session_data$fasta_metadata
            cat(sprintf(
              "*** LOAD: FASTA metadata restored (format: %s, sequences: %d) ***\n",
              session_data$fasta_metadata$fasta_format,
              session_data$fasta_metadata$num_sequences
            ))
          } else {
            cat("*** LOAD: No FASTA metadata in session data ***\n")
          }

          # [OK] NEW: Restore accession cleanup results
          cat("*** LOAD: Restoring accession cleanup results ***\n")
          if (!is.null(session_data$accession_cleanup_results)) {
            workflow_data$accession_cleanup_results <- session_data$accession_cleanup_results
            cat(sprintf(
              "*** LOAD: Accession cleanup results restored (applied: %s, method: %s) ***\n",
              session_data$accession_cleanup_results$cleanup_applied,
              session_data$accession_cleanup_results$aggregation_method
            ))
          } else {
            cat("*** LOAD: No accession cleanup results in session data ***\n")
          }

          # [OK] NEW: Restore complete RUV optimization results
          cat("*** LOAD: Restoring complete RUV optimization results ***\n")
          if (!is.null(session_data$ruv_optimization_result)) {
            workflow_data$ruv_optimization_result <- session_data$ruv_optimization_result
            cat(sprintf(
              "*** LOAD: RUV optimization results restored (k: %d, percentage: %.1f%%) ***\n",
              session_data$ruv_optimization_result$best_k,
              session_data$ruv_optimization_result$best_percentage
            ))
          } else {
            cat("*** LOAD: No complete RUV optimization results in session data ***\n")
          }

          # [OK] NEW: Restore QC parameters
          cat("*** LOAD: Restoring QC parameters ***\n")
          if (!is.null(session_data$qc_params)) {
            workflow_data$qc_params <- session_data$qc_params
            param_count <- length(unlist(session_data$qc_params, recursive = FALSE))
            cat(sprintf("*** LOAD: QC parameters restored (%d parameter groups) ***\n", param_count))
          } else {
            cat("*** LOAD: No QC parameters in session data ***\n")
          }

          # [OK] NEW: Restore mixed species FASTA analysis metadata for enrichment filtering
          cat("*** LOAD: Restoring mixed species analysis metadata ***\n")
          if (!is.null(session_data$mixed_species_analysis)) {
            workflow_data$mixed_species_analysis <- session_data$mixed_species_analysis

            if (isTRUE(session_data$mixed_species_analysis$enabled)) {
              cat(sprintf(
                "*** LOAD: Mixed species analysis restored (enabled: TRUE, organism: %s, taxon: %s) ***\n",
                session_data$mixed_species_analysis$selected_organism,
                session_data$mixed_species_analysis$selected_taxon_id
              ))

              shiny::showNotification(
                sprintf(
                  "Mixed species FASTA detected: %s. Enrichment filtering will be auto-enabled.",
                  session_data$mixed_species_analysis$selected_organism
                ),
                type = "message",
                duration = 5
              )
            } else {
              cat("*** LOAD: Mixed species analysis restored (enabled: FALSE) ***\n")
            }
          } else {
            cat("*** LOAD: No mixed species analysis metadata in session data ***\n")
          }

          # Update tab status - must replace entire list to trigger reactivity
          updated_status <- workflow_data$tab_status
          updated_status$normalization <- "complete"
          updated_status$differential_expression <- "pending"
          workflow_data$tab_status <- updated_status

          # Trigger state update for other modules
          workflow_data$state_update_trigger <- Sys.time()
        })

        # Create summary message
        summary_msg <- sprintf(
          "Session loaded successfully!\n\nData Summary:\n- Proteins: %d\n- Samples: %d\n- Contrasts: %d\n- State: %s\n- Export time: %s\n\nReady for differential abundance analysis.",
          session_data$final_protein_count,
          session_data$final_sample_count,
          ifelse(is.null(session_data$contrasts_tbl), 0, nrow(session_data$contrasts_tbl)),
          session_data$r6_current_state_name,
          format(session_data$export_timestamp, "%Y-%m-%d %H:%M:%S")
        )

        shiny::showNotification(
          summary_msg,
          type = "message",
          duration = 10
        )

        cat("=== LOAD FILTERED SESSION COMPLETED SUCCESSFULLY ===\n")
      },
      error = function(e) {
        # CRITICAL FIX: Use paste() for logger calls in error handlers to avoid interpolation bug
        error_msg <- paste("Error loading session:", e$message)
        cat(paste("***", error_msg, "\n"))
        logger::log_error(error_msg)

        shiny::showNotification(
          error_msg,
          type = "error",
          duration = 10
        )
      }
    )
  })
}
