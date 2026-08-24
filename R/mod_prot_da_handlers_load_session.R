#' Internal session loading handler for Proteomics DA module
#' @keywords internal
da_server_load_session_handler <- function(
  input,
  output,
  session,
  da_data,
  workflow_data,
  experiment_paths
) {
  shiny::observeEvent(input$load_filtered_session, {
    cat("=== LOAD FILTERED SESSION BUTTON CLICKED (DA) ===\n")

    tryCatch(
      {
        bundle <- resolveProtDaSessionBundle(experiment_paths)
        source_dir <- experiment_paths$source_dir
        if (!is.null(bundle$artifact_error)) {
          logger::log_warn(paste(
            "Portable artifact session was rejected; loaded legacy RDS:",
            conditionMessage(bundle$artifact_error)
          ))
        }
        ensureProtDaCompatibilitySession(bundle, source_dir)

        shiny::withProgress(
          message = "Loading filtered session data...",
          value = 0,
          {
            shiny::incProgress(0.3, detail = "Reading session file...")
            session_data <- bundle$session_data
            validateProtDaFilteredSession(
              session_data,
              source_dir = source_dir
            )

            cat("*** LOAD: Session data loaded successfully ***\n")
            cat(sprintf(
              "*** LOAD: Export timestamp: %s ***\n",
              session_data$export_timestamp
            ))
            cat(sprintf(
              "*** LOAD: R6 current state: %s ***\n",
              session_data$r6_current_state_name
            ))
            cat(sprintf(
              "*** LOAD: Protein count: %d ***\n",
              session_data$final_protein_count
            ))
            cat(sprintf(
              "*** LOAD: Sample count: %d ***\n",
              session_data$final_sample_count
            ))

            shiny::incProgress(0.3, detail = "Restoring R6 state...")
            applied <- applyProtDaSessionBundle(
              workflow_data,
              da_data,
              bundle
            )
            artifact_index <- restoreProtDaArtifactIndex(workflow_data)
            if (!is.null(artifact_index)) {
              da_data$da_results_list <- artifact_index
              da_data$analysis_complete <- TRUE
              workflow_data$da_analysis_results_list <- artifact_index
              restored_status <- workflow_data$tab_status
              restored_status$differential_expression <- "complete"
              restored_status$differential_abundance <- "complete"
              restored_status$enrichment_analysis <- "pending"
              workflow_data$tab_status <- restored_status
              restored_choices <- vapply(
                artifact_index$contrasts,
                `[[`,
                character(1),
                "friendly_name"
              )
              da_data$contrasts_available <- restored_choices
              shiny::updateSelectInput(
                session,
                "volcano_contrast",
                choices = stats::setNames(restored_choices, restored_choices)
              )
              shiny::updateSelectInput(
                session,
                "heatmap_contrast",
                choices = stats::setNames(restored_choices, restored_choices)
              )
              shiny::updateSelectInput(
                session,
                "table_contrast",
                choices = stats::setNames(restored_choices, restored_choices)
              )
            }
            state_snapshot <- applied$state_snapshot
            cat("*** LOAD DA: R6 state manager completely restored ***\n")
            cat(sprintf(
              "*** LOAD DA: Restored %d states ***\n",
              length(state_snapshot$r6_complete_states)
            ))
            cat(sprintf(
              "*** LOAD DA: State history: %s ***\n",
              paste(
                unlist(state_snapshot$r6_state_history),
                collapse = " -> "
              )
            ))
            cat(sprintf(
              "*** LOAD DA: Current state set to: %s ***\n",
              state_snapshot$r6_current_state_name
            ))

            shiny::incProgress(0.2, detail = "Restoring workflow data...")
            if (!is.null(session_data$contrasts_tbl)) {
              cat("*** LOAD DA: Restored contrasts_tbl to global environment ***\n")
            }
            if (!is.null(session_data$config_list)) {
              cat("*** LOAD: Restored config_list to global environment ***\n")
            }

            shiny::incProgress(0.15, detail = "Updating DA module...")
            if (!is.null(applied$formula)) {
              shiny::updateTextAreaInput(
                session,
                "formula_string",
                value = applied$formula
              )
            }

            shiny::incProgress(
              0.05,
              detail = "Loading UniProt annotations..."
            )
            cat("*** LOAD: Checking for UniProt annotations ***\n")
            if (!is.null(applied$annotations)) {
              annotation_count <- nrow(applied$annotations)
              cat(sprintf(
                "*** LOAD: Successfully loaded %d UniProt annotations ***\n",
                annotation_count
              ))
              logger::log_info(paste(
                "Loaded UniProt annotations from",
                bundle$session_path
              ))
              shiny::showNotification(
                sprintf(
                  paste(
                    "UniProt annotations loaded: %d protein annotations",
                    "available for enrichment analysis"
                  ),
                  annotation_count
                ),
                type = "message",
                duration = 5
              )
            } else {
              cat("*** LOAD: No UniProt annotations found ***\n")
              logger::log_info(paste(
                "No UniProt annotations file found - enrichment analysis",
                "may have limited functionality"
              ))
              shiny::showNotification(
                paste(
                  "No UniProt annotations found in session. Enrichment",
                  "analysis may be limited."
                ),
                type = "warning",
                duration = 5
              )
            }

            if (!is.null(session_data$fasta_metadata)) {
              cat(sprintf(
                paste0(
                  "*** LOAD: FASTA metadata restored (format: %s, ",
                  "sequences: %d) ***\n"
                ),
                session_data$fasta_metadata$fasta_format,
                session_data$fasta_metadata$num_sequences
              ))
            }
            if (!is.null(session_data$accession_cleanup_results)) {
              cat(sprintf(
                paste0(
                  "*** LOAD: Accession cleanup results restored ",
                  "(applied: %s, method: %s) ***\n"
                ),
                session_data$accession_cleanup_results$cleanup_applied,
                session_data$accession_cleanup_results$aggregation_method
              ))
            }
            if (!is.null(session_data$ruv_optimization_result)) {
              cat(sprintf(
                paste0(
                  "*** LOAD: RUV optimization results restored ",
                  "(k: %d, percentage: %.1f%%) ***\n"
                ),
                session_data$ruv_optimization_result$best_k,
                session_data$ruv_optimization_result$best_percentage
              ))
            }
            if (!is.null(session_data$qc_params)) {
              param_count <- length(unlist(
                session_data$qc_params,
                recursive = FALSE
              ))
              cat(sprintf(
                "*** LOAD: QC parameters restored (%d parameter groups) ***\n",
                param_count
              ))
            }
            if (!is.null(session_data$mixed_species_analysis)) {
              enabled <- isTRUE(session_data$mixed_species_analysis$enabled)
              cat(sprintf(
                "*** LOAD: Mixed species analysis restored (enabled: %s) ***\n",
                enabled
              ))
              if (enabled) {
                shiny::showNotification(
                  sprintf(
                    paste(
                      "Mixed species FASTA detected: %s. Enrichment",
                      "filtering will be auto-enabled."
                    ),
                    session_data$mixed_species_analysis$selected_organism
                  ),
                  type = "message",
                  duration = 5
                )
              }
            }
          }
        )

        summary_msg <- sprintf(
          paste0(
            "Session loaded successfully!\n\nData Summary:\n",
            "- Proteins: %d\n- Samples: %d\n- Contrasts: %d\n",
            "- State: %s\n- Export time: %s\n\n",
            "Ready for differential abundance analysis."
          ),
          session_data$final_protein_count,
          session_data$final_sample_count,
          ifelse(
            is.null(session_data$contrasts_tbl),
            0,
            nrow(session_data$contrasts_tbl)
          ),
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
