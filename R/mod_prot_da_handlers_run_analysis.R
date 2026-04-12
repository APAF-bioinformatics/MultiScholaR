# Handler 3: DA Analysis Execution
da_server_run_analysis_handler <- function(input, output, session, ns, da_data, workflow_data, experiment_paths) {
  shiny::observeEvent(input$run_da_analysis, {
    cat("=== STARTING DIFFERENTIAL EXPRESSION ANALYSIS ===\n")

    shiny::req(da_data$current_s4_object, da_data$contrasts_available)

    shiny::showNotification("Running differential abundance analysis...", id = "da_working", duration = NULL)

    tryCatch(
      {
        shiny::withProgress(message = "Running DE analysis...", value = 0, {
          # Step 1: Prepare contrasts table for analysis
          shiny::incProgress(0.2, detail = "Preparing contrasts for analysis...")

          # Get contrasts from global environment or create basic ones
          contrasts_tbl <- NULL
          if (exists("contrasts_tbl", envir = .GlobalEnv)) {
            contrasts_tbl <- get("contrasts_tbl", envir = .GlobalEnv)
            cat("   DA ANALYSIS Step: Using existing contrasts_tbl from global environment\n")
            cat("   DA ANALYSIS Step: contrasts_tbl structure:\n")
            str(contrasts_tbl)
          } else {
            # Create contrasts table that matches the original format
            # The original expects a data frame with contrasts in the first column
            cat("   DA ANALYSIS Step: Creating contrasts_tbl from da_data$contrasts_available\n")
            cat("   DA ANALYSIS Step: Available contrasts:\n")
            print(da_data$contrasts_available)

            contrasts_tbl <- data.frame(
              contrasts = da_data$contrasts_available,
              stringsAsFactors = FALSE
            )
            cat("   DA ANALYSIS Step: Created contrasts_tbl:\n")
            str(contrasts_tbl)
            print(contrasts_tbl)
          }

          # CRITICAL FIX: Use the correct column for contrast strings
          # The downstream functions expect "comparison=expression" format (from full_format column)
          # NOT just the raw contrast expression (from contrasts column)
          if ("full_format" %in% names(contrasts_tbl)) {
            contrast_strings_to_use <- contrasts_tbl$full_format # Use full_format column
            cat("   DA ANALYSIS Step: Using full_format column for contrast strings\n")
          } else {
            cat("   DA ANALYSIS Step: No full_format column found, auto-generating from raw contrasts\n")
            # Auto-generate full_format column from raw contrasts
            raw_contrasts <- contrasts_tbl[, 1]

            # Generate friendly names and full format
            full_format_strings <- sapply(raw_contrasts, function(contrast_string) {
              # Remove "group" prefixes if present for friendly name
              clean_string <- gsub("^group", "", contrast_string)
              clean_string <- gsub("-group", "-", clean_string)

              # Create friendly name by replacing - with _vs_
              friendly_name <- gsub("-", "_vs_", clean_string)

              # Create full format: friendly_name=original_contrast_string
              paste0(friendly_name, "=", contrast_string)
            })

            contrast_strings_to_use <- full_format_strings
            cat("   DA ANALYSIS Step: Auto-generated full_format strings:\n")
            print(contrast_strings_to_use)
          }
          cat("   DA ANALYSIS Step: contrast_strings_to_use (what goes to runTestsContrasts):\n")
          print(contrast_strings_to_use)
          cat("   DA ANALYSIS Step: Length of contrast_strings_to_use:\n")
          print(length(contrast_strings_to_use))

          # Step 2: Run actual differential abundance analysis
          shiny::incProgress(0.4, detail = "Running limma analysis...")

          # CRITICAL FIX: Follow the original RMarkdown approach - run DA analysis for ONE contrast at a time
          cat("   DA ANALYSIS Step: Running DA analysis for each contrast separately (like original workflow)\n")

          # Use original contrast names as keys to match UI dropdown expectations
          # Don't transform them - keep them exactly as they appear in contrasts_tbl
          contrast_names <- contrasts_tbl$contrasts

          cat(sprintf("   DA ANALYSIS Step: Generated contrast names: %s\n", paste(contrast_names, collapse = ", ")))

          # Run DA analysis for each contrast separately using purrr::map (like original)
          da_results_list <- seq_len(nrow(contrasts_tbl)) |>
            purrr::set_names(contrast_names) |>
            purrr::map(\(contrast_idx) {
              cat(sprintf("   DA ANALYSIS Step: Processing contrast %d: %s\n", contrast_idx, contrasts_tbl$contrasts[contrast_idx]))

              # Get single contrast (one row) like original workflow
              single_contrast_tbl <- contrasts_tbl |> dplyr::slice(contrast_idx)

              # CRITICAL FIX: Ensure the single contrast table uses the correct format
              # Create a modified version that has the full_format in the first column
              if ("full_format" %in% names(single_contrast_tbl)) {
                # Create a temporary table with full_format as the first column for runTestsContrasts
                single_contrast_for_analysis <- data.frame(
                  contrasts = single_contrast_tbl$full_format,
                  stringsAsFactors = FALSE
                )
                cat("   DA ANALYSIS Step: Using full_format for individual contrast analysis\n")
              } else {
                cat("   DA ANALYSIS Step: No full_format available, auto-generating from raw contrast\n")
                # Auto-generate full_format for this single contrast
                raw_contrast <- single_contrast_tbl$contrasts[1]

                # Remove "group" prefixes if present for friendly name
                clean_string <- gsub("^group", "", raw_contrast)
                clean_string <- gsub("-group", "-", clean_string)

                # Create friendly name by replacing - with _vs_
                friendly_name <- gsub("-", "_vs_", clean_string)

                # Create full format: friendly_name=original_contrast_string
                full_format_string <- paste0(friendly_name, "=", raw_contrast)

                single_contrast_for_analysis <- data.frame(
                  contrasts = full_format_string,
                  stringsAsFactors = FALSE
                )
                cat(sprintf("   DA ANALYSIS Step: Auto-generated full format: %s\n", full_format_string))
              }
              cat(sprintf("   DA ANALYSIS Step: Single contrast table for index %d:\n", contrast_idx))
              cat("   DA ANALYSIS Step: Structure of single_contrast_for_analysis:\n")
              str(single_contrast_for_analysis)
              cat("   DA ANALYSIS Step: Content of single_contrast_for_analysis:\n")
              print(single_contrast_for_analysis)
              cat(sprintf("   DA ANALYSIS Step: Contrast string being passed: '%s'\n", single_contrast_for_analysis$contrasts[1]))

              # Use the modular DA analysis function for this single contrast
              cat("   DA ANALYSIS Step: About to call differentialAbundanceAnalysis...\n")
              cat(sprintf("   DA ANALYSIS Step: S4 object class: %s\n", class(da_data$current_s4_object)))
              cat(sprintf("   DA ANALYSIS Step: Formula: %s\n", input$formula_string))

              # CRITICAL: Ensure S4 object has the required parameters in @args
              # The differentialAbundanceAnalysis function expects parameters in @args$differentialAbundanceAnalysis
              if (is.null(da_data$current_s4_object@args$differentialAbundanceAnalysis)) {
                cat("   DA ANALYSIS Step: Adding differentialAbundanceAnalysis parameters to S4 @args\n")
                da_data$current_s4_object@args$differentialAbundanceAnalysis <- list(
                  contrasts_tbl = single_contrast_for_analysis,
                  formula_string = input$formula_string,
                  da_q_val_thresh = input$da_q_val_thresh,
                  treat_lfc_cutoff = input$treat_lfc_cutoff,
                  eBayes_trend = TRUE,
                  eBayes_robust = TRUE,
                  args_group_pattern = "(\\d+)",
                  args_row_id = da_data$current_s4_object@protein_id_column,
                  qvalue_column = "fdr_qvalue",
                  raw_pvalue_column = "raw_pvalue"
                )
              }

              # [OK] NEW: Store UI parameters in @args for session summary
              cat("   DA ANALYSIS Step: Storing UI parameters in S4 @args\n")
              if (is.null(da_data$current_s4_object@args$daAnalysisUI)) {
                da_data$current_s4_object@args$daAnalysisUI <- list()
              }

              da_data$current_s4_object@args$daAnalysisUI <- list(
                q_value_threshold = input$da_q_val_thresh,
                log_fold_change_cutoff = input$treat_lfc_cutoff,
                formula_string = input$formula_string,
                timestamp = Sys.time()
              )

              # [OK] NEW: Also store UI parameters in workflow_data for sessionSummary
              workflow_data$da_ui_params <- list(
                q_value_threshold = input$da_q_val_thresh,
                log_fold_change_cutoff = input$treat_lfc_cutoff,
                treat_enabled = input$treat_lfc_cutoff > 0,
                formula_string = input$formula_string,
                timestamp = Sys.time()
              )
              cat("   DA ANALYSIS Step: Stored UI parameters in workflow_data for sessionSummary\n")

              # [OK] NEW: Update R6 state manager with UI parameters
              cat("   DA ANALYSIS Step: Updating R6 state with DA UI parameters\n")
              tryCatch(
                {
                  # Find the current data state and update it
                  current_data_states <- c("correlation_filtered", "normalized", "ruv_corrected", "protein_replicate_filtered")
                  available_states <- workflow_data$state_manager$getHistory()
                  current_data_state <- purrr::detect(current_data_states, ~ .x %in% available_states)

                  if (!is.null(current_data_state)) {
                    # CRITICAL FIX: Use the correct 'saveState' method instead of non-existent 'updateState'
                    # Retrieve existing config and description to preserve them
                    existing_state_info <- workflow_data$state_manager$states[[current_data_state]]

                    workflow_data$state_manager$saveState(
                      state_name = current_data_state,
                      s4_data_object = da_data$current_s4_object,
                      config_object = existing_state_info$config,
                      description = existing_state_info$description
                    )
                    cat(sprintf("   DA ANALYSIS Step: Updated state '%s' with DA UI parameters\n", current_data_state))
                  }
                },
                error = function(e) {
                  # CRITICAL FIX: Use paste() for logger calls in error handlers
                  cat(paste("   DA ANALYSIS Step: Warning - could not update state with UI parameters:", e$message, "\n"))
                }
              )

              # Also add for the helper function
              if (is.null(da_data$current_s4_object@args$differentialAbundanceAnalysisHelper)) {
                cat("   DA ANALYSIS Step: Adding differentialAbundanceAnalysisHelper parameters to S4 @args\n")
                da_data$current_s4_object@args$differentialAbundanceAnalysisHelper <- list(
                  contrasts_tbl = single_contrast_for_analysis,
                  formula_string = input$formula_string,
                  da_q_val_thresh = input$da_q_val_thresh,
                  treat_lfc_cutoff = input$treat_lfc_cutoff,
                  eBayes_trend = TRUE,
                  eBayes_robust = TRUE,
                  args_group_pattern = "(\\d+)",
                  args_row_id = da_data$current_s4_object@protein_id_column,
                  qvalue_column = "fdr_qvalue",
                  raw_pvalue_column = "raw_pvalue"
                )
              }

              result <- tryCatch(
                {
                  # DEBUG: Check if function exists
                  cat("   DA ANALYSIS Step: Checking if differentialAbundanceAnalysis function exists...\n")
                  if (!exists("differentialAbundanceAnalysis")) {
                    cat("   DA ANALYSIS Step: ERROR - differentialAbundanceAnalysis function not found!\n")
                    cat("   DA ANALYSIS Step: Attempting to use MultiScholaR namespace...\n")

                    # Try with explicit namespace
                    MultiScholaR::differentialAbundanceAnalysis(
                      theObject = da_data$current_s4_object,
                      contrasts_tbl = single_contrast_for_analysis,
                      formula_string = input$formula_string,
                      da_q_val_thresh = input$da_q_val_thresh,
                      treat_lfc_cutoff = input$treat_lfc_cutoff,
                      qvalue_column = "fdr_qvalue",
                      raw_pvalue_column = "raw_pvalue"
                    )
                  } else {
                    cat("   DA ANALYSIS Step: differentialAbundanceAnalysis function found\n")

                    # Verify S4 method dispatch
                    cat(sprintf("   DA ANALYSIS Step: S4 object class: %s\n", class(da_data$current_s4_object)))
                    cat(sprintf("   DA ANALYSIS Step: Available methods for differentialAbundanceAnalysis:\n"))
                    print(methods::showMethods("differentialAbundanceAnalysis", printTo = FALSE))

                    differentialAbundanceAnalysis(
                      theObject = da_data$current_s4_object,
                      contrasts_tbl = single_contrast_for_analysis,
                      formula_string = input$formula_string,
                      da_q_val_thresh = input$da_q_val_thresh,
                      treat_lfc_cutoff = input$treat_lfc_cutoff,
                      qvalue_column = "fdr_qvalue",
                      raw_pvalue_column = "raw_pvalue"
                    )
                  }
                },
                error = function(e) {
                  cat(paste("*** ERROR in differentialAbundanceAnalysis call:", e$message, "\n"))
                  cat(paste("*** ERROR object: \n"))
                  print(str(e))
                  cat("*** Attempting detailed error diagnosis:\n")

                  # Check parameter types
                  cat(sprintf("   - da_data$current_s4_object class: %s\n", class(da_data$current_s4_object)))
                  cat(sprintf("   - single_contrast_tbl class: %s\n", class(single_contrast_tbl)))
                  cat(sprintf("   - input$formula_string: %s\n", input$formula_string))
                  cat(sprintf("   - input$da_q_val_thresh: %s\n", input$da_q_val_thresh))
                  cat(sprintf("   - input$treat_lfc_cutoff: %s\n", input$treat_lfc_cutoff))

                  stop(e)
                }
              )

              cat(sprintf("   DA ANALYSIS Step: Successfully completed contrast %d\n", contrast_idx))
              result
            })

          cat(sprintf("   DA ANALYSIS Step: Completed DA analysis for %d contrasts\n", length(da_results_list)))
          cat(sprintf("   DA ANALYSIS Step: DA results list names: %s\n", paste(names(da_results_list), collapse = ", ")))

          # Step 3: Combine results into expected format for UI components
          shiny::incProgress(0.8, detail = "Processing results...")

          # Update UI dropdowns to use friendly names that match the data's comparison column
          cat("   DA ANALYSIS Step: Updating UI dropdowns to match data comparison column...\n")

          # Use friendly names from contrasts_tbl since they match da_proteins_long$comparison
          if ("friendly_names" %in% names(contrasts_tbl)) {
            friendly_names <- contrasts_tbl$friendly_names
            contrast_choices <- setNames(friendly_names, friendly_names)
            cat(sprintf("   DA ANALYSIS Step: Using friendly_names from contrasts_tbl: %s\n", paste(friendly_names, collapse = ", ")))
          } else {
            # Fallback: extract friendly names from full_format (part before =)
            full_format_strings <- contrasts_tbl$full_format
            friendly_names <- stringr::str_extract(full_format_strings, "^[^=]+")
            contrast_choices <- setNames(friendly_names, friendly_names)
            cat(sprintf("   DA ANALYSIS Step: Extracted friendly names from full_format: %s\n", paste(friendly_names, collapse = ", ")))
          }

          shiny::updateSelectInput(session, "volcano_contrast", choices = contrast_choices)
          shiny::updateSelectInput(session, "heatmap_contrast", choices = contrast_choices)
          shiny::updateSelectInput(session, "table_contrast", choices = contrast_choices)
          cat(sprintf("   DA ANALYSIS Step: Updated UI dropdowns with friendly names: %s\n", paste(names(contrast_choices), collapse = ", ")))
          cat(sprintf("   DA ANALYSIS Step: (These match the comparison column in da_proteins_long)\n"))

          # Combine all da_proteins_long results into a single dataframe
          cat("   DA ANALYSIS Step: Combining da_proteins_long results from all contrasts...\n")
          combined_da_proteins_long <- da_results_list |>
            purrr::keep(~ !is.null(.x$da_proteins_long)) |>
            purrr::map(~ .x$da_proteins_long) |>
            purrr::list_rbind()

          cat(sprintf("   DA ANALYSIS Step: Combined da_proteins_long has %d total rows\n", nrow(combined_da_proteins_long)))

          # Create the expected structure for UI components
          combined_results <- list(
            da_proteins_long = combined_da_proteins_long,
            individual_contrasts = da_results_list # Keep individual results for other purposes
          )

          # Add other shared elements from the first contrast result (they should be the same)
          if (length(da_results_list) > 0) {
            first_result <- da_results_list[[1]]
            if (!is.null(first_result$theObject)) {
              combined_results$theObject <- first_result$theObject
            }
          }

          da_data$da_results_list <- combined_results
          da_data$analysis_complete <- TRUE
          
          # --- TESTTHAT CHECKPOINT CP07 (see test-prot-07-da-analysis.R) ---
          .capture_checkpoint(combined_results, "cp07", "da_results")
          # --- END CP07 ---

          # Check for qvalue() failures and show prominent warning notification
          all_qvalue_warnings <- da_results_list |>
            purrr::map(~ .x$qvalue_warnings) |>
            purrr::compact() |>
            purrr::list_flatten()

          if (length(all_qvalue_warnings) > 0) {
            # Map contrast names to friendly names for user-friendly message
            failed_contrasts <- names(all_qvalue_warnings)
            if (exists("contrasts_tbl", envir = .GlobalEnv)) {
              contrasts_tbl <- get("contrasts_tbl", envir = .GlobalEnv)
              friendly_failed_names <- purrr::map_chr(failed_contrasts, function(failed_contrast) {
                # Extract the part before = if it exists
                contrast_base <- stringr::str_extract(failed_contrast, "^[^=]+")
                if (is.na(contrast_base)) contrast_base <- failed_contrast

                # Try to find matching friendly name
                if ("friendly_names" %in% names(contrasts_tbl) && "contrasts" %in% names(contrasts_tbl)) {
                  match_idx <- which(contrasts_tbl$contrasts == contrast_base)
                  if (length(match_idx) > 0) {
                    return(contrasts_tbl$friendly_names[match_idx[ match_idx[1] ]])
                  }
                }
                return(contrast_base)
              })
            } else {
              friendly_failed_names <- failed_contrasts
            }

            # Create prominent warning modal dialog
            warning_title <- paste0("[WARNING] IMPORTANT: Statistical Analysis Warning")
            warning_body <- paste0(
              "<div style='font-size: 14px; line-height: 1.6;'>",
              "<p><strong>The q-value calculation failed for ", length(all_qvalue_warnings), " contrast(s):</strong></p>",
              "<p style='margin-left: 20px;'>", paste(friendly_failed_names, collapse = ", "), "</p>",
              "<p>The analysis used <strong>Benjamini-Hochberg FDR correction (p.adjust)</strong> instead.</p>",
              "<hr>",
              "<p><strong>[WARNING] INTERPRET RESULTS WITH CAUTION:</strong></p>",
              "<ul style='margin-left: 20px;'>",
              "<li>The p-value distribution may be problematic</li>",
              "<li>Results may be less reliable than normal</li>",
              "<li>Consider reviewing your experimental design and data quality</li>",
              "<li>Check diagnostic messages in the console for details</li>",
              "</ul>",
              "</div>"
            )

            # Show prominent modal dialog
            shiny::showModal(
              shiny::modalDialog(
                title = shiny::tags$div(shiny::tags$strong(warning_title), style = "color: #d9534f; font-size: 18px;"),
                shiny::HTML(warning_body),
                size = "l", # Large modal
                easyClose = FALSE, # User must click button to close
                footer = shiny::tagList(
                  shiny::actionButton(ns("acknowledge_qvalue_warning"),
                    "I Understand - Continue",
                    class = "btn-warning",
                    style = "font-weight: bold;"
                  )
                )
              )
            )

            # Handle acknowledgment button
            shiny::observeEvent(input$acknowledge_qvalue_warning,
              {
                shiny::removeModal()
              },
              once = TRUE
            )

            cat(sprintf("   DA ANALYSIS Step: [WARNING] qvalue() failed for %d contrast(s) - user notification shown\n", length(all_qvalue_warnings)))
          }

          # Update workflow data
          workflow_data$da_analysis_results_list <- da_results_list
          # Must replace entire list to trigger reactivity
          updated_status <- workflow_data$tab_status
          updated_status$differential_abundance <- "complete"
          updated_status$enrichment_analysis <- "pending"
          workflow_data$tab_status <- updated_status

          # [OK] FIXED: Write DA results to disk using new S4 method
          shiny::incProgress(0.9, detail = "Writing results to disk...")

          cat("   DA ANALYSIS Step: Writing DA results to disk for enrichment analysis...\n")

          tryCatch(
            {
              # Get UniProt annotations for output
              uniprot_tbl <- NULL
              if (exists("uniprot_dat_cln", envir = .GlobalEnv)) {
                uniprot_tbl <- get("uniprot_dat_cln", envir = .GlobalEnv)
                cat("   DA ANALYSIS Step: Found uniprot_dat_cln for annotations\n")
              } else {
                cat("   DA ANALYSIS Step: No uniprot_dat_cln found - proceeding without annotations\n")
              }

              # Get the S4 object from the first contrast (they should all be the same)
              first_result <- da_results_list[[1]]
              s4_object <- first_result$theObject

              # Use the properly constructed da_output_dir from experiment_paths
              da_output_dir <- experiment_paths$da_output_dir
              cat(sprintf("   DA ANALYSIS Step: Using da_output_dir from experiment_paths: %s\n", da_output_dir))

              # Update S4 object parameters for file writing
              s4_object@args$outputDaResultsAllContrasts <- list(
                uniprot_tbl = uniprot_tbl,
                da_output_dir = da_output_dir,
                publication_graphs_dir = experiment_paths$publication_graphs_dir,
                file_prefix = "da_proteins",
                args_row_id = s4_object@protein_id_column,
                gene_names_column = "gene_names",
                uniprot_id_column = "Entry",
                da_q_val_thresh = input$da_q_val_thresh,
                fdr_column = "fdr_qvalue",
                log2fc_column = "log2FC"
              )

              # Call the new S4 method that writes all contrasts properly
              cat(sprintf("   DA ANALYSIS Step: Calling outputDaResultsAllContrasts for %d contrasts\n", length(da_results_list)))

              success <- outputDaResultsAllContrasts(
                theObject = s4_object,
                da_results_list_all_contrasts = da_results_list,
                uniprot_tbl = uniprot_tbl,
                da_output_dir = da_output_dir,
                publication_graphs_dir = experiment_paths$publication_graphs_dir,
                file_prefix = "da_proteins",
                args_row_id = s4_object@protein_id_column,
                gene_names_column = "gene_names",
                uniprot_id_column = "Entry"
              )

              if (success) {
                cat("   DA ANALYSIS Step: All DA results written to disk successfully\n")
              } else {
                cat("   DA ANALYSIS Step: ERROR - outputDaResultsAllContrasts returned FALSE\n")
              }

              # [OK] REMOVED: Per-contrast outputDeAnalysisResults calls - not needed!
              # outputDaResultsAllContrasts already handles everything including volcano plots
              cat("   DA ANALYSIS Step: All output handled by outputDaResultsAllContrasts\n")
            },
            error = function(e) {
              cat(sprintf("   DA ANALYSIS Step: Error writing results to disk: %s\n", e$message))
              # Don't fail the entire analysis for file writing issues
              # NOTE: Removed log_warn call to avoid logger interpolation bug in error handlers
              message(paste("Could not write DA results to disk:", e$message))
            }
          )

          shiny::incProgress(1.0, detail = "Complete!")
        })

        shiny::showNotification(
          "Differential abundance analysis completed successfully!",
          type = "message",
          duration = 5
        )

        cat("=== DIFFERENTIAL ABUNDANCE ANALYSIS COMPLETED ===\n")
      },
      error = function(e) {
        cat(paste("*** ERROR in DA analysis:", e$message, "\n"))
        shiny::showNotification(
          paste("Error in DA analysis:", e$message),
          type = "error",
          duration = 10
        )
      }
    )

    shiny::removeNotification("da_working")
  })
}

