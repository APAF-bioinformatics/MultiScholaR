#' Internal initialization handler for Proteomics DA module
#' @keywords internal
da_server_init_handlers <- function(input, output, session, da_data, workflow_data, selected_tab) {
  # CRITICAL FIX: Add tab selection observer like normalization module
  if (!is.null(selected_tab)) {
    cat("   mod_prot_da_server Step: Setting up tab selection observer\n")
    shiny::observeEvent(selected_tab(),
      {
        cat("--- Entering tab selection observer ---\n")
        cat(sprintf("   tab_observer Step: Selected tab = %s\n", selected_tab()))

        # Only trigger if DA tab is selected
        # Only trigger if DA tab is selected (recognize both 'da' and 'de' IDs)
        if (!is.null(selected_tab()) && selected_tab() == "da") {
          cat("=== DA TAB CLICKED ===\n")
          cat(sprintf("   DA TAB Step: workflow_data$state_manager is NULL = %s\n", is.null(workflow_data$state_manager)))

          if (!is.null(workflow_data$state_manager)) {
            current_state <- workflow_data$state_manager$current_state

            # Define valid states where this tab can be active
            valid_states_for_de_tab <- c("normalized", "ruv_corrected", "correlation_filtered")

            cat(sprintf("   DA TAB Step: Current state = '%s'\n", current_state))
            cat(sprintf("   DA TAB Step: Valid states for DE = %s\n", paste(valid_states_for_de_tab, collapse = ", ")))

            # Auto-trigger only fires if we've completed correlation filtering or at least normalization
            if (current_state %in% valid_states_for_de_tab) {
              cat("*** AUTO-TRIGGERING DE INITIALIZATION (correlation-filtered state found) ***\n")

              tryCatch(
                {
                  # Initialize DE analysis setup
                  cat("   DA TAB Step: Getting S4 object from state manager...\n")
                  current_s4 <- workflow_data$state_manager$getState(current_state)

                  if (!is.null(current_s4)) {
                    cat(sprintf("   DA TAB Step: S4 object retrieved, class = %s\n", class(current_s4)))
                    da_data$current_s4_object <- current_s4

                    # Check for contrasts_tbl in global environment FIRST
                    cat("   DA TAB Step: Checking for contrasts_tbl in global environment...\n")
                    cat(sprintf("   DA TAB Step: contrasts_tbl exists in global env: %s\n", exists("contrasts_tbl", envir = .GlobalEnv)))
                    if (exists("contrasts_tbl", envir = .GlobalEnv)) {
                      contrasts_tbl <- get("contrasts_tbl", envir = .GlobalEnv)
                      cat("   DA TAB Step: Found contrasts_tbl in global environment\n")
                      cat("   DA TAB Step: contrasts_tbl structure:\n")
                      str(contrasts_tbl)
                      cat("   DA TAB Step: contrasts_tbl content:\n")
                      print(contrasts_tbl)

                      # Validate contrasts_tbl has content
                      if (is.null(contrasts_tbl) || nrow(contrasts_tbl) == 0) {
                        cat("   DA TAB Step: contrasts_tbl exists but is empty, falling back to auto-generation\n")
                        da_data$contrasts_available <- NULL
                      } else if ("comparison" %in% names(contrasts_tbl)) {
                        da_data$contrasts_available <- contrasts_tbl$comparison
                        cat(sprintf("   DA TAB Step: Set contrasts from comparison column: %s\n", paste(da_data$contrasts_available, collapse = ", ")))
                      } else if ("contrasts" %in% names(contrasts_tbl)) {
                        da_data$contrasts_available <- contrasts_tbl$contrasts
                        cat(sprintf("   DA TAB Step: Set contrasts from contrasts column: %s\n", paste(da_data$contrasts_available, collapse = ", ")))
                      } else {
                        cat("   DA TAB Step: contrasts_tbl found but no recognized column names\n")
                        cat("   DA TAB Step: Available column names:\n")
                        print(names(contrasts_tbl))
                        # Try first column if it has content
                        if (ncol(contrasts_tbl) > 0) {
                          da_data$contrasts_available <- contrasts_tbl[[1]]
                          cat(sprintf("   DA TAB Step: Using first column: %s\n", paste(da_data$contrasts_available, collapse = ", ")))
                        } else {
                          cat("   DA TAB Step: contrasts_tbl has no usable content\n")
                          da_data$contrasts_available <- NULL
                        }
                      }
                    } else {
                      cat("   DA TAB Step: No contrasts_tbl in global environment. Will attempt auto-generation.\n")
                      da_data$contrasts_available <- NULL
                    }

                    # Only auto-generate contrasts if user-specified ones weren't found or were empty
                    if (is.null(da_data$contrasts_available) || length(da_data$contrasts_available) == 0) {
                      cat("   DA TAB Step: No valid user-specified contrasts found, creating basic contrasts\n")
                      shiny::showNotification("No user-specified contrasts found. Auto-generating all possible contrasts.", type = "warning", duration = 5)

                      if (!is.null(current_s4@design_matrix)) {
                        groups <- unique(current_s4@design_matrix$group)
                        if (length(groups) >= 2) {
                          group_prefixed <- paste0("group", groups)
                          basic_contrasts <- combn(group_prefixed, 2, function(x) paste0(x[1], "-", x[2]), simplify = TRUE)
                          da_data$contrasts_available <- basic_contrasts
                          cat(sprintf("   DA TAB Step: Created basic contrasts: %s\n", paste(basic_contrasts, collapse = ", ")))
                        }
                      }
                    }

                    # Extract formula from S4 @args
                    cat("   DA TAB Step: Checking for formula in S4 @args...\n")
                    if ("deAnalysisParameters" %in% names(current_s4@args)) {
                      if ("formula_string" %in% names(current_s4@args$deAnalysisParameters)) {
                        formula_from_s4 <- current_s4@args$deAnalysisParameters$formula_string
                        da_data$formula_from_s4 <- formula_from_s4
                        cat(sprintf("   DA TAB Step: Formula from S4 = %s\n", formula_from_s4))

                        # Update UI with formula from S4
                        shiny::updateTextAreaInput(
                          session,
                          "formula_string",
                          value = formula_from_s4
                        )
                      } else {
                        cat("   DA TAB Step: No formula_string in deAnalysisParameters\n")
                      }
                    } else {
                      cat("   DA TAB Step: No deAnalysisParameters in S4 @args\n")
                    }

                    # CRITICAL FIX: Validate and fix contrast format to match design matrix columns
                    if (!is.null(da_data$contrasts_available) && !is.null(current_s4@design_matrix)) {
                      cat("   DA TAB Step: Validating contrast format against design matrix...\n")
                      groups <- unique(current_s4@design_matrix$group)
                      cat(sprintf("   DA TAB Step: Available groups: %s\n", paste(groups, collapse = ", ")))

                      # Check if contrasts need "group" prefix (for formula ~ 0 + group)
                      current_contrasts <- da_data$contrasts_available
                      cat(sprintf("   DA TAB Step: Current contrasts before validation: %s\n", paste(current_contrasts, collapse = ", ")))

                      # If contrasts don't start with "group" but should (based on formula), fix them
                      if (any(grepl("~ 0 \\+ group", input$formula_string))) {
                        cat("   DA TAB Step: Formula uses '~ 0 + group', checking if contrasts need group prefix...\n")

                        # Check if any contrast references raw group names without prefix
                        needs_prefix <- any(sapply(groups, function(g) any(grepl(paste0("\\b", g, "\\b"), current_contrasts))))

                        if (needs_prefix) {
                          cat("   DA TAB Step: Contrasts need group prefix, fixing...\n")
                          fixed_contrasts <- purrr::reduce(groups, function(res, group) {
                            gsub(paste0("\\b", group, "\\b"), paste0("group", group), res)
                          }, .init = current_contrasts)
                          da_data$contrasts_available <- fixed_contrasts
                          cat(sprintf("   DA TAB Step: Fixed contrasts: %s\n", paste(fixed_contrasts, collapse = ", ")))
                        } else {
                          cat("   DA TAB Step: Contrasts already have correct format\n")
                        }
                      }
                    }
                  } else {
                    cat("   DA TAB Step: S4 object is NULL\n")
                  }

                  cat("*** DE INITIALIZATION COMPLETED SUCCESSFULLY ***\n")
                },
                error = function(e) {
                  cat(paste("*** ERROR in DE initialization:", e$message, "\n"))
                  shiny::showNotification(
                    paste("Error initializing DE analysis:", e$message),
                    type = "error",
                    duration = 10
                  )
                }
              )
            } else {
              cat(sprintf("*** State '%s' is not valid for DE analysis. User needs to complete normalization (with or without RUV) and correlation filtering. ***\n", current_state))
              shiny::showNotification(
                "Please complete the normalization (with or without RUV) and correlation filtering steps before accessing differential abundance analysis.",
                type = "warning",
                duration = 5
              )
            }
          } else {
            cat("*** workflow_data$state_manager is NULL - cannot check state ***\n")
          }
        } else {
          cat(sprintf("   tab_observer Step: Tab '%s' is not DA tab, ignoring\n", selected_tab()))
        }

        cat("--- Exiting tab selection observer ---\n")
      },
      ignoreInit = TRUE
    )
  } else {
    cat("   mod_prot_da_server Step: No selected_tab parameter provided - tab selection observer NOT set up\n")
  }

  # BACKUP: Fetch formula and contrasts from S4 object when state updates (fallback method)
  shiny::observeEvent(workflow_data$state_update_trigger,
    {
      cat("--- Entering state update trigger observer ---\n")
      cat("\n\n=== DA TAB TRIGGERED VIA STATE UPDATE: Checking for S4 object and contrasts ===\n")

      # Get current S4 object from state manager
      if (!is.null(workflow_data$state_manager)) {
        current_state <- workflow_data$state_manager$current_state
        cat(sprintf("   DA TAB Step: Current state = %s\n", current_state))

        # Should be getting the correlation-filtered protein object (final state before DE)
        # Also accept normalized or ruv_corrected states
        if (current_state %in% c("normalized", "ruv_corrected", "correlation_filtered")) {
          cat(sprintf("   DA TAB Step: State is valid for DE analysis (%s state found)\n", current_state))
          current_s4 <- workflow_data$state_manager$getState(current_state)

          if (!is.null(current_s4)) {
            cat(sprintf("   DA TAB Step: S4 object retrieved, class = %s\n", class(current_s4)))
            da_data$current_s4_object <- current_s4

            # Extract formula from S4 @args
            if ("deAnalysisParameters" %in% names(current_s4@args)) {
              if ("formula_string" %in% names(current_s4@args$deAnalysisParameters)) {
                formula_from_s4 <- current_s4@args$deAnalysisParameters$formula_string
                da_data$formula_from_s4 <- formula_from_s4
                cat(sprintf("   DA TAB Step: Formula from S4 = %s\n", formula_from_s4))

                # Update UI with formula from S4
                shiny::updateTextAreaInput(
                  session,
                  "formula_string",
                  value = formula_from_s4
                )
              } else {
                cat("   DA TAB Step: No formula_string in deAnalysisParameters\n")
              }
            } else {
              cat("   DA TAB Step: No deAnalysisParameters in S4 @args\n")
            }

            # Get contrasts from design matrix or S4 args
            if (!is.null(current_s4@design_matrix)) {
              cat("   DA TAB Step: Design matrix found in S4 object\n")
              cat(sprintf("   DA TAB Step: Design matrix dims = %d rows, %d cols\n", nrow(current_s4@design_matrix), ncol(current_s4@design_matrix)))

              # Check for contrasts_tbl in global environment
              cat(sprintf("   DA TAB Step: contrasts_tbl exists in global env: %s\n", exists("contrasts_tbl", envir = .GlobalEnv)))
              if (exists("contrasts_tbl", envir = .GlobalEnv)) {
                contrasts_tbl <- get("contrasts_tbl", envir = .GlobalEnv)
                cat("   DA TAB Step: Found contrasts_tbl in global environment\n")
                cat("   DA TAB Step: contrasts_tbl structure:\n")
                str(contrasts_tbl)
                cat("   DA TAB Step: contrasts_tbl content:\n")
                print(contrasts_tbl)

                # Validate contrasts_tbl has content
                if (is.null(contrasts_tbl) || nrow(contrasts_tbl) == 0) {
                  cat("   DA TAB Step: contrasts_tbl exists but is empty, falling back to auto-generation\n")
                  da_data$contrasts_available <- NULL
                } else if ("comparison" %in% names(contrasts_tbl)) {
                  da_data$contrasts_available <- contrasts_tbl$comparison
                  cat(sprintf("   DA TAB Step: Set contrasts from comparison column: %s\n", paste(da_data$contrasts_available, collapse = ", ")))
                } else if ("contrasts" %in% names(contrasts_tbl)) {
                  da_data$contrasts_available <- contrasts_tbl$contrasts
                  cat(sprintf("   DA TAB Step: Set contrasts from contrasts column: %s\n", paste(da_data$contrasts_available, collapse = ", ")))
                } else {
                  cat("   DA TAB Step: contrasts_tbl found but no recognized column names\n")
                  cat("   DA TAB Step: Available column names:\n")
                  print(names(contrasts_tbl))
                  # Try first column if it has content
                  if (ncol(contrasts_tbl) > 0) {
                    da_data$contrasts_available <- contrasts_tbl[[1]]
                    cat(sprintf("   DA TAB Step: Using first column: %s\n", paste(da_data$contrasts_available, collapse = ", ")))
                  } else {
                    cat("   DA TAB Step: contrasts_tbl has no usable content\n")
                    da_data$contrasts_available <- NULL
                  }
                }
              } else {
                cat("   DA TAB Step: No contrasts_tbl in global environment. Will attempt auto-generation.\n")
                da_data$contrasts_available <- NULL
              }

              # Only auto-generate contrasts if user-specified ones weren't found or were empty
              if (is.null(da_data$contrasts_available) || length(da_data$contrasts_available) == 0) {
                cat("   DA TAB Step: No valid user-specified contrasts found, creating basic contrasts\n")
                shiny::showNotification("No user-specified contrasts found. Auto-generating all possible contrasts.", type = "warning", duration = 5)

                if (!is.null(current_s4@design_matrix)) {
                  groups <- unique(current_s4@design_matrix$group)
                  if (length(groups) >= 2) {
                    group_prefixed <- paste0("group", groups)
                    basic_contrasts <- combn(group_prefixed, 2, function(x) paste0(x[1], "-", x[2]), simplify = TRUE)
                    da_data$contrasts_available <- basic_contrasts
                    cat(sprintf("   DA TAB Step: Created basic contrasts: %s\n", paste(basic_contrasts, collapse = ", ")))
                  }
                }
              }
            } else {
              cat("   DA TAB Step: No design matrix found in S4 object\n")
            }
          } else {
            cat("   DA TAB Step: S4 object is NULL\n")
          }
        } else {
          cat(sprintf("   DA TAB Step: State '%s' not valid for DA analysis (expecting: normalized, ruv_corrected, or correlation_filtered)\n", current_state))
        }
      } else {
        cat("   DA TAB Step: workflow_data$state_manager is NULL\n")
      }

      cat("=== DA TAB: Contrast detection complete ===\n")
      cat("--- Exiting state update trigger observer ---\n")
    },
    ignoreInit = TRUE
  )

  # Update contrast dropdowns when contrasts become available
  shiny::observe({
    if (!is.null(da_data$contrasts_available)) {
      contrast_choices <- setNames(da_data$contrasts_available, da_data$contrasts_available)

      # Update all contrast selectors
      shiny::updateSelectInput(session, "volcano_contrast", choices = contrast_choices)
      shiny::updateSelectInput(session, "heatmap_contrast", choices = contrast_choices)
      shiny::updateSelectInput(session, "table_contrast", choices = contrast_choices)
    }
  })

  # Display available contrasts
  output$contrasts_display <- shiny::renderText({
    if (!is.null(da_data$contrasts_available)) {
      paste(da_data$contrasts_available, collapse = "\n")
    } else {
      "No contrasts available.\nComplete normalization and\ncorrelation filtering first."
    }
  })

  # Display analysis status
  output$da_status <- shiny::renderText({
    if (da_data$analysis_complete) {
      paste(
        "[OK] Analysis Complete\n",
        sprintf("Contrasts analyzed: %d\n", length(da_data$contrasts_available)),
        "Results available in all tabs"
      )
    } else {
      "[WAITING] Waiting for analysis...\nClick 'Run DA Analysis' to start"
    }
  })
}
