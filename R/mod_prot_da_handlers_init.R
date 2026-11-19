#' Internal initialization handler for Proteomics DA module
#' @keywords internal
da_server_init_handlers <- function(input, output, session, da_data, workflow_data, selected_tab) {
  # CRITICAL FIX: Add tab selection observer like normalization module
  if (!is.null(selected_tab)) {
    shiny::observeEvent(selected_tab(),
      {
        withDIAStateScope(getDIAActiveManager(), function() {

        # Only trigger if DA tab is selected
        # Only trigger if DA tab is selected (recognize both 'da' and 'de' IDs)
        if (!is.null(selected_tab()) && selected_tab() == "da") {
          cat("=== DA TAB CLICKED ===\n")

          if (!is.null(workflow_data$state_manager)) {
            current_state <- workflow_data$state_manager$current_state

            # Define valid states where this tab can be active
            valid_states_for_de_tab <- c("ruv_covariates", "normalized", "ruv_corrected", "correlation_filtered")


            # Auto-trigger only fires if we've completed correlation filtering or at least normalization
            if (current_state %in% valid_states_for_de_tab) {
              cat("*** AUTO-TRIGGERING DE INITIALIZATION (correlation-filtered state found) ***\n")

              tryCatch(
                {
                  # Initialize DE analysis setup
                  current_s4 <- readDIAStateGuiMetadata(workflow_data$state_manager, current_state)

                  if (!is.null(current_s4)) {
                    da_data$current_s4_object <- current_s4$reference %||% current_s4$legacy_object

                    # Check for contrasts_tbl on the workflow bus FIRST
                    contrasts_tbl <- workflow_data$contrasts_tbl
                    if (!is.null(contrasts_tbl)) {
                      # Validate contrasts_tbl has content
                      if (is.null(contrasts_tbl) || nrow(contrasts_tbl) == 0) {
                        da_data$contrasts_available <- NULL
                      } else if ("comparison" %in% names(contrasts_tbl)) {
                        da_data$contrasts_available <- contrasts_tbl$comparison
                      } else if ("contrasts" %in% names(contrasts_tbl)) {
                        da_data$contrasts_available <- contrasts_tbl$contrasts
                      } else {
                        # Try first column if it has content
                        if (ncol(contrasts_tbl) > 0) {
                          da_data$contrasts_available <- contrasts_tbl[[1]]
                        } else {
                          da_data$contrasts_available <- NULL
                        }
                      }
                    } else {
                      da_data$contrasts_available <- NULL
                    }

                    # Only auto-generate contrasts if user-specified ones weren't found or were empty
                    if (is.null(da_data$contrasts_available) || length(da_data$contrasts_available) == 0) {
                      shiny::showNotification("No user-specified contrasts found. Auto-generating all possible contrasts.", type = "warning", duration = 5)

                      if (!is.null(current_s4$design_matrix)) {
                        groups <- unique(current_s4$design_matrix$group)
                        if (length(groups) >= 2) {
                          group_prefixed <- paste0("group", groups)
                          basic_contrasts <- combn(group_prefixed, 2, function(x) paste0(x[1], "-", x[2]), simplify = TRUE)
                          da_data$contrasts_available <- basic_contrasts
                        }
                      }
                    }

                    # Extract formula from S4 @args
                    if ("deAnalysisParameters" %in% names(current_s4$args)) {
                      if ("formula_string" %in% names(current_s4$args$deAnalysisParameters)) {
                        formula_from_s4 <- current_s4$args$deAnalysisParameters$formula_string
                        da_data$formula_from_s4 <- formula_from_s4

                        # Update UI with formula from S4
                        shiny::updateTextAreaInput(
                          session,
                          "formula_string",
                          value = formula_from_s4
                        )
                      }
                    }

                    # CRITICAL FIX: Validate and fix contrast format to match design matrix columns
                    if (!is.null(da_data$contrasts_available) && !is.null(current_s4$design_matrix)) {
                      groups <- unique(current_s4$design_matrix$group)

                      # Check if contrasts need "group" prefix (for formula ~ 0 + group)
                      current_contrasts <- da_data$contrasts_available

                      # If contrasts don't start with "group" but should (based on formula), fix them
                      if (any(grepl("~ 0 \\+ group", input$formula_string))) {

                        # Check if any contrast references raw group names without prefix
                        needs_prefix <- any(sapply(groups, function(g) any(grepl(paste0("\\b", g, "\\b"), current_contrasts))))

                        if (needs_prefix) {
                          fixed_contrasts <- purrr::reduce(groups, function(res, group) {
                            gsub(paste0("\\b", group, "\\b"), paste0("group", group), res)
                          }, .init = current_contrasts)
                          da_data$contrasts_available <- fixed_contrasts
                        }
                      }
                    }
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
        }

        })
      },
      ignoreInit = TRUE
    )
  }

  # BACKUP: Fetch formula and contrasts from S4 object when state updates (fallback method)
  shiny::observeEvent(workflow_data$state_update_trigger,
    {
      withDIAStateScope(getDIAActiveManager(), function() {
      cat("\n\n=== DA TAB TRIGGERED VIA STATE UPDATE: Checking for S4 object and contrasts ===\n")

      # Get current S4 object from state manager
      if (!is.null(workflow_data$state_manager)) {
        current_state <- workflow_data$state_manager$current_state
        invalidateDIAVisibleResults(da_data, workflow_data)

        # Should be getting the correlation-filtered protein object (final state before DE)
        # Also accept normalized or ruv_corrected states
        if (current_state %in% c("ruv_covariates", "normalized", "ruv_corrected", "correlation_filtered")) {
          current_s4 <- readDIAStateGuiMetadata(workflow_data$state_manager, current_state)

          if (!is.null(current_s4)) {
            da_data$current_s4_object <- current_s4$reference %||% current_s4$legacy_object

            # Extract formula from S4 @args
            if ("deAnalysisParameters" %in% names(current_s4$args)) {
              if ("formula_string" %in% names(current_s4$args$deAnalysisParameters)) {
                formula_from_s4 <- current_s4$args$deAnalysisParameters$formula_string
                da_data$formula_from_s4 <- formula_from_s4

                # Update UI with formula from S4
                shiny::updateTextAreaInput(
                  session,
                  "formula_string",
                  value = formula_from_s4
                )
              }
            }

            # Get contrasts from design matrix or S4 args
            if (!is.null(current_s4$design_matrix)) {

              # Check for contrasts_tbl on the workflow bus
              contrasts_tbl <- workflow_data$contrasts_tbl
              if (!is.null(contrasts_tbl)) {
                # Validate contrasts_tbl has content
                if (is.null(contrasts_tbl) || nrow(contrasts_tbl) == 0) {
                  da_data$contrasts_available <- NULL
                } else if ("comparison" %in% names(contrasts_tbl)) {
                  da_data$contrasts_available <- contrasts_tbl$comparison
                } else if ("contrasts" %in% names(contrasts_tbl)) {
                  da_data$contrasts_available <- contrasts_tbl$contrasts
                } else {
                  # Try first column if it has content
                  if (ncol(contrasts_tbl) > 0) {
                    da_data$contrasts_available <- contrasts_tbl[[1]]
                  } else {
                    da_data$contrasts_available <- NULL
                  }
                }
              } else {
                da_data$contrasts_available <- NULL
              }

              # Only auto-generate contrasts if user-specified ones weren't found or were empty
              if (is.null(da_data$contrasts_available) || length(da_data$contrasts_available) == 0) {
                shiny::showNotification("No user-specified contrasts found. Auto-generating all possible contrasts.", type = "warning", duration = 5)

                if (!is.null(current_s4$design_matrix)) {
                  groups <- unique(current_s4$design_matrix$group)
                  if (length(groups) >= 2) {
                    group_prefixed <- paste0("group", groups)
                    basic_contrasts <- combn(group_prefixed, 2, function(x) paste0(x[1], "-", x[2]), simplify = TRUE)
                    da_data$contrasts_available <- basic_contrasts
                  }
                }
              }
            }
          }
        }
      }

      cat("=== DA TAB: Contrast detection complete ===\n")
      })
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

  # KIT-106 W5 (D5): the contrasts and status outputs are the shared DA frame's (registerOmicDaPluginTextOutputs).
}
