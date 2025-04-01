  # Confirm reset handler with scoped reset functionality
  observeEvent(input$confirm_reset, {
    # Perform reset based on selected scope
    scope <- input$reset_scope
    
    if(scope == "all" || scope == "sample_names") {
        if(scope == "all") {
            # Reset entire matrices if scope is "all" using non-reactive initial_state
            design_matrix(initial_state$design_matrix)
            data_cln_reactive(initial_state$data_cln)
        } else { # scope == "sample_names"
            # Get the current reactive data frames
            current_matrix <- design_matrix()
            current_data <- data_cln_reactive()

            # Capture current Run names from design_matrix BEFORE reset
            current_run_names_in_matrix <- current_matrix$Run 
            initial_run_names_in_matrix <- initial_state$design_matrix$Run

            # --- Reset design_matrix$Run (assuming row order preserved) ---
            if (nrow(current_matrix) == nrow(initial_state$design_matrix)) {
                current_matrix$Run <- initial_run_names_in_matrix
            } else {
                # Fallback if rows don't match (shouldn't happen ideally)
                showNotification("Error: Row mismatch during design matrix reset. Resetting entire matrix.", type = "error", duration = 5)
                current_matrix <- initial_state$design_matrix 
            }

            # --- Reset data_cln$Run using mapping based on design_matrix changes ---
            # Create a map: current name -> initial name
            if (length(current_run_names_in_matrix) == length(initial_run_names_in_matrix)) {
                # Check for duplicate current names, which makes the map ambiguous
                if(any(duplicated(current_run_names_in_matrix))) {
                    showNotification("Warning: Duplicate sample names detected. Resetting data_cln based on initial state.", type = "warning", duration = 5)
                    # Fallback: Reset data_cln Run names based on initial_state$data_cln
                    if (nrow(current_data) == nrow(initial_state$data_cln)) {
                        current_data$Run <- initial_state$data_cln$Run
                    } else {
                        current_data <- initial_state$data_cln # Full reset if rows mismatch
                    }
                } else {
                    # Create the map only if names are unique
                    run_map_current_to_original <- setNames(initial_run_names_in_matrix, current_run_names_in_matrix)
                    
                    # Apply the mapping to data_cln
                    current_run_names_in_data <- current_data$Run
                    original_run_names_for_data <- run_map_current_to_original[current_run_names_in_data]
                    
                    # Only replace if a mapping was found, otherwise keep current name
                    current_data$Run <- ifelse(!is.na(original_run_names_for_data), 
                                               original_run_names_for_data, 
                                               current_run_names_in_data)
                }
            } else {
                 showNotification("Error: Name vector length mismatch during data reset. Resetting data_cln.", type = "error", duration = 5)
                 current_data <- initial_state$data_cln # Fallback
            }

            # Update the reactive values
            design_matrix(current_matrix)
            data_cln_reactive(current_data)
        }
    }

    if(scope == "all" || scope == "factors") {
      # Reset factors using the non-reactive initial_state list
      factors(initial_state$factors) 
      
      if(scope == "factors") {
        # When resetting just factors, keep existing assignments in design matrix
        # but remove any factors that no longer exist
        current_matrix <- design_matrix()
        # Use factors from the non-reactive initial_state list
        valid_factors <- initial_state$factors 
        
        current_matrix$factor1[!current_matrix$factor1 %in% valid_factors] <- NA
        current_matrix$factor2[!current_matrix$factor2 %in% valid_factors] <- NA
        
        design_matrix(current_matrix)
      }
    }
    
    if(scope == "all" || scope == "contrasts") {
      contrasts(data.frame(
        contrast_name = character(),
        numerator = character(),
        denominator = character(),
        factor1_num = character(),
        factor2_num = character(),
        factor1_den = character(),
        factor2_den = character(),
        stringsAsFactors = FALSE
      ))
    }
    
    if(scope == "all" || scope == "formula") {
      updateTextInput(session, "formula_string", 
                      # Reset formula using the non-reactive initial_state list
                      value = initial_state$formula) 
    }
    
    # Update UI elements
    if(scope == "all" || scope == "sample_names") {
       # Get the freshly reset runs (or initial runs if scope==all)
       # If scope was 'sample_names', design_matrix() now holds the reset matrix.
       dm_to_sort <- if(scope == "all") initial_state$design_matrix else design_matrix()
       sorted_runs <- gtools::mixedsort(dm_to_sort$Run) 
       updateSelectizeInput(session, "sample_to_rename", choices = sorted_runs, selected = "") # Clear selection
       updateSelectizeInput(session, "selected_runs", choices = sorted_runs, selected = "") # Clear selection
       updateSelectizeInput(session, "samples_to_transform", choices = sorted_runs, selected = "") # Clear selection
    }
    
    if(scope == "all") {
      # Reset groups using the non-reactive initial_state list
      groups(initial_state$groups) 
      updateTabsetPanel(session, "main_tabs", selected = "Rename Samples")
      # Also reset factor dropdowns if scope is all
       updateSelectInput(session, "factor1_select", choices = c("", initial_state$factors), selected = "")
       updateSelectInput(session, "factor2_select", choices = c("", initial_state$factors), selected = "")
       updateSelectInput(session, "contrast_group1", choices = c("", initial_state$groups), selected = "")
       updateSelectInput(session, "contrast_group2", choices = c("", initial_state$groups), selected = "")
    }
    
    # Close the modal
    removeModal()
    
    # Show success message
    scope_desc <- switch(input$reset_scope,
                       "all" = "all changes",
                       "sample_names" = "sample names",
                       "factors" = "factor definitions",
                       "contrasts" = "contrast definitions",
                       "formula" = "formula")
    showNotification(paste("Reset of", scope_desc, "completed successfully"), 
                   type = "message", 
                   duration = 3)
  }) 