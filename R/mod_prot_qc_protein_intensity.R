#' @title Protein Intensity Filter Module
#'
#' @description A Shiny module for applying protein intensity and missing value filters.
#'
#' @name mod_prot_qc_protein_intensity
NULL

#' @rdname mod_prot_qc_protein_intensity
#' @export
#' @importFrom shiny NS tagList tabPanel br fluidRow column wellPanel h4 p hr checkboxInput helpText conditionalPanel h5 numericInput div actionButton verbatimTextOutput plotOutput strong
mod_prot_qc_protein_intensity_ui <- function(id) {
  ns <- shiny::NS(id)
  
  shiny::tagList(
    shiny::tabPanel(
    "Protein Intensity Filter",
    shiny::br(),
    shiny::fluidRow(
      shiny::column(4,
        shiny::wellPanel(
          shiny::h4("Missing Value Parameters & Intensity Filter"),
          shiny::p("Configure missing value thresholds and filter proteins on intensity and missing values."),
          shiny::hr(),
          
          # Strict Mode Toggle
          shiny::checkboxInput(ns("use_strict_mode"), 
            "Strict Mode (No Missing Values Allowed)", 
            value = FALSE,
            width = "100%"
          ),
          shiny::helpText("When enabled, removes any protein with missing values in ANY sample across all groups. Overrides flexible threshold settings below."),
          
          shiny::hr(),
          
          # Flexible mode parameters (shown when strict mode is OFF)
          shiny::conditionalPanel(
            condition = paste0("!input['", ns("use_strict_mode"), "']"),
            
            # Missing value parameters (chunk 20)
            shiny::h5("Missing Value Parameters (Flexible Mode)"),
            shiny::numericInput(ns("min_reps_per_group"), 
              "Min Replicates per Group", 
              value = 2, min = 1, max = 10, step = 1,
              width = "100%"
            ),
            shiny::helpText("Minimum valid measurements required within each group (default: 2)"),
            
            shiny::numericInput(ns("min_groups"), 
              "Min Groups Required", 
              value = 2, min = 1, max = 10, step = 1,
              width = "100%"
            ),
            shiny::helpText("Minimum groups that must meet the replicate threshold (default: 2)"),
            
            shiny::hr(),
            
            # Calculated percentages (read-only display)
            shiny::h5("Calculated Thresholds"),
            shiny::p(style = "background-color: #f0f0f0; padding: 10px; border-radius: 5px;",
              shiny::strong("These values are automatically calculated from your replicate settings above:"),
              shiny::br(),
              shiny::textOutput(ns("calculated_groupwise_percent"), inline = TRUE),
              shiny::br(),
              shiny::textOutput(ns("calculated_max_groups_percent"), inline = TRUE)
            ),
            
            shiny::hr()
          ),
          
          # Intensity cutoff (always shown)
          shiny::h5("Intensity Threshold"),
          shiny::numericInput(ns("proteins_intensity_cutoff_percentile"), 
            "Protein Intensity Cutoff Percentile (%)", 
            value = 1, min = 0.1, max = 10, step = 0.1,
            width = "100%"
          ),
          shiny::helpText("Intensity threshold percentile (default: 1%)"),
          
          shiny::hr(),
          shiny::div(
            shiny::actionButton(ns("apply_protein_intensity_filter"), "Apply Filter", 
              class = "btn-primary", width = "48%"),
            shiny::actionButton(ns("revert_protein_intensity_filter"), "Revert", 
              class = "btn-warning", width = "48%", style = "margin-left: 4%")
          )
        )
      ),
      shiny::column(8,
        shiny::verbatimTextOutput(ns("protein_intensity_filter_results")),
        shiny::br(),
        shinyjqui::jqui_resizable(
          shiny::plotOutput(ns("protein_intensity_filter_plot"), height = "800px", width = "100%")
        )
      )
    )
  )
  )
}

#' @rdname mod_prot_qc_protein_intensity
#' @export
#' @importFrom shiny moduleServer reactiveVal observeEvent req showNotification removeNotification renderText renderPlot
#' @importFrom logger log_info log_error
#' @importFrom grid grid.draw
mod_prot_qc_protein_intensity_server <- function(id, workflow_data, omic_type, experiment_label) {
  shiny::moduleServer(id, function(input, output, session) {
    
    protein_intensity_filter_plot <- shiny::reactiveVal(NULL)
    
    # Step 3: Protein Intensity Filter (chunks 20+21 combined)
    shiny::observeEvent(input$apply_protein_intensity_filter, {
      shiny::req(workflow_data$state_manager)
      
      shiny::showNotification("Applying protein intensity filter...", id = "protein_intensity_filter_working", duration = NULL)
      
      tryCatch({
        # Get current ProteinQuantitativeData S4 object from the active state
        current_s4 <- workflow_data$state_manager$getState()
        shiny::req(current_s4)
        
        # Check if using strict mode
        use_strict_mode <- isTRUE(input$use_strict_mode)
        
        if (use_strict_mode) {
          # STRICT MODE: Set parameters to enforce zero tolerance
          logger::log_info("Protein Processing: Using STRICT MODE (no missing values allowed)")
          
          # Set groupwise_percentage_cutoff to 0 (no missing allowed in any group)
          current_s4 <- updateConfigParameter(
            theObject = current_s4,
            function_name = "removeRowsWithMissingValuesPercent",
            parameter_name = "groupwise_percentage_cutoff",
            new_value = 0
          )
          
          # Set max_groups_percentage_cutoff to 0 (all groups must pass)
          current_s4 <- updateConfigParameter(
            theObject = current_s4,
            function_name = "removeRowsWithMissingValuesPercent",
            parameter_name = "max_groups_percentage_cutoff",
            new_value = 0
          )
          
        } else {
          # FLEXIBLE MODE: Use existing adaptive logic
          logger::log_info("Protein Processing: Using FLEXIBLE MODE (adaptive thresholds)")
          
          # First: Update missing value parameters (chunk 20) using NEW function signature
          current_s4 <- updateMissingValueParameters(
            theObject = current_s4,
            min_reps_per_group = input$min_reps_per_group,
            min_groups = input$min_groups
          )
        }
        
        # Always update the intensity cutoff percentile
        current_s4 <- updateConfigParameter(
          theObject = current_s4,
          function_name = "removeRowsWithMissingValuesPercent",
          parameter_name = "proteins_intensity_cutoff_percentile",
          new_value = input$proteins_intensity_cutoff_percentile
        )
        
        # Apply S4 transformation (EXISTING S4 CODE - UNCHANGED)
        filtered_s4 <- removeRowsWithMissingValuesPercent(current_s4)
        
        # Track QC parameters in workflow_data
        if (is.null(workflow_data$qc_params)) {
          workflow_data$qc_params <- list()
        }
        if (is.null(workflow_data$qc_params$protein_qc)) {
          workflow_data$qc_params$protein_qc <- list()
        }
        
        workflow_data$qc_params$protein_qc$intensity_filter <- list(
          strict_mode = use_strict_mode,
          min_reps_per_group = if(!use_strict_mode) input$min_reps_per_group else NA,
          min_groups = if(!use_strict_mode) input$min_groups else NA,
          groupwise_percentage_cutoff = current_s4@args$removeRowsWithMissingValuesPercent$groupwise_percentage_cutoff,
          max_groups_percentage_cutoff = current_s4@args$removeRowsWithMissingValuesPercent$max_groups_percentage_cutoff,
          proteins_intensity_cutoff_percentile = input$proteins_intensity_cutoff_percentile,
          timestamp = Sys.time()
        )
        
        # Save new state
        workflow_data$state_manager$saveState(
          state_name = "protein_intensity_filtered",
          s4_data_object = filtered_s4,
          config_object = list(
            strict_mode = use_strict_mode,
            min_reps_per_group = if(!use_strict_mode) input$min_reps_per_group else NA,
            min_groups = if(!use_strict_mode) input$min_groups else NA,
            groupwise_percentage_cutoff = current_s4@args$removeRowsWithMissingValuesPercent$groupwise_percentage_cutoff,
            max_groups_percentage_cutoff = current_s4@args$removeRowsWithMissingValuesPercent$max_groups_percentage_cutoff,
            proteins_intensity_cutoff_percentile = input$proteins_intensity_cutoff_percentile
          ),
          description = if(use_strict_mode) "Applied STRICT protein intensity filter (no missing values)" else "Applied FLEXIBLE protein intensity filter (adaptive thresholds)"
        )
        
        # Generate summary
        protein_count <- filtered_s4@protein_quant_table |>
          dplyr::distinct(Protein.Ids) |>
          nrow()
        
        result_text <- if (use_strict_mode) {
          paste(
            "Protein Intensity Filter Applied Successfully\n",
            "============================================\n",
            "Mode: STRICT (No Missing Values)\n",
            sprintf("Proteins remaining: %d\n", protein_count),
            "Groupwise % cutoff: 0.000% (strict - no missing allowed)\n",
            "Max groups % cutoff: 0.000% (strict - all groups must pass)\n",
            sprintf("Intensity cutoff percentile: %.1f%%\n", input$proteins_intensity_cutoff_percentile),
            "State saved as: 'protein_intensity_filtered'\n"
          )
        } else {
          paste(
            "Protein Intensity Filter Applied Successfully\n",
            "============================================\n",
            "Mode: FLEXIBLE (Adaptive Thresholds)\n",
            sprintf("Proteins remaining: %d\n", protein_count),
            sprintf("Min replicates per group: %d\n", input$min_reps_per_group),
            sprintf("Min groups required: %d\n", input$min_groups),
            sprintf("Groupwise %% cutoff: %.3f%% (calculated)\n", current_s4@args$removeRowsWithMissingValuesPercent$groupwise_percentage_cutoff),
            sprintf("Max groups %% cutoff: %.3f%% (calculated)\n", current_s4@args$removeRowsWithMissingValuesPercent$max_groups_percentage_cutoff),
            sprintf("Intensity cutoff percentile: %.1f%%\n", input$proteins_intensity_cutoff_percentile),
            "State saved as: 'protein_intensity_filtered'\n"
          )
        }
        
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
    shiny::observeEvent(input$revert_protein_intensity_filter, {
      tryCatch({
        # Revert to the previous state in the history
        history <- workflow_data$state_manager$getHistory()
        if (length(history) > 1) {
          prev_state_name <- history[length(history) - 1]
          reverted_s4 <- workflow_data$state_manager$revertToState(prev_state_name)
          output$protein_intensity_filter_results <- shiny::renderText(paste("Reverted to previous state:", prev_state_name))
          logger::log_info(paste("Reverted protein intensity filter to", prev_state_name))
          shiny::showNotification("Reverted successfully", type = "message")
        } else {
          stop("No previous state to revert to.")
        }
        
      }, error = function(e) {
        msg <- paste("Error reverting:", e$message)
        logger::log_error(msg)
        shiny::showNotification(msg, type = "error")
      })
    })
    
    # Render protein intensity filter plot
    output$protein_intensity_filter_plot <- shiny::renderPlot({
      shiny::req(protein_intensity_filter_plot())
      grid::grid.draw(protein_intensity_filter_plot())
    })
  })
}

