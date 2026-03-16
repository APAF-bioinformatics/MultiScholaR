# MultiScholaR: Interactive Multi-Omics Analysis
# Copyright (C) 2024-2026 Ignatius Pang, William Klare, and APAF-bioinformatics
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

#' @title Peptide Intensity Filter Module
#'
#' @description A Shiny module for applying peptide intensity filtering.
#'
#' @name mod_prot_qc_peptide_intensity
NULL

#' @rdname mod_prot_qc_peptide_intensity
#' @export
#' @importFrom shiny NS tagList tabPanel br fluidRow column wellPanel h4 p hr numericInput helpText div actionButton verbatimTextOutput plotOutput
mod_prot_qc_peptide_intensity_ui <- function(id) {
  ns <- shiny::NS(id)
  
  shiny::tabPanel(
    "Intensity Filter",
    shiny::br(),
    shiny::fluidRow(
      shiny::column(4,
        shiny::wellPanel(
          shiny::h4("Missing Value Parameters & Intensity Filter"),
          shiny::p("Configure missing value thresholds and filter peptides on intensity and missing values."),
          shiny::hr(),
          
          # Strict Mode Toggle
          shiny::checkboxInput(ns("use_strict_mode"), 
            "Strict Mode (No Missing Values Allowed)", 
            value = FALSE,
            width = "100%"
          ),
          shiny::helpText("When enabled, removes any peptide with missing values in ANY sample across all groups. Overrides flexible threshold settings below."),
          
          shiny::hr(),
          
          # Flexible mode parameters
          shiny::conditionalPanel(
            condition = paste0("!input['", ns("use_strict_mode"), "']"),
            
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
            shiny::p(style = "background-color: #f0f0f0; padding: 10px; border-radius: 5px; color: #333;",
              shiny::strong("These values are automatically calculated from your replicate settings above:"),
              shiny::br(),
              shiny::textOutput(ns("calculated_groupwise_percent"), inline = TRUE),
              shiny::br(),
              shiny::textOutput(ns("calculated_max_groups_percent"), inline = TRUE)
            ),
            
            shiny::hr()
          ),
          
          # Intensity cutoff
          shiny::h5("Intensity Threshold"),
          shiny::numericInput(ns("intensity_cutoff_percentile"), 
            "Peptide Intensity Cutoff Percentile (%)", 
            value = 1, min = 0.1, max = 10, step = 0.1,
            width = "100%"
          ),
          shiny::helpText("Intensity threshold percentile (default: 1%)"),
          
          shiny::hr(),
          shiny::div(
            shiny::actionButton(ns("apply_intensity_filter"), "Apply Filter", 
              class = "btn-primary", width = "48%"),
            shiny::actionButton(ns("revert_intensity"), "Revert", 
              class = "btn-warning", width = "48%", style = "margin-left: 4%")
          )
        )
      ),
      shiny::column(8,
        shiny::verbatimTextOutput(ns("intensity_results")),
        shiny::br(),
        shinyjqui::jqui_resizable(
          shiny::plotOutput(ns("intensity_plot"), height = "800px", width = "100%")
        )
      )
    )
  )
}

#' @rdname mod_prot_qc_peptide_intensity
#' @export
#' @importFrom shiny moduleServer reactiveVal observeEvent req showNotification removeNotification renderText renderPlot
#' @importFrom logger log_info log_error
#' @importFrom grid grid.draw
mod_prot_qc_peptide_intensity_server <- function(id, workflow_data, omic_type, experiment_label) {
  shiny::moduleServer(id, function(input, output, session) {
    
    intensity_plot <- shiny::reactiveVal(NULL)
    ns <- session$ns
    
    # Calculate percentages for UI display
    output$calculated_groupwise_percent <- shiny::renderText({
      shiny::req(workflow_data$state_manager)
      current_s4 <- workflow_data$state_manager$getState()
      shiny::req(current_s4)
      
      # Use the utility to calculate temporarily
      temp_s4 <- updateMissingValueParameters(
        theObject = current_s4,
        min_reps_per_group = input$min_reps_per_group,
        min_groups = input$min_groups,
        function_name = "peptideIntensityFiltering",
        grouping_variable = "group" # Default for DIA-NN
      )
      
      sprintf("Groupwise %% cutoff: %.3f%%", temp_s4@args$peptideIntensityFiltering$groupwise_percentage_cutoff)
    })
    
    output$calculated_max_groups_percent <- shiny::renderText({
      shiny::req(workflow_data$state_manager)
      current_s4 <- workflow_data$state_manager$getState()
      shiny::req(current_s4)
      
      temp_s4 <- updateMissingValueParameters(
        theObject = current_s4,
        min_reps_per_group = input$min_reps_per_group,
        min_groups = input$min_groups,
        function_name = "peptideIntensityFiltering",
        grouping_variable = "group"
      )
      
      sprintf("Max groups %% cutoff: %.3f%%", temp_s4@args$peptideIntensityFiltering$max_groups_percentage_cutoff)
    })

    # Step: Intensity Filter
    shiny::observeEvent(input$apply_intensity_filter, {
      shiny::req(workflow_data$state_manager)
      
      shiny::showNotification("Applying intensity filter...", id = "intensity_working", duration = NULL)
      
      tryCatch({
        # Get current S4 object
        current_s4 <- workflow_data$state_manager$getState()
        shiny::req(current_s4)
        
        # Mode Handling
        use_strict_mode <- isTRUE(input$use_strict_mode)
        
        if (use_strict_mode) {
          logger::log_info("Peptide Processing: Using STRICT MODE")
          
          current_s4 <- updateConfigParameter(
            theObject = current_s4,
            function_name = "peptideIntensityFiltering",
            parameter_name = "groupwise_percentage_cutoff",
            new_value = 0
          )
          current_s4 <- updateConfigParameter(
            theObject = current_s4,
            function_name = "peptideIntensityFiltering",
            parameter_name = "max_groups_percentage_cutoff",
            new_value = 0
          )
        } else {
          logger::log_info("Peptide Processing: Using FLEXIBLE MODE")
          
          current_s4 <- updateMissingValueParameters(
            theObject = current_s4,
            min_reps_per_group = input$min_reps_per_group,
            min_groups = input$min_groups,
            function_name = "peptideIntensityFiltering",
            grouping_variable = "group"
          )
        }
        
        # Always update intensity cutoff
        current_s4 <- updateConfigParameter(
          theObject = current_s4,
          function_name = "peptideIntensityFiltering",
          parameter_name = "peptides_intensity_cutoff_percentile",
          new_value = input$intensity_cutoff_percentile
        )
        
        # Apply S4 transformation
        filtered_s4 <- peptideIntensityFiltering(theObject = current_s4)
        
        # Track QC parameters
        if (is.null(workflow_data$qc_params)) workflow_data$qc_params <- list()
        if (is.null(workflow_data$qc_params$peptide_qc)) workflow_data$qc_params$peptide_qc <- list()
        
        workflow_data$qc_params$peptide_qc$intensity_filter <- list(
          strict_mode = use_strict_mode,
          min_reps_per_group = if(!use_strict_mode) input$min_reps_per_group else NA,
          min_groups = if(!use_strict_mode) input$min_groups else NA,
          intensity_cutoff_percentile = input$intensity_cutoff_percentile,
          timestamp = Sys.time()
        )
        
        # Save state
        workflow_data$state_manager$saveState(
          state_name = "intensity_filtered",
          s4_data_object = filtered_s4,
          config_object = list(
            strict_mode = use_strict_mode,
            intensity_cutoff_percentile = input$intensity_cutoff_percentile
          ),
          description = if(use_strict_mode) "Applied STRICT peptide intensity filter" else "Applied FLEXIBLE peptide intensity filter"
        )
        
        # Summary
        protein_count <- filtered_s4@peptide_data |> dplyr::distinct(Protein.Ids) |> nrow()
        
        result_text <- paste(
          "Intensity Filter Applied Successfully\n",
          "====================================\n",
          sprintf("Mode: %s\n", if(use_strict_mode) "STRICT" else "FLEXIBLE"),
          sprintf("Proteins remaining: %d\n", protein_count),
          sprintf("Intensity cutoff percentile: %.1f%%\n", input$intensity_cutoff_percentile),
          sprintf("Groupwise %% cutoff: %.3f%%\n", current_s4@args$peptideIntensityFiltering$groupwise_percentage_cutoff),
          sprintf("Max groups %% cutoff: %.3f%%\n", current_s4@args$peptideIntensityFiltering$max_groups_percentage_cutoff),
          "State saved as: 'intensity_filtered'\n"
        )
        
        output$intensity_results <- shiny::renderText(result_text)
        
        # Visualization
        plot_grid <- updateProteinFiltering(
          data = filtered_s4@peptide_data,
          step_name = "4_intensity_filtered",
          omic_type = omic_type,
          experiment_label = experiment_label,
          return_grid = TRUE,
          overwrite = TRUE
        )
        intensity_plot(plot_grid)
        
        shiny::removeNotification("intensity_working")
        shiny::showNotification("Intensity filter applied successfully", type = "message")
        
      }, error = function(e) {
        msg <- paste("Error applying intensity filter:", e$message)
        logger::log_error(msg)
        shiny::showNotification(msg, type = "error", duration = 15)
        shiny::removeNotification("intensity_working")
      })
    })
    
    # Revert
    shiny::observeEvent(input$revert_intensity, {
      tryCatch({
        history <- workflow_data$state_manager$getHistory()
        if (length(history) > 1) {
          prev_state_name <- history[length(history) - 1]
          reverted_s4 <- workflow_data$state_manager$revertToState(prev_state_name)
          output$intensity_results <- shiny::renderText(paste("Reverted to previous state:", prev_state_name))
          shiny::showNotification("Reverted successfully", type = "message")
        } else {
          stop("No previous state to revert to.")
        }
      }, error = function(e) {
        shiny::showNotification(paste("Error reverting:", e$message), type = "error")
      })
    })
    
    # Render plot
    output$intensity_plot <- shiny::renderPlot({
      shiny::req(intensity_plot())
      grid::grid.draw(intensity_plot())
    })
  })
}

