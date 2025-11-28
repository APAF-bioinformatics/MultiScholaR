#' @title Peptide Imputation Module
#'
#' @description A Shiny module for applying missing value imputation.
#'
#' @name mod_prot_qc_peptide_impute
NULL

#' @rdname mod_prot_qc_peptide_impute
#' @export
#' @importFrom shiny NS tagList tabPanel br fluidRow column wellPanel h4 p hr numericInput helpText div actionButton verbatimTextOutput plotOutput
mod_prot_qc_peptide_impute_ui <- function(id) {
  ns <- shiny::NS(id)
  
  shiny::tabPanel(
    "Imputation",
    shiny::br(),
    shiny::fluidRow(
      shiny::column(4,
        shiny::wellPanel(
          shiny::h4("Missing Value Imputation"),
          shiny::p("Impute missing values using technical replicate averages."),
          shiny::hr(),
          
          shiny::numericInput(ns("proportion_missing_values"), 
            "Max Proportion Missing", 
            value = 0.5, min = 0.1, max = 1.0, step = 0.1,
            width = "100%"
          ),
          shiny::helpText("Lower = more stringent (less imputation), higher = less stringent (more imputation) (default: 0.5)"),
          
          shiny::hr(),
          shiny::div(
            shiny::actionButton(ns("apply_imputation"), "Apply Imputation", 
              class = "btn-primary", width = "48%"),
            shiny::actionButton(ns("revert_imputation"), "Revert", 
              class = "btn-warning", width = "48%", style = "margin-left: 4%")
          )
        )
      ),
      shiny::column(8,
        shiny::verbatimTextOutput(ns("imputation_results")),
        shiny::br(),
        shinyjqui::jqui_resizable(
          shiny::plotOutput(ns("imputation_plot"), height = "800px", width = "100%")
        )
      )
    )
  )
}

#' @rdname mod_prot_qc_peptide_impute
#' @export
#' @importFrom shiny moduleServer reactiveVal observeEvent req showNotification removeNotification renderText renderPlot
#' @importFrom logger log_info log_error
#' @importFrom grid grid.draw
mod_prot_qc_peptide_impute_server <- function(id, workflow_data, omic_type, experiment_label) {
  shiny::moduleServer(id, function(input, output, session) {
    
    imputation_plot <- shiny::reactiveVal(NULL)
    
    # Step 7: Missing Value Imputation (chunk 16)
    shiny::observeEvent(input$apply_imputation, {
      shiny::req(workflow_data$state_manager)
      
      shiny::showNotification("Applying missing value imputation...", id = "imputation_working", duration = NULL)
      
      tryCatch({
        # Get current S4 object from the active state
        current_s4 <- workflow_data$state_manager$getState()
        shiny::req(current_s4)
        
        logger::log_info(sprintf("QC Step: Applying missing value imputation (proportion: %s)", input$proportion_missing_values))
        
        # ✅ FIXED: Use updateConfigParameter to sync S4 object AND global config_list
        current_s4 <- updateConfigParameter(
          theObject = current_s4,
          function_name = "peptideMissingValueImputation",
          parameter_name = "proportion_missing_values",
          new_value = input$proportion_missing_values
        )
        
        # Apply S4 transformation (EXISTING S4 CODE - UNCHANGED)
        imputed_s4 <- peptideMissingValueImputation(theObject = current_s4)
        
        # Track QC parameters in workflow_data
        if (is.null(workflow_data$qc_params)) {
          workflow_data$qc_params <- list()
        }
        if (is.null(workflow_data$qc_params$peptide_qc)) {
          workflow_data$qc_params$peptide_qc <- list()
        }
        
        workflow_data$qc_params$peptide_qc$imputation <- list(
          proportion_missing_values = input$proportion_missing_values,
          timestamp = Sys.time()
        )
        
        # Save QC parameters to file for persistence (final peptide QC step)
        tryCatch({
          if (exists("experiment_paths") && !is.null(experiment_paths$source_dir)) {
            qc_params_file <- file.path(experiment_paths$source_dir, "qc_params.RDS")
            saveRDS(workflow_data$qc_params, qc_params_file)
            logger::log_info(sprintf("Saved QC parameters to: %s", qc_params_file))
          }
        }, error = function(e) {
          logger::log_warn(sprintf("Could not save QC parameters file: %s", e$message))
        })
        
        # Save new state in R6 manager
        workflow_data$state_manager$saveState(
          state_name = "imputed",
          s4_data_object = imputed_s4,
          config_object = list(
            proportion_missing_values = input$proportion_missing_values
          ),
          description = "Applied missing value imputation using technical replicates"
        )
        
        # Generate summary
        protein_count <- imputed_s4@peptide_data |>
          dplyr::distinct(Protein.Ids) |>
          nrow()
        
        result_text <- paste(
          "Missing Value Imputation Applied Successfully\n",
          "============================================\n",
          sprintf("Proteins remaining: %d\n", protein_count),
          sprintf("Max proportion missing: %.1f\n", input$proportion_missing_values),
          "State saved as: 'imputed'\n",
          "\nReady for peptide-to-protein rollup step."
        )
        
        output$imputation_results <- shiny::renderText(result_text)
        
        # Update filtering visualization and capture plot
        plot_grid <- updateProteinFiltering(
          data = imputed_s4@peptide_data,
          step_name = "8_imputed",
          omic_type = omic_type,
          experiment_label = experiment_label,
          return_grid = TRUE,
          overwrite = TRUE
        )
        imputation_plot(plot_grid)
        
        logger::log_info("Missing value imputation applied successfully")
        shiny::removeNotification("imputation_working")
        shiny::showNotification("Missing value imputation applied successfully", type = "message")
        
      }, error = function(e) {
        msg <- paste("Error applying missing value imputation:", e$message)
        logger::log_error(msg)
        shiny::showNotification(msg, type = "error", duration = 15)
        shiny::removeNotification("imputation_working")
      })
    })
    
    # Revert Imputation
    shiny::observeEvent(input$revert_imputation, {
      tryCatch({
        # Revert to the previous state in the history
        history <- workflow_data$state_manager$getHistory()
        if (length(history) > 1) {
          prev_state_name <- history[length(history) - 1]
          reverted_s4 <- workflow_data$state_manager$revertToState(prev_state_name)
          output$imputation_results <- shiny::renderText(paste("Reverted to previous state:", prev_state_name))
          logger::log_info(paste("Reverted imputation to", prev_state_name))
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
    
    # Render imputation plot
    output$imputation_plot <- shiny::renderPlot({
      shiny::req(imputation_plot())
      grid::grid.draw(imputation_plot())
    })
  })
}

