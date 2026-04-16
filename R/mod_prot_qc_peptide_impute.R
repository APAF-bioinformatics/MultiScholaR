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
runPeptideImputationStep <- function(workflow_data, proportionMissingValues) {
  current_s4 <- workflow_data$state_manager$getState()
  shiny::req(current_s4)

  logger::log_info(sprintf(
    "QC Step: Applying missing value imputation (proportion: %s)",
    proportionMissingValues
  ))

  current_s4 <- updateConfigParameter(
    theObject = current_s4,
    function_name = "peptideMissingValueImputation",
    parameter_name = "proportion_missing_values",
    new_value = proportionMissingValues
  )

  imputed_s4 <- peptideMissingValueImputation(theObject = current_s4)

  if (is.null(workflow_data$qc_params)) {
    workflow_data$qc_params <- list()
  }
  if (is.null(workflow_data$qc_params$peptide_qc)) {
    workflow_data$qc_params$peptide_qc <- list()
  }

  workflow_data$qc_params$peptide_qc$imputation <- list(
    proportion_missing_values = proportionMissingValues,
    timestamp = Sys.time()
  )

  tryCatch({
    if (exists("experiment_paths") && !is.null(experiment_paths$source_dir)) {
      qc_params_file <- file.path(experiment_paths$source_dir, "qc_params.RDS")
      saveRDS(workflow_data$qc_params, qc_params_file)
      logger::log_info(sprintf("Saved QC parameters to: %s", qc_params_file))
    }
  }, error = function(e) {
    logger::log_warn(sprintf("Could not save QC parameters file: %s", e$message))
  })

  workflow_data$state_manager$saveState(
    state_name = "imputed",
    s4_data_object = imputed_s4,
    config_object = list(
      proportion_missing_values = proportionMissingValues
    ),
    description = "Applied missing value imputation using technical replicates"
  )

  # --- TESTTHAT CHECKPOINT CP02 (see test-prot-02-qc-filtering.R) ---
  .capture_checkpoint(imputed_s4, "cp02", "qc_filtered_peptide")
  # --- END CP02 ---

  protein_count <- imputed_s4@peptide_data |>
    dplyr::distinct(Protein.Ids) |>
    nrow()

  result_text <- paste(
    "Missing Value Imputation Applied Successfully\n",
    "============================================\n",
    sprintf("Proteins remaining: %d\n", protein_count),
    sprintf("Max proportion missing: %.1f\n", proportionMissingValues),
    "State saved as: 'imputed'\n",
    "\nReady for peptide-to-protein rollup step."
  )

  list(
    imputedS4 = imputed_s4,
    resultText = result_text
  )
}

runPeptideImputationRevertStep <- function(workflow_data) {
  history <- workflow_data$state_manager$getHistory()
  if (length(history) <= 1) {
    stop("No previous state to revert to.")
  }

  previousState <- history[length(history) - 1]
  revertedS4 <- workflow_data$state_manager$revertToState(previousState)
  logger::log_info(paste("Reverted imputation to", previousState))

  list(
    previousState = previousState,
    revertedS4 = revertedS4,
    resultText = paste("Reverted to previous state:", previousState)
  )
}

updatePeptideImputationOutputs <- function(output, imputationPlot, imputationResult, omicType, experimentLabel) {
  output$imputation_results <- renderText(imputationResult$resultText)

  plotGrid <- updateProteinFiltering(
    data = imputationResult$imputedS4@peptide_data,
    step_name = "8_imputed",
    omic_type = omicType,
    experiment_label = experimentLabel,
    return_grid = TRUE,
    overwrite = TRUE
  )
  imputationPlot(plotGrid)

  invisible(plotGrid)
}

runPeptideImputationApplyObserver <- function(workflow_data,
                                              proportionMissingValues,
                                              output,
                                              imputationPlot,
                                              omicType,
                                              experimentLabel,
                                              runImputationStepFn = runPeptideImputationStep,
                                              updateOutputsFn = updatePeptideImputationOutputs,
                                              showNotificationFn = shiny::showNotification,
                                              removeNotificationFn = shiny::removeNotification,
                                              logInfoFn = logger::log_info,
                                              logErrorFn = logger::log_error) {
  showNotificationFn(
    "Applying missing value imputation...",
    id = "imputation_working",
    duration = NULL
  )

  tryCatch({
    imputationResult <- runImputationStepFn(
      workflow_data = workflow_data,
      proportionMissingValues = proportionMissingValues
    )

    plotGrid <- updateOutputsFn(
      output = output,
      imputationPlot = imputationPlot,
      imputationResult = imputationResult,
      omicType = omicType,
      experimentLabel = experimentLabel
    )

    logInfoFn("Missing value imputation applied successfully")
    removeNotificationFn("imputation_working")
    showNotificationFn(
      "Missing value imputation applied successfully",
      type = "message"
    )

    list(
      status = "success",
      imputationResult = imputationResult,
      plotGrid = plotGrid
    )
  }, error = function(e) {
    errorMessage <- paste("Error applying missing value imputation:", e$message)
    logErrorFn(errorMessage)
    showNotificationFn(errorMessage, type = "error", duration = 15)
    removeNotificationFn("imputation_working")

    list(
      status = "error",
      errorMessage = errorMessage
    )
  })
}

runPeptideImputationRevertObserver <- function(workflow_data,
                                               output,
                                               runRevertStepFn = runPeptideImputationRevertStep,
                                               renderTextFn = shiny::renderText,
                                               showNotificationFn = shiny::showNotification,
                                               logErrorFn = logger::log_error) {
  tryCatch({
    revertResult <- runRevertStepFn(workflow_data = workflow_data)
    output$imputation_results <- renderTextFn(revertResult$resultText)
    showNotificationFn("Reverted successfully", type = "message")

    list(
      status = "success",
      revertResult = revertResult
    )
  }, error = function(e) {
    errorMessage <- paste("Error reverting:", e$message)
    logErrorFn(errorMessage)
    showNotificationFn(errorMessage, type = "error")

    list(
      status = "error",
      errorMessage = errorMessage
    )
  })
}

bindPeptideImputationPlot <- function(output, imputationPlot) {
  output$imputation_plot <- renderPlot({
    req(imputationPlot())
    grid.draw(imputationPlot())
  })
}

mod_prot_qc_peptide_impute_server <- function(id, workflow_data, omic_type, experiment_label) {
  shiny::moduleServer(id, function(input, output, session) {
    
    imputation_plot <- shiny::reactiveVal(NULL)
    
    # Step 7: Missing Value Imputation (chunk 16)
    shiny::observeEvent(input$apply_imputation, {
      shiny::req(workflow_data$state_manager)

      runPeptideImputationApplyObserver(
        workflow_data = workflow_data,
        proportionMissingValues = input$proportion_missing_values,
        output = output,
        imputationPlot = imputation_plot,
        omicType = omic_type,
        experimentLabel = experiment_label
      )
    })
    
    # Revert Imputation
    shiny::observeEvent(input$revert_imputation, {
      runPeptideImputationRevertObserver(
        workflow_data = workflow_data,
        output = output
      )
    })
    
    bindPeptideImputationPlot(output = output, imputationPlot = imputation_plot)
  })
}
