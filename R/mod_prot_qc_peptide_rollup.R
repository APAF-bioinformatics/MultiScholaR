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

#' @title Peptide Precursor Rollup Module
#'
#' @description A Shiny module for applying precursor to peptide rollup.
#'
#' @name mod_prot_qc_peptide_rollup
NULL

#' @rdname mod_prot_qc_peptide_rollup
#' @export
#' @importFrom shiny NS tagList tabPanel br fluidRow column wellPanel h4 p hr div actionButton verbatimTextOutput plotOutput
mod_prot_qc_peptide_rollup_ui <- function(id) {
  ns <- shiny::NS(id)
  
  shiny::tabPanel(
    "Precursor Rollup",
    shiny::br(),
    shiny::fluidRow(
      shiny::column(4,
        shiny::wellPanel(
          shiny::h4("Precursor to Peptide Rollup"),
          shiny::p("Aggregate intensity measurements from multiple precursor ions to peptide level."),
          shiny::hr(),
          shiny::p("This step has no configurable parameters - it uses statistical methods to combine precursor measurements."),
          shiny::hr(),
          shiny::div(
            shiny::actionButton(ns("apply_rollup"), "Apply Rollup", 
              class = "btn-primary", width = "48%"),
            shiny::actionButton(ns("revert_rollup"), "Revert", 
              class = "btn-warning", width = "48%", style = "margin-left: 4%")
          )
        )
      ),
      shiny::column(8,
        shiny::verbatimTextOutput(ns("rollup_results")),
        shiny::br(),
        shinyjqui::jqui_resizable(
          shiny::plotOutput(ns("rollup_plot"), height = "800px", width = "100%")
        )
      )
    )
  )
}

#' @rdname mod_prot_qc_peptide_rollup
#' @export
#' @importFrom shiny moduleServer reactiveVal observeEvent req showNotification removeNotification renderText renderPlot
#' @importFrom logger log_info log_error
#' @importFrom grid grid.draw
runPeptideRollupApplyStep <- function(workflowData,
                                      rollupFn = rollUpPrecursorToPeptide,
                                      logInfoFn = logger::log_info,
                                      nowFn = Sys.time) {
  shiny::req(workflowData$state_manager)
  currentS4 <- workflowData$state_manager$getState()
  shiny::req(currentS4)

  logInfoFn("QC Step: Applying precursor rollup")

  rolledUpS4 <- rollupFn(currentS4)

  if (is.null(workflowData$qc_params)) {
    workflowData$qc_params <- list()
  }
  if (is.null(workflowData$qc_params$peptide_qc)) {
    workflowData$qc_params$peptide_qc <- list()
  }

  workflowData$qc_params$peptide_qc$precursor_rollup <- list(
    applied = TRUE,
    timestamp = nowFn()
  )

  workflowData$state_manager$saveState(
    state_name = "precursor_rollup",
    s4_data_object = rolledUpS4,
    config_object = list(),
    description = "Applied precursor to peptide rollup"
  )

  proteinCount <- rolledUpS4@peptide_data |>
    dplyr::distinct(Protein.Ids) |>
    nrow()

  resultText <- paste(
    "Precursor Rollup Applied Successfully\n",
    "====================================\n",
    sprintf("Proteins remaining: %d\n", proteinCount),
    "State saved as: 'precursor_rollup'\n"
  )

  list(
    rolledUpS4 = rolledUpS4,
    resultText = resultText
  )
}

updatePeptideRollupOutputs <- function(output,
                                       rollupPlot,
                                       rollupResult,
                                       omicType,
                                       experimentLabel,
                                       renderTextFn = shiny::renderText,
                                       updateProteinFilteringFn = updateProteinFiltering) {
  output$rollup_results <- renderTextFn(rollupResult$resultText)

  plotGrid <- updateProteinFilteringFn(
    data = rollupResult$rolledUpS4@peptide_data,
    step_name = "3_precursor_rollup",
    omic_type = omicType,
    experiment_label = experimentLabel,
    return_grid = TRUE,
    overwrite = TRUE
  )
  rollupPlot(plotGrid)

  invisible(plotGrid)
}

runPeptideRollupApplyObserver <- function(workflowData,
                                          output,
                                          rollupPlot,
                                          omicType,
                                          experimentLabel,
                                          runApplyStepFn = runPeptideRollupApplyStep,
                                          updateOutputsFn = updatePeptideRollupOutputs,
                                          showNotificationFn = shiny::showNotification,
                                          removeNotificationFn = shiny::removeNotification,
                                          logInfoFn = logger::log_info,
                                          logErrorFn = logger::log_error) {
  showNotificationFn(
    "Applying precursor rollup...",
    id = "rollup_working",
    duration = NULL
  )

  tryCatch({
    rollupResult <- runApplyStepFn(workflowData = workflowData)

    plotGrid <- updateOutputsFn(
      output = output,
      rollupPlot = rollupPlot,
      rollupResult = rollupResult,
      omicType = omicType,
      experimentLabel = experimentLabel
    )

    logInfoFn("Precursor rollup applied successfully")
    removeNotificationFn("rollup_working")
    showNotificationFn("Precursor rollup applied successfully", type = "message")

    list(
      status = "success",
      rollupResult = rollupResult,
      plotGrid = plotGrid
    )
  }, error = function(e) {
    errorMessage <- paste("Error applying precursor rollup:", e$message)
    logErrorFn(errorMessage)
    showNotificationFn(errorMessage, type = "error", duration = 15)
    removeNotificationFn("rollup_working")

    list(
      status = "error",
      errorMessage = errorMessage
    )
  })
}

runPeptideRollupRevertStep <- function(workflowData,
                                       logInfoFn = logger::log_info) {
  history <- workflowData$state_manager$getHistory()
  if (length(history) <= 1) {
    stop("No previous state to revert to.")
  }

  previousState <- history[length(history) - 1]
  revertedS4 <- workflowData$state_manager$revertToState(previousState)
  logInfoFn(paste("Reverted precursor rollup to", previousState))

  list(
    previousState = previousState,
    revertedS4 = revertedS4,
    resultText = paste("Reverted to previous state:", previousState)
  )
}

runPeptideRollupRevertObserver <- function(workflowData,
                                           output,
                                           runRevertStepFn = runPeptideRollupRevertStep,
                                           renderTextFn = shiny::renderText,
                                           showNotificationFn = shiny::showNotification,
                                           logErrorFn = logger::log_error) {
  tryCatch({
    revertResult <- runRevertStepFn(workflowData = workflowData)
    output$rollup_results <- renderTextFn(revertResult$resultText)
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

mod_prot_qc_peptide_rollup_server <- function(id, workflow_data, omic_type, experiment_label) {
  shiny::moduleServer(id, function(input, output, session) {
    
    rollup_plot <- shiny::reactiveVal(NULL)
    
    # Step 2: Precursor Rollup (chunk 11)
    shiny::observeEvent(input$apply_rollup, {
      runPeptideRollupApplyObserver(
        workflowData = workflow_data,
        output = output,
        rollupPlot = rollup_plot,
        omicType = omic_type,
        experimentLabel = experiment_label
      )
    })
    
    # Revert Precursor Rollup
    shiny::observeEvent(input$revert_rollup, {
      runPeptideRollupRevertObserver(
        workflowData = workflow_data,
        output = output
      )
    })
    
    # Render rollup plot
    output$rollup_plot <- shiny::renderPlot({
      shiny::req(rollup_plot())
      grid::grid.draw(rollup_plot())
    })
  })
}
