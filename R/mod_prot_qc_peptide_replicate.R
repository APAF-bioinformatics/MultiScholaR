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

#' @title Peptide Replicate Filter Module
#'
#' @description A Shiny module for applying the replicate filter.
#'
#' @name mod_prot_qc_peptide_replicate
NULL

#' @rdname mod_prot_qc_peptide_replicate
#' @export
#' @importFrom shiny NS tagList tabPanel br fluidRow column wellPanel h4 p hr textInput helpText div actionButton verbatimTextOutput plotOutput
mod_prot_qc_peptide_replicate_ui <- function(id) {
  ns <- shiny::NS(id)
  
  shiny::tabPanel(
    "Replicate Filter",
    shiny::br(),
    shiny::fluidRow(
      shiny::column(4,
        shiny::wellPanel(
          shiny::h4("Remove Single-Replicate Peptides"),
          shiny::p("Remove peptides that appear in only one replicate across all groups."),
          shiny::hr(),
          
          shiny::textInput(ns("replicate_group_column"), 
            "Replicate Group Column", 
            value = "replicates",
            width = "100%"
          ),
          shiny::helpText("Column name for grouping replicates (default: 'replicates')"),
          
          shiny::hr(),
          shiny::div(
            shiny::actionButton(ns("apply_replicate_filter"), "Apply Filter", 
              class = "btn-primary", width = "48%"),
            shiny::actionButton(ns("revert_replicate"), "Revert", 
              class = "btn-warning", width = "48%", style = "margin-left: 4%")
          )
        )
      ),
      shiny::column(8,
        shiny::verbatimTextOutput(ns("replicate_results")),
        shiny::br(),
        shinyjqui::jqui_resizable(
          shiny::plotOutput(ns("replicate_plot"), height = "800px", width = "100%")
        )
      )
    )
  )
}

#' @rdname mod_prot_qc_peptide_replicate
#' @export
#' @importFrom shiny moduleServer reactiveVal observeEvent req showNotification removeNotification renderText renderPlot
#' @importFrom logger log_info log_error
#' @importFrom grid grid.draw
runPeptideReplicateRevertStep <- function(workflowData,
                                          logInfoFn = logger::log_info) {
  history <- workflowData$state_manager$getHistory()
  if (length(history) <= 1) {
    stop("No previous state to revert to.")
  }

  previousState <- history[length(history) - 1]
  revertedS4 <- workflowData$state_manager$revertToState(previousState)
  logInfoFn(paste("Reverted replicate filter to", previousState))

  list(
    previousState = previousState,
    revertedS4 = revertedS4,
    resultText = paste("Reverted to previous state:", previousState)
  )
}

runPeptideReplicateRevertObserver <- function(workflowData,
                                              output,
                                              runRevertStepFn = runPeptideReplicateRevertStep,
                                              renderTextFn = shiny::renderText,
                                              showNotificationFn = shiny::showNotification,
                                              logErrorFn = logger::log_error) {
  tryCatch({
    revertResult <- runRevertStepFn(workflowData = workflowData)
    output$replicate_results <- renderTextFn(revertResult$resultText)
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

runPeptideReplicateApplyStep <- function(workflowData,
                                         replicateGroupColumn,
                                         removePeptidesWithOnlyOneReplicateFn = removePeptidesWithOnlyOneReplicate,
                                         logInfoFn = logger::log_info,
                                         nowFn = Sys.time) {
  shiny::req(workflowData$state_manager)
  currentS4 <- workflowData$state_manager$getState()
  shiny::req(currentS4)

  logInfoFn(sprintf(
    "QC Step: Applying replicate filter (column: %s)",
    replicateGroupColumn
  ))

  filteredS4 <- removePeptidesWithOnlyOneReplicateFn(
    currentS4,
    replicate_group_column = replicateGroupColumn
  )

  if (is.null(workflowData$qc_params)) {
    workflowData$qc_params <- list()
  }
  if (is.null(workflowData$qc_params$peptide_qc)) {
    workflowData$qc_params$peptide_qc <- list()
  }

  workflowData$qc_params$peptide_qc$replicate_filter <- list(
    replicate_group_column = replicateGroupColumn,
    timestamp = nowFn()
  )

  workflowData$state_manager$saveState(
    state_name = "replicate_filtered",
    s4_data_object = filteredS4,
    config_object = list(
      replicate_group_column = replicateGroupColumn
    ),
    description = "Applied replicate filter (removed single-replicate peptides)"
  )

  proteinCount <- filteredS4@peptide_data |>
    dplyr::distinct(Protein.Ids) |>
    nrow()

  resultText <- paste(
    "Replicate Filter Applied Successfully\n",
    "====================================\n",
    sprintf("Proteins remaining: %d\n", proteinCount),
    sprintf("Replicate group column: %s\n", replicateGroupColumn),
    "State saved as: 'replicate_filtered'\n"
  )

  list(
    filteredS4 = filteredS4,
    resultText = resultText
  )
}

updatePeptideReplicateOutputs <- function(output,
                                          replicatePlot,
                                          replicateResult,
                                          omicType,
                                          experimentLabel,
                                          renderTextFn = shiny::renderText,
                                          updateProteinFilteringFn = updateProteinFiltering) {
  output$replicate_results <- renderTextFn(replicateResult$resultText)

  plotGrid <- updateProteinFilteringFn(
    data = replicateResult$filteredS4@peptide_data,
    step_name = "7_replicate_filtered",
    omic_type = omicType,
    experiment_label = experimentLabel,
    return_grid = TRUE,
    overwrite = TRUE
  )
  replicatePlot(plotGrid)

  invisible(plotGrid)
}

runPeptideReplicateApplyObserver <- function(workflowData,
                                             replicateGroupColumn,
                                             output,
                                             replicatePlot,
                                             omicType,
                                             experimentLabel,
                                             runApplyStepFn = runPeptideReplicateApplyStep,
                                             updateOutputsFn = updatePeptideReplicateOutputs,
                                             showNotificationFn = shiny::showNotification,
                                             removeNotificationFn = shiny::removeNotification,
                                             logInfoFn = logger::log_info,
                                             logErrorFn = logger::log_error) {
  showNotificationFn(
    "Applying replicate filter...",
    id = "replicate_working",
    duration = NULL
  )

  tryCatch({
    replicateResult <- runApplyStepFn(
      workflowData = workflowData,
      replicateGroupColumn = replicateGroupColumn
    )

    plotGrid <- updateOutputsFn(
      output = output,
      replicatePlot = replicatePlot,
      replicateResult = replicateResult,
      omicType = omicType,
      experimentLabel = experimentLabel
    )

    logInfoFn("Replicate filter applied successfully")
    removeNotificationFn("replicate_working")
    showNotificationFn("Replicate filter applied successfully", type = "message")

    list(
      status = "success",
      replicateResult = replicateResult,
      plotGrid = plotGrid
    )
  }, error = function(e) {
    errorMessage <- paste("Error applying replicate filter:", e$message)
    logErrorFn(errorMessage)
    showNotificationFn(errorMessage, type = "error", duration = 15)
    removeNotificationFn("replicate_working")

    list(
      status = "error",
      errorMessage = errorMessage
    )
  })
}

mod_prot_qc_peptide_replicate_server <- function(id, workflow_data, omic_type, experiment_label) {
  shiny::moduleServer(id, function(input, output, session) {
    
    replicate_plot <- shiny::reactiveVal(NULL)
    
    # Step 6: Replicate Filter (chunk 15)
    shiny::observeEvent(input$apply_replicate_filter, {
      runPeptideReplicateApplyObserver(
        workflowData = workflow_data,
        replicateGroupColumn = input$replicate_group_column,
        output = output,
        replicatePlot = replicate_plot,
        omicType = omic_type,
        experimentLabel = experiment_label
      )
    })
    
    # Revert Replicate Filter
    shiny::observeEvent(input$revert_replicate, {
      runPeptideReplicateRevertObserver(
        workflowData = workflow_data,
        output = output
      )
    })
    
    # Render replicate filter plot
    output$replicate_plot <- shiny::renderPlot({
      shiny::req(replicate_plot())
      grid::grid.draw(replicate_plot())
    })
  })
}
