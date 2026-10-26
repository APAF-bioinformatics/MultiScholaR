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

#' @title Log2 Peptide Matrix Module (DIA with limpa)
#'
#' @description The DIA_LIMPA route's log2 matrix step (KIT-047, the limpa script's chunk 13): the peptide
#' matrix is calculated from the rolled-up peptide table and log2 transformed (values <= 0 become NA). On an
#' on-disk owner the step runs the registered `peptide_log2_matrix` declaration, on the managed worker when the
#' session carries a managed application and in process otherwise.
#'
#' @name mod_prot_qc_peptide_log2
#' @param id,workflow_data,omic_type,experiment_label,experiment_paths Runtime inputs used by this function; see the usage section for accepted values.
NULL

#' @rdname mod_prot_qc_peptide_log2
#' @export
#' @importFrom shiny NS tabPanel br fluidRow column wellPanel h4 p hr div actionButton verbatimTextOutput plotOutput
mod_prot_qc_peptide_log2_ui <- function(id) {
  ns <- shiny::NS(id)

  shiny::tabPanel(
    "Log2 Matrix",
    shiny::br(),
    shiny::fluidRow(
      shiny::column(4,
        shiny::wellPanel(
          shiny::h4("Log2 Peptide Matrix"),
          shiny::p("Pivot the rolled-up peptide quantities into a peptide-by-sample matrix and log2 transform it. Values at or below zero become missing; nothing is imputed."),
          shiny::hr(),
          shiny::div(
            shiny::actionButton(ns("apply_log2_matrix"), "Apply Log2 Matrix", class = "btn-primary", width = "48%"),
            shiny::actionButton(ns("revert_log2_matrix"), "Revert", class = "btn-warning", width = "48%", style = "margin-left: 4%")
          )
        )
      ),
      shiny::column(8,
        shiny::verbatimTextOutput(ns("log2_results")),
        shiny::br(),
        shinyjqui::jqui_resizable(
          shiny::plotOutput(ns("log2_plot"), height = "800px", width = "100%")
        )
      )
    )
  )
}

runPeptideLog2ApplyStep <- function(workflowData,
                                    calcPeptideMatrixFn = calcPeptideMatrix,
                                    log2TransformPeptideMatrixFn = log2TransformPeptideMatrix,
                                    logInfoFn = logger::log_info,
                                    nowFn = Sys.time) {
runDIAReferencePeptideLog2Matrix(workflowData)
}

updatePeptideLog2Outputs <- function(output, log2Plot, log2Result, omicType, experimentLabel,
                                     renderTextFn = shiny::renderText,
                                     updateProteinFilteringFn = updateProteinFiltering, projectDirs = NULL) {
  output$log2_results <- renderPeptideQcResultText(log2Result$resultText, renderTextFn)
  plot_arguments <- list(
    data = log2Result$filteringSummary %||% log2Result$transformedS4@peptide_data,
    step_name = "4_log2_matrix", omic_type = omicType, experiment_label = experimentLabel,
    return_grid = TRUE, overwrite = TRUE)
  if (!is.null(projectDirs)) plot_arguments$project_dirs <- projectDirs
  plotGrid <- if (!is.null(log2Result$stateReference)) {
    publishDIAQvalueOutput(getDIAActiveManager(), log2Result, plot_arguments, updateProteinFilteringFn)
  } else do.call(updateProteinFilteringFn, plot_arguments)
  log2Plot(plotGrid)
  invisible(plotGrid)
}

runPeptideLog2ApplyObserver <- function(workflowData, output, log2Plot, omicType, experimentLabel,
                                        runApplyStepFn = runPeptideLog2ApplyStep,
                                        updateOutputsFn = updatePeptideLog2Outputs,
                                        showNotificationFn = shiny::showNotification,
                                        removeNotificationFn = shiny::removeNotification,
                                        logInfoFn = logger::log_info,
                                        logErrorFn = logger::log_error, experimentPaths = NULL) {
  withDIAStateScope(getDIAActiveManager(), function() {
  showNotificationFn("Calculating the log2 peptide matrix...", id = "log2_matrix_working", duration = NULL)
  tryCatch({
    log2Result <- runApplyStepFn(workflowData = workflowData)
    output_arguments <- list(output = output, log2Plot = log2Plot, log2Result = log2Result,
      omicType = omicType, experimentLabel = experimentLabel)
    if (!is.null(experimentPaths)) {
      output_arguments$projectDirs <- stats::setNames(list(experimentPaths), paste(omicType, experimentLabel, sep = "_"))
    }
    plotGrid <- do.call(updateOutputsFn, output_arguments)
    logInfoFn("Log2 peptide matrix calculated")
    removeNotificationFn("log2_matrix_working")
    showNotificationFn("Log2 peptide matrix calculated", type = "message")
    list(status = "success", log2Result = log2Result, plotGrid = plotGrid)
  }, error = function(e) {
    errorMessage <- paste("Error calculating the log2 peptide matrix:", e$message)
    logErrorFn(errorMessage)
    showNotificationFn(errorMessage, type = "error", duration = 15)
    removeNotificationFn("log2_matrix_working")
    list(status = "error", errorMessage = errorMessage)
  })
  })
}

runPeptideLog2RevertObserver <- function(workflowData, output,
                                         renderTextFn = shiny::renderText,
                                         showNotificationFn = shiny::showNotification,
                                         logErrorFn = logger::log_error) {
  withDIAStateScope(getDIAActiveManager(), function() {
  tryCatch({
    history <- diaStateHistory(workflowData$state_manager)
    if (length(history) <= 1) stop("No previous state to revert to.")
    previousState <- history[length(history) - 1]
    selectDIAStateReference(workflowData$state_manager, previousState)
    output$log2_results <- renderPeptideQcResultText(paste("Reverted to previous state:", previousState), renderTextFn)
    showNotificationFn("Reverted successfully", type = "message")
    list(status = "success", previousState = previousState)
  }, error = function(e) {
    errorMessage <- paste("Error reverting:", e$message)
    logErrorFn(errorMessage)
    showNotificationFn(errorMessage, type = "error")
    list(status = "error", errorMessage = errorMessage)
  })
  })
}

#' @rdname mod_prot_qc_peptide_log2
#' @export
mod_prot_qc_peptide_log2_server <- function(id, workflow_data, omic_type, experiment_label, experiment_paths = NULL) {
  shiny::moduleServer(id, function(input, output, session) {
    log2_plot <- shiny::reactiveVal(NULL)
    log2_plot_identity <- shiny::reactiveVal("Log2 peptide matrix summary")

    managed_step <- installWorkflowManagedStep(session, input, workflow_data, "peptide_log2_matrix",
      publication = list(experiment_label = experiment_label, experiment_paths = experiment_paths, omic_type = omic_type),
      ui = list(text = "log2_results", plot = "log2_plot", ack = "workflow_ack", notification = "log2_matrix_working",
        label = "log2 peptide matrix", step_label = "4_log2_matrix", plot_identity = log2_plot_identity),
      deliver = function(result) {
        arguments <- list(output = output, log2Plot = log2_plot, log2Result = result,
          omicType = omic_type, experimentLabel = experiment_label)
        if (!is.null(experiment_paths)) arguments$projectDirs <- stats::setNames(list(experiment_paths),
          paste(omic_type, experiment_label, sep = "_"))
        do.call(updatePeptideLog2Outputs, arguments)
      })

    shiny::observeEvent(input$apply_log2_matrix, {
      if (!is.null(workflowSessionApplication(session))) {
        managed_step$apply(list())
        return(invisible(NULL))
      }
      apply_arguments <- list(workflowData = workflow_data, output = output, log2Plot = log2_plot,
        omicType = omic_type, experimentLabel = experiment_label)
      if (!is.null(experiment_paths)) apply_arguments$experimentPaths <- experiment_paths
      do.call(runPeptideLog2ApplyObserver, apply_arguments)
    })

    shiny::observeEvent(input$revert_log2_matrix, {
      managed_step$cancel("reverted")
      runPeptideLog2RevertObserver(workflowData = workflow_data, output = output)
    })

    output$log2_plot <- shiny::renderPlot({
      shiny::req(log2_plot())
      grid::grid.draw(log2_plot())
    }, alt = function() log2_plot_identity())
  })
}
