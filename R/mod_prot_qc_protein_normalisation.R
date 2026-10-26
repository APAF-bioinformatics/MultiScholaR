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

#' @title Between-Sample Normalisation Module (DIA with limpa)
#'
#' @description The DIA_LIMPA route's optional between-sample normalisation of the quantified protein matrix
#' (KIT-041, the limpa script's chunk 17b; the script skips it, the shipped default is `none`). On an on-disk
#' owner the step runs the registered `protein_between_sample_normalisation` declaration (managed worker or in
#' process).
#'
#' @name mod_prot_qc_protein_normalisation
#' @param id,workflow_data,experiment_paths,omic_type,experiment_label Runtime inputs used by this function; see the usage section for accepted values.
NULL

#' @rdname mod_prot_qc_protein_normalisation
#' @export
#' @importFrom shiny NS tabPanel br fluidRow column wellPanel h4 p hr div actionButton verbatimTextOutput plotOutput selectInput helpText
mod_prot_qc_protein_normalisation_ui <- function(id) {
  ns <- shiny::NS(id)

  shiny::tabPanel(
    "Between-Sample Normalisation",
    shiny::br(),
    shiny::fluidRow(
      shiny::column(4,
        shiny::wellPanel(
          shiny::h4("Between-Sample Normalisation"),
          shiny::p("Optionally normalise the DPC-quantified log2 protein matrix between samples. The limpa workflow carries the quantified matrix unchanged (none)."),
          shiny::hr(),
          shiny::selectInput(ns("norm_method"), "Normalisation method",
            choices = c("None (limpa default)" = "none", "Cyclic loess" = "cyclicloess", "Quantile" = "quantile", "Scale (median absolute value)" = "scale"),
            selected = "none", width = "100%"),
          shiny::helpText("The method and the resulting state are recorded; none commits the quantified matrix unchanged."),
          shiny::hr(),
          shiny::div(
            shiny::actionButton(ns("apply_protein_normalisation"), "Apply", class = "btn-primary", width = "48%"),
            shiny::actionButton(ns("revert_protein_normalisation"), "Revert", class = "btn-warning", width = "48%", style = "margin-left: 4%")
          )
        )
      ),
      shiny::column(8,
        shiny::verbatimTextOutput(ns("protein_normalisation_results")),
        shiny::br(),
        shinyjqui::jqui_resizable(
          shiny::plotOutput(ns("protein_normalisation_plot"), height = "800px", width = "100%")
        )
      )
    )
  )
}

runProteinNormalisationStep <- function(workflowData, normMethod,
                                        normaliseBetweenSamplesFn = normaliseBetweenSamples,
                                        logInfoFn = logger::log_info,
                                        nowFn = Sys.time) {
  runDIAReferenceProteinNormalisation(workflowData, parameters)
}

updateProteinNormalisationOutputs <- function(output, normalisationPlot, normalisationResult, omicType, experimentLabel,
                                              renderTextFn = shiny::renderText,
                                              updateProteinFilteringFn = updateProteinFiltering) {
  output$protein_normalisation_results <- renderPeptideQcResultText(normalisationResult$resultText, renderTextFn)
  plotGrid <- publishDIAProteinProgress(normalisationResult, "10_protein_normalised", omicType, experimentLabel,
    updateProteinFilteringFn, normalisationResult$normalisedS4)
  normalisationPlot(plotGrid)
  invisible(plotGrid)
}

runProteinNormalisationApplyObserver <- function(workflowData, normMethod, output, normalisationPlot, omicType, experimentLabel,
                                                 runApplyStepFn = runProteinNormalisationStep,
                                                 updateOutputsFn = updateProteinNormalisationOutputs,
                                                 showNotificationFn = shiny::showNotification,
                                                 removeNotificationFn = shiny::removeNotification,
                                                 logInfoFn = logger::log_info,
                                                 logErrorFn = logger::log_error) {
  withDIAStateScope(getDIAActiveManager(), function() {
  showNotificationFn("Applying between-sample normalisation...", id = "protein_normalisation_working", duration = NULL)
  tryCatch({
    normalisationResult <- runApplyStepFn(workflowData = workflowData, normMethod = normMethod)
    plotGrid <- updateOutputsFn(output = output, normalisationPlot = normalisationPlot, normalisationResult = normalisationResult,
      omicType = omicType, experimentLabel = experimentLabel)
    logInfoFn("Between-sample normalisation step completed")
    removeNotificationFn("protein_normalisation_working")
    showNotificationFn("Between-sample normalisation step completed", type = "message")
    list(status = "success", normalisationResult = normalisationResult, plotGrid = plotGrid)
  }, error = function(e) {
    errorMessage <- paste("Error in between-sample normalisation:", e$message)
    logErrorFn(errorMessage)
    showNotificationFn(errorMessage, type = "error", duration = 15)
    removeNotificationFn("protein_normalisation_working")
    list(status = "error", errorMessage = errorMessage)
  })
  })
}

runProteinNormalisationRevertObserver <- function(workflowData, output,
                                                  renderTextFn = shiny::renderText,
                                                  showNotificationFn = shiny::showNotification,
                                                  logErrorFn = logger::log_error) {
  withDIAStateScope(getDIAActiveManager(), function() {
  tryCatch({
    history <- diaStateHistory(workflowData$state_manager)
    if (length(history) <= 1) stop("No previous state to revert to.")
    previousState <- history[length(history) - 1]
    selectDIAStateReference(workflowData$state_manager, previousState)
    output$protein_normalisation_results <- renderPeptideQcResultText(paste("Reverted to previous state:", previousState), renderTextFn)
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

#' @rdname mod_prot_qc_protein_normalisation
#' @export
mod_prot_qc_protein_normalisation_server <- function(id, workflow_data, experiment_paths = NULL, omic_type, experiment_label) {
  shiny::moduleServer(id, function(input, output, session) {
    normalisation_plot <- shiny::reactiveVal(NULL)
    normalisation_plot_identity <- shiny::reactiveVal("Between-sample normalisation summary")

    managed_step <- installWorkflowManagedStep(session, input, workflow_data, "protein_between_sample_normalisation",
      publication = list(experiment_label = experiment_label, experiment_paths = experiment_paths, omic_type = omic_type),
      ui = list(text = "protein_normalisation_results", plot = "protein_normalisation_plot", ack = "workflow_ack",
        notification = "protein_normalisation_working", label = "between-sample normalisation",
        step_label = "10_protein_normalised", plot_identity = normalisation_plot_identity),
      deliver = function(result) updateProteinNormalisationOutputs(output = output, normalisationPlot = normalisation_plot,
        normalisationResult = result, omicType = omic_type, experimentLabel = experiment_label),
      read_result = readManagedProteinResult)

    shiny::observeEvent(input$apply_protein_normalisation, {
      if (!is.null(workflowSessionApplication(session))) {
        managed_step$apply(list(norm_method = input$norm_method))
        return(invisible(NULL))
      }
      runProteinNormalisationApplyObserver(workflowData = workflow_data, normMethod = input$norm_method, output = output,
        normalisationPlot = normalisation_plot, omicType = omic_type, experimentLabel = experiment_label)
    })

    shiny::observeEvent(input$revert_protein_normalisation, {
      managed_step$cancel("reverted")
      runProteinNormalisationRevertObserver(workflowData = workflow_data, output = output)
    })

    output$protein_normalisation_plot <- shiny::renderPlot({
      shiny::req(normalisation_plot())
      grid::grid.draw(normalisation_plot())
    }, alt = function() normalisation_plot_identity())
  })
}
