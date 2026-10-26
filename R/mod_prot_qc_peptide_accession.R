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

#' @title Peptide Accession Cleanup Module (DIA with limpa)
#'
#' @description The DIA_LIMPA route's peptide-level accession cleanup (KIT-047, the limpa script's chunk 16):
#' the peptide method of `chooseBestProteinAccession` resolves one accession per protein group from the FASTA
#' metadata and aggregates duplicate peptides. On an on-disk owner the step runs the registered
#' `peptide_accession_cleanup` declaration (managed worker or in process).
#'
#' @name mod_prot_qc_peptide_accession
#' @param id,workflow_data,omic_type,experiment_label,experiment_paths Runtime inputs used by this function; see the usage section for accepted values.
NULL

#' @rdname mod_prot_qc_peptide_accession
#' @export
#' @importFrom shiny NS tabPanel br fluidRow column wellPanel h4 p hr div actionButton verbatimTextOutput plotOutput textInput helpText selectInput
mod_prot_qc_peptide_accession_ui <- function(id) {
  ns <- shiny::NS(id)

  shiny::tabPanel(
    "Peptide Accession Cleanup",
    shiny::br(),
    shiny::fluidRow(
      shiny::column(4,
        shiny::wellPanel(
          shiny::h4("Best Protein Accession per Peptide Group"),
          shiny::p("Choose the best protein accession for each protein group using the FASTA metadata, then aggregate duplicate peptides within the resolved group."),
          shiny::hr(),
          shiny::textInput(ns("delimiter"), "Protein ID Delimiter", value = ";", width = "100%"),
          shiny::helpText("Character separating multiple protein IDs (default: ';')"),
          shiny::selectInput(ns("aggregation_method"), "Aggregation Method",
            choices = c("mean", "sum", "median"), selected = "mean", width = "100%"),
          shiny::helpText("Method for combining duplicate peptide quantities within a resolved group (the script's default: mean)"),
          shiny::hr(),
          shiny::div(
            shiny::actionButton(ns("apply_peptide_accession"), "Apply Cleanup", class = "btn-primary", width = "48%"),
            shiny::actionButton(ns("revert_peptide_accession"), "Revert", class = "btn-warning", width = "48%", style = "margin-left: 4%")
          )
        )
      ),
      shiny::column(8,
        shiny::verbatimTextOutput(ns("peptide_accession_results")),
        shiny::br(),
        shinyjqui::jqui_resizable(
          shiny::plotOutput(ns("peptide_accession_plot"), height = "800px", width = "100%")
        )
      )
    )
  )
}

runPeptideAccessionApplyStep <- function(workflowData, delimiter, aggregationMethod,
                                         chooseBestProteinAccessionFn = chooseBestProteinAccession,
                                         logInfoFn = logger::log_info,
                                         nowFn = Sys.time) {
  parameters <- list(delimiter = delimiter, aggregation_method = aggregationMethod)
  runDIAReferencePeptideAccessionCleanup(workflowData, parameters)
}

updatePeptideAccessionOutputs <- function(output, accessionPlot, accessionResult, omicType, experimentLabel,
                                          renderTextFn = shiny::renderText,
                                          updateProteinFilteringFn = updateProteinFiltering, projectDirs = NULL) {
  output$peptide_accession_results <- renderPeptideQcResultText(accessionResult$resultText, renderTextFn)
  plot_arguments <- list(
    data = accessionResult$filteringSummary %||% accessionResult$cleanedS4@peptide_data,
    step_name = "7_annotation_cleanup", omic_type = omicType, experiment_label = experimentLabel,
    return_grid = TRUE, overwrite = TRUE)
  if (!is.null(projectDirs)) plot_arguments$project_dirs <- projectDirs
  plotGrid <- if (!is.null(accessionResult$stateReference)) {
    publishDIAQvalueOutput(getDIAActiveManager(), accessionResult, plot_arguments, updateProteinFilteringFn)
  } else do.call(updateProteinFilteringFn, plot_arguments)
  accessionPlot(plotGrid)
  invisible(plotGrid)
}

runPeptideAccessionApplyObserver <- function(workflowData, delimiter, aggregationMethod, output, accessionPlot,
                                             omicType, experimentLabel,
                                             runApplyStepFn = runPeptideAccessionApplyStep,
                                             updateOutputsFn = updatePeptideAccessionOutputs,
                                             showNotificationFn = shiny::showNotification,
                                             removeNotificationFn = shiny::removeNotification,
                                             logInfoFn = logger::log_info,
                                             logErrorFn = logger::log_error, experimentPaths = NULL) {
  withDIAStateScope(getDIAActiveManager(), function() {
  showNotificationFn("Applying peptide accession cleanup...", id = "peptide_accession_working", duration = NULL)
  tryCatch({
    accessionResult <- runApplyStepFn(workflowData = workflowData, delimiter = delimiter, aggregationMethod = aggregationMethod)
    output_arguments <- list(output = output, accessionPlot = accessionPlot, accessionResult = accessionResult,
      omicType = omicType, experimentLabel = experimentLabel)
    if (!is.null(experimentPaths)) {
      output_arguments$projectDirs <- stats::setNames(list(experimentPaths), paste(omicType, experimentLabel, sep = "_"))
    }
    plotGrid <- do.call(updateOutputsFn, output_arguments)
    logInfoFn("Peptide accession cleanup completed")
    removeNotificationFn("peptide_accession_working")
    showNotificationFn("Peptide accession cleanup completed", type = "message")
    list(status = "success", accessionResult = accessionResult, plotGrid = plotGrid)
  }, error = function(e) {
    errorMessage <- paste("Error in peptide accession cleanup:", e$message)
    logErrorFn(errorMessage)
    showNotificationFn(errorMessage, type = "error", duration = 15)
    removeNotificationFn("peptide_accession_working")
    list(status = "error", errorMessage = errorMessage)
  })
  })
}

runPeptideAccessionRevertObserver <- function(workflowData, output,
                                              renderTextFn = shiny::renderText,
                                              showNotificationFn = shiny::showNotification,
                                              logErrorFn = logger::log_error) {
  withDIAStateScope(getDIAActiveManager(), function() {
  tryCatch({
    history <- diaStateHistory(workflowData$state_manager)
    if (length(history) <= 1) stop("No previous state to revert to.")
    previousState <- history[length(history) - 1]
    selectDIAStateReference(workflowData$state_manager, previousState)
    output$peptide_accession_results <- renderPeptideQcResultText(paste("Reverted to previous state:", previousState), renderTextFn)
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

#' @rdname mod_prot_qc_peptide_accession
#' @export
mod_prot_qc_peptide_accession_server <- function(id, workflow_data, omic_type, experiment_label, experiment_paths = NULL) {
  shiny::moduleServer(id, function(input, output, session) {
    accession_plot <- shiny::reactiveVal(NULL)
    accession_plot_identity <- shiny::reactiveVal("Peptide accession cleanup summary")

    managed_step <- installWorkflowManagedStep(session, input, workflow_data, "peptide_accession_cleanup",
      publication = list(experiment_label = experiment_label, experiment_paths = experiment_paths, omic_type = omic_type),
      ui = list(text = "peptide_accession_results", plot = "peptide_accession_plot", ack = "workflow_ack",
        notification = "peptide_accession_working", label = "peptide accession cleanup",
        step_label = "7_annotation_cleanup", plot_identity = accession_plot_identity),
      deliver = function(result) {
        arguments <- list(output = output, accessionPlot = accession_plot, accessionResult = result,
          omicType = omic_type, experimentLabel = experiment_label)
        if (!is.null(experiment_paths)) arguments$projectDirs <- stats::setNames(list(experiment_paths),
          paste(omic_type, experiment_label, sep = "_"))
        do.call(updatePeptideAccessionOutputs, arguments)
      })

    shiny::observeEvent(input$apply_peptide_accession, {
      if (!is.null(workflowSessionApplication(session))) {
        managed_step$apply(list(delimiter = input$delimiter, aggregation_method = input$aggregation_method))
        return(invisible(NULL))
      }
      apply_arguments <- list(workflowData = workflow_data, delimiter = input$delimiter,
        aggregationMethod = input$aggregation_method, output = output, accessionPlot = accession_plot,
        omicType = omic_type, experimentLabel = experiment_label)
      if (!is.null(experiment_paths)) apply_arguments$experimentPaths <- experiment_paths
      do.call(runPeptideAccessionApplyObserver, apply_arguments)
    })

    shiny::observeEvent(input$revert_peptide_accession, {
      managed_step$cancel("reverted")
      runPeptideAccessionRevertObserver(workflowData = workflow_data, output = output)
    })

    output$peptide_accession_plot <- shiny::renderPlot({
      shiny::req(accession_plot())
      grid::grid.draw(accession_plot())
    }, alt = function() accession_plot_identity())
  })
}
