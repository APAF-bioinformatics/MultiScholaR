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

#' @title Peptide Count Filter Module (DIA with limpa)
#'
#' @description The DIA_LIMPA route's minimum-peptide-evidence filter on the DPC-quantified protein object
#' (KIT-041, the limpa script's chunk 18): proteins with fewer unique peptides or peptidoforms than the thresholds
#' are removed from the matrix and the quantified EList. On an on-disk owner the step runs the registered
#' `protein_peptide_count_filter` declaration (managed worker or in process).
#'
#' @name mod_prot_qc_protein_peptide_count
#' @param id,workflow_data,experiment_paths,omic_type,experiment_label Runtime inputs used by this function; see the usage section for accepted values.
NULL

#' @rdname mod_prot_qc_protein_peptide_count
#' @export
#' @importFrom shiny NS tabPanel br fluidRow column wellPanel h4 p hr div actionButton verbatimTextOutput plotOutput numericInput helpText
mod_prot_qc_protein_peptide_count_ui <- function(id) {
  ns <- shiny::NS(id)

  shiny::tabPanel(
    "Peptide Count Filter",
    shiny::br(),
    shiny::fluidRow(
      shiny::column(4,
        shiny::wellPanel(
          shiny::h4("Minimum Peptide Evidence per Protein"),
          shiny::p("Keep proteins quantified from at least the given numbers of unique peptides and peptidoforms (counts from the DPC quantification step)."),
          shiny::hr(),
          shiny::numericInput(ns("num_peptides_per_protein_thresh"), "Minimum unique peptides per protein", value = 1, min = 0, step = 1, width = "100%"),
          shiny::numericInput(ns("num_peptidoforms_per_protein_thresh"), "Minimum peptidoforms per protein", value = 2, min = 0, step = 1, width = "100%"),
          shiny::helpText("The limpa workflow's defaults: 1 peptide and 2 peptidoforms."),
          shiny::hr(),
          shiny::div(
            shiny::actionButton(ns("apply_peptide_count"), "Apply Filter", class = "btn-primary", width = "48%"),
            shiny::actionButton(ns("revert_peptide_count"), "Revert", class = "btn-warning", width = "48%", style = "margin-left: 4%")
          )
        )
      ),
      shiny::column(8,
        shiny::verbatimTextOutput(ns("peptide_count_results")),
        shiny::br(),
        shinyjqui::jqui_resizable(
          shiny::plotOutput(ns("peptide_count_plot"), height = "800px", width = "100%")
        )
      )
    )
  )
}

runProteinPeptideCountStep <- function(workflowData, minPeptides, minPeptidoforms,
                                       filterMinNumPeptidesPerProteinFn = filterMinNumPeptidesPerProtein,
                                       logInfoFn = logger::log_info,
                                       nowFn = Sys.time) {
  parameters <- list(num_peptides_per_protein_thresh = as.numeric(minPeptides), num_peptidoforms_per_protein_thresh = as.numeric(minPeptidoforms))
  if (usesDIAReferenceStates(workflowData$state_manager)) return(runDIAReferenceProteinPeptideCountFilter(workflowData, parameters))

  shiny::req(workflowData$state_manager)
  manager <- workflowData$state_manager
  predecessor <- getDIACurrentReference(manager)
  currentS4 <- getDIAExecutionState(manager, "filterMinNumPeptidesPerProtein")
  shiny::req(currentS4)
  logInfoFn(sprintf("Protein Processing: peptide count filter (%s peptides, %s peptidoforms)", minPeptides, minPeptidoforms))
  proteinsBefore <- nrow(currentS4@protein_quant_table)
  filteredS4 <- suppressMessages(filterMinNumPeptidesPerProteinFn(currentS4,
    num_peptides_per_protein_thresh = parameters$num_peptides_per_protein_thresh,
    num_peptidoforms_per_protein_thresh = parameters$num_peptidoforms_per_protein_thresh, verbose = FALSE))
  proteinsAfter <- nrow(filteredS4@protein_quant_table)
  config <- list(num_peptides_per_protein_thresh = parameters$num_peptides_per_protein_thresh,
    num_peptidoforms_per_protein_thresh = parameters$num_peptidoforms_per_protein_thresh,
    proteins_before = proteinsBefore, proteins_after = proteinsAfter)
  if (is.null(workflowData$qc_params)) workflowData$qc_params <- list()
  if (is.null(workflowData$qc_params$protein_qc)) workflowData$qc_params$protein_qc <- list()
  workflowData$qc_params$protein_qc$peptide_count_filter <- c(config, list(timestamp = nowFn()))
  filteredS4 <- inheritDIAExecutionSource(filteredS4, manager, predecessor)
  manager$saveState(state_name = "peptide_count_filtered", s4_data_object = filteredS4, config_object = config,
    description = "Applied minimum peptide evidence filter on the quantified protein object")
  releaseDIAWorkingReference(manager, predecessor)
  currentS4 <- NULL
  resultText <- paste(
    "Peptide Count Filter Applied Successfully\n",
    "=========================================\n",
    sprintf("Proteins remaining: %d (removed %d)\n", proteinsAfter, proteinsBefore - proteinsAfter),
    sprintf("Min unique peptides per protein: %d\n", as.integer(minPeptides)),
    sprintf("Min peptidoforms per protein: %d\n", as.integer(minPeptidoforms)),
    "State saved as: 'peptide_count_filtered'\n")
  finishDIAProteinAction(workflowData, list(filteredS4 = filteredS4, resultText = resultText), "filteredS4")
}

updateProteinPeptideCountOutputs <- function(output, peptideCountPlot, filterResult, omicType, experimentLabel,
                                             renderTextFn = shiny::renderText,
                                             updateProteinFilteringFn = updateProteinFiltering) {
  output$peptide_count_results <- renderPeptideQcResultText(filterResult$resultText, renderTextFn)
  plotGrid <- publishDIAProteinProgress(filterResult, "11_peptide_count_filtered", omicType, experimentLabel,
    updateProteinFilteringFn, filterResult$filteredS4)
  peptideCountPlot(plotGrid)
  invisible(plotGrid)
}

runProteinPeptideCountApplyObserver <- function(workflowData, minPeptides, minPeptidoforms, output, peptideCountPlot, omicType, experimentLabel,
                                                runApplyStepFn = runProteinPeptideCountStep,
                                                updateOutputsFn = updateProteinPeptideCountOutputs,
                                                showNotificationFn = shiny::showNotification,
                                                removeNotificationFn = shiny::removeNotification,
                                                logInfoFn = logger::log_info,
                                                logErrorFn = logger::log_error) {
  withDIAStateScope(getDIAActiveManager(), function() {
  showNotificationFn("Applying peptide count filter...", id = "peptide_count_working", duration = NULL)
  tryCatch({
    filterResult <- runApplyStepFn(workflowData = workflowData, minPeptides = minPeptides, minPeptidoforms = minPeptidoforms)
    plotGrid <- updateOutputsFn(output = output, peptideCountPlot = peptideCountPlot, filterResult = filterResult,
      omicType = omicType, experimentLabel = experimentLabel)
    logInfoFn("Peptide count filter completed")
    removeNotificationFn("peptide_count_working")
    showNotificationFn("Peptide count filter completed", type = "message")
    list(status = "success", filterResult = filterResult, plotGrid = plotGrid)
  }, error = function(e) {
    errorMessage <- paste("Error in peptide count filter:", e$message)
    logErrorFn(errorMessage)
    showNotificationFn(errorMessage, type = "error", duration = 15)
    removeNotificationFn("peptide_count_working")
    list(status = "error", errorMessage = errorMessage)
  })
  })
}

runProteinPeptideCountRevertObserver <- function(workflowData, output,
                                                 renderTextFn = shiny::renderText,
                                                 showNotificationFn = shiny::showNotification,
                                                 logErrorFn = logger::log_error) {
  withDIAStateScope(getDIAActiveManager(), function() {
  tryCatch({
    history <- workflowData$state_manager$getHistory()
    if (length(history) <= 1) stop("No previous state to revert to.")
    previousState <- history[length(history) - 1]
    selectDIAStateReference(workflowData$state_manager, previousState)
    output$peptide_count_results <- renderPeptideQcResultText(paste("Reverted to previous state:", previousState), renderTextFn)
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

#' @rdname mod_prot_qc_protein_peptide_count
#' @export
mod_prot_qc_protein_peptide_count_server <- function(id, workflow_data, experiment_paths = NULL, omic_type, experiment_label) {
  shiny::moduleServer(id, function(input, output, session) {
    peptide_count_plot <- shiny::reactiveVal(NULL)
    peptide_count_plot_identity <- shiny::reactiveVal("Peptide count filter summary")

    managed_step <- installWorkflowManagedStep(session, input, workflow_data, "protein_peptide_count_filter",
      publication = list(experiment_label = experiment_label, experiment_paths = experiment_paths, omic_type = omic_type),
      ui = list(text = "peptide_count_results", plot = "peptide_count_plot", ack = "workflow_ack",
        notification = "peptide_count_working", label = "peptide count filter",
        step_label = "11_peptide_count_filtered", plot_identity = peptide_count_plot_identity),
      deliver = function(result) updateProteinPeptideCountOutputs(output = output, peptideCountPlot = peptide_count_plot,
        filterResult = result, omicType = omic_type, experimentLabel = experiment_label),
      read_result = readManagedProteinResult)

    shiny::observeEvent(input$apply_peptide_count, {
      parameters <- list(num_peptides_per_protein_thresh = as.numeric(input$num_peptides_per_protein_thresh),
        num_peptidoforms_per_protein_thresh = as.numeric(input$num_peptidoforms_per_protein_thresh))
      if (!is.null(workflowSessionApplication(session)) && usesDIAReferenceStates(workflow_data$state_manager)) {
        managed_step$apply(parameters)
        return(invisible(NULL))
      }
      runProteinPeptideCountApplyObserver(workflowData = workflow_data, minPeptides = parameters$num_peptides_per_protein_thresh,
        minPeptidoforms = parameters$num_peptidoforms_per_protein_thresh, output = output, peptideCountPlot = peptide_count_plot,
        omicType = omic_type, experimentLabel = experiment_label)
    })

    shiny::observeEvent(input$revert_peptide_count, {
      managed_step$cancel("reverted")
      runProteinPeptideCountRevertObserver(workflowData = workflow_data, output = output)
    })

    output$peptide_count_plot <- shiny::renderPlot({
      shiny::req(peptide_count_plot())
      grid::grid.draw(peptide_count_plot())
    }, alt = function() peptide_count_plot_identity())
  })
}
