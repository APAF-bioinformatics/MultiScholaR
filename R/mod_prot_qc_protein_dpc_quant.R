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

#' @title DPC Quantification Module (DIA with limpa)
#'
#' @description The DIA_LIMPA route's DPC imputation and protein quantification (KIT-047, the limpa script's
#' chunk 17): the peptide method of `proteinMissingValueImputationLimpa` on the declared engine (`limpa` by
#' default, `fastlimpadpc` selectable) produces the protein object; its Limpa results travel with the
#' state as an owned value. On an on-disk owner the step runs the registered `protein_dpc_quant` declaration
#' (managed worker or in process).
#'
#' @name mod_prot_qc_protein_dpc_quant
#' @param id,workflow_data,experiment_paths,omic_type,experiment_label Runtime inputs used by this function; see the usage section for accepted values.
NULL

#' The engines the DPC tab offers (KIT-047): the installed ones, the shipped default (limpa) first.
#' @noRd
limpaDpcEngineChoices <- function(available = function(package) requireNamespace(package, quietly = TRUE)) {
  choices <- c("limpa (reference, default)" = "limpa", "fastlimpadpc (Rust)" = "fastlimpadpc")
  choices[vapply(choices, function(engine) isTRUE(available(engine)), logical(1))]
}

#' @rdname mod_prot_qc_protein_dpc_quant
#' @export
#' @importFrom shiny NS tabPanel br fluidRow column wellPanel h4 p hr div actionButton verbatimTextOutput plotOutput radioButtons numericInput helpText
mod_prot_qc_protein_dpc_quant_ui <- function(id) {
  ns <- shiny::NS(id)
  engines <- limpaDpcEngineChoices()

  shiny::tabPanel(
    "DPC Quantification",
    shiny::br(),
    shiny::fluidRow(
      shiny::column(4,
        shiny::wellPanel(
          shiny::h4("Limpa DPC Imputation and Quantification"),
          shiny::p("Fit the detection probability curve on the log2 peptide matrix and quantify proteins with DPC-Quant (Li et al. 2025). The protein object carries the quantified matrix, its standard errors and the peptide counts per protein."),
          shiny::hr(),
          if (length(engines)) shiny::radioButtons(ns("engine"), "Engine", choices = engines, selected = engines[[1L]]) else
            shiny::p("Neither fastlimpadpc nor limpa is installed; the step cannot run."),
          shiny::helpText("The engine and its version are recorded with the step."),
          shiny::numericInput(ns("dpc_slope"), "DPC slope assumed when the fit fails", value = 0.8, min = 0.01, step = 0.1, width = "100%"),
          shiny::numericInput(ns("chunk"), "Proteins per quantification chunk", value = 1000, min = 1, step = 100, width = "100%"),
          shiny::hr(),
          shiny::div(
            shiny::actionButton(ns("apply_dpc_quant"), "Apply DPC Quantification", class = "btn-primary", width = "48%"),
            shiny::actionButton(ns("revert_dpc_quant"), "Revert", class = "btn-warning", width = "48%", style = "margin-left: 4%")
          )
        )
      ),
      shiny::column(8,
        shiny::verbatimTextOutput(ns("dpc_quant_results")),
        shiny::br(),
        shinyjqui::jqui_resizable(
          shiny::plotOutput(ns("dpc_quant_plot"), height = "800px", width = "100%")
        )
      )
    )
  )
}

runProteinDpcQuantStep <- function(workflowData, engine, dpcSlope, chunk,
                                   proteinMissingValueImputationLimpaFn = proteinMissingValueImputationLimpa,
                                   logInfoFn = logger::log_info,
                                   nowFn = Sys.time) {
  parameters <- list(engine = engine, dpc_slope = as.numeric(dpcSlope), chunk = as.numeric(chunk))
  if (usesDIAReferenceStates(workflowData$state_manager)) return(runDIAReferenceProteinDpcQuant(workflowData, parameters))

  shiny::req(workflowData$state_manager)
  manager <- workflowData$state_manager
  predecessor <- getDIACurrentReference(manager)
  currentS4 <- getDIAExecutionState(manager, "proteinMissingValueImputationLimpa")
  shiny::req(currentS4)
  if (!length(currentS4@peptide_matrix)) stop("DPC quantification needs the log2 peptide matrix state.", call. = FALSE)

  logInfoFn(sprintf("Protein Processing: Limpa DPC quantification on engine %s", engine))
  quantifiedS4 <- proteinMissingValueImputationLimpaFn(currentS4, dpc_slope = parameters$dpc_slope, verbose = FALSE,
    chunk = parameters$chunk, engine = engine)
  results <- quantifiedS4@args$limpa_dpc_quant_results
  dpc <- unname(as.numeric(results$dpc_parameters_used))
  config <- list(engine = results$engine, engine_version = results$engine_version, dpc_slope = parameters$dpc_slope,
    chunk = parameters$chunk, dpc_intercept = dpc[[1L]], dpc_slope_realised = dpc[[2L]],
    slope_interpretation = results$slope_interpretation, proteins_quantified = results$total_proteins_quantified,
    peptides_used = results$total_peptides_used, missing_percentage_peptides = results$missing_percentage_peptides,
    protein_id_column = quantifiedS4@protein_id_column)
  quantifiedS4@args$limpaDpcQuant <- config[c("engine", "engine_version", "dpc_slope", "chunk")]
  if (is.null(workflowData$qc_params)) workflowData$qc_params <- list()
  if (is.null(workflowData$qc_params$protein_qc)) workflowData$qc_params$protein_qc <- list()
  workflowData$qc_params$protein_qc$dpc_quant <- c(config, list(timestamp = nowFn()))

  quantifiedS4 <- inheritDIAExecutionSource(quantifiedS4, manager, predecessor)
  manager$saveState(state_name = "protein_dpc_quant", s4_data_object = quantifiedS4, config_object = config,
    description = "Limpa DPC imputation and protein quantification completed (ProteinQuantitativeData created)")
  releaseDIAWorkingReference(manager, predecessor)
  currentS4 <- NULL
  proteinCount <- length(unique(quantifiedS4@protein_quant_table[[quantifiedS4@protein_id_column]]))
  resultText <- paste(
    "Limpa DPC Quantification Completed Successfully\n",
    "===============================================\n",
    sprintf("Proteins quantified: %d\n", proteinCount),
    sprintf("Samples: %d\n", nrow(quantifiedS4@design_matrix)),
    sprintf("Engine: %s %s\n", config$engine, config$engine_version),
    sprintf("DPC intercept / slope: %.4f / %.4f (%s)\n", config$dpc_intercept, config$dpc_slope_realised, config$slope_interpretation),
    "S4 Class: ProteinQuantitativeData\n",
    "State saved as: 'protein_dpc_quant'\n")

  finishDIAProteinAction(workflowData, list(quantifiedS4 = quantifiedS4, resultText = resultText), "quantifiedS4")
}

updateProteinDpcQuantOutputs <- function(output, dpcQuantPlot, dpcResult, omicType, experimentLabel,
                                         renderTextFn = shiny::renderText,
                                         updateProteinFilteringFn = updateProteinFiltering) {
  output$dpc_quant_results <- renderPeptideQcResultText(dpcResult$resultText, renderTextFn)
  plotGrid <- publishDIAProteinProgress(dpcResult, "9_protein_dpc_quant", omicType, experimentLabel,
    updateProteinFilteringFn, dpcResult$quantifiedS4)
  dpcQuantPlot(plotGrid)
  invisible(plotGrid)
}

runProteinDpcQuantApplyObserver <- function(workflowData, engine, dpcSlope, chunk, output, dpcQuantPlot,
                                            omicType, experimentLabel,
                                            runApplyStepFn = runProteinDpcQuantStep,
                                            updateOutputsFn = updateProteinDpcQuantOutputs,
                                            showNotificationFn = shiny::showNotification,
                                            removeNotificationFn = shiny::removeNotification,
                                            logInfoFn = logger::log_info,
                                            logErrorFn = logger::log_error) {
  withDIAStateScope(getDIAActiveManager(), function() {
  showNotificationFn("Running Limpa DPC quantification...", id = "dpc_quant_working", duration = NULL)
  tryCatch({
    dpcResult <- runApplyStepFn(workflowData = workflowData, engine = engine, dpcSlope = dpcSlope, chunk = chunk)
    plotGrid <- updateOutputsFn(output = output, dpcQuantPlot = dpcQuantPlot, dpcResult = dpcResult,
      omicType = omicType, experimentLabel = experimentLabel)
    logInfoFn("Limpa DPC quantification completed")
    removeNotificationFn("dpc_quant_working")
    showNotificationFn("Limpa DPC quantification completed", type = "message")
    list(status = "success", dpcResult = dpcResult, plotGrid = plotGrid)
  }, error = function(e) {
    errorMessage <- paste("Error in Limpa DPC quantification:", e$message)
    logErrorFn(errorMessage)
    showNotificationFn(errorMessage, type = "error", duration = 15)
    removeNotificationFn("dpc_quant_working")
    list(status = "error", errorMessage = errorMessage)
  })
  })
}

runProteinDpcQuantRevertObserver <- function(workflowData, output,
                                             renderTextFn = shiny::renderText,
                                             showNotificationFn = shiny::showNotification,
                                             logErrorFn = logger::log_error) {
  withDIAStateScope(getDIAActiveManager(), function() {
  tryCatch({
    history <- workflowData$state_manager$getHistory()
    if (length(history) <= 1) stop("No previous state to revert to.")
    previousState <- history[length(history) - 1]
    selectDIAStateReference(workflowData$state_manager, previousState)
    output$dpc_quant_results <- renderPeptideQcResultText(paste("Reverted to previous state:", previousState), renderTextFn)
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

#' @rdname mod_prot_qc_protein_dpc_quant
#' @export
mod_prot_qc_protein_dpc_quant_server <- function(id, workflow_data, experiment_paths = NULL, omic_type, experiment_label) {
  shiny::moduleServer(id, function(input, output, session) {
    dpc_quant_plot <- shiny::reactiveVal(NULL)
    dpc_quant_plot_identity <- shiny::reactiveVal("Limpa DPC quantification summary")

    managed_step <- installWorkflowManagedStep(session, input, workflow_data, "protein_dpc_quant",
      publication = list(experiment_label = experiment_label, experiment_paths = experiment_paths, omic_type = omic_type),
      ui = list(text = "dpc_quant_results", plot = "dpc_quant_plot", ack = "workflow_ack", notification = "dpc_quant_working",
        label = "Limpa DPC quantification", step_label = "9_protein_dpc_quant", plot_identity = dpc_quant_plot_identity),
      deliver = function(result) updateProteinDpcQuantOutputs(output = output, dpcQuantPlot = dpc_quant_plot,
        dpcResult = result, omicType = omic_type, experimentLabel = experiment_label),
      read_result = readManagedProteinResult)

    shiny::observeEvent(input$apply_dpc_quant, {
      parameters <- list(engine = input$engine, dpc_slope = as.numeric(input$dpc_slope), chunk = as.numeric(input$chunk))
      if (!is.null(workflowSessionApplication(session)) && usesDIAReferenceStates(workflow_data$state_manager)) {
        managed_step$apply(parameters)
        return(invisible(NULL))
      }
      runProteinDpcQuantApplyObserver(workflowData = workflow_data, engine = parameters$engine,
        dpcSlope = parameters$dpc_slope, chunk = parameters$chunk, output = output, dpcQuantPlot = dpc_quant_plot,
        omicType = omic_type, experimentLabel = experiment_label)
    })

    shiny::observeEvent(input$revert_dpc_quant, {
      managed_step$cancel("reverted")
      runProteinDpcQuantRevertObserver(workflowData = workflow_data, output = output)
    })

    output$dpc_quant_plot <- shiny::renderPlot({
      shiny::req(dpc_quant_plot())
      grid::grid.draw(dpc_quant_plot())
    }, alt = function() dpc_quant_plot_identity())
  })
}
