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

#' @title IQ Protein Rollup Module
#'
#' @description A Shiny module for performing peptide-to-protein rollup using the IQ tool.
#'
#' @name mod_prot_qc_protein_rollup
#' @param id,workflowData,experimentPaths,createProteinDataFn,writeTsvFn,processLongFormatFn,readTsvFn,captureCheckpointFn,showNotificationFn,logInfoFn,logWarnFn,sleepFn,maxWait Runtime inputs used by this function; see the usage section for accepted values.
NULL

#' @rdname mod_prot_qc_protein_rollup
#' @export
#' @importFrom shiny NS tagList tabPanel br fluidRow column wellPanel h4 p hr div actionButton verbatimTextOutput plotOutput
mod_prot_qc_protein_rollup_ui <- function(id) {
  ns <- shiny::NS(id)
  
  shiny::tabPanel(
    "IQ Protein Rollup",
    shiny::br(),
    shiny::fluidRow(
      shiny::column(4,
        shiny::wellPanel(
          shiny::h4("Peptide-to-Protein Rollup"),
          shiny::p("Aggregate peptide-level data to protein-level quantification using the IQ algorithm (MaxLFQ)."),
          shiny::hr(),
          shiny::p("This step uses the IQ tool to implement the MaxLFQ algorithm for protein quantification, then automatically creates a ProteinQuantitativeData S4 object."),
          shiny::hr(),
          shiny::div(
            shiny::actionButton(ns("apply_iq_rollup"), "Run IQ Rollup & Create S4 Object", 
              class = "btn-primary", width = "48%"),
            shiny::actionButton(ns("revert_iq_rollup"), "Revert", 
              class = "btn-warning", width = "48%", style = "margin-left: 4%")
          )
        )
      ),
      shiny::column(8,
        shiny::verbatimTextOutput(ns("iq_rollup_results")),
        shiny::br(),
        shinyjqui::jqui_resizable(
          shiny::plotOutput(ns("iq_rollup_plot"), height = "800px", width = "100%")
        )
      )
    )
  )
}

# KIT-046: the in-memory DPC rollup (test-mode object and apply step) is gone; that route is the DIA_LIMPA workflow type.

runProteinIqRollupApplyStep <- function(workflowData,
                                        experimentPaths,
                                        createProteinDataFn = ProteinQuantitativeData,
                                        writeTsvFn = vroom::vroom_write,
                                        processLongFormatFn = iq::process_long_format,
                                        readTsvFn = vroom::vroom,
                                        captureCheckpointFn = .capture_checkpoint,
                                        showNotificationFn = shiny::showNotification,
                                        logInfoFn = logger::log_info,
                                        logWarnFn = logger::log_warn,
                                        sleepFn = Sys.sleep,
                                        maxWait = 30) {
  shiny::req(workflowData$state_manager)
  pending <- resumeDIAProteinAction(workflowData, "protein_s4_created")
  if (!is.null(pending)) return(pending)
  assertDIAExportReady(workflowData$state_manager)
  science_request <- workflowExportPaths(experimentPaths$source_dir, basename(tempfile("iq-science-")))
  experimentPaths$peptide_qc_dir <- file.path(science_request$request_dir, "peptide")
  experimentPaths$protein_qc_dir <- file.path(science_request$request_dir, "protein")
  dir.create(experimentPaths$peptide_qc_dir)
  dir.create(experimentPaths$protein_qc_dir)
  # KIT-015: an unmanaged launch runs the registered declaration in process (the IQ files under the science
  # request directories, published by the finishing step as before).
  entry <- workflowActionSpec("protein_iq_rollup")
  value <- computeDIAReferenceIQRollup(workflowData$state_manager, getDIAStateReference(workflowData$state_manager),
    list(peptide_qc_dir = experimentPaths$peptide_qc_dir, protein_qc_dir = experimentPaths$protein_qc_dir),
    variant = entry$default_variant, process_long_format = processLongFormatFn, write_fn = writeTsvFn, read_fn = readTsvFn)
  if (length(value$dropped_samples)) {
    droppedList <- paste(value$dropped_samples, collapse = ", ")
    logWarnFn(paste0("Protein Processing: ", length(value$dropped_samples), " samples dropped during IQ rollup: ", droppedList))
    showNotificationFn(paste("Warning: The following samples were dropped during protein rollup due to low peptide quality:",
      droppedList), type = "warning", duration = 10)
  }
  return(finishDIAReferenceProteinStep(workflowData, value))
}

updateProteinIqRollupOutputs <- function(output,
                                         iqRollupPlot,
                                         iqRollupResult,
                                         omicType,
                                         experimentLabel,
                                         renderTextFn = shiny::renderText,
                                         updateProteinFilteringFn = updateProteinFiltering) {
  output$iq_rollup_results <- renderPeptideQcResultText(iqRollupResult$resultText, renderTextFn)

  plotGrid <- publishDIAProteinProgress(iqRollupResult, "9_protein_s4_created", omicType,
    experimentLabel, updateProteinFilteringFn, iqRollupResult$proteinObj)
  iqRollupPlot(plotGrid)

  invisible(plotGrid)
}

runProteinIqRollupApplyObserver <- function(workflowData,
                                            experimentPaths,
                                            output,
                                            iqRollupPlot,
                                            omicType,
                                            experimentLabel,
                                            runApplyStepFn = runProteinIqRollupApplyStep,
                                            updateOutputsFn = updateProteinIqRollupOutputs,
                                            showNotificationFn = shiny::showNotification,
                                            removeNotificationFn = shiny::removeNotification,
                                            logInfoFn = logger::log_info,
                                            logErrorFn = logger::log_error) {
  withDIAStateScope(getDIAActiveManager(), function() {
  showNotificationFn(
    "Running IQ protein rollup & creating S4 object...",
    id = "iq_rollup_working",
    duration = NULL
  )

  tryCatch({
    iqRollupResult <- runApplyStepFn(
      workflowData = workflowData,
      experimentPaths = experimentPaths
    )

    plotGrid <- updateOutputsFn(
      output = output,
      iqRollupPlot = iqRollupPlot,
      iqRollupResult = iqRollupResult,
      omicType = omicType,
      experimentLabel = experimentLabel
    )

    logInfoFn("IQ protein rollup and S4 object creation completed successfully")
    removeNotificationFn("iq_rollup_working")
    showNotificationFn(
      "IQ protein rollup & S4 object creation completed successfully",
      type = "message"
    )

    list(
      status = "success",
      iqRollupResult = iqRollupResult,
      plotGrid = plotGrid
    )
  }, error = function(e) {
    errorMessage <- paste("Error in IQ protein rollup & S4 creation:", e$message)
    logErrorFn(errorMessage)
    showNotificationFn(errorMessage, type = "error", duration = 15)
    removeNotificationFn("iq_rollup_working")

    list(
      status = "error",
      errorMessage = errorMessage
    )
  })
  })
}

runProteinIqRollupRevertStep <- function(workflowData) {
  shiny::req(workflowData$state_manager)
  history <- diaStateHistory(workflowData$state_manager)
  peptideStates <- c(
    "imputed",
    "replicate_filtered",
    "sample_filtered",
    "protein_peptide_filtered",
    "intensity_filtered",
    "precursor_rollup",
    "qvalue_filtered",
    "raw_data_s4"
  )
  previousState <- intersect(rev(history), peptideStates)[1]

  if (length(previousState) == 0 || is.na(previousState)) {
    stop("No previous peptide state to revert to.")
  }

  revertedS4 <- selectDIAStateReference(workflowData$state_manager, previousState)

  list(
    previousState = previousState,
    revertedS4 = revertedS4,
    resultText = paste("Reverted to", previousState, "state")
  )
}

runProteinIqRollupRevertObserver <- function(workflowData,
                                             output,
                                             runRevertStepFn = runProteinIqRollupRevertStep,
                                             renderTextFn = shiny::renderText,
                                             showNotificationFn = shiny::showNotification,
                                             logInfoFn = logger::log_info,
                                             logErrorFn = logger::log_error) {
  withDIAStateScope(getDIAActiveManager(), function() {
  tryCatch({
    revertResult <- runRevertStepFn(workflowData = workflowData)
    output$iq_rollup_results <- renderPeptideQcResultText(revertResult$resultText, renderTextFn)
    logInfoFn(paste("Reverted IQ rollup to", revertResult$previousState))
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
  })
}

bindProteinIqRollupPlot <- function(output, iqRollupPlot, plotIdentity = NULL) {
  output$iq_rollup_plot <- if (is.null(plotIdentity)) renderPlot({
    req(iqRollupPlot())
    grid.draw(iqRollupPlot())
  }) else renderPlot({
    req(iqRollupPlot())
    grid.draw(iqRollupPlot())
  }, alt = function() plotIdentity())
}

mod_prot_qc_protein_rollup_server <- function(id, workflow_data, experiment_paths, omic_type, experiment_label) {
  shiny::moduleServer(id, function(input, output, session) {
    
    iq_rollup_plot <- shiny::reactiveVal(NULL)
    iq_rollup_plot_identity <- shiny::reactiveVal("IQ protein rollup summary")
    iq_rollup_plot_value <- iq_rollup_plot

    # KIT-015: the managed worker when the session carries a managed application (KIT-008 pattern), else the
    # same registered declaration in process.
    managed_step <- installWorkflowManagedStep(session, input, workflow_data, "protein_iq_rollup",
      publication = list(experiment_label = experiment_label, experiment_paths = experiment_paths, omic_type = omic_type),
      ui = list(text = "iq_rollup_results", plot = "iq_rollup_plot", ack = "workflow_ack", notification = "iq_rollup_working",
        label = "IQ protein rollup", step_label = "9_protein_s4_created", plot_identity = iq_rollup_plot_identity),
      deliver = function(result) updateProteinIqRollupOutputs(output = output, iqRollupPlot = iq_rollup_plot_value, iqRollupResult = result,
        omicType = omic_type, experimentLabel = experiment_label),
      read_result = readManagedProteinResult)
    # Step 1: Protein Rollup (chunk 17). KIT-046: the DIA route rolls up with IQ only; limpa is its own workflow type.
    shiny::observeEvent(input$apply_iq_rollup, {
      applyStepFn <- runProteinIqRollupApplyStep
      if (!is.null(workflowSessionApplication(session))) {
        managed_step$apply(list())
        return(invisible(NULL))
      }

      runProteinIqRollupApplyObserver(
        workflowData = workflow_data,
        experimentPaths = experiment_paths,
        output = output,
        iqRollupPlot = iq_rollup_plot,
        omicType = omic_type,
        experimentLabel = experiment_label,
        runApplyStepFn = applyStepFn
      )
    })
    
    # Revert IQ Rollup
    shiny::observeEvent(input$revert_iq_rollup, {
      managed_step$cancel("reverted")
      runProteinIqRollupRevertObserver(
        workflowData = workflow_data,
        output = output
      )
    })

    bindProteinIqRollupPlot(
      output = output,
      iqRollupPlot = iq_rollup_plot,
      plotIdentity = iq_rollup_plot_identity
    )
  })
}
