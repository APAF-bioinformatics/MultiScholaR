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
#' @param id,object Runtime inputs used by this function; see the usage section for accepted values.
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
          shiny::p("Sum linear precursor intensities across charge states and modified/unmodified forms to stripped-peptide level."),
          shiny::hr(),
          shiny::p("This step has no configurable parameters. Duplicate precursor identities and logged or invalid quantities are rejected before aggregation; all-missing outputs remain missing."),
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
.peptideRollupSummaryFromObject <- function(object) {
  if (!base::isS4(object) || !"args" %in% methods::slotNames(object)) {
    return(NULL)
  }
  object@args$rollUpPrecursorToPeptide$rollup_summary
}

.formatPeptideRollupSummary <- function(summary) {
  if (is.null(summary)) {
    return("Rollup audit: unavailable for this legacy result")
  }

  paste(
    sprintf("Aggregation: %s", summary$aggregation_rule),
    sprintf("Active protein key: %s", summary$active_protein_key),
    sprintf("Precursor identity: %s", summary$precursor_identity_mode),
    sprintf("Input precursor rows: %d", summary$input_precursor_rows),
    sprintf("Unique precursor identities: %d", summary$unique_precursor_identities),
    sprintf("Output stripped-peptide identities: %d", summary$output_stripped_peptide_identities),
    sprintf(
      "Raw partial/all-missing outputs: %d / %d",
      summary$raw_partial_missing_outputs,
      summary$raw_all_missing_outputs
    ),
    sprintf(
      "Normalised partial/all-missing outputs: %d / %d",
      summary$normalised_partial_missing_outputs,
      summary$normalised_all_missing_outputs
    ),
    sprintf("Frozen identification evidence: %s", summary$frozen_identification_evidence),
    sprintf("Protein.Ids provenance: %s", summary$protein_ids_provenance),
    sep = "\n"
  )
}

runPeptideRollupApplyStep <- function(workflowData,
                                      rollupFn = rollUpPrecursorToPeptide,
                                      logInfoFn = logger::log_info,
                                      nowFn = Sys.time) {
  runDIAReferencePrecursorRollup(workflowData)
}

updatePeptideRollupOutputs <- function(output,
                                       rollupPlot,
                                       rollupResult,
                                       omicType,
                                       experimentLabel,
                                       renderTextFn = shiny::renderText,
                                       updateProteinFilteringFn = updateProteinFiltering, projectDirs = NULL) {
  output$rollup_results <- renderPeptideQcResultText(rollupResult$resultText, renderTextFn)

  plot_arguments <- list(
    data = rollupResult$filteringSummary %||% rollupResult$rolledUpS4@peptide_data,
    step_name = "4_precursor_rollup",
    omic_type = omicType,
    experiment_label = experimentLabel,
    return_grid = TRUE,
    overwrite = TRUE
  )
  if (!is.null(projectDirs)) plot_arguments$project_dirs <- projectDirs
  plotGrid <- if (!is.null(rollupResult$stateReference)) {
    publishDIAQvalueOutput(getDIAActiveManager(), rollupResult, plot_arguments, updateProteinFilteringFn)
  } else do.call(updateProteinFilteringFn, plot_arguments)
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
                                          logErrorFn = logger::log_error, experimentPaths = NULL) {
  withDIAStateScope(getDIAActiveManager(), function() {
  showNotificationFn(
    "Applying precursor rollup...",
    id = "rollup_working",
    duration = NULL
  )

  tryCatch({
    rollupResult <- runApplyStepFn(workflowData = workflowData)

    output_arguments <- list(
      output = output,
      rollupPlot = rollupPlot,
      rollupResult = rollupResult,
      omicType = omicType,
      experimentLabel = experimentLabel
    )
    if (!is.null(experimentPaths)) {
      output_arguments$projectDirs <- stats::setNames(
        list(experimentPaths), paste(omicType, experimentLabel, sep = "_")
      )
    }
    plotGrid <- do.call(updateOutputsFn, output_arguments)

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
  })
}

runPeptideRollupRevertStep <- function(workflowData,
                                       logInfoFn = logger::log_info) {
  history <- diaStateHistory(workflowData$state_manager)
  if (length(history) <= 1) {
    stop("No previous state to revert to.")
  }

  previousState <- history[length(history) - 1]
  revertedS4 <- selectDIAStateReference(workflowData$state_manager, previousState)
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
  withDIAStateScope(getDIAActiveManager(), function() {
  tryCatch({
    revertResult <- runRevertStepFn(workflowData = workflowData)
    output$rollup_results <- renderPeptideQcResultText(revertResult$resultText, renderTextFn)
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

mod_prot_qc_peptide_rollup_server <- function(id, workflow_data, omic_type, experiment_label, experiment_paths = NULL) {
  shiny::moduleServer(id, function(input, output, session) {
    
    rollup_plot <- shiny::reactiveVal(NULL)
    rollup_plot_identity <- shiny::reactiveVal("Precursor rollup summary")

    # KIT-008: the P14 rollup runs on the managed worker (admission cap, recovery) when the session carries a
    # managed application; an unmanaged launch runs the same declaration in process.
    managed_step <- installWorkflowManagedStep(session, input, workflow_data, "precursor_rollup",
      publication = list(experiment_label = experiment_label, experiment_paths = experiment_paths, omic_type = omic_type),
      ui = list(text = "rollup_results", plot = "rollup_plot", ack = "workflow_ack", notification = "rollup_working",
        label = "precursor rollup", step_label = "4_precursor_rollup", plot_identity = rollup_plot_identity),
      deliver = function(result) {
        arguments <- list(output = output, rollupPlot = rollup_plot, rollupResult = result,
          omicType = omic_type, experimentLabel = experiment_label)
        if (!is.null(experiment_paths)) arguments$projectDirs <- stats::setNames(list(experiment_paths),
          paste(omic_type, experiment_label, sep = "_"))
        do.call(updatePeptideRollupOutputs, arguments)
      })

    # Step 2: Precursor Rollup (chunk 11)
    shiny::observeEvent(input$apply_rollup, {
      if (!is.null(workflowSessionApplication(session))) {
        managed_step$apply(list())
        return(invisible(NULL))
      }
      apply_arguments <- list(
        workflowData = workflow_data,
        output = output,
        rollupPlot = rollup_plot,
        omicType = omic_type,
        experimentLabel = experiment_label
      )
      if (!is.null(experiment_paths)) apply_arguments$experimentPaths <- experiment_paths
      do.call(runPeptideRollupApplyObserver, apply_arguments)
    })
    
    # Revert Precursor Rollup
    shiny::observeEvent(input$revert_rollup, {
      managed_step$cancel("reverted")
      runPeptideRollupRevertObserver(
        workflowData = workflow_data,
        output = output
      )
    })
    
    # Render rollup plot
    output$rollup_plot <- shiny::renderPlot({
      shiny::req(rollup_plot())
      grid::grid.draw(rollup_plot())
    }, alt = function() rollup_plot_identity())
  })
}
