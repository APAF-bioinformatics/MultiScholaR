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

#' @title Peptide Q-Value Filter Module
#'
#' @description A Shiny module for applying Q-value and proteotypic peptide filters.
#'
#' @name mod_prot_qc_peptide_qvalue
#' @param id,workflowData,qvalueThreshold,globalQvalueThreshold,proteotypicOnly,updateConfigParameterFn,qvalueFilterFn,logInfoFn,logWarnFn,nowFn,globalPGQvalueThreshold,confidenceProvenanceMode Runtime inputs used by this function; see the usage section for accepted values.
NULL

#' @rdname mod_prot_qc_peptide_qvalue
#' @export
#' @importFrom shiny NS tagList tabPanel br fluidRow column wellPanel h4 p hr helpText div actionButton verbatimTextOutput plotOutput
mod_prot_qc_peptide_qvalue_ui <- function(id) {
  ns <- shiny::NS(id)
  
  shiny::tabPanel(
    "Q-Value Filter",
    shiny::br(),
    shiny::fluidRow(
      shiny::column(4,
        shiny::wellPanel(
          shiny::h4("Q-Value & Proteotypic Peptide Filter"),
          shiny::p("Filter peptides based on statistical confidence (q-value) and ensure only proteotypic peptides are retained."),
          shiny::hr(),
          
          shiny::tags$dl(
            shiny::tags$dt("Run-specific precursor Q.Value"),
            shiny::tags$dd("\u2264 0.01 (fixed)"),
            shiny::tags$dt("Global precursor Global.Q.Value"),
            shiny::tags$dd("\u2264 0.01 (fixed)"),
            shiny::tags$dt("Global protein-group Global.PG.Q.Value"),
            shiny::tags$dd("\u2264 0.01 (fixed)"),
            shiny::tags$dt("Proteotypic"),
            shiny::tags$dd("= 1 (fixed)")
          ),
          shiny::helpText("These identification-confidence criteria are fixed for reproducible DIA-NN evidence filtering."),
          
          shiny::hr(),
          shiny::div(
            shiny::actionButton(ns("apply_qvalue_filter"), "Apply Filter", 
              class = "btn-primary", width = "48%"),
            shiny::actionButton(ns("revert_qvalue"), "Revert", 
              class = "btn-warning", width = "48%", style = "margin-left: 4%")
          )
        )
      ),
      shiny::column(8,
        shiny::verbatimTextOutput(ns("qvalue_results")),
        shiny::br(),
        shinyjqui::jqui_resizable(
          shiny::plotOutput(ns("qvalue_plot"), height = "800px", width = "100%")
        )
      )
    )
  )
}

#' @rdname mod_prot_qc_peptide_qvalue
#' @export
#' @importFrom shiny moduleServer reactiveVal observeEvent req showNotification removeNotification renderText renderPlot
#' @importFrom logger log_info log_error
#' @importFrom grid grid.draw
runPeptideQvalueApplyStep <- function(workflowData,
                                      qvalueThreshold,
                                      globalQvalueThreshold,
                                      proteotypicOnly,
                                      updateConfigParameterFn = updateConfigParameter,
                                      qvalueFilterFn = srlQvalueProteotypicPeptideClean,
                                      logInfoFn = logger::log_info,
                                      logWarnFn = logger::log_warn,
                                      nowFn = Sys.time,
                                      globalPGQvalueThreshold = .diaNnIdentificationQvalueCutoff,
                                      confidenceProvenanceMode = "diann_main_report") {
  shiny::req(workflowData$state_manager)
  runDIAReferenceQvalue(workflowData, qvalueThreshold, globalQvalueThreshold,
    proteotypicOnly, globalPGQvalueThreshold, confidenceProvenanceMode, nowFn)
}

updatePeptideQvalueOutputs <- function(output,
                                       qvaluePlot,
                                       qvalueResult,
                                       omicType,
                                       experimentLabel,
                                       renderTextFn = shiny::renderText,
                                       updateProteinFilteringFn = updateProteinFiltering, projectDirs = NULL) {
  plot_arguments <- list(
    data = qvalueResult$filteringSummary %||% qvalueResult$filteredS4@peptide_data,
    step_name = "2_qval_filtered",
    omic_type = omicType,
    experiment_label = experimentLabel,
    return_grid = TRUE,
    overwrite = TRUE
  )
  if (!is.null(projectDirs)) plot_arguments$project_dirs <- projectDirs
  plotGrid <- if (!is.null(qvalueResult$stateReference) && !is.null(projectDirs)) {
    publishDIAQvalueOutput(getDIAActiveManager(), qvalueResult, plot_arguments,
      updateProteinFilteringFn)
  } else do.call(updateProteinFilteringFn, plot_arguments)
  output$qvalue_results <- renderPeptideQcResultText(qvalueResult$resultText, renderTextFn)
  qvaluePlot(plotGrid)

  invisible(plotGrid)
}

runPeptideQvalueApplyObserver <- function(workflowData,
                                          qvalueThreshold,
                                          globalQvalueThreshold,
                                          proteotypicOnly,
                                          output,
                                          qvaluePlot,
                                          omicType,
                                          experimentLabel,
                                          runApplyStepFn = runPeptideQvalueApplyStep,
                                          updateOutputsFn = updatePeptideQvalueOutputs,
                                          showNotificationFn = shiny::showNotification,
                                          removeNotificationFn = shiny::removeNotification,
                                          logInfoFn = logger::log_info,
                                          logErrorFn = logger::log_error, experimentPaths = NULL) {
  withDIAStateScope(getDIAActiveManager(), function() {
  showNotificationFn(
    "Applying Q-value filter...",
    id = "qvalue_working",
    duration = NULL
  )

  tryCatch({
    qvalueResult <- runApplyStepFn(
      workflowData = workflowData,
      qvalueThreshold = qvalueThreshold,
      globalQvalueThreshold = globalQvalueThreshold,
      proteotypicOnly = proteotypicOnly
    )

    output_arguments <- list(
      output = output,
      qvaluePlot = qvaluePlot,
      qvalueResult = qvalueResult,
      omicType = omicType,
      experimentLabel = experimentLabel
    )
    if (!is.null(experimentPaths)) {
      output_arguments$projectDirs <- stats::setNames(
        list(experimentPaths), paste(omicType, experimentLabel, sep = "_")
      )
    }
    plotGrid <- do.call(updateOutputsFn, output_arguments)

    logInfoFn("Q-value filter applied successfully")
    removeNotificationFn("qvalue_working")
    showNotificationFn("Q-value filter applied successfully", type = "message")

    list(
      status = "success",
      qvalueResult = qvalueResult,
      plotGrid = plotGrid
    )
  }, error = function(e) {
    errorMessage <- paste("Error applying Q-value filter:", e$message)
    logErrorFn(errorMessage)
    showNotificationFn(errorMessage, type = "error", duration = 15)
    removeNotificationFn("qvalue_working")

    list(
      status = "error",
      errorMessage = errorMessage
    )
  })
  })
}

runPeptideQvalueRevertStep <- function(workflowData,
                                       revertStateName = "raw_data_s4",
                                       logInfoFn = logger::log_info) {
  history <- diaStateHistory(workflowData$state_manager)
  if (!(revertStateName %in% history)) {
    stop(sprintf("Cannot revert: '%s' state not found in history.", revertStateName))
  }

  revertedS4 <- selectDIAStateReference(workflowData$state_manager, revertStateName)
  logInfoFn("Reverted to raw data state")

  list(
    revertedS4 = revertedS4,
    revertStateName = revertStateName,
    resultText = "Reverted to raw data state"
  )
}

runPeptideQvalueRevertObserver <- function(workflowData,
                                           output,
                                           runRevertStepFn = runPeptideQvalueRevertStep,
                                           renderTextFn = shiny::renderText,
                                           showNotificationFn = shiny::showNotification,
                                           logErrorFn = logger::log_error) {
  withDIAStateScope(getDIAActiveManager(), function() {
  tryCatch({
    revertResult <- runRevertStepFn(workflowData = workflowData)
    output$qvalue_results <- renderPeptideQcResultText(revertResult$resultText, renderTextFn)
    showNotificationFn(revertResult$resultText, type = "message")

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

mod_prot_qc_peptide_qvalue_server <- function(id, workflow_data, omic_type, experiment_label, experiment_paths = NULL) {
  shiny::moduleServer(id, function(input, output, session) {
    
    qvalue_plot <- shiny::reactiveVal(NULL)
    qvalue_plot_identity <- shiny::reactiveVal("Q-value filtering summary")

    # KIT-014: the managed worker when the session carries a managed application (KIT-008 pattern), else the
    # same registered declaration in process.
    managed_step <- installWorkflowManagedStep(session, input, workflow_data, "peptide_qvalue_filter",
      publication = list(experiment_label = experiment_label, experiment_paths = experiment_paths, omic_type = omic_type),
      ui = list(text = "qvalue_results", plot = "qvalue_plot", ack = "workflow_ack", notification = "qvalue_working",
        label = "Q-value filter", step_label = "2_qval_filtered", plot_identity = qvalue_plot_identity),
      deliver = function(result) {
        arguments <- list(output = output, qvaluePlot = qvalue_plot, qvalueResult = result, omicType = omic_type, experimentLabel = experiment_label)
        if (!is.null(experiment_paths)) arguments$projectDirs <- stats::setNames(list(experiment_paths),
          paste(omic_type, experiment_label, sep = "_"))
        do.call(updatePeptideQvalueOutputs, arguments)
      })
    
    # Step 1: Q-Value Filter (chunk 10)
    shiny::observeEvent(input$apply_qvalue_filter, {
      if (!is.null(workflowSessionApplication(session))) {
        managed_step$apply(list(qvalue_threshold = .diaNnIdentificationQvalueCutoff,
          global_qvalue_threshold = .diaNnIdentificationQvalueCutoff,
          global_pg_qvalue_threshold = .diaNnIdentificationQvalueCutoff,
          proteotypic_only = TRUE, confidence_provenance_mode = "diann_main_report"))
        return(invisible(NULL))
      }
      apply_arguments <- list(
        workflowData = workflow_data,
        qvalueThreshold = .diaNnIdentificationQvalueCutoff,
        globalQvalueThreshold = .diaNnIdentificationQvalueCutoff,
        proteotypicOnly = TRUE,
        output = output,
        qvaluePlot = qvalue_plot,
        omicType = omic_type,
        experimentLabel = experiment_label
      )
      if (!is.null(experiment_paths)) apply_arguments$experimentPaths <- experiment_paths
      do.call(runPeptideQvalueApplyObserver, apply_arguments)
    })
    
    # Revert Q-Value Filter
    shiny::observeEvent(input$revert_qvalue, {
      managed_step$cancel("reverted")
      runPeptideQvalueRevertObserver(
        workflowData = workflow_data,
        output = output
      )
    })
    
    # Render Q-value filter plot
    output$qvalue_plot <- shiny::renderPlot({
      shiny::req(qvalue_plot())
      grid::grid.draw(qvalue_plot())
    }, alt = function() qvalue_plot_identity())
  })
}
