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
#' @param id,workflowData,logInfoFn Runtime inputs used by this function; see the usage section for accepted values.
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
          shiny::p("Retain a protein-group/stripped-peptide identity when any declared replicate group contains at least two distinct observed runs."),
          shiny::hr(),
          
          shiny::textInput(ns("replicate_group_column"), 
            "Replicate Group Column", 
            value = "replicates",
            width = "100%"
          ),
          shiny::helpText("Singleton observations in other groups remain when the feature is supported in at least one group (default column: 'replicates')."),
          
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
  revertedS4 <- selectDIAStateReference(workflowData$state_manager, previousState)
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
  withDIAStateScope(getDIAActiveManager(), function() {
  tryCatch({
    revertResult <- runRevertStepFn(workflowData = workflowData)
    output$replicate_results <- renderPeptideQcResultText(revertResult$resultText, renderTextFn)
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

.peptideReplicateSummaryFromObject <- function(object) {
  if (!base::isS4(object) || !"args" %in% methods::slotNames(object)) {
    return(NULL)
  }
  object@args$removePeptidesWithOnlyOneReplicate$filter_summary
}

runPeptideReplicateApplyStep <- function(workflowData,
                                         replicateGroupColumn,
                                         removePeptidesWithOnlyOneReplicateFn = removePeptidesWithOnlyOneReplicate,
                                         logInfoFn = logger::log_info,
                                         nowFn = Sys.time) {
  if (usesDIAReferenceStates(workflowData$state_manager)) {
    return(runDIAReferencePeptideGroups(workflowData, "removePeptidesWithOnlyOneReplicate",
      list(replicate_group_column = replicateGroupColumn), "replicate_filtered", "Applied distinct-run replicate support filter (supported in any group)"))
  }

  shiny::req(workflowData$state_manager)
  currentS4 <- getDIAExecutionState(workflowData$state_manager, "removePeptidesWithOnlyOneReplicate")
  predecessor <- getDIACurrentReference(workflowData$state_manager)
  shiny::req(currentS4)

  logInfoFn(sprintf(
    "QC Step: Applying replicate filter (column: %s)",
    replicateGroupColumn
  ))

  filteredS4 <- removePeptidesWithOnlyOneReplicateFn(
    currentS4,
    replicate_group_column = replicateGroupColumn
  )
  filterSummary <- .peptideReplicateSummaryFromObject(filteredS4)

  if (is.null(workflowData$qc_params)) {
    workflowData$qc_params <- list()
  }
  if (is.null(workflowData$qc_params$peptide_qc)) {
    workflowData$qc_params$peptide_qc <- list()
  }

  workflowData$qc_params$peptide_qc$replicate_filter <- list(
    replicate_group_column = replicateGroupColumn,
    minimum_distinct_runs = 2L,
    retention_policy = "supported_in_any_group",
    filter_summary = filterSummary,
    timestamp = nowFn()
  )

  filteredS4 <- .savePeptideQcState(
    state_manager = workflowData$state_manager,
    before = currentS4,
    after = filteredS4,
    stage_id = "replicate_filter",
    state_name = "replicate_filtered",
    config_object = list(
      replicate_group_column = replicateGroupColumn,
      minimum_distinct_runs = 2L,
      retention_policy = "supported_in_any_group",
      filter_summary = filterSummary
    ),
    description = "Applied distinct-run replicate support filter (supported in any group)",
    now = nowFn()
  )

  releaseDIAWorkingReference(workflowData$state_manager, predecessor)
  currentS4 <- NULL
  proteinCount <- .countPeptideProteinGroups(filteredS4)

  resultText <- paste(
    "Replicate Filter Applied Successfully\n",
    "====================================\n",
    sprintf("Proteins remaining: %d\n", proteinCount),
    sprintf("Replicate group column: %s\n", replicateGroupColumn),
    "Rule: >= 2 distinct runs in any replicate group\n",
    "Retention: globally keep supported features, including singleton observations in other groups\n",
    if (is.null(filterSummary)) "" else sprintf(
      "Peptide identities retained/removed: %d / %d\n",
      filterSummary$retained_feature_count,
      filterSummary$removed_feature_count
    ),
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
                                          updateProteinFilteringFn = updateProteinFiltering, projectDirs = NULL) {
  output$replicate_results <- renderPeptideQcResultText(replicateResult$resultText, renderTextFn)

  plot_arguments <- list(
    data = replicateResult$filteringSummary %||% replicateResult$filteredS4@peptide_data,
    step_name = "7_replicate_filtered",
    omic_type = omicType,
    experiment_label = experimentLabel,
    return_grid = TRUE,
    overwrite = TRUE
  )
  if (!is.null(projectDirs)) plot_arguments$project_dirs <- projectDirs
  plotGrid <- if (!is.null(replicateResult$stateReference)) {
    publishDIAQvalueOutput(getDIAActiveManager(), replicateResult, plot_arguments, updateProteinFilteringFn)
  } else do.call(updateProteinFilteringFn, plot_arguments)
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
                                             logErrorFn = logger::log_error, experimentPaths = NULL) {
  withDIAStateScope(getDIAActiveManager(), function() {
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

    output_arguments <- list(
      output = output,
      replicatePlot = replicatePlot,
      replicateResult = replicateResult,
      omicType = omicType,
      experimentLabel = experimentLabel
    )
    if (!is.null(experimentPaths)) {
      output_arguments$projectDirs <- stats::setNames(
        list(experimentPaths), paste(omicType, experimentLabel, sep = "_")
      )
    }
    plotGrid <- do.call(updateOutputsFn, output_arguments)

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
  })
}

mod_prot_qc_peptide_replicate_server <- function(id, workflow_data, omic_type, experiment_label, experiment_paths = NULL) {
  shiny::moduleServer(id, function(input, output, session) {
    
    replicate_plot <- shiny::reactiveVal(NULL)
    replicate_plot_identity <- shiny::reactiveVal("Peptide replicate support filtering summary")

    # KIT-014: the managed worker when the session carries a managed application (KIT-008 pattern), else the
    # same registered declaration in process.
    managed_step <- installWorkflowManagedStep(session, input, workflow_data, "peptide_replicate_filter",
      publication = list(experiment_label = experiment_label, experiment_paths = experiment_paths, omic_type = omic_type),
      ui = list(text = "replicate_results", plot = "replicate_plot", ack = "workflow_ack", notification = "replicate_working",
        label = "peptide replicate filter", step_label = "7_replicate_filtered", plot_identity = replicate_plot_identity),
      deliver = function(result) {
        arguments <- list(output = output, replicatePlot = replicate_plot, replicateResult = result, omicType = omic_type, experimentLabel = experiment_label)
        if (!is.null(experiment_paths)) arguments$projectDirs <- stats::setNames(list(experiment_paths),
          paste(omic_type, experiment_label, sep = "_"))
        do.call(updatePeptideReplicateOutputs, arguments)
      })
    
    # Step 6: Replicate Filter (chunk 15)
    shiny::observeEvent(input$apply_replicate_filter, {
      if (!is.null(workflowSessionApplication(session)) && usesDIAReferenceStates(workflow_data$state_manager)) {
        managed_step$apply(list(replicate_group_column = input$replicate_group_column))
        return(invisible(NULL))
      }
      apply_arguments <- list(
        workflowData = workflow_data,
        replicateGroupColumn = input$replicate_group_column,
        output = output,
        replicatePlot = replicate_plot,
        omicType = omic_type,
        experimentLabel = experiment_label
      )
      if (!is.null(experiment_paths)) apply_arguments$experimentPaths <- experiment_paths
      do.call(runPeptideReplicateApplyObserver, apply_arguments)
    })
    
    # Revert Replicate Filter
    shiny::observeEvent(input$revert_replicate, {
      managed_step$cancel("reverted")
      runPeptideReplicateRevertObserver(
        workflowData = workflow_data,
        output = output
      )
    })
    
    # Render replicate filter plot
    output$replicate_plot <- shiny::renderPlot({
      shiny::req(replicate_plot())
      grid::grid.draw(replicate_plot())
    }, alt = function() replicate_plot_identity())
  })
}
