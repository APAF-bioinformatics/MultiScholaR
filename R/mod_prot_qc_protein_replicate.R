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

#' @title Protein Replicate Filter Module
#'
#' @description A Shiny module for applying the protein replicate filter.
#'
#' @name mod_prot_qc_protein_replicate
#' @param id,workflowData,experimentPaths,groupingVariable,parallelCores,removeProteinsWithOnlyOneReplicateFn,writeTsvFn,saveRdsFn,logInfoFn,logWarnFn,newClusterFn,isTestModeFn,nowFn Runtime inputs used by this function; see the usage section for accepted values.
NULL

#' @rdname mod_prot_qc_protein_replicate
#' @export
#' @importFrom shiny NS tagList tabPanel br fluidRow column wellPanel h4 p hr textInput helpText numericInput div actionButton verbatimTextOutput plotOutput
mod_prot_qc_protein_replicate_ui <- function(id) {
  ns <- shiny::NS(id)
  
  shiny::tabPanel(
    "Protein Replicate Filter",
    shiny::br(),
    shiny::fluidRow(
      shiny::column(4,
        shiny::wellPanel(
          shiny::h4("Remove Single-Replicate Proteins"),
          shiny::p("Remove proteins detected in only one replicate across experimental groups."),
          shiny::hr(),
          
          shiny::textInput(ns("protein_grouping_variable"), 
            "Grouping Variable", 
            value = "group",
            width = "100%"
          ),
          shiny::helpText("Column name for experimental grouping (default: 'group')"),
          
          shiny::numericInput(ns("parallel_cores"), 
            "Parallel Processing Cores", 
            value = 4, min = 1, max = 8, step = 1,
            width = "100%"
          ),
          shiny::helpText("Number of cores for parallel processing (default: 4)"),
          
          shiny::hr(),
          shiny::div(
            shiny::actionButton(ns("apply_protein_replicate_filter"), "Apply Filter", 
              class = "btn-primary", width = "48%"),
            shiny::actionButton(ns("revert_protein_replicate_filter"), "Revert", 
              class = "btn-warning", width = "48%", style = "margin-left: 4%")
          )
        )
      ),
      shiny::column(8,
        shiny::verbatimTextOutput(ns("protein_replicate_filter_results")),
        shiny::br(),
        shinyjqui::jqui_resizable(
          shiny::plotOutput(ns("protein_replicate_filter_plot"), height = "800px", width = "100%")
        )
      )
    )
  )
}

#' @rdname mod_prot_qc_protein_replicate
#' @export
#' @importFrom shiny moduleServer reactiveVal observeEvent req showNotification removeNotification renderText renderPlot
#' @importFrom logger log_info log_error
#' @importFrom grid grid.draw
runProteinReplicateFilterApplyStep <- function(workflowData,
                                               experimentPaths,
                                               groupingVariable,
                                               parallelCores,
                                               removeProteinsWithOnlyOneReplicateFn = removeProteinsWithOnlyOneReplicate,
                                               writeTsvFn = vroom::vroom_write,
                                               saveRdsFn = saveRDS,
                                               logInfoFn = logger::log_info,
                                               logWarnFn = logger::log_warn,
                                               newClusterFn = if (exists("new_cluster")) get("new_cluster") else NULL,
                                               isTestModeFn = is_test_mode,
                                               nowFn = Sys.time) {
  pending <- if (identical(workflowData$state_manager$current_state, "protein_replicate_filtered"))
    resumeDIAProteinAction(workflowData, "protein_replicate_filtered") else NULL
  if (!is.null(pending)) return(pending)
  # KIT-015: an unmanaged launch runs the KIT-010 declaration in process with its production variant (`rows`).
  entry <- workflowActionSpec("protein_replicate_filter")
  parameters <- list(grouping_variable = groupingVariable, parallel_cores = parallelCores)
  validateWorkflowActionParameters(entry$name, parameters, workflowExecutionPolicy(FALSE), resolved = FALSE)
  value <- computeDIAReferenceProteinReplicateFilter(workflowData$state_manager, getDIAStateReference(workflowData$state_manager),
    parameters, variant = entry$default_variant)
  qc_params <- shiny::isolate(workflowData$qc_params) %||% list()
  if (is.null(qc_params$protein_qc)) qc_params$protein_qc <- list()
  qc_params$protein_qc$replicate_filter <- list(grouping_variable = groupingVariable, parallel_cores = parallelCores,
    timestamp = nowFn(), effective_parallel_cores = 1L)
  workflowData$qc_params <- qc_params
  return(completeDIAProteinAction(workflowData, value, publish = TRUE))
}

updateProteinReplicateFilterOutputs <- function(output,
                                                proteinReplicateFilterPlot,
                                                replicateFilterResult,
                                                omicType,
                                                experimentLabel,
                                                renderTextFn = shiny::renderText,
                                                updateProteinFilteringFn = updateProteinFiltering) {
  output$protein_replicate_filter_results <- renderPeptideQcResultText(replicateFilterResult$resultText, renderTextFn)

  plotGrid <- publishDIAProteinProgress(replicateFilterResult, "13_protein_replicate_filtered", omicType,
    experimentLabel, updateProteinFilteringFn, replicateFilterResult$filteredS4)
  proteinReplicateFilterPlot(plotGrid)

  invisible(plotGrid)
}

runProteinReplicateFilterApplyObserver <- function(workflowData,
                                                   experimentPaths,
                                                   groupingVariable,
                                                   parallelCores,
                                                   output,
                                                   proteinReplicateFilterPlot,
                                                   omicType,
                                                   experimentLabel,
                                                   runApplyStepFn = runProteinReplicateFilterApplyStep,
                                                   updateOutputsFn = updateProteinReplicateFilterOutputs,
                                                   completeQcStatusFn = completeProtQcWorkflowStatus,
                                                   showNotificationFn = shiny::showNotification,
                                                   removeNotificationFn = shiny::removeNotification,
                                                   logInfoFn = logger::log_info,
                                                   logErrorFn = logger::log_error) {
  withDIAStateScope(getDIAActiveManager(), function() {
  currentState <- tryCatch(workflowData$state_manager$current_state, error = function(e) NULL)
  if (identical(currentState, "protein_replicate_filtered")) {
    completeQcStatusFn(workflowData = workflowData)
    return(invisible(list(status = "skipped", reason = "already_complete")))
  }
  if (isTRUE(workflowData$protein_replicate_filter_in_progress)) {
    return(invisible(list(status = "skipped", reason = "in_progress")))
  }
  workflowData$protein_replicate_filter_in_progress <- TRUE
  on.exit({
    workflowData$protein_replicate_filter_in_progress <- FALSE
  }, add = TRUE)

  showNotificationFn(
    "Applying protein replicate filter...",
    id = "protein_replicate_filter_working",
    duration = NULL
  )

  tryCatch({
    replicateFilterResult <- runApplyStepFn(
      workflowData = workflowData,
      experimentPaths = experimentPaths,
      groupingVariable = groupingVariable,
      parallelCores = parallelCores
    )

    plotGrid <- updateOutputsFn(
      output = output,
      proteinReplicateFilterPlot = proteinReplicateFilterPlot,
      replicateFilterResult = replicateFilterResult,
      omicType = omicType,
      experimentLabel = experimentLabel
    )

    completeQcStatusFn(workflowData = workflowData)

    logInfoFn("Protein replicate filter applied successfully")
    removeNotificationFn("protein_replicate_filter_working")
    showNotificationFn(
      "Protein replicate filter applied successfully. Pipeline complete!",
      type = "message"
    )

    list(
      status = "success",
      replicateFilterResult = replicateFilterResult,
      plotGrid = plotGrid
    )
  }, error = function(e) {
    errorMessage <- paste("Error applying protein replicate filter:", e$message)
    logErrorFn(errorMessage)
    showNotificationFn(errorMessage, type = "error", duration = 15)
    removeNotificationFn("protein_replicate_filter_working")

    list(
      status = "error",
      errorMessage = errorMessage
    )
  })
  })
}

runProteinReplicateFilterRevertStep <- function(workflowData) {
  shiny::req(workflowData$state_manager)
  history <- diaStateHistory(workflowData$state_manager)

  if (length(history) <= 1) {
    stop("No previous state to revert to.")
  }

  previousState <- history[length(history) - 1]
  revertedS4 <- selectDIAStateReference(workflowData$state_manager, previousState)

  list(
    previousState = previousState,
    revertedS4 = revertedS4,
    resultText = paste("Reverted to previous state:", previousState)
  )
}

runProteinReplicateFilterRevertObserver <- function(workflowData,
                                                    output,
                                                    runRevertStepFn = runProteinReplicateFilterRevertStep,
                                                    renderTextFn = shiny::renderText,
                                                    showNotificationFn = shiny::showNotification,
                                                    logInfoFn = logger::log_info,
                                                    logErrorFn = logger::log_error) {
  withDIAStateScope(getDIAActiveManager(), function() {
  tryCatch({
    revertResult <- runRevertStepFn(workflowData = workflowData)
    output$protein_replicate_filter_results <- renderPeptideQcResultText(revertResult$resultText, renderTextFn)
    logInfoFn(paste("Reverted protein replicate filter to", revertResult$previousState))
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

bindProteinReplicateFilterPlot <- function(output, proteinReplicateFilterPlot, plotIdentity = NULL) {
  output$protein_replicate_filter_plot <- if (is.null(plotIdentity)) renderPlot({
    req(proteinReplicateFilterPlot())
    grid.draw(proteinReplicateFilterPlot())
  }) else renderPlot({
    req(proteinReplicateFilterPlot())
    grid.draw(proteinReplicateFilterPlot())
  }, alt = function() plotIdentity())
}

mod_prot_qc_protein_replicate_server <- function(id, workflow_data, experiment_paths, omic_type, experiment_label) {
  shiny::moduleServer(id, function(input, output, session) {
    logger::log_info(sprintf(
      "Protein replicate filter server registered for input: %s",
      session$ns("apply_protein_replicate_filter")
    ))
    
    protein_replicate_filter_plot <- shiny::reactiveVal(NULL)
    protein_replicate_filter_plot_identity <- shiny::reactiveVal("Protein replicate filtering summary")
    protein_replicate_filter_plot_value <- protein_replicate_filter_plot

    # KIT-015: the managed worker when the session carries a managed application (KIT-008 pattern), else the
    # same registered declaration in process.
    managed_step <- installWorkflowManagedStep(session, input, workflow_data, "protein_replicate_filter",
      publication = list(experiment_label = experiment_label, experiment_paths = experiment_paths, omic_type = omic_type),
      ui = list(text = "protein_replicate_filter_results", plot = "protein_replicate_filter_plot", ack = "workflow_ack", notification = "protein_replicate_filter_working",
        label = "protein replicate filter", step_label = "13_protein_replicate_filtered", plot_identity = protein_replicate_filter_plot_identity),
      deliver = function(result) updateProteinReplicateFilterOutputs(output = output, proteinReplicateFilterPlot = protein_replicate_filter_plot_value, replicateFilterResult = result,
        omicType = omic_type, experimentLabel = experiment_label),
      read_result = readManagedProteinResult)
    runApplyObserver <- function(groupingVariable, parallelCores) {
      if (!is.null(workflowSessionApplication(session))) {
        managed_step$apply(list(grouping_variable = groupingVariable %||% "group", parallel_cores = parallelCores %||% 4))
        return(invisible(NULL))
      }
      runProteinReplicateFilterApplyObserver(
        workflowData = workflow_data,
        experimentPaths = experiment_paths,
        groupingVariable = groupingVariable %||% "group",
        parallelCores = parallelCores %||% 4,
        output = output,
        proteinReplicateFilterPlot = protein_replicate_filter_plot,
        omicType = omic_type,
        experimentLabel = experiment_label
      )
    }

    # Step 5: Protein Replicate Filter (chunk 23)
    shiny::observeEvent(input$apply_protein_replicate_filter, {
      logger::log_info("Protein replicate filter apply observed via module input")
      runApplyObserver(
        groupingVariable = input$protein_grouping_variable,
        parallelCores = input$parallel_cores
      )
    })

    rootScopeFn <- session$rootScope
    rootSession <- if (is.function(rootScopeFn)) rootScopeFn() else session
    if (!identical(rootSession, session)) {
      shiny::observeEvent(rootSession$input[[session$ns("apply_protein_replicate_filter")]], {
        logger::log_info("Protein replicate filter apply observed via root input")
        runApplyObserver(
          groupingVariable = rootSession$input[[session$ns("protein_grouping_variable")]] %||% input$protein_grouping_variable,
          parallelCores = rootSession$input[[session$ns("parallel_cores")]] %||% input$parallel_cores
        )
      }, ignoreNULL = TRUE)
    }
    
    # Revert Protein Replicate Filter
    shiny::observeEvent(input$revert_protein_replicate_filter, {
      managed_step$cancel("reverted")
      runProteinReplicateFilterRevertObserver(
        workflowData = workflow_data,
        output = output
      )
    })

    bindProteinReplicateFilterPlot(
      output = output,
      proteinReplicateFilterPlot = protein_replicate_filter_plot,
      plotIdentity = protein_replicate_filter_plot_identity
    )
  })
}
