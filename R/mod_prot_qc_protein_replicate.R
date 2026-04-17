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
                                               nowFn = Sys.time) {
  shiny::req(workflowData$state_manager)
  currentS4 <- workflowData$state_manager$getState()
  shiny::req(currentS4)

  logInfoFn(sprintf(
    "Protein Processing: Applying protein replicate filter with %d cores",
    parallelCores
  ))

  coreUtilisation <- if (is.function(newClusterFn)) {
    newClusterFn(parallelCores)
  } else {
    NULL
  }

  filteredS4 <- removeProteinsWithOnlyOneReplicateFn(
    currentS4,
    coreUtilisation,
    grouping_variable = groupingVariable
  )

  outputFile <- if (!is.null(experimentPaths$protein_qc_dir)) {
    file.path(
      experimentPaths$protein_qc_dir,
      "remove_proteins_with_only_one_rep.tsv"
    )
  } else {
    "remove_proteins_with_only_one_rep.tsv"
  }
  writeTsvFn(filteredS4@protein_quant_table, outputFile)

  if (is.null(workflowData$qc_params)) {
    workflowData$qc_params <- list()
  }
  if (is.null(workflowData$qc_params$protein_qc)) {
    workflowData$qc_params$protein_qc <- list()
  }

  workflowData$qc_params$protein_qc$replicate_filter <- list(
    grouping_variable = groupingVariable,
    parallel_cores = parallelCores,
    timestamp = nowFn()
  )

  tryCatch({
    if (!is.null(experimentPaths$source_dir)) {
      qcParamsFile <- file.path(experimentPaths$source_dir, "qc_params.RDS")
      saveRdsFn(workflowData$qc_params, qcParamsFile)
      logInfoFn(sprintf("Saved QC parameters to: %s", qcParamsFile))
    }
  }, error = function(e) {
    logWarnFn(sprintf("Could not save QC parameters file: %s", e$message))
  })

  workflowData$state_manager$saveState(
    state_name = "protein_replicate_filtered",
    s4_data_object = filteredS4,
    config_object = list(
      grouping_variable = groupingVariable,
      parallel_cores = parallelCores,
      output_file = outputFile
    ),
    description = "Applied protein replicate filter (removed single-replicate proteins)"
  )

  proteinCount <- filteredS4@protein_quant_table |>
    dplyr::distinct(Protein.Ids) |>
    nrow()

  if (is.null(workflowData$protein_counts)) {
    workflowData$protein_counts <- list()
  }
  workflowData$protein_counts$after_qc_filtering <- proteinCount
  logInfoFn(sprintf("Tracked protein count after QC filtering: %d", proteinCount))

  resultText <- paste(
    "Protein Replicate Filter Applied Successfully\n",
    "============================================\n",
    sprintf("Proteins remaining: %d\n", proteinCount),
    sprintf("Grouping variable: %s\n", groupingVariable),
    sprintf("Parallel cores used: %d\n", parallelCores),
    sprintf("Output file: %s\n", basename(outputFile)),
    "State saved as: 'protein_replicate_filtered'\n",
    "\nProtein filtering pipeline complete!"
  )

  list(
    filteredS4 = filteredS4,
    outputFile = outputFile,
    proteinCount = proteinCount,
    resultText = resultText
  )
}

updateProteinReplicateFilterOutputs <- function(output,
                                                proteinReplicateFilterPlot,
                                                replicateFilterResult,
                                                omicType,
                                                experimentLabel,
                                                renderTextFn = shiny::renderText,
                                                updateProteinFilteringFn = updateProteinFiltering) {
  output$protein_replicate_filter_results <- renderTextFn(replicateFilterResult$resultText)

  plotGrid <- updateProteinFilteringFn(
    data = replicateFilterResult$filteredS4@protein_quant_table,
    step_name = "13_protein_replicate_filtered",
    omic_type = omicType,
    experiment_label = experimentLabel,
    return_grid = TRUE,
    overwrite = TRUE
  )
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
                                                   showNotificationFn = shiny::showNotification,
                                                   removeNotificationFn = shiny::removeNotification,
                                                   logInfoFn = logger::log_info,
                                                   logErrorFn = logger::log_error) {
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
}

runProteinReplicateFilterRevertStep <- function(workflowData) {
  shiny::req(workflowData$state_manager)
  history <- workflowData$state_manager$getHistory()

  if (length(history) <= 1) {
    stop("No previous state to revert to.")
  }

  previousState <- history[length(history) - 1]
  revertedS4 <- workflowData$state_manager$revertToState(previousState)

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
  tryCatch({
    revertResult <- runRevertStepFn(workflowData = workflowData)
    output$protein_replicate_filter_results <- renderTextFn(revertResult$resultText)
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
}

bindProteinReplicateFilterPlot <- function(output, proteinReplicateFilterPlot) {
  output$protein_replicate_filter_plot <- renderPlot({
    req(proteinReplicateFilterPlot())
    grid.draw(proteinReplicateFilterPlot())
  })
}

mod_prot_qc_protein_replicate_server <- function(id, workflow_data, experiment_paths, omic_type, experiment_label) {
  shiny::moduleServer(id, function(input, output, session) {
    
    protein_replicate_filter_plot <- shiny::reactiveVal(NULL)
    
    # Step 5: Protein Replicate Filter (chunk 23)
    shiny::observeEvent(input$apply_protein_replicate_filter, {
      runProteinReplicateFilterApplyObserver(
        workflowData = workflow_data,
        experimentPaths = experiment_paths,
        groupingVariable = input$protein_grouping_variable,
        parallelCores = input$parallel_cores,
        output = output,
        proteinReplicateFilterPlot = protein_replicate_filter_plot,
        omicType = omic_type,
        experimentLabel = experiment_label
      )
    })
    
    # Revert Protein Replicate Filter
    shiny::observeEvent(input$revert_protein_replicate_filter, {
      runProteinReplicateFilterRevertObserver(
        workflowData = workflow_data,
        output = output
      )
    })

    bindProteinReplicateFilterPlot(
      output = output,
      proteinReplicateFilterPlot = protein_replicate_filter_plot
    )
  })
}
