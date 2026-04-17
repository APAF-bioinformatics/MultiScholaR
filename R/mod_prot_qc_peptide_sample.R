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

#' @title Peptide Sample Filter Module
#'
#' @description A Shiny module for applying the sample quality filter.
#'
#' @name mod_prot_qc_peptide_sample
NULL

#' @rdname mod_prot_qc_peptide_sample
#' @export
#' @importFrom shiny NS tagList tabPanel br fluidRow column wellPanel h4 p hr numericInput helpText div actionButton verbatimTextOutput plotOutput
mod_prot_qc_peptide_sample_ui <- function(id) {
  ns <- shiny::NS(id)
  
  shiny::tabPanel(
    "Sample Quality",
    shiny::br(),
    shiny::fluidRow(
      shiny::column(4,
        shiny::wellPanel(
          shiny::h4("Minimum Peptides per Sample"),
          shiny::p("Remove samples with insufficient peptide counts (poor sample performance)."),
          shiny::hr(),
          
          shiny::numericInput(ns("min_peptides_per_sample"), 
            "Min Peptides per Sample", 
            value = 500, min = 100, max = 5000, step = 100,
            width = "100%"
          ),
          shiny::helpText("Higher = stricter sample quality filter (default: 500)"),
          
          shiny::hr(),
          shiny::div(
            shiny::actionButton(ns("apply_sample_filter"), "Apply Filter", 
              class = "btn-primary", width = "48%"),
            shiny::actionButton(ns("revert_sample"), "Revert", 
              class = "btn-warning", width = "48%", style = "margin-left: 4%")
          )
        )
      ),
      shiny::column(8,
        shiny::verbatimTextOutput(ns("sample_results")),
        shiny::br(),
        shinyjqui::jqui_resizable(
          shiny::plotOutput(ns("sample_plot"), height = "800px", width = "100%")
        )
      )
    )
  )
}

#' @rdname mod_prot_qc_peptide_sample
#' @export
#' @importFrom shiny moduleServer reactiveVal observeEvent req showNotification removeNotification renderText renderPlot
#' @importFrom logger log_info log_error
#' @importFrom grid grid.draw
runPeptideSampleApplyStep <- function(workflowData,
                                      minPeptidesPerSample,
                                      updateConfigParameterFn = updateConfigParameter,
                                      filterMinNumPeptidesPerSampleFn = filterMinNumPeptidesPerSample,
                                      logInfoFn = logger::log_info,
                                      nowFn = Sys.time) {
  shiny::req(workflowData$state_manager)
  currentS4 <- workflowData$state_manager$getState()
  shiny::req(currentS4)

  logInfoFn(sprintf(
    "QC Step: Applying sample quality filter (min: %s)",
    minPeptidesPerSample
  ))

  currentS4 <- updateConfigParameterFn(
    theObject = currentS4,
    function_name = "filterMinNumPeptidesPerSample",
    parameter_name = "peptides_per_sample_cutoff",
    new_value = minPeptidesPerSample
  )

  samplesBefore <- currentS4@peptide_data |>
    dplyr::distinct(!!rlang::sym(currentS4@sample_id)) |>
    dplyr::pull(!!rlang::sym(currentS4@sample_id))

  filteredS4 <- filterMinNumPeptidesPerSampleFn(theObject = currentS4)

  samplesAfter <- filteredS4@peptide_data |>
    dplyr::distinct(!!rlang::sym(filteredS4@sample_id)) |>
    dplyr::pull(!!rlang::sym(filteredS4@sample_id))

  samplesRemoved <- setdiff(samplesBefore, samplesAfter)
  samplesRemovedCount <- length(samplesRemoved)

  if (is.null(filteredS4@args$filterMinNumPeptidesPerSample)) {
    filteredS4@args$filterMinNumPeptidesPerSample <- list()
  }
  filteredS4@args$filterMinNumPeptidesPerSample$samples_removed <- samplesRemoved
  filteredS4@args$filterMinNumPeptidesPerSample$samples_removed_count <- samplesRemovedCount
  filteredS4@args$filterMinNumPeptidesPerSample$samples_before_count <- length(samplesBefore)
  filteredS4@args$filterMinNumPeptidesPerSample$samples_after_count <- length(samplesAfter)

  if (is.null(workflowData$qc_params)) {
    workflowData$qc_params <- list()
  }
  if (is.null(workflowData$qc_params$peptide_qc)) {
    workflowData$qc_params$peptide_qc <- list()
  }

  workflowData$qc_params$peptide_qc$sample_filter <- list(
    min_peptides_per_sample = minPeptidesPerSample,
    samples_removed = samplesRemoved,
    samples_removed_count = samplesRemovedCount,
    samples_before_count = length(samplesBefore),
    samples_after_count = length(samplesAfter),
    timestamp = nowFn()
  )

  workflowData$state_manager$saveState(
    state_name = "sample_filtered",
    s4_data_object = filteredS4,
    config_object = list(
      min_peptides_per_sample = minPeptidesPerSample,
      samples_removed = samplesRemoved,
      samples_removed_count = samplesRemovedCount
    ),
    description = "Applied minimum peptides per sample filter"
  )

  proteinCount <- filteredS4@peptide_data |>
    dplyr::distinct(Protein.Ids) |>
    nrow()
  runCount <- filteredS4@peptide_data |>
    dplyr::distinct(Run) |>
    nrow()

  resultText <- paste(
    "Sample Quality Filter Applied Successfully\n",
    "=========================================\n",
    sprintf("Proteins remaining: %d\n", proteinCount),
    sprintf("Samples remaining: %d\n", runCount),
    sprintf("Samples removed: %d\n", samplesRemovedCount),
    sprintf("Min peptides per sample: %d\n", minPeptidesPerSample),
    "State saved as: 'sample_filtered'\n"
  )

  if (samplesRemovedCount > 0) {
    resultText <- paste0(
      resultText,
      "\nRemoved samples:\n",
      paste(samplesRemoved, collapse = ", ")
    )
  }

  list(
    filteredS4 = filteredS4,
    resultText = resultText
  )
}

updatePeptideSampleOutputs <- function(output,
                                       samplePlot,
                                       sampleResult,
                                       omicType,
                                       experimentLabel,
                                       renderTextFn = shiny::renderText,
                                       updateProteinFilteringFn = updateProteinFiltering) {
  output$sample_results <- renderTextFn(sampleResult$resultText)

  plotGrid <- updateProteinFilteringFn(
    data = sampleResult$filteredS4@peptide_data,
    step_name = "6_sample_filtered",
    omic_type = omicType,
    experiment_label = experimentLabel,
    return_grid = TRUE,
    overwrite = TRUE
  )
  samplePlot(plotGrid)

  invisible(plotGrid)
}

runPeptideSampleApplyObserver <- function(workflowData,
                                          minPeptidesPerSample,
                                          output,
                                          samplePlot,
                                          omicType,
                                          experimentLabel,
                                          runApplyStepFn = runPeptideSampleApplyStep,
                                          updateOutputsFn = updatePeptideSampleOutputs,
                                          showNotificationFn = shiny::showNotification,
                                          removeNotificationFn = shiny::removeNotification,
                                          logInfoFn = logger::log_info,
                                          logErrorFn = logger::log_error) {
  showNotificationFn(
    "Applying sample quality filter...",
    id = "sample_working",
    duration = NULL
  )

  tryCatch({
    sampleResult <- runApplyStepFn(
      workflowData = workflowData,
      minPeptidesPerSample = minPeptidesPerSample
    )

    plotGrid <- updateOutputsFn(
      output = output,
      samplePlot = samplePlot,
      sampleResult = sampleResult,
      omicType = omicType,
      experimentLabel = experimentLabel
    )

    logInfoFn("Sample quality filter applied successfully")
    removeNotificationFn("sample_working")
    showNotificationFn("Sample quality filter applied successfully", type = "message")

    list(
      status = "success",
      sampleResult = sampleResult,
      plotGrid = plotGrid
    )
  }, error = function(e) {
    errorMessage <- paste("Error applying sample quality filter:", e$message)
    logErrorFn(errorMessage)
    showNotificationFn(errorMessage, type = "error", duration = 15)
    removeNotificationFn("sample_working")

    list(
      status = "error",
      errorMessage = errorMessage
    )
  })
}

runPeptideSampleRevertStep <- function(workflowData,
                                       logInfoFn = logger::log_info) {
  history <- workflowData$state_manager$getHistory()
  if (length(history) <= 1) {
    stop("No previous state to revert to.")
  }

  previousState <- history[length(history) - 1]
  revertedS4 <- workflowData$state_manager$revertToState(previousState)
  logInfoFn(paste("Reverted sample filter to", previousState))

  list(
    previousState = previousState,
    revertedS4 = revertedS4,
    resultText = paste("Reverted to previous state:", previousState)
  )
}

runPeptideSampleRevertObserver <- function(workflowData,
                                           output,
                                           runRevertStepFn = runPeptideSampleRevertStep,
                                           renderTextFn = shiny::renderText,
                                           showNotificationFn = shiny::showNotification,
                                           logErrorFn = logger::log_error) {
  tryCatch({
    revertResult <- runRevertStepFn(workflowData = workflowData)
    output$sample_results <- renderTextFn(revertResult$resultText)
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

mod_prot_qc_peptide_sample_server <- function(id, workflow_data, omic_type, experiment_label) {
  shiny::moduleServer(id, function(input, output, session) {
    
    sample_plot <- shiny::reactiveVal(NULL)
    
    # Step 5: Sample Quality Filter (chunk 14)
    shiny::observeEvent(input$apply_sample_filter, {
      runPeptideSampleApplyObserver(
        workflowData = workflow_data,
        minPeptidesPerSample = input$min_peptides_per_sample,
        output = output,
        samplePlot = sample_plot,
        omicType = omic_type,
        experimentLabel = experiment_label
      )
    })
    
    # Revert Sample Filter
    shiny::observeEvent(input$revert_sample, {
      runPeptideSampleRevertObserver(
        workflowData = workflow_data,
        output = output
      )
    })
    
    # Render sample filter plot
    output$sample_plot <- shiny::renderPlot({
      shiny::req(sample_plot())
      grid::grid.draw(sample_plot())
    })
  })
}
