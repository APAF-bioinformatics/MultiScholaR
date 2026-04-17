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
NULL

#' @rdname mod_prot_qc_peptide_qvalue
#' @export
#' @importFrom shiny NS tagList tabPanel br fluidRow column wellPanel h4 p hr numericInput helpText checkboxInput div actionButton verbatimTextOutput plotOutput
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
          
          shiny::numericInput(ns("qvalue_threshold"), 
            "Q-Value Threshold", 
            value = 0.01, min = 0.001, max = 0.1, step = 0.001,
            width = "100%"
          ),
          shiny::helpText("Lower = stricter peptide ID confidence (default: 0.01)"),
          
          shiny::numericInput(ns("global_qvalue_threshold"), 
            "Global Q-Value Threshold", 
            value = 0.01, min = 0.001, max = 0.1, step = 0.001,
            width = "100%"
          ),
          shiny::helpText("Lower = stricter protein group ID confidence (default: 0.01)"),
          
          shiny::checkboxInput(ns("proteotypic_only"), 
            "Keep only proteotypic peptides", 
            value = TRUE
          ),
          shiny::helpText("TRUE = unique peptides only, FALSE = allow shared peptides (default: TRUE)"),
          
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
                                      nowFn = Sys.time) {
  shiny::req(workflowData$state_manager)
  currentS4 <- workflowData$state_manager$getState()
  shiny::req(currentS4)

  logInfoFn(sprintf("Q-value filter: S4 class = %s", class(currentS4)[1]))
  logInfoFn(sprintf(
    "Q-value filter: S4 @peptide_data has %d rows, %d columns",
    nrow(currentS4@peptide_data),
    ncol(currentS4@peptide_data)
  ))

  idsVector <- currentS4@args$srlQvalueProteotypicPeptideClean$input_matrix_column_ids
  if (!is.null(idsVector)) {
    logInfoFn(sprintf(
      "Q-value filter: input_matrix_column_ids length = %d",
      length(idsVector)
    ))
    logInfoFn(sprintf(
      "Q-value filter: input_matrix_column_ids values: [%s]",
      paste(sapply(idsVector, function(x) paste0("'", x, "'")), collapse = ", ")
    ))

    hasWhitespace <- any(grepl("^\\s|\\s$", idsVector))
    hasEmpty <- any(idsVector == "")
    if (hasWhitespace) {
      logWarnFn("Q-value filter: input_matrix_column_ids contains values with leading/trailing whitespace!")
    }
    if (hasEmpty) {
      logWarnFn("Q-value filter: input_matrix_column_ids contains empty strings!")
    }
  }

  logInfoFn(sprintf(
    "QC Step: Applying Q-value filter with thresholds %s, %s",
    qvalueThreshold,
    globalQvalueThreshold
  ))

  currentS4 <- updateConfigParameterFn(
    theObject = currentS4,
    function_name = "srlQvalueProteotypicPeptideClean",
    parameter_name = "qvalue_threshold",
    new_value = qvalueThreshold
  )

  currentS4 <- updateConfigParameterFn(
    theObject = currentS4,
    function_name = "srlQvalueProteotypicPeptideClean",
    parameter_name = "global_qvalue_threshold",
    new_value = globalQvalueThreshold
  )

  currentS4 <- updateConfigParameterFn(
    theObject = currentS4,
    function_name = "srlQvalueProteotypicPeptideClean",
    parameter_name = "choose_only_proteotypic_peptide",
    new_value = as.numeric(proteotypicOnly)
  )

  filteredS4 <- qvalueFilterFn(theObject = currentS4)

  if (is.null(workflowData$qc_params)) {
    workflowData$qc_params <- list()
  }
  if (is.null(workflowData$qc_params$peptide_qc)) {
    workflowData$qc_params$peptide_qc <- list()
  }

  workflowData$qc_params$peptide_qc$qvalue_filter <- list(
    qvalue_threshold = qvalueThreshold,
    global_qvalue_threshold = globalQvalueThreshold,
    proteotypic_only = proteotypicOnly,
    timestamp = nowFn()
  )

  workflowData$state_manager$saveState(
    state_name = "qvalue_filtered",
    s4_data_object = filteredS4,
    config_object = list(
      qvalue_threshold = qvalueThreshold,
      global_qvalue_threshold = globalQvalueThreshold,
      proteotypic_only = proteotypicOnly
    ),
    description = "Applied Q-value and proteotypic peptide filter"
  )

  proteinCount <- filteredS4@peptide_data |>
    dplyr::distinct(Protein.Ids) |>
    nrow()

  resultText <- paste(
    "Q-Value Filter Applied Successfully\n",
    "================================\n",
    sprintf("Proteins remaining: %d\n", proteinCount),
    sprintf("Q-value threshold: %g\n", qvalueThreshold),
    sprintf("Global Q-value threshold: %g\n", globalQvalueThreshold),
    sprintf("Proteotypic only: %s\n", proteotypicOnly),
    "State saved as: 'qvalue_filtered'\n"
  )

  list(
    filteredS4 = filteredS4,
    resultText = resultText
  )
}

updatePeptideQvalueOutputs <- function(output,
                                       qvaluePlot,
                                       qvalueResult,
                                       omicType,
                                       experimentLabel,
                                       renderTextFn = shiny::renderText,
                                       updateProteinFilteringFn = updateProteinFiltering) {
  output$qvalue_results <- renderTextFn(qvalueResult$resultText)

  plotGrid <- updateProteinFilteringFn(
    data = qvalueResult$filteredS4@peptide_data,
    step_name = "2_qval_filtered",
    omic_type = omicType,
    experiment_label = experimentLabel,
    return_grid = TRUE,
    overwrite = TRUE
  )
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
                                          logErrorFn = logger::log_error) {
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

    plotGrid <- updateOutputsFn(
      output = output,
      qvaluePlot = qvaluePlot,
      qvalueResult = qvalueResult,
      omicType = omicType,
      experimentLabel = experimentLabel
    )

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
}

runPeptideQvalueRevertStep <- function(workflowData,
                                       revertStateName = "raw_data_s4",
                                       logInfoFn = logger::log_info) {
  history <- workflowData$state_manager$getHistory()
  if (!(revertStateName %in% history)) {
    stop(sprintf("Cannot revert: '%s' state not found in history.", revertStateName))
  }

  revertedS4 <- workflowData$state_manager$revertToState(revertStateName)
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
  tryCatch({
    revertResult <- runRevertStepFn(workflowData = workflowData)
    output$qvalue_results <- renderTextFn(revertResult$resultText)
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
}

mod_prot_qc_peptide_qvalue_server <- function(id, workflow_data, omic_type, experiment_label) {
  shiny::moduleServer(id, function(input, output, session) {
    
    qvalue_plot <- shiny::reactiveVal(NULL)
    
    # Step 1: Q-Value Filter (chunk 10)
    shiny::observeEvent(input$apply_qvalue_filter, {
      runPeptideQvalueApplyObserver(
        workflowData = workflow_data,
        qvalueThreshold = input$qvalue_threshold,
        globalQvalueThreshold = input$global_qvalue_threshold,
        proteotypicOnly = input$proteotypic_only,
        output = output,
        qvaluePlot = qvalue_plot,
        omicType = omic_type,
        experimentLabel = experiment_label
      )
    })
    
    # Revert Q-Value Filter
    shiny::observeEvent(input$revert_qvalue, {
      runPeptideQvalueRevertObserver(
        workflowData = workflow_data,
        output = output
      )
    })
    
    # Render Q-value filter plot
    output$qvalue_plot <- shiny::renderPlot({
      shiny::req(qvalue_plot())
      grid::grid.draw(qvalue_plot())
    })
  })
}
