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

#' @title Create Protein S4 Module
#'
#' @description UI component for creating a ProteinQuantitativeData S4 object
#' from protein-level data. This provides a user-facing button to trigger
#' the process for workflows like TMT.
#'
#' @name mod_prot_qc_protein_s4
NULL

#' @rdname mod_prot_qc_protein_s4
#' @export
#' @importFrom shiny NS tagList fluidPage wellPanel h4 p actionButton hr h5 verbatimTextOutput
mod_prot_qc_protein_s4_ui <- function(id) {
  ns <- shiny::NS(id)
  
  shiny::tagList(
    shiny::fluidPage(
    shiny::wellPanel(
      shiny::h4("Finalise Protein Data"),
      shiny::p("This workflow starts with protein-level data. Click the button below to format the data into an S4 object for downstream analysis."),
      shiny::actionButton(ns("create_protein_s4"), "Create Protein S4 Object", class = "btn-primary"),
      shiny::hr(),
      shiny::h5("Results"),
      shiny::verbatimTextOutput(ns("s4_creation_results"))
    )
  )
  )
}

runProteinS4CreationStep <- function(workflowData,
                                     proteinQuantitativeDataFn = ProteinQuantitativeData,
                                     logInfoFn = logger::log_info) {
  shiny::req(
    workflowData$data_cln,
    workflowData$design_matrix,
    workflowData$column_mapping,
    workflowData$config_list
  )

  logInfoFn("Protein S4 Creation: Starting process for protein-level workflow (e.g., TMT)")

  proteinIdCol <- workflowData$column_mapping$protein_col
  if (!proteinIdCol %in% names(workflowData$data_cln)) {
    stop(paste("Protein ID column", proteinIdCol, "not found in the data."))
  }

  proteinObj <- proteinQuantitativeDataFn(
    protein_quant_table = workflowData$data_cln,
    protein_id_column = proteinIdCol,
    protein_id_table = workflowData$data_cln |>
      dplyr::distinct(!!rlang::sym(proteinIdCol)),
    design_matrix = workflowData$design_matrix,
    sample_id = "Run",
    group_id = "group",
    technical_replicate_id = "replicates",
    args = workflowData$config_list
  )

  workflowData$state_manager$saveState(
    state_name = "protein_s4_initial",
    s4_data_object = proteinObj,
    config_object = list(
      s4_class = "ProteinQuantitativeData",
      protein_id_column = proteinIdCol
    ),
    description = "Created initial ProteinQuantitativeData S4 object from protein-level data."
  )

  proteinCount <- proteinObj@protein_quant_table |>
    dplyr::distinct(!!rlang::sym(proteinIdCol)) |>
    nrow()

  resultText <- paste(
    "Protein S4 Object Created Successfully\n",
    "======================================\n",
    sprintf("Proteins loaded: %d\n", proteinCount),
    sprintf("S4 Class: %s\n", class(proteinObj)[1]),
    "State saved as: 'protein_s4_initial'\n",
    "\nReady for protein accession cleanup."
  )

  list(
    proteinObj = proteinObj,
    proteinCount = proteinCount,
    proteinIdCol = proteinIdCol,
    resultText = resultText
  )
}

runProteinS4CreationObserver <- function(workflowData,
                                         output,
                                         runCreationStepFn = runProteinS4CreationStep,
                                         renderTextFn = shiny::renderText,
                                         showNotificationFn = shiny::showNotification,
                                         removeNotificationFn = shiny::removeNotification,
                                         logInfoFn = logger::log_info,
                                         logErrorFn = logger::log_error) {
  showNotificationFn(
    "Creating Protein S4 object from imported data...",
    id = "s4_creation_working",
    duration = NULL
  )

  tryCatch({
    creationResult <- runCreationStepFn(workflowData = workflowData)
    output$s4_creation_results <- renderTextFn(creationResult$resultText)

    logInfoFn("Protein S4 object creation from protein data completed successfully")
    removeNotificationFn("s4_creation_working")
    showNotificationFn("Protein S4 object created successfully", type = "message")

    list(
      status = "success",
      creationResult = creationResult
    )
  }, error = function(e) {
    msg <- paste("Error creating Protein S4 object:", e$message)
    logErrorFn(msg)
    showNotificationFn(msg, type = "error", duration = 15)
    removeNotificationFn("s4_creation_working")

    list(
      status = "error",
      errorMessage = msg
    )
  })
}

runProteinS4RevertStep <- function(workflowData,
                                   logInfoFn = logger::log_info) {
  shiny::req(workflowData$state_manager)

  revertedState <- workflowData$state_manager$revertToState("initial")
  logInfoFn("Reverted S4 object creation.")

  list(
    revertedState = revertedState,
    stateName = "initial",
    resultText = "Reverted to initial empty state. You may need to re-run previous steps."
  )
}

runProteinS4RevertObserver <- function(workflowData,
                                       output,
                                       runRevertStepFn = runProteinS4RevertStep,
                                       renderTextFn = shiny::renderText,
                                       showNotificationFn = shiny::showNotification,
                                       logErrorFn = logger::log_error) {
  tryCatch({
    revertResult <- runRevertStepFn(workflowData = workflowData)
    output$s4_creation_results <- renderTextFn(revertResult$resultText)
    showNotificationFn("Reverted successfully", type = "message")

    list(
      status = "success",
      revertResult = revertResult
    )
  }, error = function(e) {
    msg <- paste("Error reverting:", e$message)
    logErrorFn(msg)
    showNotificationFn(msg, type = "error")

    list(
      status = "error",
      errorMessage = msg
    )
  })
}

#' @rdname mod_prot_qc_protein_s4
#' @export
#' @importFrom shiny moduleServer observeEvent req showNotification removeNotification renderText
#' @importFrom logger log_info log_error
mod_prot_qc_protein_s4_server <- function(id, workflow_data, omic_type, experiment_label) {
  shiny::moduleServer(id, function(input, output, session) {
    
    # This module is triggered by an action button in the UI.
    shiny::observeEvent(input$create_protein_s4, {
      runProteinS4CreationObserver(
        workflowData = workflow_data,
        output = output
      )
    })
    
    # Revert S4 Creation (if needed)
    # Note: The original code had this, although there might not be a previous state to revert to if this is the first step.
    # We'll keep it but make it revert to initial.
    shiny::observeEvent(input$revert_s4_creation, {
      runProteinS4RevertObserver(
        workflowData = workflow_data,
        output = output
      )
    })
    
  })
}
