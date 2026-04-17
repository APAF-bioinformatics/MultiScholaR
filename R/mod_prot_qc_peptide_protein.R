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

#' @title Protein Peptide Count Filter Module
#'
#' @description A Shiny module for applying the protein peptide count filter.
#'
#' @name mod_prot_qc_peptide_protein
NULL

#' @rdname mod_prot_qc_peptide_protein
#' @export
#' @importFrom shiny NS tagList tabPanel br fluidRow column wellPanel h4 p hr numericInput helpText div actionButton verbatimTextOutput plotOutput
mod_prot_qc_peptide_protein_ui <- function(id) {
  ns <- shiny::NS(id)
  
  shiny::tabPanel(
    "Protein Peptides",
    shiny::br(),
    shiny::fluidRow(
      shiny::column(4,
        shiny::wellPanel(
          shiny::h4("Minimum Peptides per Protein"),
          shiny::p("Keep proteins only if they have sufficient peptide evidence (two-peptide rule)."),
          shiny::hr(),
          
          shiny::numericInput(ns("min_peptides_per_protein"), 
            "Min Peptides per Protein", 
            value = 2, min = 1, max = 5, step = 1,
            width = "100%"
          ),
          shiny::helpText("Higher = requires more unique peptides per protein (default: 2)"),
          
          shiny::numericInput(ns("min_peptidoforms_per_protein"), 
            "Min Peptidoforms per Protein", 
            value = 2, min = 1, max = 5, step = 1,
            width = "100%"
          ),
          shiny::helpText("Higher = requires more peptide forms per protein (default: 2)"),
          
          shiny::hr(),
          shiny::div(
            shiny::actionButton(ns("apply_protein_peptide_filter"), "Apply Filter", 
              class = "btn-primary", width = "48%"),
            shiny::actionButton(ns("revert_protein_peptide"), "Revert", 
              class = "btn-warning", width = "48%", style = "margin-left: 4%")
          )
        )
      ),
      shiny::column(8,
        shiny::verbatimTextOutput(ns("protein_peptida_results")),
        shiny::br(),
        shinyjqui::jqui_resizable(
          shiny::plotOutput(ns("protein_peptide_plot"), height = "800px", width = "100%")
        )
      )
    )
  )
}

#' @rdname mod_prot_qc_peptide_protein
#' @export
#' @importFrom shiny moduleServer reactiveVal observeEvent req showNotification removeNotification renderText renderPlot
#' @importFrom logger log_info log_error
#' @importFrom grid grid.draw
runProteinPeptideApplyStep <- function(workflowData,
                                       minPeptidesPerProtein,
                                       minPeptidoformsPerProtein,
                                       updateConfigParameterFn = updateConfigParameter,
                                       filterMinNumPeptidesPerProteinFn = filterMinNumPeptidesPerProtein,
                                       logInfoFn = logger::log_info,
                                       nowFn = Sys.time) {
  shiny::req(workflowData$state_manager)
  currentS4 <- workflowData$state_manager$getState()
  shiny::req(currentS4)

  logInfoFn(sprintf(
    "QC Step: Applying protein peptide count filter (min: %s)",
    minPeptidesPerProtein
  ))

  currentS4 <- updateConfigParameterFn(
    theObject = currentS4,
    function_name = "filterMinNumPeptidesPerProtein",
    parameter_name = "peptides_per_protein_cutoff",
    new_value = minPeptidesPerProtein
  )

  currentS4 <- updateConfigParameterFn(
    theObject = currentS4,
    function_name = "filterMinNumPeptidesPerProtein",
    parameter_name = "peptidoforms_per_protein_cutoff",
    new_value = minPeptidoformsPerProtein
  )

  filteredS4 <- filterMinNumPeptidesPerProteinFn(theObject = currentS4)

  if (is.null(workflowData$qc_params)) {
    workflowData$qc_params <- list()
  }
  if (is.null(workflowData$qc_params$peptide_qc)) {
    workflowData$qc_params$peptide_qc <- list()
  }

  workflowData$qc_params$peptide_qc$protein_peptide_filter <- list(
    min_peptides_per_protein = minPeptidesPerProtein,
    min_peptidoforms_per_protein = minPeptidoformsPerProtein,
    timestamp = nowFn()
  )

  workflowData$state_manager$saveState(
    state_name = "protein_peptide_filtered",
    s4_data_object = filteredS4,
    config_object = list(
      min_peptides_per_protein = minPeptidesPerProtein,
      min_peptidoforms_per_protein = minPeptidoformsPerProtein
    ),
    description = "Applied minimum peptides per protein filter"
  )

  proteinCount <- filteredS4@peptide_data |>
    dplyr::distinct(Protein.Ids) |>
    nrow()

  resultText <- paste(
    "Protein Peptide Count Filter Applied Successfully\n",
    "===============================================\n",
    sprintf("Proteins remaining: %d\n", proteinCount),
    sprintf("Min peptides per protein: %d\n", minPeptidesPerProtein),
    sprintf("Min peptidoforms per protein: %d\n", minPeptidoformsPerProtein),
    "State saved as: 'protein_peptide_filtered'\n"
  )

  list(
    filteredS4 = filteredS4,
    resultText = resultText
  )
}

updateProteinPeptideOutputs <- function(output,
                                        proteinPeptidePlot,
                                        proteinPeptideResult,
                                        omicType,
                                        experimentLabel,
                                        renderTextFn = renderText,
                                        updateProteinFilteringFn = updateProteinFiltering) {
  output$protein_peptida_results <- renderTextFn(proteinPeptideResult$resultText)

  plotGrid <- updateProteinFilteringFn(
    data = proteinPeptideResult$filteredS4@peptide_data,
    step_name = "5_protein_peptide_filtered",
    omic_type = omicType,
    experiment_label = experimentLabel,
    return_grid = TRUE,
    overwrite = TRUE
  )
  proteinPeptidePlot(plotGrid)

  invisible(plotGrid)
}

runProteinPeptideApplyObserver <- function(workflowData,
                                           minPeptidesPerProtein,
                                           minPeptidoformsPerProtein,
                                           output,
                                           proteinPeptidePlot,
                                           omicType,
                                           experimentLabel,
                                           runApplyStepFn = runProteinPeptideApplyStep,
                                           updateOutputsFn = updateProteinPeptideOutputs,
                                           showNotificationFn = shiny::showNotification,
                                           removeNotificationFn = shiny::removeNotification,
                                           logInfoFn = logger::log_info,
                                           logErrorFn = logger::log_error) {
  showNotificationFn(
    "Applying protein peptide count filter...",
    id = "protein_peptide_working",
    duration = NULL
  )

  tryCatch({
    proteinPeptideResult <- runApplyStepFn(
      workflowData = workflowData,
      minPeptidesPerProtein = minPeptidesPerProtein,
      minPeptidoformsPerProtein = minPeptidoformsPerProtein
    )

    plotGrid <- updateOutputsFn(
      output = output,
      proteinPeptidePlot = proteinPeptidePlot,
      proteinPeptideResult = proteinPeptideResult,
      omicType = omicType,
      experimentLabel = experimentLabel
    )

    logInfoFn("Protein peptide count filter applied successfully")
    removeNotificationFn("protein_peptide_working")
    showNotificationFn(
      "Protein peptide count filter applied successfully",
      type = "message"
    )

    list(
      status = "success",
      proteinPeptideResult = proteinPeptideResult,
      plotGrid = plotGrid
    )
  }, error = function(e) {
    errorMessage <- paste("Error applying protein peptide count filter:", e$message)
    logErrorFn(errorMessage)
    showNotificationFn(errorMessage, type = "error", duration = 15)
    removeNotificationFn("protein_peptide_working")

    list(
      status = "error",
      errorMessage = errorMessage
    )
  })
}

runProteinPeptideRevertStep <- function(workflowData,
                                        logInfoFn = logger::log_info) {
  history <- workflowData$state_manager$getHistory()
  if (length(history) <= 1) {
    stop("No previous state to revert to.")
  }

  previousState <- history[length(history) - 1]
  revertedS4 <- workflowData$state_manager$revertToState(previousState)
  logInfoFn(paste("Reverted protein peptide filter to", previousState))

  list(
    previousState = previousState,
    revertedS4 = revertedS4,
    resultText = paste("Reverted to previous state:", previousState)
  )
}

runProteinPeptideRevertObserver <- function(workflowData,
                                            output,
                                            runRevertStepFn = runProteinPeptideRevertStep,
                                            renderTextFn = shiny::renderText,
                                            showNotificationFn = shiny::showNotification,
                                            logErrorFn = logger::log_error) {
  tryCatch({
    revertResult <- runRevertStepFn(workflowData = workflowData)
    output$protein_peptida_results <- renderTextFn(revertResult$resultText)
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

mod_prot_qc_peptide_protein_server <- function(id, workflow_data, omic_type, experiment_label) {
  shiny::moduleServer(id, function(input, output, session) {
    
    protein_peptide_plot <- shiny::reactiveVal(NULL)
    
    # Step 4: Protein Peptide Count Filter (chunk 13)
    shiny::observeEvent(input$apply_protein_peptide_filter, {
      runProteinPeptideApplyObserver(
        workflowData = workflow_data,
        minPeptidesPerProtein = input$min_peptides_per_protein,
        minPeptidoformsPerProtein = input$min_peptidoforms_per_protein,
        output = output,
        proteinPeptidePlot = protein_peptide_plot,
        omicType = omic_type,
        experimentLabel = experiment_label
      )
    })
    
    # Revert Protein Peptide Filter
    shiny::observeEvent(input$revert_protein_peptide, {
      runProteinPeptideRevertObserver(
        workflowData = workflow_data,
        output = output
      )
    })
    
    # Render protein peptide filter plot
    output$protein_peptide_plot <- shiny::renderPlot({
      shiny::req(protein_peptide_plot())
      grid::grid.draw(protein_peptide_plot())
    })
  })
}
