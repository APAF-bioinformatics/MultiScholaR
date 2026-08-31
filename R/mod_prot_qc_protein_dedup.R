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

#' @title Protein Duplicate Removal Module
#'
#' @description A Shiny module for removing duplicate proteins.
#'
#' @name mod_prot_qc_protein_dedup
#' @param id,workflowData,aggregationMethod,aggregationResolverFn,logInfoFn Runtime inputs used by this function; see the usage section for accepted values.
NULL

#' @rdname mod_prot_qc_protein_dedup
#' @export
#' @importFrom shiny NS tagList tabPanel br fluidRow column wellPanel h4 p hr selectInput helpText div actionButton verbatimTextOutput plotOutput
mod_prot_qc_protein_dedup_ui <- function(id) {
  ns <- shiny::NS(id)
  
  shiny::tabPanel(
    "Duplicate Removal",
    shiny::br(),
    shiny::fluidRow(
      shiny::column(4,
        shiny::wellPanel(
          shiny::h4("Remove Duplicate Proteins"),
          shiny::p("Aggregate duplicate protein entries by taking the mean across matching proteins."),
          shiny::hr(),
          
          shiny::selectInput(ns("duplicate_aggregation_method"), 
            "Aggregation Method", 
            choices = c("mean", "median", "max"),
            selected = "mean",
            width = "100%"
          ),
          shiny::helpText("Method for combining duplicate protein measurements (default: mean)"),
          
          shiny::hr(),
          shiny::div(
            shiny::actionButton(ns("apply_duplicate_removal"), "Remove Duplicates", 
              class = "btn-primary", width = "48%"),
            shiny::actionButton(ns("revert_duplicate_removal"), "Revert", 
              class = "btn-warning", width = "48%", style = "margin-left: 4%")
          )
        )
      ),
      shiny::column(8,
        shiny::verbatimTextOutput(ns("duplicate_removal_results")),
        shiny::br(),
        shinyjqui::jqui_resizable(
          shiny::plotOutput(ns("duplicate_removal_plot"), height = "800px", width = "100%")
        )
      )
    )
  )
}

#' @rdname mod_prot_qc_protein_dedup
#' @export
#' @importFrom shiny moduleServer reactiveVal observeEvent req showNotification removeNotification renderText renderPlot
#' @importFrom logger log_info log_error
#' @importFrom grid grid.draw
runProteinDuplicateRemovalStep <- function(workflowData,
                                          aggregationMethod,
                                          aggregationResolverFn = base::get,
                                          logInfoFn = logger::log_info) {
  shiny::req(workflowData$state_manager)
  if (protDiaQcWorkerEligible(workflowData)) {
    worker <- runProtDiaQcProcess(protDiaQcWorkerSpec(
      workflowData,
      "protein_duplicate_aggregation",
      list(aggregationMethod = aggregationMethod)
    ))
    result <- applyProtDiaQcWorkerResult(workflowData, worker)
    result$plot_png <- worker$plot_png
    return(result)
  }

  selectedS4 <- workflowData$state_manager$getState()
  shiny::req(selectedS4)
  currentS4 <- selectedS4

  logInfoFn(sprintf(
    "Protein Processing: Removing duplicate proteins using %s",
    aggregationMethod
  ))

  proteinIdColumn <- resolveProteinQuantIdentityColumn(currentS4)
  proteinTable <- currentS4@protein_quant_table

  duplicates <- proteinTable |>
    dplyr::group_by(!!rlang::sym(proteinIdColumn)) |>
    dplyr::filter(dplyr::n() > 1) |>
    dplyr::distinct(!!rlang::sym(proteinIdColumn)) |>
    dplyr::pull(!!rlang::sym(proteinIdColumn))

  aggregationFn <- aggregationResolverFn(aggregationMethod)
  quantColumns <- setdiff(colnames(proteinTable), proteinIdColumn)
  numericQuantColumns <- quantColumns[vapply(
    proteinTable[quantColumns],
    is.numeric,
    logical(1)
  )]
  if (length(numericQuantColumns) == 0L) {
    stop("No numeric protein intensity columns available for duplicate removal.")
  }

  currentS4@protein_quant_table <- currentS4@protein_quant_table |>
    dplyr::group_by(!!rlang::sym(proteinIdColumn)) |>
    dplyr::summarise(
      dplyr::across(dplyr::all_of(numericQuantColumns), ~ aggregationFn(.x, na.rm = TRUE)),
      .groups = "drop"
    ) |>
    dplyr::ungroup()

  stateConfig <- list(
    aggregation_method = aggregationMethod,
    active_protein_key = proteinIdColumn,
    duplicates_found = duplicates,
    num_duplicates = length(duplicates)
  )
  currentS4 <- saveProtProteinQcState(
    workflow_data = workflowData,
    state_manager = workflowData$state_manager,
    before = selectedS4,
    after = currentS4,
    stage_id = "protein_duplicate_aggregation",
    state_name = "duplicates_removed",
    config_object = stateConfig,
    audit_parameters = stateConfig,
    description = "Removed duplicate proteins by aggregation",
    transformation_type = "aggregation"
  )

  proteinCount <- countDistinctProteinQuantIdentities(currentS4)

  resultText <- paste(
    "Duplicate Protein Removal Completed Successfully\n",
    "===============================================\n",
    sprintf("Proteins remaining: %d\n", proteinCount),
    sprintf("Duplicates found: %d\n", length(duplicates)),
    sprintf("Aggregation method: %s\n", aggregationMethod),
    "State saved as: 'duplicates_removed'\n"
  )

  list(
    deduplicatedS4 = currentS4,
    duplicates = duplicates,
    resultText = resultText
  )
}

updateProteinDuplicateRemovalOutputs <- function(output,
                                                 duplicateRemovalPlot,
                                                 duplicateRemovalResult,
                                                 omicType,
                                                 experimentLabel,
                                                 renderTextFn = shiny::renderText,
                                                 updateProteinFilteringFn = updateProteinFiltering,
                                                 workflowData = NULL) {
  output$duplicate_removal_results <- renderTextFn(duplicateRemovalResult$resultText)

  plotGrid <- if (is.raw(duplicateRemovalResult$plot_png)) {
    protDiaQcPlotGrid(duplicateRemovalResult$plot_png)
  } else {
    protDiaProteinQcUpdateFiltering(
      update_fn = updateProteinFilteringFn,
      workflow_data = workflowData,
      data = duplicateRemovalResult$deduplicatedS4@protein_quant_table,
      step_name = "12_duplicates_removed",
      omic_type = omicType,
      experiment_label = experimentLabel,
      return_grid = TRUE,
      overwrite = TRUE
    )
  }
  duplicateRemovalPlot(plotGrid)

  invisible(plotGrid)
}

runProteinDuplicateRemovalApplyObserver <- function(workflowData,
                                                    aggregationMethod,
                                                    output,
                                                    duplicateRemovalPlot,
                                                    omicType,
                                                    experimentLabel,
                                                    runApplyStepFn = runProteinDuplicateRemovalStep,
                                                    updateOutputsFn = updateProteinDuplicateRemovalOutputs,
                                                    showNotificationFn = shiny::showNotification,
                                                    removeNotificationFn = shiny::removeNotification,
                                                    logInfoFn = logger::log_info,
                                                    logErrorFn = logger::log_error) {
  showNotificationFn(
    "Removing duplicate proteins...",
    id = "duplicate_removal_working",
    duration = NULL
  )

  tryCatch({
    duplicateRemovalResult <- runApplyStepFn(
      workflowData = workflowData,
      aggregationMethod = aggregationMethod
    )

    plotGrid <- protDiaProteinQcInvokeOutputs(
      update_outputs_fn = updateOutputsFn,
      args = list(
        output = output,
        duplicateRemovalPlot = duplicateRemovalPlot,
        duplicateRemovalResult = duplicateRemovalResult,
        omicType = omicType,
        experimentLabel = experimentLabel
      ),
      workflow_data = workflowData
    )

    logInfoFn("Duplicate protein removal completed successfully")
    removeNotificationFn("duplicate_removal_working")
    showNotificationFn(
      "Duplicate protein removal completed successfully",
      type = "message"
    )

    list(
      status = "success",
      duplicateRemovalResult = duplicateRemovalResult,
      plotGrid = plotGrid
    )
  }, error = function(e) {
    errorMessage <- paste("Error removing duplicate proteins:", e$message)
    logErrorFn(errorMessage)
    showNotificationFn(errorMessage, type = "error", duration = 15)
    removeNotificationFn("duplicate_removal_working")

    list(
      status = "error",
      errorMessage = errorMessage
    )
  })
}

runProteinDuplicateRemovalRevertStep <- function(workflowData) {
  shiny::req(workflowData$state_manager)
  history <- workflowData$state_manager$getHistory()

  if (length(history) <= 1) {
    stop("No previous state to revert to.")
  }

  previousState <- history[length(history) - 1]
  if (protDiaQcWorkerEligible(workflowData)) {
    return(runProtDiaQcRevert(workflowData, previousState))
  }
  revertedS4 <- revertProtDiaProteinQcState(workflowData, previousState)

  list(
    previousState = previousState,
    revertedS4 = revertedS4,
    resultText = paste("Reverted to previous state:", previousState)
  )
}

runProteinDuplicateRemovalRevertObserver <- function(workflowData,
                                                     output,
                                                     runRevertStepFn = runProteinDuplicateRemovalRevertStep,
                                                     renderTextFn = shiny::renderText,
                                                     showNotificationFn = shiny::showNotification,
                                                     logInfoFn = logger::log_info,
                                                     logErrorFn = logger::log_error) {
  tryCatch({
    revertResult <- runRevertStepFn(workflowData = workflowData)
    output$duplicate_removal_results <- renderTextFn(revertResult$resultText)
    logInfoFn(paste("Reverted duplicate removal to", revertResult$previousState))
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

bindProteinDuplicateRemovalPlot <- function(output, duplicateRemovalPlot) {
  output$duplicate_removal_plot <- renderPlot({
    req(duplicateRemovalPlot())
    grid.draw(duplicateRemovalPlot())
  })
}

mod_prot_qc_protein_dedup_server <- function(id, workflow_data, omic_type, experiment_label) {
  shiny::moduleServer(id, function(input, output, session) {
    
    duplicate_removal_plot <- shiny::reactiveVal(NULL)
    
    # Step 4: Duplicate Protein Removal (chunk 22)
    shiny::observeEvent(input$apply_duplicate_removal, {
      runProteinDuplicateRemovalApplyObserver(
        workflowData = workflow_data,
        aggregationMethod = input$duplicate_aggregation_method,
        output = output,
        duplicateRemovalPlot = duplicate_removal_plot,
        omicType = omic_type,
        experimentLabel = experiment_label
      )
    })
    
    # Revert Duplicate Removal
    shiny::observeEvent(input$revert_duplicate_removal, {
      runProteinDuplicateRemovalRevertObserver(
        workflowData = workflow_data,
        output = output
      )
    })
    
    bindProteinDuplicateRemovalPlot(
      output = output,
      duplicateRemovalPlot = duplicate_removal_plot
    )
  })
}
