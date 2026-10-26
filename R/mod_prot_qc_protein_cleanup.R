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

#' @title Protein Accession Cleanup Module
#'
#' @description A Shiny module for performing protein accession cleanup.
#'
#' @name mod_prot_qc_protein_cleanup
#' @param id,proteinObject Runtime inputs used by this function; see the usage section for accepted values.
NULL

#' @rdname mod_prot_qc_protein_cleanup
#' @export
#' @importFrom shiny NS tagList tabPanel br fluidRow column wellPanel h4 p hr textInput helpText div icon strong textOutput verbatimTextOutput selectInput actionButton plotOutput
mod_prot_qc_protein_cleanup_ui <- function(id) {
  ns <- shiny::NS(id)
  
  shiny::tabPanel(
    "Accession Cleanup",
    shiny::br(),
    shiny::fluidRow(
      shiny::column(4,
        shiny::wellPanel(
          shiny::h4("Best Protein Accession Selection"),
          shiny::p("Choose the best protein accession for ambiguous mappings using FASTA metadata."),
          shiny::hr(),
          
          shiny::textInput(ns("delimiter"), 
            "Protein ID Delimiter", 
            value = ";",
            width = "100%"
          ),
          shiny::helpText("Character separating multiple protein IDs (default: ';')"),
          
          # Detected delimiters info box
          shiny::div(
            style = "background-color: #e7f3ff; padding: 10px; border-left: 4px solid #2196F3; border-radius: 4px; margin-top: 10px;",
            shiny::icon("info-circle", style = "color: #2196F3;"),
            shiny::strong(" Detected Delimiters in Your Data:"),
            shiny::br(),
            shiny::textOutput(ns("detected_delimiters"), inline = TRUE)
          ),
          
          # Sample protein IDs display
          shiny::div(
            style = "background-color: #f5f5f5; padding: 10px; border-left: 4px solid #9E9E9E; border-radius: 4px; margin-top: 10px; font-family: monospace; font-size: 11px;",
            shiny::icon("list", style = "color: #9E9E9E;"),
            shiny::strong(" Sample Protein IDs:"),
            shiny::br(),
            shiny::verbatimTextOutput(ns("sample_protein_ids"))
          ),
          
          shiny::selectInput(ns("aggregation_method"), 
            "Aggregation Method", 
            choices = c("sum", "mean", "median"),
            selected = "mean",
            width = "100%"
          ),
          shiny::helpText("Method for combining duplicate protein measurements (S4 default: mean)"),
          
          shiny::hr(),
          shiny::div(
            shiny::actionButton(ns("apply_accession_cleanup"), "Apply Cleanup", 
              class = "btn-primary", width = "48%"),
            shiny::actionButton(ns("revert_accession_cleanup"), "Revert", 
              class = "btn-warning", width = "48%", style = "margin-left: 4%")
          )
        )
      ),
      shiny::column(8,
        shiny::verbatimTextOutput(ns("accession_cleanup_results")),
        shiny::br(),
        shinyjqui::jqui_resizable(
          shiny::plotOutput(ns("accession_cleanup_plot"), height = "800px", width = "100%")
        )
      )
    )
  )
}

#' @rdname mod_prot_qc_protein_cleanup
#' @export
#' @importFrom shiny moduleServer reactiveVal observeEvent req showNotification removeNotification renderText renderPlot textOutput
#' @importFrom logger log_info log_error
#' @importFrom grid grid.draw
.proteinCleanupActiveKey <- function(proteinObject) {
  declaredKey <- tryCatch(
    proteinObject@protein_id_column,
    error = function(...) character()
  )
  if (length(declaredKey) == 1L &&
      declaredKey %in% names(proteinObject@protein_quant_table)) {
    return(declaredKey)
  }

  candidates <- c("Protein.Group", "Protein.Ids", "Protein.IDs", "protein_id")
  inferredKey <- candidates[candidates %in% names(proteinObject@protein_quant_table)][1]
  if (length(inferredKey) == 0L || is.na(inferredKey)) {
    stop("No active protein identity column is available.", call. = FALSE)
  }
  inferredKey
}

runProteinAccessionCleanupStep <- function(workflowData,
                                           delimiter,
                                           aggregationMethod,
                                           chooseBestProteinAccessionFn = chooseBestProteinAccession,
                                           nowFn = Sys.time,
                                           logInfoFn = logger::log_info,
                                           logWarnFn = logger::log_warn,
                                           saveRdsFn = saveRDS,
                                           existsFn = exists,
                                           getFn = get) {
  pending <- if (identical(workflowData$state_manager$current_state, "protein_accession_cleaned"))
    resumeDIAProteinAction(workflowData, "protein_accession_cleaned") else NULL
  if (!is.null(pending)) return(pending)
  # KIT-015: an unmanaged launch runs the registered declaration in process with its production variant.
  entry <- workflowActionSpec("protein_accession_cleanup")
  parameters <- list(delimiter = delimiter, aggregation_method = aggregationMethod)
  validateWorkflowActionParameters(entry$name, parameters, workflowExecutionPolicy(FALSE), resolved = FALSE)
  value <- computeDIAReferenceAccessionCleanup(workflowData$state_manager, getDIAStateReference(workflowData$state_manager),
    parameters, variant = entry$default_variant, seqinr = shiny::isolate(workflowData$aa_seq_tbl_final),
    fasta_metadata = shiny::isolate(workflowData$fasta_metadata), now_fn = nowFn)
  return(finishDIAReferenceProteinStep(workflowData, value))
}

updateProteinAccessionCleanupOutputs <- function(output,
                                                 accessionCleanupPlot,
                                                 cleanupResult,
                                                 omicType,
                                                 experimentLabel,
                                                 renderTextFn = shiny::renderText,
                                                 updateProteinFilteringFn = updateProteinFiltering) {
  output$accession_cleanup_results <- renderPeptideQcResultText(cleanupResult$resultText, renderTextFn)

  plotGrid <- publishDIAProteinProgress(cleanupResult, "10_protein_accession_cleaned", omicType,
    experimentLabel, updateProteinFilteringFn, cleanupResult$cleanedS4)
  accessionCleanupPlot(plotGrid)

  invisible(plotGrid)
}

runProteinAccessionCleanupApplyObserver <- function(workflowData,
                                                    delimiter,
                                                    aggregationMethod,
                                                    output,
                                                    accessionCleanupPlot,
                                                    omicType,
                                                    experimentLabel,
                                                    runApplyStepFn = runProteinAccessionCleanupStep,
                                                    updateOutputsFn = updateProteinAccessionCleanupOutputs,
                                                    showNotificationFn = shiny::showNotification,
                                                    removeNotificationFn = shiny::removeNotification,
                                                    logInfoFn = logger::log_info,
                                                    logErrorFn = logger::log_error) {
  withDIAStateScope(getDIAActiveManager(), function() {
  showNotificationFn(
    "Applying protein accession cleanup...",
    id = "accession_cleanup_working",
    duration = NULL
  )

  tryCatch({
    cleanupResult <- runApplyStepFn(
      workflowData = workflowData,
      delimiter = delimiter,
      aggregationMethod = aggregationMethod
    )

    plotGrid <- updateOutputsFn(
      output = output,
      accessionCleanupPlot = accessionCleanupPlot,
      cleanupResult = cleanupResult,
      omicType = omicType,
      experimentLabel = experimentLabel
    )

    logInfoFn("Protein accession cleanup completed")
    removeNotificationFn("accession_cleanup_working")
    showNotificationFn("Protein accession cleanup completed", type = "message")

    list(
      status = "success",
      cleanupResult = cleanupResult,
      plotGrid = plotGrid
    )
  }, error = function(e) {
    errorMessage <- paste("Error in protein accession cleanup:", e$message)
    logErrorFn(errorMessage)
    showNotificationFn(errorMessage, type = "error", duration = 15)
    removeNotificationFn("accession_cleanup_working")

    list(
      status = "error",
      errorMessage = errorMessage
    )
  })
  })
}

runProteinAccessionCleanupRevertStep <- function(workflowData) {
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

runProteinAccessionCleanupRevertObserver <- function(workflowData,
                                                     output,
                                                     runRevertStepFn = runProteinAccessionCleanupRevertStep,
                                                     renderTextFn = shiny::renderText,
                                                     showNotificationFn = shiny::showNotification,
                                                     logInfoFn = logger::log_info,
                                                     logErrorFn = logger::log_error) {
  withDIAStateScope(getDIAActiveManager(), function() {
  tryCatch({
    revertResult <- runRevertStepFn(workflowData = workflowData)
    output$accession_cleanup_results <- renderPeptideQcResultText(revertResult$resultText, renderTextFn)
    logInfoFn(paste("Reverted accession cleanup to", revertResult$previousState))
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

bindProteinAccessionCleanupPlot <- function(output, accessionCleanupPlot, plotIdentity = NULL) {
  output$accession_cleanup_plot <- if (is.null(plotIdentity)) renderPlot({
    req(accessionCleanupPlot())
    grid.draw(accessionCleanupPlot())
  }) else renderPlot({
    req(accessionCleanupPlot())
    grid.draw(accessionCleanupPlot())
  }, alt = function() plotIdentity())
}

mod_prot_qc_protein_cleanup_server <- function(id, workflow_data, experiment_paths = NULL, omic_type, experiment_label) {
  shiny::moduleServer(id, function(input, output, session) {
    
    accession_cleanup_plot <- shiny::reactiveVal(NULL)
    accession_cleanup_plot_identity <- shiny::reactiveVal("Protein accession cleanup summary")
    accession_cleanup_plot_value <- accession_cleanup_plot

    # KIT-015: the managed worker when the session carries a managed application (KIT-008 pattern), else the
    # same registered declaration in process.
    managed_step <- installWorkflowManagedStep(session, input, workflow_data, "protein_accession_cleanup",
      publication = list(experiment_label = experiment_label, experiment_paths = experiment_paths, omic_type = omic_type),
      ui = list(text = "accession_cleanup_results", plot = "accession_cleanup_plot", ack = "workflow_ack", notification = "accession_cleanup_working",
        label = "protein accession cleanup", step_label = "10_protein_accession_cleaned", plot_identity = accession_cleanup_plot_identity),
      deliver = function(result) updateProteinAccessionCleanupOutputs(output = output, accessionCleanupPlot = accession_cleanup_plot_value, cleanupResult = result,
        omicType = omic_type, experimentLabel = experiment_label),
      read_result = readManagedProteinResult)
    
    # Step 2: Protein Accession Cleanup (chunk 19)
    shiny::observeEvent(input$apply_accession_cleanup, {
      if (!is.null(workflowSessionApplication(session))) {
        managed_step$apply(list(delimiter = input$delimiter, aggregation_method = input$aggregation_method))
        return(invisible(NULL))
      }
      runProteinAccessionCleanupApplyObserver(
        workflowData = workflow_data,
        delimiter = input$delimiter,
        aggregationMethod = input$aggregation_method,
        output = output,
        accessionCleanupPlot = accession_cleanup_plot,
        omicType = omic_type,
        experimentLabel = experiment_label
      )
    })
    
    # Revert Accession Cleanup
    shiny::observeEvent(input$revert_accession_cleanup, {
      managed_step$cancel("reverted")
      runProteinAccessionCleanupRevertObserver(
        workflowData = workflow_data,
        output = output
      )
    })

    bindProteinAccessionCleanupPlot(
      output = output,
      accessionCleanupPlot = accession_cleanup_plot,
      plotIdentity = accession_cleanup_plot_identity
    )
  })
}
