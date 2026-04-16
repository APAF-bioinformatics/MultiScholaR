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
  shiny::req(workflowData$state_manager)
  currentS4 <- workflowData$state_manager$getState()
  shiny::req(currentS4)

  logInfoFn(sprintf(
    "Protein Processing: Applying accession cleanup with delimiter: %s",
    delimiter
  ))

  proteinsBefore <- currentS4@protein_quant_table |>
    dplyr::distinct(Protein.Ids) |>
    nrow()

  if (existsFn("aa_seq_tbl_final", envir = .GlobalEnv, inherits = FALSE)) {
    aaSeqTblFinal <- getFn("aa_seq_tbl_final", envir = .GlobalEnv) |>
      dplyr::rename(uniprot_acc = database_id)

    cleanedS4 <- chooseBestProteinAccessionFn(
      theObject = currentS4,
      delim = delimiter,
      seqinr_obj = aaSeqTblFinal,
      seqinr_accession_column = "uniprot_acc",
      replace_zero_with_na = TRUE,
      aggregation_method = aggregationMethod
    )
    cleanupApplied <- TRUE
  } else {
    cleanedS4 <- currentS4
    cleanupApplied <- FALSE
  }

  proteinsAfter <- cleanedS4@protein_quant_table |>
    dplyr::distinct(Protein.Ids) |>
    nrow()

  workflowData$accession_cleanup_results <- list(
    cleanup_applied = cleanupApplied,
    delimiter_used = delimiter,
    aggregation_method = aggregationMethod,
    proteins_before = proteinsBefore,
    proteins_after = proteinsAfter,
    had_full_metadata = if (!is.null(workflowData$fasta_metadata)) {
      workflowData$fasta_metadata$has_protein_evidence &&
        workflowData$fasta_metadata$has_gene_names
    } else {
      FALSE
    },
    timestamp = nowFn()
  )

  if (is.null(workflowData$qc_params)) {
    workflowData$qc_params <- list()
  }
  if (is.null(workflowData$qc_params$protein_qc)) {
    workflowData$qc_params$protein_qc <- list()
  }
  workflowData$qc_params$protein_qc$accession_cleanup <-
    workflowData$accession_cleanup_results

  tryCatch({
    if (existsFn("experiment_paths", envir = .GlobalEnv, inherits = FALSE)) {
      experimentPaths <- getFn("experiment_paths", envir = .GlobalEnv)
      if (!is.null(experimentPaths$source_dir)) {
        accessionCleanupFile <- file.path(
          experimentPaths$source_dir,
          "accession_cleanup_results.RDS"
        )
        saveRdsFn(workflowData$accession_cleanup_results, accessionCleanupFile)
        logInfoFn(sprintf(
          "Saved accession cleanup results to: %s",
          accessionCleanupFile
        ))
      }
    }
  }, error = function(e) {
    logWarnFn(sprintf(
      "Could not save accession cleanup results file: %s",
      e$message
    ))
  })

  logInfoFn(sprintf(
    "Accession cleanup results tracked: %d -> %d proteins",
    proteinsBefore,
    proteinsAfter
  ))

  workflowData$state_manager$saveState(
    state_name = "protein_accession_cleaned",
    s4_data_object = cleanedS4,
    config_object = list(
      delimiter = delimiter,
      aggregation_method = aggregationMethod,
      cleanup_applied = cleanupApplied
    ),
    description = if (cleanupApplied) {
      "Applied protein accession cleanup"
    } else {
      "Skipped accession cleanup (no FASTA data)"
    }
  )

  proteinCount <- cleanedS4@protein_quant_table |>
    dplyr::distinct(Protein.Ids) |>
    nrow()

  resultText <- if (cleanupApplied) {
    paste(
      "Protein Accession Cleanup Applied Successfully\n",
      "==============================================\n",
      sprintf("Proteins remaining: %d\n", proteinCount),
      sprintf("Delimiter: %s\n", delimiter),
      sprintf("Aggregation method: %s\n", aggregationMethod),
      "State saved as: 'protein_accession_cleaned'\n"
    )
  } else {
    paste(
      "Protein Accession Cleanup Skipped\n",
      "=================================\n",
      sprintf("Proteins remaining: %d\n", proteinCount),
      "Reason: No FASTA data available for cleanup\n",
      "State saved as: 'protein_accession_cleaned'\n"
    )
  }

  list(
    cleanedS4 = cleanedS4,
    cleanupApplied = cleanupApplied,
    resultText = resultText
  )
}

updateProteinAccessionCleanupOutputs <- function(output,
                                                 accessionCleanupPlot,
                                                 cleanupResult,
                                                 omicType,
                                                 experimentLabel,
                                                 renderTextFn = shiny::renderText,
                                                 updateProteinFilteringFn = updateProteinFiltering) {
  output$accession_cleanup_results <- renderTextFn(cleanupResult$resultText)

  plotGrid <- updateProteinFilteringFn(
    data = cleanupResult$cleanedS4@protein_quant_table,
    step_name = "10_protein_accession_cleaned",
    omic_type = omicType,
    experiment_label = experimentLabel,
    return_grid = TRUE,
    overwrite = TRUE
  )
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
}

runProteinAccessionCleanupRevertStep <- function(workflowData) {
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

runProteinAccessionCleanupRevertObserver <- function(workflowData,
                                                     output,
                                                     runRevertStepFn = runProteinAccessionCleanupRevertStep,
                                                     renderTextFn = shiny::renderText,
                                                     showNotificationFn = shiny::showNotification,
                                                     logInfoFn = logger::log_info,
                                                     logErrorFn = logger::log_error) {
  tryCatch({
    revertResult <- runRevertStepFn(workflowData = workflowData)
    output$accession_cleanup_results <- renderTextFn(revertResult$resultText)
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
}

bindProteinAccessionCleanupPlot <- function(output, accessionCleanupPlot) {
  output$accession_cleanup_plot <- renderPlot({
    req(accessionCleanupPlot())
    grid.draw(accessionCleanupPlot())
  })
}

mod_prot_qc_protein_cleanup_server <- function(id, workflow_data, omic_type, experiment_label) {
  shiny::moduleServer(id, function(input, output, session) {
    
    accession_cleanup_plot <- shiny::reactiveVal(NULL)
    
    # Step 2: Protein Accession Cleanup (chunk 19)
    shiny::observeEvent(input$apply_accession_cleanup, {
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
      runProteinAccessionCleanupRevertObserver(
        workflowData = workflow_data,
        output = output
      )
    })

    bindProteinAccessionCleanupPlot(
      output = output,
      accessionCleanupPlot = accession_cleanup_plot
    )
  })
}
