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

#' @title IQ Protein Rollup Module
#'
#' @description A Shiny module for performing peptide-to-protein rollup using the IQ tool.
#'
#' @name mod_prot_qc_protein_rollup
NULL

#' @rdname mod_prot_qc_protein_rollup
#' @export
#' @importFrom shiny NS tagList tabPanel br fluidRow column wellPanel h4 p hr div actionButton verbatimTextOutput plotOutput
mod_prot_qc_protein_rollup_ui <- function(id) {
  ns <- shiny::NS(id)
  
  shiny::tabPanel(
    "IQ Protein Rollup",
    shiny::br(),
    shiny::fluidRow(
      shiny::column(4,
        shiny::wellPanel(
          shiny::h4("Peptide-to-Protein Rollup"),
          shiny::p("Aggregate peptide-level data to protein-level quantification using the IQ algorithm (MaxLFQ)."),
          shiny::hr(),
          shiny::p("This step uses the IQ tool to implement the MaxLFQ algorithm for protein quantification, then automatically creates a ProteinQuantitativeData S4 object."),
          shiny::hr(),
          shiny::div(
            shiny::actionButton(ns("apply_iq_rollup"), "Run IQ Rollup & Create S4 Object", 
              class = "btn-primary", width = "48%"),
            shiny::actionButton(ns("revert_iq_rollup"), "Revert", 
              class = "btn-warning", width = "48%", style = "margin-left: 4%")
          )
        )
      ),
      shiny::column(8,
        shiny::verbatimTextOutput(ns("iq_rollup_results")),
        shiny::br(),
        shinyjqui::jqui_resizable(
          shiny::plotOutput(ns("iq_rollup_plot"), height = "800px", width = "100%")
        )
      )
    )
  )
}

#' @rdname mod_prot_qc_protein_rollup
#' @export
#' @importFrom shiny moduleServer reactiveVal observeEvent req showNotification removeNotification renderText renderPlot
#' @importFrom logger log_info log_error log_warn
#' @importFrom grid grid.draw
#' @importFrom tibble tibble
#' @importFrom purrr map_chr
runProteinIqRollupApplyStep <- function(workflowData,
                                        experimentPaths,
                                        createProteinDataFn = ProteinQuantitativeData,
                                        writeTsvFn = vroom::vroom_write,
                                        processLongFormatFn = iq::process_long_format,
                                        readTsvFn = vroom::vroom,
                                        captureCheckpointFn = .capture_checkpoint,
                                        showNotificationFn = shiny::showNotification,
                                        logInfoFn = logger::log_info,
                                        logWarnFn = logger::log_warn,
                                        sleepFn = Sys.sleep,
                                        maxWait = 30) {
  shiny::req(workflowData$state_manager)

  currentState <- workflowData$state_manager$current_state
  peptideS4 <- workflowData$state_manager$getState(currentState)
  shiny::req(peptideS4)

  logInfoFn("Protein Processing: Starting IQ rollup from peptide state")

  peptideValuesImputedFile <- file.path(
    experimentPaths$peptide_qc_dir,
    "peptide_values_imputed.tsv"
  )

  originalSamples <- unique(peptideS4@design_matrix$Run)
  sampleMapping <- tibble::tibble(
    Original = originalSamples,
    Alias = paste0("S_", sprintf("%03d", seq_along(originalSamples)))
  )

  peptideDataForIq <- peptideS4@peptide_data |>
    dplyr::mutate(
      Q.Value = 0.0009,
      PG.Q.Value = 0.009,
      Peptide.Imputed = ifelse(is.na(Peptide.Imputed), 0, Peptide.Imputed)
    ) |>
    dplyr::left_join(sampleMapping, by = c("Run" = "Original")) |>
    dplyr::mutate(Run = Alias) |>
    dplyr::select(-Alias)

  writeTsvFn(peptideDataForIq, peptideValuesImputedFile)

  iqOutputFile <- file.path(experimentPaths$protein_qc_dir, "iq_output_file.txt")

  processLongFormatFn(
    peptideValuesImputedFile,
    output_filename = iqOutputFile,
    sample_id = "Run",
    primary_id = "Protein.Ids",
    secondary_id = "Stripped.Sequence",
    intensity_col = "Peptide.Imputed",
    filter_double_less = c("Q.Value" = "0.01", "PG.Q.Value" = "0.01"),
    normalization = "none"
  )

  waitCount <- 0
  while (!file.exists(iqOutputFile) && waitCount < maxWait) {
    sleepFn(1)
    waitCount <- waitCount + 1
  }

  if (!file.exists(iqOutputFile)) {
    stop("IQ output file not created within timeout period")
  }

  proteinLog2QuantAliased <- readTsvFn(iqOutputFile, .name_repair = "minimal")
  currentCols <- colnames(proteinLog2QuantAliased)
  restoredCols <- purrr::map_chr(currentCols, function(col) {
    if (col == "Protein.Ids") {
      return(col)
    }

    match <- sampleMapping$Original[sampleMapping$Alias == col]
    if (length(match) > 0) {
      return(match[[1]])
    }

    col
  })
  colnames(proteinLog2QuantAliased) <- restoredCols
  proteinLog2Quant <- proteinLog2QuantAliased

  survivingSamples <- setdiff(colnames(proteinLog2Quant), "Protein.Ids")
  finalDesignMatrix <- peptideS4@design_matrix |>
    dplyr::filter(Run %in% survivingSamples)

  droppedSamples <- setdiff(originalSamples, survivingSamples)
  if (length(droppedSamples) > 0) {
    droppedList <- paste(droppedSamples, collapse = ", ")
    logWarnFn(
      paste0(
        "Protein Processing: ",
        length(droppedSamples),
        " samples dropped during IQ rollup: ",
        droppedList
      )
    )
    showNotificationFn(
      paste(
        "Warning: The following samples were dropped during protein rollup due to low peptide quality:",
        droppedList
      ),
      type = "warning",
      duration = 10
    )
  }

  logInfoFn("Protein Processing: Creating ProteinQuantitativeData S4 object")

  proteinObj <- createProteinDataFn(
    protein_quant_table = proteinLog2Quant,
    protein_id_column = "Protein.Ids",
    protein_id_table = proteinLog2Quant |> dplyr::distinct(Protein.Ids),
    design_matrix = finalDesignMatrix,
    sample_id = "Run",
    group_id = "group",
    technical_replicate_id = "replicates",
    args = peptideS4@args
  )

  workflowData$state_manager$saveState(
    state_name = "protein_s4_created",
    s4_data_object = proteinObj,
    config_object = list(
      iq_output_file = iqOutputFile,
      peptide_input_file = peptideValuesImputedFile,
      s4_class = "ProteinQuantitativeData",
      protein_id_column = "Protein.Ids"
    ),
    description = "IQ protein rollup completed and ProteinQuantitativeData S4 object created"
  )

  captureCheckpointFn(proteinObj, "cp03", "rolled_up_protein")

  proteinCount <- proteinObj@protein_quant_table |>
    dplyr::distinct(Protein.Ids) |>
    nrow()

  resultText <- paste(
    "IQ Protein Rollup & S4 Object Creation Completed Successfully\n",
    "============================================================\n",
    sprintf("Proteins quantified: %d\n", proteinCount),
    sprintf("Samples: %d\n", ncol(proteinObj@protein_quant_table) - 1),
    "Algorithm: MaxLFQ (via IQ tool)\n",
    sprintf("S4 Class: %s\n", class(proteinObj)[1]),
    sprintf("Design matrix: %s\n", paste(colnames(proteinObj@design_matrix), collapse = ", ")),
    sprintf("Output file: %s\n", basename(iqOutputFile)),
    "State saved as: 'protein_s4_created'\n",
    "\nReady for protein accession cleanup."
  )

  list(
    proteinObj = proteinObj,
    resultText = resultText
  )
}

updateProteinIqRollupOutputs <- function(output,
                                         iqRollupPlot,
                                         iqRollupResult,
                                         omicType,
                                         experimentLabel,
                                         renderTextFn = shiny::renderText,
                                         updateProteinFilteringFn = updateProteinFiltering) {
  output$iq_rollup_results <- renderTextFn(iqRollupResult$resultText)

  plotGrid <- updateProteinFilteringFn(
    data = iqRollupResult$proteinObj@protein_quant_table,
    step_name = "9_protein_s4_created",
    omic_type = omicType,
    experiment_label = experimentLabel,
    return_grid = TRUE,
    overwrite = TRUE
  )
  iqRollupPlot(plotGrid)

  invisible(plotGrid)
}

runProteinIqRollupApplyObserver <- function(workflowData,
                                            experimentPaths,
                                            output,
                                            iqRollupPlot,
                                            omicType,
                                            experimentLabel,
                                            runApplyStepFn = runProteinIqRollupApplyStep,
                                            updateOutputsFn = updateProteinIqRollupOutputs,
                                            showNotificationFn = shiny::showNotification,
                                            removeNotificationFn = shiny::removeNotification,
                                            logInfoFn = logger::log_info,
                                            logErrorFn = logger::log_error) {
  showNotificationFn(
    "Running IQ protein rollup & creating S4 object...",
    id = "iq_rollup_working",
    duration = NULL
  )

  tryCatch({
    iqRollupResult <- runApplyStepFn(
      workflowData = workflowData,
      experimentPaths = experimentPaths
    )

    plotGrid <- updateOutputsFn(
      output = output,
      iqRollupPlot = iqRollupPlot,
      iqRollupResult = iqRollupResult,
      omicType = omicType,
      experimentLabel = experimentLabel
    )

    logInfoFn("IQ protein rollup and S4 object creation completed successfully")
    removeNotificationFn("iq_rollup_working")
    showNotificationFn(
      "IQ protein rollup & S4 object creation completed successfully",
      type = "message"
    )

    list(
      status = "success",
      iqRollupResult = iqRollupResult,
      plotGrid = plotGrid
    )
  }, error = function(e) {
    errorMessage <- paste("Error in IQ protein rollup & S4 creation:", e$message)
    logErrorFn(errorMessage)
    showNotificationFn(errorMessage, type = "error", duration = 15)
    removeNotificationFn("iq_rollup_working")

    list(
      status = "error",
      errorMessage = errorMessage
    )
  })
}

runProteinIqRollupRevertStep <- function(workflowData) {
  shiny::req(workflowData$state_manager)
  history <- workflowData$state_manager$getHistory()
  peptideStates <- c(
    "imputed",
    "replicate_filtered",
    "sample_filtered",
    "protein_peptide_filtered",
    "intensity_filtered",
    "precursor_rollup",
    "qvalue_filtered",
    "raw_data_s4"
  )
  previousState <- intersect(rev(history), peptideStates)[1]

  if (length(previousState) == 0 || is.na(previousState)) {
    stop("No previous peptide state to revert to.")
  }

  revertedS4 <- workflowData$state_manager$revertToState(previousState)

  list(
    previousState = previousState,
    revertedS4 = revertedS4,
    resultText = paste("Reverted to", previousState, "state")
  )
}

runProteinIqRollupRevertObserver <- function(workflowData,
                                             output,
                                             runRevertStepFn = runProteinIqRollupRevertStep,
                                             renderTextFn = shiny::renderText,
                                             showNotificationFn = shiny::showNotification,
                                             logInfoFn = logger::log_info,
                                             logErrorFn = logger::log_error) {
  tryCatch({
    revertResult <- runRevertStepFn(workflowData = workflowData)
    output$iq_rollup_results <- renderTextFn(revertResult$resultText)
    logInfoFn(paste("Reverted IQ rollup to", revertResult$previousState))
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

bindProteinIqRollupPlot <- function(output, iqRollupPlot) {
  output$iq_rollup_plot <- renderPlot({
    req(iqRollupPlot())
    grid.draw(iqRollupPlot())
  })
}

mod_prot_qc_protein_rollup_server <- function(id, workflow_data, experiment_paths, omic_type, experiment_label) {
  shiny::moduleServer(id, function(input, output, session) {
    
    iq_rollup_plot <- shiny::reactiveVal(NULL)
    
    # Step 1: IQ Protein Rollup (chunk 17)
    shiny::observeEvent(input$apply_iq_rollup, {
      runProteinIqRollupApplyObserver(
        workflowData = workflow_data,
        experimentPaths = experiment_paths,
        output = output,
        iqRollupPlot = iq_rollup_plot,
        omicType = omic_type,
        experimentLabel = experiment_label
      )
    })
    
    # Revert IQ Rollup
    shiny::observeEvent(input$revert_iq_rollup, {
      runProteinIqRollupRevertObserver(
        workflowData = workflow_data,
        output = output
      )
    })

    bindProteinIqRollupPlot(
      output = output,
      iqRollupPlot = iq_rollup_plot
    )
  })
}
