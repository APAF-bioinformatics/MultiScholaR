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
#' @param id,workflowData,experimentPaths,createProteinDataFn,writeTsvFn,processLongFormatFn,readTsvFn,captureCheckpointFn,showNotificationFn,logInfoFn,logWarnFn,sleepFn,maxWait Runtime inputs used by this function; see the usage section for accepted values.
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
          shiny::radioButtons(
            ns("rollup_method"),
            "Protein rollup method:",
            choices = c(
              "IQ MaxLFQ" = "iq",
              "limpa DPC-Quant" = "limpa"
            ),
            selected = "iq"
          ),
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

buildProtTestModeLimpaProteinObject <- function(peptideS4,
                                                dpcSlope = 0.8,
                                                createProteinDataFn = ProteinQuantitativeData) {
  peptideData <- peptideS4@peptide_data
  proteinCol <- peptideS4@protein_id_column
  sampleCol <- peptideS4@sample_id
  peptideCol <- peptideS4@peptide_sequence_column
  quantityCol <- peptideS4@norm_quantity_column

  if (is.null(proteinCol) || !proteinCol %in% names(peptideData)) {
    proteinCandidates <- c("Protein.Ids", "Protein.Group", "Protein.IDs", "protein_id")
    proteinCol <- proteinCandidates[proteinCandidates %in% names(peptideData)][1]
  }
  if (is.na(proteinCol) || is.null(proteinCol)) {
    stop("No protein identifier column available for test-mode limpa quantification.")
  }

  if (is.null(sampleCol) || !sampleCol %in% names(peptideData)) {
    sampleCandidates <- c("Run", "Sample", "sample_id")
    sampleCol <- sampleCandidates[sampleCandidates %in% names(peptideData)][1]
  }
  if (is.na(sampleCol) || is.null(sampleCol)) {
    stop("No sample identifier column available for test-mode limpa quantification.")
  }

  if (is.null(peptideCol) || !peptideCol %in% names(peptideData)) {
    peptideCandidates <- c("Stripped.Sequence", "Modified.Sequence", "Precursor.Id", "Peptide.Sequence")
    peptideCol <- peptideCandidates[peptideCandidates %in% names(peptideData)][1]
  }
  if (is.na(peptideCol) || is.null(peptideCol)) {
    peptideCol <- proteinCol
  }

  if (is.null(quantityCol) || !quantityCol %in% names(peptideData)) {
    quantityCandidates <- c("Peptide.Imputed", "Precursor.Quantity", "Precursor.Normalised")
    quantityCol <- quantityCandidates[quantityCandidates %in% names(peptideData)][1]
  }
  if (is.na(quantityCol) || is.null(quantityCol)) {
    stop("No peptide quantity column available for test-mode limpa quantification.")
  }

  proteinLong <- peptideData |>
    dplyr::group_by(!!rlang::sym(proteinCol), !!rlang::sym(sampleCol)) |>
    dplyr::summarise(
      value = mean(as.numeric(!!rlang::sym(quantityCol)), na.rm = TRUE),
      observations = sum(!is.na(!!rlang::sym(quantityCol))),
      .groups = "drop"
    ) |>
    dplyr::mutate(value = ifelse(is.nan(value), NA_real_, value))

  proteinWide <- proteinLong |>
    dplyr::select(
      !!rlang::sym(proteinCol),
      !!rlang::sym(sampleCol),
      value
    ) |>
    tidyr::pivot_wider(
      names_from = !!rlang::sym(sampleCol),
      values_from = value
    )

  sampleNames <- unique(as.character(peptideS4@design_matrix[[sampleCol]]))
  if (length(sampleNames) == 0L || all(is.na(sampleNames))) {
    sampleNames <- unique(as.character(peptideData[[sampleCol]]))
  }
  missingSamples <- setdiff(sampleNames, colnames(proteinWide))
  for (sampleName in missingSamples) {
    proteinWide[[sampleName]] <- NA_real_
  }
  proteinWide <- proteinWide[, c(proteinCol, sampleNames), drop = FALSE]

  valueCols <- setdiff(colnames(proteinWide), proteinCol)
  for (rowIdx in seq_len(nrow(proteinWide))) {
    values <- as.numeric(proteinWide[rowIdx, valueCols, drop = TRUE])
    replacement <- if (all(is.na(values))) 0 else mean(values, na.rm = TRUE)
    values[is.na(values)] <- replacement
    proteinWide[rowIdx, valueCols] <- as.list(values)
  }

  peptideSummary <- peptideData |>
    dplyr::group_by(!!rlang::sym(proteinCol)) |>
    dplyr::summarise(
      peptide_count = dplyr::n_distinct(!!rlang::sym(peptideCol)),
      peptidoform_count = dplyr::n(),
      .groups = "drop"
    )

  proteinMatrix <- as.matrix(proteinWide[, valueCols, drop = FALSE])
  rownames(proteinMatrix) <- proteinWide[[proteinCol]]
  quantifiedElist <- list(
    E = proteinMatrix,
    genes = data.frame(
      protein.id = proteinWide[[proteinCol]],
      stringsAsFactors = FALSE
    ),
    other = list(
      standard.error = matrix(0, nrow = nrow(proteinMatrix), ncol = ncol(proteinMatrix)),
      n.observations = matrix(1, nrow = nrow(proteinMatrix), ncol = ncol(proteinMatrix))
    )
  )
  class(quantifiedElist) <- "EList"

  args <- peptideS4@args
  if (is.null(args$globalParameters)) {
    args$globalParameters <- list()
  }
  peptideMatrix <- tryCatch(peptideS4@peptide_matrix, error = function(e) NULL)
  if (is.null(peptideMatrix)) {
    peptideMatrix <- as.matrix(proteinWide[, valueCols, drop = FALSE])
    rownames(peptideMatrix) <- proteinWide[[proteinCol]]
  }
  args$globalParameters$use_limpa <- TRUE
  args$globalParameters$report_template <- "DIANN_limpa_report.rmd"
  args$proteinMissingValueImputationLimpa <- list(
    dpc_slope = dpcSlope,
    verbose = FALSE,
    quantified_protein_column = "Protein.Quantified.Limpa",
    test_mode = TRUE
  )
  args$limpa_dpc_quant_results <- list(
    dpc_parameters_used = c(NA_real_, dpcSlope),
    dpc_object_used = NULL,
    slope_interpretation = "test-mode deterministic fallback",
    y_peptide_for_dpc = peptideMatrix,
    quantified_elist = quantifiedElist,
    standard_errors = quantifiedElist$other$standard.error,
    n_observations = quantifiedElist$other$n.observations,
    peptide_counts_per_protein = peptideSummary,
    missing_percentage_peptides = round(100 * mean(is.na(peptideMatrix)), 1),
    missing_percentage_proteins = 0,
    quantification_method = "limpa_dpc_quant_test_mode",
    total_proteins_quantified = nrow(proteinMatrix),
    total_peptides_used = nrow(peptideMatrix),
    dpc_method = "limpa_dpc_quant"
  )

  createProteinDataFn(
    protein_quant_table = proteinWide,
    protein_id_column = proteinCol,
    protein_id_table = .proteinIdTableFromPeptideLineage(
      peptideS4,
      proteinWide[[proteinCol]]
    ),
    design_matrix = peptideS4@design_matrix,
    sample_id = sampleCol,
    group_id = peptideS4@group_id,
    technical_replicate_id = peptideS4@technical_replicate_id,
    args = args
  )
}

ensurePeptideMatrixForLimpaRollup <- function(peptideS4,
                                              calcPeptideMatrixFn = calcPeptideMatrix) {
  if (!methods::is(peptideS4, "PeptideQuantitativeData")) {
    return(peptideS4)
  }

  if (length(peptideS4@peptide_matrix) == 0L) {
    peptideS4 <- calcPeptideMatrixFn(peptideS4)
  }

  peptideS4
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
  if (protDiaQcWorkerEligible(workflowData)) {
    worker <- runProtDiaQcProcess(protDiaQcWorkerSpec(
      workflowData,
      "protein_iq_rollup",
      list(experimentPaths = experimentPaths, maxWait = maxWait)
    ))
    result <- applyProtDiaQcWorkerResult(workflowData, worker)
    result$plot_png <- worker$plot_png
    return(result)
  }

  currentState <- workflowStateCurrentName(workflowData$state_manager)
  peptideS4 <- workflowData$state_manager$getState(currentState)
  shiny::req(peptideS4)

  logInfoFn("Protein Processing: Starting IQ rollup from peptide state")

  proteinColumn <- peptideS4@protein_id_column
  peptideColumn <- peptideS4@peptide_sequence_column
  sampleColumn <- peptideS4@sample_id
  groupColumn <- peptideS4@group_id
  replicateColumn <- peptideS4@technical_replicate_id
  intensityColumn <- "Peptide.Imputed"
  requiredPeptideColumns <- c(
    proteinColumn,
    peptideColumn,
    sampleColumn,
    intensityColumn
  )
  missingPeptideColumns <- setdiff(
    requiredPeptideColumns,
    names(peptideS4@peptide_data)
  )
  if (length(missingPeptideColumns) > 0L) {
    stop(
      "IQ rollup input is missing required column(s): ",
      paste(missingPeptideColumns, collapse = ", "),
      call. = FALSE
    )
  }
  if (!sampleColumn %in% names(peptideS4@design_matrix)) {
    stop("The declared sample column is absent from the design matrix.", call. = FALSE)
  }

  peptideValuesImputedFile <- file.path(
    experimentPaths$peptide_qc_dir,
    "peptide_values_imputed.tsv"
  )

  originalSamples <- unique(as.character(peptideS4@design_matrix[[sampleColumn]]))
  sampleMapping <- tibble::tibble(
    Original = originalSamples,
    Alias = paste0("S_", sprintf("%03d", seq_along(originalSamples)))
  )

  peptideDataForIq <- peptideS4@peptide_data
  peptideDataForIq$Q.Value <- 0.0009
  peptideDataForIq$PG.Q.Value <- 0.009
  peptideDataForIq[[intensityColumn]][is.na(peptideDataForIq[[intensityColumn]])] <- 0
  sampleIndex <- match(as.character(peptideDataForIq[[sampleColumn]]), sampleMapping$Original)
  if (anyNA(sampleIndex)) {
    stop("Peptide data contains runs absent from the design matrix.", call. = FALSE)
  }
  peptideDataForIq[[sampleColumn]] <- sampleMapping$Alias[sampleIndex]

  writeTsvFn(peptideDataForIq, peptideValuesImputedFile)

  iqOutputFile <- file.path(experimentPaths$protein_qc_dir, "iq_output_file.txt")

  processLongFormatFn(
    peptideValuesImputedFile,
    output_filename = iqOutputFile,
    sample_id = sampleColumn,
    primary_id = proteinColumn,
    secondary_id = peptideColumn,
    intensity_col = intensityColumn,
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
  iqOutputWasParserTable <- inherits(proteinLog2QuantAliased, "spec_tbl_df")
  proteinLog2QuantAliased <- protDiaArtifactPortableTable(
    proteinLog2QuantAliased,
    "IQ MaxLFQ output"
  )
  currentCols <- colnames(proteinLog2QuantAliased)
  restoredCols <- purrr::map_chr(currentCols, function(col) {
    if (col == proteinColumn) {
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

  if (!proteinColumn %in% names(proteinLog2Quant)) {
    stop(
      sprintf("IQ output is missing the active protein key `%s`.", proteinColumn),
      call. = FALSE
    )
  }

  survivingSamples <- intersect(
    originalSamples,
    setdiff(colnames(proteinLog2Quant), proteinColumn)
  )
  finalDesignMatrix <- peptideS4@design_matrix |>
    dplyr::filter(!!rlang::sym(sampleColumn) %in% survivingSamples)

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
    protein_id_column = proteinColumn,
    protein_id_table = .proteinIdTableFromPeptideLineage(
      peptideS4,
      proteinLog2Quant[[proteinColumn]]
    ),
    design_matrix = finalDesignMatrix,
    sample_id = sampleColumn,
    group_id = groupColumn,
    technical_replicate_id = replicateColumn,
    args = peptideS4@args
  )

  iqStateConfig <- list(
      iq_output_file = iqOutputFile,
      peptide_input_file = peptideValuesImputedFile,
      s4_class = "ProteinQuantitativeData",
      protein_id_column = proteinColumn
  )
  proteinObj <- saveProtProteinQcState(
    workflow_data = workflowData,
    state_manager = workflowData$state_manager,
    before = peptideS4,
    after = proteinObj,
    stage_id = "protein_rollup",
    state_name = "protein_s4_created",
    config_object = iqStateConfig,
    audit_parameters = list(
      rollup_method = "iq_maxlfq",
      primary_id = proteinColumn,
      secondary_id = peptideColumn,
      sample_id = sampleColumn,
      intensity_column = intensityColumn,
      q_value_compatibility_value = 0.0009,
      pg_q_value_compatibility_value = 0.009,
      missing_intensity_compatibility_value = 0,
      filter_double_less = c(Q.Value = 0.01, PG.Q.Value = 0.01),
      normalization = "none",
      parser_metadata_removed = iqOutputWasParserTable,
      dropped_samples = droppedSamples
    ),
    description = "IQ protein rollup completed and ProteinQuantitativeData S4 object created",
    transformation_type = "aggregation"
  )

  captureCheckpointFn(proteinObj, "cp03", "rolled_up_protein")

  proteinCount <- proteinObj@protein_quant_table |>
    dplyr::distinct(!!rlang::sym(proteinColumn)) |>
    nrow()

  resultText <- paste(
    "IQ Protein Rollup & S4 Object Creation Completed Successfully\n",
    "============================================================\n",
    sprintf("Proteins quantified: %d\n", proteinCount),
    sprintf("Samples: %d\n", ncol(proteinObj@protein_quant_table) - 1),
    sprintf("Active protein key: %s\n", proteinColumn),
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

runProteinLimpaRollupApplyStep <- function(workflowData,
                                           experimentPaths,
                                           proteinLimpaFn = proteinMissingValueImputationLimpa,
                                           testModeFallbackFn = buildProtTestModeLimpaProteinObject,
                                           preparePeptideMatrixFn = ensurePeptideMatrixForLimpaRollup,
                                           requireNamespaceFn = requireNamespace,
                                           isTestModeFn = is_test_mode,
                                           captureCheckpointFn = .capture_checkpoint,
                                           logInfoFn = logger::log_info) {
  shiny::req(workflowData$state_manager)
  if (protDiaQcWorkerEligible(workflowData)) {
    worker <- runProtDiaQcProcess(protDiaQcWorkerSpec(
      workflowData,
      "protein_limpa_rollup",
      list(experimentPaths = experimentPaths)
    ))
    result <- applyProtDiaQcWorkerResult(workflowData, worker)
    result$plot_png <- worker$plot_png
    return(result)
  }

  currentState <- workflowStateCurrentName(workflowData$state_manager)
  selectedPeptideS4 <- workflowData$state_manager$getState(currentState)
  shiny::req(selectedPeptideS4)
  peptideS4 <- preparePeptideMatrixFn(selectedPeptideS4)

  logInfoFn("Protein Processing: Starting limpa DPC-Quant rollup from peptide state")

  proteinObj <- if (isTRUE(requireNamespaceFn("limpa", quietly = TRUE))) {
    proteinLimpaFn(peptideS4, verbose = TRUE)
  } else if (isTRUE(isTestModeFn())) {
    testModeFallbackFn(peptideS4)
  } else {
    stop("limpa package is required for limpa DPC-Quant rollup. Install it with BiocManager::install('limpa').")
  }

  if (is.null(workflowData$config_list)) {
    workflowData$config_list <- list()
  }
  if (is.null(workflowData$config_list$globalParameters)) {
    workflowData$config_list$globalParameters <- list()
  }
  workflowData$config_list$globalParameters$use_limpa <- TRUE
  workflowData$config_list$globalParameters$report_template <- "DIANN_limpa_report.rmd"
  if (is.null(proteinObj@args$globalParameters)) {
    proteinObj@args$globalParameters <- list()
  }
  proteinObj@args$globalParameters$use_limpa <- TRUE
  proteinObj@args$globalParameters$report_template <- "DIANN_limpa_report.rmd"

  limpaConfig <- list(
      rollup_method = "limpa_dpc_quant",
      s4_class = "ProteinQuantitativeData",
      protein_id_column = proteinObj@protein_id_column,
      report_template = "DIANN_limpa_report.rmd"
  )
  limpaResults <- proteinObj@args$limpa_dpc_quant_results
  limpaParameters <- proteinObj@args$proteinMissingValueImputationLimpa
  proteinObj <- saveProtProteinQcState(
    workflow_data = workflowData,
    state_manager = workflowData$state_manager,
    before = selectedPeptideS4,
    after = proteinObj,
    stage_id = "protein_rollup",
    state_name = "protein_s4_created",
    config_object = limpaConfig,
    audit_parameters = list(
      rollup_method = "limpa_dpc_quant",
      active_protein_key = resolveProteinQuantIdentityColumn(proteinObj),
      peptide_matrix_prepared = !identical(selectedPeptideS4, peptideS4),
      peptide_matrix_dimensions = dim(peptideS4@peptide_matrix),
      dpc_slope = limpaParameters$dpc_slope,
      dpc_parameters_used = limpaResults$dpc_parameters_used,
      quantification_method = limpaResults$quantification_method,
      dpc_method = limpaResults$dpc_method,
      missing_percentage_peptides = limpaResults$missing_percentage_peptides,
      missing_percentage_proteins = limpaResults$missing_percentage_proteins,
      total_proteins_quantified = limpaResults$total_proteins_quantified,
      total_peptides_used = limpaResults$total_peptides_used,
      report_template = "DIANN_limpa_report.rmd",
      test_mode = isTRUE(limpaParameters$test_mode)
    ),
    description = "limpa DPC-Quant protein rollup completed and ProteinQuantitativeData S4 object created",
    transformation_type = "aggregation"
  )

  captureCheckpointFn(proteinObj, "cp03", "rolled_up_protein_limpa")

  proteinCount <- proteinObj@protein_quant_table |>
    dplyr::distinct(!!rlang::sym(proteinObj@protein_id_column)) |>
    nrow()

  resultText <- paste(
    "limpa DPC-Quant Protein Rollup Completed Successfully\n",
    "=====================================================\n",
    sprintf("Proteins quantified: %d\n", proteinCount),
    sprintf("Samples: %d\n", ncol(proteinObj@protein_quant_table) - 1),
    "Algorithm: limpa DPC-Quant\n",
    sprintf("S4 Class: %s\n", class(proteinObj)[1]),
    "Report template: DIANN_limpa_report.rmd\n",
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
                                         updateProteinFilteringFn = updateProteinFiltering,
                                         workflowData = NULL) {
  output$iq_rollup_results <- renderTextFn(iqRollupResult$resultText)

  plotGrid <- if (is.raw(iqRollupResult$plot_png)) {
    protDiaQcPlotGrid(iqRollupResult$plot_png)
  } else {
    protDiaProteinQcUpdateFiltering(
      update_fn = updateProteinFilteringFn,
      workflow_data = workflowData,
      data = iqRollupResult$proteinObj@protein_quant_table,
      step_name = "9_protein_s4_created",
      omic_type = omicType,
      experiment_label = experimentLabel,
      return_grid = TRUE,
      overwrite = TRUE
    )
  }
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
                                            renderTextFn = shiny::renderText,
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

    plotGrid <- protDiaProteinQcInvokeOutputs(
      update_outputs_fn = updateOutputsFn,
      args = list(
        output = output,
        iqRollupPlot = iqRollupPlot,
        iqRollupResult = iqRollupResult,
        omicType = omicType,
        experimentLabel = experimentLabel
      ),
      workflow_data = workflowData
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
    output$iq_rollup_results <- renderTextFn(errorMessage)
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
  if (protDiaQcWorkerEligible(workflowData)) {
    return(runProtDiaQcRevert(workflowData, previousState))
  }

  revertedS4 <- revertProtDiaProteinQcState(workflowData, previousState)

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
    
    # Step 1: Protein Rollup (chunk 17)
    shiny::observeEvent(input$apply_iq_rollup, {
      rollupMethod <- input$rollup_method %||% "iq"
      applyStepFn <- if (identical(rollupMethod, "limpa")) {
        runProteinLimpaRollupApplyStep
      } else {
        runProteinIqRollupApplyStep
      }

      runProteinIqRollupApplyObserver(
        workflowData = workflow_data,
        experimentPaths = experiment_paths,
        output = output,
        iqRollupPlot = iq_rollup_plot,
        omicType = omic_type,
        experimentLabel = experiment_label,
        runApplyStepFn = applyStepFn
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
