# ----------------------------------------------------------------------------
# srlQvalueProteotypicPeptideCleanHelper
# ----------------------------------------------------------------------------
.diaNnIdentificationQvalueCutoff <- 0.01

.validateDiaNnIdentificationContract <- function(
    qvalue_threshold,
    global_qvalue_threshold,
    global_pg_qvalue_threshold,
    choose_only_proteotypic_peptide) {
  threshold_values <- suppressWarnings(as.numeric(c(
    qvalue_threshold,
    global_qvalue_threshold,
    global_pg_qvalue_threshold
  )))
  scalar_inputs <- all(vapply(
    list(
      qvalue_threshold,
      global_qvalue_threshold,
      global_pg_qvalue_threshold
    ),
    length,
    integer(1)
  ) == 1L)
  fixed_thresholds <- scalar_inputs &&
    length(threshold_values) == 3L &&
    all(is.finite(threshold_values)) &&
    all(abs(threshold_values - .diaNnIdentificationQvalueCutoff) < .Machine$double.eps^0.5)

  if (!fixed_thresholds) {
    stop(
      paste0(
        "DIA-NN identification q-value thresholds are fixed at 0.01 for ",
        "Q.Value, Global.Q.Value, and Global.PG.Q.Value."
      ),
      call. = FALSE
    )
  }

  proteotypic_value <- suppressWarnings(as.numeric(choose_only_proteotypic_peptide))
  if (length(proteotypic_value) != 1L ||
      !is.finite(proteotypic_value) ||
      proteotypic_value != 1) {
    stop(
      "DIA-NN identification evidence requires Proteotypic == 1.",
      call. = FALSE
    )
  }

  invisible(TRUE)
}

#' @title Filter Peptides by Q-Value and Proteotypic Status
#' @description Keep spectrum-peptide matches that are within the q-value
#'   thresholds and proteotypic, then freeze experiment-wide peptide and
#'   peptidoform identity counts for each protein group.
#' @param input_table Peptide quantities table in long format
#' @param qvalue_threshold Maximum precursor q-value
#' @param global_qvalue_threshold Maximum global q-value
#' @param global_pg_qvalue_threshold Maximum global protein-group q-value
#' @param choose_only_proteotypic_peptide Whether to retain only proteotypic
#'   identifications
#' @param confidence_provenance_mode Either `diann_main_report`, which requires
#'   all confidence fields, or the explicit `upstream_prefiltered` compatibility
#'   mode, which records any unavailable fields as unverified.
#' @param input_matrix_column_ids Columns retained in the filtered table
#' @param protein_id_column Protein-group identity column
#' @param peptide_sequence_column Stripped peptide identity column
#' @param modified_peptide_sequence_column Modified peptidoform identity column
#' @param q_value_column Precursor q-value column
#' @param global_q_value_column Global precursor q-value column
#' @param global_pg_q_value_column Global protein-group q-value column
#' @param proteotypic_peptide_sequence_column Proteotypic indicator column
#' @param exclude_decoys Exclude confidently classified decoys from the
#'   biological analysis view.
#' @param exclude_contaminants Exclude confidently classified contaminants from
#'   the biological analysis view.
#' @param contaminant_manifest Optional user-supplied contaminant accessions.
#' @param protein_ids_column Accession-provenance column used for classification.
#' @param return_exclusion_result Return classification/ledger metadata together
#'   with the filtered data. The default preserves the historical data-frame API.
#' @export
srlQvalueProteotypicPeptideCleanHelper <- function(input_table
                                             , qvalue_threshold = 0.01
                                             , global_qvalue_threshold = 0.01
                                             , choose_only_proteotypic_peptide = 1
                                             ,   input_matrix_column_ids = c("Run"
                                                                       , "Precursor.Id"
                                                                       , "Protein.Ids"
                                                                       , "Stripped.Sequence"
                                                                       , "Modified.Sequence"
                                                                       , "Precursor.Charge"
                                                                       , "Precursor.Quantity"
                                                                       , "Precursor.Normalised")
                                             , protein_id_column = Protein.Ids
                                             , peptide_sequence_column = Stripped.Sequence
                                             , modified_peptide_sequence_column = Modified.Sequence
                                             , q_value_column = Q.Value
                                             , global_q_value_column = Global.Q.Value
                                             , proteotypic_peptide_sequence_column = Proteotypic
                                             , global_pg_qvalue_threshold = 0.01
                                             , confidence_provenance_mode = "diann_main_report"
                                             , global_pg_q_value_column = Global.PG.Q.Value
                                             , exclude_decoys = TRUE
                                             , exclude_contaminants = TRUE
                                             , contaminant_manifest = NULL
                                             , protein_ids_column = "Protein.Ids"
                                             , return_exclusion_result = FALSE) {

  confidence_provenance_mode <- match.arg(
    confidence_provenance_mode,
    c("diann_main_report", "upstream_prefiltered")
  )
  .validateDiaNnIdentificationContract(
    qvalue_threshold = qvalue_threshold,
    global_qvalue_threshold = global_qvalue_threshold,
    global_pg_qvalue_threshold = global_pg_qvalue_threshold,
    choose_only_proteotypic_peptide = choose_only_proteotypic_peptide
  )

  qvalue_threshold <- .diaNnIdentificationQvalueCutoff
  global_qvalue_threshold <- .diaNnIdentificationQvalueCutoff
  global_pg_qvalue_threshold <- .diaNnIdentificationQvalueCutoff

  protein_id_name <- resolvePeptideQcColumnArgument(
    substitute(protein_id_column),
    protein_id_column,
    names(input_table),
    environment()
  )
  peptide_sequence_name <- resolvePeptideQcColumnArgument(
    substitute(peptide_sequence_column),
    peptide_sequence_column,
    names(input_table),
    environment()
  )
  modified_sequence_name <- resolvePeptideQcColumnArgument(
    substitute(modified_peptide_sequence_column),
    modified_peptide_sequence_column,
    names(input_table),
    environment()
  )

  q_val_name <- resolvePeptideQcColumnArgument(
    substitute(q_value_column),
    q_value_column,
    names(input_table),
    environment()
  )
  global_q_val_name <- resolvePeptideQcColumnArgument(
    substitute(global_q_value_column),
    global_q_value_column,
    names(input_table),
    environment()
  )
  global_pg_q_val_name <- resolvePeptideQcColumnArgument(
    substitute(global_pg_q_value_column),
    global_pg_q_value_column,
    names(input_table),
    environment()
  )
  proteotypic_name <- resolvePeptideQcColumnArgument(
    substitute(proteotypic_peptide_sequence_column),
    proteotypic_peptide_sequence_column,
    names(input_table),
    environment()
  )
  confidence_columns <- c(
    q_value = q_val_name,
    global_q_value = global_q_val_name,
    global_pg_q_value = global_pg_q_val_name,
    proteotypic = proteotypic_name
  )

  # [OK] DIAGNOSTIC + DEFENSIVE: Check output column availability
  missing_cols <- input_matrix_column_ids[!input_matrix_column_ids %in% names(input_table)]
  if (identical(confidence_provenance_mode, "upstream_prefiltered")) {
    missing_cols <- setdiff(missing_cols, unname(confidence_columns))
  }
  
  if (length(missing_cols) > 0) {
    error_msg <- paste0(
      "Q-value filter error: Required output columns not found in data.\n",
      "Missing columns: ", paste(missing_cols, collapse = ", "), "\n",
      "Required columns: ", paste(input_matrix_column_ids, collapse = ", "), "\n",
      "Available columns: ", paste(names(input_table), collapse = ", "), "\n\n",
      "This may be caused by:\n",
      "1. Whitespace in column names from config.ini parsing\n",
      "2. Column names with special characters or encoding issues\n",
      "3. Importing data from a different workflow stage"
    )
    logger::log_error(error_msg)
    stop(error_msg)
  }
  
  # [OK] ALSO CHECK: Filter columns exist
  filter_cols <- unname(confidence_columns)
  missing_filter_cols <- filter_cols[!filter_cols %in% names(input_table)]

  if (identical(confidence_provenance_mode, "diann_main_report") &&
      length(missing_filter_cols) > 0) {
    error_msg <- paste0(
      "Q-value filter error: Required filter columns not found in data.\n",
      "Missing filter columns: ", paste(missing_filter_cols, collapse = ", "), "\n",
      "Available columns: ", paste(names(input_table), collapse = ", ")
    )
    logger::log_error(error_msg)
    stop(error_msg)
  }

  if (identical(confidence_provenance_mode, "upstream_prefiltered") &&
      length(missing_filter_cols) > 0) {
    warning(
      paste0(
        "Explicit upstream_prefiltered mode: the following confidence fields ",
        "are unavailable and were not verified: ",
        paste(missing_filter_cols, collapse = ", "),
        "."
      ),
      call. = FALSE
    )
  }

  qvalue_filter <- rep(TRUE, nrow(input_table))
  confidence_thresholds <- c(
    q_value = qvalue_threshold,
    global_q_value = global_qvalue_threshold,
    global_pg_q_value = global_pg_qvalue_threshold
  )
  for (metric_name in names(confidence_thresholds)) {
    metric_column <- confidence_columns[[metric_name]]
    if (metric_column %in% names(input_table)) {
      qvalue_filter <- qvalue_filter &
        input_table[[metric_column]] <= confidence_thresholds[[metric_name]]
    }
  }
  if (proteotypic_name %in% names(input_table)) {
    qvalue_filter <- qvalue_filter & input_table[[proteotypic_name]] == 1
  }

  q_valid_rows <- input_table |>
    dplyr::mutate(.source_row_id = dplyr::row_number()) |>
    dplyr::filter(qvalue_filter) |>
    dplyr::select(dplyr::all_of(unique(c(
      ".source_row_id",
      input_matrix_column_ids[input_matrix_column_ids %in% names(input_table)],
      filter_cols[filter_cols %in% names(input_table)]
    ))))

  exclusion_result <- classifyPeptideBiologicalExclusions(
    input_table = q_valid_rows,
    protein_id_column = protein_id_name,
    protein_ids_column = protein_ids_column,
    contaminant_manifest = contaminant_manifest,
    exclude_decoys = .asPeptideExclusionFlag(exclude_decoys, "exclude_decoys"),
    exclude_contaminants = .asPeptideExclusionFlag(
      exclude_contaminants,
      "exclude_contaminants"
    )
  )

  search_srl_quant_cln <- .annotateProteinIdentificationEvidence(
    input_table = exclusion_result$analysis_data,
    protein_id_column = protein_id_name,
    peptide_sequence_column = peptide_sequence_name,
    modified_peptide_sequence_column = modified_sequence_name
  )

  if (isTRUE(return_exclusion_result)) {
    return(list(
      data = search_srl_quant_cln,
      q_valid_classified_data = exclusion_result$classified_data,
      exclusion_ledger = exclusion_result$exclusion_ledger,
      exclusion_summary = exclusion_result$summary,
      contaminant_manifest = exclusion_result$manifest
    ))
  }

  search_srl_quant_cln

}

