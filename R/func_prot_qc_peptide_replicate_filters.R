# ----------------------------------------------------------------------------
# removePeptidesWithOnlyOneReplicateHelper
# ----------------------------------------------------------------------------
#---------------------------------------------------------------------------------------
#' @title Remove Peptides with Single Replicate
#' @description Remove peptides that only have data for one technical replicate for all sample.
#' This can be repurposed for removing peptides that only have one biological replicates for all experimental groups.
#' This function can be repurposed for filtering on proteins as well (we just have to create a dummy variable for peptide_sequence_column)
#' @export
removePeptidesWithOnlyOneReplicateHelper <- function(input_table
                                               , samples_id_tbl
                                               , input_table_sample_id_column = Run
                                               , sample_id_tbl_sample_id_column  =  ms_filename
                                               , replicate_group_column = general_sample_info
                                               , protein_id_column = Protein.Ids
                                               , peptide_sequence_column = Stripped.Sequence
                                               , core_utilisation ) {

  input_sample_str <- resolvePeptideQcColumnArgument(substitute(input_table_sample_id_column), input_table_sample_id_column, names(input_table), environment())
  sample_tbl_sample_str <- resolvePeptideQcColumnArgument(substitute(sample_id_tbl_sample_id_column), sample_id_tbl_sample_id_column, names(samples_id_tbl), environment())
  replicate_group_str <- resolvePeptideQcColumnArgument(substitute(replicate_group_column), replicate_group_column, names(samples_id_tbl), environment())
  protein_id_str <- resolvePeptideQcColumnArgument(substitute(protein_id_column), protein_id_column, names(input_table), environment())
  peptide_seq_str <- resolvePeptideQcColumnArgument(substitute(peptide_sequence_column), peptide_sequence_column, names(input_table), environment())
  join_cols <- stats::setNames(sample_tbl_sample_str, input_sample_str)

  # Count the number of technical replicates per sample and peptide combination
  num_tech_reps_per_sample_and_peptide <- NA
  if (any(is.na(core_utilisation))) {
    num_tech_reps_per_sample_and_peptide <- input_table |>
      dplyr::left_join(samples_id_tbl, by = join_cols) |>
      dplyr::group_by(
        dplyr::across(dplyr::all_of(c(replicate_group_str, protein_id_str, peptide_seq_str)))
      ) |>
      #partition(core_utilisation) |>
      dplyr::summarise(counts = dplyr::n(), .groups = "drop") |>
      #collect() |>
      dplyr::ungroup()
  } else {
    num_tech_reps_per_sample_and_peptide <- input_table |>
      dplyr::left_join(samples_id_tbl, by = join_cols) |>
      dplyr::group_by(
        dplyr::across(dplyr::all_of(c(replicate_group_str, protein_id_str, peptide_seq_str)))
      ) |>
      partition(core_utilisation) |>
      dplyr::summarise(counts = dplyr::n(), .groups = "drop") |>
      collect() |>
      dplyr::ungroup()
  }

  # Any peptides found in more than one replicates in any patient will be kept for analysis
  removed_peptides_with_only_one_replicate <- input_table |>
    dplyr::inner_join( num_tech_reps_per_sample_and_peptide |>
                  dplyr::filter( counts >  1) |>
                  dplyr::select(-counts, -dplyr::all_of(replicate_group_str)) |>
                  dplyr::distinct()
                , by = c(protein_id_str, peptide_seq_str) )  |>
    dplyr::distinct()

  removed_peptides_with_only_one_replicate
}

# ----------------------------------------------------------------------------
# filterMinNumPeptidesPerProteinHelper
# ----------------------------------------------------------------------------
# Identification evidence must be calculated while the original stripped and
# modified sequence identities are still available.  In particular, neither a
# row count nor a sum of per-run peptidoform counts is an experiment-wide
# identity count.
.proteinIdentificationEvidenceColumns <- c(
  "identification_peptide_count",
  "identification_peptidoform_count"
)

.calculateProteinIdentificationEvidence <- function(
    input_table,
    protein_id_column = Protein.Group,
    peptide_sequence_column = Stripped.Sequence,
    modified_peptide_sequence_column = Modified.Sequence) {
  input_columns <- names(input_table)
  protein_id_str <- resolvePeptideQcColumnArgument(
    substitute(protein_id_column),
    protein_id_column,
    input_columns,
    environment()
  )
  peptide_sequence_str <- resolvePeptideQcColumnArgument(
    substitute(peptide_sequence_column),
    peptide_sequence_column,
    input_columns,
    environment()
  )
  modified_sequence_str <- resolvePeptideQcColumnArgument(
    substitute(modified_peptide_sequence_column),
    modified_peptide_sequence_column,
    input_columns,
    environment()
  )

  input_table |>
    dplyr::filter(
      !is.na(.data[[protein_id_str]]),
      trimws(as.character(.data[[protein_id_str]])) != ""
    ) |>
    dplyr::mutate(
      .identification_peptide_identity = trimws(
        as.character(.data[[peptide_sequence_str]])
      ),
      .identification_peptidoform_identity = trimws(
        as.character(.data[[modified_sequence_str]])
      )
    ) |>
    dplyr::group_by(.data[[protein_id_str]]) |>
    dplyr::summarise(
      identification_peptide_count = dplyr::n_distinct(
        .data$.identification_peptide_identity[
          !is.na(.data$.identification_peptide_identity) &
            nzchar(.data$.identification_peptide_identity)
        ]
      ),
      identification_peptidoform_count = dplyr::n_distinct(
        .data$.identification_peptidoform_identity[
          !is.na(.data$.identification_peptidoform_identity) &
            nzchar(.data$.identification_peptidoform_identity)
        ]
      ),
      .groups = "drop"
    )
}

.annotateProteinIdentificationEvidence <- function(
    input_table,
    protein_id_column = Protein.Group,
    peptide_sequence_column = Stripped.Sequence,
    modified_peptide_sequence_column = Modified.Sequence) {
  protein_id_str <- resolvePeptideQcColumnArgument(
    substitute(protein_id_column),
    protein_id_column,
    names(input_table),
    environment()
  )
  peptide_sequence_str <- resolvePeptideQcColumnArgument(
    substitute(peptide_sequence_column),
    peptide_sequence_column,
    names(input_table),
    environment()
  )
  modified_sequence_str <- resolvePeptideQcColumnArgument(
    substitute(modified_peptide_sequence_column),
    modified_peptide_sequence_column,
    names(input_table),
    environment()
  )

  evidence <- .calculateProteinIdentificationEvidence(
    input_table = input_table,
    protein_id_column = protein_id_str,
    peptide_sequence_column = peptide_sequence_str,
    modified_peptide_sequence_column = modified_sequence_str
  )

  input_table |>
    dplyr::select(-dplyr::any_of(.proteinIdentificationEvidenceColumns)) |>
    dplyr::left_join(evidence, by = protein_id_str)
}

#' @export
#' @title Filter Proteins by Minimum Distinct Peptide Evidence
#' @description Keep protein groups only when they meet the configured numbers
#'   of distinct stripped peptide sequences and distinct modified peptidoforms.
#'   Frozen identification counts created at the q-value filtering stage are
#'   preferred so later quantitative filtering cannot change identification
#'   support.
#' @param input_table Peptide quantities table in long format
#' @param num_peptides_per_protein_thresh Minimum number of peptides per protein
#' @param num_peptidoforms_per_protein_thresh Minimum number of peptidoforms per protein
#' @param protein_id_column Protein ID column name as string
#' @param peptide_sequence_column Peptide sequence column name
#' @param modified_peptide_sequence_column Modified peptide sequence column name
#' @param peptidoform_ids_column List-column containing the modified peptide
#'   identities retained by an older precursor-to-peptide rollup. This is kept
#'   only for compatibility with already-created in-memory objects.
#' @param core_utilisation core_utilisation to use for parallel processing
filterMinNumPeptidesPerProteinHelper <- function( input_table
          , num_peptides_per_protein_thresh = 1
          , num_peptidoforms_per_protein_thresh = 2
          , protein_id_column = Protein.Ids
          , peptide_sequence_column = Stripped.Sequence
          , modified_peptide_sequence_column = Modified.Sequence
          , peptidoform_ids_column = "peptidoform_ids"
          , core_utilisation) {

  protein_id_str <- resolvePeptideQcColumnArgument(substitute(protein_id_column), protein_id_column, names(input_table), environment())
  peptide_sequence_str <- resolvePeptideQcColumnArgument(substitute(peptide_sequence_column), peptide_sequence_column, names(input_table), environment())
  modified_sequence_str <- resolvePeptideQcColumnArgument(substitute(modified_peptide_sequence_column), modified_peptide_sequence_column, names(input_table), environment())

  frozen_columns_present <- .proteinIdentificationEvidenceColumns %in% names(input_table)
  if (any(frozen_columns_present) && !all(frozen_columns_present)) {
    stop(
      "filterMinNumPeptidesPerProtein: frozen identification evidence is incomplete.",
      call. = FALSE
    )
  }

  if (all(frozen_columns_present)) {
    frozen_evidence <- input_table |>
      dplyr::select(
        dplyr::all_of(protein_id_str),
        dplyr::all_of(.proteinIdentificationEvidenceColumns)
      ) |>
      dplyr::distinct()

    conflicting_evidence <- frozen_evidence |>
      dplyr::count(.data[[protein_id_str]], name = ".evidence_rows") |>
      dplyr::filter(.data$.evidence_rows != 1L)

    if (nrow(conflicting_evidence) > 0L) {
      stop(
        "filterMinNumPeptidesPerProtein: conflicting frozen identification evidence for one or more protein groups.",
        call. = FALSE
      )
    }

    num_peptides_per_protein <- frozen_evidence |>
      dplyr::transmute(
        !!protein_id_str := .data[[protein_id_str]],
        peptides_for_protein_count = as.integer(
          .data$identification_peptide_count
        ),
        peptidoforms_for_protein_count = as.integer(
          .data$identification_peptidoform_count
        )
      )
  } else if (modified_sequence_str %in% names(input_table)) {
    exact_evidence <- .calculateProteinIdentificationEvidence(
      input_table = input_table,
      protein_id_column = protein_id_str,
      peptide_sequence_column = peptide_sequence_str,
      modified_peptide_sequence_column = modified_sequence_str
    )

    num_peptides_per_protein <- exact_evidence |>
      dplyr::transmute(
        !!protein_id_str := .data[[protein_id_str]],
        peptides_for_protein_count = as.integer(
          .data$identification_peptide_count
        ),
        peptidoforms_for_protein_count = as.integer(
          .data$identification_peptidoform_count
        )
      )
  } else if (peptidoform_ids_column %in% names(input_table)) {
    warning(
      "Using legacy peptidoform_ids provenance. Re-run q-value filtering to create frozen identification evidence.",
      call. = FALSE
    )
    peptide_counts <- input_table |>
      dplyr::distinct(.data[[protein_id_str]], .data[[peptide_sequence_str]]) |>
      dplyr::count(
        .data[[protein_id_str]],
        name = "peptides_for_protein_count"
      )

    peptidoform_counts <- input_table |>
      dplyr::select(
        dplyr::all_of(protein_id_str),
        dplyr::all_of(peptidoform_ids_column)
      ) |>
      tidyr::unnest_longer(
        dplyr::all_of(peptidoform_ids_column),
        values_to = ".peptidoform_id",
        keep_empty = FALSE
      ) |>
      dplyr::filter(
        !is.na(.data$.peptidoform_id),
        trimws(as.character(.data$.peptidoform_id)) != ""
      ) |>
      dplyr::distinct(.data[[protein_id_str]], .data$.peptidoform_id) |>
      dplyr::count(
        .data[[protein_id_str]],
        name = "peptidoforms_for_protein_count"
      )

    num_peptides_per_protein <- peptide_counts |>
      dplyr::left_join(peptidoform_counts, by = protein_id_str)
  } else {
    stop(
      paste0(
        "filterMinNumPeptidesPerProtein: exact experiment-wide peptide evidence is unavailable. ",
        "Run q-value filtering before precursor rollup, or supply Modified.Sequence."
      ),
      call. = FALSE
    )
  }

  num_peptides_per_protein <- num_peptides_per_protein |>
    dplyr::mutate(
      peptidoforms_for_protein_count = dplyr::coalesce(
        .data$peptidoforms_for_protein_count,
        0L
      )
    )

  input_table_for_join <- input_table |>
    dplyr::select(
      -dplyr::any_of(c(
        "peptides_for_protein_count",
        "peptidoforms_for_protein_count"
      ))
    )

  protein_peptide_cln <- NA
  valid_thresholds <- length(num_peptides_per_protein_thresh) == 1L &&
    length(num_peptidoforms_per_protein_thresh) == 1L &&
    is.finite(num_peptides_per_protein_thresh) &&
    is.finite(num_peptidoforms_per_protein_thresh) &&
    num_peptides_per_protein_thresh >= 1 &&
    num_peptidoforms_per_protein_thresh >= 1

  if (valid_thresholds) {

    protein_peptide_cln <- input_table_for_join |>
      inner_join( num_peptides_per_protein
                  , by = protein_id_str) |>
      dplyr::filter(   peptidoforms_for_protein_count >= num_peptidoforms_per_protein_thresh
                      ,
                      peptides_for_protein_count >= num_peptides_per_protein_thresh
                     )
  } else {
    stop(
      paste0(
        "filterMinNumPeptidesPerProtein: num_peptides_per_protein_thresh and ",
        "num_peptidoforms_per_protein_thresh must be provided as finite values >= 1."
      )
    )
  }

  protein_peptide_cln
}

# ----------------------------------------------------------------------------
# filterMinNumPeptidesPerSampleHelper
# ----------------------------------------------------------------------------
#' @export
#' @title Filter Samples by Minimum Number of Peptides
#' @description Remove sample if it has less than a certain number of peptides identified
#' @param List of samples to keep regardless of how many peptides it has because it is am important sample
filterMinNumPeptidesPerSampleHelper <- function ( input_table
                                            , peptides_per_sample_cutoff = 5000
                                            , sample_id_column = Run
                                            , core_utilisation
                                            , inclusion_list = c()) {

  sample_id_str <- resolvePeptideQcColumnArgument(substitute(sample_id_column), sample_id_column, names(input_table), environment())
  samples_passing_filter <- NA
  if (any(is.na(core_utilisation))) {
    samples_passing_filter <- input_table |>
      dplyr::group_by(.data[[sample_id_str]]) |>
      #partition(core_utilisation) |>
      summarise( counts = n()) |>
      #collect() |>
      ungroup() |>
      dplyr::filter( counts >= peptides_per_sample_cutoff |
                       ( .data[[sample_id_str]] %in% inclusion_list)  ) |>
      dplyr::select(-counts)

  } else {
    samples_passing_filter <- input_table |>
      dplyr::group_by(.data[[sample_id_str]]) |>
      partition(core_utilisation) |>
      summarise( counts = n()) |>
      collect() |>
      ungroup() |>
      dplyr::filter( counts >= peptides_per_sample_cutoff |
                       ( .data[[sample_id_str]] %in% inclusion_list)  ) |>
      dplyr::select(-counts)
  }

  filtered_table <- input_table |>
    inner_join( samples_passing_filter, by = sample_id_str)

  filtered_table
}

# ----------------------------------------------------------------------------
# srlQvalueProteotypicPeptideCleanHelper
# ----------------------------------------------------------------------------
#' @title Filter Peptides by Q-Value and Proteotypic Status
#' @description Keep spectrum-peptide matches that are within the q-value
#'   thresholds and proteotypic, then freeze experiment-wide peptide and
#'   peptidoform identity counts for each protein group.
#' @param input_table Peptide quantities table in long format
#' @param qvalue_threshold Maximum precursor q-value
#' @param global_qvalue_threshold Maximum global q-value
#' @param choose_only_proteotypic_peptide Whether to retain only proteotypic
#'   identifications
#' @param input_matrix_column_ids Columns retained in the filtered table
#' @param protein_id_column Protein-group identity column
#' @param peptide_sequence_column Stripped peptide identity column
#' @param modified_peptide_sequence_column Modified peptidoform identity column
#' @param q_value_column Precursor q-value column
#' @param global_q_value_column Global q-value column
#' @param proteotypic_peptide_sequence_column Proteotypic indicator column
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
                                             , proteotypic_peptide_sequence_column = Proteotypic) {

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

  # [OK] DIAGNOSTIC + DEFENSIVE: Check output column availability
  missing_cols <- input_matrix_column_ids[!input_matrix_column_ids %in% names(input_table)]
  
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
  q_val_name <- resolvePeptideQcColumnArgument(substitute(q_value_column), q_value_column, names(input_table), environment())
  global_q_val_name <- resolvePeptideQcColumnArgument(substitute(global_q_value_column), global_q_value_column, names(input_table), environment())
  proteotypic_name <- resolvePeptideQcColumnArgument(substitute(proteotypic_peptide_sequence_column), proteotypic_peptide_sequence_column, names(input_table), environment())
  
  use_proteotypic_filter <- isTRUE(
    suppressWarnings(as.numeric(choose_only_proteotypic_peptide)) == 1
  )
  filter_cols <- c(
    q_val_name,
    global_q_val_name,
    if (use_proteotypic_filter) proteotypic_name else character()
  )
  missing_filter_cols <- filter_cols[!filter_cols %in% names(input_table)]
  
  if (length(missing_filter_cols) > 0) {
    error_msg <- paste0(
      "Q-value filter error: Required filter columns not found in data.\n",
      "Missing filter columns: ", paste(missing_filter_cols, collapse = ", "), "\n",
      "Available columns: ", paste(names(input_table), collapse = ", ")
    )
    logger::log_error(error_msg)
    stop(error_msg)
  }

  qvalue_filter <- input_table[[q_val_name]] <= qvalue_threshold &
    input_table[[global_q_val_name]] <= global_qvalue_threshold
  if (use_proteotypic_filter) {
    qvalue_filter <- qvalue_filter &
      input_table[[proteotypic_name]] == choose_only_proteotypic_peptide
  }

  search_srl_quant_cln <- input_table |>
    dplyr::filter(qvalue_filter) |>
    dplyr::select(all_of(unique(c(input_matrix_column_ids, filter_cols))))

  search_srl_quant_cln <- .annotateProteinIdentificationEvidence(
    input_table = search_srl_quant_cln,
    protein_id_column = protein_id_name,
    peptide_sequence_column = peptide_sequence_name,
    modified_peptide_sequence_column = modified_sequence_name
  )

  search_srl_quant_cln

}
