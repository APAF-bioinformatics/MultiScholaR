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
#' @export
#' @title Filter Proteins by Minimum Number of Peptides
#' @description Keep the proteins only if they have two or more peptides.
#' @param input_table Peptide quantities table in long format
#' @param num_peptides_per_protein_thresh Minimum number of peptides per protein
#' @param num_peptidoforms_per_protein_thresh Minimum number of peptidoforms per protein
#' @param protein_id_column Protein ID column name as string
#' @param core_utilisation core_utilisation to use for parallel processing
filterMinNumPeptidesPerProteinHelper <- function( input_table
          , num_peptides_per_protein_thresh = 1
          , num_peptidoforms_per_protein_thresh = 2
          , protein_id_column = Protein.Ids
          , core_utilisation) {

  protein_id_str <- resolvePeptideQcColumnArgument(substitute(protein_id_column), protein_id_column, names(input_table), environment())
  num_peptides_per_protein <- NA
  if (any(is.na(core_utilisation))) {
    num_peptides_per_protein <- input_table |>
      dplyr::group_by(.data[[protein_id_str]]) |>
      dplyr::summarise( peptides_for_protein_count = n()
                 , peptidoforms_for_protein_count = sum( peptidoform_count, na.rm=TRUE)) |>
      ungroup()
  } else {
    num_peptides_per_protein <- input_table |>
      dplyr::group_by(.data[[protein_id_str]]) |>
      partition(core_utilisation) |>
      dplyr::summarise( peptides_for_protein_count = n()
                 , peptidoforms_for_protein_count = sum( peptidoform_count, na.rm=TRUE)) |>
      collect() |>
      ungroup()
  }

  protein_peptide_cln <- NA
  if ( !is.na(num_peptides_per_protein_thresh) &
       !is.na(num_peptidoforms_per_protein_thresh )  ) {

    print(num_peptides_per_protein)

    protein_peptide_cln <- input_table |>
      inner_join( num_peptides_per_protein
                  , by = protein_id_str) |>
      dplyr::filter(   peptidoforms_for_protein_count >= num_peptidoforms_per_protein_thresh
                      ,
                      peptides_for_protein_count >= num_peptides_per_protein_thresh
                     )
  } else {
    stop("filterMinNumPeptidesPerProtein: num_peptides_per_protein_thresh and num_peptidoforms_per_protein_thresh must be provided.")
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
#' @description Keep spectrum-peptide matches that is within q-value threshold and are proteotypic
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
                                             , q_value_column = Q.Value
                                             , global_q_value_column = Global.Q.Value
                                             , proteotypic_peptide_sequence_column = Proteotypic) {

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

  search_srl_quant_cln

}
