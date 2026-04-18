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

  # Count the number of technical replicates per sample and peptide combination
  num_tech_reps_per_sample_and_peptide <- NA
  if( length(which(is.na(core_utilisation))) == 0 ) {
    num_tech_reps_per_sample_and_peptide <- input_table |>
      left_join( samples_id_tbl, by=join_by( {{input_table_sample_id_column}} == {{sample_id_tbl_sample_id_column}} ) ) |>
      group_by( {{replicate_group_column}}, {{protein_id_column}}, {{peptide_sequence_column}}) |>
      #partition(core_utilisation) |>
      summarise(counts = n() ) |>
      #collect() |>
      ungroup()
  } else {
    num_tech_reps_per_sample_and_peptide <- input_table |>
      left_join( samples_id_tbl, by=join_by( {{input_table_sample_id_column}} == {{sample_id_tbl_sample_id_column}} ) ) |>
      group_by( {{replicate_group_column}}, {{protein_id_column}}, {{peptide_sequence_column}}) |>
      partition(core_utilisation) |>
      summarise(counts = n() ) |>
      collect() |>
      ungroup()
  }

  # Any peptides found in more than one replicates in any patient will be kept for analysis
  removed_peptides_with_only_one_replicate <- input_table |>
    inner_join( num_tech_reps_per_sample_and_peptide |>
                  dplyr::filter( counts >  1) |>
                  dplyr::select(-counts, -{{replicate_group_column}}) |>
                  distinct()
                , by=join_by( {{protein_id_column}},
                              {{peptide_sequence_column}}) )  |>
    distinct()

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

  num_peptides_per_protein <- NA
  if( length(which(is.na(core_utilisation))) == 0 ) {
    num_peptides_per_protein <- input_table |>
      group_by( {{protein_id_column}} ) |>
      dplyr::summarise( peptides_for_protein_count = n()
                 , peptidoforms_for_protein_count = sum( peptidoform_count, na.rm=TRUE)) |>
      ungroup()
  } else {
    num_peptides_per_protein <- input_table |>
      group_by( {{protein_id_column}} ) |>
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
                  , by = join_by({{protein_id_column}})) |>
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

  samples_passing_filter <- NA
  if( length(which(is.na(core_utilisation))) == 0 ) {
    samples_passing_filter <- input_table |>
      group_by( {{sample_id_column}} ) |>
      #partition(core_utilisation) |>
      summarise( counts = n()) |>
      #collect() |>
      ungroup() |>
      dplyr::filter( counts >= peptides_per_sample_cutoff |
                       ( {{sample_id_column}} %in% inclusion_list)  ) |>
      dplyr::select(-counts)

  } else {
    samples_passing_filter <- input_table |>
      group_by( {{sample_id_column}} ) |>
      partition(core_utilisation) |>
      summarise( counts = n()) |>
      collect() |>
      ungroup() |>
      dplyr::filter( counts >= peptides_per_sample_cutoff |
                       ( {{sample_id_column}} %in% inclusion_list)  ) |>
      dplyr::select(-counts)
  }

  filtered_table <- input_table |>
    inner_join( samples_passing_filter, by = join_by({{sample_id_column}}))

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
  q_val_name <- rlang::as_name(rlang::ensym(q_value_column))
  global_q_val_name <- rlang::as_name(rlang::ensym(global_q_value_column))
  proteotypic_name <- rlang::as_name(rlang::ensym(proteotypic_peptide_sequence_column))
  
  filter_cols <- c(q_val_name, global_q_val_name, proteotypic_name)
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

  search_srl_quant_cln <- input_table |>
    dplyr::filter( {{q_value_column}} < qvalue_threshold &
                     {{global_q_value_column}} < global_qvalue_threshold &
                     {{proteotypic_peptide_sequence_column}} == choose_only_proteotypic_peptide ) |>
    dplyr::select(all_of(unique(c(input_matrix_column_ids, filter_cols))))

  search_srl_quant_cln

}

