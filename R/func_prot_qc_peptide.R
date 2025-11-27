# ============================================================================
# func_prot_qc_peptide.R
# ============================================================================
# Purpose: Peptide-level quality control and filtering functions
# 
# This file contains functions for peptide-level QC filtering, including
# intensity filtering, missing value filtering, replicate filtering, and
# q-value filtering. Functions in this file are used by mod_prot_qc_peptide.R
# and related QC modules.
#
# NOTE: Peptides are part of the proteomics workflow, hence "prot" prefix.
# This file contains peptide-specific QC functions within the proteomics context.
#
# Functions to extract here:
# - peptideIntensityFiltering(): S4 method for peptide intensity filtering
# - peptideIntensityFilteringHelper(): Helper for peptide intensity filtering
# - removePeptidesWithMissingValuesPercent(): S4 method for missing value filtering
# - removePeptidesWithMissingValuesPercentHelper(): Helper for missing value filtering
# - removePeptidesWithOnlyOneReplicate(): S4 method for replicate filtering
# - removePeptidesWithOnlyOneReplicateHelper(): Helper for replicate filtering
# - filterMinNumPeptidesPerProtein(): S4 method for filtering by peptides per protein
# - filterMinNumPeptidesPerSample(): S4 method for filtering by peptides per sample
# - srlQvalueProteotypicPeptideClean(): S4 method for q-value filtering
# - Additional peptide QC helper functions
#
# Dependencies:
# - dplyr, tidyr
# - func_general_helpers.R (for utility functions)
# ============================================================================

# TODO: Extract the following functions from their current locations:

# Function 1: peptideIntensityFiltering()
# Current location: R/peptideVsSamplesS4Objects.R
# Type: S4 method (exportMethods)
# Description: Filters peptides based on intensity thresholds
# setMethod(f = "peptideIntensityFiltering", ...) {
#   # Extract from R/peptideVsSamplesS4Objects.R
# }

# Function 2: peptideIntensityFilteringHelper()
# Current location: R/peptideVsSamplesS4Objects.R
# Description: Helper function for peptide intensity filtering
# peptideIntensityFilteringHelper <- function(...) {
#   # Extract from R/peptideVsSamplesS4Objects.R
# }

# Function 3: removePeptidesWithMissingValuesPercent()
# Current location: R/peptideVsSamplesS4Objects.R
# Type: S4 method (exportMethods)
# Description: Removes peptides with excessive missing values
# setMethod(f = "removePeptidesWithMissingValuesPercent", ...) {
#   # Extract from R/peptideVsSamplesS4Objects.R
# }

# Function 4: removePeptidesWithMissingValuesPercentHelper()
# Current location: R/peptideVsSamplesS4Objects.R
# Description: Helper for removing peptides with missing values
# removePeptidesWithMissingValuesPercentHelper <- function(...) {
#   # Extract from R/peptideVsSamplesS4Objects.R
# }

# Function 5: removePeptidesWithOnlyOneReplicate()
# Current location: R/peptideVsSamplesS4Objects.R
# Type: S4 method (exportMethods)
# Description: Removes peptides present in only one replicate
# setMethod(f = "removePeptidesWithOnlyOneReplicate", ...) {
#   # Extract from R/peptideVsSamplesS4Objects.R
# }

# Function 6: removePeptidesWithOnlyOneReplicateHelper()
# Current location: R/peptideVsSamplesS4Objects.R
# Description: Helper for removing peptides with only one replicate
# removePeptidesWithOnlyOneReplicateHelper <- function(...) {
#   # Extract from R/peptideVsSamplesS4Objects.R
# }

# Function 7: filterMinNumPeptidesPerProtein()
# Current location: R/peptideVsSamplesS4Objects.R
# Type: S4 method (exportMethods)
# Description: Filters proteins with insufficient peptides
# setMethod(f = "filterMinNumPeptidesPerProtein", ...) {
#   # Extract from R/peptideVsSamplesS4Objects.R
# }

# Function 8: filterMinNumPeptidesPerProteinHelper()
# Current location: R/peptideVsSamplesS4Objects.R
# Description: Helper for filtering by peptides per protein
# filterMinNumPeptidesPerProteinHelper <- function(...) {
#   # Extract from R/peptideVsSamplesS4Objects.R
# }

# Function 9: filterMinNumPeptidesPerSample()
# Current location: R/peptideVsSamplesS4Objects.R
# Type: S4 method (exportMethods)
# Description: Filters samples with insufficient peptides
# setMethod(f = "filterMinNumPeptidesPerSample", ...) {
#   # Extract from R/peptideVsSamplesS4Objects.R
# }

# Function 10: filterMinNumPeptidesPerSampleHelper()
# Current location: R/peptideVsSamplesS4Objects.R
# Description: Helper for filtering by peptides per sample
# filterMinNumPeptidesPerSampleHelper <- function(...) {
#   # Extract from R/peptideVsSamplesS4Objects.R
# }

# Function 11: srlQvalueProteotypicPeptideClean()
# Current location: R/peptideVsSamplesS4Objects.R
# Type: S4 method (exportMethods)
# Description: Filters peptides based on q-value and proteotypic status
# setMethod(f = "srlQvalueProteotypicPeptideClean", ...) {
#   # Extract from R/peptideVsSamplesS4Objects.R
# }

# Function 12: srlQvalueProteotypicPeptideCleanHelper()
# Current location: R/peptideVsSamplesS4Objects.R
# Description: Helper for q-value and proteotypic filtering
# srlQvalueProteotypicPeptideCleanHelper <- function(...) {
#   # Extract from R/peptideVsSamplesS4Objects.R
# }

# Function 13: checkPeptideNAPercentages()
# Current location: R/de_proteins_functions.R
# Description: Checks peptide NA percentages
# checkPeptideNAPercentages <- function(...) {
#   # Extract from R/de_proteins_functions.R
# }

# Function 14: removePeptidesOnlyInHek293()
# Current location: R/de_proteins_functions.R
# Description: Removes peptides only found in HEK293 cells
# removePeptidesOnlyInHek293 <- function(...) {
#   # Extract from R/de_proteins_functions.R
# }

# Function 15: removePeptidesWithoutAbundances()
# Current location: R/de_proteins_functions.R
# Description: Removes peptides without abundance values
# removePeptidesWithoutAbundances <- function(...) {
#   # Extract from R/de_proteins_functions.R
# }

# Function 16: filterByScoreAndGetSimilarPeptides()
# Current location: R/enrichment_functions.R
# Description: Filters peptides by score and finds similar peptides
# filterByScoreAndGetSimilarPeptides <- function(...) {
#   # Extract from R/enrichment_functions.R
# }

# Function 17: filterPeptideAndExtractProbabilities()
# Current location: R/enrichment_functions.R
# Description: Filters peptides and extracts probabilities
# filterPeptideAndExtractProbabilities <- function(...) {
#   # Extract from R/enrichment_functions.R
# }

# Function 18: groupParalogPeptides()
# Current location: R/enrichment_functions.R
# Description: Groups paralog peptides
# groupParalogPeptides <- function(...) {
#   # Extract from R/enrichment_functions.R
# }


# ----------------------------------------------------------------------------
# peptideIntensityFilteringHelper
# ----------------------------------------------------------------------------
#' @export
#' @title Filter Peptides by Intensity and Proportion
#' @description Remove peptide based on the intensity threshold and the proportion of samples below the threshold
peptideIntensityFilteringHelper <- function(input_table
                                      , min_peptide_intensity_threshold = 15
                                      , peptides_proportion_of_samples_below_cutoff = 1
                                      , protein_id_column = Protein.Ids
                                      , peptide_sequence_column = Stripped.Sequence
                                      , peptide_quantity_column = Peptide.Normalised
                                      , core_utilisation) {
  print(">>> Entering peptideIntensityFilteringHelper <<<")
  
  print(sprintf("      peptideIntensityFilteringHelper Arg: min_peptide_intensity_threshold = %g", min_peptide_intensity_threshold))
  print(sprintf("      peptideIntensityFilteringHelper Arg: peptides_proportion_of_samples_below_cutoff = %g", peptides_proportion_of_samples_below_cutoff))
  print(sprintf("      peptideIntensityFilteringHelper Arg: protein_id_column = %s", deparse(substitute(protein_id_column))))
  print(sprintf("      peptideIntensityFilteringHelper Arg: peptide_sequence_column = %s", deparse(substitute(peptide_sequence_column))))
  print(sprintf("      peptideIntensityFilteringHelper Arg: peptide_quantity_column = %s", deparse(substitute(peptide_quantity_column))))
  
  print(sprintf("      Data State (input_table): Dims = %d rows, %d cols", nrow(input_table), ncol(input_table)))
  print("      Data State (input_table) Structure:")
  utils::str(input_table)
  print("      Data State (input_table) Head:")
  print(head(input_table))

  num_values_per_peptide <- NA

  print("      peptideIntensityFilteringHelper Step: Checking core utilisation...")
  print(sprintf("      peptideIntensityFilteringHelper: is.na(core_utilisation) = %s", is.na(core_utilisation)))

  if( length(which(is.na(core_utilisation))) == 0 ) {
    print("      peptideIntensityFilteringHelper Step: Processing WITHOUT parallelization...")
    num_values_per_peptide <- input_table |>
      mutate(  below_intensity_threshold = case_when( {{peptide_quantity_column}} < min_peptide_intensity_threshold ~ 1,
                                                      TRUE ~ 0) ) |>
      group_by( {{protein_id_column}}, {{peptide_sequence_column}}) |>
      #partition(core_utilisation) |>
      summarise (samples_counts = n(),
                 num_below_intesnity_treshold = sum(below_intensity_threshold)) |>
      #collect() |>
      ungroup() |>
      dplyr::filter( num_below_intesnity_treshold/samples_counts < peptides_proportion_of_samples_below_cutoff )
  } else {
    print("      peptideIntensityFilteringHelper Step: Processing WITH parallelization...")
    num_values_per_peptide <- input_table |>
      mutate(  below_intensity_threshold = case_when( {{peptide_quantity_column}} < min_peptide_intensity_threshold ~ 1,
                                                      TRUE ~ 0) ) |>
      group_by( {{protein_id_column}}, {{peptide_sequence_column}}) |>
      partition(core_utilisation) |>
      summarise (samples_counts = n(),
                 num_below_intesnity_treshold = sum(below_intensity_threshold)) |>
      collect() |>
      ungroup() |>
      dplyr::filter( num_below_intesnity_treshold/samples_counts < peptides_proportion_of_samples_below_cutoff )

  }

  print("      peptideIntensityFilteringHelper Step: Filtering calculations completed.")
  print(sprintf("      Data State (num_values_per_peptide): Dims = %d rows, %d cols", nrow(num_values_per_peptide), ncol(num_values_per_peptide)))
  print("      Data State (num_values_per_peptide) Head:")
  print(head(num_values_per_peptide, 10))

  print("      peptideIntensityFilteringHelper Step: Applying inner_join to keep only passing peptides...")
  peptide_normalised_pif_cln <- input_table |>
    inner_join ( num_values_per_peptide |>
                   dplyr::select( -num_below_intesnity_treshold, -samples_counts)
                 , by = join_by( {{protein_id_column}}, {{peptide_sequence_column}} ) )

  print(sprintf("      peptideIntensityFilteringHelper: Original table had %d rows", nrow(input_table)))
  print(sprintf("      peptideIntensityFilteringHelper: Filtered table has %d rows", nrow(peptide_normalised_pif_cln)))
  print(sprintf("      peptideIntensityFilteringHelper: ACTUALLY REMOVED: %d rows", nrow(input_table) - nrow(peptide_normalised_pif_cln)))

  # Count proteins before and after
  proteins_before <- input_table |> dplyr::distinct({{protein_id_column}}) |> nrow()
  proteins_after <- peptide_normalised_pif_cln |> dplyr::distinct({{protein_id_column}}) |> nrow()
  print(sprintf("      peptideIntensityFilteringHelper: Proteins before: %d, after: %d, removed: %d", proteins_before, proteins_after, proteins_before - proteins_after))

  print("<<< Exiting peptideIntensityFilteringHelper <<<")
  return(peptide_normalised_pif_cln)

}


# ----------------------------------------------------------------------------
# removePeptidesWithMissingValuesPercentHelper
# ----------------------------------------------------------------------------
#' @title Remove Peptides with Missing Values
#'@param input_table An input table with a column containing the row ID and the rest of the columns representing abundance values for each sample.
#'@param cols A tidyselect command to select the columns. This includes the functions starts_with(), ends_with(), contains(), matches(), and num_range()
#'@param design_matrix A data frame with a column containing the sample ID (as per the sample_id param) and the experimental group (as per the group param). Each row as the sample ID as row name in the data frame.
#'@param sample_id The name of the column in design_matrix table that has the sample ID.
#'@param row_id A unique ID for each row of the 'input_table' variable.
#'@param group_column The name of the column in design_matrix table that has the experimental group.
#'@param groupwise_percentage_cutoff The maximum percentage of values below threshold allow in each group for a peptide .
#'@param max_groups_percentage_cutoff The maximum percentage of groups allowed with too many samples with peptide abundance values below threshold.
#'@param abundance_threshold Abundance threshold in which the protein in the sample must be above for it to be considered for inclusion into data analysis.
#'@param temporary_abundance_column The name of a temporary column to keep the abundance value you want to filter upon
#'@return A list, the name of each element is the sample ID and each element is a vector containing the protein accessions (e.g. row_id) with enough number of values.
#'@export
removePeptidesWithMissingValuesPercentHelper <- function(input_table
                                               , design_matrix
                                               , sample_id
                                               , protein_id_column
                                               , peptide_sequence_column
                                               , grouping_variable
                                               , groupwise_percentage_cutoff = 1
                                               , max_groups_percentage_cutoff = 50
                                               , abundance_threshold
                                               , abundance_column = "Abundance") {

  abundance_long <- input_table |>
    mutate( row_id = purrr::map2_chr( {{protein_id_column}}
                                     , {{peptide_sequence_column}}
                                     , \(x,y)paste(x , y, sep="_")) ) |>
    mutate( {{sample_id}} := purrr::map_chr(   {{sample_id}}  , as.character)   ) |>
    left_join(  design_matrix |>
                mutate(  {{sample_id}} := purrr::map_chr( {{sample_id}} , as.character ))
                , by = join_by({{sample_id}} ) )

  count_values_per_group <- abundance_long |>
    distinct( {{sample_id}} , {{ grouping_variable }} ) |>
    group_by( {{ grouping_variable }} ) |>
    summarise(  num_per_group = n()) |>
    ungroup()

  count_values_missing_per_group <- abundance_long |>
    mutate(is_missing = ifelse( !is.na( !!sym( abundance_column ))
                                & !!sym( abundance_column ) > abundance_threshold
                                , 0, 1)) |>
    group_by( row_id, {{ grouping_variable }} ) |>
    summarise( num_missing_per_group = sum(is_missing)) |>
    ungroup()

  count_percent_missing_per_group <- count_values_missing_per_group |>
    full_join( count_values_per_group,
               by = join_by( {{ grouping_variable }} )) |>
    mutate(  perc_below_thresh_per_group = num_missing_per_group / num_per_group * 100 )

  total_num_of_groups <- count_values_per_group |> nrow()

  remove_rows_temp <- count_percent_missing_per_group |>
    dplyr::filter(groupwise_percentage_cutoff <  perc_below_thresh_per_group) |>
    group_by( row_id ) |>
    summarise( percent  = n()/total_num_of_groups*100 ) |>
    ungroup() |>
    dplyr::filter(percent > max_groups_percentage_cutoff)

  print(nrow(remove_rows_temp))

  filtered_tbl <- input_table |>
    mutate( row_id = purrr::map2_chr( {{protein_id_column}}
                                     , {{peptide_sequence_column}}
                                     , \(x,y)paste(x , y, sep="_")) ) |>
    dplyr::anti_join(remove_rows_temp, by = join_by(row_id)) |>
    dplyr::select(-row_id)

  return(filtered_tbl)

}


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

  # ✅ DIAGNOSTIC + DEFENSIVE: Check output column availability
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
  
  # ✅ ALSO CHECK: Filter columns exist
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
    dplyr::select(all_of( input_matrix_column_ids))

  search_srl_quant_cln

}


# ----------------------------------------------------------------------------
# checkPeptideNAPercentages
# ----------------------------------------------------------------------------
#' Check Missing Value Percentages in Peptide Data
#' 
#' @description Calculate and report the percentage of missing values (NAs) in peptide data
#' at different levels: total dataset, per sample, and per group.
#' 
#' @param peptide_obj A PeptideQuantitativeData S4 object
#' @param verbose Logical, whether to print detailed results (default: TRUE)
#' 
#' @return A list containing:
#' \itemize{
#'   \item total_na_percent: Overall percentage of NAs in the dataset
#'   \item per_sample_na: Data frame with NA percentages per sample
#'   \item per_group_na: Data frame with NA percentages per group
#'   \item summary_stats: Summary statistics of NA distribution
#' }
#' 
#' @export
checkPeptideNAPercentages <- function(peptide_obj, verbose = TRUE) {
  
  # Validate input
  if (!is(peptide_obj, "PeptideQuantitativeData")) {
    stop("Input must be a PeptideQuantitativeData S4 object")
  }
  
  # Extract data from S4 object
  peptide_matrix <- peptide_obj@peptide_matrix
  design_matrix <- peptide_obj@design_matrix
  sample_id_col <- peptide_obj@sample_id
  group_id_col <- peptide_obj@group_id
  
  # Validate that matrix and design matrix are compatible
  if (ncol(peptide_matrix) != nrow(design_matrix)) {
    stop("Number of samples in peptide_matrix doesn't match design_matrix rows")
  }
  
  # Calculate total NA percentage
  total_values <- length(peptide_matrix)
  total_nas <- sum(is.na(peptide_matrix))
  total_na_percent <- (total_nas / total_values) * 100
  
  # Calculate per-sample NA percentages
  sample_na_counts <- apply(peptide_matrix, 2, function(x) sum(is.na(x)))
  sample_na_percentages <- (sample_na_counts / nrow(peptide_matrix)) * 100
  
  per_sample_na <- data.frame(
    sample = colnames(peptide_matrix),
    na_count = sample_na_counts,
    na_percentage = sample_na_percentages,
    stringsAsFactors = FALSE
  )
  
  # Add group information to per-sample results
  per_sample_na <- merge(per_sample_na, design_matrix, 
                        by.x = "sample", by.y = sample_id_col, all.x = TRUE)
  
  # Calculate per-group NA percentages
  per_group_na <- per_sample_na %>%
    group_by(!!sym(group_id_col)) %>%
    summarise(
      num_samples = n(),
      mean_na_percentage = mean(na_percentage, na.rm = TRUE),
      median_na_percentage = median(na_percentage, na.rm = TRUE),
      min_na_percentage = min(na_percentage, na.rm = TRUE),
      max_na_percentage = max(na_percentage, na.rm = TRUE),
      sd_na_percentage = sd(na_percentage, na.rm = TRUE),
      .groups = 'drop'
    ) %>%
    arrange(mean_na_percentage)
  
  # Calculate summary statistics
  summary_stats <- list(
    total_peptides = nrow(peptide_matrix),
    total_samples = ncol(peptide_matrix),
    total_groups = length(unique(design_matrix[[group_id_col]])),
    total_values = total_values,
    total_nas = total_nas,
    mean_na_per_sample = mean(sample_na_percentages),
    median_na_per_sample = median(sample_na_percentages),
    min_na_per_sample = min(sample_na_percentages),
    max_na_per_sample = max(sample_na_percentages)
  )
  
  # Print results if verbose
  if (verbose) {
    cat("\n=== Peptide Data Missing Value Analysis ===\n")
    cat(sprintf("Dataset dimensions: %d peptides × %d samples\n", 
                nrow(peptide_matrix), ncol(peptide_matrix)))
    cat(sprintf("Number of groups: %d\n", summary_stats$total_groups))
    cat(sprintf("Total missing values: %s out of %s (%.2f%%)\n", 
                format(total_nas, big.mark = ","),
                format(total_values, big.mark = ","),
                total_na_percent))
    
    cat("\n--- Per-Sample Missing Value Summary ---\n")
    cat(sprintf("Mean NA%% per sample: %.2f%%\n", summary_stats$mean_na_per_sample))
    cat(sprintf("Median NA%% per sample: %.2f%%\n", summary_stats$median_na_per_sample))
    cat(sprintf("Range: %.2f%% - %.2f%%\n", 
                summary_stats$min_na_per_sample, summary_stats$max_na_per_sample))
    
    cat("\n--- Per-Group Missing Value Summary ---\n")
    print(per_group_na)
    
    cat("\n--- Samples with Highest Missing Values ---\n")
    top_missing_samples <- per_sample_na %>%
      arrange(desc(na_percentage)) %>%
      head(min(5, nrow(per_sample_na)))
    print(top_missing_samples[, c("sample", group_id_col, "na_percentage")])
    
    cat("\n--- Samples with Lowest Missing Values ---\n")
    bottom_missing_samples <- per_sample_na %>%
      arrange(na_percentage) %>%
      head(min(5, nrow(per_sample_na)))
    print(bottom_missing_samples[, c("sample", group_id_col, "na_percentage")])
  }
  
  # Return results
  results <- list(
    total_na_percent = total_na_percent,
    per_sample_na = per_sample_na,
    per_group_na = per_group_na,
    summary_stats = summary_stats
  )
  
  return(invisible(results))
}


# ----------------------------------------------------------------------------
# removePeptidesOnlyInHek293
# ----------------------------------------------------------------------------
#---------------------------------------------------------------------------------------
#' @export
removePeptidesOnlyInHek293 <- function( input_table
                                        , metadata_table
                                        , input_table_sample_id_column = "Run"
                                        , sample_id_tbl_sample_id_column  =  "ms_filename"
                                        , protein_id_column = Protein.Ids
                                        , peptide_sequence_column = Stripped.Sequence
                                        , hek_string = "HEK"
                                        , general_sample_info = general_sample_info
                                        , core_utilisation= core_utilisation) {


  peptides_found_in_hek_samples_only <- NA
  if( length(which(is.na(core_utilisation))) == 0 ) {

    peptides_found_in_hek_samples_only <- input_table |>
      left_join( metadata_table |>
                   dplyr::distinct( {{sample_id_tbl_sample_id_column}}, {{general_sample_info}})
                 , by=join_by({{input_table_sample_id_column}} == {{sample_id_tbl_sample_id_column}}) ) |>
      mutate( sample_type = case_when ( str_detect( {{general_sample_info}}, hek_string) ~ hek_string,
                                        TRUE ~ "Cohort_Sample")) |>
      group_by( {{protein_id_column}}, {{peptide_sequence_column}}, sample_type ) |>
      partition(core_utilisation) |>
      summarise( counts = n()) |>
      collect() |>
      ungroup() |>
      pivot_wider ( names_from = sample_type
                    , values_from = counts ) |>
      dplyr::filter( !is.na( !!sym(hek_string)) & is.na( Cohort_Sample) ) |>
      dplyr::distinct( {{protein_id_column}}, {{peptide_sequence_column}} )

  } else {
    peptides_found_in_hek_samples_only <- input_table |>
      left_join( metadata_table |>
                   dplyr::distinct( {{sample_id_tbl_sample_id_column}}, {{general_sample_info}})
                 , by=join_by({{input_table_sample_id_column}} == {{sample_id_tbl_sample_id_column}}) ) |>
      mutate( sample_type = case_when ( str_detect( {{general_sample_info}}, hek_string) ~ hek_string,
                                        TRUE ~ "Cohort_Sample")) |>
      group_by( {{protein_id_column}}, {{peptide_sequence_column}}, sample_type ) |>
      #partition(core_utilisation) |>
      summarise( counts = n()) |>
      #collect() |>
      ungroup() |>
      pivot_wider ( names_from = sample_type
                    , values_from = counts ) |>
      dplyr::filter( !is.na( !!sym(hek_string)) & is.na( Cohort_Sample) ) |>
      dplyr::distinct( {{protein_id_column}}, {{peptide_sequence_column}} )

  }

  removed_peptides_only_in_hek_samples <- input_table |>
    anti_join( peptides_found_in_hek_samples_only
               , by = join_by(  {{protein_id_column}}, {{peptide_sequence_column}} ) )

}


# ----------------------------------------------------------------------------
# compareTwoPeptideDataObjects
# ----------------------------------------------------------------------------
# I want to input two peptide data objects and compare them,
# to see how the number of proteins and peptides changes and how the number of samples changed
# Use set diff or set intersect to compare the peptides, proteins, samples in the two objects
#'@export
compareTwoPeptideDataObjects <- function( object_a, object_b) {

  object_a_peptides <- object_a@peptide_data |>
    distinct(!!sym(object_a@protein_id_column), !!sym(object_a@peptide_sequence_column))

  object_b_peptides <- object_b@peptide_data |>
    distinct(!!sym(object_b@protein_id_column), !!sym(object_b@peptide_sequence_column))

  object_a_proteins <- object_a@peptide_data |>
    distinct(!!sym(object_a@protein_id_column)) |>
    dplyr::pull(!!sym(object_a@protein_id_column))

  object_b_proteins <- object_b@peptide_data |>
    distinct(!!sym(object_b@protein_id_column)) |>
    dplyr::pull(!!sym(object_b@protein_id_column))

  object_a_samples <- object_a@design_matrix |>
    distinct(!!sym(object_a@sample_id)) |>
    dplyr::pull(!!sym(object_a@sample_id))

  object_b_samples <- object_b@design_matrix |>
    distinct(!!sym(object_b@sample_id)) |>
    dplyr::pull(!!sym(object_b@sample_id))


  peptides_in_a_not_b <- nrow( dplyr::setdiff( object_a_peptides, object_b_peptides) )
  peptides_intersect_a_and_b <- nrow( dplyr::intersect( object_a_peptides, object_b_peptides) )
  peptides_in_b_not_a <- nrow(  dplyr::setdiff( object_b_peptides, object_a_peptides) )

  proteins_in_a_not_b <- length( setdiff( object_a_proteins, object_b_proteins) )
  proteins_intersect_a_and_b <- length( intersect( object_a_proteins, object_b_proteins) )
  proteins_in_b_not_a <- length( setdiff( object_b_proteins, object_a_proteins) )


  samples_in_a_not_b <- length( setdiff( object_a_samples, object_b_samples) )
  samples_intersect_a_and_b <- length( intersect( object_a_samples, object_b_samples) )
  samples_in_b_not_a <- length( setdiff( object_b_samples, object_a_samples) )

  comparisons_list <- list(  peptides = list( in_a_not_b = peptides_in_a_not_b
                                             , intersect_a_and_b = peptides_intersect_a_and_b
                                             , in_b_not_a = peptides_in_b_not_a)
                            , proteins = list( in_a_not_b = proteins_in_a_not_b
                                               , intersect_a_and_b = proteins_intersect_a_and_b
                                               , in_b_not_a = proteins_in_b_not_a)
                            , samples = list( in_a_not_b = samples_in_a_not_b
                                              , intersect_a_and_b = samples_intersect_a_and_b
                                              , in_b_not_a = samples_in_b_not_a)
  )

  comparison_tibble <- comparisons_list |>
    purrr::map_df( tibble::as_tibble) |>
    add_column( Levels = c("peptides", "proteins", "samples")) |>
    relocate( Levels, .before="in_a_not_b")

  comparison_tibble


}


# ----------------------------------------------------------------------------
# summarisePeptideObject
# ----------------------------------------------------------------------------
#'@export
summarisePeptideObject <- function(theObject) {

  num_peptides <- theObject@peptide_data |>
    distinct(!!sym(theObject@protein_id_column), !!sym(theObject@peptide_sequence_column))

  num_proteins <- theObject@peptide_data |>
    distinct(!!sym(theObject@protein_id_column)) |>
    dplyr::pull(!!sym(theObject@protein_id_column))

  num_samples <- theObject@design_matrix |>
    distinct(!!sym(theObject@sample_id)) |>
    dplyr::pull(!!sym(theObject@sample_id))

  summary_list <- list( num_peptides = nrow(num_peptides)
                       , num_proteins = length(num_proteins)
                       , num_samples = length(num_samples))

  summary_list


}


# ----------------------------------------------------------------------------
# calculatePeptidePearsonCorrelation
# ----------------------------------------------------------------------------
calculatePeptidePearsonCorrelation <- function(temp_obj, tech_rep_remove_regex, correlation_group) {
  data_table <- temp_obj$data_table
  id_column <- temp_obj$id_column
  design_matrix <- temp_obj$design_matrix
  sample_id <- temp_obj$sample_id
  
  # Get sample columns (exclude ID column)
  sample_columns <- setdiff(colnames(data_table), id_column)
  
  # Filter out technical replicates if regex provided
  if(!is.null(tech_rep_remove_regex) && tech_rep_remove_regex != "") {
    sample_columns <- sample_columns[!grepl(tech_rep_remove_regex, sample_columns)]
  }
  
  # Create correlation matrix
  peptide_matrix_for_corr <- data_table |>
    column_to_rownames(id_column) |>
    select(all_of(sample_columns)) |>
    as.matrix()
  
  # Calculate correlations between all sample pairs
  sample_correlations <- cor(peptide_matrix_for_corr, use = "pairwise.complete.obs")
  
  # Extract upper triangle (avoid duplicate pairs and self-correlations)
  upper_tri_indices <- which(upper.tri(sample_correlations), arr.ind = TRUE)
  
  correlation_results <- data.frame(
    sample1 = rownames(sample_correlations)[upper_tri_indices[,1]],
    sample2 = colnames(sample_correlations)[upper_tri_indices[,2]],
    pearson_correlation = sample_correlations[upper_tri_indices],
    stringsAsFactors = FALSE
  )
  
  # Remove NA correlations
  correlation_results <- correlation_results[!is.na(correlation_results$pearson_correlation), ]
  
  return(correlation_results)
}


# ----------------------------------------------------------------------------
# peptideIntensityFiltering
# ----------------------------------------------------------------------------
#'@export
setMethod( f="peptideIntensityFiltering"
           , signature="PeptideQuantitativeData"
           , definition = function( theObject, peptides_intensity_cutoff_percentile = NULL, peptides_proportion_of_samples_below_cutoff = NULL, core_utilisation = NULL) {
             print("--- Entering peptideIntensityFiltering S4 Method ---")
             
             peptide_data <- theObject@peptide_data
             raw_quantity_column <- theObject@raw_quantity_column
             norm_quantity_column <- theObject@norm_quantity_column

             print("   peptideIntensityFiltering: Extracting input data...")
             print(sprintf("      Arg: raw_quantity_column = %s", raw_quantity_column))
             print(sprintf("      Arg: norm_quantity_column = %s", norm_quantity_column))
             print(sprintf("      Data State (peptide_data): Dims = %d rows, %d cols", nrow(peptide_data), ncol(peptide_data)))

             print("   peptideIntensityFiltering: Resolving parameters with checkParamsObjectFunctionSimplify...")
             peptides_intensity_cutoff_percentile <- checkParamsObjectFunctionSimplify( theObject
                                                                                    , "peptides_intensity_cutoff_percentile")
             print(sprintf("      Resolved peptides_intensity_cutoff_percentile = %g", peptides_intensity_cutoff_percentile))

             peptides_proportion_of_samples_below_cutoff <- checkParamsObjectFunctionSimplify( theObject
                                                                                                , "peptides_proportion_of_samples_below_cutoff")
             print(sprintf("      Resolved peptides_proportion_of_samples_below_cutoff = %g", peptides_proportion_of_samples_below_cutoff))

             core_utilisation <- checkParamsObjectFunctionSimplify( theObject, "core_utilisation", NA)
             print(sprintf("      Resolved core_utilisation = %s", ifelse(is.na(core_utilisation), "NA", as.character(core_utilisation))))

             print("   peptideIntensityFiltering: Updating parameters in S4 object...")
             theObject <- updateParamInObject(theObject, "peptides_intensity_cutoff_percentile")
             theObject <- updateParamInObject(theObject, "peptides_proportion_of_samples_below_cutoff")
             theObject <- updateParamInObject(theObject, "core_utilisation")

             print("   peptideIntensityFiltering: Calculating intensity threshold...")
             # Get non-missing values for threshold calculation
             valid_values <- peptide_data |> dplyr::pull(!!sym(raw_quantity_column))
             valid_values <- valid_values[!is.na(valid_values) & !is.nan(valid_values) & !is.infinite(valid_values)]
             
             print(sprintf("      peptideIntensityFiltering: Found %d valid intensity values", length(valid_values)))
             print(sprintf("      peptideIntensityFiltering: Valid values range: %g to %g", min(valid_values, na.rm=TRUE), max(valid_values, na.rm=TRUE)))
             
             min_peptide_intensity_threshold <- ceiling( quantile( peptide_data |> dplyr::pull(!!sym(raw_quantity_column)), na.rm=TRUE, probs = c(peptides_intensity_cutoff_percentile/100) ))[1]
             print(sprintf("      peptideIntensityFiltering: Calculated min_peptide_intensity_threshold = %g (percentile %g%%)", 
                          min_peptide_intensity_threshold, peptides_intensity_cutoff_percentile))

             print("   peptideIntensityFiltering: About to call helper function...")
             print(sprintf("      Helper Args: min_peptide_intensity_threshold = %g", min_peptide_intensity_threshold))
             print(sprintf("      Helper Args: peptides_proportion_of_samples_below_cutoff = %g", peptides_proportion_of_samples_below_cutoff))
             print(sprintf("      Helper Args: protein_id_column = %s", theObject@protein_id_column))
             print(sprintf("      Helper Args: peptide_sequence_column = %s", theObject@peptide_sequence_column))
             print(sprintf("      Helper Args: peptide_quantity_column = %s", raw_quantity_column))

             peptide_normalised_pif_cln <- peptideIntensityFilteringHelper( peptide_data
                                                                      , min_peptide_intensity_threshold = min_peptide_intensity_threshold
                                                                      , peptides_proportion_of_samples_below_cutoff = peptides_proportion_of_samples_below_cutoff
                                                                      , protein_id_column = !!sym( theObject@protein_id_column)
                                                                      , peptide_sequence_column = !!sym(theObject@peptide_sequence_column)
                                                                      , peptide_quantity_column = !!sym(raw_quantity_column)
                                                                      , core_utilisation = core_utilisation)

             print(sprintf("   peptideIntensityFiltering: Helper function returned. New dims = %d rows, %d cols", 
                          nrow(peptide_normalised_pif_cln), ncol(peptide_normalised_pif_cln)))

             theObject@peptide_data <- peptide_normalised_pif_cln

             print("   peptideIntensityFiltering: Cleaning design matrix...")
             theObject <- cleanDesignMatrixPeptide(theObject)

             print("--- Exiting peptideIntensityFiltering S4 Method ---")
             return(theObject)
           })


# ----------------------------------------------------------------------------
# removePeptidesWithMissingValuesPercent
# ----------------------------------------------------------------------------
#'@export
setMethod( f = "removePeptidesWithMissingValuesPercent"
           , signature="PeptideQuantitativeData"
           , definition=function( theObject
                                  , grouping_variable = NULL
                                  , groupwise_percentage_cutoff = NULL
                                  , max_groups_percentage_cutoff = NULL
                                  , peptides_intensity_cutoff_percentile = NULL) {

             peptide_data <- theObject@peptide_data
             protein_id_column <- theObject@protein_id_column
             peptide_sequence_column <- theObject@peptide_sequence_column
             raw_quantity_column <- theObject@raw_quantity_column
             norm_quantity_column <- theObject@norm_quantity_column
             sample_id <- theObject@sample_id

             design_matrix <- theObject@design_matrix

             grouping_variable <- checkParamsObjectFunctionSimplify( theObject
                                                                   , "grouping_variable"
                                                                   , NULL)
             groupwise_percentage_cutoff <- checkParamsObjectFunctionSimplify( theObject
                                                                                   , "groupwise_percentage_cutoff"
                                                                                   , 50)
             max_groups_percentage_cutoff <- checkParamsObjectFunctionSimplify( theObject
                                                                                   , "max_groups_percentage_cutoff"
                                                                                   , 50)
             peptides_intensity_cutoff_percentile <- checkParamsObjectFunctionSimplify( theObject
                                                                                    , "peptides_intensity_cutoff_percentile"
                                                                                    , 50)

             theObject <- updateParamInObject(theObject, "grouping_variable")
             theObject <- updateParamInObject(theObject, "groupwise_percentage_cutoff")
             theObject <- updateParamInObject(theObject, "max_groups_percentage_cutoff")
             theObject <- updateParamInObject(theObject, "peptides_intensity_cutoff_percentile")

             min_protein_intensity_threshold <- ceiling( quantile( peptide_data |>
                                                                     dplyr::filter( !is.nan(!!sym(norm_quantity_column)) & !is.infinite(!!sym(norm_quantity_column))) |>
                                                                     dplyr::pull(!!sym(norm_quantity_column))
                                                                   , na.rm=TRUE
                                                                   , probs = c(peptides_intensity_cutoff_percentile/100) ))[1]

             # print(min_protein_intensity_threshold )

             theObject@peptide_data <- removePeptidesWithMissingValuesPercentHelper( peptide_data
                                                                               , design_matrix = design_matrix
                                                                               , sample_id = !!sym(sample_id)
                                                                               , protein_id_column = !!sym(protein_id_column)
                                                                               , peptide_sequence_column = !!sym(peptide_sequence_column)
                                                                               , grouping_variable = !!sym(grouping_variable)
                                                                               , groupwise_percentage_cutoff = groupwise_percentage_cutoff
                                                                               , max_groups_percentage_cutoff = max_groups_percentage_cutoff
                                                                               , abundance_threshold = peptides_intensity_cutoff_percentile
                                                                               , abundance_column =  norm_quantity_column )


             theObject <- cleanDesignMatrixPeptide(theObject)

             return(theObject)

           })


# ----------------------------------------------------------------------------
# removePeptidesWithOnlyOneReplicate
# ----------------------------------------------------------------------------
#'@export
setMethod( f="removePeptidesWithOnlyOneReplicate"
           , signature="PeptideQuantitativeData"
           , definition = function( theObject, replicate_group_column = NULL, core_utilisation = NULL) {

             peptide_data <- theObject@peptide_data
             sample_id_column <- theObject@sample_id
             design_matrix <- theObject@design_matrix


             grouping_variable <- checkParamsObjectFunctionSimplifyAcceptNull( theObject
                                                                       , "replicate_group_column"
                                                                       , NULL)

             core_utilisation <- checkParamsObjectFunctionSimplify( theObject
                                                           , "core_utilisation"
                                                           , NA)

             theObject <- updateParamInObject(theObject, "replicate_group_column")
             theObject <- updateParamInObject(theObject, "core_utilisation")

             theObject@peptide_data <- removePeptidesWithOnlyOneReplicateHelper( input_table = peptide_data
                                                                                             , samples_id_tbl = design_matrix
                                                                                             , input_table_sample_id_column = !!sym(sample_id_column)
                                                                                             , sample_id_tbl_sample_id_column  = !!sym(sample_id_column)
                                                                                             , replicate_group_column = !!sym(replicate_group_column)
                                                                                             , core_utilisation = core_utilisation)
             theObject <- cleanDesignMatrixPeptide(theObject)

             theObject
           })


# ----------------------------------------------------------------------------
# filterMinNumPeptidesPerProtein
# ----------------------------------------------------------------------------
#'@title Filter the proteins based on the number of peptides and peptidoforms
#'@description Keep the proteins only if they have two or more peptides.
#'@param theObject Object of class PeptideQuantitativeData
#'@param num_peptides_per_protein_thresh Minimum number of peptides per protein
#'@param num_peptidoforms_per_protein_thresh Minimum number of peptidoforms per protein
#'@param core_utilisation core_utilisation to use for parallel processing
#'@export
setMethod( f="filterMinNumPeptidesPerProtein"
           , signature="PeptideQuantitativeData"
           , definition = function( theObject, ... ) {
             
             # Extract specific parameters from ...
             args <- list(...)
             num_peptides_per_protein_thresh <- args$num_peptides_per_protein_thresh
             num_peptidoforms_per_protein_thresh <- args$num_peptidoforms_per_protein_thresh
             core_utilisation <- args$core_utilisation
             peptide_data <- theObject@peptide_data
             protein_id_column <- theObject@protein_id_column

             num_peptides_per_protein_thresh <- checkParamsObjectFunctionSimplify( theObject
                                                                                   , "num_peptides_per_protein_thresh"
                                                                                   , 1)

             num_peptidoforms_per_protein_thresh <- checkParamsObjectFunctionSimplify( theObject
                                                                                       , "num_peptidoforms_per_protein_thresh"
                                                                                       , 2)

             core_utilisation <- checkParamsObjectFunctionSimplify( theObject, "core_utilisation", NA)


             theObject <- updateParamInObject(theObject, "num_peptides_per_protein_thresh")
             theObject <- updateParamInObject(theObject, "num_peptidoforms_per_protein_thresh")
             theObject <- updateParamInObject(theObject, "core_utilisation")

             theObject@peptide_data <- filterMinNumPeptidesPerProteinHelper ( input_table = peptide_data
                                                                        , num_peptides_per_protein_thresh = num_peptides_per_protein_thresh
                                                                        , num_peptidoforms_per_protein_thresh = num_peptidoforms_per_protein_thresh
                                                                        , protein_id_column = !!sym(protein_id_column)
                                                                        , core_utilisation = core_utilisation)

             theObject <- cleanDesignMatrixPeptide(theObject)

             theObject
           })


# ----------------------------------------------------------------------------
# filterMinNumPeptidesPerSample
# ----------------------------------------------------------------------------
#'@export
setMethod( f="filterMinNumPeptidesPerSample"
           , signature="PeptideQuantitativeData"
           , definition = function( theObject
                                    , peptides_per_sample_cutoff = NULL
                                    , core_utilisation = NULL
                                    , inclusion_list = NULL) {

             peptide_data <- theObject@peptide_data
             sample_id_column <- theObject@sample_id

             peptides_per_sample_cutoff <- checkParamsObjectFunctionSimplify( theObject
                                                                              , "peptides_per_sample_cutoff"
                                                                              , 5000)

             inclusion_list <- checkParamsObjectFunctionSimplifyAcceptNull( theObject
                                                                            , "inclusion_list"
                                                                            , NULL)

             core_utilisation <- checkParamsObjectFunctionSimplify( theObject, "core_utilisation", NA)

             theObject <- updateParamInObject(theObject, "peptides_per_sample_cutoff")
             theObject <- updateParamInObject(theObject, "inclusion_list")
             theObject <- updateParamInObject(theObject, "core_utilisation")

             theObject@peptide_data <- filterMinNumPeptidesPerSampleHelper( peptide_data
                                            , peptides_per_sample_cutoff = peptides_per_sample_cutoff
                                            , sample_id_column = !!sym(sample_id_column)
                                            , core_utilisation
                                            , inclusion_list = inclusion_list )

             theObject <- cleanDesignMatrixPeptide(theObject)

             theObject
           })


# ----------------------------------------------------------------------------
# srlQvalueProteotypicPeptideClean
# ----------------------------------------------------------------------------
#'@export
setMethod( f ="srlQvalueProteotypicPeptideClean"
           , signature="PeptideQuantitativeData"
           , definition=function ( theObject
                                  , qvalue_threshold = NULL
                                  , global_qvalue_threshold = NULL
                                  , choose_only_proteotypic_peptide = NULL
                                  , input_matrix_column_ids =  NULL
                                  ) {
             peptide_data <- theObject@peptide_data
             protein_id_column <- theObject@protein_id_column
             q_value_column <- theObject@q_value_column
             global_q_value_column <- theObject@global_q_value_column
             peptide_sequence_column <- theObject@peptide_sequence_column
             proteotypic_peptide_sequence_column <- theObject@proteotypic_peptide_sequence_column
             raw_quantity_column <- theObject@raw_quantity_column
             norm_quantity_column <- theObject@norm_quantity_column

             qvalue_threshold <- checkParamsObjectFunctionSimplify( theObject, "qvalue_threshold", 0.01)

             global_qvalue_threshold <- checkParamsObjectFunctionSimplify( theObject, "global_qvalue_threshold", 0.01)

             choose_only_proteotypic_peptide <- checkParamsObjectFunctionSimplify( theObject
                                                                                   , "choose_only_proteotypic_peptide"
                                                                                   , 1 )

             input_matrix_column_ids <- checkParamsObjectFunctionSimplify( theObject
                                                                           , "input_matrix_column_ids" )

             theObject <- updateParamInObject(theObject, "qvalue_threshold")
             theObject <- updateParamInObject(theObject, "global_qvalue_threshold")
             theObject <- updateParamInObject(theObject, "choose_only_proteotypic_peptide")

             dia_nn_default_columns <- c("Protein.Ids"
                                        , "Stripped.Sequence"
                                        , "Q.Value"
                                        , "Global.Q.Value"
                                        , "Precursor.Quantity"
                                        , "Precursor.Normalised")

             theObject <- updateParamInObject(theObject, "input_matrix_column_ids")

             # print( paste("qvalue_threshold: ", qvalue_threshold))
             search_srl_quant_cln <- srlQvalueProteotypicPeptideCleanHelper( input_table = peptide_data
                                                                       , input_matrix_column_ids = unique(c(input_matrix_column_ids
                                                                                                      , protein_id_column
                                                                                                      , peptide_sequence_column
                                                                                                      , peptide_sequence_column))
                                                                       , protein_id_column = !!sym(protein_id_column)
                                                                       , q_value_column = !!sym(q_value_column)
                                                                       , global_q_value_column = !!sym(global_q_value_column)
                                                                       , global_qvalue_threshold = global_qvalue_threshold
                                                                       , qvalue_threshold = qvalue_threshold
                                                                       , choose_only_proteotypic_peptide = choose_only_proteotypic_peptide)

             theObject@peptide_data <- search_srl_quant_cln

             theObject <- cleanDesignMatrixPeptide(theObject)

             return(theObject)
           })

