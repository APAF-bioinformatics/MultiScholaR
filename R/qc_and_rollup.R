#---------------------------------------------------------------------------------------
#' @title Save Computation Time Record
#' @description This function records the time elapsed for a specific computation step
#' into a pre-initialized data frame.
#'
#' @param compute_time_record A data frame with columns "step_name" and "time_elapsed"
#'   to store the records.
#' @param step_name A character string naming the computation step.
#' @param toc_output The output from `tictoc::toc()`, which contains the callback message
#'   with the elapsed time.
#'
#' @return The updated `compute_time_record` data frame.
#' @export
saveTimeRecord <- function ( compute_time_record, step_name, toc_output ) {

  if(length(which(compute_time_record[,"step_name"] == "")) == 0) {
    stop("saveTimeRecord: No more room on compute time record table.")
  }
  # find_first_empty_slot
  first_empty_slot <- which(compute_time_record[,"step_name"] == "")[1]

  print(first_empty_slot)

  compute_time_record[first_empty_slot, "step_name"] <- step_name
  compute_time_record[first_empty_slot, "time_elapsed"] <- toc_output$callback_msg

  return(compute_time_record)
}

#---------------------------------------------------------------------------------------
# Peptide intensity filtering helper functions


#' @title Count the Number of Samples
#' @description Counts the number of unique samples in a given data table.
#'
#' @param input_table The input data frame.
#' @param sample_id_column The unquoted column name that contains the sample identifiers.
#'
#' @return An integer representing the number of unique samples.
#' @export
count_num_samples <- function( input_table
                               , sample_id_column = Run) {
  num_samples <- input_table |>
    distinct( {{sample_id_column}}) |>
    count()

  num_samples[[1,1]]
}



#' @title Count the Number of Proteins
#' @description Counts the number of unique proteins in a given data table.
#'
#' @param input_table The input data frame.
#' @param protein_id_column The unquoted column name that contains the protein identifiers.
#'
#' @return An integer representing the number of unique proteins.
#' @export
count_num_proteins <- function( input_table
                                , protein_id_column = Protein.Ids) {
  num_proteins <- input_table |>
    distinct( {{protein_id_column}}) |>
    count()

  num_proteins[[1,1]]
}


#' @title Count the Number of Peptides
#' @description Counts the number of unique peptide-protein combinations in a given data table.
#'
#' @param input_table The input data frame.
#' @param protein_id_column The unquoted column name that contains the protein identifiers.
#' @param peptide_sequence_column The unquoted column name that contains the peptide sequences.
#'
#' @return An integer representing the number of unique peptides.
#' @export
count_num_peptides <- function( input_table
                                , protein_id_column = Protein.Ids
                                , peptide_sequence_column = Stripped.Sequence ) {
  num_peptides <- input_table |>
    distinct( {{protein_id_column}}, {{peptide_sequence_column}}) |>
    count()

  num_peptides[[1,1]]
}





#' @title Plot Peptide and Protein Counts Per Sample
#' @description Generates a plot showing the number of unique proteins and unique
#' peptides identified in each sample.
#'
#' @param input_table The input data frame, typically peptide-level data.
#' @param intensity_column The unquoted column name for intensity values, used to filter out non-quantified entries.
#' @param protein_id_column The unquoted column name for protein identifiers.
#' @param peptide_id_column The unquoted column name for peptide identifiers.
#' @param sample_id_column The unquoted column name for sample identifiers.
#' @param peptide_sequence_column The unquoted column name for peptide sequences (legacy, same as peptide_id_column).
#'
#' @return A ggplot object showing protein and peptide counts per sample.
#' @export
plotPeptidesProteinsCountsPerSampleHelper <- function( input_table
                                                 , intensity_column = Peptide.RawQuantity
                                                 , protein_id_column = Protein.Ids
                                                 , peptide_id_column = Stripped.Sequence
                                                 , sample_id_column = Run
                                                 , peptide_sequence_column = Stripped.Sequence ) {

  num_proteins_per_sample <- input_table |>
    dplyr::filter( !is.na(  {{intensity_column}} )) |>
    distinct( {{sample_id_column}}, {{protein_id_column}} ) |>
    group_by( {{sample_id_column}} ) |>
    summarise( count = n()  ) |>
    ungroup()

  num_peptides_per_sample <- input_table |>
    dplyr::filter( !is.na(  {{intensity_column}} )) |>
    distinct( {{sample_id_column}}, {{protein_id_column}}, {{peptide_id_column}} ) |>
    group_by( {{sample_id_column}}) |>
    summarise( count = n()  ) |>
    ungroup()

  combined_counts <- num_proteins_per_sample |>
    mutate( type = "Protein" ) |>
    bind_rows( num_peptides_per_sample |>
                 mutate( type = "Peptide")) |>
    pivot_wider( id_cols = {{sample_id_column}}
                 , names_from = type
                 , values_from = count)

  output_plot <- combined_counts |>
    ggplot( aes( reorder({{sample_id_column}}, Peptide) )) +
    geom_point(aes(y = Peptide/10, shape="Peptide" ),  show.legend = TRUE) +
    geom_point(aes(y = Protein, shape="Protein"  ),  show.legend = TRUE) +
    scale_y_continuous(name = "Protein",
                       sec.axis = sec_axis(\(x) { x*10 }, name =  "Peptide")) +
    apafTheme() +
    theme( axis.text.x=element_blank()
           , axis.ticks.x=element_blank()
           , panel.grid.major.x = element_blank() ) +
    xlab("Samples") +
    scale_shape_manual ( values = c("Peptide" = 1
                                    , "Protein" = 2) ) +
    labs( shape = "Category")

  output_plot
}


#' @title Filter Peptides by Q-value and Proteotypic Status
#' @description This helper function filters a peptide table to keep only high-confidence,
#' proteotypic peptides based on specified q-value thresholds.
#'
#' @param input_table The input data frame of peptide-spectrum matches.
#' @param qvalue_threshold The maximum q-value for a peptide to be retained.
#' @param global_qvalue_threshold The maximum global q-value to be retained.
#' @param choose_only_proteotypic_peptide A flag (typically 1 or TRUE) to select only proteotypic peptides.
#' @param input_matrix_column_ids A character vector of column names to keep in the output table.
#' @param protein_id_column The unquoted column name for protein identifiers.
#' @param q_value_column The unquoted column name for the peptide q-value.
#' @param global_q_value_column The unquoted column name for the global q-value.
#' @param proteotypic_peptide_sequence_column The unquoted column name indicating if a peptide is proteotypic.
#'
#' @return A data frame containing the filtered peptides.
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



  search_srl_quant_cln <- input_table |>
    dplyr::filter( {{q_value_column}} < qvalue_threshold &
                     {{global_q_value_column}} < global_qvalue_threshold &
                     {{proteotypic_peptide_sequence_column}} == choose_only_proteotypic_peptide ) |>
    dplyr::select(all_of( input_matrix_column_ids))

  search_srl_quant_cln

}

#' @title Roll Up Precursor Quantities to Peptide Level
#' @description This helper function aggregates precursor-level quantities to the
#' peptide level. It sums the quantities of different precursor charge states and
#' modifications for each unique peptide sequence.
#'
#' @param input_table The input data frame of precursor-level data.
#' @param sample_id_column The unquoted column name for sample identifiers.
#' @param protein_id_column The unquoted column name for protein identifiers.
#' @param peptide_sequence_column The unquoted column name for the (stripped) peptide sequence.
#' @param modified_peptide_sequence_column The unquoted column name for the modified peptide sequence.
#' @param precursor_quantity_column The unquoted column name for raw precursor quantity.
#' @param precursor_normalised_column The unquoted column name for normalized precursor quantity.
#' @param core_utilisation The number of cores for parallel processing.
#'
#' @return A data frame with aggregated peptide-level quantities.
#' @export
rollUpPrecursorToPeptideHelper <- function( input_table
                                      , sample_id_column = Run
                                      , protein_id_column = Protein.Ids
                                      , peptide_sequence_column = Stripped.Sequence
                                      , modified_peptide_sequence_column = Modified.Sequence
                                      , precursor_quantity_column = Precursor.Quantity
                                      , precursor_normalised_column = Precursor.Normalised
                                      , core_utilisation) {

  peptide_normalised_tbl <- NA
  if( length(which(is.na(core_utilisation))) == 0 ) {

    peptide_normalised_tbl <- input_table  |>
      group_by( {{sample_id_column}}, {{protein_id_column}}, {{peptide_sequence_column}}, {{modified_peptide_sequence_column}} ) |>
      summarise( Peptide.RawQuantity = sum( {{precursor_quantity_column}} )
                 ,  Peptide.Normalised = sum( {{precursor_normalised_column}} ) ) |>
      ungroup() |>
      group_by( {{sample_id_column}}, {{protein_id_column}}, {{peptide_sequence_column}} ) |>
      summarise( Peptide.RawQuantity = sum( Peptide.RawQuantity )
                 ,  Peptide.Normalised = sum( Peptide.Normalised )
                 ,  peptidoform_count = n()) |>
      ungroup()

  } else {
    peptide_normalised_tbl <- input_table  |>

      group_by( {{sample_id_column}}, {{protein_id_column}}, {{peptide_sequence_column}}, {{modified_peptide_sequence_column}} ) |>
      partition(core_utilisation) |>
      summarise( Peptide.RawQuantity = sum( {{precursor_quantity_column}} )
                 ,  Peptide.Normalised = sum( {{precursor_normalised_column}} ) ) |>
      collect() |>
      ungroup() |>

      group_by( {{sample_id_column}}, {{protein_id_column}}, {{peptide_sequence_column}} ) |>
      partition(core_utilisation) |>
      summarise( Peptide.RawQuantity = sum( Peptide.RawQuantity )
                 ,  Peptide.Normalised = sum( Peptide.Normalised )
                 , peptidoform_count = n() ) |>
      collect() |>
      ungroup()

  }

  peptide_normalised_tbl
}

#' @title Filter Peptides by Intensity Threshold
#' @description This helper filters out peptides that have low abundance across a
#' large proportion of samples.
#'
#' @param input_table The input data frame of peptide-level data.
#' @param min_peptide_intensity_threshold The minimum intensity for a peptide to be considered "present".
#' @param peptides_proportion_of_samples_below_cutoff The maximum proportion of samples where a
#'   peptide can be below the threshold. If a peptide exceeds this, it is removed.
#' @param protein_id_column The unquoted column name for protein identifiers.
#' @param peptide_sequence_column The unquoted column name for peptide sequences.
#' @param peptide_quantity_column The unquoted column name for peptide quantity.
#' @param core_utilisation The number of cores for parallel processing.
#'
#' @return A data frame with low-intensity peptides filtered out.
#' @export
peptideIntensityFilteringHelper <- function(input_table
                                      , min_peptide_intensity_threshold = 15
                                      , peptides_proportion_of_samples_below_cutoff = 1
                                      , protein_id_column = Protein.Ids
                                      , peptide_sequence_column = Stripped.Sequence
                                      , peptide_quantity_column = Peptide.Normalised
                                      , core_utilisation) {
  num_values_per_peptide <- NA

  if( length(which(is.na(core_utilisation))) == 0 ) {
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

  peptide_normalised_pif_cln <- input_table |>
    inner_join ( num_values_per_peptide |>
                   dplyr::select( -num_below_intesnity_treshold, -samples_counts)
                 , by = join_by( {{protein_id_column}}, {{peptide_sequence_column}} ) )


  peptide_normalised_pif_cln


}



#' @title Filter Peptides by Percentage of Missing Values in Groups
#' @description This helper function removes peptides that have a high percentage of
#' missing values within experimental groups. A peptide is removed if it fails to
#' meet the `groupwise_percentage_cutoff` in more than `max_groups_percentage_cutoff`
#' of the groups.
#'
#' @param input_table The input data frame of peptide-level data.
#' @param design_matrix A data frame describing the experimental design.
#' @param sample_id The unquoted column name for sample identifiers in both tables.
#' @param protein_id_column The unquoted column name for protein identifiers.
#' @param peptide_sequence_column The unquoted column name for peptide sequences.
#' @param grouping_variable The unquoted column name in the design matrix for grouping samples.
#' @param groupwise_percentage_cutoff The maximum percentage of missing values allowed within a single group.
#' @param max_groups_percentage_cutoff The maximum percentage of groups that can fail the `groupwise_percentage_cutoff`.
#' @param abundance_threshold An intensity threshold below which values are considered missing.
#' @param abundance_column The name of the abundance column as a string.
#'
#' @return A data frame with the filtered peptides.
#' @export
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

#' @title Filter Proteins by Minimum Number of Peptides
#' @description This helper function filters out proteins that are not supported by a
#' minimum number of unique peptides or peptidoforms.
#'
#' @param input_table The input data frame of peptide-level data.
#' @param num_peptides_per_protein_thresh The minimum number of unique peptides required per protein.
#' @param num_peptidoforms_per_protein_thresh The minimum number of unique peptidoforms (e.g., modified sequences) required per protein.
#' @param protein_id_column The unquoted column name for protein identifiers.
#' @param core_utilisation The number of cores for parallel processing.
#'
#' @return A data frame with filtered protein-peptide data.
#' @export
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


#' @title Filter Samples by Minimum Number of Peptides
#' @description This helper function removes samples that have fewer than a specified
#' number of identified peptides.
#'
#' @param input_table The input data frame of peptide-level data.
#' @param peptides_per_sample_cutoff The minimum number of peptides required for a sample to be kept.
#' @param sample_id_column The unquoted column name for sample identifiers.
#' @param core_utilisation The number of cores for parallel processing.
#' @param inclusion_list A character vector of sample IDs to keep regardless of their peptide count.
#'
#' @return A data frame with low-peptide-count samples removed.
#' @export
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

#--------------------------------------------------------------------------------------------------------------------------

#' @title Log2 Transform with Pseudo-Count
#' @description This function applies a log2 transformation to a numeric matrix. It adds
#' a small pseudo-count to non-zero values to avoid issues with zero values before
#' transformation.
#'
#' @param input_matrix A numeric matrix of abundance data.
#'
#' @return The log2-transformed matrix.
#' @export
log2Transformation <- function(input_matrix) {

  pseudo_count <- min( input_matrix[input_matrix> 0] , na.rm=TRUE)/100
  input_matrix[input_matrix> 0 & !is.na(input_matrix)] <- input_matrix[input_matrix> 0 & !is.na(input_matrix)] + pseudo_count
  input_matrix <- log2(input_matrix)

  return(input_matrix )
}


#--------------------------------------------------------------------------------------------------------------------------

#' @title Get Pairs of Samples for Comparison
#' @description Creates a data frame of all unique pairs of samples within the same
#' technical replicate group, for use in correlation analysis.
#'
#' @param input_table A data frame containing sample IDs and their corresponding replicate group.
#' @param run_id_column The name of the sample ID column as a string.
#' @param replicate_group_column The name of the replicate group column as a string.
#'
#' @return A data frame with three columns: the replicate group, sample ID for the first sample in a pair, and sample ID for the second.
#' @export
getPairsOfSamplesTable <- function ( input_table
                                     , run_id_column
                                     , replicate_group_column) {

  pairs_for_comparison <- input_table |>
    inner_join( input_table, by = join_by(!!rlang::sym( replicate_group_column) )) |>
    dplyr::filter( !!rlang::sym( paste0( run_id_column, ".x")) >  !!rlang::sym(paste0( run_id_column, ".y")) ) |>
    arrange( !!rlang::sym( replicate_group_column) ) |>
    relocate( !!rlang::sym( replicate_group_column)
              , .before=paste0( run_id_column, ".x"))

  pairs_for_comparison
}



#' @title Calculate Pearson Correlation Between Two Samples
#' @description Calculates the Pearson correlation coefficient for peptide abundances
#' between a specified pair of samples.
#'
#' @param ms_filename_x The identifier for the first sample.
#' @param ms_filename_y The identifier for the second sample.
#' @param input_table A data frame containing peptide-level data for multiple samples.
#' @param sample_id_column The unquoted column name for sample identifiers.
#' @param protein_id_column The unquoted column name for protein identifiers.
#' @param peptide_sequence_column The unquoted column name for peptide sequences.
#' @param peptide_normalised_column The name of the normalized abundance column as a string.
#'
#' @return The Pearson correlation value, or `NA` if no common peptides are found.
#' @export
calulatePearsonCorrelation <- function( ms_filename_x, ms_filename_y, input_table
                                        , sample_id_column = Run
                                        , protein_id_column = Protein.Ids
                                        , peptide_sequence_column = Stripped.Sequence
                                        , peptide_normalised_column = "Peptide.Normalised")  {

  tab_x <- input_table |>
    dplyr::filter( {{sample_id_column}} == ms_filename_x )

  tab_y <- input_table |>
    dplyr::filter( {{sample_id_column}} == ms_filename_y )

  merged_tbl <- tab_x |>
    inner_join( tab_y, by=join_by( {{protein_id_column}}, {{peptide_sequence_column}}) )

  # merged_tbl |>
  #   dplyr::filter(!is.na( !!sym(paste0(peptide_normalised_column, ".x")) ) & !is.na( !!sym(paste0(peptide_normalised_column, ".x")))) |>
  #   head() |> print()

  # print( paste(ms_filename_x, ms_filename_y))
  input_x <-  merged_tbl[[ paste0(peptide_normalised_column, ".x") ]]
  input_y <- merged_tbl[[paste0(peptide_normalised_column, ".y")]]
  if( length(input_x) > 0 & length(input_y) >0  ) {
    cor_result <- cor( input_x
                       , input_y
                       , use="pairwise.complete.obs")

    cor_result
  } else {
    return( NA )
  }

}


#' @title Calculate Pearson Correlation for All Sample Pairs
#' @description A helper function that iterates through all pairs of samples within
#' replicate groups and calculates their Pearson correlation using parallel processing.
#'
#' @param samples_id_tbl A data frame with sample IDs and replicate group information.
#' @param run_id_column The name of the sample ID column in `samples_id_tbl` as a string.
#' @param replicate_group_column The name of the replicate group column in `samples_id_tbl` as a string.
#' @param input_table The main data frame with peptide abundance data.
#' @param num_of_cores The number of cores to use for parallel computation.
#' @param sample_id_column The unquoted column name for sample IDs in `input_table`.
#' @param protein_id_column The unquoted column name for protein IDs in `input_table`.
#' @param peptide_sequence_column The unquoted column name for peptide sequences in `input_table`.
#' @param peptide_normalised_column The name of the abundance column in `input_table` as a string.
#'
#' @return A data frame with the Pearson correlation for each pair of samples.
#' @export
calulatePearsonCorrelationForSamplePairsHelper <- function( samples_id_tbl
                                                      , run_id_column = "ms_filename"
                                                      , replicate_group_column = "general_sample_info"
                                                      , input_table
                                                      , num_of_cores = 1
                                                      , sample_id_column = Run
                                                      , protein_id_column = Protein.Ids
                                                      , peptide_sequence_column = Stripped.Sequence
                                                      , peptide_normalised_column = "Peptide.Normalised") {


  pairs_for_comparison <- getPairsOfSamplesTable(samples_id_tbl
                                                 , run_id_column = run_id_column
                                                 , replicate_group_column = replicate_group_column)

  plan(multisession, workers = num_of_cores)

  pearson_correlation_per_pair <- pairs_for_comparison |>
    mutate( pearson_correlation = furrr::future_map2_dbl( !!rlang::sym( paste0( run_id_column, ".x"))
                                                          , !!rlang::sym(paste0( run_id_column, ".y"))
                                                          , \(x,y){
                                                            input_table_filt <- input_table |>
                                                              dplyr::filter( {{sample_id_column}} == x | {{sample_id_column}} == y)

                                                            calulatePearsonCorrelation( ms_filename_x = x
                                                                                        , ms_filename_y = y
                                                                                        , input_table = input_table_filt
                                                                                        , sample_id_column = {{sample_id_column}}
                                                                                        , protein_id_column = {{protein_id_column}}
                                                                                        , peptide_sequence_column = {{peptide_sequence_column}}
                                                                                        , peptide_normalised_column = {{peptide_normalised_column}}) }))

  pearson_correlation_per_pair

}


#' @title Filter Samples by Peptide Correlation Threshold
#' @description Removes samples that are poorly correlated with their technical replicates,
#' based on a pre-calculated table of pairwise correlations.
#'
#' @param pearson_correlation_per_pair A data frame of pairwise sample correlations.
#' @param peptide_keep_samples_with_min_num_peptides The main peptide data frame to be filtered.
#' @param min_pearson_correlation_threshold The minimum correlation for a sample to be kept.
#' @param filename_column_x The unquoted column name for the first sample in a pair.
#' @param filename_column_y The unquoted column name for the second sample in a pair.
#' @param correlation_column The unquoted column name for the correlation value.
#' @param filename_id_column The name of the sample ID column in the main data frame as a string.
#'
#' @return A filtered data frame containing only the high-correlation samples.
#' @export
filterSamplesByPeptideCorrelationThreshold <- function(pearson_correlation_per_pair
                                                , peptide_keep_samples_with_min_num_peptides
                                                , min_pearson_correlation_threshold = 0.75
                                                , filename_column_x = ms_filename.x
                                                , filename_column_y = ms_filename.y
                                                , correlation_column = pearson_correlation
                                                , filename_id_column = "Run" ) {
  # Samples to keep include all those pairs of samples with correlation score passing threshold
  samples_to_keep <-  pearson_correlation_per_pair |>
    dplyr::filter( {{correlation_column}} >= min_pearson_correlation_threshold) |>
    pivot_longer( cols =c({{filename_column_x}}, {{filename_column_y}})
                  , values_to = filename_id_column ) |>
    dplyr::distinct( !!rlang::sym(filename_id_column ) )

  samples_above_correlation_theshold <- peptide_keep_samples_with_min_num_peptides |>
    inner_join( samples_to_keep
               , by=join_by(  !!rlang::sym(filename_id_column ) == !!rlang::sym(filename_id_column ) )) |>
    distinct()

  samples_above_correlation_theshold

}

#' @title Find Sample Pairs Below Correlation Threshold
#' @description Identifies pairs of samples that fall below a specified Pearson
#' correlation threshold.
#'
#' @param pearson_correlation_per_pair A data frame of pairwise sample correlations.
#' @param peptide_keep_samples_with_min_num_peptides The main peptide data frame (used to define the universe of samples).
#' @param min_pearson_correlation_threshold The correlation threshold.
#' @param filename_column_x The unquoted column name for the first sample in a pair.
#' @param filename_column_y The unquoted column name for the second sample in a pair.
#' @param correlation_column The unquoted column name for the correlation value.
#' @param filename_id_column The name of the sample ID column in the main data frame as a string.
#'
#' @return A data frame listing the sample pairs that are below the correlation threshold.
#' @export
findSamplesPairBelowPeptideCorrelationThreshold <- function(pearson_correlation_per_pair
                                                     , peptide_keep_samples_with_min_num_peptides
                                                     , min_pearson_correlation_threshold = 0.75
                                                     , filename_column_x = ms_filename.x
                                                     , filename_column_y = ms_filename.y
                                                     , correlation_column = pearson_correlation
                                                     , filename_id_column = "Run" ) {

  samples_to_keep <-  pearson_correlation_per_pair |>
    dplyr::filter( {{correlation_column}} >= min_pearson_correlation_threshold) |>
    pivot_longer( cols =c({{filename_column_x}}, {{filename_column_y}})
                  , values_to = filename_id_column ) |>
    dplyr::distinct( !!rlang::sym(filename_id_column ) )

  samples_below_correlation_theshold <- pearson_correlation_per_pair |>
    pivot_longer( cols =c({{filename_column_x}}, {{filename_column_y}})
                  , values_to = filename_id_column ) |>
    dplyr::distinct( !!rlang::sym(filename_id_column ) ) |>
    innner_join( samples_to_keep
               , by= join_by( !!rlang::sym(filename_id_column ) == !!rlang::sym(filename_id_column ) ) )

  samples_below_correlation_theshold

}


#---------------------------------------------------------------------------------------


#' @title Helper to Filter Samples by Protein Correlation
#' @description A helper function that filters a protein intensity matrix, keeping only
#' the columns (samples) that meet a certain correlation threshold with their replicates.
#'
#' @param pearson_correlation_per_pair A data frame of pairwise sample correlations.
#' @param protein_intensity_table A wide-format data frame with proteins as rows and samples as columns.
#' @param min_pearson_correlation_threshold The minimum correlation for a sample to be kept.
#' @param filename_column_x The unquoted column name for the first sample in a pair.
#' @param filename_column_y The unquoted column name for the second sample in a pair.
#' @param protein_id_column The unquoted column name for protein identifiers.
#' @param correlation_column The unquoted column name for the correlation value.
#'
#' @return A filtered protein intensity data frame.
#' @export
filterSamplesByProteinCorrelationThresholdHelper <- function(pearson_correlation_per_pair
                                                       , protein_intensity_table
                                                       , min_pearson_correlation_threshold = 0.75
                                                       , filename_column_x = ms_filename.x
                                                       , filename_column_y = ms_filename.y
                                                       , protein_id_column = Protein.Ids
                                                       , correlation_column = pearson_correlation ) {

  # All Samples
  all_samples <-  pearson_correlation_per_pair |>
    pivot_longer( cols =c({{filename_column_x}}, {{filename_column_y}})
                  , values_to = "temp_column" ) |>
    dplyr::distinct( temp_column )

  # Samples to keep include all those pairs of samples with correlation score passing threshold
  samples_to_keep <-  pearson_correlation_per_pair |>
    dplyr::filter( {{correlation_column}} >= min_pearson_correlation_threshold) |>
    pivot_longer( cols =c({{filename_column_x}}, {{filename_column_y}})
                  , values_to = "temp_column" ) |>
    dplyr::distinct( temp_column )

  # Samples to keep anyway
  samples_to_keep_anyway <-setdiff( setdiff(colnames(protein_intensity_table), (all_samples |> dplyr::pull( temp_column )))
                                    ,  as_string({{protein_id_column}})  )

  print( samples_to_keep_anyway)

  # Samples in the table to keep
  samples_to_keep_subset <- colnames(protein_intensity_table)[colnames(protein_intensity_table) %in% (samples_to_keep |> dplyr::pull( temp_column ))]

  samples_above_correlation_theshold <- protein_intensity_table |>
    dplyr::select( {{protein_id_column}}, all_of( c(samples_to_keep_anyway, samples_to_keep_subset)))

  samples_above_correlation_theshold

}

#---------------------------------------------------------------------------------------
#' @title Remove Peptides with Only One Replicate
#' @description This helper function removes peptides that are observed in only one
#' technical replicate within any sample group. This helps to filter out less
#' reliable identifications.
#'
#' @param input_table The input data frame of peptide-level data.
#' @param samples_id_tbl A data frame with sample IDs and their corresponding replicate groups.
#' @param input_table_sample_id_column The unquoted column name for sample IDs in `input_table`.
#' @param sample_id_tbl_sample_id_column The unquoted column name for sample IDs in `samples_id_tbl`.
#' @param replicate_group_column The unquoted column name for replicate groups in `samples_id_tbl`.
#' @param protein_id_column The unquoted column name for protein identifiers.
#' @param peptide_sequence_column The unquoted column name for peptide sequences.
#' @param core_utilisation The number of cores for parallel processing.
#'
#' @return A data frame with the filtered peptides.
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

#---------------------------------------------------------------------------------------
#' @title Remove Proteins with Only One Replicate
#' @description This function removes proteins that are observed in only one technical
#' replicate within any sample group.
#'
#' @param input_table The input data frame of protein-level data.
#' @param samples_id_tbl A data frame with sample IDs and their corresponding replicate groups.
#' @param input_table_sample_id_column The unquoted column name for sample IDs in `input_table`.
#' @param sample_id_tbl_sample_id_column The unquoted column name for sample IDs in `samples_id_tbl`.
#' @param replicate_group_column The unquoted column name for replicate groups in `samples_id_tbl`.
#' @param protein_id_column The unquoted column name for protein identifiers.
#' @param core_utilisation The number of cores for parallel processing.
#'
#' @return A data frame with the filtered proteins.
#' @export
removeProteinWithOnlyOneReplicate <- function(input_table
                                               , samples_id_tbl
                                               , input_table_sample_id_column = Run
                                               , sample_id_tbl_sample_id_column  =  ms_filename
                                               , replicate_group_column = general_sample_info
                                               , protein_id_column = Protein.Ids
                                               , core_utilisation ) {

  # Count the number of technical replicates per sample and protein combination
  num_tech_reps_per_sample_and_protein <- NA
  if( length(which(is.na(core_utilisation))) == 0 ) {
    num_tech_reps_per_sample_and_protein <- input_table |>
      left_join( samples_id_tbl, by=join_by( {{input_table_sample_id_column}} == {{sample_id_tbl_sample_id_column}} ) ) |>
      group_by( {{replicate_group_column}}, {{protein_id_column}}) |>
      #partition(core_utilisation) |>
      summarise(counts = n() ) |>
      #collect() |>
      ungroup()
  } else {
    num_tech_reps_per_sample_and_protein <- input_table |>
      left_join( samples_id_tbl, by=join_by( {{input_table_sample_id_column}} == {{sample_id_tbl_sample_id_column}} ) ) |>
      group_by( {{replicate_group_column}}, {{protein_id_column}}) |>
      partition(core_utilisation) |>
      summarise(counts = n() ) |>
      collect() |>
      ungroup()
  }

  # Any proteins found in more than one replicates in any patient will be kept for analysis
  removed_proteins_with_only_one_replicate <- input_table |>
    inner_join( num_tech_reps_per_sample_and_protein |>
                  dplyr::filter( counts >  1) |>
                  dplyr::select(-counts, -{{replicate_group_column}}) |>
                  distinct()
                , by=join_by( {{protein_id_column}}) )  |>
    distinct()

  removed_proteins_with_only_one_replicate
}

##-----------------------------------------------------------------------------------------

#' @title Impute Missing Peptide-Level Values
#' @description This helper function performs missing value imputation for peptide-level
#' data. It imputes missing values within a replicate group by using the average
#' of the observed values, but only if the proportion of missing values is below
#' a certain threshold.
#'
#' @param input_table The input data frame of peptide-level data.
#' @param metadata_table A data frame with sample metadata, including replicate group information.
#' @param input_table_sample_id_column The unquoted column name for sample IDs in `input_table`.
#' @param sample_id_tbl_sample_id_column The unquoted column name for sample IDs in `metadata_table`.
#' @param replicate_group_column The unquoted column name for replicate groups in `metadata_table`.
#' @param protein_id_column The unquoted column name for protein identifiers.
#' @param peptide_sequence_column The unquoted column name for peptide sequences.
#' @param quantity_to_impute_column The unquoted column name for the abundance values to be imputed.
#' @param imputed_value_column The unquoted column name for the new column that will hold the imputed values.
#' @param hek_string A character string to identify control samples (e.g., HEK) to be excluded from imputation.
#' @param proportion_missing_values The maximum proportion of missing values within a group to allow for imputation.
#' @param core_utilisation The number of cores for parallel processing.
#'
#' @return A data frame with missing values imputed.
#' @export
peptideMissingValueImputationHelper <- function( input_table
                                           , metadata_table
                                           , input_table_sample_id_column = Run
                                           , sample_id_tbl_sample_id_column  =  ms_filename
                                           , replicate_group_column = general_sample_info
                                           , protein_id_column = Protein.Ids
                                           , peptide_sequence_column = Stripped.Sequence
                                           , quantity_to_impute_column = Peptide.Normalised
                                           , imputed_value_column = Peptide.Imputed
                                           , hek_string = "HEK"
                                           , proportion_missing_values = 0.50
                                           , core_utilisation ) {

  # Max number of technical replicates per group
  num_tech_rep_per_sample <-  metadata_table  |>
    dplyr::filter( !str_detect( {{replicate_group_column}}, hek_string))  |>
    distinct( {{sample_id_tbl_sample_id_column}}, {{replicate_group_column}}) |>
    group_by( {{replicate_group_column}}) |>
    summarise(total_num_tech_rep = n()) |>
    ungroup()

  # Count the number of technical replicates per sample and peptide combination
  num_tech_reps_per_sample_and_peptide <- NA

  if( length(which(is.na(core_utilisation))) == 0 ) {
    num_tech_reps_per_sample_and_peptide <- input_table |>
      left_join( metadata_table
                 , by=join_by( {{input_table_sample_id_column}} == {{sample_id_tbl_sample_id_column}} ) ) |>
      dplyr::filter( !str_detect( {{replicate_group_column}}, hek_string))  |>
      dplyr::filter( !is.na({{quantity_to_impute_column}}) ) |>
      distinct( {{replicate_group_column}}, {{protein_id_column}}, {{peptide_sequence_column}}, {{quantity_to_impute_column}}) |>
      group_by( {{replicate_group_column}}, {{protein_id_column}}, {{peptide_sequence_column}}) |>
      #partition(core_utilisation) |>
      summarise( num_tech_rep = n()
                 , average_value = mean({{quantity_to_impute_column}}, na.rm=TRUE )) |>
      #collect() |>
      ungroup()
  } else {
    num_tech_reps_per_sample_and_peptide <- input_table |>
      left_join( metadata_table
                 , by=join_by( {{input_table_sample_id_column}} == {{sample_id_tbl_sample_id_column}} ) ) |>
      dplyr::filter( !str_detect( {{replicate_group_column}}, hek_string))  |>
      dplyr::filter( !is.na({{quantity_to_impute_column}}) ) |>
      distinct( {{replicate_group_column}}, {{protein_id_column}}, {{peptide_sequence_column}}, {{quantity_to_impute_column}}) |>
      group_by( {{replicate_group_column}}, {{protein_id_column}}, {{peptide_sequence_column}}) |>
      partition(core_utilisation) |>
      summarise( num_tech_rep = n()
                 , average_value = mean({{quantity_to_impute_column}}, na.rm=TRUE)) |>
      collect() |>
      ungroup()

  }

  ## Calculate proportion of replicates in a group that is missing
  rows_needing_imputation_temp <-  num_tech_reps_per_sample_and_peptide |>
    left_join( num_tech_rep_per_sample
               , by = join_by( {{replicate_group_column}} ) )


  print(rows_needing_imputation_temp)

  rows_needing_imputation <-   rows_needing_imputation_temp |>
    dplyr::filter(    (1 - num_tech_rep / total_num_tech_rep ) < proportion_missing_values )

  get_combinations_part_1 <- metadata_table |>
    distinct( {{sample_id_tbl_sample_id_column}}, {{replicate_group_column}}) |>
    left_join(  input_table |>
                  distinct( {{input_table_sample_id_column}}, {{protein_id_column}}, {{peptide_sequence_column}})
                , by =join_by( {{sample_id_tbl_sample_id_column}} == {{input_table_sample_id_column}}))

  all_peptides_combination <- metadata_table |>
    distinct( {{sample_id_tbl_sample_id_column}}, {{replicate_group_column}}) |>
    group_by({{replicate_group_column}} ) |>
    nest( data = c({{sample_id_tbl_sample_id_column}}) ) |>
    left_join( get_combinations_part_1 |>
                 dplyr::select( -{{sample_id_tbl_sample_id_column}}) |>
                 dplyr::distinct( {{replicate_group_column}}, {{protein_id_column}}, {{peptide_sequence_column}})
               , by = join_by( {{replicate_group_column}}))  |>
    unnest( data ) |>
    ungroup({{replicate_group_column}})


  make_imputation <- all_peptides_combination |>
    left_join( input_table
               , by = join_by( {{sample_id_tbl_sample_id_column}} == {{input_table_sample_id_column}}
                               , {{protein_id_column}} == {{protein_id_column}}
                               , {{peptide_sequence_column}} == {{peptide_sequence_column}} ) ) |>
    left_join(rows_needing_imputation
              , by = join_by( {{replicate_group_column}}
                              , {{protein_id_column}}
                              , {{peptide_sequence_column}}  ))  |>
    dplyr::filter(!is.na({{protein_id_column}}) & !is.na( {{peptide_sequence_column}} )) |>
    mutate( is_imputed = case_when (is.na({{quantity_to_impute_column}})
                                    & !is.na(average_value)  ~ TRUE
                                    , TRUE ~ FALSE) ) |>
    mutate ( {{imputed_value_column}} := case_when (is.na({{quantity_to_impute_column}})
                                                    & !is.na(average_value)  ~ average_value
                                                    , TRUE ~ {{quantity_to_impute_column}} ) ) |>
    dplyr::select( -num_tech_rep
                   , - average_value
                   , - total_num_tech_rep
                   , - {{replicate_group_column}} ) |>
    dplyr::rename( {{input_table_sample_id_column}} := {{sample_id_tbl_sample_id_column}})

  make_imputation
}


#---------------------------------------------------------------------------------------

#' @title Calculate Percent Missing Peptides per Replicate
#' @description Calculates the percentage of missing peptide identifications for each
#' sample (replicate) and merges this information with sample metadata.
#'
#' @param input_table The input data frame of peptide-level data.
#' @param metadata_table A data frame with sample metadata.
#' @param protein_id_column The unquoted column name for protein identifiers.
#' @param intensity_column The unquoted column name for intensity values.
#' @param replicate_id_column The unquoted column name for sample/replicate identifiers.
#' @param peptide_sequence_column The unquoted column name for peptide sequences.
#'
#' @return A data frame with metadata and a new `percent_missing` column.
#' @export
calculatePercentMissingPeptidePerReplicate <- function( input_table
                                                        , metadata_table
                                                        , protein_id_column = Protein.Ids
                                                        , intensity_column = Peptide.Normalised
                                                        , replicate_id_column = Run
                                                        , peptide_sequence_column = Stripped.Sequence ) {

  # Total number of peptides with values per run
  total_num_of_peptides_with_values_per_run <- input_table |>
    left_join( metadata_table, by=join_by({{replicate_id_column}})) |>
    dplyr::filter( !is.na( {{intensity_column}} )) |>
    group_by( {{replicate_id_column}}) |>
    summarise(counts = n()) |>
    ungroup()

  # Total number of peptides
  total_num_of_peptides <- input_table  |>
    left_join( metadata_table, by=join_by({{replicate_id_column}}) ) |>
    distinct( {{protein_id_column}}, {{peptide_sequence_column}} ) |>
    nrow()

  percent_missing_per_run <- total_num_of_peptides_with_values_per_run |>
    mutate( percent_missing = (1 - (counts / total_num_of_peptides)) * 100 ) |>
    left_join( metadata_table, by=join_by({{replicate_id_column}}))

  return( percent_missing_per_run )
}



#' @title Calculate Percent Missing Proteins per Replicate
#' @description Calculates the percentage of missing protein quantifications for each
#' sample (replicate) and merges this information with sample metadata.
#'
#' @param input_table The input data frame of protein-level data.
#' @param metadata_table A data frame with sample metadata.
#' @param protein_id_column The unquoted column name for protein identifiers.
#' @param intensity_column The unquoted column name for intensity values.
#' @param replicate_id_column The unquoted column name for sample/replicate identifiers.
#'
#' @return A data frame with metadata and a new `percent_missing` column.
#' @export
calculatePercentMissingProteinPerReplicate <- function( input_table
                                                        , metadata_table
                                                        , protein_id_column = Protein.Ids
                                                        , intensity_column = Log2.Protein.Imputed
                                                        , replicate_id_column = Run ) {

  # Total number of peptides with values per run
  total_num_of_proteins_with_values_per_run <- input_table |>
    left_join( metadata_table, by=join_by({{replicate_id_column}})) |>
    dplyr::filter( !is.na( {{intensity_column}} )) |>
    group_by( {{replicate_id_column}}) |>
    summarise(num_proteins_with_values = n()) |>
    ungroup()

  # Total number of peptides
  total_num_of_proteins <- input_table  |>
    left_join( metadata_table, by=join_by({{replicate_id_column}}) ) |>
    distinct( {{protein_id_column}} ) |>
    nrow()

  percent_missing_per_run <- total_num_of_proteins_with_values_per_run |>
    mutate( percent_missing = (1 - (num_proteins_with_values / total_num_of_proteins)) * 100 ) |>
    left_join( metadata_table, by=join_by( {{replicate_id_column}}))

  return( percent_missing_per_run )
}



#' @title Plot Histogram of Percent Missing Values
#' @description Generates a histogram to visualize the distribution of the percentage
#' of missing values across samples.
#'
#' @param percent_missing_table A data frame containing a column with the percentage of missing values.
#' @param percent_missing_column The unquoted column name for the percent missing values.
#'
#' @return A ggplot object representing the histogram.
#' @export
plotHistogramOfPercentMissingPerIndvidual <- function( percent_missing_table
                                                       , percent_missing_column = percent_missing) {
  percent_missing_table |>
    ggplot( aes( {{percent_missing_column}})) +
    geom_histogram() +
    apafTheme() +
    xlab("Percent Missing") +
    ylab("Count")


}


#---------------------------------------------------------------------------------------
#' @title Remove Peptides Found Only in Control Samples
#' @description This function filters out peptides that are identified only in control
#' samples (e.g., HEK293) and not in any of the experimental cohort samples.
#'
#' @param input_table The input data frame of peptide-level data.
#' @param metadata_table A data frame with sample metadata.
#' @param input_table_sample_id_column The name of the sample ID column in `input_table` as a string.
#' @param sample_id_tbl_sample_id_column The name of the sample ID column in `metadata_table` as a string.
#' @param protein_id_column The unquoted column name for protein identifiers.
#' @param peptide_sequence_column The unquoted column name for peptide sequences.
#' @param hek_string A character string to identify control samples.
#' @param general_sample_info The unquoted column name in `metadata_table` that contains the sample group information.
#' @param core_utilisation The number of cores for parallel processing.
#'
#' @return A data frame with control-only peptides removed.
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


#------------------------------------------------------------------------------------------------

#' @title Get Data for a Single RLE Plot
#' @description Calculates the data required to generate a Relative Log Expression (RLE)
#' plot from a matrix of abundance data.
#'
#' @param input_matrix A numeric matrix with samples as columns and features (proteins/peptides) as rows.
#'
#' @return A data frame in long format suitable for plotting with ggplot2.
#' @export
getOneRlePlotData <- function( input_matrix ) {


  #if(!( length(which( is.na(input_matrix[, 1]) | is.nan(input_matrix[, 1]) | is.infinite(input_matrix[, 1]) )) > 0 )){

  input_matrix[is.infinite(input_matrix)  | is.nan(input_matrix) ] <- NA

  deviations <- input_matrix - Biobase::rowMedians(input_matrix, na.rm=TRUE)

  stats <-  graphics::boxplot(
    deviations,
    outcol="lightgray",
    cex=0.1,
    cex.axis=0.7,
    las=2,
    outline=FALSE)

  rownames(stats$stats) <- c("lower whisker", "lower hinge", "median", "upper hinge", "upper whisker")
  colnames(stats$stats ) <- colnames(deviations )

  results <- stats$stats |>
    as.data.frame() |>
    rownames_to_column("Quantiles") |>
    pivot_longer( cols = !contains("Quantiles")) |>
    mutate( Quantiles = factor( Quantiles, levels=rev(c( "lower whisker", "lower hinge", "median", "upper hinge", "upper whisker" ))))

  # print(head( results) )



  return(results)
  #}

}

#' @title Plot RLE QC Plot
#' @description Generates a Relative Log Expression (RLE) plot from pre-calculated data.
#'
#' @param input_table A data frame in long format, typically from `getOneRlePlotData`.
#' @param x_value The unquoted column name for the x-axis (sample names).
#' @param y_value The unquoted column name for the y-axis (deviations).
#' @param quantiles_column The unquoted column name for the quantile groups.
#'
#' @return A ggplot object representing the RLE plot.
#' @export
plotRleQc <- function( input_table
                       , x_value = name
                       , y_value = value
                       , quantiles_column = Quantiles ) {

  rle_results <- input_table |>
    ggplot(aes( x={{x_value}}, y={{y_value}},  group={{quantiles_column}}, col={{quantiles_column}})) +
    geom_line()  +
    apafTheme() +
    theme(axis.text.x = element_blank()) +
    xlab("Samples") +
    ylab("Relative log expression") +
    labs(col = "Boxplot features")

  rle_results
}

#------------------------------------------------------------------------------------------------
#' @title Scale, Center, and Fill Missing Values
#' @description This function scales and centers a numeric matrix and then fills any
#' remaining missing values with a value derived from the minimum of the scaled data.
#'
#' @param input_matrix A numeric matrix of abundance data.
#'
#' @return A numeric matrix that has been scaled, centered, and has missing values filled.
#' @export
scaleCenterAndFillMissing <- function( input_matrix) {

  input_matrix_scaled <- scale(input_matrix, center = TRUE, scale = TRUE)


  min_data_point <- min(input_matrix_scaled, na.rm=TRUE)


  input_matrix_scaled_fill_missing <-  input_matrix_scaled
  input_matrix_scaled_fill_missing[is.na(input_matrix_scaled_fill_missing)] <- min_data_point*2

  input_matrix_scaled_fill_missing
}

#------------------------------------------------------------------------------------------------
#' @title Compare UMAP Components in a Pairs Plot
#' @description Generates a pairs plot (using `GGally::ggpairs`) to visualize the
#' relationships between the first few UMAP components, colored by a covariate.
#'
#' @param input_table A data frame containing UMAP components (e.g., V1, V2, ...).
#' @param columns A character vector of the UMAP component columns to plot.
#' @param covariate The unquoted column name of the covariate to use for coloring.
#'
#' @return A `ggpairs` plot object.
#' @export
compareUmapComponentsPairs <- function(input_table, columns = c("V1", "V2","V3","V4"), covariate) {

  pm <- umap_data |>
    ggpairs( columns = columns, ggplot2::aes(colour = {{covariate}}), legend = 1)  +
    apafTheme()

  pm

}


#------------------------------------------------------------------------------------------------
#' @title Create a UMAP Factor Plot
#' @description Generates a scatter plot of two UMAP components, with points colored
#' by a factor variable using a specified color rule.
#'
#' @param input_data A data frame containing UMAP components and the coloring variable.
#' @param header The unquoted column name of the factor variable to use for coloring.
#' @param legend_label The label for the plot legend.
#' @param x The unquoted column name for the x-axis component.
#' @param y The unquoted column name for the y-axis component.
#' @param colour_rule A named vector specifying the colors for each level of the factor.
#'
#' @return A ggplot object.
#' @export
umap_factor_plot <- function(input_data, header, legend_label, x = V1, y = V2, colour_rule) {


  input_data |>
    mutate( !!sym( {{header}}) := factor( !!sym( {{header}})) ) |>
    ggplot(aes( {{x}}, {{y}}, color = !!sym( {{header}}) )) +
    geom_point() +
    scale_colour_manual( name = legend_label
                         , values=colour_rule
                         , breaks=names( colour_rule)) +
    apafTheme()


}

#' @title Save a List of Plots to a Single PDF
#' @description Iterates through a list of plot objects and saves them to a single
#' multi-page PDF file.
#'
#' @param list A list of plot objects that can be printed.
#' @param filename The path and name of the output PDF file.
#'
#' @return Invisibly returns `NULL`. This function is used for its side effect of creating a file.
#' @export
saveListOfPdfs <- function(list, filename) {
  #start pdf
  cairo_pdf(filename)

  #loop
  #purrr::walk( list, print)
  for (p in list) {
    print(p)
  }

  #end pdf
  dev.off()

  invisible(NULL)
}

#------------------------------------------------------------------------------------------------


#' @title Get Sample Correlation Matrix
#' @description Calculates a Pearson correlation matrix for all non-control samples.
#'
#' @param input_table A wide-format data frame with features as rows and samples as columns.
#' @param metadata_tbl A data frame with sample metadata.
#' @param is_HEK_column The unquoted column name in `metadata_tbl` indicating control samples.
#' @param use The `use` argument for the `cor` function (e.g., "pairwise.complete.obs").
#' @param method The `method` argument for the `cor` function (e.g., "pearson").
#'
#' @return A numeric matrix of sample-sample correlations.
#' @export
getSamplesCorrelationMatrix <- function(input_table
                                        , metadata_tbl
                                        , is_HEK_column = is_HEK
                                        , use ="pairwise.complete.obs"
                                        , method = "pearson") {

  without_hek_samples <- metadata_tbl |>
    dplyr::filter( {{is_HEK_column}} == FALSE) |>
    pull(Run)

  correlation_samples_to_use <- intersect( colnames(input_table), without_hek_samples) |> sort()

  correlation_between_samples <-  cor(input_table[, correlation_samples_to_use], use = use, method=method)
  which(is.na(correlation_between_samples))
  correlation_between_samples[is.na(correlation_between_samples)] <- 0

  return( correlation_between_samples)
}

#------------------------------------------------------------------------------------------------
#' @title Get a Large Categorical Color Palette
#' @description Combines several `RColorBrewer` palettes to create a large, diverse
#' color palette for categorical variables with many levels.
#'
#' @return A character vector of color hex codes.
#' @export
getCategoricalColourPalette <- function() {
  set1_colour <- brewer.pal(9,'Set1')
  set2_colour <- brewer.pal(8,'Set2')
  set3_colour <- brewer.pal(12,'Set3')
  pastel1_colour <- brewer.pal(9,'Pastel1')
  pastel2_colour <- brewer.pal(8,'Pastel2')
  dark2_colour <- brewer.pal(8,'Dark2')
  accent_colour <- brewer.pal(8,'Accent')
  paired_colour <- brewer.pal(12,'Paired')

  set1_2_3_colour <- c( set1_colour, set2_colour, set3_colour
                        , pastel1_colour, pastel2_colour, dark2_colour
                        , accent_colour, paired_colour )

  return(set1_2_3_colour)
}



#' @title Get a Continuous Color Palette
#' @description Creates a color palette for a continuous variable, mapping colors to
#' the unique values present in the data.
#'
#' @param metadata_tbl The metadata table containing the continuous variable.
#' @param column_name The name of the continuous variable column as a string.
#' @param palette_name The name of the `RColorBrewer` palette to use.
#' @param na_colour The color to use for NA values.
#'
#' @return A named character vector of color hex codes.
# getOneContinousPalette <- function(metadata_tbl, column_name, palette_name, num_colours=9) {
#
#   list_of_values <-  metadata_tbl |>
#     dplyr::select( all_of(column_name)  ) |>
#     distinct()  |>
#     dplyr::filter(!is.na(!!sym(column_name))) |>
#     arrange( !!sym( column_name)  ) |>
#     pull()
#
#   min_value <- min(list_of_values)
#   max_value <-  max(list_of_values)
#
#   if(min_value > 1) {
#     min_value <- floor(min_value)
#   }
#
#   if(max_value > 1) {
#     max_value <- ceiling(max_value)
#   }
#
#   list_of_names <- levels(cut(list_of_values, breaks=seq( min_value, max_value, length.out=num_colours) ))
#
#   na_name <- metadata_tbl |>
#     dplyr::select( all_of(column_name)  ) |>
#     distinct()  |>
#     dplyr::filter(is.na(!!sym(column_name))) |>
#     pull()
#
#   list_of_colours <- brewer.pal(num_colours, palette_name)
#   names(list_of_colours) <- list_of_names
#
#   if(length(na_name) == 1) {
#      new_list_of_colours <- c(list_of_colours, NA)
#      names(new_list_of_colours) <- c(names(list_of_colours), "NA")
#      return( new_list_of_colours )
#   }
#
#   return( list_of_colours )
# }

#' @export
getOneContinousPalette <- function(metadata_tbl, column_name, palette_name, na_colour = "white") {
  number_of_values <- metadata_tbl |>
    dplyr::select( all_of(column_name)  ) |>
    dplyr::filter(!is.na(!!sym(column_name))) |>
    distinct()  |>
    arrange( !!sym( column_name)  ) |>
    pull() |>
    length()

  list_of_names <- metadata_tbl |>
    dplyr::select( all_of(column_name)  )  |>
    dplyr::filter(!is.na(!!sym(column_name))) |>
    distinct()  |>
    arrange( !!sym( column_name)  )  |>
    pull()

  na_name <- metadata_tbl |>
    dplyr::select( all_of(column_name)  ) |>
    distinct()  |>
    dplyr::filter(is.na(!!sym(column_name))) |>
    pull()

  list_of_colours <- brewer.pal(number_of_values, palette_name)
  names(list_of_colours) <- list_of_names

  if(length(na_name) == 1) {
    new_list_of_colours <- c(list_of_colours, na_colour )
    names(new_list_of_colours) <- c(names(list_of_colours), "NA")
    return( new_list_of_colours )
  }

  return( list_of_colours )
}

#' @title Get Color Rules for Multiple Continuous Variables
#' @description Generates a list of color palettes, one for each specified continuous
#' variable in the metadata.
#'
#' @param metadata_tbl The metadata table.
#' @param metadata_column_labels A named vector for relabeling columns for display.
#' @param metadata_column_selected The selected columns to process.
#' @param continous_scale_columns A character vector of the names of continuous columns.
#' @param na_colour The color to use for NA values.
#'
#' @return A named list of color palettes.
#' @export
getContinousColourRules <- function( metadata_tbl
                                     , metadata_column_labels
                                     , metadata_column_selected
                                     , continous_scale_columns
                                     , na_colour = "white" ) {

  metadata_column_labels_copy <- metadata_column_labels
  names( metadata_column_labels_copy) <- metadata_column_selected

  list_of_continuous_colour_palette <- c( "Greys", "Blues", "Greens", "Purples", "Reds", "Oranges", "BuGn"
                                          , "BuPu", "GnBu", "OrRd", "PuBu", "PuBuGn", "PuRd"
                                          ,  "RdPu", "YlGn", "YlGnBu", "YlOrBr", "YlOrRd" )

  if( length(continous_scale_columns) > length( list_of_continuous_colour_palette )) {
    list_of_continuous_colour_palette <- rep( list_of_continuous_colour_palette, length.out = length(continous_scale_columns))
  }

  list_of_continous_colour_rules <- purrr::map2(continous_scale_columns
                                                ,  list_of_continuous_colour_palette[seq_along(continous_scale_columns)]
                                                , \(column, palette_name) { getOneContinousPalette(metadata_tbl
                                                                                                   , column
                                                                                                   , palette_name
                                                                                                   , na_colour = na_colour) } )

  names( list_of_continous_colour_rules) <- metadata_column_labels_copy[continous_scale_columns]

  return(list_of_continous_colour_rules)

}

#' @title Get Combined Color Rules for Annotations
#' @description A wrapper function that generates a combined list of color rules for
#' both categorical and continuous variables, suitable for heatmap annotations.
#'
#' @param metadata_tbl The metadata table.
#' @param metadata_column_labels A named vector for relabeling columns.
#' @param metadata_column_selected A character vector of the columns to process.
#' @param categorical_columns A character vector of categorical column names.
#' @param continous_scale_columns A character vector of continuous column names.
#' @param ms_machine_column The name of the mass spec machine column (treated specially).
#' @param sample_id_column The unquoted column name for sample IDs.
#' @param columns_to_exclude A character vector of columns to exclude from the final rules.
#' @param na_colour The color to use for NA values.
#'
#' @return A named list of color palettes.
#' @export
getCategoricalAndContinuousColourRules <- function( metadata_tbl
                                                    , metadata_column_labels
                                                    , metadata_column_selected
                                                    , categorical_columns
                                                    , continous_scale_columns
                                                    , ms_machine_column
                                                    , sample_id_column = Run
                                                    , columns_to_exclude
                                                    , na_colour = "white" ) {

  metadata_column_labels_copy <- metadata_column_labels
  names( metadata_column_labels_copy) <- metadata_column_selected

  if( ms_machine_column %in% columns_to_exclude ) {
    metadata_column_selected <- setdiff( metadata_column_selected, ms_machine_column)
  }

  cln_meatadata_tbl <- metadata_tbl |>
    column_to_rownames(as_name( enquo(sample_id_column ))) |>
    dplyr::select( all_of( c(metadata_column_selected) ) )

  colour_rules <- getCategoricalColourRules( metadata_tbl =  cln_meatadata_tbl
                                             , metadata_column_labels = metadata_column_labels
                                             , metadata_column_selected = metadata_column_selected
                                             , categorical_columns = categorical_columns
                                             , ms_machine_column = ms_machine_column
                                             , columns_to_exclude = columns_to_exclude
                                             , na_colour = na_colour)

  print("Add column annotation")
  colnames(cln_meatadata_tbl) <-  metadata_column_labels_copy[metadata_column_selected]


  continous_colour_list <- getContinousColourRules( metadata_tbl
                                                    , metadata_column_labels
                                                    , metadata_column_selected
                                                    , continous_scale_columns
                                                    , na_colour = na_colour)

  categorical_and_continuous_colour_rules <- c( colour_rules, continous_colour_list)

  columns_to_use <- setdiff(names(categorical_and_continuous_colour_rules), metadata_column_labels_copy[columns_to_exclude])
  categorical_and_continuous_colour_rules_filt <- categorical_and_continuous_colour_rules[columns_to_use]

  return(categorical_and_continuous_colour_rules_filt)
}



#' @title Convert a Continuous Variable to Categorical
#' @description Bins a continuous variable into a specified number of categories.
#'
#' @param metadata_tbl The metadata table.
#' @param column_name The name of the continuous column to convert, as a string.
#' @param num_colours The number of bins (categories) to create.
#'
#' @return A factor representing the binned variable.
#' @export
changeToCategorical <- function(metadata_tbl, column_name, num_colours=9) {

  list_of_values <-  metadata_tbl |>
    dplyr::select( all_of(column_name)  ) |>
    #dplyr::filter(!is.na(!!sym(column_name))) |>
    pull()

  min_value <- min(list_of_values, na.rm =TRUE)
  max_value <-  max(list_of_values, na.rm =TRUE)

  if(min_value > 1) {
    min_value <- floor(min_value)
  }

  if(max_value > 1) {
    max_value <- ceiling(max_value)
  }

  formatted_list_of_values <- cut(list_of_values, breaks=seq( min_value, max_value, length.out=num_colours) )

  formatted_list_of_values
}

#------------------------------------------------------------------------------------------------

#' @title Generate a Sample Correlation Heatmap
#' @description Creates a complex heatmap of sample-sample correlations, with annotations
#' for sample metadata.
#'
#' @param correlation_matrix A sample-sample correlation matrix.
#' @param metadata_tbl The sample metadata table.
#' @param is_HEK_column The unquoted column in `metadata_tbl` indicating control samples.
#' @param metadata_column_labels A named vector for relabeling columns.
#' @param metadata_column_selected A character vector of columns to use for annotation.
#' @param colour_rules A named list of color palettes for annotations.
#' @param columns_to_exclude A character vector of annotation columns to exclude.
#' @param sample_id_column The unquoted column name for sample IDs.
#' @param use_raster A boolean indicating whether to render the heatmap as a raster image.
#' @param raster_device The graphics device for rasterization (e.g., "CairoPDF").
#' @param heatmap_legend_param A list of parameters for the heatmap legend.
#' @param heatmap_width The width of the heatmap.
#' @param heatmap_height The height of the heatmap.
#'
#' @return A list containing the `Heatmap` object and a list of `Legend` objects.
#' @export
getSamplesCorrelationHeatMap <- function(correlation_matrix
                                         , metadata_tbl
                                         , is_HEK_column = is_HEK
                                         , metadata_column_labels
                                         , metadata_column_selected
                                         , colour_rules
                                         , columns_to_exclude
                                         , sample_id_column = Run
                                         , use_raster = TRUE
                                         , raster_device = "CairoPDF"
                                         , heatmap_legend_param = list(title = "Correlation")
                                         , heatmap_width = ncol(correlation_matrix)*unit(0.05, "cm")
                                         , heatmap_height = nrow(correlation_matrix)*unit(0.05, "cm")
) {

  names( metadata_column_labels) <- metadata_column_selected

  without_hek_samples <- metadata_tbl |>
    dplyr::filter( {{is_HEK_column}} == FALSE) |>
    pull({{sample_id_column}})

  correlation_samples_to_use <- intersect( colnames(correlation_matrix), without_hek_samples) |> sort()

  cln_meatadata_orig_col_name <- metadata_tbl |>
    dplyr::filter( {{is_HEK_column}} == FALSE) |>
    dplyr::filter( {{sample_id_column}} %in% correlation_samples_to_use) |>
    arrange( {{sample_id_column}}) |>
    column_to_rownames( as_name( enquo( sample_id_column)  ))|>
    dplyr::select( all_of( setdiff( metadata_column_selected, columns_to_exclude) ) )

  columns_to_use <- setdiff(names(colour_rules), metadata_column_labels[columns_to_exclude])
  colour_rules_filt <- colour_rules[columns_to_use]

  print("Add column annotation")
  cln_meatadata_tbl <- cln_meatadata_orig_col_name
  colnames(cln_meatadata_tbl) <-  metadata_column_labels[colnames(cln_meatadata_orig_col_name)]

  top_annotation <- HeatmapAnnotation(df = cln_meatadata_tbl |>
                                        dplyr::select( - any_of(metadata_column_labels[columns_to_exclude] ))
                                      , col = colour_rules_filt
                                      , show_legend = FALSE )

  print("Add row annotation")
  row_ha <- rowAnnotation( df = cln_meatadata_tbl |>
                             dplyr::select( - any_of(metadata_column_labels[columns_to_exclude] ))
                           , col = colour_rules_filt
                           , show_legend = FALSE)

  output_heatmap <- Heatmap(correlation_matrix[correlation_samples_to_use, correlation_samples_to_use]
                            , name="Correlation"
                            , left_annotation = row_ha
                            , top_annotation = top_annotation
                            , show_row_names=FALSE
                            , show_column_names = FALSE
                            , use_raster = use_raster
                            , raster_device = raster_device
                            , row_title_gp = gpar(fontsize = 13.2)
                            , column_title_gp = gpar(fontsize = 13.2)
                            , row_names_gp = gpar(fontsize = 12, fontfamily = "sans")
                            , column_names_gp = gpar(fontsize = 12)
                            , heatmap_legend_param = heatmap_legend_param
                            , heatmap_height = heatmap_height
                            , heatmap_width = heatmap_width
  )

  output_legends <- purrr::map2 (colour_rules_filt
                                 , names(colour_rules_filt)
                                 , \(rule, title) { Legend(labels = names(rule), title = title,
                                                           legend_gp = gpar(fill = rule))} )

  # output_legends <- packLegend(list = list_of_legends)

  return( list( heatmap = output_heatmap
                , legend = output_legends ))


}


#' @title Calculate Heatmap Size
#' @description Calculates the final width and height of a `ComplexHeatmap` object
#' after it has been drawn.
#'
#' @param ht A `Heatmap` or `HeatmapList` object.
#' @param unit The unit for the output dimensions (e.g., "inch").
#'
#' @return A numeric vector containing the width and height.
#' @export
calcHtSize = function(ht, unit = "inch") {
  pdf(NULL)
  ht = draw(ht)
  w = ComplexHeatmap:::width(ht)
  w = convertX(w, unit, valueOnly = TRUE)
  h = ComplexHeatmap:::height(ht)
  h = convertY(h, unit, valueOnly = TRUE)
  dev.off()

  c(w, h)
}

#------------------------------------------------------------------------------------------------

#' @title Pivot Peptide Intensity Matrix to Long Format
#' @description Converts a wide-format peptide intensity matrix (features x samples)
#' to a long-format data frame.
#'
#' @param input_matrix The wide-format numeric matrix.
#' @param sample_id_column The desired name for the sample ID column in the output table.
#' @param sequence_column The desired name for the peptide sequence column.
#' @param protein_id_column The name of the protein ID column (present in rownames).
#' @param quantity_column The desired name for the abundance column.
#' @param unlog_data A boolean indicating whether to apply a 2^x transformation to the data.
#'
#' @return A long-format data frame.
#' @export
peptidesIntensityMatrixPivotLonger <- function( input_matrix
                                                , sample_id_column
                                                , sequence_column
                                                , protein_id_column
                                                , quantity_column
                                                , unlog_data = TRUE) {

  output_matrix <- input_matrix |>
    as.data.frame() |>
    rownames_to_column(protein_id_column) |>
    pivot_longer(cols=!contains(protein_id_column)
                 , names_to = sample_id_column
                 , values_to =  quantity_column ) |>
    separate( col = protein_id_column
              , into=c(protein_id_column, "Stripped.Sequence"), sep="_")

  if ( unlog_data == TRUE) {
    output_matrix <- output_matrix|>
      mutate( {{quantity_column}} := 2^(!!sym(quantity_column)))

  }

  output_matrix
}

#------------------------------------------------------------------------------------------------

#' @title Pivot Protein Intensity Matrix to Long Format
#' @description Converts a wide-format protein intensity matrix (proteins x samples)
#' to a long-format data frame.
#'
#' @param input_matrix The wide-format numeric matrix.
#' @param sample_id_column The desired name for the sample ID column.
#' @param protein_id_column The name of the protein ID column (present in rownames).
#' @param quantity_column The desired name for the abundance column.
#'
#' @return A long-format data frame.
#' @export
proteinIntensityMatrixPivotLonger <- function( input_matrix
                                               , sample_id_column
                                               , protein_id_column
                                               , quantity_column) {

  output_matrix <- input_matrix |>
    as.data.frame() |>
    rownames_to_column(protein_id_column) |>
    pivot_longer(cols=!contains(protein_id_column)
                 , names_to = sample_id_column
                 , values_to =  quantity_column )

  output_matrix
}


#------------------------------------------------------------------------------------------------

#' @title Helper to Remove Proteins with Only One Replicate
#' @description A helper function that removes proteins observed in only one replicate
#' across all sample groups. A protein must be seen in >1 replicate in at least
#' >1 group to be retained.
#'
#' @param input_table The input data frame of protein-level data.
#' @param samples_id_tbl A data frame with sample IDs and replicate group information.
#' @param input_table_sample_id_column The unquoted column name for sample IDs in `input_table`.
#' @param sample_id_tbl_sample_id_column The unquoted column name for sample IDs in `samples_id_tbl`.
#' @param replicate_group_column The unquoted column name for replicate groups.
#' @param protein_id_column The unquoted column name for protein identifiers.
#' @param quantity_column The unquoted column name for abundance values.
#' @param core_utilisation The number of cores for parallel processing.
#'
#' @return A data frame with the filtered proteins.
#' @export
removeProteinsWithOnlyOneReplicateHelper <- function(input_table
                                               , samples_id_tbl
                                               , input_table_sample_id_column = Run
                                               , sample_id_tbl_sample_id_column  =  ms_filename
                                               , replicate_group_column = general_sample_info
                                               , protein_id_column = Protein.Ids
                                               , quantity_column = Protein.Normalised
                                               , core_utilisation ) {

  # Count the number of technical replicates per sample and peptide combination
  num_tech_reps_per_sample_and_protein <- NA
  if( length(which(is.na(core_utilisation))) == 0 ) {
    num_tech_reps_per_sample_and_protein <- input_table |>
      left_join( samples_id_tbl, by=join_by( {{input_table_sample_id_column}} == {{sample_id_tbl_sample_id_column}} ) ) |>
      dplyr::filter( !is.na( {{quantity_column}}))  |>
      group_by( {{replicate_group_column}}, {{protein_id_column}} ) |>
      #partition(core_utilisation) |>
      summarise(counts = n() ) |>
      #collect() |>
      ungroup()
  } else {
    num_tech_reps_per_sample_and_protein <- input_table |>
      left_join( samples_id_tbl, by=join_by( {{input_table_sample_id_column}} == {{sample_id_tbl_sample_id_column}} ) ) |>
      dplyr::filter( !is.na( {{quantity_column}}))  |>
      group_by( {{replicate_group_column}}, {{protein_id_column}}) |>
      partition(core_utilisation) |>
      summarise(counts = n() ) |>
      collect() |>
      ungroup()
  }

  ## Need to have two or more replicates in at least two groups to be included
  proteins_in_two_or_more_groups_with_two_or_more_replicates <- num_tech_reps_per_sample_and_protein |>
    dplyr::filter(counts > 1) |>
    group_by({ { protein_id_column } }) |>
    summarise(num_groups = n()) |>
    ungroup() |>
    dplyr::filter(num_groups > 1) |>
    dplyr::select(-num_groups) |>
    distinct()


  removed_proteins_with_only_one_replicate <- input_table |>
    inner_join(proteins_in_two_or_more_groups_with_two_or_more_replicates
               , by = join_by({ { protein_id_column } }))  |>
    distinct()

  removed_proteins_with_only_one_replicate
}


##-----------------------------------------------------------------------------------------

#' @title Impute Missing Protein-Level Values
#' @description Performs missing value imputation for protein-level data. It imputes
#' missing values within a replicate group using the average of observed values,
#' but only if more than one replicate has a value.
#'
#' @param input_table The input data frame of protein-level data.
#' @param metadata_table A data frame with sample metadata.
#' @param input_table_sample_id_column The unquoted column name for sample IDs in `input_table`.
#' @param sample_id_tbl_sample_id_column The unquoted column name for sample IDs in `metadata_table`.
#' @param replicate_group_column The unquoted column name for replicate groups.
#' @param protein_id_column The unquoted column name for protein identifiers.
#' @param quantity_to_impute_column The unquoted column name for abundance values to impute.
#' @param imputed_value_column The unquoted column name for the new column with imputed values.
#' @param hek_string A character string to identify control samples to exclude.
#' @param core_utilisation The number of cores for parallel processing.
#'
#' @return A data frame with missing values imputed.
#' @export
proteinMissingValueImputation <- function( input_table
                                           , metadata_table
                                           , input_table_sample_id_column = Run
                                           , sample_id_tbl_sample_id_column  =  ms_filename
                                           , replicate_group_column = general_sample_info
                                           , protein_id_column = Protein.Ids
                                           , quantity_to_impute_column = Protein.Normalised
                                           , imputed_value_column = Protein.Imputed
                                           , hek_string = "HEK"
                                           , core_utilisation ) {

  # Max number of technical replicates
  num_tech_rep_per_sample <-  metadata_table  |>
    dplyr::filter( !str_detect( {{replicate_group_column}}, hek_string))  |>
    distinct( {{sample_id_tbl_sample_id_column}}, {{replicate_group_column}}) |>
    group_by( {{replicate_group_column}}) |>
    summarise(total_num_tech_rep = n()) |>
    ungroup()

  # Count the number of technical replicates per sample and protein combination
  num_tech_reps_per_sample_and_protein <- NA

  if( length(which(is.na(core_utilisation))) == 0 ) {
    num_tech_reps_per_sample_and_protein <- input_table |>
      left_join( metadata_table
                 , by=join_by( {{input_table_sample_id_column}} == {{sample_id_tbl_sample_id_column}} ) ) |>
      dplyr::filter( !str_detect( {{replicate_group_column}}, hek_string))  |>
      dplyr::filter( !is.na( {{quantity_to_impute_column}}))  |>
      distinct( {{replicate_group_column}}, {{protein_id_column}}, {{quantity_to_impute_column}}) |>
      group_by( {{replicate_group_column}}, {{protein_id_column}} ) |>
      #partition(core_utilisation) |>
      summarise( num_tech_rep = n()
                 , average_value = mean({{quantity_to_impute_column}}, na.rm=TRUE )) |>
      #collect() |>
      ungroup()
  } else {
    num_tech_reps_per_sample_and_protein <- input_table |>
      left_join( metadata_table
                 , by=join_by( {{input_table_sample_id_column}} == {{sample_id_tbl_sample_id_column}} ) ) |>
      dplyr::filter( !str_detect( {{replicate_group_column}}, hek_string))  |>
      dplyr::filter( !is.na( {{quantity_to_impute_column}}))  |>
      distinct( {{replicate_group_column}}, {{protein_id_column}}, {{quantity_to_impute_column}}) |>
      group_by( {{replicate_group_column}}, {{protein_id_column}}) |>
      partition(core_utilisation) |>
      summarise( num_tech_rep = n()
                 , average_value = mean({{quantity_to_impute_column}}, na.rm=TRUE)) |>
      collect() |>
      ungroup()

  }

  # total number of tech replicates > actual number technical replicates with data > 1
  rows_needing_imputation <-  num_tech_reps_per_sample_and_protein |>
    left_join( num_tech_rep_per_sample
               , by = join_by( {{replicate_group_column}} ) ) |>
    dplyr::filter( total_num_tech_rep > num_tech_rep &
                     num_tech_rep > 1)

  get_combinations_part_1 <- metadata_table |>
    distinct( {{sample_id_tbl_sample_id_column}}, {{replicate_group_column}} ) |>
    left_join(  input_table |>
                  distinct( {{input_table_sample_id_column}}, {{protein_id_column}} )
                , by =join_by( {{sample_id_tbl_sample_id_column}} == {{input_table_sample_id_column}}) )

  all_proteins_combination <- metadata_table |>
    distinct( {{sample_id_tbl_sample_id_column}}, {{replicate_group_column}}) |>
    group_by({{replicate_group_column}} ) |>
    nest( data = {{sample_id_tbl_sample_id_column}} )  |>
    left_join( get_combinations_part_1 |>
                 dplyr::select( -{{sample_id_tbl_sample_id_column}}) |>
                 dplyr::distinct( {{replicate_group_column}}, {{protein_id_column}})
               , by = join_by( {{replicate_group_column}}))  |>
    unnest( data ) |>
    ungroup({{replicate_group_column}})


  make_imputation <- all_proteins_combination |>
    left_join( input_table
               , by = join_by( {{sample_id_tbl_sample_id_column}} == {{input_table_sample_id_column}}
                               , {{protein_id_column}} == {{protein_id_column}} ) ) |>
    left_join(rows_needing_imputation
              , by = join_by( {{replicate_group_column}}
                              , {{protein_id_column}} ))  |>
    dplyr::filter(!is.na({{protein_id_column}})  ) |>
    mutate( is_imputed = case_when (is.na({{quantity_to_impute_column}})
                                    & !is.na(average_value)  ~ TRUE
                                    , TRUE ~ FALSE) ) |>
    mutate ( {{imputed_value_column}} := case_when (is.na({{quantity_to_impute_column}})
                                                    & !is.na(average_value)  ~ average_value
                                                    , TRUE ~ {{quantity_to_impute_column}} ) ) |>
    dplyr::select( -num_tech_rep
                   , - average_value
                   , - total_num_tech_rep
                   , - {{replicate_group_column}} ) |>
    dplyr::rename( {{input_table_sample_id_column}} := {{sample_id_tbl_sample_id_column}})

  make_imputation
}


#------------------------------------------------------------------------------------------------

#' @title Average Protein Intensity Across Replicates
#' @description Calculates the average protein intensity across technical replicates
#' for each protein within each sample group.
#'
#' @param input_table The input data frame of protein-level data.
#' @param metadata_table A data frame with sample metadata.
#' @param protein_id_column The unquoted column name for protein identifiers.
#' @param input_table_sample_id_column The unquoted column name for sample IDs in `input_table`.
#' @param sample_id_tbl_sample_id_column The unquoted column name for sample IDs in `metadata_table`.
#' @param replicate_group_column The unquoted column name for replicate groups.
#' @param quantity_column The unquoted column name for abundance values.
#' @param avg_quantity_column The desired name for the new column with averaged abundances.
#'
#' @return A data frame with averaged protein intensities.
#' @export
avgReplicateProteinIntensity <- function( input_table
                                          , metadata_table
                                          , protein_id_column = protein_id_column
                                          , input_table_sample_id_column = Run
                                          , sample_id_tbl_sample_id_column  =  Run
                                          , replicate_group_column = collaborator_patient_id
                                          , quantity_column = Log2.Protein.Imputed
                                          , avg_quantity_column = Avg.Log2.Protein.Imputed) {

  avg_log2_protein_intensity_imputed <- input_table |>
    inner_join( metadata_table
                , by=join_by( {{input_table_sample_id_column}} == {{sample_id_tbl_sample_id_column}})) |>
    group_by( {{protein_id_column}},  {{replicate_group_column}} ) |>
    summarise ( {{avg_quantity_column}} := mean({{quantity_column}}, na.rm=TRUE))  |>
    ungroup()

  avg_log2_protein_intensity_imputed
}

#------------------------------------------------------------------------------------------------

#' @title Plot Density of Protein Intensity
#' @description Generates a density plot of protein intensities, faceted by whether
#' proteins were quantified from single or multiple peptides.
#'
#' @param protein_intensity_long_tbl A long-format data frame of protein intensities.
#' @param number_of_peptides_per_protein_per_sample A data frame mapping protein-sample pairs to the number of peptides.
#' @param protein_id_column The unquoted column name for protein identifiers.
#' @param sample_id_column The unquoted column name for sample identifiers.
#' @param num_peptides_column The unquoted column name for the number of peptides.
#' @param protein_intensity_column The unquoted column name for protein intensity values.
#'
#' @return A ggplot object.
#' @export
plotDensityOfProteinIntensityPerSample <- function( protein_intensity_long_tbl
                                                    , number_of_peptides_per_protein_per_sample
                                                    , protein_id_column = Protein.Ids
                                                    , sample_id_column = Run
                                                    , num_peptides_column = num_peptides_after_impute
                                                    , protein_intensity_column = Log2.Protein.Imputed) {

  protein_intensity_vs_num_peptides_for_replicates <- protein_intensity_long_tbl |>
    left_join( number_of_peptides_per_protein_per_sample
               , by = join_by( {{protein_id_column}}, {{sample_id_column}} )) |>
    mutate( peptides_status = ifelse( {{num_peptides_column}} == 1, "Multiple Peptides"
                                      , "Single Peptide")) |>
    ggplot(aes( {{protein_intensity_column}}, group=peptides_status, fill= peptides_status, alpha=0.5 )) +
    geom_density()  +
    scale_alpha(guide = 'none') +
    apafTheme()  +
    xlab("log2 Protein Intensity") +
    ylab("Density") +
    labs( fill = "Peptide") +
    scale_y_continuous( expand = expansion(  mult=c(0, 0.1)))

}



#------------------------------------------------------------------------------------------------

#' @title Plot Protein Quantification Across Samples
#' @description Creates a bar plot showing the number of proteins quantified in
#' different percentages of samples, faceted by single vs. multiple peptide evidence.
#'
#' @param protein_intensity_long_tbl A long-format data frame of protein intensities.
#' @param number_of_peptides_per_protein_per_sample A data frame mapping protein-sample pairs to the number of peptides.
#' @param protein_id_column The unquoted column name for protein identifiers.
#' @param sample_id_column The unquoted column name for sample identifiers.
#' @param num_peptides_column The unquoted column name for the number of peptides.
#' @param protein_intensity_column The unquoted column name for protein intensity values.
#'
#' @return A ggplot object.
#' @export
plotPercentSamplesVsProteinQuantified <- function ( protein_intensity_long_tbl = frozen_protein_table
                                                    , number_of_peptides_per_protein_per_sample = number_of_peptides_per_protein_per_sample
                                                    , protein_id_column = Protein.Ids
                                                    , sample_id_column = Run
                                                    , num_peptides_column = num_peptides_after_impute
                                                    , protein_intensity_column = Log2.Protein.Imputed) {

  samples_vs_intensity <-  protein_intensity_long_tbl  |>
    left_join( number_of_peptides_per_protein_per_sample
               , by = join_by(  {{protein_id_column}}, {{sample_id_column}} )) |>
    mutate( peptides_status = ifelse( num_peptides_after_impute == 1, "Multiple Peptides"
                                      , "Single Peptide"))

  total_num_samples <- samples_vs_intensity |>
    distinct( {{sample_id_column}}) |>
    nrow()

  summarise_peptide_status <- function ( input_vector) {

    if( "Multiple Peptides" %in% input_vector  ) {
      return ( "Multiple Peptides" )
    } else {
      return ( "Single Peptide")
    }
  }

  num_samples_per_protein <- samples_vs_intensity |>
    dplyr::filter(!is.na({{protein_intensity_column}})) |>
    group_by( {{protein_id_column}} ) |>
    summarise( num_values = n()
               , peptides_status = summarise_peptide_status(peptides_status )) |>
    ungroup()  |>
    mutate ( percentage = num_values / total_num_samples * 100 )

  num_samples_per_protein |>
    mutate( percentage_bin = cut( percentage, breaks = c(0, 20, 40, 60, 80, 100) )) |>
    ggplot( aes( percentage_bin, fill = peptides_status, group = peptides_status)) +
    geom_bar( position = "dodge" ) +
    apafTheme()  +
    xlab("Percentage of Samples") +
    ylab("Num. Quantified Proteins ") +
    labs( fill = "Peptide") +
    scale_y_continuous( expand = expansion(  mult=c(0, 0.1)))

}

#------------------------------------------------------------------------------------------------
# Codes to format experimental design table for pairwise comparison of groups
#' @title Create "Each vs. All" Dummy Columns in Design Matrix
#' @description Modifies a design matrix by adding new dummy-coded columns for each
#' level of a specified factor, facilitating "each vs. all" comparisons.
#'
#' @param input_table The design matrix data frame.
#' @param id_cols The unquoted column name for the unique sample/row identifiers.
#' @param column The unquoted column name of the factor to be dummy-coded.
#'
#' @return The modified design matrix with new dummy columns.
#' @export
cleanDesignMatrixCreateEachVersusAllColumns <- function(input_table, id_cols, column ) {

  id_col_name <-  as_string(as_name(enquo(id_cols)))
  column_string <- as_string(as_name(enquo(column)))

  new_columns_tab <-  input_table |>
    mutate( my_value = TRUE ) |>
    pivot_wider( id_cols = {{id_cols}}
                 , names_from = {{column}}
                 , values_from = my_value
                 , values_fill = FALSE
                 , names_prefix = paste0(column_string, ".")  )

  new_column_names <- base::setdiff( colnames( new_columns_tab), id_col_name )

  # print (id_col_name)
  # print(new_column_names)

  return_table <- input_table |>
    left_join( new_columns_tab
               , by = join_by( {{id_cols}} )) |>
    relocate( all_of(new_column_names), .after={{column}} )

  return_table

}

#' @title Clean Category Names
#' @description Cleans special characters from category names to make them valid
#' for use as column names or in formulas.
#'
#' @param x A character string to be cleaned.
#'
#' @return The cleaned character string.
#' @export
cleanDesignMatrixCleanCategories <- function(x ) {

  str_replace_all(x, ">=", "ge") |>
    str_replace_all( "<=", "le") |>
    str_replace_all( ">", "gt") |>
    str_replace_all( "<", "lt") |>
    str_replace_all( "\\+", ".POS") |>
    str_replace_all( "\\-", ".NEG") |>
    str_replace_all( " ", "\\.") |>
    str_replace_all("&", "and") |>
    str_replace_all("/", "_") |>
    str_replace_all("\\:", ".")
}

# clean_categories("<= 5")
# clean_categories("> 3")
# clean_categories(">& /3")

#' @title Map Category Cleaning Function over a Column
#' @description Applies the `cleanDesignMatrixCleanCategories` function to each element
#' of a specified column in a data frame.
#'
#' @param input_table The data frame.
#' @param column The unquoted column name to be cleaned.
#'
#' @return The data frame with the specified column cleaned.
#' @export
cleanDesignMatrixCleanCategoriesMap <- function( input_table, column ) {
  input_table |>
    mutate( {{column}} := purrr::map_chr( {{column}}, cleanDesignMatrixCleanCategories ) )

}

#------------------------------------------------------------------------------------------------

#' @title Generate a Protein Intensity Heatmap
#' @description Creates a complex heatmap of protein intensities across samples, with
#' annotations for sample metadata.
#'
#' @param protein_matrix A wide-format matrix of protein intensities (proteins x samples).
#' @param metadata_tbl The sample metadata table.
#' @param is_HEK_column The unquoted column in `metadata_tbl` indicating control samples.
#' @param metadata_column_selected A character vector of columns to use for annotation.
#' @param metadata_column_labels A named vector for relabeling columns.
#' @param colour_rules A named list of color palettes for annotations.
#' @param columns_to_exclude A character vector of annotation columns to exclude.
#' @param core_utilisation_samples A boolean passed to `ComplexHeatmap` for parallelization.
#' @param sort_by_sample_id A boolean indicating whether to sort samples by ID.
#' @param sample_id_column The unquoted column name for sample IDs.
#' @param use_raster A boolean indicating whether to render the heatmap as a raster image.
#' @param raster_device The graphics device for rasterization.
#' @param heatmap_legend_param A list of parameters for the heatmap legend.
#'
#' @return A list containing the `Heatmap` object and a list of `Legend` objects.
#' @export
getProteinsHeatMap <- function( protein_matrix
                                , metadata_tbl
                                , is_HEK_column = is_HEK
                                , metadata_column_selected
                                , metadata_column_labels
                                # , categorical_columns
                                # , continous_scale_columns
                                # , ms_machine_column
                                , colour_rules
                                , columns_to_exclude
                                , core_utilisation_samples = TRUE
                                , sort_by_sample_id = TRUE
                                , sample_id_column = Run
                                , use_raster = TRUE
                                , raster_device = "CairoTIFF"
                                , heatmap_legend_param = list(title = "Intensity")) {

  metadata_column_labels_copy <- metadata_column_labels
  names( metadata_column_labels_copy) <- metadata_column_selected

  print("Without HEK samples")
  without_hek_samples <- metadata_tbl |>
    dplyr::filter( {{is_HEK_column}} == FALSE) |>
    pull({{sample_id_column}})

  samples_to_use <- intersect( colnames(protein_matrix), without_hek_samples)

  if( sort_by_sample_id == TRUE) {
    samples_to_use <- samples_to_use |>
      sort()
  }

  cln_meatadata_tbl_orig_col_names <- metadata_tbl |>
    dplyr::filter( {{is_HEK_column}} == FALSE)  |>
    dplyr::filter( {{sample_id_column}} %in% samples_to_use) |>
    arrange( {{sample_id_column}}) |>
    dplyr::select( {{sample_id_column}}, all_of( setdiff(metadata_column_selected, columns_to_exclude) ) ) |>
    distinct()  |>
    column_to_rownames( as_label( enquo(sample_id_column) ) )

  print("Add column annotation")
  cln_meatadata_tbl <- cln_meatadata_tbl_orig_col_names[ samples_to_use
                                                         , setdiff(metadata_column_selected, columns_to_exclude)]
  colnames(cln_meatadata_tbl) <- metadata_column_labels_copy[setdiff(metadata_column_selected, columns_to_exclude)]
  colour_rules_filt <- colour_rules[metadata_column_labels_copy[setdiff(metadata_column_selected, columns_to_exclude)] ]

  # print(colour_rules)
  # print(metadata_column_labels_copy[setdiff(metadata_column_selected, columns_to_exclude)] )
  # print(colour_rules)
  # print(colour_rules_filt)

  # print(colour_rules_filt)

  print("Set top located annotation")
  top_annotation <- HeatmapAnnotation( df = cln_meatadata_tbl
                                       , col = colour_rules_filt
                                       , show_legend = FALSE
                                       , annotation_name_side = "left" )

  print("Print Heatmap")

  heatmap <- Heatmap( protein_matrix[, samples_to_use]
                      , name="Intensity"
                      , top_annotation = top_annotation
                      , show_row_names=FALSE
                      , show_column_names = FALSE
                      , use_raster = use_raster
                      , raster_device = raster_device
                      , row_title_gp = gpar(fontsize = 13.2)
                      , column_title_gp = gpar(fontsize = 13.2)
                      , row_names_gp = gpar(fontsize = 12, fontfamily = "sans")
                      , column_names_gp = gpar(fontsize = 12)
                      , heatmap_legend_param = heatmap_legend_param
                      , core_utilisation_columns = core_utilisation_samples )

  output_legends <- purrr::map2 (colour_rules_filt
                                 , names(colour_rules_filt)
                                 , \(rule, title) { Legend(labels = names(rule), title = title,
                                                           legend_gp = gpar(fill = rule))} )

  return(list( heatmap = heatmap
               , legend = output_legends ))

}


#------------------------------------------------------------------------------------------------

#' @title Calculate Percent Missing Values per Protein Category
#' @description For each protein, calculates the percentage of missing values for
#' each level of each categorical variable in the experimental design.
#'
#' @param intensity_wide_table A wide-format data frame of protein intensities.
#' @param protein_id The name of the protein ID column as a string.
#' @param pattern A tidyselect pattern to select the sample columns.
#' @param experimental_design_table The experimental design data frame.
#' @param names_to The desired name for the sample ID column after pivoting.
#' @param values_to The desired name for the abundance column after pivoting.
#' @param is_missing_column The unquoted name for the new boolean column indicating missingness.
#'
#' @return A data frame summarizing the missing value percentages.
#' @export
calculatePercentMissingPerProtein <- function( intensity_wide_table
                                               , protein_id = "uniprot_acc"
                                               , pattern = ! tidyselect::matches( protein_id )
                                               , experimental_design_table
                                               , names_to = "sample_collaborator_sample_id"
                                               , values_to = "Avg.Log2.Protein.Imputed"
                                               , is_missing_column = is_missing ) {

  # print(deparse1(substitute(!!sym({{protein_id}}) )) )

  intensity_long_table <- intensity_wide_table |>
    pivot_longer( cols= {{pattern}}
                  , names_to = names_to
                  , values_to = values_to )


  intensity_vs_design_matrix <- intensity_long_table |>
    mutate( {{names_to}}:= purrr::map_chr(!!rlang::sym(names_to), as.character) ) |>
    left_join( experimental_design_table
               , by = join_by({{names_to}}))


  list_of_columns_to_pivot <- setdiff( colnames( experimental_design_table)
                                       , c( protein_id
                                            , names_to
                                            , values_to
                                            , as_string(as_name(enquo(is_missing)))  ))

  intensity_vs_design_matrix_cln <- intensity_vs_design_matrix |>
    mutate( {{is_missing_column}} := case_when( is.nan( !!sym(values_to)) |
                                                  is.na(Avg.Log2.Protein.Imputed)  ~ TRUE
                                                , TRUE ~ FALSE)) |>
    relocate({{is_missing_column}}, .after=!!sym(values_to)  ) |>
    pivot_longer( cols = all_of(list_of_columns_to_pivot)
                  , names_to = "parameter_name"
                  , values_to = "values" )

  missing_value_per_category <- intensity_vs_design_matrix_cln |>
    group_by( !!sym( protein_id)
              , parameter_name
              , values ) |>
    summarise( num_values = n()
               , num_missing =  sum( is_missing)  ) |>
    ungroup( ) |>
    mutate( perc_missing = num_missing/num_values * 100 ) |>
    mutate( num_present = num_values - num_missing  ) |>
    mutate( perc_present = 100  - perc_missing ) |>
    dplyr::select( uniprot_acc
                   , parameter_name
                   , values
                   , num_missing
                   , num_present
                   , num_values
                   , perc_missing
                   , perc_present ) |>
    mutate( compare_column = paste0( parameter_name, as.character(values)))

  missing_value_per_category
}

#' @title Perform Fisher's Exact Test for Missing Values
#' @description For each contrast, performs a Fisher's exact test to determine if the
#' proportion of missing values for each protein is significantly different
#' between the two groups in the contrast.
#'
#' @param contrasts_table A data frame defining the contrasts to be tested.
#' @param missing_value_per_category A data frame from `calculatePercentMissingPerProtein`.
#'
#' @return A data frame with the results of the Fisher's exact tests, including p-values and FDR.
#' @export
calculateMissingValuesPerProteinFishersTest <- function( contrasts_table, missing_value_per_category) {

  contrasts_table_separated <- contrasts_table |>
    separate( col=contrasts, sep = "[=-]", into=c("contrast_name", "left", "right"))

  runFisherTest <- function( a1, b1, a2, b2) {
    fisher.test( matrix( c( a1, b1, a2, b2), 2, 2, byrow = TRUE))$p.value
  }

  plan(multisession, workers = 8)


  contasts_missing_counts_tbl <- contrasts_table_separated |>
    left_join( missing_value_per_category
               , by=join_by( left == compare_column)  ) |>
    left_join( missing_value_per_category
               , by=join_by( right == compare_column
                             , uniprot_acc == uniprot_acc )
               , suffix = c(".left", ".right")) |>
    dplyr::filter( !( is.na(num_missing.left)
                      & is.na(num_present.left)
                      & is.na(num_missing.right)
                      & is.na(num_present.right ))) |>
    mutate( fisher_test = furrr::future_pmap_dbl ( list( a1 = num_missing.left
                                                         , b1 = num_present.left
                                                         , a2 = num_missing.right
                                                         , b2 = num_present.right )
                                                   , \(a1,a2,b1, b2){ runFisherTest( a1=a1, b1=b1, a2=a2, b2=b2)} )    )

  # fisher.test(matrix( c(8, 19, 27, 73), 2,2, byrow=TRUE))
  # fisher.test(matrix( c(8, 19, 27, 73), 2,2, byrow=FALSE))

  contasts_missing_fdr_tbl <- contasts_missing_counts_tbl |>
    nest(.by=contrast_name, .key="tables" ) |>
    dplyr::mutate( updated_tables = purrr::map(tables
                                               , \(x){ x |> bind_cols( data.frame(fdr=p.adjust(x$fisher_test, method= "fdr"))) }) ) |>
    dplyr::select(-tables) |>
    unnest( updated_tables )

  contasts_missing_fdr_tbl

}

#------------------------------------------------------------------------------------------------


#' @title APAF ggplot2 Theme
#' @description A custom ggplot2 theme for APAF (Australian Proteome Analysis Facility)
#' style plots.
#'
#' @return A ggplot2 theme object.
#' @export
apafTheme <- function() {
  theme(
    # Set font family and size
    text = element_text(family = "Arial", size = 12),
    # Add rectangular box around the plot
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
    # Add grid lines
    panel.grid.major = element_line(color = "gray", linetype = "dashed"),
    panel.grid.minor = element_blank(),
    # Set plot background color
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA),
    # Set axis line and tick colors
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black"),
    # Set axis label colors and sizes
    axis.text.x = element_text(color = "black", size = 10),
    axis.text.y = element_text(color = "black", size = 10),
    # Set legend title and label colors and sizes
    legend.title = element_text(color = "black", size = 12),
    legend.text = element_text(color = "black", size = 10),
    # Set plot title and subtitle colors and sizes
    plot.title = element_text(color = "black", size = 14),
    plot.subtitle = element_text(color = "black", size = 12),
    # Set plot margin sizes
    plot.margin = unit(c(1, 1, 1, 1), "cm")

  )
}

##################################################################################################################

#' FilteringProgress Class
#' 
#' @description
#' An S4 class to track and store the progress of protein filtering steps in
#' proteomics data analysis. This class maintains records of protein and peptide
#' counts at each filtering stage.
#' 
#' @slot steps Character vector storing names of filtering steps
#' @slot proteins Numeric vector storing protein counts for each step
#' @slot total_peptides Numeric vector storing total peptide counts for each step
#' @slot peptides_per_protein List storing peptides per protein distributions for each step
#' @slot proteins_per_run List storing proteins per run counts for each step
#' @slot peptides_per_run List storing peptides per run counts for each step
#' 
#' @export
setClass("FilteringProgress",
  slots = list(
    steps = "character",
    proteins = "numeric",
    total_peptides = "numeric",
    peptides_per_protein = "list",
    proteins_per_run = "list",
    peptides_per_run = "list"
  )
)

##################################################################################################################

#' Initialize a new FilteringProgress object
#' 
#' @description
#' Creates a new FilteringProgress object with empty slots to track protein
#' filtering progress.
#' 
#' @return A new FilteringProgress object
#' 
#' @examples
#' filtering_progress <- new("FilteringProgress",
#'   steps = character(),
#'   proteins = numeric(),
#'   total_peptides = numeric(),
#'   peptides_per_protein = list(),
#'   proteins_per_run = list(),
#'   peptides_per_run = list()
#' )
#' 
#' @export
filtering_progress <- new("FilteringProgress",
  steps = character(),
  proteins = numeric(),
  total_peptides = numeric(),
  peptides_per_protein = list(),
  proteins_per_run = list(),
  peptides_per_run = list()
)

##################################################################################################################

#' Generate a color palette
#' 
#' @param n Number of colors needed
#' @param base_color Base color to use
#' @return Vector of colors
#' @export 
get_color_palette <- function(n, base_color) {
  colorRampPalette(c(base_color, "black"))(n)
}


#' Count unique proteins in peptide or protein data
#' 
#' @param data A data frame or S4 object containing protein/peptide data
#' @return Integer count of unique proteins
#' @export
countUniqueProteins <- function(data) {
  if (isS4(data)) {
    if ("peptide_data" %in% slotNames(data)) {
      return(data@peptide_data |> 
             distinct(Protein.Ids) |> 
             nrow())
    }
    if ("protein_quant_table" %in% slotNames(data)) {
      return(nrow(data@protein_quant_table))
    }
  }
  
  # For regular dataframes
  if ("Protein.Ids" %in% names(data)) {
    # Check if it's a protein quantification table
    if (all(sapply(data[setdiff(names(data), "Protein.Ids")], is.numeric))) {
      return(nrow(data))  # Each row is a unique protein
    }
    return(distinct(data, Protein.Ids) |> nrow())
  }
  stop("No Protein.Ids column found")
}

#' Count proteins per run
#' 
#' @param data A data frame or S4 object containing protein/peptide data
#' @return Data frame with run IDs and protein counts
#' @export
countProteinsPerRun <- function(data) {
  if (isS4(data)) {
    if ("peptide_data" %in% slotNames(data)) {
      return(data@peptide_data |>
             group_by(Run) |>
             summarise(n_proteins = n_distinct(Protein.Ids), 
                      .groups = "drop") |>
             arrange(Run))
    }
    if ("protein_quant_table" %in% slotNames(data)) {
      data <- data@protein_quant_table
      run_cols <- setdiff(names(data), "Protein.Ids")
      
      # For each run (column), count non-NA values
      result <- data.frame(
        Run = run_cols,
        n_proteins = sapply(run_cols, function(col) {
          sum(!is.na(data[[col]]))
        })
      ) |> arrange(Run)
      
      return(result)
    }
  }
  
  # For regular dataframes
  if ("Protein.Ids" %in% names(data)) {
    # Check if it's a protein quantification table
    if (all(sapply(data[setdiff(names(data), "Protein.Ids")], is.numeric))) {
      run_cols <- setdiff(names(data), "Protein.Ids")
      
      # For each run (column), count non-NA values
      result <- data.frame(
        Run = run_cols,
        n_proteins = sapply(run_cols, function(col) {
          sum(!is.na(data[[col]]))
        })
      ) |> arrange(Run)
      
      return(result)
    }
    
    # For peptide data
    if ("Run" %in% names(data)) {
      return(data |>
             group_by(Run) |>
             summarise(n_proteins = n_distinct(Protein.Ids), 
                      .groups = "drop") |>
             arrange(Run))
    }
  }
  stop("Required columns not found")
}

#' Calculate total unique peptides
#' 
#' @param data A data frame or S4 object containing protein/peptide data
#' @return Integer count of unique peptide-protein combinations
#' @export
calcTotalPeptides <- function(data) {
  # For protein quantification data, return NA
  if (isS4(data)) {
    if ("protein_quant_table" %in% slotNames(data)) {
      return(NA_integer_)
    }
    if ("peptide_data" %in% slotNames(data)) {
      return(data@peptide_data |>
             distinct(Protein.Ids, Stripped.Sequence) |>
             nrow())
    }
  }
  
  # For regular dataframes, check if it's protein quantification data
  if ("Protein.Ids" %in% names(data)) {
    if (all(sapply(data[setdiff(names(data), "Protein.Ids")], is.numeric))) {
      return(NA_integer_)
    }
    
    if ("Stripped.Sequence" %in% names(data)) {
      return(distinct(data, Protein.Ids, Stripped.Sequence) |> nrow())
    }
  }
  stop("Required columns not found")
}

#' Calculate peptides per protein
#' 
#' @param data A data frame or S4 object containing protein/peptide data
#' @return Data frame with protein IDs and peptide counts
#' @export
calcPeptidesPerProtein <- function(data) {
  # For protein quantification data, return empty data frame
  if (isS4(data)) {
    if ("protein_quant_table" %in% slotNames(data)) {
      return(data.frame(Protein.Ids = character(), 
                       n_peptides = integer()))
    }
    if ("peptide_data" %in% slotNames(data)) {
      return(data@peptide_data |>
             group_by(Protein.Ids) |>
             summarise(n_peptides = n_distinct(Stripped.Sequence), 
                      .groups = "drop"))
    }
  }
  
  # For regular dataframes, check if it's protein quantification data
  if ("Protein.Ids" %in% names(data)) {
    if (all(sapply(data[setdiff(names(data), "Protein.Ids")], is.numeric))) {
      return(data.frame(Protein.Ids = character(), 
                       n_peptides = integer()))
    }
    
    if ("Stripped.Sequence" %in% names(data)) {
      return(data |>
             group_by(Protein.Ids) |>
             summarise(n_peptides = n_distinct(Stripped.Sequence), 
                      .groups = "drop"))
    }
  }
  stop("Required columns not found")
}

#' Count peptides per run
#' 
#' @param data A data frame or S4 object containing protein/peptide data
#' @return Data frame with run IDs and peptide counts
#' @export
countPeptidesPerRun <- function(data) {
  # For protein quantification data, return empty data frame
  if (isS4(data)) {
    if ("protein_quant_table" %in% slotNames(data)) {
      return(data.frame(Run = character(), 
                       n_peptides = integer()))
    }
    if ("peptide_data" %in% slotNames(data)) {
      return(data@peptide_data |>
             group_by(Run) |>
             summarise(n_peptides = n_distinct(Stripped.Sequence), 
                      .groups = "drop") |>
             arrange(Run))
    }
  }
  
  # For regular dataframes, check if it's protein quantification data
  if ("Protein.Ids" %in% names(data)) {
    if (all(sapply(data[setdiff(names(data), "Protein.Ids")], is.numeric))) {
      return(data.frame(Run = character(), 
                       n_peptides = integer()))
    }
    
    if (all(c("Run", "Stripped.Sequence") %in% names(data))) {
      return(data |>
             group_by(Run) |>
             summarise(n_peptides = n_distinct(Stripped.Sequence), 
                      .groups = "drop") |>
             arrange(Run))
    }
  }
  stop("Required columns not found")
}

#' @title Update and Visualize Filtering Progress
#' @description Tracks and visualizes the impact of filtering steps on peptide 
#'   and protein counts. Updates a global `FilteringProgress` object and optionally 
#'   saves plots summarizing the changes. Handles both peptide-level and 
#'   protein-level data inputs.
#' 
#' @details 
#' This function acts as a central hub for monitoring data reduction throughout 
#' a filtering workflow. It performs the following actions:
#' \itemize{
#'   \item Initializes or retrieves a global S4 object named `filtering_progress` 
#'     of class `FilteringProgress`.
#'   \item Calculates key metrics (total unique proteins, proteins per run, 
#'     total unique peptides, peptides per protein distribution, peptides per run) 
#'     based on the input `data`. Peptide metrics are only calculated or updated 
#'     if `data` is identified as peptide-level data. For protein-level data, 
#'     peptide metrics from the last peptide step (if any) are carried forward or 
#'     initialized as empty/NA.
#'   \item Adds or updates these metrics in the `filtering_progress` object 
#'     under the specified `step_name`.
#'   \item Generates summary plots using `ggplot2`:
#'     \itemize{
#'       \item Bar plot of total unique proteins per step.
#'       \item Bar plot of total unique peptides per step (or placeholder if only protein data).
#'       \item Box plot of peptides per protein distribution per step (or placeholder).
#'       \item Line plot of proteins per run across steps.
#'       \item Line plot of peptides per run across steps (or placeholder).
#'     }
#'   \item If `omic_type` and `experiment_label` are provided and valid paths can be 
#'     derived from the global `project_dirs` object, the generated plots are saved 
#'     as PNG files into the derived `time_dir`. Warnings are issued if paths cannot be 
#'     derived or if `project_dirs` is not found.
#'   \item If `return_grid` is `TRUE`, arranges the plots into a single grid using 
#'     `gridExtra` and returns the grid object (grob). Also saves this combined grid 
#'     if plot saving is enabled.
#'   \item If `return_grid` is `FALSE` (default), prints each plot individually 
#'     and returns the list of plot objects invisibly.
#' }
#' 
#' **Important:** This function relies on and modifies a global variable named 
#' `filtering_progress`. For saving plots, it depends on the global `project_dirs` 
#' object (expected to be populated by `setupDirectories()`) and the successful 
#' derivation of `time_dir` from it using `omic_type` and `experiment_label`.
#' 
#' @param data The input data object. Can be a data frame (expected to conform 
#'   to typical peptide or protein quantification structures) or an S4 object 
#'   containing relevant slots (e.g., inheriting from `SummarizedExperiment`). 
#'   The function attempts to automatically detect if it\'s peptide or protein data.
#' @param step_name A character string uniquely identifying the current filtering 
#'   step (e.g., "InitialData", "FilteredByQuality", "Normalized"). This name is 
#'   used for tracking in the `filtering_progress` object and plot labels.
#' @param omic_type Optional character string. The type of omics data 
#'   (e.g., "proteomics", "metabolomics"). Used with `experiment_label` to 
#'   derive save paths from the global `project_dirs` object. If `NULL` (default) 
#'   or `experiment_label` is `NULL`, plots are not saved.
#' @param experiment_label Optional character string. The specific experiment 
#'   label (e.g., "workshop_data"). Used with `omic_type` to derive save paths 
#'   from the global `project_dirs` object. If `NULL` (default) or `omic_type` 
#'   is `NULL`, plots are not saved.
#' @param overwrite Logical. If `TRUE`, allows overwriting an existing entry for 
#'   `step_name` in the `filtering_progress` object. If `FALSE` (default) and 
#'   `step_name` already exists, the function will stop with an error.
#' @param return_grid Logical. If `TRUE`, returns a single combined plot grid 
#'   object created with `gridExtra::grid.arrange()`. If `FALSE` (default), prints 
#'   individual plots and returns an invisible list of the ggplot objects.
#' 
#' @return If `return_grid` is `TRUE`, returns a `grob` object (a grid graphical object). 
#'   If `return_grid` is `FALSE`, returns an invisible list containing the individual 
#'   `ggplot` objects (`proteins_total`, `proteins_per_run`, `peptides_total`, 
#'   `peptides_per_protein`, `peptides_per_run`). Has side effects: modifies the 
#'   global `filtering_progress` object and potentially saves plots to disk if 
#'   `omic_type` and `experiment_label` are provided and paths are valid.
#'   
#' @importFrom ggplot2 ggplot aes geom_bar geom_text labs theme_minimal theme element_text panel_grid_major_x geom_line geom_point scale_color_manual annotate theme_void geom_boxplot coord_cartesian quantile ggsave
#' @importFrom dplyr bind_rows mutate group_by ungroup mean %>%
#' @importFrom forcats fct_reorder
#' @importFrom gridExtra arrangeGrob grid.arrange
#' @importFrom methods isS4 slotNames new
#' @importFrom stats quantile
#' 
#' @export
updateProteinFiltering <- function(data, step_name, 
                                 omic_type = NULL, experiment_label = NULL,
                                 overwrite = FALSE, return_grid = FALSE) {
    
    # Initialize filtering_progress if it doesn\'t exist
    if (!exists("filtering_progress", envir = .GlobalEnv)) {
        filtering_progress <- new("FilteringProgress")
        assign("filtering_progress", filtering_progress, envir = .GlobalEnv)
    }
    
    # Get the current filtering_progress object
    filtering_progress <- get("filtering_progress", envir = .GlobalEnv)
    
    # Path derivation and save_plots logic
    derived_time_dir <- NULL
    save_plots <- FALSE

    if (!is.null(omic_type) && !is.null(experiment_label)) {
        if (!exists("project_dirs", envir = .GlobalEnv)) {
            warning("Global object \'project_dirs\' not found. Plots will not be saved. Ensure \'setupDirectories()\' has been run.")
        } else {
            project_dirs_global <- get("project_dirs", envir = .GlobalEnv)
            omic_project_key <- paste0(omic_type, "_", experiment_label)

            if (!omic_project_key %in% names(project_dirs_global)) {
                warning(paste0("Entry for \'", omic_project_key, "\' not found in global \'project_dirs\'. Plots will not be saved."))
            } else {
                current_project_paths <- project_dirs_global[[omic_project_key]]
                if (is.null(current_project_paths)) {
                    warning(paste0("Entry for \'", omic_project_key, "\' in global \'project_dirs\' is NULL. Plots will not be saved."))
                } else {
                    derived_publication_graphs_dir <- current_project_paths$publication_graphs_dir
                    temp_time_dir <- current_project_paths$time_dir

                    if (is.null(temp_time_dir) || !is.character(temp_time_dir) || length(temp_time_dir) != 1 ||
                        is.null(derived_publication_graphs_dir) || !is.character(derived_publication_graphs_dir) || length(derived_publication_graphs_dir) != 1) {
                        warning(paste0("\'time_dir\' or \'publication_graphs_dir\' is missing, not a character string, or not a single path for \'", omic_project_key,
                                       "\' in global \'project_dirs\'. Plots will not be saved."))
                    } else {
                        if (!dir.exists(temp_time_dir)) {
                            warning(paste0("The derived \'time_dir\' (", temp_time_dir, ") for \'", omic_project_key,
                                           "\' does not exist. Plots will not be saved. Ensure directories are created via setupDirectories()."))
                        } else {
                            derived_time_dir <- temp_time_dir
                            save_plots <- TRUE
                            message(paste0("Plots will be saved to: ", derived_time_dir))
                        }
                    }
                }
            }
        }
    } else {
        # Message if omic_type/label are missing and saving might have been expected
        if (return_grid && (is.null(omic_type) || is.null(experiment_label))) {
             message("omic_type and/or experiment_label not provided. Plots will not be saved.")
        }
    }
    
    # Determine if we\'re working with protein_quant_table
    is_protein_quant <- if (isS4(data)) {
        "protein_quant_table" %in% slotNames(data)
    } else {
        # For data frames, check if it looks like a protein quant table
        if ("Protein.Ids" %in% names(data)) {
            all(sapply(data[setdiff(names(data), "Protein.Ids")], is.numeric))
        } else {
            FALSE
        }
    }
    
    # Calculate protein metrics (always done)
    protein_count <- countUniqueProteins(data)
    proteins_per_run <- countProteinsPerRun(data)
    
    # Ensure consistent data types in proteins_per_run
    proteins_per_run$Run <- as.character(proteins_per_run$Run)
    proteins_per_run$n_proteins <- as.numeric(proteins_per_run$n_proteins)
    
    # Update filtering progress based on data type
    if (step_name %in% filtering_progress@steps) {
        if (!overwrite) {
            stop("Step name \'", step_name, "\' already exists. Use overwrite = TRUE to replace it.")
        }
        idx <- which(filtering_progress@steps == step_name)
        
        # Always update protein metrics
        filtering_progress@proteins[idx] <- protein_count
        filtering_progress@proteins_per_run[[idx]] <- proteins_per_run
        
        if (!is_protein_quant) {
            # Update peptide metrics only for peptide data
            filtering_progress@total_peptides[idx] <- calcTotalPeptides(data)
            peptides_per_protein <- calcPeptidesPerProtein(data)
            peptides_per_run <- countPeptidesPerRun(data)
            
            # Ensure consistent data types
            peptides_per_protein$Protein.Ids <- as.character(peptides_per_protein$Protein.Ids)
            peptides_per_protein$n_peptides <- as.numeric(peptides_per_protein$n_peptides)
            
            peptides_per_run$Run <- as.character(peptides_per_run$Run)
            peptides_per_run$n_peptides <- as.numeric(peptides_per_run$n_peptides)
            
            filtering_progress@peptides_per_protein[[idx]] <- peptides_per_protein
            filtering_progress@peptides_per_run[[idx]] <- peptides_per_run
        }
    } else {
        filtering_progress@steps <- c(filtering_progress@steps, step_name)
        filtering_progress@proteins <- c(filtering_progress@proteins, protein_count)
        filtering_progress@proteins_per_run <- c(filtering_progress@proteins_per_run, 
                                               list(proteins_per_run))
        
        if (!is_protein_quant) {
            # Add peptide metrics only for peptide data
            filtering_progress@total_peptides <- c(filtering_progress@total_peptides, 
                                                 calcTotalPeptides(data))
            
            peptides_per_protein <- calcPeptidesPerProtein(data)
            peptides_per_run <- countPeptidesPerRun(data)
            
            # Ensure consistent data types
            peptides_per_protein$Protein.Ids <- as.character(peptides_per_protein$Protein.Ids)
            peptides_per_protein$n_peptides <- as.numeric(peptides_per_protein$n_peptides)
            
            peptides_per_run$Run <- as.character(peptides_per_run$Run)
            peptides_per_run$n_peptides <- as.numeric(peptides_per_run$n_peptides)
            
            filtering_progress@peptides_per_protein <- c(filtering_progress@peptides_per_protein, 
                                                       list(peptides_per_protein))
            filtering_progress@peptides_per_run <- c(filtering_progress@peptides_per_run, 
                                                   list(peptides_per_run))
        } else {
            # For protein data, maintain existing peptide metrics or add NA/empty entries
            if (length(filtering_progress@total_peptides) > 0) {
                filtering_progress@total_peptides <- c(filtering_progress@total_peptides, 
                                                     filtering_progress@total_peptides[length(filtering_progress@total_peptides)])
                filtering_progress@peptides_per_protein <- c(filtering_progress@peptides_per_protein, 
                                                           filtering_progress@peptides_per_protein[length(filtering_progress@peptides_per_protein)])
                filtering_progress@peptides_per_run <- c(filtering_progress@peptides_per_run, 
                                                       filtering_progress@peptides_per_run[length(filtering_progress@peptides_per_run)])
            } else {
                filtering_progress@total_peptides <- c(filtering_progress@total_peptides, NA_integer_)
                filtering_progress@peptides_per_protein <- c(filtering_progress@peptides_per_protein, 
                                                           list(data.frame(Protein.Ids = character(), 
                                                                         n_peptides = integer())))
                filtering_progress@peptides_per_run <- c(filtering_progress@peptides_per_run, 
                                                       list(data.frame(Run = character(), 
                                                                     n_peptides = integer())))
            }
        }
    }
    
    # Update the global filtering_progress object
    assign("filtering_progress", filtering_progress, envir = .GlobalEnv)
    
    # Create base protein count plot (always shown)
    p1 <- ggplot(data.frame(
        step = factor(filtering_progress@steps, levels = filtering_progress@steps),
        proteins = filtering_progress@proteins
    ), aes(x = step, y = proteins)) +
        geom_bar(stat = "identity", fill = "steelblue", width = 0.7) +
        geom_text(aes(label = proteins), 
                  vjust = -0.5, 
                  size = 4) +
        labs(
            title = "Number of Proteins",
            x = "Filtering Step",
            y = "Unique Proteins"
        ) +
        theme_minimal() +
        theme(
            axis.text.x = element_text(angle = 45, hjust = 1),
            panel.grid.major.x = element_blank()
        )
    
    # Create proteins per run plot (always shown)
    # First ensure all data frames in the list have consistent column types
    proteins_per_run_list <- lapply(filtering_progress@proteins_per_run, function(df) {
        df$Run <- as.character(df$Run)
        df$n_proteins <- as.numeric(df$n_proteins)
        return(df)
    })
    
    p4 <- bind_rows(proteins_per_run_list, .id = "step") |>
        mutate(step = filtering_progress@steps[as.numeric(step)]) |>
        group_by(Run) |>
        mutate(avg_proteins = mean(n_proteins)) |>
        ungroup() |>
        # Run is already character from our preprocessing
        mutate(Run = fct_reorder(Run, avg_proteins)) |>
        ggplot(aes(x = Run, y = n_proteins, 
                  group = step, 
                  color = factor(step, levels = filtering_progress@steps))) +
        geom_line() +
        geom_point() +
        labs(
            title = "Proteins per Run",
            x = "Run ID (ordered by average protein count)",
            y = "Number of Proteins",
            color = "Step"
        ) +
        theme_minimal() +
        theme(
            axis.text.x = element_text(angle = 45, hjust = 1),
            panel.grid.major.x = element_blank()
        ) +
        scale_color_manual(values = get_color_palette(length(filtering_progress@steps), "steelblue"))
    
    # Initialize peptide plots
    if (is_protein_quant) {
        # For protein data, create empty placeholder plots if no peptide data exists
        if (all(is.na(filtering_progress@total_peptides))) {
            p2 <- p3 <- p5 <- ggplot() + 
                annotate("text", x = 0.5, y = 0.5, 
                        label = "No peptide data available for protein quantification data") +
                theme_void()
        } else {
            # If peptide data exists from previous steps, create plots with existing data
            p2 <- ggplot(data.frame(
                step = factor(filtering_progress@steps, levels = filtering_progress@steps),
                total_peptides = filtering_progress@total_peptides
            ), aes(x = step, y = total_peptides)) +
                geom_bar(stat = "identity", fill = "forestgreen", width = 0.7) +
                geom_text(aes(label = total_peptides), 
                          vjust = -0.5, 
                          size = 4) +
                labs(
                    title = "Total Unique Peptides (from last peptide data)",
                    x = "Filtering Step",
                    y = "Unique Peptides"
                ) +
                theme_minimal() +
                theme(
                    axis.text.x = element_text(angle = 45, hjust = 1),
                    panel.grid.major.x = element_blank()
                )
            
            # Ensure consistent data types in peptides_per_protein list
            peptides_per_protein_list <- lapply(filtering_progress@peptides_per_protein, function(df) {
                if (nrow(df) > 0) {
                    df$Protein.Ids <- as.character(df$Protein.Ids)
                    df$n_peptides <- as.numeric(df$n_peptides)
                }
                return(df)
            })
            
            p3 <- ggplot() +
                geom_boxplot(data = bind_rows(peptides_per_protein_list, .id = "step") |>
                             mutate(step = filtering_progress@steps[as.numeric(step)]),
                           aes(x = factor(step, levels = filtering_progress@steps), 
                               y = n_peptides),
                           fill = "darkred",
                           alpha = 0.5,
                           outlier.shape = NA) +
                labs(
                    title = "Peptides per Protein Distribution (from last peptide data)",
                    x = "Filtering Step",
                    y = "Number of Peptides"
                ) +
                theme_minimal() +
                theme(
                    axis.text.x = element_text(angle = 45, hjust = 1),
                    panel.grid.major.x = element_blank()
                ) +
                coord_cartesian(
                    ylim = c(0, 
                             quantile(bind_rows(peptides_per_protein_list)$n_peptides, 0.95))
                )
            
            # Ensure consistent data types in peptides_per_run list
            peptides_per_run_list <- lapply(filtering_progress@peptides_per_run, function(df) {
                if (nrow(df) > 0) {
                    df$Run <- as.character(df$Run)
                    df$n_peptides <- as.numeric(df$n_peptides)
                }
                return(df)
            })
            
            p5 <- bind_rows(peptides_per_run_list, .id = "step") |>
                mutate(step = filtering_progress@steps[as.numeric(step)]) |>
                group_by(Run) |>
                mutate(avg_peptides = mean(n_peptides)) |>
                ungroup() |>
                # Run is already character from our preprocessing
                mutate(Run = fct_reorder(Run, avg_peptides)) |>
                ggplot(aes(x = Run, y = n_peptides, 
                          group = step, 
                          color = factor(step, levels = filtering_progress@steps))) +
                geom_line() +
                geom_point() +
                labs(
                    title = "Peptides per Run (from last peptide data)",
                    x = "Run ID (ordered by average peptide count)",
                    y = "Number of Peptides",
                    color = "Step"
                ) +
                theme_minimal() +
                theme(
                    axis.text.x = element_text(angle = 45, hjust = 1),
                    panel.grid.major.x = element_blank()
                ) +
                scale_color_manual(values = get_color_palette(length(filtering_progress@steps), "forestgreen"))
        }
    } else {
        # For peptide data, create normal plots
        p2 <- ggplot(data.frame(
            step = factor(filtering_progress@steps, levels = filtering_progress@steps),
            total_peptides = filtering_progress@total_peptides
        ), aes(x = step, y = total_peptides)) +
            geom_bar(stat = "identity", fill = "forestgreen", width = 0.7) +
            geom_text(aes(label = total_peptides), 
                      vjust = -0.5, 
                      size = 4) +
            labs(
                title = "Total Unique Peptides",
                x = "Filtering Step",
                y = "Unique Peptides"
            ) +
            theme_minimal() +
            theme(
                axis.text.x = element_text(angle = 45, hjust = 1),
                panel.grid.major.x = element_blank()
            )
        
        # Ensure consistent data types in peptides_per_protein list
        peptides_per_protein_list <- lapply(filtering_progress@peptides_per_protein, function(df) {
            if (nrow(df) > 0) {
                df$Protein.Ids <- as.character(df$Protein.Ids)
                df$n_peptides <- as.numeric(df$n_peptides)
            }
            return(df)
        })
        
        p3 <- ggplot() +
            geom_boxplot(data = bind_rows(peptides_per_protein_list, .id = "step") |>
                         mutate(step = filtering_progress@steps[as.numeric(step)]),
                       aes(x = factor(step, levels = filtering_progress@steps), 
                           y = n_peptides),
                       fill = "darkred",
                       alpha = 0.5,
                       outlier.shape = NA) +
            labs(
                title = "Peptides per Protein Distribution",
                x = "Filtering Step",
                y = "Number of Peptides"
            ) +
            theme_minimal() +
            theme(
                axis.text.x = element_text(angle = 45, hjust = 1),
                panel.grid.major.x = element_blank()
            ) +
            coord_cartesian(
                ylim = c(0, 
                         quantile(bind_rows(peptides_per_protein_list)$n_peptides, 0.95))
            )
        
        # Ensure consistent data types in peptides_per_run list
        peptides_per_run_list <- lapply(filtering_progress@peptides_per_run, function(df) {
            if (nrow(df) > 0) {
                df$Run <- as.character(df$Run)
                df$n_peptides <- as.numeric(df$n_peptides)
            }
            return(df)
        })
        
        p5 <- bind_rows(peptides_per_run_list, .id = "step") |>
            mutate(step = filtering_progress@steps[as.numeric(step)]) |>
            group_by(Run) |>
            mutate(avg_peptides = mean(n_peptides)) |>
            ungroup() |>
            # Run is already character from our preprocessing
            mutate(Run = fct_reorder(Run, avg_peptides)) |>
            ggplot(aes(x = Run, y = n_peptides, 
                      group = step, 
                      color = factor(step, levels = filtering_progress@steps))) +
            geom_line() +
            geom_point() +
            labs(
                title = "Peptides per Run",
                x = "Run ID (ordered by average peptide count)",
                y = "Number of Peptides",
                color = "Step"
            ) +
            theme_minimal() +
            theme(
                axis.text.x = element_text(angle = 45, hjust = 1),
                panel.grid.major.x = element_blank()
            ) +
            scale_color_manual(values = get_color_palette(length(filtering_progress@steps), "forestgreen"))
    }
    
    # Create plot list based on data type
    plot_list <- list(
        proteins_total = p1,
        proteins_per_run = p4,
        peptides_total = p2,
        peptides_per_protein = p3,
        peptides_per_run = p5
    )
    
    # Save plots if derived_time_dir is valid and save_plots is TRUE
    if (save_plots) {
        for (plot_name in names(plot_list)) {
            filename <- file.path(derived_time_dir,
                                sprintf("%s_%s.png", step_name, plot_name))
            ggsave(filename, 
                   plot = plot_list[[plot_name]], 
                   width = 10, 
                   height = 8, 
                   dpi = 300)
        }
    }
    
    # Return/display plots based on return_grid parameter
    if(return_grid) {
        if (!is_protein_quant || !all(is.na(filtering_progress@total_peptides))) {
            # Create full grid with all plots if peptide data exists
            grid1 <- gridExtra::arrangeGrob(p1, p2, p3, ncol = 3)
            grid2 <- gridExtra::arrangeGrob(p4, ncol = 1)
            grid3 <- gridExtra::arrangeGrob(p5, ncol = 1)
            
            grid_plot <- gridExtra::grid.arrange(
                grid1,
                grid2,
                grid3,
                heights = c(1, 1, 1)
            )
        } else {
            # For protein_quant_table without peptide data, only show protein plots
            grid_plot <- gridExtra::grid.arrange(
                p1,
                p4,
                ncol = 1,
                heights = c(1, 1)
            )
        }
        
        # Save the grid if derived_time_dir is valid and save_plots is TRUE
        if (save_plots) {
            filename <- file.path(derived_time_dir,
                                sprintf("%s_combined_plots.png", step_name))
            ggsave(filename, 
                   plot = grid_plot, 
                   width = 15, 
                   height = if (!is_protein_quant || !all(is.na(filtering_progress@total_peptides))) 18 else 12,
                   dpi = 300)
        }
        
        return(grid_plot)
    } else {
        # Print each plot individually
        for(plot_obj in plot_list) { # Changed loop variable to avoid conflict with base::plot
            print(plot_obj)
        }
        # Return the list invisibly
        invisible(plot_list)
    }
}

