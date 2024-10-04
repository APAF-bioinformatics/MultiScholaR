#---------------------------------------------------------------------------------------
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


# Count the number of samples in the input table
#' @export
count_num_samples <- function( input_table
                               , sample_id_column = Run) {
  num_samples <- input_table |>
    distinct( {{sample_id_column}}) |>
    count()

  num_samples[[1,1]]
}



# Count the number of peptides in the input table
#' @export
count_num_proteins <- function( input_table
                                , protein_id_column = Protein.Ids) {
  num_proteins <- input_table |>
    distinct( {{protein_id_column}}) |>
    count()

  num_proteins[[1,1]]
}


# Count the number of peptides in the input table
#' @export
count_num_peptides <- function( input_table
                                , protein_id_column = Protein.Ids
                                , peptide_sequence_column = Stripped.Sequence ) {
  num_peptides <- input_table |>
    distinct( {{protein_id_column}}, {{peptide_sequence_column}}) |>
    count()

  num_peptides[[1,1]]
}


#' @description
#' Plot the number of replicate samples analyzed by each mass spec. machine
#'
#' @export
plotNumSamplesPerMachine <- function( metadata_table
                                      , machine_id_column = Machine
                                      , machine_colour_values = getCmriMachineColour()) {

  counts_tbl <- metadata_table |>
    group_by( {{machine_id_column}} ) |>
    summarise(counts  = n() ) |>
    ungroup()

  results_plot <-  counts_tbl |>
    ggplot( aes( {{machine_id_column}}, counts, fill={{machine_id_column}})) +
    geom_bar(stat = "identity") +
    apafTheme() +
    ylab( "Number of sample runs") +
    ggtitle("Sample Distribution") +
    scale_y_continuous( expand = expansion(  mult=c(0, 0.1)))

  #  xlim( 0, x_max_value)

  if(!is.null(machine_colour_values)) {
    results_plot <- results_plot +
      scale_fill_manual( values= machine_colour_values)
  }

  results_plot

}


#' plotPeptidesProteinsCountsPerSampleHelper
#' @description Plot the number of proteins and peptides identified per sample
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



  search_srl_quant_cln <- input_table |>
    dplyr::filter( {{q_value_column}} < qvalue_threshold &
                     {{global_q_value_column}} < global_qvalue_threshold &
                     {{proteotypic_peptide_sequence_column}} == choose_only_proteotypic_peptide ) |>
    dplyr::select(all_of( input_matrix_column_ids))

  search_srl_quant_cln

}

#' @description  Peptides of with charges and modifications are rolled up (summed) together
#' @export
rollUpPrecursorToPeptideHelper <- function( input_table
                                      , sample_id_column = Run
                                      , protein_id_column = Protein.Ids
                                      , peptide_sequence_column = Stripped.Sequence
                                      , modified_peptide_sequence_column = Modified.Sequence
                                      , precursor_quantity_column = Precursor.Quantity
                                      , precursor_normalized_column = Precursor.Normalised
                                      , core_utilisation) {

  peptide_normalized_tbl <- NA
  if( length(which(is.na(core_utilisation))) == 0 ) {

    peptide_normalized_tbl <- input_table  |>
      group_by( {{sample_id_column}}, {{protein_id_column}}, {{peptide_sequence_column}}, {{modified_peptide_sequence_column}} ) |>
      summarise( Peptide.RawQuantity = sum( {{precursor_quantity_column}} )
                 ,  Peptide.Normalized = sum( {{precursor_normalized_column}} ) ) |>
      ungroup() |>
      group_by( {{sample_id_column}}, {{protein_id_column}}, {{peptide_sequence_column}} ) |>
      summarise( Peptide.RawQuantity = sum( Peptide.RawQuantity )
                 ,  Peptide.Normalized = sum( Peptide.Normalized )
                 ,  peptidoform_count = n()) |>
      ungroup()

  } else {
    peptide_normalized_tbl <- input_table  |>

      group_by( {{sample_id_column}}, {{protein_id_column}}, {{peptide_sequence_column}}, {{modified_peptide_sequence_column}} ) |>
      partition(core_utilisation) |>
      summarise( Peptide.RawQuantity = sum( {{precursor_quantity_column}} )
                 ,  Peptide.Normalized = sum( {{precursor_normalized_column}} ) ) |>
      collect() |>
      ungroup() |>

      group_by( {{sample_id_column}}, {{protein_id_column}}, {{peptide_sequence_column}} ) |>
      partition(core_utilisation) |>
      summarise( Peptide.RawQuantity = sum( Peptide.RawQuantity )
                 ,  Peptide.Normalized = sum( Peptide.Normalized )
                 , peptidoform_count = n() ) |>
      collect() |>
      ungroup()

  }

  peptide_normalized_tbl
}

#' @export
#' @description Remove peptide based on the intensity threshold and the proportion of samples below the threshold
peptideIntensityFilteringHelper <- function(input_table
                                      , min_peptide_intensity_threshold = 15
                                      , peptides_proportion_of_samples_below_cutoff = 1
                                      , protein_id_column = Protein.Ids
                                      , peptide_sequence_column = Stripped.Sequence
                                      , peptide_quantity_column = Peptide.Normalized
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

  peptide_normalized_pif_cln <- input_table |>
    inner_join ( num_values_per_peptide |>
                   dplyr::select( -num_below_intesnity_treshold, -samples_counts)
                 , by = join_by( {{protein_id_column}}, {{peptide_sequence_column}} ) )


  peptide_normalized_pif_cln


}



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
                                               , cols
                                               , design_matrix
                                               , sample_id
                                               , protein_id_column
                                               , peptide_sequence_column
                                               , group_column
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
    distinct( {{sample_id}} , {{ group_column }} ) |>
    group_by( {{ group_column }} ) |>
    summarise(  num_per_group = n()) |>
    ungroup()

  count_values_missing_per_group <- abundance_long |>
    mutate(is_missing = ifelse( !is.na( !!sym( abundance_column ))
                                & !!sym( abundance_column ) > abundance_threshold
                                , 0, 1)) |>
    group_by( row_id, {{ group_column }} ) |>
    summarise( num_missing_per_group = sum(is_missing)) |>
    ungroup()

  count_percent_missing_per_group <- count_values_missing_per_group |>
    full_join( count_values_per_group,
               by = join_by( {{ group_column }} )) |>
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

#' @export
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


#' @export
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

#--------------------------------------------------------------------------------------------------------------------------

#' @export
#' @description Log 2 transformation with pseudo count
log2Transformation <- function(input_matrix) {

  pseudo_count <- min( input_matrix[input_matrix> 0] , na.rm=TRUE)/100
  input_matrix[input_matrix> 0 & !is.na(input_matrix)] <- input_matrix[input_matrix> 0 & !is.na(input_matrix)] + pseudo_count
  input_matrix <- log2(input_matrix)

  return(input_matrix )
}


#--------------------------------------------------------------------------------------------------------------------------

#'@param input_table A table with two columns, the Run ID column and the technical replicate group column
#'@param run_id_column A string representing the name of the column with the Run ID (or sample ID)
#'@param replicate_group_column A string representing the name of the column with the technical replicate group
#'@return A table with three columns, the technical replicate group column, run ID X column, and run ID Y column
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



#'@description Calculate the Pearson correlation of the abundances of peptides between two samples X and Y.
#'@param ms_filename_x A string representing the sample file name X (for a pair of sample in the same technical replicate group) for correlation score calculation.
#'@param ms_filename_y A string representing the sample file name Y (for a pair of sample in the same technical replicate group) for correlation score calculation.
#'@param input_table A data frame with the following columns: 1. Sample file name or Run name, 2. Protein IDs, 3. Stripped peptide sequence, 4. Normalized peptide abundances
#'@param sample_id_column Sample ID column, tidyverse format (default = Run).
#'@param protein_id_column Protein accession column, tidyverse format (default = Protein.Ids).
#'@param peptide_sequence_column Peptide sequence column, tidyverse fromat (default =  Stripped.Sequence).
#'@param peptide_normalized_column Normalized peptide abundance column name as string (default = "Peptide.Normalized").
#'@return The pearson correlation value of the abundances of peptides between two samples X and Y.
#' @export
calulatePearsonCorrelation <- function( ms_filename_x, ms_filename_y, input_table
                                        , sample_id_column = Run
                                        , protein_id_column = Protein.Ids
                                        , peptide_sequence_column = Stripped.Sequence
                                        , peptide_normalized_column = "Peptide.Normalized")  {

  tab_x <- input_table |>
    dplyr::filter( {{sample_id_column}} == ms_filename_x )

  tab_y <- input_table |>
    dplyr::filter( {{sample_id_column}} == ms_filename_y )

  merged_tbl <- tab_x |>
    inner_join( tab_y, by=join_by( {{protein_id_column}}, {{peptide_sequence_column}}) )

  # merged_tbl |>
  #   dplyr::filter(!is.na( !!sym(paste0(peptide_normalized_column, ".x")) ) & !is.na( !!sym(paste0(peptide_normalized_column, ".x")))) |>
  #   head() |> print()

  # print( paste(ms_filename_x, ms_filename_y))
  input_x <-  merged_tbl[[ paste0(peptide_normalized_column, ".x") ]]
  input_y <- merged_tbl[[paste0(peptide_normalized_column, ".y")]]
  if( length(input_x) > 0 & length(input_y) >0  ) {
    cor_result <- cor( input_x
                       , input_y
                       , use="pairwise.complete.obs")

    cor_result
  } else {
    return( NA )
  }

}


#' @export
calulatePearsonCorrelationForSamplePairsHelper <- function( samples_id_tbl
                                                      , run_id_column = "ms_filename"
                                                      , replicate_group_column = "general_sample_info"
                                                      , input_table
                                                      , num_of_cores = 1
                                                      , sample_id_column = Run
                                                      , protein_id_column = Protein.Ids
                                                      , peptide_sequence_column = Stripped.Sequence
                                                      , peptide_normalized_column = "Peptide.Normalized") {


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
                                                                                        , peptide_normalized_column = {{peptide_normalized_column}}) }))

  pearson_correlation_per_pair

}


#'@description Remove samples which is correlated with any technical replicate samples
#'@param pearson_correlation_per_pair A data frame with the following columns: 1. ID of technical replicate group, 2. sample file name X, 3. sample file name Y, 4. Pearson correlation of the abundances of peptides between sample X and Y.
#'@param peptide_keep_samples_with_min_num_peptides A data frame with the following columns: 1. Sample file name or Run name, 2. Protein IDs, 3. Stripped peptide sequence, 4. Normalized peptide abundances
#'@param min_pearson_correlation_threshold Minimum pearson correlation for a pair of files to be considered to be consistent and kept for further analysis
#'@param filename_column_x Name of column containing the sample file name X (for a pair of sample in the same technical replicate group). Tidyverse column header format, not a string.
#'@param filename_column_y Name of column containing the sample file name Y (for a pair of sample in the same technical replicate group). Tidyverse column header format, not a string.
#'@param correlation_column Name of column containing the Pearson's correlation score between Sample X and Y. Tidyverse column header format, not a string.
#'@param filename_id_column A string indicating the name of column that contains the sample ID or Run ID in the data frame `peptide_keep_samples_with_min_num_peptides`.
#'@return A table without samples that are poorly correlated with the rest of the samples in the technical replicate group. Contains the following columns: 1. Sample file name or Run name, 2. Protein IDs, 3. Stripped peptide sequence, 4. Normalized peptide abundances
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


#'@description Remove samples which is correlated with any technical replicate samples
#'@param protein_intensity_table A data frame with the following columns: 1. ID of technical replicate group, 2. sample file name X, 3. sample file name Y, 4. Pearson correlation of the abundances of peptides between sample X and Y.
#'@param peptide_keep_samples_with_min_num_peptides A data frame with the proteins as rows and samples ID as columns.
#'@param min_pearson_correlation_threshold Minimum pearson correlation for a pair of files to be considered to be consistent and kept for further analysis
#'@param filename_column_x Name of column containing the sample file name X (for a pair of sample in the same technical replicate group). Tidyverse column header format, not a string.
#'@param filename_column_y Name of column containing the sample file name Y (for a pair of sample in the same technical replicate group). Tidyverse column header format, not a string.
#'@param protein_id_column Name of column containing the protein ID. Tidyverse column header format, not a string.
#'@param correlation_column Name of column containing the Pearson's correlation score between Sample X and Y. Tidyverse column header format, not a string.
#'@param filename_id_column A string indicating the name of column that contains the sample ID or Run ID in the data frame `peptide_keep_samples_with_min_num_peptides`.
#'@return A table without samples that are poorly correlated with the rest of the samples in the technical replicate group. Contains the following columns: 1. Sample file name or Run name, 2. Protein IDs, 3. Stripped peptide sequence, 4. Normalized peptide abundances
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
  samples_to_keep_anyway <-setdiff( setdiff(colnames(protein_intensity_table), (all_samples |> pull( temp_column )))
                                    ,  as_string({{protein_id_column}})  )

  print( samples_to_keep_anyway)

  # Samples in the table to keep
  samples_to_keep_subset <- colnames(protein_intensity_table)[colnames(protein_intensity_table) %in% (samples_to_keep |> pull( temp_column ))]

  samples_above_correlation_theshold <- protein_intensity_table |>
    dplyr::select( {{protein_id_column}}, all_of( c(samples_to_keep_anyway, samples_to_keep_subset)))

  samples_above_correlation_theshold

}

#---------------------------------------------------------------------------------------
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

#---------------------------------------------------------------------------------------
#' @description Remove proteins that only have data for one technical replicate for all sample.
#' This can be repurposed for removing proteins that only have one biological replicates for all experimental groups.
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

#' peptideMissingValueImputationHelper
#' @description Perform peptide level missing value imputation
#'@param input_table A data frame with the following columns: 1. Sample file name or Run name, 2. Protein IDs, 3. Stripped peptide sequence, 4. Normalized peptide abundances
#'@param metadata_table A data table with the following columns: 1. the sample file name or run name (as per parameter sample_id_tbl_sample_id_column), 2. The replicate group ID (as per parameter replicate_group_column)
#'@param input_table_sample_id_column The name of the column in the input_table that contained the run information or sample file name as per the input_table parameter (default: Run)
#'@param sample_id_tbl_sample_id_column The name of the column in the input_table that contained the run information or sample file name as per the metadata_table parameter (default: ms_filename)
#'@param replicate_group_column (default: general_sample_info)
#'@param protein_id_column Protein accession column, tidyverse format (default = Protein.Ids).
#'@param peptide_sequence_column Peptide sequence column, tidyverse fromat (default =  Stripped.Sequence).
#'@param quantity_to_impute_column Name of column containing the peptide abundance that needs to be normalized in tidyverse format (default: Peptide.RawQuantity)
#'@param hek_string The string denoting samples that are controls using HEK cells (default: "HEK")
#'@param proportion_missing_values The proportion of sample replicates in a group that is missing below which the peptide intensity will be imputed (default: 0.50)
#'@export
peptideMissingValueImputationHelper <- function( input_table
                                           , metadata_table
                                           , input_table_sample_id_column = Run
                                           , sample_id_tbl_sample_id_column  =  ms_filename
                                           , replicate_group_column = general_sample_info
                                           , protein_id_column = Protein.Ids
                                           , peptide_sequence_column = Stripped.Sequence
                                           , quantity_to_impute_column = Peptide.Normalized
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

#' calculatePercentMissingPeptidePerReplicate
#' @description Calculate percentage of peptides from each sample that is missing and merge with metadata
#' @export
calculatePercentMissingPeptidePerReplicate <- function( input_table
                                                        , metadata_table
                                                        , protein_id_column = Protein.Ids
                                                        , intensity_column = Peptide.Normalized
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



#' calculatePercentMissingProteinPerReplicate
#' @description Calculate percentage of proteins from each sample that is missing and merge with metadata
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


#' @description
#' Get a standard colour for each mass spect machine at CMRI
#' @export
getCmriMachineColour <- function() {

  colour_definition <- RColorBrewer::brewer.pal(8, "Set2")
  names(colour_definition) <- paste0("M0", 1:8 )
  colour_definition
}

#' @description
#' A density plot of percent missing per individual per machine. Machine with standardized colour for each machine.
#'  Depends on the output of the function `calculatePercentMissingPerReplicate`
#' @export
plotDensityOfPercentMissingPerIndvidual <- function( percent_missing_table
                                                     , percent_missing_column = percent_missing
                                                     , machine_id_column = Machine
                                                     , colour_definition = getCmriMachineColour() ) {
  percent_missing_table |>
    ggplot( aes( {{percent_missing_column}}
                 , group = {{machine_id_column}}
                 , fill={{machine_id_column}}
                 , alpha=0.5
    )) +
    geom_density() +
    scale_fill_manual( values = colour_definition) +
    scale_alpha(guide = 'none') +
    apafTheme()  +
    scale_y_continuous( expand = expansion(  mult=c(0, 0.1)))

}

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

#' plotBarplotMissingnessPerReplicate
#' @description Plot missing rate for each replicate sample in a bar plot. Colour the bar according to the 'fill_column'.
#'  Depends on the output of the function `calculatePercentMissingPerReplicate`
#' @export
plotBarplotMissingnessPerReplicate <- function( input_table
                                                , percent_missing_column = percent_missing
                                                , fill_column = Machine
                                                , fill_colour_values = getCmriMachineColour()) {

  summary_data <- summary( input_table |> pull( {{percent_missing_column}} ) )

  #print(summary_data)

  result_plot <- input_table |>
    ggplot( aes( reorder(Run, {{percent_missing_column}})
                 , {{percent_missing_column}}
                 , fill={{fill_column}})) +
    geom_bar(stat = "identity") +
    scale_y_continuous(breaks = seq(0, summary_data[6], by=10), limits=c(0,summary_data[6])) +
    xlab("Samples") +
    ylab("Missing Rate (%)") +
    apafTheme() +
    theme(    panel.grid.major.y = element_line(color = "gray", linetype = "solid")
              , panel.grid.major.x = element_blank()
              , axis.text.x=element_blank(),
              axis.ticks.x=element_blank()) +
    geom_hline(yintercept = c(summary_data[c(2,3,5)]), linetype = "dashed", color = "grey") +
    scale_y_continuous( expand = expansion(  mult=c(0, 0.1)))

  if ( !is.null( fill_colour_values )) {
    result_plot <- result_plot +
      scale_fill_manual( values =  fill_colour_values )
  }

  result_plot

}

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


#------------------------------------------------------------------------------------------------

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
#'@param input_matrix Samples are columns, rows are proteins or peptides
#'@export
scaleCenterAndFillMissing <- function( input_matrix) {

  input_matrix_scaled <- scale(input_matrix, center = TRUE, scale = TRUE)


  min_data_point <- min(input_matrix_scaled, na.rm=TRUE)


  input_matrix_scaled_fill_missing <-  input_matrix_scaled
  input_matrix_scaled_fill_missing[is.na(input_matrix_scaled_fill_missing)] <- min_data_point*2

  input_matrix_scaled_fill_missing
}

#------------------------------------------------------------------------------------------------
#' @export
compareUmapComponentsPairs <- function(input_table, columns = c("V1", "V2","V3","V4"), covariate) {

  pm <- umap_data |>
    ggpairs( columns = columns, ggplot2::aes(colour = {{covariate}}), legend = 1)  +
    apafTheme()

  pm

}


#------------------------------------------------------------------------------------------------
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


#' getSamplesCorrelationMatrix
#' @description Calculate the Pearson's correlation score between sample
#' @param input_table Table with samples as columns and peptides as rows. Contains the log peptide intensity values.
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



#' getCategoricalColourRules
#' @export
getCategoricalColourRules <- function( metadata_tbl
                                       , metadata_column_labels
                                       , metadata_column_selected
                                       , categorical_columns
                                       , ms_machine_column
                                       , columns_to_exclude
                                       , colour_palette = getCategoricalColourPalette()
                                       , na_colour = "white" ) {

  names(metadata_column_labels) <- metadata_column_selected

  if(!is.na(ms_machine_column) ) {
    colour_palette <- setdiff(colour_palette, getCmriMachineColour())
  }

  categorical_columns_tbl <- purrr::map( setdiff( categorical_columns, ms_machine_column)
                                         , \(x)  { metadata_tbl |>
                                             dplyr::distinct(!!rlang::sym(x) ) |>
                                             dplyr::mutate( values =  purrr::map_chr( !!rlang::sym(x) , as.character ) ) |>
                                             mutate( column_name = x)  |>
                                             dplyr::select(-!!rlang::sym(x )) }  ) |>
    bind_rows()

  categorical_colour_tbl_part_1 <- categorical_columns_tbl |>
    dplyr::filter( !is.na(values) ) |>
    mutate( colours =  colour_palette[seq_len( nrow( categorical_columns_tbl|>
                                                       dplyr::filter( !is.na(values) )) ) ])

  print(categorical_colour_tbl_part_1)

  categorical_colour_tbl_part_2 <- categorical_columns_tbl |>
    dplyr::filter( is.na(values) )    |>
    mutate( values = "NA") |>
    mutate( colours = na_colour )

  categorical_colour_tbl <- categorical_colour_tbl_part_1 |>
    bind_rows( categorical_colour_tbl_part_2 ) |>
    arrange( column_name) |>
    group_by(column_name) |>
    summarise( values_list = list(values)
               , colours_list = list(colours)) |>
    ungroup() |>
    mutate( colours_values_list = purrr::map2(values_list, colours_list, \(x,y){  names(y) <- x ; return(y)  }  ) )

  colour_rules <- categorical_colour_tbl |> pull( colours_values_list )
  names( colour_rules ) <-  metadata_column_labels[categorical_colour_tbl |> pull( column_name )]

  if( ! ms_machine_column %in% columns_to_exclude) {
    # Deal with the machine colour
    list_of_machines_used <- metadata_tbl |>
      distinct( !!(sym(ms_machine_column)) ) |>
      arrange( !!(sym(ms_machine_column)) ) |>
      pull(!!(sym(ms_machine_column)))

    # get the machines that we need and format as list with the ms_machine_column as the top-level name
    # e.g. list( Machine = list( M01 = ""#8DA0CB", M02 = "#FFD92F", M04 = "#E78AC3", M05= "#66C2A5" ))
    machine_colour_rule <- list( getCmriMachineColour()[list_of_machines_used] )
    names(machine_colour_rule) <- metadata_column_labels[ms_machine_column]

    colour_rules <- c( colour_rules, machine_colour_rule )
  }

  return( colour_rules)
}

#' getOneContinousPalette
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

#' getContinousColourRules
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

#' getCategoricalAndContinuousColourRules
#' @param metadata_tbl This is the table containing sample ID and other columns containing clinical variables / metadata
#' @param metadata_column_labels This is the nice
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



#' getOneContinousPalette
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

#' getSamplesCorrelationHeatMap
#' @description get the
#' @param correlation_matrix Output from the `getSamplesCorrelationMatrix` function
#' @param metadata_tbl This is the table containing sample ID and other columns containing clinical variables / metadata
#' @param is_HEK_column A logical column in the metadata table that indicates if the sample is a HEK sample
#' @param metadata_column_selected A list of column names in string selected from the metadata tbl
#' @param metadata_column_labels A list of column names in string to rename each of the columns selected in the param `metadata_column_selected`
#' @param categorical_columns A vector of string with all the names of the categorical data column  present in the `metadata_tbl` table
#' @param continous_scale_columns  A vector of string with all the names of the continuous data column  present in the `metadata_tbl` table
#' @param ms_machine_column A string of the column name describing the mass spectrometer machine used to analyze each sample
#' @param sample_id_column A string describing the column name of the sample ID column
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

#' @description
#'  Pivot peptide intensity matrix into long format table.
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

#' @description
#' Pivot protein intensity matrix into long format
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

#' @description
#' Remove proteins that have been identified in only one replicate across all patients
#' (e.g. identified no more than one relpicate in any patient)
#' @export
removeProteinsWithOnlyOneReplicateHelper <- function(input_table
                                               , samples_id_tbl
                                               , input_table_sample_id_column = Run
                                               , sample_id_tbl_sample_id_column  =  ms_filename
                                               , replicate_group_column = general_sample_info
                                               , protein_id_column = Protein.Ids
                                               , quantity_column = Protein.Normalized
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

#' proteinMissingValueImputation
#' @description Perform protein level missing value imputation
#'@param input_table A data frame with the following columns: 1. Sample file name or Run name, 2. Protein IDs, 3. Normalized protein abundances
#'@param metadata_table A data table with the following columns: 1. the sample file name or run name (as per parameter sample_id_tbl_sample_id_column), 2. The replicate group ID (as per parameter replicate_group_column)
#'@param input_table_sample_id_column The name of the column in the input_table that contained the run information or sample file name as per the input_table parameter (default: Run)
#'@param sample_id_tbl_sample_id_column The name of the column in the input_table that contained the run information or sample file name as per the metadata_table parameter (default: ms_filename)
#'@param replicate_group_column (default: general_sample_info)
#'@param protein_id_column Protein accession column, tidyverse format (default = Protein.Ids).
#'@param quantity_to_impute_column Name of column containing the peptide abundance that needs to be normalized in tidyverse format (default: Peptide.RawQuantity)
#'@param hek_string The string denoting samples that are controls using HEK cells (default: "HEK")
#'@export
proteinMissingValueImputation <- function( input_table
                                           , metadata_table
                                           , input_table_sample_id_column = Run
                                           , sample_id_tbl_sample_id_column  =  ms_filename
                                           , replicate_group_column = general_sample_info
                                           , protein_id_column = Protein.Ids
                                           , quantity_to_impute_column = Protein.Normalized
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

#' @description
#' Protein average values from replicate samples
#' #'@export
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

#' @export
cleanDesignMatrixCleanCategoriesMap <- function( input_table, column ) {
  input_table |>
    mutate( {{column}} := purrr::map_chr( {{column}}, cleanDesignMatrixCleanCategories ) )

}

#------------------------------------------------------------------------------------------------

#' @description
#' get protein intensity heatmap
#'@export
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


#' @description
#' ProCan ggplot2 theme. Rectangle box around each plot.
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

