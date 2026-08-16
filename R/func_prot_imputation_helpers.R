# ----------------------------------------------------------------------------
# proteinMissingValueImputation
# ----------------------------------------------------------------------------
#' proteinMissingValueImputation
#' @description Perform protein level missing value imputation
#'@param input_table A data frame with the following columns: 1. Sample file name or Run name, 2. Protein IDs, 3. Normalised protein abundances
#'@param metadata_table A data table with the following columns: 1. the sample file name or run name (as per parameter sample_id_tbl_sample_id_column), 2. The replicate group ID (as per parameter replicate_group_column)
#'@param input_table_sample_id_column The name of the column in the input_table that contained the run information or sample file name as per the input_table parameter (default: Run)
#'@param sample_id_tbl_sample_id_column The name of the column in the input_table that contained the run information or sample file name as per the metadata_table parameter (default: ms_filename)
#'@param replicate_group_column (default: general_sample_info)
#'@param protein_id_column Protein accession column, tidyverse format (default = Protein.Ids).
#'@param quantity_to_impute_column Name of column containing the peptide abundance that needs to be normalised in tidyverse format (default: Peptide.RawQuantity)
#'@param hek_string The string denoting samples that are controls using HEK cells (default: "HEK")
#'@export
#' @param imputed_value_column,core_utilisation Runtime inputs used by this function; see the usage section for accepted values.
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

  if( length(which(is.na(core_utilisation))) > 0 ) {
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

