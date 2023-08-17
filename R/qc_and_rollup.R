#' @export
# Remove peptide based on proportion of samples below intensity threshold
peptideIntensityFiltering <- function(input_table
                                      , min_peptide_intensity_threshold = 15
                                      , proportion_samples_below_intensity_threshold = 1
                                      , protein_id_column = Protein.Ids
                                      , peptide_sequence_column = Stripped.Sequence
                                      , quantity_column = Peptide.RawQuantity
                                      , cluster) {
  num_values_per_peptide <- NA
  
  if( length(which(is.na(cluster))) == 0 ) {
    num_values_per_peptide <- input_table |>
      mutate(  below_intensity_threshold = case_when( {{quantity_column}} < min_peptide_intensity_threshold ~ 1,
                                                      TRUE ~ 0) ) |>
      group_by( {{protein_id_column}}, {{peptide_sequence_column}}) |> 
      #partition(cluster) |>
      summarise (samples_counts = n(),
                 num_below_intesnity_treshold = sum(below_intensity_threshold)) |>
      #collect() |>
      ungroup() |>
      dplyr::filter( num_below_intesnity_treshold/samples_counts < proportion_samples_below_intensity_threshold )
  } else {
    num_values_per_peptide <- input_table |>
      mutate(  below_intensity_threshold = case_when( {{quantity_column}} < min_peptide_intensity_threshold ~ 1,
                                                      TRUE ~ 0) ) |>
      group_by( {{protein_id_column}}, {{peptide_sequence_column}}) |> 
      partition(cluster) |>
      summarise (samples_counts = n(),
                 num_below_intesnity_treshold = sum(below_intensity_threshold)) |>
      collect() |>
      ungroup() |>
      dplyr::filter( num_below_intesnity_treshold/samples_counts < proportion_samples_below_intensity_threshold )    
    
  }
  
  peptide_normalized_pif_cln <- input_table |>
    inner_join ( num_values_per_peptide |>
                   dplyr::select( -num_below_intesnity_treshold, -samples_counts)
                 , by = join_by( {{protein_id_column}}, {{peptide_sequence_column}} ) ) 
  
  
  peptide_normalized_pif_cln
  
  
}

