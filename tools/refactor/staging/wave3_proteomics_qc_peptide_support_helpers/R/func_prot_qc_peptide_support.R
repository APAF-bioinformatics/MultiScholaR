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
    cat(sprintf("Dataset dimensions: %d peptides x %d samples\n", 
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

