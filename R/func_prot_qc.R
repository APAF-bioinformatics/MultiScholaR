# ============================================================================
# func_prot_qc.R
# ============================================================================
# Purpose: Protein-level quality control and filtering functions
# 
# This file contains functions for protein-level QC filtering, including
# intensity filtering, missing value filtering, replicate filtering, and
# protein cleanup operations. Functions in this file are used by 
# mod_prot_qc_protein.R and related QC modules.
#
# NOTE: This file handles protein-level QC. Peptide-level QC is in
# func_prot_qc_peptide.R. The "protein" suffix distinguishes it from
# peptide QC functions.
#
# Functions to extract here:
# - proteinIntensityFiltering(): S4 method for protein intensity filtering
# - removeProteinsWithOnlyOneReplicate(): S4 method for replicate filtering
# - removeRowsWithMissingValuesPercent(): S4 method for missing value filtering
# - removeEmptyRows(): Remove rows with all zero/NA values
# - calculatePercentMissingPerProtein(): Calculate missing value percentages
# - Additional protein QC helper functions
#
# Dependencies:
# - dplyr, tidyr
# - func_general_helpers.R (for utility functions)
# ============================================================================

# TODO: Extract the following functions from their current locations:

# Function 1: proteinIntensityFiltering()
# Current location: R/proteinVsSamplesS4Objects.R
# Type: S4 method (exportMethods)
# Description: Filters proteins based on intensity thresholds
# setMethod(f = "proteinIntensityFiltering", ...) {
#   # Extract from R/proteinVsSamplesS4Objects.R
# }

# Function 2: proteinIntensityFilteringHelper()
# Current location: R/proteinVsSamplesS4Objects.R
# Description: Helper function for protein intensity filtering
# proteinIntensityFilteringHelper <- function(...) {
#   # Extract from R/proteinVsSamplesS4Objects.R
# }

# Function 3: removeProteinsWithOnlyOneReplicate()
# Current location: R/proteinVsSamplesS4Objects.R
# Type: S4 method (exportMethods)
# Description: Removes proteins present in only one replicate
# setMethod(f = "removeProteinsWithOnlyOneReplicate", ...) {
#   # Extract from R/proteinVsSamplesS4Objects.R
# }

# Function 4: removeProteinWithOnlyOneReplicate()
# Current location: R/de_proteins_functions.R
# Description: Removes proteins with only one replicate (non-S4 version)
# removeProteinWithOnlyOneReplicate <- function(...) {
#   # Extract from R/de_proteins_functions.R
# }

# Function 5: removeProteinsWithOnlyOneReplicateHelper()
# Current location: R/proteinVsSamplesS4Objects.R
# Description: Helper for removing proteins with only one replicate
# removeProteinsWithOnlyOneReplicateHelper <- function(...) {
#   # Extract from R/proteinVsSamplesS4Objects.R
# }

# Function 6: removeRowsWithMissingValuesPercent()
# Current location: R/proteinVsSamplesS4Objects.R
# Type: S4 method (exportMethods)
# Description: Removes proteins with excessive missing values
# setMethod(f = "removeRowsWithMissingValuesPercent", ...) {
#   # Extract from R/proteinVsSamplesS4Objects.R
# }

# Function 7: removeRowsWithMissingValuesPercentHelper()
# Current location: R/proteinVsSamplesS4Objects.R
# Description: Helper for removing rows with missing values
# removeRowsWithMissingValuesPercentHelper <- function(...) {
#   # Extract from R/proteinVsSamplesS4Objects.R
# }

# Function 8: removeEmptyRows()
# Current location: R/de_proteins_functions.R
# Description: Removes rows where all sample columns are zero or NA
# removeEmptyRows <- function(input_table, col_pattern, row_id) {
#   # Extract from R/de_proteins_functions.R
# }

# Function 9: removeRowsWithMissingValues()
# Current location: R/helper_functions.R
# Description: Removes rows with missing values
# removeRowsWithMissingValues <- function(...) {
#   # Extract from R/helper_functions.R
# }

# Function 10: calculatePercentMissingPerProtein()
# Current location: R/de_proteins_functions.R
# Description: Calculates percentage of missing values per protein
# calculatePercentMissingPerProtein <- function(...) {
#   # Extract from R/de_proteins_functions.R
# }

# Function 11: calculatePercentMissingProteinPerReplicate()
# Current location: R/de_proteins_functions.R
# Description: Calculates percentage of missing values per replicate
# calculatePercentMissingProteinPerReplicate <- function(...) {
#   # Extract from R/de_proteins_functions.R
# }

# Function 12: calculatePercentMissingPeptidePerReplicate()
# Current location: R/de_proteins_functions.R
# Description: Calculates percentage of missing peptides per replicate
# calculatePercentMissingPeptidePerReplicate <- function(...) {
#   # Extract from R/de_proteins_functions.R
# }

# Function 13: calculateMissingValuesPerProteinFishersTest()
# Current location: R/de_proteins_functions.R
# Description: Performs Fisher's test on missing values per protein
# calculateMissingValuesPerProteinFishersTest <- function(...) {
#   # Extract from R/de_proteins_functions.R
# }

# Function 14: checkProteinNAPercentages()
# Current location: R/de_proteins_functions.R
# Description: Checks protein NA percentages
# checkProteinNAPercentages <- function(...) {
#   # Extract from R/de_proteins_functions.R
# }

# Function 15: getProteinNARecommendations()
# Current location: R/de_proteins_functions.R
# Description: Gets recommendations for protein NA filtering
# getProteinNARecommendations <- function(...) {
#   # Extract from R/de_proteins_functions.R
# }


# ----------------------------------------------------------------------------
# removeEmptyRows
# ----------------------------------------------------------------------------
#'Remove rows in the table where the columns specified by the column regular expression pattern are all zero or NA value.
#'@param input_table Input table with columns recording protein abundances for each sample. The name of these columns matches a regular expression pattern, defined by 'col_pattern'. Remove rows with all samples having no protein abundance.
#'@param col_pattern String representing regular expression pattern that matches the name of columns containing the protein abundance values.
#'@param row_id The column name with the row_id, tidyverse style name.
#'@return A data frame with the rows without abundance values removed.
#'@export
removeEmptyRows <- function(input_table, col_pattern, row_id) {

  temp_col_name <- paste0("temp_", as_string(as_name(enquo(row_id))))

  temp_input_table <- input_table |>
    dplyr::mutate(!!rlang::sym(temp_col_name) := row_number())

  sites_to_accept <- temp_input_table |>
    mutate(across(matches(col_pattern, perl = TRUE), \(x){ (is.na(x) | x == 0) })) |>
    dplyr::filter(!if_all(matches(col_pattern, perl = TRUE), \(x){ x == TRUE })) |>
    dplyr::select({ { temp_col_name } })

  ## Removing entries where all the "Reporter intensity corrected" rows are zero
  filtered_table <- temp_input_table |>
    inner_join(sites_to_accept, by = temp_col_name) |>
    dplyr::select(-temp_col_name)

  return(filtered_table)
}


# ----------------------------------------------------------------------------
# removeProteinsWithOnlyOneReplicateHelper
# ----------------------------------------------------------------------------
#' @title Remove Proteins with Single Replicate
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


# ----------------------------------------------------------------------------
# removeRowsWithMissingValues
# ----------------------------------------------------------------------------
## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#' For each experimental group, identify proteins that have more than accepted number of missing values per group.
#'@param input_table An input table with a column containing the row ID and the rest of the columns representing abundance values for each sample.
#'@param cols A tidyselect command to select the columns. This includes the functions starts_with(), ends_with(), contains(), matches(), and num_range()
#'@param design_matrix A data frame with a column containing the sample ID (as per the sample_id param) and the experimental group (as per the group param). Each row as the sample ID as row name in the data frame.
#'@param sample_id The name of the column in design_matrix table that has the sample ID.
#'@param row_id A unique ID for each row of the 'input_table' variable.
#'@param grouping_variable The name of the column in design_matrix table that has the experimental group.
#'@param max_num_samples_miss_per_group An integer representing the maximum number of samples with missing values per group.
#'@param abundance_threshold Abundance threshold in which the protein in the sample must be above for it to be considered for inclusion into data analysis.
#'@param temporary_abundance_column The name of a temporary column, as a string, to keep the abundance value you want to filter upon
#'@return A list, the name of each element is the sample ID and each element is a vector containing the protein accessions (e.g. row_id) with enough number of values.
#'@export
removeRowsWithMissingValues <- function(input_table, cols, design_matrix, sample_id, row_id, grouping_variable, max_num_samples_miss_per_group, abundance_threshold
                                        , temporary_abundance_column = "Abundance") {

  abundance_long <- input_table |>
    pivot_longer(cols = { { cols } },
                 names_to =  as_string(as_name(enquo(sample_id))) ,
                 values_to = temporary_abundance_column  ) |>
    mutate( {{sample_id}} := purrr::map_chr(   {{sample_id}}  , as.character)   ) |>
    left_join(design_matrix |>
                mutate(  {{sample_id}} := purrr::map_chr(    {{sample_id}} , as.character)   )
              , by = as_string(as_name(enquo(sample_id))))

  count_missing_values_per_group <- abundance_long |>
    mutate(is_missing = ifelse(!is.na( !!sym(temporary_abundance_column)) & !!sym(temporary_abundance_column) > abundance_threshold, 0, 1)) |>
    group_by( {{ row_id }}, {{ grouping_variable }} ) |>
    summarise(num_missing_values = sum(is_missing)) |>
    ungroup()

  remove_rows_temp <- count_missing_values_per_group |>
    dplyr::filter(max_num_samples_miss_per_group < num_missing_values) |>
    dplyr::select(-num_missing_values, -{ { grouping_variable } }) |>
    distinct({ { row_id } })

  filtered_tbl <- input_table |>
    dplyr::anti_join(remove_rows_temp, by = as_string(as_name(enquo(row_id))))

  return(filtered_tbl)

}


# ----------------------------------------------------------------------------
# removeRowsWithMissingValuesPercentHelper
# ----------------------------------------------------------------------------
## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#'@title Remove rows with missing values
#'@param input_table An input table with a column containing the row ID and the rest of the columns representing abundance values for each sample.
#'@param cols A tidyselect command to select the columns. This includes the functions starts_with(), ends_with(), contains(), matches(), and num_range()
#'@param design_matrix A data frame with a column containing the sample ID (as per the sample_id param) and the experimental group (as per the group param). Each row as the sample ID as row name in the data frame.
#'@param sample_id The name of the column in design_matrix table that has the sample ID.
#'@param row_id A unique ID for each row of the 'input_table' variable.
#'@param grouping_variable The name of the column in design_matrix table that has the experimental group.
#'@param groupwise_percentage_cutoff The maximum percentage of values below threshold allow in each group for a protein .
#'@param max_groups_percentage_cutoff The maximum percentage of groups allowed with too many samples with protein abundance values below threshold.
#'@param temporary_abundance_column The name of a temporary column to keep the abundance value you want to filter upon
#'@param proteins_intensity_cutoff_percentile The percentile of the protein intensity values to be used as the minimum threshold for protein intensity.
#'@return A list, the name of each element is the sample ID and each element is a vector containing the protein accessions (e.g. row_id) with enough number of values.
#'@export
removeRowsWithMissingValuesPercentHelper <- function(input_table
                                                     , cols
                                                     , design_matrix
                                                     , sample_id # symbol, e.g. Run
                                                     , row_id    # symbol
                                                     , grouping_variable # symbol
                                                     , groupwise_percentage_cutoff = 1
                                                     , max_groups_percentage_cutoff = 50
                                                     , proteins_intensity_cutoff_percentile = 1
                                                     , temporary_abundance_column = "Abundance") {

  message("╔═══════════════════════════════════════════════════════════════════════════╗")
  message("║  DEBUG66: Entering removeRowsWithMissingValuesPercentHelper (OPTIMIZED)   ║")
  message("╚═══════════════════════════════════════════════════════════════════════════╝")
  
  # 1. Setup strings and symbols
  sample_id_str <- rlang::as_string(rlang::ensym(sample_id))
  row_id_str <- rlang::as_string(rlang::ensym(row_id))
  group_var_str <- rlang::as_string(rlang::ensym(grouping_variable))
  
  message(sprintf("   DEBUG66: Pivoting %d rows x %d cols. Memory: %.1f MB", 
                  nrow(input_table), ncol(input_table), sum(gc()[,2])))
  
  # 2. Prepare minimal design matrix (Select only needed columns to save memory in join)
  design_matrix_minimal <- design_matrix |>
    dplyr::select(!!rlang::ensym(sample_id), !!rlang::ensym(grouping_variable)) |>
    dplyr::mutate(!!rlang::sym(sample_id_str) := as.character(!!rlang::sym(sample_id_str)))
  
  # 3. Pivot and Join in one pipeline to reduce intermediate copies
  # Use .data[[]] for robust column access within dplyr
  abundance_long <- input_table |>
    tidyr::pivot_longer(
      cols = !all_of(cols),
      names_to = sample_id_str,
      values_to = temporary_abundance_column
    ) |>
    dplyr::mutate(
      !!rlang::sym(sample_id_str) := as.character(!!rlang::sym(sample_id_str)),
      !!rlang::sym(temporary_abundance_column) := dplyr::if_else(
        is.nan(!!rlang::sym(temporary_abundance_column)), 
        NA_real_, 
        !!rlang::sym(temporary_abundance_column)
      )
    ) |>
    dplyr::left_join(design_matrix_minimal, by = sample_id_str)

  message(sprintf("   DEBUG66: Pivot/Join complete. Long table: %d rows. Memory: %.1f MB", 
                  nrow(abundance_long), sum(gc()[,2])))
  
  # Force garbage collection after big join
  gc()

  # 4. Calculate Threshold
  # Extract vector directly to avoid data.frame overhead for simple stats
  valid_values <- abundance_long[[temporary_abundance_column]]
  valid_values <- valid_values[!is.na(valid_values) & !is.nan(valid_values) & !is.infinite(valid_values)]
  
  min_protein_intensity_threshold <- if(length(valid_values) > 0) {
    ceiling(quantile(valid_values, probs = proteins_intensity_cutoff_percentile/100, na.rm = TRUE))[1]
  } else {
    0
  }
  message(sprintf("   DEBUG66: Threshold calculated: %g", min_protein_intensity_threshold))

  # 5. Count samples per group
  count_values_per_group <- abundance_long |>
    dplyr::distinct(!!rlang::sym(sample_id_str), !!rlang::sym(group_var_str)) |>
    dplyr::count(!!rlang::sym(group_var_str), name = "num_per_group")

  # 6. Calculate missing stats
  # Perform aggregation directly
  count_percent_missing_per_group <- abundance_long |>
    dplyr::mutate(is_missing = dplyr::if_else(
      !is.na(!!rlang::sym(temporary_abundance_column)) & 
      !!rlang::sym(temporary_abundance_column) > min_protein_intensity_threshold, 
      0, 1
    )) |>
    dplyr::group_by(!!rlang::sym(row_id_str), !!rlang::sym(group_var_str)) |>
    dplyr::summarise(num_missing_per_group = sum(is_missing), .groups = "drop") |>
    dplyr::left_join(count_values_per_group, by = group_var_str) |>
    dplyr::mutate(perc_missing_per_group = num_missing_per_group / num_per_group * 100)

  # 7. Identify proteins to remove
  total_num_of_groups <- nrow(count_values_per_group)
  
  proteins_to_remove <- count_percent_missing_per_group |>
    dplyr::filter(perc_missing_per_group > groupwise_percentage_cutoff) |>
    dplyr::group_by(!!rlang::sym(row_id_str)) |>
    dplyr::summarise(percent_groups_failed = n() / total_num_of_groups * 100, .groups = "drop") |>
    dplyr::filter(percent_groups_failed > max_groups_percentage_cutoff)
    
  num_removed <- nrow(proteins_to_remove)
  message(sprintf("   DEBUG66: Proteins to remove: %d", num_removed))

  # 8. Filter original table
  if (num_removed > 0) {
    filtered_tbl <- input_table |>
      dplyr::anti_join(proteins_to_remove, by = row_id_str)
  } else {
    filtered_tbl <- input_table
  }

  message(sprintf("   DEBUG66: Final table: %d rows. Memory: %.1f MB", 
                  nrow(filtered_tbl), sum(gc()[,2])))
  
  # Final GC
  gc()
  
  message("╔═══════════════════════════════════════════════════════════════════════════╗")
  message("║  DEBUG66: Exiting removeRowsWithMissingValuesPercentHelper                ║")
  message("╚═══════════════════════════════════════════════════════════════════════════╝")
  
  return(filtered_tbl)

}


# ----------------------------------------------------------------------------
# checkProteinNAPercentages
# ----------------------------------------------------------------------------
#' Check Missing Value Percentages in Protein Data
#' 
#' @description Calculate and report the percentage of missing values (NAs) in protein data
#' at different levels: total dataset, per sample, and per group.
#' 
#' @param protein_obj A ProteinQuantitativeData S4 object
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
checkProteinNAPercentages <- function(protein_obj, verbose = TRUE) {
  
  # Validate input
  if (!is(protein_obj, "ProteinQuantitativeData")) {
    stop("Input must be a ProteinQuantitativeData S4 object")
  }
  
  # Extract data from S4 object
  protein_quant_table <- protein_obj@protein_quant_table
  design_matrix <- protein_obj@design_matrix
  sample_id_col <- protein_obj@sample_id
  group_id_col <- protein_obj@group_id
  protein_id_col <- protein_obj@protein_id_column
  
  # Identify sample columns (exclude protein ID column)
  sample_columns <- setdiff(colnames(protein_quant_table), protein_id_col)
  
  # Validate that sample columns match design matrix
  if (length(sample_columns) != nrow(design_matrix)) {
    stop("Number of sample columns doesn't match design_matrix rows")
  }
  
  # Extract quantitative data matrix (samples only)
  protein_matrix <- as.matrix(protein_quant_table[, sample_columns])
  rownames(protein_matrix) <- protein_quant_table[[protein_id_col]]
  
  # Calculate total NA percentage
  total_values <- length(protein_matrix)
  total_nas <- sum(is.na(protein_matrix))
  total_na_percent <- (total_nas / total_values) * 100
  
  # Calculate per-sample NA percentages
  sample_na_counts <- apply(protein_matrix, 2, function(x) sum(is.na(x)))
  sample_na_percentages <- (sample_na_counts / nrow(protein_matrix)) * 100
  
  per_sample_na <- data.frame(
    sample = names(sample_na_counts),
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
    total_proteins = nrow(protein_matrix),
    total_samples = ncol(protein_matrix),
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
    cat("\n=== Protein Data Missing Value Analysis ===\n")
    cat(sprintf("Dataset dimensions: %d proteins × %d samples\n", 
                nrow(protein_matrix), ncol(protein_matrix)))
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
# getProteinNARecommendations
# ----------------------------------------------------------------------------
#' Get Recommendations for Handling Protein-Level Missing Values
#' 
#' @description Provides specific recommendations for dealing with missing values
#' in protein data based on the percentage and distribution of NAs.
#' 
#' @param protein_obj A ProteinQuantitativeData S4 object
#' @param include_code Logical, whether to include example R code (default: TRUE)
#' 
#' @return Prints recommendations and invisibly returns a list of strategies
#' 
#' @export
getProteinNARecommendations <- function(protein_obj, include_code = TRUE) {
  
  # Get NA analysis
  na_results <- checkProteinNAPercentages(protein_obj, verbose = FALSE)
  na_percent <- na_results$total_na_percent
  
  cat("\n=== PROTEIN NA HANDLING RECOMMENDATIONS ===\n")
  cat(sprintf("Your data: %.1f%% NAs across %d proteins\n\n", 
              na_percent, na_results$summary_stats$total_proteins))
  
  if (na_percent < 15) {
    cat("🎯 RECOMMENDATION: Complete Case Analysis\n")
    cat("• Your data has excellent protein coverage\n")
    cat("• Can proceed with standard analysis on proteins with complete data\n")
    if (include_code) {
      cat("\n📝 Example code:\n")
      cat("complete_proteins <- protein_obj@protein_quant_table[complete.cases(protein_obj@protein_quant_table), ]\n")
    }
    
  } else if (na_percent >= 15 && na_percent < 40) {
    cat("🎯 RECOMMENDATION: Consider Protein-Level Imputation\n")
    cat("• Moderate missing values - imputation could be beneficial\n")
    cat("• Options: KNN, minimum value, or mixed imputation strategies\n")
    cat("• Alternative: Filter to proteins detected in ≥X samples per group\n")
    if (include_code) {
      cat("\n📝 Example filtering code:\n")
      cat("# Keep proteins detected in ≥50% of samples per group\n")
      cat("filtered_proteins <- filterProteinsByGroupDetection(protein_obj, min_detection_rate = 0.5)\n")
    }
    
  } else if (na_percent >= 40 && na_percent < 60) {
    cat("🎯 RECOMMENDATION: Strict Filtering + Targeted Imputation\n")
    cat("• High missing values suggest challenging sample/detection conditions\n")
    cat("• Focus on well-detected proteins (present in majority of samples)\n")
    cat("• Consider group-wise detection requirements\n")
    if (include_code) {
      cat("\n📝 Example approach:\n")
      cat("# Keep proteins detected in ≥70% of samples in at least one group\n")
      cat("robust_proteins <- filterProteinsByGroupwise(protein_obj, min_group_detection = 0.7)\n")
    }
    
  } else {
    cat("⚠️  RECOMMENDATION: Review Data Quality\n")
    cat("• Very high missing values (>60%) suggest potential issues\n")
    cat("• Check: sample quality, peptide identification, rollup parameters\n")
    cat("• Consider more stringent protein identification criteria\n")
    cat("• May need to focus only on highly abundant/well-detected proteins\n")
  }
  
  cat("\n📚 STRATEGIES SUMMARY:\n")
  cat("1. Complete Case: Use only proteins with no NAs\n")
  cat("2. Filtering: Remove proteins with >X% missing values\n")
  cat("3. Group-wise: Require detection in ≥Y% samples per group\n")
  cat("4. Imputation: Fill NAs with estimated values (KNN, minimum, etc.)\n")
  cat("5. Hybrid: Combine filtering + imputation\n")
  
  cat("\n💡 TIP: Protein NAs ≠ Data Quality Issues\n")
  cat("Missing proteins often reflect:\n")
  cat("• Low abundance proteins below detection limit\n")
  cat("• Sample-specific biology (some proteins not expressed)\n")
  cat("• Normal variation in complex proteomes\n\n")
  
  strategies <- list(
    na_percent = na_percent,
    primary_recommendation = if (na_percent < 15) "complete_case" 
                            else if (na_percent < 40) "imputation_or_filtering"
                            else if (na_percent < 60) "strict_filtering"
                            else "data_quality_review",
    alternative_strategies = c("complete_case", "group_wise_filtering", "imputation", "hybrid")
  )
  
  return(invisible(strategies))
}


# ----------------------------------------------------------------------------
# validatePostImputationProteinData
# ----------------------------------------------------------------------------
#' Validate Post-Imputation Protein Data
#' 
#' @description A simple wrapper to validate protein data after imputation,
#' specifically checking if imputation was successful.
#' 
#' @param protein_obj A ProteinQuantitativeData S4 object (post-imputation)
#' @param expected_na_percent Expected NA percentage (default: varies based on protein data)
#' @param tolerance Tolerance for expected percentage (default: 10%)
#' 
#' @return Logical indicating if validation passed, with detailed output
#' 
#' @export
validatePostImputationProteinData <- function(protein_obj, expected_na_percent = NULL, tolerance = 10) {
  
  cat("\n=== POST-IMPUTATION PROTEIN DATA VALIDATION ===\n")
  cat("Note: Protein-level NAs occur even after peptide imputation because:\n")
  cat("• Proteins need ≥1 detected peptide to get a quantification\n")
  cat("• Some proteins detected only in subset of samples\n")
  cat("• This is normal proteomics data behavior!\n\n")
  
  # Run the full NA analysis
  na_results <- checkProteinNAPercentages(protein_obj, verbose = TRUE)
  
  # Set expected NA percentage if not provided (proteins often have some NAs)
  if (is.null(expected_na_percent)) {
    # For protein data, NAs are very common due to missing peptides/proteins
    # Typical ranges: 20-50% depending on sample complexity and detection method
    expected_na_percent <- 35  # Realistic expectation for protein data
    cat(sprintf("Note: Using default expected NA%% of %.1f%% for protein data\n", expected_na_percent))
    cat("(Protein-level NAs are normal due to incomplete protein detection across samples)\n")
  }
  
  # Check if validation passes
  actual_na_percent <- na_results$total_na_percent
  is_valid <- abs(actual_na_percent - expected_na_percent) <= tolerance
  
  cat("\n--- VALIDATION RESULT ---\n")
  cat(sprintf("Expected NA%%: %.2f%% (± %.2f%%)\n", expected_na_percent, tolerance))
  cat(sprintf("Actual NA%%: %.2f%%\n", actual_na_percent))
  
  if (is_valid) {
    cat("✓ VALIDATION PASSED: Protein data NA levels are within expected range!\n")
  } else {
    cat("✗ VALIDATION FAILED: Unexpected NA percentage detected!\n")
    if (actual_na_percent > expected_na_percent + tolerance) {
      cat("  → Issue: More NAs than expected. Check for missing proteins/peptides.\n")
    } else {
      cat("  → Issue: Fewer NAs than expected. Possible over-imputation.\n")
    }
  }
  
  # Additional warnings for common issues
  if (actual_na_percent > 50) {
    cat("⚠ WARNING: Very high NA percentage (>50%) suggests data quality issues!\n")
  }
  
  if (actual_na_percent < 10) {
    cat("ℹ INFO: Very low NA percentage (<10%) - excellent protein coverage!\n")
  }
  
  # Educational information about protein NAs
  if (actual_na_percent > 20 && actual_na_percent < 50) {
    cat("ℹ INFO: NA percentage is typical for protein-level data\n")
    cat("  → This reflects biological reality: not all proteins detected in all samples\n")
    cat("  → Consider: protein-level imputation OR complete-case analysis\n")
  }
  
  if (na_results$summary_stats$max_na_per_sample > actual_na_percent + 10) {
    cat("⚠ WARNING: Large variation in NA% between samples detected!\n")
    cat("  → Some samples may have much lower protein coverage.\n")
  }
  
  # Check for problematic samples (>80% missing)
  high_missing_samples <- na_results$per_sample_na[na_results$per_sample_na$na_percentage > 80, ]
  if (nrow(high_missing_samples) > 0) {
    cat("⚠ WARNING: Samples with >80% missing proteins detected:\n")
    print(high_missing_samples[, c("sample", "na_percentage")])
  }
  
  cat("\n")
  return(invisible(list(
    is_valid = is_valid,
    actual_na_percent = actual_na_percent,
    expected_na_percent = expected_na_percent,
    full_results = na_results
  )))
}


# ----------------------------------------------------------------------------
# getPairsOfSamplesTable
# ----------------------------------------------------------------------------
#' @title Get Pairs of Samples for Comparison
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


# ----------------------------------------------------------------------------
# calulatePearsonCorrelation
# ----------------------------------------------------------------------------
#' @title Calculate Pearson Correlation
#'@description Calculate the Pearson correlation of the abundances of peptides between two samples X and Y.
#'@param ms_filename_x A string representing the sample file name X (for a pair of sample in the same technical replicate group) for correlation score calculation.
#'@param ms_filename_y A string representing the sample file name Y (for a pair of sample in the same technical replicate group) for correlation score calculation.
#'@param input_table A data frame with the following columns: 1. Sample file name or Run name, 2. Protein IDs, 3. Stripped peptide sequence, 4. Normalised peptide abundances
#'@param sample_id_column Sample ID column name as string (default = "Run").
#'@param protein_id_column Protein accession column name as string (default = "Protein.Ids").
#'@param peptide_sequence_column Peptide sequence column name as string (default = "Stripped.Sequence").
#'@param peptide_normalised_column Normalised peptide abundance column name as string (default = "Peptide.Normalised").
#'@return The pearson correlation value of the abundances of peptides between two samples X and Y.
#' @export
calulatePearsonCorrelation <- function( ms_filename_x, ms_filename_y, input_table
                                        , sample_id_column = "Run"
                                        , protein_id_column = "Protein.Ids"
                                        , peptide_sequence_column = "Stripped.Sequence"
                                        , peptide_normalised_column = "Peptide.Normalised")  {

  tab_x <- input_table |>
    dplyr::filter( !!rlang::sym(sample_id_column) == ms_filename_x )

  tab_y <- input_table |>
    dplyr::filter( !!rlang::sym(sample_id_column) == ms_filename_y )

  merged_tbl <- tab_x |>
    inner_join( tab_y, by=join_by( !!rlang::sym(protein_id_column), !!rlang::sym(peptide_sequence_column)) )

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


# ----------------------------------------------------------------------------
# calculatePearsonCorrelationMatrix
# ----------------------------------------------------------------------------
#' @title Calculate Pearson Correlation Matrix (Optimized)
#' @description Computes all pairwise sample correlations using a single matrix operation.
#' Much faster than iterative approach for large datasets. Pivots data to wide format
#' (proteins x samples) and computes the full correlation matrix in one vectorized call.
#' @param input_table Long-format data table with sample, protein, and value columns
#' @param sample_id_column String name of sample ID column
#' @param protein_id_column String name of protein ID column
#' @param peptide_normalised_column String name of normalized value column
#' @return Correlation matrix (samples x samples)
#' @importFrom tidyr pivot_wider
#' @export
calculatePearsonCorrelationMatrix <- function(input_table, 
                                               sample_id_column,
                                               protein_id_column,
                                               peptide_normalised_column) {
  
  message("*** PEARSON MATRIX: Pivoting data to wide format... ***")
  pivot_start <- Sys.time()
  

  # Pivot to wide: rows = proteins, columns = samples
  wide_matrix <- input_table |>
    tidyr::pivot_wider(
      id_cols = dplyr::all_of(protein_id_column),
      names_from = dplyr::all_of(sample_id_column),
      values_from = dplyr::all_of(peptide_normalised_column),
      values_fn = mean  # Handle any duplicates
    ) |>
    dplyr::select(-dplyr::all_of(protein_id_column)) |>
    as.matrix()
  
  pivot_elapsed <- as.numeric(difftime(Sys.time(), pivot_start, units = "secs"))
  message(sprintf("*** PEARSON MATRIX: Pivot completed in %.2f seconds (%d proteins x %d samples) ***", 
                  pivot_elapsed, nrow(wide_matrix), ncol(wide_matrix)))
  
  message("*** PEARSON MATRIX: Computing correlation matrix... ***")
  cor_start <- Sys.time()
  
  # Single cor() call computes ALL pairwise correlations
  cor_matrix <- cor(wide_matrix, use = "pairwise.complete.obs")
  
  cor_elapsed <- as.numeric(difftime(Sys.time(), cor_start, units = "secs"))
  message(sprintf("*** PEARSON MATRIX: Correlation matrix computed in %.2f seconds ***", cor_elapsed))
  
  cor_matrix
}


# ----------------------------------------------------------------------------
# calculatePearsonCorrelationOptimized
# ----------------------------------------------------------------------------
#' @title Calculate Pearson Correlation (Optimized for pre-filtered data)
#' @description Fast correlation calculation for pre-filtered sample data.
#' This function is optimized to work with pre-partitioned data and avoids
#' redundant filtering operations.
#' @param data_x Pre-filtered data for sample X
#' @param data_y Pre-filtered data for sample Y  
#' @param protein_id_column String name of protein ID column
#' @param peptide_sequence_column String name of peptide sequence column
#' @param peptide_normalised_column String name of normalized peptide column
#' @return Numeric Pearson correlation value, or NA_real_ if insufficient data
#' @export
calculatePearsonCorrelationOptimized <- function(data_x, data_y, 
                                                  protein_id_column,
                                                  peptide_sequence_column, 
                                                  peptide_normalised_column) {
  
  # Quick validation
  if (is.null(data_x) || is.null(data_y) || nrow(data_x) == 0 || nrow(data_y) == 0) {
    return(NA_real_)
  }
  
  # Inner join on protein + peptide (already filtered, just need matching)
  merged_tbl <- data_x |>
    dplyr::inner_join(
      data_y,
      by = c(protein_id_column, peptide_sequence_column),
      suffix = c(".x", ".y")
    )
  
  if (nrow(merged_tbl) == 0) {
    return(NA_real_)
  }
  
  # Extract values
  values_x <- merged_tbl[[paste0(peptide_normalised_column, ".x")]]
  values_y <- merged_tbl[[paste0(peptide_normalised_column, ".y")]]
  
  # Calculate correlation
  if (length(values_x) > 0 && length(values_y) > 0) {
    cor(values_x, values_y, use = "pairwise.complete.obs")
  } else {
    NA_real_
  }
}


# ----------------------------------------------------------------------------
# calulatePearsonCorrelationForSamplePairsHelper
# ----------------------------------------------------------------------------
#' @title Calculate Pearson Correlation for Sample Pairs
#' @description Helper function to calculate pairwise Pearson correlations between samples
#' using optimized sequential processing with pre-partitioned data for fast performance.
#' Data is pre-partitioned by sample ID to eliminate redundant filtering operations.
#' @param samples_id_tbl Table of sample IDs
#' @param run_id_column String name of run ID column (default: "ms_filename")
#' @param replicate_group_column String name of replicate group column (default: "general_sample_info")
#' @param input_table Input data table
#' @param num_of_cores Number of CPU cores (parameter retained for compatibility but not used)
#' @param sample_id_column String name of sample ID column (default: "Run")
#' @param protein_id_column String name of protein ID column (default: "Protein.Ids")
#' @param peptide_sequence_column String name of peptide sequence column (default: "Stripped.Sequence")
#' @param peptide_normalised_column String name of normalized peptide column (default: "Peptide.Normalised")
#' @return Data frame with pairwise correlations
#' @importFrom purrr map2_dbl map_chr
#' @export
calulatePearsonCorrelationForSamplePairsHelper <- function( samples_id_tbl
                                                      , run_id_column = "ms_filename"
                                                      , replicate_group_column = "general_sample_info"
                                                      , input_table
                                                      , num_of_cores = 1
                                                      , sample_id_column = "Run"
                                                      , protein_id_column = "Protein.Ids"
                                                      , peptide_sequence_column = "Stripped.Sequence"
                                                      , peptide_normalised_column = "Peptide.Normalised") {


  pairs_for_comparison <- getPairsOfSamplesTable(samples_id_tbl
                                                 , run_id_column = run_id_column
                                                 , replicate_group_column = replicate_group_column)
  
  # Log pair generation
  num_pairs <- nrow(pairs_for_comparison)
  message(sprintf("*** PEARSON HELPER: Generated %d sample pairs for correlation analysis ***", num_pairs))

  message(sprintf("--- DEBUG66 [calulatePearsonCorrelationForSamplePairsHelper]: Processing %d pairs ---", num_pairs))

  # Track total calculation time
  total_start_time <- Sys.time()
  
  # MATRIX-BASED OPTIMIZATION: Compute all correlations in one vectorized operation
  # This replaces the slow group_split + iterative approach
  message("*** PEARSON HELPER: Using matrix-based correlation (fast vectorized approach) ***")
  
  # Compute full correlation matrix in one call
  message(sprintf("   [calulatePearsonCorrelationForSamplePairsHelper] Step: Computing correlation matrix..."))
  mat_start <- Sys.time()
  cor_matrix <- calculatePearsonCorrelationMatrix(
    input_table = input_table,
    sample_id_column = sample_id_column,
    protein_id_column = protein_id_column,
    peptide_normalised_column = peptide_normalised_column
  )
  mat_end <- Sys.time()
  message(sprintf("   [calulatePearsonCorrelationForSamplePairsHelper] Matrix computed. Dim: %d x %d. Duration: %.2f secs", 
                  nrow(cor_matrix), ncol(cor_matrix), as.numeric(difftime(mat_end, mat_start, units = "secs"))))
  
  # Extract correlations for the specific pairs we need
  message("*** PEARSON HELPER: Extracting correlations for specified pairs... ***")
  extract_start <- Sys.time()
  
  sample_pairs_x <- pairs_for_comparison[[paste0(run_id_column, ".x")]]
  sample_pairs_y <- pairs_for_comparison[[paste0(run_id_column, ".y")]]
  
  message(sprintf("   [calulatePearsonCorrelationForSamplePairsHelper] Extracting %d values from matrix...", length(sample_pairs_x)))

  # Vectorized extraction from correlation matrix
  correlations <- purrr::map2_dbl(
    sample_pairs_x,
    sample_pairs_y,
    \(x, y) {
      # Direct lookup from correlation matrix - O(1) per pair
      if (x %in% colnames(cor_matrix) && y %in% colnames(cor_matrix)) {
        cor_matrix[x, y]
      } else {
        NA_real_
      }
    }
  )
  
  extract_elapsed <- as.numeric(difftime(Sys.time(), extract_start, units = "secs"))
  message(sprintf("*** PEARSON HELPER: Pair extraction completed in %.2f seconds ***", extract_elapsed))
  
  # Add results back to dataframe
  pearson_correlation_per_pair <- pairs_for_comparison |>
    mutate(pearson_correlation = correlations)
  
  # Log completion with timing statistics
  total_elapsed <- as.numeric(difftime(Sys.time(), total_start_time, units = "secs"))
  message(sprintf("*** PEARSON HELPER: Completed %d correlations in %.1f seconds (matrix approach) ***", 
                  num_pairs, total_elapsed))

  pearson_correlation_per_pair

}


# ----------------------------------------------------------------------------
# filterSamplesByPeptideCorrelationThreshold
# ----------------------------------------------------------------------------
#' @title Filter Samples by Peptide Correlation Threshold
#'@description Remove samples which is correlated with any technical replicate samples
#'@param pearson_correlation_per_pair A data frame with the following columns: 1. ID of technical replicate group, 2. sample file name X, 3. sample file name Y, 4. Pearson correlation of the abundances of peptides between sample X and Y.
#'@param peptide_keep_samples_with_min_num_peptides A data frame with the following columns: 1. Sample file name or Run name, 2. Protein IDs, 3. Stripped peptide sequence, 4. Normalised peptide abundances
#'@param min_pearson_correlation_threshold Minimum pearson correlation for a pair of files to be considered to be consistent and kept for further analysis
#'@param filename_column_x Name of column containing the sample file name X (for a pair of sample in the same technical replicate group). Tidyverse column header format, not a string.
#'@param filename_column_y Name of column containing the sample file name Y (for a pair of sample in the same technical replicate group). Tidyverse column header format, not a string.
#'@param correlation_column Name of column containing the Pearson's correlation score between Sample X and Y. Tidyverse column header format, not a string.
#'@param filename_id_column A string indicating the name of column that contains the sample ID or Run ID in the data frame `peptide_keep_samples_with_min_num_peptides`.
#'@return A table without samples that are poorly correlated with the rest of the samples in the technical replicate group. Contains the following columns: 1. Sample file name or Run name, 2. Protein IDs, 3. Stripped peptide sequence, 4. Normalised peptide abundances
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


# ----------------------------------------------------------------------------
# findSamplesPairBelowPeptideCorrelationThreshold
# ----------------------------------------------------------------------------
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


# ----------------------------------------------------------------------------
# filterSamplesByProteinCorrelationThresholdHelper
# ----------------------------------------------------------------------------
#' @title Filter Samples by Protein Correlation
#'@description Remove samples which is correlated with any technical replicate samples
#'@param protein_intensity_table A data frame with the following columns: 1. ID of technical replicate group, 2. sample file name X, 3. sample file name Y, 4. Pearson correlation of the abundances of peptides between sample X and Y.
#'@param peptide_keep_samples_with_min_num_peptides A data frame with the proteins as rows and samples ID as columns.
#'@param min_pearson_correlation_threshold Minimum pearson correlation for a pair of files to be considered to be consistent and kept for further analysis
#'@param filename_column_x Name of column containing the sample file name X (for a pair of sample in the same technical replicate group). Tidyverse column header format, not a string.
#'@param filename_column_y Name of column containing the sample file name Y (for a pair of sample in the same technical replicate group). Tidyverse column header format, not a string.
#'@param protein_id_column Name of column containing the protein ID. Tidyverse column header format, not a string.
#'@param correlation_column Name of column containing the Pearson's correlation score between Sample X and Y. Tidyverse column header format, not a string.
#'@param filename_id_column A string indicating the name of column that contains the sample ID or Run ID in the data frame `peptide_keep_samples_with_min_num_peptides`.
#'@return A table without samples that are poorly correlated with the rest of the samples in the technical replicate group. Contains the following columns: 1. Sample file name or Run name, 2. Protein IDs, 3. Stripped peptide sequence, 4. Normalised peptide abundances
#' @export
filterSamplesByProteinCorrelationThresholdHelper <- function(pearson_correlation_per_pair
                                                       , protein_intensity_table
                                                       , min_pearson_correlation_threshold = 0.75
                                                       , filename_column_x = ms_filename.x
                                                       , filename_column_y = ms_filename.y
                                                       , protein_id_column = Protein.Ids
                                                       , correlation_column = pearson_correlation ) {

  message("--- DEBUG66 [filterSamplesByProteinCorrelationThresholdHelper]: Entry ---")
  message(sprintf("   [filterSamplesByProteinCorrelationThresholdHelper] protein_intensity_table dims: %d x %d", 
                  nrow(protein_intensity_table), ncol(protein_intensity_table)))
  message(sprintf("   [filterSamplesByProteinCorrelationThresholdHelper] Threshold: %.2f", min_pearson_correlation_threshold))

  # All Samples
  all_samples <-  pearson_correlation_per_pair |>
    pivot_longer( cols =c({{filename_column_x}}, {{filename_column_y}})
                  , values_to = "temp_column" ) |>
    dplyr::distinct( temp_column )
  
  message(sprintf("   [filterSamplesByProteinCorrelationThresholdHelper] Total unique samples in pairs: %d", nrow(all_samples)))

  # Samples to keep include all those pairs of samples with correlation score passing threshold
  samples_to_keep <-  pearson_correlation_per_pair |>
    dplyr::filter( {{correlation_column}} >= min_pearson_correlation_threshold) |>
    pivot_longer( cols =c({{filename_column_x}}, {{filename_column_y}})
                  , values_to = "temp_column" ) |>
    dplyr::distinct( temp_column )
  
  message(sprintf("   [filterSamplesByProteinCorrelationThresholdHelper] Samples passing correlation threshold: %d", nrow(samples_to_keep)))

  # Samples to keep anyway
  # Use robust method to get protein ID column name as string to prevent bad_alloc on setdiff
  pid_col_name <- tryCatch({
    if (is.character(protein_id_column)) protein_id_column else rlang::as_string(rlang::ensym(protein_id_column))
  }, error = function(e) rlang::as_string(rlang::ensym(protein_id_column)))
  
  message(sprintf("   [filterSamplesByProteinCorrelationThresholdHelper] Protein ID Column Name resolved to: '%s'", pid_col_name))
  
  # DEBUG: Trace setdiff inputs
  all_cols <- colnames(protein_intensity_table)
  pair_samples <- all_samples |> dplyr::pull(temp_column)
  
  message(sprintf("   [filterSamplesByProteinCorrelationThresholdHelper] Calculating 'samples_to_keep_anyway'..."))
  message(sprintf("      - Total columns: %d", length(all_cols)))
  message(sprintf("      - Pair samples: %d", length(pair_samples)))
  
  cols_not_in_pairs <- setdiff(all_cols, pair_samples)
  message(sprintf("      - Cols NOT in pairs: %d", length(cols_not_in_pairs)))
  
  samples_to_keep_anyway <- setdiff(cols_not_in_pairs, pid_col_name)
  message(sprintf("      - Final 'samples_to_keep_anyway' (excluding ID col): %d", length(samples_to_keep_anyway)))
  
  print(samples_to_keep_anyway)

  # Samples in the table to keep
  samples_to_keep_subset <- colnames(protein_intensity_table)[colnames(protein_intensity_table) %in% (samples_to_keep |> dplyr::pull( temp_column ))]
  message(sprintf("   [filterSamplesByProteinCorrelationThresholdHelper] Samples from pairs to keep: %d", length(samples_to_keep_subset)))

  samples_above_correlation_threshold <- protein_intensity_table |>
    dplyr::select( {{protein_id_column}}, all_of( c(samples_to_keep_anyway, samples_to_keep_subset)))
  
  message(sprintf("   [filterSamplesByProteinCorrelationThresholdHelper] Result dims: %d x %d", 
                  nrow(samples_above_correlation_threshold), ncol(samples_above_correlation_threshold)))

  samples_above_correlation_threshold

}


# ----------------------------------------------------------------------------
# removeProteinWithOnlyOneReplicate
# ----------------------------------------------------------------------------
#---------------------------------------------------------------------------------------
#' @title Remove Proteins with Single Replicate
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


# ----------------------------------------------------------------------------
# calculatePercentMissingPeptidePerReplicate
# ----------------------------------------------------------------------------
#' calculatePercentMissingPeptidePerReplicate
#' @description Calculate percentage of peptides from each sample that is missing and merge with metadata
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


# ----------------------------------------------------------------------------
# calculatePercentMissingProteinPerReplicate
# ----------------------------------------------------------------------------
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


# ----------------------------------------------------------------------------
# getSamplesCorrelationMatrix
# ----------------------------------------------------------------------------
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


# ----------------------------------------------------------------------------
# avgReplicateProteinIntensity
# ----------------------------------------------------------------------------
#' @title Average Replicate Protein Intensity
#' @description
#' Protein average values from replicate samples
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


# ----------------------------------------------------------------------------
# calculatePercentMissingPerProtein
# ----------------------------------------------------------------------------
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


# ----------------------------------------------------------------------------
# calculateMissingValuesPerProteinFishersTest
# ----------------------------------------------------------------------------
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


# ----------------------------------------------------------------------------
# getRowsToKeepList
# ----------------------------------------------------------------------------
## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#' For each experimental group, identify proteins that does have enough number of samples with abundance values.
#'@param input_table An input table with a column containing the row ID and the rest of the columns representing abundance values for each sample.
#'@param cols A tidyselect command to select the columns. This includes the functions starts_with(), ends_with(), contains(), matches(), and num_range()
#'@param design_matrix A data frame with a column containing the sample ID (as per the sample_id param) and the experimental group (as per the group param). Each row as the sample ID as row name in the data frame.
#'@param sample_id The name of the column in design_matrix table that has the sample ID.
#'@param row_id A unique ID for each row of the 'input_table' variable.
#'@param grouping_variable The name of the column in design_matrix table that has the experimental group.
#'@param min_num_samples_per_group An integer representing the minimum number of samples per group.
#'@param abundance_threshold Abundance threshold in which the protein in the sample must be above for it to be considered for inclusion into data analysis.
#'@return A list, the name of each element is the sample ID and each element is a vector containing the protein accessions (e.g. row_id) with enough number of values.
#'@export
getRowsToKeepList <- function(input_table, cols, design_matrix, sample_id, row_id, grouping_variable, min_num_samples_per_group, abundance_threshold) {

  abundance_long <- input_table |>
    pivot_longer(cols = { { cols } },
                 names_to = as_string(as_name(enquo(sample_id))),
                 values_to = "Abundance") |>
    left_join(design_matrix, by = as_string(as_name(enquo(sample_id))))


  count_values_per_group <- abundance_long |>
    mutate(has_value = ifelse(!is.na(Abundance) & Abundance > abundance_threshold, 1, 0)) |>
    group_by({ { row_id } }, { { grouping_variable } }) |>
    summarise(num_values = sum(has_value)) |>
    ungroup()


  kept_rows_temp <- count_values_per_group |>
    dplyr::filter(num_values >= min_num_samples_per_group) |>
    dplyr::select(-num_values) |>
    group_by({ { grouping_variable } }) |>
    nest(data = c({ { row_id } })) |>
    ungroup() |>
    mutate(data = purrr::map(data, \(x){ x[, as_name(enquo(row_id))][[1]] }))


  sample_rows_lists <- kept_rows_temp$data
  names(sample_rows_lists) <- kept_rows_temp[, as_name(enquo(grouping_variable))][[1]]

  return(sample_rows_lists)

}


# ----------------------------------------------------------------------------
# averageValuesFromReplicates
# ----------------------------------------------------------------------------
## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#'@title Average values from replicates
#'@param design_matrix Contains the sample_id column and the average_replicates_id column
#'@export
averageValuesFromReplicates <- function(input_table, design_matrix, group_pattern, row_id, sample_id, average_replicates_id) {

  output_table <- input_table |>
    as.data.frame() |>
    rownames_to_column(  row_id   ) |>
    pivot_longer( cols=matches(group_pattern),
                  names_to = sample_id,
                  values_to = "value") |>
    left_join( design_matrix, by = sample_id ) |>
    group_by( !!rlang::sym(average_replicates_id) ,  !!rlang::sym(row_id) ) |>
    summarise( value = mean(value, na.rm = TRUE)) |>
    ungroup() |>
    pivot_wider( names_from = !!rlang::sym(average_replicates_id),
                 values_from = "value") |>
    column_to_rownames(row_id) |>
    as.matrix()

  return( output_table )
}


# ----------------------------------------------------------------------------
# proteinTechRepCorrelationHelper
# ----------------------------------------------------------------------------
## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#' @title Calculate protein technical replicate correlation
#' @param design_matrix_tech_rep: design matrix with the technical replicates
#' @param data_matrix: input data matrix
#' @param sample_id_column: column name of the sample ID. This is the unique identifier for each sample.
#' @param tech_rep_column: column name of the technical replicates. Technical replicates of the same sample will have the same value.
#' @param tech_rep_num_column: column name of the technical replicate number. This is a unique number for each technical replicate for each sample.
#' @export
proteinTechRepCorrelationHelper <- function( design_matrix_tech_rep, data_matrix
                                             , protein_id_column = "Protein.Ids"
                                             , sample_id_column="Sample_ID", tech_rep_column = "replicates", tech_rep_num_column = "tech_rep_num", tech_rep_remove_regex = "pool" ) {

  tech_reps_list <- design_matrix_tech_rep |> dplyr::pull( !!sym(tech_rep_num_column )) |> unique()

  frozen_protein_matrix_tech_rep <- data_matrix  |>
    as.data.frame() |>
    rownames_to_column(protein_id_column) |>
    pivot_longer( cols=!matches( protein_id_column)
                  , values_to = "log2_intensity"
                  , names_to = sample_id_column) |>
    left_join( design_matrix_tech_rep
               , by = join_by( !!sym(sample_id_column) == !!sym(sample_id_column))) |>
    dplyr::filter( !str_detect(  !!sym(tech_rep_column) , tech_rep_remove_regex ) ) |>
    dplyr::select(!!sym( protein_id_column), !!sym(tech_rep_column), log2_intensity, !!sym(tech_rep_num_column)) |>
    dplyr::filter( !!sym(tech_rep_num_column ) %in% tech_reps_list ) |>
    pivot_wider( id_cols = c(!!sym( protein_id_column), !!sym(tech_rep_column))
                 , names_from = !!sym(tech_rep_num_column)
                 , values_from = log2_intensity) |>
    nest( data=!matches(protein_id_column)) |>
    mutate( data = purrr::map( data, \(x){ x |> column_to_rownames(tech_rep_column)} ) ) |>
    mutate( pearson = purrr::map_dbl( data, \(x){  if( length(which(!is.na(x[,1]))) > 0 & length(which(!is.na(x[,2]))) > 0) { cor(x, use="pairwise.complete.obs")[1,2] } else { NA_real_ }   })) |>
    mutate( spearman = purrr::map_dbl( data, \(x){  if( length(which(!is.na(x[,1]))) > 0 & length(which(!is.na(x[,2]))) > 0) { cor(x, use="pairwise.complete.obs", method="spearman")[1,2] } else { NA_real_ }  }))

  frozen_protein_matrix_tech_rep
}


# ----------------------------------------------------------------------------
# updateProteinFiltering
# ----------------------------------------------------------------------------
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
#' @importFrom ggplot2 ggplot aes geom_bar geom_text labs theme_minimal theme element_text element_blank geom_line geom_point scale_color_manual annotate theme_void geom_boxplot coord_cartesian ggsave
#' @importFrom dplyr bind_rows mutate group_by ungroup %>%
#' @importFrom forcats fct_reorder
#' @importFrom gridExtra arrangeGrob grid.arrange
#' @importFrom methods slotNames new is
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

    # DEBUG66: Memory check
    message("--- DEBUG66 [updateProteinFiltering]: Entry ---")
    message(sprintf("   [updateProteinFiltering] Step Name: %s", step_name))
    message(sprintf("   [updateProteinFiltering] filtering_progress size: %s", format(object.size(filtering_progress), units = "auto")))
    gc()
    
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
    is_protein_quant <- if (methods::is(data, "S4")) {
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
    message("   [updateProteinFiltering] Generating P1 (Protein Count Bar)...")
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
    message("   [updateProteinFiltering] Generating P4 (Proteins Per Run Line)...")
    # First ensure all data frames in the list have consistent column types
    proteins_per_run_list <- lapply(filtering_progress@proteins_per_run, function(df) {
        df$Run <- as.character(df$Run)
        df$n_proteins <- as.numeric(df$n_proteins)
        return(df)
    })
    
    message(sprintf("      [P4] Binding %d data frames...", length(proteins_per_run_list)))
    p4_data <- bind_rows(proteins_per_run_list, .id = "step") 
    message(sprintf("      [P4] Combined data rows: %d", nrow(p4_data)))
    
    p4 <- p4_data |>
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
    
    message("   [updateProteinFiltering] P4 Generated.")
    
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
        message(sprintf("   [updateProteinFiltering] Saving individual plots to %s...", derived_time_dir))
        for (plot_name in names(plot_list)) {
            filename <- file.path(derived_time_dir,
                                sprintf("%s_%s.png", step_name, plot_name))
            message(sprintf("      Saving %s...", plot_name))
            tryCatch({
                ggsave(filename, 
                       plot = plot_list[[plot_name]], 
                       width = 10, 
                       height = 8, 
                       dpi = 300)
            }, error = function(e) message(sprintf("Warning: Failed to save %s: %s", plot_name, e$message)))
        }
    }
    
    # Return/display plots based on return_grid parameter
    if(return_grid) {
        message("   [updateProteinFiltering] Generating final grid...")
        message(sprintf("      Memory before grid: %s", format(sum(gc()[,2]), units="auto")))
        
        tryCatch({
            if (!is_protein_quant || !all(is.na(filtering_progress@total_peptides))) {
                # Create full grid with all plots if peptide data exists
                message("      Combining all 5 plots...")
                grid1 <- gridExtra::arrangeGrob(p1, p2, p3, ncol = 3)
                grid2 <- gridExtra::arrangeGrob(p4, ncol = 1)
                grid3 <- gridExtra::arrangeGrob(p5, ncol = 1)
                
                # Use arrangeGrob to prevent immediate drawing, which might double-render
                grid_plot <- gridExtra::arrangeGrob(
                    grid1,
                    grid2,
                    grid3,
                    heights = c(1, 1, 1)
                )
            } else {
                # For protein_quant_table without peptide data, only show protein plots
                message("      Combining protein plots (p1, p4)...")
                grid_plot <- gridExtra::arrangeGrob(
                    p1,
                    p4,
                    ncol = 1,
                    heights = c(1, 1)
                )
            }
            
            message(sprintf("      Memory after grid creation: %s", format(sum(gc()[,2]), units="auto")))
            
            # Save the grid if derived_time_dir is valid and save_plots is TRUE
            if (save_plots) {
                message("      Saving combined grid plot...")
                filename <- file.path(derived_time_dir,
                                    sprintf("%s_combined_plots.png", step_name))
                ggsave(filename, 
                       plot = grid_plot, 
                       width = 15, 
                       height = if (!is_protein_quant || !all(is.na(filtering_progress@total_peptides))) 18 else 12,
                       dpi = 300)
            }
            
            return(grid_plot)
            
        }, error = function(e) {
            message(sprintf("ERROR in grid generation: %s", e$message))
            return(NULL)
        })
        
    } else {
        # Print each plot individually
        message("   [updateProteinFiltering] Printing individual plots...")
        for(plot_obj in plot_list) { 
            print(plot_obj)
        }
        # Return the list invisibly
        invisible(plot_list)
    }
}

