# ============================================================================
# func_peptide_qc_imputation.R
# ============================================================================
# Purpose: Peptide and protein missing value imputation functions
# 
# This file contains functions for imputing missing values in peptide and
# protein data, including limpa-based imputation and validation functions.
# Functions in this file are used by mod_prot_qc_peptide_impute.R and
# related imputation modules.
#
# Functions to extract here:
# - peptideMissingValueImputation(): S4 method for peptide imputation
# - peptideMissingValueImputationLimpa(): S4 method for limpa peptide imputation
# - peptideMissingValueImputationHelper(): Helper for peptide imputation
# - proteinMissingValueImputationLimpa(): S4 method for limpa protein imputation
# - imputePerCol(): Impute missing values per column
# - validatePostImputationData(): Validate imputed peptide data
# - validatePostImputationProteinData(): Validate imputed protein data
# - Additional imputation helper functions
#
# Dependencies:
# - limpa package
# - func_general_helpers.R (for utility functions)
# ============================================================================

# TODO: Extract the following functions from their current locations:

# Function 1: peptideMissingValueImputation()
# Current location: R/peptideVsSamplesS4Objects.R
# Type: S4 method (exportMethods)
# Description: Imputes missing values in peptide data
# setMethod(f = "peptideMissingValueImputation", ...) {
#   # Extract from R/peptideVsSamplesS4Objects.R
# }

# Function 2: peptideMissingValueImputationLimpa()
# Current location: R/peptideVsSamplesS4Objects.R
# Type: S4 method (exportMethods)
# Description: Imputes missing values using limpa package
# setMethod(f = "peptideMissingValueImputationLimpa", ...) {
#   # Extract from R/peptideVsSamplesS4Objects.R
# }

# Function 3: peptideMissingValueImputationHelper()
# Current location: R/peptideVsSamplesS4Objects.R
# Description: Helper function for peptide imputation
# peptideMissingValueImputationHelper <- function(...) {
#   # Extract from R/peptideVsSamplesS4Objects.R
# }

# Function 4: proteinMissingValueImputationLimpa()
# Current location: R/limpa_functions.R
# Type: S4 method (exportMethods)
# Description: Imputes missing values in protein data using limpa
# setMethod(f = "proteinMissingValueImputationLimpa", signature = "ProteinQuantitativeData", ...) {
#   # Extract from R/limpa_functions.R
# }

# Function 5: proteinMissingValueImputation()
# Current location: R/proteinVsSamplesS4Objects.R
# Description: Imputes missing values in protein data (non-limpa)
# proteinMissingValueImputation <- function(...) {
#   # Extract from R/proteinVsSamplesS4Objects.R
# }

# Function 6: imputePerCol()
# Current location: R/helper_functions.R
# Description: Imputes missing values per column
# imputePerCol <- function(...) {
#   # Extract from R/helper_functions.R
# }

# Function 7: validatePostImputationData()
# Current location: R/peptideVsSamplesS4Objects.R
# Description: Validates peptide data after imputation
# validatePostImputationData <- function(...) {
#   # Extract from R/peptideVsSamplesS4Objects.R
# }

# Function 8: validatePostImputationProteinData()
# Current location: R/proteinVsSamplesS4Objects.R
# Description: Validates protein data after imputation
# validatePostImputationProteinData <- function(...) {
#   # Extract from R/proteinVsSamplesS4Objects.R
# }


# ----------------------------------------------------------------------------
# imputePerCol
# ----------------------------------------------------------------------------
## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#' Data imputation function
#'@param df Data matrix
#'@param width Adjustment factor to the observed standard deviation
#'@param downshift Downshift the mean value by this downshift factor multiplied by the observed standard deviation.
#'@return Data matrix with the missing values from each column replaced with a value randomly sampled from the normal distribution with adjusted mean and standard deviation. The normal distribution parameters are based on the observed distribution of the same column.
#'@export
imputePerCol <- function(temp, width = 0.3, downshift = 1.8) {

  temp[!is.finite(temp)] <- NA

  temp.sd <- width * sd(temp, na.rm = TRUE)   # shrink sd width
  temp.mean <- mean(temp, na.rm = TRUE) -
    downshift * sd(temp, na.rm = TRUE)   # shift mean of imputed values

  n.missing <- sum(is.na(temp))
  temp[is.na(temp)] <- rnorm(n.missing, mean = temp.mean, sd = temp.sd)
  return(temp)
}


# ----------------------------------------------------------------------------
# validatePostImputationData
# ----------------------------------------------------------------------------
#' Validate Post-Imputation Peptide Data
#' 
#' @description A simple wrapper to validate peptide data after imputation,
#' specifically checking if imputation was successful (should show 0% NAs).
#' 
#' @param peptide_obj A PeptideQuantitativeData S4 object (post-imputation)
#' @param expected_na_percent Expected NA percentage (default: 0 for post-imputation)
#' @param tolerance Tolerance for expected percentage (default: 0.1%)
#' 
#' @return Logical indicating if validation passed, with detailed output
#' 
#' @export
validatePostImputationData <- function(peptide_obj, expected_na_percent = 0, tolerance = 0.1) {
  
  cat("\n=== POST-IMPUTATION VALIDATION ===\n")
  
  # Run the full NA analysis
  na_results <- checkPeptideNAPercentages(peptide_obj, verbose = TRUE)
  
  # Check if imputation was successful
  actual_na_percent <- na_results$total_na_percent
  is_valid <- abs(actual_na_percent - expected_na_percent) <= tolerance
  
  cat("\n--- VALIDATION RESULT ---\n")
  cat(sprintf("Expected NA%%: %.2f%% (± %.2f%%)\n", expected_na_percent, tolerance))
  cat(sprintf("Actual NA%%: %.2f%%\n", actual_na_percent))
  
  if (is_valid) {
    cat("✓ VALIDATION PASSED: Imputation appears successful!\n")
  } else {
    cat("✗ VALIDATION FAILED: Unexpected NA percentage detected!\n")
    if (actual_na_percent > expected_na_percent + tolerance) {
      cat("  → Issue: More NAs than expected. Imputation may have failed.\n")
    } else {
      cat("  → Issue: Fewer NAs than expected. Check data integrity.\n")
    }
  }
  
  # Additional warnings for common issues
  if (actual_na_percent > 10) {
    cat("⚠ WARNING: High NA percentage suggests imputation problems!\n")
  }
  
  if (na_results$summary_stats$max_na_per_sample > actual_na_percent + 5) {
    cat("⚠ WARNING: Large variation in NA% between samples detected!\n")
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
# peptideMissingValueImputationHelper
# ----------------------------------------------------------------------------
#' peptideMissingValueImputationHelper
#' @description Perform peptide level missing value imputation
#'@param input_table A data frame with the following columns: 1. Sample file name or Run name, 2. Protein IDs, 3. Stripped peptide sequence, 4. Normalised peptide abundances
#'@param metadata_table A data table with the following columns: 1. the sample file name or run name (as per parameter sample_id_tbl_sample_id_column), 2. The replicate group ID (as per parameter replicate_group_column)
#'@param input_table_sample_id_column The name of the column in the input_table that contained the run information or sample file name as per the input_table parameter (default: Run)
#'@param sample_id_tbl_sample_id_column The name of the column in the input_table that contained the run information or sample file name as per the metadata_table parameter (default: ms_filename)
#'@param replicate_group_column (default: general_sample_info)
#'@param protein_id_column Protein accession column, tidyverse format (default = Protein.Ids).
#'@param peptide_sequence_column Peptide sequence column, tidyverse fromat (default =  Stripped.Sequence).
#'@param quantity_to_impute_column Name of column containing the peptide abundance that needs to be normalised in tidyverse format (default: Peptide.RawQuantity)
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


# ----------------------------------------------------------------------------
# peptideMissingValueImputation
# ----------------------------------------------------------------------------
#'@export
setMethod( f="peptideMissingValueImputation"
           , signature="PeptideQuantitativeData"
           , definition = function( theObject,  imputed_value_column = NULL, proportion_missing_values = NULL, core_utilisation = NULL) {
             peptide_data <- theObject@peptide_data
             raw_quantity_column <- theObject@raw_quantity_column
             sample_id_column <- theObject@sample_id
             replicate_group_column <- theObject@technical_replicate_id
             design_matrix <- theObject@design_matrix


             imputed_value_column <- checkParamsObjectFunctionSimplifyAcceptNull( theObject
                                                           , "imputed_value_column"
                                                           , NULL)

             proportion_missing_values <- checkParamsObjectFunctionSimplifyAcceptNull( theObject
                                                           , "proportion_missing_values"
                                                           , NULL)

             core_utilisation <- checkParamsObjectFunctionSimplify( theObject
                                                           , "core_utilisation"
                                                           , NA)

             theObject <- updateParamInObject(theObject, "imputed_value_column")
             theObject <- updateParamInObject(theObject, "proportion_missing_values")
             theObject <- updateParamInObject(theObject, "core_utilisation")

             peptide_values_imputed <- peptideMissingValueImputationHelper( input_table = peptide_data
                                                                      , metadata_table = design_matrix
                                                                      , quantity_to_impute_column = !!sym( raw_quantity_column )
                                                                      , imputed_value_column = !!sym(imputed_value_column)
                                                                      , core_utilisation = core_utilisation
                                                                      , input_table_sample_id_column = !!sym( sample_id_column)
                                                                      , sample_id_tbl_sample_id_column = !!sym( sample_id_column)
                                                                      , replicate_group_column = !!sym( replicate_group_column)
                                                                      , proportion_missing_values = proportion_missing_values )

             theObject@peptide_data <- peptide_values_imputed

             theObject <- cleanDesignMatrixPeptide(theObject)

             theObject
           })


# ----------------------------------------------------------------------------
# peptideMissingValueImputationLimpa
# ----------------------------------------------------------------------------
#' @export
setMethod(f="peptideMissingValueImputationLimpa"
          , signature="PeptideQuantitativeData"
          , definition = function(theObject, 
                                  imputed_value_column = NULL, 
                                  use_log2_transform = TRUE,
                                  verbose = TRUE,
                                  ensure_matrix = TRUE) {
            
            # Load required packages
            if (!requireNamespace("limpa", quietly = TRUE)) {
              stop("limpa package is required but not installed. Please install it using: BiocManager::install('limpa')")
            }
            
            # Parameter validation and defaults
            imputed_value_column <- checkParamsObjectFunctionSimplifyAcceptNull(
              theObject, "imputed_value_column", "Peptide.Imputed.Limpa"
            )
            
            use_log2_transform <- checkParamsObjectFunctionSimplify(
              theObject, "use_log2_transform", TRUE
            )
            
            verbose <- checkParamsObjectFunctionSimplify(
              theObject, "verbose", TRUE
            )
            
            # Update parameters in object
            theObject <- updateParamInObject(theObject, "imputed_value_column")
            theObject <- updateParamInObject(theObject, "use_log2_transform")
            theObject <- updateParamInObject(theObject, "verbose")
            
            # Ensure peptide matrix is calculated if requested
            if (ensure_matrix && (!"peptide_matrix" %in% slotNames(theObject) || 
                                  is.null(theObject@peptide_matrix) || 
                                  length(theObject@peptide_matrix) == 0)) {
              if (verbose) {
                log_info("Peptide matrix not found. Calculating peptide matrix...")
              }
              theObject <- calcPeptideMatrix(theObject)
            }
            
            # Extract data
            peptide_data <- theObject@peptide_data
            peptide_matrix <- theObject@peptide_matrix
            raw_quantity_column <- theObject@raw_quantity_column
            sample_id_column <- theObject@sample_id
            design_matrix <- theObject@design_matrix
            
            if (verbose) {
              log_info("Starting limpa-based missing value imputation...")
              log_info("Data dimensions: {nrow(peptide_matrix)} peptides x {ncol(peptide_matrix)} samples")
              log_info("Missing value percentage: {round(100 * mean(is.na(peptide_matrix)), 1)}%")
            }
            
            # Prepare data for limpa (peptides as rows, samples as columns)
            # limpa expects log2-transformed data
            y_peptide <- peptide_matrix
            
            # Transform to log2 if requested and data is not already log-transformed
            if (use_log2_transform && !theObject@is_logged_data) {
              if (verbose) {
                log_info("Applying log2 transformation...")
              }
              # Add small constant to avoid log(0)
              y_peptide <- log2(y_peptide + 1)
            } else if (use_log2_transform && theObject@is_logged_data) {
              if (verbose) {
                log_warn("Data already log2 transformed, skipping additional transformation")
              }
              # Data already log2, use as-is
            } else if (!use_log2_transform && !theObject@is_logged_data) {
              if (verbose) {
                log_info("Converting raw intensities to log2 scale for limpa...")
              }
              # limpa expects log2 data, so transform raw data
              y_peptide <- log2(y_peptide + 1)
            } else {
              # !use_log2_transform && theObject@is_logged_data
              if (verbose) {
                log_info("Using existing log2 transformed data (no additional transformation)")
              }
              # Data already log2, use as-is - this is the correct case!
            }
            
            # Check for infinite or NaN values
            if (any(is.infinite(y_peptide) | is.nan(y_peptide), na.rm = TRUE)) {
              if (verbose) {
                log_warn("Infinite or NaN values detected. Replacing with NA...")
              }
              y_peptide[is.infinite(y_peptide) | is.nan(y_peptide)] <- NA
            }
            
            # Estimate Detection Probability Curve
            if (verbose) {
              log_info("Estimating detection probability curve...")
            }
            
            tryCatch({
              dpcfit <- limpa::dpc(y_peptide)
              
              if (verbose) {
                log_info("DPC parameters estimated:")
                log_info("  beta0 (intercept): {round(dpcfit$dpc[1], 4)}")
                log_info("  beta1 (slope): {round(dpcfit$dpc[2], 4)}")
                
                # Interpret the slope
                slope_interpretation <- if (dpcfit$dpc[2] < 0.3) {
                  "nearly random missing"
                } else if (dpcfit$dpc[2] < 0.7) {
                  "moderate intensity-dependent missing"
                } else if (dpcfit$dpc[2] < 1.2) {
                  "strong intensity-dependent missing"
                } else {
                  "very strong intensity-dependent missing (approaching left-censoring)"
                }
                log_info("  Interpretation: {slope_interpretation}")
              }
              
              # Perform row-wise imputation using limpa
              if (verbose) {
                log_info("Performing row-wise imputation using DPC model...")
              }
              
              y_imputed <- limpa::dpcImpute(y_peptide, dpc = dpcfit)
              
              if (verbose) {
                log_info("Imputation completed successfully")
                log_info("No missing values remaining: {!any(is.na(y_imputed$E))}")
              }
              
              # Extract the imputed matrix
              imputed_matrix <- y_imputed$E
              
              # Transform back to original scale if necessary
              if (use_log2_transform && !theObject@is_logged_data) {
                if (verbose) {
                  log_info("Converting back from log2 scale...")
                }
                imputed_matrix <- 2^imputed_matrix - 1
                # Ensure no negative values
                imputed_matrix[imputed_matrix < 0] <- 0
              }
              
              # Convert back to long format and merge with original data
              if (verbose) {
                log_info("Converting imputed data back to original format...")
              }
              
              # Create peptide IDs that match the matrix rownames
              peptide_ids <- rownames(imputed_matrix)
              
              # Convert imputed matrix to long format
              imputed_long <- imputed_matrix |>
                as.data.frame() |>
                tibble::rownames_to_column("peptide_id") |>
                tidyr::pivot_longer(cols = -peptide_id, 
                                   names_to = sample_id_column, 
                                   values_to = imputed_value_column) |>
                tidyr::separate(peptide_id, 
                               into = c(theObject@protein_id_column, theObject@peptide_sequence_column), 
                               sep = "%")
              
              # Merge with original peptide data
              updated_peptide_data <- peptide_data |>
                dplyr::left_join(imputed_long, 
                                by = c(theObject@protein_id_column, 
                                      theObject@peptide_sequence_column,
                                      sample_id_column))
              
              # Update the object
              theObject@peptide_data <- updated_peptide_data
              theObject@peptide_matrix <- imputed_matrix
              
              # Update norm_quantity_column to point to the new imputed column
              # This ensures plotting functions use the final imputed data
              theObject@norm_quantity_column <- imputed_value_column
              
              # Store DPC results in the object for future reference
              if (is.null(theObject@args)) {
                theObject@args <- list()
              }
              theObject@args$limpa_dpc_results <- list(
                dpc_parameters = dpcfit$dpc,  # Numeric vector c(intercept, slope)
                dpc_object = dpcfit,          # Full DPC object (preferred for dpcQuant)
                missing_percentage_before = round(100 * mean(is.na(y_peptide)), 1),
                missing_percentage_after = round(100 * mean(is.na(imputed_matrix)), 1),
                slope_interpretation = slope_interpretation,
                dpc_method = "limpa_dpc",
                # Store the original y_peptide data for recreating DPC plot
                y_peptide_for_dpc = y_peptide
              )
              
              # Clean design matrix
              theObject <- cleanDesignMatrixPeptide(theObject)
              
              if (verbose) {
                log_info("limpa-based imputation completed successfully!")
                log_info("New imputed column: {imputed_value_column}")
                log_info("DPC parameters stored in object@args$limpa_dpc_results")
              }
              
              return(theObject)
              
            }, error = function(e) {
              log_error("Error during limpa imputation: {e$message}")
              stop(paste("limpa imputation failed:", e$message))
            })
          })


# ----------------------------------------------------------------------------
# proteinMissingValueImputationLimpa
# ----------------------------------------------------------------------------
#' Protein-Level Missing Value Imputation using limpa Package
#'
#' This function applies limpa's DPC-based missing value imputation directly to 
#' protein-level quantification data. This is useful when you already have protein
#' quantification data with missing values that need to be handled.
#'
#' @param theObject A ProteinQuantitativeData object with protein-level data
#' @param dpc_results DPC results to use. If NULL, will estimate using dpc_slope
#' @param dpc_slope Default DPC slope to use if no DPC results available (default: 0.8)
#' @param quantified_protein_column Name for the column containing quantified protein values
#' @param verbose Whether to print progress messages. Default is TRUE
#' @param chunk When verbose=TRUE, how often to output progress information (default: 1000)
#'
#' @details
#' This method treats each protein as a separate "feature" and applies DPC-based
#' imputation using limpa's dpcImpute function. This is appropriate when you have
#' protein-level data with missing values that follow intensity-dependent patterns.
#'
#' @return Updated ProteinQuantitativeData object with imputed protein values
#'
#' @export
setMethod(f="proteinMissingValueImputationLimpa"
          , signature="ProteinQuantitativeData"
          , definition = function(theObject, 
                                  dpc_results = NULL,
                                  dpc_slope = 0.8,
                                  quantified_protein_column = NULL,
                                  verbose = TRUE,
                                  chunk = 1000) {
            
            # Load required packages
            if (!requireNamespace("limpa", quietly = TRUE)) {
              stop("limpa package is required but not installed. Please install it using: BiocManager::install('limpa')")
            }
            
            # Parameter validation and defaults
            quantified_protein_column <- if (is.null(quantified_protein_column)) {
              "Protein.Imputed.Limpa"
            } else {
              quantified_protein_column
            }
            
            # Extract data from protein object
            protein_quant_table <- theObject@protein_quant_table
            protein_id_column <- theObject@protein_id_column
            sample_id_column <- theObject@sample_id
            design_matrix <- theObject@design_matrix
            
            if (verbose) {
              log_info("Starting limpa-based protein-level missing value imputation...")
            }
            
            # Convert to matrix format (proteins as rows, samples as columns)
            # First, identify sample columns (exclude protein ID and other metadata)
            sample_columns <- setdiff(colnames(protein_quant_table), 
                                     c(protein_id_column, "description", "gene_name", 
                                       "protein_name", "organism", "length"))
            
            if (verbose) {
              log_info("Converting protein data to matrix format...")
              log_info("Found {length(sample_columns)} sample columns")
            }
            
            # Create protein matrix
            protein_matrix <- protein_quant_table |>
              dplyr::select(all_of(c(protein_id_column, sample_columns))) |>
              tibble::column_to_rownames(protein_id_column) |>
              as.matrix()
            
            if (verbose) {
              log_info("Protein matrix dimensions: {nrow(protein_matrix)} proteins x {ncol(protein_matrix)} samples")
              log_info("Missing value percentage: {round(100 * mean(is.na(protein_matrix)), 1)}%")
            }
            
            # Check if we need log2 transformation
            # Assume if max value > 50, data is not log2 transformed
            max_val <- max(protein_matrix, na.rm = TRUE)
            needs_log_transform <- max_val > 50
            
            if (needs_log_transform) {
              if (verbose) {
                log_info("Converting to log2 scale for limpa (max value: {round(max_val, 2)})...")
              }
              protein_matrix <- log2(protein_matrix + 1)
            } else {
              if (verbose) {
                log_info("Data appears to be log2-scale already (max value: {round(max_val, 2)})")
              }
            }
            
            # Handle infinite or NaN values
            if (any(is.infinite(protein_matrix) | is.nan(protein_matrix), na.rm = TRUE)) {
              if (verbose) {
                log_warn("Infinite or NaN values detected. Replacing with NA...")
              }
              protein_matrix[is.infinite(protein_matrix) | is.nan(protein_matrix)] <- NA
            }
            
            # Get or estimate DPC parameters
            dpc_params <- NULL
            if (!is.null(dpc_results)) {
              dpc_params <- dpc_results
              if (verbose) {
                log_info("Using provided DPC results")
              }
            } else {
              if (verbose) {
                log_info("Estimating DPC from protein data...")
              }
              tryCatch({
                dpcfit <- limpa::dpc(protein_matrix)
                dpc_params <- dpcfit
                if (verbose) {
                  log_info("DPC parameters estimated:")
                  log_info("  beta0 (intercept): {round(dpcfit$dpc[1], 4)}")
                  log_info("  beta1 (slope): {round(dpcfit$dpc[2], 4)}")
                }
              }, error = function(e) {
                if (verbose) {
                  log_warn("DPC estimation failed, using default slope: {dpc_slope}")
                }
                dpc_params <- NULL
              })
            }
            
            # Apply missing value imputation
            if (verbose) {
              log_info("Applying protein-level missing value imputation...")
            }
            
            tryCatch({
              # Apply dpcImpute to protein matrix
              if (!is.null(dpc_params)) {
                imputed_result <- limpa::dpcImpute(protein_matrix, dpc = dpc_params, verbose = verbose, chunk = chunk)
              } else {
                imputed_result <- limpa::dpcImpute(protein_matrix, dpc.slope = dpc_slope, verbose = verbose, chunk = chunk)
              }
              
              if (verbose) {
                log_info("Protein-level imputation completed successfully")
                log_info("No missing values remaining: {!any(is.na(imputed_result$E))}")
              }
              
              # Extract imputed matrix
              imputed_matrix <- imputed_result$E
              
              # Transform back to original scale if necessary
              if (needs_log_transform) {
                if (verbose) {
                  log_info("Converting back from log2 scale...")
                }
                imputed_matrix <- 2^imputed_matrix - 1
                # Ensure no negative values
                imputed_matrix[imputed_matrix < 0] <- 0
              }
              
              # Convert back to long format
              if (verbose) {
                log_info("Converting imputed data back to original format...")
              }
              
              imputed_long <- imputed_matrix |>
                as.data.frame() |>
                tibble::rownames_to_column(protein_id_column) |>
                tidyr::pivot_longer(cols = -all_of(protein_id_column), 
                                   names_to = sample_id_column, 
                                   values_to = quantified_protein_column)
              
              # Merge with original protein data
              updated_protein_data <- protein_quant_table |>
                dplyr::left_join(imputed_long, by = c(protein_id_column, sample_id_column))
              
              # Update the object
              theObject@protein_quant_table <- updated_protein_data
              
              # Store DPC results for future reference
              if (is.null(theObject@args)) {
                theObject@args <- list()
              }
              
              theObject@args$limpa_protein_imputation_results <- list(
                dpc_parameters_used = if (!is.null(dpc_params)) {
                  if (is.list(dpc_params) && !is.null(dpc_params$dpc)) {
                    dpc_params$dpc  # Extract parameters from DPC object
                  } else if (is.numeric(dpc_params)) {
                    dpc_params  # Already numeric parameters  
                  } else {
                    c(NA, dpc_slope)
                  }
                } else {
                  c(NA, dpc_slope)
                },
                dpc_object_used = if (is.list(dpc_params) && !is.null(dpc_params$dpc)) dpc_params else NULL,
                quantified_protein_column = quantified_protein_column,
                missing_percentage_before = round(100 * mean(is.na(protein_matrix)), 1),
                missing_percentage_after = 0,  # DPC imputation produces complete data
                imputation_method = "limpa_dpc_protein_imputation",
                total_proteins_imputed = nrow(imputed_matrix)
              )
              
              if (verbose) {
                log_info("limpa protein-level imputation completed successfully!")
                log_info("New imputed column: {quantified_protein_column}")
                log_info("DPC results stored in object@args$limpa_protein_imputation_results")
              }
              
              return(theObject)
              
            }, error = function(e) {
              log_error(paste("Error during limpa protein imputation:", e$message))
              stop(paste("limpa protein imputation failed:", e$message))
            })
          })

