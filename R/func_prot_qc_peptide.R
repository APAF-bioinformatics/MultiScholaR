# MultiScholaR: Interactive Multi-Omics Analysis
# Copyright (C) 2024-2026 Ignatius Pang, William Klare, and APAF-bioinformatics
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

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
# Current location: R/da_proteins_functions.R
# Description: Checks peptide NA percentages
# checkPeptideNAPercentages <- function(...) {
#   # Extract from R/da_proteins_functions.R
# }

# Function 14: removePeptidesOnlyInHek293()
# Current location: R/da_proteins_functions.R
# Description: Removes peptides only found in HEK293 cells
# removePeptidesOnlyInHek293 <- function(...) {
#   # Extract from R/da_proteins_functions.R
# }

# Function 15: removePeptidesWithoutAbundances()
# Current location: R/da_proteins_functions.R
# Description: Removes peptides without abundance values
# removePeptidesWithoutAbundances <- function(...) {
#   # Extract from R/da_proteins_functions.R
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
# check_case_collision_columns
# ----------------------------------------------------------------------------
#' @title Check for Case-Collision Column Names
#' @description Detects columns in a data frame that differ only by capitalisation
#'   (e.g., both "Group" and "group"). Emits a warning for each collision found.
#' @param df A data frame to check
#' @param df_label A label for the data frame (used in warning messages)
#' @return Invisible NULL. Called for its side-effect (warnings).
#' @keywords internal
#' @export
check_case_collision_columns <- function(df, df_label = "data frame") {
  col_names <- colnames(df)
  lower_names <- tolower(col_names)
  
  # Find groups of columns that have the same lowercase name
  dup_groups <- split(col_names, lower_names)
  collisions <- dup_groups[lengths(dup_groups) > 1]
  
  if (length(collisions) > 0) {
    lapply(collisions, function(collision_set) {
      warning(sprintf(
        "Case-collision detected in %s: columns %s differ only by capitalisation. This may cause unexpected behaviour.",
        df_label, paste(sprintf("'%s'", collision_set), collapse = " and ")
      ), call. = FALSE)
    })
  }
  
  invisible(NULL)
}

# ----------------------------------------------------------------------------
# peptideIntensityFilteringHelper
# ----------------------------------------------------------------------------
#' @export
#' @title Filter Peptides by Intensity and Proportion
#' @description Remove peptide based on the intensity threshold and the proportion of samples below the threshold
peptideIntensityFilteringHelper <- function(input_table = NULL
                                      , design_matrix = NULL
                                      , min_peptide_intensity_threshold = 15
                                      , sample_id_column = "Run"
                                      , grouping_variable = "group"
                                      , groupwise_percentage_cutoff = 1
                                      , max_groups_percentage_cutoff = 50
                                      , protein_id_column = Protein.Ids
                                      , peptide_sequence_column = Stripped.Sequence
                                      , peptide_quantity_column = Peptide.Normalised
                                      , core_utilisation = NA
                                      , ...) {
  
  message("+===========================================================================+")
  message("|  DEBUG66: Entering peptideIntensityFilteringHelper (OPTIMIZED)            |")
  message("+===========================================================================+")
  
  # Handle aliases and extra arguments
  extra_args <- list(...)
  
  # Alias for input_table
  if (is.null(input_table)) {
    if ("input_data" %in% names(extra_args)) {
      input_table <- extra_args[["input_data"]]
    } else {
      stop("peptideIntensityFilteringHelper: input_table or input_data must be provided.")
    }
  }

  # Alias for peptide_quantity_column
  if ("raw_quantity_column" %in% names(extra_args)) {
      peptide_quantity_column <- extra_args[["raw_quantity_column"]]
  }

  # Robust resolution of column names (handles strings, symbols, or rlang patterns)
  resolve_col <- function(x) {
    if (is.null(x)) return(NULL)
    
    # Try capturing expression as name first (handles symbols like Protein.Ids)
    captured <- tryCatch({
      rlang::as_name(rlang::enquo(x))
    }, error = function(e) NULL)
    
    if (!is.null(captured) && nchar(captured) > 0) {
      # Check if the captured name is actually the variable holding the string
      val <- tryCatch(eval(rlang::enquo(x), envir = parent.frame()), error = function(e) NULL)
      if (is.character(val) && length(val) == 1) return(val)
      return(captured)
    }
    
    # Fallback to string evaluation
    val <- tryCatch(eval(rlang::enquo(x), envir = parent.frame()), error = function(e) NULL)
    if (is.character(val) && length(val) == 1) return(val)
    
    return(rlang::as_label(rlang::enquo(x)))
  }

  sample_id_str <- resolve_col(sample_id_column)
  group_var_str <- resolve_col(grouping_variable)
  peptide_quant_str <- resolve_col(peptide_quantity_column)
  protein_id_str <- resolve_col(protein_id_column)
  peptide_seq_str <- resolve_col(peptide_sequence_column)

  message(sprintf("   DEBUG66: Resolved columns: SampleID=%s, Group=%s, Quant=%s, ProteinID=%s, PeptSeq=%s",
                  sample_id_str, group_var_str, peptide_quant_str, protein_id_str, peptide_seq_str))
                  
  message(sprintf("   DEBUG66: Design Matrix Cols: %s", paste(colnames(design_matrix), collapse = ", ")))
  check_case_collision_columns(design_matrix, "design_matrix")
  message(sprintf("   DEBUG66: Input Table Cols: %s", paste(colnames(input_table), collapse = ", ")))

  message(sprintf("   DEBUG66: Args: threshold=%g, group_cutoff=%g%%, max_groups_fail=%g%%", 
                  min_peptide_intensity_threshold, groupwise_percentage_cutoff, max_groups_percentage_cutoff))

  # Case-insensitive column resolution for design matrix
  dm_cols <- colnames(design_matrix)
  if (!group_var_str %in% dm_cols) {
    match_idx <- which(tolower(dm_cols) == tolower(group_var_str))
    if (length(match_idx) == 1) {
      message(sprintf("   DEBUG66: Column '%s' not found in design matrix, using '%s' (case-insensitive match)", 
                      group_var_str, dm_cols[match_idx]))
      group_var_str <- dm_cols[match_idx]
    } else {
      stop(sprintf("Column '%s' not found in design matrix. Available: %s", 
                    group_var_str, paste(dm_cols, collapse = ", ")))
    }
  }
  if (!sample_id_str %in% dm_cols) {
    match_idx <- which(tolower(dm_cols) == tolower(sample_id_str))
    if (length(match_idx) == 1) {
      message(sprintf("   DEBUG66: Column '%s' not found in design matrix, using '%s' (case-insensitive match)", 
                      sample_id_str, dm_cols[match_idx]))
      sample_id_str <- dm_cols[match_idx]
    }
  }

  # Prepare minimal design matrix
  design_matrix_minimal <- design_matrix |>
    dplyr::select(!!rlang::sym(sample_id_str), !!rlang::sym(group_var_str)) |>
    dplyr::mutate(!!rlang::sym(sample_id_str) := as.character(!!rlang::sym(sample_id_str)))

  # Handle long vs wide input
  is_long <- sample_id_str %in% colnames(input_table)
  
  if (is_long) {
    message("   DEBUG66: Input table appears to be in LONG format. Skipping pivot.")
    abundance_long <- input_table |>
      dplyr::mutate(
        !!rlang::sym(sample_id_str) := as.character(!!rlang::sym(sample_id_str)),
        !!rlang::sym(peptide_quant_str) := dplyr::if_else(
          is.nan(!!rlang::sym(peptide_quant_str)), 
          NA_real_, 
          !!rlang::sym(peptide_quant_str)
        )
      ) |>
      dplyr::left_join(design_matrix_minimal, by = sample_id_str)
  } else {
    message("   DEBUG66: Input table appears to be in WIDE format. Pivoting long.")
    abundance_long <- input_table |>
      tidyr::pivot_longer(
        cols = !c(!!rlang::sym(protein_id_str), !!rlang::sym(peptide_seq_str)),
        names_to = sample_id_str,
        values_to = peptide_quant_str
      ) |>
      dplyr::mutate(
        !!rlang::sym(sample_id_str) := as.character(!!rlang::sym(sample_id_str)),
        !!rlang::sym(peptide_quant_str) := dplyr::if_else(
          is.nan(!!rlang::sym(peptide_quant_str)), 
          NA_real_, 
          !!rlang::sym(peptide_quant_str)
        )
      ) |>
      dplyr::left_join(design_matrix_minimal, by = sample_id_str)
  }

  gc()

  # Count samples per group
  count_values_per_group <- abundance_long |>
    dplyr::distinct(!!rlang::sym(sample_id_str), !!rlang::sym(group_var_str)) |>
    dplyr::count(!!rlang::sym(group_var_str), name = "num_per_group")

  # Calculate missing stats per group
  count_percent_missing_per_group <- abundance_long |>
    dplyr::mutate(is_below = dplyr::if_else(
      !is.na(!!rlang::sym(peptide_quant_str)) & 
      !!rlang::sym(peptide_quant_str) >= min_peptide_intensity_threshold, 
      0, 1
    )) |>
    dplyr::group_by(!!rlang::sym(protein_id_str), !!rlang::sym(peptide_seq_str), !!rlang::sym(group_var_str)) |>
    dplyr::summarise(num_observed_above = sum(1 - is_below, na.rm = TRUE), .groups = "drop") |>
    # CRITICAL: Ensure all peptide-group combinations exist. 
    # If a peptide has 0 rows for a group, complete() will add it with num_observed_above = 0.
    tidyr::complete(
      tidyr::nesting(!!rlang::sym(protein_id_str), !!rlang::sym(peptide_seq_str)),
      !!rlang::sym(group_var_str),
      fill = list(num_observed_above = 0)
    ) |>
    dplyr::left_join(count_values_per_group, by = group_var_str) |>
    dplyr::mutate(num_below_per_group = num_per_group - num_observed_above) |>
    dplyr::mutate(perc_below_per_group = num_below_per_group / num_per_group * 100)

  # Identify peptides to remove
  total_num_of_groups <- nrow(count_values_per_group)
  
  peptides_to_remove <- count_percent_missing_per_group |>
    dplyr::group_by(!!rlang::sym(protein_id_str), !!rlang::sym(peptide_seq_str)) |>
    dplyr::summarise(
      num_groups_failed = sum(perc_below_per_group > groupwise_percentage_cutoff, na.rm = TRUE),
      percent_groups_failed = (num_groups_failed / total_num_of_groups) * 100,
      .groups = "drop"
    ) |>
    dplyr::filter(percent_groups_failed > max_groups_percentage_cutoff)
    
  num_removed <- nrow(peptides_to_remove)
  message(sprintf("   Results: Identified %d peptides to remove.", num_removed))

  # Filter original table
  if (num_removed > 0) {
    peptide_normalised_pif_cln <- input_table |>
      dplyr::anti_join(peptides_to_remove, by = c(protein_id_str, peptide_seq_str))
  } else {
    peptide_normalised_pif_cln <- input_table
  }

  # Statistics for logging
  proteins_before <- dplyr::n_distinct(input_table[[protein_id_str]])
  proteins_after <- dplyr::n_distinct(peptide_normalised_pif_cln[[protein_id_str]])
  message(sprintf("   Results: Peptides: %d -> %d. Proteins: %d -> %d.", 
                  nrow(input_table), nrow(peptide_normalised_pif_cln),
                  proteins_before, proteins_after))

  message("+===========================================================================+")
  message("|  DEBUG66: Exiting peptideIntensityFilteringHelper                         |")
  message("+===========================================================================+")
  
  gc()
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
#'@param protein_id_column Protein ID column name
#'@param peptide_sequence_column Peptide sequence column name
#'@param grouping_variable The name of the column in design_matrix table that has the experimental group.
#'@param groupwise_percentage_cutoff The maximum percentage of values below threshold allow in each group for a peptide .
#'@param max_groups_percentage_cutoff The maximum percentage of groups allowed with too many samples with peptide abundance values below threshold.
#'@param abundance_threshold Abundance threshold in which the protein in the sample must be above for it to be considered for inclusion into data analysis.
#'@param abundance_column The name of the column containing the abundance values.
#'@return A filtered data.frame
#'@export
removePeptidesWithMissingValuesPercentHelper <- function(input_table
                                               , design_matrix
                                               , sample_id
                                               , protein_id_column
                                               , peptide_sequence_column
                                               , grouping_variable
                                               , groupwise_percentage_cutoff = 1
                                               , max_groups_percentage_cutoff = 50
                                               , abundance_threshold = 0
                                               , abundance_column = "Abundance") {

  message("+===========================================================================+")
  message("|  DEBUG66: Entering removePeptidesWithMissingValuesPercentHelper (OPTIMIZED)|")
  message("+===========================================================================+")

  # Setup strings and symbols
  resolve_col <- function(x) {
    if (is.null(x)) return(NULL)
    if (is.character(x)) return(x)
    return(rlang::as_name(rlang::enquo(x)))
  }

  sample_id_str <- resolve_col(sample_id)
  group_var_str <- resolve_col(grouping_variable)
  protein_id_str <- resolve_col(protein_id_column)
  peptide_seq_str <- resolve_col(peptide_sequence_column)
  abundance_col_str <- resolve_col(abundance_column)

  message(sprintf("   DEBUG66: Args: threshold=%g, group_cutoff=%g%%, max_groups_fail=%g%%", 
                  abundance_threshold, groupwise_percentage_cutoff, max_groups_percentage_cutoff))
  check_case_collision_columns(design_matrix, "design_matrix")

  # Case-insensitive column resolution for design matrix
  dm_cols <- colnames(design_matrix)
  if (!group_var_str %in% dm_cols) {
    match_idx <- which(tolower(dm_cols) == tolower(group_var_str))
    if (length(match_idx) == 1) {
      message(sprintf("   DEBUG66: Column '%s' not found in design matrix, using '%s' (case-insensitive match)", 
                      group_var_str, dm_cols[match_idx]))
      group_var_str <- dm_cols[match_idx]
    } else {
      stop(sprintf("Column '%s' not found in design matrix. Available: %s", 
                    group_var_str, paste(dm_cols, collapse = ", ")))
    }
  }
  if (!sample_id_str %in% dm_cols) {
    match_idx <- which(tolower(dm_cols) == tolower(sample_id_str))
    if (length(match_idx) == 1) {
      message(sprintf("   DEBUG66: Column '%s' not found in design matrix, using '%s' (case-insensitive match)", 
                      sample_id_str, dm_cols[match_idx]))
      sample_id_str <- dm_cols[match_idx]
    }
  }

  # Prepare minimal design matrix
  design_matrix_minimal <- design_matrix |>
    dplyr::select(!!rlang::sym(sample_id_str), !!rlang::sym(group_var_str)) |>
    dplyr::mutate(!!rlang::sym(sample_id_str) := as.character(!!rlang::sym(sample_id_str)))

  # Handle long vs wide input
  is_long <- sample_id_str %in% colnames(input_table)
  
  if (is_long) {
    message("   DEBUG66: Input table appears to be in LONG format. Skipping pivot.")
    abundance_long <- input_table |>
      dplyr::mutate(
        !!rlang::sym(sample_id_str) := as.character(!!rlang::sym(sample_id_str)),
        !!rlang::sym(abundance_col_str) := dplyr::if_else(
          is.nan(!!rlang::sym(abundance_col_str)), 
          NA_real_, 
          !!rlang::sym(abundance_col_str)
        )
      ) |>
      dplyr::left_join(design_matrix_minimal, by = sample_id_str)
  } else {
    message("   DEBUG66: Input table appears to be in WIDE format. Pivoting long.")
    abundance_long <- input_table |>
      tidyr::pivot_longer(
        cols = !c(!!rlang::sym(protein_id_str), !!rlang::sym(peptide_seq_str)),
        names_to = sample_id_str,
        values_to = abundance_col_str
      ) |>
      dplyr::mutate(
        !!rlang::sym(sample_id_str) := as.character(!!rlang::sym(sample_id_str)),
        !!rlang::sym(abundance_col_str) := dplyr::if_else(
          is.nan(!!rlang::sym(abundance_col_str)), 
          NA_real_, 
          !!rlang::sym(abundance_col_str)
        )
      ) |>
      dplyr::left_join(design_matrix_minimal, by = sample_id_str)
  }

  gc()

  # Count samples per group
  count_values_per_group <- abundance_long |>
    dplyr::distinct(!!rlang::sym(sample_id_str), !!rlang::sym(group_var_str)) |>
    dplyr::count(!!rlang::sym(group_var_str), name = "num_per_group")

  # Calculate missing stats per group
  count_percent_missing_per_group <- abundance_long |>
    dplyr::mutate(is_below = dplyr::if_else(
      !is.na(!!rlang::sym(abundance_col_str)) & 
      !!rlang::sym(abundance_col_str) >= abundance_threshold, 
      0, 1
    )) |>
    dplyr::group_by(!!rlang::sym(protein_id_str), !!rlang::sym(peptide_seq_str), !!rlang::sym(group_var_str)) |>
    dplyr::summarise(num_observed_above = sum(1 - is_below, na.rm = TRUE), .groups = "drop") |>
    # CRITICAL: Ensure all peptide-group combinations exist. 
    # If a peptide has 0 rows for a group, complete() will add it with num_observed_above = 0.
    tidyr::complete(
      tidyr::nesting(!!rlang::sym(protein_id_str), !!rlang::sym(peptide_seq_str)),
      !!rlang::sym(group_var_str),
      fill = list(num_observed_above = 0)
    ) |>
    dplyr::left_join(count_values_per_group, by = group_var_str) |>
    dplyr::mutate(num_below_per_group = num_per_group - num_observed_above) |>
    dplyr::mutate(perc_below_per_group = num_below_per_group / num_per_group * 100)

  # Identify peptides to remove
  total_num_of_groups <- nrow(count_values_per_group)
  
  peptides_to_remove <- count_percent_missing_per_group |>
    dplyr::group_by(!!rlang::sym(protein_id_str), !!rlang::sym(peptide_seq_str)) |>
    dplyr::summarise(
      num_groups_failed = sum(perc_below_per_group > groupwise_percentage_cutoff, na.rm = TRUE),
      percent_groups_failed = (num_groups_failed / total_num_of_groups) * 100,
      .groups = "drop"
    ) |>
    dplyr::filter(percent_groups_failed > max_groups_percentage_cutoff)
    
  num_removed <- nrow(peptides_to_remove)
  message(sprintf("   Results: Identified %d peptides to remove.", num_removed))

  # Filter original table
  if (num_removed > 0) {
    if (is_long) {
      filtered_tbl <- input_table |>
        dplyr::anti_join(peptides_to_remove, by = c(protein_id_str, peptide_seq_str))
    } else {
      # For wide format, convert peptides_to_remove back to IDs to filter
      filtered_tbl <- input_table |>
        dplyr::anti_join(peptides_to_remove, by = c(protein_id_str, peptide_seq_str))
    }
  } else {
    filtered_tbl <- input_table
  }

  message(sprintf("   Results: Peptides: %d -> %d", nrow(input_table), nrow(filtered_tbl)))
  
  message("+===========================================================================+")
  message("|  Exiting removePeptidesWithMissingValuesPercentHelper                     |")
  message("+===========================================================================+")

  gc()
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


# ----------------------------------------------------------------------------
# peptideIntensityFiltering
# ----------------------------------------------------------------------------
#'@export
setMethod( f="peptideIntensityFiltering"
           , signature="PeptideQuantitativeData"
           , definition = function( theObject, 
                                    grouping_variable = NULL, 
                                    groupwise_percentage_cutoff = NULL, 
                                    max_groups_percentage_cutoff = NULL, 
                                    peptides_intensity_cutoff_percentile = NULL, 
                                    core_utilisation = NULL) {
             message("--- Entering peptideIntensityFiltering S4 Method ---")
             
             peptide_data <- theObject@peptide_data
             raw_quantity_column <- theObject@raw_quantity_column
             norm_quantity_column <- theObject@norm_quantity_column

             message("   peptideIntensityFiltering: Extracting input data...")
             message(sprintf("      Arg: raw_quantity_column = %s", raw_quantity_column))
             message(sprintf("      Arg: norm_quantity_column = %s", norm_quantity_column))
             message(sprintf("      Data State (peptide_data): Dims = %d rows, %d cols", nrow(peptide_data), ncol(peptide_data)))
             message(sprintf("      Columns: %s", paste(colnames(peptide_data), collapse = ", ")))

             message("   peptideIntensityFiltering: Resolving parameters...")
             grouping_variable <- checkParamsObjectFunctionSimplify( theObject, "grouping_variable", "group")
             groupwise_percentage_cutoff <- checkParamsObjectFunctionSimplify( theObject, "groupwise_percentage_cutoff", 1)
             max_groups_percentage_cutoff <- checkParamsObjectFunctionSimplify( theObject, "max_groups_percentage_cutoff", 50)
             peptides_intensity_cutoff_percentile <- checkParamsObjectFunctionSimplify( theObject, "peptides_intensity_cutoff_percentile", 1)
             core_utilisation <- checkParamsObjectFunctionSimplify( theObject, "core_utilisation", NA)

             message(sprintf("      Resolved: grouping_variable = %s", grouping_variable))
             message(sprintf("      Resolved: groupwise_cutoff = %g%%, max_groups_fail = %g%%", groupwise_percentage_cutoff, max_groups_percentage_cutoff))
             message(sprintf("      Resolved: intensity_percentile = %g%%", peptides_intensity_cutoff_percentile))

             message("   peptideIntensityFiltering: Updating parameters in S4 object...")
             theObject <- updateParamInObject(theObject, "grouping_variable")
             theObject <- updateParamInObject(theObject, "groupwise_percentage_cutoff")
             theObject <- updateParamInObject(theObject, "max_groups_percentage_cutoff")
             theObject <- updateParamInObject(theObject, "peptides_intensity_cutoff_percentile")
             theObject <- updateParamInObject(theObject, "core_utilisation")

             message("   peptideIntensityFiltering: Calculating intensity threshold...")
             # Get non-missing values for threshold calculation
             valid_values <- peptide_data |> dplyr::pull(!!sym(raw_quantity_column))
             valid_values <- valid_values[!is.na(valid_values) & !is.nan(valid_values) & !is.infinite(valid_values)]
             
             if (length(valid_values) == 0) {
               message("      WARNING: No valid intensity values found for threshold calculation.")
               min_peptide_intensity_threshold <- 0
             } else {
               min_peptide_intensity_threshold <- ceiling( quantile( valid_values, na.rm=TRUE, probs = c(peptides_intensity_cutoff_percentile/100) ))[1]
             }
             message(sprintf("      Calculated min_peptide_intensity_threshold = %g", min_peptide_intensity_threshold))

             message("   peptideIntensityFiltering: Calling helper function...")
             peptide_normalised_pif_cln <- peptideIntensityFilteringHelper( 
                                              input_table = peptide_data
                                              , design_matrix = theObject@design_matrix
                                              , min_peptide_intensity_threshold = min_peptide_intensity_threshold
                                              , sample_id_column = theObject@sample_id
                                              , grouping_variable = grouping_variable
                                              , groupwise_percentage_cutoff = groupwise_percentage_cutoff
                                              , max_groups_percentage_cutoff = max_groups_percentage_cutoff
                                              , protein_id_column = theObject@protein_id_column
                                              , peptide_sequence_column = theObject@peptide_sequence_column
                                              , peptide_quantity_column = raw_quantity_column
                                              , core_utilisation = core_utilisation)

             message(sprintf("   peptideIntensityFiltering: Helper returned %d rows", nrow(peptide_normalised_pif_cln)))

             theObject@peptide_data <- peptide_normalised_pif_cln

             message("   peptideIntensityFiltering: Cleaning design matrix...")
             theObject <- cleanDesignMatrixPeptide(theObject)

             message("--- Exiting peptideIntensityFiltering S4 Method ---")
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
             
             message("--- Entering removePeptidesWithMissingValuesPercent S4 Method ---")
             
             peptide_data <- theObject@peptide_data
             protein_id_column <- theObject@protein_id_column
             peptide_sequence_column <- theObject@peptide_sequence_column
             raw_quantity_column <- theObject@raw_quantity_column
             norm_quantity_column <- theObject@norm_quantity_column
             sample_id <- theObject@sample_id
             design_matrix <- theObject@design_matrix

             message("   removePeptidesWithMissingValuesPercent: Resolving parameters...")
             grouping_variable <- checkParamsObjectFunctionSimplify( theObject, "grouping_variable", "group")
             groupwise_percentage_cutoff <- checkParamsObjectFunctionSimplify( theObject, "groupwise_percentage_cutoff", 50)
             max_groups_percentage_cutoff <- checkParamsObjectFunctionSimplify( theObject, "max_groups_percentage_cutoff", 50)
             peptides_intensity_cutoff_percentile <- checkParamsObjectFunctionSimplify( theObject, "peptides_intensity_cutoff_percentile", 1)

             message(sprintf("      Resolved: grouping_variable = %s", grouping_variable))
             message(sprintf("      Resolved: groupwise_cutoff = %g%%, max_groups_fail = %g%%", groupwise_percentage_cutoff, max_groups_percentage_cutoff))
             message(sprintf("      Resolved: intensity_percentile = %g%%", peptides_intensity_cutoff_percentile))

             theObject <- updateParamInObject(theObject, "grouping_variable")
             theObject <- updateParamInObject(theObject, "groupwise_percentage_cutoff")
             theObject <- updateParamInObject(theObject, "max_groups_percentage_cutoff")
             theObject <- updateParamInObject(theObject, "peptides_intensity_cutoff_percentile")

             message("   removePeptidesWithMissingValuesPercent: Calculating intensity threshold...")
             # Filter out non-numeric/invalid values for threshold calculation
             valid_values <- peptide_data |> 
               dplyr::pull(!!sym(norm_quantity_column))
             valid_values <- valid_values[!is.na(valid_values) & !is.nan(valid_values) & !is.infinite(valid_values)]
             
             if (length(valid_values) == 0) {
               message("      WARNING: No valid intensity values found for threshold calculation.")
               min_peptide_intensity_threshold <- 0
             } else {
               min_peptide_intensity_threshold <- ceiling( quantile( valid_values, na.rm=TRUE, probs = c(peptides_intensity_cutoff_percentile/100) ))[1]
             }
             message(sprintf("      Calculated min_peptide_intensity_threshold = %g", min_peptide_intensity_threshold))

             message("   removePeptidesWithMissingValuesPercent: Calling helper function...")
             theObject@peptide_data <- removePeptidesWithMissingValuesPercentHelper( 
                                                 input_table = peptide_data
                                               , design_matrix = design_matrix
                                               , sample_id = !!sym(sample_id)
                                               , protein_id_column = !!sym(protein_id_column)
                                               , peptide_sequence_column = !!sym(peptide_sequence_column)
                                               , grouping_variable = !!sym(grouping_variable)
                                               , groupwise_percentage_cutoff = groupwise_percentage_cutoff
                                               , max_groups_percentage_cutoff = max_groups_percentage_cutoff
                                               , abundance_threshold = min_peptide_intensity_threshold
                                               , abundance_column = !!sym(norm_quantity_column) )


             message("   removePeptidesWithMissingValuesPercent: Cleaning design matrix...")
             theObject <- cleanDesignMatrixPeptide(theObject)

             message("--- Exiting removePeptidesWithMissingValuesPercent S4 Method ---")
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

