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
  cat(sprintf("Expected NA%%: %.2f%% (+/- %.2f%%)\n", expected_na_percent, tolerance))
  cat(sprintf("Actual NA%%: %.2f%%\n", actual_na_percent))
  
  if (is_valid) {
    cat("[OK] VALIDATION PASSED: Imputation appears successful!\n")
  } else {
    cat("[FAIL] VALIDATION FAILED: Unexpected NA percentage detected!\n")
    if (actual_na_percent > expected_na_percent + tolerance) {
      cat("  -> Issue: More NAs than expected. Imputation may have failed.\n")
    } else {
      cat("  -> Issue: Fewer NAs than expected. Check data integrity.\n")
    }
  }
  
  # Additional warnings for common issues
  if (actual_na_percent > 10) {
    cat("[WARNING] WARNING: High NA percentage suggests imputation problems!\n")
  }
  
  if (na_results$summary_stats$max_na_per_sample > actual_na_percent + 5) {
    cat("[WARNING] WARNING: Large variation in NA% between samples detected!\n")
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
#' @description Impute an eligible missing peptide measurement with the mean of
#'   the observed raw quantities from the same declared technical-replicate
#'   group. Replicate availability is counted by distinct run identity.
#' @param input_table Long-format peptide measurements.
#' @param metadata_table Experimental design containing one row per run.
#' @param input_table_sample_id_column Run column in `input_table`.
#' @param sample_id_tbl_sample_id_column Run column in `metadata_table`.
#' @param replicate_group_column Technical-replicate group column in the design.
#' @param protein_id_column Active protein-group identity column.
#' @param peptide_sequence_column Stripped peptide sequence column.
#' @param quantity_to_impute_column Raw peptide quantity used for the mean.
#' @param imputed_value_column Output quantity column.
#' @param hek_string Deprecated opt-in regular expression. When supplied, runs
#'   whose run ID or replicate-group ID matches are excluded from imputation.
#'   The default is `NULL`; names containing `HEK` are not implicitly excluded.
#' @param proportion_missing_values Maximum eligible missing fraction, inclusive.
#' @param core_utilisation Optional historical parallel backend. Validation and
#'   imputation are performed deterministically in-process.
#' @param exclusion_column Optional design column containing an explicit logical
#'   (or unambiguous 0/1, true/false) run-exclusion flag.
#' @param return_imputation_result Return data plus support and audit metadata
#'   when `TRUE`; retain the historical data-frame return when `FALSE`.
#' @export
peptideMissingValueImputationHelper <- function(input_table,
                                                 metadata_table,
                                                 input_table_sample_id_column = Run,
                                                 sample_id_tbl_sample_id_column = ms_filename,
                                                 replicate_group_column = general_sample_info,
                                                 protein_id_column = Protein.Ids,
                                                 peptide_sequence_column = Stripped.Sequence,
                                                 quantity_to_impute_column = Peptide.Normalised,
                                                 imputed_value_column = Peptide.Imputed,
                                                 hek_string = NULL,
                                                 proportion_missing_values = 0.50,
                                                 core_utilisation = NA,
                                                 exclusion_column = NULL,
                                                 return_imputation_result = FALSE) {
  input_sample_column <- rlang::as_name(rlang::enquo(input_table_sample_id_column))
  design_sample_column <- rlang::as_name(rlang::enquo(sample_id_tbl_sample_id_column))
  replicate_column <- rlang::as_name(rlang::enquo(replicate_group_column))
  protein_column <- rlang::as_name(rlang::enquo(protein_id_column))
  peptide_column <- rlang::as_name(rlang::enquo(peptide_sequence_column))
  quantity_column <- rlang::as_name(rlang::enquo(quantity_to_impute_column))
  imputed_column <- rlang::as_name(rlang::enquo(imputed_value_column))

  if (!is.data.frame(input_table) || !is.data.frame(metadata_table)) {
    stop("peptideMissingValueImputationHelper: input and metadata must be data frames.", call. = FALSE)
  }
  required_input <- c(input_sample_column, protein_column, peptide_column, quantity_column)
  required_design <- c(design_sample_column, replicate_column)
  missing_input <- setdiff(required_input, names(input_table))
  missing_design <- setdiff(required_design, names(metadata_table))
  if (length(missing_input) > 0L || length(missing_design) > 0L) {
    stop(
      paste0(
        "peptideMissingValueImputationHelper: missing required column(s): ",
        paste(c(missing_input, missing_design), collapse = ", "),
        "."
      ),
      call. = FALSE
    )
  }
  if (identical(quantity_column, imputed_column)) {
    stop(
      "peptideMissingValueImputationHelper: raw and imputed quantity columns must be different.",
      call. = FALSE
    )
  }
  maximum_missing <- suppressWarnings(as.numeric(proportion_missing_values))
  if (length(maximum_missing) != 1L || is.na(maximum_missing) ||
      !is.finite(maximum_missing) || maximum_missing < 0 || maximum_missing > 1) {
    stop(
      "peptideMissingValueImputationHelper: `proportion_missing_values` must be between 0 and 1.",
      call. = FALSE
    )
  }
  if (!is.numeric(input_table[[quantity_column]])) {
    stop(
      sprintf("peptideMissingValueImputationHelper: quantity column `%s` must be numeric.", quantity_column),
      call. = FALSE
    )
  }
  non_finite_rows <- which(!is.na(input_table[[quantity_column]]) &
                             !is.finite(input_table[[quantity_column]]))
  if (length(non_finite_rows) > 0L) {
    stop(
      sprintf(
        "peptideMissingValueImputationHelper: quantity column `%s` contains non-finite observed values at row(s) %s.",
        quantity_column,
        paste(utils::head(non_finite_rows, 10L), collapse = ", ")
      ),
      call. = FALSE
    )
  }

  internal_columns <- c(
    ".impute_run", ".impute_group", ".impute_protein", ".impute_peptide",
    ".impute_quantity", ".impute_source_order", ".impute_design_order",
    ".impute_excluded", ".impute_group_size"
  )
  collisions <- intersect(internal_columns, union(names(input_table), names(metadata_table)))
  if (length(collisions) > 0L) {
    stop(
      paste0(
        "peptideMissingValueImputationHelper: reserved internal column(s) present: ",
        paste(collisions, collapse = ", "),
        "."
      ),
      call. = FALSE
    )
  }

  finish_result <- function(data, support_table, summary) {
    result <- list(data = data, support_table = support_table, summary = summary)
    if (isTRUE(return_imputation_result)) result else data
  }

  annotate_without_imputation <- function(data) {
    data[[imputed_column]] <- data[[quantity_column]]
    data$is_imputed <- rep(FALSE, nrow(data))
    data
  }

  if (nrow(input_table) == 0L) {
    empty_imputation <- annotate_without_imputation(input_table)
    return(finish_result(
      empty_imputation,
      data.frame(),
      list(
        status = "skipped",
        skip_reason = "empty_input",
        quantity_column = quantity_column,
        imputed_value_column = imputed_column,
        eligibility_denominator = "distinct_design_runs_in_technical_replicate_group",
        maximum_missing_fraction = maximum_missing,
        eligibility_operator = "<=",
        exclusion_source = "none",
        input_rows = 0L,
        output_rows = 0L,
        imputed_rows = 0L
      )
    ))
  }

  input_runs <- as.character(input_table[[input_sample_column]])
  design_runs <- as.character(metadata_table[[design_sample_column]])
  design_groups <- as.character(metadata_table[[replicate_column]])
  missing_run_rows <- which(is.na(input_runs) | !nzchar(trimws(input_runs)))
  missing_design_run_rows <- which(is.na(design_runs) | !nzchar(trimws(design_runs)))
  missing_group_rows <- which(is.na(design_groups) | !nzchar(trimws(design_groups)))
  if (length(missing_run_rows) > 0L || length(missing_design_run_rows) > 0L ||
      length(missing_group_rows) > 0L) {
    stop(
      paste0(
        "peptideMissingValueImputationHelper: missing run/group identity at row(s): input=",
        paste(utils::head(missing_run_rows, 10L), collapse = ","),
        "; design_run=", paste(utils::head(missing_design_run_rows, 10L), collapse = ","),
        "; design_group=", paste(utils::head(missing_group_rows, 10L), collapse = ","),
        "."
      ),
      call. = FALSE
    )
  }
  duplicate_design_runs <- unique(design_runs[
    duplicated(design_runs) | duplicated(design_runs, fromLast = TRUE)
  ])
  if (length(duplicate_design_runs) > 0L) {
    stop(
      sprintf(
        "peptideMissingValueImputationHelper: design must contain one row per run; duplicate run ID(s): %s.",
        paste(utils::head(duplicate_design_runs, 10L), collapse = ", ")
      ),
      call. = FALSE
    )
  }
  unmatched_runs <- setdiff(unique(input_runs), design_runs)
  if (length(unmatched_runs) > 0L) {
    stop(
      sprintf(
        "peptideMissingValueImputationHelper: unmapped input run ID(s): %s.",
        paste(utils::head(unmatched_runs, 10L), collapse = ", ")
      ),
      call. = FALSE
    )
  }

  missing_feature_rows <- which(
    is.na(input_table[[protein_column]]) |
      !nzchar(trimws(as.character(input_table[[protein_column]]))) |
      is.na(input_table[[peptide_column]]) |
      !nzchar(trimws(as.character(input_table[[peptide_column]])))
  )
  if (length(missing_feature_rows) > 0L) {
    stop(
      sprintf(
        "peptideMissingValueImputationHelper: protein/peptide identity is missing at row(s) %s.",
        paste(utils::head(missing_feature_rows, 10L), collapse = ", ")
      ),
      call. = FALSE
    )
  }
  input_keys <- data.frame(
    run = input_runs,
    protein = as.character(input_table[[protein_column]]),
    peptide = as.character(input_table[[peptide_column]]),
    stringsAsFactors = FALSE
  )
  duplicate_keys <- duplicated(input_keys) | duplicated(input_keys, fromLast = TRUE)
  if (any(duplicate_keys)) {
    offending <- unique(input_keys[duplicate_keys, , drop = FALSE])
    offending_labels <- paste(offending$run, offending$protein, offending$peptide, sep = "/")
    stop(
      sprintf(
        paste0(
          "peptideMissingValueImputationHelper: duplicate run-feature observation(s) ",
          "must be resolved before imputation: %s."
        ),
        paste(utils::head(offending_labels, 10L), collapse = ", ")
      ),
      call. = FALSE
    )
  }

  coerce_exclusion_flag <- function(values, column_name) {
    if (is.logical(values)) {
      parsed <- values
    } else if (is.numeric(values)) {
      parsed <- ifelse(values %in% c(0, 1), as.logical(values), NA)
    } else {
      normalised <- tolower(trimws(as.character(values)))
      parsed <- rep(NA, length(normalised))
      parsed[normalised %in% c("true", "t", "yes", "y", "1")] <- TRUE
      parsed[normalised %in% c("false", "f", "no", "n", "0")] <- FALSE
    }
    if (any(is.na(parsed))) {
      stop(
        sprintf(
          "peptideMissingValueImputationHelper: exclusion column `%s` must contain only explicit true/false values.",
          column_name
        ),
        call. = FALSE
      )
    }
    as.logical(parsed)
  }

  exclusion_sources <- character(0)
  excluded <- rep(FALSE, nrow(metadata_table))
  if (!is.null(exclusion_column) && length(exclusion_column) > 0L &&
      !is.na(exclusion_column[[1L]]) && nzchar(trimws(as.character(exclusion_column[[1L]])))) {
    exclusion_column <- trimws(as.character(exclusion_column[[1L]]))
    if (!exclusion_column %in% names(metadata_table)) {
      stop(
        sprintf("peptideMissingValueImputationHelper: exclusion column `%s` is absent from the design.", exclusion_column),
        call. = FALSE
      )
    }
    excluded <- excluded | coerce_exclusion_flag(metadata_table[[exclusion_column]], exclusion_column)
    exclusion_sources <- c(exclusion_sources, paste0("design_column:", exclusion_column))
  } else {
    exclusion_column <- NULL
  }
  if (!is.null(hek_string) && length(hek_string) > 0L &&
      !is.na(hek_string[[1L]]) && nzchar(as.character(hek_string[[1L]]))) {
    legacy_pattern <- as.character(hek_string[[1L]])
    warning(
      paste0(
        "`hek_string` is deprecated; use an explicit design `exclusion_column`. ",
        "The supplied pattern is being applied to run and replicate-group IDs."
      ),
      call. = FALSE
    )
    legacy_excluded <- tryCatch(
      grepl(legacy_pattern, design_runs) | grepl(legacy_pattern, design_groups),
      error = function(error) {
        stop(
          paste0("peptideMissingValueImputationHelper: invalid `hek_string` regex: ", error$message),
          call. = FALSE
        )
      }
    )
    excluded <- excluded | legacy_excluded
    exclusion_sources <- c(exclusion_sources, "deprecated_hek_string")
  }
  exclusion_source <- if (length(exclusion_sources) == 0L) {
    "none"
  } else {
    paste(exclusion_sources, collapse = "+")
  }

  design_map <- data.frame(
    .impute_run = design_runs,
    .impute_group = design_groups,
    .impute_excluded = excluded,
    .impute_design_order = seq_len(nrow(metadata_table)),
    stringsAsFactors = FALSE
  )
  included_design <- design_map[!design_map$.impute_excluded, , drop = FALSE]
  included_design <- included_design |>
    dplyr::group_by(.impute_group) |>
    dplyr::mutate(.impute_group_size = dplyr::n_distinct(.impute_run)) |>
    dplyr::ungroup()

  working <- input_table
  working$.impute_run <- input_runs
  working$.impute_protein <- input_keys$protein
  working$.impute_peptide <- input_keys$peptide
  working$.impute_quantity <- working[[quantity_column]]
  working$.impute_source_order <- seq_len(nrow(working))

  if (nrow(included_design) == 0L || !any(included_design$.impute_group_size >= 2L)) {
    skipped_data <- annotate_without_imputation(input_table)
    return(finish_result(
      skipped_data,
      data.frame(),
      list(
        status = "skipped",
        skip_reason = "no_eligible_technical_replicate_groups",
        quantity_column = quantity_column,
        imputed_value_column = imputed_column,
        protein_identity_column = protein_column,
        peptide_identity_column = peptide_column,
        sample_id_column = input_sample_column,
        replicate_group_column = replicate_column,
        eligibility_denominator = "distinct_design_runs_in_technical_replicate_group",
        maximum_missing_fraction = maximum_missing,
        eligibility_operator = "<=",
        exclusion_source = exclusion_source,
        excluded_runs = design_map$.impute_run[design_map$.impute_excluded],
        input_rows = nrow(input_table),
        output_rows = nrow(skipped_data),
        imputed_rows = 0L,
        technical_replicate_groups = 0L
      )
    ))
  }

  if (inherits(core_utilisation, "multidplyr_cluster")) {
    warning(
      paste0(
        "peptideMissingValueImputationHelper: validated run-feature imputation ",
        "is evaluated in-process for deterministic duplicate and provenance semantics."
      ),
      call. = FALSE
    )
  }

  mapped_input <- dplyr::inner_join(
    working,
    included_design,
    by = ".impute_run"
  )
  feature_groups <- mapped_input |>
    dplyr::distinct(.impute_group, .impute_protein, .impute_peptide)
  support_grid <- feature_groups |>
    dplyr::inner_join(
      included_design,
      by = ".impute_group",
      relationship = "many-to-many"
    ) |>
    dplyr::left_join(
      working,
      by = c(".impute_run", ".impute_protein", ".impute_peptide")
    )

  support_internal <- support_grid |>
    dplyr::group_by(.impute_group, .impute_protein, .impute_peptide) |>
    dplyr::summarise(
      total_distinct_runs = dplyr::n_distinct(.impute_run),
      observed_distinct_runs = dplyr::n_distinct(.impute_run[!is.na(.impute_quantity)]),
      missing_distinct_runs = total_distinct_runs - observed_distinct_runs,
      missing_fraction = missing_distinct_runs / total_distinct_runs,
      mean_observed_raw_quantity = if (observed_distinct_runs > 0L) {
        mean(.impute_quantity, na.rm = TRUE)
      } else {
        NA_real_
      },
      group_has_technical_replicates = total_distinct_runs >= 2L,
      eligible_for_imputation = group_has_technical_replicates &
        observed_distinct_runs > 0L & missing_distinct_runs > 0L &
        missing_fraction <= maximum_missing,
      .groups = "drop"
    )

  completed_rows <- support_grid |>
    dplyr::filter(.impute_group_size >= 2L) |>
    dplyr::left_join(
      support_internal,
      by = c(".impute_group", ".impute_protein", ".impute_peptide")
    )
  completed_rows[[input_sample_column]] <- completed_rows$.impute_run
  completed_rows[[protein_column]] <- completed_rows$.impute_protein
  completed_rows[[peptide_column]] <- completed_rows$.impute_peptide
  completed_rows$is_imputed <- is.na(completed_rows$.impute_quantity) &
    completed_rows$eligible_for_imputation
  completed_rows[[imputed_column]] <- completed_rows$.impute_quantity
  completed_rows[[imputed_column]][completed_rows$is_imputed] <-
    completed_rows$mean_observed_raw_quantity[completed_rows$is_imputed]

  technical_run_ids <- included_design$.impute_run[included_design$.impute_group_size >= 2L]
  passthrough_rows <- dplyr::inner_join(working, design_map, by = ".impute_run") |>
    dplyr::filter(.impute_excluded | !.impute_run %in% technical_run_ids)
  passthrough_rows[[imputed_column]] <- passthrough_rows$.impute_quantity
  passthrough_rows$is_imputed <- rep(FALSE, nrow(passthrough_rows))

  output_columns <- unique(c(names(input_table), imputed_column, "is_imputed"))
  ordering_columns <- c(".impute_source_order", ".impute_design_order")
  completed_output <- completed_rows[, c(ordering_columns, output_columns), drop = FALSE]
  passthrough_output <- passthrough_rows[, c(ordering_columns, output_columns), drop = FALSE]
  existing_output <- dplyr::bind_rows(
    completed_output[!is.na(completed_output$.impute_source_order), , drop = FALSE],
    passthrough_output
  ) |>
    dplyr::arrange(.impute_source_order)
  generated_output <- completed_output[is.na(completed_output$.impute_source_order), , drop = FALSE] |>
    dplyr::arrange(.impute_design_order)
  output <- dplyr::bind_rows(existing_output, generated_output) |>
    dplyr::select(-dplyr::all_of(ordering_columns))

  support_table <- support_internal |>
    dplyr::transmute(
      technical_replicate_group = .impute_group,
      protein_identity = .impute_protein,
      peptide_identity = .impute_peptide,
      total_distinct_runs,
      observed_distinct_runs,
      missing_distinct_runs,
      missing_fraction,
      mean_observed_raw_quantity,
      group_has_technical_replicates,
      eligible_for_imputation
    )
  imputed_rows <- sum(output$is_imputed)
  summary <- list(
    status = "applied",
    skip_reason = NULL,
    quantity_column = quantity_column,
    imputed_value_column = imputed_column,
    protein_identity_column = protein_column,
    peptide_identity_column = peptide_column,
    sample_id_column = input_sample_column,
    replicate_group_column = replicate_column,
    eligibility_denominator = "distinct_design_runs_in_technical_replicate_group",
    maximum_missing_fraction = maximum_missing,
    eligibility_operator = "<=",
    exclusion_source = exclusion_source,
    exclusion_column = exclusion_column,
    excluded_runs = design_map$.impute_run[design_map$.impute_excluded],
    input_rows = nrow(input_table),
    output_rows = nrow(output),
    observed_rows = sum(!is.na(output[[quantity_column]])),
    missing_not_imputed_rows = sum(is.na(output[[imputed_column]])),
    generated_rows = nrow(output) - nrow(input_table),
    imputed_rows = as.integer(imputed_rows),
    technical_replicate_groups = dplyr::n_distinct(
      included_design$.impute_group[included_design$.impute_group_size >= 2L]
    ),
    eligible_group_features = sum(support_table$eligible_for_imputation)
  )

  finish_result(output, support_table, summary)
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


# ----------------------------------------------------------------------------
# peptideMissingValueImputation
# ----------------------------------------------------------------------------
#'@export
setMethod( f="peptideMissingValueImputation"
           , signature="PeptideQuantitativeData"
           , definition = function(theObject,
                                   imputed_value_column = NULL,
                                   proportion_missing_values = NULL,
                                   core_utilisation = NULL,
                                   exclusion_column = NULL,
                                   hek_string = NULL) {
             section <- theObject@args$peptideMissingValueImputation
             if (!is.list(section)) {
               section <- list()
             }
             resolve_section_value <- function(explicit_value, parameter_name, default_value = NULL) {
               if (!is.null(explicit_value)) {
                 return(explicit_value)
               }
               if (parameter_name %in% names(section) && !is.null(section[[parameter_name]])) {
                 return(section[[parameter_name]])
               }
               default_value
             }

             imputed_value_column <- resolve_section_value(
               imputed_value_column,
               "imputed_value_column",
               "Peptide.Imputed"
             )
             proportion_missing_values <- resolve_section_value(
               proportion_missing_values,
               "proportion_missing_values",
               0.50
             )
             core_utilisation <- resolve_section_value(core_utilisation, "core_utilisation", NA)
             exclusion_column <- resolve_section_value(exclusion_column, "exclusion_column", NULL)
             hek_string <- resolve_section_value(hek_string, "hek_string", NULL)
             if (length(exclusion_column) == 0L || is.na(exclusion_column[[1L]]) ||
                 !nzchar(trimws(as.character(exclusion_column[[1L]])))) {
               exclusion_column <- NULL
             } else {
               exclusion_column <- trimws(as.character(exclusion_column[[1L]]))
             }

             replicate_group_column <- theObject@technical_replicate_id
             has_declared_replicates <- length(replicate_group_column) == 1L &&
               !is.na(replicate_group_column) && nzchar(trimws(replicate_group_column)) &&
               replicate_group_column %in% names(theObject@design_matrix)

             if (!has_declared_replicates) {
               peptide_values_imputed <- theObject@peptide_data
               peptide_values_imputed[[imputed_value_column]] <-
                 peptide_values_imputed[[theObject@raw_quantity_column]]
               peptide_values_imputed$is_imputed <- rep(FALSE, nrow(peptide_values_imputed))
               filter_result <- list(
                 data = peptide_values_imputed,
                 support_table = data.frame(),
                 summary = list(
                   status = "skipped",
                   skip_reason = "technical_replicate_column_not_declared",
                   quantity_column = theObject@raw_quantity_column,
                   imputed_value_column = imputed_value_column,
                   eligibility_denominator = "distinct_design_runs_in_technical_replicate_group",
                   maximum_missing_fraction = as.numeric(proportion_missing_values),
                   eligibility_operator = "<=",
                   exclusion_source = if (is.null(exclusion_column)) "none" else paste0("design_column:", exclusion_column),
                   input_rows = nrow(theObject@peptide_data),
                   output_rows = nrow(peptide_values_imputed),
                   imputed_rows = 0L
                 )
               )
             } else {
               filter_result <- peptideMissingValueImputationHelper(
                 input_table = theObject@peptide_data,
                 metadata_table = theObject@design_matrix,
                 quantity_to_impute_column = !!rlang::sym(theObject@raw_quantity_column),
                 imputed_value_column = !!rlang::sym(imputed_value_column),
                 core_utilisation = core_utilisation,
                 input_table_sample_id_column = !!rlang::sym(theObject@sample_id),
                 sample_id_tbl_sample_id_column = !!rlang::sym(theObject@sample_id),
                 replicate_group_column = !!rlang::sym(replicate_group_column),
                 protein_id_column = !!rlang::sym(theObject@protein_id_column),
                 peptide_sequence_column = !!rlang::sym(theObject@peptide_sequence_column),
                 hek_string = hek_string,
                 proportion_missing_values = proportion_missing_values,
                 exclusion_column = exclusion_column,
                 return_imputation_result = TRUE
               )
               if (is.data.frame(filter_result)) {
                 filter_result <- list(
                   data = filter_result,
                   support_table = NULL,
                   summary = list(status = "applied", skip_reason = NULL)
                 )
               }
             }

             section$imputed_value_column <- imputed_value_column
             section$proportion_missing_values <- as.numeric(proportion_missing_values)
             section$core_utilisation <- core_utilisation
             section$exclusion_column <- exclusion_column
             section$hek_string <- hek_string
             section$quantity_column <- theObject@raw_quantity_column
             section$support_table <- filter_result$support_table
             section$imputation_summary <- filter_result$summary
             theObject@args$peptideMissingValueImputation <- section
             theObject@peptide_data <- filter_result$data
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

              slope_interpretation <- if (dpcfit$dpc[2] < 0.3) {
                "nearly random missing"
              } else if (dpcfit$dpc[2] < 0.7) {
                "moderate intensity-dependent missing"
              } else if (dpcfit$dpc[2] < 1.2) {
                "strong intensity-dependent missing"
              } else {
                "very strong intensity-dependent missing (approaching left-censoring)"
              }
              
              if (verbose) {
                log_info("DPC parameters estimated:")
                log_info("  beta0 (intercept): {round(dpcfit$dpc[1], 4)}")
                log_info("  beta1 (slope): {round(dpcfit$dpc[2], 4)}")
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
              
              # Convert with the authoritative feature map; matrix row labels
              # are never parsed as biological identifiers.
              imputed_identity <- .peptideMatrixToIdentityLong(
                theObject,
                imputed_matrix,
                imputed_value_column
              )
              theObject <- imputed_identity$theObject
              imputed_long <- imputed_identity$data
              
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
