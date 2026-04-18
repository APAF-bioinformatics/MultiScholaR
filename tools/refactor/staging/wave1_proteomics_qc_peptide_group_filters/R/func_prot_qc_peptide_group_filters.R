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

