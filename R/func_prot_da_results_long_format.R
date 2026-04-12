# createDaResultsLongFormat
# ----------------------------------------------------------------------------
## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Create the de_protein_long and de_phos_long tables
#' @export
createDaResultsLongFormat <- function(lfc_qval_tbl,
                                      norm_counts_input_tbl,
                                      raw_counts_input_tbl,
                                      row_id,
                                      sample_id,
                                      group_id,
                                      group_pattern,
                                      design_matrix_norm,
                                      design_matrix_raw,
                                      expression_column = log_intensity,
                                      protein_id_table) {
  message("   DEBUG66: createDaResultsLongFormat - Starting norm_counts processing")
  message(sprintf("      DEBUG66: norm_counts_input_tbl dims = %d x %d", nrow(norm_counts_input_tbl), ncol(norm_counts_input_tbl)))
  message(sprintf("      DEBUG66: group_pattern = %s", group_pattern))
  message(sprintf("      DEBUG66: row_id = %s", row_id))
  message(sprintf("      DEBUG66: sample_id = %s", sample_id))
  message(sprintf("      DEBUG66: group_id = %s", group_id))

  norm_counts <- norm_counts_input_tbl |>
    as.data.frame() |>
    rownames_to_column(row_id) |>
    pivot_longer(
      cols = matches(group_pattern),
      names_to = sample_id,
      values_to = "log2norm"
    ) |>
    left_join(design_matrix_norm, by = sample_id) |>
    group_by(!!sym(row_id), !!sym(group_id)) |>
    arrange(!!sym(row_id), !!sym(group_id), !!sym(sample_id)) |>
    mutate(replicate_number = paste0("log2norm.", row_number())) |>
    ungroup() |>
    pivot_wider(
      id_cols = c(!!sym(row_id), !!sym(group_id)),
      names_from = replicate_number,
      values_from = log2norm
    ) |>
    mutate({{ group_id }} := purrr::map_chr(!!sym(group_id), as.character))

  message("   DEBUG66: norm_counts processing completed")


  # print(head(norm_counts))

  raw_counts <- raw_counts_input_tbl |>
    as.data.frame() |>
    rownames_to_column(row_id) |>
    pivot_longer(
      cols = matches(group_pattern),
      names_to = sample_id,
      values_to = "raw"
    ) |>
    left_join(design_matrix_raw, by = sample_id) |>
    group_by(!!sym(row_id), !!sym(group_id)) |>
    arrange(!!sym(row_id), !!sym(group_id), !!sym(sample_id)) |>
    mutate(replicate_number = paste0("raw.", row_number())) |>
    ungroup() |>
    pivot_wider(
      id_cols = c(!!sym(row_id), !!sym(group_id)),
      names_from = replicate_number,
      values_from = raw
    ) |>
    mutate({{ group_id }} := purrr::map_chr(!!sym(group_id), as.character))

  # print(head(raw_counts))

  left_join_columns <- rlang::set_names(
    c(row_id, group_id),
    c(row_id, "left_group")
  )

  right_join_columns <- rlang::set_names(
    c(row_id, group_id),
    c(row_id, "right_group")
  )

  # print(head(lfc_qval_tbl))

  # DEBUG66: Commented out print statements that were causing confusion
  # print( row_id)
  # print(colnames( protein_id_table)[1])

  da_proteins_long <- lfc_qval_tbl |>
    dplyr::select(-lqm, -colour, -analysis_type) |>
    dplyr::mutate({{ expression_column }} := str_replace_all({{ expression_column }}, group_id, "")) |>
    separate_wider_delim({{ expression_column }}, delim = "-", names = c("left_group", "right_group")) |>
    left_join(norm_counts, by = left_join_columns) |>
    left_join(norm_counts,
      by = right_join_columns,
      suffix = c(".left", ".right")
    ) |>
    left_join(raw_counts, by = left_join_columns) |>
    left_join(raw_counts,
      by = right_join_columns,
      suffix = c(".left", ".right")
    ) |>
    left_join(protein_id_table,
      by = join_by(!!sym(row_id) == !!sym(colnames(protein_id_table)[1]))
    ) |>
    arrange(comparison, fdr_qvalue, log2FC) |>
    distinct()

  # --- NEW: Rename columns to use sample IDs if single contrast ---
  # Only perform this renaming if we have a single comparison, to ensure unique mapping
  if (length(unique(da_proteins_long$comparison)) == 1) {
    # Get the groups involved
    this_left_group <- unique(da_proteins_long$left_group)
    this_right_group <- unique(da_proteins_long$right_group)

    # Ensure we have exactly one left and one right group
    if (length(this_left_group) == 1 && length(this_right_group) == 1) {
      message(sprintf("   createDaResultsLongFormat: Renaming columns for contrast %s vs %s", this_left_group, this_right_group))

      # Helper to get sorted sample IDs for a group
      get_samples_for_group <- function(dm, grp) {
        dm |>
          dplyr::filter(!!sym(group_id) == grp) |>
          dplyr::arrange(!!sym(sample_id)) |>
          dplyr::pull(!!sym(sample_id))
      }

      left_samples <- get_samples_for_group(design_matrix_norm, this_left_group)
      right_samples <- get_samples_for_group(design_matrix_norm, this_right_group)

      # Helper to generate rename mapping using vectorized operations
      # Returns: named vector c(new_name = old_name) for dplyr::rename
      generate_rename_map <- function(df, prefix, suffix, samples, group_name) {
        indices <- seq_along(samples)
        old_cols <- paste0(prefix, ".", indices, suffix)
        new_cols <- paste0(prefix, ".", samples, ".", group_name)

        # Only include columns that exist in the dataframe
        valid_idx <- old_cols %in% colnames(df)

        if (any(valid_idx)) {
          return(setNames(old_cols[valid_idx], new_cols[valid_idx]))
        } else {
          return(character(0))
        }
      }

      # Generate mappings for all 4 sets of columns
      map1 <- generate_rename_map(da_proteins_long, "log2norm", ".left", left_samples, this_left_group)
      map2 <- generate_rename_map(da_proteins_long, "raw", ".left", left_samples, this_left_group)
      map3 <- generate_rename_map(da_proteins_long, "log2norm", ".right", right_samples, this_right_group)
      map4 <- generate_rename_map(da_proteins_long, "raw", ".right", right_samples, this_right_group)

      # Combine all mappings
      all_mappings <- c(map1, map2, map3, map4)

      # Apply renaming in a single vectorized step
      if (length(all_mappings) > 0) {
        da_proteins_long <- da_proteins_long |> dplyr::rename(!!!all_mappings)
        message(sprintf("   createDaResultsLongFormat: Renamed %d columns", length(all_mappings)))
      }
    }
  }
  # --- END NEW ---

  # Rename group columns to numerator/denominator for clarity
  da_proteins_long <- da_proteins_long |>
    dplyr::rename(numerator = left_group, denominator = right_group)

  da_proteins_long
}
# ----------------------------------------------------------------------------
# getTypeOfGrouping
# ----------------------------------------------------------------------------
#' Assign experimental group list
#' @param design_matrix A data frame representing the design matrix.
#' @param group_id A string representing the name of the group ID column used in the design matrix.
#' @param sample_id A string representing the name of the sample ID column used in the design matrix.
#' @return A list where each element name is the name of a treatment group and each element is a vector containing the sample IDs within the treatment group.
#' @export
getTypeOfGrouping <- function(design_matrix, group_id, sample_id) {
  temp_type_of_grouping <- design_matrix |>
    dplyr::select(!!rlang::sym(group_id), !!rlang::sym(sample_id)) |>
    group_by(!!rlang::sym(group_id)) |>
    summarise(!!rlang::sym(sample_id) := list(!!rlang::sym(sample_id))) |>
    ungroup()

  type_of_grouping <- temp_type_of_grouping |> dplyr::pull(!!rlang::sym(sample_id))
  names(type_of_grouping) <- temp_type_of_grouping |> dplyr::pull(!!rlang::sym(group_id))

  return(type_of_grouping)
}

# ----------------------------------------------------------------------------
# extractResults
# ----------------------------------------------------------------------------
#' @export
extractResults <- function(results_list) {
  extracted <- purrr::map(results_list, \(x){
    x$results
  })

  names(extracted) <- names(results_list)

  return(extracted)
}
