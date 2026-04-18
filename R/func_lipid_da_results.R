# ----------------------------------------------------------------------------
# createLipidDaResultsLongFormat
# ----------------------------------------------------------------------------
#' Create lipidomics DA results in long format with sample intensity columns
#'
#' @description Creates a long-format DA results table that includes individual
#'   sample intensity values, mirroring the proteomics `createDaResultsLongFormat()`.
#'   The output includes columns for each sample's intensity value, named as
#'   `intensity.{sample_id}.{group}`.
#'
#' @param lfc_qval_tbl Data frame with DE statistics (logFC, pvalue, fdr_qvalue).
#'   Must have a lipid_id column and a comparison column containing the
#'   contrast string in format "groupA-groupB".
#' @param expr_matrix Expression matrix with lipids as rows, samples as columns.
#'   Rownames must match lipid_id column in lfc_qval_tbl.
#' @param design_matrix Design matrix with sample and group assignments.
#' @param sample_id_col Name of sample ID column in design_matrix.
#' @param group_id_col Name of group ID column in design_matrix.
#' @param lipid_id_col Name of lipid ID column.
#'
#' @return Data frame with DE results plus intensity columns for each sample.
#'
#' @importFrom dplyr left_join mutate select rename arrange distinct filter pull
#' @importFrom tidyr pivot_longer pivot_wider separate_wider_delim
#' @importFrom tibble rownames_to_column
#' @importFrom rlang sym set_names
#' @importFrom stringr str_replace_all
#' @export
createLipidDaResultsLongFormat <- function(
  lfc_qval_tbl,
  expr_matrix,
  design_matrix,
  sample_id_col = "sample_id",
  group_id_col = "group",
  lipid_id_col = "lipid_id"
) {
    logger::log_info("--- Entering createLipidDaResultsLongFormat ---")
    logger::log_info(sprintf("   lfc_qval_tbl: %d rows", nrow(lfc_qval_tbl)))
    logger::log_info(sprintf(
        "   expr_matrix: %d lipids x %d samples",
        nrow(expr_matrix), ncol(expr_matrix)
    ))

    # Convert expression matrix to long format with intensity values
    intensity_long <- expr_matrix |>
        as.data.frame() |>
        tibble::rownames_to_column(lipid_id_col) |>
        tidyr::pivot_longer(
            cols = -!!rlang::sym(lipid_id_col),
            names_to = sample_id_col,
            values_to = "intensity"
        ) |>
        dplyr::left_join(design_matrix, by = sample_id_col)

    logger::log_info(sprintf("   intensity_long: %d rows", nrow(intensity_long)))

    # Create intensity columns per group
    # Format: intensity.{sample_id}.{group}
    intensity_wide <- intensity_long |>
        dplyr::mutate(
            col_name = paste0("intensity.", !!rlang::sym(sample_id_col), ".", !!rlang::sym(group_id_col))
        ) |>
        dplyr::select(!!rlang::sym(lipid_id_col), col_name, intensity) |>
        tidyr::pivot_wider(
            id_cols = !!rlang::sym(lipid_id_col),
            names_from = col_name,
            values_from = intensity
        )

    logger::log_info(sprintf(
        "   intensity_wide: %d lipids x %d columns",
        nrow(intensity_wide), ncol(intensity_wide)
    ))

    # Parse contrast to get numerator/denominator groups
    # Contrast format is "groupA-groupB" where groupA is numerator, groupB is denominator
    if ("comparison" %in% colnames(lfc_qval_tbl)) {
        # Extract unique comparisons
        contrasts <- unique(lfc_qval_tbl$comparison)
        logger::log_info(sprintf(
            "   Processing %d contrasts: %s",
            length(contrasts), paste(contrasts, collapse = ", ")
        ))

        # Parse first contrast to get group prefix pattern
        # Assumes format like "groupTreatment-groupControl"
        first_contrast <- contrasts[1]
        parts <- strsplit(first_contrast, "-")[[1]]

        if (length(parts) == 2) {
            # Try to detect group prefix (e.g., "group" from "groupTreatment")
            left_part <- parts[1]
            right_part <- parts[2]

            # Add numerator/denominator columns by parsing comparison
            tryCatch(
                {
                    lfc_qval_tbl <- lfc_qval_tbl |>
                        tidyr::separate_wider_delim(
                            comparison,
                            delim = "-",
                            names = c("numerator", "denominator"),
                            cols_remove = FALSE
                        )
                },
                error = function(e) {
                    cat(sprintf("[WARN] Could not parse contrast groups: %s\n", e$message))
                    lfc_qval_tbl$numerator <<- NA_character_
                    lfc_qval_tbl$denominator <<- NA_character_
                }
            )
        }
    }

    # Join DE results with intensity values
    da_results_long <- lfc_qval_tbl |>
        dplyr::left_join(intensity_wide, by = lipid_id_col) |>
        dplyr::arrange(comparison, fdr_qvalue, logFC) |>
        dplyr::distinct()

    logger::log_info(sprintf(
        "   Final da_results_long: %d rows x %d columns",
        nrow(da_results_long), ncol(da_results_long)
    ))

    logger::log_info("--- Exiting createLipidDaResultsLongFormat ---")

    return(da_results_long)
}

