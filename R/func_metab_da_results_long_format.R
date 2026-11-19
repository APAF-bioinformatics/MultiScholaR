# ----------------------------------------------------------------------------
# createMetabDaResultsLongFormat
# ----------------------------------------------------------------------------
#' Create metabolomics DA results in long format with sample intensity columns
#'
#' @description Creates a long-format DA results table that includes individual
#'   sample intensity values, mirroring the proteomics `createDaResultsLongFormat()`.
#'   The output includes columns for each sample's intensity value, named as
#'   `intensity.{sample_id}.{group}`.
#'
#' @param lfc_qval_tbl Data frame with DA statistics (logFC, pvalue, fdr_qvalue).
#'   Must have a metabolite_id column and a comparison column containing the
#'   contrast string in format "groupA-groupB".
#' @param expr_matrix Expression matrix with metabolites as rows, samples as columns.
#'   Rownames must match metabolite_id column in lfc_qval_tbl.
#' @param design_matrix Design matrix with sample and group assignments.
#' @param sample_id_col Name of sample ID column in design_matrix.
#' @param group_id_col Name of group ID column in design_matrix.
#' @param metabolite_id_col Name of metabolite ID column.
#'
#' @return Data frame with DA results plus intensity columns for each sample.
#'
#' @importFrom dplyr left_join mutate select rename arrange distinct filter pull
#' @importFrom tidyr pivot_longer pivot_wider separate_wider_delim
#' @importFrom tibble rownames_to_column
#' @importFrom rlang sym set_names
#' @importFrom stringr str_replace_all
#' @export
createMetabDaResultsLongFormat <- function(
  lfc_qval_tbl,
  expr_matrix,
  design_matrix,
  sample_id_col = "sample_id",
  group_id_col = "group",
  metabolite_id_col = "metabolite_id"
) {

    logger::log_info(sprintf("   lfc_qval_tbl: %d rows", nrow(lfc_qval_tbl)))
    logger::log_info(sprintf(
        "   expr_matrix: %d metabolites x %d samples",
        nrow(expr_matrix), ncol(expr_matrix)
    ))

    # Convert expression matrix to long format with intensity values
    intensity_long <- tryCatch(
        {
            expr_matrix |>
                as.data.frame() |>
                tibble::rownames_to_column(metabolite_id_col) |>
                tidyr::pivot_longer(
                    cols = -!!rlang::sym(metabolite_id_col),
                    names_to = sample_id_col,
                    values_to = "intensity"
                ) |>
                dplyr::left_join(design_matrix, by = sample_id_col)
        },
        error = function(e) {
            cat(sprintf("[ERROR] pivot_longer/left_join failed: %s\n", e$message))
            stop(e)
        }
    )


    # Create intensity columns per group
    # Format: intensity.{sample_id}.{group}
    intensity_wide <- tryCatch(
        {
            intensity_long |>
                dplyr::mutate(
                    col_name = paste0("intensity.", !!rlang::sym(sample_id_col), ".", !!rlang::sym(group_id_col))
                ) |>
                dplyr::select(!!rlang::sym(metabolite_id_col), col_name, intensity) |>
                tidyr::pivot_wider(
                    id_cols = !!rlang::sym(metabolite_id_col),
                    names_from = col_name,
                    values_from = intensity
                )
        },
        error = function(e) {
            cat(sprintf("[ERROR] pivot_wider failed: %s\n", e$message))
            stop(e)
        }
    )


    # Parse contrast to get numerator/denominator groups
    if ("comparison" %in% colnames(lfc_qval_tbl)) {
        # Try to parse numerator/denominator
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
                lfc_qval_tbl$numerator <- NA_character_
                lfc_qval_tbl$denominator <- NA_character_
            }
        )
    }

    # Join DE results with intensity values
    da_results_long <- tryCatch(
        {
            lfc_qval_tbl |>
                dplyr::left_join(intensity_wide, by = metabolite_id_col) |>
                dplyr::arrange(comparison, fdr_qvalue, logFC) |>
                dplyr::distinct()
        },
        error = function(e) {
            cat(sprintf("[ERROR] final join failed: %s\n", e$message))
            stop(e)
        }
    )

    return(da_results_long)
}

