
# ----------------------------------------------------------------------------
# getLipidQuantData
# ----------------------------------------------------------------------------
#' @title Extract Quantitative Data and Sample Names from Assay Tibble
#' @description Separates annotation columns from quantitative data columns
#'              in a lipidomics assay tibble.
#'
#' @param assay_data A tibble/data.frame representing one lipidomics assay,
#'                   with lipid annotations and sample intensity columns.
#' @param sample_columns Optional character vector of explicit sample column names.
#'                       When provided, these columns are used directly instead of
#'                       guessing based on numeric type. This is the preferred method
#'                       as MS-DIAL data contains many numeric annotation columns
#'                       (scores, m/z, RT) that are NOT sample data.
#'
#' @return A list containing:
#'         - `quant_data`: A data frame with only the quantitative (sample) columns.
#'         - `sample_names`: A character vector of the sample column names.
#'         - `annotation_data`: A data frame with the non-sample annotation columns.
#'
#' @importFrom dplyr select where
#' @keywords internal
#' @noRd
#' @export
getLipidQuantData <- function(assay_data, sample_columns = NULL) {
    if (!is.null(sample_columns) && length(sample_columns) > 0) {
        # Use explicit sample columns (preferred - avoids including numeric annotation cols)
        valid_cols <- intersect(sample_columns, colnames(assay_data))
        if (length(valid_cols) == 0) {
            warning("None of the provided sample_columns exist in assay_data. Falling back to numeric detection.")
            quant_cols <- sapply(assay_data, is.numeric)
            quant_data <- assay_data[, quant_cols, drop = FALSE]
            sample_names <- colnames(quant_data)
        } else {
            quant_data <- assay_data[, valid_cols, drop = FALSE]
            sample_names <- valid_cols
        }
    } else {
        # Fallback: guess by numeric type (unreliable for MS-DIAL data with numeric scores)
        quant_cols <- sapply(assay_data, is.numeric)
        quant_data <- assay_data[, quant_cols, drop = FALSE]
        sample_names <- colnames(quant_data)
    }

    # Identify annotation columns (everything NOT in sample_names)
    annotation_data <- assay_data[, setdiff(colnames(assay_data), sample_names), drop = FALSE]

    return(list(
        quant_data = as.data.frame(quant_data),
        sample_names = sample_names,
        annotation_data = as.data.frame(annotation_data)
    ))
}

