# ----------------------------------------------------------------------------
# getMetaboliteQuantData
# ----------------------------------------------------------------------------
#' Extract quantitative data columns from an assay data frame
#'
#' @description Helper function to separate metabolite quantitative data
#'   (sample columns) from metadata columns (ID, annotation, etc.).
#'
#' @param assay_df Data frame containing metabolite data for one assay.
#' @param metabolite_id_col Name of the metabolite ID column.
#' @param annotation_col Name of the annotation column.
#' @param additional_meta_cols Additional columns to exclude from quant data.
#'
#' @return A list with:
#'   - quant_data: Data frame with only sample columns
#'   - meta_data: Data frame with only metadata columns
#'   - sample_cols: Names of sample columns
#'
#' @export
getMetaboliteQuantData <- function(
  assay_df,
  metabolite_id_col = "Alignment ID",
  annotation_col = "Metabolite name",
  additional_meta_cols = NULL
) {
    # Identify metadata columns to exclude
    meta_cols <- c(metabolite_id_col, annotation_col)
    if (!is.null(additional_meta_cols)) {
        meta_cols <- c(meta_cols, additional_meta_cols)
    }

    # Get sample columns (everything that's not metadata)
    all_cols <- colnames(assay_df)
    sample_cols <- setdiff(all_cols, meta_cols)

    # Also exclude any obviously non-numeric columns
    sample_cols <- sample_cols[sapply(assay_df[, sample_cols, drop = FALSE], is.numeric)]

    list(
        quant_data = assay_df[, sample_cols, drop = FALSE],
        meta_data = assay_df[, intersect(meta_cols, all_cols), drop = FALSE],
        sample_cols = sample_cols
    )
}

