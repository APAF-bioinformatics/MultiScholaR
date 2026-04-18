#' @export
setMethod(
  f = "plotRle",
  signature = "ProteinQuantitativeData",
  definition = function(theObject, grouping_variable, yaxis_limit = c(), sample_label = NULL) {
    protein_quant_table <- theObject@protein_quant_table
    protein_id_column <- theObject@protein_id_column
    design_matrix <- theObject@design_matrix
    sample_id <- theObject@sample_id

    frozen_protein_matrix <- protein_quant_table |>
      column_to_rownames(protein_id_column) |>
      as.matrix()

    design_matrix <- as.data.frame(design_matrix)

    if (!is.null(sample_label)) {
      if (sample_label %in% colnames(design_matrix)) {
        rownames(design_matrix) <- design_matrix[, sample_label]
        colnames(frozen_protein_matrix) <- design_matrix[, sample_label]
      }
    } else {
      rownames(design_matrix) <- design_matrix[, sample_id]
    }

    # print( design_matrix)

    rowinfo_vector <- NA
    if (!is.na(grouping_variable)) {
      rowinfo_vector <- design_matrix[colnames(frozen_protein_matrix), grouping_variable]
    }

    print(rownames(design_matrix))
    print(colnames(frozen_protein_matrix))
    print(rowinfo_vector)
    # Handle missing/non-finite values
    working_matrix <- frozen_protein_matrix
    working_matrix[!is.finite(working_matrix)] <- NA

    rle_plot_before_cyclic_loess <- plotRleHelper(t(working_matrix),
      rowinfo = rowinfo_vector,
      yaxis_limit = yaxis_limit
    )

    return(rle_plot_before_cyclic_loess)
  }
)

