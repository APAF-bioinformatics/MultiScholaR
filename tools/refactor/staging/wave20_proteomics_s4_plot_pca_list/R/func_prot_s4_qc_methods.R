#' @export
setMethod(
  f = "plotPcaList",
  signature = "ProteinQuantitativeData",
  definition = function(theObject, grouping_variables_list, label_column, title, font_size = 8, cv_percentile = 0.90) {
    protein_quant_table <- theObject@protein_quant_table
    protein_id_column <- theObject@protein_id_column
    design_matrix <- theObject@design_matrix
    sample_id <- theObject@sample_id

    frozen_protein_matrix <- protein_quant_table |>
      column_to_rownames(protein_id_column) |>
      as.matrix()

    frozen_protein_matrix_pca <- frozen_protein_matrix
    frozen_protein_matrix_pca[!is.finite(frozen_protein_matrix_pca)] <- NA

    if (is.na(label_column) || label_column == "") {
      label_column <- ""
    }

    pca_plots_list <- plotPcaListHelper(frozen_protein_matrix_pca,
      design_matrix,
      sample_id_column = sample_id,
      grouping_variables_list = grouping_variables_list,
      label_column = label_column,
      title = title,
      geom.text.size = font_size
    )

    return(pca_plots_list)
  }
)

