#' @export
setMethod(
  f = "plotRleList",
  signature = "ProteinQuantitativeData",
  definition = function(theObject, list_of_columns, yaxis_limit = c()) {
    protein_quant_table <- theObject@protein_quant_table
    protein_id_column <- theObject@protein_id_column
    design_matrix <- theObject@design_matrix
    sample_id <- theObject@sample_id

    frozen_protein_matrix <- protein_quant_table |>
      column_to_rownames(protein_id_column) |>
      as.matrix()

    design_matrix <- as.data.frame(design_matrix)
    rownames(design_matrix) <- design_matrix[, sample_id]

    # print( design_matrix)

    runOneRle <- function(column_name) {
      rowinfo_vector <- NA

      if (column_name %in% colnames(design_matrix)) {
        rowinfo_vector <- design_matrix[colnames(frozen_protein_matrix), column_name]
      }

      rle_plot_before_cyclic_loess <- plotRleHelper(t(frozen_protein_matrix),
        rowinfo = rowinfo_vector,
        yaxis_limit = yaxis_limit
      )

      return(rle_plot_before_cyclic_loess)
    }

    list_of_rle_plots <- purrr::map(list_of_columns, runOneRle)

    names(list_of_rle_plots) <- list_of_columns

    return(list_of_rle_plots)
  }
)

