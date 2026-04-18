#' @export
setMethod(
  f = "getPcaMatrix",
  signature = "ProteinQuantitativeData",
  definition = function(theObject) {
    protein_quant_table <- theObject@protein_quant_table
    protein_id_column <- theObject@protein_id_column
    design_matrix <- theObject@design_matrix
    sample_id <- theObject@sample_id


    frozen_protein_matrix <- protein_quant_table |>
      column_to_rownames(protein_id_column) |>
      as.matrix()

    frozen_protein_matrix_pca <- frozen_protein_matrix
    frozen_protein_matrix_pca[!is.finite(frozen_protein_matrix_pca)] <- NA


    pca_mixomics_before_cyclic_loess <- mixOmics::pca(t(as.matrix(frozen_protein_matrix_pca)))$variates$X |>
      as.data.frame() |>
      rownames_to_column(sample_id) |>
      left_join(design_matrix, by = sample_id)


    return(pca_mixomics_before_cyclic_loess)
  }
)

