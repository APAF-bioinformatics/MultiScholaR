## ----------------------------------------------------------------------------------------------------------------------------------------------------------------------
#' @title Normalise between Arrays
#' @export
#' @param theObject Object of class ProteinQuantitativeData
#' @param normalisation_method Method to use for normalisation. Options are cyclicloess, quantile, scale, none
setMethod(
  f = "normaliseBetweenSamples",
  signature = "ProteinQuantitativeData",
  definition = function(theObject, normalisation_method = NULL) {
    protein_quant_table <- theObject@protein_quant_table
    protein_id_column <- theObject@protein_id_column
    design_matrix <- theObject@design_matrix
    sample_id <- theObject@sample_id

    normalisation_method <- checkParamsObjectFunctionSimplify(
      theObject,
      "normalisation_method",
      "cyclicloess"
    )

    theObject <- updateParamInObject(theObject, "normalisation_method")

    frozen_protein_matrix <- protein_quant_table |>
      column_to_rownames(protein_id_column) |>
      as.matrix()

    frozen_protein_matrix[!is.finite(frozen_protein_matrix)] <- NA

    normalised_frozen_protein_matrix <- frozen_protein_matrix

    print(paste0("normalisation_method = ", normalisation_method))

    switch(normalisation_method,
      cyclicloess = {
        normalised_frozen_protein_matrix <- normalizeCyclicLoess(frozen_protein_matrix)
      },
      quantile = {
        normalised_frozen_protein_matrix <- normalizeQuantiles(frozen_protein_matrix)
      },
      scale = {
        normalised_frozen_protein_matrix <- normalizeMedianAbsValues(frozen_protein_matrix)
      },
      none = {
        normalised_frozen_protein_matrix <- frozen_protein_matrix
      }
    )

    normalised_frozen_protein_matrix[!is.finite(normalised_frozen_protein_matrix)] <- NA

    # normalised_frozen_protein_matrix_filt <- as.data.frame( normalised_frozen_protein_matrix ) |>
    #   dplyr::filter( if_all( everything(), \(x) { !is.na(x) } ) ) |>
    #   as.matrix()

    theObject@protein_quant_table <- normalised_frozen_protein_matrix |>
      as.data.frame() |>
      rownames_to_column(protein_id_column)

    theObject <- cleanDesignMatrix(theObject)

    updated_object <- theObject

    return(updated_object)
  }
)

