## ----------------------------------------------------------------------------------------------------------------------------------------------------------------------
#' Average Technical Replicates
#'
#' Averages protein intensities across technical replicates, collapsing the
#' design matrix accordingly.
#'
#' @param theObject A ProteinQuantitativeData object
#' @param design_matrix_columns Additional columns to keep in the design matrix
#' @return The modified ProteinQuantitativeData object with averaged technical replicates
#' @rdname averageTechReps
#' @export
setMethod(
  f = "averageTechReps",
  signature = "ProteinQuantitativeData",
  definition = function(theObject, design_matrix_columns = c()) {
    protein_quant_table <- theObject@protein_quant_table
    protein_id_column <- theObject@protein_id_column
    design_matrix <- theObject@design_matrix
    group_id <- theObject@group_id
    sample_id <- theObject@sample_id
    replicate_group_column <- theObject@technical_replicate_id

    theObject@protein_quant_table <- protein_quant_table |>
      pivot_longer(
        cols = !matches(protein_id_column),
        names_to = sample_id,
        values_to = "Log2.Protein.Imputed"
      ) |>
      left_join(design_matrix,
        by = join_by(!!sym(sample_id) == !!sym(sample_id))
      ) |>
      group_by(!!sym(protein_id_column), !!sym(replicate_group_column)) |>
      summarise(Log2.Protein.Imputed = mean(Log2.Protein.Imputed, na.rm = TRUE)) |>
      ungroup() |>
      pivot_wider(
        names_from = !!sym(replicate_group_column),
        values_from = Log2.Protein.Imputed
      )

    theObject@sample_id <- theObject@technical_replicate_id

    theObject@design_matrix <- design_matrix |>
      dplyr::select(-!!sym(sample_id)) |>
      dplyr::select(all_of(unique(c(replicate_group_column, group_id, design_matrix_columns)))) |>
      distinct()

    theObject@sample_id <- replicate_group_column
    theObject@technical_replicate_id <- NA_character_

    theObject <- cleanDesignMatrix(theObject)

    theObject
  }
)

