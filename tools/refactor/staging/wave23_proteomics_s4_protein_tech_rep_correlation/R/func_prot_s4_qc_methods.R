#' @export
setMethod(
  f = "proteinTechRepCorrelation",
  signature = "ProteinQuantitativeData",
  definition = function(theObject, tech_rep_num_column = NULL, tech_rep_remove_regex = NULL) {
    protein_quant_table <- theObject@protein_quant_table
    protein_id_column <- theObject@protein_id_column
    design_matrix <- theObject@design_matrix
    sample_id <- theObject@sample_id
    tech_rep_column <- theObject@technical_replicate_id

    tech_rep_num_column <- checkParamsObjectFunctionSimplifyAcceptNull(theObject, "tech_rep_num_column", NULL)
    tech_rep_remove_regex <- checkParamsObjectFunctionSimplifyAcceptNull(theObject, "tech_rep_remove_regex", NULL)

    theObject <- updateParamInObject(theObject, "tech_rep_num_column")
    theObject <- updateParamInObject(theObject, "tech_rep_remove_regex")

    frozen_protein_matrix <- protein_quant_table |>
      column_to_rownames(protein_id_column) |>
      as.matrix()

    frozen_protein_matrix_pca <- frozen_protein_matrix
    frozen_protein_matrix_pca[!is.finite(frozen_protein_matrix_pca)] <- NA

    protein_matrix_tech_rep <- proteinTechRepCorrelationHelper(design_matrix, frozen_protein_matrix_pca,
      protein_id_column = protein_id_column,
      sample_id_column = sample_id,
      tech_rep_column = tech_rep_column,
      tech_rep_num_column = tech_rep_num_column,
      tech_rep_remove_regex = tech_rep_remove_regex
    )

    return(protein_matrix_tech_rep)
  }
)

