#' @export
setMethod(
  f = "getRuvIIIReplicateMatrix",
  signature = "ProteinQuantitativeData",
  definition = function(theObject, ruv_grouping_variable = NULL) {
    protein_quant_table <- theObject@protein_quant_table
    protein_id_column <- theObject@protein_id_column
    design_matrix <- theObject@design_matrix
    group_id <- theObject@group_id
    sample_id <- theObject@sample_id
    replicate_group_column <- theObject@technical_replicate_id

    ruv_grouping_variable <- checkParamsObjectFunctionSimplify(theObject, "ruv_grouping_variable", NULL)

    theObject <- updateParamInObject(theObject, "ruv_grouping_variable")

    ruvIII_replicates_matrix <- getRuvIIIReplicateMatrixHelper(
      design_matrix,
      !!sym(sample_id),
      !!sym(ruv_grouping_variable)
    )
    return(ruvIII_replicates_matrix)
  }
)

