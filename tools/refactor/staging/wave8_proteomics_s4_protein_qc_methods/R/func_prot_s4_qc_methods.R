#' @export
setMethod(
  f = "removeProteinsWithOnlyOneReplicate",
  definition = function(theObject, core_utilisation = NULL, grouping_variable = NULL) {
    protein_quant_table <- theObject@protein_quant_table
    samples_id_tbl <- theObject@design_matrix
    sample_id_tbl_sample_id_column <- theObject@sample_id
    # replicate_group_column <- theObject@technical_replicate_id
    protein_id_column <- theObject@protein_id_column

    input_table_sample_id_column <- theObject@sample_id
    quantity_column <- "log_values"

    grouping_variable <- checkParamsObjectFunctionSimplifyAcceptNull(
      theObject,
      "grouping_variable",
      NULL
    )

    core_utilisation <- checkParamsObjectFunctionSimplify(
      theObject,
      "core_utilisation",
      NA
    )

    theObject <- updateParamInObject(theObject, "grouping_variable")
    theObject <- updateParamInObject(theObject, "core_utilisation")

    data_long_cln <- protein_quant_table |>
      pivot_longer(
        cols = !matches(protein_id_column),
        names_to = input_table_sample_id_column,
        values_to = quantity_column
      )

    protein_quant_table <- removeProteinsWithOnlyOneReplicateHelper(
      input_table = data_long_cln,
      samples_id_tbl = samples_id_tbl,
      input_table_sample_id_column = !!sym(input_table_sample_id_column),
      sample_id_tbl_sample_id_column = !!sym(sample_id_tbl_sample_id_column),
      replicate_group_column = !!sym(grouping_variable),
      protein_id_column = !!sym(protein_id_column),
      quantity_column = !!sym(quantity_column),
      core_utilisation = core_utilisation
    )


    theObject@protein_quant_table <- protein_quant_table |>
      pivot_wider(
        id_cols = !!sym(protein_id_column),
        names_from = !!sym(input_table_sample_id_column),
        values_from = !!sym(quantity_column)
      )

    theObject <- cleanDesignMatrix(theObject)

    updated_object <- theObject

    return(updated_object)
  }
)

