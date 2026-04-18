## ----------------------------------------------------------------------------------------------------------------------------------------------------------------------
#' Get Low Coefficient of Variation Proteins
#'
#' Identifies proteins with the lowest coefficient of variation for use as
#' negative controls in RUV normalization.
#'
#' @param theObject A ProteinQuantitativeData object
#' @param percentage_as_neg_ctrl Percentage of proteins to use as negative controls (default: 10)
#' @param num_neg_ctrl Number of negative control proteins (overrides percentage if set)
#' @return A logical vector indicating which proteins are selected as negative controls
#' @export
setMethod(
  f = "getLowCoefficientOfVariationProteins",
  signature = "ProteinQuantitativeData",
  definition = function(
    theObject,
    percentage_as_neg_ctrl = NULL,
    num_neg_ctrl = NULL
  ) {
    percentage_as_neg_ctrl <- checkParamsObjectFunctionSimplify(theObject, "percentage_as_neg_ctrl", 10)
    num_neg_ctrl <- checkParamsObjectFunctionSimplify(
      theObject,
      "num_neg_ctrl",
      round(nrow(theObject@protein_quant_table) * percentage_as_neg_ctrl / 100, 0)
    )

    theObject <- updateParamInObject(theObject, "percentage_as_neg_ctrl")
    theObject <- updateParamInObject(theObject, "num_neg_ctrl")

    list_of_control_genes <- theObject@protein_quant_table |>
      column_to_rownames(theObject@protein_id_column) |>
      t() |>
      as.data.frame() |>
      summarise(across(everything(), ~ sd(.) / mean(.))) |>
      t() |>
      as.data.frame() |>
      dplyr::rename(coefficient_of_variation = "V1") |>
      tibble::rownames_to_column(theObject@protein_id_column) |>
      arrange(coefficient_of_variation) |>
      head(num_neg_ctrl)

    control_gene_index_helper <- theObject@protein_quant_table |>
      dplyr::select(theObject@protein_id_column) |>
      mutate(index = row_number()) |>
      left_join(list_of_control_genes, by = theObject@protein_id_column) |>
      mutate(is_selected = case_when(
        is.na(coefficient_of_variation) ~ FALSE,
        TRUE ~ TRUE
      )) |>
      arrange(index) |>
      dplyr::select(!!sym(theObject@protein_id_column), is_selected) |>
      column_to_rownames(theObject@protein_id_column) |>
      t()

    control_gene_index <- control_gene_index_helper[1, ]

    control_gene_index
  }
)

