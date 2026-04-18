#' @title Choose Best Protein Accession Sum Duplicates
#' @name chooseBestProteinAccessionSumDuplicates,ProteinQuantitativeData-method
#' @export
setMethod(
  f = "chooseBestProteinAccessionSumDuplicates",
  signature = "ProteinQuantitativeData",
  definition = function(theObject, delim = ";", quant_columns_pattern = "\\d+", islogged = TRUE) {
    protein_quant_table <- theObject@protein_quant_table
    protein_id_column <- theObject@protein_id_column

    protein_log2_quant_cln <- protein_quant_table |>
      mutate(!!sym(protein_id_column) := str_split_i(!!sym(protein_id_column), delim, 1)) |>
      group_by(!!sym(protein_id_column)) |>
      summarise(across(
        matches(quant_columns_pattern),
        \(x){
          if (islogged == TRUE) {
            log2(sum(2^x, na.rm = TRUE))
          } else {
            sum(x, na.rm = TRUE)
          }
        }
      )) |>
      ungroup()

    theObject@protein_quant_table <- protein_log2_quant_cln

    theObject
  }
)

