#' @export
setMethod(
  f = "setProteinData",
  signature = "ProteinQuantitativeData",
  definition = function(theObject, protein_quant_table, protein_id_column) {
    theObject@protein_quant_table <- protein_quant_table
    theObject@protein_id_column <- protein_id_column

    return(theObject)
  }
)

