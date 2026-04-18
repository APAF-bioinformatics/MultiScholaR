# Helper function to get counts table
getCountsTable <- function(obj) {
    if (inherits(obj, "MetaboliteAssayData")) {
        message(sprintf("   Getting counts table for object of class: %s", class(obj)[1]))
        message(sprintf(
            "   Returning metabolite_data with dimensions: %d rows, %d cols",
            nrow(obj@metabolite_data), ncol(obj@metabolite_data)
        ))
        obj@metabolite_data
    } else if (inherits(obj, "ProteinQuantitativeData")) {
        message(sprintf(
            "   Returning protein_quant_table with dimensions: %d rows, %d cols",
            nrow(obj@protein_quant_table), ncol(obj@protein_quant_table)
        ))
        obj@protein_quant_table
    } else {
        message(sprintf("   ERROR: Unsupported object type: %s", class(obj)[1]))
        stop("Unsupported object type")
    }
}

