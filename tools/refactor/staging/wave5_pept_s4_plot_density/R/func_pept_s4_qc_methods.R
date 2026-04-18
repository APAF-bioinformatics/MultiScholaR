# Density plot for PeptideQuantitativeData

#' @export
setMethod(
  f = "plotDensity",
  signature = "PeptideQuantitativeData",
  definition = function(theObject, grouping_variable, title = "", font_size = 8) {
    peptide_data <- theObject@peptide_data
    quant_col <- theObject@norm_quantity_column
    sample_id <- theObject@sample_id
    
    peptide_data |>
      filter(!is.na(!!sym(quant_col))) |>
      ggplot(aes(x = !!sym(quant_col), color = !!sym(sample_id))) +
      geom_density() +
      apafTheme() +
      labs(title = title, x = "Log2 Intensity", y = "Density") +
      theme(legend.position = "none")
  }
)

## ----------------------------------------------------------------------------------------------------------------------------------------------------------------------

