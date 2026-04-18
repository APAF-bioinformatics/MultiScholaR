#' @export
setMethod(
  f = "proteinIntensityFiltering",
  signature = "ProteinQuantitativeData",
  definition = function(
    theObject,
    proteins_intensity_cutoff_percentile = NULL,
    proteins_proportion_of_samples_below_cutoff = NULL,
    core_utilisation = NULL
  ) {
    protein_quant_table <- theObject@protein_quant_table

    proteins_intensity_cutoff_percentile <- checkParamsObjectFunctionSimplify(
      theObject,
      "proteins_intensity_cutoff_percentile",
      NULL
    )
    proteins_proportion_of_samples_below_cutoff <- checkParamsObjectFunctionSimplify(
      theObject,
      "proteins_proportion_of_samples_below_cutoff",
      NULL
    )
    core_utilisation <- checkParamsObjectFunctionSimplify(
      theObject,
      "core_utilisation",
      NA
    )

    theObject <- updateParamInObject(theObject, "proteins_intensity_cutoff_percentile")
    theObject <- updateParamInObject(theObject, "proteins_proportion_of_samples_below_cutoff")
    theObject <- updateParamInObject(theObject, "core_utilisation")


    data_long_cln <- protein_quant_table |>
      pivot_longer(
        cols = !matches(theObject@protein_id_column),
        names_to = theObject@sample_id,
        values_to = "log_values"
      ) |>
      mutate(temp = "")

    min_peptide_intensity_threshold <- ceiling(quantile(data_long_cln$log_values, na.rm = TRUE, probs = c(proteins_intensity_cutoff_percentile)))[1]

    peptide_normalised_pif_cln <- peptideIntensityFilteringHelper(data_long_cln,
      min_peptide_intensity_threshold = min_peptide_intensity_threshold,
      proteins_proportion_of_samples_below_cutoff = proteins_proportion_of_samples_below_cutoff,
      protein_id_column = !!sym(theObject@protein_id_column),
      peptide_sequence_column = temp,
      peptide_quantity_column = log_values,
      core_utilisation = core_utilisation
    )


    theObject@protein_quant_table <- peptide_normalised_pif_cln |>
      dplyr::select(-temp) |>
      pivot_wider(id_cols = theObject@protein_id_column, names_from = !!sym(theObject@sample_id), values_from = log_values)

    theObject <- cleanDesignMatrix(theObject)

    updated_object <- theObject

    return(updated_object)
  }
)

