#' @export
setMethod(
  f = "filterSamplesByProteinCorrelationThreshold",
  signature = "ProteinQuantitativeData",
  definition = function(theObject, pearson_correlation_per_pair = NULL, min_pearson_correlation_threshold = NULL) {
    message("+===========================================================================+")
    message("|  DEBUG66: Entering filterSamplesByProteinCorrelationThreshold             |")
    message("+===========================================================================+")

    # Memory tracking - Entry
    entry_mem <- checkMemoryBoth("Entry", context = "filterSamplesByProteinCorrelationThreshold")

    pearson_correlation_per_pair <- checkParamsObjectFunctionSimplify(theObject,
      "pearson_correlation_per_pair",
      default_value = NULL
    )
    min_pearson_correlation_threshold <- checkParamsObjectFunctionSimplify(theObject,
      "min_pearson_correlation_threshold",
      default_value = 0.75
    )

    theObject <- updateParamInObject(theObject, "pearson_correlation_per_pair")
    theObject <- updateParamInObject(theObject, "min_pearson_correlation_threshold")

    # Memory tracking - Before helper
    pre_helper_mem <- checkMemoryBoth("Before helper", context = "filterSamplesByProteinCorrelationThreshold")

    filtered_table <- filterSamplesByProteinCorrelationThresholdHelper(
      pearson_correlation_per_pair,
      protein_intensity_table = theObject@protein_quant_table,
      min_pearson_correlation_threshold = min_pearson_correlation_threshold,
      filename_column_x = !!sym(paste0(theObject@sample_id, ".x")),
      filename_column_y = !!sym(paste0(theObject@sample_id, ".y")),
      protein_id_column = theObject@protein_id_column,
      correlation_column = pearson_correlation
    )

    # Memory tracking - After helper
    reportMemoryDelta(pre_helper_mem, "helper function", context = "filterSamplesByProteinCorrelationThreshold")

    theObject@protein_quant_table <- filtered_table

    # Memory tracking - Before cleanDesignMatrix
    pre_clean_mem <- checkMemoryBoth("Before cleanDesignMatrix", context = "filterSamplesByProteinCorrelationThreshold")

    theObject <- cleanDesignMatrix(theObject)

    # Memory tracking - After cleanDesignMatrix
    reportMemoryDelta(pre_clean_mem, "cleanDesignMatrix", context = "filterSamplesByProteinCorrelationThreshold")

    # Memory tracking - Exit
    reportMemoryDelta(entry_mem, "TOTAL filterSamplesByProteinCorrelationThreshold", context = "filterSamplesByProteinCorrelationThreshold")

    message("+===========================================================================+")
    message("|  DEBUG66: Exiting filterSamplesByProteinCorrelationThreshold              |")
    message("+===========================================================================+")

    theObject
  }
)

