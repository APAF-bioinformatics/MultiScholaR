#' @title Remove Rows with High Missing Value Percentages
#' @description Filters rows (features) based on the percentage of missing values
#'   within specified groups.
#'
#' @param theObject The quantitative data object.
#' @param ruv_grouping_variable Character string. The column name in the design matrix
#'   used for grouping samples to assess missingness. Defaults to the value
#'   in the object's arguments or NULL.
#' @param groupwise_percentage_cutoff Numeric (0-100). The maximum percentage of
#'   missing values allowed for a feature within any single group defined by
#'   `ruv_grouping_variable`. Features exceeding this in any group might be removed
#'   based on `max_groups_percentage_cutoff`. Defaults to the value in the
#'   object's arguments or 50.
#' @param max_groups_percentage_cutoff Numeric (0-100). The maximum percentage of
#'   groups in which a feature can exceed the `groupwise_percentage_cutoff`
#'   before being removed. Defaults to the value in the object's arguments or 50.
#' @param intensity_cutoff_percentile Numeric (0-100). An intensity percentile used by the
#'   underlying helper function, potentially to define low-abundance values that
#'   might be treated as missing. The specific config parameter name (e.g.,
#'   `proteins_intensity_cutoff_percentile`) is looked up internally based on the
#'   object's class. Defaults to the value in the object's arguments or 1.
#' @param ... Additional arguments (currently unused).
#'
#' @return The input object (`theObject`) with rows filtered.
#'
#' @export
#' @importFrom methods setGeneric
setGeneric(
    name = "removeRowsWithMissingValuesPercent",
    def = function(
        theObject,
        ruv_grouping_variable = NULL,
        groupwise_percentage_cutoff = NULL,
        max_groups_percentage_cutoff = NULL
        # Use a generic parameter name here in the generic definition
        , intensity_cutoff_percentile = NULL,
        ...) {
        standardGeneric("removeRowsWithMissingValuesPercent")
    }
    # Signature only needs the main object for dispatch usually
    , signature = c("theObject")
)


#' @describeIn removeRowsWithMissingValuesPercent Method for ProteinQuantitativeData
#' @export
#' @importFrom methods setMethod slot slot<- is
#' @importFrom rlang sym !! matches
setMethod(
    f = "removeRowsWithMissingValuesPercent",
    signature = "ProteinQuantitativeData",
    definition = function(theObject,
                          ruv_grouping_variable = NULL,
                          groupwise_percentage_cutoff = NULL,
                          max_groups_percentage_cutoff = NULL,
                          intensity_cutoff_percentile = NULL) { # Match generic arg name

        protein_quant_table <- slot(theObject, "protein_quant_table")
        protein_id_column <- slot(theObject, "protein_id_column")
        design_matrix <- slot(theObject, "design_matrix")
        sample_id <- slot(theObject, "sample_id")
        # replicate_group_column <- slot(theObject, "technical_replicate_id") # If needed

        # --- Parameter Resolution ---
        # Use specific config names for lookup inside the method
        config_ruv_grouping <- "ruv_grouping_variable"
        config_groupwise_cutoff <- "groupwise_percentage_cutoff"
        config_max_groups_cutoff <- "max_groups_percentage_cutoff"
        config_intensity_percentile <- "proteins_intensity_cutoff_percentile" # Specific to Protein

        ruv_grouping_variable_final <- checkParamsObjectFunctionSimplify(theObject, config_ruv_grouping, ruv_grouping_variable)
        groupwise_percentage_cutoff_final <- checkParamsObjectFunctionSimplify(theObject, config_groupwise_cutoff, groupwise_percentage_cutoff, 50)
        max_groups_percentage_cutoff_final <- checkParamsObjectFunctionSimplify(theObject, config_max_groups_cutoff, max_groups_percentage_cutoff, 50)
        intensity_cutoff_percentile_final <- checkParamsObjectFunctionSimplify(theObject, config_intensity_percentile, intensity_cutoff_percentile, 1)

        # --- Update Object Params (Optional) ---
        # theObject <- updateParamInObject(theObject, config_ruv_grouping, ruv_grouping_variable_final)
        # theObject <- updateParamInObject(theObject, config_groupwise_cutoff, groupwise_percentage_cutoff_final)
        # ... etc ...

        # --- Filtering ---
        slot(theObject, "protein_quant_table") <- removeRowsWithMissingValuesPercentHelper(
            protein_quant_table,
            cols = !matches(protein_id_column),
            design_matrix = design_matrix,
            sample_id = !!sym(sample_id),
            row_id = !!sym(protein_id_column),
            grouping_variable = !!sym(ruv_grouping_variable_final),
            groupwise_percentage_cutoff = groupwise_percentage_cutoff_final,
            max_groups_percentage_cutoff = max_groups_percentage_cutoff_final
            # *** Pass the resolved value using the GENERIC parameter name ***
            , intensity_cutoff_percentile = intensity_cutoff_percentile_final,
            temporary_abundance_column = "Log_Abundance"
        ) # Adjust name if needed

        theObject <- cleanDesignMatrix(theObject)
        return(theObject)
    }
)

# Helper functions assumed to exist elsewhere:
# checkParamsObjectFunctionSimplify <- function(theObject, config_param_name, function_arg_value, default = NULL) { ... }
# updateParamInObject <- function(theObject, config_param_name, value_to_store) { ... }
# removeRowsWithMissingValuesPercentHelper <- function(quant_table, cols, design_matrix, sample_id, row_id, grouping_variable, groupwise_percentage_cutoff, max_groups_percentage_cutoff, proteins_intensity_cutoff_percentile, temporary_abundance_column) { ... }
# cleanDesignMatrix <- function(theObject) { ... }

#' @describeIn removeRowsWithMissingValuesPercent Method for MetaboliteQuantitativeData
#' @export
#' @importFrom methods setMethod slot slot<- is
#' @importFrom rlang sym !! matches
setMethod(
    f = "removeRowsWithMissingValuesPercent",
    signature = "MetaboliteQuantitativeData" # Changed Signature
    , definition = function(theObject,
                            ruv_grouping_variable = NULL,
                            groupwise_percentage_cutoff = NULL,
                            max_groups_percentage_cutoff = NULL,
                            intensity_cutoff_percentile = NULL) { # Match generic arg name

        metabolite_quant_table <- slot(theObject, "metabolite_quant_table") # Changed Slot
        metabolite_id_column <- slot(theObject, "metabolite_id_column") # Changed Slot
        design_matrix <- slot(theObject, "design_matrix")
        sample_id <- slot(theObject, "sample_id")
        # replicate_group_column <- slot(theObject, "technical_replicate_id") # If needed

        # --- Parameter Resolution ---
        config_ruv_grouping <- "ruv_grouping_variable"
        config_groupwise_cutoff <- "groupwise_percentage_cutoff"
        config_max_groups_cutoff <- "max_groups_percentage_cutoff"
        # *** Use the correct config name for metabolites ***
        config_intensity_percentile <- "metabolites_intensity_cutoff_percentile"

        ruv_grouping_variable_final <- checkParamsObjectFunctionSimplify(theObject, config_ruv_grouping, ruv_grouping_variable)
        groupwise_percentage_cutoff_final <- checkParamsObjectFunctionSimplify(theObject, config_groupwise_cutoff, groupwise_percentage_cutoff, 50)
        max_groups_percentage_cutoff_final <- checkParamsObjectFunctionSimplify(theObject, config_max_groups_cutoff, max_groups_percentage_cutoff, 50)
        intensity_cutoff_percentile_final <- checkParamsObjectFunctionSimplify(theObject, config_intensity_percentile, intensity_cutoff_percentile, 1) # Default might differ

        # --- Update Object Params (Optional) ---
        # ...

        # --- Filtering ---
        # *** Update the slot being assigned to ***
        slot(theObject, "metabolite_quant_table") <- removeRowsWithMissingValuesPercentHelper(
            metabolite_quant_table # Changed table
            ,
            cols = !matches(metabolite_id_column) # Changed ID col
            , design_matrix = design_matrix,
            sample_id = !!sym(sample_id),
            row_id = !!sym(metabolite_id_column) # Changed ID col
            , grouping_variable = !!sym(ruv_grouping_variable_final),
            groupwise_percentage_cutoff = groupwise_percentage_cutoff_final,
            max_groups_percentage_cutoff = max_groups_percentage_cutoff_final
            # *** Pass the resolved value using the GENERIC parameter name ***
            , intensity_cutoff_percentile = intensity_cutoff_percentile_final,
            temporary_abundance_column = "Log_Abundance"
        ) # Adjust name if needed

        theObject <- cleanDesignMatrix(theObject)
        return(theObject)
    }
)


#' @describeIn removeRowsWithMissingValuesPercent Method for TranscriptQuantitativeData
#' @export
#' @importFrom methods setMethod slot slot<- is
#' @importFrom rlang sym !! matches
setMethod(
    f = "removeRowsWithMissingValuesPercent",
    signature = "TranscriptQuantitativeData" # Changed Signature
    , definition = function(theObject,
                            ruv_grouping_variable = NULL,
                            groupwise_percentage_cutoff = NULL,
                            max_groups_percentage_cutoff = NULL,
                            intensity_cutoff_percentile = NULL) { # Match generic arg name

        transcript_quant_table <- slot(theObject, "transcript_quant_table") # Changed Slot
        transcript_id_column <- slot(theObject, "transcript_id_column") # Changed Slot
        design_matrix <- slot(theObject, "design_matrix")
        sample_id <- slot(theObject, "sample_id")
        # replicate_group_column <- slot(theObject, "technical_replicate_id") # If needed

        # --- Parameter Resolution ---
        config_ruv_grouping <- "ruv_grouping_variable"
        config_groupwise_cutoff <- "groupwise_percentage_cutoff"
        config_max_groups_cutoff <- "max_groups_percentage_cutoff"
        # *** Use the correct config name for transcripts ***
        config_intensity_percentile <- "transcripts_intensity_cutoff_percentile"

        ruv_grouping_variable_final <- checkParamsObjectFunctionSimplify(theObject, config_ruv_grouping, ruv_grouping_variable)
        groupwise_percentage_cutoff_final <- checkParamsObjectFunctionSimplify(theObject, config_groupwise_cutoff, groupwise_percentage_cutoff, 50)
        max_groups_percentage_cutoff_final <- checkParamsObjectFunctionSimplify(theObject, config_max_groups_cutoff, max_groups_percentage_cutoff, 50)
        intensity_cutoff_percentile_final <- checkParamsObjectFunctionSimplify(theObject, config_intensity_percentile, intensity_cutoff_percentile, 1) # Default might differ

        # --- Update Object Params (Optional) ---
        # ...

        # --- Filtering ---
        # *** Update the slot being assigned to ***
        slot(theObject, "transcript_quant_table") <- removeRowsWithMissingValuesPercentHelper(
            transcript_quant_table # Changed table
            ,
            cols = !matches(transcript_id_column) # Changed ID col
            , design_matrix = design_matrix,
            sample_id = !!sym(sample_id),
            row_id = !!sym(transcript_id_column) # Changed ID col
            , grouping_variable = !!sym(ruv_grouping_variable_final),
            groupwise_percentage_cutoff = groupwise_percentage_cutoff_final,
            max_groups_percentage_cutoff = max_groups_percentage_cutoff_final
            # *** Pass the resolved value using the GENERIC parameter name ***
            , intensity_cutoff_percentile = intensity_cutoff_percentile_final,
            temporary_abundance_column = "Log_Abundance"
        ) # Adjust name if needed

        theObject <- cleanDesignMatrix(theObject)
        return(theObject)
    }
)
