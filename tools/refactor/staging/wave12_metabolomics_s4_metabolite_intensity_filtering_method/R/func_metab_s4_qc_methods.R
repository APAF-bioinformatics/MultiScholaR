#' @title Metabolite Intensity Filtering Method for MetaboliteAssayData
#'
#' @description
#' Filters metabolites in *all* assays of a MetaboliteAssayData object.
#' It removes metabolites that have intensities below a certain percentile threshold
#' in a proportion of samples exceeding a defined cutoff. The threshold is calculated
#' independently for each assay.
#'
#' @describeIn metaboliteIntensityFiltering Method for MetaboliteAssayData
#'
#' @param theObject A MetaboliteAssayData object.
#' @param metabolites_intensity_cutoff_percentile See generic definition.
#' @param metabolites_proportion_of_samples_below_cutoff See generic definition.
#'
#' @importFrom dplyr pull select all_of across
#' @importFrom rlang sym
#' @importFrom stats quantile
#'
#' @return An updated MetaboliteAssayData object.
#' @export
setMethod(
    f = "metaboliteIntensityFiltering",
    signature = "MetaboliteAssayData",
    definition = function(theObject, metabolites_intensity_cutoff_percentile = NULL, metabolites_proportion_of_samples_below_cutoff = NULL) {
        # --- Parameter Resolution (Done once) ---
        config_intensity_percentile <- "metabolites_intensity_cutoff_percentile"
        raw_intensity_percentile <- checkParamsObjectFunctionSimplify(
            theObject,
            config_intensity_percentile,
            metabolites_intensity_cutoff_percentile
        )
        message("Raw intensity percentile from config/param: ", raw_intensity_percentile)
        cleaned_intensity_percentile <- trimws(sub("#.*$", "", raw_intensity_percentile))
        intensity_cutoff_percentile_final <- as.numeric(cleaned_intensity_percentile)

        config_proportion_cutoff <- "metabolites_proportion_of_samples_below_cutoff"
        raw_proportion_cutoff <- checkParamsObjectFunctionSimplify(
            theObject,
            config_proportion_cutoff,
            metabolites_proportion_of_samples_below_cutoff
        )
        message("Raw proportion cutoff from config/param: ", raw_proportion_cutoff)
        cleaned_proportion_cutoff <- trimws(sub("#.*$", "", raw_proportion_cutoff))
        proportion_of_samples_below_cutoff_final <- as.numeric(cleaned_proportion_cutoff)

        if (is.na(intensity_cutoff_percentile_final)) {
            stop("Failed to convert cleaned metabolites_intensity_cutoff_percentile ('", cleaned_intensity_percentile, "' from raw '", raw_intensity_percentile, "') to numeric. Check config.ini or parameter value.")
        }
        if (is.na(proportion_of_samples_below_cutoff_final)) {
            stop("Failed to convert cleaned metabolites_proportion_of_samples_below_cutoff ('", cleaned_proportion_cutoff, "' from raw '", raw_proportion_cutoff, "') to numeric. Check config.ini or parameter value.")
        }

        # --- Update Object Parameters (Done once) ---
        theObject <- updateParamInObject(theObject, config_intensity_percentile)
        theObject <- updateParamInObject(theObject, config_proportion_cutoff)

        # --- Process Each Assay in the List ---
        metabolite_id_col <- theObject@metabolite_id_column
        original_assay_list <- theObject@metabolite_data
        original_assay_names <- names(original_assay_list)

        if (length(original_assay_list) == 0) {
            warning("MetaboliteAssayData object has no assays in 'metabolite_data' slot. No filtering performed.")
            return(theObject)
        }

        # Iterate using indices
        filtered_assay_list <- lapply(seq_along(original_assay_list), function(i) {
            assay_table <- original_assay_list[[i]]
            # Determine assay name for messages (use index if no name)
            assay_name_for_msg <- if (!is.null(original_assay_names) && nzchar(original_assay_names[i])) {
                original_assay_names[i]
            } else {
                as.character(i) # Use index as fallback name
            }
            message("\nProcessing assay: ", assay_name_for_msg)

            if (!(metabolite_id_col %in% names(assay_table))) {
                warning("Metabolite ID column '", metabolite_id_col, "' not found in assay '", assay_name_for_msg, "'. Skipping this assay.")
                return(assay_table) # Return the original table if ID is missing
            }

            # Identify numeric sample columns for this assay
            sample_cols <- names(assay_table)[sapply(assay_table, is.numeric)]

            if (length(sample_cols) == 0) {
                warning("No numeric sample columns found in assay '", assay_name_for_msg, "'. Skipping filtering for this assay.")
                return(assay_table)
            }

            # Extract intensity values for this assay
            all_intensity_values <- assay_table |>
                dplyr::select(all_of(sample_cols)) |>
                unlist()

            if (length(all_intensity_values) == 0 || all(is.na(all_intensity_values))) {
                warning("No valid intensity values found in assay '", assay_name_for_msg, "' to calculate threshold. Skipping filtering for this assay.")
                return(assay_table)
            }

            # Calculate threshold specifically for this assay
            min_metabolite_intensity_threshold <- ceiling(quantile(all_intensity_values,
                na.rm = TRUE,
                probs = c(intensity_cutoff_percentile_final / 100)
            ))[1]

            message("Calculated minimum intensity threshold for assay '", assay_name_for_msg, "': ", min_metabolite_intensity_threshold)

            # Filter using Helper
            filtered_assay <- metaboliteIntensityFilteringHelper(
                assay_table = assay_table,
                min_metabolite_intensity_threshold = min_metabolite_intensity_threshold,
                metabolites_proportion_of_samples_below_cutoff = proportion_of_samples_below_cutoff_final,
                metabolite_id_column = metabolite_id_col
            )

            message("Filtered assay '", assay_name_for_msg, "'. Original rows: ", nrow(assay_table), ", Filtered rows: ", nrow(filtered_assay))
            return(filtered_assay)
        })

        # Restore original names if they existed
        if (!is.null(original_assay_names)) {
            names(filtered_assay_list) <- original_assay_names
        }

        # Assign the list of filtered assays back to the object
        theObject@metabolite_data <- filtered_assay_list

        # Optional: Call a generic cleanup/design matrix function if applicable
        # theObject <- cleanDesignMatrix(theObject) # If a generic method exists

        return(theObject)
    }
)

