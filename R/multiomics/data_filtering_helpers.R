# Helper functions for data filtering across omics types

#' Remove rows with high missing value percentages (Helper Function)
#'
#' Internal helper function performing the core logic for missing value filtering.
#'
#' @param input_table An input table with a column containing the row ID and the rest of the columns representing abundance values for each sample.
#' @param cols A tidyselect command to select the columns. This includes the functions starts_with(), ends_with(), contains(), matches(), and num_range()
#' @param design_matrix A data frame with a column containing the sample ID (as per the sample_id param) and the experimental group (as per the group param). Each row as the sample ID as row name in the data frame.
#' @param sample_id The name of the column in design_matrix table that has the sample ID.
#' @param row_id A unique ID for each row of the 'input_table' variable.
#' @param grouping_variable The name of the column in design_matrix table that has the experimental group.
#' @param groupwise_percentage_cutoff The maximum percentage of values below threshold allow in each group for a protein .
#' @param max_groups_percentage_cutoff The maximum percentage of groups allowed with too many samples with protein abundance values below threshold.
#' @param intensity_cutoff_percentile ***Generic Parameter Name*** The percentile of the intensity values to be used as the minimum threshold.
#' @param temporary_abundance_column The name of a temporary column to keep the abundance value you want to filter upon
#' @return A filtered tibble.
#'
#' @importFrom tidyr pivot_longer pivot_wider nest
#' @importFrom dplyr mutate left_join distinct group_by summarise ungroup anti_join full_join filter pull arrange join_by
#' @importFrom rlang sym !! enquo as_name as_string
#' @importFrom purrr map_chr
#' @importFrom stats quantile
#' @export
removeRowsWithMissingValuesPercentHelper <- function(
    input_table,
    cols,
    design_matrix,
    sample_id,
    row_id,
    grouping_variable,
    groupwise_percentage_cutoff = 1,
    max_groups_percentage_cutoff = 50
    # Renamed parameter:
    , intensity_cutoff_percentile = 1,
    temporary_abundance_column = "Abundance") {
    abundance_long <- input_table |>
        pivot_longer(
            cols = {{ cols }},
            names_to = as_string(as_name(enquo(sample_id))),
            values_to = temporary_abundance_column
        ) |>
        mutate({{ sample_id }} := purrr::map_chr({{ sample_id }}, as.character)) |>
        mutate(!!sym(temporary_abundance_column) := case_when(
            is.nan(!!sym(temporary_abundance_column)) ~ NA_real_,
            TRUE ~ !!sym(temporary_abundance_column)
        )) |>
        left_join(
            design_matrix |>
                mutate({{ sample_id }} := purrr::map_chr({{ sample_id }}, as.character)),
            by = join_by({{ sample_id }})
        )

    # Use the renamed parameter here:
    min_intensity_threshold <- ceiling(quantile(
        abundance_long |>
            dplyr::filter(!is.nan(!!sym(temporary_abundance_column)) & !is.infinite(!!sym(temporary_abundance_column))) |>
            dplyr::pull(!!sym(temporary_abundance_column)),
        na.rm = TRUE
        # Use the generic parameter name here:
        , probs = c(intensity_cutoff_percentile / 100)
    ))[1]

    count_values_per_group <- abundance_long |>
        distinct({{ sample_id }}, {{ grouping_variable }}) |>
        group_by({{ grouping_variable }}) |>
        summarise(num_per_group = n()) |>
        ungroup()

    count_values_missing_per_group <- abundance_long |>
        mutate(is_missing = ifelse(!is.na(!!sym(temporary_abundance_column))
        # Use the calculated threshold here:
        & !!sym(temporary_abundance_column) > min_intensity_threshold,
        0, 1
        )) |>
        group_by({{ row_id }}, {{ grouping_variable }}) |>
        summarise(num_missing_per_group = sum(is_missing)) |>
        ungroup()

    count_percent_missing_per_group <- count_values_missing_per_group |>
        full_join(count_values_per_group,
            by = join_by({{ grouping_variable }})
        ) |>
        mutate(perc_missing_per_group = num_missing_per_group / num_per_group * 100)

    total_num_of_groups <- count_values_per_group |> nrow()

    remove_rows_temp <- count_percent_missing_per_group |>
        dplyr::filter(groupwise_percentage_cutoff < perc_missing_per_group) |>
        group_by({{ row_id }}) |>
        summarise(percent = n() / total_num_of_groups * 100) |>
        ungroup() |>
        dplyr::filter(percent > max_groups_percentage_cutoff)

    filtered_tbl <- input_table |>
        dplyr::anti_join(remove_rows_temp, by = join_by({{ row_id }}))

    return(filtered_tbl)
}

# =============================================================================

#' Clean Design Matrix
#'
#' Filters the design matrix to retain only metadata for samples present
#' in the quantitative data and ensures consistent sample order.
#'
#' @param theObject The quantitative data object (e.g., ProteinQuantitativeData).
#'
#' @return The object with its design_matrix slot updated.
#'
#' @export
#' @importFrom methods setGeneric
setGeneric(
    name = "cleanDesignMatrix",
    def = function(theObject) {
        standardGeneric("cleanDesignMatrix")
    }
)

#' @describeIn cleanDesignMatrix Method for ProteinQuantitativeData
#' @export
#' @importFrom methods setMethod slot
#' @importFrom dplyr inner_join rename filter
#' @importFrom rlang sym !!
setMethod(
    f = "cleanDesignMatrix",
    signature = "ProteinQuantitativeData",
    definition = function(theObject) {
        # Get sample IDs from the quantitative data columns (excluding the feature ID column)
        samples_id_vector <- setdiff(colnames(slot(theObject, "protein_quant_table")), slot(theObject, "protein_id_column"))
        sample_id_col_name <- slot(theObject, "sample_id")

        slot(theObject, "design_matrix") <- data.frame(temp_sample_id = samples_id_vector) |>
            inner_join(slot(theObject, "design_matrix"),
                by = join_by(temp_sample_id == !!sym(sample_id_col_name))
            ) |>
            dplyr::rename(!!sym(sample_id_col_name) := "temp_sample_id") |>
            # Ensure only samples present in the vector remain
            dplyr::filter(!!sym(sample_id_col_name) %in% samples_id_vector)

        return(theObject)
    }
)

#' @describeIn cleanDesignMatrix Method for MetaboliteQuantitativeData
#' @export
#' @importFrom methods setMethod slot
#' @importFrom dplyr inner_join rename filter
#' @importFrom rlang sym !!
setMethod(
    f = "cleanDesignMatrix",
    signature = "MetaboliteQuantitativeData",
    definition = function(theObject) {
        # Get sample IDs from the quantitative data columns (excluding the feature ID column)
        samples_id_vector <- setdiff(colnames(slot(theObject, "metabolite_quant_table")), slot(theObject, "metabolite_id_column"))
        sample_id_col_name <- slot(theObject, "sample_id")

        slot(theObject, "design_matrix") <- data.frame(temp_sample_id = samples_id_vector) |>
            inner_join(slot(theObject, "design_matrix"),
                by = join_by(temp_sample_id == !!sym(sample_id_col_name))
            ) |>
            dplyr::rename(!!sym(sample_id_col_name) := "temp_sample_id") |>
            # Ensure only samples present in the vector remain
            dplyr::filter(!!sym(sample_id_col_name) %in% samples_id_vector)

        return(theObject)
    }
)

#' @describeIn cleanDesignMatrix Method for TranscriptQuantitativeData
#' @export
#' @importFrom methods setMethod slot
#' @importFrom dplyr inner_join rename filter
#' @importFrom rlang sym !!
setMethod(
    f = "cleanDesignMatrix",
    signature = "TranscriptQuantitativeData",
    definition = function(theObject) {
        # Get sample IDs from the quantitative data columns (excluding the feature ID column)
        samples_id_vector <- setdiff(colnames(slot(theObject, "transcript_quant_table")), slot(theObject, "transcript_id_column"))
        sample_id_col_name <- slot(theObject, "sample_id")

        slot(theObject, "design_matrix") <- data.frame(temp_sample_id = samples_id_vector) |>
            inner_join(slot(theObject, "design_matrix"),
                by = join_by(temp_sample_id == !!sym(sample_id_col_name))
            ) |>
            dplyr::rename(!!sym(sample_id_col_name) := "temp_sample_id") |>
            # Ensure only samples present in the vector remain
            dplyr::filter(!!sym(sample_id_col_name) %in% samples_id_vector)

        return(theObject)
    }
)
