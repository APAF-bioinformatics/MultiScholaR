# ----------------------------------------------------------------------------
# countUniqueMetabolites
# ----------------------------------------------------------------------------
#' @title Count Unique Metabolites in an Assay
#' @description Counts the number of unique, non-NA identifiers in the specified
#'              metabolite ID column of an assay tibble.
#'
#' @param assay_data A tibble/data.frame for one assay.
#' @param metabolite_id_col Character string, the name of the column containing
#'                         metabolite identifiers.
#'
#' @return A single numeric value: the count of unique metabolites.
#'         Returns 0 if the column doesn't exist or has no non-NA values.
#'
#' @importFrom dplyr pull distinct n
#' @keywords internal
#' @noRd
#' @export
countUniqueMetabolites <- function(assay_data, metabolite_id_col) {
    if (!metabolite_id_col %in% colnames(assay_data)) {
        warning("Metabolite ID column '", metabolite_id_col, "' not found in assay data.")
        return(0)
    }

    # Count unique non-NA metabolite IDs - use base R subsetting (dplyr::pull doesn't support .data[[var]])
    metabolite_ids <- assay_data[[metabolite_id_col]]
    metabolite_ids <- stats::na.omit(metabolite_ids)
    unique_ids <- dplyr::distinct(data.frame(id = metabolite_ids))
    return(nrow(unique_ids))
}

# ----------------------------------------------------------------------------
# countMetabolitesPerSample
# ----------------------------------------------------------------------------
#' @title Count Detected Metabolites per Sample
#' @description For each sample column in an assay's quantitative data,
#'              counts the number of non-missing, non-zero intensity values.
#'
#' @param assay_data A tibble/data.frame for one assay.
#' @param sample_id_col ***Currently unused.***
#' @param metabolite_id_col ***Currently unused.***
#' @param sample_columns Character vector of explicit sample column names from design matrix.
#'
#' @return A data frame with columns 'Run' (sample name) and 'n_detected'.
#'
#' @importFrom tidyr pivot_longer
#' @importFrom dplyr group_by summarise n filter
#' @keywords internal
#' @noRd
#' @export
countMetabolitesPerSample <- function(assay_data, sample_id_col, metabolite_id_col, sample_columns = NULL) {
    quant_info <- getMetaboliteQuantData(assay_data, sample_columns = sample_columns)
    quant_data <- quant_info$quant_data
    sample_names <- quant_info$sample_names

    if (length(sample_names) == 0) {
        return(data.frame(Run = character(), n_detected = integer()))
    }

    # Add a temporary row identifier if no metabolite ID column exists or is non-unique
    # This ensures pivot_longer works correctly even without a unique ID col.
    quant_data$..temp_row_id.. <- seq_len(nrow(quant_data))

    quant_data |>
        tidyr::pivot_longer(
            cols = dplyr::all_of(sample_names),
            names_to = "Run",
            values_to = "Intensity"
        ) |>
        # Define 'detected' as non-NA and greater than 0 (adjust if needed)
        dplyr::filter(!is.na(.data$Intensity) & .data$Intensity > 0) |>
        dplyr::group_by(.data$Run) |>
        dplyr::summarise(n_detected = dplyr::n(), .groups = "drop")
}

# ----------------------------------------------------------------------------
# calculateMissingness
# ----------------------------------------------------------------------------
#' @title Calculate Overall Missingness Percentage
#' @description Calculates the percentage of missing values (NA or zero)
#'              in the quantitative portion of an assay tibble.
#'
#' @param assay_data A tibble/data.frame for one assay.
#' @param sample_id_col ***Currently unused.***
#' @param sample_columns Character vector of explicit sample column names from design matrix.
#'
#' @return A single numeric value: the percentage of missing data points.
#'
#' @keywords internal
#' @noRd
#' @export
calculateMissingness <- function(assay_data, sample_id_col, sample_columns = NULL) {
    # Get only the quantitative columns (sample columns)
    quant_info <- getMetaboliteQuantData(assay_data, sample_columns = sample_columns)
    quant_data <- quant_info$quant_data
    sample_names <- quant_info$sample_names

    # Validate input
    if (nrow(quant_data) == 0 || ncol(quant_data) == 0 || length(sample_names) == 0) {
        warning("No valid data for missingness calculation")
        return(NA_real_)
    }

    # Exclude any non-sample columns
    if ("..temp_row_id.." %in% colnames(quant_data)) {
        quant_data <- quant_data[, setdiff(colnames(quant_data), "..temp_row_id.."), drop = FALSE]
    }

    # Count missing (NA or zero) values in all sample columns
    total_cells <- nrow(quant_data) * length(sample_names)
    missing_values <- 0

    for (col in sample_names) {
        missing_values <- missing_values + sum(is.na(quant_data[[col]]) | quant_data[[col]] == 0)
    }

    # Calculate percentage
    missing_pct <- (missing_values / total_cells) * 100

    # Debug output to verify calculation
    message(
        "DEBUG: Missing values: ", missing_values, ", Total cells: ", total_cells,
        ", Percentage: ", missing_pct, "%"
    )

    return(missing_pct)
}

# ----------------------------------------------------------------------------
# calculateSumIntensityPerSample
# ----------------------------------------------------------------------------
#' @title Calculate Sum Intensity per Sample (TIC Proxy)
#' @description Calculates the sum of all intensities for each sample column
#'              in the quantitative portion of an assay tibble.
#'
#' @param assay_data A tibble/data.frame for one assay.
#' @param sample_id_col ***Currently unused.***
#' @param sample_columns Character vector of explicit sample column names from design matrix.
#'
#' @return A data frame with columns 'Run' (sample name) and 'sum_intensity'.
#'
#' @importFrom tidyr pivot_longer
#' @importFrom dplyr group_by summarise
#' @keywords internal
#' @noRd
#' @export
calculateSumIntensityPerSample <- function(assay_data, sample_id_col, sample_columns = NULL) {
    quant_info <- getMetaboliteQuantData(assay_data, sample_columns = sample_columns)
    quant_data <- quant_info$quant_data
    sample_names <- quant_info$sample_names

    if (length(sample_names) == 0) {
        return(data.frame(Run = character(), sum_intensity = numeric()))
    }

    # Add a temporary row identifier
    quant_data$..temp_row_id.. <- seq_len(nrow(quant_data))

    quant_data |>
        tidyr::pivot_longer(
            cols = dplyr::all_of(sample_names),
            names_to = "Run",
            values_to = "Intensity"
        ) |>
        dplyr::group_by(.data$Run) |>
        # Sum intensities, replacing NA with 0 for the sum
        dplyr::summarise(sum_intensity = sum(.data$Intensity, na.rm = TRUE), .groups = "drop")
}

# ----------------------------------------------------------------------------
# calculateTotalUniqueMetabolitesAcrossAssays
# ----------------------------------------------------------------------------
#' @title Calculate Total Unique Metabolites Across Assays
#' @description Finds all unique metabolite IDs present across a list of assay tibbles.
#'
#' @param assay_list A list of tibbles/data.frames, one for each assay.
#' @param metabolite_id_col Character string, the name of the column containing
#'                         metabolite identifiers in each assay tibble.
#'
#' @return A single numeric value: the total count of unique metabolite IDs
#'         across all provided assays.
#'
#' @importFrom purrr map
#' @importFrom dplyr distinct n pull
#' @keywords internal
#' @noRd
#' @export
calculateTotalUniqueMetabolitesAcrossAssays <- function(assay_list, metabolite_id_col) {
    if (length(assay_list) == 0) {
        return(0)
    }

    # Extract the metabolite ID column from each assay, handling missing columns
    # Use base R subsetting - dplyr::pull doesn't support {{ var }} for string column names
    all_ids <- purrr::map(assay_list, ~ {
        if (metabolite_id_col %in% colnames(.x)) {
            .x[[metabolite_id_col]]
        } else {
            NULL # Return NULL if column is missing
        }
    }) |>
        unlist() # Combine all IDs into a single vector

    # Remove NAs
    all_ids <- stats::na.omit(all_ids)

    # If no valid IDs found, return 0
    if (length(all_ids) == 0) {
        return(0)
    }

    # Count unique IDs using a more explicit approach
    unique_ids <- dplyr::distinct(data.frame(id = all_ids))
    return(nrow(unique_ids))
}

