# ----------------------------------------------------------------------------
# getFilteringProgressLipidomics
# ----------------------------------------------------------------------------
#' @title Initialize or Retrieve Global Lipidomics Filtering Progress Object
#' @description Checks for a global object named `filtering_progress_lipidomics`
#'              of class `FilteringProgressLipidomics`. If it doesn't exist,
#'              it creates and assigns a new one to the global environment.
#'
#' @return The global `FilteringProgressLipidomics` object.
#' @keywords internal
#' @noRd
#' @export
getFilteringProgressLipidomics <- function() {
    if (!exists("filtering_progress_lipidomics", envir = .GlobalEnv)) {
        filtering_progress_lipidomics <- new("FilteringProgressLipidomics")
        assign("filtering_progress_lipidomics", filtering_progress_lipidomics, envir = .GlobalEnv)
        message("Initialized global 'filtering_progress_lipidomics' object.") # Optional message
    }
    get("filtering_progress_lipidomics", envir = .GlobalEnv)
}

# ----------------------------------------------------------------------------
# updateFilteringProgressLipidomics
# ----------------------------------------------------------------------------
#' @title Update the Global Lipidomics Filtering Progress Object
#' @description Modifies the global `filtering_progress_lipidomics` object
#'              by adding or overwriting data for a specific step.
#'
#' @param prog_met The `FilteringProgressLipidomics` object (retrieved globally).
#' @param step_name The name of the step to add or update.
#' @param current_assay_names Character vector of assay names for this step.
#' @param metrics_list A nested list containing the calculated metrics for each assay
#'                      for the current step.
#' @param total_lipids The total unique lipids calculated across assays for this step.
#' @param overwrite Logical, whether to overwrite if `step_name` exists.
#'
#' @return The updated `FilteringProgressLipidomics` object (invisibly).
#'         Has the side effect of updating the global object.
#' @keywords internal
#' @noRd
#' @export
updateFilteringProgressLipidomics <- function(prog_met,
                                              step_name,
                                              current_assay_names,
                                              metrics_list,
                                              total_lipids,
                                              overwrite = FALSE) {
    if (step_name %in% prog_met@steps) {
        if (!overwrite) {
            stop("Step name '", step_name, "' already exists in filtering_progress_lipidomics. Use overwrite = TRUE to replace it.")
        }
        idx <- which(prog_met@steps == step_name)

        # Overwrite existing data
        prog_met@assay_names[[idx]] <- current_assay_names
        prog_met@n_lipids_per_assay[[idx]] <- lapply(metrics_list, `[[`, "n_lipids")
        prog_met@n_lipids_total[idx] <- total_lipids
        prog_met@detected_per_sample[[idx]] <- lapply(metrics_list, `[[`, "detected_per_sample")
        prog_met@missingness_per_assay[[idx]] <- lapply(metrics_list, `[[`, "missingness")
        prog_met@sum_intensity_per_sample[[idx]] <- lapply(metrics_list, `[[`, "sum_intensity_per_sample")
        prog_met@cv_distribution_per_assay[[idx]] <- lapply(metrics_list, `[[`, "cv_distribution")
        prog_met@is_metrics_per_assay[[idx]] <- lapply(metrics_list, `[[`, "is_metrics")
    } else {
        # Append new data
        prog_met@steps <- c(prog_met@steps, step_name)
        prog_met@assay_names <- c(prog_met@assay_names, list(current_assay_names))
        prog_met@n_lipids_per_assay <- c(prog_met@n_lipids_per_assay, list(lapply(metrics_list, `[[`, "n_lipids")))
        prog_met@n_lipids_total <- c(prog_met@n_lipids_total, total_lipids)
        prog_met@detected_per_sample <- c(prog_met@detected_per_sample, list(lapply(metrics_list, `[[`, "detected_per_sample")))
        prog_met@missingness_per_assay <- c(prog_met@missingness_per_assay, list(lapply(metrics_list, `[[`, "missingness")))
        prog_met@sum_intensity_per_sample <- c(prog_met@sum_intensity_per_sample, list(lapply(metrics_list, `[[`, "sum_intensity_per_sample")))
        prog_met@cv_distribution_per_assay <- c(prog_met@cv_distribution_per_assay, list(lapply(metrics_list, `[[`, "cv_distribution")))
        prog_met@is_metrics_per_assay <- c(prog_met@is_metrics_per_assay, list(lapply(metrics_list, `[[`, "is_metrics")))
    }

    # Update the global object
    assign("filtering_progress_lipidomics", prog_met, envir = .GlobalEnv)
    invisible(prog_met)
}

# ----------------------------------------------------------------------------
# countUniqueLipids
# ----------------------------------------------------------------------------
#' @title Count Unique Lipids in an Assay
#' @description Counts the number of unique, non-NA identifiers in the specified
#'              lipid ID column of an assay tibble.
#'
#' @param assay_data A tibble/data.frame for one assay.
#' @param lipid_id_col Character string, the name of the column containing
#'                         lipid identifiers.
#'
#' @return A single numeric value: the count of unique lipids.
#'         Returns 0 if the column doesn't exist or has no non-NA values.
#'
#' @importFrom dplyr pull distinct n
#' @keywords internal
#' @noRd
#' @export
countUniqueLipids <- function(assay_data, lipid_id_col) {
    if (!lipid_id_col %in% colnames(assay_data)) {
        warning("Lipid ID column '", lipid_id_col, "' not found in assay data.")
        return(0)
    }

    # Count unique non-NA lipid IDs - use base R subsetting (dplyr::pull doesn't support .data[[var]])
    lipid_ids <- assay_data[[lipid_id_col]]
    lipid_ids <- stats::na.omit(lipid_ids)
    unique_ids <- dplyr::distinct(data.frame(id = lipid_ids))
    return(nrow(unique_ids))
}

# ----------------------------------------------------------------------------
# countLipidsPerSample
# ----------------------------------------------------------------------------
#' @title Count Detected Lipids per Sample
#' @description For each sample column in an assay's quantitative data,
#'              counts the number of non-missing, non-zero intensity values.
#'
#' @param assay_data A tibble/data.frame for one assay.
#' @param sample_id_col ***Currently unused.***
#' @param lipid_id_col ***Currently unused.***
#' @param sample_columns Character vector of explicit sample column names from design matrix.
#'
#' @return A data frame with columns 'Run' (sample name) and 'n_detected'.
#'
#' @importFrom tidyr pivot_longer
#' @importFrom dplyr group_by summarise n filter
#' @keywords internal
#' @noRd
#' @export
countLipidsPerSample <- function(assay_data, sample_id_col, lipid_id_col, sample_columns = NULL) {
    quant_info <- getLipidQuantData(assay_data, sample_columns = sample_columns)
    quant_data <- quant_info$quant_data
    sample_names <- quant_info$sample_names

    if (length(sample_names) == 0) {
        return(data.frame(Run = character(), n_detected = integer()))
    }

    # Add a temporary row identifier if no lipid ID column exists or is non-unique
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
# calculateLipidMissingness
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
calculateLipidMissingness <- function(assay_data, sample_id_col, sample_columns = NULL) {
    # Get only the quantitative columns (sample columns)
    quant_info <- getLipidQuantData(assay_data, sample_columns = sample_columns)
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
# calculateLipidSumIntensityPerSample
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
calculateLipidSumIntensityPerSample <- function(assay_data, sample_id_col, sample_columns = NULL) {
    quant_info <- getLipidQuantData(assay_data, sample_columns = sample_columns)
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
# calculateTotalUniqueLipidsAcrossAssays
# ----------------------------------------------------------------------------
#' @title Calculate Total Unique Lipids Across Assays
#' @description Finds all unique lipid IDs present across a list of assay tibbles.
#'
#' @param assay_list A list of tibbles/data.frames, one for each assay.
#' @param lipid_id_col Character string, the name of the column containing
#'                         lipid identifiers in each assay tibble.
#'
#' @return A single numeric value: the total count of unique lipid IDs
#'         across all provided assays.
#'
#' @importFrom purrr map
#' @importFrom dplyr distinct n pull
#' @keywords internal
#' @noRd
#' @export
calculateTotalUniqueLipidsAcrossAssays <- function(assay_list, lipid_id_col) {
    if (length(assay_list) == 0) {
        return(0)
    }

    # Extract the lipid ID column from each assay, handling missing columns
    # Use base R subsetting - dplyr::pull doesn't support {{ var }} for string column names
    all_ids <- purrr::map(assay_list, ~ {
        if (lipid_id_col %in% colnames(.x)) {
            .x[[lipid_id_col]]
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

