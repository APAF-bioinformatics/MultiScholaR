# ----------------------------------------------------------------------------
# calculateMetaboliteCVs
# ----------------------------------------------------------------------------
#' @title Calculate Within-Group Metabolite CVs
#' @description Calculates the Coefficient of Variation (CV) for each metabolite
#'              across replicate samples within each experimental group.
#'
#' @param assay_data A tibble/data.frame for one assay.
#' @param design_matrix A data frame linking samples to experimental design
#'                      (must contain sample_id_col and group_id_col).
#' @param group_id_col Character string, the name of the grouping column in design_matrix.
#' @param replicate_id_col Character string, the name of the replicate identifier column
#'                         in design_matrix (used indirectly to ensure grouping works).
#' @param sample_id_col Character string, the name of the sample ID column in design_matrix
#'                      (matching column names in assay_data).
#' @param metabolite_id_col Character string, the name of the metabolite ID column in assay_data.
#' @param sample_columns Character vector of explicit sample column names from design matrix.
#'
#' @return A data frame with columns 'metabolite_id', 'group', and 'cv'.
#'         Returns an empty data frame if inputs are invalid or calculations fail.
#'
#' @importFrom tidyr pivot_longer
#' @importFrom dplyr left_join group_by summarise filter select rename all_of
#' @importFrom rlang sym
#' @importFrom stats sd na.omit
#' @keywords internal
#' @noRd
#' @export
calculateMetaboliteCVs <- function(assay_data,
                                   design_matrix,
                                   group_id_col,
                                   replicate_id_col,
                                   sample_id_col,
                                   metabolite_id_col,
                                   sample_columns = NULL) {
    # --- Input Validation --- #
    required_design_cols <- c(sample_id_col, group_id_col)
    if (!all(required_design_cols %in% colnames(design_matrix))) {
        warning("Design matrix missing required columns: ", paste(setdiff(required_design_cols, colnames(design_matrix)), collapse = ", "))
        return(data.frame(metabolite_id = character(), group = character(), cv = numeric()))
    }
    if (!metabolite_id_col %in% colnames(assay_data)) {
        warning("Metabolite ID column '", metabolite_id_col, "' not found in assay data.")
        return(data.frame(metabolite_id = character(), group = character(), cv = numeric()))
    }

    # --- Data Preparation --- #
    quant_info <- getMetaboliteQuantData(assay_data, sample_columns = sample_columns)
    quant_data <- quant_info$quant_data
    sample_names <- quant_info$sample_names
    annotation_data <- quant_info$annotation_data

    if (length(sample_names) == 0 || nrow(quant_data) == 0) {
        warning("No sample columns or no data rows found for CV calculation")
        return(data.frame(metabolite_id = character(), group = character(), cv = numeric()))
    }

    # Ensure metabolite IDs are present for joining
    quant_data[[metabolite_id_col]] <- annotation_data[[metabolite_id_col]]

    # Select relevant columns from design matrix and ensure consistent types
    # Use !!rlang::sym() for string column names - {{ }} is for defused symbols only
    design_subset <- design_matrix |>
        dplyr::select(dplyr::all_of(required_design_cols)) |>
        dplyr::rename(Run = !!rlang::sym(sample_id_col), group = !!rlang::sym(group_id_col)) |>
        dplyr::mutate(Run = as.character(Run)) # Convert Run to character to match pivoted data

    # --- Debug output to check groups --- #
    message("Groups in design matrix: ", paste(unique(design_subset$group), collapse = ", "))
    message("Number of unique samples per group:")
    for (g in unique(design_subset$group)) {
        n_samples <- sum(design_subset$group == g)
        message("  - ", g, ": ", n_samples, " samples")
    }

    # --- CV Calculation --- #
    # First pivot the data to long format
    long_data <- quant_data |>
        tidyr::pivot_longer(
            cols = dplyr::all_of(sample_names),
            names_to = "Run",
            values_to = "Intensity"
        ) |>
        # Ensure Run column is character for joining
        dplyr::mutate(Run = as.character(.data$Run))

    # Join with design info
    long_data_with_groups <- dplyr::left_join(long_data, design_subset, by = "Run")

    # Check if joining worked correctly
    if (sum(is.na(long_data_with_groups$group)) > 0) {
        message("Warning: Some samples couldn't be matched to groups. Check sample IDs in design matrix.")
        message("Unmatched samples: ", paste(unique(long_data_with_groups$Run[is.na(long_data_with_groups$group)]), collapse = ", "))
    }

    # Remove rows where joining failed or Intensity is NA or 0
    filtered_data <- stats::na.omit(long_data_with_groups) |>
        dplyr::filter(.data$Intensity > 0) # Filter out zero values which can inflate CV

    # Group by metabolite and experimental group
    cv_data <- filtered_data |>
        dplyr::group_by(dplyr::across(dplyr::all_of(c(metabolite_id_col, "group")))) |>
        # Calculate SD and Mean, requiring at least 2 data points for SD
        dplyr::summarise(
            n_samples = dplyr::n(),
            mean_intensity = mean(.data$Intensity, na.rm = TRUE),
            sd_intensity = if (dplyr::n() >= 2) stats::sd(.data$Intensity, na.rm = TRUE) else NA_real_,
            .groups = "drop"
        ) |>
        # Calculate CV, handle mean close to zero or NA sd
        dplyr::mutate(
            cv = dplyr::case_when(
                n_samples < 2 ~ NA_real_, # CV is NA if fewer than 2 samples in group
                is.na(.data$sd_intensity) ~ NA_real_, # CV is NA if sd is NA
                abs(.data$mean_intensity) < .Machine$double.eps ~ NA_real_, # Avoid division by tiny number
                TRUE ~ (.data$sd_intensity / .data$mean_intensity) * 100
            ),
            # Cap CV at a reasonable maximum to avoid extreme outliers
            cv = pmin(cv, 200) # Cap at 200% to avoid extreme outliers affecting visualization
        ) |>
        # Keep only relevant columns and rename metabolite ID column back
        # Use !!rlang::sym() for string column names - {{ }} is for defused symbols only
        dplyr::select(dplyr::all_of(c(metabolite_id_col, "group", "cv", "n_samples"))) |>
        dplyr::rename(metabolite_id = !!rlang::sym(metabolite_id_col))

    # Summary statistics for debugging
    message("CV calculation complete")
    message("CV summary statistics per group:")
    for (g in unique(cv_data$group)) {
        group_cvs <- cv_data$cv[cv_data$group == g & !is.na(cv_data$cv)]
        if (length(group_cvs) > 0) {
            group_stats <- summary(group_cvs)
            message("  - Group ", g, ":")
            message(
                "    Min: ", round(group_stats[1], 1), "%, 1st Qu: ", round(group_stats[2], 1),
                "%, Median: ", round(group_stats[3], 1), "%, Mean: ", round(group_stats[4], 1),
                "%, 3rd Qu: ", round(group_stats[5], 1), "%, Max: ", round(group_stats[6], 1), "%"
            )
            message("    Number of metabolites with CV: ", length(group_cvs))
        } else {
            message("  - Group ", g, ": No valid CV values")
        }
    }

    return(cv_data)
}

# ----------------------------------------------------------------------------
# getInternalStandardMetrics
# ----------------------------------------------------------------------------
#' @title Calculate Internal Standard Metrics
#' @description Identifies internal standards (IS) based on a regex pattern and
#'              calculates their mean intensity and CV across all samples.
#'
#' @param assay_data A tibble/data.frame for one assay.
#' @param is_pattern Character string, regex pattern to identify IS in the metabolite_id_col.
#'                  If NULL, NA, or empty, the function returns an empty data frame.
#' @param metabolite_id_col Character string, the name of the metabolite ID column.
#' @param sample_id_col ***Currently unused.***
#' @param sample_columns Character vector of explicit sample column names from design matrix.
#'
#' @return A data frame with columns 'is_id', 'mean_intensity', 'cv'.
#'         Returns an empty data frame if no IS are found or pattern is invalid.
#'
#' @importFrom tidyr pivot_longer
#' @importFrom dplyr filter group_by summarise mutate select rename all_of
#' @importFrom stringr str_detect
#' @importFrom stats sd na.omit
#' @keywords internal
#' @noRd
#' @export
getInternalStandardMetrics <- function(assay_data,
                                       is_pattern,
                                       metabolite_id_col,
                                       sample_id_col,
                                       sample_columns = NULL) {
    # --- Input Validation --- #
    if (is.null(is_pattern) || is.na(is_pattern) || is_pattern == "") {
        # message("Internal standard pattern is missing or invalid. Skipping IS metrics.")
        return(data.frame(is_id = character(), mean_intensity = numeric(), cv = numeric()))
    }
    if (!metabolite_id_col %in% colnames(assay_data)) {
        warning("Metabolite ID column '", metabolite_id_col, "' not found in assay data for IS metrics.")
        return(data.frame(is_id = character(), mean_intensity = numeric(), cv = numeric()))
    }

    # --- Data Preparation --- #
    quant_info <- getMetaboliteQuantData(assay_data, sample_columns = sample_columns)
    quant_data <- quant_info$quant_data
    sample_names <- quant_info$sample_names
    annotation_data <- quant_info$annotation_data

    if (length(sample_names) == 0 || nrow(quant_data) == 0) {
        return(data.frame(is_id = character(), mean_intensity = numeric(), cv = numeric()))
    }

    # Ensure metabolite_id_col is in annotation_data (may be classified as numeric if integer IDs)
    if (!metabolite_id_col %in% colnames(annotation_data)) {
        if (metabolite_id_col %in% colnames(assay_data)) {
            annotation_data[[metabolite_id_col]] <- as.character(assay_data[[metabolite_id_col]])
        } else {
            warning("Metabolite ID column '", metabolite_id_col, "' not found in annotation data for IS metrics.")
            return(data.frame(is_id = character(), mean_intensity = numeric(), cv = numeric()))
        }
    }

    # Ensure metabolite IDs are present in quant_data for downstream filtering
    quant_data[[metabolite_id_col]] <- annotation_data[[metabolite_id_col]]

    # --- Identify and Filter IS --- #
    # Convert ID column to character for reliable regex matching
    id_values <- as.character(annotation_data[[metabolite_id_col]])

    # Try user-provided pattern first, then fallback to common ITSD patterns
    is_rows <- character(0)

    if (!is.null(is_pattern) && nzchar(is_pattern)) {
        # User-provided pattern
        tryCatch(
            {
                matches <- stringr::str_detect(id_values, is_pattern)
                is_rows <- id_values[matches]
            },
            error = function(e) {
                warning("Invalid regex pattern '", is_pattern, "': ", e$message)
            }
        )
    }

    # Fallback: try common ITSD naming conventions if no matches found
    if (length(is_rows) == 0) {
        common_is_patterns <- c(
            "^IS[_-]" # IS_ or IS- prefix (e.g., IS_Caffeine, IS-Leucine)
            , "^ITSD[_-]" # ITSD_ prefix (e.g., ITSD_A, ITSD_001)
            , "ISTD" # ISTD anywhere (e.g., Caffeine-ISTD, ISTD_001)
            , "ITSD" # ITSD anywhere (alternate spelling)
            , "[_-]d\\d+$" # Deuterated suffix (e.g., Caffeine_d3, Leucine-d10)
            , "[_-]d\\d+[_-]" # Deuterated mid-name (e.g., Caffeine_d3_IS)
            , "\\(d\\d+\\)" # Deuterated in parens (e.g., Caffeine(d3))
            , "(?i)internal.?standard" # "Internal Standard" / "Internal_Standard"
            , "(?i)^is\\d+$" # IS followed by numbers (e.g., IS1, IS23)
            , "13C[_-]?labeled" # 13C labeled (e.g., 13C-labeled_reference)
            , "15N[_-]?labeled" # 15N labeled
            , "-13C\\d*-" # 13C in middle (e.g., Glucose-13C6-IS)
            , "-15N\\d*-" # 15N in middle
        )
        combined_pattern <- paste(common_is_patterns, collapse = "|")

        tryCatch(
            {
                matches <- stringr::str_detect(id_values, combined_pattern)
                is_rows <- id_values[matches]
                if (length(is_rows) > 0) {
                    message("No IS found with user pattern. Using fallback patterns, found ", length(is_rows), " candidates.")
                } else {
                    # Debug: show what we're actually searching
                    sample_ids <- head(id_values, 5)
                    message(
                        "ITSD detection: No matches found in column '", metabolite_id_col,
                        "'. Sample values: ", paste(sample_ids, collapse = ", ")
                    )
                }
            },
            error = function(e) {
                warning("Fallback IS pattern matching failed: ", e$message)
            }
        )
    }

    is_rows <- annotation_data[[metabolite_id_col]][annotation_data[[metabolite_id_col]] %in% is_rows]

    if (length(is_rows) == 0) {
        # message("No internal standards found matching pattern: ", is_pattern)
        return(data.frame(is_id = character(), mean_intensity = numeric(), cv = numeric()))
    }

    is_data <- quant_data |>
        dplyr::filter(.data[[metabolite_id_col]] %in% is_rows)

    # --- Calculate Metrics --- #
    is_metrics <- is_data |>
        tidyr::pivot_longer(
            cols = dplyr::all_of(sample_names),
            names_to = "Run",
            values_to = "Intensity"
        ) |>
        stats::na.omit() |> # Remove missing values before calculating mean/sd
        dplyr::group_by(dplyr::across(dplyr::all_of(metabolite_id_col))) |>
        dplyr::summarise(
            mean_intensity = mean(.data$Intensity, na.rm = TRUE),
            sd_intensity = if (dplyr::n() >= 2) stats::sd(.data$Intensity, na.rm = TRUE) else NA_real_,
            n_samples = dplyr::n(), # Keep track of how many samples contributed
            .groups = "drop"
        ) |>
        dplyr::mutate(
            cv = dplyr::case_when(
                is.na(.data$sd_intensity) ~ NA_real_,
                abs(.data$mean_intensity) < .Machine$double.eps ~ NA_real_,
                TRUE ~ (.data$sd_intensity / .data$mean_intensity) * 100
            )
        ) |>
        dplyr::select(dplyr::all_of(c(metabolite_id_col, "mean_intensity", "cv"))) |>
        dplyr::rename(is_id = dplyr::all_of(metabolite_id_col))

    return(is_metrics)
}

