# ----------------------------------------------------------------------------
# calculateLipidCVs
# ----------------------------------------------------------------------------
#' @title Calculate Within-Group Lipid CVs
#' @description Calculates the Coefficient of Variation (CV) for each lipid
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
#' @param lipid_id_col Character string, the name of the lipid ID column in assay_data.
#' @param sample_columns Character vector of explicit sample column names from design matrix.
#'
#' @return A data frame with columns 'lipid_id', 'group', and 'cv'.
#'         Returns an empty data frame if inputs are invalid or calculations fail.
#'
#' @importFrom tidyr pivot_longer
#' @importFrom dplyr left_join group_by summarise filter select rename all_of
#' @importFrom rlang sym
#' @importFrom stats sd na.omit
#' @keywords internal
#' @noRd
#' @export
calculateLipidCVs <- function(assay_data,
                              design_matrix,
                              group_id_col,
                              replicate_id_col,
                              sample_id_col,
                              lipid_id_col,
                              sample_columns = NULL) {
    # --- Input Validation --- #
    required_design_cols <- c(sample_id_col, group_id_col)
    if (!all(required_design_cols %in% colnames(design_matrix))) {
        warning("Design matrix missing required columns: ", paste(setdiff(required_design_cols, colnames(design_matrix)), collapse = ", "))
        return(data.frame(lipid_id = character(), group = character(), cv = numeric()))
    }
    if (!lipid_id_col %in% colnames(assay_data)) {
        warning("Lipid ID column '", lipid_id_col, "' not found in assay data.")
        return(data.frame(lipid_id = character(), group = character(), cv = numeric()))
    }

    # --- Data Preparation --- #
    quant_info <- getLipidQuantData(assay_data, sample_columns = sample_columns)
    quant_data <- quant_info$quant_data
    sample_names <- quant_info$sample_names
    annotation_data <- quant_info$annotation_data

    if (length(sample_names) == 0 || nrow(quant_data) == 0) {
        warning("No sample columns or no data rows found for CV calculation")
        return(data.frame(lipid_id = character(), group = character(), cv = numeric()))
    }

    # Ensure lipid IDs are present for joining
    quant_data[[lipid_id_col]] <- annotation_data[[lipid_id_col]]

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

    # Group by lipid and experimental group
    cv_data <- filtered_data |>
        dplyr::group_by(dplyr::across(dplyr::all_of(c(lipid_id_col, "group")))) |>
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
        # Keep only relevant columns and rename lipid ID column back
        # Use !!rlang::sym() for string column names - {{ }} is for defused symbols only
        dplyr::select(dplyr::all_of(c(lipid_id_col, "group", "cv", "n_samples"))) |>
        dplyr::rename(lipid_id = !!rlang::sym(lipid_id_col))

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
            message("    Number of lipids with CV: ", length(group_cvs))
        } else {
            message("  - Group ", g, ": No valid CV values")
        }
    }

    return(cv_data)
}

# ----------------------------------------------------------------------------
# getLipidInternalStandardMetrics
# ----------------------------------------------------------------------------
#' @title Calculate Internal Standard Metrics
#' @description Identifies internal standards (IS) based on a regex pattern and
#'              calculates their mean intensity and CV across all samples.
#'
#' @param assay_data A tibble/data.frame for one assay.
#' @param is_pattern Character string, regex pattern to identify IS in the lipid_id_col.
#'                  If NULL, NA, or empty, the function returns an empty data frame.
#' @param lipid_id_col Character string, the name of the lipid ID column.
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
getLipidInternalStandardMetrics <- function(assay_data,
                                       is_pattern,
                                       lipid_id_col,
                                       sample_id_col,
                                       sample_columns = NULL) {
    # --- Input Validation --- #
    if (is.null(is_pattern) || is.na(is_pattern) || is_pattern == "") {
        # message("Internal standard pattern is missing or invalid. Skipping IS metrics.")
        return(data.frame(is_id = character(), mean_intensity = numeric(), cv = numeric()))
    }
    if (!lipid_id_col %in% colnames(assay_data)) {
        warning("Lipid ID column '", lipid_id_col, "' not found in assay data for IS metrics.")
        return(data.frame(is_id = character(), mean_intensity = numeric(), cv = numeric()))
    }

    # --- Data Preparation --- #
    quant_info <- getLipidQuantData(assay_data, sample_columns = sample_columns)
    quant_data <- quant_info$quant_data
    sample_names <- quant_info$sample_names
    annotation_data <- quant_info$annotation_data

    if (length(sample_names) == 0 || nrow(quant_data) == 0) {
        return(data.frame(is_id = character(), mean_intensity = numeric(), cv = numeric()))
    }

    # Ensure lipid_id_col is in annotation_data (may be classified as numeric if integer IDs)
    if (!lipid_id_col %in% colnames(annotation_data)) {
        if (lipid_id_col %in% colnames(assay_data)) {
            annotation_data[[lipid_id_col]] <- as.character(assay_data[[lipid_id_col]])
        } else {
            warning("Lipid ID column '", lipid_id_col, "' not found in annotation data for IS metrics.")
            return(data.frame(is_id = character(), mean_intensity = numeric(), cv = numeric()))
        }
    }

    # Ensure lipid IDs are present in quant_data for downstream filtering
    quant_data[[lipid_id_col]] <- annotation_data[[lipid_id_col]]

    # --- Identify and Filter IS --- #
    # Convert ID column to character for reliable regex matching
    id_values <- as.character(annotation_data[[lipid_id_col]])

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
                        "ITSD detection: No matches found in column '", lipid_id_col,
                        "'. Sample values: ", paste(sample_ids, collapse = ", ")
                    )
                }
            },
            error = function(e) {
                warning("Fallback IS pattern matching failed: ", e$message)
            }
        )
    }

    is_rows <- annotation_data[[lipid_id_col]][annotation_data[[lipid_id_col]] %in% is_rows]

    if (length(is_rows) == 0) {
        # message("No internal standards found matching pattern: ", is_pattern)
        return(data.frame(is_id = character(), mean_intensity = numeric(), cv = numeric()))
    }

    is_data <- quant_data |>
        dplyr::filter(.data[[lipid_id_col]] %in% is_rows)

    # --- Calculate Metrics --- #
    is_metrics <- is_data |>
        tidyr::pivot_longer(
            cols = dplyr::all_of(sample_names),
            names_to = "Run",
            values_to = "Intensity"
        ) |>
        stats::na.omit() |> # Remove missing values before calculating mean/sd
        dplyr::group_by(dplyr::across(dplyr::all_of(lipid_id_col))) |>
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
        dplyr::select(dplyr::all_of(c(lipid_id_col, "mean_intensity", "cv"))) |>
        dplyr::rename(is_id = dplyr::all_of(lipid_id_col))

    return(is_metrics)
}

# ----------------------------------------------------------------------------
# generateLipidFilteringPlots
# ----------------------------------------------------------------------------
#' @title Generate Lipidomics Filtering Progress Plots
#' @description Creates a set of quality control plots based on the metrics
#'              stored in the global `FilteringProgressLipidomics` object.
#'
#' @details
#' Generates visualizations for key lipidomics QC metrics across processing steps:
#' \itemize{
#'   \item Total unique lipids across assays per step (bar chart)
#'   \item Lipids per assay per step (bar chart)
#'   \item Detected lipids per sample per step (line chart)
#'   \item Missing value percentage per assay per step (bar chart)
#'   \item Total intensity per sample per step (line chart)
#'   \item Coefficient of variation distribution per step (box plot)
#'   \item Internal standard metrics (CV and intensity) per step (if available)
#' }
#'
#' @param prog_met A `FilteringProgressLipidomics` object containing the tracked
#'                 metrics across various processing steps. If not provided, the
#'                 function retrieves the global object.
#'
#' @return A list of ggplot objects, one for each visualization type.
#'
#' @importFrom ggplot2 ggplot aes geom_bar geom_text labs theme_minimal theme element_text element_blank geom_line geom_point scale_color_manual annotate theme_void geom_boxplot coord_cartesian facet_wrap geom_violin scale_fill_brewer position_dodge
#' @importFrom dplyr bind_rows mutate group_by ungroup
#' @importFrom forcats fct_reorder
#' @importFrom tidyr pivot_longer
#' @keywords internal
#' @noRd
#' @export
generateLipidFilteringPlots <- function(prog_met = NULL) {
    # Get the global object if not provided
    if (is.null(prog_met)) {
        prog_met <- getFilteringProgressLipidomics()
    }


    # Return empty list if no steps have been tracked
    if (length(prog_met@steps) == 0) {
        message("No lipidomics filtering steps have been tracked yet.")
        return(list())
    }

    plot_list <- list()

    # --- 1. Total Lipids Plot (All Assays Combined) --- #
    plot_list$total_lipids <- tryCatch(
        {
            ggplot(data.frame(
                step = factor(prog_met@steps, levels = prog_met@steps),
                total_lipids = prog_met@n_lipids_total
            ), aes(x = step, y = total_lipids)) +
                geom_bar(stat = "identity", fill = "steelblue", width = 0.7) +
                geom_text(aes(label = total_lipids),
                    vjust = -0.5,
                    size = 4
                ) +
                labs(
                    title = "Total Unique Lipids (All Assays)",
                    x = "Filtering Step",
                    y = "Unique Lipids"
                ) +
                theme_minimal() +
                theme(
                    axis.text.x = element_text(angle = 45, hjust = 1),
                    panel.grid.major.x = element_blank()
                )
        },
        error = function(e) {
            NULL
        }
    )

    # --- 2. Lipids Per Assay Plot --- #
    # First create a data frame with lipids per assay per step - use purrr
    lipids_per_assay_df <- purrr::map_dfr(seq_along(prog_met@steps), function(step_idx) {
        step <- prog_met@steps[step_idx]
        assay_names <- prog_met@assay_names[[step_idx]]
        lipid_counts <- unlist(prog_met@n_lipids_per_assay[[step_idx]])

        if (length(lipid_counts) > 0) {
            return(data.frame(
                step = step,
                assay = names(lipid_counts),
                n_lipids = as.numeric(lipid_counts)
            ))
        }
        return(NULL) # Return NULL if no lipid counts (will be filtered by map_dfr)
    })

    if (nrow(lipids_per_assay_df) > 0) {
        plot_list$lipids_per_assay <- ggplot(
            lipids_per_assay_df,
            aes(
                x = factor(step, levels = prog_met@steps),
                y = n_lipids,
                fill = assay
            )
        ) +
            geom_bar(stat = "identity", position = "dodge", width = 0.7) +
            geom_text(aes(label = n_lipids),
                position = position_dodge(width = 0.7),
                vjust = -0.5,
                size = 3
            ) +
            labs(
                title = "Lipids per Assay",
                x = "Filtering Step",
                y = "Unique Lipids",
                fill = "Assay"
            ) +
            theme_minimal() +
            theme(
                axis.text.x = element_text(angle = 45, hjust = 1),
                panel.grid.major.x = element_blank()
            ) +
            scale_fill_brewer(palette = "Set1")
    }


    # --- 3. Detected Lipids Per Sample Plot --- #
    # Create a data frame with detected lipids per sample per step - use purrr
    detected_per_sample_df <- purrr::map_dfr(seq_along(prog_met@steps), function(step_idx) {
        step <- prog_met@steps[step_idx]
        assay_names <- prog_met@assay_names[[step_idx]]

        # Use purrr::map_dfr to combine results from all assays
        purrr::map_dfr(assay_names, function(assay_name) {
            if (assay_name %in% names(prog_met@detected_per_sample[[step_idx]])) {
                detected_df <- prog_met@detected_per_sample[[step_idx]][[assay_name]]
                if (nrow(detected_df) > 0) {
                    detected_df$step <- step
                    detected_df$assay <- assay_name
                    return(detected_df)
                }
            }
            return(NULL) # Return NULL if no data (will be filtered by map_dfr)
        })
    })

    if (nrow(detected_per_sample_df) > 0) {
        # Ensure consistent data types
        detected_per_sample_df$Run <- as.character(detected_per_sample_df$Run)
        detected_per_sample_df$n_detected <- as.numeric(detected_per_sample_df$n_detected)
        detected_per_sample_df$step <- factor(detected_per_sample_df$step,
            levels = prog_met@steps
        )

        plot_list$detected_per_sample <- detected_per_sample_df |>
            group_by(Run, assay) |>
            mutate(avg_detected = mean(n_detected)) |>
            ungroup() |>
            mutate(Run = fct_reorder(Run, avg_detected)) |>
            ggplot(aes(
                x = Run, y = n_detected,
                group = interaction(step, assay),
                color = step,
                linetype = assay
            )) +
            geom_line() +
            geom_point() +
            labs(
                title = "Detected Lipids per Sample",
                x = "Sample ID (ordered by average detected lipids)",
                y = "Number of Detected Lipids",
                color = "Step",
                linetype = "Assay"
            ) +
            theme_minimal() +
            theme(
                axis.text.x = element_text(angle = 45, hjust = 1),
                panel.grid.major.x = element_blank()
            ) +
            scale_color_brewer(palette = "Set2")
    }


    # --- 4. Missingness Percentage Per Assay Plot --- #
    # Create a data frame with missingness percentage per assay per step - use purrr
    missingness_df <- purrr::map_dfr(seq_along(prog_met@steps), function(step_idx) {
        step <- prog_met@steps[step_idx]
        assay_names <- prog_met@assay_names[[step_idx]]

        # Use purrr::map_dfr to combine results from all assays
        purrr::map_dfr(assay_names, function(assay_name) {
            if (assay_name %in% names(prog_met@missingness_per_assay[[step_idx]])) {
                missingness <- prog_met@missingness_per_assay[[step_idx]][[assay_name]]
                if (!is.null(missingness) && !is.na(missingness)) {
                    return(data.frame(
                        step = step,
                        assay = assay_name,
                        missingness = as.numeric(missingness)
                    ))
                }
            }
            return(NULL) # Return NULL if no data (will be filtered by map_dfr)
        })
    })

    if (nrow(missingness_df) > 0) {
        plot_list$missingness <- ggplot(
            missingness_df,
            aes(
                x = factor(step, levels = prog_met@steps),
                y = missingness,
                fill = assay
            )
        ) +
            geom_bar(stat = "identity", position = "dodge", width = 0.7) +
            geom_text(aes(label = sprintf("%.1f%%", missingness)),
                position = position_dodge(width = 0.7),
                vjust = -0.5,
                size = 3
            ) +
            labs(
                title = "Missing Values Percentage per Assay",
                x = "Filtering Step",
                y = "Missing Values (%)",
                fill = "Assay"
            ) +
            theme_minimal() +
            theme(
                axis.text.x = element_text(angle = 45, hjust = 1),
                panel.grid.major.x = element_blank()
            ) +
            scale_fill_brewer(palette = "Set1")
    }


    # --- 5. Sum Intensity Per Sample Plot --- #
    # Create a data frame with sum intensity per sample per step - use purrr
    sum_intensity_df <- purrr::map_dfr(seq_along(prog_met@steps), function(step_idx) {
        step <- prog_met@steps[step_idx]
        assay_names <- prog_met@assay_names[[step_idx]]

        # Use purrr::map_dfr to combine results from all assays
        purrr::map_dfr(assay_names, function(assay_name) {
            if (assay_name %in% names(prog_met@sum_intensity_per_sample[[step_idx]])) {
                intensity_df <- prog_met@sum_intensity_per_sample[[step_idx]][[assay_name]]
                if (nrow(intensity_df) > 0) {
                    intensity_df$step <- step
                    intensity_df$assay <- assay_name
                    return(intensity_df)
                }
            }
            return(NULL) # Return NULL if no data (will be filtered by map_dfr)
        })
    })

    if (nrow(sum_intensity_df) > 0) {
        # Ensure consistent data types
        sum_intensity_df$Run <- as.character(sum_intensity_df$Run)
        sum_intensity_df$sum_intensity <- as.numeric(sum_intensity_df$sum_intensity)
        sum_intensity_df$step <- factor(sum_intensity_df$step,
            levels = prog_met@steps
        )

        plot_list$sum_intensity <- sum_intensity_df |>
            group_by(Run, assay) |>
            mutate(avg_intensity = mean(sum_intensity)) |>
            ungroup() |>
            mutate(
                Run = fct_reorder(Run, avg_intensity),
                # Apply log2 transformation directly to the values
                log2_intensity = log2(pmax(sum_intensity, 1)) # Avoid log(0) with pmax
            ) |>
            ggplot(aes(
                x = Run, y = log2_intensity, # Plot the transformed values
                group = interaction(step, assay),
                color = step,
                linetype = assay
            )) +
            geom_line() +
            geom_point() +
            labs(
                title = "Total Signal Intensity per Sample",
                x = "Sample ID (ordered by average intensity)",
                y = "Sum Intensity (log2 scale)",
                color = "Step",
                linetype = "Assay"
            ) +
            theme_minimal() +
            theme(
                axis.text.x = element_text(angle = 45, hjust = 1),
                panel.grid.major.x = element_blank()
            ) +
            scale_color_brewer(palette = "Set2")
    }


    # --- 6. CV Distribution Plot --- #
    # Create a data frame with CV distribution per step - use purrr
    cv_df <- purrr::map_dfr(seq_along(prog_met@steps), function(step_idx) {
        step <- prog_met@steps[step_idx]
        assay_names <- prog_met@assay_names[[step_idx]]

        # Use purrr::map_dfr to combine results from all assays
        purrr::map_dfr(assay_names, function(assay_name) {
            if (assay_name %in% names(prog_met@cv_distribution_per_assay[[step_idx]])) {
                step_cv_df <- prog_met@cv_distribution_per_assay[[step_idx]][[assay_name]]
                if (nrow(step_cv_df) > 0) {
                    step_cv_df$step <- step
                    step_cv_df$assay <- assay_name
                    return(step_cv_df)
                }
            }
            return(NULL) # Return NULL if no data (will be filtered by map_dfr)
        })
    })

    if (nrow(cv_df) > 0) {
        # Ensure consistent data types and remove outliers
        cv_df$cv <- as.numeric(cv_df$cv)
        cv_df$step <- factor(cv_df$step, levels = prog_met@steps)

        # Calculate 95th percentile for y-axis limit
        q95 <- quantile(cv_df$cv, 0.95, na.rm = TRUE)

        plot_list$cv_distribution <- ggplot(
            cv_df,
            aes(
                x = step,
                y = cv,
                fill = assay
            )
        ) +
            geom_boxplot(alpha = 0.7, outlier.shape = NA) +
            labs(
                title = "CV Distribution by Group",
                x = "Filtering Step",
                y = "Coefficient of Variation (%)",
                fill = "Assay"
            ) +
            theme_minimal() +
            theme(
                axis.text.x = element_text(angle = 45, hjust = 1),
                panel.grid.major.x = element_blank()
            ) +
            coord_cartesian(ylim = c(0, q95)) +
            scale_fill_brewer(palette = "Set1") +
            facet_wrap(~group, scales = "free_y")
    }

    # --- 7. Internal Standards Metrics Plot --- #
    # Create a data frame with internal standard metrics per step - use purrr
    is_df <- purrr::map_dfr(seq_along(prog_met@steps), function(step_idx) {
        step <- prog_met@steps[step_idx]
        assay_names <- prog_met@assay_names[[step_idx]]

        # Use purrr::map_dfr to combine results from all assays
        purrr::map_dfr(assay_names, function(assay_name) {
            if (assay_name %in% names(prog_met@is_metrics_per_assay[[step_idx]])) {
                step_is_df <- prog_met@is_metrics_per_assay[[step_idx]][[assay_name]]
                if (nrow(step_is_df) > 0) {
                    step_is_df$step <- step
                    step_is_df$assay <- assay_name
                    return(step_is_df)
                }
            }
            return(NULL) # Return NULL if no data (will be filtered by map_dfr)
        })
    })

    if (nrow(is_df) > 0) {
        # Ensure consistent data types
        is_df$cv <- as.numeric(is_df$cv)
        is_df$mean_intensity <- as.numeric(is_df$mean_intensity)
        is_df$step <- factor(is_df$step, levels = prog_met@steps)

        # Create CV plot for internal standards
        plot_list$is_cv <- ggplot(
            is_df,
            aes(
                x = step,
                y = cv,
                fill = assay
            )
        ) +
            geom_violin(alpha = 0.7, trim = FALSE, scale = "width") +
            geom_boxplot(width = 0.1, fill = "white", alpha = 0.7) +
            labs(
                title = "Internal Standards CV",
                x = "Filtering Step",
                y = "Coefficient of Variation (%)",
                fill = "Assay"
            ) +
            theme_minimal() +
            theme(
                axis.text.x = element_text(angle = 45, hjust = 1),
                panel.grid.major.x = element_blank()
            ) +
            scale_fill_brewer(palette = "Set1")

        # Create mean intensity plot for internal standards
        plot_list$is_intensity <- is_df |>
            # Apply log2 transformation directly to the data
            mutate(log2_intensity = log2(pmax(mean_intensity, 1))) |> # Avoid log(0) with pmax
            ggplot(aes(
                x = step,
                y = log2_intensity, # Plot transformed values
                fill = assay
            )) +
            geom_violin(alpha = 0.7, trim = FALSE, scale = "width") +
            geom_boxplot(width = 0.1, fill = "white", alpha = 0.7) +
            labs(
                title = "Internal Standards Mean Intensity",
                x = "Filtering Step",
                y = "Mean Intensity (log2 scale)",
                fill = "Assay"
            ) +
            theme_minimal() +
            theme(
                axis.text.x = element_text(angle = 45, hjust = 1),
                panel.grid.major.x = element_blank()
            ) +
            scale_fill_brewer(palette = "Set1")
    }


    return(plot_list)
}

