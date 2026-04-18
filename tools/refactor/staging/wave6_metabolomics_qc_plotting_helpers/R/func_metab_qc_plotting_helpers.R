# ----------------------------------------------------------------------------
# generateMetaboliteFilteringPlots
# ----------------------------------------------------------------------------
#' @title Generate Metabolomics Filtering Progress Plots
#' @description Creates a set of quality control plots based on the metrics
#'              stored in the global `FilteringProgressMetabolomics` object.
#'
#' @details
#' Generates visualizations for key metabolomics QC metrics across processing steps:
#' \itemize{
#'   \item Total unique metabolites across assays per step (bar chart)
#'   \item Metabolites per assay per step (bar chart)
#'   \item Detected metabolites per sample per step (line chart)
#'   \item Missing value percentage per assay per step (bar chart)
#'   \item Total intensity per sample per step (line chart)
#'   \item Coefficient of variation distribution per step (box plot)
#'   \item Internal standard metrics (CV and intensity) per step (if available)
#' }
#'
#' @param prog_met A `FilteringProgressMetabolomics` object containing the tracked
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
generateMetaboliteFilteringPlots <- function(prog_met = NULL) {
    # Get the global object if not provided
    if (is.null(prog_met)) {
        prog_met <- getFilteringProgressMetabolomics()
    }


    # Return empty list if no steps have been tracked
    if (length(prog_met@steps) == 0) {
        message("No metabolomics filtering steps have been tracked yet.")
        return(list())
    }

    plot_list <- list()

    # --- 1. Total Metabolites Plot (All Assays Combined) --- #
    plot_list$total_metabolites <- tryCatch(
        {
            ggplot(data.frame(
                step = factor(prog_met@steps, levels = prog_met@steps),
                total_metabolites = prog_met@n_metabolites_total
            ), aes(x = step, y = total_metabolites)) +
                geom_bar(stat = "identity", fill = "steelblue", width = 0.7) +
                geom_text(aes(label = total_metabolites),
                    vjust = -0.5,
                    size = 4
                ) +
                labs(
                    title = "Total Unique Metabolites (All Assays)",
                    x = "Filtering Step",
                    y = "Unique Metabolites"
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

    # --- 2. Metabolites Per Assay Plot --- #
    # First create a data frame with metabolites per assay per step - use purrr
    metabolites_per_assay_df <- purrr::map_dfr(seq_along(prog_met@steps), function(step_idx) {
        step <- prog_met@steps[step_idx]
        assay_names <- prog_met@assay_names[[step_idx]]
        metabolite_counts <- unlist(prog_met@n_metabolites_per_assay[[step_idx]])

        if (length(metabolite_counts) > 0) {
            return(data.frame(
                step = step,
                assay = names(metabolite_counts),
                n_metabolites = as.numeric(metabolite_counts)
            ))
        }
        return(NULL) # Return NULL if no metabolite counts (will be filtered by map_dfr)
    })

    if (nrow(metabolites_per_assay_df) > 0) {
        plot_list$metabolites_per_assay <- ggplot(
            metabolites_per_assay_df,
            aes(
                x = factor(step, levels = prog_met@steps),
                y = n_metabolites,
                fill = assay
            )
        ) +
            geom_bar(stat = "identity", position = "dodge", width = 0.7) +
            geom_text(aes(label = n_metabolites),
                position = position_dodge(width = 0.7),
                vjust = -0.5,
                size = 3
            ) +
            labs(
                title = "Metabolites per Assay",
                x = "Filtering Step",
                y = "Unique Metabolites",
                fill = "Assay"
            ) +
            theme_minimal() +
            theme(
                axis.text.x = element_text(angle = 45, hjust = 1),
                panel.grid.major.x = element_blank()
            ) +
            scale_fill_brewer(palette = "Set1")
    }


    # --- 3. Detected Metabolites Per Sample Plot --- #
    # Create a data frame with detected metabolites per sample per step - use purrr
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
                title = "Detected Metabolites per Sample",
                x = "Sample ID (ordered by average detected metabolites)",
                y = "Number of Detected Metabolites",
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

