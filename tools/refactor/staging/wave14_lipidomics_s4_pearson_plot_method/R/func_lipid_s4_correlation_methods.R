#' @title Plot Pearson Correlation
#' @name plotPearson,LipidomicsAssayData-method
#' @importFrom purrr map set_names
#' @importFrom ggplot2 ggplot aes geom_histogram scale_y_continuous xlab ylab theme element_blank
#' @export
setMethod(
    f = "plotPearson",
    signature = "LipidomicsAssayData",
    definition = function(theObject, tech_rep_remove_regex = "pool", correlation_group = NA) {
        # Get the list of correlation tibbles (one per assay)
        # tech_rep_remove_regex and correlation_group are passed down
        correlation_list <- pearsonCorForSamplePairs(theObject,
            tech_rep_remove_regex = tech_rep_remove_regex,
            correlation_group = correlation_group
        )

        if (length(correlation_list) == 0) {
            warning("No correlation results generated (likely no valid assays). Returning empty list.")
            return(list())
        }

        # Ensure list is named (pearsonCorForSamplePairs should have handled this, but double-check)
        if (is.null(names(correlation_list))) {
            names(correlation_list) <- paste0("Assay_", seq_along(correlation_list))
        }


        # --- Plotting Logic per Assay's Correlation Results ---
        pearson_plots_list <- purrr::map(seq_along(correlation_list), function(i) {
            assay_name <- names(correlation_list)[i]
            correlation_vec <- correlation_list[[i]]

            # Check if the correlation data is valid
            if (is.null(correlation_vec) || nrow(correlation_vec) == 0 || !"pearson_correlation" %in% colnames(correlation_vec)) {
                warning(sprintf("Assay '%s': Invalid or empty correlation data provided. Skipping Pearson plot.", assay_name))
                return(NULL)
            }

            # Check for all NA values
            if (all(is.na(correlation_vec$pearson_correlation))) {
                warning(sprintf("Assay '%s': All Pearson correlation values are NA. Skipping plot.", assay_name))
                return(NULL)
            }

            # Calculate breaks carefully, handling potential NAs and edge cases
            min_cor <- min(correlation_vec$pearson_correlation, na.rm = TRUE)
            # Ensure min_cor is finite; default if not
            if (!is.finite(min_cor)) min_cor <- 0

            # Use finer breaks, similar to protein version, clamped to [0, 1]
            # Note: Protein version uses 0.001 step, using 0.01 here for potentially better visibility first.
            hist_breaks <- seq(0, 1, 0.01)

            # --- Create Plot ---
            tryCatch(
                {
                    pearson_plot <- correlation_vec |>
                        ggplot(aes(pearson_correlation)) +
                        geom_histogram(breaks = hist_breaks, na.rm = TRUE) +
                        # Set x-axis limits and breaks
                        scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.1), expand = c(0, 0)) +
                        # # Set fixed y-axis scale, similar to protein version
                        # scale_y_continuous(breaks = seq(0, 4, 1), limits = c(0, 4), expand = c(0, 0)) +
                        xlab("Pearson Correlation") +
                        ylab("Counts") +
                        # ggtitle(paste(assay_name)) +
                        theme_bw() +
                        theme(
                            panel.grid.major = element_blank(),
                            panel.grid.minor = element_blank()
                        )

                    return(pearson_plot)
                },
                error = function(e) {
                    warning(sprintf("Assay '%s': Error creating Pearson histogram: %s. Skipping.", assay_name, e$message))
                    return(NULL)
                }
            )
        })

        # Set names for the list of plots
        names(pearson_plots_list) <- names(correlation_list)

        # Remove NULL elements (skipped assays)
        pearson_plots_list <- pearson_plots_list[!sapply(pearson_plots_list, is.null)]

        return(pearson_plots_list)
    }
)

