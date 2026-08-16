#' @title Plot Density for list of ggplot objects
#' @name plotDensity,list-method
#' @importFrom purrr map set_names
#' @importFrom tibble as_tibble
#' @importFrom ggplot2 ggplot aes geom_boxplot theme_bw labs theme element_blank element_text margin
#' @importFrom patchwork plot_layout plot_annotation
#' @importFrom rlang sym !!
#' @export
#' @param theObject,grouping_variable,title,font_size Runtime inputs used by this function; see the usage section for accepted values.
setMethod(
    f = "plotDensity",
    signature = c(theObject = "list"), # Explicitly define signature argument
    definition = function(theObject, grouping_variable, title = "", font_size = 8) {
        # --- Input is a list: Assume list of ggplot objects (Use existing PCA data) ---

        # Basic validation
        if (!all(sapply(theObject, function(x) inherits(x, "ggplot")))) {
            stop("If 'theObject' is a list, all its elements must be ggplot objects.")
        }
        if (!is.character(grouping_variable) || length(grouping_variable) != 1 || is.na(grouping_variable)) {
            stop("`grouping_variable` must be a single non-NA character string.")
        }

        pca_plots_list <- theObject # Rename for clarity

        if (length(pca_plots_list) == 0) {
            warning("Input list of ggplot objects is empty. Returning empty list.")
            return(list())
        }

        # Ensure list is named, provide default names if not
        if (is.null(names(pca_plots_list))) {
            names(pca_plots_list) <- paste0("Plot_", seq_along(pca_plots_list))
            warning("Input ggplot list was unnamed. Using default names (Plot_1, Plot_2, ...).")
        }

        # --- Plotting Logic per Input ggplot ---
        density_plots_list <- purrr::map(seq_along(pca_plots_list), function(i) {
            pca_plot <- pca_plots_list[[i]]
            plot_name <- names(pca_plots_list)[i]
            # Use title override if provided, otherwise use original plot title or name
            plot_title_final <- if (!is.null(title) && title != "") paste(title, "-", plot_name) else tryCatch(pca_plot$labels$title, error = function(e) plot_name)

            # --- Extract PCA data from the ggplot object ---
            pca_data_for_plot <- NULL
            if (!is.null(pca_plot$data) && is.data.frame(pca_plot$data) && all(c("PC1", "PC2", grouping_variable) %in% colnames(pca_plot$data))) {
                pca_data_for_plot <- tibble::as_tibble(pca_plot$data)
            } else {
                warning(sprintf("Plot '%s': Could not reliably extract required data (PC1, PC2, %s) from the ggplot object's internal structure. Skipping density plot generation.", plot_name, grouping_variable))
                return(NULL)
            }

            if (!grouping_variable %in% colnames(pca_data_for_plot)) {
                warning(sprintf("Plot '%s': Grouping variable '%s' not found in extracted data. Skipping.", plot_name, grouping_variable))
                return(NULL)
            }

            # --- Create Density/Box Plots (using extracted data) ---
            tryCatch(
                {
                    # Create PC1 boxplot
                    pc1_box <- ggplot(pca_data_for_plot, aes(x = !!rlang::sym(grouping_variable), y = PC1, fill = !!rlang::sym(grouping_variable))) +
                        geom_boxplot(notch = FALSE) +
                        theme_bw() +
                        labs(
                            title = plot_title_final, # Use final title
                            x = "",
                            y = "PC1"
                        ) +
                        theme(
                            legend.position = "none",
                            axis.text.x = element_blank(),
                            axis.ticks.x = element_blank(),
                            text = element_text(size = font_size),
                            plot.margin = margin(b = 0, t = 5, l = 5, r = 5),
                            panel.grid.major = element_blank(),
                            panel.grid.minor = element_blank(),
                            panel.background = element_blank()
                        )

                    # Create PC2 boxplot
                    pc2_box <- ggplot(pca_data_for_plot, aes(x = !!rlang::sym(grouping_variable), y = PC2, fill = !!rlang::sym(grouping_variable))) +
                        geom_boxplot(notch = FALSE) +
                        theme_bw() +
                        labs(
                            x = "",
                            y = "PC2"
                        ) +
                        theme(
                            legend.position = "none",
                            axis.text.x = element_blank(),
                            axis.ticks.x = element_blank(),
                            text = element_text(size = font_size),
                            plot.margin = margin(t = 0, b = 5, l = 5, r = 5),
                            panel.grid.major = element_blank(),
                            panel.grid.minor = element_blank(),
                            panel.background = element_blank()
                        )

                    # Combine plots with minimal spacing
                    combined_plot <- pc1_box / pc2_box +
                        patchwork::plot_layout(heights = c(1, 1)) +
                        patchwork::plot_annotation(theme = theme(plot.margin = margin(0, 0, 0, 0)))

                    return(combined_plot)
                },
                error = function(e) {
                    warning(sprintf("Plot '%s': Error creating density boxplots from ggplot input: %s. Skipping.", plot_name, e$message))
                    return(NULL)
                }
            )
        })

        # Set names for the list of plots
        names(density_plots_list) <- names(pca_plots_list)

        # Remove NULL elements (skipped plots)
        density_plots_list <- density_plots_list[!sapply(density_plots_list, is.null)]

        return(density_plots_list)
    }
)
