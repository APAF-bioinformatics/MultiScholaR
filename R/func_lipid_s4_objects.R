# ============================================================================
# func_lipid_s4_objects.R
# ============================================================================
# Purpose: Lipidomics S4 class definitions and methods
#
# This file contains S4 class definitions and methods specific to lipidomics,
# including LipidomicsAssayData, LipidomicsDifferentialAbundanceResults classes,
# and all lipidomics-specific S4 methods.
#
# Consolidated from:
# - lipidVsSamplesS4Objects.R (LipidomicsAssayData class + 16 methods)
# - lipid_da_analysis_wrapper.R (LipidomicsDifferentialAbundanceResults class + 2 methods)
# - lipid_normalization.R (logTransformAssays, normaliseUntransformedData)
# - lipid_qc.R (lipidIntensityFilteringHelper)
# - QC_visualisation.R (FilteringProgressLipidomics class)
#
# Dependencies:
# - methods package
# - func_general_s4_generics.R (for generic definitions)
# ============================================================================




## -----------------------------------------------------------------------------
## Plotting Methods for LipidomicsAssayData
## -----------------------------------------------------------------------------












## Lipid pair-correlation helper now lives in
## func_lipid_qc_support_helpers.R.



## -----------------------------------------------------------------------------
## Negative Control Selection Methods for LipidomicsAssayData
## -----------------------------------------------------------------------------








# Get the differential expression results in wide format





# htmlwidgets::saveWidget( widget = theObject@interactive_volcano_plot
# , file = file.path( output_dir
#                     , paste0(colnames(r_obj$coefficients)[coef], ".html"))  #the path & file name
# , selfcontained = TRUE #creates a single html file
# )


## ----------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Create a QC composite figure
# COMMENTED OUT: GridPlotData class does not exist, method not in use
#
# #' @export
# setMethod(f = "createGridQCLipidomics",
#           signature = "GridPlotData",
#           definition = function(theObject, pca_titles = NULL, density_titles = NULL, rle_titles = NULL, pearson_titles = NULL, save_path = NULL, file_name = "pca_density_rle_pearson_corr_plots_merged") {
#
#             # --- Defensive check on input titles ---
#             stopifnot("pca_titles must be a character vector or NULL" = is.null(pca_titles) || is.character(pca_titles))
#             stopifnot("density_titles must be a character vector or NULL" = is.null(density_titles) || is.character(density_titles))
#             stopifnot("rle_titles must be a character vector or NULL" = is.null(rle_titles) || is.character(rle_titles))
#             stopifnot("pearson_titles must be a character vector or NULL" = is.null(pearson_titles) || is.character(pearson_titles))
#
#             # --- Identify all unique assay names from the plot lists ---
#             all_plot_lists <- c(theObject@pca_plots, theObject@density_plots, theObject@rle_plots, theObject@pearson_plots)
#             all_plot_lists <- all_plot_lists[sapply(all_plot_lists, function(x) is.list(x) && length(x) > 0)]
#
#             assay_names <- if (length(all_plot_lists) > 0) unique(unlist(lapply(all_plot_lists, names))) else character(0)
#
#             if (length(assay_names) == 0 || all(sapply(assay_names, is.null))) {
#                 warning("No assays with named plots found. Cannot generate composite QC plot.", immediate. = TRUE)
#                 return(list())
#             }
#
#             # --- Loop over each assay to create a composite plot ---
#             composite_plots_list <- purrr::map(assay_names, function(current_assay_name) {
#                 message(sprintf("--- Generating composite QC plot for assay: %s ---", current_assay_name))
#
#                 # Helper to extract and prepare plots for the current assay
#                 prepare_plot_row <- function(plot_groups_list) {
#                     plots <- purrr::map(plot_groups_list, ~ .x[[current_assay_name]])
#                     # Replace any NULLs with a blank plot to maintain grid alignment
#                     lapply(plots, function(p) if(is.null(p)) ggplot() + theme_void() else p)
#                 }
#
#                 # Extract and prepare plots for the current assay
#                 pca_plots_assay <- prepare_plot_row(theObject@pca_plots)
#                 density_plots_assay <- prepare_plot_row(theObject@density_plots)
#                 rle_plots_assay <- prepare_plot_row(theObject@rle_plots)
#                 pearson_plots_assay <- prepare_plot_row(theObject@pearson_plots)
#
#                 # --- Plot creation helper functions ---
#                 createLabelPlot <- function(title) {
#                   ggplot() +
#                     annotate("text", x = 0, y = 0.5, label = title, size = 5, hjust = 0) +
#                     xlim(0, 1) +
#                     theme_void() +
#                     theme(plot.margin = margin(5, 5, 5, 5), panel.background = element_blank())
#                 }
#                 # These functions now only apply theme, as titles are handled by labels
#                 applyTheme <- function(plot) {
#                     if (inherits(plot, "ggplot") && !is.null(plot$data)) { # Check if it's not an empty plot
#                          plot <- plot + theme(text = element_text(size = 15),
#                                    panel.grid.major = element_blank(),
#                                    panel.grid.minor = element_blank(),
#                                    panel.background = element_blank())
#                          if("patchwork" %in% class(plot)) { # for density plots
#                              plot <- plot & theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
#                          }
#                     }
#                     plot
#                 }
#
#
#                 # --- Generate and combine plots for the current assay ---
#                 plots_to_combine <- list()
#                 add_plot_row <- function(plots, titles, create_fn) {
#                     # Only add row if there are titles provided for it
#                     if (!is.null(titles) && length(titles) > 0 && length(plots) > 0) {
#                         # Number of columns is now determined by the number of plots for THIS row
#                         num_cols_row <- length(plots)
#                         # Ensure the number of titles matches the number of plots
#                         if(length(titles) != num_cols_row){
#                              warning(sprintf("Mismatch between number of titles (%d) and number of plots (%d). Adjusting titles to match plots.", length(titles), num_cols_row), immediate. = TRUE)
#                              # Create a character vector of the correct length, filled with blanks
#                              adjusted_titles <- character(num_cols_row)
#                              # Copy the provided titles into the new vector
#                              n_to_copy <- min(length(titles), num_cols_row)
#                              if (n_to_copy > 0) {
#                                 adjusted_titles[1:n_to_copy] <- titles[1:n_to_copy]
#                              }
#                              titles <- adjusted_titles
#                         }
#                         list(wrap_plots(lapply(titles, createLabelPlot), ncol = num_cols_row),
#                              wrap_plots(lapply(plots, create_fn), ncol = num_cols_row))
#                     } else {
#                         list()
#                     }
#                 }
#
#                 plots_to_combine <- c(plots_to_combine, add_plot_row(pca_plots_assay, pca_titles, applyTheme))
#                 plots_to_combine <- c(plots_to_combine, add_plot_row(density_plots_assay, density_titles, applyTheme))
#                 plots_to_combine <- c(plots_to_combine, add_plot_row(rle_plots_assay, rle_titles, applyTheme))
#                 plots_to_combine <- c(plots_to_combine, add_plot_row(pearson_plots_assay, pearson_titles, applyTheme))
#
#                 if (length(plots_to_combine) == 0) {
#                     warning(paste("No plots to combine for assay:", current_assay_name))
#                     return(NULL)
#                 }
#
#                 num_rows <- length(plots_to_combine) / 2
#                 layout_heights <- rep(c(0.1, 1), num_rows)
#
#                 # Determine overall width by the row with the maximum number of plots
#                 max_cols <- max(
#                     length(pca_plots_assay),
#                     length(density_plots_assay),
#                     length(rle_plots_assay),
#                     length(pearson_plots_assay)
#                 )
#
#                 combined_plot <- wrap_plots(plots_to_combine, ncol = 1) + plot_layout(heights = layout_heights)
#
#                 if (!is.null(save_path)) {
#                   assay_file_name <- paste0(file_name, "_", current_assay_name)
#                   sapply(c("png", "pdf", "svg"), function(ext) {
#                     ggsave(
#                       plot = combined_plot,
#                       filename = file.path(save_path, paste0(assay_file_name, ".", ext)),
#                       width = 3.5 * max_cols,
#                       height = 4 * num_rows
#                     )
#                   })
#                   message(paste("Plots saved for assay '", current_assay_name, "' in", save_path))
#                 }
#
#                 return(combined_plot)
#             })
#
#             names(composite_plots_list) <- assay_names
#             composite_plots_list[!sapply(composite_plots_list, is.null)]
#           })



## ----------------------------------------------------------------------------------------------------------------------------------------------------------------------


## ----------------------------------------------------------------------------------------------------------------------------------------------------------------------


# ==========================================
# Content from lipid_normalization.R
# ==========================================
## Log-transform and normalization methods now live in
## func_lipid_s4_normalization_methods.R.



# Optional: Add other normalization methods (e.g., PQN, Median) here later
# setMethod(f = "normaliseUntransformedData",
#           signature = signature(theObject = "LipidomicsAssayData", method = "character"),
#           definition = function(theObject, method = "PQN", ...) { ... }
# )
# ==========================================
# Content from lipid_qc.R
# ==========================================
## Lipid intensity-filtering helper now lives in
## func_lipid_qc.R.
## The public S4 method shell now lives in
## func_lipid_s4_normalization_methods.R.

#-------------------------------------------------------------------------------

## Duplicate-resolution helpers now live in func_lipid_s4_duplicate_helpers.R.
## Keep the public S4 method shell in this wrapper file.



#' @title Resolve Duplicate Features for LipidomicsAssayData
#' @name resolveDuplicateFeatures,LipidomicsAssayData-method
#' @export
setMethod("resolveDuplicateFeatures",
    signature = "LipidomicsAssayData",
    definition = function(theObject, itsd_pattern_columns = NULL) {
        resolveDuplicateFeaturesForLipidObject(
            theObject = theObject,
            itsd_pattern_columns = itsd_pattern_columns
        )
    }
)


# ==========================================
# FilteringProgressLipidomics class from QC_visualisation.R
# ==========================================

## FilteringProgressLipidomics helpers now live in
## func_lipid_s4_progress_helpers.R.

# ============================================================================
# Filter Samples by Lipid Correlation Threshold
# ============================================================================
