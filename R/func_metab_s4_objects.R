# ============================================================================
# func_metab_s4_objects.R
# ============================================================================
# Purpose: Metabolomics S4 class definitions and methods
#
# This file contains S4 class definitions and methods specific to metabolomics,
# including MetaboliteAssayData, MetabolomicsDifferentialAbundanceResults classes,
# and all metabolomics-specific S4 methods.
#
# Consolidated from:
# - metaboliteVsSamplesS4Objects.R (MetaboliteAssayData class + 16 methods)
# - metabolite_da_analysis_wrapper.R (MetabolomicsDifferentialAbundanceResults class + 2 methods)
# - metabolite_normalization.R (logTransformAssays, normaliseUntransformedData)
# - metabolite_qc.R (metaboliteIntensityFiltering)
# - QC_visualisation.R (FilteringProgressMetabolomics class)
#
# Dependencies:
# - methods package
# - func_general_s4_generics.R (for generic definitions)
# ============================================================================

# ==========================================
# Content from metaboliteVsSamplesS4Objects.R
# ==========================================
#' Metabolite Assay Data S4 Class
#'
#' An S4 class to store and manage multiple metabolomics quantitative datasets
#' derived from different platforms or ionization modes, along with associated
#' experimental design and metadata.
#'
#' @slot metabolite_data A named list of data frames. Each data frame contains
#'   quantitative data for one assay (e.g., LCMS-Pos). Rows should represent
#'   metabolites/features, and columns should represent samples. Each data frame
#'   MUST contain the column specified by `metabolite_id_column`. The names of
#'   the list elements should describe the assay source (e.g., "LCMS_Pos").
#' @slot metabolite_id_column Character string. The name of the column within
#'   each assay data frame that contains the **primary feature identifier** (used
#'   for uniquely identifying rows in quantitative analysis).
#' @slot annotation_id_column Character string. The name of the column (often
#'   in feature metadata, but could be in assay tables) that contains the
#'   **biological or database annotation identifier** (e.g., HMDB ID, KEGG ID,
#'   chemical name, CHEBI ID). This is the ID typically used for downstream
#'   biological mapping.
#' @slot database_identifier_type Character string. Describes the type or format
#'   of identifiers found in the column specified by `annotation_id_column`
#'   (e.g., "HMDB", "KEGG", "CHEBI", "Mixed_CHEBI_Unknown", "InternalName").
#' @slot internal_standard_regex Character string. A regular expression used to
#'   identify features that are internal standards based on their identifier in
#'   the `metabolite_id_column` or `annotation_id_column`. Set to `NA_character_` or `""`
#'   if not applicable or no internal standards are used.
#' @slot design_matrix A data frame containing the experimental design. Must
#'   include the column specified by `sample_id`.
#' @slot sample_id Character string. The name of the column in `design_matrix`
#'   that contains unique sample identifiers. These identifiers must correspond
#'   to the sample column names in the assay data frames.
#' @slot group_id Character string. The name of the column in `design_matrix`
#'   that defines the primary experimental groups or conditions.
#' @slot technical_replicate_id Character string. The name of the column in
#'   `design_matrix` that identifies technical replicates, if applicable. Use
#'   `NA_character_` if there are no technical replicates.
#' @slot args A list, typically populated from a configuration file, holding
#'   parameters used during processing.
#'
#' @importFrom methods setClass slotNames slot new
#' @importFrom dplyr pull distinct
#' @importFrom rlang sym !!
#' @exportClass MetaboliteAssayData
MetaboliteAssayData <- setClass("MetaboliteAssayData",
    slots = c(
        metabolite_data = "list",
        metabolite_id_column = "character",
        annotation_id_column = "character",
        database_identifier_type = "character",
        internal_standard_regex = "character",
        design_matrix = "data.frame",
        sample_id = "character",
        group_id = "character",
        technical_replicate_id = "character",
        args = "list"
    ),
    prototype = list(
        metabolite_data = list(),
        metabolite_id_column = "database_identifier",
        annotation_id_column = "metabolite_identification",
        database_identifier_type = "Unknown",
        internal_standard_regex = NA_character_,
        design_matrix = data.frame(),
        sample_id = "Sample_ID",
        group_id = "group",
        technical_replicate_id = NA_character_,
        args = list()
    ),
    validity = function(object) {
        errors <- character()
        # --- Get required info ---
        sample_id_col <- object@sample_id
        metabolite_id_col <- object@metabolite_id_column
        design_matrix <- object@design_matrix
        metabolite_data <- object@metabolite_data

        # --- Basic slot type checks (as before) ---
        if (!is.list(metabolite_data)) {
            errors <- c(errors, "`metabolite_data` must be a list.")
        } else if (length(metabolite_data) > 0 && !all(sapply(metabolite_data, is.data.frame))) {
            errors <- c(errors, "All elements in `metabolite_data` must be data frames.")
        }
        if (!is.character(object@metabolite_id_column) || length(object@metabolite_id_column) != 1) {
            errors <- c(errors, "`metabolite_id_column` must be a single character string.")
        }
        if (!is.character(object@annotation_id_column) || length(object@annotation_id_column) != 1) {
            errors <- c(errors, "`annotation_id_column` must be a single character string.")
        }
        if (!is.character(object@database_identifier_type) || length(object@database_identifier_type) != 1) {
            errors <- c(errors, "`database_identifier_type` must be a single character string.")
        }
        if (!is.character(object@internal_standard_regex) || length(object@internal_standard_regex) != 1) {
            errors <- c(errors, "`internal_standard_regex` must be a single character string (can be NA_character_).")
        }
        if (!is.data.frame(object@design_matrix)) {
            errors <- c(errors, "`design_matrix` must be a data frame.")
        }
        if (!is.character(object@sample_id) || length(object@sample_id) != 1) {
            errors <- c(errors, "`sample_id` must be a single character string.")
        }
        if (!is.character(object@group_id) || length(object@group_id) != 1) {
            errors <- c(errors, "`group_id` must be a single character string.")
        }
        if (!is.character(object@technical_replicate_id) || length(object@technical_replicate_id) != 1) {
            errors <- c(errors, "`technical_replicate_id` must be a single character string (can be NA_character_).")
        }
        if (!is.list(object@args)) {
            errors <- c(errors, "`args` must be a list.")
        }

        # --- Content Checks ---
        # Check design matrix first
        if (!is.data.frame(design_matrix)) {
            # Error already added by basic checks, but prevent further processing
        } else if (!sample_id_col %in% colnames(design_matrix)) {
            errors <- c(errors, paste0("`sample_id` column ('", sample_id_col, "') not found in `design_matrix`."))
        } else {
            # Get unique, sorted sample IDs from design matrix (ensure character)
            samples_in_design <- tryCatch(
                design_matrix[[sample_id_col]] |> as.character() |> unique() |> sort(),
                error = function(e) {
                    errors <- c(errors, "Error extracting sample IDs from design matrix.")
                    character(0)
                }
            )

            if (length(samples_in_design) == 0 && length(errors) == 0) {
                errors <- c(errors, "No valid sample IDs found in design matrix.")
            }

            # Proceed with assay checks only if design matrix looks okay so far
            if (length(metabolite_data) > 0 && length(errors) == 0) {
                assay_names_vec <- names(metabolite_data)
                if (is.null(assay_names_vec)) assay_names_vec <- paste0("Assay_", seq_along(metabolite_data))
                names(metabolite_data) <- assay_names_vec # Ensure the list is named for lapply output

                # Use lapply to check each assay and collect results/errors
                assay_check_results <- lapply(assay_names_vec, function(assay_name) {
                    assay_df <- metabolite_data[[assay_name]]
                    assay_errors <- character()

                    # Check metabolite ID column exists
                    if (!metabolite_id_col %in% colnames(assay_df)) {
                        assay_errors <- c(assay_errors, paste0("Assay '", assay_name, "': `metabolite_id_column` ('", metabolite_id_col, "') not found."))
                    }

                    # Identify actual sample columns in the assay
                    assay_colnames <- colnames(assay_df)
                    actual_sample_cols_in_assay <- intersect(assay_colnames, samples_in_design) |> sort()

                    # Store results for later checks
                    list(
                        errors = assay_errors,
                        sample_cols = actual_sample_cols_in_assay
                    )
                })

                # Aggregate errors from individual assay checks
                all_assay_errors <- unlist(lapply(assay_check_results, `[[`, "errors"))
                errors <- c(errors, all_assay_errors)

                # Perform cross-assay consistency checks if no individual errors found yet
                if (length(errors) == 0 && length(assay_check_results) > 1) {
                    first_assay_samples <- assay_check_results[[1]]$sample_cols
                    consistency_check <- sapply(assay_check_results[-1], function(res) {
                        identical(res$sample_cols, first_assay_samples)
                    })
                    if (!all(consistency_check)) {
                        mismatched_assays <- assay_names_vec[c(FALSE, !consistency_check)] # Get names of inconsistent assays
                        errors <- c(errors, paste0("Actual sample columns differ between assays. First mismatch found in: ", mismatched_assays[1]))
                    }
                }

                # Perform comparison with design matrix if no errors found yet
                if (length(errors) == 0 && length(assay_check_results) >= 1) {
                    first_assay_samples <- assay_check_results[[1]]$sample_cols # Get samples from first (or only) assay
                    if (!identical(first_assay_samples, samples_in_design)) {
                        errors <- c(errors, paste0("Sample columns in assays do not exactly match unique sample IDs ('", sample_id_col, "') in `design_matrix`."))
                        # Add more detail:
                        missing_in_assay <- setdiff(samples_in_design, first_assay_samples)
                        extra_in_assay <- setdiff(first_assay_samples, samples_in_design)
                        if (length(missing_in_assay) > 0) errors <- c(errors, paste0("   Samples in design_matrix missing from assay columns: ", paste(utils::head(missing_in_assay, 10), collapse = ", "), ifelse(length(missing_in_assay) > 10, "...", "")))
                        if (length(extra_in_assay) > 0) errors <- c(errors, paste0("   Sample columns in assay not found in design_matrix: ", paste(utils::head(extra_in_assay, 10), collapse = ", "), ifelse(length(extra_in_assay) > 10, "...", "")))
                    }
                }
            }
        }

        # --- Final Check ---
        if (length(errors) == 0) TRUE else errors
    }
)


#' Create MetaboliteAssayData Object
#'
#' Constructor function for the MetaboliteAssayData class.
#'
#' @param metabolite_data Named list of data frames (assays).
#' @param design_matrix Experimental design data frame.
#' @param metabolite_id_column Name of the **primary feature ID** column within assays (e.g., `"database_identifier"`).
#' @param annotation_id_column Name of the **annotation ID** column (e.g., `"metabolite_identification"`).
#' @param sample_id Name of the sample ID column in design_matrix and assays.
#' @param group_id Name of the group column in design_matrix.
#' @param technical_replicate_id Name of the technical replicate column in design_matrix (use NA_character_ if none).
#' @param database_identifier_type Type of identifier in the `annotation_id_column` (e.g., `"Mixed_CHEBI_Unknown"`).
#' @param internal_standard_regex Regex to identify internal standards. Use `NA_character_` or `""` if none.
#' @param args List of arguments (e.g., from config).
#'
#' @return A MetaboliteAssayData object.
#' @export
#' @examples
#' \dontrun{
#' # Assuming lcms_pos_df, lcms_neg_df, gcms_df are data frames
#' # with 'Metabolite' as ID column and samples as other columns
#' # Assuming design_df has 'SampleID', 'Group', 'Replicate' columns
#' assays_list <- list(
#'     LCMS_Pos = lcms_pos_df,
#'     LCMS_Neg = lcms_neg_df,
#'     GCMS = gcms_df
#' )
#' config <- list(...) # Your config list
#'
#' met_assay_obj <- createMetaboliteAssayData(
#'     metabolite_data = assays_list,
#'     design_matrix = design_df,
#'     metabolite_id_column = "Metabolite",
#'     sample_id = "SampleID",
#'     group_id = "Group",
#'     technical_replicate_id = "Replicate",
#'     database_identifier_type = "InternalName",
#'     internal_standard_regex = "^IS_",
#'     args = config
#' )
#' }
createMetaboliteAssayData <- function(
  metabolite_data,
  design_matrix,
  metabolite_id_column = "database_identifier",
  annotation_id_column = "metabolite_identification",
  sample_id = "Sample_ID",
  group_id = "group",
  technical_replicate_id = NA_character_,
  database_identifier_type = "Unknown",
  internal_standard_regex = NA_character_,
  args = list()
) {
    # Perform basic checks before creating the object
    stopifnot(is.list(metabolite_data))
    stopifnot(all(sapply(metabolite_data, is.data.frame)))
    stopifnot(is.data.frame(design_matrix))
    # Add more checks as needed...

    obj <- new("MetaboliteAssayData",
        metabolite_data = metabolite_data,
        metabolite_id_column = metabolite_id_column,
        annotation_id_column = annotation_id_column,
        database_identifier_type = database_identifier_type,
        internal_standard_regex = internal_standard_regex,
        design_matrix = design_matrix,
        sample_id = sample_id,
        group_id = group_id,
        technical_replicate_id = technical_replicate_id,
        args = args
    )
    # Validity check is automatically called by 'new'
    return(obj)
}

## -----------------------------------------------------------------------------
## Plotting Methods for MetaboliteAssayData
## -----------------------------------------------------------------------------









# ------------------------------------------------------- #




## -----------------------------------------------------------------------------
## Negative Control Selection Methods for MetaboliteAssayData
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
# setMethod(f = "createGridQCMetabolomics",
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
# Content from metabolite_normalization.R
# ==========================================





#-------------------------------------------------------------------------------


#' Find Duplicate Feature IDs in MetaboliteAssayData Assays
#'
#' Iterates through each assay tibble in a MetaboliteAssayData object
#' and identifies any duplicated values in the specified feature ID column.
#'
#' @param theObject A MetaboliteAssayData object.
#'
#' @return A named list. Each element corresponds to an assay in the input object.
#'   If duplicates are found in an assay, the element will be a tibble showing
#'   the duplicated identifiers and their counts. If no duplicates are found,
#'   the element will be NULL.



## -----------------------------------------------------------------------------
## Helper Function to Resolve Duplicate Features by Intensity
## -----------------------------------------------------------------------------


