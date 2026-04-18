# ----------------------------------------------------------------------------
# updateLipidFiltering
# ----------------------------------------------------------------------------
#' @title Update and Visualize Lipidomics Filtering Progress
#' @description Tracks and visualizes the impact of filtering steps on lipidomics
#'              data, updating a global `FilteringProgressLipidomics` object.
#'              Generates QC plots summarizing the changes.
#'
#' @details
#' This function serves as the primary interface for tracking lipidomics QC.
#' It performs the following:
#' \itemize{
#'   \item Initializes or retrieves the global `FilteringProgressLipidomics` object.
#'   \item Takes a `LipidomicsAssayData` object (or similar S4 containing assays list
#'         and design matrix) as input for the current processing step.
#'   \item Extracts the list of assay data frames/tibbles and the design matrix.
#'   \item For each assay, calls helper functions (`countUniqueLipids`,
#'         `countLipidsPerSample`, `calculateLipidMissingness`,
#'         `calculateLipidSumIntensityPerSample`, `calculateLipidCVs`,
#'         `getLipidInternalStandardMetrics`) to calculate QC metrics.
#'   \item Calls `calculateTotalUniqueLipidsAcrossAssays` for an overall count.
#'   \item Updates the global `FilteringProgressLipidomics` object with the
#'         calculated metrics for the given `step_name`.
#'   \item Generates summary plots visualizing the tracked metrics across steps.
#'   \item Optionally saves plots to disk if `publication_graphs_dir` and a time
#'         directory are available.
#'   \item Returns either a combined grid plot or an invisible list of plots.
#' }
#'
#' **Important:** Relies on and modifies the global `filtering_progress_lipidomics` object.
#' Requires helper functions (defined previously in this file) to be available.
#' For plot saving, requires either `time_dir` in the global environment or access to
#' the appropriate directory via `project_dirs$omics_type$time_dir`.
#'
#' @param theObject A S4 object (e.g., `LipidomicsAssayData`, `SummarizedExperiment`,
#'                  `MultiAssayExperiment`) containing lipidomics data. Must provide
#'                  access to a list of assays (data frames/tibbles with lipid rows
#'                  and sample columns) and a colData/design matrix linking samples to groups.
#' @param step_name Character string uniquely identifying the current filtering step.
#' @param publication_graphs_dir Optional path for saving plots. If provided, the function
#'                              will try to find the corresponding time_dir.
#' @param omics_type Optional character string specifying the omics type (e.g., "lipidomics").
#'                  If provided and project_dirs exists in the global environment, will use
#'                  `project_dirs[[omics_type]]$time_dir` for plot saving.
#' @param time_dir Optional explicit path to the time directory. If provided, this overrides
#'                other methods of finding the time directory.
#' @param overwrite Logical, whether to overwrite existing data for `step_name`.
#' @param return_grid Logical, whether to return a `gridExtra` combined plot.
#' @param group_id_col Character, name of the column in `colData(theObject)` specifying groups.
#' @param sample_id_col Character, name of the sample ID column in `colData(theObject)`.
#' @keywords internal
#' @noRd
calculateLipidFilteringAssayMetrics <- function(current_assay_data,
                                                current_assay_name,
                                                design_matrix,
                                                group_id_col,
                                                sample_id_col,
                                                lipid_id_col,
                                                is_pattern = NULL,
                                                sample_columns = NULL) {
    if (!is.data.frame(current_assay_data)) {
        warning("Assay ", current_assay_name, " is not a data frame. Skipping metrics calculation.")
        return(list(
            n_lipids = 0,
            detected_per_sample = data.frame(Run = character(), n_detected = integer()),
            missingness = NA_real_,
            sum_intensity_per_sample = data.frame(Run = character(), sum_intensity = numeric()),
            cv_distribution = data.frame(lipid_id = character(), group = character(), cv = numeric()),
            is_metrics = data.frame(is_id = character(), mean_intensity = numeric(), cv = numeric())
        ))
    }

    n_met <- tryCatch(
        {
            countUniqueLipids(current_assay_data, lipid_id_col)
        },
        error = function(e) {
            stop(e)
        }
    )

    det_per_sample <- tryCatch(
        {
            countLipidsPerSample(current_assay_data, sample_id_col, lipid_id_col, sample_columns = sample_columns)
        },
        error = function(e) {
            stop(e)
        }
    )

    miss <- tryCatch(
        {
            calculateLipidMissingness(current_assay_data, sample_id_col, sample_columns = sample_columns)
        },
        error = function(e) {
            stop(e)
        }
    )

    sum_int <- tryCatch(
        {
            calculateLipidSumIntensityPerSample(current_assay_data, sample_id_col, sample_columns = sample_columns)
        },
        error = function(e) {
            stop(e)
        }
    )

    cv_dist <- tryCatch(
        {
            calculateLipidCVs(current_assay_data, design_matrix, group_id_col, NULL, sample_id_col, lipid_id_col, sample_columns = sample_columns)
        },
        error = function(e) {
            stop(e)
        }
    )

    is_met <- tryCatch(
        {
            getLipidInternalStandardMetrics(current_assay_data, is_pattern, lipid_id_col, sample_id_col, sample_columns = sample_columns)
        },
        error = function(e) {
            stop(e)
        }
    )

    list(
        n_lipids = n_met,
        detected_per_sample = det_per_sample,
        missingness = miss,
        sum_intensity_per_sample = sum_int,
        cv_distribution = cv_dist,
        is_metrics = is_met
    )
}

prepareLipidFilteringContext <- function(theObject,
                                         group_id_col = NULL,
                                         sample_id_col = NULL,
                                         lipid_id_col = NULL,
                                         is_pattern = NULL) {
    if (!isS4(theObject)) {
        stop("`theObject` must be an S4 object.")
    }

    if (inherits(theObject, "LipidomicsAssayData")) {
        design_matrix <- theObject@design_matrix
        assay_list <- theObject@lipid_data
        assay_names <- names(assay_list)
        if (is.null(assay_names)) assay_names <- paste0("Assay_", seq_along(assay_list))
        names(assay_list) <- assay_names
    } else {
        if (!requireNamespace("SummarizedExperiment", quietly = TRUE)) {
            stop("Package 'SummarizedExperiment' needed for this function to work.")
        }
        if (!("assays" %in% methods::slotNames(theObject) || canCoerce(theObject, "SummarizedExperiment"))) {
            stop("`theObject` must have an accessible `assays` method or slot.")
        }
        if (!("colData" %in% methods::slotNames(theObject) || canCoerce(theObject, "SummarizedExperiment"))) {
            stop("`theObject` must have an accessible `colData` method or slot.")
        }

        design_matrix <- SummarizedExperiment::colData(theObject)
        assay_list <- SummarizedExperiment::assays(theObject)
        assay_list <- lapply(assay_list, as.data.frame)
        if (is.null(names(assay_list))) assay_names <- paste0("Assay_", seq_along(assay_list)) else assay_names <- names(assay_list)
        names(assay_list) <- assay_names
    }

    if (is.null(group_id_col) && "group_id" %in% slotNames(theObject)) group_id_col <- theObject@group_id
    if (is.null(sample_id_col) && "sample_id" %in% slotNames(theObject)) sample_id_col <- theObject@sample_id
    if (is.null(lipid_id_col) && "lipid_id_column" %in% slotNames(theObject)) lipid_id_col <- theObject@lipid_id_column
    if (is.null(is_pattern) && "internal_standard_regex" %in% slotNames(theObject)) is_pattern <- theObject@internal_standard_regex

    if (is.null(group_id_col)) stop("group_id_col must be provided or accessible via theObject@group_id")
    if (is.null(sample_id_col)) stop("sample_id_col must be provided or accessible via theObject@sample_id")
    if (is.null(lipid_id_col)) stop("lipid_id_col must be provided or accessible via theObject@lipid_id_column")

    if (!sample_id_col %in% colnames(design_matrix)) {
        if (identical(rownames(design_matrix), as.character(design_matrix[[sample_id_col]]))) {
            warning("Sample ID column '", sample_id_col, "' not found, but rownames might match? Check object structure.")
        } else if (!is.null(rownames(design_matrix)) && sample_id_col == "Run") {
            message("Moving rownames of design matrix to column: ", sample_id_col)
            design_matrix <- as.data.frame(design_matrix)
            design_matrix[[sample_id_col]] <- rownames(design_matrix)
            rownames(design_matrix) <- NULL
        } else {
            warning("Sample ID column '", sample_id_col, "' not found in design matrix and rownames don't seem to match or weren't checked.")
        }
    }

    design_matrix <- as.data.frame(design_matrix)

    list(
        assay_list = assay_list,
        assay_names = assay_names,
        design_matrix = design_matrix,
        group_id_col = group_id_col,
        sample_id_col = sample_id_col,
        lipid_id_col = lipid_id_col,
        is_pattern = is_pattern,
        sample_columns = as.character(design_matrix[[sample_id_col]])
    )
}

finalizeLipidFilteringStep <- function(prog_met,
                                       step_name,
                                       assay_names,
                                       metrics_list_this_step,
                                       assay_list,
                                       lipid_id_col,
                                       overwrite = FALSE) {
    total_lipids <- tryCatch(
        {
            calculateTotalUniqueLipidsAcrossAssays(assay_list, lipid_id_col)
        },
        error = function(e) {
            stop(e)
        }
    )

    tryCatch(
        {
            updateFilteringProgressLipidomics(
                prog_met,
                step_name,
                assay_names,
                metrics_list_this_step,
                total_lipids,
                overwrite
            )
        },
        error = function(e) {
            stop(e)
        }
    )

    plot_list <- tryCatch(
        {
            generateLipidFilteringPlots(getFilteringProgressLipidomics())
        },
        error = function(e) {
            stop(e)
        }
    )

    list(
        total_lipids = total_lipids,
        plot_list = plot_list
    )
}

resolveLipidFilteringPlotSaveDir <- function(publication_graphs_dir = NULL,
                                             omics_type = NULL,
                                             time_dir = NULL) {
    actual_save_dir <- NULL

    message("--- Plot Saving Diagnostics ---")
    message(sprintf("Value of publication_graphs_dir argument in function call: %s", ifelse(is.null(publication_graphs_dir), "NULL", publication_graphs_dir)))
    message(sprintf("Value of omics_type argument in function call: %s", ifelse(is.null(omics_type), "NULL", omics_type)))
    message(sprintf("Value of time_dir argument in function call: %s", ifelse(is.null(time_dir), "NULL", time_dir)))
    message("Attempting to determine save directory using global project_dirs, omic_type, and experiment_label...")

    if (exists("project_dirs", envir = .GlobalEnv) &&
        exists("omic_type", envir = .GlobalEnv) &&
        exists("experiment_label", envir = .GlobalEnv)) {
        message("Global variables project_dirs, omic_type, and experiment_label found.")
        local_project_dirs <- get("project_dirs", envir = .GlobalEnv)
        local_omic_type <- get("omic_type", envir = .GlobalEnv)
        local_experiment_label <- get("experiment_label", envir = .GlobalEnv)
        message(sprintf("Global omic_type value: '%s', Global experiment_label value: '%s'", local_omic_type, local_experiment_label))

        omics_key <- paste0(local_omic_type, "_", local_experiment_label)
        message(sprintf("Constructed omics_key for project_dirs: '%s'", omics_key))

        if (omics_key %in% names(local_project_dirs) &&
            !is.null(local_project_dirs[[omics_key]]) &&
            "time_dir" %in% names(local_project_dirs[[omics_key]])) {
            message(sprintf("omics_key '%s' found in project_dirs and has a 'time_dir' entry.", omics_key))
            retrieved_time_dir <- local_project_dirs[[omics_key]]$time_dir
            message(sprintf("Retrieved time_dir from project_dirs: %s", ifelse(is.null(retrieved_time_dir), "NULL", retrieved_time_dir)))

            if (is.null(retrieved_time_dir) || !is.character(retrieved_time_dir) || !nzchar(retrieved_time_dir)) {
                warning(sprintf("project_dirs[['%s']]$time_dir is NULL, not a character string, or empty. Plots will not be saved.", omics_key))
                message("Reason: retrieved_time_dir is invalid.")
            } else {
                actual_save_dir <- retrieved_time_dir
                message(sprintf("Successfully set actual_save_dir for plot saving to: '%s'", actual_save_dir))
            }
        } else {
            warning(sprintf("Could not find omics_key '%s' in project_dirs, or it lacks a 'time_dir' entry. Plots will not be saved.", omics_key))
            message(sprintf(
                "Details: omics_key '%s' in names(project_dirs): %s. project_dirs[['%s']] is NULL: %s. 'time_dir' in names(project_dirs[['%s']]): %s.",
                omics_key, omics_key %in% names(local_project_dirs),
                omics_key, is.null(local_project_dirs[[omics_key]]),
                omics_key, if (omics_key %in% names(local_project_dirs) && !is.null(local_project_dirs[[omics_key]])) "time_dir" %in% names(local_project_dirs[[omics_key]]) else NA
            ))
            message("Available keys in project_dirs: ", paste(names(local_project_dirs), collapse = ", "))
        }
    } else {
        warning("One or more global variables ('project_dirs', 'omic_type', 'experiment_label') not found. Plots will not be saved.")
        message(sprintf(
            "Exists project_dirs: %s, Exists omic_type: %s, Exists experiment_label: %s",
            exists("project_dirs", envir = .GlobalEnv),
            exists("omic_type", envir = .GlobalEnv),
            exists("experiment_label", envir = .GlobalEnv)
        ))
    }

    actual_save_dir
}

saveLipidFilteringPlots <- function(plot_list,
                                    step_name,
                                    actual_save_dir = NULL,
                                    return_grid = FALSE,
                                    publication_graphs_dir = NULL) {
    if (!is.null(actual_save_dir)) {
        message(sprintf("Proceeding to save plots to derived directory: %s", actual_save_dir))
        if (!dir.exists(actual_save_dir)) {
            dir.create(actual_save_dir, recursive = TRUE)
            message("Created directory for QC plots: ", actual_save_dir)
        }

        purrr::iwalk(plot_list, function(plot, plot_name) {
            filename <- file.path(actual_save_dir, sprintf("%s_%s.png", step_name, plot_name))
            message(sprintf("Saving plot: %s", filename))
            ggsave(filename,
                plot = plot,
                width = 10,
                height = 8,
                dpi = 300
            )
        })

        if (return_grid && length(plot_list) > 0 && !is.null(plot_list[[1]]) && inherits(plot_list[[1]], "ggplot")) {
            pdf(NULL)
            grid_plot_obj <- do.call(gridExtra::arrangeGrob, c(plot_list, ncol = 2))
            invisible(dev.off())
            filename_grid <- file.path(actual_save_dir, sprintf("%s_combined_plots.png", step_name))
            message(sprintf("Saving combined grid plot: %s", filename_grid))
            ggsave(filename_grid, plot = grid_plot_obj, width = 15, height = 15, dpi = 300)
        }
        message("Lipidomics QC plots saved to: ", actual_save_dir)
    } else {
        message("No valid save directory determined from global project_dirs. Plots will not be saved.")
        if (!is.null(publication_graphs_dir)) {
            warning("Function was called with a publication_graphs_dir path, but plot saving still failed because a valid time_dir could not be derived from global project_dirs.")
        }
    }
    message("--- End of Plot Saving Diagnostics ---")

    invisible(actual_save_dir)
}

returnLipidFilteringPlots <- function(plot_list,
                                      return_grid = FALSE) {
    if (return_grid) {
        has_plots <- length(plot_list) > 0
        has_first_plot <- has_plots && !is.null(plot_list[[1]])
        can_build_grid <- has_first_plot && inherits(plot_list[[1]], "ggplot")

        if (!can_build_grid) {
            return(NULL)
        }

        grid_plot_obj <- tryCatch(
            {
                pdf(NULL)
                result <- do.call(gridExtra::arrangeGrob, c(plot_list, ncol = 2))
                invisible(dev.off())
                result
            },
            error = function(e) {
                tryCatch(invisible(dev.off()), error = function(e2) NULL)
                NULL
            }
        )

        return(grid_plot_obj)
    }

    if (length(plot_list) > 0) {
        message("Printing plots individually as return_grid is FALSE or grid could not be formed.")
        purrr::walk(plot_list, function(plot) {
            if (inherits(plot, "ggplot")) {
                print(plot)
            } else {
                message("Encountered a non-ggplot object in plot_list when trying to print individually.")
            }
        })
    } else {
        message("Plot list is empty, nothing to print.")
    }

    invisible(plot_list)
}

