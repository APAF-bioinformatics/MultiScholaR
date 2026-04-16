# ----------------------------------------------------------------------------
# getFilteringProgressMetabolomics
# ----------------------------------------------------------------------------
#' @title Initialize or Retrieve Global Metabolomics Filtering Progress Object
#' @description Checks for a global object named `filtering_progress_metabolomics`
#'              of class `FilteringProgressMetabolomics`. If it doesn't exist,
#'              it creates and assigns a new one to the global environment.
#'
#' @return The global `FilteringProgressMetabolomics` object.
#' @keywords internal
#' @noRd
#' @export
getFilteringProgressMetabolomics <- function() {
    if (!exists("filtering_progress_metabolomics", envir = .GlobalEnv)) {
        filtering_progress_metabolomics <- new("FilteringProgressMetabolomics")
        assign("filtering_progress_metabolomics", filtering_progress_metabolomics, envir = .GlobalEnv)
        message("Initialized global 'filtering_progress_metabolomics' object.") # Optional message
    }
    get("filtering_progress_metabolomics", envir = .GlobalEnv)
}

# ----------------------------------------------------------------------------
# updateFilteringProgressMetabolomics
# ----------------------------------------------------------------------------
#' @title Update the Global Metabolomics Filtering Progress Object
#' @description Modifies the global `filtering_progress_metabolomics` object
#'              by adding or overwriting data for a specific step.
#'
#' @param prog_met The `FilteringProgressMetabolomics` object (retrieved globally).
#' @param step_name The name of the step to add or update.
#' @param current_assay_names Character vector of assay names for this step.
#' @param metrics_list A nested list containing the calculated metrics for each assay
#'                      for the current step.
#' @param total_metabolites The total unique metabolites calculated across assays for this step.
#' @param overwrite Logical, whether to overwrite if `step_name` exists.
#'
#' @return The updated `FilteringProgressMetabolomics` object (invisibly).
#'         Has the side effect of updating the global object.
#' @keywords internal
#' @noRd
#' @export
updateFilteringProgressMetabolomics <- function(prog_met,
                                                step_name,
                                                current_assay_names,
                                                metrics_list,
                                                total_metabolites,
                                                overwrite = FALSE) {
    if (step_name %in% prog_met@steps) {
        if (!overwrite) {
            stop("Step name '", step_name, "' already exists in filtering_progress_metabolomics. Use overwrite = TRUE to replace it.")
        }
        idx <- which(prog_met@steps == step_name)

        # Overwrite existing data
        prog_met@assay_names[[idx]] <- current_assay_names
        prog_met@n_metabolites_per_assay[[idx]] <- lapply(metrics_list, `[[`, "n_metabolites")
        prog_met@n_metabolites_total[idx] <- total_metabolites
        prog_met@detected_per_sample[[idx]] <- lapply(metrics_list, `[[`, "detected_per_sample")
        prog_met@missingness_per_assay[[idx]] <- lapply(metrics_list, `[[`, "missingness")
        prog_met@sum_intensity_per_sample[[idx]] <- lapply(metrics_list, `[[`, "sum_intensity_per_sample")
        prog_met@cv_distribution_per_assay[[idx]] <- lapply(metrics_list, `[[`, "cv_distribution")
        prog_met@is_metrics_per_assay[[idx]] <- lapply(metrics_list, `[[`, "is_metrics")
    } else {
        # Append new data
        prog_met@steps <- c(prog_met@steps, step_name)
        prog_met@assay_names <- c(prog_met@assay_names, list(current_assay_names))
        prog_met@n_metabolites_per_assay <- c(prog_met@n_metabolites_per_assay, list(lapply(metrics_list, `[[`, "n_metabolites")))
        prog_met@n_metabolites_total <- c(prog_met@n_metabolites_total, total_metabolites)
        prog_met@detected_per_sample <- c(prog_met@detected_per_sample, list(lapply(metrics_list, `[[`, "detected_per_sample")))
        prog_met@missingness_per_assay <- c(prog_met@missingness_per_assay, list(lapply(metrics_list, `[[`, "missingness")))
        prog_met@sum_intensity_per_sample <- c(prog_met@sum_intensity_per_sample, list(lapply(metrics_list, `[[`, "sum_intensity_per_sample")))
        prog_met@cv_distribution_per_assay <- c(prog_met@cv_distribution_per_assay, list(lapply(metrics_list, `[[`, "cv_distribution")))
        prog_met@is_metrics_per_assay <- c(prog_met@is_metrics_per_assay, list(lapply(metrics_list, `[[`, "is_metrics")))
    }

    # Update the global object
    assign("filtering_progress_metabolomics", prog_met, envir = .GlobalEnv)
    invisible(prog_met)
}

