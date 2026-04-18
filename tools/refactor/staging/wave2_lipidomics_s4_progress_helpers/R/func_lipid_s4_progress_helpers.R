#' FilteringProgressLipidomics S4 Class
#'
#' @description
#' An S4 class to track lipidomics data filtering progress through QC steps.
#'
#' @export
setClass("FilteringProgressLipidomics",
    slots = c(
        steps = "character",
        assay_names = "list",
        n_lipids_per_assay = "list",
        n_lipids_total = "numeric",
        detected_per_sample = "list",
        missingness_per_assay = "list",
        sum_intensity_per_sample = "list",
        cv_distribution_per_assay = "list",
        is_metrics_per_assay = "list"
    )
)

#' @title Initialize or Retrieve Global Lipidomics Filtering Progress Object
#' @description Checks for a global object named
#' filtering_progress_lipidomics
#'              of class FilteringProgressLipidomics. If it doesn't exist,
#'              it creates and assigns a new one to the global environment.
#'
#' @return The global FilteringProgressLipidomics object.
#' @keywords internal
#' @noRd
#' @export
getFilteringProgressLipidomics <- function() {
    if (!exists("filtering_progress_lipidomics", envir = .GlobalEnv)) {
        filtering_progress_lipidomics <- new("FilteringProgressLipidomics")
        assign("filtering_progress_lipidomics", filtering_progress_lipidomics, envir = .GlobalEnv)
        message("Initialized global 'filtering_progress_lipidomics' object.")
    }
    get("filtering_progress_lipidomics", envir = .GlobalEnv)
}

