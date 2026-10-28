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
getFilteringProgressLipidomics <- function(owner = workflowOmicsProgressOwner("lipidomics")) {
    # KIT-030: the owner is the lipidomics bus registered on the Shiny session, else the global
    # environment (scripts and workbooks unchanged).
    value <- shiny::isolate(owner$filtering_progress_lipidomics)
    if (is.null(value)) {
        value <- new("FilteringProgressLipidomics")
        owner$filtering_progress_lipidomics <- value
        message("Initialized global 'filtering_progress_lipidomics' object.")
    }
    value
}

