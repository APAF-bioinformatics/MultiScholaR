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
#' @param owner Optional workflow/session environment that owns progress state.
#' @return The FilteringProgressLipidomics object for the selected owner.
#' @keywords internal
#' @noRd
#' @export
getFilteringProgressLipidomics <- function(owner = NULL) {
    if (!is.null(owner)) {
        progress <- owner$filtering_progress_lipidomics
        if (!methods::is(progress, "FilteringProgressLipidomics")) {
            progress <- methods::new("FilteringProgressLipidomics")
            owner$filtering_progress_lipidomics <- progress
        }
        return(progress)
    }
    if (!exists("filtering_progress_lipidomics", envir = .GlobalEnv)) {
        filtering_progress_lipidomics <- methods::new(
            "FilteringProgressLipidomics"
        )
        assign("filtering_progress_lipidomics", filtering_progress_lipidomics, envir = .GlobalEnv)
        message("Initialized global 'filtering_progress_lipidomics' object.")
    }
    get("filtering_progress_lipidomics", envir = .GlobalEnv)
}

releaseFilteringProgressLipidomics <- function(owner = NULL) {
    if (!is.null(owner)) owner$filtering_progress_lipidomics <- NULL
    invisible(TRUE)
}
