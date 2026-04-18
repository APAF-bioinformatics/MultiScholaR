#' FilteringProgressMetabolomics S4 Class
#'
#' @description
#' An S4 class to track metabolomics data filtering progress through QC steps.
#'
#' @export
setClass("FilteringProgressMetabolomics",
    slots = c(
        steps = "character",
        assay_names = "list",
        n_metabolites_per_assay = "list",
        n_metabolites_total = "numeric",
        detected_per_sample = "list",
        missingness_per_assay = "list",
        sum_intensity_per_sample = "list",
        cv_distribution_per_assay = "list",
        is_metrics_per_assay = "list"
    )
)

