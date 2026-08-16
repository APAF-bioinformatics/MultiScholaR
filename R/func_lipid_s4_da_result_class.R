# ==========================================
# Content from lipid_da_analysis_wrapper.R
# ==========================================
#' LipidomicsDifferentialAbundanceResults S4 Class
#'
#' @description
#' S4 class to store essential results from lipidomics differential abundance analysis.
#' This class contains the original data object, fitted model, and results table.
#'
#' @slot theObject The original LipidomicsAssayData object used for analysis
#' @slot fit.eb The fitted eBayes model from limma analysis (ANY type to allow limma MArrayLM when available)
#' @slot contrasts_results_table Data frame with differential abundance statistics
#' @slot num_sig_diff_exp_bar_plot Per-contrast significant-feature count plots
#' @slot num_sig_diff_table Summary table of significant-feature counts
#' @slot volcano_plot Per-contrast static volcano plots
#' @slot interactive_volcano_plot Per-contrast interactive volcano plots
#' @slot p_value_dist_plot Per-contrast p-value distribution plots
#' @slot results_table_long Differential abundance results in long format
#' @slot results_table_wide Differential abundance results in wide format
#'
#' @export
setClass("LipidomicsDifferentialAbundanceResults",
    slots = c(
        theObject = "LipidomicsAssayData",
        fit.eb = "ANY",
        contrasts_results_table = "list",
        num_sig_diff_exp_bar_plot = "list",
        num_sig_diff_table = "data.frame",
        volcano_plot = "list",
        interactive_volcano_plot = "list",
        p_value_dist_plot = "list",
        results_table_long = "data.frame",
        results_table_wide = "data.frame"
    ),
    prototype = list(
        theObject = NULL,
        fit.eb = NULL,
        contrasts_results_table = list(),
        num_sig_diff_exp_bar_plot = list(),
        num_sig_diff_table = data.frame(),
        volcano_plot = list(),
        interactive_volcano_plot = list(),
        p_value_dist_plot = list(),
        results_table_long = data.frame(),
        results_table_wide = data.frame()
    )
)
