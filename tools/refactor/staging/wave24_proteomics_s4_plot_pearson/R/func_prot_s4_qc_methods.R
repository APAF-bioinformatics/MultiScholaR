## ----------------------------------------------------------------------------------------------------------------------------------------------------------------------
#' @title Plot Pearson Correlation
#' @param theObject is an object of the type ProteinQuantitativeData
#' @param tech_rep_remove_regex DEPRECATED - use exclude_pool_samples instead
#' @param correlation_group is the group where every pair of samples are compared
#' @param exclude_pool_samples Logical. If TRUE (default), automatically exclude samples from groups containing "Pool" or "QC" in their name from correlation analysis.
#' @export
setMethod(
  f = "plotPearson",
  signature = "ProteinQuantitativeData",
  definition = function(theObject, tech_rep_remove_regex = NULL, correlation_group = NA, exclude_pool_samples = TRUE) {
    correlation_group_to_use <- correlation_group

    if (is.na(correlation_group)) {
      correlation_group_to_use <- theObject@technical_replicate_id
    }

    # Handle deprecated parameter (backward compatibility)
    if (!is.null(tech_rep_remove_regex)) {
      message("*** plotPearson: WARNING - tech_rep_remove_regex is deprecated, use exclude_pool_samples instead ***")
    }

    correlation_vec <- pearsonCorForSamplePairs(theObject,
      tech_rep_remove_regex = tech_rep_remove_regex,
      correlation_group = correlation_group_to_use,
      exclude_pool_samples = exclude_pool_samples
    )

    pearson_plot <- correlation_vec |>
      ggplot(aes(pearson_correlation)) +
      geom_histogram(breaks = seq(min(round(correlation_vec$pearson_correlation - 0.5, 2), na.rm = TRUE), 1, 0.001)) +
      scale_y_continuous(breaks = seq(0, 4, 1), limits = c(0, 4)) +
      xlab("Pearson Correlation") +
      ylab("Counts") +
      theme(
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank()
      )

    return(pearson_plot)
  }
)

