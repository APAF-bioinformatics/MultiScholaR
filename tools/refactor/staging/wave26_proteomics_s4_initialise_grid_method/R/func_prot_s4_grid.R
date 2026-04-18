#' @export
setMethod(
  "InitialiseGrid",
  signature(dummy = "ANY"),
  function(dummy = NULL) {
    new("GridPlotData",
      pca_plots = list(),
      density_plots = list(),
      rle_plots = list(),
      pearson_plots = list(),
      cancor_plots = list(),
      limpa_plots = list(),
      pca_titles = list(),
      density_titles = list(),
      rle_titles = list(),
      pearson_titles = list(),
      cancor_titles = list(),
      limpa_titles = list()
    )
  }
)

