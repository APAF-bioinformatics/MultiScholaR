#' Plot Top MOFA Weights for a Specific View and Factor
#'
#' @description
#' This function visualizes the top 20 features (e.g., genes, proteins)
#' with the largest absolute weights for a given factor and view from a
#' Multi-Omics Factor Analysis (MOFA) model. The weights are presented
#' as a bar plot, colored by their direction (positive or negative).
#' The gene/protein names are cleaned by removing "_transcriptome" or
#' "_proteome" suffixes.
#'
#' The plot is also saved to a file using the `savePlot` helper function,
#' which expects `project_dirs`, `omic_type`, and `experiment_label` to be
#' defined in the calling environment to determine the output path.
#'
#' @param model A trained MOFA model object, typically the output of
#'   `MOFA2::run_mofa()` or a similar MOFA model fitting function. It is
#'   expected to be compatible with `MOFA2::get_weights()`.
#' @param view Character string. The name of the view (e.g., "proteomics",
#'   "transcriptomics") within the MOFA model for which to extract and plot weights.
#' @param factor_level An unquoted expression specifying the factor column name
#'   (e.g., `Factor1`, `Factor2`) from the weights matrix of the specified view.
#'   Defaults to `Factor1`. This is used with `dplyr`'s non-standard evaluation.
#'
#' @return A `ggplot` object representing the bar plot of the top MOFA weights.
#'
#' @importFrom dplyr mutate as_data_frame case_when arrange desc row_number filter
#' @importFrom stringr str_replace
#' @importFrom tibble rownames_to_column
#' @importFrom ggplot2 ggplot aes geom_bar scale_fill_discrete theme_bw theme element_text xlab ylab coord_flip ggtitle
#' @importFrom MOFA2 get_weights
#' @export
#'
#' @examples
#' \dontrun{
#' # This example assumes a trained MOFA model ('dummy_mofa_model') exists,
#' # and necessary helper functions/variables for 'savePlot' are available
#' # (e.g., 'project_dirs', 'omic_type', 'experiment_label', 'savePlot' function).
#' # It also assumes MOFA2 and tidyverse packages are loaded.
#'
#' # Minimal setup for example:
#' # library(MOFA2) # if not loaded
#' # library(dplyr)
#' # library(ggplot2)
#'
#' # # Create a dummy MOFA model (replace with your actual model)
#' # set.seed(123)
#' # view1_weights <- matrix(rnorm(30*2), ncol=2, dimnames=list(paste0("Feat",1:30,"_transcriptome"), c("Factor1","Factor2")))
#' # dummy_mofa_model <- list(weights = list(transcriptome = view1_weights))
#' # if (!exists("get_weights", mode = "function") & requireNamespace("MOFA2", quietly=TRUE)) {
#' #   get_weights <- MOFA2::get_weights # Ensure get_weights is available
#' # } else if (!exists("get_weights", mode="function")) {
#' #   get_weights <- function(model, views) model$weights[[views]] # Simple mock
#' # }
#' #
#' # # Mock project_dirs and other globals for savePlot
#' # project_dirs <- list(); omic_type <- "integration"; experiment_label <- "example"
#' # project_dirs[[paste0(omic_type, "_", experiment_label)]] <- list(mofa_plots_dir = tempdir())
#' # savePlot <- function(plot, base_path, plot_name, formats) {
#' #   ggsave(file.path(base_path, paste0(plot_name, ".", formats[1])), plot)
#' #   message("Plot saved to: ", file.path(base_path, paste0(plot_name, ".", formats[1])))
#' # }
#'
#' # mofa_weights_plot <- plotMofaWeights(
#' #   model = dummy_mofa_model,
#' #   view = "transcriptome",
#' #   factor_level = Factor1
#' # )
#' # print(mofa_weights_plot)
#' }
plotMofaWeights <- function( model, view, factor_level = Factor1) {
  
  input_table <- get_weights(model)[[view]] 
  
  view_weights_helper <- input_table |>
    as.data.frame () |>
    dplyr::mutate( abs_factor1 = abs({{factor_level}})) |>
    mutate( Direction = case_when( sign({{factor_level}}) == 1 ~ "Positive"
                                   , TRUE ~ "Negative") ) |>
    mutate( Direction = factor( Direction, levels=c("Positive", "Negative"))) |>
    arrange( desc(abs_factor1) ) |> 
    mutate( rank = row_number()) |>
    
    dplyr::filter (rank <= 20) |>
    rownames_to_column("gene_symbol") |>
    mutate( gene_symbol = str_replace( gene_symbol, "_transcriptome$", "")) |>
    mutate( gene_symbol = str_replace( gene_symbol, "_proteome$", "")) 
  
  
  plotMofaWeightsAll <- view_weights_helper |>
    mutate( gene_symbol = factor( gene_symbol, levels=rev(unique( view_weights_helper$gene_symbol) )) ) |>
    ggplot(aes( gene_symbol, abs(Factor1), fill= Direction ) ) +
    geom_bar(stat = "identity") +
    scale_fill_discrete( type= c("red", "blue")) +
    theme_bw() +
    theme (axis.text.x = element_text (angle = 90, vjust = 1)) +
    xlab("Molecules") +
    ylab( "Factor 1 Weights") + 
    coord_flip() +
    ggtitle( view )
  
  savePlot(
    plot = plotMofaWeightsAll,
    base_path = project_dirs[[paste0(omic_type, "_", experiment_label)]]$mofa_plots_dir,
    plot_name = paste0(view, "_plot_top_weights_cln"),
    formats = c("png", "pdf")
  )
  
  return(plotMofaWeightsAll)
}