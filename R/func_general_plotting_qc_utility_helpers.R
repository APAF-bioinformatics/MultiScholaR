# ----------------------------------------------------------------------------
# plotPeptidesProteinsCountsPerSampleHelper
# ----------------------------------------------------------------------------
#' plotPeptidesProteinsCountsPerSampleHelper
#' @description Plot the number of proteins and peptides identified per sample
#' @importFrom purrr map map2 map_chr pmap
#' @importFrom dplyr filter select mutate pull rename distinct
#' @importFrom ggplot2 ggplot aes
#' @importFrom ggrepel geom_text_repel
#' @importFrom RColorBrewer brewer.pal
#' @importFrom grid grid.draw convertX convertY gpar unit
#' @importFrom gridExtra arrangeGrob
#' @importFrom graphics boxplot
#' @importFrom rlang base_env
#' @importFrom utils str
#' @importFrom htmlwidgets saveWidget
#' @importFrom GGally ggpairs
#' @importFrom ggpubr ggarrange
#' @export
plotPeptidesProteinsCountsPerSampleHelper <- function(
  input_table,
  intensity_column = Peptide.RawQuantity,
  protein_id_column = Protein.Ids,
  peptide_id_column = Stripped.Sequence,
  sample_id_column = Run,
  peptide_sequence_column = Stripped.Sequence
) {
  num_proteins_per_sample <- input_table |>
    dplyr::filter(!is.na({{ intensity_column }})) |>
    distinct({{ sample_id_column }}, {{ protein_id_column }}) |>
    group_by({{ sample_id_column }}) |>
    summarise(count = n()) |>
    ungroup()

  num_peptides_per_sample <- input_table |>
    dplyr::filter(!is.na({{ intensity_column }})) |>
    distinct({{ sample_id_column }}, {{ protein_id_column }}, {{ peptide_id_column }}) |>
    group_by({{ sample_id_column }}) |>
    summarise(count = n()) |>
    ungroup()

  combined_counts <- num_proteins_per_sample |>
    mutate(type = "Protein") |>
    bind_rows(num_peptides_per_sample |>
      mutate(type = "Peptide")) |>
    pivot_wider(
      id_cols = {{ sample_id_column }},
      names_from = type,
      values_from = count
    )

  output_plot <- combined_counts |>
    ggplot(aes(reorder({{ sample_id_column }}, Peptide))) +
    geom_point(aes(y = Peptide / 10, shape = "Peptide"), show.legend = TRUE) +
    geom_point(aes(y = Protein, shape = "Protein"), show.legend = TRUE) +
    scale_y_continuous(
      name = "Protein",
      sec.axis = sec_axis(\(x) {
        x * 10
      }, name = "Peptide")
    ) +
    apafTheme() +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      panel.grid.major.x = element_blank()
    ) +
    xlab("Samples") +
    scale_shape_manual(values = c(
      "Peptide" = 1,
      "Protein" = 2
    )) +
    labs(shape = "Category")

  # Strip captured environment to prevent memory bloat
  output_plot$plot_env <- rlang::base_env()
  output_plot
}

# ----------------------------------------------------------------------------
# plotHistogramOfPercentMissingPerIndvidual
# ----------------------------------------------------------------------------
#' @export
plotHistogramOfPercentMissingPerIndvidual <- function(
  percent_missing_table,
  percent_missing_column = percent_missing
) {
  percent_missing_table |>
    ggplot(aes({{ percent_missing_column }})) +
    geom_histogram() +
    apafTheme() +
    xlab("Percent Missing") +
    ylab("Count")
}

# ----------------------------------------------------------------------------
# getOneRlePlotData
# ----------------------------------------------------------------------------
#' @export
getOneRlePlotData <- function(input_matrix) {
  # if(!( length(which( is.na(input_matrix[, 1]) | is.nan(input_matrix[, 1]) | is.infinite(input_matrix[, 1]) )) > 0 )){

  input_matrix[is.infinite(input_matrix) | is.nan(input_matrix)] <- NA

  deviations <- input_matrix - Biobase::rowMedians(input_matrix, na.rm = TRUE)

  stats <- graphics::boxplot(
    deviations,
    outcol = "lightgray",
    cex = 0.1,
    cex.axis = 0.7,
    las = 2,
    outline = FALSE
  )

  rownames(stats$stats) <- c("lower whisker", "lower hinge", "median", "upper hinge", "upper whisker")
  colnames(stats$stats) <- colnames(deviations)

  results <- stats$stats |>
    as.data.frame() |>
    rownames_to_column("Quantiles") |>
    pivot_longer(cols = !contains("Quantiles")) |>
    mutate(Quantiles = factor(Quantiles, levels = rev(c("lower whisker", "lower hinge", "median", "upper hinge", "upper whisker"))))

  # print(head( results) )


  return(results)
  # }
}

# ----------------------------------------------------------------------------
# plotRleQc
# ----------------------------------------------------------------------------
#' @export
plotRleQc <- function(
  input_table,
  x_value = name,
  y_value = value,
  quantiles_column = Quantiles
) {
  rle_results <- input_table |>
    ggplot(aes(x = {{ x_value }}, y = {{ y_value }}, group = {{ quantiles_column }}, col = {{ quantiles_column }})) +
    geom_line() +
    apafTheme() +
    theme(axis.text.x = element_blank()) +
    xlab("Samples") +
    ylab("Relative log expression") +
    labs(col = "Boxplot features")

  # Strip captured environment to prevent memory bloat
  rle_results$plot_env <- rlang::base_env()
  rle_results
}

# ----------------------------------------------------------------------------
# compareUmapComponentsPairs
# ----------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------
#' @export
compareUmapComponentsPairs <- function(input_table, columns = c("V1", "V2", "V3", "V4"), covariate) {
  pm <- umap_data |>
    ggpairs(columns = columns, ggplot2::aes(colour = {{ covariate }}), legend = 1) +
    apafTheme()

  pm
}

# ----------------------------------------------------------------------------
# umap_factor_plot
# ----------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------
#' @export
umap_factor_plot <- function(input_data, header, legend_label, x = V1, y = V2, colour_rule) {
  input_data |>
    mutate(!!sym({{ header }}) := factor(!!sym({{ header }}))) |>
    ggplot(aes({{ x }}, {{ y }}, color = !!sym({{ header }}))) +
    geom_point() +
    scale_colour_manual(
      name = legend_label,
      values = colour_rule,
      breaks = names(colour_rule)
    ) +
    apafTheme()
}

