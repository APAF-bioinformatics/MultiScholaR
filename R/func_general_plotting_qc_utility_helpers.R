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
#' @export
#' @param input_table,intensity_column,protein_id_column,peptide_id_column,sample_id_column,peptide_sequence_column Runtime inputs used by this function; see the usage section for accepted values.
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
  requireOptionalPackage("Biobase", "RLE plot data")
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
  requireOptionalPackage("GGally", "UMAP component pair plots")
  pm <- input_table |>
    GGally::ggpairs(columns = columns, ggplot2::aes(colour = {{ covariate }}), legend = 1) +
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

# ----------------------------------------------------------------------------
# saveSelfContainedWidget (KIT-036)
# ----------------------------------------------------------------------------
#' Save an htmlwidget as one self-contained HTML file without pandoc
#'
#' `htmlwidgets::saveWidget(selfcontained = TRUE)` hands the rendered page to pandoc to embed the widget's
#' script and style dependencies; on a page that already carries a large data payload (the DA volcano page is
#' about 15 MB for 300 samples) pandoc re-parses the whole document and peaks at about 3 GiB resident for
#' several seconds, in whichever process saves the page (KIT-017 post-hoc audit PH-01). The dependencies
#' themselves are under 1 MB of plain script and style files, so this writer inlines them in R: each script as
#' a `<script>` element, each stylesheet as a `<style>` element, and the widget's own head and body as
#' `htmltools` renders them. The page is one file that opens without network access, as before.
#'
#' A dependency with attachments, or a file whose contents would close its own element (`</script` or
#' `</style`), cannot be inlined this way; the writer then falls back to `htmlwidgets::saveWidget`, as it does
#' when `getOption("multischolar.widget_pandoc")` is `TRUE` (the pandoc path kept for one release).
#' @param widget An htmlwidget.
#' @param file Output path.
#' @param title Page title (the widget's class by default).
#' @param background Body background colour, as `saveWidget` sets it.
#' @return `file`, invisibly, with attribute `writer` (`"inline"` or `"pandoc"`).
#' @noRd
saveSelfContainedWidget <- function(widget, file, title = NULL, background = "white") {
  # Anything that is not an htmlwidget takes the writer it always had (which refuses it, or a test's stand-in).
  if (!inherits(widget, "htmlwidget") || isTRUE(getOption("multischolar.widget_pandoc", FALSE))) {
    htmlwidgets::saveWidget(widget, file, selfcontained = TRUE)
    return(invisible(structure(file, writer = "pandoc")))
  }
  rendered <- htmltools::renderTags(htmltools::as.tags(widget, standalone = TRUE))
  dependencies <- htmltools::resolveDependencies(rendered$dependencies)
  read_file <- function(path) {
    text <- readChar(path, file.size(path), useBytes = TRUE)
    Encoding(text) <- "UTF-8"
    text
  }
  inline <- character()
  for (dependency in dependencies) {
    root <- dependency$src$file
    if (is.null(root) || length(dependency$attachment)) return(saveSelfContainedWidgetFallback(widget, file, title, background))
    if (!is.null(dependency$package)) root <- system.file(root, package = dependency$package)
    for (name in dependency$stylesheet) {
      css <- read_file(file.path(root, name))
      if (grepl("</style", css, fixed = TRUE, useBytes = TRUE)) return(saveSelfContainedWidgetFallback(widget, file, title, background))
      inline <- c(inline, "<style>", css, "</style>")
    }
    for (name in dependency$script) {
      path <- if (is.list(name)) name$src else name
      js <- read_file(file.path(root, path))
      if (grepl("</script", js, fixed = TRUE, useBytes = TRUE)) return(saveSelfContainedWidgetFallback(widget, file, title, background))
      inline <- c(inline, "<script>", js, "</script>")
    }
    if (length(dependency$head)) inline <- c(inline, dependency$head)
  }
  page <- c("<!DOCTYPE html>", "<html>", "<head>", "<meta charset=\"utf-8\"/>",
    paste0("<title>", htmltools::htmlEscape(title %||% class(widget)[[1L]]), "</title>"),
    inline, as.character(rendered$head), "</head>",
    paste0("<body style=\"background-color: ", htmltools::htmlEscape(background), ";\">"),
    as.character(rendered$html), "</body>", "</html>")
  connection <- file(file, open = "wb")
  on.exit(close(connection), add = TRUE)
  writeLines(enc2utf8(page), connection, useBytes = TRUE)
  invisible(structure(file, writer = "inline"))
}

#' @noRd
saveSelfContainedWidgetFallback <- function(widget, file, title, background) {
  # The call the page saves made before KIT-036, so a test's stand-in for saveWidget sees the same arguments.
  htmlwidgets::saveWidget(widget, file, selfcontained = TRUE)
  invisible(structure(file, writer = "pandoc"))
}
