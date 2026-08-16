# ----------------------------------------------------------------------------
# getCategoricalColourPalette
# ----------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------
#' @export
getCategoricalColourPalette <- function() {
  set1_colour <- brewer.pal(9, "Set1")
  set2_colour <- brewer.pal(8, "Set2")
  set3_colour <- brewer.pal(12, "Set3")
  pastel1_colour <- brewer.pal(9, "Pastel1")
  pastel2_colour <- brewer.pal(8, "Pastel2")
  dark2_colour <- brewer.pal(8, "Dark2")
  accent_colour <- brewer.pal(8, "Accent")
  paired_colour <- brewer.pal(12, "Paired")

  set1_2_3_colour <- c(
    set1_colour, set2_colour, set3_colour,
    pastel1_colour, pastel2_colour, dark2_colour,
    accent_colour, paired_colour
  )

  return(set1_2_3_colour)
}

# ----------------------------------------------------------------------------
# getCmriMachineColour
# ----------------------------------------------------------------------------
#' Standard colours for CMRI mass spectrometers
#'
#' @return A named character vector of colours for machines M01 through M08.
#' @export
getCmriMachineColour <- function() {
  colours <- RColorBrewer::brewer.pal(8, "Set2")
  names(colours) <- sprintf("M%02d", seq_along(colours))
  colours
}

# ----------------------------------------------------------------------------
# getCategoricalColourRules
# ----------------------------------------------------------------------------
#' Build colour rules for categorical sample metadata
#'
#' @param metadata_tbl Sample metadata containing the selected columns.
#' @param metadata_column_labels Display labels corresponding to
#'   `metadata_column_selected`.
#' @param metadata_column_selected Metadata columns available for annotation.
#' @param categorical_columns Columns that should use categorical colours.
#' @param ms_machine_column Optional mass-spectrometer column.
#' @param columns_to_exclude Columns that should not receive colour rules.
#' @param colour_palette Character vector of colours to allocate.
#' @param na_colour Colour used for missing categorical values.
#'
#' @return A named list of named colour vectors.
#' @export
getCategoricalColourRules <- function(
  metadata_tbl,
  metadata_column_labels,
  metadata_column_selected,
  categorical_columns,
  ms_machine_column,
  columns_to_exclude,
  colour_palette = getCategoricalColourPalette(),
  na_colour = "white"
) {
  stopifnot(length(metadata_column_labels) == length(metadata_column_selected))
  label_by_column <- stats::setNames(metadata_column_labels, metadata_column_selected)
  excluded <- intersect(columns_to_exclude, names(metadata_tbl))
  machine_enabled <- length(ms_machine_column) == 1L &&
    !is.na(ms_machine_column) &&
    nzchar(ms_machine_column) &&
    ms_machine_column %in% names(metadata_tbl) &&
    !ms_machine_column %in% excluded

  if (machine_enabled) {
    colour_palette <- setdiff(colour_palette, unname(getCmriMachineColour()))
  }

  categorical_columns <- intersect(
    setdiff(categorical_columns, c(ms_machine_column, excluded)),
    names(metadata_tbl)
  )
  category_values <- lapply(categorical_columns, function(column) {
    sort(unique(as.character(metadata_tbl[[column]][!is.na(metadata_tbl[[column]])])))
  })
  required_colours <- sum(lengths(category_values))
  if (required_colours > length(colour_palette)) {
    colour_palette <- grDevices::colorRampPalette(colour_palette)(required_colours)
  }

  next_colour <- 1L
  colour_rules <- lapply(seq_along(categorical_columns), function(index) {
    column <- categorical_columns[[index]]
    values <- category_values[[index]]
    colours <- colour_palette[seq.int(next_colour, length.out = length(values))]
    next_colour <<- next_colour + length(values)
    names(colours) <- values
    if (anyNA(metadata_tbl[[column]])) {
      colours <- c(colours, "NA" = na_colour)
    }
    colours
  })
  names(colour_rules) <- unname(label_by_column[categorical_columns])

  if (machine_enabled) {
    machine_values <- sort(unique(as.character(metadata_tbl[[ms_machine_column]])))
    machine_colours <- getCmriMachineColour()[machine_values]
    machine_colours[is.na(machine_colours)] <- na_colour
    names(machine_colours) <- machine_values
    colour_rules[[unname(label_by_column[[ms_machine_column]])]] <- machine_colours
  }

  colour_rules
}

# ----------------------------------------------------------------------------
# getOneContinousPalette
# ----------------------------------------------------------------------------
#' @export
getOneContinousPalette <- function(metadata_tbl, column_name, palette_name, na_colour = "white") {
  number_of_values <- metadata_tbl |>
    dplyr::select(all_of(column_name)) |>
    dplyr::filter(!is.na(!!sym(column_name))) |>
    distinct() |>
    arrange(!!sym(column_name)) |>
    pull() |>
    length()

  list_of_names <- metadata_tbl |>
    dplyr::select(all_of(column_name)) |>
    dplyr::filter(!is.na(!!sym(column_name))) |>
    distinct() |>
    arrange(!!sym(column_name)) |>
    pull()

  na_name <- metadata_tbl |>
    dplyr::select(all_of(column_name)) |>
    distinct() |>
    dplyr::filter(is.na(!!sym(column_name))) |>
    pull()

  list_of_colours <- brewer.pal(number_of_values, palette_name)
  names(list_of_colours) <- list_of_names

  if (length(na_name) == 1) {
    new_list_of_colours <- c(list_of_colours, na_colour)
    names(new_list_of_colours) <- c(names(list_of_colours), "NA")
    return(new_list_of_colours)
  }

  return(list_of_colours)
}

# ----------------------------------------------------------------------------
# getContinousColourRules
# ----------------------------------------------------------------------------
#' getContinousColourRules
#' @export
#' @param metadata_tbl,metadata_column_labels,metadata_column_selected,continous_scale_columns,na_colour Runtime inputs used by this function; see the usage section for accepted values.
getContinousColourRules <- function(
  metadata_tbl,
  metadata_column_labels,
  metadata_column_selected,
  continous_scale_columns,
  na_colour = "white"
) {
  metadata_column_labels_copy <- metadata_column_labels
  names(metadata_column_labels_copy) <- metadata_column_selected

  list_of_continuous_colour_palette <- c(
    "Greys", "Blues", "Greens", "Purples", "Reds", "Oranges", "BuGn",
    "BuPu", "GnBu", "OrRd", "PuBu", "PuBuGn", "PuRd",
    "RdPu", "YlGn", "YlGnBu", "YlOrBr", "YlOrRd"
  )

  if (length(continous_scale_columns) > length(list_of_continuous_colour_palette)) {
    list_of_continuous_colour_palette <- rep(list_of_continuous_colour_palette, length.out = length(continous_scale_columns))
  }

  list_of_continous_colour_rules <- purrr::map2(
    continous_scale_columns,
    list_of_continuous_colour_palette[seq_along(continous_scale_columns)],
    \(column, palette_name) {
      getOneContinousPalette(metadata_tbl,
        column,
        palette_name,
        na_colour = na_colour
      )
    }
  )

  names(list_of_continous_colour_rules) <- metadata_column_labels_copy[continous_scale_columns]

  return(list_of_continous_colour_rules)
}

# ----------------------------------------------------------------------------
# getCategoricalAndContinuousColourRules
# ----------------------------------------------------------------------------
#' getCategoricalAndContinuousColourRules
#' @param metadata_tbl This is the table containing sample ID and other columns containing clinical variables / metadata
#' @param metadata_column_labels This is the nice
#' @export
#' @param metadata_column_selected,categorical_columns,continous_scale_columns,ms_machine_column,sample_id_column,columns_to_exclude,na_colour Runtime inputs used by this function; see the usage section for accepted values.
getCategoricalAndContinuousColourRules <- function(
  metadata_tbl,
  metadata_column_labels,
  metadata_column_selected,
  categorical_columns,
  continous_scale_columns,
  ms_machine_column,
  sample_id_column = Run,
  columns_to_exclude,
  na_colour = "white"
) {
  metadata_column_labels_copy <- metadata_column_labels
  names(metadata_column_labels_copy) <- metadata_column_selected

  if (ms_machine_column %in% columns_to_exclude) {
    metadata_column_selected <- setdiff(metadata_column_selected, ms_machine_column)
  }

  cln_meatadata_tbl <- metadata_tbl |>
    column_to_rownames(as_name(enquo(sample_id_column))) |>
    dplyr::select(all_of(c(metadata_column_selected)))

  colour_rules <- getCategoricalColourRules(
    metadata_tbl = cln_meatadata_tbl,
    metadata_column_labels = unname(metadata_column_labels_copy[metadata_column_selected]),
    metadata_column_selected = metadata_column_selected,
    categorical_columns = categorical_columns,
    ms_machine_column = ms_machine_column,
    columns_to_exclude = columns_to_exclude,
    na_colour = na_colour
  )

  print("Add column annotation")
  colnames(cln_meatadata_tbl) <- metadata_column_labels_copy[metadata_column_selected]


  continous_colour_list <- getContinousColourRules(metadata_tbl,
    metadata_column_labels,
    metadata_column_selected,
    continous_scale_columns,
    na_colour = na_colour
  )

  categorical_and_continuous_colour_rules <- c(colour_rules, continous_colour_list)

  columns_to_use <- setdiff(names(categorical_and_continuous_colour_rules), metadata_column_labels_copy[columns_to_exclude])
  categorical_and_continuous_colour_rules_filt <- categorical_and_continuous_colour_rules[columns_to_use]

  return(categorical_and_continuous_colour_rules_filt)
}

# ----------------------------------------------------------------------------
# apafTheme
# ----------------------------------------------------------------------------
#' @title ProCan ggplot2 theme. Rectangle box around each plot.
#' @description Standard ggplot2 theme for ProCan plots
#' @export
apafTheme <- function() {
  theme(
    # Set font family and size
    text = element_text(family = "sans", size = 12),
    # Add rectangular box around the plot
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
    # Add grid lines
    panel.grid.major = element_line(color = "gray", linetype = "dashed"),
    panel.grid.minor = element_blank(),
    # Set plot background color
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA),
    # Set axis line and tick colors
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black"),
    # Set axis label colors and sizes
    axis.text.x = element_text(color = "black", size = 10),
    axis.text.y = element_text(color = "black", size = 10),
    # Set legend title and label colors and sizes
    legend.title = element_text(color = "black", size = 12),
    legend.text = element_text(color = "black", size = 10),
    # Set plot title and subtitle colors and sizes
    plot.title = element_text(color = "black", size = 14),
    plot.subtitle = element_text(color = "black", size = 12),
    # Set plot margin sizes
    plot.margin = unit(c(1, 1, 1, 1), "cm")
  )
}

# ----------------------------------------------------------------------------
# get_color_palette
# ----------------------------------------------------------------------------
#' Generate a color palette
#'
#' @param n Number of colors needed
#' @param base_color Base color to use
#' @return Vector of colors
#' @export
get_color_palette <- function(n, base_color) {
  colorRampPalette(c(base_color, "black"))(n)
}
