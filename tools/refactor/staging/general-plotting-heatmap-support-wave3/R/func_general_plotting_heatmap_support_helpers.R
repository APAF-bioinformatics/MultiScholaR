# ----------------------------------------------------------------------------
# getSamplesCorrelationHeatMap
# ----------------------------------------------------------------------------
#' getSamplesCorrelationHeatMap
#' @description get the
#' @param correlation_matrix Output from the `getSamplesCorrelationMatrix` function
#' @param metadata_tbl This is the table containing sample ID and other columns containing clinical variables / metadata
#' @param is_HEK_column A logical column in the metadata table that indicates if the sample is a HEK sample
#' @param metadata_column_selected A list of column names in string selected from the metadata tbl
#' @param metadata_column_labels A list of column names in string to rename each of the columns selected in the param `metadata_column_selected`
#' @param categorical_columns A vector of string with all the names of the categorical data column  present in the `metadata_tbl` table
#' @param continous_scale_columns  A vector of string with all the names of the continuous data column  present in the `metadata_tbl` table
#' @param ms_machine_column A string of the column name describing the mass spectrometer machine used to analyze each sample
#' @param sample_id_column A string describing the column name of the sample ID column
#' @export
getSamplesCorrelationHeatMap <- function(
  correlation_matrix,
  metadata_tbl,
  is_HEK_column = is_HEK,
  metadata_column_labels,
  metadata_column_selected,
  colour_rules,
  columns_to_exclude,
  sample_id_column = Run,
  use_raster = TRUE,
  raster_device = "CairoPDF",
  heatmap_legend_param = list(title = "Correlation"),
  heatmap_width = ncol(correlation_matrix) * unit(0.05, "cm"),
  heatmap_height = nrow(correlation_matrix) * unit(0.05, "cm")
) {
  names(metadata_column_labels) <- metadata_column_selected

  without_hek_samples <- metadata_tbl |>
    dplyr::filter({{ is_HEK_column }} == FALSE) |>
    pull({{ sample_id_column }})

  correlation_samples_to_use <- intersect(colnames(correlation_matrix), without_hek_samples) |> sort()

  cln_meatadata_orig_col_name <- metadata_tbl |>
    dplyr::filter({{ is_HEK_column }} == FALSE) |>
    dplyr::filter({{ sample_id_column }} %in% correlation_samples_to_use) |>
    arrange({{ sample_id_column }}) |>
    column_to_rownames(as_name(enquo(sample_id_column))) |>
    dplyr::select(all_of(setdiff(metadata_column_selected, columns_to_exclude)))

  columns_to_use <- setdiff(names(colour_rules), metadata_column_labels[columns_to_exclude])
  colour_rules_filt <- colour_rules[columns_to_use]

  print("Add column annotation")
  cln_meatadata_tbl <- cln_meatadata_orig_col_name
  colnames(cln_meatadata_tbl) <- metadata_column_labels[colnames(cln_meatadata_orig_col_name)]

  top_annotation <- HeatmapAnnotation(
    df = cln_meatadata_tbl |>
      dplyr::select(-any_of(metadata_column_labels[columns_to_exclude])),
    col = colour_rules_filt,
    show_legend = FALSE
  )

  print("Add row annotation")
  row_ha <- rowAnnotation(
    df = cln_meatadata_tbl |>
      dplyr::select(-any_of(metadata_column_labels[columns_to_exclude])),
    col = colour_rules_filt,
    show_legend = FALSE
  )

  output_heatmap <- Heatmap(correlation_matrix[correlation_samples_to_use, correlation_samples_to_use],
    name = "Correlation",
    left_annotation = row_ha,
    top_annotation = top_annotation,
    show_row_names = FALSE,
    show_column_names = FALSE,
    use_raster = use_raster,
    raster_device = raster_device,
    row_title_gp = gpar(fontsize = 13.2),
    column_title_gp = gpar(fontsize = 13.2),
    row_names_gp = gpar(fontsize = 12, fontfamily = "sans"),
    column_names_gp = gpar(fontsize = 12),
    heatmap_legend_param = heatmap_legend_param,
    heatmap_height = heatmap_height,
    heatmap_width = heatmap_width
  )

  output_legends <- purrr::map2(
    colour_rules_filt,
    names(colour_rules_filt),
    \(rule, title) {
      Legend(
        labels = names(rule), title = title,
        legend_gp = gpar(fill = rule)
      )
    }
  )

  # output_legends <- packLegend(list = list_of_legends)

  return(list(
    heatmap = output_heatmap,
    legend = output_legends
  ))
}

# ----------------------------------------------------------------------------
# getProteinsHeatMap
# ----------------------------------------------------------------------------
#' @title get protein intensity heatmap
#' @description Generates a heatmap of protein intensities
#' @export
getProteinsHeatMap <- function(
  protein_matrix,
  metadata_tbl,
  is_HEK_column = is_HEK,
  metadata_column_selected,
  metadata_column_labels
  # , categorical_columns
  # , continous_scale_columns
  # , ms_machine_column
  , colour_rules,
  columns_to_exclude,
  core_utilisation_samples = TRUE,
  sort_by_sample_id = TRUE,
  sample_id_column = Run,
  use_raster = TRUE,
  raster_device = "CairoTIFF",
  heatmap_legend_param = list(title = "Intensity")
) {
  metadata_column_labels_copy <- metadata_column_labels
  names(metadata_column_labels_copy) <- metadata_column_selected

  print("Without HEK samples")
  without_hek_samples <- metadata_tbl |>
    dplyr::filter({{ is_HEK_column }} == FALSE) |>
    pull({{ sample_id_column }})

  samples_to_use <- intersect(colnames(protein_matrix), without_hek_samples)

  if (sort_by_sample_id == TRUE) {
    samples_to_use <- samples_to_use |>
      sort()
  }

  cln_meatadata_tbl_orig_col_names <- metadata_tbl |>
    dplyr::filter({{ is_HEK_column }} == FALSE) |>
    dplyr::filter({{ sample_id_column }} %in% samples_to_use) |>
    arrange({{ sample_id_column }}) |>
    dplyr::select({{ sample_id_column }}, all_of(setdiff(metadata_column_selected, columns_to_exclude))) |>
    distinct() |>
    column_to_rownames(as_label(enquo(sample_id_column)))

  print("Add column annotation")
  cln_meatadata_tbl <- cln_meatadata_tbl_orig_col_names[
    samples_to_use,
    setdiff(metadata_column_selected, columns_to_exclude)
  ]
  colnames(cln_meatadata_tbl) <- metadata_column_labels_copy[setdiff(metadata_column_selected, columns_to_exclude)]
  colour_rules_filt <- colour_rules[metadata_column_labels_copy[setdiff(metadata_column_selected, columns_to_exclude)]]

  # print(colour_rules)
  # print(metadata_column_labels_copy[setdiff(metadata_column_selected, columns_to_exclude)] )
  # print(colour_rules)
  # print(colour_rules_filt)

  # print(colour_rules_filt)

  print("Set top located annotation")
  top_annotation <- HeatmapAnnotation(
    df = cln_meatadata_tbl,
    col = colour_rules_filt,
    show_legend = FALSE,
    annotation_name_side = "left"
  )

  print("Print Heatmap")

  heatmap <- Heatmap(protein_matrix[, samples_to_use],
    name = "Intensity",
    top_annotation = top_annotation,
    show_row_names = FALSE,
    show_column_names = FALSE,
    use_raster = use_raster,
    raster_device = raster_device,
    row_title_gp = gpar(fontsize = 13.2),
    column_title_gp = gpar(fontsize = 13.2),
    row_names_gp = gpar(fontsize = 12, fontfamily = "sans"),
    column_names_gp = gpar(fontsize = 12),
    heatmap_legend_param = heatmap_legend_param,
    core_utilisation_columns = core_utilisation_samples
  )

  output_legends <- purrr::map2(
    colour_rules_filt,
    names(colour_rules_filt),
    \(rule, title) {
      Legend(
        labels = names(rule), title = title,
        legend_gp = gpar(fill = rule)
      )
    }
  )

  return(list(
    heatmap = heatmap,
    legend = output_legends
  ))
}

# ----------------------------------------------------------------------------
# save_heatmap_products
# ----------------------------------------------------------------------------
#' Save Heatmap Products
#'
#' @description Saves a heatmap plot to PNG and PDF, writes generation parameters to a Markdown file,
#' and saves cluster assignments to a CSV file.
#'
#' @param heatmap_obj A ComplexHeatmap object.
#' @param row_clusters A named vector of cluster assignments (names are molecule IDs, values are cluster IDs). Can be NULL.
#' @param params_list A list of parameters used to generate the heatmap (for documentation).
#' @param output_dir The directory where the "Heatmap" subdirectory will be created/used.
#' @param file_prefix A string to prefix the filenames with (e.g., "prot_contrast1").
#'
#' @importFrom grDevices png pdf dev.off
#' @importFrom ComplexHeatmap draw
#' @importFrom utils write.csv
#' @export
save_heatmap_products <- function(heatmap_obj, row_clusters, params_list, output_dir, file_prefix) {
  logger::log_info("--- Entering save_heatmap_products ---")

  # 1. Create directory
  target_dir <- file.path(output_dir, "Heatmap")
  if (!dir.exists(target_dir)) {
    dir.create(target_dir, recursive = TRUE, showWarnings = FALSE)
    logger::log_info(sprintf("   Created directory: %s", target_dir))
  }

  # 2. Save Heatmap Images (PNG & PDF)
  tryCatch(
    {
      # PNG
      png_path <- file.path(target_dir, paste0(file_prefix, "_heatmap.png"))
      grDevices::png(png_path, width = 10, height = 8, units = "in", res = 300)
      ComplexHeatmap::draw(heatmap_obj)
      grDevices::dev.off()
      logger::log_info(sprintf("   Saved heatmap PNG: %s", basename(png_path)))

      # PDF
      pdf_path <- file.path(target_dir, paste0(file_prefix, "_heatmap.pdf"))
      grDevices::pdf(pdf_path, width = 10, height = 8)
      ComplexHeatmap::draw(heatmap_obj)
      grDevices::dev.off()
      logger::log_info(sprintf("   Saved heatmap PDF: %s", basename(pdf_path)))
    },
    error = function(e) {
      logger::log_error(sprintf("   Error saving heatmap images: %s", e$message))
    }
  )

  # 3. Save Methods Markdown
  tryCatch(
    {
      md_path <- file.path(target_dir, paste0(file_prefix, "_heatmap_methods.md"))

      # Construct Markdown content
      md_lines <- c(
        "# Heatmap Generation Parameters",
        "",
        "The following parameters were used to generate this heatmap:",
        "",
        "| Parameter | Value |",
        "| :--- | :--- |"
      )

      for (param_name in names(params_list)) {
        val <- params_list[[param_name]]
        # Convert logicals to text
        if (is.logical(val)) val <- ifelse(val, "TRUE", "FALSE")
        # Handle vectors (comma separated)
        if (length(val) > 1) val <- paste(val, collapse = ", ")
        # Escape pipes for markdown table
        val <- gsub("|", "\\|", as.character(val), fixed = TRUE)

        md_lines <- c(md_lines, sprintf("| %s | %s |", param_name, val))
      }

      md_lines <- c(md_lines, "", sprintf("Generated on: %s", Sys.time()))

      writeLines(md_lines, md_path)
      logger::log_info(sprintf("   Saved methods MD: %s", basename(md_path)))
    },
    error = function(e) {
      logger::log_error(sprintf("   Error saving methods markdown: %s", e$message))
    }
  )

  # 4. Save Cluster CSV
  if (!is.null(row_clusters)) {
    tryCatch(
      {
        csv_path <- file.path(target_dir, paste0(file_prefix, "_clusters.csv"))

        cluster_df <- data.frame(
          molecule_id = names(row_clusters),
          cluster_id = as.vector(row_clusters),
          stringsAsFactors = FALSE
        )

        utils::write.csv(cluster_df, csv_path, row.names = FALSE)
        logger::log_info(sprintf("   Saved clusters CSV: %s", basename(csv_path)))
      },
      error = function(e) {
        logger::log_error(sprintf("   Error saving clusters CSV: %s", e$message))
      }
    )
  } else {
    logger::log_info("   No cluster information to save.")
  }

  logger::log_info("--- Exiting save_heatmap_products ---")
}

