# ----------------------------------------------------------------------------
# plotOneVolcano
# ----------------------------------------------------------------------------
## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#' Draw the volcano plot, used in publication graphs
#' @param input_data Input data with the `log_q_value_column`, the `log_fc_column`, the `points_type_label` and the `points_color` columns.
#' @param log_q_value_column The name of the column representing the log q-value.
#' @param log_fc_column The name of the column representing the log fold-change.
#' @param points_type_label A column in input table with the type of points based on log fold-change and q-value (e.g. "Not sig., logFC >= 1" = "orange" , "Sig., logFC >= 1" = "purple" , "Sig., logFC < 1" = "blue" , "Not sig." )
#' @param points_color A column in input table with the colour of the points corresponding to each type of points (e.g. orange, purple, blue black, )
#' @param q_val_thresh A numerical value specifying the q-value threshold for statistically significant proteins.
#' @param log2FC_thresh A numerical value specifying the log fold-change threshold to draw a vertical line
#' @param formula_string The formula string used in the facet_grid command for the ggplot scatter plot.
#' @export
plotOneVolcano <- function(input_data, input_title,
                           log_q_value_column = lqm,
                           log_fc_column = logFC,
                           points_type_label = label,
                           points_color = colour,
                           q_val_thresh = 0.05,
                           log2FC_thresh = 1) {
  colour_tbl <- input_data |>
    distinct({{ points_type_label }}, {{ points_color }})

  # print(colour_tbl)

  colour_map <- colour_tbl |>
    dplyr::pull({{ points_color }}) |>
    as.vector()

  names(colour_map) <- colour_tbl |>
    dplyr::pull({{ points_type_label }})

  avail_labels <- input_data |>
    distinct({{ points_type_label }}) |>
    dplyr::pull({{ points_type_label }})

  avail_colours <- colour_map[avail_labels]

  # print(avail_labels)
  # print(avail_colours)

  volcano_plot <- input_data |>
    ggplot(aes(
      y = {{ log_q_value_column }},
      x = {{ log_fc_column }}
    )) +
    geom_point(aes(col = label)) +
    scale_colour_manual(values = avail_colours) +
    geom_hline(yintercept = -log10(q_val_thresh)) +
    theme_bw() +
    xlab(expression(Log[2](`fold-change`))) +
    ylab(expression(-log[10](`q-value`))) +
    labs(title = input_title) + # Remove legend title
    theme(legend.title = element_blank()) +
    # theme(legend.position = "none")  +
    theme(axis.text.x = element_text(size = 13)) +
    theme(axis.text.y = element_text(size = 13)) +
    theme(axis.title.x = element_text(size = 12)) +
    theme(axis.title.y = element_text(size = 12)) +
    theme(plot.title = element_text(size = 12)) +
    theme(legend.text = element_text(size = 12)) # +
  # theme(legend.title = element_text(size = 12))

  if (!is.na(log2FC_thresh)) {
    volcano_plot <- volcano_plot +
      geom_vline(xintercept = log2FC_thresh, colour = "black", linewidth = 0.2) +
      geom_vline(xintercept = -log2FC_thresh, colour = "black", linewidth = 0.2)
  }

  volcano_plot
}

# ----------------------------------------------------------------------------
# plotOneVolcanoNoVerticalLines
# ----------------------------------------------------------------------------
## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#' Draw the volcano plot, used in publication graphs
#' @param input_data Input data with the `log_q_value_column`, the `log_fc_column`, the `points_type_label` and the `points_color` columns.
#' @param log_q_value_column The name of the column representing the log q-value.
#' @param log_fc_column The name of the column representing the log fold-change.
#' @param points_type_label A column in input table with the type of points based on log fold-change and q-value (e.g. "Not sig., logFC >= 1" = "orange" , "Sig., logFC >= 1" = "purple" , "Sig., logFC < 1" = "blue" , "Not sig." )
#' @param points_color A column in input table with the colour of the points corresponding to each type of points (e.g. orange, purple, blue black, )
#' @param q_val_thresh A numerical value specifying the q-value threshold for statistically significant proteins.
#' @param gene_name The column representing the gene name
#' @param formula_string The formula string used in the facet_grid command for the ggplot scatter plot.
#' @export
plotOneVolcanoNoVerticalLines <- function(input_data, input_title,
                                          log_q_value_column = lqm,
                                          log_fc_column = logFC,
                                          points_type_label = label,
                                          points_color = colour,
                                          q_val_thresh = 0.05) {
  colour_tbl <- input_data |>
    distinct({{ points_type_label }}, {{ points_color }})

  colour_map <- colour_tbl |>
    dplyr::pull({{ points_color }}) |>
    as.vector()

  names(colour_map) <- colour_tbl |>
    dplyr::pull({{ points_type_label }})

  avail_labels <- input_data |>
    distinct({{ points_type_label }}) |>
    dplyr::pull({{ points_type_label }})

  avail_colours <- colour_map[avail_labels]

  volcano_plot <- input_data |>
    ggplot(aes(
      y = {{ log_q_value_column }},
      x = {{ log_fc_column }},
      col = {{ points_type_label }}
    )) +
    geom_point()

  volcano_plot <- volcano_plot +
    scale_colour_manual(values = avail_colours) +
    # geom_vline(xintercept = 1, colour = "black", size = 0.2) +
    # geom_vline(xintercept = -1, colour = "black", size = 0.2) +
    geom_hline(yintercept = -log10(q_val_thresh)) +
    theme_bw() +
    xlab(expression(Log[2](`fold-change`))) +
    ylab(expression(-log[10](FDR))) +
    labs(title = input_title) + # Remove legend title
    theme(legend.title = element_blank()) +
    # theme(legend.position = "none")  +
    theme(axis.text.x = element_text(size = 13)) +
    theme(axis.text.y = element_text(size = 13)) +
    theme(axis.title.x = element_text(size = 12)) +
    theme(axis.title.y = element_text(size = 12)) +
    theme(plot.title = element_text(size = 12)) +
    theme(legend.text = element_text(size = 12)) # +
  # theme(legend.title = element_text(size = 12))

  volcano_plot
}

# ----------------------------------------------------------------------------
# printOneVolcanoPlotWithProteinLabel
# ----------------------------------------------------------------------------
#'  This function creates a volcano plot with protein labels
#' @param input_table The input table to be used for the volcano plot, contains the protein_id_column, fdr_column and log2FC_column
#' @param uniprot_table The uniprot table to be used for the volcano plot, contains the uniprot_protein_id_column and gene_name_column
#' @param protein_id_column The column name in the input_table that contains the protein ids (tidyverse format)
#' @param uniprot_protein_id_column The column name in the uniprot_table that contains the uniprot protein ids (tidyverse format)
#' @param gene_name_column The column name in the uniprot_table that contains the gene names (tidyverse format)
#' @param number_of_genes Increasing P-value rank for the number of proteins to display on the volcano plot, default is 100`
#' @param fdr_threshold The FDR threshold to use for the volcano plot, default is 0.05
#' @param fdr_column The column name in the input_table that contains the FDR values (tidyverse format)
#' @param log2FC_column The column name in the input_table that contains the log2FC values (tidyverse format)
#' @param input_title The title to use for the volcano plot
#' @param max.overlaps The maximum number of overlaps to allow for the protein labels using ggrepel (default is 20)
#' @return A ggplot object with the volcano plot and protein labels
printOneVolcanoPlotWithProteinLabel <- function(
  input_table,
  uniprot_table,
  protein_id_column = Protein.Ids,
  uniprot_protein_id_column = Entry,
  gene_name_column = gene_name,
  number_of_genes = 100,
  fdr_threshold = 0.05,
  fdr_column = fdr_qvalue,
  log2FC_column = log2FC,
  input_title = "Proteomics",
  include_protein_label = TRUE,
  max.overlaps = 20
) {
  proteomics_volcano_tbl <- prepareDataForVolcanoPlot(input_table,
    protein_id_column = {{ protein_id_column }},
    uniprot_table = uniprot_table,
    uniprot_protein_id_column = {{ uniprot_protein_id_column }},
    gene_name_column = {{ gene_name_column }},
    number_of_genes = number_of_genes,
    fdr_threshold = fdr_threshold,
    fdr_column = {{ fdr_column }},
    log2FC_column = {{ log2FC_column }}
  )

  proteomics_volcano_plot <- plotOneVolcanoNoVerticalLines(proteomics_volcano_tbl,
    input_title = input_title,
    log_q_value_column = lqm,
    log_fc_column = log2FC,
    points_type_label = label,
    points_color = colour,
    q_val_thresh = fdr_threshold
  )

  proteomics_volcano_plot_with_proteins_label <- proteomics_volcano_plot

  if (include_protein_label == TRUE) {
    proteomics_volcano_plot_with_proteins_label <- proteomics_volcano_plot +
      ggrepel::geom_text_repel(aes(label = gene_name_significant),
        show.legend = FALSE,
        max.overlaps = max.overlaps
      )
  }

  proteomics_volcano_plot_with_proteins_label
}

# ----------------------------------------------------------------------------
# getGlimmaVolcanoProteomics
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
# getGlimmaVolcanoProteomics
# ----------------------------------------------------------------------------
getGlimmaVolcanoProteomics <- function(
  volcano_plot_tab,
  uniprot_column = best_uniprot_acc,
  gene_name_column = gene_name,
  display_columns = c("PROTEIN_NAMES"),
  additional_annotations = NULL,
  additional_annotations_join_column = NULL,
  counts_tbl = NULL,
  groups = NULL,
  output_dir,
  fdr_column = "fdr_qvalue",
  log2fc_column = "log2FC",
  da_q_val_thresh = 0.05,
  contrast_name = "Contrast",
  ...
) {
  
  # 1. Clean and prepare basic annotation
  volcano_plot_tab_cln <- volcano_plot_tab |>
    dplyr::select(
      {{ uniprot_column }},
      {{ gene_name_column }}, any_of(display_columns)
    ) |>
    distinct()

  if (!is.null(additional_annotations) &
    !is.null(additional_annotations_join_column)) {
    volcano_plot_tab_cln <- volcano_plot_tab_cln |>
      left_join(additional_annotations,
        by = join_by({{ uniprot_column }} == {{ additional_annotations_join_column }})
      ) |>
      dplyr::select(
        {{ uniprot_column }},
        {{ gene_name_column }},
        any_of(display_columns)
      )
  }

  anno_tbl <- volcano_plot_tab_cln |>
    mutate(gene_name = case_when(
      is.na({{ gene_name_column }}) | {{ gene_name_column }} == "" ~ as.character({{ uniprot_column }}),
      TRUE ~ as.character({{ gene_name_column }})
    ))

  # CRITICAL FIX: Glimma Best Practice: Replace dots with underscores in column names
  # and coerce all to character, replacing NAs with empty strings.
  anno_tbl <- anno_tbl |>
    dplyr::rename_with(~ gsub("\\.", "_", .)) |>
    dplyr::mutate(across(everything(), as.character)) |>
    dplyr::mutate(across(everything(), ~ tidyr::replace_na(., ""))) |>
    as.data.frame()

  # 2. Extract plot data (x, y) from the main table
  # Since this function is called from outputDaAnalysisResults, volcano_plot_tab already has the metrics.
  plot_data <- volcano_plot_tab |>
    dplyr::mutate(
      FDR = !!sym(fdr_column),
      FDR = ifelse(FDR == 0, min(FDR[FDR > 0], na.rm = TRUE) * 0.1, FDR),
      negLog10FDR = -log10(FDR),
      logFC = !!sym(log2fc_column)
    ) |>
    dplyr::filter(!is.infinite(negLog10FDR), !is.na(negLog10FDR), !is.na(logFC))

  if (nrow(plot_data) == 0) return(NULL)
  
  # Add status for coloring
  plot_data <- plot_data |>
    dplyr::mutate(
      Status = case_when(
        logFC >= 1 & FDR < as.double(da_q_val_thresh) ~ 1,
        logFC <= -1 & FDR < as.double(da_q_val_thresh) ~ -1,
        TRUE ~ 0
      )
    )

  # 3. Match anno_tbl rows to plot_data exactly
  # Retrieve the actual column name for uniprot_column to use in matching
  id_col_str <- rlang::as_label(rlang::enquo(uniprot_column))
  id_col_cln <- gsub("\\.", "_", id_col_str)
  
  match_idx <- match(plot_data[[id_col_str]], anno_tbl[[id_col_cln]])
  anno_tbl <- anno_tbl[match_idx, , drop = FALSE]
  
  # CRITICAL: Use the same ID for rownames as will be used for counts
  rownames(anno_tbl) <- as.character(plot_data[[id_col_str]])

  # 4. Prepare Counts matrix
  if (!is.null(counts_tbl)) {
    matching_col <- NULL
    
    # Determine which column in plot_data maps to rownames(counts_tbl)
    id_col_str <- rlang::as_label(rlang::enquo(uniprot_column))
    potential_cols <- unique(c(id_col_str, "Protein.Ids", "Protein.ID", "uniprot_acc", "Entry", names(plot_data)))
    
    for (col in potential_cols) {
      if (col %in% names(plot_data)) {
        if (sum(plot_data[[col]] %in% rownames(counts_tbl)) > 0) {
          matching_col <- col
          break
        }
      }
    }
    
    if (!is.null(matching_col)) {
      # Ensure plot_data only includes items present in the counts table
      has_counts <- plot_data[[matching_col]] %in% rownames(counts_tbl)
      
      if (sum(has_counts) < nrow(plot_data)) {
         plot_data <- plot_data[has_counts, , drop = FALSE]
         anno_tbl <- anno_tbl[has_counts, , drop = FALSE]
         gene_names <- anno_tbl$gene_name
      }
      
      counts_filtered <- counts_tbl[plot_data[[matching_col]], , drop = FALSE]
      rownames(counts_filtered) <- gene_names
    } else if (nrow(counts_tbl) == nrow(plot_data)) {
      counts_filtered <- counts_tbl
      rownames(counts_filtered) <- gene_names
    } else {
      counts_filtered <- NULL
    }
  } else {
    counts_filtered <- NULL
  }

  # 5. Generate glimmaXY Plot with temporary file workaround
  temp_html_name <- paste0("volcano_", gsub("[^A-Za-z0-9]", "_", contrast_name), "_temp.html")
  
  Glimma::glimmaXY(
    x = plot_data$logFC,
    y = plot_data$negLog10FDR,
    xlab = "logFC",
    ylab = "negLog10FDR",
    counts = counts_filtered,
    groups = groups,
    status = plot_data$Status,
    anno = anno_tbl,
    display.columns = colnames(anno_tbl),
    transform.counts = "none",
    status.cols = c("#1052bd", "silver", "#cc212f"),
    sample.cols = if (!is.null(groups)) {
      unique_groups <- unique(groups)
      group_colors <- stats::setNames(
        grDevices::hcl.colors(length(unique_groups), "Set2"), 
        unique_groups
      )
      group_colors[groups]
    } else if (!is.null(counts_filtered)) {
      rep("#1f77b4", ncol(counts_filtered))
    } else {
      NULL
    },
    main = paste("Interactive Volcano Plot:", contrast_name),
    html = temp_html_name
  )

  # 6. Post-process to inject CSS and move to final destination
  final_html_path <- file.path(output_dir, paste0(contrast_name, ".html"))
  
  if (file.exists(temp_html_name)) {
    html_lines <- readLines(temp_html_name)
    css_fix <- "<style> .dataTables_wrapper { overflow-x: auto !important; } .dataTables_scroll { display: table !important; width: auto !important; min-width: 100% !important; } .dataTables_scrollHead, .dataTables_scrollBody, .dataTables_scrollHeadInner, .dataTables_scrollHeadInner > table, .dataTables_scrollBody > table { display: contents !important; } .dataTables_scrollBody thead { display: none !important; } .dataTable th, .dataTable td { white-space: nowrap !important; } </style></head>"
    html_lines <- sub("</head>", css_fix, html_lines)
    writeLines(html_lines, temp_html_name)

    file.rename(temp_html_name, final_html_path)
  }
}

# ----------------------------------------------------------------------------
# getGlimmaVolcanoProteomicsWidget
# ----------------------------------------------------------------------------
#' @export
getGlimmaVolcanoProteomicsWidget <- function(
  r_obj,
  coef,
  volcano_plot_tab,
  uniprot_column = best_uniprot_acc,
  gene_name_column = gene_name,
  display_columns = c("PROTEIN_NAMES"),
  additional_annotations = NULL,
  additional_annotations_join_column = NULL,
  counts_tbl = NULL,
  groups = NULL,
  ...
) {
  message("--- Entering getGlimmaVolcanoProteomicsWidget ---")
  message(sprintf("   getGlimmaVolcanoProteomicsWidget Arg: coef = %d", coef))
  message(sprintf("   getGlimmaVolcanoProteomicsWidget Arg: r_obj class = %s", class(r_obj)))
  message(sprintf("   getGlimmaVolcanoProteomicsWidget Arg: volcano_plot_tab dims = %d rows, %d cols", nrow(volcano_plot_tab), ncol(volcano_plot_tab)))
  
  if (coef <= ncol(r_obj$coefficients)) {
    message("   getGlimmaVolcanoProteomicsWidget Step: Coefficient validation passed...")

    # Extract best_uniprot_acc from rownames
    best_uniprot_acc <- str_split(rownames(r_obj@.Data[[1]]), " |:") |>
      purrr::map_chr(1)

    volcano_plot_tab_cln <- volcano_plot_tab |>
      dplyr::select(
        {{ uniprot_column }},
        {{ gene_name_column }}, any_of(display_columns)
      ) |>
      distinct()

    if (!is.null(additional_annotations) &
      !is.null(additional_annotations_join_column)) {
      volcano_plot_tab_cln <- volcano_plot_tab_cln |>
        left_join(additional_annotations,
          by = join_by({{ uniprot_column }} == {{ additional_annotations_join_column }})
        ) |>
        dplyr::select(
          {{ uniprot_column }},
          {{ gene_name_column }},
          any_of(display_columns)
        )
    }

    # [OK] REFACTORED: Use centralized UniProt normalization
    best_uniprot_acc_base <- normalizeUniprotAccession(best_uniprot_acc, remove_isoform = TRUE)

    anno_tbl <- data.frame(
      uniprot_acc = rownames(r_obj@.Data[[1]]),
      temp_column = best_uniprot_acc,
      temp_column_base = best_uniprot_acc_base
    ) |>
      dplyr::rename({{ uniprot_column }} := temp_column)

    # First try exact match
    anno_tbl <- anno_tbl |>
      left_join(volcano_plot_tab_cln,
        by = join_by({{ uniprot_column }} == {{ uniprot_column }})
      )

    # CRITICAL FIX: For entries with no gene_name (NA or empty string), try base accession match
    missing_gene_mask <- is.na(anno_tbl$gene_name) | anno_tbl$gene_name == ""

    if (any(missing_gene_mask)) {
      message(sprintf("   getGlimmaVolcanoProteomicsWidget: %d proteins missing gene names, trying base accession match...", sum(missing_gene_mask)))

      base_match_tbl <- data.frame(temp_column_base = anno_tbl$temp_column_base[missing_gene_mask]) |>
        left_join(
          volcano_plot_tab_cln |>
            dplyr::mutate(temp_base = normalizeUniprotAccession({{ uniprot_column }}, remove_isoform = TRUE)) |>
            dplyr::select(temp_base, gene_name_alt = {{ gene_name_column }}) |>
            dplyr::distinct(temp_base, .keep_all = TRUE),
          by = join_by(temp_column_base == temp_base)
        )

      found_alt <- !is.na(base_match_tbl$gene_name_alt) & base_match_tbl$gene_name_alt != ""
      if (any(found_alt)) {
        which_missing <- which(missing_gene_mask)
        anno_tbl$gene_name[which_missing[found_alt]] <- base_match_tbl$gene_name_alt[found_alt]
      }
    }

    # Remove temporary base column
    anno_tbl <- anno_tbl |> dplyr::select(-temp_column_base)

    # Fallback - use accession as gene name if still NA or empty
    anno_tbl <- anno_tbl |>
      mutate(gene_name = case_when(
        is.na(gene_name) | gene_name == "" ~ {{ uniprot_column }},
        TRUE ~ gene_name
      ))

    # Glimma Best Practice: Replace dots with underscores in column names
    colnames(anno_tbl) <- gsub("\\.", "_", colnames(anno_tbl))
    
    # Glimma Best Practice: Replace NAs with empty strings and coerce to character
    anno_tbl <- anno_tbl |>
      mutate(across(everything(), ~ as.character(ifelse(is.na(.), "", .))))

    # CRITICAL: Use the original IDs for rownames to match counts matrix
    rownames(anno_tbl) <- as.character(rownames(r_obj@.Data[[1]]))

    # Update rownames in r_obj to gene names for plot labels if desired,
    # but Glimma v2 glimmaVolcano uses anno for everything.
    # We will keep r_obj names as is to match anno rownames.
    # rownames(r_obj@.Data[[1]]) <- gene_names # (DON'T DO THIS IF IT BREAKS MATCHING)

    # Transform p-values to q-values
    # Check if qvalue calculation is possible
    valid_p <- r_obj$p.value[, coef][!is.na(r_obj$p.value[, coef])]
    if (length(valid_p) > 0 && any(valid_p < 1)) {
       tryCatch({
         r_obj$p.value[, coef] <- qvalue(r_obj$p.value[, coef])$qvalues
       }, error = function(e) {
         message(paste("   getGlimmaVolcanoProteomicsWidget: qvalue failed, using raw p-values:", e$message))
       })
    }

    # Prepare status
    status_result <- decideTests(r_obj, adjust.method = "none")
    if (any(is.na(status_result))) {
      status_result[is.na(status_result)] <- 0
    }

    # Call glimmaVolcano
    widget <- glimmaVolcano(r_obj,
      coef = coef,
      counts = counts_tbl,
      groups = groups,
      anno = anno_tbl,
      display.columns = colnames(anno_tbl),
      status = status_result,
      p.adj.method = "none",
      transform.counts = "none",
      sample.cols = if (!is.null(groups)) {
        unique_groups <- unique(groups)
        group_colors <- stats::setNames(
          grDevices::hcl.colors(length(unique_groups), "Set2"), 
          unique_groups
        )
        group_colors[groups]
      } else if (!is.null(counts_tbl)) {
        rep("#1f77b4", ncol(counts_tbl))
      } else {
        NULL
      },
      ...
    )

    message("--- Exiting getGlimmaVolcanoProteomicsWidget ---")
    return(widget)
  }
}

# ----------------------------------------------------------------------------
# getGlimmaVolcanoPhosphoproteomics
# ----------------------------------------------------------------------------
getGlimmaVolcanoPhosphoproteomics <- function(
  r_obj,
  coef,
  volcano_plot_tab,
  sites_id_column = sites_id,
  sites_id_display_column = sites_id_short,
  display_columns = c("sequence", "PROTEIN_NAMES"),
  additional_annotations = NULL,
  additional_annotations_join_column = NULL,
  counts_tbl = NULL,
  output_dir
) {
  if (coef <= ncol(r_obj$coefficients)) {
    volcano_plot_tab_cln <- volcano_plot_tab |>
      dplyr::distinct(
        {{ sites_id_column }},
        {{ sites_id_display_column }},
        any_of(display_columns)
      )

    if (!is.null(additional_annotations) &
      !is.null(additional_annotations_join_column)) {
      volcano_plot_tab_cln <- volcano_plot_tab_cln |>
        left_join(additional_annotations,
          by = join_by({{ sites_id_column }} == {{ additional_annotations_join_column }})
        ) |>
        dplyr::select(any_of(display_columns))
    }

    anno_tbl <- data.frame(sites_id = rownames(r_obj@.Data[[1]])) |>
      left_join(volcano_plot_tab_cln,
        by = join_by(sites_id == {{ sites_id_column }})
      )

    sites_id_short_list <- anno_tbl |>
      dplyr::pull(sites_id_short)

    rownames(r_obj@.Data[[1]]) <- sites_id_short_list

    # coef <- seq_len( ncol(r_obj$coefficients))[1]

    r_obj$p.value[, coef] <- qvalue(r_obj$p.value[, coef])$qvalues


    htmlwidgets::saveWidget(
      widget = glimmaVolcano(r_obj,
        coef = coef,
        counts = counts_tbl,
        anno = anno_tbl,
        display.columns = display_columns,
        p.adj.method = "none",
        transform.counts = "none"
      ) # the plotly object
      , file = file.path(
        output_dir,
        paste0(colnames(r_obj$coefficients)[coef], ".html")
      ) # the path & file name
      , selfcontained = TRUE # creates a single html file
    )
  }
}

# ----------------------------------------------------------------------------
# plotVolcano
# ----------------------------------------------------------------------------
#' Draw the volcano plot.
#' @param selected_data A table that is generated by running the function \code{\link{get_significant_data}}.
#' @param log_q_value_column The name of the column representing the log q-value.
#' @param log_fc_column The name of the column representing the log fold-change.
#' @param q_val_thresh A numerical value specifying the q-value threshold for statistically significant proteins.
#' @param formula_string The formula string used in the facet_grid command for the ggplot scatter plot.
#' @export
plotVolcano <- function(
  selected_data,
  log_q_value_column = lqm,
  log_fc_column = logFC,
  q_val_thresh = 0.05,
  formula_string = "analysis_type ~ comparison"
) {
  volplot_gg.all <- selected_data |>
    ggplot(aes(y = {{ log_q_value_column }}, x = {{ log_fc_column }})) +
    geom_point(aes(col = colour)) +
    scale_colour_manual(
      values = c(levels(selected_data$colour)),
      labels = c(
        paste0(
          "Not significant, logFC > ",
          1
        ),
        paste0(
          "Significant, logFC >= ",
          1
        ),
        paste0(
          "Significant, logFC <",
          1
        ),
        "Not Significant"
      )
    ) +
    geom_vline(xintercept = 1, colour = "black", linewidth = 0.2) +
    geom_vline(xintercept = -1, colour = "black", linewidth = 0.2) +
    geom_hline(yintercept = -log10(q_val_thresh)) +
    theme_bw() +
    xlab("Log fold changes") +
    ylab("-log10 q-value") +
    theme(legend.position = "none")

  volplot_gg.plot <- volplot_gg.all
  if (!is.na(formula_string) | formula_string != "") {
    volplot_gg.plot <- volplot_gg.all +
      facet_grid(as.formula(formula_string),
        labeller = labeller(facet_category = label_wrap_gen(width = 10))
      )
  }

  volplot_gg.plot
}

