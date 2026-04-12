# ----------------------------------------------------------------------------
# writeInteractiveVolcanoPlotProteomics
# ----------------------------------------------------------------------------
## Create proteomics interactive volcano plot
#' @export
# da_analysis_results_list$contrasts_results$fit.eb
# No full stops in the nme of columns of interactive table in glimma plot. It won't display column with full stop in the column name.
writeInteractiveVolcanoPlotProteomics <- function(
  da_proteins_long,
  uniprot_tbl,
  fit.eb = NULL, # No longer strictly required but kept for signature compatibility
  publication_graphs_dir,
  args_row_id = "uniprot_acc",
  fdr_column = "fdr_qvalue",
  raw_p_value_column = "raw_pvalue",
  log2fc_column = "log2FC",
  da_q_val_thresh = 0.05,
  counts_tbl = NULL,
  groups = NULL,
  uniprot_id_column = "Entry",
  gene_names_column = "gene_names",
  display_columns = c("best_uniprot_acc")
) {

  logger::log_info("--- Entering writeInteractiveVolcanoPlotProteomics ---")

  # Dynamically identify row ID column
  if (!(args_row_id %in% names(da_proteins_long))) {
    potential_ids <- c("Protein.Ids", "Protein.ID", "Entry", "uniprot_acc", "sites_id")
    found_id <- names(da_proteins_long)[names(da_proteins_long) %in% potential_ids][1]
    if (!is.na(found_id)) args_row_id <- found_id
  }

  volcano_plot_tab <- da_proteins_long |>
    dplyr::mutate(
      best_uniprot_acc = stringr::str_extract(!!sym(args_row_id), "^[^|:]+"),
      best_uniprot_acc_base = gsub("-\\d+$", "", best_uniprot_acc)
    )

  if (!is.null(uniprot_tbl)) {
    id_col_uniprot <- intersect(c("Entry", "UNIPROTKB", "Protein.Ids"), names(uniprot_tbl))[1]
    gene_col_uniprot <- intersect(c("gene_names", "GENENAME", "Gene.Name", "GeneName"), names(uniprot_tbl))[1]

    if (!is.na(id_col_uniprot) && !is.na(gene_col_uniprot)) {
      mapping_df <- uniprot_tbl |>
        dplyr::select(all_of(c(id_col_uniprot, gene_col_uniprot))) |>
        dplyr::rename(Entry = !!sym(id_col_uniprot), GeneSymbol = !!sym(gene_col_uniprot)) |>
        dplyr::mutate(GeneSymbol = stringr::str_extract(GeneSymbol, "^[^ ;:]+")) |>
        dplyr::distinct(Entry, .keep_all = TRUE)

      volcano_plot_tab <- volcano_plot_tab |>
        dplyr::left_join(mapping_df, by = c("best_uniprot_acc_base" = "Entry")) |>
        dplyr::mutate(
          gene_name = dplyr::coalesce(GeneSymbol, best_uniprot_acc)
        )
    } else {
      volcano_plot_tab$gene_name <- volcano_plot_tab$best_uniprot_acc
    }
  } else {
    volcano_plot_tab$gene_name <- volcano_plot_tab$best_uniprot_acc
  }


  output_dir <- file.path(
    publication_graphs_dir,
    "Interactive_Volcano_Plots"
  )

  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  # Iterate over uniquely defined contrasts instead of relying on coefficients index
  unique_contrasts <- unique(da_proteins_long$comparison)
  logger::log_info(sprintf("   Found %d unique contrasts to plot", length(unique_contrasts)))

  purrr::walk(
    unique_contrasts,
    \(contrast_name) {
      logger::log_info(sprintf("   Generating static interactive plot for: %s", contrast_name))

      contrast_data <- volcano_plot_tab |>
        dplyr::filter(comparison == contrast_name)

      # Filter counts to only the IDs in the contrast results to avoid length mismatches
      counts_mat_filtered <- NULL
      if (!is.null(counts_tbl)) {
        ids_in_contrast <- contrast_data[[args_row_id]]
        valid_ids <- intersect(ids_in_contrast, rownames(counts_tbl))
        if(length(valid_ids) > 0) {
           counts_mat_filtered <- counts_tbl[valid_ids, , drop = FALSE]
        }
      }

      getGlimmaVolcanoProteomics(
        volcano_plot_tab = contrast_data,
        uniprot_column = best_uniprot_acc,
        gene_name_column = gene_name,
        display_columns = display_columns,
        counts_tbl = counts_mat_filtered,
        groups = groups,
        output_dir = output_dir,
        fdr_column = fdr_column,
        log2fc_column = log2fc_column,
        da_q_val_thresh = da_q_val_thresh,
        contrast_name = contrast_name
      )
    }
  )
  logger::log_info("--- Exiting writeInteractiveVolcanoPlotProteomics ---")
}

# ----------------------------------------------------------------------------
# writeInteractiveVolcanoPlotProteomicsWidget
# ----------------------------------------------------------------------------
#' @export
# da_analysis_results_list$contrasts_results$fit.eb
# No full stops in the nme of columns of interactive table in glimma plot. It won't display column with full stop in the column name.
writeInteractiveVolcanoPlotProteomicsWidget <- function(
  da_proteins_long,
  uniprot_tbl,
  fit.eb,
  args_row_id = "uniprot_acc",
  fdr_column = "fdr_qvalue",
  raw_p_value_column = "raw_pvalue",
  log2fc_column = "log2FC",
  da_q_val_thresh = 0.05,
  counts_tbl = NULL,
  groups = NULL,
  uniprot_id_column = "Entry",
  gene_names_column = "gene_names",
  display_columns = c("best_uniprot_acc")
) {
  volcano_plot_tab <- da_proteins_long |>
    left_join(uniprot_tbl, by = join_by(!!sym(args_row_id) == !!sym(uniprot_id_column))) |>
    dplyr::rename(UNIPROT_GENENAME = gene_names_column) |>
    mutate(UNIPROT_GENENAME = purrr::map_chr(UNIPROT_GENENAME, \(x){
      tryCatch(
        {
          if (is.na(x) || is.null(x) || x == "") {
            ""
          } else {
            split_result <- str_split(x, " ")[[1]]
            if (length(split_result) > 0) split_result[1] else ""
          }
        },
        error = function(e) ""
      )
    })) |>
    mutate(lqm = -log10(!!sym(fdr_column))) |>
    dplyr::mutate(label = case_when(
      abs(!!sym(log2fc_column)) >= 1 & !!sym(fdr_column) >= da_q_val_thresh ~ "Not sig., logFC >= 1",
      abs(!!sym(log2fc_column)) >= 1 & !!sym(fdr_column) < da_q_val_thresh ~ "Sig., logFC >= 1",
      abs(!!sym(log2fc_column)) < 1 & !!sym(fdr_column) < da_q_val_thresh ~ "Sig., logFC < 1",
      TRUE ~ "Not sig."
    )) |>
    dplyr::mutate(colour = case_when(
      abs(!!sym(log2fc_column)) >= 1 & !!sym(fdr_column) >= da_q_val_thresh ~ "orange",
      abs(!!sym(log2fc_column)) >= 1 & !!sym(fdr_column) < da_q_val_thresh ~ "purple",
      abs(!!sym(log2fc_column)) < 1 & !!sym(fdr_column) < da_q_val_thresh ~ "blue",
      TRUE ~ "black"
    )) |>
    dplyr::mutate(gene_name = str_split(UNIPROT_GENENAME, " |:") |> purrr::map_chr(1)) |>
    dplyr::mutate(best_uniprot_acc = str_split(!!sym(args_row_id), ":") |> purrr::map_chr(1)) |>
    dplyr::mutate(analysis_type = comparison) |>
    dplyr::select(
      best_uniprot_acc, lqm, !!sym(fdr_column), !!sym(raw_p_value_column), !!sym(log2fc_column), comparison, label, colour, gene_name,
      any_of(display_columns)
    ) |>
    dplyr::mutate(my_alpha = case_when(
      gene_name != "" ~ 1,
      TRUE ~ 0.5
    ))

  interactive_volcano_plots <- purrr::map(
    seq_len(ncol(fit.eb$coefficients)),
    \(coef) {
      # print(paste0( "coef = ", coef))
      getGlimmaVolcanoProteomicsWidget(fit.eb,
        coef = coef,
        volcano_plot_tab = volcano_plot_tab,
        uniprot_column = best_uniprot_acc,
        gene_name_column = gene_name,
        display_columns = display_columns,
        counts_tbl = counts_tbl,
        groups = groups
      )
    }
  )

  interactive_volcano_plots
}
