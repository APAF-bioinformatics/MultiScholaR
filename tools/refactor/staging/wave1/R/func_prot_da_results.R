# ----------------------------------------------------------------------------
# outputDaAnalysisResults
# ----------------------------------------------------------------------------
#' @export
outputDaAnalysisResults <- function(
  da_analysis_results_list,
  theObject,
  uniprot_tbl,
  da_output_dir = NULL,
  publication_graphs_dir = NULL,
  file_prefix = NULL,
  plots_format = NULL,
  args_row_id = NULL,
  da_q_val_thresh = NULL,
  gene_names_column = NULL,
  fdr_column = NULL,
  raw_p_value_column = NULL,
  log2fc_column = NULL,
  uniprot_id_column = NULL,
  display_columns = NULL
) {
  cat("*** ENTERING outputDaAnalysisResults ***\n")
  cat("DEBUG: file_prefix =", file_prefix, "\n")
  cat("DEBUG: da_output_dir =", da_output_dir, "\n")
  cat("DEBUG: uniprot_tbl is null:", is.null(uniprot_tbl), "\n")
  cat("DEBUG: da_analysis_results_list names:", paste(names(da_analysis_results_list), collapse = ", "), "\n")


  uniprot_tbl <- checkParamsObjectFunctionSimplify(theObject, "uniprot_tbl", NULL)
  da_output_dir <- checkParamsObjectFunctionSimplify(theObject, "da_output_dir", NULL)
  publication_graphs_dir <- checkParamsObjectFunctionSimplify(theObject, "publication_graphs_dir", NULL)
  file_prefix <- checkParamsObjectFunctionSimplify(theObject, "file_prefix", "da_proteins")
  plots_format <- checkParamsObjectFunctionSimplify(theObject, "plots_format", c("pdf", "png"))
  args_row_id <- checkParamsObjectFunctionSimplify(theObject, "args_row_id", "uniprot_acc")
  da_q_val_thresh <- checkParamsObjectFunctionSimplify(theObject, "da_q_val_thresh", 0.05)
  gene_names_column <- checkParamsObjectFunctionSimplify(theObject, "gene_names_column", "gene_names")
  fdr_column <- checkParamsObjectFunctionSimplify(theObject, "fdr_column", "fdr_qvalue")
  raw_p_value_column <- checkParamsObjectFunctionSimplify(theObject, "raw_p_value_column", "raw_pvalue")
  log2fc_column <- checkParamsObjectFunctionSimplify(theObject, "log2fc_column", "log2FC")
  uniprot_id_column <- checkParamsObjectFunctionSimplify(theObject, "uniprot_id_column", "Entry")
  display_columns <- checkParamsObjectFunctionSimplify(theObject, "display_columns", c("best_uniprot_acc"))

  theObject <- updateParamInObject(theObject, "uniprot_tbl")
  theObject <- updateParamInObject(theObject, "da_output_dir")
  theObject <- updateParamInObject(theObject, "publication_graphs_dir")
  theObject <- updateParamInObject(theObject, "file_prefix")
  theObject <- updateParamInObject(theObject, "plots_format")
  theObject <- updateParamInObject(theObject, "args_row_id")
  theObject <- updateParamInObject(theObject, "da_q_val_thresh")
  theObject <- updateParamInObject(theObject, "gene_names_column")
  theObject <- updateParamInObject(theObject, "fdr_column")
  theObject <- updateParamInObject(theObject, "raw_p_value_column")
  theObject <- updateParamInObject(theObject, "log2fc_column")
  theObject <- updateParamInObject(theObject, "uniprot_id_column")
  theObject <- updateParamInObject(theObject, "display_columns")

  cat("DEBUG: Finished parameter setup, starting plot creation\n")

  ## PCA plot
  cat("DEBUG: Starting PCA plot section\n")
  plot_pca_plot <- da_analysis_results_list$pca_plot

  dir.create(file.path(publication_graphs_dir, "PCA"),
    recursive = TRUE,
    showWarnings = FALSE
  )
  cat("DEBUG: PCA plot section completed\n")

  for (format_ext in plots_format) {
    file_name <- file.path(publication_graphs_dir, "PCA", paste0("PCA_plot.", format_ext))
    ggsave(filename = file_name, plot = plot_pca_plot, limitsize = FALSE)
  }

  plot_pca_plot_with_labels <- da_analysis_results_list$pca_plot_with_labels
  for (format_ext in plots_format) {
    file_name <- file.path(publication_graphs_dir, "PCA", paste0("PCA_plot_with_sample_ids.", format_ext))
    ggsave(filename = file_name, plot = plot_pca_plot_with_labels, limitsize = FALSE)
  }

  ## RLE plot
  plot_rle_plot <- da_analysis_results_list$rle_plot

  dir.create(file.path(publication_graphs_dir, "RLE"),
    recursive = TRUE,
    showWarnings = FALSE
  )

  for (format_ext in plots_format) {
    file_name <- file.path(publication_graphs_dir, "RLE", paste0("RLE_plot.", format_ext))
    ggsave(filename = file_name, plot = plot_rle_plot, limitsize = FALSE)
  }

  ## Save the number of values graph
  plot_num_of_values <- da_analysis_results_list$plot_num_of_values

  for (format_ext in plots_format) {
    file_name <- file.path(da_output_dir, paste0("num_of_values.", format_ext))
    ggsave(filename = file_name, plot = plot_num_of_values, limitsize = FALSE)
  }

  ## Contrasts results
  ## This plot is used to check the mean-variance relationship of the expression data, after fitting a linear model.
  pdf(file.path(da_output_dir, "plotSA_after_ruvIII.pdf"))
  plotSA(da_analysis_results_list$contrasts_results$fit.eb)
  dev.off()

  png(file.path(da_output_dir, "plotSA_after_ruvIII.png"))
  plotSA(da_analysis_results_list$contrasts_results$fit.eb)
  dev.off()

  saveRDS(
    da_analysis_results_list$contrasts_results$fit.eb,
    file.path(da_output_dir, "fit.eb.RDS")
  )

  ## Values for volcano plts

  ## Write all the results in one single table
  significant_rows <- da_analysis_results_list$significant_rows

  significant_rows |>
    dplyr::select(-colour, -lqm) |>
    vroom::vroom_write(file.path(da_output_dir, "lfc_qval_long.tsv"))

  significant_rows |>
    dplyr::select(-colour, -lqm) |>
    writexl::write_xlsx(file.path(da_output_dir, "lfc_qval_long.xlsx"))

  ## Print Volcano plot
  volplot_plot <- da_analysis_results_list$volplot_plot

  for (format_ext in plots_format) {
    file_name <- file.path(da_output_dir, paste0("volplot_gg_all.", format_ext))
    ggsave(filename = file_name, plot = volplot_plot, width = 7.29, height = 6)
  }


  ## Number of values graph
  plot_num_of_values <- da_analysis_results_list$plot_num_of_values

  for (format_ext in plots_format) {
    file_name <- file.path(da_output_dir, paste0("num_of_values.", format_ext))
    ggsave(filename = file_name, plot = plot_num_of_values, limitsize = FALSE)
  }

  ## Contrasts results
  ## This plot is used to check the mean-variance relationship of the expression data, after fitting a linear model.
  contrasts_results <- da_analysis_results_list$contrasts_results
  for (format_ext in plots_format) {
    file_name <- file.path(da_output_dir, paste0("plotSA_after_ruvIII", format_ext))

    if (format_ext == "pdf") {
      pdf(file_name)
    } else if (format_ext == "png") {
      png(file_name)
    }

    plotSA(contrasts_results$fit.eb)
    dev.off()
  }

  saveRDS(
    contrasts_results$fit.eb,
    file.path(da_output_dir, "fit.eb.RDS")
  )

  ## Values for volcano plts

  ## Write all the results in one single table
  significant_rows <- da_analysis_results_list$significant_rows
  significant_rows |>
    dplyr::select(-colour, -lqm) |>
    vroom::vroom_write(file.path(da_output_dir, "lfc_qval_long.tsv"))

  significant_rows |>
    dplyr::select(-colour, -lqm) |>
    writexl::write_xlsx(file.path(da_output_dir, "lfc_qval_long.xlsx"))


  ## Count the number of up or down significnat differentially expressed proteins.
  if (!is.null(da_analysis_results_list$num_sig_da_genes_barplot_only_significant)) {
    num_sig_da_genes_barplot_only_significant <- da_analysis_results_list$num_sig_da_genes_barplot_only_significant
    num_of_comparison_only_significant <- da_analysis_results_list$num_of_comparison_only_significant

    savePlot(num_sig_da_genes_barplot_only_significant,
      base_path = da_output_dir,
      plot_name = paste0(file_prefix, "_num_sda_entities_barplot_only_significant"),
      formats = c("pdf", "png", "svg"),
      width = (num_of_comparison_only_significant + 2) * 7 / 6,
      height = 6
    )
  }


  ## Count the number of up or down significnat differentially expressed proteins.
  num_sig_da_molecules_first_go <- da_analysis_results_list$num_sig_da_molecules_first_go
  vroom::vroom_write(
    num_sig_da_molecules_first_go$table,
    file.path(
      da_output_dir,
      paste0(file_prefix, "_num_significant_differentially_abundant_all.tab")
    )
  )

  writexl::write_xlsx(
    num_sig_da_molecules_first_go$table,
    file.path(
      da_output_dir,
      paste0(file_prefix, "_num_significant_differentially_abundant_all.xlsx")
    )
  )


  ## Print p-values distribution figure
  pvalhist <- da_analysis_results_list$pvalhist
  for (format_ext in plots_format) {
    file_name <- file.path(da_output_dir, paste0(file_prefix, "_p_values_distn.", format_ext))
    ggsave(
      filename = file_name,
      plot = pvalhist,
      height = 10,
      width = 7
    )
  }


  ## Create wide format output file
  cat("DEBUG: Starting wide format output section\n")
  da_proteins_wide <- da_analysis_results_list$da_proteins_wide
  cat("DEBUG: da_proteins_wide extracted successfully\n")

  vroom::vroom_write(
    da_proteins_wide,
    file.path(
      da_output_dir,
      paste0(file_prefix, "_wide.tsv")
    )
  )

  writexl::write_xlsx(
    da_proteins_wide,
    file.path(
      da_output_dir,
      paste0(file_prefix, "_wide.xlsx")
    )
  )
  cat("DEBUG: Wide format files written\n")

  cat("DEBUG: Starting wide_annot creation\n")

  tryCatch(
    {
      da_proteins_wide_annot <- da_proteins_wide |>
        mutate(uniprot_acc_cleaned = str_split(!!sym(args_row_id), "-") |>
          purrr::map_chr(1)) |>
        left_join(uniprot_tbl, by = join_by(uniprot_acc_cleaned == Entry)) |>
        dplyr::select(-uniprot_acc_cleaned) |>
        mutate(gene_name = ifelse(
          !is.na(gene_names) & gene_names != "",
          sapply(gene_names, function(x) {
            if (is.na(x) || x == "") {
              return("")
            }
            split_result <- strsplit(as.character(x), " |:")[[1]]
            if (length(split_result) > 0) split_result[1] else ""
          }),
          ""
        )) |>
        relocate(gene_name, .after = !!sym(args_row_id))

      cat("DEBUG: wide_annot creation successful\n")
    },
    error = function(e) {
      cat("DEBUG: ERROR in wide_annot creation:", e$message, "\n")
      # Create a fallback version without annotations
      da_proteins_wide_annot <- da_proteins_wide |>
        mutate(gene_name = "")
      cat("DEBUG: Created fallback wide_annot without annotations\n")
    }
  )

  vroom::vroom_write(
    da_proteins_wide_annot,
    file.path(
      da_output_dir,
      paste0(file_prefix, "_wide_annot.tsv")
    )
  )

  writexl::write_xlsx(
    da_proteins_wide_annot,
    file.path(
      da_output_dir,
      paste0(file_prefix, "_wide_annot.xlsx")
    )
  )

  cat("DEBUG: Wide_annot files written, proceeding to long format\n")

  ## Create long format output file
  da_proteins_long <- da_analysis_results_list$da_proteins_long

  cat("DEBUG: da_proteins_long exists:", !is.null(da_proteins_long), "\n")
  if (!is.null(da_proteins_long)) {
    cat("DEBUG: da_proteins_long dimensions:", dim(da_proteins_long), "\n")
  } else {
    cat("DEBUG: da_proteins_long is NULL - cannot create long_annot\n")
    return(NULL)
  }

  vroom::vroom_write(
    da_proteins_long,
    file.path(
      da_output_dir,
      paste0(file_prefix, "_long.tsv")
    )
  )

  writexl::write_xlsx(
    da_proteins_long,
    file.path(
      da_output_dir,
      paste0(file_prefix, "_long.xlsx")
    )
  )

  cat("DEBUG: Starting long_annot creation\n")
  cat("DEBUG: da_proteins_long dimensions:", dim(da_proteins_long), "\n")
  cat("DEBUG: uniprot_tbl is null:", is.null(uniprot_tbl), "\n")
  if (!is.null(uniprot_tbl)) {
    cat("DEBUG: uniprot_tbl dimensions:", dim(uniprot_tbl), "\n")
  }
  cat("DEBUG: args_row_id:", args_row_id, "\n")

  da_proteins_long_annot <- da_proteins_long |>
    mutate(uniprot_acc_cleaned = str_split(!!sym(args_row_id), "-") |>
      purrr::map_chr(1)) |>
    left_join(uniprot_tbl, by = join_by(uniprot_acc_cleaned == Entry)) |>
    dplyr::select(-uniprot_acc_cleaned) |>
    mutate(gene_name = ifelse(
      !is.na(gene_names) & gene_names != "",
      sapply(gene_names, function(x) {
        if (is.na(x) || x == "") {
          return("")
        }
        split_result <- strsplit(as.character(x), " |:")[[1]]
        if (length(split_result) > 0) split_result[1] else ""
      }),
      ""
    )) |>
    relocate(gene_name, .after = !!sym(args_row_id))

  cat("DEBUG: da_proteins_long_annot dimensions:", dim(da_proteins_long_annot), "\n")
  cat("DEBUG: da_proteins_long_annot is null:", is.null(da_proteins_long_annot), "\n")
  if (!is.null(da_proteins_long_annot) && nrow(da_proteins_long_annot) > 0) {
    cat("DEBUG: long_annot columns:", paste(names(da_proteins_long_annot), collapse = ", "), "\n")
  } else {
    cat("DEBUG: long_annot is empty or null - annotation pipeline failed\n")
  }

  long_annot_file_path <- file.path(da_output_dir, paste0(file_prefix, "_long_annot.tsv"))
  cat("DEBUG: Attempting to write long_annot to:", long_annot_file_path, "\n")

  vroom::vroom_write(da_proteins_long_annot, long_annot_file_path)

  cat("DEBUG: long_annot file written, checking if exists:", file.exists(long_annot_file_path), "\n")

  writexl::write_xlsx(
    da_proteins_long_annot,
    file.path(
      da_output_dir,
      paste0(file_prefix, "_long_annot.xlsx")
    )
  )

  ## Static volcano plots
  dir.create(file.path(publication_graphs_dir, "Volcano_Plots"),
    recursive = TRUE,
    showWarnings = FALSE
  )

  list_of_volcano_plots <- da_analysis_results_list$list_of_volcano_plots

  # Print diagnostic info about the volcano plots
  message(sprintf("Number of volcano plots: %d", nrow(list_of_volcano_plots)))

  purrr::walk2(
    list_of_volcano_plots %>% dplyr::pull(title),
    list_of_volcano_plots %>% dplyr::pull(plot),
    \(x, y){
      # gg_save_logging ( .y, file_name_part, plots_format)

      savePlot(y,
        base_path = file.path(publication_graphs_dir, "Volcano_Plots"),
        plot_name = x,
        formats = plots_format, width = 7, height = 7
      )
    }
  )

  # Generate a multi-page PDF with all volcano plots
  volcano_plots_list <- list_of_volcano_plots %>% dplyr::pull(plot)

  # Generate combined PDF with all plots, one per page
  pdf_file <- file.path(publication_graphs_dir, "Volcano_Plots", "list_of_volcano_plots.pdf")
  pdf(file = pdf_file, width = 7, height = 7, onefile = TRUE)
  purrr::walk(volcano_plots_list, print)
  invisible(dev.off())

  # Verify the PDF was created with the right number of pages
  message(sprintf("Created multi-page PDF at %s", pdf_file))

  list_of_volcano_plots_with_gene_names <- da_analysis_results_list$list_of_volcano_plots_with_gene_names

  # Print diagnostic info about the labeled volcano plots
  message(sprintf("Number of labeled volcano plots: %d", nrow(list_of_volcano_plots_with_gene_names)))

  purrr::walk2(
    list_of_volcano_plots_with_gene_names %>% dplyr::pull(title),
    list_of_volcano_plots_with_gene_names %>% dplyr::pull(plot),
    \(x, y) {
      savePlot(
        x,
        file.path(publication_graphs_dir, "Volcano_Plots"),
        paste0(y, "_with_protein_labels")
      )
    }
  )

  # Generate a multi-page PDF with all labeled volcano plots
  volcano_plots_with_genes_list <- list_of_volcano_plots_with_gene_names %>% dplyr::pull(plot)

  # Generate combined PDF with all labeled plots, one per page
  pdf_file_with_genes <- file.path(publication_graphs_dir, "Volcano_Plots", "list_of_volcano_plots_with_gene_names.pdf")
  pdf(file = pdf_file_with_genes, width = 7, height = 7, onefile = TRUE)
  purrr::walk(volcano_plots_with_genes_list, print)
  invisible(dev.off())

  # Verify the labeled PDF was created with the right number of pages
  message(sprintf("Created multi-page labeled PDF at %s", pdf_file_with_genes))

  ## Number of significant molecules
  createDirIfNotExists(file.path(publication_graphs_dir, "NumSigDaMolecules"))
  vroom::vroom_write(
    da_analysis_results_list$num_sig_da_molecules,
    file.path(publication_graphs_dir, "NumSigDaMolecules", paste0(file_prefix, "_num_sig_da_molecules.tab"))
  )


  if (!is.null(da_analysis_results_list$num_sig_da_genes_barplot_only_significant)) {
    num_sig_da_genes_barplot_only_significant <- da_analysis_results_list$num_sig_da_genes_barplot_only_significant
    num_of_comparison_only_significant <- da_analysis_results_list$num_of_comparison_only_significant

    savePlot(num_sig_da_genes_barplot_only_significant,
      base_path = file.path(publication_graphs_dir, "NumSigDaMolecules"),
      plot_name = paste0(file_prefix, "_num_sig_da_molecules."),
      formats = plots_format,
      width = (num_of_comparison_only_significant + 2) * 7 / 6,
      height = 6
    )
  }


  if (!is.null(da_analysis_results_list$num_sig_da_genes_barplot_with_not_significant)) {
    num_sig_da_genes_barplot_with_not_significant <- da_analysis_results_list$num_sig_da_genes_barplot_with_not_significant
    num_of_comparison_with_not_significant <- da_analysis_results_list$num_of_comparison_with_not_significant

    print("print bar plot")

    savePlot(num_sig_da_genes_barplot_with_not_significant,
      base_path = file.path(publication_graphs_dir, "NumSigDaMolecules"),
      plot_name = paste0(file_prefix, "_num_sig_da_molecules_with_not_significant"),
      formats = plots_format,
      width = (num_of_comparison_with_not_significant + 2) * 7 / 6,
      height = 6
    )
  }

  ## Write interactive volcano plot
  counts_mat <- (da_analysis_results_list$theObject)@protein_quant_table |>
    column_to_rownames((da_analysis_results_list$theObject)@protein_id_column) |>
    as.matrix()

  this_design_matrix <- da_analysis_results_list$theObject@design_matrix

  rownames(this_design_matrix) <- this_design_matrix[, da_analysis_results_list$theObject@sample_id]

  this_groups <- this_design_matrix[colnames(counts_mat), "group"]

  writeInteractiveVolcanoPlotProteomics(da_proteins_long,
    uniprot_tbl = uniprot_tbl,
    fit.eb = contrasts_results$fit.eb,
    publication_graphs_dir = publication_graphs_dir,
    args_row_id = args_row_id,
    fdr_column = fdr_column,
    raw_p_value_column = raw_p_value_column,
    log2fc_column = log2fc_column,
    da_q_val_thresh = da_q_val_thresh,
    counts_tbl = counts_mat,
    groups = this_groups,
    uniprot_id_column = uniprot_id_column,
    gene_names_column = gene_names_column,
    display_columns = display_columns
  )
}

# ----------------------------------------------------------------------------
# saveDaProteinList
# ----------------------------------------------------------------------------
#' Save the list of output tables from differential expression analysis of proteins or phosphopeptides into a file and in a specific directory.
#' @param list_of_da_tables A list, each element is a table of log fold-change and q-values from differential expression analysis of proteins / phosphopeptides. Each element in the list has a name, usually the name of the pairwise comparison.
#' @param row_id Add row ID to the output table based on the name (protein or phosphopeptid ID) of each row
#' @param sort_by_column Each table in the list_of_da_tables is sorted in ascending order
#' @param results_dir The results directory to store the output file
#' @param file_suffix The file suffix string to aadd to the name of each comparison from the list_of_da_tables.
#' @export
saveDaProteinList <- function(list_of_da_tables, row_id, sort_by_column = fdr_qvalue, results_dir, file_suffix) {
  purrr::walk2(
    list_of_da_tables, names(list_of_da_tables),
    \(.x, .y) {
      vroom::vroom_write(
        .x |>
          rownames_to_column(row_id) |>
          arrange({{ sort_by_column }}),
        path = file.path(results_dir, paste0(.y, file_suffix))
      )
    }
  )
}

# ----------------------------------------------------------------------------
# createDaResultsLongFormat
# ----------------------------------------------------------------------------
## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Create the de_protein_long and de_phos_long tables
#' @export
createDaResultsLongFormat <- function(lfc_qval_tbl,
                                      norm_counts_input_tbl,
                                      raw_counts_input_tbl,
                                      row_id,
                                      sample_id,
                                      group_id,
                                      group_pattern,
                                      design_matrix_norm,
                                      design_matrix_raw,
                                      expression_column = log_intensity,
                                      protein_id_table) {
  message("   DEBUG66: createDaResultsLongFormat - Starting norm_counts processing")
  message(sprintf("      DEBUG66: norm_counts_input_tbl dims = %d x %d", nrow(norm_counts_input_tbl), ncol(norm_counts_input_tbl)))
  message(sprintf("      DEBUG66: group_pattern = %s", group_pattern))
  message(sprintf("      DEBUG66: row_id = %s", row_id))
  message(sprintf("      DEBUG66: sample_id = %s", sample_id))
  message(sprintf("      DEBUG66: group_id = %s", group_id))

  norm_counts <- norm_counts_input_tbl |>
    as.data.frame() |>
    rownames_to_column(row_id) |>
    pivot_longer(
      cols = matches(group_pattern),
      names_to = sample_id,
      values_to = "log2norm"
    ) |>
    left_join(design_matrix_norm, by = sample_id) |>
    group_by(!!sym(row_id), !!sym(group_id)) |>
    arrange(!!sym(row_id), !!sym(group_id), !!sym(sample_id)) |>
    mutate(replicate_number = paste0("log2norm.", row_number())) |>
    ungroup() |>
    pivot_wider(
      id_cols = c(!!sym(row_id), !!sym(group_id)),
      names_from = replicate_number,
      values_from = log2norm
    ) |>
    mutate({{ group_id }} := purrr::map_chr(!!sym(group_id), as.character))

  message("   DEBUG66: norm_counts processing completed")


  # print(head(norm_counts))

  raw_counts <- raw_counts_input_tbl |>
    as.data.frame() |>
    rownames_to_column(row_id) |>
    pivot_longer(
      cols = matches(group_pattern),
      names_to = sample_id,
      values_to = "raw"
    ) |>
    left_join(design_matrix_raw, by = sample_id) |>
    group_by(!!sym(row_id), !!sym(group_id)) |>
    arrange(!!sym(row_id), !!sym(group_id), !!sym(sample_id)) |>
    mutate(replicate_number = paste0("raw.", row_number())) |>
    ungroup() |>
    pivot_wider(
      id_cols = c(!!sym(row_id), !!sym(group_id)),
      names_from = replicate_number,
      values_from = raw
    ) |>
    mutate({{ group_id }} := purrr::map_chr(!!sym(group_id), as.character))

  # print(head(raw_counts))

  left_join_columns <- rlang::set_names(
    c(row_id, group_id),
    c(row_id, "left_group")
  )

  right_join_columns <- rlang::set_names(
    c(row_id, group_id),
    c(row_id, "right_group")
  )

  # print(head(lfc_qval_tbl))

  # DEBUG66: Commented out print statements that were causing confusion
  # print( row_id)
  # print(colnames( protein_id_table)[1])

  da_proteins_long <- lfc_qval_tbl |>
    dplyr::select(-lqm, -colour, -analysis_type) |>
    dplyr::mutate({{ expression_column }} := str_replace_all({{ expression_column }}, group_id, "")) |>
    separate_wider_delim({{ expression_column }}, delim = "-", names = c("left_group", "right_group")) |>
    left_join(norm_counts, by = left_join_columns) |>
    left_join(norm_counts,
      by = right_join_columns,
      suffix = c(".left", ".right")
    ) |>
    left_join(raw_counts, by = left_join_columns) |>
    left_join(raw_counts,
      by = right_join_columns,
      suffix = c(".left", ".right")
    ) |>
    left_join(protein_id_table,
      by = join_by(!!sym(row_id) == !!sym(colnames(protein_id_table)[1]))
    ) |>
    arrange(comparison, fdr_qvalue, log2FC) |>
    distinct()

  # --- NEW: Rename columns to use sample IDs if single contrast ---
  # Only perform this renaming if we have a single comparison, to ensure unique mapping
  if (length(unique(da_proteins_long$comparison)) == 1) {
    # Get the groups involved
    this_left_group <- unique(da_proteins_long$left_group)
    this_right_group <- unique(da_proteins_long$right_group)

    # Ensure we have exactly one left and one right group
    if (length(this_left_group) == 1 && length(this_right_group) == 1) {
      message(sprintf("   createDaResultsLongFormat: Renaming columns for contrast %s vs %s", this_left_group, this_right_group))

      # Helper to get sorted sample IDs for a group
      get_samples_for_group <- function(dm, grp) {
        dm |>
          dplyr::filter(!!sym(group_id) == grp) |>
          dplyr::arrange(!!sym(sample_id)) |>
          dplyr::pull(!!sym(sample_id))
      }

      left_samples <- get_samples_for_group(design_matrix_norm, this_left_group)
      right_samples <- get_samples_for_group(design_matrix_norm, this_right_group)

      # Helper to generate rename mapping using vectorized operations
      # Returns: named vector c(new_name = old_name) for dplyr::rename
      generate_rename_map <- function(df, prefix, suffix, samples, group_name) {
        indices <- seq_along(samples)
        old_cols <- paste0(prefix, ".", indices, suffix)
        new_cols <- paste0(prefix, ".", samples, ".", group_name)

        # Only include columns that exist in the dataframe
        valid_idx <- old_cols %in% colnames(df)

        if (any(valid_idx)) {
          return(setNames(old_cols[valid_idx], new_cols[valid_idx]))
        } else {
          return(character(0))
        }
      }

      # Generate mappings for all 4 sets of columns
      map1 <- generate_rename_map(da_proteins_long, "log2norm", ".left", left_samples, this_left_group)
      map2 <- generate_rename_map(da_proteins_long, "raw", ".left", left_samples, this_left_group)
      map3 <- generate_rename_map(da_proteins_long, "log2norm", ".right", right_samples, this_right_group)
      map4 <- generate_rename_map(da_proteins_long, "raw", ".right", right_samples, this_right_group)

      # Combine all mappings
      all_mappings <- c(map1, map2, map3, map4)

      # Apply renaming in a single vectorized step
      if (length(all_mappings) > 0) {
        da_proteins_long <- da_proteins_long |> dplyr::rename(!!!all_mappings)
        message(sprintf("   createDaResultsLongFormat: Renamed %d columns", length(all_mappings)))
      }
    }
  }
  # --- END NEW ---

  # Rename group columns to numerator/denominator for clarity
  da_proteins_long <- da_proteins_long |>
    dplyr::rename(numerator = left_group, denominator = right_group)

  da_proteins_long
}

# ----------------------------------------------------------------------------
# countStatDaGenes
# ----------------------------------------------------------------------------
#' Count the number of statistically significant differentially abundant proteins (according to user-defined threshold)
#' @param lfc_thresh A numerical value specifying the log fold-change threhold (absolute value) for calling statistically significant proteins.
#' @param q_val_thresh A numerical value specifying the q-value threshold for statistically significant proteins.
#' @param log_fc_column The name of the log fold-change column (tidyverse style).
#' @param q_value_column The name of the q-value column (tidyverse style).
#' @return A table with the following columns:
#' status: The status could be Significant and Up, Significant and Down or Not significant
#' counts: The number of proteins wit this status
#' @export
countStatDaGenes <- function(data,
                             lfc_thresh = 0,
                             q_val_thresh = 0.05,
                             log_fc_column = log2FC,
                             q_value_column = fdr_qvalue) {
  # comparison <- as.data.frame(data) |>
  #   distinct(comparison) |>
  #   pull(comparison)

  selected_data <- data |>
    dplyr::mutate(status = case_when(
      {{ q_value_column }} >= q_val_thresh ~ "Not significant",
      {{ log_fc_column }} >= lfc_thresh & {{ q_value_column }} < q_val_thresh ~ "Significant and Up",
      {{ log_fc_column }} < lfc_thresh & {{ q_value_column }} < q_val_thresh ~ "Significant and Down",
      TRUE ~ "Not significant"
    ))

  counts <- selected_data |>
    group_by(status) |>
    summarise(counts = n()) |>
    ungroup()

  all_possible_status <- data.frame(status = c("Not significant", "Significant and Up", "Significant and Down"))

  results <- all_possible_status |>
    left_join(counts, by = c("status" = "status")) |>
    mutate(counts = ifelse(is.na(counts), 0, counts))

  return(results)
}

# ----------------------------------------------------------------------------
# countStatDaGenesHelper
# ----------------------------------------------------------------------------
## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#' @export
countStatDaGenesHelper <- function(
  da_table,
  description,
  facet_column = analysis_type,
  comparison_column = "comparison",
  expression_column = "expression",
  lfc_thresh = 0,
  q_val_thresh = 0.05,
  log_fc_column = logFC,
  q_value_column = fdr_qvalue
) {
  message("--- Entering countStatDaGenesHelper (DEBUG66) ---")
  message(paste("   countStatDaGenesHelper: da_table class =", class(da_table)))
  message(paste("   countStatDaGenesHelper: da_table is list =", is.list(da_table)))
  message(paste("   countStatDaGenesHelper: da_table length =", length(da_table)))
  message("   countStatDaGenesHelper: da_table names:")
  print(names(da_table))

  # CRITICAL FIX 1: If da_table is a list of tables from limma, we need to process each one
  da_table_updated <- if (is.list(da_table) && !is.data.frame(da_table)) {
    purrr::imap(da_table, \(x, n) {
      if (is.data.frame(x)) {
        # Ensure comparison column is set if missing
        if (!"comparison" %in% colnames(x)) {
          x <- x |> dplyr::mutate(comparison = n)
        }
      }
      countStatDaGenes(x,
                       lfc_thresh = lfc_thresh,
                       q_val_thresh = q_val_thresh,
                       log_fc_column = {{ log_fc_column }},
                       q_value_column = {{ q_value_column }})
    })
  } else {
    list(countStatDaGenes(da_table,
                          lfc_thresh = lfc_thresh,
                       q_val_thresh = q_val_thresh,
                       log_fc_column = {{ log_fc_column }},
                       q_value_column = {{ q_value_column }}))
  }

  message("   countStatDaGenesHelper: da_table_updated created")
  message(paste("   countStatDaGenesHelper: da_table_updated length =", length(da_table_updated)))
  message("   countStatDaGenesHelper: da_table_updated names:")
  print(names(da_table_updated))

  list_of_tables <- purrr::map2(
    da_table_updated,
    names(da_table_updated),
    \(.x, .y){
      message(paste("      [map2] Processing element with name:", .y))
      .x |> mutate(!!sym(comparison_column) := .y)
    }
  )

  message("   countStatDaGenesHelper: list_of_tables created")
  message("   countStatDaGenesHelper: About to bind_rows...")

  bound_tables <- list_of_tables |> bind_rows()
  message(paste("   countStatDaGenesHelper: bound_tables dims =", nrow(bound_tables), "x", ncol(bound_tables)))

  faceted_tables <- bound_tables |> mutate({{ facet_column }} := description)
  message("   countStatDaGenesHelper: facet column added")
  message("   countStatDaGenesHelper: Checking comparison column values:")
  print(unique(faceted_tables[[comparison_column]]))

  message(paste("   countStatDaGenesHelper: About to separate_wider_delim on column:", comparison_column))
  message(paste("   countStatDaGenesHelper: Looking for delimiter: ="))

  # Check if any values contain "="
  has_delimiter <- any(grepl("=", faceted_tables[[comparison_column]]))
  message(paste("   countStatDaGenesHelper: Any values contain '=' ?", has_delimiter))

  if (has_delimiter) {
    merged_tables <- faceted_tables |>
      separate_wider_delim(!!sym(comparison_column),
        delim = "=",
        names = c(
          comparison_column,
          expression_column
        ),
        too_few = "align_start"
      )
    message("   countStatDaGenesHelper: Separation complete")
  } else {
    message("   countStatDaGenesHelper: No '=' found, skipping separation and using original comparison names")
    # Ensure expression column exists even if no separation
    merged_tables <- faceted_tables
    if (!expression_column %in% colnames(merged_tables)) {
      merged_tables[[expression_column]] <- merged_tables[[comparison_column]]
    }
  }

  message("--- Exiting countStatDaGenesHelper (DEBUG66) ---")
  merged_tables
}

# ----------------------------------------------------------------------------
# printCountDaGenesTable
# ----------------------------------------------------------------------------
## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#' Format results table for use in volcano plots, counting number of significant proteins, p-values distribution histogram.
#' @param list_of_da_tables A list with each element being a results table with log fold-change and q-value per protein.
#' @param list_of_descriptions  A list of strings describing the parameters used to generate the result table.
#' @param formula_string The formula string used in the facet_grid command for the ggplot scatter plot.
#' @param facet_column The name of the column describing the type of analysis or parameters used to generate the result table (tidyverse style). This is related to the \code{list_of_descriptions} parameter above.
#' @param comparison_column The name of the column describing the contrasts or comparison between groups (tidyverse style).
#' @param expression_column The name of the column that will contain the formula expressions of the contrasts.
#' @export
printCountDaGenesTable <- function(
  list_of_da_tables,
  list_of_descriptions,
  formula_string = "analysis_type ~ comparison",
  facet_column = analysis_type,
  comparison_column = "comparison",
  expression_column = "expression"
) {
  num_significant_da_genes_all <- purrr::map2(
    list_of_da_tables,
    list_of_descriptions,
    function(a, b) {
      countStatDaGenesHelper(
        da_table = a,
        description = b,
        facet_column = {{ facet_column }},
        comparison_column = comparison_column,
        expression_column = expression_column
      )
    }
  ) |>
    bind_rows()

  num_sig_da_genes_barplot <- num_significant_da_genes_all |>
    dplyr::filter(status != "Not significant") |>
    ggplot(aes(x = status, y = counts)) +
    geom_bar(stat = "identity") +
    geom_text(stat = "identity", aes(label = counts), vjust = -0.5) +
    theme(axis.text.x = element_text(angle = 90))

  # print(head(num_sig_da_genes_barplot))

  if (!is.na(formula_string)) {
    num_sig_da_genes_barplot <- num_sig_da_genes_barplot +
      facet_grid(as.formula(formula_string))
  }


  return(list(plot = num_sig_da_genes_barplot, table = num_significant_da_genes_all))
}

# ----------------------------------------------------------------------------
# getSignificantData
# ----------------------------------------------------------------------------
## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#' Format results table for use in volcano plots, counting number of significant proteins, p-values distribution histogram.
#' @param list_of_da_tables A list with each element being a results table with log fold-change and q-value per protein.
#' @param list_of_descriptions  A list of strings describing the parameters used to generate the result table.
#' @param row_id The name of the row ID column (tidyverse style).
#' @param p_value_column The name of the raw p-value column (tidyverse style).
#' @param q_value_column The name of the q-value column (tidyverse style).
#' @param log_q_value_column The name of the log q-value column (tidyverse style).
#' @param log_fc_column The name of the log fold-change column (tidyverse style).
#' @param comparison_column The name of the column describing the contrasts or comparison between groups (tidyverse style).
#' @param expression_column The name of the column that will contain the formula expressions of the contrasts.
#' @param facet_column The name of the column describing the type of analysis or parameters used to generate the result table (tidyverse style). This is related to the \code{list_of_descriptions} parameter above.
#' @param q_val_thresh A numerical value specifying the q-value threshold for statistically significant proteins.
#' @return A table with the following columns:
#' row_id:  The protein ID, this column is derived from the input to the row_id column.
#' log_q_value_column: The log (base 10) q-value, this column name is derived from the input to the log_q_value_column.
#' q_value_column:  The q-value, this column name is derived from the input to the q_value_column.
#'   p_value_column: The p-value, this column name is derived from the input to p_value_column.
#'   log_fc_column: The log fold-change, this column name is derived from the input to log_fc_column.
#'   comparison_column: The comparison, this column name is derived from the input to comparison_column.
#'   expression_column: The formula expression for the contrasts, this column name is derived from the input to expression_column.
#'   facet_column: The analysis type, this column name is derived from the input to facet_column.
#'   colour The colour of the dots used in the volcano plot.
#'   orange = Absolute Log fold-change >= 1 and q-value >= threshold
#'   purple = Absolute Log fold-change >= 1 and q-value < threshold
#'   blue = Absolute Log fold-change < 1 and q-value < threshold
#'   black = all other values
#' @export
getSignificantData <- function(
  list_of_da_tables,
  list_of_descriptions,
  row_id = uniprot_acc,
  p_value_column = raw_pvalue,
  q_value_column = fdr_qvalue,
  fdr_value_column = fdr_value_bh_adjustment,
  log_q_value_column = lqm,
  log_fc_column = logFC,
  comparison_column = "comparison",
  expression_column = "log_intensity",
  facet_column = analysis_type,
  q_val_thresh = 0.05
) {
  message("--- Entering getSignificantData (DEBUG66) ---")
  message(paste("   getSignificantData: list_of_da_tables class =", class(list_of_da_tables)))
  message(paste("   getSignificantData: list_of_da_tables length =", length(list_of_da_tables)))
  message("   getSignificantData: list_of_da_tables structure:")
  str(list_of_da_tables, max.level = 2)

  get_row_binded_table <- function(de_table_list, description) {
    message("   --- Entering get_row_binded_table (DEBUG66) ---")
    message(paste("      get_row_binded_table: de_table_list class =", class(de_table_list)))
    message(paste("      get_row_binded_table: de_table_list length =", length(de_table_list)))
    message("      get_row_binded_table: de_table_list names:")
    print(names(de_table_list))

    # This internal helper is now more robust. It checks if the table
    # already has the row_id as a column. If not, it converts rownames.
    # This handles both old and new data structures.

    processed_list <- purrr::map(de_table_list, function(tbl) {
      row_id_col <- as_string(as_name(enquo(row_id)))
      message(paste("         [map] Processing table, row_id_col =", row_id_col))
      message(paste("         [map] Table class =", class(tbl)))
      message(paste("         [map] Table is data.frame =", is.data.frame(tbl)))

      if (!row_id_col %in% colnames(tbl)) {
        # If row_id is not a column, convert rownames
        message(paste("         [map] row_id not in columns, converting rownames"))
        tbl <- tbl |> rownames_to_column(var = row_id_col)
      } else {
        message(paste("         [map] row_id already in columns"))
      }

      return(tbl)
    })

    message("      get_row_binded_table: processed_list created")
    message(paste("      get_row_binded_table: processed_list length =", length(processed_list)))

    output <- processed_list |>
      purrr::map2(names(processed_list), \(.x, .y){
        message(paste("         [map2] Processing element with name:", .y))
        message(paste("         [map2] comparison_column already exists:", comparison_column %in% colnames(.x)))
        # If the 'comparison' column doesn't already exist, create it from the list name
        if (!comparison_column %in% colnames(.x)) {
          message(paste("         [map2] Adding comparison column with value:", .y))
          .x <- .x |> mutate({{ comparison_column }} := .y)
        }
        .x
      }) |>
      bind_rows() |>
      mutate({{ facet_column }} := description)

    message("      get_row_binded_table: output table created")
    message(paste("      get_row_binded_table: output dims =", nrow(output), "x", ncol(output)))
    message("      get_row_binded_table: Checking comparison column values:")
    print(unique(output[[comparison_column]]))

    # This separator logic assumes a specific format like 'ContrastName=ExpressionType'
    # in the comparison column. We need to make this conditional as well.
    # Check if any values in the comparison column contain '=' before trying to separate.
    has_delimiter <- any(grepl("=", output[[comparison_column]]))
    message(paste("      get_row_binded_table: Any values contain '=' ?", has_delimiter))

    if (has_delimiter) {
      message("      get_row_binded_table: Separating on '=' delimiter")
      output <- output |>
        separate_wider_delim(
          {{ comparison_column }},
          delim = "=",
          names = c(comparison_column, expression_column)
        )
      message("      get_row_binded_table: Separation complete")
    } else {
      message("      get_row_binded_table: No '=' delimiter found, skipping separation")
    }

    message("   --- Exiting get_row_binded_table (DEBUG66) ---")
    return(output)
  }

  message("   getSignificantData: About to call get_row_binded_table for each list element...")

  logfc_tbl_all <- purrr::map2(
    list_of_da_tables, list_of_descriptions,
    function(a, b) {
      get_row_binded_table(de_table_list = a, description = b)
    }
  ) |>
    bind_rows()

  message("   getSignificantData: All tables bound")
  message(paste("   getSignificantData: logfc_tbl_all dims =", nrow(logfc_tbl_all), "x", ncol(logfc_tbl_all)))

  # CRITICAL FIX: Ensure expression_column exists if it doesn't
  expr_col_name <- expression_column
  if (!expr_col_name %in% colnames(logfc_tbl_all)) {
    message(sprintf("   getSignificantData: expression_column '%s' not found, creating dummy", expr_col_name))
    logfc_tbl_all[[expr_col_name]] <- NA_real_
  }

  selected_data <- logfc_tbl_all |>
    mutate({{ log_q_value_column }} := -log10(fdr_qvalue)) |>
    dplyr::select(
      {{ row_id }}, {{ log_q_value_column }}, {{ q_value_column }}, {{ p_value_column }}, {{ log_fc_column }},
      {{ comparison_column }}, {{ expression_column }},
      {{ facet_column }}
    ) |>
    dplyr::mutate(colour = case_when(
      abs({{ log_fc_column }}) >= 1 & {{ q_value_column }} >= q_val_thresh ~ "orange",
      abs({{ log_fc_column }}) >= 1 & {{ q_value_column }} < q_val_thresh ~ "purple",
      abs({{ log_fc_column }}) < 1 & {{ q_value_column }} < q_val_thresh ~ "blue",
      TRUE ~ "black"
    )) |>
    dplyr::mutate(colour = factor(colour, levels = c("black", "orange", "blue", "purple")))

  message("--- Exiting getSignificantData (DEBUG66) ---")
  selected_data
}

# ----------------------------------------------------------------------------
# getTypeOfGrouping
# ----------------------------------------------------------------------------
#' Assign experimental group list
#' @param design_matrix A data frame representing the design matrix.
#' @param group_id A string representing the name of the group ID column used in the design matrix.
#' @param sample_id A string representing the name of the sample ID column used in the design matrix.
#' @return A list where each element name is the name of a treatment group and each element is a vector containing the sample IDs within the treatment group.
#' @export
getTypeOfGrouping <- function(design_matrix, group_id, sample_id) {
  temp_type_of_grouping <- design_matrix |>
    dplyr::select(!!rlang::sym(group_id), !!rlang::sym(sample_id)) |>
    group_by(!!rlang::sym(group_id)) |>
    summarise(!!rlang::sym(sample_id) := list(!!rlang::sym(sample_id))) |>
    ungroup()

  type_of_grouping <- temp_type_of_grouping |> dplyr::pull(!!rlang::sym(sample_id))
  names(type_of_grouping) <- temp_type_of_grouping |> dplyr::pull(!!rlang::sym(group_id))

  return(type_of_grouping)
}

# ----------------------------------------------------------------------------
# extractResults
# ----------------------------------------------------------------------------
#' @export
extractResults <- function(results_list) {
  extracted <- purrr::map(results_list, \(x){
    x$results
  })

  names(extracted) <- names(results_list)

  return(extracted)
}

# ----------------------------------------------------------------------------
# outputDaResultsAllContrasts
# ----------------------------------------------------------------------------
#' @export
setMethod(
  f = "outputDaResultsAllContrasts",
  signature = "ProteinQuantitativeData",
  definition = function(theObject,
                        da_results_list_all_contrasts = NULL,
                        uniprot_tbl = NULL,
                        da_output_dir = NULL,
                        publication_graphs_dir = NULL,
                        file_prefix = "da_proteins",
                        args_row_id = NULL,
                        gene_names_column = "gene_names",
                        uniprot_id_column = "Entry") {
    message("--- Entering outputDaResultsAllContrasts ---")
    message(sprintf("   outputDaResultsAllContrasts: da_output_dir = %s", da_output_dir))
    message(sprintf("   outputDaResultsAllContrasts: file_prefix = %s", file_prefix))

    # Extract parameters from S4 object with fallbacks
    uniprot_tbl <- checkParamsObjectFunctionSimplify(theObject, "uniprot_tbl", uniprot_tbl)
    da_output_dir <- checkParamsObjectFunctionSimplify(theObject, "da_output_dir", da_output_dir)
    publication_graphs_dir <- checkParamsObjectFunctionSimplify(theObject, "publication_graphs_dir", publication_graphs_dir)
    args_row_id <- checkParamsObjectFunctionSimplify(theObject, "args_row_id", args_row_id)
    gene_names_column <- checkParamsObjectFunctionSimplify(theObject, "gene_names_column", gene_names_column)
    uniprot_id_column <- checkParamsObjectFunctionSimplify(theObject, "uniprot_id_column", uniprot_id_column)
 
    # CRITICAL FIX: Normalize paths for Windows (fixes C:// double slash issue)
    # Only normalize if the path is not NULL and doesn't already exist
    if (!is.null(da_output_dir)) {
      # Use normalizePath with mustWork=FALSE to handle non-existent dirs
      da_output_dir <- gsub("//+", "/", da_output_dir) # Remove double slashes first
      da_output_dir <- normalizePath(da_output_dir, winslash = "/", mustWork = FALSE)
      message(sprintf("   outputDaResultsAllContrasts: Normalized da_output_dir = %s", da_output_dir))
    }

    if (!is.null(publication_graphs_dir)) {
      publication_graphs_dir <- gsub("//+", "/", publication_graphs_dir) # Remove double slashes first
      publication_graphs_dir <- normalizePath(publication_graphs_dir, winslash = "/", mustWork = FALSE)
      message(sprintf("   outputDaResultsAllContrasts: Normalized publication_graphs_dir = %s", publication_graphs_dir))
    }

    # Ensure output directory exists (CRITICAL FIX!)
    if (!dir.exists(da_output_dir)) {
      dir.create(da_output_dir, recursive = TRUE, showWarnings = FALSE)
      message(sprintf("   outputDaResultsAllContrasts: Created output directory: %s", da_output_dir))
    }

    message(sprintf("   outputDaResultsAllContrasts: Processing %d contrasts", length(da_results_list_all_contrasts)))

    # Write results for each contrast with UNIQUE filenames
    for (contrast_name in names(da_results_list_all_contrasts)) {
      message(sprintf("   outputDaResultsAllContrasts: Writing files for contrast: %s", contrast_name))

      contrast_result <- da_results_list_all_contrasts[[contrast_name]]

      if (!is.null(contrast_result$da_proteins_long)) {
        # Clean contrast name for safe filename
        safe_contrast_name <- gsub("[^A-Za-z0-9_-]", "_", contrast_name)

        # Create annotated version (FIXED: proper conditional logic)
        da_proteins_long_annot <- contrast_result$da_proteins_long |>
          mutate(uniprot_acc_cleaned = str_split(!!sym(args_row_id), "-") |>
            purrr::map_chr(1))

        # Add UniProt annotations if available
        if (!is.null(uniprot_tbl) && nrow(uniprot_tbl) > 0) {
          da_proteins_long_annot <- da_proteins_long_annot |>
            left_join(uniprot_tbl, by = join_by(uniprot_acc_cleaned == !!sym(uniprot_id_column))) |>
            dplyr::select(-uniprot_acc_cleaned)
        } else {
          da_proteins_long_annot <- da_proteins_long_annot |>
            dplyr::select(-uniprot_acc_cleaned)
        }

        # Add gene_name column with proper conditional logic
        if ("gene_names" %in% names(da_proteins_long_annot)) {
          da_proteins_long_annot <- da_proteins_long_annot |>
            mutate(gene_name = purrr::map_chr(gene_names, \(x){
              tryCatch(
                {
                  if (is.na(x) || is.null(x) || x == "") {
                    ""
                  } else {
                    split_result <- str_split(x, " |:")[[1]]
                    if (length(split_result) > 0) split_result[1] else ""
                  }
                },
                error = function(e) ""
              )
            }))
        } else if (gene_names_column %in% names(da_proteins_long_annot)) {
          da_proteins_long_annot <- da_proteins_long_annot |>
            mutate(gene_name = purrr::map_chr(.data[[gene_names_column]], \(x){
              tryCatch(
                {
                  if (is.na(x) || is.null(x) || x == "") {
                    ""
                  } else {
                    split_result <- str_split(x, " |:")[[1]]
                    if (length(split_result) > 0) split_result[1] else ""
                  }
                },
                error = function(e) ""
              )
            }))
        } else {
          da_proteins_long_annot <- da_proteins_long_annot |>
            mutate(gene_name = "")
        }

        # Relocate gene_name column
        da_proteins_long_annot <- da_proteins_long_annot |>
          relocate(gene_name, .after = !!sym(args_row_id))

        # CRITICAL FIX: Use contrast-specific filename!
        contrast_filename <- paste0(file_prefix, "_", safe_contrast_name, "_long_annot.tsv")
        output_path <- file.path(da_output_dir, contrast_filename)

        message(sprintf("   outputDaResultsAllContrasts: Writing %s", output_path))

        # Write the file
        vroom::vroom_write(da_proteins_long_annot, output_path)

        # Verify file was written
        if (file.exists(output_path)) {
          file_size <- file.size(output_path)
          message(sprintf(
            "   outputDaResultsAllContrasts: SUCCESS - File written: %s (%d bytes)",
            contrast_filename, file_size
          ))
        } else {
          message(sprintf("   outputDaResultsAllContrasts: ERROR - File NOT written: %s", output_path))
        }

        # Also write Excel version
        excel_filename <- paste0(file_prefix, "_", safe_contrast_name, "_long_annot.xlsx")
        excel_path <- file.path(da_output_dir, excel_filename)
        writexl::write_xlsx(da_proteins_long_annot, excel_path)

        message(sprintf("   outputDaResultsAllContrasts: Also wrote Excel: %s", excel_filename))

        # [OK] NEW: Generate volcano plot for this contrast
        message(sprintf("   outputDaResultsAllContrasts: Generating volcano plot for contrast: %s", contrast_name))

        # Extract parameters - try direct access first, then use checkParams
        da_q_val_thresh <- if (!is.null(theObject@args$outputDaResultsAllContrasts$da_q_val_thresh)) {
          theObject@args$outputDaResultsAllContrasts$da_q_val_thresh
        } else {
          checkParamsObjectFunctionSimplify(theObject, "da_q_val_thresh", 0.05)
        }

        fdr_column <- if (!is.null(theObject@args$outputDaResultsAllContrasts$fdr_column)) {
          theObject@args$outputDaResultsAllContrasts$fdr_column
        } else {
          checkParamsObjectFunctionSimplify(theObject, "fdr_column", "fdr_qvalue")
        }

        log2fc_column <- if (!is.null(theObject@args$outputDaResultsAllContrasts$log2fc_column)) {
          theObject@args$outputDaResultsAllContrasts$log2fc_column
        } else {
          checkParamsObjectFunctionSimplify(theObject, "log2fc_column", "log2FC")
        }

        # Prepare data for volcano plot (same logic as in differentialAbundanceAnalysisHelper)
        volcano_data <- contrast_result$da_proteins_long |>
          mutate(lqm = -log10(!!sym(fdr_column))) |>
          dplyr::mutate(label = case_when(
            !!sym(fdr_column) < da_q_val_thresh ~ "Significant",
            TRUE ~ "Not sig."
          )) |>
          dplyr::mutate(colour = case_when(
            !!sym(fdr_column) < da_q_val_thresh ~ "purple",
            TRUE ~ "black"
          )) |>
          dplyr::mutate(colour = factor(colour, levels = c("black", "purple")))

        # Generate the volcano plot
        volcano_plot <- plotOneVolcanoNoVerticalLines(
          volcano_data,
          contrast_name,
          log_q_value_column = lqm,
          log_fc_column = !!sym(log2fc_column)
        )

        # Create Volcano_Plots directory if it doesn't exist
        volcano_dir <- file.path(publication_graphs_dir, "Volcano_Plots")
        if (!dir.exists(volcano_dir)) {
          dir.create(volcano_dir, recursive = TRUE, showWarnings = FALSE)
          message(sprintf("   outputDaResultsAllContrasts: Created volcano plots directory: %s", volcano_dir))
        }

        # Save volcano plot with contrast-specific filename
        volcano_png_path <- file.path(volcano_dir, paste0(safe_contrast_name, ".png"))
        volcano_pdf_path <- file.path(volcano_dir, paste0(safe_contrast_name, ".pdf"))

        # Save as PNG
        tryCatch(
          {
            ggplot2::ggsave(volcano_png_path, volcano_plot, width = 7, height = 7, dpi = 300)
            message(sprintf("   outputDaResultsAllContrasts: Saved volcano plot PNG: %s", basename(volcano_png_path)))
          },
          error = function(e) {
            message(sprintf("   outputDaResultsAllContrasts: ERROR saving PNG: %s", e$message))
          }
        )

        # Save as PDF
        tryCatch(
          {
            ggplot2::ggsave(volcano_pdf_path, volcano_plot, width = 7, height = 7)
            message(sprintf("   outputDaResultsAllContrasts: Saved volcano plot PDF: %s", basename(volcano_pdf_path)))
          },
          error = function(e) {
            message(sprintf("   outputDaResultsAllContrasts: ERROR saving PDF: %s", e$message))
          }
        )
      }
    }

    # [OK] NEW: Generate combined multi-page PDF with all volcano plots
    message("   outputDaResultsAllContrasts: Creating combined volcano plots PDF...")

    volcano_dir <- file.path(publication_graphs_dir, "Volcano_Plots")
    combined_pdf_path <- file.path(volcano_dir, "all_volcano_plots_combined.pdf")

    tryCatch(
      {
        # Re-extract parameters for combined PDF generation (in case they were modified)
        da_q_val_thresh <- if (!is.null(theObject@args$outputDaResultsAllContrasts$da_q_val_thresh)) {
          theObject@args$outputDaResultsAllContrasts$da_q_val_thresh
        } else {
          0.05
        }

        fdr_column <- if (!is.null(theObject@args$outputDaResultsAllContrasts$fdr_column)) {
          theObject@args$outputDaResultsAllContrasts$fdr_column
        } else {
          "fdr_qvalue"
        }

        log2fc_column <- if (!is.null(theObject@args$outputDaResultsAllContrasts$log2fc_column)) {
          theObject@args$outputDaResultsAllContrasts$log2fc_column
        } else {
          "log2FC"
        }

        # Collect all volcano plots
        all_volcano_plots <- list()

        for (contrast_name in names(da_results_list_all_contrasts)) {
          contrast_result <- da_results_list_all_contrasts[[contrast_name]]

          if (!is.null(contrast_result$da_proteins_long)) {
            # Recreate the volcano plot
            volcano_data <- contrast_result$da_proteins_long |>
              mutate(lqm = -log10(!!sym(fdr_column))) |>
              dplyr::mutate(label = case_when(
                !!sym(fdr_column) < da_q_val_thresh ~ "Significant",
                TRUE ~ "Not sig."
              )) |>
              dplyr::mutate(colour = case_when(
                !!sym(fdr_column) < da_q_val_thresh ~ "purple",
                TRUE ~ "black"
              )) |>
              dplyr::mutate(colour = factor(colour, levels = c("black", "purple")))

            volcano_plot <- plotOneVolcanoNoVerticalLines(
              volcano_data,
              contrast_name,
              log_q_value_column = lqm,
              log_fc_column = !!sym(log2fc_column)
            )

            all_volcano_plots[[contrast_name]] <- volcano_plot
          }
        }

        # Generate multi-page PDF
        if (length(all_volcano_plots) > 0) {
          pdf(file = combined_pdf_path, width = 7, height = 7, onefile = TRUE)
          purrr::walk(all_volcano_plots, print)
          invisible(dev.off())
          message(sprintf(
            "   outputDaResultsAllContrasts: Created combined PDF with %d volcano plots: %s",
            length(all_volcano_plots), basename(combined_pdf_path)
          ))
        }
      },
      error = function(e) {
        message(sprintf("   outputDaResultsAllContrasts: ERROR creating combined PDF: %s", e$message))
      }
    )

    # ============================================================================
    # NEW: Aggregate and save Num Sig DE Molecules from all contrasts
    # ============================================================================
    message("   outputDaResultsAllContrasts: Processing NumSigDaMolecules figures...")

    tryCatch(
      {
        # Create NumSigDaMolecules directory
        numsigde_dir <- file.path(publication_graphs_dir, "NumSigDaMolecules")
        if (!dir.exists(numsigde_dir)) {
          dir.create(numsigde_dir, recursive = TRUE, showWarnings = TRUE)
          if (!dir.exists(numsigde_dir)) {
            stop("Failed to create NumSigDaMolecules directory: ", numsigde_dir)
          }
          message(sprintf("   outputDaResultsAllContrasts: Created NumSigDaMolecules directory: %s", numsigde_dir))
        }

        # Collect num_sig_da_molecules tables from all contrasts
        all_numsig_tables <- list()

        for (contrast_name in names(da_results_list_all_contrasts)) {
          contrast_result <- da_results_list_all_contrasts[[contrast_name]]

          # Check for num_sig_da_molecules_first_go in the result
          if (!is.null(contrast_result$num_sig_da_molecules_first_go)) {
            numsig_data <- contrast_result$num_sig_da_molecules_first_go

            # Extract the table component
            if (!is.null(numsig_data$table)) {
              all_numsig_tables[[contrast_name]] <- numsig_data$table
              message(sprintf("   outputDaResultsAllContrasts: Found NumSigDE table for contrast: %s", contrast_name))
            }
          }
        }

        # If we have any tables, aggregate and save them
        if (length(all_numsig_tables) > 0) {
          # Combine all tables
          combined_numsig_table <- dplyr::bind_rows(all_numsig_tables)

          # Write the combined table
          table_path <- file.path(numsigde_dir, paste0(file_prefix, "_num_sig_da_molecules.tab"))
          vroom::vroom_write(combined_numsig_table, table_path)
          message(sprintf("   outputDaResultsAllContrasts: Wrote NumSigDE table: %s", basename(table_path)))

          # Also write as Excel
          excel_path <- file.path(numsigde_dir, paste0(file_prefix, "_num_sig_da_molecules.xlsx"))
          writexl::write_xlsx(combined_numsig_table, excel_path)
          message(sprintf("   outputDaResultsAllContrasts: Wrote NumSigDE Excel: %s", basename(excel_path)))

          # Generate combined barplot from the aggregated data
          # Filter to only significant results
          sig_only <- combined_numsig_table |>
            dplyr::filter(status != "Not significant")

          if (nrow(sig_only) > 0) {
            num_sig_de_barplot <- sig_only |>
              ggplot2::ggplot(ggplot2::aes(x = status, y = counts)) +
              ggplot2::geom_bar(stat = "identity", fill = "steelblue") +
              ggplot2::geom_text(stat = "identity", ggplot2::aes(label = counts), vjust = -0.5) +
              ggplot2::theme_bw() +
              ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)) +
              ggplot2::labs(
                title = "Number of Significant DE Molecules",
                x = "Direction",
                y = "Count"
              )

            # Add faceting if we have comparison column
            if ("comparison" %in% names(sig_only)) {
              num_sig_de_barplot <- num_sig_de_barplot +
                ggplot2::facet_wrap(~comparison, scales = "free_x")
            }

            # Calculate appropriate width based on number of comparisons
            num_comparisons <- length(unique(sig_only$comparison))
            plot_width <- max(7, (num_comparisons + 2) * 7 / 6)

            # Save as PNG
            png_path <- file.path(numsigde_dir, paste0(file_prefix, "_num_sig_da_molecules.png"))
            ggplot2::ggsave(png_path, num_sig_de_barplot, width = plot_width, height = 6, dpi = 300)
            message(sprintf("   outputDaResultsAllContrasts: Saved NumSigDE barplot PNG: %s", basename(png_path)))

            # Save as PDF
            pdf_path <- file.path(numsigde_dir, paste0(file_prefix, "_num_sig_da_molecules.pdf"))
            ggplot2::ggsave(pdf_path, num_sig_de_barplot, width = plot_width, height = 6)
            message(sprintf("   outputDaResultsAllContrasts: Saved NumSigDE barplot PDF: %s", basename(pdf_path)))
          } else {
            message("   outputDaResultsAllContrasts: No significant DE molecules found - skipping barplot")
          }
        } else {
          message("   outputDaResultsAllContrasts: No NumSigDE data found in contrast results - skipping")
        }
      },
      error = function(e) {
        message(sprintf("   outputDaResultsAllContrasts: ERROR processing NumSigDaMolecules: %s", e$message))
      }
    )

    message("--- Exiting outputDaResultsAllContrasts ---")
    return(TRUE)
  }
)

# ----------------------------------------------------------------------------
# getDaResultsLongFormat
# ----------------------------------------------------------------------------
# Get the differential abundance results in wide format
#' @export
setMethod(
  f = "getDaResultsLongFormat",
  signature = "list",
  definition = function(objectsList) {
    return_object_list <- purrr::map(objectsList, function(object) {
      # Correctly access the metabolite data from the nested 'theObject' slot.
      # This defensively handles cases where the slot might hold a list of
      # data frames (correct) or a single data frame (incorrect but handled).
      counts_data_slot <- object@theObject@metabolite_data
      counts_table_to_use <- if (is.list(counts_data_slot) && !is.data.frame(counts_data_slot)) {
        counts_data_slot[[1]]
      } else {
        counts_data_slot
      }

      id_col_name <- object@theObject@metabolite_id_column

      # Bind the list of data frames into a single tidy data frame
      tidy_results <- object@contrasts_results_table |>
        dplyr::bind_rows(.id = "comparison") |>
        dplyr::mutate(comparision_short = str_split_i(comparison, "=", 1))

      # Pivot the tidy data frame to a wide format using the correct ID column
      long_results <- tidy_results |>
        # Join with original counts using the correct ID column.
        dplyr::left_join(counts_table_to_use, by = join_by(!!sym(id_col_name) == !!sym(id_col_name))) |>
        dplyr::distinct()

      print(head(long_results))

      # Assign to the correct slot and return the object
      object@results_table_long <- long_results
      return(object)
    })

    return(return_object_list)
  }
)

# ----------------------------------------------------------------------------
# getDaResultsWideFormat
# ----------------------------------------------------------------------------
#' @export
setMethod(
  f = "getDaResultsWideFormat",
  signature = "list",
  definition = function(
    objectsList,
    qvalue_column = "fdr_qvalue",
    raw_pvalue_column = "raw_pvalue",
    log2fc_column = "logFC"
  ) {
    return_object_list <- purrr::map(objectsList, function(object) {
      # Correctly access the metabolite data from the nested 'theObject' slot.
      # This defensively handles cases where the slot might hold a list of
      # data frames (correct) or a single data frame (incorrect but handled).
      counts_data_slot <- object@theObject@metabolite_data
      counts_table_to_use <- if (is.list(counts_data_slot) && !is.data.frame(counts_data_slot)) {
        counts_data_slot[[1]]
      } else {
        counts_data_slot
      }

      id_col_name <- object@theObject@metabolite_id_column

      # Bind the list of data frames into a single tidy data frame
      tidy_results <- object@contrasts_results_table |>
        dplyr::bind_rows(.id = "comparison") |>
        dplyr::mutate(comparision_short = str_split_i(comparison, "=", 1))

      # Pivot the tidy data frame to a wide format using the correct ID column
      wida_results <- tidy_results |>
        tidyr::pivot_wider(
          id_cols = c(!!sym(id_col_name)),
          names_from = c(comparision_short),
          names_sep = ":",
          values_from = c(!!sym(log2fc_column), !!sym(qvalue_column), !!sym(raw_pvalue_column))
        ) |>
        # Join with original counts using the correct ID column.
        dplyr::left_join(counts_table_to_use, by = join_by(!!sym(id_col_name) == !!sym(id_col_name))) |>
        dplyr::arrange(dplyr::across(matches(qvalue_column))) |>
        dplyr::distinct()

      # Assign to the correct slot and return the object
      object@results_table_wide <- wida_results
      return(object)
    })

    return(return_object_list)
  }
)

