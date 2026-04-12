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
