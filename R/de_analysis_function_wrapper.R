#' @title Wrapper for Differential Expression Analysis
#' @description A comprehensive wrapper function that orchestrates a full differential
#' expression (DE) analysis workflow. It handles parameter retrieval, data
#' preprocessing, QC plot generation (PCA, RLE), limma-based DE analysis,
#' and generation of results tables and plots.
#'
#' @details
#' This function acts as a high-level pipeline for DE analysis. The workflow is as follows:
#' 1.  It retrieves analysis parameters (e.g., `contrasts_tbl`, `formula_string`)
#'     from the input S4 object (`theObject`) or uses function arguments as overrides.
#' 2.  It preprocesses group names in the design matrix to ensure they are valid R names
#'     (e.g., by prepending "grp_" to names starting with numbers).
#' 3.  It generates QC plots: a Relative Log Expression (RLE) plot and Principal
#'     Component Analysis (PCA) plots (with and without sample labels).
#' 4.  It performs DE analysis using `limma` by calling `runTestsContrasts` with the
#'     specified formula and contrasts.
#' 5.  It processes the results to identify significant DE features and prepares data
#'     for volcano plots.
#' 6.  It generates various outputs and stores them in a list, including static
#'     volcano plots, p-value distribution histograms, and tables of DE feature counts.
#' 7.  It creates wide and long format data frames of the DE results.
#'
#' @param theObject An S4 object containing the quantification data (e.g., in
#'   `protein_quant_table`), design matrix, and analysis parameters in its slots.
#' @param contrasts_tbl A data frame or tibble defining the contrasts to be tested.
#'   If `NULL`, it's retrieved from `theObject`.
#' @param formula_string A character string representing the formula for the design
#'   matrix (e.g., "~ 0 + group"). If `NULL`, retrieved from `theObject`.
#' @param group_id The name of the column in the design matrix that defines the groups.
#'   If `NULL`, retrieved from `theObject`.
#' @param de_q_val_thresh The q-value (FDR) threshold for determining significance.
#'   If `NULL`, retrieved from `theObject`.
#' @param treat_lfc_cutoff The log-fold change threshold for `limma::treat`.
#'   A value of 0 implies a standard eBayes test. If `NULL`, retrieved from `theObject`.
#' @param eBayes_trend Logical. Whether to allow for a mean-variance trend in `limma::eBayes`.
#'   If `NULL`, retrieved from `theObject`.
#' @param eBayes_robust Logical. Whether to use robust estimation in `limma::eBayes`.
#'   If `NULL`, retrieved from `theObject`.
#' @param args_group_pattern A regex pattern to extract group information from sample names.
#'   If `NULL`, retrieved from `theObject`.
#' @param args_row_id The column name for feature/row identifiers (e.g., "uniprot_acc").
#'   If `NULL`, retrieved from `theObject`.
#' @param qvalue_column The name of the column containing q-values in the results.
#'   Defaults to "fdr_qvalue".
#' @param raw_pvalue_colum The name of the column containing raw p-values.
#'   Defaults to "raw_pvalue".
#'
#' @return A list containing the following elements:
#' \itemize{
#'   \item `theObject`: The updated S4 object.
#'   \item `rle_plot`: A ggplot object for the RLE plot.
#'   \item `pca_plot`: A ggplot object for the PCA plot without labels.
#'   \item `pca_plot_with_labels`: A ggplot object for the PCA plot with labels.
#'   \item `plot_num_of_values`: A plot showing the number of quantified values.
#'   \item `contrasts_results`: The raw results object from `runTestsContrasts`.
#'   \item `contrasts_results_table`: A data frame of the DE results.
#'   \item `significant_rows`: A data frame of significant results, formatted for plotting.
#'   \item `volplot_plot`: A combined ggplot object of all volcano plots.
#'   \item `num_sig_de_molecules_first_go`: A table counting significant features.
#'   \item `pvalhist`: A ggplot object of the p-value distribution histogram.
#'   \item `norm_counts`: A data frame of the normalized data.
#'   \item `de_proteins_wide`: A wide-format data frame of DE results.
#'   \item `de_proteins_long`: A long-format data frame of DE results.
#'   \item `list_of_volcano_plots`: A list of individual static volcano plots.
#'   \item `list_of_volcano_plots_with_gene_names`: A list of labeled static volcano plots.
#'   \item `num_sig_de_molecules`: A summary table of significant feature counts.
#'   \item `num_sig_de_genes_barplot_only_significant`: A bar plot of significant counts.
#'   \item `num_sig_de_genes_barplot_with_not_significant`: A bar plot of all feature counts.
#' }
#' @export
deAnalysisWrapperFunction <- function( theObject
                                       , contrasts_tbl = NULL
                                       , formula_string = NULL
                                       , group_id = NULL
                                       , de_q_val_thresh = NULL
                                       , treat_lfc_cutoff = NULL
                                       , eBayes_trend = NULL
                                       , eBayes_robust = NULL
                                       , args_group_pattern = NULL
                                       , args_row_id = NULL
                                       , qvalue_column = "fdr_qvalue"
                                       , raw_pvalue_colum = "raw_pvalue") {

  contrasts_tbl <- checkParamsObjectFunctionSimplify( theObject, "contrasts_tbl", NULL)
  formula_string <- checkParamsObjectFunctionSimplify( theObject, "formula_string", " ~ 0 + group")
  group_id <- checkParamsObjectFunctionSimplify( theObject, "group_id", "group")
  de_q_val_thresh <- checkParamsObjectFunctionSimplify( theObject, "de_q_val_thresh", 0.05)
  treat_lfc_cutoff <- checkParamsObjectFunctionSimplify( theObject, "treat_lfc_cutoff", 0)
  eBayes_trend <- checkParamsObjectFunctionSimplify( theObject, "eBayes_trend", TRUE)
  eBayes_robust <- checkParamsObjectFunctionSimplify( theObject, "eBayes_robust", TRUE)
  args_group_pattern <- checkParamsObjectFunctionSimplify( theObject, "args_group_pattern", "(\\d+)")
  args_row_id <- checkParamsObjectFunctionSimplify( theObject, "args_row_id", "uniprot_acc")

  # Add preprocessing for group names that start with numbers
  design_matrix <- theObject@design_matrix
  group_col <- design_matrix[["group"]]
  
  # Check if any group names start with numbers and create mapping
  starts_with_number <- grepl("^[0-9]", group_col)
  if(any(starts_with_number)) {
    original_groups <- unique(group_col)
    safe_groups <- purrr::map_chr(original_groups, \(x) {
      if(grepl("^[0-9]", x)) paste0("grp_", x) else x
    })
    group_mapping <- setNames(original_groups, safe_groups)
    
    # Update design matrix with safe names
    design_matrix[["group"]] <- purrr::map_chr(group_col, \(x) {
      if(grepl("^[0-9]", x)) paste0("grp_", x) else x
    })
    
    # Update contrasts table if it exists
    if(!is.null(contrasts_tbl)) {
      contrasts_tbl[[1]] <- purrr::map_chr(contrasts_tbl[[1]], \(x) {
        for(orig in names(group_mapping)) {
          x <- gsub(group_mapping[orig], orig, x, fixed = TRUE)
        }
        x
      })
    }
    
    theObject@design_matrix <- design_matrix
  }

  theObject <- updateParamInObject(theObject, "contrasts_tbl")
  theObject <- updateParamInObject(theObject, "formula_string")
  theObject <- updateParamInObject(theObject, "group_id")
  theObject <- updateParamInObject(theObject, "de_q_val_thresh")
  theObject <- updateParamInObject(theObject, "treat_lfc_cutoff")
  theObject <- updateParamInObject(theObject, "eBayes_trend")
  theObject <- updateParamInObject(theObject, "eBayes_robust")
  theObject <- updateParamInObject(theObject, "args_group_pattern")
  theObject <- updateParamInObject(theObject, "args_row_id")

  return_list <- list()
  return_list$theObject <- theObject

  ## plot RLE plot
  rle_plot <-   plotRle(theObject = theObject, theObject@group_id  ) +
    theme(axis.text.x = element_text(size = 13))   +
    theme(axis.text.y = element_text(size = 13))  +
    theme(axis.title.x = element_text(size = 12))  +
    theme(axis.title.y = element_text(size = 12))  +
    theme(plot.title = element_text(size = 12)) +
    theme(legend.text = element_text(size = 12)) +
    theme(legend.title = element_text(size = 12)) +
    xlab("Samples")

  return_list$rle_plot <- rle_plot

  ## plot PCA plot

  pca_plot <-  plotPca( theObject
                        , grouping_variable = theObject@group_id
                        , label_column = ""
                        , title = ""
                        , font_size = 8) +
    theme_bw() +
    theme(axis.text.x = element_text(size = 12)) +
    theme(axis.text.y = element_text(size = 12)) +
    theme(axis.title.x = element_text(size = 12)) +
    theme(axis.title.y = element_text(size = 12)) +
    theme(plot.title = element_text(size = 12)) +
    theme(legend.text = element_text(size = 12)) +
    theme(legend.title = element_text(size = 12))

  return_list$pca_plot <- pca_plot


  pca_plot_with_labels <-  plotPca( theObject
                                    , grouping_variable = theObject@group_id
                                    , label_column = theObject@sample_id
                                    , title = ""
                                    , font_size = 8) +
    theme_bw() +
    theme(axis.text.x = element_text(size = 12)) +
    theme(axis.text.y = element_text(size = 12)) +
    theme(axis.title.x = element_text(size = 12)) +
    theme(axis.title.y = element_text(size = 12)) +
    theme(plot.title = element_text(size = 12)) +
    theme(legend.text = element_text(size = 12)) +
    theme(legend.title = element_text(size = 12))

  return_list$pca_plot_with_labels <- pca_plot_with_labels

  ## Count the number of values
  return_list$plot_num_of_values <- plotNumOfValuesNoLog(theObject@protein_quant_table)

  ## Compare the different experimental groups and obtain lists of differentially expressed proteins.")

  rownames( theObject@design_matrix ) <- theObject@design_matrix |> dplyr::pull( one_of(theObject@sample_id ))


  contrasts_results <- runTestsContrasts(as.matrix(column_to_rownames(theObject@protein_quant_table, theObject@protein_id_column)),
                                         contrast_strings = contrasts_tbl[, 1][[1]],
                                         design_matrix = theObject@design_matrix,
                                         formula_string = formula_string,
                                         weights = NA,
                                         treat_lfc_cutoff = as.double(treat_lfc_cutoff),
                                         eBayes_trend = as.logical(eBayes_trend),
                                         eBayes_robust = as.logical(eBayes_robust))

  # Map back to original group names in results if needed
  if(exists("group_mapping")) {
    contrasts_results_table <- contrasts_results$results |>
      dplyr::mutate(comparison = purrr::map_chr(comparison, \(x) {
        result <- x
        for(safe_name in names(group_mapping)) {
          result <- gsub(safe_name, group_mapping[safe_name], result, fixed = TRUE)
        }
        result
      }))
  } else {
    contrasts_results_table <- contrasts_results$results
  }

  return_list$contrasts_results <- contrasts_results
  return_list$contrasts_results_table <- contrasts_results_table

  ## Prepare data for drawing the volcano plots

  significant_rows <- getSignificantData(list_of_de_tables = list(contrasts_results_table),
                                         list_of_descriptions = list("RUV applied"),
                                         row_id = !!sym(args_row_id),
                                         p_value_column = !!sym(raw_pvalue_colum),
                                         q_value_column = !!sym(qvalue_column),
                                         fdr_value_column = fdr_value_bh_adjustment,
                                         log_q_value_column = lqm,
                                         log_fc_column = logFC,
                                         comparison_column = "comparison",
                                         expression_column = "log_intensity",
                                         facet_column = analysis_type,
                                         q_val_thresh = de_q_val_thresh) |>
    dplyr::rename(log2FC = "logFC")

  return_list$significant_rows <- significant_rows

  # Print the volcano plots
  volplot_plot <- plotVolcano(significant_rows,
                              log_q_value_column = lqm,
                              log_fc_column = log2FC,
                              q_val_thresh = de_q_val_thresh,
                              formula_string = "analysis_type ~ comparison")


  return_list$volplot_plot <- volplot_plot

  ## Count the number of up or down significant differentially expressed proteins.
  num_sig_de_molecules_first_go <- printCountDeGenesTable(list_of_de_tables = list(contrasts_results_table),
                                                          list_of_descriptions = list( "RUV applied"),
                                                          formula_string = "analysis_type ~ comparison")

  return_list$num_sig_de_molecules_first_go <- num_sig_de_molecules_first_go

  ## Print p-values distribution figure
  pvalhist <- printPValuesDistribution(significant_rows,
                                       p_value_column = !!sym(raw_pvalue_colum),
                                       formula_string = "analysis_type ~ comparison")

  return_list$pvalhist <- pvalhist

  ## Create wide format output file
  norm_counts <- NA

  counts_table_to_use <- theObject@protein_quant_table

  norm_counts <- counts_table_to_use |>
    as.data.frame() |>
    column_to_rownames(args_row_id) |>
    set_colnames(paste0(colnames(counts_table_to_use[-1]), ".log2norm")) |>
    rownames_to_column(args_row_id)

  print(head( norm_counts))

  return_list$norm_counts <- norm_counts

  de_proteins_wide <- significant_rows |>
    dplyr::filter(analysis_type == "RUV applied") |>
    dplyr::select(-lqm, -colour, -analysis_type) |>
    pivot_wider(id_cols = c(!!sym(args_row_id)),
                names_from = c(comparison),
                names_sep = ":",
                values_from = c(log2FC, !!sym(qvalue_column), !!sym(raw_pvalue_colum))) |>
    left_join(counts_table_to_use, by = join_by( !!sym(args_row_id)  == !!sym(theObject@protein_id_column)  )   ) |>
    left_join(theObject@protein_id_table, by = join_by( !!sym(args_row_id) == !!sym(theObject@protein_id_column) )  ) |>
    dplyr::arrange(across(matches("!!sym(qvalue_column)"))) |>
    distinct()


  return_list$de_proteins_wide <- de_proteins_wide


  ## Create long format output file

  de_proteins_long <- createDeResultsLongFormat( lfc_qval_tbl = significant_rows |>
                                                   dplyr::filter(analysis_type == "RUV applied") ,
                                                 norm_counts_input_tbl = as.matrix(column_to_rownames(theObject@protein_quant_table, theObject@protein_id_column)),
                                                 raw_counts_input_tbl = as.matrix(column_to_rownames(theObject@protein_quant_table, theObject@protein_id_column)),
                                                 row_id = args_row_id,
                                                 sample_id = theObject@sample_id,
                                                 group_id = group_id,
                                                 group_pattern = args_group_pattern,
                                                 design_matrix_norm = theObject@design_matrix,
                                                 design_matrix_raw =  theObject@design_matrix,
                                                 ##POTENTIAL ISSUE
                                                 protein_id_table = theObject@protein_id_table)

  return_list$de_proteins_long <- de_proteins_long


  ## Plot static volcano plot
  static_volcano_plot_data <- de_proteins_long |>
    mutate( lqm = -log10(!!sym(qvalue_column) ) ) |>
    dplyr::mutate(label = case_when( !!sym(qvalue_column) < de_q_val_thresh ~ "Significant",
                                     TRUE ~ "Not sig.")) |>
    dplyr::mutate(colour = case_when( !!sym(qvalue_column) < de_q_val_thresh ~ "purple",
                                      TRUE ~ "black")) |>
    dplyr::mutate(colour = factor(colour, levels = c("black", "purple")))

  list_of_volcano_plots <- static_volcano_plot_data %>%
    group_by( comparison) %>%
    nest() %>%
    ungroup() %>%
    mutate( title = paste( comparison)) %>%
    mutate( plot = purrr::map2( data, title, \(x,y) { plotOneVolcanoNoVerticalLines(x, y
                                                                                     , log_q_value_column = lqm
                                                                                     , log_fc_column = log2FC) } ) )

  return_list$list_of_volcano_plots <- list_of_volcano_plots


  list_of_volcano_plots_with_gene_names <-  static_volcano_plot_data %>%
    group_by( comparison) %>%
    nest() %>%
    ungroup() %>%
    mutate( title = paste( comparison)) %>%
    mutate( plot = purrr::map( data, \(x) {
      printOneVolcanoPlotWithProteinLabel( input_table=  x
                                           , uniprot_table = uniprot_dat_cln |>
                                             mutate( gene_name = purrr::map_chr( gene_names
                                                                                 , \(x) str_split(x, "; ")[[1]][1])) )
       } ) )

  return_list$list_of_volcano_plots_with_gene_names <- list_of_volcano_plots_with_gene_names


  ## Return the number of significant molecules
  num_sig_de_molecules <- significant_rows %>%
    dplyr::mutate(status = case_when(!!sym(qvalue_column)  >= de_q_val_thresh ~ "Not significant",
                                     log2FC >= 0 & !!sym(qvalue_column) < de_q_val_thresh ~ "Significant and Up",
                                     log2FC < 0 &  !!sym(qvalue_column) < de_q_val_thresh ~ "Significant and Down",
                                     TRUE ~ "Not significant")) %>%
    group_by( comparison,  status) %>% # expression, analysis_type,
    summarise(counts = n()) %>%
    ungroup()

  formula_string <- ". ~ comparison"

  return_list$num_sig_de_molecules <- num_sig_de_molecules

  if (num_sig_de_molecules %>%
      dplyr::filter(status != "Not significant") |>
      nrow() > 0 ) {

    num_sig_de_genes_barplot_only_significant <- num_sig_de_molecules %>%
      dplyr::filter(status != "Not significant") %>%
      ggplot(aes(x = status, y = counts)) +
      geom_bar(stat = "identity") +
      geom_text(stat = 'identity', aes(label = counts), vjust = -0.5) +
      theme(axis.text.x = element_text(angle = 90))  +
      facet_wrap(as.formula(formula_string))

    num_of_comparison_only_significant <- num_sig_de_molecules |>
      distinct(comparison) |>
      nrow()

    return_list$num_sig_de_genes_barplot_only_significant <- num_sig_de_genes_barplot_only_significant
    return_list$num_of_comparison_only_significant <- num_of_comparison_only_significant
  }

  if (num_sig_de_molecules %>%
      dplyr::filter(status != "Not significant") |>
      nrow() > 0 ) {

    num_sig_de_genes_barplot_with_not_significant <- num_sig_de_molecules %>%
      ggplot(aes(x = status, y = counts)) +
      geom_bar(stat = "identity") +
      geom_text(stat = 'identity', aes(label = counts), vjust = -0.5) +
      theme(axis.text.x = element_text(angle = 90))  +
      facet_wrap(as.formula(formula_string))

    num_of_comparison_with_not_significant <- num_sig_de_molecules |>
      distinct(comparison) |>
      nrow()

    return_list$num_sig_de_genes_barplot_with_not_significant <- num_sig_de_genes_barplot_with_not_significant
    return_list$num_of_comparison_with_not_significant <- num_of_comparison_with_not_significant

  }

  return_list

}



## Create proteomics interactive volcano plot
#' @title Create and Save an Interactive Proteomics Volcano Plot
#' @description This function generates an interactive volcano plot using the Glimma
#' package and saves it as an HTML file. It's designed for proteomics data,
#' integrating differential expression results with UniProt annotations.
#'
#' @details The function joins the long-format DE results with a UniProt annotation
#' table. It then calculates values needed for the volcano plot (`-log10(FDR)`)
#' and assigns labels and colors based on significance and fold-change. Finally,
#' it iterates through each contrast (coefficient) in the `limma` fit object and
#' calls `getGlimmaVolcanoProteomics` to generate and save a self-contained
#' interactive HTML file for each one.
#'
#' @param de_proteins_long A data frame in long format containing DE results. Must
#'   include columns for log-fold change, FDR/q-value, and feature IDs.
#' @param uniprot_tbl A data frame with UniProt annotations, including at least
#'   a protein/uniprot accession column and a gene name column.
#' @param fit.eb An `MArrayLM` object returned by `limma::eBayes`.
#' @param publication_graphs_dir The base directory where the output plots will be
#'   saved, inside a subdirectory named "Interactive_Volcano_Plots".
#' @param args_row_id The column name in `de_proteins_long` that contains the
#'   primary feature identifiers to join with `uniprot_tbl`. Defaults to "uniprot_acc".
#' @param fdr_column The name of the column containing FDR or q-values. Defaults to "fdr_qvalue".
#' @param raw_p_value_column The name of the column with raw p-values. Defaults to "raw_pvalue".
#' @param log2fc_column The name of the column with log2 fold changes. Defaults to "log2FC".
#' @param de_q_val_thresh The significance threshold for FDR/q-value. Defaults to 0.05.
#' @param counts_tbl A matrix or data frame of normalized counts/intensities, with
#'   features as rows and samples as columns.
#' @param groups A vector specifying the group membership for each sample (column) in `counts_tbl`.
#' @param uniprot_id_column The column name in `uniprot_tbl` containing the UniProt accession
#'   numbers. Defaults to "Entry".
#' @param gene_names_column The column name in `uniprot_tbl` containing gene names.
#'   Defaults to "Gene Names".
#' @param display_columns A character vector of additional columns from the joined data
#'   to display in the interactive table. Defaults to "best_uniprot_acc".
#'
#' @return This function does not return a value. It is called for its side effect of
#'   writing HTML files to disk.
#' @export
writeInteractiveVolcanoPlotProteomics <- function( de_proteins_long
                                                   , uniprot_tbl
                                                   , fit.eb
                                                   , publication_graphs_dir
                                                   , args_row_id = "uniprot_acc"
                                                   , fdr_column = "fdr_qvalue"
                                                   , raw_p_value_column = "raw_pvalue"
                                                   , log2fc_column = "log2FC"
                                                   , de_q_val_thresh = 0.05
                                                   , counts_tbl = NULL
                                                   , groups = NULL
                                                   , uniprot_id_column = "Entry"
                                                   , gene_names_column = "Gene Names"
                                                   , display_columns = c( "best_uniprot_acc" )) {


  volcano_plot_tab <- de_proteins_long  |>
    dplyr::mutate(best_uniprot_acc = str_split(!!sym(args_row_id), ":" ) |> purrr::map_chr(1)  ) |>
    left_join(uniprot_tbl, by = join_by( best_uniprot_acc == !!sym( uniprot_id_column) ) ) |>
    dplyr::rename( UNIPROT_GENENAME = gene_names_column ) |>
    mutate( UNIPROT_GENENAME = purrr::map_chr( UNIPROT_GENENAME, \(x){str_split(x, " |:")[[1]][1]})) |>
    dplyr::mutate(gene_name = UNIPROT_GENENAME ) |>
    mutate( lqm = -log10(!!sym(fdr_column)))  |>
    dplyr::mutate(label = case_when(abs(!!sym(log2fc_column)) >= 1 & !!sym(fdr_column) >= de_q_val_thresh ~ "Not sig., logFC >= 1",
                                    abs(!!sym(log2fc_column)) >= 1 & !!sym(fdr_column) < de_q_val_thresh ~ "Sig., logFC >= 1",
                                    abs(!!sym(log2fc_column)) < 1 & !!sym(fdr_column) < de_q_val_thresh ~ "Sig., logFC < 1",
                                    TRUE ~ "Not sig.")) |>
    dplyr::mutate(colour = case_when(abs(!!sym(log2fc_column)) >= 1 & !!sym(fdr_column) >= de_q_val_thresh ~ "orange",
                                     abs(!!sym(log2fc_column)) >= 1 & !!sym(fdr_column) < de_q_val_thresh ~ "purple",
                                     abs(!!sym(log2fc_column)) < 1 & !!sym(fdr_column) < de_q_val_thresh ~ "blue",
                                     TRUE ~ "black")) |>
    dplyr::mutate(analysis_type = comparison)  |>
    dplyr::select( best_uniprot_acc, lqm, !!sym(fdr_column), !!sym(raw_p_value_column), !!sym(log2fc_column), comparison, label, colour,  gene_name
                   , any_of( display_columns ))   |>
    dplyr::mutate( my_alpha = case_when ( gene_name !=  "" ~ 1
                                          , TRUE ~ 0.5))

  output_dir <- file.path( publication_graphs_dir
                           ,  "Interactive_Volcano_Plots")

  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)


  purrr::walk( seq_len( ncol(fit.eb$coefficients))
               , \(coef) { # print(coef)

                 print( paste("run volcano plot", coef) )

                 #                    print(head( counts_tbl))
                 #                    print(groups)
                 #                    print( output_dir)

                 getGlimmaVolcanoProteomics( fit.eb
                                             , coef = coef
                                             , volcano_plot_tab  = volcano_plot_tab
                                             , uniprot_column = best_uniprot_acc
                                             , gene_name_column = gene_name
                                             , display_columns = display_columns
                                             , counts_tbl = counts_tbl
                                             , groups = groups
                                             , output_dir = output_dir ) } )

}




#' @title Create Interactive Proteomics Volcano Plot as a Widget
#' @description This function generates an interactive volcano plot using Glimma and
#' returns it as an HTML widget, suitable for embedding in R Markdown documents or
#' Shiny applications.
#'
#' @details This function is similar to `writeInteractiveVolcanoPlotProteomics` but
#' instead of writing HTML files to disk, it returns a list of widget objects.
#' It prepares the data by joining DE results with UniProt annotations and then
#' calls `getGlimmaVolcanoProteomicsWidget` for each contrast to generate the
#' interactive plots.
#'
#' @inheritParams writeInteractiveVolcanoPlotProteomics
#'
#' @return A list of HTML widget objects, where each element is an interactive
#'   volcano plot for one of the contrasts tested.
#' @export
writeInteractiveVolcanoPlotProteomicsWidget <- function( de_proteins_long
                                                         , uniprot_tbl
                                                         , fit.eb
                                                         , args_row_id = "uniprot_acc"
                                                         , fdr_column = "fdr_qvalue"
                                                         , raw_p_value_column = "raw_pvalue"
                                                         , log2fc_column = "log2FC"
                                                         , de_q_val_thresh = 0.05
                                                         , counts_tbl = NULL
                                                         , groups = NULL
                                                         , uniprot_id_column = "Entry"
                                                         , gene_names_column = "Gene Names"
                                                         , display_columns = c( "best_uniprot_acc" )) {


  volcano_plot_tab <- de_proteins_long  |>
    left_join(uniprot_tbl, by = join_by( !!sym(args_row_id) == !!sym( uniprot_id_column) ) ) |>
    dplyr::rename( UNIPROT_GENENAME = gene_names_column ) |>
    mutate( UNIPROT_GENENAME = purrr::map_chr( UNIPROT_GENENAME, \(x){str_split(x, " ")[[1]][1]})) |>
    mutate( lqm = -log10(!!sym(fdr_column)))  |>
    dplyr::mutate(label = case_when(abs(!!sym(log2fc_column)) >= 1 & !!sym(fdr_column) >= de_q_val_thresh ~ "Not sig., logFC >= 1",
                                    abs(!!sym(log2fc_column)) >= 1 & !!sym(fdr_column) < de_q_val_thresh ~ "Sig., logFC >= 1",
                                    abs(!!sym(log2fc_column)) < 1 & !!sym(fdr_column) < de_q_val_thresh ~ "Sig., logFC < 1",
                                    TRUE ~ "Not sig.")) |>
    dplyr::mutate(colour = case_when(abs(!!sym(log2fc_column)) >= 1 & !!sym(fdr_column) >= de_q_val_thresh ~ "orange",
                                     abs(!!sym(log2fc_column)) >= 1 & !!sym(fdr_column) < de_q_val_thresh ~ "purple",
                                     abs(!!sym(log2fc_column)) < 1 & !!sym(fdr_column) < de_q_val_thresh ~ "blue",
                                     TRUE ~ "black")) |>
    dplyr::mutate(gene_name = str_split(UNIPROT_GENENAME, " |:" ) |> purrr::map_chr(1)  ) |>
    dplyr::mutate(best_uniprot_acc = str_split(!!sym(args_row_id), ":" ) |> purrr::map_chr(1)  ) |>
    dplyr::mutate(analysis_type = comparison)  |>
    dplyr::select( best_uniprot_acc, lqm, !!sym(fdr_column), !!sym(raw_p_value_column), !!sym(log2fc_column), comparison, label, colour,  gene_name
                   , any_of( display_columns ))   |>
    dplyr::mutate( my_alpha = case_when ( gene_name !=  "" ~ 1
                                          , TRUE ~ 0.5))

  interactive_volcano_plots <- purrr::map(seq_len( ncol(fit.eb$coefficients))
                                          , \(coef) {
                                            # print(paste0( "coef = ", coef))
                                            getGlimmaVolcanoProteomicsWidget( fit.eb
                                                                              , coef = coef
                                                                              , volcano_plot_tab  = volcano_plot_tab
                                                                              , uniprot_column = best_uniprot_acc
                                                                              , gene_name_column = gene_name
                                                                              , display_columns = display_columns
                                                                              , counts_tbl = counts_tbl
                                                                              , groups = groups ) } )

  interactive_volcano_plots
}

#' @title Output and Save All Differential Expression Analysis Results
#' @description A wrapper function that takes the list of results generated by
#' `deAnalysisWrapperFunction` and writes all outputs—plots and data tables—to
#' disk in specified formats.
#'
#' @details This function is the primary outputter for the DE analysis pipeline.
#' It systematically saves:
#' - QC plots (PCA, RLE)
#' - `limma::plotSA` diagnostic plot
#' - The `MArrayLM` fit object from `eBayes`
#' - DE results in long and wide formats (as both .tsv and .xlsx)
#' - Annotated versions of the DE results tables
#' - Static volcano plots for each contrast, including a multi-page PDF
#' - P-value distribution histograms
#' - Bar plots summarizing the number of significant DE features
#' - Interactive volcano plots generated by `writeInteractiveVolcanoPlotProteomics`.
#'
#' @param de_analysis_results_list The list object returned by `deAnalysisWrapperFunction`.
#' @param theObject The main S4 object, used to retrieve file paths and parameters if not
#'   provided directly.
#' @param uniprot_tbl A data frame with UniProt annotations for annotating results.
#' @param de_output_dir Directory to save tabular results.
#' @param publication_graphs_dir Directory to save plot outputs.
#' @param file_prefix A prefix for all output filenames.
#' @param plots_format A character vector of file formats to save plots in (e.g., c("pdf", "png")).
#' @param args_row_id The column name for feature identifiers.
#' @param de_q_val_thresh The q-value significance threshold.
#' @param gene_names_column The column name for gene names in `uniprot_tbl`.
#' @param fdr_column The column name for FDR/q-values.
#' @param raw_p_value_column The column name for raw p-values.
#' @param log2fc_column The column name for log2 fold changes.
#' @param uniprot_id_column The column name for UniProt accessions in `uniprot_tbl`.
#' @param display_columns Additional columns to include in interactive plots.
#'
#' @return The function is called for its side effects and does not return a value.
#' @export
outputDeAnalysisResults <- function(de_analysis_results_list
                                    , theObject
                                    , uniprot_tbl
                                    , de_output_dir = NULL
                                    , publication_graphs_dir = NULL
                                    , file_prefix = NULL
                                    , plots_format = NULL
                                    , args_row_id = NULL
                                    , de_q_val_thresh = NULL
                                    , gene_names_column = NULL
                                    , fdr_column = NULL
                                    , raw_p_value_column = NULL
                                    , log2fc_column = NULL
                                    , uniprot_id_column = NULL
                                    , display_columns = NULL
) {


  uniprot_tbl <- checkParamsObjectFunctionSimplify(theObject, "uniprot_tbl", NULL)
  de_output_dir <- checkParamsObjectFunctionSimplify(theObject, "de_output_dir", NULL)
  publication_graphs_dir <- checkParamsObjectFunctionSimplify(theObject, "publication_graphs_dir", NULL)
  file_prefix <- checkParamsObjectFunctionSimplify(theObject, "file_prefix", "de_proteins")
  plots_format <- checkParamsObjectFunctionSimplify(theObject, "plots_format", c("pdf", "png"))
  args_row_id <- checkParamsObjectFunctionSimplify(theObject, "args_row_id", "uniprot_acc")
  de_q_val_thresh <- checkParamsObjectFunctionSimplify(theObject, "de_q_val_thresh", 0.05)
  gene_names_column <- checkParamsObjectFunctionSimplify(theObject, "gene_names_column", "Gene Names")
  fdr_column <- checkParamsObjectFunctionSimplify(theObject, "fdr_column", "fdr_qvalue")
  raw_p_value_column <- checkParamsObjectFunctionSimplify(theObject, "raw_p_value_column", "raw_pvalue")
  log2fc_column <- checkParamsObjectFunctionSimplify(theObject, "log2fc_column", "log2FC")
  uniprot_id_column <- checkParamsObjectFunctionSimplify(theObject, "uniprot_id_column", "Entry")
  display_columns <- checkParamsObjectFunctionSimplify(theObject, "display_columns", c( "best_uniprot_acc" ))

  theObject <- updateParamInObject(theObject, "uniprot_tbl")
  theObject <- updateParamInObject(theObject, "de_output_dir")
  theObject <- updateParamInObject(theObject, "publication_graphs_dir")
  theObject <- updateParamInObject(theObject, "file_prefix")
  theObject <- updateParamInObject(theObject, "plots_format")
  theObject <- updateParamInObject(theObject, "args_row_id")
  theObject <- updateParamInObject(theObject, "de_q_val_thresh")
  theObject <- updateParamInObject(theObject, "gene_names_column")
  theObject <- updateParamInObject(theObject, "fdr_column")
  theObject <- updateParamInObject(theObject, "raw_p_value_column")
  theObject <- updateParamInObject(theObject, "log2fc_column")
  theObject <- updateParamInObject(theObject, "uniprot_id_column")
  theObject <- updateParamInObject(theObject, "display_columns")

  ## PCA plot
  plot_pca_plot <- de_analysis_results_list$pca_plot

  dir.create(file.path( publication_graphs_dir, "PCA")
             , recursive = TRUE
             , showWarnings = FALSE)

  for( format_ext in plots_format) {
    file_name <- file.path( publication_graphs_dir, "PCA", paste0("PCA_plot.",format_ext))
    ggsave(filename = file_name, plot = plot_pca_plot, limitsize = FALSE)
  }

  plot_pca_plot_with_labels <- de_analysis_results_list$pca_plot_with_labels
  for( format_ext in plots_format) {
    file_name <- file.path( publication_graphs_dir, "PCA", paste0("PCA_plot_with_sample_ids.",format_ext))
    ggsave(filename = file_name, plot = plot_pca_plot_with_labels, limitsize = FALSE)
  }

  ## RLE plot
  plot_rle_plot <- de_analysis_results_list$rle_plot

  dir.create(file.path( publication_graphs_dir, "RLE")
             , recursive = TRUE
             , showWarnings = FALSE)

  for( format_ext in plots_format) {
    file_name <- file.path( publication_graphs_dir, "RLE", paste0("RLE_plot.",format_ext))
    ggsave(filename = file_name, plot = plot_rle_plot, limitsize = FALSE)
  }

  ## Save the number of values graph
  plot_num_of_values <- de_analysis_results_list$plot_num_of_values

  for( format_ext in plots_format) {
    file_name <- file.path(de_output_dir, paste0("num_of_values.",format_ext))
    ggsave(filename = file_name, plot = plot_num_of_values, limitsize = FALSE)
  }

  ## Contrasts results
  ## This plot is used to check the mean-variance relationship of the expression data, after fitting a linear model.
  pdf(file.path(de_output_dir, "plotSA_after_ruvIII.pdf" ))
  plotSA(de_analysis_results_list$contrasts_results$fit.eb)
  dev.off()

  png(file.path(de_output_dir, "plotSA_after_ruvIII.png" ))
  plotSA(de_analysis_results_list$contrasts_results$fit.eb)
  dev.off()

  saveRDS( de_analysis_results_list$contrasts_results$fit.eb,
           file.path(de_output_dir, "fit.eb.RDS" ) )

  ## Values for volcano plts

  ## Write all the results in one single table
  significant_rows <- de_analysis_results_list$significant_rows

  significant_rows |>
    dplyr::select(-colour, -lqm) |>
    vroom::vroom_write(file.path(de_output_dir, "lfc_qval_long.tsv"))

  significant_rows |>
    dplyr::select(-colour, -lqm) |>
    writexl::write_xlsx(file.path(de_output_dir, "lfc_qval_long.xlsx"))

  ## Print Volcano plot
  volplot_plot <- de_analysis_results_list$volplot_plot

  for( format_ext in plots_format) {
    file_name <- file.path(de_output_dir, paste0("volplot_gg_all.",format_ext))
    ggsave(filename = file_name, plot = volplot_plot, width = 7.29, height = 6)
  }


  ## Number of values graph
  plot_num_of_values <- de_analysis_results_list$plot_num_of_values

  for( format_ext in plots_format) {
    file_name <- file.path(de_output_dir, paste0("num_of_values.",format_ext))
    ggsave(filename = file_name, plot = plot_num_of_values, limitsize = FALSE)
  }

  ## Contrasts results
  ## This plot is used to check the mean-variance relationship of the expression data, after fitting a linear model.
  contrasts_results <- de_analysis_results_list$contrasts_results
  for( format_ext in plots_format) {
    file_name <- file.path(de_output_dir, paste0("plotSA_after_ruvIII",format_ext))

    if( format_ext == "pdf") {
      pdf(file_name)
    } else if(format_ext == "png") {
      png(file_name)
    }

    plotSA(contrasts_results$fit.eb)
    dev.off()

  }

  saveRDS( contrasts_results$fit.eb,
           file.path(de_output_dir, "fit.eb.RDS" ) )

  ## Values for volcano plts

  ## Write all the results in one single table
  significant_rows <- de_analysis_results_list$significant_rows
  significant_rows |>
    dplyr:::select(-colour, -lqm) |>
    vroom::vroom_write(file.path(de_output_dir, "lfc_qval_long.tsv"))

  significant_rows |>
    dplyr:::select(-colour, -lqm) |>
    writexl::write_xlsx(file.path(de_output_dir, "lfc_qval_long.xlsx"))


  ## Count the number of up or down significnat differentially expressed proteins.
  if( !is.null(de_analysis_results_list$num_sig_de_genes_barplot_only_significant)) {
    num_sig_de_genes_barplot_only_significant <- de_analysis_results_list$num_sig_de_genes_barplot_only_significant
    num_of_comparison_only_significant <- de_analysis_results_list$num_of_comparison_only_significant

    savePlot(num_sig_de_genes_barplot_only_significant,
             base_path = de_output_dir,
             plot_name = paste0(file_prefix, "_num_sda_entities_barplot_only_significant"),
             formats =  c("pdf", "png", "svg"),
             width = (num_of_comparison_only_significant + 2) *7/6,
             height = 6)



  }


  ## Count the number of up or down significnat differentially expressed proteins.
  num_sig_de_molecules_first_go <- de_analysis_results_list$num_sig_de_molecules_first_go
  vroom::vroom_write(num_sig_de_molecules_first_go$table,
                     file.path(de_output_dir,
                               paste0(file_prefix, "_num_significant_differentially_abundant_all.tab") ))

  writexl::write_xlsx(num_sig_de_molecules_first_go$table,
                      file.path(de_output_dir,
                                paste0(file_prefix, "_num_significant_differentially_abundant_all.xlsx")) )



  ## Print p-values distribution figure
  pvalhist <- de_analysis_results_list$pvalhist
  for( format_ext in plots_format) {
    file_name<-file.path(de_output_dir,paste0(file_prefix, "_p_values_distn.",format_ext))
    ggsave(filename = file_name,
           plot = pvalhist,
           height = 10,
           width = 7)

  }



  ## Create wide format output file
  de_proteins_wide <- de_analysis_results_list$de_proteins_wide
  vroom::vroom_write( de_proteins_wide,
                      file.path( de_output_dir,
                                 paste0(file_prefix, "_wide.tsv")))

  writexl::write_xlsx( de_proteins_wide,
                       file.path( de_output_dir,
                                  paste0(file_prefix, "_wide.xlsx")))

  de_proteins_wide_annot <- de_proteins_wide |>
    mutate( uniprot_acc_cleaned = str_split( !!sym(args_row_id), "-" )  |>
              purrr::map_chr(1) ) |>
    left_join( uniprot_tbl, by = join_by( uniprot_acc_cleaned == Entry ) ) |>
    dplyr::select( -uniprot_acc_cleaned)  |>
    mutate( gene_name = purrr::map_chr( !!sym(gene_names_column), \(x){str_split(x, " |:")[[1]][1]})) |>
    relocate( gene_name, .after = !!sym(args_row_id))

  vroom::vroom_write( de_proteins_wide_annot,
                      file.path( de_output_dir,
                                 paste0(file_prefix, "_wide_annot.tsv")))

  writexl::write_xlsx( de_proteins_wide_annot,
                       file.path( de_output_dir,
                                  paste0(file_prefix, "_wide_annot.xlsx")))

  ## Create long format output file
  de_proteins_long <- de_analysis_results_list$de_proteins_long
  vroom::vroom_write( de_proteins_long,
                      file.path( de_output_dir,
                                 paste0(file_prefix, "_long.tsv")))

  writexl::write_xlsx( de_proteins_long,
                       file.path( de_output_dir,
                                  paste0(file_prefix, "_long.xlsx")))

  de_proteins_long_annot <- de_proteins_long |>
    mutate( uniprot_acc_cleaned = str_split( !!sym(args_row_id), "-" )  |>
              purrr::map_chr(1) )|>
    left_join( uniprot_tbl, by = join_by( uniprot_acc_cleaned == Entry ) )  |>
    dplyr::select( -uniprot_acc_cleaned)  |>
    mutate( gene_name = purrr::map_chr( !!sym(gene_names_column), \(x){str_split(x, " |:")[[1]][1]})) |>
    relocate( gene_name, .after = !!sym(args_row_id))

  vroom::vroom_write( de_proteins_long_annot,
                      file.path( de_output_dir,
                                 paste0(file_prefix, "_long_annot.tsv")))

  writexl::write_xlsx( de_proteins_long_annot,
                       file.path( de_output_dir,
                                  paste0(file_prefix, "_long_annot.xlsx")))

  ## Static volcano plots
  dir.create(file.path( publication_graphs_dir, "Volcano_Plots")
             , recursive = TRUE
             , showWarnings = FALSE)

  list_of_volcano_plots <- de_analysis_results_list$list_of_volcano_plots
  
  # Print diagnostic info about the volcano plots
  message(sprintf("Number of volcano plots: %d", nrow(list_of_volcano_plots)))

  purrr::walk2( list_of_volcano_plots %>% dplyr::pull(title),
                list_of_volcano_plots %>% dplyr::pull(plot),
                \(x,y){
                # gg_save_logging ( .y, file_name_part, plots_format)

                savePlot( y
                          , base_path = file.path( publication_graphs_dir, "Volcano_Plots")
                          , plot_name =  x
                          , formats = plots_format, width = 7, height = 7)

                })

  # Generate a multi-page PDF with all volcano plots
  volcano_plots_list <- list_of_volcano_plots %>% dplyr::pull(plot)
  
  # Generate combined PDF with all plots, one per page
  pdf_file <- file.path(publication_graphs_dir, "Volcano_Plots", "list_of_volcano_plots.pdf")
  pdf(file = pdf_file, width = 7, height = 7, onefile = TRUE)
  purrr::walk(volcano_plots_list, print)
  invisible(dev.off())
  
  # Verify the PDF was created with the right number of pages
  message(sprintf("Created multi-page PDF at %s", pdf_file))

  list_of_volcano_plots_with_gene_names <- de_analysis_results_list$list_of_volcano_plots_with_gene_names
  
  # Print diagnostic info about the labeled volcano plots
  message(sprintf("Number of labeled volcano plots: %d", nrow(list_of_volcano_plots_with_gene_names)))

  purrr::walk2( list_of_volcano_plots_with_gene_names %>% dplyr::pull(title)
                , list_of_volcano_plots_with_gene_names %>% dplyr::pull(plot)
                , \(x, y) {

                  savePlot( x
                            , file.path( publication_graphs_dir, "Volcano_Plots")
                            , paste0( y,"_with_protein_labels"))
                })

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
  createDirIfNotExists(file.path(publication_graphs_dir, "NumSigDeMolecules"))
  vroom::vroom_write( de_analysis_results_list$num_sig_de_molecules,
                      file.path(publication_graphs_dir, "NumSigDeMolecules", paste0(file_prefix, "_num_sig_de_molecules.tab" ) ))


  if( !is.null(de_analysis_results_list$num_sig_de_genes_barplot_only_significant)) {
    num_sig_de_genes_barplot_only_significant <- de_analysis_results_list$num_sig_de_genes_barplot_only_significant
    num_of_comparison_only_significant <- de_analysis_results_list$num_of_comparison_only_significant

    savePlot(plot = num_sig_de_genes_barplot_only_significant,
             base_path = file.path(publication_graphs_dir, "NumSigDeMolecules"),
             plot_name = paste0(file_prefix, "_num_sig_de_molecules."),
             formats = plots_format,
             width = (num_of_comparison_only_significant + 2) *7/6,
             height = 6)


  }


  if(!is.null( de_analysis_results_list$num_sig_de_genes_barplot_with_not_significant )) {
    num_sig_de_genes_barplot_with_not_significant <- de_analysis_results_list$num_sig_de_genes_barplot_with_not_significant
    num_of_comparison_with_not_significant <- de_analysis_results_list$num_of_comparison_with_not_significant

    print("print bar plot")

    savePlot(num_sig_de_genes_barplot_with_not_significant,
             base_path = file.path(publication_graphs_dir, "NumSigDeMolecules"),
             plot_name = paste0(file_prefix, "_num_sig_de_molecules_with_not_significant"),
             formats = plots_format,
             width = (num_of_comparison_with_not_significant + 2) *7/6,
             height = 6)
  }

  ## Write interactive volcano plot
  counts_mat <- (de_analysis_results_list$theObject)@protein_quant_table |>
    column_to_rownames((de_analysis_results_list$theObject)@protein_id_column  ) |>
    as.matrix()

  this_design_matrix <- de_analysis_results_list$theObject@design_matrix

  rownames( this_design_matrix ) <- this_design_matrix[,de_analysis_results_list$theObject@sample_id]

  this_groups <- this_design_matrix[colnames( counts_mat), "group"]

  writeInteractiveVolcanoPlotProteomics( de_proteins_long
                                         , uniprot_tbl = uniprot_tbl
                                         , fit.eb = contrasts_results$fit.eb
                                         , publication_graphs_dir= publication_graphs_dir
                                         , args_row_id = args_row_id
                                         , fdr_column = fdr_column
                                         , raw_p_value_column = raw_p_value_column
                                         , log2fc_column = log2fc_column
                                         , de_q_val_thresh = de_q_val_thresh
                                         , counts_tbl = counts_mat

                                         , groups = this_groups
                                         , uniprot_id_column = uniprot_id_column
                                         , gene_names_column = gene_names_column
                                         , display_columns = display_columns )

}



#' @title Main Wrapper to Write Interactive Volcano Plot
#' @description A convenience wrapper around `writeInteractiveVolcanoPlotProteomics`.
#' Its primary purpose is to extract all necessary parameters from `theObject` and
#' the `de_analysis_results_list` before calling the main interactive plotting function.
#'
#' @details This function serves as a simplified entry point for generating interactive
#' volcano plots. It uses `checkParamsObjectFunctionSimplify` to pull most of its
#' arguments from slots within `theObject`, reducing the number of explicit arguments
#' needed in the function call. It's particularly useful in a sequential pipeline
#' where `theObject` has been populated with all necessary parameters.
#'
#' @inheritParams outputDeAnalysisResults
#'
#' @return The function is called for its side effects and does not return a value.
#' @export
writeInteractiveVolcanoPlotProteomicsMain <- function(de_analysis_results_list
                                                      , theObject
                                                      , uniprot_tbl
                                                      , publication_graphs_dir = NULL
                                                      , file_prefix = NULL
                                                      , plots_format = NULL
                                                      , args_row_id = NULL
                                                      , de_q_val_thresh = NULL
                                                      , gene_names_column = NULL
                                                      , fdr_column = NULL
                                                      , raw_p_value_column = NULL
                                                      , log2fc_column = NULL
                                                      , uniprot_id_column = NULL
                                                      , display_columns = NULL
) {


  uniprot_tbl <- checkParamsObjectFunctionSimplify(theObject, "uniprot_tbl", NULL)
  publication_graphs_dir <- checkParamsObjectFunctionSimplify(theObject, "publication_graphs_dir", NULL)
  args_row_id <- checkParamsObjectFunctionSimplify(theObject, "args_row_id", "uniprot_acc")
  de_q_val_thresh <- checkParamsObjectFunctionSimplify(theObject, "de_q_val_thresh", 0.05)
  gene_names_column <- checkParamsObjectFunctionSimplify(theObject, "gene_names_column", "Gene Names")
  fdr_column <- checkParamsObjectFunctionSimplify(theObject, "fdr_column", "fdr_qvalue")
  raw_p_value_column <- checkParamsObjectFunctionSimplify(theObject, "raw_p_value_column", "raw_pvalue")
  log2fc_column <- checkParamsObjectFunctionSimplify(theObject, "log2fc_column", "log2FC")
  uniprot_id_column <- checkParamsObjectFunctionSimplify(theObject, "uniprot_id_column", "Entry")
  display_columns <- checkParamsObjectFunctionSimplify(theObject, "display_columns", c( "best_uniprot_acc" ))

  theObject <- updateParamInObject(theObject, "uniprot_tbl")
  theObject <- updateParamInObject(theObject, "publication_graphs_dir")
  theObject <- updateParamInObject(theObject, "args_row_id")
  theObject <- updateParamInObject(theObject, "de_q_val_thresh")
  theObject <- updateParamInObject(theObject, "gene_names_column")
  theObject <- updateParamInObject(theObject, "fdr_column")
  theObject <- updateParamInObject(theObject, "raw_p_value_column")
  theObject <- updateParamInObject(theObject, "log2fc_column")
  theObject <- updateParamInObject(theObject, "uniprot_id_column")
  theObject <- updateParamInObject(theObject, "display_columns")


  ## Write interactive volcano plot

  de_proteins_long <- de_analysis_results_list$de_proteins_long
  contrasts_results <- de_analysis_results_list$contrasts_results

  counts_mat <- (de_analysis_results_list$theObject)@protein_quant_table |>
    column_to_rownames((de_analysis_results_list$theObject)@protein_id_column  ) |>
    as.matrix()

  this_design_matrix <- de_analysis_results_list$theObject@design_matrix

  rownames( this_design_matrix ) <- this_design_matrix[,de_analysis_results_list$theObject@sample_id]

  this_groups <- this_design_matrix[colnames( counts_mat), "group"]

  writeInteractiveVolcanoPlotProteomics( de_proteins_long
                                         , uniprot_tbl = uniprot_tbl
                                         , fit.eb = contrasts_results$fit.eb
                                         , publication_graphs_dir= publication_graphs_dir
                                         , args_row_id = args_row_id
                                         , fdr_column = fdr_column
                                         , raw_p_value_column = raw_p_value_column
                                         , log2fc_column = log2fc_column
                                         , de_q_val_thresh = de_q_val_thresh
                                         , counts_tbl = counts_mat

                                         , groups = this_groups
                                         , uniprot_id_column = uniprot_id_column
                                         , gene_names_column = gene_names_column
                                         , display_columns = display_columns )

}
