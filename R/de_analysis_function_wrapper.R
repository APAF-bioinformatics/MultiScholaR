


#' @export
deAnalysisWrapperFunction <- function( theObject
                                       , contrasts_tbl
                                       , formula_string = " ~ 0 + group"
                                       , group_id = "group"
                                       , de_q_val_thresh = 0.05
                                       , treat_lfc_cutoff = 0
                                       , eBayes_trend = TRUE
                                       , eBayes_robust = TRUE
                                       , args_group_pattern = "(\\d+)"
                                       , args_row_id = "uniprot_acc" ) {

  return_list <- list()

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
                           , group_column = theObject@group_id
                           , label_column = ""
                           , title = ""
                           , geom_text_size = 8) +
    theme_bw() +
    theme(axis.text.x = element_text(size = 12)) +
    theme(axis.text.y = element_text(size = 12)) +
    theme(axis.title.x = element_text(size = 12)) +
    theme(axis.title.y = element_text(size = 12)) +
    theme(plot.title = element_text(size = 12)) +
    theme(legend.text = element_text(size = 12)) +
    theme(legend.title = element_text(size = 12))

  return_list$pca_plot <- pca_plot

  ## Count the number of values
  return_list$plot_num_of_values <- plotNumOfValuesNoLog(theObject@protein_data)

  ## Compare the different experimental groups and obtain lists of differentially expressed proteins.")

  rownames( theObject@design_matrix ) <- theObject@design_matrix |> pull( one_of(theObject@sample_id ))

  # requires statmod library
  contrasts_results <- runTestsContrasts(theObject@protein_data |> column_to_rownames(theObject@protein_id_column  ) |> as.matrix(),
                                         contrast_strings = contrasts_tbl[, 1][[1]],
                                         design_matrix = theObject@design_matrix,
                                         formula_string = formula_string,
                                         weights = NA,
                                         treat_lfc_cutoff = as.double(treat_lfc_cutoff),
                                         eBayes_trend = as.logical(eBayes_trend),
                                         eBayes_robust = as.logical(eBayes_robust))

  contrasts_results_table <- contrasts_results$results

  return_list$contrasts_results <- contrasts_results
  return_list$contrasts_results_table <- contrasts_results_table

  ## Prepare data for drawing the volcano plots

  significant_rows <- getSignificantData(list_of_de_tables = list(contrasts_results_table),
                                         list_of_descriptions = list("RUV applied"),
                                         row_id = !!sym(args_row_id),
                                         p_value_column = p.mod,
                                         q_value_column = q.mod,
                                         fdr_value_column = fdr.mod,
                                         log_q_value_column = lqm,
                                         log_fc_column = logFC,
                                         comparison_column = comparison,
                                         expression_column = log_intensity,
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

  ## Count the number of up or down significnat differentially expressed proteins.
  num_sig_de_molecules_first_go <- printCountDeGenesTable(list_of_de_tables = list(contrasts_results_table),
                                                 list_of_descriptions = list( "RUV applied"),
                                                 formula_string = "analysis_type ~ comparison")

  return_list$num_sig_de_molecules_first_go <- num_sig_de_molecules_first_go

  ## Print p-values distribution figure
  pvalhist <- printPValuesDistribution(significant_rows,
                                       p_value_column = p.mod,
                                       formula_string = "analysis_type ~ comparison")

  return_list$pvalhist <- pvalhist

  ## Create wide format output file
  norm_counts <- NA

  counts_table_to_use <- theObject@protein_data

  norm_counts <- counts_table_to_use |>
    as.data.frame() |>
    set_colnames(paste0(colnames(counts_table_to_use), ".log2norm")) |>
    rownames_to_column(args_row_id)

  return_list$norm_counts <- norm_counts

  de_proteins_wide <- significant_rows |>
    dplyr::filter(analysis_type == "RUV applied") |>
    dplyr::select(-lqm, -colour, -analysis_type) |>
    pivot_wider(id_cols = c(!!sym(args_row_id)),
                names_from = c(comparison),
                names_sep = ":",
                values_from = c(log2FC, q.mod, p.mod)) |>
    left_join(counts_table_to_use, by = join_by( !!sym(args_row_id)  == !!sym(theObject@protein_id_column)  )   ) |>
    dplyr::arrange(across(matches("q.mod"))) |>
    distinct()


  return_list$de_proteins_wide <- de_proteins_wide


  ## Create long format output file

  de_proteins_long <- createDeResultsLongFormat( lfc_qval_tbl = significant_rows |>
                                                   dplyr::filter(analysis_type == "RUV applied") ,
                                                 norm_counts_input_tbl = counts_table_to_use |> column_to_rownames( theObject@protein_id_column) |> as.matrix(),
                                                 raw_counts_input_tbl = counts_table_to_use |> column_to_rownames(theObject@protein_id_column) |> as.matrix(),
                                                 row_id = args_row_id,
                                                 sample_id = theObject@sample_id,
                                                 group_id = group_id,
                                                 group_pattern = args_group_pattern,
                                                 design_matrix_norm = theObject@design_matrix,
                                                 design_matrix_raw =  theObject@design_matrix )

  return_list$de_proteins_long <- de_proteins_long


  ## Plot static volcano plot
  static_volcano_plot_data <- de_proteins_long |>
    mutate( lqm = -log10(q.mod ) ) |>
    dplyr::mutate(label = case_when( q.mod < de_q_val_thresh ~ "Significant",
                                     TRUE ~ "Not sig.")) |>
    dplyr::mutate(colour = case_when( q.mod < de_q_val_thresh ~ "purple",
                                      TRUE ~ "black")) |>
    dplyr::mutate(colour = factor(colour, levels = c("black", "purple")))

  list_of_volcano_plots <- static_volcano_plot_data %>%
    group_by( comparison) %>%
    nest() %>%
    ungroup() %>%
    mutate( title = paste( comparison)) %>%
    mutate( plot = purrr:::map2( data, title, \(x,y) { plotOneVolcanoNoVerticalLines(x, y
                                                                                     , log_q_value_column = lqm
                                                                                     , log_fc_column = log2FC) } ) )

  return_list$list_of_volcano_plots <- list_of_volcano_plots


  ## Return the number of significant molecules
  num_sig_de_molecules <- significant_rows %>%
    dplyr::mutate(status = case_when(q.mod  >= de_q_val_thresh ~ "Not significant",
                                       log2FC >= 0 & q.mod < de_q_val_thresh ~ "Significant and Up",
                                     log2FC < 0 &  q.mod < de_q_val_thresh ~ "Significant and Down",
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

#' @export
# de_analysis_results_list$contrasts_results$fit.eb

writeInteractiveVolcanoPlotProteomics <- function( de_proteins_long, uniprot_tbl, fit.eb, args_row_id = "uniprot_acc", publication_graphs_dir, de_q_val_thresh = 0.05) {


    volcano_plot_tab <- de_proteins_long  |>
      left_join(uniprot_tbl, by = join_by( !!sym(args_row_id) == Entry ) ) |>
      dplyr::rename( UNIPROT_GENENAME = "Gene Names" ) |>
      dplyr::rename(PROTEIN_NAMES = "Protein names") |>
      mutate( UNIPROT_GENENAME = purrr::map_chr( UNIPROT_GENENAME, \(x){str_split(x, " ")[[1]][1]})) |>
      mutate( lqm = -log10(q.mod))  |>
      dplyr::mutate(label = case_when(abs(log2FC) >= 1 & q.mod >= de_q_val_thresh ~ "Not sig., logFC >= 1",
                                      abs(log2FC) >= 1 & q.mod < de_q_val_thresh ~ "Sig., logFC >= 1",
                                      abs(log2FC) < 1 & q.mod < de_q_val_thresh ~ "Sig., logFC < 1",
                                      TRUE ~ "Not sig.")) |>
      dplyr::mutate(colour = case_when(abs(log2FC) >= 1 & q.mod >= de_q_val_thresh ~ "orange",
                                       abs(log2FC) >= 1 & q.mod < de_q_val_thresh ~ "purple",
                                       abs(log2FC) < 1 & q.mod < de_q_val_thresh ~ "blue",
                                       TRUE ~ "black")) |>
      dplyr::mutate(gene_name = str_split(UNIPROT_GENENAME, " |:" ) |> purrr::map_chr(1)  ) |>
      dplyr::mutate(best_uniprot_acc = str_split(!!sym(args_row_id), ":" ) |> purrr::map_chr(1)  ) |>
      dplyr::mutate(analysis_type = comparison)  |>
      dplyr::select( best_uniprot_acc, lqm, q.mod, p.mod, log2FC, comparison, label, colour,  gene_name, `PROTEIN_NAMES`)   |>
      dplyr::mutate( my_alpha = case_when ( gene_name !=  "" ~ 1
                                            , TRUE ~ 0.5))

    print(head(volcano_plot_tab))

    ncol(fit.eb$coefficients)
    colnames(fit.eb$coefficients)

    output_dir <- file.path( publication_graphs_dir
                             ,  "Interactive_Volcano_Plots")

    dir.create(output_dir, recursive = TRUE)

    purrr::walk( seq_len( ncol(fit.eb$coefficients))
                 , \(coef) { # print(coef)
                   ProteomeRiver::getGlimmaVolcanoProteomics( fit.eb
                                                              , coef = coef
                                                              , volcano_plot_tab  = volcano_plot_tab
                                                              , uniprot_column = best_uniprot_acc
                                                              , gene_name_column = gene_name
                                                              , display_columns = c( "best_uniprot_acc",  "PROTEIN_NAMES"   )
                                                              , output_dir = output_dir ) } )

}


#' @export
outputDeAnalysisResults <- function(de_analysis_results_list, uniprot_tbl
                                    , de_output_dir
                                    , publication_graphs_dir
                                    , file_prefix
                                    , plots_format
                                    , args_row_id = "uniprot_acc"
                                    , de_q_val_thresh = 0.05) {

  ## PCA plot
  plot_pca_plot <- de_analysis_results_list$pca_plot

  for( format_ext in plots_format) {
    file_name <- file.path( publication_graphs_dir, "PCA", paste0("PCA_plot.",format_ext))
    ggsave(filename = file_name, plot = plot_pca_plot, limitsize = FALSE)
  }

  ## RLE plot
  plot_rle_plot <- de_analysis_results_list$rle_plot

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

  ## Count the number of up or down significnat differentially expressed proteins.
  num_sig_de_molecules <- de_analysis_results_list$num_sig_de_molecules

  for( format_ext in plots_format) {
    file_name<-file.path(de_output_dir, paste0("num_sda_entities_barplot.",format_ext))
    ggsave(filename = file_name,
           plot = num_sig_de_molecules$plot,
           height = 10,
           width = 7)
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


  ## Print Volcano plot

  volcano_plot <- de_analysis_results_list$volcano_plot
  for( format_ext in plots_format) {
    file_name <- file.path(de_output_dir,paste0("volplot_gg_all.",format_ext))
    ggsave(filename = file_name, plot = volplot_plot, width = 7.29, height = 6)
  }

  ## Count the number of up or down significnat differentially expressed proteins.
  num_sig_de_molecules <- de_analysis_results_list$num_sig_de_molecules
  for( format_ext in plots_format) {
    file_name<-file.path(de_output_dir,paste0("num_sda_entities_barplot.",format_ext))
    ggsave(filename = file_name,
           plot = num_sig_de_molecules$plot,
           height = 10,
           width = 7)

  }



  ## Count the number of up or down significnat differentially expressed proteins.
  num_sig_de_molecules_first_go <- de_analysis_results_list$num_sig_de_molecules_first_go
  vroom::vroom_write(num_sig_de_molecules_first_go$table,
                     file.path(de_output_dir,
                               "num_significant_differentially_abundant_all.tab"))

  writexl::write_xlsx(num_sig_de_molecules_first_go$table,
                      file.path(de_output_dir,
                                "num_significant_differentially_abundant_all.xlsx"))



  ## Print p-values distribution figure
  pvalhist <- de_analysis_results_list$pvalhist
  for( format_ext in plots_format) {
    file_name<-file.path(de_output_dir,paste0("p_values_distn.",format_ext))
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



  ## Create long format output file
  de_proteins_long <- de_analysis_results_list$de_proteins_long
  vroom::vroom_write( de_proteins_long,
                      file.path( de_output_dir,
                                 paste0(file_prefix, "_long.tsv")))

  writexl::write_xlsx( de_proteins_long,
                       file.path( de_output_dir,
                                  paste0(file_prefix, "_long.xlsx")))


  ## Static volcano plots
  dir.create(file.path( publication_graphs_dir, "Volcano_Plots"), recursive = TRUE)

  list_of_volcano_plots <- de_analysis_results_list$list_of_volcano_plots

  purrr::walk2( list_of_volcano_plots %>% pull(title),
                list_of_volcano_plots %>% pull(plot),
                ~{file_name_part <- file.path( publication_graphs_dir, "Volcano_Plots", paste0(.x, "."))
                # gg_save_logging ( .y, file_name_part, plots_format)
                for( format_ext in plots_format) {
                  file_name <- paste0(file_name_part, format_ext)
                    ggsave(plot=.y
                           , filename = file_name
                           , width=7
                           , height=7 )
                }
                } )

  ggsave(
    filename = file.path(publication_graphs_dir, "Volcano_Plots", "list_of_volcano_plots.pdf" ),
    plot = gridExtra::marrangeGrob( (list_of_volcano_plots  %>% pull(plot)), nrow=1, ncol=1),
    width = 7, height = 7
  )

  ## Number of significant molecules
  createDirIfNotExists(file.path(publication_graphs_dir, "NumSigDeMolecules"))
  vroom::vroom_write( de_analysis_results_list$num_sig_de_molecules,
                      file.path(publication_graphs_dir, "NumSigDeMolecules", "num_sig_de_molecules.tab" ) )


  if( !is.null(de_analysis_results_list$num_sig_de_genes_barplot_only_significant)) {
    num_sig_de_genes_barplot_only_significant <- de_analysis_results_list$num_sig_de_genes_barplot_only_significant
    num_of_comparison_only_significant <- de_analysis_results_list$num_of_comparison_only_significant


    ggsave(filename = file.path(publication_graphs_dir, "NumSigDeMolecules", "num_sig_de_genes_barplot.png" ),
           plot = num_sig_de_genes_barplot_only_significant,
           height = 6,
           width = (num_of_comparison_only_significant + 2) *7/6 )

    ggsave(filename = file.path(publication_graphs_dir, "NumSigDeMolecules", "num_sig_de_genes_barplot.pdf" ),
           plot = num_sig_de_genes_barplot_only_significant,
           height = 6,
           width = (num_of_comparison_only_significant + 2) *7/6 )

  }


  if(!is.null( de_analysis_results_list$num_sig_de_genes_barplot_with_not_significant )) {
    num_sig_de_genes_barplot_with_not_significant <- de_analysis_results_list$num_sig_de_genes_barplot_with_not_significant
    num_of_comparison_with_not_significant <- de_analysis_results_list$num_of_comparison_with_not_significant

    ggsave(filename = file.path(publication_graphs_dir, "NumSigDeMolecules", "num_sig_de_molecules_with_not_significant.png" ),
           plot = num_sig_de_genes_barplot_with_not_significant,
           height = 6,
           width = (num_of_comparison_with_not_significant + 2) *7/6 )

    ggsave(filename = file.path(publication_graphs_dir, "NumSigDeMolecules", "num_sig_de_molecules_with_not_significant.pdf" ),
           plot = num_sig_de_genes_barplot_with_not_significant,
           height = 6,
           width = (num_of_comparison_with_not_significant + 2) *7/6 )


  }

  ## Write interactive volcano plot
  writeInteractiveVolcanoPlotProteomics( de_proteins_long
                                         , uniprot_tbl
                                         , contrasts_results$fit.eb
                                         , args_row_id = args_row_id
                                         , publication_graphs_dir
                                         , de_q_val_thresh = de_q_val_thresh)


}








