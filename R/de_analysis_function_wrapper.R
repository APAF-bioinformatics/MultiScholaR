


#' @export
deAnalysisWrapperFunction <- function( theObject
                                       , contrasts_tbl
                                       , formula_string = " ~ 0 + group"
                                       , de_q_val_thresh = 0.05
                                       , treat_lfc_cutoff = 0
                                       , eBayes_trend = TRUE
                                       , eBayes_robust = TRUE
                                       , args_group_pattern = "(\\d+)"
                                       , args_row_id = "uniprot_acc" ) {
  
  return_list <- list()
  
  ## Count the number of values 
  return_list$plot_num_of_values <- plotNumOfValuesNoLog(theObject@protein_data)
  
  ## Compare the different experimental groups and obtain lists of differentially expressed proteins.")
  
  rownames( theObject@design_matrix ) <- theObject@design_matrix |> pull( one_of("replicates"))
  
  # requires statmod library
  contrasts_results <- runTestsContrasts(theObject@protein_data,
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
                                         expression_column = expression,
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
  num_sig_de_molecules <- printCountDeGenesTable(list_of_de_tables = list(contrasts_results_table),
                                                 list_of_descriptions = list( "RUV applied"),
                                                 formula_string = "analysis_type ~ comparison")
  
  return_list$num_sig_de_molecules <- num_sig_de_molecules
  
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
    left_join(norm_counts, by = args_row_id) |>
    dplyr::arrange(across(matches("q.mod"))) |>
    distinct()
  
  
  return_list$de_proteins_wide <- de_proteins_wide
  
  
  ## Create long format output file
  
  de_proteins_long <- createDeResultsLongFormat( lfc_qval_tbl = significant_rows |>
                                                   dplyr::filter(analysis_type == "RUV applied") ,
                                                 norm_counts_input_tbl = counts_table_to_use,
                                                 raw_counts_input_tbl = counts_table_to_use,
                                                 row_id = args_row_id,
                                                 sample_id = theObject@sample_id,
                                                 group_id = theObject@group_id,
                                                 group_pattern = args_group_pattern,
                                                 design_matrix_norm = theObject@design_matrix,
                                                 design_matrix_raw =  theObject@design_matrix )
  
  return_list$de_proteins_long <- de_proteins_long
  
  
  return_list
  
}

#' @export

outputDeAnalysisResults <- function(de_analysis_results_list, de_output_dir, plots_format ) {
  
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
    file_name <- file.path(de_output_dir,paste0("volplot_gg_all.",format_ext))
    ggsave(filename = file_name, plot = volplot_plot, width = 7.29, height = 6)
  }
  
  ## Count the number of up or down significnat differentially expressed proteins.
  num_sig_de_molecules <- de_analysis_results_list$num_sig_de_molecules
  
  
  for( format_ext in args$plots_format) {
    file_name<-file.path(args$output_dir,paste0("num_sda_entities_barplot.",format_ext))
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
  num_sig_de_molecules <- de_analysis_results_list$num_sig_de_molecules
  vroom::vroom_write(num_sig_de_molecules$table,
                     file.path(de_output_dir,
                               "num_significant_differentially_abundant_all.tab"))
  
  writexl::write_xlsx(num_sig_de_molecules$table,
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
                                 paste0(args$file_prefix, "_wide.tsv")))
  
  writexl::write_xlsx( de_proteins_wide,
                       file.path( de_output_dir,
                                  paste0(args$file_prefix, "_wide.xlsx")))
  
  
  
  ## Create long format output file
  de_proteins_long <- de_analysis_results_list$de_proteins_long
  vroom::vroom_write( de_proteins_long,
                      file.path( de_output_dir,
                                 paste0(args$file_prefix, "_long.tsv")))
  
  writexl::write_xlsx( de_proteins_long,
                       file.path( de_output_dir,
                                  paste0(args$file_prefix, "_long.xlsx")))
  
  
}


