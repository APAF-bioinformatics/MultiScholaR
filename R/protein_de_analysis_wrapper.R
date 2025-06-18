##----------------------------------------------------------------------------------------------------------------------------------------------------------------------
#' Differential Expression Analysis for Proteomics
#' 
#' This file contains modular functions for proteomics differential expression analysis,
#' breaking up the monolithic deAnalysisWrapperFunction into focused components.
#' 
#' Functions:
#' - differentialExpressionAnalysis: Main S4 generic for DE analysis
#' - differentialExpressionAnalysisHelper: Core limma-based analysis
#' - generateVolcanoPlotGlimma: Interactive volcano plot generation
#' - generateDEHeatmap: Heatmap visualization with clustering
#' 
##----------------------------------------------------------------------------------------------------------------------------------------------------------------------

#' @export
setGeneric( name ="differentialExpressionAnalysis"
            , def=function(theObject
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
                           , raw_pvalue_column = "raw_pvalue") {
              standardGeneric("differentialExpressionAnalysis")
            })

#' @export
setMethod( f ="differentialExpressionAnalysis"
           , signature = "ProteinQuantitativeData"
           , definition=function( theObject
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
                                  , raw_pvalue_column = "raw_pvalue" ) {

  message("--- Entering differentialExpressionAnalysis ---")
  message(sprintf("   differentialExpressionAnalysis: theObject class = %s", class(theObject)))

  # Run DE analysis using the helper function
  results_list <- differentialExpressionAnalysisHelper(  theObject
                                                        , contrasts_tbl = contrasts_tbl
                                                        , formula_string = formula_string
                                                        , group_id = group_id
                                                        , de_q_val_thresh = de_q_val_thresh
                                                        , treat_lfc_cutoff = treat_lfc_cutoff
                                                        , eBayes_trend = eBayes_trend
                                                        , eBayes_robust = eBayes_robust
                                                        , args_group_pattern = args_group_pattern
                                                        , args_row_id = args_row_id
                                                        , qvalue_column = qvalue_column
                                                        , raw_pvalue_column = raw_pvalue_column
                                                        )

  message("--- Exiting differentialExpressionAnalysis ---")
  return(results_list)

})

##----------------------------------------------------------------------------------------------------------------------------------------------------------------------
#' Core Differential Expression Analysis Helper
#' 
#' @export
setGeneric( name ="differentialExpressionAnalysisHelper"
            , def=function(theObject
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
                           , raw_pvalue_column = "raw_pvalue") {
              standardGeneric("differentialExpressionAnalysisHelper")
            })

#' @export
setMethod( f ="differentialExpressionAnalysisHelper"
           , signature = "ProteinQuantitativeData"
           , definition=function( theObject
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
                                  , raw_pvalue_column = "raw_pvalue" ) {

  message("--- Entering differentialExpressionAnalysisHelper ---")

  # Extract parameters from S4 object with fallbacks
  contrasts_tbl <- checkParamsObjectFunctionSimplify( theObject, "contrasts_tbl", contrasts_tbl)
  formula_string <- checkParamsObjectFunctionSimplify( theObject, "formula_string", "~ 0 + group")
  group_id <- checkParamsObjectFunctionSimplify( theObject, "group_id", "group")
  de_q_val_thresh <- checkParamsObjectFunctionSimplify( theObject, "de_q_val_thresh", 0.05)
  treat_lfc_cutoff <- checkParamsObjectFunctionSimplify( theObject, "treat_lfc_cutoff", 0)
  eBayes_trend <- checkParamsObjectFunctionSimplify( theObject, "eBayes_trend", TRUE)
  eBayes_robust <- checkParamsObjectFunctionSimplify( theObject, "eBayes_robust", TRUE)
  args_group_pattern <- checkParamsObjectFunctionSimplify( theObject, "args_group_pattern", "(\\d+)")
  args_row_id <- checkParamsObjectFunctionSimplify( theObject, "args_row_id", "uniprot_acc")

  message(sprintf("   differentialExpressionAnalysisHelper: formula_string = %s", formula_string))
  message(sprintf("   differentialExpressionAnalysisHelper: group_id = %s", group_id))
  message(sprintf("   differentialExpressionAnalysisHelper: de_q_val_thresh = %f", de_q_val_thresh))

  # Handle group names that start with numbers (same pattern as original wrapper)
  design_matrix <- theObject@design_matrix
  group_col <- design_matrix[[group_id]]
  
  # Check if any group names start with numbers and create mapping
  starts_with_number <- grepl("^[0-9]", group_col)
  if(any(starts_with_number)) {
    message("   differentialExpressionAnalysisHelper Step: Fixing group names that start with numbers...")
    original_groups <- unique(group_col)
    safe_groups <- purrr::map_chr(original_groups, \(x) {
      if(grepl("^[0-9]", x)) paste0("grp_", x) else x
    })
    group_mapping <- setNames(original_groups, safe_groups)
    
    # Update design matrix with safe names
    design_matrix[[group_id]] <- purrr::map_chr(group_col, \(x) {
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

  # Update object with validated parameters
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

  message("   differentialExpressionAnalysisHelper Step: Generating QC plots...")

  # Generate QC plots (RLE, PCA, density)
  rle_plot <- plotRle(theObject = theObject, theObject@group_id) +
    theme(axis.text.x = element_text(size = 13)) +
    theme(axis.text.y = element_text(size = 13)) +
    theme(axis.title.x = element_text(size = 12)) +
    theme(axis.title.y = element_text(size = 12)) +
    theme(plot.title = element_text(size = 12)) +
    theme(legend.text = element_text(size = 12)) +
    theme(legend.title = element_text(size = 12)) +
    xlab("Samples")

  return_list$rle_plot <- rle_plot

  pca_plot <- plotPca( theObject
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

  pca_plot_with_labels <- plotPca( theObject
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

  # Count the number of values
  return_list$plot_num_of_values <- plotNumOfValuesNoLog(theObject@protein_quant_table)

  message("   differentialExpressionAnalysisHelper Step: Running limma contrasts analysis...")

  # Prepare design matrix for limma
  rownames( theObject@design_matrix ) <- theObject@design_matrix |> dplyr::pull( one_of(theObject@sample_id ))

  # Run the core limma analysis using existing function
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

  message("   differentialExpressionAnalysisHelper Step: Preparing data for visualization...")

  # Prepare data for volcano plots
  significant_rows <- getSignificantData(list_of_de_tables = list(contrasts_results_table),
                                         list_of_descriptions = list("RUV applied"),
                                         row_id = !!sym(args_row_id),
                                         p_value_column = !!sym(raw_pvalue_column),
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

  # Generate volcano plots
  volplot_plot <- plotVolcano(significant_rows,
                              log_q_value_column = lqm,
                              log_fc_column = log2FC,
                              q_val_thresh = de_q_val_thresh,
                              formula_string = "analysis_type ~ comparison")

  return_list$volplot_plot <- volplot_plot

  # Count significant molecules
  num_sig_de_molecules <- printCountDeGenesTable(list_of_de_tables = list(contrasts_results_table),
                                                 list_of_descriptions = list("RUV applied"),
                                                 formula_string = "analysis_type ~ comparison")

  return_list$num_sig_de_molecules_first_go <- num_sig_de_molecules

  # P-values distribution plot
  pvalhist <- printPValuesDistribution(significant_rows,
                                       p_value_column = !!sym(raw_pvalue_column),
                                       formula_string = "analysis_type ~ comparison")

  return_list$pvalhist <- pvalhist

  message("   differentialExpressionAnalysisHelper Step: Creating results tables...")

  # Create wide format output
  counts_table_to_use <- theObject@protein_quant_table

  norm_counts <- counts_table_to_use |>
    as.data.frame() |>
    column_to_rownames(args_row_id) |>
    set_colnames(paste0(colnames(counts_table_to_use[-1]), ".log2norm")) |>
    rownames_to_column(args_row_id)

  return_list$norm_counts <- norm_counts

  de_proteins_wide <- significant_rows |>
    dplyr::filter(analysis_type == "RUV applied") |>
    dplyr::select(-lqm, -colour, -analysis_type) |>
    pivot_wider(id_cols = c(!!sym(args_row_id)),
                names_from = c(comparison),
                names_sep = ":",
                values_from = c(log2FC, !!sym(qvalue_column), !!sym(raw_pvalue_column))) |>
    left_join(counts_table_to_use, by = join_by( !!sym(args_row_id) == !!sym(theObject@protein_id_column))) |>
    left_join(theObject@protein_id_table, by = join_by( !!sym(args_row_id) == !!sym(theObject@protein_id_column))) |>
    dplyr::arrange(across(matches(paste0("!!sym(", qvalue_column, ")")))) |>
    distinct()

  return_list$de_proteins_wide <- de_proteins_wide

  # Create long format output
  de_proteins_long <- createDeResultsLongFormat( lfc_qval_tbl = significant_rows |>
                                                   dplyr::filter(analysis_type == "RUV applied"),
                                                 norm_counts_input_tbl = as.matrix(column_to_rownames(theObject@protein_quant_table, theObject@protein_id_column)),
                                                 raw_counts_input_tbl = as.matrix(column_to_rownames(theObject@protein_quant_table, theObject@protein_id_column)),
                                                 row_id = args_row_id,
                                                 sample_id = theObject@sample_id,
                                                 group_id = group_id,
                                                 group_pattern = args_group_pattern,
                                                 design_matrix_norm = theObject@design_matrix,
                                                 design_matrix_raw = theObject@design_matrix,
                                                 protein_id_table = theObject@protein_id_table)

  return_list$de_proteins_long <- de_proteins_long

  # Static volcano plots with gene names
  static_volcano_plot_data <- de_proteins_long |>
    mutate( lqm = -log10(!!sym(qvalue_column))) |>
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

  # Additional summary statistics
  num_sig_de_molecules <- significant_rows %>%
    dplyr::mutate(status = case_when(!!sym(qvalue_column) >= de_q_val_thresh ~ "Not significant",
                                     log2FC >= 0 & !!sym(qvalue_column) < de_q_val_thresh ~ "Significant and Up",
                                     log2FC < 0 & !!sym(qvalue_column) < de_q_val_thresh ~ "Significant and Down",
                                     TRUE ~ "Not significant")) %>%
    group_by( comparison, status) %>%
    summarise(counts = n()) %>%
    ungroup()

  return_list$num_sig_de_molecules <- num_sig_de_molecules

  # Create barplots if significant results exist
  if (num_sig_de_molecules %>%
      dplyr::filter(status != "Not significant") |>
      nrow() > 0 ) {

    num_sig_de_genes_barplot_only_significant <- num_sig_de_molecules %>%
      dplyr::filter(status != "Not significant") %>%
      ggplot(aes(x = status, y = counts)) +
      geom_bar(stat = "identity") +
      geom_text(stat = 'identity', aes(label = counts), vjust = -0.5) +
      theme(axis.text.x = element_text(angle = 90)) +
      facet_wrap(~ comparison)

    return_list$num_sig_de_genes_barplot_only_significant <- num_sig_de_genes_barplot_only_significant

    num_sig_de_genes_barplot_with_not_significant <- num_sig_de_molecules %>%
      ggplot(aes(x = status, y = counts)) +
      geom_bar(stat = "identity") +
      geom_text(stat = 'identity', aes(label = counts), vjust = -0.5) +
      theme(axis.text.x = element_text(angle = 90)) +
      facet_wrap(~ comparison)

    return_list$num_sig_de_genes_barplot_with_not_significant <- num_sig_de_genes_barplot_with_not_significant
  }

  message("--- Exiting differentialExpressionAnalysisHelper ---")
  return(return_list)

})

##----------------------------------------------------------------------------------------------------------------------------------------------------------------------
#' Generate Interactive Volcano Plot using Glimma
#' 
#' @export
generateVolcanoPlotGlimma <- function( de_results_list
                                       , selected_contrast = NULL
                                       , uniprot_tbl = NULL
                                       , args_row_id = "uniprot_acc"
                                       , fdr_column = "fdr_qvalue"
                                       , raw_p_value_column = "raw_pvalue"
                                       , log2fc_column = "log2FC"
                                       , de_q_val_thresh = 0.05
                                       , uniprot_id_column = "Entry"
                                       , gene_names_column = "Gene Names"
                                       , display_columns = c( "best_uniprot_acc" )) {

  message("--- Entering generateVolcanoPlotGlimma ---")
  message(sprintf("   generateVolcanoPlotGlimma Arg: selected_contrast = %s", selected_contrast))

  if (is.null(de_results_list) || is.null(de_results_list$de_proteins_long)) {
    message("   generateVolcanoPlotGlimma: No DE results available")
    return(NULL)
  }

  if (is.null(selected_contrast)) {
    message("   generateVolcanoPlotGlimma: No contrast selected")
    return(NULL)
  }

  # Get the contrast-specific results
  de_proteins_long <- de_results_list$de_proteins_long
  contrasts_results <- de_results_list$contrasts_results

  # Filter for selected contrast
  contrast_data <- de_proteins_long |>
    dplyr::filter(comparison == selected_contrast)

  if (nrow(contrast_data) == 0) {
    message(sprintf("   generateVolcanoPlotGlimma: No data found for contrast %s", selected_contrast))
    return(NULL)
  }

  message("   generateVolcanoPlotGlimma Step: Preparing volcano plot data...")

  # Prepare volcano plot table
  volcano_plot_tab <- contrast_data |>
    dplyr::mutate(best_uniprot_acc = str_split(!!sym(args_row_id), ":") |> purrr::map_chr(1))

  # Add UniProt annotations if available
  if (!is.null(uniprot_tbl)) {
    volcano_plot_tab <- volcano_plot_tab |>
      left_join(uniprot_tbl, by = join_by( best_uniprot_acc == !!sym( uniprot_id_column) )) |>
      dplyr::rename( UNIPROT_GENENAME = !!sym(gene_names_column) ) |>
      mutate( UNIPROT_GENENAME = purrr::map_chr( UNIPROT_GENENAME, \(x){str_split(x, " |:")[[1]][1]})) |>
      dplyr::mutate(gene_name = UNIPROT_GENENAME )
  } else {
    volcano_plot_tab <- volcano_plot_tab |>
      dplyr::mutate(gene_name = best_uniprot_acc)
  }

  # Add visualization columns
  volcano_plot_tab <- volcano_plot_tab |>
    mutate( lqm = -log10(!!sym(fdr_column))) |>
    dplyr::mutate(label = case_when(abs(!!sym(log2fc_column)) >= 1 & !!sym(fdr_column) >= de_q_val_thresh ~ "Not sig., logFC >= 1",
                                    abs(!!sym(log2fc_column)) >= 1 & !!sym(fdr_column) < de_q_val_thresh ~ "Sig., logFC >= 1",
                                    abs(!!sym(log2fc_column)) < 1 & !!sym(fdr_column) < de_q_val_thresh ~ "Sig., logFC < 1",
                                    TRUE ~ "Not sig.")) |>
    dplyr::mutate(colour = case_when(abs(!!sym(log2fc_column)) >= 1 & !!sym(fdr_column) >= de_q_val_thresh ~ "orange",
                                     abs(!!sym(log2fc_column)) >= 1 & !!sym(fdr_column) < de_q_val_thresh ~ "purple",
                                     abs(!!sym(log2fc_column)) < 1 & !!sym(fdr_column) < de_q_val_thresh ~ "blue",
                                     TRUE ~ "black")) |>
    dplyr::mutate(analysis_type = comparison) |>
    dplyr::select( best_uniprot_acc, lqm, !!sym(fdr_column), !!sym(raw_p_value_column), !!sym(log2fc_column), comparison, label, colour, gene_name
                   , any_of( display_columns )) |>
    dplyr::mutate( my_alpha = case_when ( gene_name != "" ~ 1
                                          , TRUE ~ 0.5))

  message("   generateVolcanoPlotGlimma Step: Generating interactive plot widget...")

  # Find the coefficient index for this contrast
  coef_names <- colnames(contrasts_results$fit.eb$coefficients)
  coef_index <- which(coef_names == selected_contrast)

  if (length(coef_index) == 0) {
    message(sprintf("   generateVolcanoPlotGlimma: Contrast %s not found in coefficients", selected_contrast))
    return(NULL)
  }

  # Prepare counts matrix and groups
  counts_mat <- de_results_list$theObject@protein_quant_table |>
    column_to_rownames(de_results_list$theObject@protein_id_column) |>
    as.matrix()

  this_design_matrix <- de_results_list$theObject@design_matrix
  rownames( this_design_matrix ) <- this_design_matrix[, de_results_list$theObject@sample_id]
  this_groups <- this_design_matrix[colnames( counts_mat), "group"]

  # Generate the Glimma widget
  glimma_widget <- getGlimmaVolcanoProteomicsWidget( contrasts_results$fit.eb
                                                    , coef = coef_index
                                                    , volcano_plot_tab = volcano_plot_tab
                                                    , uniprot_column = best_uniprot_acc
                                                    , gene_name_column = gene_name
                                                    , display_columns = display_columns
                                                    , counts_tbl = counts_mat
                                                    , groups = this_groups )

  message("--- Exiting generateVolcanoPlotGlimma ---")
  return(glimma_widget)
}

##----------------------------------------------------------------------------------------------------------------------------------------------------------------------
#' Generate DE Results Heatmap with Advanced Clustering
#' 
#' @export
generateDEHeatmap <- function( de_results_list
                              , selected_contrast = NULL
                              , top_n_genes = 50
                              , clustering_method = "ward.D2"
                              , distance_method = "euclidean"
                              , cluster_rows = TRUE
                              , cluster_cols = TRUE
                              , scale_data = "row"
                              , color_scheme = "RdBu"
                              , show_gene_names = FALSE
                              , de_q_val_thresh = 0.05
                              , qvalue_column = "fdr_qvalue"
                              , log2fc_column = "log2FC") {

  message("--- Entering generateDEHeatmap ---")
  message(sprintf("   generateDEHeatmap Arg: selected_contrast = %s", selected_contrast))
  message(sprintf("   generateDEHeatmap Arg: top_n_genes = %d", top_n_genes))

  if (is.null(de_results_list) || is.null(de_results_list$de_proteins_long)) {
    message("   generateDEHeatmap: No DE results available")
    return(NULL)
  }

  if (is.null(selected_contrast)) {
    message("   generateDEHeatmap: No contrast selected")
    return(NULL)
  }

  message("   generateDEHeatmap Step: Filtering data for selected contrast...")

  # Get the contrast-specific results
  de_proteins_long <- de_results_list$de_proteins_long
  
  # Filter for selected contrast and significant results
  contrast_data <- de_proteins_long |>
    dplyr::filter(comparison == selected_contrast) |>
    dplyr::filter(!!sym(qvalue_column) < de_q_val_thresh)

  if (nrow(contrast_data) == 0) {
    message(sprintf("   generateDEHeatmap: No significant genes found for contrast %s", selected_contrast))
    return(NULL)
  }

  message(sprintf("   generateDEHeatmap Step: Found %d significant genes", nrow(contrast_data)))

  # Order by absolute fold change and take top N
  top_genes_data <- contrast_data |>
    dplyr::arrange(desc(abs(!!sym(log2fc_column)))) |>
    head(top_n_genes)

  message(sprintf("   generateDEHeatmap Step: Selected top %d genes for heatmap", nrow(top_genes_data)))

  # Extract expression matrix for these genes
  protein_ids <- top_genes_data |> dplyr::pull(!!sym(names(top_genes_data)[1])) # First column should be protein ID

  # Get expression data from the S4 object
  expr_matrix <- de_results_list$theObject@protein_quant_table |>
    dplyr::filter(!!sym(de_results_list$theObject@protein_id_column) %in% protein_ids) |>
    column_to_rownames(de_results_list$theObject@protein_id_column) |>
    as.matrix()

  if (nrow(expr_matrix) == 0) {
    message("   generateDEHeatmap: No expression data found for selected genes")
    return(NULL)
  }

  message(sprintf("   generateDEHeatmap Step: Expression matrix dims = %d rows, %d cols", nrow(expr_matrix), ncol(expr_matrix)))

  # Apply scaling
  if (scale_data == "row") {
    message("   generateDEHeatmap Step: Applying row scaling...")
    expr_matrix <- t(scale(t(expr_matrix)))
  } else if (scale_data == "column") {
    message("   generateDEHeatmap Step: Applying column scaling...")
    expr_matrix <- scale(expr_matrix)
  } else if (scale_data == "both") {
    message("   generateDEHeatmap Step: Applying both row and column scaling...")
    expr_matrix <- scale(t(scale(t(expr_matrix))))
  } else {
    message("   generateDEHeatmap Step: No scaling applied")
  }

  message("   generateDEHeatmap Step: Setting up clustering parameters...")

  # Set up clustering
  row_clust <- NULL
  col_clust <- NULL

  if (cluster_rows || cluster_cols) {
    # Calculate distance matrices
    if (distance_method %in% c("pearson", "spearman")) {
      if (cluster_rows) {
        if (distance_method == "pearson") {
          row_dist <- as.dist(1 - cor(t(expr_matrix), method = "pearson", use = "complete.obs"))
        } else {
          row_dist <- as.dist(1 - cor(t(expr_matrix), method = "spearman", use = "complete.obs"))
        }
        row_clust <- hclust(row_dist, method = clustering_method)
      }
      
      if (cluster_cols) {
        if (distance_method == "pearson") {
          col_dist <- as.dist(1 - cor(expr_matrix, method = "pearson", use = "complete.obs"))
        } else {
          col_dist <- as.dist(1 - cor(expr_matrix, method = "spearman", use = "complete.obs"))
        }
        col_clust <- hclust(col_dist, method = clustering_method)
      }
    } else {
      # Standard distance metrics
      if (cluster_rows) {
        row_dist <- dist(expr_matrix, method = distance_method)
        row_clust <- hclust(row_dist, method = clustering_method)
      }
      
      if (cluster_cols) {
        col_dist <- dist(t(expr_matrix), method = distance_method)
        col_clust <- hclust(col_dist, method = clustering_method)
      }
    }
  }

  message("   generateDEHeatmap Step: Setting up color scheme...")

  # Choose color palette
  colors <- switch(color_scheme,
    "RdBu" = colorRampPalette(c("red", "white", "blue"))(100),
    "RdYlBu" = colorRampPalette(c("red", "yellow", "blue"))(100),
    "coolwarm" = colorRampPalette(c("blue", "white", "red"))(100),
    "viridis" = viridis::viridis(100),
    "plasma" = viridis::plasma(100),
    "inferno" = viridis::inferno(100),
    colorRampPalette(c("red", "white", "blue"))(100)
  )

  message("   generateDEHeatmap Step: Generating heatmap...")

  # Generate the heatmap
  if (requireNamespace("ComplexHeatmap", quietly = TRUE)) {
    message("   generateDEHeatmap: Using ComplexHeatmap for advanced visualization")
    
    heatmap_plot <- ComplexHeatmap::Heatmap(
      expr_matrix,
      name = "Expression",
      col = colors,
      cluster_rows = if (cluster_rows) row_clust else FALSE,
      cluster_columns = if (cluster_cols) col_clust else FALSE,
      show_row_names = show_gene_names,
      show_column_names = TRUE,
      show_row_dend = cluster_rows,
      show_column_dend = cluster_cols,
      column_title = paste("Heatmap - Top", nrow(expr_matrix), "genes\nContrast:", selected_contrast),
      heatmap_legend_param = list(title = "Scaled\nExpression")
    )
  } else {
    message("   generateDEHeatmap: Using base R heatmap")
    
    # Create a plot function that returns the heatmap
    heatmap_plot <- function() {
      heatmap(
        expr_matrix,
        main = paste("Heatmap - Top", nrow(expr_matrix), "genes\nContrast:", selected_contrast),
        col = colors,
        Rowv = if (cluster_rows) as.dendrogram(row_clust) else NA,
        Colv = if (cluster_cols) as.dendrogram(col_clust) else NA,
        labRow = if (show_gene_names) rownames(expr_matrix) else rep("", nrow(expr_matrix)),
        cexRow = if (show_gene_names) 0.8 else 0.1
      )
    }
  }

  message("--- Exiting generateDEHeatmap ---")
  return(heatmap_plot)
} 