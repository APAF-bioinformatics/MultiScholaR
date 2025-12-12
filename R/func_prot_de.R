# ============================================================================
# func_prot_de.R
# ============================================================================
# Purpose: Protein differential expression analysis functions
# 
# This file contains functions for protein differential expression analysis,
# including limma-based analysis, result formatting, and visualization.
# Functions in this file are used by mod_prot_de.R and related DE modules.
#
# Functions to extract here:
# - differentialExpressionAnalysis(): S4 method for DE analysis
# - differentialExpressionAnalysisHelper(): Helper for DE analysis
# - deAnalysisWrapperFunction(): Wrapper function for DE analysis
# - outputDeResultsAllContrasts(): S4 method for outputting DE results
# - generateVolcanoPlotGlimma(): Generate interactive volcano plots
# - generateDEHeatmap(): Generate DE heatmaps
# - createDeResultsLongFormat(): Create long format DE results
# - getDeResultsLongFormat(): S4 method to get long format results
# - getDeResultsWideFormat(): S4 method to get wide format results
# - Additional DE helper functions
#
# Dependencies:
# - limma, edgeR
# - func_general_plotting.R (for visualization)
# - func_general_helpers.R (for utility functions)
# ============================================================================

# TODO: Extract the following functions from their current locations:

# Function 1: differentialExpressionAnalysis() (protein method)
# Current location: R/protein_de_analysis_wrapper.R
# Type: S4 method (exportMethods)
# Description: Performs differential expression analysis on proteins
# setMethod(f = "differentialExpressionAnalysis", signature = "ProteinQuantitativeData", ...) {
#   # Extract from R/protein_de_analysis_wrapper.R
# }

# Function 2: differentialExpressionAnalysisHelper()
# Current location: R/protein_de_analysis_wrapper.R
# Type: S4 method (exportMethods)
# Description: Helper function for DE analysis
# setMethod(f = "differentialExpressionAnalysisHelper", ...) {
#   # Extract from R/protein_de_analysis_wrapper.R
# }

# Function 3: deAnalysisWrapperFunction()
# Current location: R/de_analysis_function_wrapper.R
# Description: Wrapper function for DE analysis
# deAnalysisWrapperFunction <- function(...) {
#   # Extract from R/de_analysis_function_wrapper.R
# }

# Function 4: outputDeResultsAllContrasts()
# Current location: R/protein_de_analysis_wrapper.R
# Type: S4 method (exportMethods)
# Description: Outputs DE results for all contrasts
# setMethod(f = "outputDeResultsAllContrasts", ...) {
#   # Extract from R/protein_de_analysis_wrapper.R
# }

# Function 5: outputDeAnalysisResults()
# Current location: R/de_analysis_function_wrapper.R
# Description: Outputs DE analysis results
# outputDeAnalysisResults <- function(...) {
#   # Extract from R/de_analysis_function_wrapper.R
# }

# Function 6: generateVolcanoPlotGlimma()
# Current location: R/protein_de_analysis_wrapper.R
# Description: Generates interactive volcano plots using Glimma
# generateVolcanoPlotGlimma <- function(...) {
#   # Extract from R/protein_de_analysis_wrapper.R
# }

# Function 7: generateDEHeatmap()
# Current location: R/protein_de_analysis_wrapper.R
# Description: Generates heatmaps for DE results
# generateDEHeatmap <- function(...) {
#   # Extract from R/protein_de_analysis_wrapper.R
# }

# Function 8: createDeResultsLongFormat()
# Current location: R/de_analysis_function_wrapper.R
# Description: Creates long format DE results
# createDeResultsLongFormat <- function(...) {
#   # Extract from R/de_analysis_function_wrapper.R
# }

# Function 9: getDeResultsLongFormat()
# Current location: R/metaboliteVsSamplesS4Objects.R
# Type: S4 method (exportMethods)
# Description: Gets DE results in long format
# setMethod(f = "getDeResultsLongFormat", ...) {
#   # Extract from R/metaboliteVsSamplesS4Objects.R
# }

# Function 10: getDeResultsWideFormat()
# Current location: R/metaboliteVsSamplesS4Objects.R
# Type: S4 method (exportMethods)
# Description: Gets DE results in wide format
# setMethod(f = "getDeResultsWideFormat", ...) {
#   # Extract from R/metaboliteVsSamplesS4Objects.R
# }

# Function 11: prepareDataForVolcanoPlot()
# Current location: R/de_proteins_functions.R
# Description: Prepares data for volcano plot visualization
# prepareDataForVolcanoPlot <- function(...) {
#   # Extract from R/de_proteins_functions.R
# }

# Function 12: ebFit()
# Current location: R/de_proteins_functions.R
# Description: Empirical Bayes fitting for DE analysis
# ebFit <- function(...) {
#   # Extract from R/de_proteins_functions.R
# }

# Function 13: runTest()
# Current location: R/de_proteins_functions.R
# Description: Runs statistical test for DE
# runTest <- function(...) {
#   # Extract from R/de_proteins_functions.R
# }

# Function 14: runTests()
# Current location: R/de_proteins_functions.R
# Description: Runs multiple statistical tests
# runTests <- function(...) {
#   # Extract from R/de_proteins_functions.R
# }

# Function 15: runTestsContrasts()
# Current location: R/de_proteins_functions.R
# Description: Runs tests for multiple contrasts
# runTestsContrasts <- function(...) {
#   # Extract from R/de_proteins_functions.R
# }

# Function 16: saveDeProteinList()
# Current location: R/de_proteins_functions.R
# Description: Saves list of DE proteins
# saveDeProteinList <- function(...) {
#   # Extract from R/de_proteins_functions.R
# }


# ----------------------------------------------------------------------------
# generateVolcanoPlotGlimma
# ----------------------------------------------------------------------------
##----------------------------------------------------------------------------------------------------------------------------------------------------------------------
#' Generate Interactive Volcano Plot using Glimma
#' 
#' @importFrom purrr map map2 walk walk2 compact imap
#' @importFrom stringr str_extract str_detect
#' @importFrom dplyr filter select mutate bind_rows arrange distinct relocate
#' @importFrom tibble rownames_to_column
#' @importFrom ggplot2 ggplot aes geom_bar geom_text theme_bw element_text labs facet_wrap ggsave
#' @importFrom viridis viridis plasma inferno
#' @importFrom writexl write_xlsx
#' @importFrom vroom vroom_write
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
                                       , gene_names_column = "gene_names"
                                       , display_columns = c( "best_uniprot_acc" )) {

  message("--- Entering generateVolcanoPlotGlimma ---")
  message(paste("   generateVolcanoPlotGlimma Arg: selected_contrast =", selected_contrast))

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

  # CRITICAL FIX: Extract comparison name from selected_contrast if it contains "="
  # The selected_contrast might be the full name like "GA_Elevated.minus.GA_Control=groupGA_Elevated-groupGA_Control"
  # But the comparison column in de_proteins_long contains only "GA_Elevated.minus.GA_Control"
  comparison_to_search <- stringr::str_extract(selected_contrast, "^[^=]+")
  message(paste("   generateVolcanoPlotGlimma: Searching for comparison =", comparison_to_search))
  message(paste("   generateVolcanoPlotGlimma: Original selected_contrast =", selected_contrast))

  # Filter for selected contrast using the extracted comparison name
  contrast_data <- de_proteins_long |>
    dplyr::filter(comparison == comparison_to_search)

  if (nrow(contrast_data) == 0) {
    message(paste("   generateVolcanoPlotGlimma: No data found for contrast", comparison_to_search))
    # DEBUG: Show what comparisons are actually available
    available_comparisons <- unique(de_proteins_long$comparison)
    message(paste("   generateVolcanoPlotGlimma: Available comparisons:", paste(available_comparisons, collapse = ", ")))
    return(NULL)
  }

  message(paste("   generateVolcanoPlotGlimma: Found", nrow(contrast_data), "rows for contrast", comparison_to_search))
  message("   generateVolcanoPlotGlimma Step: Preparing volcano plot data...")

  # Prepare volcano plot table
  volcano_plot_tab <- contrast_data |>
    dplyr::mutate(best_uniprot_acc = str_split(!!sym(args_row_id), ":") |> purrr::map_chr(1))

  # Add UniProt annotations if available
  if (!is.null(uniprot_tbl)) {
    volcano_plot_tab <- volcano_plot_tab |>
      left_join(uniprot_tbl, by = join_by( best_uniprot_acc == !!sym( uniprot_id_column) )) |>
      dplyr::rename( UNIPROT_GENENAME = !!sym(gene_names_column) ) |>
      mutate( UNIPROT_GENENAME = purrr::map_chr( UNIPROT_GENENAME, \(x){
        tryCatch({
          if(is.na(x) || is.null(x) || x == "") {
            ""
          } else {
            split_result <- str_split(x, " |:")[[1]]
            if(length(split_result) > 0) split_result[1] else ""
          }
        }, error = function(e) "")
      })) |>
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
  # CRITICAL FIX: Use grep for reliable pattern matching instead of stringr
  coef_names <- colnames(contrasts_results$fit.eb$coefficients)
  message(paste("   generateVolcanoPlotGlimma: Available coefficient names:", paste(coef_names, collapse = ", ")))
  message(paste("   generateVolcanoPlotGlimma: Looking for pattern:", paste0("^", comparison_to_search, "=")))
  
  # Use grep to find coefficient that starts with friendly name followed by "="
  coef_index <- grep(paste0("^", comparison_to_search, "="), coef_names)
  message(paste("   generateVolcanoPlotGlimma: Pattern match found indices:", paste(coef_index, collapse = ", ")))
  
  # Fallback: try exact match with the friendly name (for backwards compatibility)
  if (length(coef_index) == 0) {
    message("   generateVolcanoPlotGlimma: Pattern match failed, trying exact match with friendly name")
    coef_index <- which(coef_names == comparison_to_search)
  }
  
  # Last resort: try exact match with selected_contrast
  if (length(coef_index) == 0) {
    message("   generateVolcanoPlotGlimma: Exact friendly name match failed, trying selected_contrast")
    coef_index <- which(coef_names == selected_contrast)
  }

  if (length(coef_index) == 0) {
    message(paste("   generateVolcanoPlotGlimma: FINAL FAILURE - No coefficient found for:", selected_contrast))
    message(paste("   generateVolcanoPlotGlimma: Also tried:", comparison_to_search))
    message("   generateVolcanoPlotGlimma: Available coefficient patterns:")
    for (i in seq_along(coef_names)) {
      message(sprintf("     [%d] %s", i, coef_names[i]))
    }
    message(paste("   generateVolcanoPlotGlimma: Pattern attempted:", paste0("^", comparison_to_search, "=")))
    return(NULL)
  }

  message(paste("   generateVolcanoPlotGlimma: Using coefficient index", coef_index, "for contrast", coef_names[coef_index]))

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


# ----------------------------------------------------------------------------
# generateDEHeatmap
# ----------------------------------------------------------------------------
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


# ----------------------------------------------------------------------------
# deAnalysisWrapperFunction
# ----------------------------------------------------------------------------
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
  group_col <- design_matrix[[theObject@group_id]]

  # Check if any group names start with numbers and create mapping
  starts_with_number <- grepl("^[0-9]", group_col)
  if(any(starts_with_number)) {
    original_groups <- unique(group_col)
    safe_groups <- purrr::map_chr(original_groups, \(x) {
      if(grepl("^[0-9]", x)) paste0("grp_", x) else x
    })
    group_mapping <- setNames(original_groups, safe_groups)

    # Update design matrix with safe names
    design_matrix[[theObject@group_id]] <- purrr::map_chr(group_col, \(x) {
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
  if (inherits(theObject, "MetaboliteAssayData")) {
    # Convert metabolite data to matrix format
    metabolite_matrix <- getDataMatrix(theObject)
    return_list$plot_num_of_values <- plotNumOfValuesNoLog(metabolite_matrix)
  } else if(inherits(theObject, "ProteinQuantitativeData") ) {
    return_list$plot_num_of_values <- plotNumOfValuesNoLog(theObject@protein_quant_table)
  }

  ## Compare the different experimental groups and obtain lists of differentially expressed proteins.")

  rownames( theObject@design_matrix ) <- theObject@design_matrix |> dplyr::pull( one_of(theObject@sample_id ))

  # Check if object contains DPC-Quant results and use limpa dpc if available
  use_dpc_de <- FALSE
  dpc_quant_results <- NULL
  
  if (!is.null(theObject@args$limpa_dpc_quant_results)) {
    dpc_quant_results <- theObject@args$limpa_dpc_quant_results
    use_dpc_de <- TRUE
    cat("   DE ANALYSIS Step: Detected DPC-Quant results - using limpa dpcDE for uncertainty-weighted analysis\n")
    cat("   DE ANALYSIS Step: DPC parameters used:", paste(dpc_quant_results$dpc_parameters_used, collapse=", "), "\n")
  } else {
    cat("   DE ANALYSIS Step: No DPC-Quant results found - using standard limma analysis\n")
  }

  # CRITICAL FIX: Use the correct column for contrast strings
  # The downstream functions expect "comparison=expression" format (from full_format column)
  # NOT just the raw contrast expression (from contrasts column)
  if ("full_format" %in% names(contrasts_tbl)) {
    contrast_strings_to_use <- contrasts_tbl$full_format  # Use full_format column
    cat("   DE ANALYSIS Step: Using full_format column for contrast strings\n")
  } else {
    cat("   DE ANALYSIS Step: No full_format column found, auto-generating from raw contrasts\n")
    # Auto-generate full_format column from raw contrasts
    raw_contrasts <- contrasts_tbl[, 1][[1]]
    
    # Generate friendly names and full format
    full_format_strings <- sapply(raw_contrasts, function(contrast_string) {
      # Remove "group" prefixes if present for friendly name
      clean_string <- gsub("^group", "", contrast_string)
      clean_string <- gsub("-group", "-", clean_string)
      
      # Create friendly name by replacing - with _vs_
      friendly_name <- gsub("-", "_vs_", clean_string)
      
      # Create full format: friendly_name=original_contrast_string
      paste0(friendly_name, "=", contrast_string)
    })
    
    contrast_strings_to_use <- full_format_strings
    cat("   DE ANALYSIS Step: Auto-generated full_format strings:\n")
    print(contrast_strings_to_use)
  }

  # Run differential expression analysis
  if (use_dpc_de && requireNamespace("limpa", quietly = TRUE)) {
    cat("   DE ANALYSIS Step: Running limpa dpcDE analysis with uncertainty weights\n")
    
    # The EList is now pre-filtered and stored, so we can use it directly
    y_elist_filtered <- dpc_quant_results$quantified_elist
    
    if(is.null(y_elist_filtered)) {
      stop("FATAL: The quantified_elist is missing from the object's args. It should have been created by proteinMissingValueImputationLimpa.")
    }
    
    # Create design matrix for dpcDE
    design_matrix_for_dpcde <- model.matrix(as.formula(formula_string), theObject@design_matrix)
    
    cat("   DE ANALYSIS Step: Calling limpa::dpcDE\n")
    cat("   DE ANALYSIS Step: Protein matrix dims:", nrow(y_elist_filtered$E), "x", ncol(y_elist_filtered$E), "\n")
    cat("   DE ANALYSIS Step: Design matrix dims:", nrow(design_matrix_for_dpcde), "x", ncol(design_matrix_for_dpcde), "\n")
    
    # Run dpcDE using the synchronized EList
    dpc_fit <- limpa::dpcDE(y_elist_filtered, design_matrix_for_dpcde, plot = FALSE)
    
    # Convert dpcDE results to format compatible with runTestsContrasts
    contrasts_results <- convertDpcDEToStandardFormat(
      dpc_fit = dpc_fit,
      contrast_strings = contrast_strings_to_use,
      design_matrix = design_matrix_for_dpcde,
      eBayes_trend = as.logical(eBayes_trend),
      eBayes_robust = as.logical(eBayes_robust)
    )
    
    cat("   DE ANALYSIS Step: dpcDE analysis completed successfully\n")
    
  } else {
    # Standard limma analysis (existing code)
    cat("   DE ANALYSIS Step: Running standard limma analysis\n")
    
    data_matrix <- getDataMatrix(theObject)

    contrasts_results <- runTestsContrasts(data_matrix,
                                           contrast_strings = contrast_strings_to_use,
                                           design_matrix = theObject@design_matrix,
                                           formula_string = formula_string,
                                           weights = NA,
                                           treat_lfc_cutoff = as.double(treat_lfc_cutoff),
                                           eBayes_trend = as.logical(eBayes_trend),
                                           eBayes_robust = as.logical(eBayes_robust))
  }

  # Extract the contrast name (contains "=" delimiter needed for downstream functions)
  contrast_name <- names(contrasts_results$results)[1]
  message(paste("   DEBUG66: contrast_name extracted =", contrast_name))

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

  # --- NEW: Generate P-value distribution plots for diagnostics ---
  raw_pval_hist_list <- purrr::imap(contrasts_results_table, \(de_tbl, contrast_name) {
    
    # Check if raw_pvalue column exists
    if (!"raw_pvalue" %in% colnames(de_tbl)) {
      warning(paste("raw_pvalue column not found for contrast:", contrast_name))
      return(NULL)
    }
    
    # Create the histogram
    p <- ggplot(de_tbl, aes(x = raw_pvalue)) +
      geom_histogram(aes(y = after_stat(density)), bins = 30, boundary = 0, color = "black", fill = "lightblue") +
      labs(
        title = paste("Raw P-value Distribution for:", contrast_name),
        subtitle = "A uniform distribution (red line) is expected under the null hypothesis.",
        x = "Raw P-value",
        y = "Density"
      ) +
      theme_bw() +
      geom_hline(yintercept = 1, linetype = "dashed", color = "red")
      
    return(p)
  })
  
  # Remove any NULLs from list if a plot failed and add to the main return list
  return_list$raw_pval_histograms <- purrr::compact(raw_pval_hist_list)
  # --- END NEW ---

  ## Prepare data for drawing the volcano plots
  
  # CRITICAL FIX: getSignificantData expects list_of_de_tables to be a list where each element 
  # is itself a list of data.frames. Since we're processing one contrast at a time, we need to 
  # wrap the single data.frame in an additional list layer.
  # CRITICAL FIX 2: The list element name MUST contain "=" delimiter because countStatDeGenesHelper
  # expects to split on "=" to separate comparison name from expression. Use contrast_name which
  # contains the full format like "H4_vs_WT=groupH4-groupWT"
  nested_list <- list(contrasts_results_table)
  names(nested_list) <- contrast_name  # Use actual contrast name with "=" delimiter
  significant_rows <- getSignificantData(list_of_de_tables = list(nested_list),
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
  # CRITICAL FIX: Same as above - wrap in additional list layer and use contrast_name
  nested_list_for_count <- list(contrasts_results_table)
  names(nested_list_for_count) <- contrast_name  # Use actual contrast name with "=" delimiter
  num_sig_de_molecules_first_go <- printCountDeGenesTable(list_of_de_tables = list(nested_list_for_count),
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

  message("--- Checking object type and data access ---")
  message(sprintf("   Object class: %s", class(theObject)))
  
  counts_table_to_use <- getCountsTable(theObject)

  norm_counts <- counts_table_to_use |>
    as.data.frame() |>
    column_to_rownames(args_row_id) |>
    set_colnames(paste0(colnames(counts_table_to_use[-1]), ".log2norm")) |>
    rownames_to_column(args_row_id)

  print(head( norm_counts))

  return_list$norm_counts <- norm_counts

  # Helper to handle protein ID joining based on object type
  dealWithProeinIdJoining <- function(df) {
      if (inherits(theObject, "MetaboliteAssayData")) {
        # For metabolites, we don't need to join with protein_id_table
        df
      } else if (inherits(theObject, "ProteinQuantitativeData")) {
        # For proteins, join with protein_id_table
        df |>
          left_join(theObject@protein_id_table,
                   by = join_by(!!sym(args_row_id) == !!sym(theObject@protein_id_column)))
      } else {
        stop(sprintf("Unsupported object type: %s", class(theObject)))
      }
    }

  de_proteins_wide <- significant_rows |>
    dplyr::filter(analysis_type == "RUV applied") |>
    dplyr::select(-lqm, -colour, -analysis_type) |>
    pivot_wider(id_cols = c(!!sym(args_row_id)),
                names_from = c(comparison),
                names_sep = ":",
                values_from = c(log2FC, !!sym(qvalue_column), !!sym(raw_pvalue_colum))) |>
    left_join(counts_table_to_use, by = join_by( !!sym(args_row_id)  == !!sym(args_row_id)  )   ) |>
    dealWithProeinIdJoining() |>
    dplyr::arrange(across(matches("!!sym(qvalue_column)"))) |>
    distinct()


  return_list$de_proteins_wide <- de_proteins_wide


  ## Create long format output file

  # Create appropriate ID table based on object type
  id_table <- if (inherits(theObject, "MetaboliteAssayData")) {
    message("   Using metabolite data as ID table for MetaboliteAssayData")
    # For metabolites, create a simple ID table from the metabolite data
    theObject@metabolite_data |>
      dplyr::select(!!sym(args_row_id)) |>
      dplyr::distinct()
  } else if (inherits(theObject, "ProteinQuantitativeData")) {
    message("   Using protein_id_table for ProteinQuantitativeData")
    theObject@protein_id_table
  } else {
    stop(sprintf("Unsupported object type: %s", class(theObject)))
  }

  de_proteins_long <- createDeResultsLongFormat( lfc_qval_tbl = significant_rows |>
                                                   dplyr::filter(analysis_type == "RUV applied") ,
                                                 norm_counts_input_tbl = getDataMatrix(theObject),
                                                 raw_counts_input_tbl = getDataMatrix(theObject),
                                                 row_id = args_row_id,
                                                 sample_id = theObject@sample_id,
                                                 group_id = group_id,
                                                 group_pattern = args_group_pattern,
                                                 design_matrix_norm = theObject@design_matrix,
                                                 design_matrix_raw =  theObject@design_matrix,
                                                 protein_id_table = id_table)

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


  # Generate volcano plots with gene names
  cat("DEBUG_66: About to generate volcano plots with gene names\n")
  cat(sprintf("DEBUG_66: static_volcano_plot_data has %d rows\n", nrow(static_volcano_plot_data)))
  cat("DEBUG_66: static_volcano_plot_data structure:\n")
  str(static_volcano_plot_data)
  
  list_of_volcano_plots_with_gene_names <- tryCatch({
    cat("DEBUG_66: Starting volcano plot generation pipeline\n")
    
    # Step 1: Group by comparison
    grouped_data <- static_volcano_plot_data %>%
      group_by(comparison)
    cat(sprintf("DEBUG_66: Grouped data has %d groups\n", n_groups(grouped_data)))
    
    # Step 2: Nest
    nested_data <- grouped_data %>%
      nest()
    cat(sprintf("DEBUG_66: Nested data has %d rows\n", nrow(nested_data)))
    cat("DEBUG_66: Nested data structure:\n")
    str(nested_data)
    
    # Step 3: Ungroup
    ungrouped_data <- nested_data %>%
      ungroup()
    
    # Step 4: Add title
    with_title <- ungrouped_data %>%
      mutate(title = paste(comparison))
    cat("DEBUG_66: Added titles successfully\n")
    
    # Step 5: Add plots with detailed error handling
    with_plots <- with_title %>%
      mutate(plot = purrr::map(data, \(x) {
        cat(sprintf("DEBUG_66: Processing plot for data with %d rows\n", nrow(x)))
        
        # Check uniprot_tbl structure
        cat(sprintf("DEBUG_66: uniprot_tbl has %d rows\n", nrow(uniprot_tbl)))
        cat("DEBUG_66: uniprot_tbl columns:\n")
        print(colnames(uniprot_tbl))
        
        # Check if gene_names column exists
        if (!"gene_names" %in% colnames(uniprot_tbl)) {
          cat("DEBUG_66: ERROR - gene_names column not found in uniprot_tbl!\n")
          cat("DEBUG_66: Available columns in uniprot_tbl:\n")
          print(colnames(uniprot_tbl))
          stop("gene_names column missing from uniprot_tbl")
        }
        
        # Process gene names with detailed debugging
        cat("DEBUG_66: Processing gene names\n")
        uniprot_with_gene_names <- tryCatch({
          uniprot_tbl |>
            mutate(gene_name = purrr::map_chr(gene_names, \(x) {
              tryCatch({
                if(is.na(x) || is.null(x) || x == "") {
                  ""
                } else {
                  split_result <- str_split(x, "; ")[[1]]
                  if(length(split_result) > 0) split_result[1] else ""
                }
              }, error = function(e) {
                cat(sprintf("DEBUG_66: Error processing gene name: %s\n", e$message))
                ""
              })
            }))
        }, error = function(e) {
          cat(sprintf("DEBUG_66: ERROR in gene name processing: %s\n", e$message))
          cat("DEBUG_66: Error details:\n")
          print(e)
          stop(e)
        })
        
        cat("DEBUG_66: Gene names processed successfully\n")
        
        # Call plot function
        printOneVolcanoPlotWithProteinLabel(
          input_table = x,
          uniprot_table = uniprot_with_gene_names,
          protein_id_column = protein_id_column,
          input_title = "",
          fdr_threshold = de_q_val_thresh,
          number_of_genes = 10
        )
      }))
    
    cat("DEBUG_66: All plots generated successfully\n")
    with_plots
    
  }, error = function(e) {
    cat(sprintf("DEBUG_66: CRITICAL ERROR in volcano plot generation: %s\n", e$message))
    cat("DEBUG_66: Full error object:\n")
    print(e)
    cat("DEBUG_66: Traceback:\n")
    print(traceback())
    NULL
  })

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


# ----------------------------------------------------------------------------
# outputDeAnalysisResults
# ----------------------------------------------------------------------------
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

  cat("*** ENTERING outputDeAnalysisResults ***\n")
  cat("DEBUG: file_prefix =", file_prefix, "\n")
  cat("DEBUG: de_output_dir =", de_output_dir, "\n")
  cat("DEBUG: uniprot_tbl is null:", is.null(uniprot_tbl), "\n")
  cat("DEBUG: de_analysis_results_list names:", paste(names(de_analysis_results_list), collapse=", "), "\n")


  uniprot_tbl <- checkParamsObjectFunctionSimplify(theObject, "uniprot_tbl", NULL)
  de_output_dir <- checkParamsObjectFunctionSimplify(theObject, "de_output_dir", NULL)
  publication_graphs_dir <- checkParamsObjectFunctionSimplify(theObject, "publication_graphs_dir", NULL)
  file_prefix <- checkParamsObjectFunctionSimplify(theObject, "file_prefix", "de_proteins")
  plots_format <- checkParamsObjectFunctionSimplify(theObject, "plots_format", c("pdf", "png"))
  args_row_id <- checkParamsObjectFunctionSimplify(theObject, "args_row_id", "uniprot_acc")
  de_q_val_thresh <- checkParamsObjectFunctionSimplify(theObject, "de_q_val_thresh", 0.05)
  gene_names_column <- checkParamsObjectFunctionSimplify(theObject, "gene_names_column", "gene_names")
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

  cat("DEBUG: Finished parameter setup, starting plot creation\n")

  ## PCA plot
  cat("DEBUG: Starting PCA plot section\n")
  plot_pca_plot <- de_analysis_results_list$pca_plot

  dir.create(file.path( publication_graphs_dir, "PCA")
             , recursive = TRUE
             , showWarnings = FALSE)
  cat("DEBUG: PCA plot section completed\n")

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
    dplyr::select(-colour, -lqm) |>
    vroom::vroom_write(file.path(de_output_dir, "lfc_qval_long.tsv"))

  significant_rows |>
    dplyr::select(-colour, -lqm) |>
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
  cat("DEBUG: Starting wide format output section\n")
  de_proteins_wide <- de_analysis_results_list$de_proteins_wide
  cat("DEBUG: de_proteins_wide extracted successfully\n")
  
  vroom::vroom_write( de_proteins_wide,
                      file.path( de_output_dir,
                                 paste0(file_prefix, "_wide.tsv")))

  writexl::write_xlsx( de_proteins_wide,
                       file.path( de_output_dir,
                                  paste0(file_prefix, "_wide.xlsx")))
  cat("DEBUG: Wide format files written\n")

  cat("DEBUG: Starting wide_annot creation\n")
  
  tryCatch({
    de_proteins_wide_annot <- de_proteins_wide |>
      mutate( uniprot_acc_cleaned = str_split( !!sym(args_row_id), "-" )  |>
                purrr::map_chr(1) ) |>
      left_join( uniprot_tbl, by = join_by( uniprot_acc_cleaned == Entry ) ) |>
      dplyr::select( -uniprot_acc_cleaned)  |>
      mutate( gene_name = ifelse(
        !is.na(gene_names) & gene_names != "", 
        sapply(gene_names, function(x) {
          if(is.na(x) || x == "") return("")
          split_result <- strsplit(as.character(x), " |:")[[1]]
          if(length(split_result) > 0) split_result[1] else ""
        }),
        ""
      )) |>
      relocate( gene_name, .after = !!sym(args_row_id))
    
    cat("DEBUG: wide_annot creation successful\n")
    
  }, error = function(e) {
    cat("DEBUG: ERROR in wide_annot creation:", e$message, "\n")
    # Create a fallback version without annotations
    de_proteins_wide_annot <- de_proteins_wide |>
      mutate(gene_name = "")
    cat("DEBUG: Created fallback wide_annot without annotations\n")
  })

  vroom::vroom_write( de_proteins_wide_annot,
                      file.path( de_output_dir,
                                 paste0(file_prefix, "_wide_annot.tsv")))

  writexl::write_xlsx( de_proteins_wide_annot,
                       file.path( de_output_dir,
                                  paste0(file_prefix, "_wide_annot.xlsx")))

  cat("DEBUG: Wide_annot files written, proceeding to long format\n")

    ## Create long format output file
  de_proteins_long <- de_analysis_results_list$de_proteins_long
  
  cat("DEBUG: de_proteins_long exists:", !is.null(de_proteins_long), "\n")
  if(!is.null(de_proteins_long)) {
    cat("DEBUG: de_proteins_long dimensions:", dim(de_proteins_long), "\n")
  } else {
    cat("DEBUG: de_proteins_long is NULL - cannot create long_annot\n")
    return(NULL)
  }
  
  vroom::vroom_write( de_proteins_long,
                      file.path( de_output_dir,
                                 paste0(file_prefix, "_long.tsv")))

  writexl::write_xlsx( de_proteins_long,
                       file.path( de_output_dir,
                                 paste0(file_prefix, "_long.xlsx")))

  cat("DEBUG: Starting long_annot creation\n")
  cat("DEBUG: de_proteins_long dimensions:", dim(de_proteins_long), "\n")
  cat("DEBUG: uniprot_tbl is null:", is.null(uniprot_tbl), "\n")
  if(!is.null(uniprot_tbl)) {
    cat("DEBUG: uniprot_tbl dimensions:", dim(uniprot_tbl), "\n")
  }
  cat("DEBUG: args_row_id:", args_row_id, "\n")

  de_proteins_long_annot <- de_proteins_long |>
    mutate( uniprot_acc_cleaned = str_split( !!sym(args_row_id), "-" )  |>
              purrr::map_chr(1) )|>
    left_join( uniprot_tbl, by = join_by( uniprot_acc_cleaned == Entry ) )  |>
    dplyr::select( -uniprot_acc_cleaned)  |>
    mutate( gene_name = ifelse(
      !is.na(gene_names) & gene_names != "", 
      sapply(gene_names, function(x) {
        if(is.na(x) || x == "") return("")
        split_result <- strsplit(as.character(x), " |:")[[1]]
        if(length(split_result) > 0) split_result[1] else ""
      }),
      ""
    )) |>
    relocate( gene_name, .after = !!sym(args_row_id))

  cat("DEBUG: de_proteins_long_annot dimensions:", dim(de_proteins_long_annot), "\n")
  cat("DEBUG: de_proteins_long_annot is null:", is.null(de_proteins_long_annot), "\n")
  if(!is.null(de_proteins_long_annot) && nrow(de_proteins_long_annot) > 0) {
    cat("DEBUG: long_annot columns:", paste(names(de_proteins_long_annot), collapse=", "), "\n")
  } else {
    cat("DEBUG: long_annot is empty or null - annotation pipeline failed\n")
  }

  long_annot_file_path <- file.path( de_output_dir, paste0(file_prefix, "_long_annot.tsv"))
  cat("DEBUG: Attempting to write long_annot to:", long_annot_file_path, "\n")

  vroom::vroom_write( de_proteins_long_annot, long_annot_file_path)
  
  cat("DEBUG: long_annot file written, checking if exists:", file.exists(long_annot_file_path), "\n")

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

    savePlot(num_sig_de_genes_barplot_only_significant,
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


# ----------------------------------------------------------------------------
# prepareDataForVolcanoPlot
# ----------------------------------------------------------------------------
## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#' @title Prepare data for volcano plot
#' @param input_table The input table with the log fold-change and q-value columns.
#' @param protein_id_column The name of the column representing the protein ID (tidyverse style).
#' @param uniprot_table The uniprot table with the imporatnt info on each protein
#' @param uniprot_protein_id_column The name of the column representing the protein ID in the uniprot table (tidyverse style).
#' @param number_of_genes The number of genes to show in the volcano plot.
#' @param fdr_threshold The FDR threshold for the volcano plot.
#' @param fdr_column The name of the column representing the FDR value (tidyverse style).
#' @param log2FC_column The name of the column representing the log fold-change (tidyverse style).
#' @return A table with the following columns:
#' label  The label of the significant proteins.
#' log2FC The log2 fold-change of the significant proteins.
#' lqm  The -log10 of the q-value.
#' colour The colour of the significant proteins.
#' rank_positive: The rank of the positive fold-change values.
#' rank_negative: The rank of the negative fold-change values.
#' gene_name_significant  The gene name of the significant proteins.
#'
#' @export
prepareDataForVolcanoPlot <- function(input_table
                                      , protein_id_column = uniprot_acc
                                      , uniprot_table
                                      , uniprot_protein_id_column = uniprot_acc_first
                                      , gene_name_column = gene_name
                                      , number_of_genes = 3000
                                      , fdr_threshold = 0.05
                                      , fdr_column = q.mod
                                      , log2FC_column = log2FC){

  temp_col_name <-  as_string(as_name(enquo(protein_id_column)))

  # Check if we're dealing with metabolite data (no protein IDs)
  # This handles the case where the input table is for metabolites and doesn't have Protein.Ids or uniprot_acc
  is_metabolite_data <- !any(grepl("Protein.Ids|uniprot_acc", names(input_table)))

  if (is_metabolite_data) {
    message("   Processing metabolite data (no UniProt ID mapping)...")
    proteomics_volcano_tbl <- input_table |>
      dplyr::select({{protein_id_column}}, {{fdr_column}}, {{log2FC_column}}) |>
      mutate(
        colour = case_when(
          {{fdr_column}} < fdr_threshold & {{log2FC_column}} > 0 ~ "red",
          {{fdr_column}} < fdr_threshold & {{log2FC_column}} < 0 ~ "blue",
          TRUE ~ "grey"
        ),
        lqm = -log10({{fdr_column}}),
        label = case_when(
          {{fdr_column}} < fdr_threshold & {{log2FC_column}} > 0 ~ "Significant Increase",
          {{fdr_column}} < fdr_threshold & {{log2FC_column}} < 0 ~ "Significant Decrease",
          TRUE ~ "Not significant"
        ),
        label = factor(label, levels = c(
          "Significant Increase",
          "Significant Decrease",
          "Not significant"
        )),
        rank_positive = case_when(
          {{log2FC_column}} > 0 ~ {{fdr_column}},
          TRUE ~ NA_real_
        ) |> rank(),
        rank_negative = case_when(
          {{log2FC_column}} < 0 ~ {{fdr_column}},
          TRUE ~ NA_real_
        ) |> rank(),
        # For metabolites, we use the protein_id_column (usually Name) as the "gene name"
        gene_name_significant = case_when(
          {{fdr_column}} < fdr_threshold &
            (rank_positive <= number_of_genes |
               rank_negative <= number_of_genes) ~ as.character({{protein_id_column}}),
          TRUE ~ NA_character_
        )
      )
  } else {
    message("   Processing protein data (performing UniProt ID mapping)...")
    proteomics_volcano_tbl <- input_table |>
      dplyr::mutate( uniprot_acc_first  = purrr::map_chr( {{protein_id_column}}, \(x) { str_split( x, ":")[[1]][1] } ) ) |>
      dplyr::relocate( uniprot_acc_first, .after=temp_col_name) |>
      dplyr::select (uniprot_acc_first, {{fdr_column}}, {{log2FC_column}}) |>
      left_join(uniprot_table
                , by = join_by( uniprot_acc_first == {{uniprot_protein_id_column}})) |>
      mutate( colour = case_when ( {{fdr_column}} < fdr_threshold & {{log2FC_column}} > 0 ~ "red"
                                   , {{fdr_column}} < fdr_threshold & {{log2FC_column}} < 0 ~ "blue"
                                   , TRUE ~ "grey" )) |>
      mutate( lqm = -log10({{fdr_column}})) |>
      mutate( label = case_when (  {{fdr_column}} < fdr_threshold & {{log2FC_column}} > 0 ~ "Significant Increase"
                                   , {{fdr_column}} < fdr_threshold & {{log2FC_column}} < 0 ~ "Significant Decrease"
                                   , TRUE ~ "Not significant" )) |>
      mutate( label = factor( label, levels = c( "Significant Increase"
                                                 , "Significant Decrease"
                                                 , "Not significant" ))) |>
      mutate ( rank_positive = case_when( {{log2FC_column}} > 0 ~ {{fdr_column}}
                                          , TRUE ~ NA_real_) |> rank() ) |>
      mutate ( rank_negative = case_when( {{log2FC_column}} < 0 ~ {{fdr_column}}
                                          , TRUE ~ NA_real_) |> rank() ) |>
      mutate ( {{gene_name_column}} := purrr::map_chr( {{gene_name_column}}, \(x) { str_split( x, " ")[[1]][1] } ) ) |>
      mutate( gene_name_significant = case_when( {{fdr_column}} < fdr_threshold &
                                                   ( rank_positive <= number_of_genes |
                                                     rank_negative <= number_of_genes ) ~  {{gene_name_column}}
                                                 , TRUE ~ NA ) )
  }

  proteomics_volcano_tbl
}


# ----------------------------------------------------------------------------
# ebFit
# ----------------------------------------------------------------------------
#' Run the Empircal Bayes Statistics for Differential Expression in the limma package
#' @param ID List of protein accessions / row names.
#' @param design Output from running the function \code{\link{model.matrix}}.
#' @param contr.matrix Output from the function \code{\link{makeContrasts}}.
#' @seealso \code{\link{model.matrix}}
#' @seealso \code{\link{makeContrasts}}
#' @export
ebFit <- function(data, design, contr.matrix)
{
  fit <- lmFit(data, design)
  fit.c <- contrasts.fit(fit, contrasts = contr.matrix)

  fit.eb <- suppressWarnings(eBayes(fit.c))

  logFC <- fit.eb$coefficients[, 1]
  df.r <- fit.eb$df.residual
  df.0 <- rep(fit.eb$df.prior, dim(data)[1])
  s2.0 <- rep(fit.eb$s2.prior, dim(data)[1])
  s2 <- (fit.eb$sigma)^2
  s2.post <- fit.eb$s2.post
  t.ord <- fit.eb$coefficients[, 1] /
    fit.eb$sigma /
    fit.eb$stdev.unscaled[, 1]
  t.mod <- fit.eb$t[, 1]
  p.ord <- 2 * pt(-abs(t.ord), fit.eb$df.residual)
  raw_pvalue <- fit.eb$p.value[, 1]
  q.ord <- qvalue(p.ord)$q
  fdr_qvalue <- qvalue(raw_pvalue)$q

  return(list(table = data.frame(logFC, t.ord, t.mod, p.ord, raw_pvalue, q.ord, fdr_qvalue, df.r, df.0, s2.0, s2, s2.post),
              fit.eb = fit.eb))
}


# ----------------------------------------------------------------------------
# runTest
# ----------------------------------------------------------------------------
## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#' Analyse one contrast (e.g. compare a pair of experimental groups) and output the q-values per protein.
#'@param ID List of protein accessions / row names.
#'@param A String representing the name of experimental group A for pairwise comparison of B - A.
#'@param B String representing the name of experimental group B for pairwise comparison of B - A.
#'@param group_A Names of all the columns / samples that are in experimental group A.
#'@param group_B Names of all the columns / samples that are in experimental group B.
#'@param design_matrix A data frame with a column containing the sample ID (as per the sample_id param) and the experimental group (as per the group param). Each row as the sample ID as row name in the data frame.
#'@param formula_string A formula string representing the experimental design. e.g. ("~ 0 + group")
#'@param contrast_variable String representing the contrast variable, which is also used in the formula string. (e.g. "group")
#'@param weights Numeric matrix for adjusting each sample and gene.
#'@return A data frame with the following columns:
#' row.names = the protein accessions
#' comparison A string showing log({group B's name}) minus log({group A's name})
#' meanA     mean of the normalised log abundance value of the gene across samples from experimental group A
#' meanB     mean of the normalised log abundance value of the gene across samples from experimental group B
#' logFC     log fold-change
#' tstats    t-test statistics
#' tmod      moderated t-test statistics
#' pval      t-test p-value
#' raw_pvalue      moderated t-test p-value
#' qval      t-test q-value
#' fdr_qvalue     moderated t-test q-value
#'@export
runTest <- function(ID, A, B, group_A, group_B, design_matrix, formula_string,
                    contrast_variable = "group",
                    weights = NA) {


  ff <- as.formula(formula_string)
  mod_frame <- model.frame(ff, design_matrix)
  design_m <- model.matrix(ff, mod_frame)


  # print("My design matrix")
  # print(design_m)
  # print( paste( "nrow(weights)", nrow(weights), "nrow(design_m)", nrow(design_m)))

  if (!is.na(weights)) {
    if (nrow(weights) == nrow(design_m)) {
      design_m <- cbind(design_m, weights)
    } else {
      stop("Stop: nrow(weights) should be equal to nrow(design_m)")
    }

  }

  # print(paste("group_A = ", group_A))
  # print(paste("group_B = ", group_B))

  contr.matrix <- makeContrasts(contrasts = paste0(group_B, "vs", group_A, "=", contrast_variable, group_B, "-", contrast_variable, group_A),
                                levels = colnames(design_m))

  eb_fit_list <- ebFit(cbind(A, B), design_m, contr.matrix = contr.matrix)

  r <- eb_fit_list$table
  fit.eb <- eb_fit_list$fit.eb

  return(list(table = data.frame(row.names = row.names(r),
                                 comparison = paste("log(", group_B, ") minus log(", group_A, ")", sep = ""),
                                 meanA = rowMeans(A),
                                 meanB = rowMeans(B),
                                 logFC = r$logFC,
                                 tstats = r$t.ord,
                                 tmod = r$t.mod,
                                 pval = r$p.ord,
                                 raw_pvalue = r$raw_pvalue,
                                 qval = r$q.ord,
                                 fdr_qvalue = r$fdr_qvalue),
              fit.eb = fit.eb))
}


# ----------------------------------------------------------------------------
# runTests
# ----------------------------------------------------------------------------
#' Compare a pair of experimental groups and output the log fold-change and q-values per protein.
#'@param ID List of protein accessions / row names.
#'@param data Data frame containing the log (base 2) protein abundance values where each column represents a sample and each row represents a protein group, and proteins as rows. The data is preferably median-scaled, with missing values imputed, and batch-effects removed.
#'@param test_pairs Input file with a table listing all the pairs of experimental groups to compare. First column represents group A and second column represents group B. Linear model comparisons (e.g. Contrasts) would be group B minus group A.
#'@param sample_columns A vector of column names (e.g. strings) representing samples which would be used in the statistical tests. Each column contains protein abundance values.
#'@param sample_rows_list A list, the name of each element is the sample ID and each element is a vector containing the protein accessions (e.g. row_id) with enough number of values. It is usually the output from the function \code{get_rows_to_keep_list}.
#'@param type_of_grouping A list where each element name is the name of a treatment group and each element is a vector containing the sample IDs within the treatment group. It is usually the output from the function \code{get_type_of_grouping}.
#'@param design_matrix A data frame with a column containing the sample ID (as per the sample_id param) and the experimental group (as per the group param). Each row as the sample ID as row name in the data frame.
#'@param formula_string A formula string representing the experimental design. e.g. ("~ 0 + group")
#'@param contrast_variable String representing the contrast variable, which is also used in the formula string. (e.g. "group")
#'@param weights Numeric matrix for adjusting each sample and gene.
#'@return A list of data frames, the name of each element represents each pairwise comparison. Each data frame has the following columns:
#' row.names = the protein accessions
#' comparison A string showing log({group B's name}) minus log({group A's name})
#' meanA     mean of the normalised log abundance value of the gene across samples from experimental group A
#' meanB     mean of the normalised log abundance value of the gene across samples from experimental group B
#' logFC     log fold-change
#' tstats    t-test statistics
#' tmod      moderated t-test statistics
#' pval      t-test p-value
#' raw_pvalue      moderated t-test p-value
#' qval      t-test q-value
#' fdr_qvalue     moderated t-test q-value
#' @seealso \code{\link{get_rows_to_keep_list}}
#' @seealso \code{\link{get_type_of_grouping}}
#'@export
runTests <- function(ID, data, test_pairs, sample_columns, sample_rows_list = NA, type_of_grouping, design_matrix, formula_string, contrast_variable = "group", weights = NA) {
  r <- list()
  for (i in 1:nrow(test_pairs)) {

    rows_to_keep <- rownames(data)


    if (length(sample_rows_list) > 0) {
      if (!is.na(sample_rows_list) &
        #  Check that sample group exists as names inside sample_rows_list
        length(which(c(test_pairs[i, "A"], test_pairs[i, "B"]) %in% names(sample_rows_list))) > 0) {

        rows_to_keep <- unique(sample_rows_list[[ test_pairs[[i, "A"]]]],
                               sample_rows_list[[ test_pairs[[i, "B"]]]])
      }
    }

    tmp <- data[rows_to_keep, sample_columns]
    rep <- colnames(tmp)

    # print( paste( test_pairs[i,]$A, test_pairs[i,]$B) )
    A <- tmp[, type_of_grouping[test_pairs[i,]$A][[1]]]
    B <- tmp[, type_of_grouping[test_pairs[i,]$B][[1]]]

    subset_weights <- NA

    if (!is.na(weights)) {
      subset_weights <- weights[c(colnames(A), colnames(B)),]
    }

    # print(colnames(A))
    # print(colnames(B))
    tmp <- unname(cbind(A, B))
    Aname <- paste(test_pairs[i,]$A, 1:max(1, ncol(A)), sep = "_")
    Bname <- paste(test_pairs[i,]$B, 1:max(1, ncol(B)), sep = "_")
    colnames(tmp) <- c(Aname, Bname)

    selected_sample_ids <- c(type_of_grouping[test_pairs[i,]$A][[1]], type_of_grouping[test_pairs[i,]$B][[1]])
    design_matrix_subset <- design_matrix[selected_sample_ids, , drop = FALSE]

    # print("My design matrix 1")
    # print( selected_sample_ids)
    # print(design_matrix)
    # print( dim(design_matrix))

    group_A <- test_pairs[i,]$A
    group_B <- test_pairs[i,]$B

    x <- runTest(ID, A, B, group_A, group_B, design_matrix = design_matrix_subset,
                 formula_string = formula_string, contrast_variable = contrast_variable,
                 weights = subset_weights)

    comparison <- paste(group_B, " vs ", group_A, sep = "")

    r[[comparison]] <- list(results = x$table, counts = t(cbind(A, B)), fit.eb = x$fit.eb)
  }
  r
}


# ----------------------------------------------------------------------------
# runTestsContrasts
# ----------------------------------------------------------------------------
#' Run the linear model fitting and statistical tests for a set of contrasts, then adjust with Empirical Bayes function
#'@param data Data frame containing the log (base 2) protein abundance values where each column represents a sample and each row represents a protein group, and proteins as rows. The data is preferably median-scaled, with missing values imputed, and batch-effects removed.
#'@param contrast_strings Input file with a table listing all the experimental contrasts to analyse. It will be in the format required for the function \code{makeContrasts} in the limma package.
#'The contrast string consists of variable that each consist of concatenating the column name (e.g. group) and the string representing the group type (e.g. A) in the design matrix.
#'@param design_matrix A data frame with a column containing the sample ID (as per the sample_id param) and the experimental group (as per the group param). Each row as the sample ID as row name in the data frame.
#'@param formula_string A formula string representing the experimental design. e.g. ("~ 0 + group")
#' @param p_value_column The name of the raw p-value column (tidyverse style).
#' @param q_value_column The name of the q-value column (tidyverse style).
#' @param fdr_value_column The name of the fdr-value column (tidyverse style).
#'@return A list containing two elements. $results returns a list of tables containing logFC and q-values. $fit.eb returns the Empiracle Bayes output object.
#'@export
runTestsContrasts <- function(data,
                              contrast_strings,
                              design_matrix,
                              formula_string,
                              p_value_column = raw_pvalue,
                              q_value_column = fdr_qvalue,
                              fdr_value_column = fdr_value_bh_adjustment,
                              weights = NA,
                              treat_lfc_cutoff = NA,
                              eBayes_trend = FALSE,
                              eBayes_robust = FALSE) {

  message("--- Entering runTestsContrasts ---")
  message(sprintf("   runTestsContrasts: data dims = %d x %d, %d contrasts", nrow(data), ncol(data), length(contrast_strings)))
  message(sprintf("   runTestsContrasts: contrasts = %s", paste(contrast_strings, collapse=", ")))
  message(sprintf("   runTestsContrasts: treat_lfc_cutoff = %s", treat_lfc_cutoff))

  # Create formula and design matrix
  ff <- as.formula(formula_string)
  mod_frame <- model.frame(ff, design_matrix)
  design_m <- model.matrix(ff, mod_frame)
  message(sprintf("   runTestsContrasts: design_m dims = %d x %d", nrow(design_m), ncol(design_m)))

  # Subset data to match design matrix
  data_subset <- data[, rownames( design_m)]
  message(sprintf("   runTestsContrasts: data_subset dims = %d x %d", nrow(data_subset), ncol(data_subset)))

  # Create contrast matrix
  message("   runTestsContrasts: Creating contrast matrix...")
  contr.matrix <- makeContrasts(contrasts = contrast_strings,
                                levels = colnames(design_m))
  message(sprintf("   runTestsContrasts: contr.matrix dims = %d x %d", nrow(contr.matrix), ncol(contr.matrix)))

  # Check weights
  if (!is.na(weights)) {
    message("   runTestsContrasts: Attaching weights...")
    if (nrow(weights) == nrow(design_m)) {
      design_m <- cbind(design_m, weights)
    } else {
      stop("Stop: nrow(weights) should be equal to nrow(design_m)")
    }
  }

  # Run limma analysis
  message("   runTestsContrasts: Running lmFit...")
  fit <- lmFit(data_subset, design = design_m)
  
  message("   runTestsContrasts: Running contrasts.fit...")
  cfit <- contrasts.fit(fit, contrasts = contr.matrix)
  
  message("   runTestsContrasts: Running eBayes...")
  eb.fit <- eBayes( cfit, trend = eBayes_trend, robust = eBayes_robust )

  # Run treat or standard analysis
  t.fit <- NA
  result_tables <- NA
  if( !is.na( treat_lfc_cutoff)) {
    message("   runTestsContrasts: Running treat analysis...")
    t.fit <- treat(eb.fit, lfc=as.double(treat_lfc_cutoff))
    
    message(sprintf("   runTestsContrasts: Processing %d contrasts with topTreat...", length(contrast_strings)))
    # Track qvalue failures for user notification
    qvalue_failures <- list()
    
    result_tables <- purrr::map(contrast_strings,
                                function(contrast) {
                                  message(sprintf("      [map] Processing contrast: %s", contrast))
                                  qvalue_failed <- FALSE
                                  
                                  tryCatch({
                                    message(sprintf("      About to call topTreat with coef = %s", contrast))
                                    de_tbl <- topTreat(t.fit, coef = contrast, n = Inf)
                                    message(sprintf("      [map] topTreat success: %d rows", nrow(de_tbl)))
                                    
                                    message("      Adding qvalue column...")
                                    # Safe qvalue computation: handle invalid p-values (NA, Inf, NaN)
                                    valid_p_idx <- which(!is.na(de_tbl$P.Value) & is.finite(de_tbl$P.Value))
                                    if (length(valid_p_idx) > 0) {
                                      # Compute q-values only for valid p-values
                                      valid_p_values <- de_tbl$P.Value[valid_p_idx]
                                      
                                      # Diagnostic: Log p-value distribution statistics
                                      message(sprintf("      Diagnostic: Valid p-values: %d of %d total", length(valid_p_idx), nrow(de_tbl)))
                                      message(sprintf("      Diagnostic: P-value range: [%.6f, %.6f]", min(valid_p_values), max(valid_p_values)))
                                      message(sprintf("      Diagnostic: P-value mean: %.6f, median: %.6f", mean(valid_p_values), median(valid_p_values)))
                                      
                                      # Edge case checks that might cause qvalue() to fail
                                      all_zeros <- all(valid_p_values == 0)
                                      all_ones <- all(valid_p_values == 1)
                                      too_few <- length(valid_p_values) < 3
                                      
                                      if (all_zeros) {
                                        message("      Warning: All p-values are 0 - qvalue() cannot compute, using p.adjust()")
                                        use_qvalue <- FALSE
                                      } else if (all_ones) {
                                        message("      Warning: All p-values are 1 - qvalue() may fail, using p.adjust()")
                                        use_qvalue <- FALSE
                                      } else if (too_few) {
                                        message(sprintf("      Warning: Too few p-values (%d < 3) for qvalue() estimation, using p.adjust()", length(valid_p_values)))
                                        use_qvalue <- FALSE
                                      } else {
                                        use_qvalue <- TRUE
                                      }
                                      
                                      q_values_all <- rep(NA_real_, nrow(de_tbl))
                                      if (use_qvalue) {
                                        tryCatch({
                                          q_values_valid <- qvalue(valid_p_values)$q
                                          q_values_all[valid_p_idx] <- q_values_valid
                                          message("      qvalue() computation successful")
                                        }, error = function(e) {
                                          qvalue_failed <<- TRUE
                                          message(sprintf("      Warning: qvalue() failed during computation: %s", e$message))
                                          message(sprintf("      Diagnostic: P-value distribution may be problematic for qvalue smoothing algorithm"))
                                          message(sprintf("      Diagnostic: Falling back to p.adjust() method='BH' (Benjamini-Hochberg FDR)"))
                                          # Fallback to p.adjust if qvalue fails
                                          q_values_all[valid_p_idx] <- p.adjust(valid_p_values, method = "BH")
                                          message(sprintf("      Diagnostic: Assigned %d p.adjust() values to q-value column", length(valid_p_idx)))
                                        })
                                      } else {
                                        # Use p.adjust directly for edge cases
                                        qvalue_failed <<- TRUE
                                        q_values_all[valid_p_idx] <- p.adjust(valid_p_values, method = "BH")
                                        message("      Using p.adjust() due to edge case detection")
                                        message(sprintf("      Diagnostic: Assigned %d p.adjust() values to q-value column", length(valid_p_idx)))
                                      }
                                    } else {
                                      # All p-values are invalid, set all q-values to NA
                                      q_values_all <- rep(NA_real_, nrow(de_tbl))
                                      message("      Warning: All p-values are invalid (NA, Inf, or NaN), setting q-values to NA")
                                    }
                                    
                                    # Debug: Verify assignment before mutate
                                    message(sprintf("      Diagnostic: q_values_all has %d non-NA values before assignment", sum(!is.na(q_values_all))))
                                    
                                    de_tbl <- de_tbl |>
                                      mutate({ { q_value_column } } := q_values_all)
                                    
                                    # Debug: Verify assignment after mutate
                                    assigned_col <- de_tbl[[rlang::as_name(rlang::ensym(q_value_column))]]
                                    message(sprintf("      Diagnostic: Assigned column has %d non-NA values after mutate", sum(!is.na(assigned_col))))
                                    message("      qvalue column added")
                                    
                                    message("      Adding FDR column...")
                                    # Use the same safe logic for FDR column - only compute for valid p-values
                                    fdr_values_all <- rep(NA_real_, nrow(de_tbl))
                                    if (length(valid_p_idx) > 0) {
                                      fdr_values_all[valid_p_idx] <- p.adjust(valid_p_values, method = "BH")
                                    }
                                    de_tbl <- de_tbl |>
                                      mutate({ { fdr_value_column } } := fdr_values_all)
                                    message("      FDR column added")
                                    
                                    message("      Renaming P.Value column...")
                                    de_tbl <- de_tbl |>
                                      dplyr::rename({ { p_value_column } } := P.Value)
                                    message("      P.Value column renamed")
                                    
                                    message(sprintf("   [map] Completed processing contrast: %s", contrast))
                                    if (qvalue_failed) {
                                      qvalue_failures[[contrast]] <<- TRUE
                                    }
                                    return(de_tbl)
                                    
                                  }, error = function(e) {
                                    message(sprintf("      [map] ERROR in contrast %s: %s", contrast, e$message))
                                    message(sprintf("      [map] ERROR call stack: %s", capture.output(traceback())))
                                    stop(e)
                                  })
                                }
    )
  } else {
    message("   runTestsContrasts: Running standard analysis...")
    t.fit <- eb.fit
    
    message(sprintf("   runTestsContrasts: Processing %d contrasts with topTable...", length(contrast_strings)))
    # Track qvalue failures for user notification
    qvalue_failures <- list()
    
    result_tables <- purrr::map(contrast_strings,
                                function(contrast) {
                                  message(sprintf("      [map] Processing contrast: %s", contrast))
                                  qvalue_failed <- FALSE
                                  
                                  tryCatch({
                                    message(sprintf("      About to call topTable with coef = %s", contrast))
                                    de_tbl <- topTable(t.fit, coef = contrast, n = Inf)
                                    message(sprintf("      [map] topTable success: %d rows", nrow(de_tbl)))
                                    
                                    message("      Adding qvalue column...")
                                    # Safe qvalue computation: handle invalid p-values (NA, Inf, NaN)
                                    valid_p_idx <- which(!is.na(de_tbl$P.Value) & is.finite(de_tbl$P.Value))
                                    if (length(valid_p_idx) > 0) {
                                      # Compute q-values only for valid p-values
                                      valid_p_values <- de_tbl$P.Value[valid_p_idx]
                                      
                                      # Diagnostic: Log p-value distribution statistics
                                      message(sprintf("      Diagnostic: Valid p-values: %d of %d total", length(valid_p_idx), nrow(de_tbl)))
                                      message(sprintf("      Diagnostic: P-value range: [%.6f, %.6f]", min(valid_p_values), max(valid_p_values)))
                                      message(sprintf("      Diagnostic: P-value mean: %.6f, median: %.6f", mean(valid_p_values), median(valid_p_values)))
                                      
                                      # Edge case checks that might cause qvalue() to fail
                                      all_zeros <- all(valid_p_values == 0)
                                      all_ones <- all(valid_p_values == 1)
                                      too_few <- length(valid_p_values) < 3
                                      
                                      if (all_zeros) {
                                        message("      Warning: All p-values are 0 - qvalue() cannot compute, using p.adjust()")
                                        use_qvalue <- FALSE
                                      } else if (all_ones) {
                                        message("      Warning: All p-values are 1 - qvalue() may fail, using p.adjust()")
                                        use_qvalue <- FALSE
                                      } else if (too_few) {
                                        message(sprintf("      Warning: Too few p-values (%d < 3) for qvalue() estimation, using p.adjust()", length(valid_p_values)))
                                        use_qvalue <- FALSE
                                      } else {
                                        use_qvalue <- TRUE
                                      }
                                      
                                      q_values_all <- rep(NA_real_, nrow(de_tbl))
                                      if (use_qvalue) {
                                        tryCatch({
                                          q_values_valid <- qvalue(valid_p_values)$q
                                          q_values_all[valid_p_idx] <- q_values_valid
                                          message("      qvalue() computation successful")
                                        }, error = function(e) {
                                          qvalue_failed <<- TRUE
                                          message(sprintf("      Warning: qvalue() failed during computation: %s", e$message))
                                          message(sprintf("      Diagnostic: P-value distribution may be problematic for qvalue smoothing algorithm"))
                                          message(sprintf("      Diagnostic: Falling back to p.adjust() method='BH' (Benjamini-Hochberg FDR)"))
                                          # Fallback to p.adjust if qvalue fails
                                          q_values_all[valid_p_idx] <- p.adjust(valid_p_values, method = "BH")
                                          message(sprintf("      Diagnostic: Assigned %d p.adjust() values to q-value column", length(valid_p_idx)))
                                        })
                                      } else {
                                        # Use p.adjust directly for edge cases
                                        qvalue_failed <<- TRUE
                                        q_values_all[valid_p_idx] <- p.adjust(valid_p_values, method = "BH")
                                        message("      Using p.adjust() due to edge case detection")
                                        message(sprintf("      Diagnostic: Assigned %d p.adjust() values to q-value column", length(valid_p_idx)))
                                      }
                                    } else {
                                      # All p-values are invalid, set all q-values to NA
                                      q_values_all <- rep(NA_real_, nrow(de_tbl))
                                      message("      Warning: All p-values are invalid (NA, Inf, or NaN), setting q-values to NA")
                                    }
                                    
                                    # Debug: Verify assignment before mutate
                                    message(sprintf("      Diagnostic: q_values_all has %d non-NA values before assignment", sum(!is.na(q_values_all))))
                                    
                                    de_tbl <- de_tbl |>
                                      mutate({ { q_value_column } } := q_values_all)
                                    
                                    # Debug: Verify assignment after mutate
                                    assigned_col <- de_tbl[[rlang::as_name(rlang::ensym(q_value_column))]]
                                    message(sprintf("      Diagnostic: Assigned column has %d non-NA values after mutate", sum(!is.na(assigned_col))))
                                    message("      qvalue column added")
                                    
                                    message("      Adding FDR column...")
                                    # Use the same safe logic for FDR column - only compute for valid p-values
                                    fdr_values_all <- rep(NA_real_, nrow(de_tbl))
                                    if (length(valid_p_idx) > 0) {
                                      fdr_values_all[valid_p_idx] <- p.adjust(valid_p_values, method = "BH")
                                    }
                                    de_tbl <- de_tbl |>
                                      mutate({ { fdr_value_column } } := fdr_values_all)
                                    message("      FDR column added")
                                    
                                    message("      Renaming P.Value column...")
                                    de_tbl <- de_tbl |>
                                      dplyr::rename({ { p_value_column } } := P.Value)
                                    message("      P.Value column renamed")
                                    
                                    message(sprintf("   [map] Completed processing contrast: %s", contrast))
                                    if (qvalue_failed) {
                                      qvalue_failures[[contrast]] <<- TRUE
                                    }
                                    return(de_tbl)
                                    
                                  }, error = function(e) {
                                    message(sprintf("      [map] ERROR in contrast %s: %s", contrast, e$message))
                                    message(sprintf("      [map] ERROR call stack: %s", capture.output(traceback())))
                                    stop(e)
                                  })
                                }
    )
  }

  names(result_tables) <- contrast_strings
  message("--- Exiting runTestsContrasts ---")
  return(list(results = result_tables, fit.eb = t.fit, qvalue_warnings = qvalue_failures))
}


# ----------------------------------------------------------------------------
# saveDeProteinList
# ----------------------------------------------------------------------------
#'Save the list of output tables from differential expression analysis of proteins or phosphopeptides into a file and in a specific directory.
#'@param list_of_de_tables A list, each element is a table of log fold-change and q-values from differential expression analysis of proteins / phosphopeptides. Each element in the list has a name, usually the name of the pairwise comparison.
#'@param row_id Add row ID to the output table based on the name (protein or phosphopeptid ID) of each row
#'@param sort_by_column Each table in the list_of_de_tables is sorted in ascending order
#'@param results_dir The results directory to store the output file
#'@param file_suffix The file suffix string to aadd to the name of each comparison from the list_of_de_tables.
#'@export
saveDeProteinList <- function(list_of_de_tables, row_id, sort_by_column = fdr_qvalue, results_dir, file_suffix) {

  purrr::walk2(list_of_de_tables, names(list_of_de_tables),

               \(.x) {vroom::vroom_write(.x |>
                                     rownames_to_column(row_id) |>
                                     arrange({ { sort_by_column } }),
                                   path = file.path(results_dir, paste0(.y, file_suffix)))} )

}


# ----------------------------------------------------------------------------
# createDeResultsLongFormat
# ----------------------------------------------------------------------------
## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Create the de_protein_long and de_phos_long tables
#'@export
createDeResultsLongFormat <- function( lfc_qval_tbl,
                                       norm_counts_input_tbl,
                                       raw_counts_input_tbl,
                                       row_id,
                                       sample_id,
                                       group_id,
                                       group_pattern,
                                       design_matrix_norm,
                                       design_matrix_raw,
                                       expression_column = log_intensity,
                                       protein_id_table
) {

  message("   DEBUG66: createDeResultsLongFormat - Starting norm_counts processing")
  message(sprintf("      DEBUG66: norm_counts_input_tbl dims = %d x %d", nrow(norm_counts_input_tbl), ncol(norm_counts_input_tbl)))
  message(sprintf("      DEBUG66: group_pattern = %s", group_pattern))
  message(sprintf("      DEBUG66: row_id = %s", row_id))
  message(sprintf("      DEBUG66: sample_id = %s", sample_id))
  message(sprintf("      DEBUG66: group_id = %s", group_id))
  
  norm_counts <- norm_counts_input_tbl |>
    as.data.frame() |>
    rownames_to_column(row_id) |>
    pivot_longer(cols = matches(group_pattern),
                 names_to = sample_id,
                 values_to = "log2norm") |>
    left_join(design_matrix_norm, by = sample_id) |>
    group_by(!!sym(row_id), !!sym(group_id)) |>
    arrange(!!sym(row_id), !!sym(group_id), !!sym(sample_id)) |>
    mutate(replicate_number = paste0("log2norm.", row_number())) |>
    ungroup() |>
    pivot_wider(id_cols = c(!!sym(row_id), !!sym(group_id)),
                names_from = replicate_number,
                values_from = log2norm) |>
    mutate( {{group_id}} := purrr::map_chr( !!sym(group_id), as.character))
  
  message("   DEBUG66: norm_counts processing completed")


  # print(head(norm_counts))

  raw_counts <- raw_counts_input_tbl |>
    as.data.frame() |>
    rownames_to_column(row_id) |>
    pivot_longer(cols = matches(group_pattern),
                 names_to = sample_id,
                 values_to = "raw") |>
    left_join(design_matrix_raw, by = sample_id) |>
    group_by(!!sym(row_id), !!sym(group_id)) |>
    arrange(!!sym(row_id), !!sym(group_id), !!sym(sample_id)) |>
    mutate(replicate_number = paste0("raw.", row_number())) |>
    ungroup() |>
    pivot_wider(id_cols = c(!!sym(row_id), !!sym(group_id)),
                names_from = replicate_number,
                values_from = raw)  |>
    mutate( {{group_id}} := purrr::map_chr( !!sym(group_id), as.character))

  # print(head(raw_counts))

  left_join_columns <- rlang::set_names(c(row_id, group_id ),
                                        c(row_id, "left_group"))

  right_join_columns <- rlang::set_names(c(row_id, group_id ),
                                         c(row_id, "right_group"))

   # print(head(lfc_qval_tbl))

  # DEBUG66: Commented out print statements that were causing confusion
  # print( row_id)
  # print(colnames( protein_id_table)[1])

  de_proteins_long <- lfc_qval_tbl |>
    dplyr::select(-lqm, -colour, -analysis_type) |>
    dplyr::mutate( {{expression_column}} := str_replace_all({{expression_column}}, group_id, "")) |>
    separate_wider_delim( {{expression_column}}, delim = "-", names = c("left_group", "right_group"))  |>
   left_join(norm_counts, by = left_join_columns) |>
    left_join(norm_counts, by = right_join_columns,
              suffix = c(".left", ".right")) |>
    left_join(raw_counts, by = left_join_columns) |>
    left_join(raw_counts, by = right_join_columns,
              suffix = c(".left", ".right")) |>
    left_join( protein_id_table
               , by = join_by( !!sym(row_id) == !!sym( colnames( protein_id_table)[1]))) |>
    arrange( comparison, fdr_qvalue, log2FC) |>
    distinct()

  # --- NEW: Rename columns to use sample IDs if single contrast ---
  # Only perform this renaming if we have a single comparison, to ensure unique mapping
  if (length(unique(de_proteins_long$comparison)) == 1) {

    # Get the groups involved
    this_left_group <- unique(de_proteins_long$left_group)
    this_right_group <- unique(de_proteins_long$right_group)

    # Ensure we have exactly one left and one right group
    if (length(this_left_group) == 1 && length(this_right_group) == 1) {

      message(sprintf("   createDeResultsLongFormat: Renaming columns for contrast %s vs %s", this_left_group, this_right_group))

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
      map1 <- generate_rename_map(de_proteins_long, "log2norm", ".left", left_samples, this_left_group)
      map2 <- generate_rename_map(de_proteins_long, "raw", ".left", left_samples, this_left_group)
      map3 <- generate_rename_map(de_proteins_long, "log2norm", ".right", right_samples, this_right_group)
      map4 <- generate_rename_map(de_proteins_long, "raw", ".right", right_samples, this_right_group)
      
      # Combine all mappings
      all_mappings <- c(map1, map2, map3, map4)
      
      # Apply renaming in a single vectorized step
      if (length(all_mappings) > 0) {
        de_proteins_long <- de_proteins_long |> dplyr::rename(!!!all_mappings)
        message(sprintf("   createDeResultsLongFormat: Renamed %d columns", length(all_mappings)))
      }
    }
  }
  # --- END NEW ---

  # Rename group columns to numerator/denominator for clarity
  de_proteins_long <- de_proteins_long |>
    dplyr::rename(numerator = left_group, denominator = right_group)

  de_proteins_long
}


# ----------------------------------------------------------------------------
# countStatDeGenes
# ----------------------------------------------------------------------------
#' Count the number of statistically significant differentially expressed proteins (according to user-defined threshold)
#' @param lfc_thresh A numerical value specifying the log fold-change threhold (absolute value) for calling statistically significant proteins.
#' @param q_val_thresh A numerical value specifying the q-value threshold for statistically significant proteins.
#' @param log_fc_column The name of the log fold-change column (tidyverse style).
#' @param q_value_column The name of the q-value column (tidyverse style).
#' @return A table with the following columns:
#' status: The status could be Significant and Up, Significant and Down or Not significant
#' counts: The number of proteins wit this status
#'@export
countStatDeGenes <- function(data,
                             lfc_thresh = 0,
                             q_val_thresh = 0.05,
                             log_fc_column = log2FC,
                             q_value_column = fdr_qvalue) {

  # comparison <- as.data.frame(data) |>
  #   distinct(comparison) |>
  #   pull(comparison)

  selected_data <- data |>
    dplyr::mutate(status = case_when( { { q_value_column } } >= q_val_thresh ~ "Not significant"
                                      , { { log_fc_column } } >= lfc_thresh & { { q_value_column } } < q_val_thresh ~ "Significant and Up"
                                      , { { log_fc_column } } < lfc_thresh & { { q_value_column } } < q_val_thresh ~ "Significant and Down"
                                      , TRUE ~ "Not significant"))

  counts <- selected_data |>
    group_by(status) |>
    summarise(counts = n()) |>
    ungroup()

  all_possible_status <- data.frame( status = c( "Not significant", "Significant and Up", "Significant and Down"))

  results <- all_possible_status |>
    left_join( counts, by=c("status" = "status")) |>
    mutate( counts = ifelse(is.na(counts), 0, counts))

  return(results)
}


# ----------------------------------------------------------------------------
# countStatDeGenesHelper
# ----------------------------------------------------------------------------
## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#' @export
countStatDeGenesHelper <- function(de_table
                                   , description
                                   , facet_column = analysis_type
                                   , comparison_column = "comparison"
                                   , expression_column = "expression") {

  message("--- Entering countStatDeGenesHelper (DEBUG66) ---")
  message(paste("   countStatDeGenesHelper: de_table class =", class(de_table)))
  message(paste("   countStatDeGenesHelper: de_table is list =", is.list(de_table)))
  message(paste("   countStatDeGenesHelper: de_table length =", length(de_table)))
  message("   countStatDeGenesHelper: de_table names:")
  print(names(de_table))

  de_table_updated <- purrr::map(de_table, \(x){  countStatDeGenes(x,
                                                                   lfc_thresh = 0,
                                                                   q_val_thresh = 0.05,
                                                                   log_fc_column = logFC,
                                                                   q_value_column = fdr_qvalue)})

  message("   countStatDeGenesHelper: de_table_updated created")
  message(paste("   countStatDeGenesHelper: de_table_updated length =", length(de_table_updated)))
  message("   countStatDeGenesHelper: de_table_updated names:")
  print(names(de_table_updated))

  list_of_tables <- purrr::map2(de_table_updated
                                ,names(de_table_updated)
                                ,\(.x, .y){ 
                                  message(paste("      [map2] Processing element with name:", .y))
                                  .x |> mutate(!!sym(comparison_column) := .y) 
                                })

  message("   countStatDeGenesHelper: list_of_tables created")
  message("   countStatDeGenesHelper: About to bind_rows...")
  
  bound_tables <- list_of_tables |> bind_rows()
  message(paste("   countStatDeGenesHelper: bound_tables dims =", nrow(bound_tables), "x", ncol(bound_tables)))
  
  faceted_tables <- bound_tables |> mutate({ { facet_column } } := description)
  message("   countStatDeGenesHelper: facet column added")
  message("   countStatDeGenesHelper: Checking comparison column values:")
  print(unique(faceted_tables[[comparison_column]]))
  
  message(paste("   countStatDeGenesHelper: About to separate_wider_delim on column:", comparison_column))
  message(paste("   countStatDeGenesHelper: Looking for delimiter: ="))
  
  # Check if any values contain "="
  has_delimiter <- any(grepl("=", faceted_tables[[comparison_column]]))
  message(paste("   countStatDeGenesHelper: Any values contain '=' ?", has_delimiter))
  
  if (!has_delimiter) {
    message("   countStatDeGenesHelper: WARNING - No '=' found in comparison column values!")
    message("   countStatDeGenesHelper: This will cause separate_wider_delim to fail!")
    message("   countStatDeGenesHelper: Comparison column values are:")
    print(faceted_tables[[comparison_column]])
    stop("countStatDeGenesHelper: comparison column values do not contain '=' delimiter. Check list element naming in calling function.")
  }

  merged_tables <- faceted_tables |>
    separate_wider_delim( !!sym(comparison_column ),
                          delim = "=",
                          names = c(comparison_column,
                                    expression_column))

  message("--- Exiting countStatDeGenesHelper (DEBUG66) ---")
  merged_tables
}


# ----------------------------------------------------------------------------
# printCountDeGenesTable
# ----------------------------------------------------------------------------
## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#' Format results table for use in volcano plots, counting number of significant proteins, p-values distribution histogram.
#' @param list_of_de_tables A list with each element being a results table with log fold-change and q-value per protein.
#' @param list_of_descriptions  A list of strings describing the parameters used to generate the result table.
#' @param formula_string The formula string used in the facet_grid command for the ggplot scatter plot.
#' @param facet_column The name of the column describing the type of analysis or parameters used to generate the result table (tidyverse style). This is related to the \code{list_of_descriptions} parameter above.
#' @param comparison_column The name of the column describing the contrasts or comparison between groups (tidyverse style).
#' @param expression_column The name of the column that will contain the formula expressions of the contrasts.
#'@export
printCountDeGenesTable <- function(  list_of_de_tables
                                     , list_of_descriptions
                                     , formula_string = "analysis_type ~ comparison"
                                     , facet_column = analysis_type
                                     , comparison_column = "comparison"
                                     , expression_column = "expression") {




  num_significant_de_genes_all <- purrr::map2(list_of_de_tables,
                                              list_of_descriptions,
                                              function(a, b) { countStatDeGenesHelper( de_table = a
                                                                                       , description = b
                                                                                       , facet_column = {{facet_column}}
                                                                                       , comparison_column = comparison_column
                                                                                       , expression_column = expression_column) }) |>
    bind_rows()

  num_sig_de_genes_barplot <- num_significant_de_genes_all |>
    dplyr::filter(status != "Not significant") |>
    ggplot(aes(x = status, y = counts)) +
    geom_bar(stat = "identity") +
    geom_text(stat = 'identity', aes(label = counts), vjust = -0.5) +
    theme(axis.text.x = element_text(angle = 90))

  # print(head(num_sig_de_genes_barplot))

  if (!is.na(formula_string)) {

    num_sig_de_genes_barplot <- num_sig_de_genes_barplot +
      facet_grid(as.formula(formula_string))
  }


  return(list(plot = num_sig_de_genes_barplot, table = num_significant_de_genes_all))
}


# ----------------------------------------------------------------------------
# getSignificantData
# ----------------------------------------------------------------------------
## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#' Format results table for use in volcano plots, counting number of significant proteins, p-values distribution histogram.
#' @param list_of_de_tables A list with each element being a results table with log fold-change and q-value per protein.
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
getSignificantData <- function( list_of_de_tables
                                , list_of_descriptions
                                , row_id = uniprot_acc
                                , p_value_column = raw_pvalue
                                , q_value_column = fdr_qvalue
                                , fdr_value_column = fdr_value_bh_adjustment
                                , log_q_value_column = lqm
                                , log_fc_column = logFC
                               , comparison_column = "comparison"
                               , expression_column = "log_intensity"
                               , facet_column = analysis_type
                               , q_val_thresh = 0.05) {

  message("--- Entering getSignificantData (DEBUG66) ---")
  message(paste("   getSignificantData: list_of_de_tables class =", class(list_of_de_tables)))
  message(paste("   getSignificantData: list_of_de_tables length =", length(list_of_de_tables)))
  message("   getSignificantData: list_of_de_tables structure:")
  str(list_of_de_tables, max.level = 2)

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
          .x <- .x |> mutate({ { comparison_column } } := .y)
        }
        .x
        }) |>
      bind_rows() |>
      mutate({ { facet_column } } := description)

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
            { { comparison_column } },
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
  
  logfc_tbl_all <- purrr::map2(list_of_de_tables, list_of_descriptions,
                               function(a, b) { get_row_binded_table(de_table_list = a, description = b) }) |>
    bind_rows()

  message("   getSignificantData: All tables bound")
  message(paste("   getSignificantData: logfc_tbl_all dims =", nrow(logfc_tbl_all), "x", ncol(logfc_tbl_all)))

  selected_data <- logfc_tbl_all |>
    mutate({ { log_q_value_column } } := -log10(fdr_qvalue)) |>
    dplyr::select({ { row_id } }, { { log_q_value_column } }, { { q_value_column } }, { { p_value_column } }, { { log_fc_column } },
                  { { comparison_column } }, { { expression_column } },
                  { { facet_column } }) |>
    dplyr::mutate(colour = case_when(abs({ { log_fc_column } }) >= 1 & { { q_value_column } } >= q_val_thresh ~ "orange",
                                     abs({ { log_fc_column } }) >= 1 & { { q_value_column } } < q_val_thresh ~ "purple",
                                     abs({ { log_fc_column } }) < 1 & { { q_value_column } } < q_val_thresh ~ "blue",
                                     TRUE ~ "black")) |>
    dplyr::mutate(colour = factor(colour, levels = c("black", "orange", "blue", "purple")))

  message("--- Exiting getSignificantData (DEBUG66) ---")
  selected_data

}


# ----------------------------------------------------------------------------
# getTypeOfGrouping
# ----------------------------------------------------------------------------
#'Assign experimental group list
#'@param design_matrix A data frame representing the design matrix.
#'@param group_id A string representing the name of the group ID column used in the design matrix.
#'@param sample_id A string representing the name of the sample ID column used in the design matrix.
#'@return A list where each element name is the name of a treatment group and each element is a vector containing the sample IDs within the treatment group.
#'@export
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
#'@export
extractResults <- function(results_list) {

  extracted <- purrr::map(results_list, \(x){ x$results })

  names(extracted) <- names(results_list)

  return(extracted)
}


# ----------------------------------------------------------------------------
# writeInteractiveVolcanoPlotProteomics
# ----------------------------------------------------------------------------
## Create proteomics interactive volcano plot
#' @export
# de_analysis_results_list$contrasts_results$fit.eb
# No full stops in the nme of columns of interactive table in glimma plot. It won't display column with full stop in the column name.
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
                                                   , gene_names_column = "gene_names"
                                                   , display_columns = c( "best_uniprot_acc" )) {


  volcano_plot_tab <- de_proteins_long  |>
    dplyr::mutate(best_uniprot_acc = str_split(!!sym(args_row_id), ":" ) |> purrr::map_chr(1)  ) |>
    left_join(uniprot_tbl, by = join_by( best_uniprot_acc == !!sym( uniprot_id_column) ) ) |>
    dplyr::rename( UNIPROT_GENENAME = gene_names_column ) |>
    mutate( UNIPROT_GENENAME = purrr::map_chr( UNIPROT_GENENAME, \(x){
      tryCatch({
        if(is.na(x) || is.null(x) || x == "") {
          ""
        } else {
          split_result <- str_split(x, " |:")[[1]]
          if(length(split_result) > 0) split_result[1] else ""
        }
      }, error = function(e) "")
    })) |>
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


# ----------------------------------------------------------------------------
# writeInteractiveVolcanoPlotProteomicsWidget
# ----------------------------------------------------------------------------
#' @export
# de_analysis_results_list$contrasts_results$fit.eb
# No full stops in the nme of columns of interactive table in glimma plot. It won't display column with full stop in the column name.
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
                                                         , gene_names_column = "gene_names"
                                                         , display_columns = c( "best_uniprot_acc" )) {


  volcano_plot_tab <- de_proteins_long  |>
    left_join(uniprot_tbl, by = join_by( !!sym(args_row_id) == !!sym( uniprot_id_column) ) ) |>
    dplyr::rename( UNIPROT_GENENAME = gene_names_column ) |>
    mutate( UNIPROT_GENENAME = purrr::map_chr( UNIPROT_GENENAME, \(x){
      tryCatch({
        if(is.na(x) || is.null(x) || x == "") {
          ""
        } else {
          split_result <- str_split(x, " ")[[1]]
          if(length(split_result) > 0) split_result[1] else ""
        }
      }, error = function(e) "")
    })) |>
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


# ----------------------------------------------------------------------------
# writeInteractiveVolcanoPlotProteomicsMain
# ----------------------------------------------------------------------------
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
  gene_names_column <- checkParamsObjectFunctionSimplify(theObject, "gene_names_column", "gene_names")
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

  # Use helper to extract counts table instead of direct slot access
  counts_mat <- getCountsTable(de_analysis_results_list$theObject) |>
    column_to_rownames(args_row_id) |>
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


# ----------------------------------------------------------------------------
# getDataMatrix
# ----------------------------------------------------------------------------
# Helper function to get data matrix
getDataMatrix <- function(obj) {
    if (inherits(obj, "MetaboliteAssayData")) {
        message(sprintf("   Getting data matrix for object of class: %s", class(obj)[1]))
        message(sprintf("   Processing MetaboliteAssayData"))
        message(sprintf("   Metabolite data dimensions: %d rows, %d cols",
                      nrow(obj@metabolite_data), ncol(obj@metabolite_data)))
        matrix_data <- as.matrix(obj@metabolite_data[, -1]) # Exclude Name column
        colnames(matrix_data) <- colnames(obj@metabolite_data)[-1]
        rownames(matrix_data) <- obj@metabolite_data$Name
        message(sprintf("   Created matrix with dimensions: %d rows, %d cols",
                      nrow(matrix_data), ncol(matrix_data)))
        matrix_data
    } else if (inherits(obj, "ProteinQuantitativeData")) {
        message(sprintf("   Processing ProteinQuantitativeData"))
        message(sprintf("   Protein quant table dimensions: %d rows, %d cols",
                      nrow(obj@protein_quant_table), ncol(obj@protein_quant_table)))
        result <- as.matrix(column_to_rownames(obj@protein_quant_table, obj@protein_id_column))
        message(sprintf("   Created matrix with dimensions: %d rows, %d cols",
                      nrow(result), ncol(result)))
        result
    } else {
        message(sprintf("   ERROR: Unsupported object type: %s", class(obj)[1]))
        stop("Unsupported object type")
    }
}


# ----------------------------------------------------------------------------
# differentialExpressionAnalysis
# ----------------------------------------------------------------------------
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

  # IMMEDIATE ERROR CATCH - Check if we even get here
  message("*** ENTERING differentialExpressionAnalysis METHOD ***")
  message(sprintf("*** METHOD SIGNATURE MATCHED: ProteinQuantitativeData ***"))
  
  # Try to catch the index error immediately
  tryCatch({
    message("*** Testing parameter access ***")
    if (!is.null(contrasts_tbl)) {
      test <- contrasts_tbl[[1]]
      message("*** Parameter access successful ***")
    }
  }, error = function(e) {
    message(sprintf("*** IMMEDIATE ERROR: %s ***", e$message))
    message(sprintf("*** ERROR CLASS: %s ***", class(e)))
    stop(e)
  })

  message("--- Entering differentialExpressionAnalysis ---")
  message(sprintf("   differentialExpressionAnalysis: theObject class = %s", class(theObject)))
  message(sprintf("   differentialExpressionAnalysis: contrasts_tbl provided = %s", !is.null(contrasts_tbl)))
  if(!is.null(contrasts_tbl)) {
    message(sprintf("   differentialExpressionAnalysis: contrasts_tbl dims = %d x %d", nrow(contrasts_tbl), ncol(contrasts_tbl)))
    message(sprintf("   differentialExpressionAnalysis: contrasts_tbl content = %s", paste(contrasts_tbl[[1]], collapse=", ")))
  }

  # Wrap the helper function call in tryCatch to get better error info
  message("   differentialExpressionAnalysis: About to call differentialExpressionAnalysisHelper...")
  
  results_list <- tryCatch({
    differentialExpressionAnalysisHelper(  theObject
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
  }, error = function(e) {
    # CRITICAL FIX: Use paste() for logger calls in error handlers to avoid interpolation bug
    message(paste("   differentialExpressionAnalysis ERROR in helper function:", e$message))
    message(paste("   differentialExpressionAnalysis ERROR call stack:", capture.output(traceback())))
    stop(e)
  })

  message("   differentialExpressionAnalysis: Helper function completed successfully!")
  message("--- Exiting differentialExpressionAnalysis ---")
  return(results_list)

})


# ----------------------------------------------------------------------------
# differentialExpressionAnalysisHelper
# ----------------------------------------------------------------------------
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
  
  # IMMEDIATE PARAMETER VALIDATION TO CATCH INDEX ERROR
  message("   PARAMETER VALIDATION Step: Checking all input parameters...")
  message(sprintf("      theObject is NULL: %s", is.null(theObject)))
  message(sprintf("      contrasts_tbl is NULL: %s", is.null(contrasts_tbl)))
  if (!is.null(contrasts_tbl)) {
    message(sprintf("      contrasts_tbl class: %s", class(contrasts_tbl)))
    message(sprintf("      contrasts_tbl type: %s", typeof(contrasts_tbl)))
    message("      contrasts_tbl print:")
    print(contrasts_tbl)
    message("      Trying to access contrasts_tbl[[1]]...")
    tryCatch({
      first_col <- contrasts_tbl[[1]]
      message(sprintf("      SUCCESS: contrasts_tbl[[1]] class: %s", class(first_col)))
      message(sprintf("      SUCCESS: contrasts_tbl[[1]] length: %d", length(first_col)))
      message(sprintf("      SUCCESS: contrasts_tbl[[1]] content: %s", paste(first_col, collapse=", ")))
    }, error = function(e) {
      message(sprintf("      ERROR accessing contrasts_tbl[[1]]: %s", e$message))
    })
  }
  
  message("   DEBUG66: Initial parameter inspection")
  message(sprintf("      DEBUG66: theObject class = %s", class(theObject)))
  message(sprintf("      DEBUG66: contrasts_tbl param is.null = %s", is.null(contrasts_tbl)))
  if (!is.null(contrasts_tbl)) {
    message(sprintf("      DEBUG66: contrasts_tbl param class = %s", class(contrasts_tbl)))
    message("      DEBUG66: contrasts_tbl param structure:")
    str(contrasts_tbl)
  }

  message("   differentialExpressionAnalysisHelper Step: Extracting parameters...")
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
  
  message("   DEBUG66: After parameter extraction")
  message(sprintf("      DEBUG66: contrasts_tbl class = %s", class(contrasts_tbl)))
  message("      DEBUG66: contrasts_tbl structure:")
  str(contrasts_tbl)

  message(sprintf("   differentialExpressionAnalysisHelper: formula_string = %s", formula_string))
  message(sprintf("   differentialExpressionAnalysisHelper: group_id = %s", group_id))
  message(sprintf("   differentialExpressionAnalysisHelper: de_q_val_thresh = %f", de_q_val_thresh))

  message("   differentialExpressionAnalysisHelper Step: Processing group names...")
  # Handle group names that start with numbers (same pattern as original wrapper)
  design_matrix <- theObject@design_matrix
  group_col <- design_matrix[[group_id]]
  message(sprintf("   differentialExpressionAnalysisHelper: group_col length = %d", length(group_col)))
  message(sprintf("   differentialExpressionAnalysisHelper: group_col content = %s", paste(head(group_col, 10), collapse=", ")))
  
  # Check if any group names start with numbers and create mapping
  starts_with_number <- grepl("^[0-9]", group_col)
  message(sprintf("   differentialExpressionAnalysisHelper: any start with number? %s", any(starts_with_number)))
  
  if(any(starts_with_number)) {
    message("   differentialExpressionAnalysisHelper Step: Fixing group names that start with numbers...")
    
    original_groups <- unique(group_col)
    message(sprintf("   differentialExpressionAnalysisHelper: original_groups = %s", paste(original_groups, collapse=", ")))
    message(sprintf("   differentialExpressionAnalysisHelper: About to process %d original groups with purrr::map_chr", length(original_groups)))
    
    tryCatch({
      message("      DEBUG66: Before safe_groups purrr::map_chr")
      message(sprintf("      DEBUG66: original_groups class = %s", class(original_groups)))
      message(sprintf("      DEBUG66: original_groups length = %d", length(original_groups)))
      message("      DEBUG66: original_groups content:")
      print(original_groups)
      
      safe_groups <- purrr::map_chr(original_groups, \(x) {
        message(sprintf("         DEBUG66: Processing group item: '%s' (class: %s)", x, class(x)))
        result <- if(grepl("^[0-9]", x)) paste0("grp_", x) else x
        message(sprintf("         DEBUG66: Result for '%s' -> '%s'", x, result))
        result
      })
      message("   differentialExpressionAnalysisHelper: safe_groups processing SUCCESS")
      message("      DEBUG66: safe_groups result:")
      print(safe_groups)
    }, error = function(e) {
      message(sprintf("   differentialExpressionAnalysisHelper ERROR in safe_groups purrr::map_chr: %s", e$message))
      message("      DEBUG66: Error details:")
      message(sprintf("      DEBUG66: Error class: %s", class(e)))
      print(str(e))
      stop(e)
    })
    
    group_mapping <- setNames(original_groups, safe_groups)
    message("   differentialExpressionAnalysisHelper: group_mapping created")
    
    # Update design matrix with safe names
    message("   differentialExpressionAnalysisHelper: About to update design matrix with purrr::map_chr")
    tryCatch({
      design_matrix[[group_id]] <- purrr::map_chr(group_col, \(x) {
        if(grepl("^[0-9]", x)) paste0("grp_", x) else x
      })
      message("   differentialExpressionAnalysisHelper: design matrix update SUCCESS")
    }, error = function(e) {
      message(sprintf("   differentialExpressionAnalysisHelper ERROR in design matrix purrr::map_chr: %s", e$message))
      stop(e)
    })
    
    # Update contrasts table if it exists
    if(!is.null(contrasts_tbl)) {
      message("   differentialExpressionAnalysisHelper: About to update contrasts table with purrr::map_chr")
      message(sprintf("   differentialExpressionAnalysisHelper: contrasts_tbl[[1]] = %s", paste(contrasts_tbl[[1]], collapse=", ")))
      
      message("      DEBUG66: Contrasts table inspection before purrr::map_chr")
      message(sprintf("      DEBUG66: contrasts_tbl class = %s", class(contrasts_tbl)))
      message(sprintf("      DEBUG66: contrasts_tbl dim = %d x %d", nrow(contrasts_tbl), ncol(contrasts_tbl)))
      message("      DEBUG66: contrasts_tbl structure:")
      str(contrasts_tbl)
      message(sprintf("      DEBUG66: contrasts_tbl[[1]] class = %s", class(contrasts_tbl[[1]])))
      message(sprintf("      DEBUG66: contrasts_tbl[[1]] length = %d", length(contrasts_tbl[[1]])))
      message("      DEBUG66: group_mapping content:")
      print(group_mapping)
      
      tryCatch({
        contrasts_tbl[[1]] <- purrr::map_chr(contrasts_tbl[[1]], \(x) {
          message(sprintf("         DEBUG66: Processing contrast: '%s' (class: %s)", x, class(x)))
          message(sprintf("         DEBUG66: group_mapping names: %s", paste(names(group_mapping), collapse=", ")))
          result <- x
          for(orig in names(group_mapping)) {
            message(sprintf("            DEBUG66: Checking gsub: '%s' -> '%s'", group_mapping[orig], orig))
            result <- gsub(group_mapping[orig], orig, result, fixed = TRUE)
          }
          message(sprintf("         DEBUG66: Final result: '%s' -> '%s'", x, result))
          result
        })
        message("   differentialExpressionAnalysisHelper: contrasts table update SUCCESS")
        message("      DEBUG66: Updated contrasts_tbl[[1]]:")
        print(contrasts_tbl[[1]])
      }, error = function(e) {
        message(sprintf("   differentialExpressionAnalysisHelper ERROR in contrasts table purrr::map_chr: %s", e$message))
        message("      DEBUG66: Error details:")
        message(sprintf("      DEBUG66: Error class: %s", class(e)))
        print(str(e))
        stop(e)
      })
    }
    
    theObject@design_matrix <- design_matrix
  }

  message("   differentialExpressionAnalysisHelper Step: Updating S4 object parameters...")
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
  message("   DEBUG66: About to prepare design matrix rownames")
  message(paste("      DEBUG66: theObject@sample_id =", theObject@sample_id))
  message(paste("      DEBUG66: design_matrix dims before =", nrow(theObject@design_matrix), "x", ncol(theObject@design_matrix)))
  message("      DEBUG66: design_matrix column names:")
  print(colnames(theObject@design_matrix))
  
  rownames( theObject@design_matrix ) <- theObject@design_matrix |> dplyr::pull( one_of(theObject@sample_id ))
  
  message("      DEBUG66: design_matrix rownames set successfully")

  # Run the core limma analysis using existing function
  protein_quant_matrix <- as.matrix(column_to_rownames(theObject@protein_quant_table, theObject@protein_id_column))
  contrast_strings_to_use <- contrasts_tbl[, 1]
  message(paste("   differentialExpressionAnalysisHelper: About to call runTestsContrasts with", length(contrast_strings_to_use), "contrasts"))
  message(paste("   differentialExpressionAnalysisHelper: protein_quant_matrix dims =", nrow(protein_quant_matrix), "x", ncol(protein_quant_matrix)))
  
  contrasts_results <- runTestsContrasts(protein_quant_matrix,
                                         contrast_strings = contrast_strings_to_use,
                                         design_matrix = theObject@design_matrix,
                                         formula_string = formula_string,
                                         weights = NA,
                                         treat_lfc_cutoff = as.double(treat_lfc_cutoff),
                                         eBayes_trend = as.logical(eBayes_trend),
                                         eBayes_robust = as.logical(eBayes_robust))

  message("   differentialExpressionAnalysisHelper Step: runTestsContrasts completed successfully!")

  # The result from runTestsContrasts is a LIST with a 'results'
  # element, which ITSELF is a list of tables. For a single contrast run, we 
  # need to extract the first table from the nested list.
  
  # Check if the expected structure exists
  if (is.null(contrasts_results) || is.null(contrasts_results$results) || length(contrasts_results$results) == 0) {
      stop("Error: DE analysis function did not return results in the expected format.")
  }
  
  # Extract the results table (it's the first element in the nested list)
  de_results_table <- contrasts_results$results[[1]]
  
  # Get the name of the contrast from the list element name
  contrast_name <- names(contrasts_results$results)[1]
  
  # Ensure it's a data frame before proceeding
  if (!is.data.frame(de_results_table)) {
    stop("Error: DE analysis results are not in the expected data frame format.")
  }

  # CRITICAL FIX 3.0: The 'topTreat' table needs a 'comparison' column added to it,
  # containing the name of the contrast. It also needs the protein IDs from rownames.
  de_results_table <- de_results_table |>
    tibble::rownames_to_column(var = args_row_id) |>
    dplyr::mutate(comparison = contrast_name)

  # Map back to original group names in results if needed
  if(exists("group_mapping")) {
    contrasts_results_table <- de_results_table |>
      dplyr::mutate(comparison = purrr::map_chr(comparison, \(x) {
        result <- x
        for(safe_name in names(group_mapping)) {
          result <- gsub(safe_name, group_mapping[safe_name], result, fixed = TRUE)
        }
        result
      }))
  } else {
    contrasts_results_table <- de_results_table
  }

  return_list$contrasts_results <- contrasts_results
  return_list$contrasts_results_table <- contrasts_results_table
  
  # Extract qvalue warnings if present
  if (!is.null(contrasts_results$qvalue_warnings) && length(contrasts_results$qvalue_warnings) > 0) {
    return_list$qvalue_warnings <- contrasts_results$qvalue_warnings
    message(sprintf("   differentialExpressionAnalysisHelper Step: qvalue() failed for %d contrast(s)", length(contrasts_results$qvalue_warnings)))
  }

  message("   differentialExpressionAnalysisHelper Step: Preparing data for visualization...")
  message(paste("   DEBUG66: contrast_name for list naming =", contrast_name))

  # Prepare data for volcano plots
  # CRITICAL FIX: getSignificantData expects list_of_de_tables to be a list where each element 
  # is itself a list of data.frames. Since we're processing one contrast at a time, we need to 
  # wrap the single data.frame in an additional list layer.
  # CRITICAL FIX 2: The list element name MUST contain "=" delimiter because countStatDeGenesHelper
  # expects to split on "=" to separate comparison name from expression. Use contrast_name which
  # contains the full format like "H4_vs_WT=groupH4-groupWT"
  nested_list <- list(contrasts_results_table)
  names(nested_list) <- contrast_name  # Use actual contrast name with "=" delimiter
  significant_rows <- getSignificantData(list_of_de_tables = list(nested_list),
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
  # CRITICAL FIX: Same as above - wrap in additional list layer and use contrast_name
  nested_list_for_count <- list(contrasts_results_table)
  names(nested_list_for_count) <- contrast_name  # Use actual contrast name with "=" delimiter
  num_sig_de_molecules <- printCountDeGenesTable(list_of_de_tables = list(nested_list_for_count),
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


# ----------------------------------------------------------------------------
# outputDeResultsAllContrasts
# ----------------------------------------------------------------------------
#' @export
setMethod(f = "outputDeResultsAllContrasts",
          signature = "ProteinQuantitativeData",
          definition = function(theObject,
                               de_results_list_all_contrasts = NULL,
                               uniprot_tbl = NULL,
                               de_output_dir = NULL,
                               publication_graphs_dir = NULL,
                               file_prefix = "de_proteins",
                               args_row_id = NULL,
                               gene_names_column = "gene_names",
                               uniprot_id_column = "Entry") {

  message("--- Entering outputDeResultsAllContrasts ---")
  message(sprintf("   outputDeResultsAllContrasts: de_output_dir = %s", de_output_dir))
  message(sprintf("   outputDeResultsAllContrasts: file_prefix = %s", file_prefix))
  
  # Extract parameters from S4 object with fallbacks
  uniprot_tbl <- checkParamsObjectFunctionSimplify(theObject, "uniprot_tbl", uniprot_tbl)
  de_output_dir <- checkParamsObjectFunctionSimplify(theObject, "de_output_dir", de_output_dir)
  publication_graphs_dir <- checkParamsObjectFunctionSimplify(theObject, "publication_graphs_dir", publication_graphs_dir)
  args_row_id <- checkParamsObjectFunctionSimplify(theObject, "args_row_id", args_row_id)
  gene_names_column <- checkParamsObjectFunctionSimplify(theObject, "gene_names_column", gene_names_column)
  uniprot_id_column <- checkParamsObjectFunctionSimplify(theObject, "uniprot_id_column", uniprot_id_column)
  
  # CRITICAL FIX: Normalize paths for Windows (fixes C:// double slash issue)
  # Only normalize if the path is not NULL and doesn't already exist
  if (!is.null(de_output_dir)) {
    # Use normalizePath with mustWork=FALSE to handle non-existent dirs
    de_output_dir <- gsub("//+", "/", de_output_dir)  # Remove double slashes first
    de_output_dir <- normalizePath(de_output_dir, winslash = "/", mustWork = FALSE)
    message(sprintf("   outputDeResultsAllContrasts: Normalized de_output_dir = %s", de_output_dir))
  }
  
  if (!is.null(publication_graphs_dir)) {
    publication_graphs_dir <- gsub("//+", "/", publication_graphs_dir)  # Remove double slashes first
    publication_graphs_dir <- normalizePath(publication_graphs_dir, winslash = "/", mustWork = FALSE)
    message(sprintf("   outputDeResultsAllContrasts: Normalized publication_graphs_dir = %s", publication_graphs_dir))
  }
  
  # Ensure output directory exists (CRITICAL FIX!)
  if (!dir.exists(de_output_dir)) {
    dir.create(de_output_dir, recursive = TRUE, showWarnings = FALSE)
    message(sprintf("   outputDeResultsAllContrasts: Created output directory: %s", de_output_dir))
  }
  
  message(sprintf("   outputDeResultsAllContrasts: Processing %d contrasts", length(de_results_list_all_contrasts)))
  
  # Write results for each contrast with UNIQUE filenames
  for (contrast_name in names(de_results_list_all_contrasts)) {
    message(sprintf("   outputDeResultsAllContrasts: Writing files for contrast: %s", contrast_name))
    
    contrast_result <- de_results_list_all_contrasts[[contrast_name]]
    
    if (!is.null(contrast_result$de_proteins_long)) {
      # Clean contrast name for safe filename
      safe_contrast_name <- gsub("[^A-Za-z0-9_-]", "_", contrast_name)
      
             # Create annotated version (FIXED: proper conditional logic)
       de_proteins_long_annot <- contrast_result$de_proteins_long |>
         mutate(uniprot_acc_cleaned = str_split(!!sym(args_row_id), "-") |>
                  purrr::map_chr(1))
       
       # Add UniProt annotations if available
       if (!is.null(uniprot_tbl) && nrow(uniprot_tbl) > 0) {
         de_proteins_long_annot <- de_proteins_long_annot |>
           left_join(uniprot_tbl, by = join_by(uniprot_acc_cleaned == !!sym(uniprot_id_column))) |>
           dplyr::select(-uniprot_acc_cleaned)
       } else {
         de_proteins_long_annot <- de_proteins_long_annot |>
           dplyr::select(-uniprot_acc_cleaned)
       }
       
       # Add gene_name column with proper conditional logic
       if ("gene_names" %in% names(de_proteins_long_annot)) {
         de_proteins_long_annot <- de_proteins_long_annot |>
           mutate(gene_name = purrr::map_chr(gene_names, \(x){
             tryCatch({
               if(is.na(x) || is.null(x) || x == "") {
                 ""
               } else {
                 split_result <- str_split(x, " |:")[[1]]
                 if(length(split_result) > 0) split_result[1] else ""
               }
             }, error = function(e) "")
           }))
       } else if (gene_names_column %in% names(de_proteins_long_annot)) {
         de_proteins_long_annot <- de_proteins_long_annot |>
           mutate(gene_name = purrr::map_chr(.data[[gene_names_column]], \(x){
             tryCatch({
               if(is.na(x) || is.null(x) || x == "") {
                 ""
               } else {
                 split_result <- str_split(x, " |:")[[1]]
                 if(length(split_result) > 0) split_result[1] else ""
               }
             }, error = function(e) "")
           }))
       } else {
         de_proteins_long_annot <- de_proteins_long_annot |>
           mutate(gene_name = "")
       }
       
       # Relocate gene_name column
       de_proteins_long_annot <- de_proteins_long_annot |>
         relocate(gene_name, .after = !!sym(args_row_id))
      
      # CRITICAL FIX: Use contrast-specific filename!
      contrast_filename <- paste0(file_prefix, "_", safe_contrast_name, "_long_annot.tsv")
      output_path <- file.path(de_output_dir, contrast_filename)
      
      message(sprintf("   outputDeResultsAllContrasts: Writing %s", output_path))
      
      # Write the file
      vroom::vroom_write(de_proteins_long_annot, output_path)
      
      # Verify file was written
      if (file.exists(output_path)) {
        file_size <- file.size(output_path)
        message(sprintf("   outputDeResultsAllContrasts: SUCCESS - File written: %s (%d bytes)", 
                       contrast_filename, file_size))
      } else {
        message(sprintf("   outputDeResultsAllContrasts: ERROR - File NOT written: %s", output_path))
      }
      
      # Also write Excel version
      excel_filename <- paste0(file_prefix, "_", safe_contrast_name, "_long_annot.xlsx")
      excel_path <- file.path(de_output_dir, excel_filename)
      writexl::write_xlsx(de_proteins_long_annot, excel_path)
      
      message(sprintf("   outputDeResultsAllContrasts: Also wrote Excel: %s", excel_filename))
      
      # ✅ NEW: Generate volcano plot for this contrast
      message(sprintf("   outputDeResultsAllContrasts: Generating volcano plot for contrast: %s", contrast_name))
      
      # Extract parameters - try direct access first, then use checkParams
      de_q_val_thresh <- if (!is.null(theObject@args$outputDeResultsAllContrasts$de_q_val_thresh)) {
        theObject@args$outputDeResultsAllContrasts$de_q_val_thresh
      } else {
        checkParamsObjectFunctionSimplify(theObject, "de_q_val_thresh", 0.05)
      }
      
      fdr_column <- if (!is.null(theObject@args$outputDeResultsAllContrasts$fdr_column)) {
        theObject@args$outputDeResultsAllContrasts$fdr_column
      } else {
        checkParamsObjectFunctionSimplify(theObject, "fdr_column", "fdr_qvalue")
      }
      
      log2fc_column <- if (!is.null(theObject@args$outputDeResultsAllContrasts$log2fc_column)) {
        theObject@args$outputDeResultsAllContrasts$log2fc_column
      } else {
        checkParamsObjectFunctionSimplify(theObject, "log2fc_column", "log2FC")
      }
      
      # Prepare data for volcano plot (same logic as in differentialExpressionAnalysisHelper)
      volcano_data <- contrast_result$de_proteins_long |>
        mutate(lqm = -log10(!!sym(fdr_column))) |>
        dplyr::mutate(label = case_when(!!sym(fdr_column) < de_q_val_thresh ~ "Significant",
                                       TRUE ~ "Not sig.")) |>
        dplyr::mutate(colour = case_when(!!sym(fdr_column) < de_q_val_thresh ~ "purple",
                                        TRUE ~ "black")) |>
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
        message(sprintf("   outputDeResultsAllContrasts: Created volcano plots directory: %s", volcano_dir))
      }
      
      # Save volcano plot with contrast-specific filename
      volcano_png_path <- file.path(volcano_dir, paste0(safe_contrast_name, ".png"))
      volcano_pdf_path <- file.path(volcano_dir, paste0(safe_contrast_name, ".pdf"))
      
      # Save as PNG
      tryCatch({
        ggplot2::ggsave(volcano_png_path, volcano_plot, width = 7, height = 7, dpi = 300)
        message(sprintf("   outputDeResultsAllContrasts: Saved volcano plot PNG: %s", basename(volcano_png_path)))
      }, error = function(e) {
        message(sprintf("   outputDeResultsAllContrasts: ERROR saving PNG: %s", e$message))
      })
      
      # Save as PDF
      tryCatch({
        ggplot2::ggsave(volcano_pdf_path, volcano_plot, width = 7, height = 7)
        message(sprintf("   outputDeResultsAllContrasts: Saved volcano plot PDF: %s", basename(volcano_pdf_path)))
      }, error = function(e) {
        message(sprintf("   outputDeResultsAllContrasts: ERROR saving PDF: %s", e$message))
      })
    }
  }
  
  # ✅ NEW: Generate combined multi-page PDF with all volcano plots
  message("   outputDeResultsAllContrasts: Creating combined volcano plots PDF...")
  
  volcano_dir <- file.path(publication_graphs_dir, "Volcano_Plots")
  combined_pdf_path <- file.path(volcano_dir, "all_volcano_plots_combined.pdf")
  
  tryCatch({
    # Re-extract parameters for combined PDF generation (in case they were modified)
    de_q_val_thresh <- if (!is.null(theObject@args$outputDeResultsAllContrasts$de_q_val_thresh)) {
      theObject@args$outputDeResultsAllContrasts$de_q_val_thresh
    } else {
      0.05
    }
    
    fdr_column <- if (!is.null(theObject@args$outputDeResultsAllContrasts$fdr_column)) {
      theObject@args$outputDeResultsAllContrasts$fdr_column
    } else {
      "fdr_qvalue"
    }
    
    log2fc_column <- if (!is.null(theObject@args$outputDeResultsAllContrasts$log2fc_column)) {
      theObject@args$outputDeResultsAllContrasts$log2fc_column
    } else {
      "log2FC"
    }
    
    # Collect all volcano plots
    all_volcano_plots <- list()
    
    for (contrast_name in names(de_results_list_all_contrasts)) {
      contrast_result <- de_results_list_all_contrasts[[contrast_name]]
      
      if (!is.null(contrast_result$de_proteins_long)) {
        # Recreate the volcano plot
        volcano_data <- contrast_result$de_proteins_long |>
          mutate(lqm = -log10(!!sym(fdr_column))) |>
          dplyr::mutate(label = case_when(!!sym(fdr_column) < de_q_val_thresh ~ "Significant",
                                         TRUE ~ "Not sig.")) |>
          dplyr::mutate(colour = case_when(!!sym(fdr_column) < de_q_val_thresh ~ "purple",
                                          TRUE ~ "black")) |>
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
      message(sprintf("   outputDeResultsAllContrasts: Created combined PDF with %d volcano plots: %s", 
                     length(all_volcano_plots), basename(combined_pdf_path)))
    }
    
  }, error = function(e) {
    message(sprintf("   outputDeResultsAllContrasts: ERROR creating combined PDF: %s", e$message))
  })
  
  # ============================================================================
  # NEW: Aggregate and save Num Sig DE Molecules from all contrasts
  # ============================================================================
  message("   outputDeResultsAllContrasts: Processing NumSigDeMolecules figures...")
  
  tryCatch({
    # Create NumSigDeMolecules directory
    numsigde_dir <- file.path(publication_graphs_dir, "NumSigDeMolecules")
    if (!dir.exists(numsigde_dir)) {
      dir.create(numsigde_dir, recursive = TRUE, showWarnings = TRUE)
      if (!dir.exists(numsigde_dir)) {
        stop("Failed to create NumSigDeMolecules directory: ", numsigde_dir)
      }
      message(sprintf("   outputDeResultsAllContrasts: Created NumSigDeMolecules directory: %s", numsigde_dir))
    }
    
    # Collect num_sig_de_molecules tables from all contrasts
    all_numsig_tables <- list()
    
    for (contrast_name in names(de_results_list_all_contrasts)) {
      contrast_result <- de_results_list_all_contrasts[[contrast_name]]
      
      # Check for num_sig_de_molecules_first_go in the result
      if (!is.null(contrast_result$num_sig_de_molecules_first_go)) {
        numsig_data <- contrast_result$num_sig_de_molecules_first_go
        
        # Extract the table component
        if (!is.null(numsig_data$table)) {
          all_numsig_tables[[contrast_name]] <- numsig_data$table
          message(sprintf("   outputDeResultsAllContrasts: Found NumSigDE table for contrast: %s", contrast_name))
        }
      }
    }
    
    # If we have any tables, aggregate and save them
    if (length(all_numsig_tables) > 0) {
      # Combine all tables
      combined_numsig_table <- dplyr::bind_rows(all_numsig_tables)
      
      # Write the combined table
      table_path <- file.path(numsigde_dir, paste0(file_prefix, "_num_sig_de_molecules.tab"))
      vroom::vroom_write(combined_numsig_table, table_path)
      message(sprintf("   outputDeResultsAllContrasts: Wrote NumSigDE table: %s", basename(table_path)))
      
      # Also write as Excel
      excel_path <- file.path(numsigde_dir, paste0(file_prefix, "_num_sig_de_molecules.xlsx"))
      writexl::write_xlsx(combined_numsig_table, excel_path)
      message(sprintf("   outputDeResultsAllContrasts: Wrote NumSigDE Excel: %s", basename(excel_path)))
      
      # Generate combined barplot from the aggregated data
      # Filter to only significant results
      sig_only <- combined_numsig_table |>
        dplyr::filter(status != "Not significant")
      
      if (nrow(sig_only) > 0) {
        num_sig_de_barplot <- sig_only |>
          ggplot2::ggplot(ggplot2::aes(x = status, y = counts)) +
          ggplot2::geom_bar(stat = "identity", fill = "steelblue") +
          ggplot2::geom_text(stat = 'identity', ggplot2::aes(label = counts), vjust = -0.5) +
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
            ggplot2::facet_wrap(~ comparison, scales = "free_x")
        }
        
        # Calculate appropriate width based on number of comparisons
        num_comparisons <- length(unique(sig_only$comparison))
        plot_width <- max(7, (num_comparisons + 2) * 7/6)
        
        # Save as PNG
        png_path <- file.path(numsigde_dir, paste0(file_prefix, "_num_sig_de_molecules.png"))
        ggplot2::ggsave(png_path, num_sig_de_barplot, width = plot_width, height = 6, dpi = 300)
        message(sprintf("   outputDeResultsAllContrasts: Saved NumSigDE barplot PNG: %s", basename(png_path)))
        
        # Save as PDF
        pdf_path <- file.path(numsigde_dir, paste0(file_prefix, "_num_sig_de_molecules.pdf"))
        ggplot2::ggsave(pdf_path, num_sig_de_barplot, width = plot_width, height = 6)
        message(sprintf("   outputDeResultsAllContrasts: Saved NumSigDE barplot PDF: %s", basename(pdf_path)))
        
      } else {
        message("   outputDeResultsAllContrasts: No significant DE molecules found - skipping barplot")
      }
      
    } else {
      message("   outputDeResultsAllContrasts: No NumSigDE data found in contrast results - skipping")
    }
    
  }, error = function(e) {
    message(sprintf("   outputDeResultsAllContrasts: ERROR processing NumSigDeMolecules: %s", e$message))
  })
  
  message("--- Exiting outputDeResultsAllContrasts ---")
  return(TRUE)
}) 


# ----------------------------------------------------------------------------
# getDeResultsLongFormat
# ----------------------------------------------------------------------------
# Get the differential expression results in wide format
#'@export
setMethod(f = "getDeResultsLongFormat"
          , signature = "list"
          , definition = function(objectsList)  {


            return_object_list <-  purrr::map( objectsList, function(object) {

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
                dplyr::bind_rows(.id="comparison") |>
                dplyr::mutate( comparision_short = str_split_i( comparison, "=", 1))

              # Pivot the tidy data frame to a wide format using the correct ID column
              long_results <- tidy_results |>
                # Join with original counts using the correct ID column.
                dplyr::left_join(counts_table_to_use, by = join_by( !!sym(id_col_name) == !!sym(id_col_name))) |>
                dplyr::distinct()

              print(head( long_results))

              # Assign to the correct slot and return the object
              object@results_table_long <- long_results
              return(object)
            })

            return(return_object_list)
          })


# ----------------------------------------------------------------------------
# getDeResultsWideFormat
# ----------------------------------------------------------------------------
#'@export
setMethod(f = "getDeResultsWideFormat"
          , signature = "list"
          , definition = function(objectsList
                                , qvalue_column = "fdr_qvalue"
                                , raw_pvalue_column = "raw_pvalue"
                                , log2fc_column = "logFC")  {


            return_object_list <-  purrr::map( objectsList, function(object) {

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
                dplyr::bind_rows(.id="comparison") |>
                dplyr::mutate( comparision_short = str_split_i( comparison, "=", 1))

              # Pivot the tidy data frame to a wide format using the correct ID column
              wide_results <- tidy_results |>
                tidyr::pivot_wider(id_cols = c(!!sym(id_col_name)),
                                   names_from = c(comparision_short),
                                   names_sep = ":",
                                   values_from = c(!!sym(log2fc_column), !!sym(qvalue_column), !!sym(raw_pvalue_column))) |>
                # Join with original counts using the correct ID column.
                dplyr::left_join(counts_table_to_use, by = join_by( !!sym(id_col_name) == !!sym(id_col_name))) |>
                dplyr::arrange(dplyr::across(matches(qvalue_column))) |>
                dplyr::distinct()

              # Assign to the correct slot and return the object
              object@results_table_wide <- wide_results
              return(object)
            })

            return(return_object_list)
          })

