#' MetabolomicsDifferentialAbundanceResults S4 Class
#'
#' @description
#' S4 class to store essential results from metabolomics differential abundance analysis.
#' This class contains the original data object, fitted model, and results table.
#'
#' @slot theObject The original MetaboliteAssayData object used for analysis
#' @slot fit.eb The fitted eBayes model from limma analysis
#' @slot contrasts_results_table Data frame with differential abundance statistics
#'
#' @export
setClass("MetabolomicsDifferentialAbundanceResults",
         slots = c(
           theObject = "MetaboliteAssayData",
           fit.eb = "MArrayLM",
           contrasts_results_table = "list"
         ),
         prototype = list(
           theObject = NULL,
           fit.eb = NULL,
           contrasts_results_table = list()
         )
)

##----------------------------------------------------------------------------------------------------------------------------------------------------------------------
#'@export
setGeneric( name ="differentialAbundanceAnalysis"
            , def=function(objectsList
                           , contrasts_tbl = NULL
                           , formula_string = NULL
                           , de_q_val_thresh = NULL
                           , treat_lfc_cutoff = NULL
                           , eBayes_trend = NULL
                           , eBayes_robust = NULL
                           , args_group_pattern = NULL
                           , args_row_id = NULL) {
              standardGeneric("differentialAbundanceAnalysis")})



#'@export
setMethod( f ="differentialAbundanceAnalysis"
           , signature = "list"
           , definition=function( objectsList
                                  , contrasts_tbl = NULL
                                  , formula_string = NULL
                                  , de_q_val_thresh = NULL
                                  , treat_lfc_cutoff = NULL
                                  , eBayes_trend = NULL
                                  , eBayes_robust = NULL
                                  , args_group_pattern = NULL
                                  , args_row_id = NULL ) {

             # Validate that all objects in the list are MetaboliteAssayData
             if (!all(purrr::map_lgl(objectsList, ~inherits(.x, "MetaboliteAssayData")))) {
               stop("All objects in objectsList must be of class MetaboliteAssayData")
             }

             # Run DE analysis and explicitly set names
             results_list <- purrr::map(    objectsList
                                            , \( obj) {
                                              differentialAbundanceAnalysisHelper(  obj
                                                                                    , contrasts_tbl = contrasts_tbl
                                                                                    , formula_string = formula_string
                                                                                    , de_q_val_thresh = de_q_val_thresh
                                                                                    , treat_lfc_cutoff = treat_lfc_cutoff
                                                                                    , eBayes_trend = eBayes_trend
                                                                                    , eBayes_robust = eBayes_robust
                                                                                    , args_group_pattern = args_group_pattern
                                                                                    , args_row_id = args_row_id
                                              )
                                            })

             # Set names if the input list had names
             if (!is.null(names(objectsList))) {
               names(results_list) <- names(objectsList)
             }

             return(results_list)

           })

##----------------------------------------------------------------------------------------------------------------------------------------------------------------------
#'@export
setGeneric( name ="differentialAbundanceAnalysisHelper"
            , def=function(theObject
                           , contrasts_tbl = NULL
                           , formula_string = NULL
                           , de_q_val_thresh = NULL
                           , treat_lfc_cutoff = NULL
                           , eBayes_trend = NULL
                           , eBayes_robust = NULL
                           , args_group_pattern = NULL
                           , args_row_id = NULL) {
              standardGeneric("differentialAbundanceAnalysisHelper")
            })

#'@export
setMethod( f ="differentialAbundanceAnalysisHelper"
           , signature = "MetaboliteAssayData"
           , definition=function( theObject
                                  , contrasts_tbl = NULL
                                  , formula_string = NULL
                                  , de_q_val_thresh = NULL
                                  , treat_lfc_cutoff = NULL
                                  , eBayes_trend = NULL
                                  , eBayes_robust = NULL
                                  , args_group_pattern = NULL
                                  , args_row_id = NULL ) {

  contrasts_tbl <- checkParamsObjectFunctionSimplify( theObject, "contrasts_tbl", NULL)
  formula_string <- checkParamsObjectFunctionSimplify( theObject, "formula_string", " ~ 0 + group")
  de_q_val_thresh <- checkParamsObjectFunctionSimplify( theObject, "de_q_val_thresh", 0.05)
  treat_lfc_cutoff <- checkParamsObjectFunctionSimplify( theObject, "treat_lfc_cutoff", 0)
  eBayes_trend <- checkParamsObjectFunctionSimplify( theObject, "eBayes_trend", TRUE)
  eBayes_robust <- checkParamsObjectFunctionSimplify( theObject, "eBayes_robust", TRUE)
  args_group_pattern <- checkParamsObjectFunctionSimplify( theObject, "args_group_pattern", "(\\d+)")
  args_row_id <- checkParamsObjectFunctionSimplify( theObject, "args_row_id", "uniprot_acc")

  print(formula_string)

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
  theObject <- updateParamInObject(theObject, "de_q_val_thresh")
  theObject <- updateParamInObject(theObject, "treat_lfc_cutoff")
  theObject <- updateParamInObject(theObject, "eBayes_trend")
  theObject <- updateParamInObject(theObject, "eBayes_robust")
  theObject <- updateParamInObject(theObject, "args_group_pattern")
  theObject <- updateParamInObject(theObject, "args_row_id")

  return_list <- list()
  return_list$theObject <- theObject

  message("--- Entering differentialAbundanceAnalysis ---")
  message(sprintf("   differentialAbundanceAnalysis: theObject class = %s", class(theObject)))

  # Debug: Print parameter values before conversion
  message(sprintf("   eBayes_trend value: %s (class: %s)", eBayes_trend, class(eBayes_trend)))
  message(sprintf("   eBayes_robust value: %s (class: %s)", eBayes_robust, class(eBayes_robust)))

  # Ensure logical values - handle character strings that might represent logical values
  if (is.character(eBayes_trend)) {
    eBayes_trend <- as.logical(toupper(eBayes_trend) %in% c("TRUE", "T", "1", "YES"))
    message(sprintf("   eBayes_trend converted from character '%s' to logical: %s", eBayes_trend, eBayes_trend))
  } else {
    eBayes_trend <- as.logical(eBayes_trend)
  }

  if (is.character(eBayes_robust)) {
    eBayes_robust <- as.logical(toupper(eBayes_robust) %in% c("TRUE", "T", "1", "YES"))
    message(sprintf("   eBayes_robust converted from character '%s' to logical: %s", eBayes_robust, eBayes_robust))
  } else {
    eBayes_robust <- as.logical(eBayes_robust)
  }

  message(sprintf("   eBayes_trend after conversion: %s (class: %s)", eBayes_trend, class(eBayes_trend)))
  message(sprintf("   eBayes_robust after conversion: %s (class: %s)", eBayes_robust, class(eBayes_robust)))

  ## Compare the different experimental groups and obtain lists of differentially expressed metabolites

  rownames( theObject@design_matrix ) <- theObject@design_matrix |> dplyr::pull( one_of(theObject@sample_id ))

  # Prepare data matrix for DE analysis
  data_matrix <- NA

  matrix_data <- as.matrix(theObject@metabolite_data[, -1]) # Exclude Name column
  colnames(matrix_data) <- colnames(theObject@metabolite_data)[-1]
  rownames(matrix_data) <- theObject@metabolite_data$Name
  data_matrix <- matrix_data

  contrasts_results <- runTestsContrasts(
    data_matrix,
    contrast_strings = contrasts_tbl$contrasts,
    design_matrix = theObject@design_matrix,
    formula_string = formula_string,
    treat_lfc_cutoff = treat_lfc_cutoff,
    eBayes_trend = eBayes_trend,
    eBayes_robust = eBayes_robust
  )

  # Combine all contrast results into a single data frame
  contrasts_results_table <-  contrasts_results$results

  # Map back to original group names in results if needed
  if(exists("group_mapping")) {
    contrasts_results_table <- contrasts_results_table |>
      dplyr::mutate(comparison = purrr::map_chr(comparison, \(x) {
        result <- x
        for(safe_name in names(group_mapping)) {
          result <- gsub(safe_name, group_mapping[safe_name], result, fixed = TRUE)
        }
        result
      }))
  }

  return_list$fit.eb <- contrasts_results$fit.eb
  return_list$contrasts_results_table <- contrasts_results_table

  # Create and return the S4 object
  result_object <- new("MetabolomicsDifferentialAbundanceResults",
                       theObject = return_list$theObject,
                       fit.eb = return_list$fit.eb,
                       contrasts_results_table = return_list$contrasts_results_table
  )

  return(result_object)
})

#'
#'
#' ## Create proteomics interactive volcano plot
#' #' @export
#' # de_analysis_results_list$contrasts_results$fit.eb
#' # No full stops in the nme of columns of interactive table in glimma plot. It won't display column with full stop in the column name.
#' writeInteractiveVolcanoPlotProteomics <- function( de_proteins_long
#'                                                    , uniprot_tbl
#'                                                    , fit.eb
#'                                                    , publication_graphs_dir
#'                                                    , args_row_id = "uniprot_acc"
#'                                                    , fdr_column = "fdr_qvalue"
#'                                                    , raw_p_value_column = "raw_pvalue"
#'                                                    , log2fc_column = "log2FC"
#'                                                    , de_q_val_thresh = 0.05
#'                                                    , counts_tbl = NULL
#'                                                    , groups = NULL
#'                                                    , uniprot_id_column = "Entry"
#'                                                    , gene_names_column = "Gene Names"
#'                                                    , display_columns = c( "best_uniprot_acc" )) {
#'
#'
#'   volcano_plot_tab <- de_proteins_long  |>
#'     dplyr::mutate(best_uniprot_acc = str_split(!!sym(args_row_id), ":" ) |> purrr::map_chr(1)  ) |>
#'     left_join(uniprot_tbl, by = join_by( best_uniprot_acc == !!sym( uniprot_id_column) ) ) |>
#'     dplyr::rename( UNIPROT_GENENAME = gene_names_column ) |>
#'     mutate( UNIPROT_GENENAME = purrr::map_chr( UNIPROT_GENENAME, \(x){str_split(x, " |:")[[1]][1]})) |>
#'     dplyr::mutate(gene_name = UNIPROT_GENENAME ) |>
#'     mutate( lqm = -log10(!!sym(fdr_column)))  |>
#'     dplyr::mutate(label = case_when(abs(!!sym(log2fc_column)) >= 1 & !!sym(fdr_column) >= de_q_val_thresh ~ "Not sig., logFC >= 1",
#'                                     abs(!!sym(log2fc_column)) >= 1 & !!sym(fdr_column) < de_q_val_thresh ~ "Sig., logFC >= 1",
#'                                     abs(!!sym(log2fc_column)) < 1 & !!sym(fdr_column) < de_q_val_thresh ~ "Sig., logFC < 1",
#'                                     TRUE ~ "Not sig.")) |>
#'     dplyr::mutate(colour = case_when(abs(!!sym(log2fc_column)) >= 1 & !!sym(fdr_column) >= de_q_val_thresh ~ "orange",
#'                                      abs(!!sym(log2fc_column)) >= 1 & !!sym(fdr_column) < de_q_val_thresh ~ "purple",
#'                                      abs(!!sym(log2fc_column)) < 1 & !!sym(fdr_column) < de_q_val_thresh ~ "blue",
#'                                      TRUE ~ "black")) |>
#'     dplyr::mutate(analysis_type = comparison)  |>
#'     dplyr::select( best_uniprot_acc, lqm, !!sym(fdr_column), !!sym(raw_p_value_column), !!sym(log2fc_column), comparison, label, colour,  gene_name
#'                    , any_of( display_columns ))   |>
#'     dplyr::mutate( my_alpha = case_when ( gene_name !=  "" ~ 1
#'                                           , TRUE ~ 0.5))
#'
#'   output_dir <- file.path( publication_graphs_dir
#'                            ,  "Interactive_Volcano_Plots")
#'
#'   dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
#'
#'
#'   purrr::walk( seq_len( ncol(fit.eb$coefficients))
#'                , \(coef) { # print(coef)
#'
#'                  print( paste("run volcano plot", coef) )
#'
#'                  #                    print(head( counts_tbl))
#'                  #                    print(groups)
#'                  #                    print( output_dir)
#'
#'                  getGlimmaVolcanoProteomics( fit.eb
#'                                              , coef = coef
#'                                              , volcano_plot_tab  = volcano_plot_tab
#'                                              , uniprot_column = best_uniprot_acc
#'                                              , gene_name_column = gene_name
#'                                              , display_columns = display_columns
#'                                              , counts_tbl = counts_tbl
#'                                              , groups = groups
#'                                              , output_dir = output_dir ) } )
#'
#' }
#'
#'
#'
#'
#' #' @export
#' # de_analysis_results_list$contrasts_results$fit.eb
#' # No full stops in the nme of columns of interactive table in glimma plot. It won't display column with full stop in the column name.
#' writeInteractiveVolcanoPlotProteomicsWidget <- function( de_proteins_long
#'                                                          , uniprot_tbl
#'                                                          , fit.eb
#'                                                          , args_row_id = "uniprot_acc"
#'                                                          , fdr_column = "fdr_qvalue"
#'                                                          , raw_p_value_column = "raw_pvalue"
#'                                                          , log2fc_column = "log2FC"
#'                                                          , de_q_val_thresh = 0.05
#'                                                          , counts_tbl = NULL
#'                                                          , groups = NULL
#'                                                          , uniprot_id_column = "Entry"
#'                                                          , gene_names_column = "Gene Names"
#'                                                          , display_columns = c( "best_uniprot_acc" )) {
#'
#'
#'   volcano_plot_tab <- de_proteins_long  |>
#'     left_join(uniprot_tbl, by = join_by( !!sym(args_row_id) == !!sym( uniprot_id_column) ) ) |>
#'     dplyr::rename( UNIPROT_GENENAME = gene_names_column ) |>
#'     mutate( UNIPROT_GENENAME = purrr::map_chr( UNIPROT_GENENAME, \(x){str_split(x, " ")[[1]][1]})) |>
#'     mutate( lqm = -log10(!!sym(fdr_column)))  |>
#'     dplyr::mutate(label = case_when(abs(!!sym(log2fc_column)) >= 1 & !!sym(fdr_column) >= de_q_val_thresh ~ "Not sig., logFC >= 1",
#'                                     abs(!!sym(log2fc_column)) >= 1 & !!sym(fdr_column) < de_q_val_thresh ~ "Sig., logFC >= 1",
#'                                     abs(!!sym(log2fc_column)) < 1 & !!sym(fdr_column) < de_q_val_thresh ~ "Sig., logFC < 1",
#'                                     TRUE ~ "Not sig.")) |>
#'     dplyr::mutate(colour = case_when(abs(!!sym(log2fc_column)) >= 1 & !!sym(fdr_column) >= de_q_val_thresh ~ "orange",
#'                                      abs(!!sym(log2fc_column)) >= 1 & !!sym(fdr_column) < de_q_val_thresh ~ "purple",
#'                                      abs(!!sym(log2fc_column)) < 1 & !!sym(fdr_column) < de_q_val_thresh ~ "blue",
#'                                      TRUE ~ "black")) |>
#'     dplyr::mutate(gene_name = str_split(UNIPROT_GENENAME, " |:" ) |> purrr::map_chr(1)  ) |>
#'     dplyr::mutate(best_uniprot_acc = str_split(!!sym(args_row_id), ":" ) |> purrr::map_chr(1)  ) |>
#'     dplyr::mutate(analysis_type = comparison)  |>
#'     dplyr::select( best_uniprot_acc, lqm, !!sym(fdr_column), !!sym(raw_p_value_column), !!sym(log2fc_column), comparison, label, colour,  gene_name
#'                    , any_of( display_columns ))   |>
#'     dplyr::mutate( my_alpha = case_when ( gene_name !=  "" ~ 1
#'                                           , TRUE ~ 0.5))
#'
#'   interactive_volcano_plots <- purrr::map(seq_len( ncol(fit.eb$coefficients))
#'                                           , \(coef) {
#'                                             # print(paste0( "coef = ", coef))
#'                                             getGlimmaVolcanoProteomicsWidget( fit.eb
#'                                                                               , coef = coef
#'                                                                               , volcano_plot_tab  = volcano_plot_tab
#'                                                                               , uniprot_column = best_uniprot_acc
#'                                                                               , gene_name_column = gene_name
#'                                                                               , display_columns = display_columns
#'                                                                               , counts_tbl = counts_tbl
#'                                                                               , groups = groups ) } )
#'
#'   interactive_volcano_plots
#' }
#'
#' #' @export
#' outputDeAnalysisResults <- function(de_analysis_results_list
#'                                     , theObject
#'                                     , uniprot_tbl
#'                                     , de_output_dir = NULL
#'                                     , publication_graphs_dir = NULL
#'                                     , file_prefix = NULL
#'                                     , plots_format = NULL
#'                                     , args_row_id = NULL
#'                                     , de_q_val_thresh = NULL
#'                                     , gene_names_column = NULL
#'                                     , fdr_column = NULL
#'                                     , raw_p_value_column = NULL
#'                                     , log2fc_column = NULL
#'                                     , uniprot_id_column = NULL
#'                                     , display_columns = NULL
#' ) {
#'
#'
#'   uniprot_tbl <- checkParamsObjectFunctionSimplify(theObject, "uniprot_tbl", NULL)
#'   de_output_dir <- checkParamsObjectFunctionSimplify(theObject, "de_output_dir", NULL)
#'   publication_graphs_dir <- checkParamsObjectFunctionSimplify(theObject, "publication_graphs_dir", NULL)
#'   file_prefix <- checkParamsObjectFunctionSimplify(theObject, "file_prefix", "de_proteins")
#'   plots_format <- checkParamsObjectFunctionSimplify(theObject, "plots_format", c("pdf", "png"))
#'   args_row_id <- checkParamsObjectFunctionSimplify(theObject, "args_row_id", "uniprot_acc")
#'   de_q_val_thresh <- checkParamsObjectFunctionSimplify(theObject, "de_q_val_thresh", 0.05)
#'   gene_names_column <- checkParamsObjectFunctionSimplify(theObject, "gene_names_column", "Gene Names")
#'   fdr_column <- checkParamsObjectFunctionSimplify(theObject, "fdr_column", "fdr_qvalue")
#'   raw_p_value_column <- checkParamsObjectFunctionSimplify(theObject, "raw_p_value_column", "raw_pvalue")
#'   log2fc_column <- checkParamsObjectFunctionSimplify(theObject, "log2fc_column", "log2FC")
#'   uniprot_id_column <- checkParamsObjectFunctionSimplify(theObject, "uniprot_id_column", "Entry")
#'   display_columns <- checkParamsObjectFunctionSimplify(theObject, "display_columns", c( "best_uniprot_acc" ))
#'
#'   theObject <- updateParamInObject(theObject, "uniprot_tbl")
#'   theObject <- updateParamInObject(theObject, "de_output_dir")
#'   theObject <- updateParamInObject(theObject, "publication_graphs_dir")
#'   theObject <- updateParamInObject(theObject, "file_prefix")
#'   theObject <- updateParamInObject(theObject, "plots_format")
#'   theObject <- updateParamInObject(theObject, "args_row_id")
#'   theObject <- updateParamInObject(theObject, "de_q_val_thresh")
#'   theObject <- updateParamInObject(theObject, "gene_names_column")
#'   theObject <- updateParamInObject(theObject, "fdr_column")
#'   theObject <- updateParamInObject(theObject, "raw_p_value_column")
#'   theObject <- updateParamInObject(theObject, "log2fc_column")
#'   theObject <- updateParamInObject(theObject, "uniprot_id_column")
#'   theObject <- updateParamInObject(theObject, "display_columns")
#'
#'   ## PCA plot
#'   plot_pca_plot <- de_analysis_results_list$pca_plot
#'
#'   dir.create(file.path( publication_graphs_dir, "PCA")
#'              , recursive = TRUE
#'              , showWarnings = FALSE)
#'
#'   for( format_ext in plots_format) {
#'     file_name <- file.path( publication_graphs_dir, "PCA", paste0("PCA_plot.",format_ext))
#'     ggsave(filename = file_name, plot = plot_pca_plot, limitsize = FALSE)
#'   }
#'
#'   plot_pca_plot_with_labels <- de_analysis_results_list$pca_plot_with_labels
#'   for( format_ext in plots_format) {
#'     file_name <- file.path( publication_graphs_dir, "PCA", paste0("PCA_plot_with_sample_ids.",format_ext))
#'     ggsave(filename = file_name, plot = plot_pca_plot_with_labels, limitsize = FALSE)
#'   }
#'
#'   ## RLE plot
#'   plot_rle_plot <- de_analysis_results_list$rle_plot
#'
#'   dir.create(file.path( publication_graphs_dir, "RLE")
#'              , recursive = TRUE
#'              , showWarnings = FALSE)
#'
#'   for( format_ext in plots_format) {
#'     file_name <- file.path( publication_graphs_dir, "RLE", paste0("RLE_plot.",format_ext))
#'     ggsave(filename = file_name, plot = plot_rle_plot, limitsize = FALSE)
#'   }
#'
#'   ## Save the number of values graph
#'   plot_num_of_values <- de_analysis_results_list$plot_num_of_values
#'
#'   for( format_ext in plots_format) {
#'     file_name <- file.path(de_output_dir, paste0("num_of_values.",format_ext))
#'     ggsave(filename = file_name, plot = plot_num_of_values, limitsize = FALSE)
#'   }
#'
#'   ## Contrasts results
#'   ## This plot is used to check the mean-variance relationship of the expression data, after fitting a linear model.
#'   pdf(file.path(de_output_dir, "plotSA_after_ruvIII.pdf" ))
#'   plotSA(de_analysis_results_list$contrasts_results$fit.eb)
#'   dev.off()
#'
#'   png(file.path(de_output_dir, "plotSA_after_ruvIII.png" ))
#'   plotSA(de_analysis_results_list$contrasts_results$fit.eb)
#'   dev.off()
#'
#'   saveRDS( de_analysis_results_list$contrasts_results$fit.eb,
#'            file.path(de_output_dir, "fit.eb.RDS" ) )
#'
#'   ## Values for volcano plts
#'
#'   ## Write all the results in one single table
#'   significant_rows <- de_analysis_results_list$significant_rows
#'
#'   significant_rows |>
#'     dplyr::select(-colour, -lqm) |>
#'     vroom::vroom_write(file.path(de_output_dir, "lfc_qval_long.tsv"))
#'
#'   significant_rows |>
#'     dplyr::select(-colour, -lqm) |>
#'     writexl::write_xlsx(file.path(de_output_dir, "lfc_qval_long.xlsx"))
#'
#'   ## Print Volcano plot
#'   volplot_plot <- de_analysis_results_list$volplot_plot
#'
#'   for( format_ext in plots_format) {
#'     file_name <- file.path(de_output_dir, paste0("volplot_gg_all.",format_ext))
#'     ggsave(filename = file_name, plot = volplot_plot, width = 7.29, height = 6)
#'   }
#'
#'
#'   ## Number of values graph
#'   plot_num_of_values <- de_analysis_results_list$plot_num_of_values
#'
#'   for( format_ext in plots_format) {
#'     file_name <- file.path(de_output_dir, paste0("num_of_values.",format_ext))
#'     ggsave(filename = file_name, plot = plot_num_of_values, limitsize = FALSE)
#'   }
#'
#'   ## Contrasts results
#'   ## This plot is used to check the mean-variance relationship of the expression data, after fitting a linear model.
#'   contrasts_results <- de_analysis_results_list$contrasts_results
#'   for( format_ext in plots_format) {
#'     file_name <- file.path(de_output_dir, paste0("plotSA_after_ruvIII",format_ext))
#'
#'     if( format_ext == "pdf") {
#'       pdf(file_name)
#'     } else if(format_ext == "png") {
#'       png(file_name)
#'     }
#'
#'     plotSA(contrasts_results$fit.eb)
#'     dev.off()
#'
#'   }
#'
#'   saveRDS( contrasts_results$fit.eb,
#'            file.path(de_output_dir, "fit.eb.RDS" ) )
#'
#'   ## Values for volcano plts
#'
#'   ## Write all the results in one single table
#'   significant_rows <- de_analysis_results_list$significant_rows
#'   significant_rows |>
#'     dplyr:::select(-colour, -lqm) |>
#'     vroom::vroom_write(file.path(de_output_dir, "lfc_qval_long.tsv"))
#'
#'   significant_rows |>
#'     dplyr:::select(-colour, -lqm) |>
#'     writexl::write_xlsx(file.path(de_output_dir, "lfc_qval_long.xlsx"))
#'
#'
#'   ## Count the number of up or down significnat differentially expressed proteins.
#'   if( !is.null(de_analysis_results_list$num_sig_de_genes_barplot_only_significant)) {
#'     num_sig_de_genes_barplot_only_significant <- de_analysis_results_list$num_sig_de_genes_barplot_only_significant
#'     num_of_comparison_only_significant <- de_analysis_results_list$num_of_comparison_only_significant
#'
#'     savePlot(num_sig_de_genes_barplot_only_significant,
#'              base_path = de_output_dir,
#'              plot_name = paste0(file_prefix, "_num_sda_entities_barplot_only_significant"),
#'              formats =  c("pdf", "png", "svg"),
#'              width = (num_of_comparison_only_significant + 2) *7/6,
#'              height = 6)
#'
#'
#'
#'   }
#'
#'
#'   ## Count the number of up or down significnat differentially expressed proteins.
#'   num_sig_de_molecules_first_go <- de_analysis_results_list$num_sig_de_molecules_first_go
#'   vroom::vroom_write(num_sig_de_molecules_first_go$table,
#'                      file.path(de_output_dir,
#'                                paste0(file_prefix, "_num_significant_differentially_abundant_all.tab") ))
#'
#'   writexl::write_xlsx(num_sig_de_molecules_first_go$table,
#'                       file.path(de_output_dir,
#'                                 paste0(file_prefix, "_num_significant_differentially_abundant_all.xlsx")) )
#'
#'
#'
#'   ## Print p-values distribution figure
#'   pvalhist <- de_analysis_results_list$pvalhist
#'   for( format_ext in plots_format) {
#'     file_name<-file.path(de_output_dir,paste0(file_prefix, "_p_values_distn.",format_ext))
#'     ggsave(filename = file_name,
#'            plot = pvalhist,
#'            height = 10,
#'            width = 7)
#'
#'   }
#'
#'
#'
#'   ## Create wide format output file
#'   de_proteins_wide <- de_analysis_results_list$de_proteins_wide
#'   vroom::vroom_write( de_proteins_wide,
#'                       file.path( de_output_dir,
#'                                  paste0(file_prefix, "_wide.tsv")))
#'
#'   writexl::write_xlsx( de_proteins_wide,
#'                        file.path( de_output_dir,
#'                                   paste0(file_prefix, "_wide.xlsx")))
#'
#'   de_proteins_wide_annot <- de_proteins_wide |>
#'     mutate( uniprot_acc_cleaned = str_split( !!sym(args_row_id), "-" )  |>
#'               purrr::map_chr(1) ) |>
#'     left_join( uniprot_tbl, by = join_by( uniprot_acc_cleaned == Entry ) ) |>
#'     dplyr::select( -uniprot_acc_cleaned)  |>
#'     mutate( gene_name = purrr::map_chr( !!sym(gene_names_column), \(x){str_split(x, " |:")[[1]][1]})) |>
#'     relocate( gene_name, .after = !!sym(args_row_id))
#'
#'   vroom::vroom_write( de_proteins_wide_annot,
#'                       file.path( de_output_dir,
#'                                  paste0(file_prefix, "_wide_annot.tsv")))
#'
#'   writexl::write_xlsx( de_proteins_wide_annot,
#'                        file.path( de_output_dir,
#'                                   paste0(file_prefix, "_wide_annot.xlsx")))
#'
#'   ## Create long format output file
#'   de_proteins_long <- de_analysis_results_list$de_proteins_long
#'   vroom::vroom_write( de_proteins_long,
#'                       file.path( de_output_dir,
#'                                  paste0(file_prefix, "_long.tsv")))
#'
#'   writexl::write_xlsx( de_proteins_long,
#'                        file.path( de_output_dir,
#'                                   paste0(file_prefix, "_long.xlsx")))
#'
#'   de_proteins_long_annot <- de_proteins_long |>
#'     mutate( uniprot_acc_cleaned = str_split( !!sym(args_row_id), "-" )  |>
#'               purrr::map_chr(1) )|>
#'     left_join( uniprot_tbl, by = join_by( uniprot_acc_cleaned == Entry ) )  |>
#'     dplyr::select( -uniprot_acc_cleaned)  |>
#'     mutate( gene_name = purrr::map_chr( !!sym(gene_names_column), \(x){str_split(x, " |:")[[1]][1]})) |>
#'     relocate( gene_name, .after = !!sym(args_row_id))
#'
#'   vroom::vroom_write( de_proteins_long_annot,
#'                       file.path( de_output_dir,
#'                                  paste0(file_prefix, "_long_annot.tsv")))
#'
#'   writexl::write_xlsx( de_proteins_long_annot,
#'                        file.path( de_output_dir,
#'                                   paste0(file_prefix, "_long_annot.xlsx")))
#'
#'   ## Static volcano plots
#'   dir.create(file.path( publication_graphs_dir, "Volcano_Plots")
#'              , recursive = TRUE
#'              , showWarnings = FALSE)
#'
#'   list_of_volcano_plots <- de_analysis_results_list$list_of_volcano_plots
#'
#'   # Print diagnostic info about the volcano plots
#'   message(sprintf("Number of volcano plots: %d", nrow(list_of_volcano_plots)))
#'
#'   purrr::walk2( list_of_volcano_plots %>% dplyr::pull(title),
#'                 list_of_volcano_plots %>% dplyr::pull(plot),
#'                 \(x,y){
#'                   # gg_save_logging ( .y, file_name_part, plots_format)
#'
#'                   savePlot( y
#'                             , base_path = file.path( publication_graphs_dir, "Volcano_Plots")
#'                             , plot_name =  x
#'                             , formats = plots_format, width = 7, height = 7)
#'
#'                 })
#'
#'   # Generate a multi-page PDF with all volcano plots
#'   volcano_plots_list <- list_of_volcano_plots %>% dplyr::pull(plot)
#'
#'   # Generate combined PDF with all plots, one per page
#'   pdf_file <- file.path(publication_graphs_dir, "Volcano_Plots", "list_of_volcano_plots.pdf")
#'   pdf(file = pdf_file, width = 7, height = 7, onefile = TRUE)
#'   purrr::walk(volcano_plots_list, print)
#'   invisible(dev.off())
#'
#'   # Verify the PDF was created with the right number of pages
#'   message(sprintf("Created multi-page PDF at %s", pdf_file))
#'
#'   list_of_volcano_plots_with_gene_names <- de_analysis_results_list$list_of_volcano_plots_with_gene_names
#'
#'   # Print diagnostic info about the labeled volcano plots
#'   message(sprintf("Number of labeled volcano plots: %d", nrow(list_of_volcano_plots_with_gene_names)))
#'
#'   purrr::walk2( list_of_volcano_plots_with_gene_names %>% dplyr::pull(title)
#'                 , list_of_volcano_plots_with_gene_names %>% dplyr::pull(plot)
#'                 , \(x, y) {
#'
#'                   savePlot( x
#'                             , file.path( publication_graphs_dir, "Volcano_Plots")
#'                             , paste0( y,"_with_protein_labels"))
#'                 })
#'
#'   # Generate a multi-page PDF with all labeled volcano plots
#'   volcano_plots_with_genes_list <- list_of_volcano_plots_with_gene_names %>% dplyr::pull(plot)
#'
#'   # Generate combined PDF with all labeled plots, one per page
#'   pdf_file_with_genes <- file.path(publication_graphs_dir, "Volcano_Plots", "list_of_volcano_plots_with_gene_names.pdf")
#'   pdf(file = pdf_file_with_genes, width = 7, height = 7, onefile = TRUE)
#'   purrr::walk(volcano_plots_with_genes_list, print)
#'   invisible(dev.off())
#'
#'   # Verify the labeled PDF was created with the right number of pages
#'   message(sprintf("Created multi-page labeled PDF at %s", pdf_file_with_genes))
#'
#'   ## Number of significant molecules
#'   createDirIfNotExists(file.path(publication_graphs_dir, "NumSigDeMolecules"))
#'   vroom::vroom_write( de_analysis_results_list$num_sig_de_molecules,
#'                       file.path(publication_graphs_dir, "NumSigDeMolecules", paste0(file_prefix, "_num_sig_de_molecules.tab" ) ))
#'
#'
#'   if( !is.null(de_analysis_results_list$num_sig_de_genes_barplot_only_significant)) {
#'     num_sig_de_genes_barplot_only_significant <- de_analysis_results_list$num_sig_de_genes_barplot_only_significant
#'     num_of_comparison_only_significant <- de_analysis_results_list$num_of_comparison_only_significant
#'
#'     savePlot(num_sig_de_genes_barplot_only_significant,
#'              base_path = file.path(publication_graphs_dir, "NumSigDeMolecules"),
#'              plot_name = paste0(file_prefix, "_num_sig_de_molecules."),
#'              formats = plots_format,
#'              width = (num_of_comparison_only_significant + 2) *7/6,
#'              height = 6)
#'
#'
#'   }
#'
#'
#'   if(!is.null( de_analysis_results_list$num_sig_de_genes_barplot_with_not_significant )) {
#'     num_sig_de_genes_barplot_with_not_significant <- de_analysis_results_list$num_sig_de_genes_barplot_with_not_significant
#'     num_of_comparison_with_not_significant <- de_analysis_results_list$num_of_comparison_with_not_significant
#'
#'     print("print bar plot")
#'
#'     savePlot(num_sig_de_genes_barplot_with_not_significant,
#'              base_path = file.path(publication_graphs_dir, "NumSigDeMolecules"),
#'              plot_name = paste0(file_prefix, "_num_sig_de_molecules_with_not_significant"),
#'              formats = plots_format,
#'              width = (num_of_comparison_with_not_significant + 2) *7/6,
#'              height = 6)
#'   }
#'
#'   ## Write interactive volcano plot
#'   counts_mat <- getCountsTable(de_analysis_results_list$theObject) |>
#'     column_to_rownames((de_analysis_results_list$theObject)@protein_id_column  ) |>
#'     as.matrix()
#'
#'   this_design_matrix <- de_analysis_results_list$theObject@design_matrix
#'
#'   rownames( this_design_matrix ) <- this_design_matrix[,de_analysis_results_list$theObject@sample_id]
#'
#'   this_groups <- this_design_matrix[colnames( counts_mat), "group"]
#'
#'   writeInteractiveVolcanoPlotProteomics( de_proteins_long
#'                                          , uniprot_tbl = uniprot_tbl
#'                                          , fit.eb = contrasts_results$fit.eb
#'                                          , publication_graphs_dir= publication_graphs_dir
#'                                          , args_row_id = args_row_id
#'                                          , fdr_column = fdr_column
#'                                          , raw_p_value_column = raw_p_value_column
#'                                          , log2fc_column = log2fc_column
#'                                          , de_q_val_thresh = de_q_val_thresh
#'                                          , counts_tbl = counts_mat
#'
#'                                          , groups = this_groups
#'                                          , uniprot_id_column = uniprot_id_column
#'                                          , gene_names_column = gene_names_column
#'                                          , display_columns = display_columns )
#'
#' }
#'
#'
#'
#' #' @export
#' writeInteractiveVolcanoPlotProteomicsMain <- function(de_analysis_results_list
#'                                                       , theObject
#'                                                       , uniprot_tbl
#'                                                       , publication_graphs_dir = NULL
#'                                                       , file_prefix = NULL
#'                                                       , plots_format = NULL
#'                                                       , args_row_id = NULL
#'                                                       , de_q_val_thresh = NULL
#'                                                       , gene_names_column = NULL
#'                                                       , fdr_column = NULL
#'                                                       , raw_p_value_column = NULL
#'                                                       , log2fc_column = NULL
#'                                                       , uniprot_id_column = NULL
#'                                                       , display_columns = NULL
#' ) {
#'
#'
#'   uniprot_tbl <- checkParamsObjectFunctionSimplify(theObject, "uniprot_tbl", NULL)
#'   publication_graphs_dir <- checkParamsObjectFunctionSimplify(theObject, "publication_graphs_dir", NULL)
#'   args_row_id <- checkParamsObjectFunctionSimplify(theObject, "args_row_id", "uniprot_acc")
#'   de_q_val_thresh <- checkParamsObjectFunctionSimplify(theObject, "de_q_val_thresh", 0.05)
#'   gene_names_column <- checkParamsObjectFunctionSimplify(theObject, "gene_names_column", "Gene Names")
#'   fdr_column <- checkParamsObjectFunctionSimplify(theObject, "fdr_column", "fdr_qvalue")
#'   raw_p_value_column <- checkParamsObjectFunctionSimplify(theObject, "raw_p_value_column", "raw_pvalue")
#'   log2fc_column <- checkParamsObjectFunctionSimplify(theObject, "log2fc_column", "log2FC")
#'   uniprot_id_column <- checkParamsObjectFunctionSimplify(theObject, "uniprot_id_column", "Entry")
#'   display_columns <- checkParamsObjectFunctionSimplify(theObject, "display_columns", c( "best_uniprot_acc" ))
#'
#'   theObject <- updateParamInObject(theObject, "uniprot_tbl")
#'   theObject <- updateParamInObject(theObject, "publication_graphs_dir")
#'   theObject <- updateParamInObject(theObject, "args_row_id")
#'   theObject <- updateParamInObject(theObject, "de_q_val_thresh")
#'   theObject <- updateParamInObject(theObject, "gene_names_column")
#'   theObject <- updateParamInObject(theObject, "fdr_column")
#'   theObject <- updateParamInObject(theObject, "raw_p_value_column")
#'   theObject <- updateParamInObject(theObject, "log2fc_column")
#'   theObject <- updateParamInObject(theObject, "uniprot_id_column")
#'   theObject <- updateParamInObject(theObject, "display_columns")
#'
#'
#'   ## Write interactive volcano plot
#'
#'   de_proteins_long <- de_analysis_results_list$de_proteins_long
#'   contrasts_results <- de_analysis_results_list$contrasts_results
#'
#'   counts_mat <- getCountsTable(de_analysis_results_list$theObject) |>
#'     column_to_rownames((de_analysis_results_list$theObject)@protein_id_column  ) |>
#'     as.matrix()
#'
#'   this_design_matrix <- de_analysis_results_list$theObject@design_matrix
#'
#'   rownames( this_design_matrix ) <- this_design_matrix[,de_analysis_results_list$theObject@sample_id]
#'
#'   this_groups <- this_design_matrix[colnames( counts_mat), "group"]
#'
#'   writeInteractiveVolcanoPlotProteomics( de_proteins_long
#'                                          , uniprot_tbl = uniprot_tbl
#'                                          , fit.eb = contrasts_results$fit.eb
#'                                          , publication_graphs_dir= publication_graphs_dir
#'                                          , args_row_id = args_row_id
#'                                          , fdr_column = fdr_column
#'                                          , raw_p_value_column = raw_p_value_column
#'                                          , log2fc_column = log2fc_column
#'                                          , de_q_val_thresh = de_q_val_thresh
#'                                          , counts_tbl = counts_mat
#'
#'                                          , groups = this_groups
#'                                          , uniprot_id_column = uniprot_id_column
#'                                          , gene_names_column = gene_names_column
#'                                          , display_columns = display_columns )
#'
#' }
#'
#' # Helper function to get data matrix
#' getDataMatrix <- function(obj) {
#'   if (inherits(obj, "MetaboliteAssayData")) {
#'     message(sprintf("   Getting data matrix for object of class: %s", class(obj)[1]))
#'     message(sprintf("   Processing MetaboliteAssayData"))
#'     message(sprintf("   Metabolite data dimensions: %d rows, %d cols",
#'                     nrow(obj@metabolite_data), ncol(obj@metabolite_data)))
#'     matrix_data <- as.matrix(obj@metabolite_data[, -1]) # Exclude Name column
#'     colnames(matrix_data) <- colnames(obj@metabolite_data)[-1]
#'     rownames(matrix_data) <- obj@metabolite_data$Name
#'     message(sprintf("   Created matrix with dimensions: %d rows, %d cols",
#'                     nrow(matrix_data), ncol(matrix_data)))
#'     matrix_data
#'   } else if (inherits(obj, "ProteinQuantitativeData")) {
#'     message(sprintf("   Processing ProteinQuantitativeData"))
#'     message(sprintf("   Protein quant table dimensions: %d rows, %d cols",
#'                     nrow(obj@protein_quant_table), ncol(obj@protein_quant_table)))
#'     result <- as.matrix(column_to_rownames(obj@protein_quant_table, obj@protein_id_column))
#'     message(sprintf("   Created matrix with dimensions: %d rows, %d cols",
#'                     nrow(result), ncol(result)))
#'     result
#'   } else {
#'     message(sprintf("   ERROR: Unsupported object type: %s", class(obj)[1]))
#'     stop("Unsupported object type")
#'   }
#' }
#'
#' # Helper function to get counts table
#' getCountsTable <- function(obj) {
#'   if (inherits(obj, "MetaboliteAssayData")) {
#'     message(sprintf("   Getting counts table for object of class: %s", class(obj)[1]))
#'     message(sprintf("   Returning metabolite_data with dimensions: %d rows, %d cols",
#'                     nrow(obj@metabolite_data), ncol(obj@metabolite_data)))
#'     obj@metabolite_data
#'   } else if (inherits(obj, "ProteinQuantitativeData")) {
#'     message(sprintf("   Returning protein_quant_table with dimensions: %d rows, %d cols",
#'                     nrow(obj@protein_quant_table), ncol(obj@protein_quant_table)))
#'     obj@protein_quant_table
#'   } else {
#'     message(sprintf("   ERROR: Unsupported object type: %s", class(obj)[1]))
#'     stop("Unsupported object type")
#'   }
#' }
