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

# Helper function to get counts table
getCountsTable <- function(obj) {
    if (inherits(obj, "MetaboliteAssayData")) {
        message(sprintf("   Getting counts table for object of class: %s", class(obj)[1]))
        message(sprintf("   Returning metabolite_data with dimensions: %d rows, %d cols",
                      nrow(obj@metabolite_data), ncol(obj@metabolite_data)))
        obj@metabolite_data
    } else if (inherits(obj, "ProteinQuantitativeData")) {
        message(sprintf("   Returning protein_quant_table with dimensions: %d rows, %d cols",
                      nrow(obj@protein_quant_table), ncol(obj@protein_quant_table)))
        obj@protein_quant_table
    } else {
        message(sprintf("   ERROR: Unsupported object type: %s", class(obj)[1]))
        stop("Unsupported object type")
    }
}