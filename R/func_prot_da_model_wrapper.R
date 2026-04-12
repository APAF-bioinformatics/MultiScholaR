# ----------------------------------------------------------------------------
# daAnalysisWrapperFunction
# ----------------------------------------------------------------------------
#' @export
daAnalysisWrapperFunction <- function(
  theObject,
  contrasts_tbl = NULL,
  formula_string = NULL,
  group_id = NULL,
  da_q_val_thresh = NULL,
  treat_lfc_cutoff = NULL,
  eBayes_trend = NULL,
  eBayes_robust = NULL,
  args_group_pattern = NULL,
  args_row_id = NULL,
  qvalue_column = "fdr_qvalue",
  raw_pvalue_column = "raw_pvalue"
) {
  contrasts_tbl <- checkParamsObjectFunctionSimplify(theObject, "contrasts_tbl", NULL)
  formula_string <- checkParamsObjectFunctionSimplify(theObject, "formula_string", " ~ 0 + group")
  group_id <- checkParamsObjectFunctionSimplify(theObject, "group_id", "group")
  da_q_val_thresh <- checkParamsObjectFunctionSimplify(theObject, "da_q_val_thresh", 0.05)
  treat_lfc_cutoff <- checkParamsObjectFunctionSimplify(theObject, "treat_lfc_cutoff", 0)
  eBayes_trend <- checkParamsObjectFunctionSimplify(theObject, "eBayes_trend", TRUE)
  eBayes_robust <- checkParamsObjectFunctionSimplify(theObject, "eBayes_robust", TRUE)
  args_group_pattern <- checkParamsObjectFunctionSimplify(theObject, "args_group_pattern", "(\\d+)")
  args_row_id <- checkParamsObjectFunctionSimplify(theObject, "args_row_id", "uniprot_acc")

  # Add preprocessing for group names that start with numbers
  design_matrix <- theObject@design_matrix
  group_col <- design_matrix[[theObject@group_id]]

  # Check if any group names start with numbers and create mapping
  starts_with_number <- grepl("^[0-9]", group_col)
  if (any(starts_with_number)) {
    original_groups <- unique(group_col)
    safe_groups <- purrr::map_chr(original_groups, \(x) {
      if (grepl("^[0-9]", x)) paste0("grp_", x) else x
    })
    group_mapping <- setNames(original_groups, safe_groups)

    # Update design matrix with safe names
    design_matrix[[theObject@group_id]] <- purrr::map_chr(group_col, \(x) {
      if (grepl("^[0-9]", x)) paste0("grp_", x) else x
    })

    # Update contrasts table if it exists
    if (!is.null(contrasts_tbl)) {
      contrasts_tbl[[1]] <- purrr::map_chr(contrasts_tbl[[1]], \(x) {
        for (orig in names(group_mapping)) {
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
  theObject <- updateParamInObject(theObject, "da_q_val_thresh")
  theObject <- updateParamInObject(theObject, "treat_lfc_cutoff")
  theObject <- updateParamInObject(theObject, "eBayes_trend")
  theObject <- updateParamInObject(theObject, "eBayes_robust")
  theObject <- updateParamInObject(theObject, "args_group_pattern")
  theObject <- updateParamInObject(theObject, "args_row_id")

  return_list <- list()
  return_list$theObject <- theObject

  ## plot RLE plot
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

  ## plot PCA plot

  pca_plot <- plotPca(theObject,
    grouping_variable = theObject@group_id,
    label_column = "",
    title = "",
    font_size = 8
  ) +
    theme_bw() +
    theme(axis.text.x = element_text(size = 12)) +
    theme(axis.text.y = element_text(size = 12)) +
    theme(axis.title.x = element_text(size = 12)) +
    theme(axis.title.y = element_text(size = 12)) +
    theme(plot.title = element_text(size = 12)) +
    theme(legend.text = element_text(size = 12)) +
    theme(legend.title = element_text(size = 12))

  return_list$pca_plot <- pca_plot


  pca_plot_with_labels <- plotPca(theObject,
    grouping_variable = theObject@group_id,
    label_column = theObject@sample_id,
    title = "",
    font_size = 8
  ) +
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
  } else if (inherits(theObject, "ProteinQuantitativeData")) {
    return_list$plot_num_of_values <- plotNumOfValuesNoLog(theObject@protein_quant_table)
  }

  ## Compare the different experimental groups and obtain lists of differentially expressed proteins.")

  rownames(theObject@design_matrix) <- theObject@design_matrix |> dplyr::pull(one_of(theObject@sample_id))

  # Check if object contains DPC-Quant results and use limpa dpc if available
  use_dpc_de <- FALSE
  dpc_quant_results <- NULL

  if (!is.null(theObject@args$limpa_dpc_quant_results)) {
    dpc_quant_results <- theObject@args$limpa_dpc_quant_results
    use_dpc_da <- TRUE
    cat("   DA ANALYSIS Step: Detected DPC-Quant results - using limpa dpcDA for uncertainty-weighted analysis\n")
    cat("   DA ANALYSIS Step: DPC parameters used:", paste(dpc_quant_results$dpc_parameters_used, collapse = ", "), "\n")
  } else {
    cat("   DA ANALYSIS Step: No DPC-Quant results found - using standard limma analysis\n")
  }

  # CRITICAL FIX: Use the correct column for contrast strings
  # The downstream functions expect "comparison=expression" format (from full_format column)
  # NOT just the raw contrast expression (from contrasts column)
  if ("full_format" %in% names(contrasts_tbl)) {
    contrast_strings_to_use <- contrasts_tbl$full_format # Use full_format column
    cat("   DA ANALYSIS Step: Using full_format column for contrast strings\n")
  } else {
    cat("   DA ANALYSIS Step: No full_format column found, auto-generating from raw contrasts\n")
    # Auto-generate full_format column from raw contrasts
    raw_contrasts <- contrasts_tbl[, 1][[1]]

    # Generate friendly names and full format
    full_format_strings <- sapply(raw_contrasts, function(contrast_string) {
      if (grepl("=", contrast_string)) {
        return(contrast_string)
      }
      # Remove "group" prefixes if present for friendly name
      clean_string <- gsub("^group", "", contrast_string)
      clean_string <- gsub("-group", "-", clean_string)

      # Create friendly name by replacing - with _vs_
      friendly_name <- gsub("-", "_vs_", clean_string)

      # Create full format: friendly_name=original_contrast_string
      paste0(friendly_name, "=", contrast_string)
    })

    contrast_strings_to_use <- full_format_strings
    cat("   DA ANALYSIS Step: Auto-generated full_format strings:\n")
    print(contrast_strings_to_use)
  }

  # Run differential expression analysis
  if (use_dpc_de && requireNamespace("limpa", quietly = TRUE)) {
    cat("   DA ANALYSIS Step: Running limpa dpcDE analysis with uncertainty weights\n")

    # The EList is now pre-filtered and stored, so we can use it directly
    y_elist_filtered <- dpc_quant_results$quantified_elist

    if (is.null(y_elist_filtered)) {
      stop("FATAL: The quantified_elist is missing from the object's args. It should have been created by proteinMissingValueImputationLimpa.")
    }

    # Create design matrix for dpcDA
    design_matrix_for_dpcda <- model.matrix(as.formula(formula_string), theObject@design_matrix)

    cat("   DA ANALYSIS Step: Calling limpa::dpcDE\n")
    cat("   DA ANALYSIS Step: Protein matrix dims:", nrow(y_elist_filtered$E), "x", ncol(y_elist_filtered$E), "\n")
    cat("   DA ANALYSIS Step: Design matrix dims:", nrow(design_matrix_for_dpcde), "x", ncol(design_matrix_for_dpcde), "\n")

    # Run dpcDE using the synchronized EList
    dpc_fit <- limpa::dpcDE(y_elist_filtered, design_matrix_for_dpcde, plot = FALSE)

    # Convert dpcDE results to format compatible with runTestsContrasts
    contrasts_results <- convertDpcDAToStandardFormat(
      dpc_fit = dpc_fit,
      contrast_strings = contrast_strings_to_use,
      design_matrix = design_matrix_for_dpcde,
      eBayes_trend = as.logical(eBayes_trend),
      eBayes_robust = as.logical(eBayes_robust)
    )

    cat("   DA ANALYSIS Step: dpcDE analysis completed successfully\n")
  } else {
    # Standard limma analysis (existing code)
    cat("   DA ANALYSIS Step: Running standard limma analysis\n")

    data_matrix <- getDataMatrix(theObject)

    contrasts_results <- runTestsContrasts(data_matrix,
      contrast_strings = contrast_strings_to_use,
      design_matrix = theObject@design_matrix,
      formula_string = formula_string,
      weights = NA,
      treat_lfc_cutoff = as.double(treat_lfc_cutoff),
      eBayes_trend = as.logical(eBayes_trend),
      eBayes_robust = as.logical(eBayes_robust)
    )
  }

  # Extract the contrast name (contains "=" delimiter needed for downstream functions)
  contrast_name <- names(contrasts_results$results)[1]
  message(paste("   DEBUG66: contrast_name extracted =", contrast_name))

  # Map back to original group names in results if needed
  if (exists("group_mapping")) {
    contrasts_results_table <- contrasts_results$results |>
      dplyr::mutate(comparison = purrr::map_chr(comparison, \(x) {
        result <- x
        for (safe_name in names(group_mapping)) {
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
  raw_pval_hist_list <- purrr::imap(contrasts_results_table, \(da_tbl, contrast_name) {
    # Check if raw_pvalue column exists
    if (!"raw_pvalue" %in% colnames(da_tbl)) {
      warning(paste("raw_pvalue column not found for contrast:", contrast_name))
      return(NULL)
    }

    # Create the histogram
    p <- ggplot(da_tbl, aes(x = raw_pvalue)) +
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

  # CRITICAL FIX: getSignificantData expects list_of_da_tables to be a list where each element
  # is itself a list of data.frames. Since we're processing one contrast at a time, we need to
  # wrap the single data.frame in an additional list layer.
  # CRITICAL FIX 2: The list element name MUST contain "=" delimiter because countStatDaGenesHelper
  # expects to split on "=" to separate comparison name from expression. Use contrast_name which
  # contains the full format like "H4_vs_WT=groupH4-groupWT"
  nested_list <- list(contrasts_results_table)
  names(nested_list) <- contrast_name # Use actual contrast name with "=" delimiter
  significant_rows <- getSignificantData(
    list_of_da_tables = list(nested_list),
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
    q_val_thresh = da_q_val_thresh
  ) |>
    dplyr::rename(log2FC = "logFC")

  return_list$significant_rows <- significant_rows

  # Print the volcano plots
  volplot_plot <- plotVolcano(significant_rows,
    log_q_value_column = lqm,
    log_fc_column = log2FC,
    q_val_thresh = as.double(da_q_val_thresh),
    formula_string = "analysis_type ~ comparison"
  )


  return_list$volplot_plot <- volplot_plot

  ## Count the number of up or down significant differentially expressed proteins.
  # CRITICAL FIX: Same as above - wrap in additional list layer and use contrast_name
  nested_list_for_count <- list(contrasts_results_table)
  names(nested_list_for_count) <- contrast_name # Use actual contrast name with "=" delimiter
  num_sig_da_molecules_first_go <- printCountDaGenesTable(
    list_of_da_tables = list(nested_list_for_count),
    list_of_descriptions = list("RUV applied"),
    formula_string = "analysis_type ~ comparison"
  )

  return_list$num_sig_da_molecules_first_go <- num_sig_da_molecules_first_go

  ## Print p-values distribution figure
  pvalhist <- printPValuesDistribution(significant_rows,
    p_value_column = !!sym(raw_pvalue_column),
    formula_string = "analysis_type ~ comparison"
  )

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

  print(head(norm_counts))

  return_list$norm_counts <- norm_counts

  # Helper to handle protein ID joining based on object type
  dealWithProteinIdJoining <- function(df) {
    if (inherits(theObject, "MetaboliteAssayData")) {
      # For metabolites, we don't need to join with protein_id_table
      df
    } else if (inherits(theObject, "ProteinQuantitativeData")) {
      # For proteins, join with protein_id_table
      df |>
        left_join(theObject@protein_id_table,
          by = join_by(!!sym(args_row_id) == !!sym(theObject@protein_id_column))
        )
    } else {
      stop(sprintf("Unsupported object type: %s", class(theObject)))
    }
  }

  da_proteins_wide <- significant_rows |>
    dplyr::filter(analysis_type == "RUV applied") |>
    dplyr::select(-lqm, -colour, -analysis_type) |>
    pivot_wider(
      id_cols = c(!!sym(args_row_id)),
      names_from = c(comparison),
      names_sep = ":",
      values_from = c(log2FC, !!sym(qvalue_column), !!sym(raw_pvalue_column))
    ) |>
    left_join(counts_table_to_use, by = join_by(!!sym(args_row_id) == !!sym(args_row_id))) |>
    dealWithProteinIdJoining() |>
    dplyr::arrange(across(matches("!!sym(qvalue_column)"))) |>
    distinct()


  return_list$da_proteins_wide <- da_proteins_wide


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

  da_proteins_long <- createDaResultsLongFormat(
    lfc_qval_tbl = significant_rows |>
      dplyr::filter(analysis_type == "RUV applied"),
    norm_counts_input_tbl = getDataMatrix(theObject),
    raw_counts_input_tbl = getDataMatrix(theObject),
    row_id = args_row_id,
    sample_id = theObject@sample_id,
    group_id = group_id,
    group_pattern = args_group_pattern,
    design_matrix_norm = theObject@design_matrix,
    design_matrix_raw = theObject@design_matrix,
    protein_id_table = id_table
  )

  return_list$da_proteins_long <- da_proteins_long


  ## Plot static volcano plot
  static_volcano_plot_data <- da_proteins_long |>
    mutate(lqm = -log10(!!sym(qvalue_column))) |>
    dplyr::mutate(label = case_when(
      !!sym(qvalue_column) < da_q_val_thresh ~ "Significant",
      TRUE ~ "Not sig."
    )) |>
    dplyr::mutate(colour = case_when(
      !!sym(qvalue_column) < da_q_val_thresh ~ "purple",
      TRUE ~ "black"
    )) |>
    dplyr::mutate(colour = factor(colour, levels = c("black", "purple")))

  list_of_volcano_plots <- static_volcano_plot_data %>%
    group_by(comparison) %>%
    nest() %>%
    ungroup() %>%
    mutate(title = paste(comparison)) %>%
    mutate(plot = purrr::map2(data, title, \(x, y) {
      plotOneVolcanoNoVerticalLines(x, y,
        log_q_value_column = lqm,
        log_fc_column = log2FC
      )
    }))

  return_list$list_of_volcano_plots <- list_of_volcano_plots


  # Generate volcano plots with gene names
  cat("DEBUG_66: About to generate volcano plots with gene names\n")
  cat(sprintf("DEBUG_66: static_volcano_plot_data has %d rows\n", nrow(static_volcano_plot_data)))
  cat("DEBUG_66: static_volcano_plot_data structure:\n")
  str(static_volcano_plot_data)

  list_of_volcano_plots_with_gene_names <- tryCatch(
    {
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
          uniprot_with_gene_names <- tryCatch(
            {
              uniprot_tbl |>
                mutate(gene_name = purrr::map_chr(gene_names, \(x) {
                  tryCatch(
                    {
                      if (is.na(x) || is.null(x) || x == "") {
                        ""
                      } else {
                        split_result <- str_split(x, "; ")[[1]]
                        if (length(split_result) > 0) split_result[1] else ""
                      }
                    },
                    error = function(e) {
                      cat(sprintf("DEBUG_66: Error processing gene name: %s\n", e$message))
                      ""
                    }
                  )
                }))
            },
            error = function(e) {
              cat(sprintf("DEBUG_66: ERROR in gene name processing: %s\n", e$message))
              cat("DEBUG_66: Error details:\n")
              print(e)
              stop(e)
            }
          )

          cat("DEBUG_66: Gene names processed successfully\n")

          # Call plot function
          printOneVolcanoPlotWithProteinLabel(
            input_table = x,
            uniprot_table = uniprot_with_gene_names,
            protein_id_column = protein_id_column,
            input_title = "",
            fdr_threshold = da_q_val_thresh,
            number_of_genes = 10
          )
        }))

      cat("DEBUG_66: All plots generated successfully\n")
      with_plots
    },
    error = function(e) {
      cat(sprintf("DEBUG_66: CRITICAL ERROR in volcano plot generation: %s\n", e$message))
      cat("DEBUG_66: Full error object:\n")
      print(e)
      cat("DEBUG_66: Traceback:\n")
      print(traceback())
      NULL
    }
  )

  return_list$list_of_volcano_plots_with_gene_names <- list_of_volcano_plots_with_gene_names


  ## Return the number of significant molecules
  num_sig_da_molecules <- significant_rows %>%
    dplyr::mutate(status = case_when(
      !!sym(qvalue_column) >= da_q_val_thresh ~ "Not significant",
      log2FC >= 0 & !!sym(qvalue_column) < da_q_val_thresh ~ "Significant and Up",
      log2FC < 0 & !!sym(qvalue_column) < da_q_val_thresh ~ "Significant and Down",
      TRUE ~ "Not significant"
    )) %>%
    group_by(comparison, status) %>% # expression, analysis_type,
    summarise(counts = n()) %>%
    ungroup()

  formula_string <- ". ~ comparison"

  return_list$num_sig_da_molecules <- num_sig_da_molecules

  if (num_sig_da_molecules %>%
    dplyr::filter(status != "Not significant") |>
    nrow() > 0) {
    num_sig_da_genes_barplot_only_significant <- num_sig_da_molecules %>%
      dplyr::filter(status != "Not significant") %>%
      ggplot(aes(x = status, y = counts)) +
      geom_bar(stat = "identity") +
      geom_text(stat = "identity", aes(label = counts), vjust = -0.5) +
      theme(axis.text.x = element_text(angle = 90)) +
      facet_wrap(as.formula(formula_string))

    num_of_comparison_only_significant <- num_sig_da_molecules |>
      distinct(comparison) |>
      nrow()

    return_list$num_sig_da_genes_barplot_only_significant <- num_sig_da_genes_barplot_only_significant
    return_list$num_of_comparison_only_significant <- num_of_comparison_only_significant
  }

  if (num_sig_da_molecules %>%
    dplyr::filter(status != "Not significant") |>
    nrow() > 0) {
    num_sig_da_genes_barplot_with_not_significant <- num_sig_da_molecules %>%
      ggplot(aes(x = status, y = counts)) +
      geom_bar(stat = "identity") +
      geom_text(stat = "identity", aes(label = counts), vjust = -0.5) +
      theme(axis.text.x = element_text(angle = 90)) +
      facet_wrap(as.formula(formula_string))

    num_of_comparison_with_not_significant <- num_sig_da_molecules |>
      distinct(comparison) |>
      nrow()

    return_list$num_sig_da_genes_barplot_with_not_significant <- num_sig_da_genes_barplot_with_not_significant
    return_list$num_of_comparison_with_not_significant <- num_of_comparison_with_not_significant
  }

  return_list
}
