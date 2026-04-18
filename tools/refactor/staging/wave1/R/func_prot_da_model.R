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
prepareDataForVolcanoPlot <- function(
  input_table,
  protein_id_column = uniprot_acc,
  uniprot_table,
  uniprot_protein_id_column = uniprot_acc_first,
  gene_name_column = gene_name,
  number_of_genes = 3000,
  fdr_threshold = 0.05,
  fdr_column = q.mod,
  log2FC_column = log2FC
) {
  temp_col_name <- as_string(as_name(enquo(protein_id_column)))

  # Check if we're dealing with metabolite data (no protein IDs)
  # This handles the case where the input table is for metabolites and doesn't have Protein.Ids or uniprot_acc
  is_metabolite_data <- !any(grepl("Protein.Ids|uniprot_acc", names(input_table)))

  if (is_metabolite_data) {
    message("   Processing metabolite data (no UniProt ID mapping)...")
    proteomics_volcano_tbl <- input_table |>
      dplyr::select({{ protein_id_column }}, {{ fdr_column }}, {{ log2FC_column }}) |>
      mutate(
        colour = case_when(
          {{ fdr_column }} < fdr_threshold & {{ log2FC_column }} > 0 ~ "red",
          {{ fdr_column }} < fdr_threshold & {{ log2FC_column }} < 0 ~ "blue",
          TRUE ~ "grey"
        ),
        lqm = -log10({{ fdr_column }}),
        label = case_when(
          {{ fdr_column }} < fdr_threshold & {{ log2FC_column }} > 0 ~ "Significant Increase",
          {{ fdr_column }} < fdr_threshold & {{ log2FC_column }} < 0 ~ "Significant Decrease",
          TRUE ~ "Not significant"
        ),
        label = factor(label, levels = c(
          "Significant Increase",
          "Significant Decrease",
          "Not significant"
        )),
        rank_positive = case_when(
          {{ log2FC_column }} > 0 ~ {{ fdr_column }},
          TRUE ~ NA_real_
        ) |> rank(),
        rank_negative = case_when(
          {{ log2FC_column }} < 0 ~ {{ fdr_column }},
          TRUE ~ NA_real_
        ) |> rank(),
        # For metabolites, we use the protein_id_column (usually Name) as the "gene name"
        gene_name_significant = case_when(
          {{ fdr_column }} < fdr_threshold &
            (rank_positive <= number_of_genes |
              rank_negative <= number_of_genes) ~ as.character({{ protein_id_column }}),
          TRUE ~ NA_character_
        )
      )
  } else {
    message("   Processing protein data (performing UniProt ID mapping)...")
    proteomics_volcano_tbl <- input_table |>
      dplyr::mutate(uniprot_acc_first = purrr::map_chr({{ protein_id_column }}, \(x) {
        str_split(x, ":")[[1]][1]
      })) |>
      dplyr::relocate(uniprot_acc_first, .after = temp_col_name) |>
      dplyr::select(uniprot_acc_first, {{ fdr_column }}, {{ log2FC_column }}) |>
      left_join(uniprot_table,
        by = join_by(uniprot_acc_first == {{ uniprot_protein_id_column }})
      ) |>
      mutate(colour = case_when(
        {{ fdr_column }} < fdr_threshold & {{ log2FC_column }} > 0 ~ "red",
        {{ fdr_column }} < fdr_threshold & {{ log2FC_column }} < 0 ~ "blue",
        TRUE ~ "grey"
      )) |>
      mutate(lqm = -log10({{ fdr_column }})) |>
      mutate(label = case_when(
        {{ fdr_column }} < fdr_threshold & {{ log2FC_column }} > 0 ~ "Significant Increase",
        {{ fdr_column }} < fdr_threshold & {{ log2FC_column }} < 0 ~ "Significant Decrease",
        TRUE ~ "Not significant"
      )) |>
      mutate(label = factor(label, levels = c(
        "Significant Increase",
        "Significant Decrease",
        "Not significant"
      ))) |>
      mutate(rank_positive = case_when(
        {{ log2FC_column }} > 0 ~ {{ fdr_column }},
        TRUE ~ NA_real_
      ) |> rank()) |>
      mutate(rank_negative = case_when(
        {{ log2FC_column }} < 0 ~ {{ fdr_column }},
        TRUE ~ NA_real_
      ) |> rank()) |>
      mutate({{ gene_name_column }} := purrr::map_chr({{ gene_name_column }}, \(x) {
        str_split(x, " ")[[1]][1]
      })) |>
      mutate(gene_name_significant = case_when(
        {{ fdr_column }} < fdr_threshold &
          (rank_positive <= number_of_genes |
            rank_negative <= number_of_genes) ~ {{ gene_name_column }},
        TRUE ~ NA
      ))
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
ebFit <- function(data, design, contr.matrix) {
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

  return(list(
    table = data.frame(logFC, t.ord, t.mod, p.ord, raw_pvalue, q.ord, fdr_qvalue, df.r, df.0, s2.0, s2, s2.post),
    fit.eb = fit.eb
  ))
}

# ----------------------------------------------------------------------------
# runTest
# ----------------------------------------------------------------------------
## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#' Analyse one contrast (e.g. compare a pair of experimental groups) and output the q-values per protein.
#' @param ID List of protein accessions / row names.
#' @param A String representing the name of experimental group A for pairwise comparison of B - A.
#' @param B String representing the name of experimental group B for pairwise comparison of B - A.
#' @param group_A Names of all the columns / samples that are in experimental group A.
#' @param group_B Names of all the columns / samples that are in experimental group B.
#' @param design_matrix A data frame with a column containing the sample ID (as per the sample_id param) and the experimental group (as per the group param). Each row as the sample ID as row name in the data frame.
#' @param formula_string A formula string representing the experimental design. e.g. ("~ 0 + group")
#' @param contrast_variable String representing the contrast variable, which is also used in the formula string. (e.g. "group")
#' @param weights Numeric matrix for adjusting each sample and gene.
#' @return A data frame with the following columns:
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
#' @export
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

  contr.matrix <- makeContrasts(
    contrasts = paste0(group_B, "vs", group_A, "=", contrast_variable, group_B, "-", contrast_variable, group_A),
    levels = colnames(design_m)
  )

  eb_fit_list <- ebFit(cbind(A, B), design_m, contr.matrix = contr.matrix)

  r <- eb_fit_list$table
  fit.eb <- eb_fit_list$fit.eb

  return(list(
    table = data.frame(
      row.names = row.names(r),
      comparison = paste("log(", group_B, ") minus log(", group_A, ")", sep = ""),
      meanA = rowMeans(A),
      meanB = rowMeans(B),
      logFC = r$logFC,
      tstats = r$t.ord,
      tmod = r$t.mod,
      pval = r$p.ord,
      raw_pvalue = r$raw_pvalue,
      qval = r$q.ord,
      fdr_qvalue = r$fdr_qvalue
    ),
    fit.eb = fit.eb
  ))
}

# ----------------------------------------------------------------------------
# runTests
# ----------------------------------------------------------------------------
#' Compare a pair of experimental groups and output the log fold-change and q-values per protein.
#' @param ID List of protein accessions / row names.
#' @param data Data frame containing the log (base 2) protein abundance values where each column represents a sample and each row represents a protein group, and proteins as rows. The data is preferably median-scaled, with missing values imputed, and batch-effects removed.
#' @param test_pairs Input file with a table listing all the pairs of experimental groups to compare. First column represents group A and second column represents group B. Linear model comparisons (e.g. Contrasts) would be group B minus group A.
#' @param sample_columns A vector of column names (e.g. strings) representing samples which would be used in the statistical tests. Each column contains protein abundance values.
#' @param sample_rows_list A list, the name of each element is the sample ID and each element is a vector containing the protein accessions (e.g. row_id) with enough number of values. It is usually the output from the function \code{get_rows_to_keep_list}.
#' @param type_of_grouping A list where each element name is the name of a treatment group and each element is a vector containing the sample IDs within the treatment group. It is usually the output from the function \code{get_type_of_grouping}.
#' @param design_matrix A data frame with a column containing the sample ID (as per the sample_id param) and the experimental group (as per the group param). Each row as the sample ID as row name in the data frame.
#' @param formula_string A formula string representing the experimental design. e.g. ("~ 0 + group")
#' @param contrast_variable String representing the contrast variable, which is also used in the formula string. (e.g. "group")
#' @param weights Numeric matrix for adjusting each sample and gene.
#' @return A list of data frames, the name of each element represents each pairwise comparison. Each data frame has the following columns:
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
#' @export
runTests <- function(ID, data, test_pairs, sample_columns, sample_rows_list = NA, type_of_grouping, design_matrix, formula_string, contrast_variable = "group", weights = NA) {
  r <- list()
  for (i in 1:nrow(test_pairs)) {
    rows_to_keep <- rownames(data)


    if (length(sample_rows_list) > 0) {
      if (!is.na(sample_rows_list) &
        #  Check that sample group exists as names inside sample_rows_list
        length(which(c(test_pairs[i, "A"], test_pairs[i, "B"]) %in% names(sample_rows_list))) > 0) {
        rows_to_keep <- unique(
          sample_rows_list[[test_pairs[[i, "A"]]]],
          sample_rows_list[[test_pairs[[i, "B"]]]]
        )
      }
    }

    tmp <- data[rows_to_keep, sample_columns]
    rep <- colnames(tmp)

    # print( paste( test_pairs[i,]$A, test_pairs[i,]$B) )
    A <- tmp[, type_of_grouping[test_pairs[i, ]$A][[1]]]
    B <- tmp[, type_of_grouping[test_pairs[i, ]$B][[1]]]

    subset_weights <- NA

    if (!is.na(weights)) {
      subset_weights <- weights[c(colnames(A), colnames(B)), ]
    }

    # print(colnames(A))
    # print(colnames(B))
    tmp <- unname(cbind(A, B))
    Aname <- paste(test_pairs[i, ]$A, 1:max(1, ncol(A)), sep = "_")
    Bname <- paste(test_pairs[i, ]$B, 1:max(1, ncol(B)), sep = "_")
    colnames(tmp) <- c(Aname, Bname)

    selected_sample_ids <- c(type_of_grouping[test_pairs[i, ]$A][[1]], type_of_grouping[test_pairs[i, ]$B][[1]])
    design_matrix_subset <- design_matrix[selected_sample_ids, , drop = FALSE]

    # print("My design matrix 1")
    # print( selected_sample_ids)
    # print(design_matrix)
    # print( dim(design_matrix))

    group_A <- test_pairs[i, ]$A
    group_B <- test_pairs[i, ]$B

    x <- runTest(ID, A, B, group_A, group_B,
      design_matrix = design_matrix_subset,
      formula_string = formula_string, contrast_variable = contrast_variable,
      weights = subset_weights
    )

    comparison <- paste(group_B, " vs ", group_A, sep = "")

    r[[comparison]] <- list(results = x$table, counts = t(cbind(A, B)), fit.eb = x$fit.eb)
  }
  r
}

# ----------------------------------------------------------------------------
# runTestsContrasts
# ----------------------------------------------------------------------------
#' Run the linear model fitting and statistical tests for a set of contrasts, then adjust with Empirical Bayes function
#' @param data Data frame containing the log (base 2) protein abundance values where each column represents a sample and each row represents a protein group, and proteins as rows. The data is preferably median-scaled, with missing values imputed, and batch-effects removed.
#' @param contrast_strings Input file with a table listing all the experimental contrasts to analyse. It will be in the format required for the function \code{makeContrasts} in the limma package.
#' The contrast string consists of variable that each consist of concatenating the column name (e.g. group) and the string representing the group type (e.g. A) in the design matrix.
#' @param design_matrix A data frame with a column containing the sample ID (as per the sample_id param) and the experimental group (as per the group param). Each row as the sample ID as row name in the data frame.
#' @param formula_string A formula string representing the experimental design. e.g. ("~ 0 + group")
#' @param p_value_column The name of the raw p-value column (tidyverse style).
#' @param q_value_column The name of the q-value column (tidyverse style).
#' @param fdr_value_column The name of the fdr-value column (tidyverse style).
#' @return A list containing two elements. $results returns a list of tables containing logFC and q-values. $fit.eb returns the Empiracle Bayes output object.
#' @export
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
  message(sprintf("   runTestsContrasts: contrasts = %s", paste(contrast_strings, collapse = ", ")))
  message(sprintf("   runTestsContrasts: treat_lfc_cutoff = %s", treat_lfc_cutoff))

  # Create formula and design matrix
  ff <- as.formula(formula_string)
  mod_frame <- model.frame(ff, design_matrix)
  design_m <- model.matrix(ff, mod_frame)
  message(sprintf("   runTestsContrasts: design_m dims = %d x %d", nrow(design_m), ncol(design_m)))

  # Subset data to match design matrix
  data_subset <- data[, rownames(design_m)]
  message(sprintf("   runTestsContrasts: data_subset dims = %d x %d", nrow(data_subset), ncol(data_subset)))

  # Create contrast matrix
  message("   runTestsContrasts: Creating contrast matrix...")
  contr.matrix <- makeContrasts(
    contrasts = contrast_strings,
    levels = colnames(design_m)
  )
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

  # Check for technical replicates
  # Blocking factor should be constructed from replicate ID (biological replicate)
  # But we need to ensure it maps correctly to the samples in the subset
  # The design_matrix here is already subsetted/ordered match mod_frame in model.matrix construction
  # However, we need the original columns (group, replicates) which might be in the 'design_matrix' argument (which is the params data.frame)
  # BUT 'design_matrix' input argument is strictly the data.frame with metadata

  # Ensure we have the metadata for the subsetted samples
  # The samples in data_subset are columns matching rownames(design_m)
  # design_m was created from design_matrix (the metadata dataframe)
  samples_in_model <- rownames(design_m)

  # Check if we have replicates column
  has_replicates <- "replicates" %in% colnames(design_matrix)

  fit <- NULL

  if (has_replicates) {
    # Extract replicates for the samples in the model, maintaining order
    # We assume design_matrix rownames are the sample IDs (which was set in helper/wrapper functions)
    design_matrix_subset <- design_matrix[samples_in_model, ]

    # Construct blocking factor: Biological Replicate ID (e.g. "Control_1")
    # This identifies the biological unit that is technically replicated
    # Combining group + replicate number gives unique biological ID
    # Use paste0 to ensure character vector
    block <- paste(design_matrix_subset$group, design_matrix_subset$replicates, sep = "_")

    # Check if there are any technical replicates (duplicated blocks)
    if (any(duplicated(block))) {
      message("   runTestsContrasts: Detected technical replicates. Calculating duplicateCorrelation...")
      message(sprintf(
        "   runTestsContrasts: Block defined by group_replicates. %d unique blocks for %d samples.",
        length(unique(block)), length(block)
      ))

      # Calculate consensus correlation
      # Note: duplicateCorrelation can be slow for large datasets
      dup_cor <- duplicateCorrelation(data_subset, design = design_m, block = block)

      message(sprintf("   runTestsContrasts: Consensus correlation = %.4f", dup_cor$consensus.correlation))

      # Run lmFit with correlation and block
      message("   runTestsContrasts: Running lmFit with duplicateCorrelation...")
      fit <- lmFit(data_subset, design = design_m, block = block, correlation = dup_cor$consensus.correlation)
    } else {
      message("   runTestsContrasts: No technical replicates detected (unique blocks). Running standard lmFit...")
      fit <- lmFit(data_subset, design = design_m)
    }
  } else {
    message("   runTestsContrasts: 'replicates' column not found. Running standard lmFit...")
    fit <- lmFit(data_subset, design = design_m)
  }

  message("   runTestsContrasts: Running contrasts.fit...")
  cfit <- contrasts.fit(fit, contrasts = contr.matrix)

  message("   runTestsContrasts: Running eBayes...")
  eb.fit <- eBayes(cfit, trend = eBayes_trend, robust = eBayes_robust)

  # Run treat or standard analysis
  t.fit <- NA
  result_tables <- NA
  if (!is.na(treat_lfc_cutoff)) {
    message("   runTestsContrasts: Running treat analysis...")
    t.fit <- treat(eb.fit, lfc = as.double(treat_lfc_cutoff))

    message(sprintf("   runTestsContrasts: Processing %d contrasts with topTreat...", length(contrast_strings)))
    # Track qvalue failures for user notification
    qvalue_failures <- list()

    result_tables <- purrr::map(
      contrast_strings,
      function(contrast) {
        message(sprintf("      [map] Processing contrast: %s", contrast))
        qvalue_failed <- FALSE

        tryCatch(
          {
            message(sprintf("      About to call topTreat with coef = %s", contrast))
            da_tbl <- topTreat(t.fit, coef = contrast, n = Inf)
            message(sprintf("      [map] topTreat success: %d rows", nrow(da_tbl)))

            message("      Adding qvalue column...")
            # Safe qvalue computation: handle invalid p-values (NA, Inf, NaN)
            valid_p_idx <- which(!is.na(da_tbl$P.Value) & is.finite(da_tbl$P.Value))
            if (length(valid_p_idx) > 0) {
              # Compute q-values only for valid p-values
              valid_p_values <- da_tbl$P.Value[valid_p_idx]

              # Diagnostic: Log p-value distribution statistics
              message(sprintf("      Diagnostic: Valid p-values: %d of %d total", length(valid_p_idx), nrow(da_tbl)))
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

              q_values_all <- rep(NA_real_, nrow(da_tbl))
              if (use_qvalue) {
                tryCatch(
                  {
                    q_values_valid <- qvalue(valid_p_values)$q
                    q_values_all[valid_p_idx] <- q_values_valid
                    message("      qvalue() computation successful")
                  },
                  error = function(e) {
                    qvalue_failed <<- TRUE
                    message(sprintf("      Warning: qvalue() failed during computation: %s", e$message))
                    message(sprintf("      Diagnostic: P-value distribution may be problematic for qvalue smoothing algorithm"))
                    message(sprintf("      Diagnostic: Falling back to p.adjust() method='BH' (Benjamini-Hochberg FDR)"))
                    # Fallback to p.adjust if qvalue fails
                    q_values_all[valid_p_idx] <- p.adjust(valid_p_values, method = "BH")
                    message(sprintf("      Diagnostic: Assigned %d p.adjust() values to q-value column", length(valid_p_idx)))
                  }
                )
              } else {
                # Use p.adjust directly for edge cases
                qvalue_failed <<- TRUE
                q_values_all[valid_p_idx] <- p.adjust(valid_p_values, method = "BH")
                message("      Using p.adjust() due to edge case detection")
                message(sprintf("      Diagnostic: Assigned %d p.adjust() values to q-value column", length(valid_p_idx)))
              }
            } else {
              # All p-values are invalid, set all q-values to NA
              q_values_all <- rep(NA_real_, nrow(da_tbl))
              message("      Warning: All p-values are invalid (NA, Inf, or NaN), setting q-values to NA")
            }

            # Debug: Verify assignment before mutate
            message(sprintf("      Diagnostic: q_values_all has %d non-NA values before assignment", sum(!is.na(q_values_all))))

            da_tbl <- da_tbl |>
              mutate({{ q_value_column }} := q_values_all)

            # Debug: Verify assignment after mutate
            assigned_col <- da_tbl[[rlang::as_name(rlang::ensym(q_value_column))]]
            message(sprintf("      Diagnostic: Assigned column has %d non-NA values after mutate", sum(!is.na(assigned_col))))
            message("      qvalue column added")

            message("      Adding FDR column...")
            # Use the same safe logic for FDR column - only compute for valid p-values
            fdr_values_all <- rep(NA_real_, nrow(da_tbl))
            if (length(valid_p_idx) > 0) {
              fdr_values_all[valid_p_idx] <- p.adjust(valid_p_values, method = "BH")
            }
            da_tbl <- da_tbl |>
              mutate({{ fdr_value_column }} := fdr_values_all)
            message("      FDR column added")

            message("      Renaming P.Value column...")
            da_tbl <- da_tbl |>
              dplyr::rename({{ p_value_column }} := P.Value)
            message("      P.Value column renamed")

            message(sprintf("   [map] Completed processing contrast: %s", contrast))
            if (qvalue_failed) {
              qvalue_failures[[contrast]] <<- TRUE
            }
            return(da_tbl)
          },
          error = function(e) {
            message(sprintf("      [map] ERROR in contrast %s: %s", contrast, e$message))
            message(sprintf("      [map] ERROR call stack: %s", capture.output(traceback())))
            stop(e)
          }
        )
      }
    )
  } else {
    message("   runTestsContrasts: Running standard analysis...")
    t.fit <- eb.fit

    message(sprintf("   runTestsContrasts: Processing %d contrasts with topTable...", length(contrast_strings)))
    # Track qvalue failures for user notification
    qvalue_failures <- list()

    result_tables <- purrr::map(
      contrast_strings,
      function(contrast) {
        message(sprintf("      [map] Processing contrast: %s", contrast))
        qvalue_failed <- FALSE

        tryCatch(
          {
            message(sprintf("      About to call topTable with coef = %s", contrast))
            da_tbl <- topTable(t.fit, coef = contrast, n = Inf)
            message(sprintf("      [map] topTable success: %d rows", nrow(da_tbl)))

            message("      Adding qvalue column...")
            # Safe qvalue computation: handle invalid p-values (NA, Inf, NaN)
            valid_p_idx <- which(!is.na(da_tbl$P.Value) & is.finite(da_tbl$P.Value))
            if (length(valid_p_idx) > 0) {
              # Compute q-values only for valid p-values
              valid_p_values <- da_tbl$P.Value[valid_p_idx]

              # Diagnostic: Log p-value distribution statistics
              message(sprintf("      Diagnostic: Valid p-values: %d of %d total", length(valid_p_idx), nrow(da_tbl)))
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

              q_values_all <- rep(NA_real_, nrow(da_tbl))
              if (use_qvalue) {
                tryCatch(
                  {
                    q_values_valid <- qvalue(valid_p_values)$q
                    q_values_all[valid_p_idx] <- q_values_valid
                    message("      qvalue() computation successful")
                  },
                  error = function(e) {
                    qvalue_failed <<- TRUE
                    message(sprintf("      Warning: qvalue() failed during computation: %s", e$message))
                    message(sprintf("      Diagnostic: P-value distribution may be problematic for qvalue smoothing algorithm"))
                    message(sprintf("      Diagnostic: Falling back to p.adjust() method='BH' (Benjamini-Hochberg FDR)"))
                    # Fallback to p.adjust if qvalue fails
                    q_values_all[valid_p_idx] <- p.adjust(valid_p_values, method = "BH")
                    message(sprintf("      Diagnostic: Assigned %d p.adjust() values to q-value column", length(valid_p_idx)))
                  }
                )
              } else {
                # Use p.adjust directly for edge cases
                qvalue_failed <<- TRUE
                q_values_all[valid_p_idx] <- p.adjust(valid_p_values, method = "BH")
                message("      Using p.adjust() due to edge case detection")
                message(sprintf("      Diagnostic: Assigned %d p.adjust() values to q-value column", length(valid_p_idx)))
              }
            } else {
              # All p-values are invalid, set all q-values to NA
              q_values_all <- rep(NA_real_, nrow(da_tbl))
              message("      Warning: All p-values are invalid (NA, Inf, or NaN), setting q-values to NA")
            }

            # Debug: Verify assignment before mutate
            message(sprintf("      Diagnostic: q_values_all has %d non-NA values before assignment", sum(!is.na(q_values_all))))

            da_tbl <- da_tbl |>
              mutate({{ q_value_column }} := q_values_all)

            # Debug: Verify assignment after mutate
            assigned_col <- da_tbl[[rlang::as_name(rlang::ensym(q_value_column))]]
            message(sprintf("      Diagnostic: Assigned column has %d non-NA values after mutate", sum(!is.na(assigned_col))))
            message("      qvalue column added")

            message("      Adding FDR column...")
            # Use the same safe logic for FDR column - only compute for valid p-values
            fdr_values_all <- rep(NA_real_, nrow(da_tbl))
            if (length(valid_p_idx) > 0) {
              fdr_values_all[valid_p_idx] <- p.adjust(valid_p_values, method = "BH")
            }
            da_tbl <- da_tbl |>
              mutate({{ fdr_value_column }} := fdr_values_all)
            message("      FDR column added")

            message("      Renaming P.Value column...")
            da_tbl <- da_tbl |>
              dplyr::rename({{ p_value_column }} := P.Value)
            message("      P.Value column renamed")

            message(sprintf("   [map] Completed processing contrast: %s", contrast))
            if (qvalue_failed) {
              qvalue_failures[[contrast]] <<- TRUE
            }
            return(da_tbl)
          },
          error = function(e) {
            message(sprintf("      [map] ERROR in contrast %s: %s", contrast, e$message))
            message(sprintf("      [map] ERROR call stack: %s", capture.output(traceback())))
            stop(e)
          }
        )
      }
    )
  }

  names(result_tables) <- contrast_strings
  message("--- Exiting runTestsContrasts ---")
  return(list(results = result_tables, fit.eb = t.fit, qvalue_warnings = qvalue_failures))
}

# ----------------------------------------------------------------------------
# differentialAbundanceAnalysis
# ----------------------------------------------------------------------------
#' @export
setMethod(
  f = "differentialAbundanceAnalysis",
  signature = "ProteinQuantitativeData",
  definition = function(
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
    # IMMEDIATE ERROR CATCH - Check if we even get here
    message("*** ENTERING differentialAbundanceAnalysis METHOD ***")
    message(sprintf("*** METHOD SIGNATURE MATCHED: ProteinQuantitativeData ***"))

    # Try to catch the index error immediately
    tryCatch(
      {
        message("*** Testing parameter access ***")
        if (!is.null(contrasts_tbl)) {
          test <- contrasts_tbl[[1]]
          message("*** Parameter access successful ***")
        }
      },
      error = function(e) {
        message(sprintf("*** IMMEDIATE ERROR: %s ***", e$message))
        message(sprintf("*** ERROR CLASS: %s ***", class(e)))
        stop(e)
      }
    )

    message("--- Entering differentialAbundanceAnalysis ---")
    message(sprintf("   differentialAbundanceAnalysis: theObject class = %s", class(theObject)))
    message(sprintf("   differentialAbundanceAnalysis: contrasts_tbl provided = %s", !is.null(contrasts_tbl)))
    if (!is.null(contrasts_tbl)) {
      message(sprintf("   differentialAbundanceAnalysis: contrasts_tbl dims = %d x %d", nrow(contrasts_tbl), ncol(contrasts_tbl)))
      message(sprintf("   differentialAbundanceAnalysis: contrasts_tbl content = %s", paste(contrasts_tbl[[1]], collapse = ", ")))
    }

    # Wrap the helper function call in tryCatch to get better error info
    message("   differentialAbundanceAnalysis: About to call differentialAbundanceAnalysisHelper...")

    results_list <- tryCatch(
      {
        differentialAbundanceAnalysisHelper(theObject,
          contrasts_tbl = contrasts_tbl,
          formula_string = formula_string,
          group_id = group_id,
          da_q_val_thresh = da_q_val_thresh,
          treat_lfc_cutoff = treat_lfc_cutoff,
          eBayes_trend = eBayes_trend,
          eBayes_robust = eBayes_robust,
          args_group_pattern = args_group_pattern,
          args_row_id = args_row_id,
          qvalue_column = qvalue_column,
          raw_pvalue_column = raw_pvalue_column
        )
      },
      error = function(e) {
        # CRITICAL FIX: Use paste() for logger calls in error handlers to avoid interpolation bug
        message(paste("   differentialAbundanceAnalysis ERROR in helper function:", e$message))
        message(paste("   differentialAbundanceAnalysis ERROR call stack:", capture.output(traceback())))
        stop(e)
      }
    )

    message("   differentialAbundanceAnalysis: Helper function completed successfully!")
    message("--- Exiting differentialAbundanceAnalysis ---")
    return(results_list)
  }
)

# ----------------------------------------------------------------------------
# differentialAbundanceAnalysisHelper
# ----------------------------------------------------------------------------
#' @export
setMethod(
  f = "differentialAbundanceAnalysisHelper",
  signature = "ProteinQuantitativeData",
  definition = function(
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
    message("--- Entering differentialAbundanceAnalysisHelper ---")

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
      tryCatch(
        {
          first_col <- contrasts_tbl[[1]]
          message(sprintf("      SUCCESS: contrasts_tbl[[1]] class: %s", class(first_col)))
          message(sprintf("      SUCCESS: contrasts_tbl[[1]] length: %d", length(first_col)))
          message(sprintf("      SUCCESS: contrasts_tbl[[1]] content: %s", paste(first_col, collapse = ", ")))
        },
        error = function(e) {
          message(sprintf("      ERROR accessing contrasts_tbl[[1]]: %s", e$message))
        }
      )
    }

    message("   DEBUG66: Initial parameter inspection")
    message(sprintf("      DEBUG66: theObject class = %s", class(theObject)))
    message(sprintf("      DEBUG66: contrasts_tbl param is.null = %s", is.null(contrasts_tbl)))
    if (!is.null(contrasts_tbl)) {
      message(sprintf("      DEBUG66: contrasts_tbl param class = %s", class(contrasts_tbl)))
      message("      DEBUG66: contrasts_tbl param structure:")
      str(contrasts_tbl)
    }

    message("   differentialAbundanceAnalysisHelper Step: Extracting parameters...")
    # Extract parameters from S4 object with fallbacks
    contrasts_tbl <- checkParamsObjectFunctionSimplify(theObject, "contrasts_tbl", contrasts_tbl)
    formula_string <- checkParamsObjectFunctionSimplify(theObject, "formula_string", "~ 0 + group")
    group_id <- checkParamsObjectFunctionSimplify(theObject, "group_id", "group")
    da_q_val_thresh <- checkParamsObjectFunctionSimplify(theObject, "da_q_val_thresh", 0.05)
    treat_lfc_cutoff <- checkParamsObjectFunctionSimplify(theObject, "treat_lfc_cutoff", 0)
    eBayes_trend <- checkParamsObjectFunctionSimplify(theObject, "eBayes_trend", TRUE)
    eBayes_robust <- checkParamsObjectFunctionSimplify(theObject, "eBayes_robust", TRUE)
    args_group_pattern <- checkParamsObjectFunctionSimplify(theObject, "args_group_pattern", "(\\d+)")
    args_row_id <- checkParamsObjectFunctionSimplify(theObject, "args_row_id", "uniprot_acc")

    message("   DEBUG66: After parameter extraction")
    message(sprintf("      DEBUG66: contrasts_tbl class = %s", class(contrasts_tbl)))
    message("      DEBUG66: contrasts_tbl structure:")
    str(contrasts_tbl)

    message(sprintf("   differentialAbundanceAnalysisHelper: formula_string = %s", formula_string))
    message(sprintf("   differentialAbundanceAnalysisHelper: group_id = %s", group_id))
    message(sprintf("   differentialAbundanceAnalysisHelper: da_q_val_thresh = %f", da_q_val_thresh))

    message("   differentialAbundanceAnalysisHelper Step: Processing group names...")
    # Handle group names that start with numbers (same pattern as original wrapper)
    design_matrix <- theObject@design_matrix
    group_col <- design_matrix[[group_id]]
    message(sprintf("   differentialAbundanceAnalysisHelper: group_col length = %d", length(group_col)))
    message(sprintf("   differentialAbundanceAnalysisHelper: group_col content = %s", paste(head(group_col, 10), collapse = ", ")))

    # Check if any group names start with numbers and create mapping
    starts_with_number <- grepl("^[0-9]", group_col)
    message(sprintf("   differentialAbundanceAnalysisHelper: any start with number? %s", any(starts_with_number)))

    if (any(starts_with_number)) {
      message("   differentialAbundanceAnalysisHelper Step: Fixing group names that start with numbers...")

      original_groups <- unique(group_col)
      message(sprintf("   differentialAbundanceAnalysisHelper: original_groups = %s", paste(original_groups, collapse = ", ")))
      message(sprintf("   differentialAbundanceAnalysisHelper: About to process %d original groups with purrr::map_chr", length(original_groups)))

      tryCatch(
        {
          message("      DEBUG66: Before safe_groups purrr::map_chr")
          message(sprintf("      DEBUG66: original_groups class = %s", class(original_groups)))
          message(sprintf("      DEBUG66: original_groups length = %d", length(original_groups)))
          message("      DEBUG66: original_groups content:")
          print(original_groups)

          safe_groups <- purrr::map_chr(original_groups, \(x) {
            message(sprintf("         DEBUG66: Processing group item: '%s' (class: %s)", x, class(x)))
            result <- if (grepl("^[0-9]", x)) paste0("grp_", x) else x
            message(sprintf("         DEBUG66: Result for '%s' -> '%s'", x, result))
            result
          })
          message("   differentialAbundanceAnalysisHelper: safe_groups processing SUCCESS")
          message("      DEBUG66: safe_groups result:")
          print(safe_groups)
        },
        error = function(e) {
          message(sprintf("   differentialAbundanceAnalysisHelper ERROR in safe_groups purrr::map_chr: %s", e$message))
          message("      DEBUG66: Error details:")
          message(sprintf("      DEBUG66: Error class: %s", class(e)))
          print(str(e))
          stop(e)
        }
      )

      group_mapping <- setNames(original_groups, safe_groups)
      message("   differentialAbundanceAnalysisHelper: group_mapping created")

      # Update design matrix with safe names
      message("   differentialAbundanceAnalysisHelper: About to update design matrix with purrr::map_chr")
      tryCatch(
        {
          design_matrix[[group_id]] <- purrr::map_chr(group_col, \(x) {
            if (grepl("^[0-9]", x)) paste0("grp_", x) else x
          })
          message("   differentialAbundanceAnalysisHelper: design matrix update SUCCESS")
        },
        error = function(e) {
          message(sprintf("   differentialAbundanceAnalysisHelper ERROR in design matrix purrr::map_chr: %s", e$message))
          stop(e)
        }
      )

      # Update contrasts table if it exists
      if (!is.null(contrasts_tbl)) {
        message("   differentialAbundanceAnalysisHelper: About to update contrasts table with purrr::map_chr")
        message(sprintf("   differentialAbundanceAnalysisHelper: contrasts_tbl[[1]] = %s", paste(contrasts_tbl[[1]], collapse = ", ")))

        message("      DEBUG66: Contrasts table inspection before purrr::map_chr")
        message(sprintf("      DEBUG66: contrasts_tbl class = %s", class(contrasts_tbl)))
        message(sprintf("      DEBUG66: contrasts_tbl dim = %d x %d", nrow(contrasts_tbl), ncol(contrasts_tbl)))
        message("      DEBUG66: contrasts_tbl structure:")
        str(contrasts_tbl)
        message(sprintf("      DEBUG66: contrasts_tbl[[1]] class = %s", class(contrasts_tbl[[1]])))
        message(sprintf("      DEBUG66: contrasts_tbl[[1]] length = %d", length(contrasts_tbl[[1]])))
        message("      DEBUG66: group_mapping content:")
        print(group_mapping)

        tryCatch(
          {
            contrasts_tbl[[1]] <- purrr::map_chr(contrasts_tbl[[1]], \(x) {
              message(sprintf("         DEBUG66: Processing contrast: '%s' (class: %s)", x, class(x)))
              message(sprintf("         DEBUG66: group_mapping names: %s", paste(names(group_mapping), collapse = ", ")))
              result <- x
              for (orig in names(group_mapping)) {
                message(sprintf("            DEBUG66: Checking gsub: '%s' -> '%s'", group_mapping[orig], orig))
                result <- gsub(group_mapping[orig], orig, result, fixed = TRUE)
              }
              message(sprintf("         DEBUG66: Final result: '%s' -> '%s'", x, result))
              result
            })
            message("   differentialAbundanceAnalysisHelper: contrasts table update SUCCESS")
            message("      DEBUG66: Updated contrasts_tbl[[1]]:")
            print(contrasts_tbl[[1]])
          },
          error = function(e) {
            message(sprintf("   differentialAbundanceAnalysisHelper ERROR in contrasts table purrr::map_chr: %s", e$message))
            message("      DEBUG66: Error details:")
            message(sprintf("      DEBUG66: Error class: %s", class(e)))
            print(str(e))
            stop(e)
          }
        )
      }

      theObject@design_matrix <- design_matrix
    }

    message("   differentialAbundanceAnalysisHelper Step: Updating S4 object parameters...")
    # Update object with validated parameters
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

    message("   differentialAbundanceAnalysisHelper Step: Generating QC plots...")

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

    # Count the number of values
    return_list$plot_num_of_values <- plotNumOfValuesNoLog(theObject@protein_quant_table)

    message("   differentialAbundanceAnalysisHelper Step: Running limma contrasts analysis...")

    # Prepare design matrix for limma
    message("   DEBUG66: About to prepare design matrix rownames")
    message(paste("      DEBUG66: theObject@sample_id =", theObject@sample_id))
    message(paste("      DEBUG66: design_matrix dims before =", nrow(theObject@design_matrix), "x", ncol(theObject@design_matrix)))
    message("      DEBUG66: design_matrix column names:")
    print(colnames(theObject@design_matrix))

    rownames(theObject@design_matrix) <- theObject@design_matrix |> dplyr::pull(one_of(theObject@sample_id))

    message("      DEBUG66: design_matrix rownames set successfully")

    # Run the core limma analysis using existing function
    protein_quant_matrix <- as.matrix(column_to_rownames(theObject@protein_quant_table, theObject@protein_id_column))
    
    # DEBUG: Check if rownames are blank
    if (any(rownames(protein_quant_matrix) == "") || any(is.na(rownames(protein_quant_matrix)))) {
      message("   WARNING DEBUG66: protein_quant_matrix has blank or NA rownames!")
      message(paste("      Sample of rownames:", paste(head(rownames(protein_quant_matrix), 5), collapse = ", ")))
    }
    
    contrast_strings_to_use <- contrasts_tbl[, 1]
    message(paste("   differentialAbundanceAnalysisHelper: About to call runTestsContrasts with", length(contrast_strings_to_use), "contrasts"))
    message(paste("   differentialAbundanceAnalysisHelper: protein_quant_matrix dims =", nrow(protein_quant_matrix), "x", ncol(protein_quant_matrix)))

    contrasts_results <- runTestsContrasts(protein_quant_matrix,
      contrast_strings = contrast_strings_to_use,
      design_matrix = theObject@design_matrix,
      formula_string = formula_string,
      weights = NA,
      treat_lfc_cutoff = as.double(treat_lfc_cutoff),
      eBayes_trend = as.logical(eBayes_trend),
      eBayes_robust = as.logical(eBayes_robust)
    )

    message("   differentialAbundanceAnalysisHelper Step: runTestsContrasts completed successfully!")

    # The result from runTestsContrasts is a LIST with a 'results'
    # element, which ITSELF is a list of tables. For a single contrast run, we
    # need to extract the first table from the nested list.

    # Check if the expected structure exists
    if (is.null(contrasts_results) || is.null(contrasts_results$results) || length(contrasts_results$results) == 0) {
      stop("Error: DE analysis function did not return results in the expected format.")
    }

    # Extract the results table (it's the first element in the nested list)
    da_results_table <- contrasts_results$results[[1]]

    # Get the name of the contrast from the list element name
    contrast_name <- names(contrasts_results$results)[1]

    # Ensure it's a data frame before proceeding
    if (!is.data.frame(da_results_table)) {
      stop("Error: DE analysis results are not in the expected data frame format.")
    }

    # CRITICAL FIX 3.0: The 'topTreat' table needs a 'comparison' column added to it,
    # containing the name of the contrast. It also needs the protein IDs from rownames.
    da_results_table <- da_results_table |>
      tibble::rownames_to_column(var = args_row_id) |>
      dplyr::mutate(comparison = contrast_name)

    # Map back to original group names in results if needed
    if (exists("group_mapping")) {
      contrasts_results_table <- da_results_table |>
        dplyr::mutate(comparison = purrr::map_chr(comparison, \(x) {
          result <- x
          for (safe_name in names(group_mapping)) {
            result <- gsub(safe_name, group_mapping[safe_name], result, fixed = TRUE)
          }
          result
        }))
    } else {
      contrasts_results_table <- da_results_table
    }

    return_list$contrasts_results <- contrasts_results
    return_list$contrasts_results_table <- contrasts_results_table

    # Extract qvalue warnings if present
    if (!is.null(contrasts_results$qvalue_warnings) && length(contrasts_results$qvalue_warnings) > 0) {
      return_list$qvalue_warnings <- contrasts_results$qvalue_warnings
      message(sprintf("   differentialAbundanceAnalysisHelper Step: qvalue() failed for %d contrast(s)", length(contrasts_results$qvalue_warnings)))
    }

    message("   differentialAbundanceAnalysisHelper Step: Preparing data for visualization...")
    message(paste("   DEBUG66: contrast_name for list naming =", contrast_name))

    # Prepare data for volcano plots
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

    # Generate volcano plots
    volplot_plot <- plotVolcano(significant_rows,
      log_q_value_column = lqm,
      log_fc_column = log2FC,
      q_val_thresh = da_q_val_thresh,
      formula_string = "analysis_type ~ comparison"
    )

    return_list$volplot_plot <- volplot_plot

    # Count significant molecules
    # CRITICAL FIX: Same as above - wrap in additional list layer and use contrast_name
    nested_list_for_count <- list(contrasts_results_table)
    names(nested_list_for_count) <- contrast_name # Use actual contrast name with "=" delimiter
    num_sig_da_molecules <- printCountDaGenesTable(
      list_of_da_tables = list(nested_list_for_count),
      list_of_descriptions = list("RUV applied"),
      formula_string = "analysis_type ~ comparison"
    )

    return_list$num_sig_da_molecules_first_go <- num_sig_da_molecules

    # P-values distribution plot
    pvalhist <- printPValuesDistribution(significant_rows,
      p_value_column = !!sym(raw_pvalue_column),
      formula_string = "analysis_type ~ comparison"
    )

    return_list$pvalhist <- pvalhist

    message("   differentialAbundanceAnalysisHelper Step: Creating results tables...")

    # Create wide format output
    counts_table_to_use <- theObject@protein_quant_table

    norm_counts <- counts_table_to_use |>
      as.data.frame() |>
      column_to_rownames(args_row_id) |>
      set_colnames(paste0(colnames(counts_table_to_use[-1]), ".log2norm")) |>
      rownames_to_column(args_row_id)

    return_list$norm_counts <- norm_counts

    da_proteins_wide <- significant_rows |>
      dplyr::filter(analysis_type == "RUV applied") |>
      dplyr::select(-lqm, -colour, -analysis_type) |>
      pivot_wider(
        id_cols = c(!!sym(args_row_id)),
        names_from = c(comparison),
        names_sep = ":",
        values_from = c(log2FC, !!sym(qvalue_column), !!sym(raw_pvalue_column))
      ) |>
      left_join(counts_table_to_use, by = join_by(!!sym(args_row_id) == !!sym(theObject@protein_id_column))) |>
      left_join(theObject@protein_id_table, by = join_by(!!sym(args_row_id) == !!sym(theObject@protein_id_column))) |>
      dplyr::arrange(across(matches(paste0("!!sym(", qvalue_column, ")")))) |>
      distinct()

    return_list$da_proteins_wide <- da_proteins_wide

    # Create long format output
    da_proteins_long <- createDaResultsLongFormat(
      lfc_qval_tbl = significant_rows |>
        dplyr::filter(analysis_type == "RUV applied"),
      norm_counts_input_tbl = as.matrix(column_to_rownames(theObject@protein_quant_table, theObject@protein_id_column)),
      raw_counts_input_tbl = as.matrix(column_to_rownames(theObject@protein_quant_table, theObject@protein_id_column)),
      row_id = args_row_id,
      sample_id = theObject@sample_id,
      group_id = group_id,
      group_pattern = args_group_pattern,
      design_matrix_norm = theObject@design_matrix,
      design_matrix_raw = theObject@design_matrix,
      protein_id_table = theObject@protein_id_table
    )

    return_list$da_proteins_long <- da_proteins_long

    # Static volcano plots with gene names
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

    # Additional summary statistics
    num_sig_da_molecules <- significant_rows %>%
      dplyr::mutate(status = case_when(
        !!sym(qvalue_column) >= da_q_val_thresh ~ "Not significant",
        log2FC >= 0 & !!sym(qvalue_column) < da_q_val_thresh ~ "Significant and Up",
        log2FC < 0 & !!sym(qvalue_column) < da_q_val_thresh ~ "Significant and Down",
        TRUE ~ "Not significant"
      )) %>%
      group_by(comparison, status) %>%
      summarise(counts = n()) %>%
      ungroup()

    return_list$num_sig_da_molecules <- num_sig_da_molecules

    # Create barplots if significant results exist
    if (num_sig_da_molecules %>%
      dplyr::filter(status != "Not significant") |>
      nrow() > 0) {
      num_sig_da_genes_barplot_only_significant <- num_sig_da_molecules %>%
        dplyr::filter(status != "Not significant") %>%
        ggplot(aes(x = status, y = counts)) +
        geom_bar(stat = "identity") +
        geom_text(stat = "identity", aes(label = counts), vjust = -0.5) +
        theme(axis.text.x = element_text(angle = 90)) +
        facet_wrap(~comparison)

      return_list$num_sig_da_genes_barplot_only_significant <- num_sig_da_genes_barplot_only_significant

      num_sig_da_genes_barplot_with_not_significant <- num_sig_da_molecules %>%
        ggplot(aes(x = status, y = counts)) +
        geom_bar(stat = "identity") +
        geom_text(stat = "identity", aes(label = counts), vjust = -0.5) +
        theme(axis.text.x = element_text(angle = 90)) +
        facet_wrap(~comparison)

      return_list$num_sig_da_genes_barplot_with_not_significant <- num_sig_da_genes_barplot_with_not_significant
    }

    message("--- Exiting differentialAbundanceAnalysisHelper ---")
    return(return_list)
  }
)

