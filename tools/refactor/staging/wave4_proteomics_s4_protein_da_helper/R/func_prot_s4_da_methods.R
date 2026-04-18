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
    # CRITICAL FIX: getSignificantData expects list_of_de_tables to be a list where each element
    # is itself a list of data.frames. Since we're processing one contrast at a time, we need to
    # wrap the single data.frame in an additional list layer.
    # CRITICAL FIX 2: The list element name MUST contain "=" delimiter because countStatDeGenesHelper
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
      q_val_thresh = as.double(da_q_val_thresh),
      formula_string = "analysis_type ~ comparison"
    )

    return_list$volplot_plot <- volplot_plot

    # Count significant molecules
    # CRITICAL FIX: Same as above - wrap in additional list layer and use contrast_name
    nested_list_for_count <- list(contrasts_results_table)
    names(nested_list_for_count) <- contrast_name # Use actual contrast name with "=" delimiter
    num_sig_de_molecules <- printCountDaGenesTable(
      list_of_da_tables = list(nested_list_for_count),
      list_of_descriptions = list("RUV applied"),
      formula_string = "analysis_type ~ comparison"
    )

    return_list$num_sig_de_molecules_first_go <- num_sig_de_molecules

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
    num_sig_de_molecules <- significant_rows %>%
      dplyr::mutate(status = case_when(
        !!sym(qvalue_column) >= da_q_val_thresh ~ "Not significant",
        log2FC >= 0 & !!sym(qvalue_column) < da_q_val_thresh ~ "Significant and Up",
        log2FC < 0 & !!sym(qvalue_column) < da_q_val_thresh ~ "Significant and Down",
        TRUE ~ "Not significant"
      )) %>%
      group_by(comparison, status) %>%
      summarise(counts = n()) %>%
      ungroup()

    return_list$num_sig_de_molecules <- num_sig_de_molecules

    # Create barplots if significant results exist
    if (num_sig_de_molecules %>%
      dplyr::filter(status != "Not significant") |>
      nrow() > 0) {
      num_sig_da_genes_barplot_only_significant <- num_sig_de_molecules %>%
        dplyr::filter(status != "Not significant") %>%
        ggplot(aes(x = status, y = counts)) +
        geom_bar(stat = "identity") +
        geom_text(stat = "identity", aes(label = counts), vjust = -0.5) +
        theme(axis.text.x = element_text(angle = 90)) +
        facet_wrap(~comparison)

      return_list$num_sig_da_genes_barplot_only_significant <- num_sig_da_genes_barplot_only_significant

      num_sig_da_genes_barplot_with_not_significant <- num_sig_de_molecules %>%
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

