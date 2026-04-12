# outputDaResultsAllContrasts
# ----------------------------------------------------------------------------
#' @export
setMethod(
  f = "outputDaResultsAllContrasts",
  signature = "ProteinQuantitativeData",
  definition = function(theObject,
                        da_results_list_all_contrasts = NULL,
                        uniprot_tbl = NULL,
                        da_output_dir = NULL,
                        publication_graphs_dir = NULL,
                        file_prefix = "da_proteins",
                        args_row_id = NULL,
                        gene_names_column = "gene_names",
                        uniprot_id_column = "Entry") {
    message("--- Entering outputDaResultsAllContrasts ---")
    message(sprintf("   outputDaResultsAllContrasts: da_output_dir = %s", da_output_dir))
    message(sprintf("   outputDaResultsAllContrasts: file_prefix = %s", file_prefix))

    # Extract parameters from S4 object with fallbacks
    uniprot_tbl <- checkParamsObjectFunctionSimplify(theObject, "uniprot_tbl", uniprot_tbl)
    da_output_dir <- checkParamsObjectFunctionSimplify(theObject, "da_output_dir", da_output_dir)
    publication_graphs_dir <- checkParamsObjectFunctionSimplify(theObject, "publication_graphs_dir", publication_graphs_dir)
    args_row_id <- checkParamsObjectFunctionSimplify(theObject, "args_row_id", args_row_id)
    gene_names_column <- checkParamsObjectFunctionSimplify(theObject, "gene_names_column", gene_names_column)
    uniprot_id_column <- checkParamsObjectFunctionSimplify(theObject, "uniprot_id_column", uniprot_id_column)
 
    # CRITICAL FIX: Normalize paths for Windows (fixes C:// double slash issue)
    # Only normalize if the path is not NULL and doesn't already exist
    if (!is.null(da_output_dir)) {
      # Use normalizePath with mustWork=FALSE to handle non-existent dirs
      da_output_dir <- gsub("//+", "/", da_output_dir) # Remove double slashes first
      da_output_dir <- normalizePath(da_output_dir, winslash = "/", mustWork = FALSE)
      message(sprintf("   outputDaResultsAllContrasts: Normalized da_output_dir = %s", da_output_dir))
    }

    if (!is.null(publication_graphs_dir)) {
      publication_graphs_dir <- gsub("//+", "/", publication_graphs_dir) # Remove double slashes first
      publication_graphs_dir <- normalizePath(publication_graphs_dir, winslash = "/", mustWork = FALSE)
      message(sprintf("   outputDaResultsAllContrasts: Normalized publication_graphs_dir = %s", publication_graphs_dir))
    }

    # Ensure output directory exists (CRITICAL FIX!)
    if (!dir.exists(da_output_dir)) {
      dir.create(da_output_dir, recursive = TRUE, showWarnings = FALSE)
      message(sprintf("   outputDaResultsAllContrasts: Created output directory: %s", da_output_dir))
    }

    message(sprintf("   outputDaResultsAllContrasts: Processing %d contrasts", length(da_results_list_all_contrasts)))

    # Write results for each contrast with UNIQUE filenames
    for (contrast_name in names(da_results_list_all_contrasts)) {
      message(sprintf("   outputDaResultsAllContrasts: Writing files for contrast: %s", contrast_name))

      contrast_result <- da_results_list_all_contrasts[[contrast_name]]

      if (!is.null(contrast_result$da_proteins_long)) {
        # Clean contrast name for safe filename
        safe_contrast_name <- gsub("[^A-Za-z0-9_-]", "_", contrast_name)

        # Create annotated version (FIXED: proper conditional logic)
        da_proteins_long_annot <- contrast_result$da_proteins_long |>
          mutate(uniprot_acc_cleaned = str_split(!!sym(args_row_id), "-") |>
            purrr::map_chr(1))

        # Add UniProt annotations if available
        if (!is.null(uniprot_tbl) && nrow(uniprot_tbl) > 0) {
          da_proteins_long_annot <- da_proteins_long_annot |>
            left_join(uniprot_tbl, by = join_by(uniprot_acc_cleaned == !!sym(uniprot_id_column))) |>
            dplyr::select(-uniprot_acc_cleaned)
        } else {
          da_proteins_long_annot <- da_proteins_long_annot |>
            dplyr::select(-uniprot_acc_cleaned)
        }

        # Add gene_name column with proper conditional logic
        if ("gene_names" %in% names(da_proteins_long_annot)) {
          da_proteins_long_annot <- da_proteins_long_annot |>
            mutate(gene_name = purrr::map_chr(gene_names, \(x){
              tryCatch(
                {
                  if (is.na(x) || is.null(x) || x == "") {
                    ""
                  } else {
                    split_result <- str_split(x, " |:")[[1]]
                    if (length(split_result) > 0) split_result[1] else ""
                  }
                },
                error = function(e) ""
              )
            }))
        } else if (gene_names_column %in% names(da_proteins_long_annot)) {
          da_proteins_long_annot <- da_proteins_long_annot |>
            mutate(gene_name = purrr::map_chr(.data[[gene_names_column]], \(x){
              tryCatch(
                {
                  if (is.na(x) || is.null(x) || x == "") {
                    ""
                  } else {
                    split_result <- str_split(x, " |:")[[1]]
                    if (length(split_result) > 0) split_result[1] else ""
                  }
                },
                error = function(e) ""
              )
            }))
        } else {
          da_proteins_long_annot <- da_proteins_long_annot |>
            mutate(gene_name = "")
        }

        # Relocate gene_name column
        da_proteins_long_annot <- da_proteins_long_annot |>
          relocate(gene_name, .after = !!sym(args_row_id))

        # CRITICAL FIX: Use contrast-specific filename!
        contrast_filename <- paste0(file_prefix, "_", safe_contrast_name, "_long_annot.tsv")
        output_path <- file.path(da_output_dir, contrast_filename)

        message(sprintf("   outputDaResultsAllContrasts: Writing %s", output_path))

        # Write the file
        vroom::vroom_write(da_proteins_long_annot, output_path)

        # Verify file was written
        if (file.exists(output_path)) {
          file_size <- file.size(output_path)
          message(sprintf(
            "   outputDaResultsAllContrasts: SUCCESS - File written: %s (%d bytes)",
            contrast_filename, file_size
          ))
        } else {
          message(sprintf("   outputDaResultsAllContrasts: ERROR - File NOT written: %s", output_path))
        }

        # Also write Excel version
        excel_filename <- paste0(file_prefix, "_", safe_contrast_name, "_long_annot.xlsx")
        excel_path <- file.path(da_output_dir, excel_filename)
        writexl::write_xlsx(da_proteins_long_annot, excel_path)

        message(sprintf("   outputDaResultsAllContrasts: Also wrote Excel: %s", excel_filename))

        # [OK] NEW: Generate volcano plot for this contrast
        message(sprintf("   outputDaResultsAllContrasts: Generating volcano plot for contrast: %s", contrast_name))

        # Extract parameters - try direct access first, then use checkParams
        da_q_val_thresh <- if (!is.null(theObject@args$outputDaResultsAllContrasts$da_q_val_thresh)) {
          theObject@args$outputDaResultsAllContrasts$da_q_val_thresh
        } else {
          checkParamsObjectFunctionSimplify(theObject, "da_q_val_thresh", 0.05)
        }

        fdr_column <- if (!is.null(theObject@args$outputDaResultsAllContrasts$fdr_column)) {
          theObject@args$outputDaResultsAllContrasts$fdr_column
        } else {
          checkParamsObjectFunctionSimplify(theObject, "fdr_column", "fdr_qvalue")
        }

        log2fc_column <- if (!is.null(theObject@args$outputDaResultsAllContrasts$log2fc_column)) {
          theObject@args$outputDaResultsAllContrasts$log2fc_column
        } else {
          checkParamsObjectFunctionSimplify(theObject, "log2fc_column", "log2FC")
        }

        # Prepare data for volcano plot (same logic as in differentialAbundanceAnalysisHelper)
        volcano_data <- contrast_result$da_proteins_long |>
          mutate(lqm = -log10(!!sym(fdr_column))) |>
          dplyr::mutate(label = case_when(
            !!sym(fdr_column) < da_q_val_thresh ~ "Significant",
            TRUE ~ "Not sig."
          )) |>
          dplyr::mutate(colour = case_when(
            !!sym(fdr_column) < da_q_val_thresh ~ "purple",
            TRUE ~ "black"
          )) |>
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
          message(sprintf("   outputDaResultsAllContrasts: Created volcano plots directory: %s", volcano_dir))
        }

        # Save volcano plot with contrast-specific filename
        volcano_png_path <- file.path(volcano_dir, paste0(safe_contrast_name, ".png"))
        volcano_pdf_path <- file.path(volcano_dir, paste0(safe_contrast_name, ".pdf"))

        # Save as PNG
        tryCatch(
          {
            ggplot2::ggsave(volcano_png_path, volcano_plot, width = 7, height = 7, dpi = 300)
            message(sprintf("   outputDaResultsAllContrasts: Saved volcano plot PNG: %s", basename(volcano_png_path)))
          },
          error = function(e) {
            message(sprintf("   outputDaResultsAllContrasts: ERROR saving PNG: %s", e$message))
          }
        )

        # Save as PDF
        tryCatch(
          {
            ggplot2::ggsave(volcano_pdf_path, volcano_plot, width = 7, height = 7)
            message(sprintf("   outputDaResultsAllContrasts: Saved volcano plot PDF: %s", basename(volcano_pdf_path)))
          },
          error = function(e) {
            message(sprintf("   outputDaResultsAllContrasts: ERROR saving PDF: %s", e$message))
          }
        )
      }
    }

    # [OK] NEW: Generate combined multi-page PDF with all volcano plots
    message("   outputDaResultsAllContrasts: Creating combined volcano plots PDF...")

    volcano_dir <- file.path(publication_graphs_dir, "Volcano_Plots")
    combined_pdf_path <- file.path(volcano_dir, "all_volcano_plots_combined.pdf")

    tryCatch(
      {
        # Re-extract parameters for combined PDF generation (in case they were modified)
        da_q_val_thresh <- if (!is.null(theObject@args$outputDaResultsAllContrasts$da_q_val_thresh)) {
          theObject@args$outputDaResultsAllContrasts$da_q_val_thresh
        } else {
          0.05
        }

        fdr_column <- if (!is.null(theObject@args$outputDaResultsAllContrasts$fdr_column)) {
          theObject@args$outputDaResultsAllContrasts$fdr_column
        } else {
          "fdr_qvalue"
        }

        log2fc_column <- if (!is.null(theObject@args$outputDaResultsAllContrasts$log2fc_column)) {
          theObject@args$outputDaResultsAllContrasts$log2fc_column
        } else {
          "log2FC"
        }

        # Collect all volcano plots
        all_volcano_plots <- list()

        for (contrast_name in names(da_results_list_all_contrasts)) {
          contrast_result <- da_results_list_all_contrasts[[contrast_name]]

          if (!is.null(contrast_result$da_proteins_long)) {
            # Recreate the volcano plot
            volcano_data <- contrast_result$da_proteins_long |>
              mutate(lqm = -log10(!!sym(fdr_column))) |>
              dplyr::mutate(label = case_when(
                !!sym(fdr_column) < da_q_val_thresh ~ "Significant",
                TRUE ~ "Not sig."
              )) |>
              dplyr::mutate(colour = case_when(
                !!sym(fdr_column) < da_q_val_thresh ~ "purple",
                TRUE ~ "black"
              )) |>
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
          message(sprintf(
            "   outputDaResultsAllContrasts: Created combined PDF with %d volcano plots: %s",
            length(all_volcano_plots), basename(combined_pdf_path)
          ))
        }
      },
      error = function(e) {
        message(sprintf("   outputDaResultsAllContrasts: ERROR creating combined PDF: %s", e$message))
      }
    )

    # ============================================================================
    # NEW: Aggregate and save Num Sig DE Molecules from all contrasts
    # ============================================================================
    message("   outputDaResultsAllContrasts: Processing NumSigDaMolecules figures...")

    tryCatch(
      {
        # Create NumSigDaMolecules directory
        numsigde_dir <- file.path(publication_graphs_dir, "NumSigDaMolecules")
        if (!dir.exists(numsigde_dir)) {
          dir.create(numsigde_dir, recursive = TRUE, showWarnings = TRUE)
          if (!dir.exists(numsigde_dir)) {
            stop("Failed to create NumSigDaMolecules directory: ", numsigde_dir)
          }
          message(sprintf("   outputDaResultsAllContrasts: Created NumSigDaMolecules directory: %s", numsigde_dir))
        }

        # Collect num_sig_da_molecules tables from all contrasts
        all_numsig_tables <- list()

        for (contrast_name in names(da_results_list_all_contrasts)) {
          contrast_result <- da_results_list_all_contrasts[[contrast_name]]

          # Check for num_sig_da_molecules_first_go in the result
          if (!is.null(contrast_result$num_sig_da_molecules_first_go)) {
            numsig_data <- contrast_result$num_sig_da_molecules_first_go

            # Extract the table component
            if (!is.null(numsig_data$table)) {
              all_numsig_tables[[contrast_name]] <- numsig_data$table
              message(sprintf("   outputDaResultsAllContrasts: Found NumSigDE table for contrast: %s", contrast_name))
            }
          }
        }

        # If we have any tables, aggregate and save them
        if (length(all_numsig_tables) > 0) {
          # Combine all tables
          combined_numsig_table <- dplyr::bind_rows(all_numsig_tables)

          # Write the combined table
          table_path <- file.path(numsigde_dir, paste0(file_prefix, "_num_sig_da_molecules.tab"))
          vroom::vroom_write(combined_numsig_table, table_path)
          message(sprintf("   outputDaResultsAllContrasts: Wrote NumSigDE table: %s", basename(table_path)))

          # Also write as Excel
          excel_path <- file.path(numsigde_dir, paste0(file_prefix, "_num_sig_da_molecules.xlsx"))
          writexl::write_xlsx(combined_numsig_table, excel_path)
          message(sprintf("   outputDaResultsAllContrasts: Wrote NumSigDE Excel: %s", basename(excel_path)))

          # Generate combined barplot from the aggregated data
          # Filter to only significant results
          sig_only <- combined_numsig_table |>
            dplyr::filter(status != "Not significant")

          if (nrow(sig_only) > 0) {
            num_sig_de_barplot <- sig_only |>
              ggplot2::ggplot(ggplot2::aes(x = status, y = counts)) +
              ggplot2::geom_bar(stat = "identity", fill = "steelblue") +
              ggplot2::geom_text(stat = "identity", ggplot2::aes(label = counts), vjust = -0.5) +
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
                ggplot2::facet_wrap(~comparison, scales = "free_x")
            }

            # Calculate appropriate width based on number of comparisons
            num_comparisons <- length(unique(sig_only$comparison))
            plot_width <- max(7, (num_comparisons + 2) * 7 / 6)

            # Save as PNG
            png_path <- file.path(numsigde_dir, paste0(file_prefix, "_num_sig_da_molecules.png"))
            ggplot2::ggsave(png_path, num_sig_de_barplot, width = plot_width, height = 6, dpi = 300)
            message(sprintf("   outputDaResultsAllContrasts: Saved NumSigDE barplot PNG: %s", basename(png_path)))

            # Save as PDF
            pdf_path <- file.path(numsigde_dir, paste0(file_prefix, "_num_sig_da_molecules.pdf"))
            ggplot2::ggsave(pdf_path, num_sig_de_barplot, width = plot_width, height = 6)
            message(sprintf("   outputDaResultsAllContrasts: Saved NumSigDE barplot PDF: %s", basename(pdf_path)))
          } else {
            message("   outputDaResultsAllContrasts: No significant DE molecules found - skipping barplot")
          }
        } else {
          message("   outputDaResultsAllContrasts: No NumSigDE data found in contrast results - skipping")
        }
      },
      error = function(e) {
        message(sprintf("   outputDaResultsAllContrasts: ERROR processing NumSigDaMolecules: %s", e$message))
      }
    )

    message("--- Exiting outputDaResultsAllContrasts ---")
    return(TRUE)
  }
)

# ----------------------------------------------------------------------------
# getDaResultsLongFormat
# ----------------------------------------------------------------------------
# Get the differential abundance results in wide format
#' @export
setMethod(
  f = "getDaResultsLongFormat",
  signature = "list",
  definition = function(objectsList) {
    return_object_list <- purrr::map(objectsList, function(object) {
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
        dplyr::bind_rows(.id = "comparison") |>
        dplyr::mutate(comparision_short = str_split_i(comparison, "=", 1))

      # Pivot the tidy data frame to a wide format using the correct ID column
      long_results <- tidy_results |>
        # Join with original counts using the correct ID column.
        dplyr::left_join(counts_table_to_use, by = join_by(!!sym(id_col_name) == !!sym(id_col_name))) |>
        dplyr::distinct()

      print(head(long_results))

      # Assign to the correct slot and return the object
      object@results_table_long <- long_results
      return(object)
    })

    return(return_object_list)
  }
)

# ----------------------------------------------------------------------------
# getDaResultsWideFormat
# ----------------------------------------------------------------------------
#' @export
setMethod(
  f = "getDaResultsWideFormat",
  signature = "list",
  definition = function(
    objectsList,
    qvalue_column = "fdr_qvalue",
    raw_pvalue_column = "raw_pvalue",
    log2fc_column = "logFC"
  ) {
    return_object_list <- purrr::map(objectsList, function(object) {
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
        dplyr::bind_rows(.id = "comparison") |>
        dplyr::mutate(comparision_short = str_split_i(comparison, "=", 1))

      # Pivot the tidy data frame to a wide format using the correct ID column
      wida_results <- tidy_results |>
        tidyr::pivot_wider(
          id_cols = c(!!sym(id_col_name)),
          names_from = c(comparision_short),
          names_sep = ":",
          values_from = c(!!sym(log2fc_column), !!sym(qvalue_column), !!sym(raw_pvalue_column))
        ) |>
        # Join with original counts using the correct ID column.
        dplyr::left_join(counts_table_to_use, by = join_by(!!sym(id_col_name) == !!sym(id_col_name))) |>
        dplyr::arrange(dplyr::across(matches(qvalue_column))) |>
        dplyr::distinct()

      # Assign to the correct slot and return the object
      object@results_table_wide <- wida_results
      return(object)
    })

    return(return_object_list)
  }
)

