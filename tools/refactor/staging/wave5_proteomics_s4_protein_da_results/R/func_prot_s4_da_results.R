#' @export
setMethod(
  f = "outputDaResultsAllContrasts",
  signature = "ProteinQuantitativeData",
  definition = function(theObject,
                        da_results_list_all_contrasts = NULL,
                        uniprot_tbl = NULL,
                        da_output_dir = NULL,
                        publication_graphs_dir = NULL,
                        file_prefix = "de_proteins",
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

        # [NEW]: Generate volcano plot for this contrast
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

    # [NEW]: Generate combined multi-page PDF with all volcano plots
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

    message("--- Exiting outputDaResultsAllContrasts ---")
    return(TRUE)
  }
)

