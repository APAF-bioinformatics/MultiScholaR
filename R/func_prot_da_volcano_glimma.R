# ----------------------------------------------------------------------------
# generateProtDAVolcanoPlotGlimma
# ----------------------------------------------------------------------------
## ----------------------------------------------------------------------------------------------------------------------------------------------------------------------
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
generateProtDAVolcanoPlotGlimma <- function(
  da_results_list,
  selected_contrast = NULL,
  uniprot_tbl = NULL,
  args_row_id = "uniprot_acc",
  fdr_column = "fdr_qvalue",
  raw_p_value_column = "raw_pvalue",
  log2fc_column = "log2FC",
  da_q_val_thresh = 0.05,
  uniprot_id_column = "Entry",
  gene_names_column = "gene_names",
  display_columns = c("best_uniprot_acc"),
  output_dir = NULL,
  ...
) {
  logger::log_info("--- Entering generateProtDAVolcanoPlotGlimma (glimmaXY refactor) ---")

  # --- Semi-automated Test Capture Logic ---
  if (getOption("multischolar.capture_glimma_data", FALSE)) {
    # tryCatch({
    #   # Use provided output_dir or fallback to a temp directory
    #   base_dir <- if (!is.null(output_dir)) output_dir else tempdir()
    #   capture_dir <- file.path(base_dir, "glimma_fixtures")
    #
    #   if (!dir.exists(capture_dir)) dir.create(capture_dir, recursive = TRUE)
    #
    #   timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
    #   filename <- file.path(capture_dir, paste0("glimma_snapshot_", timestamp, ".rds"))
    #
    #   logger::log_info(sprintf("   Capturing test data to: %s", filename))
    #   saveRDS(list(
    #     da_results_list = da_results_list,
    #     selected_contrast = selected_contrast,
    #     uniprot_tbl = uniprot_tbl,
    #     args_row_id = args_row_id,
    #     fdr_column = fdr_column,
    #     raw_p_value_column = raw_p_value_column,
    #     log2fc_column = log2fc_column,
    #     da_q_val_thresh = da_q_val_thresh,
    #     uniprot_id_column = uniprot_id_column,
    #     gene_names_column = gene_names_column,
    #     display_columns = display_columns
    #   ), file = filename)
    # }, error = function(e) {
    #   logger::log_warn(sprintf("   Failed to capture test data: %s", conditionMessage(e)))
    # })
  }
  # ------------------------------------------

  logger::log_info(sprintf("   selected_contrast = %s", as.character(selected_contrast)[1]))

  if (is.null(da_results_list) || is.null(da_results_list$da_proteins_long)) {
    logger::log_warn("   No DA results available")
    return(NULL)
  }

  if (is.null(selected_contrast)) {
    logger::log_info("   No contrast selected")
    return(NULL)
  }

  da_proteins_long <- da_results_list$da_proteins_long
  comparison_to_search <- selected_contrast

  # Robust contrast matching
  # First, log available contrasts for debugging
  available_contrasts <- unique(da_proteins_long$comparison)
  logger::log_info(sprintf("   Glimma: Available contrasts in data: %s", paste(available_contrasts, collapse = ", ")))
  logger::log_info(sprintf("   Glimma: Target contrast: %s", as.character(selected_contrast)[1]))

  contrast_data <- da_proteins_long |>
    dplyr::filter(comparison == selected_contrast)

  if (nrow(contrast_data) == 0) {
    # Try fuzzy match if exact match fails
    # 1. Try removing characters like '=', '-', ' '
    clean_target <- gsub("[^A-Za-z0-9]", "", selected_contrast)

    match_idx <- which(sapply(available_contrasts, function(x) {
      gsub("[^A-Za-z0-9]", "", x) == clean_target
    }))

    if (length(match_idx) > 0) {
      comparison_to_search <- available_contrasts[match_idx[1]]
      logger::log_info(sprintf("   Glimma: Found fuzzy match: %s", comparison_to_search))
      contrast_data <- da_proteins_long |> dplyr::filter(comparison == comparison_to_search)
    }
  }

  if (nrow(contrast_data) == 0) {
    # Try prefix match
    potential_prefix <- stringr::str_extract(selected_contrast, "^[^=]+")
    if (is.na(potential_prefix)) potential_prefix <- stringr::str_extract(selected_contrast, "^[^-]+")

    if (!is.na(potential_prefix)) {
       match_idx <- grepl(potential_prefix, da_proteins_long$comparison, fixed = TRUE)
       if (any(match_idx)) {
         comparison_to_search <- da_proteins_long$comparison[which(match_idx)[1]]
         contrast_data <- da_proteins_long |> dplyr::filter(comparison == comparison_to_search)
         logger::log_info(sprintf("   Glimma: Found prefix match: %s", comparison_to_search))
       }
    }
  }

  if (nrow(contrast_data) == 0 && length(available_contrasts) > 0) {
    # Final fallback: just use the first available contrast
    first_contrast <- available_contrasts[1]
    logger::log_warn(sprintf("   Glimma: No data found for contrast %s. Falling back to first available: %s", selected_contrast, first_contrast))
    contrast_data <- da_proteins_long |> dplyr::filter(comparison == first_contrast)
    comparison_to_search <- first_contrast
  }

  if (nrow(contrast_data) == 0) {
    logger::log_warn(sprintf("   No data found for contrast %s or any matching variants", selected_contrast))
    return(NULL)
  }

  # Dynamically identify ID column if not already found
  if (!(args_row_id %in% names(da_proteins_long))) {
    potential_ids <- c("Protein.Ids", "Protein.ID", "Entry", "uniprot_acc", "sites_id")
    found_id <- names(da_proteins_long)[names(da_proteins_long) %in% potential_ids][1]

    if (!is.na(found_id)) {
      args_row_id <- found_id
    } else {
      args_row_id <- names(da_proteins_long)[1]
    }
  }

  # Prepare annotation data
  plot_data <- contrast_data |>
    dplyr::mutate(
      best_uniprot_acc = purrr::map_chr(as.character(Protein.Ids), \(x) {
        # Robustly handle FASTA/UniProt/NCBI headers
        # Remove common prefix like '>'
        x_cln <- gsub("^>", "", x)

        # Support both | and : as delimiters
        parts <- unlist(strsplit(x_cln, "\\||:"))

        if (length(parts) > 1) {
          # Handle sp|acc|name, tr|acc|name, lcl|acc
          if (trimws(parts[1]) %in% c("sp", "tr", "lcl")) {
             return(trimws(parts[2]))
          }
          # Handle common case where ID is just the first part
          return(trimws(parts[1]))
        }
        return(trimws(x_cln))
      }),
      best_uniprot_acc_base = gsub("-\\d+$", "", best_uniprot_acc)
    )


  # Merge UniProt annotations if available
  if (!is.null(uniprot_tbl)) {
    id_col_uniprot <- intersect(c("Entry", "UNIPROTKB", "Protein.Ids"), names(uniprot_tbl))[1]
    gene_col_uniprot <- intersect(c("gene_names", "GENENAME", "Gene.Name", "GeneName"), names(uniprot_tbl))[1]

    if (!is.na(id_col_uniprot) && !is.na(gene_col_uniprot)) {
      mapping_df <- uniprot_tbl |>
        dplyr::select(all_of(c(id_col_uniprot, gene_col_uniprot))) |>
        dplyr::rename(Entry = !!sym(id_col_uniprot), GeneSymbol = !!sym(gene_col_uniprot)) |>
        dplyr::mutate(GeneSymbol = stringr::str_extract(GeneSymbol, "^[^ ;:]+")) |>
        dplyr::distinct(Entry, .keep_all = TRUE)

      plot_data <- plot_data |>
        dplyr::left_join(mapping_df, by = c("best_uniprot_acc_base" = "Entry")) |>
        dplyr::mutate(
          gene_name = dplyr::coalesce(GeneSymbol, best_uniprot_acc)
        )
    } else {
      plot_data$gene_name <- plot_data$best_uniprot_acc
    }
  } else {
    plot_data$gene_name <- plot_data$best_uniprot_acc
  }

  # Prepare plotting metrics
  # Handle zero FDR to avoid Inf values in -log10
  plot_data <- plot_data |>
    dplyr::mutate(
      FDR = as.numeric(!!sym(fdr_column)),
      # Avoid having zero FDR which leads to Inf in -log10
      FDR = ifelse(FDR == 0, min(FDR[FDR > 0], na.rm = TRUE) * 0.1, FDR),
      negLog10FDR = -log10(FDR),
      logFC = as.numeric(!!sym(log2fc_column)),
      Status = dplyr::case_when(
        logFC >= 1 & FDR < as.double(da_q_val_thresh) ~ 1,
        logFC <= -1 & FDR < as.double(da_q_val_thresh) ~ -1,
        TRUE ~ 0
      )
    )

  # Remove rows with Inf/NA plotting metrics AND ensure Protein.Ids are unique for rownames
  plot_data <- plot_data |>
    dplyr::filter(!is.infinite(negLog10FDR), !is.na(negLog10FDR), !is.na(logFC)) |>
    # Ensure Protein.Ids is not blank
    dplyr::filter(!is.na(Protein.Ids), Protein.Ids != "") |>
    dplyr::distinct(Protein.Ids, .keep_all = TRUE)


  if (nrow(plot_data) == 0) {
    logger::log_warn("   Skipping glimmaXY - no valid rows after removing Inf/NA and duplicates")
    return(NULL)
  }

  # Prepare Counts Matrix with GUID rownames (to avoid tibble rownames crash and duplicates)
  counts_mat <- NULL
  groups_vec <- NULL
  if (!is.null(da_results_list$theObject)) {
    # CRITICAL FIX: Use check.names = FALSE to prevent sample name mangling (e.g. replacing spaces with dots)
    counts_df <- as.data.frame(da_results_list$theObject@protein_quant_table, check.names = FALSE)

    if (args_row_id %in% names(counts_df)) {
      # 1. Ensure exact match between plot_data and counts_df
      common_prots <- intersect(plot_data$Protein.Ids, counts_df[[args_row_id]])
      plot_data <- plot_data |> dplyr::filter(Protein.Ids %in% common_prots)

      counts_df <- counts_df |> dplyr::filter(!!sym(args_row_id) %in% common_prots)

      # 2. Ensure exact row order mapped to plot_data
      match_idx <- match(plot_data$Protein.Ids, counts_df[[args_row_id]])
      counts_df <- counts_df[match_idx, , drop = FALSE]

      # 3. Avoid tibble rownames crash by ensuring it's a data.frame and reset rownames first
      counts_df <- as.data.frame(counts_df, check.names = FALSE)
      rownames(counts_df) <- NULL
      counts_mat <- counts_df |> tibble::column_to_rownames(args_row_id) |> as.matrix()

      this_design <- da_results_list$theObject@design_matrix
      sample_id_col <- da_results_list$theObject@sample_id
      group_id_col <- da_results_list$theObject@group_id

      # MATCHING SAMPLES FOR GROUPING
      if (group_id_col %in% names(this_design)) {
        common_samples <- intersect(trimws(tolower(colnames(counts_mat))),
                                    trimws(tolower(as.character(this_design[[sample_id_col]]))))

        if (length(common_samples) > 0) {
          # Use original case for the continuous mapping definition
          design_mapping <- data.frame(
            orig_name = as.character(this_design[[sample_id_col]]),
            norm_name = trimws(tolower(as.character(this_design[[sample_id_col]]))),
            group = as.character(this_design[[group_id_col]]),
            stringsAsFactors = FALSE
          ) |>
          dplyr::filter(norm_name %in% common_samples)

          # Match matrix columns to design exactly to prevent ordering issues
          mat_cols_norm <- trimws(tolower(colnames(counts_mat)))
          cols_to_keep <- colnames(counts_mat)[mat_cols_norm %in% design_mapping$norm_name]

          # Subset matrix to only those samples
          counts_mat <- counts_mat[, cols_to_keep, drop = FALSE]

          # Re-evaluate samples in matrix after subsetting
          samples_in_mat <- colnames(counts_mat)
          groups_vec <- sapply(samples_in_mat, function(s) {
            match_idx <- match(trimws(tolower(s)), design_mapping$norm_name)
            if (!is.na(match_idx)) return(design_mapping$group[match_idx])
            return("Unknown")
          })

          # CRITICAL: Glimma often expects a factor for groups to correctly partition the expression plot
          groups_vec <- factor(groups_vec)

          logger::log_info(sprintf("   Glimma: Grouping expression plot by '%s' with %d samples", group_id_col, length(groups_vec)))
          logger::log_info(sprintf("   Glimma: Unique groups detected: %s", paste(unique(groups_vec), collapse = ", ")))
        } else {
          logger::log_warn("   Glimma: No common samples found (case-insensitive check failed)")
        }
      } else {
        logger::log_warn(sprintf("   Glimma: Group column '%s' not found in design matrix", group_id_col))
      }

      # Sync plot_data to available proteins in counts (defensive)
      # Coerce both to character for matching
      available_ids <- intersect(as.character(plot_data$Protein.Ids), as.character(rownames(counts_mat)))
      plot_data <- plot_data |> dplyr::filter(as.character(Protein.Ids) %in% available_ids)
      counts_mat <- counts_mat[as.character(plot_data$Protein.Ids), , drop = FALSE]

      logger::log_info(sprintf("   Glimma: Synchronized %d proteins with expression data", nrow(counts_mat)))
    }
  }

  # Build Clean Annotation Dataframe
  clean_anno <- data.frame(
    Protein_Ids = as.character(plot_data$Protein.Ids),
    gene = ifelse(is.na(plot_data$gene_name) | plot_data$gene_name == "",
                  gsub("^lcl\\|", "", as.character(plot_data$Protein.Ids)),
                  as.character(plot_data$gene_name)),
    gene_name = ifelse(is.na(plot_data$gene_name) | plot_data$gene_name == "",
                       gsub("^lcl\\|", "", as.character(plot_data$Protein.Ids)),
                       as.character(plot_data$gene_name)),
    FDR = signif(as.numeric(plot_data$FDR), 4),
    stringsAsFactors = FALSE
  )

  # Check for UniProt accession in mapping table if available
  if ("best_uniprot_acc" %in% names(plot_data)) {
    clean_anno$UniProt <- sapply(as.character(plot_data$best_uniprot_acc), function(x) {
      if (is.na(x) || x == "") return("")
      if (grepl("\\|", x)) {
        parts <- strsplit(as.character(x), "\\|")[[1]]
        if (length(parts) >= 2) return(parts[2])
      }
      return(as.character(x))
    })
  } else {
    clean_anno$UniProt <- gsub("^lcl\\|", "", as.character(plot_data$Protein.Ids))
  }

  # Add other display columns from display_columns argument
  purrr::walk(display_columns, function(col) {
    col_cln <- gsub("\\.", "_", col)
    if (col %in% names(plot_data) && !(col_cln %in% names(clean_anno))) {
      clean_anno[[col_cln]] <- as.character(plot_data[[col]])
    }
  })

  clean_anno[is.na(clean_anno)] <- ""
  rownames(clean_anno) <- as.character(plot_data$Protein.Ids)

  # Unify stats and annotations into a single display_df for Glimma
  display_df <- clean_anno |>
    dplyr::mutate(
      logFC = signif(as.numeric(plot_data$logFC), 4),
      negLog10FDR = signif(as.numeric(plot_data$negLog10FDR), 4)
    ) |>
    dplyr::rename_with(~ gsub("[\\. ]", "_", .))

  # Ensure Protein_Ids is first and IN the columns to display
  if ("Protein_Ids" %in% names(display_df)) {
    display_df <- display_df |> dplyr::relocate(Protein_Ids)
  }

  display_cols_final <- colnames(display_df)

  # FINAL CLEANING for Glimma table stability
  # Glimma table can fail if there are NAs or non-atomic types in the anno dataframe
  logger::log_info("FINAL CLEANING for Glimma table stability")
  display_df <- as.data.frame(display_df)
  for (col in names(display_df)) {
    display_df[[col]] <- as.character(display_df[[col]])
    display_df[[col]][is.na(display_df[[col]])] <- ""
  }

  # 4. Prepare Final Data for Glimma
  # EXACT Pattern from PDI Working Workflow

  # Sanitize the annotation dataframe to prevent D3.js / DataTables crashing
  # 1. Drop math columns, keeping all possible identifier/annotation columns
  # 2. Rename columns to replace dots with underscores (DataTables interprets dots as nested JSON)
  # 3. Convert to character to drop unexpected factors/complex types
  # 4. Replace NAs with empty strings which JSON serializers handle safely
  clean_anno <- plot_data |>
    dplyr::select(-any_of(c("logFC", "FDR", "PValue", "negLog10FDR", "left_group", "right_group", "Status", "log2FC", "fdr_qvalue", "raw_pvalue"))) |>
    dplyr::rename_with(~ gsub("[\\. ]", "_", .)) |>
    dplyr::mutate(dplyr::across(everything(), as.character)) |>
    dplyr::mutate(dplyr::across(everything(), ~ ifelse(is.na(.), "", .))) |>
    as.data.frame()

  # Ensure Protein_Ids is first
  if ("Protein_Ids" %in% names(clean_anno)) {
    clean_anno <- clean_anno |> dplyr::relocate(Protein_Ids)
  }

  # Sync counts_mat rows to plot_data exactly
  if (!is.null(counts_mat)) {
    counts_mat <- counts_mat[as.character(plot_data$Protein.Ids), , drop = FALSE]
  }

  logger::log_info(sprintf("   Glimma: Calling glimmaXY with %d proteins and %d samples",
                           nrow(plot_data), if (is.null(counts_mat)) 0 else ncol(counts_mat)))

  # 5. Call Glimma publicly (standardized pattern)
  # CRITICAL: status must be explicitly passed to highlight significant genes
  # 1 = Up, -1 = Down, 0 = Not Sig
  status_vec <- as.integer(plot_data$Status)

  logger::log_info(sprintf("   Glimma: Highlighting %d significant genes (%d up, %d down)",
                           sum(status_vec != 0), sum(status_vec == 1), sum(status_vec == -1)))

  widget <- Glimma::glimmaXY(
    x = plot_data$logFC,
    y = plot_data$negLog10FDR,
    xlab = "logFC",
    ylab = "negLog10FDR",
    status = status_vec,
    counts = counts_mat,
    groups = groups_vec,
    transform.counts = "none",
    anno = clean_anno,
    status.cols = c("#1052bd", "silver", "#cc212f"),
    main = paste("Interactive Volcano Plot:", comparison_to_search)
  )

  # CRITICAL: Inject CSS fix to re-unify the split DataTables layout
  # This ensures column widths are synchronized under a single scrollbar
  # Using htmlwidgets::prependContent allows this to work in Shiny and standalone
  css_fix <- "<style> .dataTables_wrapper { overflow-x: auto !important; } .dataTables_scroll { display: table !important; width: auto !important; min-width: 100% !important; } .dataTables_scrollHead, .dataTables_scrollBody, .dataTables_scrollHeadInner, .dataTables_scrollHeadInner > table, .dataTables_scrollBody > table { display: contents !important; } .dataTables_scrollBody thead { display: none !important; } .dataTable th, .dataTable td { white-space: nowrap !important; } </style>"
  widget <- htmlwidgets::prependContent(widget, htmltools::HTML(css_fix))

  logger::log_info("--- Exiting generateProtDAVolcanoPlotGlimma ---")
  return(widget)
}
