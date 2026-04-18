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

# ----------------------------------------------------------------------------
# generateProtDAVolcanoStatic
# ----------------------------------------------------------------------------
#' Generate Static Volcano Plot for Proteomics DA Results
#'
#' @description Creates a static volcano plot using ggplot2 and ggrepel.
#'
#' @param da_results_list A list containing DA results (da_proteins_long, etc.)
#' @param selected_contrast Character string of the contrast to plot.
#' @param da_q_val_thresh Numeric threshold for q-value significance (default: 0.05).
#' @param lfc_threshold Numeric threshold for log2 fold-change (default: 1).
#' @param show_labels Logical, whether to label top significant proteins.
#' @param n_labels Number of top proteins to label (default: 10).
#'
#' @return A ggplot2 object.
#'
#' @importFrom ggplot2 ggplot aes geom_point geom_hline geom_vline labs theme_minimal scale_color_manual
#' @importFrom ggrepel geom_text_repel
#' @importFrom dplyr filter arrange slice_head
#' @export
generateProtDAVolcanoStatic <- function(
  da_results_list,
  selected_contrast = NULL,
  da_q_val_thresh = 0.05,
  lfc_threshold = 1,
  show_labels = TRUE,
  n_labels = 10
) {
  if (is.null(da_results_list) || is.null(da_results_list$da_proteins_long)) {
    return(NULL)
  }

  if (is.null(selected_contrast)) {
    return(NULL)
  }

  # Extract comparison name
  comparison_to_search <- stringr::str_extract(selected_contrast, "^[^=]+")
  if (is.na(comparison_to_search)) {
    comparison_to_search <- selected_contrast
  }

  # Filter data for selected contrast
  plot_data <- da_results_list$da_proteins_long |>
    dplyr::filter(comparison == comparison_to_search)

  if (nrow(plot_data) == 0) {
    # Try fuzzy match if exact match fails
    available_comparisons <- unique(da_results_list$da_proteins_long$comparison)
    matching_key <- which(stringr::str_detect(available_comparisons, fixed(comparison_to_search)))
    if (length(matching_key) > 0) {
      comparison_to_search <- available_comparisons[matching_key[1]]
      plot_data <- da_results_list$da_proteins_long |>
        dplyr::filter(comparison == comparison_to_search)
    }
  }

  if (nrow(plot_data) == 0) {
    return(NULL)
  }

  # Add plotting columns
  # CRITICAL FIX: Look up the ID column dynamically instead of hardcoding uniprot_acc
  id_col <- intersect(c("Protein.Ids", "Protein.ID", "Entry", "uniprot_acc", "sites_id"), names(plot_data))
  if (length(id_col) > 0) {
    id_col <- id_col[1]
    message(paste("   generateProtDAVolcanoStatic: Using ID column =", id_col))
  } else {
    # Fallback to first column if none of the above match
    id_col <- names(plot_data)[1]
    message(paste("   generateProtDAVolcanoStatic: Fallback to first column as ID =", id_col))
  }

  plot_data <- plot_data |>
    dplyr::mutate(
      neg_log10_q = -log10(fdr_qvalue),
      display_name = if ("gene_name" %in% names(plot_data)) {
        ifelse(!is.na(gene_name) & gene_name != "", gene_name, !!sym(id_col))
      } else {
        !!sym(id_col)
      },
      significant_label = case_when(
        fdr_qvalue < as.double(da_q_val_thresh) & log2FC >= lfc_threshold ~ "Up",
        fdr_qvalue < as.double(da_q_val_thresh) & log2FC <= -lfc_threshold ~ "Down",
        TRUE ~ "NS"
      )
    )

  # Define colors
  volcano_colors <- c("Up" = "#d9534f", "Down" = "#5bc0de", "NS" = "grey70")

  # Create base plot
  p <- ggplot2::ggplot(plot_data, aes(x = log2FC, y = neg_log10_q, color = significant_label)) +
    ggplot2::geom_point(alpha = 0.6, size = 1.5) +
    ggplot2::geom_hline(yintercept = -log10(da_q_val_thresh), linetype = "dashed", color = "grey50") +
    ggplot2::geom_vline(xintercept = c(-lfc_threshold, lfc_threshold), linetype = "dashed", color = "grey50") +
    ggplot2::scale_color_manual(values = volcano_colors) +
    ggplot2::labs(
      title = paste("Volcano Plot:", comparison_to_search),
      x = "Log2 Fold Change",
      y = "-Log10 Q-value",
      color = "Status"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      legend.position = "bottom",
      plot.title = element_text(hjust = 0.5, face = "bold")
    )

  # Add labels for top proteins
  if (show_labels) {
    top_to_label <- plot_data |>
      dplyr::filter(significant_label != "NS") |>
      dplyr::arrange(fdr_qvalue) |>
      dplyr::slice_head(n = n_labels)

    if (nrow(top_to_label) > 0) {
      p <- p + ggrepel::geom_text_repel(
        data = top_to_label,
        aes(label = display_name),
        size = 3,
        box.padding = 0.5,
        max.overlaps = 15,
        show.legend = FALSE
      )
    }
  }

  return(p)
}

# ----------------------------------------------------------------------------
# writeInteractiveVolcanoPlotProteomics
# ----------------------------------------------------------------------------
## Create proteomics interactive volcano plot
#' @export
# da_analysis_results_list$contrasts_results$fit.eb
# No full stops in the nme of columns of interactive table in glimma plot. It won't display column with full stop in the column name.
writeInteractiveVolcanoPlotProteomics <- function(
  da_proteins_long,
  uniprot_tbl,
  fit.eb = NULL, # No longer strictly required but kept for signature compatibility
  publication_graphs_dir,
  args_row_id = "uniprot_acc",
  fdr_column = "fdr_qvalue",
  raw_p_value_column = "raw_pvalue",
  log2fc_column = "log2FC",
  da_q_val_thresh = 0.05,
  counts_tbl = NULL,
  groups = NULL,
  uniprot_id_column = "Entry",
  gene_names_column = "gene_names",
  display_columns = c("best_uniprot_acc")
) {

  logger::log_info("--- Entering writeInteractiveVolcanoPlotProteomics ---")

  # Dynamically identify row ID column 
  if (!(args_row_id %in% names(da_proteins_long))) {
    potential_ids <- c("Protein.Ids", "Protein.ID", "Entry", "uniprot_acc", "sites_id")
    found_id <- names(da_proteins_long)[names(da_proteins_long) %in% potential_ids][1]
    if (!is.na(found_id)) args_row_id <- found_id
  }

  volcano_plot_tab <- da_proteins_long |>
    dplyr::mutate(
      best_uniprot_acc = stringr::str_extract(!!sym(args_row_id), "^[^|:]+"),
      best_uniprot_acc_base = gsub("-\\d+$", "", best_uniprot_acc)
    )

  if (!is.null(uniprot_tbl)) {
    id_col_uniprot <- intersect(c("Entry", "UNIPROTKB", "Protein.Ids"), names(uniprot_tbl))[1]
    gene_col_uniprot <- intersect(c("gene_names", "GENENAME", "Gene.Name", "GeneName"), names(uniprot_tbl))[1]

    if (!is.na(id_col_uniprot) && !is.na(gene_col_uniprot)) {
      mapping_df <- uniprot_tbl |>
        dplyr::select(all_of(c(id_col_uniprot, gene_col_uniprot))) |>
        dplyr::rename(Entry = !!sym(id_col_uniprot), GeneSymbol = !!sym(gene_col_uniprot)) |>
        dplyr::mutate(GeneSymbol = stringr::str_extract(GeneSymbol, "^[^ ;:]+")) |>
        dplyr::distinct(Entry, .keep_all = TRUE)

      volcano_plot_tab <- volcano_plot_tab |>
        dplyr::left_join(mapping_df, by = c("best_uniprot_acc_base" = "Entry")) |>
        dplyr::mutate(
          gene_name = dplyr::coalesce(GeneSymbol, best_uniprot_acc)
        )
    } else {
      volcano_plot_tab$gene_name <- volcano_plot_tab$best_uniprot_acc
    }
  } else {
    volcano_plot_tab$gene_name <- volcano_plot_tab$best_uniprot_acc
  }


  output_dir <- file.path(
    publication_graphs_dir,
    "Interactive_Volcano_Plots"
  )

  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  # Iterate over uniquely defined contrasts instead of relying on coefficients index
  unique_contrasts <- unique(da_proteins_long$comparison)
  logger::log_info(sprintf("   Found %d unique contrasts to plot", length(unique_contrasts)))
  
  purrr::walk(
    unique_contrasts,
    \(contrast_name) {
      logger::log_info(sprintf("   Generating static interactive plot for: %s", contrast_name))

      contrast_data <- volcano_plot_tab |> 
        dplyr::filter(comparison == contrast_name)
      
      # Filter counts to only the IDs in the contrast results to avoid length mismatches
      counts_mat_filtered <- NULL
      if (!is.null(counts_tbl)) {
        ids_in_contrast <- contrast_data[[args_row_id]]
        valid_ids <- intersect(ids_in_contrast, rownames(counts_tbl))
        if(length(valid_ids) > 0) {
           counts_mat_filtered <- counts_tbl[valid_ids, , drop = FALSE]
        }
      }

      getGlimmaVolcanoProteomics(
        volcano_plot_tab = contrast_data,
        uniprot_column = best_uniprot_acc,
        gene_name_column = gene_name,
        display_columns = display_columns,
        counts_tbl = counts_mat_filtered,
        groups = groups,
        output_dir = output_dir,
        fdr_column = fdr_column,
        log2fc_column = log2fc_column,
        da_q_val_thresh = da_q_val_thresh,
        contrast_name = contrast_name
      )
    }
  )
  logger::log_info("--- Exiting writeInteractiveVolcanoPlotProteomics ---")
}

# ----------------------------------------------------------------------------
# writeInteractiveVolcanoPlotProteomicsWidget
# ----------------------------------------------------------------------------
#' @export
# da_analysis_results_list$contrasts_results$fit.eb
# No full stops in the nme of columns of interactive table in glimma plot. It won't display column with full stop in the column name.
writeInteractiveVolcanoPlotProteomicsWidget <- function(
  da_proteins_long,
  uniprot_tbl,
  fit.eb,
  args_row_id = "uniprot_acc",
  fdr_column = "fdr_qvalue",
  raw_p_value_column = "raw_pvalue",
  log2fc_column = "log2FC",
  da_q_val_thresh = 0.05,
  counts_tbl = NULL,
  groups = NULL,
  uniprot_id_column = "Entry",
  gene_names_column = "gene_names",
  display_columns = c("best_uniprot_acc")
) {
  volcano_plot_tab <- da_proteins_long |>
    left_join(uniprot_tbl, by = join_by(!!sym(args_row_id) == !!sym(uniprot_id_column))) |>
    dplyr::rename(UNIPROT_GENENAME = gene_names_column) |>
    mutate(UNIPROT_GENENAME = purrr::map_chr(UNIPROT_GENENAME, \(x){
      tryCatch(
        {
          if (is.na(x) || is.null(x) || x == "") {
            ""
          } else {
            split_result <- str_split(x, " ")[[1]]
            if (length(split_result) > 0) split_result[1] else ""
          }
        },
        error = function(e) ""
      )
    })) |>
    mutate(lqm = -log10(!!sym(fdr_column))) |>
    dplyr::mutate(label = case_when(
      abs(!!sym(log2fc_column)) >= 1 & !!sym(fdr_column) >= da_q_val_thresh ~ "Not sig., logFC >= 1",
      abs(!!sym(log2fc_column)) >= 1 & !!sym(fdr_column) < da_q_val_thresh ~ "Sig., logFC >= 1",
      abs(!!sym(log2fc_column)) < 1 & !!sym(fdr_column) < da_q_val_thresh ~ "Sig., logFC < 1",
      TRUE ~ "Not sig."
    )) |>
    dplyr::mutate(colour = case_when(
      abs(!!sym(log2fc_column)) >= 1 & !!sym(fdr_column) >= da_q_val_thresh ~ "orange",
      abs(!!sym(log2fc_column)) >= 1 & !!sym(fdr_column) < da_q_val_thresh ~ "purple",
      abs(!!sym(log2fc_column)) < 1 & !!sym(fdr_column) < da_q_val_thresh ~ "blue",
      TRUE ~ "black"
    )) |>
    dplyr::mutate(gene_name = str_split(UNIPROT_GENENAME, " |:") |> purrr::map_chr(1)) |>
    dplyr::mutate(best_uniprot_acc = str_split(!!sym(args_row_id), ":") |> purrr::map_chr(1)) |>
    dplyr::mutate(analysis_type = comparison) |>
    dplyr::select(
      best_uniprot_acc, lqm, !!sym(fdr_column), !!sym(raw_p_value_column), !!sym(log2fc_column), comparison, label, colour, gene_name,
      any_of(display_columns)
    ) |>
    dplyr::mutate(my_alpha = case_when(
      gene_name != "" ~ 1,
      TRUE ~ 0.5
    ))

  interactive_volcano_plots <- purrr::map(
    seq_len(ncol(fit.eb$coefficients)),
    \(coef) {
      # print(paste0( "coef = ", coef))
      getGlimmaVolcanoProteomicsWidget(fit.eb,
        coef = coef,
        volcano_plot_tab = volcano_plot_tab,
        uniprot_column = best_uniprot_acc,
        gene_name_column = gene_name,
        display_columns = display_columns,
        counts_tbl = counts_tbl,
        groups = groups
      )
    }
  )

  interactive_volcano_plots
}

# ----------------------------------------------------------------------------
# writeInteractiveVolcanoPlotProteomicsMain
# ----------------------------------------------------------------------------
#' @export
writeInteractiveVolcanoPlotProteomicsMain <- function(
  da_analysis_results_list,
  theObject,
  uniprot_tbl,
  de_analysis_results_list = NULL,
  publication_graphs_dir = NULL,
  file_prefix = NULL,
  plots_format = NULL,
  args_row_id = NULL,
  da_q_val_thresh = NULL,
  de_q_val_thresh = NULL,
  gene_names_column = NULL,
  fdr_column = NULL,
  raw_p_value_column = NULL,
  log2fc_column = NULL,
  uniprot_id_column = NULL,
  display_columns = NULL
) {
  # Accept legacy de_* aliases but normalize to the current da_* path internally.
  results_list <- da_analysis_results_list
  if (is.null(results_list) && !is.null(de_analysis_results_list)) {
    results_list <- de_analysis_results_list
  }
  if (is.null(results_list)) {
    stop("A DA/DE analysis results list must be supplied.", call. = FALSE)
  }

  safe_get_slot <- function(res_list, primary_name, legacy_name = NULL) {
    if (!is.null(res_list[[primary_name]])) {
      return(res_list[[primary_name]])
    }
    if (!is.null(legacy_name) && !is.null(res_list[[legacy_name]])) {
      return(res_list[[legacy_name]])
    }
    NULL
  }

  uniprot_tbl <- checkParamsObjectFunctionSimplify(theObject, "uniprot_tbl", NULL)
  publication_graphs_dir <- checkParamsObjectFunctionSimplify(theObject, "publication_graphs_dir", NULL)
  args_row_id <- checkParamsObjectFunctionSimplify(theObject, "args_row_id", "uniprot_acc")
  da_q_val_thresh <- if (!is.null(da_q_val_thresh)) {
    da_q_val_thresh
  } else if (!is.null(de_q_val_thresh)) {
    de_q_val_thresh
  } else {
    checkParamsObjectFunctionSimplify(theObject, "da_q_val_thresh", 0.05)
  }
  gene_names_column <- checkParamsObjectFunctionSimplify(theObject, "gene_names_column", "gene_names")
  fdr_column <- checkParamsObjectFunctionSimplify(theObject, "fdr_column", "fdr_qvalue")
  raw_p_value_column <- checkParamsObjectFunctionSimplify(theObject, "raw_p_value_column", "raw_pvalue")
  log2fc_column <- checkParamsObjectFunctionSimplify(theObject, "log2fc_column", "log2FC")
  uniprot_id_column <- checkParamsObjectFunctionSimplify(theObject, "uniprot_id_column", "Entry")
  display_columns <- checkParamsObjectFunctionSimplify(theObject, "display_columns", c("best_uniprot_acc"))

  theObject <- updateParamInObject(theObject, "uniprot_tbl")
  theObject <- updateParamInObject(theObject, "publication_graphs_dir")
  theObject <- updateParamInObject(theObject, "args_row_id")
  theObject <- updateParamInObject(theObject, "da_q_val_thresh")
  theObject <- updateParamInObject(theObject, "gene_names_column")
  theObject <- updateParamInObject(theObject, "fdr_column")
  theObject <- updateParamInObject(theObject, "raw_p_value_column")
  theObject <- updateParamInObject(theObject, "log2fc_column")
  theObject <- updateParamInObject(theObject, "uniprot_id_column")
  theObject <- updateParamInObject(theObject, "display_columns")


  ## Write interactive volcano plot

  da_proteins_long <- safe_get_slot(results_list, "da_proteins_long", "de_proteins_long")
  contrasts_results <- safe_get_slot(results_list, "contrasts_results")

  if (is.null(da_proteins_long) || is.null(contrasts_results)) {
    stop("Results list does not contain the expected DA/DE volcano inputs.", call. = FALSE)
  }

  # Use helper to extract counts table instead of direct slot access
  result_object <- safe_get_slot(results_list, "theObject")
  if (is.null(result_object)) {
    result_object <- theObject
  }

  counts_mat <- getCountsTable(result_object) |>
    column_to_rownames(args_row_id) |>
    as.matrix()

  this_design_matrix <- result_object@design_matrix

  rownames(this_design_matrix) <- this_design_matrix[, result_object@sample_id]

  this_groups <- this_design_matrix[colnames(counts_mat), "group"]

  writeInteractiveVolcanoPlotProteomics(da_proteins_long,
    uniprot_tbl = uniprot_tbl,
    fit.eb = contrasts_results$fit.eb,
    publication_graphs_dir = publication_graphs_dir,
    args_row_id = args_row_id,
    fdr_column = fdr_column,
    raw_p_value_column = raw_p_value_column,
    log2fc_column = log2fc_column,
    da_q_val_thresh = da_q_val_thresh,
    counts_tbl = counts_mat,
    groups = this_groups,
    uniprot_id_column = uniprot_id_column,
    gene_names_column = gene_names_column,
    display_columns = display_columns
  )
}
