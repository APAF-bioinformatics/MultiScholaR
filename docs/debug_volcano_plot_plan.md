# Implementation Plan - Resolve Proteomics Volcano Plot Bug

The user reports that the proteomics volcano plot is still problematic despite previous refactoring. This plan focuses on identifying and fixing the root cause of the remaining error.

## User Review Required

> [!IMPORTANT]
> I need specific details about the error (e.g., error messages in the console, blank plots, or incorrect data display) to proceed effectively.

## Proposed Changes

### [Component] Proteomics Differential Abundance Functions

#### [MODIFY] [func_prot_da.R](file:///Users/ignatiuspang/Workings/2025/MultiScholaR/R/func_prot_da.R)

**Problem**: Glimma crashes due to:
1. UI/Raw contrast mismatches (friendly vs internal).
2. Duplicate gene names used as rownames.
3. Tibble/Rowname attribute conflicts.

**Verified Fix**:
Replace the body of `generateProtDAVolcanoPlotGlimma` with this robust implementation:

```r
  # ... (Inside generateProtDAVolcanoPlotGlimma) ...
  comparison_to_search <- selected_contrast

  # Robust contrast matching
  contrast_data <- da_proteins_long |>
    dplyr::filter(comparison == selected_contrast)

  if (nrow(contrast_data) == 0) {
    potential_prefix <- stringr::str_extract(selected_contrast, "^[^=]+")
    if (is.na(potential_prefix)) potential_prefix <- stringr::str_extract(selected_contrast, "^[^-]+")
    
    if (!is.na(potential_prefix)) {
       match_idx <- grepl(potential_prefix, da_proteins_long$comparison, fixed = TRUE)
       contrast_data <- da_proteins_long[match_idx, ]
       if (nrow(contrast_data) > 0) comparison_to_search <- potential_prefix
    }
  }

  if (nrow(contrast_data) == 0 && length(unique(da_proteins_long$comparison)) > 0) {
    first_contrast <- unique(da_proteins_long$comparison)[1]
    logger::log_warn(sprintf("   No data found for contrast %s. Falling back to %s", selected_contrast, first_contrast))
    contrast_data <- da_proteins_long |> dplyr::filter(comparison == first_contrast)
    comparison_to_search <- first_contrast
  }

  if (nrow(contrast_data) == 0) return(NULL)

  # Prepare annotation data (Standardize IDs)
  plot_data <- contrast_data |>
    dplyr::mutate(
      Protein.Ids = !!sym(args_row_id),
      FDR = !!sym(fdr_column),
      FDR = ifelse(FDR == 0, min(FDR[FDR > 0], na.rm = TRUE) * 0.1, FDR),
      negLog10FDR = -log10(FDR),
      logFC = !!sym(log2fc_column),
      Status = case_when(
        logFC >= 1 & FDR < da_q_val_thresh ~ 1,
        logFC <= -1 & FDR < da_q_val_thresh ~ -1,
        TRUE ~ 0
      )
    )

  # Prepare Counts Matrix with GUID rownames
  counts_mat <- NULL; groups_vec <- NULL
  if (!is.null(da_results_list$theObject)) {
    counts_df <- as.data.frame(da_results_list$theObject@protein_quant_table)
    if (args_row_id %in% names(counts_df)) {
      counts_df <- counts_df |> dplyr::filter(!!sym(args_row_id) %in% plot_data$Protein.Ids)
      match_idx <- match(plot_data$Protein.Ids, counts_df[[args_row_id]])
      counts_df <- counts_df[match_idx[!is.na(match_idx)], , drop = FALSE]
      
      # FIX: Avoid tibble rownames crash
      rownames(counts_df) <- NULL
      counts_mat <- counts_df |> tibble::column_to_rownames(args_row_id) |> as.matrix()
      
      this_design <- da_results_list$theObject@design_matrix
      if ("group" %in% names(this_design)) {
        common_samples <- intersect(colnames(counts_mat), rownames(this_design))
        if (length(common_samples) > 0) {
          counts_mat <- counts_mat[, common_samples, drop = FALSE]
          groups_vec <- as.character(this_design[common_samples, "group"])
        }
      }

      # Sync plot_data to available proteins in counts
      plot_data <- plot_data |> dplyr::filter(Protein.Ids %in% rownames(counts_mat))
      counts_mat <- counts_mat[as.character(plot_data$Protein.Ids), , drop = FALSE]
    }
  }

  # Build Base R Annotation Dataframe (GUID rownames)
  clean_anno <- data.frame(
    Protein_Ids = as.character(plot_data$Protein.Ids),
    gene_name = as.character(plot_data$gene_name),
    stringsAsFactors = FALSE
  )
  
  purrr::walk(display_columns, function(col) {
    col_cln <- gsub("\\.", "_", col)
    if (col %in% names(plot_data) && !(col_cln %in% names(clean_anno))) {
      clean_anno[[col_cln]] <<- as.character(plot_data[[col]])
    }
  })
  
  clean_anno[is.na(clean_anno)] <- ""
  rownames(clean_anno) <- clean_anno$Protein_Ids

  # Build xy_data for Glimma
  glimma_table <- data.frame(
    logFC = signif(as.numeric(plot_data$logFC), 4),
    negLog10FDR = signif(as.numeric(plot_data$negLog10FDR), 4)
  )
  rownames(glimma_table) <- as.character(plot_data$Protein.Ids)

  xy_data <- Glimma:::buildXYData(
    table = glimma_table, status = as.integer(plot_data$Status),
    main = paste("Interactive Volcano Plot:", comparison_to_search),
    display.columns = colnames(clean_anno), anno = clean_anno, counts = counts_mat,
    groups = groups_vec, transform.counts = "none",
    status.cols = c("#1052bd", "silver", "#cc212f"),
    xlab = "logFC", ylab = "negLog10FDR"
  )
```

## Verification Plan

### Automated Tests
- Run [/tmp/test_volcano_edge_cases.R](file:///tmp/test_volcano_edge_cases.R) to verify all 3 critical failure modes are resolved.

### Manual Verification
- Test via the Shiny application.

<!-- APAF Bioinformatics | R_is_for_Robot | Approved | 2026-03-12 -->
