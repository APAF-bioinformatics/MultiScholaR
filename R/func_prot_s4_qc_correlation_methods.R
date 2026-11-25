#' @export
setMethod(
  f = "proteinTechRepCorrelation",
  signature = "ProteinQuantitativeData",
  definition = function(theObject, tech_rep_num_column = NULL, tech_rep_remove_regex = NULL) {
    protein_quant_table <- theObject@protein_quant_table
    protein_id_column <- theObject@protein_id_column
    design_matrix <- theObject@design_matrix
    sample_id <- theObject@sample_id
    tech_rep_column <- theObject@technical_replicate_id

    tech_rep_num_column <- checkParamsObjectFunctionSimplifyAcceptNull(theObject, "tech_rep_num_column", NULL)
    tech_rep_remove_regex <- checkParamsObjectFunctionSimplifyAcceptNull(theObject, "tech_rep_remove_regex", NULL)

    theObject <- updateParamInObject(theObject, "tech_rep_num_column")
    theObject <- updateParamInObject(theObject, "tech_rep_remove_regex")

    frozen_protein_matrix <- protein_quant_table |>
      column_to_rownames(protein_id_column) |>
      as.matrix()

    frozen_protein_matrix_pca <- frozen_protein_matrix
    frozen_protein_matrix_pca[!is.finite(frozen_protein_matrix_pca)] <- NA

    protein_matrix_tech_rep <- proteinTechRepCorrelationHelper(design_matrix, frozen_protein_matrix_pca,
      protein_id_column = protein_id_column,
      sample_id_column = sample_id,
      tech_rep_column = tech_rep_column,
      tech_rep_num_column = tech_rep_num_column,
      tech_rep_remove_regex = tech_rep_remove_regex
    )

    return(protein_matrix_tech_rep)
  }
)

## ----------------------------------------------------------------------------------------------------------------------------------------------------------------------
#' @title Plot Pearson Correlation
#' @param theObject is an object of the type ProteinQuantitativeData
#' @param tech_rep_remove_regex DEPRECATED - use exclude_pool_samples instead
#' @param correlation_group is the group where every pair of samples are compared
#' @param exclude_pool_samples Logical. If TRUE (default), automatically exclude samples from groups containing "Pool" or "QC" in their name from correlation analysis.
#' @export
setMethod(
  f = "plotPearson",
  signature = "ProteinQuantitativeData",
  definition = function(theObject, tech_rep_remove_regex = NULL, correlation_group = NA, exclude_pool_samples = TRUE) {
    correlation_group_to_use <- correlation_group

    if (is.na(correlation_group)) {
      correlation_group_to_use <- theObject@technical_replicate_id
    }

    # Handle deprecated parameter (backward compatibility)
    if (!is.null(tech_rep_remove_regex)) {
      message("*** plotPearson: WARNING - tech_rep_remove_regex is deprecated, use exclude_pool_samples instead ***")
    }

    correlation_vec <- pearsonCorForSamplePairs(theObject,
      tech_rep_remove_regex = tech_rep_remove_regex,
      correlation_group = correlation_group_to_use,
      exclude_pool_samples = exclude_pool_samples
    )

    # KIT-259 (A2): no within-group pair with a correlation gives a labelled figure, not an error.
    if (!any(is.finite(correlation_vec$pearson_correlation))) {
      return(buildProtPearsonUnavailablePlot(correlation_vec, correlation_group_to_use))
    }

    pearson_plot <- correlation_vec |>
      ggplot(aes(pearson_correlation)) +
      geom_histogram(breaks = seq(min(round(correlation_vec$pearson_correlation - 0.5, 2), na.rm = TRUE), 1, 0.001)) +
      protPearsonCountScale() +
      xlab("Pearson Correlation") +
      ylab("Counts") +
      theme(
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank()
      )

    return(pearson_plot)
  }
)

#' The count axis of the proteomics within-group Pearson histogram (KIT-259, A1): from 0 to the tallest bar, at least
#' 4, with whole-number ticks. The axis was fixed at 0 to 4, which removed every bar holding more than 4 pairs. A panel
#' whose bars hold at most 4 is drawn as before (ticks 0 to 4). The breaks function receives the expanded range, hence
#' the floor of its upper end. Shared by `plotPearson` and `rebuildWorkflowDiaQcPlot`, which must stay identical.
#' @return A ggplot2 continuous y scale.
#' @noRd
protPearsonCountScale <- function() {
  ggplot2::scale_y_continuous(
    limits = function(range) c(0, max(4, range[[2L]])),
    breaks = function(range) {
      top <- floor(range[[2L]])
      seq(0, top, by = max(1, ceiling(top / 4)))
    }
  )
}

#' The proteomics within-group Pearson histogram when no within-group pair has a correlation to draw (KIT-259, A2):
#' a labelled figure in place of the error at `seq(Inf, ...)`. It keeps the pairs as its data, so that the figure
#' rebuilt from the saved product says the same.
#' @param data The within-group pairs (column `pearson_correlation`), possibly none.
#' @param group The variable the pairs are grouped by ('Color by').
#' @return A ggplot with a text label and no axes.
#' @noRd
buildProtPearsonUnavailablePlot <- function(data, group) {
  if (!is.data.frame(data)) data <- data.frame()
  pairs <- if ("pearson_correlation" %in% names(data)) nrow(data) else 0L
  label <- if (pairs == 0L) {
    paste0("No two samples share a value of '", group, "' (Pool and QC groups are left out), so there are no ",
      "within-group pairs to correlate. Choose a 'Color by' with at least two samples per value.")
  } else {
    paste0("None of the ", pairs, " within-group sample pairs has a correlation (too few shared values, or a ",
      "sample with no variation), so there is nothing to draw.")
  }
  ggplot2::ggplot(data) +
    ggplot2::annotate("text", x = 0.5, y = 0.5, label = paste(strwrap(label, 70L), collapse = "\n"), size = 5) +
    ggplot2::xlim(0, 1) + ggplot2::ylim(0, 1) +
    ggplot2::theme_void()
}

#' @title Pearson Correlation for Sample Pairs
#' @param theObject is an object of the type ProteinQuantitativeData
#' @param tech_rep_remove_regex DEPRECATED - use exclude_pool_samples instead
#' @param correlation_group is the group where every pair of samples are compared
#' @param exclude_pool_samples Logical. If TRUE (default), automatically exclude samples from groups containing "Pool" or "QC" in their name from correlation analysis. Pool/QC samples are excluded from within-group correlation calculations but remain in RUV-III correction.
#' @export
setMethod(
  f = "pearsonCorForSamplePairs",
  signature = "ProteinQuantitativeData",
  definition = function(theObject, tech_rep_remove_regex = NULL, correlation_group = NA, exclude_pool_samples = TRUE) {

    # Memory tracking - Entry
    entry_mem <- checkMemoryBoth("Entry", context = "pearsonCorForSamplePairs")

    protein_quant_table <- theObject@protein_quant_table
    protein_id_column <- theObject@protein_id_column
    design_matrix <- theObject@design_matrix
    sample_id <- theObject@sample_id

    replicate_group_column <- theObject@technical_replicate_id
    if (!is.na(correlation_group)) {
      replicate_group_column <- correlation_group
    }

    # Handle deprecated parameter (backward compatibility)
    if (!is.null(tech_rep_remove_regex)) {
      message("*** PEARSON: WARNING - tech_rep_remove_regex is deprecated, use exclude_pool_samples instead ***")
    }

    exclude_pool_samples <- checkParamsObjectFunctionSimplifyAcceptNull(theObject, "exclude_pool_samples", TRUE)
    theObject <- updateParamInObject(theObject, "exclude_pool_samples")

    # --- OPTIMIZED MATRIX APPROACH ---

    # 1. Extract numeric matrix
    # Get sample columns (exclude ID)
    sample_cols <- setdiff(colnames(protein_quant_table), protein_id_column)

    # Ensure they match design matrix samples
    # Note: protein_quant_table might contain pools/controls not in design if not cleaned, but usually they should align.
    valid_samples <- intersect(sample_cols, design_matrix[[sample_id]])

    if (length(valid_samples) == 0) {
      stop("No matching samples found between protein data and design matrix")
    }

    # Subset and convert to matrix
    mat <- as.matrix(protein_quant_table[, valid_samples, drop = FALSE])

    message(sprintf("   [pearsonCorForSamplePairs] Matrix dimensions: %d proteins x %d samples", nrow(mat), ncol(mat)))

    # 2. Calculate Correlation Matrix
    # fast vectorized correlation
    start_time <- Sys.time()
    cor_mat <- cor(mat, method = "pearson", use = "pairwise.complete.obs")
    end_time <- Sys.time()
    message(sprintf("   [pearsonCorForSamplePairs] Correlation matrix computed in %.2f seconds", as.numeric(difftime(end_time, start_time, units = "secs"))))

    # 3. Convert to Long Format (Pairwise) and Filter

    # Map samples to groups using design matrix
    # Create a named vector for fast lookup
    sample_to_group <- setNames(as.character(design_matrix[[replicate_group_column]]), as.character(design_matrix[[sample_id]]))

    # Get upper triangle indices to avoid duplicates and self-correlation
    upper_tri_indices <- which(upper.tri(cor_mat), arr.ind = TRUE)

    # Extract names and values
    samples1 <- rownames(cor_mat)[upper_tri_indices[, 1]]
    samples2 <- colnames(cor_mat)[upper_tri_indices[, 2]]
    cor_values <- cor_mat[upper_tri_indices]

    # Define column names expected by filterSamplesByProteinCorrelationThreshold
    col_x <- paste0(sample_id, ".x")
    col_y <- paste0(sample_id, ".y")

    # Create dataframe
    # Pre-allocate list for speed then convert
    result_list <- list()
    result_list[[col_x]] <- samples1
    result_list[[col_y]] <- samples2
    result_list[["pearson_correlation"]] <- cor_values

    result_df <- as.data.frame(result_list, stringsAsFactors = FALSE)

    # Add groups
    result_df$Group_1 <- sample_to_group[result_df[[col_x]]]
    result_df$Group_2 <- sample_to_group[result_df[[col_y]]]

    # Filter for within-group correlations only (Group_1 == Group_2)
    # Handle NA groups just in case
    result_df <- result_df[!is.na(result_df$Group_1) & !is.na(result_df$Group_2) & result_df$Group_1 == result_df$Group_2, ]

    # Add the group column required by the S4 object / downstream logic
    # The helper seems to expect the group info might be joined or passed, but let's see.
    # The return value is just the correlation table.
    # However, for `exclude_pool_samples`, we need the group.

    result_df[[replicate_group_column]] <- result_df$Group_1

    # Clean up temporary group columns
    result_df$Group_1 <- NULL
    result_df$Group_2 <- NULL

    # 4. Exclude Pool Samples Logic
    if (exclude_pool_samples) {
      # Detect Pool/QC groups (case-insensitive detection for "pool" or "qc")
      is_pool_qc_group <- grepl("pool|qc", result_df[[replicate_group_column]], ignore.case = TRUE)
      pool_qc_count <- sum(is_pool_qc_group)

      if (pool_qc_count > 0) {
        message(sprintf("*** PEARSON: Excluded %d sample pairs from Pool/QC groups ***", pool_qc_count))
        result_df <- result_df[!is_pool_qc_group, ]
      } else {
        message("*** PEARSON: No Pool/QC groups detected to exclude ***")
      }
    }

    message(sprintf("   [pearsonCorForSamplePairs] Final pair count: %d", nrow(result_df)))

    # Memory tracking - Exit
    reportMemoryDelta(entry_mem, "TOTAL pearsonCorForSamplePairs", context = "pearsonCorForSamplePairs")


    # Ensure the returned dataframe has the expected columns for the next step
    # filterSamplesByProteinCorrelationThreshold expects:
    # - {sample_id}.x
    # - {sample_id}.y
    # - pearson_correlation
    # And it might carry over replicate_group_column (it was in the original helper output)

    return(result_df)
  }
)

