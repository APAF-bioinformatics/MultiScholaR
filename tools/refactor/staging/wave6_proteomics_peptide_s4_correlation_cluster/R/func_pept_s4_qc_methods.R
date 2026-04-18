#' @export
setMethod(
  f = "plotPearson",
  signature = "PeptideQuantitativeData",
  definition = function(theObject, tech_rep_remove_regex = NULL, correlation_group = NA, exclude_pool_samples = TRUE) {
    correlation_vec <- pearsonCorForSamplePairs(theObject,
      tech_rep_remove_regex = tech_rep_remove_regex,
      correlation_group = correlation_group,
      exclude_pool_samples = exclude_pool_samples
    )

    if (nrow(correlation_vec) == 0) {
      return(ggplot() + labs(title = "No within-group sample pairs found for correlation"))
    }

    pearson_plot <- correlation_vec |>
      ggplot(aes(pearson_correlation)) +
      geom_histogram(bins = 100) +
      xlab("Pearson Correlation") +
      ylab("Counts") +
      apafTheme()
      
    return(pearson_plot)
  }
)

#' @export
setMethod(
  f = "pearsonCorForSamplePairs",
  signature = "PeptideQuantitativeData",
  definition = function(theObject, tech_rep_remove_regex = NULL, correlation_group = NA, exclude_pool_samples = TRUE) {
    message("+===========================================================================+")
    message("|  DEBUG66: Entering pearsonCorForSamplePairs (Peptides)                   |")
    message("+===========================================================================+")

    design_matrix <- theObject@design_matrix
    sample_id <- theObject@sample_id
    peptide_matrix <- theObject@peptide_matrix

    replicate_group_column <- theObject@technical_replicate_id
    if (!is.na(correlation_group)) {
      replicate_group_column <- correlation_group
    }

    # Use the matrix directly
    mat <- peptide_matrix
    
    # Ensure they match design matrix samples
    valid_samples <- intersect(colnames(mat), design_matrix[[sample_id]])
    if (length(valid_samples) == 0) {
      stop("No matching samples found between peptide data and design matrix")
    }
    mat <- mat[, valid_samples, drop = FALSE]

    # Calculate Correlation Matrix
    cor_mat <- cor(mat, method = "pearson", use = "pairwise.complete.obs")

    # Map samples to groups
    sample_to_group <- setNames(as.character(design_matrix[[replicate_group_column]]), as.character(design_matrix[[sample_id]]))

    # Get upper triangle indices
    upper_tri_indices <- which(upper.tri(cor_mat), arr.ind = TRUE)
    samples1 <- rownames(cor_mat)[upper_tri_indices[, 1]]
    samples2 <- colnames(cor_mat)[upper_tri_indices[, 2]]
    cor_values <- cor_mat[upper_tri_indices]

    col_x <- paste0(sample_id, ".x")
    col_y <- paste0(sample_id, ".y")

    result_df <- data.frame(
      temp_s1 = samples1,
      temp_s2 = samples2,
      pearson_correlation = cor_values,
      stringsAsFactors = FALSE
    )
    colnames(result_df)[1:2] <- c(col_x, col_y)

    # Add groups
    result_df$Group_1 <- sample_to_group[result_df[[col_x]]]
    result_df$Group_2 <- sample_to_group[result_df[[col_y]]]

    # Filter for within-group correlations
    result_df <- result_df[!is.na(result_df$Group_1) & !is.na(result_df$Group_2) & result_df$Group_1 == result_df$Group_2, ]
    result_df[[replicate_group_column]] <- result_df$Group_1

    # Clean up
    result_df$Group_1 <- NULL
    result_df$Group_2 <- NULL

    if (exclude_pool_samples) {
      is_pool_qc_group <- grepl("pool|qc", result_df[[replicate_group_column]], ignore.case = TRUE)
      if (sum(is_pool_qc_group) > 0) {
        result_df <- result_df[!is_pool_qc_group, ]
      }
    }

    message("|  DEBUG66: Exiting pearsonCorForSamplePairs (Peptides)                    |")
    return(result_df)
  }
)

