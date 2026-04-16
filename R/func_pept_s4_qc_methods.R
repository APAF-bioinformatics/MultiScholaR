#'@export
#' @export
setMethod(f="plotRle"
          , signature="PeptideQuantitativeData"
          , definition=function(theObject, grouping_variable, yaxis_limit = c(), sample_label = NULL) {
            peptide_matrix <- theObject@peptide_matrix
            design_matrix <- theObject@design_matrix
            sample_id <- theObject@sample_id

            design_matrix <- as.data.frame(design_matrix)

            if(!is.null(sample_label)) {
              if (sample_label %in% colnames(design_matrix)) {
                rownames(design_matrix) <- design_matrix[,sample_label]
                colnames(peptide_matrix) <- design_matrix[,sample_label]
              }
            } else {
              rownames(design_matrix) <- design_matrix[,sample_id]
            }

            rowinfo_vector <- NA
            if(!is.na(grouping_variable)){
              rowinfo_vector <- design_matrix[colnames(peptide_matrix), grouping_variable]
            }

            # Handle missing/non-finite values
            working_matrix <- peptide_matrix
            working_matrix[!is.finite(working_matrix)] <- NA

            rle_plot <- plotRleHelper(t(working_matrix)
                                     , rowinfo = rowinfo_vector
                                     , yaxis_limit = yaxis_limit)

            return(rle_plot)
          })

# Renamed from plotPca to avoid S4 generic conflict
plotPcaDispatch <- function(theObject, grouping_variable, shape_variable = NULL, label_column, title, font_size=8, cv_percentile = 0.90) {
  # Defensive checks
  if (!is.character(grouping_variable) || length(grouping_variable) != 1) {
    stop("grouping_variable must be a single character string")
  }

  if (!is.null(shape_variable) && (!is.character(shape_variable) || length(shape_variable) != 1)) {
    stop("shape_variable must be NULL or a single character string")
  }

  if (!grouping_variable %in% colnames(theObject@design_matrix)) {
    stop(sprintf("grouping_variable '%s' not found in design matrix", grouping_variable))
  }

  if (!is.null(shape_variable) && !shape_variable %in% colnames(theObject@design_matrix)) {
    stop(sprintf("shape_variable '%s' not found in design matrix", shape_variable))
  }

  data_matrix <- NULL
  ## I want to check the class of theObject here
 if( class(theObject) == "PeptideQuantitativeData") {
   data_matrix <- theObject@peptide_matrix

 } else if( class(theObject) == "ProteinQuantitativeData") {
   data_matrix <- theObject@protein_quant_table |>
     column_to_rownames(var = "Protein.Ids") |>
     as.matrix()
 }

  design_matrix <- theObject@design_matrix
  sample_id <- theObject@sample_id

  # Prepare matrix for PCA (data should already be log2 transformed)
  data_matrix_pca <- data_matrix
  data_matrix_pca[!is.finite(data_matrix_pca)] <- NA

  if(is.na(label_column) || label_column == "") {
    label_column <- ""
  }

  pca_plot <- plotPcaHelper(data_matrix_pca
                            , design_matrix = design_matrix
                            , sample_id_column = sample_id
                            , grouping_variable = grouping_variable
                            , shape_variable = shape_variable
                            , label_column = label_column
                            , title = title
                            , geom.text.size = font_size
                            , cv_percentile = cv_percentile)

  return(pca_plot)
}

#'@export
setMethod(f="plotPca"
          , signature="PeptideQuantitativeData"
          , definition=function(theObject, grouping_variable, shape_variable = NULL, label_column, title, font_size=8, cv_percentile = 0.90) {
            # Defensive checks
            if (!is.character(grouping_variable) || length(grouping_variable) != 1) {
              stop("grouping_variable must be a single character string")
            }

            if (!is.null(shape_variable) && (!is.character(shape_variable) || length(shape_variable) != 1)) {
              stop("shape_variable must be NULL or a single character string")
            }

            if (!grouping_variable %in% colnames(theObject@design_matrix)) {
              stop(sprintf("grouping_variable '%s' not found in design matrix", grouping_variable))
            }

            if (!is.null(shape_variable) && !shape_variable %in% colnames(theObject@design_matrix)) {
              stop(sprintf("shape_variable '%s' not found in design matrix", shape_variable))
            }

            peptide_matrix <- theObject@peptide_matrix
            design_matrix <- theObject@design_matrix
            sample_id <- theObject@sample_id

            # Prepare matrix for PCA (data should already be log2 transformed)
            peptide_matrix_pca <- peptide_matrix
            peptide_matrix_pca[!is.finite(peptide_matrix_pca)] <- NA

            if(is.na(label_column) || label_column == "") {
              label_column <- ""
            }

            pca_plot <- plotPcaHelper(peptide_matrix_pca
                                     , design_matrix = design_matrix
                                     , sample_id_column = sample_id
                                     , grouping_variable = grouping_variable
                                     , shape_variable = shape_variable
                                     , label_column = label_column
                                     , title = title
                                     , geom.text.size = font_size
                                     , cv_percentile = cv_percentile)

            return(pca_plot)
          })

# PCA plot for PeptideQuantitativeData

#' @export
setMethod(
  f = "plotPca",
  signature = "PeptideQuantitativeData",
  definition = function(theObject, grouping_variable, shape_variable = NULL, label_column = NULL, title = NULL, font_size = 8, cv_percentile = 0.90) {
    # Defensive checks
    if (!is.character(grouping_variable) || length(grouping_variable) != 1) {
      stop("grouping_variable must be a single character string")
    }

    if (!is.null(shape_variable) && (!is.character(shape_variable) || length(shape_variable) != 1)) {
      stop("shape_variable must be NULL or a single character string")
    }

    if (!grouping_variable %in% colnames(theObject@design_matrix)) {
      stop(sprintf("grouping_variable '%s' not found in design matrix", grouping_variable))
    }

    if (!is.null(shape_variable) && !shape_variable %in% colnames(theObject@design_matrix)) {
      stop(sprintf("shape_variable '%s' not found in design matrix", shape_variable))
    }

    # Extract data from the matrix slot
    peptide_matrix <- theObject@peptide_matrix
    design_matrix <- theObject@design_matrix
    sample_id <- theObject@sample_id

    # Handle missing values - replace INF/NAN with NA
    working_matrix <- peptide_matrix
    working_matrix[!is.finite(working_matrix)] <- NA

    if (is.na(label_column) || label_column == "") {
      label_column <- ""
    }

    required_cols <- c(sample_id, grouping_variable)
    if (!is.null(shape_variable)) {
      required_cols <- c(required_cols, shape_variable)
    }

    tryCatch(
      {
        pca_plot <- plotPcaHelper(working_matrix,
          design_matrix,
          sample_id_column = sample_id,
          grouping_variable = grouping_variable,
          shape_variable = shape_variable,
          label_column = label_column,
          title = title,
          geom.text.size = font_size
        )
        return(pca_plot)
      },
      error = function(e) {
        stop(sprintf("Error in plotPcaHelper for PeptideQuantitativeData: %s", e$message))
      }
    )
  }
)

## ----------------------------------------------------------------------------------------------------------------------------------------------------------------------

# Density plot for PeptideQuantitativeData

#' @export
setMethod(
  f = "plotDensity",
  signature = "PeptideQuantitativeData",
  definition = function(theObject, grouping_variable, title = "", font_size = 8) {
    peptide_data <- theObject@peptide_data
    quant_col <- theObject@norm_quantity_column
    sample_id <- theObject@sample_id

    peptide_data |>
      filter(!is.na(!!sym(quant_col))) |>
      ggplot(aes(x = !!sym(quant_col), color = !!sym(sample_id))) +
      geom_density() +
      apafTheme() +
      labs(title = title, x = "Log2 Intensity", y = "Density") +
      theme(legend.position = "none")
  }
)

## ----------------------------------------------------------------------------------------------------------------------------------------------------------------------

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

setMethod(f="plotDensity"
          , signature="ggplot2::ggplot"
          , definition=function(theObject, grouping_variable, title = "", font_size = 8) {
            # First try to get data directly from the ggplot object's data element
            if (!is.null(theObject$data) && is.data.frame(theObject$data)) {
              pca_data <- as_tibble(theObject$data)
            } else {
              # Fall back to other extraction methods
              pca_data <- as_tibble(ggplot_build(theObject)$data[[1]])

              # If the data doesn't have PC1/PC2, try to extract from the plot's environment
              if (!("PC1" %in% colnames(pca_data) && "PC2" %in% colnames(pca_data))) {
                # Try to get the data from the plot's environment
                if (exists("data", envir = environment(theObject$mapping$x))) {
                  pca_data <- as_tibble(get("data", envir = environment(theObject$mapping$x)))
                } else {
                  stop("Could not extract PCA data from the ggplot object")
                }
              }
            }

            # Check if grouping variable exists in the data
            if (!grouping_variable %in% colnames(pca_data)) {
              stop(sprintf("grouping_variable '%s' not found in the data", grouping_variable))
            }

            # Create PC1 density plot
            pc1_density <- ggplot(pca_data, aes(x = PC1, fill = !!sym(grouping_variable), color = !!sym(grouping_variable))) +
              geom_density(alpha = 0.3) +
              theme_bw() +
              labs(title = title,
                   x = "PC1",
                   y = "Density") +
              theme(
                text = element_text(size = font_size),
                plot.margin = margin(b = 0, t = 5, l = 5, r = 5),
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                panel.background = element_blank()
              )

            # Create PC2 density plot
            pc2_density <- ggplot(pca_data, aes(x = PC2, fill = !!sym(grouping_variable), color = !!sym(grouping_variable))) +
              geom_density(alpha = 0.3) +
              theme_bw() +
              labs(x = "PC2",
                   y = "Density") +
              theme(
                text = element_text(size = font_size),
                plot.margin = margin(t = 0, b = 5, l = 5, r = 5),
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                panel.background = element_blank()
              )

            # Combine plots with minimal spacing
            combined_plot <- pc1_density / pc2_density +
              plot_layout(heights = c(1, 1)) +
              plot_annotation(theme = theme(plot.margin = margin(0, 0, 0, 0)))

            return(combined_plot)
          })
