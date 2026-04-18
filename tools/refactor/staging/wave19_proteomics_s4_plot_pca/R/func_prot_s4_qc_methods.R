#' @export
setMethod(
  f = "plotPca",
  signature = "ProteinQuantitativeData",
  definition = function(theObject, grouping_variable, shape_variable = NULL, label_column, title, font_size = 8, cv_percentile = 0.90) {
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

    protein_quant_table <- theObject@protein_quant_table
    protein_id_column <- theObject@protein_id_column
    design_matrix <- theObject@design_matrix
    sample_id <- theObject@sample_id

    frozen_protein_matrix <- protein_quant_table |>
      column_to_rownames(protein_id_column) |>
      as.matrix()

    frozen_protein_matrix_pca <- frozen_protein_matrix
    frozen_protein_matrix_pca[!is.finite(frozen_protein_matrix_pca)] <- NA

    if (is.na(label_column) || label_column == "") {
      label_column <- ""
    }

    required_cols <- c(sample_id, grouping_variable)
    if (!is.null(shape_variable)) {
      required_cols <- c(required_cols, shape_variable)
    }
    missing_cols <- setdiff(required_cols, colnames(design_matrix))
    if (length(missing_cols) > 0) {
      stop(sprintf("Missing columns in design matrix: %s", paste(missing_cols, collapse = ", ")))
    }

    tryCatch(
      {
        pca_plot <- plotPcaHelper(frozen_protein_matrix_pca,
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
        stop(sprintf("Error in plotPcaHelper: %s", e$message))
      }
    )
  }
)

