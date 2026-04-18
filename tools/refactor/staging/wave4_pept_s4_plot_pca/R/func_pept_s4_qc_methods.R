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

