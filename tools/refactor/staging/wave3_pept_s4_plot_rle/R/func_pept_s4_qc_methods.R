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

