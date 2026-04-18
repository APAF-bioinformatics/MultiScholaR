## normalise between Arrays
#'@export
#'@param theObject Object of class PeptideQuantitativeData
#'@param normalisation_method Method to use for normalisation. Options are cyclicloess, quantile, scale, none
#' @title Normalize Between Samples for PeptideQuantitativeData
#' @name normaliseBetweenSamples,PeptideQuantitativeData-method
#' @export
setMethod(f="normaliseBetweenSamples"
          , signature="PeptideQuantitativeData"
          , definition=function(theObject, normalisation_method = NULL) {
            peptide_data <- theObject@peptide_data
            design_matrix <- theObject@design_matrix
            sample_id <- theObject@sample_id
            current_quant_column <- theObject@norm_quantity_column

            normalisation_method <- checkParamsObjectFunctionSimplify(theObject
                                                                      , "normalisation_method"
                                                                      , "cyclicloess")

            theObject <- updateParamInObject(theObject, "normalisation_method")

            # Create matrix from peptide_data (like protein version from protein_quant_table)
            # Create unique peptide identifier for matrix rownames
            peptide_unique_id <- paste(peptide_data[[theObject@protein_id_column]], 
                                      peptide_data[[theObject@peptide_sequence_column]], 
                                      sep = "%")
            
            # Create wide format data frame then convert to matrix (exactly like protein version)
            temp_peptide_wide <- peptide_data |>
              mutate(peptide_id = peptide_unique_id) |>
              select(peptide_id, !!sym(sample_id), !!sym(current_quant_column)) |>
              pivot_wider(names_from = !!sym(sample_id), 
                         values_from = !!sym(current_quant_column),
                         values_fn = mean)

            # Convert to matrix exactly like protein version: column_to_rownames() then as.matrix()
            frozen_peptide_matrix <- temp_peptide_wide |>
              column_to_rownames("peptide_id") |>
              as.matrix()

            frozen_peptide_matrix[!is.finite(frozen_peptide_matrix)] <- NA

            normalised_frozen_peptide_matrix <- frozen_peptide_matrix

            print(paste0("normalisation_method = ", normalisation_method))

            switch(normalisation_method
                   , cyclicloess = {
                     normalised_frozen_peptide_matrix <- limma::normalizeCyclicLoess(frozen_peptide_matrix)
                   }
                   , quantile = {
                     normalised_frozen_peptide_matrix <- limma::normalizeQuantiles(frozen_peptide_matrix)
                   }
                   , scale = {
                     normalised_frozen_peptide_matrix <- limma::normalizeMedianAbsValues(frozen_peptide_matrix)
                   }
                   , none = {
                     normalised_frozen_peptide_matrix <- frozen_peptide_matrix
                   }
            )

            normalised_frozen_peptide_matrix[!is.finite(normalised_frozen_peptide_matrix)] <- NA

            # Convert back to data frame (exactly like protein version)
            normalised_peptide_table <- normalised_frozen_peptide_matrix |>
              as.data.frame() |>
              rownames_to_column("peptide_id")

            # Update peptide_data by joining with normalized values and replacing the quant column
            updated_peptide_data <- peptide_data |>
              mutate(peptide_id = peptide_unique_id) |>
              select(-!!sym(current_quant_column)) |>
              left_join(normalised_peptide_table |> 
                       pivot_longer(cols = -peptide_id, 
                                   names_to = sample_id, 
                                   values_to = current_quant_column),
                       by = c("peptide_id", sample_id)) |>
              select(-peptide_id)

            # Update both slots
            theObject@peptide_data <- updated_peptide_data
            theObject@peptide_matrix <- normalised_frozen_peptide_matrix

            theObject <- cleanDesignMatrixPeptide(theObject)

            return(theObject)
          })

#' Log2 Transform Peptide Matrix
#'
#' Transforms raw peptide intensity values to log2 scale for downstream normalization and RUV analysis.
#' This should be called after calcPeptideMatrix() and before normaliseBetweenSamples().
#' @export
#'
#' @param theObject A PeptideQuantitativeData object
#' @return PeptideQuantitativeData object with log2 transformed data
#' @export
setMethod(f="log2TransformPeptideMatrix"
          , signature="PeptideQuantitativeData"
          , definition=function(theObject) {
            
            if (theObject@is_logged_data) {
              warning("Data appears to already be log-transformed (is_logged_data = TRUE). Skipping transformation.")
              return(theObject)
            }
            
            peptide_data <- theObject@peptide_data
            peptide_matrix <- theObject@peptide_matrix
            current_quant_column <- theObject@norm_quantity_column
            
            # Log2 transform the matrix
            log2_peptide_matrix <- peptide_matrix
            
            # Handle zeros and negative values
            log2_peptide_matrix[log2_peptide_matrix <= 0] <- NA
            log2_peptide_matrix <- log2(log2_peptide_matrix)
            
            # Update matrix
            theObject@peptide_matrix <- log2_peptide_matrix
            
            # Also update peptide_data to maintain consistency
            # Create long format from log2 matrix
            log2_long <- log2_peptide_matrix |>
              as.data.frame() |>
              rownames_to_column("peptide_row_id") |>
              pivot_longer(cols = -peptide_row_id, 
                          names_to = theObject@sample_id, 
                          values_to = "log2_value")

            # Update peptide_data: match by protein%peptide ID and sample
            updated_peptide_data <- peptide_data |>
              mutate(peptide_row_id = paste(!!sym(theObject@protein_id_column), 
                                           !!sym(theObject@peptide_sequence_column), 
                                           sep = "%")) |>
              left_join(log2_long, by = c("peptide_row_id", theObject@sample_id)) |>
              mutate(!!sym(current_quant_column) := ifelse(!is.na(log2_value), 
                                                           log2_value, 
                                                           !!sym(current_quant_column))) |>
              select(-peptide_row_id, -log2_value)

            theObject@peptide_data <- updated_peptide_data
            
            # Mark as logged
            theObject@is_logged_data <- TRUE
            
            theObject <- cleanDesignMatrixPeptide(theObject)
            
            message("Peptide data successfully log2 transformed. Raw intensities converted to log2 scale.")
            message(paste("is_logged_data flag set to:", theObject@is_logged_data))
            
            return(theObject)
          })

