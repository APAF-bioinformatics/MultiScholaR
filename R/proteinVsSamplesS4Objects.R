##----------------------------------------------------------------------------------------------------------------------------------------------------------------------

#' Protein-Level Missing Value Imputation using limpa Package
#'
#' This function applies limpa's DPC-based missing value imputation directly to 
#' protein-level quantification data. This is useful when you already have protein
#' quantification data with missing values that need to be handled.
#'
#' @param theObject A ProteinQuantitativeData object with protein-level data
#' @param dpc_results DPC results to use. If NULL, will estimate using dpc_slope
#' @param dpc_slope Default DPC slope to use if no DPC results available (default: 0.8)
#' @param quantified_protein_column Name for the column containing quantified protein values
#' @param verbose Whether to print progress messages. Default is TRUE
#' @param chunk When verbose=TRUE, how often to output progress information (default: 1000)
#'
#' @details
#' This method treats each protein as a separate "feature" and applies DPC-based
#' imputation using limpa's dpcImpute function. This is appropriate when you have
#' protein-level data with missing values that follow intensity-dependent patterns.
#'
#' @return Updated ProteinQuantitativeData object with imputed protein values
#'
#' @importFrom limpa dpc dpcImpute
#' @export
setMethod(f="proteinMissingValueImputationLimpa"
          , signature="ProteinQuantitativeData"
          , definition = function(theObject, 
                                  dpc_results = NULL,
                                  dpc_slope = 0.8,
                                  quantified_protein_column = NULL,
                                  verbose = TRUE,
                                  chunk = 1000) {
            
            # Load required packages
            if (!requireNamespace("limpa", quietly = TRUE)) {
              stop("limpa package is required but not installed. Please install it using: BiocManager::install('limpa')")
            }
            
            # Parameter validation and defaults
            quantified_protein_column <- if (is.null(quantified_protein_column)) {
              "Protein.Imputed.Limpa"
            } else {
              quantified_protein_column
            }
            
            # Extract data from protein object
            protein_quant_table <- theObject@protein_quant_table
            protein_id_column <- theObject@protein_id_column
            sample_id_column <- theObject@sample_id
            design_matrix <- theObject@design_matrix
            
            if (verbose) {
              log_info("Starting limpa-based protein-level missing value imputation...")
            }
            
            # Convert to matrix format (proteins as rows, samples as columns)
            # First, identify sample columns (exclude protein ID and other metadata)
            sample_columns <- setdiff(colnames(protein_quant_table), 
                                     c(protein_id_column, "description", "gene_name", 
                                       "protein_name", "organism", "length"))
            
            if (verbose) {
              log_info("Converting protein data to matrix format...")
              log_info("Found {length(sample_columns)} sample columns")
            }
            
            # Create protein matrix
            protein_matrix <- protein_quant_table |>
              dplyr::select(all_of(c(protein_id_column, sample_columns))) |>
              tibble::column_to_rownames(protein_id_column) |>
              as.matrix()

            if (verbose) {
              log_info("Protein matrix dimensions: {nrow(protein_matrix)} proteins x {ncol(protein_matrix)} samples")
              log_info("Missing value percentage: {round(100 * mean(is.na(protein_matrix)), 1)}%")
            }
            
            # Check if we need log2 transformation
            # Assume if max value > 50, data is not log2 transformed
            max_val <- max(protein_matrix, na.rm = TRUE)
            needs_log_transform <- max_val > 50
            
            if (needs_log_transform) {
              if (verbose) {
                log_info("Converting to log2 scale for limpa (max value: {round(max_val, 2)})...")
              }
              protein_matrix <- log2(protein_matrix + 1)
            } else {
              if (verbose) {
                log_info("Data appears to be log2-scale already (max value: {round(max_val, 2)})")
              }
            }
            
            # Handle infinite or NaN values
            if (any(is.infinite(protein_matrix) | is.nan(protein_matrix), na.rm = TRUE)) {
              if (verbose) {
                log_warn("Infinite or NaN values detected. Replacing with NA...")
              }
              protein_matrix[is.infinite(protein_matrix) | is.nan(protein_matrix)] <- NA
            }
            
            # Get or estimate DPC parameters
            dpc_params <- NULL
            if (!is.null(dpc_results)) {
              dpc_params <- dpc_results
              if (verbose) {
                log_info("Using provided DPC results")
              }
            } else {
              if (verbose) {
                log_info("Estimating DPC from protein data...")
              }
            tryCatch({
                dpcfit <- limpa::dpc(protein_matrix)
                dpc_params <- dpcfit
                if (verbose) {
                  log_info("DPC parameters estimated:")
                  log_info("  beta0 (intercept): {round(dpcfit$dpc[1], 4)}")
                  log_info("  beta1 (slope): {round(dpcfit$dpc[2], 4)}")
                }
            }, error = function(e) {
                if (verbose) {
                  log_warn("DPC estimation failed, using default slope: {dpc_slope}")
                }
                dpc_params <- NULL
              })
            }
            
            # Apply missing value imputation
            if (verbose) {
              log_info("Applying protein-level missing value imputation...")
            }
            
              tryCatch({
              # Apply dpcImpute to protein matrix
              if (!is.null(dpc_params)) {
                imputed_result <- limpa::dpcImpute(protein_matrix, dpc = dpc_params, verbose = verbose, chunk = chunk)
            } else {
                imputed_result <- limpa::dpcImpute(protein_matrix, dpc.slope = dpc_slope, verbose = verbose, chunk = chunk)
              }
              
              if (verbose) {
                log_info("Protein-level imputation completed successfully")
                log_info("No missing values remaining: {!any(is.na(imputed_result$E))}")
              }
              
              # Extract imputed matrix
              imputed_matrix <- imputed_result$E
              
              # Transform back to original scale if necessary
              if (needs_log_transform) {
                if (verbose) {
                  log_info("Converting back from log2 scale...")
                }
                imputed_matrix <- 2^imputed_matrix - 1
                # Ensure no negative values
                imputed_matrix[imputed_matrix < 0] <- 0
              }
              
              # Convert back to long format
              if (verbose) {
                log_info("Converting imputed data back to original format...")
              }
              
              imputed_long <- imputed_matrix |>
                      as.data.frame() |>
                tibble::rownames_to_column(protein_id_column) |>
                tidyr::pivot_longer(cols = -all_of(protein_id_column), 
                                   names_to = sample_id_column, 
                                   values_to = quantified_protein_column)
              
              # Merge with original protein data
              updated_protein_data <- protein_quant_table |>
                dplyr::left_join(imputed_long, by = c(protein_id_column, sample_id_column))
              
              # Update the object
              theObject@protein_quant_table <- updated_protein_data
              
              # Store DPC results for future reference
              if (is.null(theObject@args)) {
                theObject@args <- list()
              }
              
              theObject@args$limpa_protein_imputation_results <- list(
                dpc_parameters_used = if (!is.null(dpc_params)) {
                  if (is.list(dpc_params) && !is.null(dpc_params$dpc)) {
                    dpc_params$dpc  # Extract parameters from DPC object
                  } else if (is.numeric(dpc_params)) {
                    dpc_params  # Already numeric parameters  
             } else {
                    c(NA, dpc_slope)
                  }
                } else {
                  c(NA, dpc_slope)
                },
                dpc_object_used = if (is.list(dpc_params) && !is.null(dpc_params$dpc)) dpc_params else NULL,
                quantified_protein_column = quantified_protein_column,
                missing_percentage_before = round(100 * mean(is.na(protein_matrix)), 1),
                missing_percentage_after = 0,  # DPC imputation produces complete data
                imputation_method = "limpa_dpc_protein_imputation",
                total_proteins_imputed = nrow(imputed_matrix)
              )
              
              if (verbose) {
                log_info("limpa protein-level imputation completed successfully!")
                log_info("New imputed column: {quantified_protein_column}")
                log_info("DPC results stored in object@args$limpa_protein_imputation_results")
              }
            
            return(theObject)
              
              }, error = function(e) {
              log_error(paste("Error during limpa protein imputation:", e$message))
              stop(paste("limpa protein imputation failed:", e$message))
            })
          })
