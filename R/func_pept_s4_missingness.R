#' @export
setMethod(f="peptideMissingValueImputationLimpa"
          , signature="PeptideQuantitativeData"
          , definition = function(theObject,
                                  imputed_value_column = NULL,
                                  use_log2_transform = TRUE,
                                  verbose = TRUE,
                                  ensure_matrix = TRUE) {

            # Load required packages
            if (!requireNamespace("limpa", quietly = TRUE)) {
              stop("limpa package is required but not installed. Please install it using: BiocManager::install('limpa')")
            }

            # Parameter validation and defaults
            imputed_value_column <- checkParamsObjectFunctionSimplifyAcceptNull(
              theObject, "imputed_value_column", "Peptide.Imputed.Limpa"
            )

            use_log2_transform <- checkParamsObjectFunctionSimplify(
              theObject, "use_log2_transform", TRUE
            )

            verbose <- checkParamsObjectFunctionSimplify(
              theObject, "verbose", TRUE
            )

            # Update parameters in object
            theObject <- updateParamInObject(theObject, "imputed_value_column")
            theObject <- updateParamInObject(theObject, "use_log2_transform")
            theObject <- updateParamInObject(theObject, "verbose")

            # Ensure peptide matrix is calculated if requested
            if (ensure_matrix && (!"peptide_matrix" %in% slotNames(theObject) ||
                                  is.null(theObject@peptide_matrix) ||
                                  length(theObject@peptide_matrix) == 0)) {
              if (verbose) {
                log_info("Peptide matrix not found. Calculating peptide matrix...")
              }
              theObject <- calcPeptideMatrix(theObject)
            }

            # Extract data
            peptide_data <- theObject@peptide_data
            peptide_matrix <- theObject@peptide_matrix
            raw_quantity_column <- theObject@raw_quantity_column
            sample_id_column <- theObject@sample_id
            design_matrix <- theObject@design_matrix

            if (verbose) {
              log_info("Starting limpa-based missing value imputation...")
              log_info("Data dimensions: {nrow(peptide_matrix)} peptides x {ncol(peptide_matrix)} samples")
              log_info("Missing value percentage: {round(100 * mean(is.na(peptide_matrix)), 1)}%")
            }

            # Prepare data for limpa (peptides as rows, samples as columns)
            # limpa expects log2-transformed data
            y_peptide <- peptide_matrix

            # Transform to log2 if requested and data is not already log-transformed
            if (use_log2_transform && !theObject@is_logged_data) {
              if (verbose) {
                log_info("Applying log2 transformation...")
              }
              # Add small constant to avoid log(0)
              y_peptide <- log2(y_peptide + 1)
            } else if (use_log2_transform && theObject@is_logged_data) {
              if (verbose) {
                log_warn("Data already log2 transformed, skipping additional transformation")
              }
              # Data already log2, use as-is
            } else if (!use_log2_transform && !theObject@is_logged_data) {
              if (verbose) {
                log_info("Converting raw intensities to log2 scale for limpa...")
              }
              # limpa expects log2 data, so transform raw data
              y_peptide <- log2(y_peptide + 1)
            } else {
              # !use_log2_transform && theObject@is_logged_data
              if (verbose) {
                log_info("Using existing log2 transformed data (no additional transformation)")
              }
              # Data already log2, use as-is - this is the correct case!
            }

            # Check for infinite or NaN values
            if (any(is.infinite(y_peptide) | is.nan(y_peptide), na.rm = TRUE)) {
              if (verbose) {
                log_warn("Infinite or NaN values detected. Replacing with NA...")
              }
              y_peptide[is.infinite(y_peptide) | is.nan(y_peptide)] <- NA
            }

            # Estimate Detection Probability Curve
            if (verbose) {
              log_info("Estimating detection probability curve...")
            }

            tryCatch({
              dpcfit <- limpa::dpc(y_peptide)

              if (verbose) {
                log_info("DPC parameters estimated:")
                log_info("  beta0 (intercept): {round(dpcfit$dpc[1], 4)}")
                log_info("  beta1 (slope): {round(dpcfit$dpc[2], 4)}")

                # Interpret the slope
                slope_interpretation <- if (dpcfit$dpc[2] < 0.3) {
                  "nearly random missing"
                } else if (dpcfit$dpc[2] < 0.7) {
                  "moderate intensity-dependent missing"
                } else if (dpcfit$dpc[2] < 1.2) {
                  "strong intensity-dependent missing"
                } else {
                  "very strong intensity-dependent missing (approaching left-censoring)"
                }
                log_info("  Interpretation: {slope_interpretation}")
              }

              # Perform row-wise imputation using limpa
              if (verbose) {
                log_info("Performing row-wise imputation using DPC model...")
              }

              y_imputed <- limpa::dpcImpute(y_peptide, dpc = dpcfit)

              if (verbose) {
                log_info("Imputation completed successfully")
                log_info("No missing values remaining: {!any(is.na(y_imputed$E))}")
              }

              # Extract the imputed matrix
              imputed_matrix <- y_imputed$E

              # Transform back to original scale if necessary
              if (use_log2_transform && !theObject@is_logged_data) {
                if (verbose) {
                  log_info("Converting back from log2 scale...")
                }
                imputed_matrix <- 2^imputed_matrix - 1
                # Ensure no negative values
                imputed_matrix[imputed_matrix < 0] <- 0
              }

              # Convert back to long format and merge with original data
              if (verbose) {
                log_info("Converting imputed data back to original format...")
              }

              # Create peptide IDs that match the matrix rownames
              peptide_ids <- rownames(imputed_matrix)

              # Convert imputed matrix to long format
              imputed_long <- imputed_matrix |>
                as.data.frame() |>
                tibble::rownames_to_column("peptide_id") |>
                tidyr::pivot_longer(cols = -peptide_id,
                                   names_to = sample_id_column,
                                   values_to = imputed_value_column) |>
                tidyr::separate(peptide_id,
                               into = c(theObject@protein_id_column, theObject@peptide_sequence_column),
                               sep = "%")

              # Merge with original peptide data
              updated_peptide_data <- peptide_data |>
                dplyr::left_join(imputed_long,
                                by = c(theObject@protein_id_column,
                                      theObject@peptide_sequence_column,
                                      sample_id_column))

              # Update the object
              theObject@peptide_data <- updated_peptide_data
              theObject@peptide_matrix <- imputed_matrix

              # Update norm_quantity_column to point to the new imputed column
              # This ensures plotting functions use the final imputed data
              theObject@norm_quantity_column <- imputed_value_column

              # Store DPC results in the object for future reference
              if (is.null(theObject@args)) {
                theObject@args <- list()
              }
              theObject@args$limpa_dpc_results <- list(
                dpc_parameters = dpcfit$dpc,  # Numeric vector c(intercept, slope)
                dpc_object = dpcfit,          # Full DPC object (preferred for dpcQuant)
                missing_percentage_before = round(100 * mean(is.na(y_peptide)), 1),
                missing_percentage_after = round(100 * mean(is.na(imputed_matrix)), 1),
                slope_interpretation = slope_interpretation,
                dpc_method = "limpa_dpc",
                # Store the original y_peptide data for recreating DPC plot
                y_peptide_for_dpc = y_peptide
              )

              # Clean design matrix
              theObject <- cleanDesignMatrixPeptide(theObject)

              if (verbose) {
                log_info("limpa-based imputation completed successfully!")
                log_info("New imputed column: {imputed_value_column}")
                log_info("DPC parameters stored in object@args$limpa_dpc_results")
              }

              return(theObject)

            }, error = function(e) {
              log_error("Error during limpa imputation: {e$message}")
              stop(paste("limpa imputation failed:", e$message))
            })
          })
