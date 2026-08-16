# ----------------------------------------------------------------------------
# peptideMissingValueImputation
# ----------------------------------------------------------------------------
#'@export
setMethod( f="peptideMissingValueImputation"
           , signature="PeptideQuantitativeData"
           , definition = function(theObject,
                                   imputed_value_column = NULL,
                                   proportion_missing_values = NULL,
                                   core_utilisation = NULL,
                                   exclusion_column = NULL,
                                   hek_string = NULL) {
             section <- theObject@args$peptideMissingValueImputation
             if (!is.list(section)) {
               section <- list()
             }
             resolve_section_value <- function(explicit_value, parameter_name, default_value = NULL) {
               if (!is.null(explicit_value)) {
                 return(explicit_value)
               }
               if (parameter_name %in% names(section) && !is.null(section[[parameter_name]])) {
                 return(section[[parameter_name]])
               }
               default_value
             }

             imputed_value_column <- resolve_section_value(
               imputed_value_column,
               "imputed_value_column",
               "Peptide.Imputed"
             )
             proportion_missing_values <- resolve_section_value(
               proportion_missing_values,
               "proportion_missing_values",
               0.50
             )
             core_utilisation <- resolve_section_value(core_utilisation, "core_utilisation", NA)
             exclusion_column <- resolve_section_value(exclusion_column, "exclusion_column", NULL)
             hek_string <- resolve_section_value(hek_string, "hek_string", NULL)
             if (length(exclusion_column) == 0L || is.na(exclusion_column[[1L]]) ||
                 !nzchar(trimws(as.character(exclusion_column[[1L]])))) {
               exclusion_column <- NULL
             } else {
               exclusion_column <- trimws(as.character(exclusion_column[[1L]]))
             }

             replicate_group_column <- theObject@technical_replicate_id
             has_declared_replicates <- length(replicate_group_column) == 1L &&
               !is.na(replicate_group_column) && nzchar(trimws(replicate_group_column)) &&
               replicate_group_column %in% names(theObject@design_matrix)

             if (!has_declared_replicates) {
               peptide_values_imputed <- theObject@peptide_data
               peptide_values_imputed[[imputed_value_column]] <-
                 peptide_values_imputed[[theObject@raw_quantity_column]]
               peptide_values_imputed$is_imputed <- rep(FALSE, nrow(peptide_values_imputed))
               filter_result <- list(
                 data = peptide_values_imputed,
                 support_table = data.frame(),
                 summary = list(
                   status = "skipped",
                   skip_reason = "technical_replicate_column_not_declared",
                   quantity_column = theObject@raw_quantity_column,
                   imputed_value_column = imputed_value_column,
                   eligibility_denominator = "distinct_design_runs_in_technical_replicate_group",
                   maximum_missing_fraction = as.numeric(proportion_missing_values),
                   eligibility_operator = "<=",
                   exclusion_source = if (is.null(exclusion_column)) "none" else paste0("design_column:", exclusion_column),
                   input_rows = nrow(theObject@peptide_data),
                   output_rows = nrow(peptide_values_imputed),
                   imputed_rows = 0L
                 )
               )
             } else {
               filter_result <- peptideMissingValueImputationHelper(
                 input_table = theObject@peptide_data,
                 metadata_table = theObject@design_matrix,
                 quantity_to_impute_column = !!rlang::sym(theObject@raw_quantity_column),
                 imputed_value_column = !!rlang::sym(imputed_value_column),
                 core_utilisation = core_utilisation,
                 input_table_sample_id_column = !!rlang::sym(theObject@sample_id),
                 sample_id_tbl_sample_id_column = !!rlang::sym(theObject@sample_id),
                 replicate_group_column = !!rlang::sym(replicate_group_column),
                 protein_id_column = !!rlang::sym(theObject@protein_id_column),
                 peptide_sequence_column = !!rlang::sym(theObject@peptide_sequence_column),
                 hek_string = hek_string,
                 proportion_missing_values = proportion_missing_values,
                 exclusion_column = exclusion_column,
                 return_imputation_result = TRUE
               )
               if (is.data.frame(filter_result)) {
                 filter_result <- list(
                   data = filter_result,
                   support_table = NULL,
                   summary = list(status = "applied", skip_reason = NULL)
                 )
               }
             }

             section$imputed_value_column <- imputed_value_column
             section$proportion_missing_values <- as.numeric(proportion_missing_values)
             section$core_utilisation <- core_utilisation
             section$exclusion_column <- exclusion_column
             section$hek_string <- hek_string
             section$quantity_column <- theObject@raw_quantity_column
             section$support_table <- filter_result$support_table
             section$imputation_summary <- filter_result$summary
             theObject@args$peptideMissingValueImputation <- section
             theObject@peptide_data <- filter_result$data
             theObject <- cleanDesignMatrixPeptide(theObject)
             theObject
           })

#' Impute Missing Peptide Values with limpa
#'
#' Uses limpa's detection-probability-curve approach to model the relationship
#' between peptide intensity and detection probability before imputing missing
#' values.
#'
#' @param theObject A `PeptideQuantitativeData` object.
#' @param imputed_value_column Name of the column that receives imputed values.
#' @param use_log2_transform Whether to transform values to log2 scale before
#'   imputation and restore the original scale afterwards.
#' @param verbose Whether to emit progress messages.
#' @param ensure_matrix Whether to calculate the peptide matrix when absent.
#'
#' @return The updated `PeptideQuantitativeData` object.
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

              slope_interpretation <- if (dpcfit$dpc[2] < 0.3) {
                "nearly random missing"
              } else if (dpcfit$dpc[2] < 0.7) {
                "moderate intensity-dependent missing"
              } else if (dpcfit$dpc[2] < 1.2) {
                "strong intensity-dependent missing"
              } else {
                "very strong intensity-dependent missing (approaching left-censoring)"
              }

              if (verbose) {
                log_info("DPC parameters estimated:")
                log_info("  beta0 (intercept): {round(dpcfit$dpc[1], 4)}")
                log_info("  beta1 (slope): {round(dpcfit$dpc[2], 4)}")
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

              # Convert with the authoritative feature map; matrix row labels
              # are never parsed as biological identifiers.
              imputed_identity <- .peptideMatrixToIdentityLong(
                theObject,
                imputed_matrix,
                imputed_value_column
              )
              theObject <- imputed_identity$theObject
              imputed_long <- imputed_identity$data

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
