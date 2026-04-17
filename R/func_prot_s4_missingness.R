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
#' @export
setMethod(
  f = "proteinMissingValueImputationLimpa",
  signature = "ProteinQuantitativeData",
  definition = function(theObject,
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
    sample_columns <- setdiff(
      colnames(protein_quant_table),
      c(
        protein_id_column, "description", "gene_name",
        "protein_name", "organism", "length"
      )
    )

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
      tryCatch(
        {
          dpcfit <- limpa::dpc(protein_matrix)
          dpc_params <- dpcfit
          if (verbose) {
            log_info("DPC parameters estimated:")
            log_info("  beta0 (intercept): {round(dpcfit$dpc[1], 4)}")
            log_info("  beta1 (slope): {round(dpcfit$dpc[2], 4)}")
          }
        },
        error = function(e) {
          if (verbose) {
            log_warn("DPC estimation failed, using default slope: {dpc_slope}")
          }
          dpc_params <- NULL
        }
      )
    }

    # Apply missing value imputation
    if (verbose) {
      log_info("Applying protein-level missing value imputation...")
    }

    tryCatch(
      {
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
          tidyr::pivot_longer(
            cols = -all_of(protein_id_column),
            names_to = sample_id_column,
            values_to = quantified_protein_column
          )

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
              dpc_params$dpc # Extract parameters from DPC object
            } else if (is.numeric(dpc_params)) {
              dpc_params # Already numeric parameters
            } else {
              c(NA, dpc_slope)
            }
          } else {
            c(NA, dpc_slope)
          },
          dpc_object_used = if (is.list(dpc_params) && !is.null(dpc_params$dpc)) dpc_params else NULL,
          quantified_protein_column = quantified_protein_column,
          missing_percentage_before = round(100 * mean(is.na(protein_matrix)), 1),
          missing_percentage_after = 0, # DPC imputation produces complete data
          imputation_method = "limpa_dpc_protein_imputation",
          total_proteins_imputed = nrow(imputed_matrix)
        )

        if (verbose) {
          log_info("limpa protein-level imputation completed successfully!")
          log_info("New imputed column: {quantified_protein_column}")
          log_info("DPC results stored in object@args$limpa_protein_imputation_results")
        }

        return(theObject)
      },
      error = function(e) {
        log_error(paste("Error during limpa protein imputation:", e$message))
        stop(paste("limpa protein imputation failed:", e$message))
      }
    )
  }
)

#' @export
setMethod(
  f = "preservePeptideNaValues",
  signature = c("PeptideQuantitativeData", "ProteinQuantitativeData"),
  definition = function(peptide_obj, protein_obj) {
    preservePeptideNaValuesHelper(peptide_obj, protein_obj)
  }
)

#' @export
preservePeptideNaValuesHelper <- function(peptide_obj, protein_obj) {
  sample_id_column <- peptide_obj@sample_id
  protein_id_column <- peptide_obj@protein_id_column

  check_peptide_value <- peptide_obj@peptide_data |>
    group_by(!!sym(sample_id_column), !!sym(protein_id_column)) |>
    summarise(
      Peptide.Normalised = sum(Peptide.Normalised, na.rm = TRUE),
      is_na = sum(is.na(Peptide.Normalised)),
      num_values = n()
    ) |>
    mutate(Peptide.Normalised = if_else(is_na == num_values, NA_real_, Peptide.Normalised)) |>
    ungroup() |>
    arrange(!!sym(sample_id_column)) |>
    pivot_wider(
      id_cols = !!sym(protein_id_column),
      names_from = !!sym(sample_id_column),
      values_from = Peptide.Normalised,
      values_fill = NA_real_
    )

  check_peptide_value_cln <- check_peptide_value[
    rownames(protein_obj@protein_quant_table),
    colnames(protein_obj@protein_quant_table)
  ]

  if (length(which(rownames(protein_obj@protein_quant_table) == rownames(check_peptide_value_cln))) != nrow(check_peptide_value_cln)) {
    stop("The rows in the protein object and the peptide object do not match")
  }

  if (length(which(colnames(protein_obj@protein_quant_table) == colnames(check_peptide_value_cln))) != ncol(check_peptide_value_cln)) {
    stop("The columns in the protein object and the peptide object do not match")
  }

  protein_obj@protein_quant_table[is.na(check_peptide_value_cln)] <- NA

  protein_obj
}

#' Protein Quantification using limpa Package DPC-Quant Method
#'
#' This function uses the limpa package's Detection Probability Curve (DPC) approach
#' combined with DPC-Quant for sophisticated protein-level quantification from peptide data.
#' This method leverages the DPC estimated at the peptide level to perform robust
#' protein quantification that accounts for missing value mechanisms.
#'
#' @param theObject A PeptideQuantitativeData object with peptide-level data
#' @param dpc_results DPC results from peptide-level analysis. If NULL, will use
#'   results stored in theObject@args$limpa_dpc_results. If those don't exist,
#'   will estimate DPC using default slope.
#' @param dpc_slope Default DPC slope to use if no DPC results available (default: 0.8)
#' @param quantified_protein_column Name for the new column containing quantified protein values.
#'   Default is "Protein.Quantified.Limpa"
#' @param verbose Whether to print progress messages. Default is TRUE
#' @param chunk When verbose=TRUE, how often to output progress information (default: 1000)
#'
#' @details
#' The DPC-Quant method represents missing values probabilistically using the Detection
#' Probability Curve and returns maximum posterior estimates for all protein log2-expression
#' values. Unlike simple imputation, this method provides:
#'
#' 1. Protein-level quantification from peptide data using hierarchical modeling
#' 2. Uncertainty quantification via standard errors for each protein estimate
#' 3. No missing values in final protein summaries
#' 4. Proper handling of peptides with sparse observations
#'
#' The process follows these steps:
#' 1. Extract or estimate Detection Probability Curve from peptide data
#' 2. Apply dpcQuant() to obtain protein-level summaries
#' 3. Extract quantified protein matrix and uncertainty estimates
#' 4. Create ProteinQuantitativeData object with results
#'
#' @return ProteinQuantitativeData object with quantified protein values and metadata
#'
#' @export
setMethod(
  f = "proteinMissingValueImputationLimpa",
  signature = "PeptideQuantitativeData",
  definition = function(theObject,
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
      "Protein.Quantified.Limpa"
    } else {
      quantified_protein_column
    }

    dpc_slope <- checkParamsObjectFunctionSimplify(theObject, "dpc_slope", dpc_slope)
    verbose <- checkParamsObjectFunctionSimplify(theObject, "verbose", verbose)
    chunk <- checkParamsObjectFunctionSimplify(theObject, "chunk", chunk)

    # Extract data from peptide object
    peptide_data <- theObject@peptide_data
    peptide_matrix <- theObject@peptide_matrix
    protein_id_column <- theObject@protein_id_column
    peptide_sequence_column <- theObject@peptide_sequence_column
    sample_id_column <- theObject@sample_id
    design_matrix <- theObject@design_matrix

    if (verbose) {
      log_info("Starting limpa DPC-Quant protein quantification...")
      log_info("Peptide data dimensions: {nrow(peptide_matrix)} peptides x {ncol(peptide_matrix)} samples")
      log_info("Missing value percentage in peptides: {round(100 * mean(is.na(peptide_matrix)), 1)}%")
    }

    # Check if peptide matrix exists and is properly formatted
    if (is.null(peptide_matrix) || length(peptide_matrix) == 0) {
      stop("Peptide matrix not found. Please run calcPeptideMatrix() first.")
    }

    # Prepare peptide data for limpa (peptides as rows, samples as columns)
    y_peptide <- peptide_matrix

    # Ensure data is log2 transformed for limpa
    if (!theObject@is_logged_data) {
      if (verbose) {
        log_info("Converting to log2 scale for limpa DPC-Quant...")
      }
      y_peptide <- log2(y_peptide + 1)
    } else {
      if (verbose) {
        log_info("Using existing log2 transformed peptide data")
      }
    }

    # Handle infinite or NaN values
    if (any(is.infinite(y_peptide) | is.nan(y_peptide), na.rm = TRUE)) {
      if (verbose) {
        log_warn("Infinite or NaN values detected. Replacing with NA...")
      }
      y_peptide[is.infinite(y_peptide) | is.nan(y_peptide)] <- NA
    }

    # NEW: Remove rows that are all NA, which cause limpa optimization to fail
    all_na_rows <- rowSums(is.na(y_peptide)) == ncol(y_peptide)
    if (any(all_na_rows)) {
      num_removed <- sum(all_na_rows)
      if (verbose) {
        log_warn("Found and removed {num_removed} peptides with no observations across all samples.")
      }
      y_peptide <- y_peptide[!all_na_rows, ]
    }

    # Get or estimate DPC parameters
    dpc_params <- NULL
    slope_interpretation <- "nearly random missing" # Default interpretation
    if (!is.null(dpc_results)) {
      # Use provided DPC results
      dpc_params <- dpc_results
      if (verbose) {
        log_info("Using provided DPC results")
        if (is.numeric(dpc_results) && length(dpc_results) == 2) {
          log_info("  beta0 (intercept): {round(dpc_results[1], 4)}")
          log_info("  beta1 (slope): {round(dpc_results[2], 4)}")
        }
      }
    } else if (!is.null(theObject@args$limpa_dpc_results)) {
      # Prefer full DPC object if available, otherwise use parameters
      if (!is.null(theObject@args$limpa_dpc_results$dpc_object)) {
        dpc_params <- theObject@args$limpa_dpc_results$dpc_object
        # Extract slope interpretation if available
        if (!is.null(theObject@args$limpa_dpc_results$slope_interpretation)) {
          slope_interpretation <- theObject@args$limpa_dpc_results$slope_interpretation
        }
        if (verbose) {
          log_info("Using full DPC object from previous peptide analysis:")
          log_info("  beta0 (intercept): {round(dpc_params$dpc[1], 4)}")
          log_info("  beta1 (slope): {round(dpc_params$dpc[2], 4)}")
          log_info("  Interpretation: {slope_interpretation}")
        }
      } else if (!is.null(theObject@args$limpa_dpc_results$dpc_parameters)) {
        dpc_params <- theObject@args$limpa_dpc_results$dpc_parameters
        # Extract slope interpretation if available
        if (!is.null(theObject@args$limpa_dpc_results$slope_interpretation)) {
          slope_interpretation <- theObject@args$limpa_dpc_results$slope_interpretation
        }
        if (verbose) {
          log_info("Using DPC parameters from previous peptide analysis:")
          log_info("  beta0 (intercept): {round(dpc_params[1], 4)}")
          log_info("  beta1 (slope): {round(dpc_params[2], 4)}")
          log_info("  Interpretation: {slope_interpretation}")
        }
      }
    } else {
      # Estimate DPC parameters on-the-fly if not available
      if (verbose) {
        log_info("No DPC results found. Estimating DPC parameters from peptide data...")
      }
      tryCatch(
        {
          # Estimate DPC using the peptide data
          dpcfit <- limpa::dpc(y_peptide)
          dpc_params <- dpcfit # Use the full DPC object

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

          if (verbose) {
            log_info("DPC parameters estimated on-the-fly:")
            log_info("  beta0 (intercept): {round(dpcfit$dpc[1], 4)}")
            log_info("  beta1 (slope): {round(dpcfit$dpc[2], 4)}")
            log_info("  Interpretation: {slope_interpretation}")
          }
        },
        error = function(e) {
          if (verbose) {
            log_warn("DPC estimation failed: {e$message}. Using default slope: {dpc_slope}")
          }
          dpc_params <<- NULL # Will let dpcQuant estimate with default slope
          slope_interpretation <<- "unable to determine (DPC estimation failed)"
        }
      )
    }

    # Create protein ID mapping from peptide data
    # Extract protein IDs from rownames (format: "ProteinID%PeptideSequence")
    rownames_split <- strsplit(rownames(y_peptide), "%")
    protein_ids <- sapply(rownames_split, function(x) x[1])

    # Calculate peptide and peptidoform counts per protein for future filtering
    # This ensures we carry this critical info to the protein level
    peptide_summary <- peptide_data |>
      dplyr::group_by(!!sym(protein_id_column)) |>
      dplyr::summarise(
        peptide_count = dplyr::n_distinct(!!sym(peptide_sequence_column)),
        peptidoform_count = dplyr::n(), # Total peptidoforms is just the number of rows per protein
        .groups = "drop"
      )

    if (verbose) {
      unique_proteins <- nrow(peptide_summary)
      log_info("Found {unique_proteins} unique proteins from {nrow(y_peptide)} peptidoforms")

      peptide_count_summary <- peptide_summary |>
        dplyr::summarise(
          min_peptides = min(peptide_count),
          max_peptides = max(peptide_count),
          median_peptides = median(peptide_count),
          min_peptidoforms = min(peptidoform_count),
          max_peptidoforms = max(peptidoform_count),
          median_peptidoforms = median(peptidoform_count),
          proteins_with_2plus_peptides = sum(peptide_count >= 2)
        )
      log_info("Unique peptide counts per protein - min: {peptide_count_summary$min_peptides}, max: {peptide_count_summary$max_peptides}, median: {round(peptide_count_summary$median_peptides, 1)}")
      log_info("Total peptidoform counts per protein - min: {peptide_count_summary$min_peptidoforms}, max: {peptide_count_summary$max_peptidoforms}, median: {round(peptide_count_summary$median_peptidoforms, 1)}")
      log_info("Proteins with >=2 unique peptides: {peptide_count_summary$proteins_with_2plus_peptides} / {unique_proteins}")
    }

    # Apply DPC-Quant for protein quantification
    if (verbose) {
      log_info("Applying DPC-Quant for protein quantification...")
    }

    tryCatch(
      {
        # Create EList-like object for limpa
        y_elist <- list(
          E = y_peptide,
          genes = data.frame(
            protein.id = protein_ids,
            peptide.id = rownames(y_peptide),
            stringsAsFactors = FALSE
          )
        )
        class(y_elist) <- "EList"

        # Apply dpcQuant
        if (!is.null(dpc_params)) {
          protein_quant_result <- limpa::dpcQuant(
            y = y_elist,
            protein.id = "protein.id",
            dpc = dpc_params,
            verbose = verbose,
            chunk = chunk
          )
        } else {
          protein_quant_result <- limpa::dpcQuant(
            y = y_elist,
            protein.id = "protein.id",
            dpc.slope = dpc_slope,
            verbose = verbose,
            chunk = chunk
          )
        }

        if (verbose) {
          log_info("DPC-Quant completed successfully")
          log_info("Quantified proteins: {nrow(protein_quant_result$E)}")
          log_info("Samples: {ncol(protein_quant_result$E)}")
          log_info("No missing values in protein quantification: {!any(is.na(protein_quant_result$E))}")
        }

        # Extract quantified protein matrix
        protein_matrix <- protein_quant_result$E
        rownames(protein_matrix) <- protein_quant_result$genes$protein.id
        colnames(protein_matrix) <- colnames(y_peptide)

        # Extract standard errors and observation counts
        standard_errors <- protein_quant_result$other$standard.error
        n_observations <- protein_quant_result$other$n.observations

        # Transform back to original scale if necessary
        if (!theObject@is_logged_data) {
          if (verbose) {
            log_info("Converting back from log2 scale...")
          }
          protein_matrix <- 2^protein_matrix - 1
          # Ensure no negative values
          protein_matrix[protein_matrix < 0] <- 0
        }

        # Convert protein matrix to wide format data frame (samples as columns)
        if (verbose) {
          log_info("Converting protein data to wide format...")
        }

        protein_wide <- protein_matrix |>
          as.data.frame() |>
          tibble::rownames_to_column(protein_id_column)

        # --- NEW: Synchronize the peptide summary with the actual quantified proteins ---
        # dpcQuant may not be able to quantify all proteins present in the input.
        # We must ensure the summary table only contains proteins that are in the final matrix.
        final_protein_ids <- protein_wide[[protein_id_column]]
        peptide_summary_synced <- peptide_summary |>
          dplyr::filter(!!sym(protein_id_column) %in% final_protein_ids)

        if (verbose) {
          log_info("Synchronized peptide summary table with quantified proteins. Kept {nrow(peptide_summary_synced)} of {nrow(peptide_summary)} entries.")
        }
        # --- END NEW ---

        # Create ProteinQuantitativeData object
        if (verbose) {
          log_info("Creating ProteinQuantitativeData object...")
        }

        protein_obj <- ProteinQuantitativeData(
          protein_quant_table = protein_wide,
          protein_id_column = protein_id_column,
          protein_id_table = protein_wide |> dplyr::distinct(!!sym(protein_id_column)),
          design_matrix = design_matrix,
          sample_id = sample_id_column,
          group_id = theObject@group_id,
          technical_replicate_id = theObject@technical_replicate_id,
          args = theObject@args
        )

        # Store DPC-Quant results in args for future reference
        if (is.null(protein_obj@args)) {
          protein_obj@args <- list()
        }

        protein_obj@args$limpa_dpc_quant_results <- list(
          dpc_parameters_used = if (!is.null(dpc_params)) {
            if (is.list(dpc_params) && !is.null(dpc_params$dpc)) {
              dpc_params$dpc # Extract parameters from DPC object
            } else if (is.numeric(dpc_params)) {
              dpc_params # Already numeric parameters
            } else {
              c(NA, dpc_slope)
            }
          } else {
            c(NA, dpc_slope)
          },
          dpc_object_used = if (is.list(dpc_params) && !is.null(dpc_params$dpc)) dpc_params else NULL,
          slope_interpretation = slope_interpretation, # Store slope interpretation for QC plots
          y_peptide_for_dpc = y_peptide, # Store peptide data for DPC plotting
          quantified_elist = protein_quant_result, # Store the entire EList object
          standard_errors = standard_errors,
          n_observations = n_observations,
          peptide_counts_per_protein = peptide_summary_synced, # Use the new SYNCED summary table
          missing_percentage_peptides = round(100 * mean(is.na(y_peptide)), 1),
          missing_percentage_proteins = 0, # DPC-Quant produces complete protein data
          quantification_method = "limpa_dpc_quant",
          total_proteins_quantified = nrow(protein_matrix),
          total_peptides_used = nrow(y_peptide),
          dpc_method = "limpa_dpc_quant" # Identify the method type
        )

        if (verbose) {
          log_info("limpa DPC-Quant protein quantification completed successfully!")
          log_info("Protein data stored in wide format with {ncol(protein_matrix)} sample columns")
          log_info("Proteins quantified: {nrow(protein_matrix)}")
          log_info("DPC-Quant results stored in object@args$limpa_dpc_quant_results")
        }

        return(protein_obj)
      },
      error = function(e) {
        log_error(paste("Error during limpa DPC-Quant:", e$message))
        stop(paste("limpa DPC-Quant failed:", e$message))
      }
    )
  }
)
