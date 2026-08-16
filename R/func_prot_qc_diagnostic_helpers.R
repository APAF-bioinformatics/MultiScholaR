# ----------------------------------------------------------------------------
# checkProteinNAPercentages
# ----------------------------------------------------------------------------
#' Check Missing Value Percentages in Protein Data
#'
#' @description Calculate and report the percentage of missing values (NAs) in protein data
#' at different levels: total dataset, per sample, and per group.
#'
#' @param protein_obj A ProteinQuantitativeData S4 object
#' @param verbose Logical, whether to print detailed results (default: TRUE)
#'
#' @return A list containing:
#' \itemize{
#'   \item total_na_percent: Overall percentage of NAs in the dataset
#'   \item per_sample_na: Data frame with NA percentages per sample
#'   \item per_group_na: Data frame with NA percentages per group
#'   \item summary_stats: Summary statistics of NA distribution
#' }
#'
#' @export
checkProteinNAPercentages <- function(protein_obj, verbose = TRUE) {
  # Validate input
  if (!is(protein_obj, "ProteinQuantitativeData")) {
    stop("Input must be a ProteinQuantitativeData S4 object")
  }

  # Extract data from S4 object
  protein_quant_table <- protein_obj@protein_quant_table
  design_matrix <- protein_obj@design_matrix
  sample_id_col <- protein_obj@sample_id
  group_id_col <- protein_obj@group_id
  protein_id_col <- protein_obj@protein_id_column

  # Identify sample columns (exclude protein ID column)
  sample_columns <- setdiff(colnames(protein_quant_table), protein_id_col)

  # Validate that sample columns match design matrix
  if (length(sample_columns) != nrow(design_matrix)) {
    stop("Number of sample columns doesn't match design_matrix rows")
  }

  # Extract quantitative data matrix (samples only)
  protein_matrix <- as.matrix(protein_quant_table[, sample_columns])
  rownames(protein_matrix) <- protein_quant_table[[protein_id_col]]

  # Calculate total NA percentage
  total_values <- length(protein_matrix)
  total_nas <- sum(is.na(protein_matrix))
  total_na_percent <- (total_nas / total_values) * 100

  # Calculate per-sample NA percentages
  sample_na_counts <- apply(protein_matrix, 2, function(x) sum(is.na(x)))
  sample_na_percentages <- (sample_na_counts / nrow(protein_matrix)) * 100

  per_sample_na <- data.frame(
    sample = names(sample_na_counts),
    na_count = sample_na_counts,
    na_percentage = sample_na_percentages,
    stringsAsFactors = FALSE
  )

  # Add group information to per-sample results
  per_sample_na <- merge(per_sample_na, design_matrix,
    by.x = "sample", by.y = sample_id_col, all.x = TRUE
  )

  # Calculate per-group NA percentages
  per_group_na <- per_sample_na %>%
    group_by(!!sym(group_id_col)) %>%
    summarise(
      num_samples = n(),
      mean_na_percentage = mean(na_percentage, na.rm = TRUE),
      median_na_percentage = median(na_percentage, na.rm = TRUE),
      min_na_percentage = min(na_percentage, na.rm = TRUE),
      max_na_percentage = max(na_percentage, na.rm = TRUE),
      sd_na_percentage = sd(na_percentage, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(mean_na_percentage)

  # Calculate summary statistics
  summary_stats <- list(
    total_proteins = nrow(protein_matrix),
    total_samples = ncol(protein_matrix),
    total_groups = length(unique(design_matrix[[group_id_col]])),
    total_values = total_values,
    total_nas = total_nas,
    mean_na_per_sample = mean(sample_na_percentages),
    median_na_per_sample = median(sample_na_percentages),
    min_na_per_sample = min(sample_na_percentages),
    max_na_per_sample = max(sample_na_percentages)
  )

  # Print results if verbose
  if (verbose) {
    cat("\n=== Protein Data Missing Value Analysis ===\n")
    cat(sprintf(
      "Dataset dimensions: %d proteins x %d samples\n",
      nrow(protein_matrix), ncol(protein_matrix)
    ))
    cat(sprintf("Number of groups: %d\n", summary_stats$total_groups))
    cat(sprintf(
      "Total missing values: %s out of %s (%.2f%%)\n",
      format(total_nas, big.mark = ","),
      format(total_values, big.mark = ","),
      total_na_percent
    ))

    cat("\n--- Per-Sample Missing Value Summary ---\n")
    cat(sprintf("Mean NA%% per sample: %.2f%%\n", summary_stats$mean_na_per_sample))
    cat(sprintf("Median NA%% per sample: %.2f%%\n", summary_stats$median_na_per_sample))
    cat(sprintf(
      "Range: %.2f%% - %.2f%%\n",
      summary_stats$min_na_per_sample, summary_stats$max_na_per_sample
    ))

    cat("\n--- Per-Group Missing Value Summary ---\n")
    print(per_group_na)

    cat("\n--- Samples with Highest Missing Values ---\n")
    top_missing_samples <- per_sample_na %>%
      arrange(desc(na_percentage)) %>%
      head(min(5, nrow(per_sample_na)))
    print(top_missing_samples[, c("sample", group_id_col, "na_percentage")])

    cat("\n--- Samples with Lowest Missing Values ---\n")
    bottom_missing_samples <- per_sample_na %>%
      arrange(na_percentage) %>%
      head(min(5, nrow(per_sample_na)))
    print(bottom_missing_samples[, c("sample", group_id_col, "na_percentage")])
  }

  # Return results
  results <- list(
    total_na_percent = total_na_percent,
    per_sample_na = per_sample_na,
    per_group_na = per_group_na,
    summary_stats = summary_stats
  )

  return(invisible(results))
}

# ----------------------------------------------------------------------------
# getProteinNARecommendations
# ----------------------------------------------------------------------------
#' Get Recommendations for Handling Protein-Level Missing Values
#'
#' @description Provides specific recommendations for dealing with missing values
#' in protein data based on the percentage and distribution of NAs.
#'
#' @param protein_obj A ProteinQuantitativeData S4 object
#' @param include_code Logical, whether to include example R code (default: TRUE)
#'
#' @return Prints recommendations and invisibly returns a list of strategies
#'
#' @export
getProteinNARecommendations <- function(protein_obj, include_code = TRUE) {
  # Get NA analysis
  na_results <- checkProteinNAPercentages(protein_obj, verbose = FALSE)
  na_percent <- na_results$total_na_percent

  cat("\n=== PROTEIN NA HANDLING RECOMMENDATIONS ===\n")
  cat(sprintf(
    "Your data: %.1f%% NAs across %d proteins\n\n",
    na_percent, na_results$summary_stats$total_proteins
  ))

  if (na_percent < 15) {
    cat("[RECOMMENDATION] RECOMMENDATION: Complete Case Analysis\n")
    cat("* Your data has excellent protein coverage\n")
    cat("* Can proceed with standard analysis on proteins with complete data\n")
    if (include_code) {
      cat("\n[NOTE] Example code:\n")
      cat("complete_proteins <- protein_obj@protein_quant_table[complete.cases(protein_obj@protein_quant_table), ]\n")
    }
  } else if (na_percent >= 15 && na_percent < 40) {
    cat("[RECOMMENDATION] RECOMMENDATION: Consider Protein-Level Imputation\n")
    cat("* Moderate missing values - imputation could be beneficial\n")
    cat("* Options: KNN, minimum value, or mixed imputation strategies\n")
    cat("* Alternative: Filter to proteins detected in >=X samples per group\n")
    if (include_code) {
      cat("\n[NOTE] Example filtering code:\n")
      cat("# Keep proteins detected in >=50% of samples per group\n")
      cat("filtered_proteins <- filterProteinsByGroupDetection(protein_obj, min_detection_rate = 0.5)\n")
    }
  } else if (na_percent >= 40 && na_percent < 60) {
    cat("[RECOMMENDATION] RECOMMENDATION: Strict Filtering + Targeted Imputation\n")
    cat("* High missing values suggest challenging sample/detection conditions\n")
    cat("* Focus on well-detected proteins (present in majority of samples)\n")
    cat("* Consider group-wise detection requirements\n")
    if (include_code) {
      cat("\n[NOTE] Example approach:\n")
      cat("# Keep proteins detected in >=70% of samples in at least one group\n")
      cat("robust_proteins <- filterProteinsByGroupwise(protein_obj, min_group_detection = 0.7)\n")
    }
  } else {
    cat("[WARNING]  RECOMMENDATION: Review Data Quality\n")
    cat("* Very high missing values (>60%) suggest potential issues\n")
    cat("* Check: sample quality, peptide identification, rollup parameters\n")
    cat("* Consider more stringent protein identification criteria\n")
    cat("* May need to focus only on highly abundant/well-detected proteins\n")
  }

  cat("\n[REFERENCE] STRATEGIES SUMMARY:\n")
  cat("1. Complete Case: Use only proteins with no NAs\n")
  cat("2. Filtering: Remove proteins with >X% missing values\n")
  cat("3. Group-wise: Require detection in >=Y% samples per group\n")
  cat("4. Imputation: Fill NAs with estimated values (KNN, minimum, etc.)\n")
  cat("5. Hybrid: Combine filtering + imputation\n")

  cat("\n[TIP] TIP: Protein NAs != Data Quality Issues\n")
  cat("Missing proteins often reflect:\n")
  cat("* Low abundance proteins below detection limit\n")
  cat("* Sample-specific biology (some proteins not expressed)\n")
  cat("* Normal variation in complex proteomes\n\n")

  strategies <- list(
    na_percent = na_percent,
    primary_recommendation = if (na_percent < 15) {
      "complete_case"
    } else if (na_percent < 40) {
      "imputation_or_filtering"
    } else if (na_percent < 60) {
      "strict_filtering"
    } else {
      "data_quality_review"
    },
    alternative_strategies = c("complete_case", "group_wise_filtering", "imputation", "hybrid")
  )

  return(invisible(strategies))
}

# ----------------------------------------------------------------------------
# validatePostImputationProteinData
# ----------------------------------------------------------------------------
#' Validate Post-Imputation Protein Data
#'
#' @description A simple wrapper to validate protein data after imputation,
#' specifically checking if imputation was successful.
#'
#' @param protein_obj A ProteinQuantitativeData S4 object (post-imputation)
#' @param expected_na_percent Expected NA percentage (default: varies based on protein data)
#' @param tolerance Tolerance for expected percentage (default: 10%)
#'
#' @return Logical indicating if validation passed, with detailed output
#'
#' @export
validatePostImputationProteinData <- function(protein_obj, expected_na_percent = NULL, tolerance = 10) {
  cat("\n=== POST-IMPUTATION PROTEIN DATA VALIDATION ===\n")
  cat("Note: Protein-level NAs occur even after peptide imputation because:\n")
  cat("* Proteins need >=1 detected peptide to get a quantification\n")
  cat("* Some proteins detected only in subset of samples\n")
  cat("* This is normal proteomics data behavior!\n\n")

  # Run the full NA analysis
  na_results <- checkProteinNAPercentages(protein_obj, verbose = TRUE)

  # Set expected NA percentage if not provided (proteins often have some NAs)
  if (is.null(expected_na_percent)) {
    # For protein data, NAs are very common due to missing peptides/proteins
    # Typical ranges: 20-50% depending on sample complexity and detection method
    expected_na_percent <- 35 # Realistic expectation for protein data
    cat(sprintf("Note: Using default expected NA%% of %.1f%% for protein data\n", expected_na_percent))
    cat("(Protein-level NAs are normal due to incomplete protein detection across samples)\n")
  }

  # Check if validation passes
  actual_na_percent <- na_results$total_na_percent
  is_valid <- abs(actual_na_percent - expected_na_percent) <= tolerance

  cat("\n--- VALIDATION RESULT ---\n")
  cat(sprintf("Expected NA%%: %.2f%% (+/- %.2f%%)\n", expected_na_percent, tolerance))
  cat(sprintf("Actual NA%%: %.2f%%\n", actual_na_percent))

  if (is_valid) {
    cat("[OK] VALIDATION PASSED: Protein data NA levels are within expected range!\n")
  } else {
    cat("[FAIL] VALIDATION FAILED: Unexpected NA percentage detected!\n")
    if (actual_na_percent > expected_na_percent + tolerance) {
      cat("  -> Issue: More NAs than expected. Check for missing proteins/peptides.\n")
    } else {
      cat("  -> Issue: Fewer NAs than expected. Possible over-imputation.\n")
    }
  }

  # Additional warnings for common issues
  if (actual_na_percent > 50) {
    cat("[WARNING] WARNING: Very high NA percentage (>50%) suggests data quality issues!\n")
  }

  if (actual_na_percent < 10) {
    cat("[INFO] INFO: Very low NA percentage (<10%) - excellent protein coverage!\n")
  }

  # Educational information about protein NAs
  if (actual_na_percent > 20 && actual_na_percent < 50) {
    cat("[INFO] INFO: NA percentage is typical for protein-level data\n")
    cat("  -> This reflects biological reality: not all proteins detected in all samples\n")
    cat("  -> Consider: protein-level imputation OR complete-case analysis\n")
  }

  if (na_results$summary_stats$max_na_per_sample > actual_na_percent + 10) {
    cat("[WARNING] WARNING: Large variation in NA% between samples detected!\n")
    cat("  -> Some samples may have much lower protein coverage.\n")
  }

  # Check for problematic samples (>80% missing)
  high_missing_samples <- na_results$per_sample_na[na_results$per_sample_na$na_percentage > 80, ]
  if (nrow(high_missing_samples) > 0) {
    cat("[WARNING] WARNING: Samples with >80% missing proteins detected:\n")
    print(high_missing_samples[, c("sample", "na_percentage")])
  }

  cat("\n")
  return(invisible(list(
    is_valid = is_valid,
    actual_na_percent = actual_na_percent,
    expected_na_percent = expected_na_percent,
    full_results = na_results
  )))
}

# ----------------------------------------------------------------------------
# getSamplesCorrelationMatrix
# ----------------------------------------------------------------------------
#' getSamplesCorrelationMatrix
#' @description Calculate the Pearson's correlation score between sample
#' @param input_table Table with samples as columns and peptides as rows. Contains the log peptide intensity values.
#' @export
#' @param metadata_tbl,is_HEK_column,use,method Runtime inputs used by this function; see the usage section for accepted values.
getSamplesCorrelationMatrix <- function(
  input_table,
  metadata_tbl,
  is_HEK_column = is_HEK,
  use = "pairwise.complete.obs",
  method = "pearson"
) {
  without_hek_samples <- metadata_tbl |>
    dplyr::filter({{ is_HEK_column }} == FALSE) |>
    pull(Run)

  correlation_samples_to_use <- intersect(colnames(input_table), without_hek_samples) |> sort()

  correlation_between_samples <- cor(input_table[, correlation_samples_to_use], use = use, method = method)
  which(is.na(correlation_between_samples))
  correlation_between_samples[is.na(correlation_between_samples)] <- 0

  return(correlation_between_samples)
}

