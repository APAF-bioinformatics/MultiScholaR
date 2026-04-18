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

# ----------------------------------------------------------------------------
# updateProteinFiltering
# ----------------------------------------------------------------------------
#' @title Update and Visualize Filtering Progress
#' @description Tracks and visualizes the impact of filtering steps on peptide
#'   and protein counts. Updates a global `FilteringProgress` object and optionally
#'   saves plots summarizing the changes. Handles both peptide-level and
#'   protein-level data inputs.
#'
#' @details
#' This function acts as a central hub for monitoring data reduction throughout
#' a filtering workflow. It performs the following actions:
#' \itemize{
#'   \item Initializes or retrieves a global S4 object named `filtering_progress`
#'     of class `FilteringProgress`.
#'   \item Calculates key metrics (total unique proteins, proteins per run,
#'     total unique peptides, peptides per protein distribution, peptides per run)
#'     based on the input `data`. Peptide metrics are only calculated or updated
#'     if `data` is identified as peptide-level data. For protein-level data,
#'     peptide metrics from the last peptide step (if any) are carried forward or
#'     initialized as empty/NA.
#'   \item Adds or updates these metrics in the `filtering_progress` object
#'     under the specified `step_name`.
#'   \item Generates summary plots using `ggplot2`:
#'     \itemize{
#'       \item Bar plot of total unique proteins per step.
#'       \item Bar plot of total unique peptides per step (or placeholder if only protein data).
#'       \item Box plot of peptides per protein distribution per step (or placeholder).
#'       \item Line plot of proteins per run across steps.
#'       \item Line plot of peptides per run across steps (or placeholder).
#'     }
#'   \item If `omic_type` and `experiment_label` are provided and valid paths can be
#'     derived from the global `project_dirs` object, the generated plots are saved
#'     as PNG files into the derived `time_dir`. Warnings are issued if paths cannot be
#'     derived or if `project_dirs` is not found.
#'   \item If `return_grid` is `TRUE`, arranges the plots into a single grid using
#'     `gridExtra` and returns the grid object (grob). Also saves this combined grid
#'     if plot saving is enabled.
#'   \item If `return_grid` is `FALSE` (default), prints each plot individually
#'     and returns the list of plot objects invisibly.
#' }
#'
#' **Important:** This function relies on and modifies a global variable named
#' `filtering_progress`. For saving plots, it depends on the global `project_dirs`
#' object (expected to be populated by `setupDirectories()`) and the successful
#' derivation of `time_dir` from it using `omic_type` and `experiment_label`.
#'
#' @param data The input data object. Can be a data frame (expected to conform
#'   to typical peptide or protein quantification structures) or an S4 object
#'   containing relevant slots (e.g., inheriting from `SummarizedExperiment`).
#'   The function attempts to automatically detect if it\'s peptide or protein data.
#' @param step_name A character string uniquely identifying the current filtering
#'   step (e.g., "InitialData", "FilteredByQuality", "Normalized"). This name is
#'   used for tracking in the `filtering_progress` object and plot labels.
#' @param omic_type Optional character string. The type of omics data
#'   (e.g., "proteomics", "metabolomics"). Used with `experiment_label` to
#'   derive save paths from the global `project_dirs` object. If `NULL` (default)
#'   or `experiment_label` is `NULL`, plots are not saved.
#' @param experiment_label Optional character string. The specific experiment
#'   label (e.g., "workshop_data"). Used with `omic_type` to derive save paths
#'   from the global `project_dirs` object. If `NULL` (default) or `omic_type`
#'   is `NULL`, plots are not saved.
#' @param overwrite Logical. If `TRUE`, allows overwriting an existing entry for
#'   `step_name` in the `filtering_progress` object. If `FALSE` (default) and
#'   `step_name` already exists, the function will stop with an error.
#' @param return_grid Logical. If `TRUE`, returns a single combined plot grid
#'   object created with `gridExtra::grid.arrange()`. If `FALSE` (default), prints
#'   individual plots and returns an invisible list of the ggplot objects.
#'
#' @return If `return_grid` is `TRUE`, returns a `grob` object (a grid graphical object).
#'   If `return_grid` is `FALSE`, returns an invisible list containing the individual
#'   `ggplot` objects (`proteins_total`, `proteins_per_run`, `peptides_total`,
#'   `peptides_per_protein`, `peptides_per_run`). Has side effects: modifies the
#'   global `filtering_progress` object and potentially saves plots to disk if
#'   `omic_type` and `experiment_label` are provided and paths are valid.
#'
#' @importFrom ggplot2 ggplot aes geom_bar geom_text labs theme_minimal theme element_text element_blank geom_line geom_point scale_color_manual annotate theme_void geom_boxplot coord_cartesian ggsave
#' @importFrom dplyr bind_rows mutate group_by ungroup %>%
#' @importFrom forcats fct_reorder
#' @importFrom gridExtra arrangeGrob grid.arrange
#' @importFrom methods slotNames new is
#' @importFrom stats quantile
#' @importFrom multidplyr partition new_cluster cluster_library
#' @importFrom future plan
#'
#' @export
updateProteinFiltering <- function(data, step_name,
                                   omic_type = NULL, experiment_label = NULL,
                                   overwrite = FALSE, return_grid = FALSE,
                                   formats = c("png", "pdf"),
                                   project_dirs = NULL) {
  # Initialize filtering_progress if it doesn\'t exist
  if (!exists("filtering_progress", envir = .GlobalEnv)) {
    filtering_progress <- new("FilteringProgress")
    assign("filtering_progress", filtering_progress, envir = .GlobalEnv)
  }

  # Get the current filtering_progress object
  filtering_progress <- get("filtering_progress", envir = .GlobalEnv)

  # DEBUG66: Memory check
  message("--- DEBUG66 [updatePeptideFiltering]: Entry ---")
  message(sprintf("   [updatePeptideFiltering] Step Name: %s", step_name))
  message(sprintf("   [updatePeptideFiltering] filtering_progress size: %s", format(object.size(filtering_progress), units = "auto")))
  gc()

  # Path derivation and save_plots logic
  derived_time_dir <- NULL
  save_plots <- FALSE

  if (!is.null(omic_type) && !is.null(experiment_label)) {
    # Check if project_dirs is provided as argument or exists in Global Env
    project_dirs_to_use <- if (!is.null(project_dirs)) {
      project_dirs
    } else if (exists("project_dirs", envir = .GlobalEnv)) {
      get("project_dirs", envir = .GlobalEnv)
    } else {
      NULL
    }

    if (is.null(project_dirs_to_use)) {
      warning("Object 'project_dirs' not found (neither as argument nor in Global Env). Plots will not be saved. Ensure 'setupDirectories()' has been run.")
    } else {
      project_dirs_global <- project_dirs_to_use
      omic_project_key <- paste0(omic_type, "_", experiment_label)

      if (!omic_project_key %in% names(project_dirs_global)) {
        warning(paste0("Entry for \'", omic_project_key, "\' not found in global \'project_dirs\'. Plots will not be saved."))
      } else {
        current_project_paths <- project_dirs_global[[omic_project_key]]
        if (is.null(current_project_paths)) {
          warning(paste0("Entry for \'", omic_project_key, "\' in global \'project_dirs\' is NULL. Plots will not be saved."))
        } else {
          derived_publication_graphs_dir <- current_project_paths$publication_graphs_dir
          temp_time_dir <- current_project_paths$time_dir

          if (is.null(temp_time_dir) || !is.character(temp_time_dir) || length(temp_time_dir) != 1 ||
            is.null(derived_publication_graphs_dir) || !is.character(derived_publication_graphs_dir) || length(derived_publication_graphs_dir) != 1) {
            warning(paste0(
              "\'time_dir\' or \'publication_graphs_dir\' is missing, not a character string, or not a single path for \'", omic_project_key,
              "\' in global \'project_dirs\'. Plots will not be saved."
            ))
          } else {
            if (!dir.exists(temp_time_dir)) {
              warning(paste0(
                "The derived \'time_dir\' (", temp_time_dir, ") for \'", omic_project_key,
                "\' does not exist. Plots will not be saved. Ensure directories are created via setupDirectories()."
              ))
            } else {
              derived_time_dir <- temp_time_dir
              save_plots <- TRUE
              message(paste0("Plots will be saved to: ", derived_time_dir))
            }
          }
        }
      }
    }
  } else {
    # Message if omic_type/label are missing and saving might have been expected
    if (return_grid && (is.null(omic_type) || is.null(experiment_label))) {
      message("omic_type and/or experiment_label not provided. Plots will not be saved.")
    }
  }

  # Determine if we\'re working with protein_quant_table
  is_protein_quant <- if (methods::is(data, "S4")) {
    "protein_quant_table" %in% slotNames(data)
  } else {
    # For data frames, check if it looks like a protein quant table
    if ("Protein.Ids" %in% names(data)) {
      all(sapply(data[setdiff(names(data), "Protein.Ids")], is.numeric))
    } else {
      FALSE
    }
  }

  # Calculate protein metrics (always done)
  protein_count <- countUniqueProteins(data)
  proteins_per_run <- countProteinsPerRun(data)

  # Ensure consistent data types in proteins_per_run
  proteins_per_run$Run <- as.character(proteins_per_run$Run)
  proteins_per_run$n_proteins <- as.numeric(proteins_per_run$n_proteins)

  # Update filtering progress based on data type
  if (step_name %in% filtering_progress@steps) {
    if (!overwrite) {
      stop("Step name \'", step_name, "\' already exists. Use overwrite = TRUE to replace it.")
    }
    idx <- which(filtering_progress@steps == step_name)

    # Always update protein metrics
    filtering_progress@proteins[idx] <- protein_count
    filtering_progress@proteins_per_run[[idx]] <- proteins_per_run

    if (!is_protein_quant) {
      # Update peptide metrics only for peptide data
      filtering_progress@total_peptides[idx] <- calcTotalPeptides(data)
      peptides_per_protein <- calcPeptidesPerProtein(data)
      peptides_per_run <- countPeptidesPerRun(data)

      # Ensure consistent data types
      peptides_per_protein$Protein.Ids <- as.character(peptides_per_protein$Protein.Ids)
      peptides_per_protein$n_peptides <- as.numeric(peptides_per_protein$n_peptides)

      peptides_per_run$Run <- as.character(peptides_per_run$Run)
      peptides_per_run$n_peptides <- as.numeric(peptides_per_run$n_peptides)

      filtering_progress@peptides_per_protein[[idx]] <- peptides_per_protein
      filtering_progress@peptides_per_run[[idx]] <- peptides_per_run
    }
  } else {
    filtering_progress@steps <- c(filtering_progress@steps, step_name)
    filtering_progress@proteins <- c(filtering_progress@proteins, protein_count)
    filtering_progress@proteins_per_run <- c(
      filtering_progress@proteins_per_run,
      list(proteins_per_run)
    )

    if (!is_protein_quant) {
      # Add peptide metrics only for peptide data
      filtering_progress@total_peptides <- c(
        filtering_progress@total_peptides,
        calcTotalPeptides(data)
      )

      peptides_per_protein <- calcPeptidesPerProtein(data)
      peptides_per_run <- countPeptidesPerRun(data)

      # Ensure consistent data types
      peptides_per_protein$Protein.Ids <- as.character(peptides_per_protein$Protein.Ids)
      peptides_per_protein$n_peptides <- as.numeric(peptides_per_protein$n_peptides)

      peptides_per_run$Run <- as.character(peptides_per_run$Run)
      peptides_per_run$n_peptides <- as.numeric(peptides_per_run$n_peptides)

      filtering_progress@peptides_per_protein <- c(
        filtering_progress@peptides_per_protein,
        list(peptides_per_protein)
      )
      filtering_progress@peptides_per_run <- c(
        filtering_progress@peptides_per_run,
        list(peptides_per_run)
      )
    } else {
      # For protein data, maintain existing peptide metrics or add NA/empty entries
      if (length(filtering_progress@total_peptides) > 0) {
        filtering_progress@total_peptides <- c(
          filtering_progress@total_peptides,
          filtering_progress@total_peptides[length(filtering_progress@total_peptides)]
        )
        filtering_progress@peptides_per_protein <- c(
          filtering_progress@peptides_per_protein,
          filtering_progress@peptides_per_protein[length(filtering_progress@peptides_per_protein)]
        )
        filtering_progress@peptides_per_run <- c(
          filtering_progress@peptides_per_run,
          filtering_progress@peptides_per_run[length(filtering_progress@peptides_per_run)]
        )
      } else {
        filtering_progress@total_peptides <- c(filtering_progress@total_peptides, NA_integer_)
        filtering_progress@peptides_per_protein <- c(
          filtering_progress@peptides_per_protein,
          list(data.frame(
            Protein.Ids = character(),
            n_peptides = integer()
          ))
        )
        filtering_progress@peptides_per_run <- c(
          filtering_progress@peptides_per_run,
          list(data.frame(
            Run = character(),
            n_peptides = integer()
          ))
        )
      }
    }
  }

  # Update the global filtering_progress object
  assign("filtering_progress", filtering_progress, envir = .GlobalEnv)

  # Create base protein count plot (always shown)
  message("   [updateProteinFiltering] Generating P1 (Protein Count Bar)...")
  p1 <- ggplot(data.frame(
    step = factor(filtering_progress@steps, levels = filtering_progress@steps),
    proteins = filtering_progress@proteins
  ), aes(x = step, y = proteins)) +
    geom_bar(stat = "identity", fill = "steelblue", width = 0.7) +
    geom_text(aes(label = proteins),
      vjust = -0.5,
      size = 4
    ) +
    labs(
      title = "Number of Proteins",
      x = "Filtering Step",
      y = "Unique Proteins"
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      panel.grid.major.x = element_blank()
    )

  # Create proteins per run plot (always shown)
  message("   [updateProteinFiltering] Generating P4 (Proteins Per Run Line)...")
  # First ensure all data frames in the list have consistent column types
  proteins_per_run_list <- lapply(filtering_progress@proteins_per_run, function(df) {
    df$Run <- as.character(df$Run)
    df$n_proteins <- as.numeric(df$n_proteins)
    return(df)
  })

  message(sprintf("      [P4] Binding %d data frames...", length(proteins_per_run_list)))
  p4_data <- bind_rows(proteins_per_run_list, .id = "step")
  message(sprintf("      [P4] Combined data rows: %d", nrow(p4_data)))

  p4 <- p4_data |>
    mutate(step = filtering_progress@steps[as.numeric(step)]) |>
    group_by(Run) |>
    mutate(avg_proteins = mean(n_proteins)) |>
    ungroup() |>
    # Run is already character from our preprocessing
    mutate(Run = fct_reorder(Run, avg_proteins)) |>
    ggplot(aes(
      x = Run, y = n_proteins,
      group = step,
      color = factor(step, levels = filtering_progress@steps)
    )) +
    geom_line() +
    geom_point() +
    labs(
      title = "Proteins per Run",
      x = "Run ID (ordered by average protein count)",
      y = "Number of Proteins",
      color = "Step"
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      panel.grid.major.x = element_blank()
    ) +
    scale_color_manual(values = get_color_palette(length(filtering_progress@steps), "steelblue"))

  message("   [updateProteinFiltering] P4 Generated.")

  # Initialize peptide plots
  if (is_protein_quant) {
    # For protein data, create empty placeholder plots if no peptide data exists
    if (all(is.na(filtering_progress@total_peptides))) {
      p2 <- p3 <- p5 <- ggplot() +
        annotate("text",
          x = 0.5, y = 0.5,
          label = "No peptide data available for protein quantification data"
        ) +
        theme_void()
    } else {
      # If peptide data exists from previous steps, create plots with existing data
      p2 <- ggplot(data.frame(
        step = factor(filtering_progress@steps, levels = filtering_progress@steps),
        total_peptides = filtering_progress@total_peptides
      ), aes(x = step, y = total_peptides)) +
        geom_bar(stat = "identity", fill = "forestgreen", width = 0.7) +
        geom_text(aes(label = total_peptides),
          vjust = -0.5,
          size = 4
        ) +
        labs(
          title = "Total Unique Peptides (from last peptide data)",
          x = "Filtering Step",
          y = "Unique Peptides"
        ) +
        theme_minimal() +
        theme(
          axis.text.x = element_text(angle = 45, hjust = 1),
          panel.grid.major.x = element_blank()
        )

      # Ensure consistent data types in peptides_per_protein list
      peptides_per_protein_list <- lapply(filtering_progress@peptides_per_protein, function(df) {
        if (nrow(df) > 0) {
          df$Protein.Ids <- as.character(df$Protein.Ids)
          df$n_peptides <- as.numeric(df$n_peptides)
        }
        return(df)
      })

      p3 <- ggplot() +
        geom_boxplot(
          data = bind_rows(peptides_per_protein_list, .id = "step") |>
            mutate(step = filtering_progress@steps[as.numeric(step)]),
          aes(
            x = factor(step, levels = filtering_progress@steps),
            y = n_peptides
          ),
          fill = "darkred",
          alpha = 0.5,
          outlier.shape = NA
        ) +
        labs(
          title = "Peptides per Protein Distribution (from last peptide data)",
          x = "Filtering Step",
          y = "Number of Peptides"
        ) +
        theme_minimal() +
        theme(
          axis.text.x = element_text(angle = 45, hjust = 1),
          panel.grid.major.x = element_blank()
        ) +
        coord_cartesian(
          ylim = c(
            0,
            quantile(bind_rows(peptides_per_protein_list)$n_peptides, 0.95)
          )
        )

      # Ensure consistent data types in peptides_per_run list
      peptides_per_run_list <- lapply(filtering_progress@peptides_per_run, function(df) {
        if (nrow(df) > 0) {
          df$Run <- as.character(df$Run)
          df$n_peptides <- as.numeric(df$n_peptides)
        }
        return(df)
      })

      p5 <- bind_rows(peptides_per_run_list, .id = "step") |>
        mutate(step = filtering_progress@steps[as.numeric(step)]) |>
        group_by(Run) |>
        mutate(avg_peptides = mean(n_peptides)) |>
        ungroup() |>
        # Run is already character from our preprocessing
        mutate(Run = fct_reorder(Run, avg_peptides)) |>
        ggplot(aes(
          x = Run, y = n_peptides,
          group = step,
          color = factor(step, levels = filtering_progress@steps)
        )) +
        geom_line() +
        geom_point() +
        labs(
          title = "Peptides per Run (from last peptide data)",
          x = "Run ID (ordered by average peptide count)",
          y = "Number of Peptides",
          color = "Step"
        ) +
        theme_minimal() +
        theme(
          axis.text.x = element_text(angle = 45, hjust = 1),
          panel.grid.major.x = element_blank()
        ) +
        scale_color_manual(values = get_color_palette(length(filtering_progress@steps), "forestgreen"))
    }
  } else {
    # For peptide data, create normal plots
    p2 <- ggplot(data.frame(
      step = factor(filtering_progress@steps, levels = filtering_progress@steps),
      total_peptides = filtering_progress@total_peptides
    ), aes(x = step, y = total_peptides)) +
      geom_bar(stat = "identity", fill = "forestgreen", width = 0.7) +
      geom_text(aes(label = total_peptides),
        vjust = -0.5,
        size = 4
      ) +
      labs(
        title = "Total Unique Peptides",
        x = "Filtering Step",
        y = "Unique Peptides"
      ) +
      theme_minimal() +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid.major.x = element_blank()
      )

    # Ensure consistent data types in peptides_per_protein list
    peptides_per_protein_list <- lapply(filtering_progress@peptides_per_protein, function(df) {
      if (nrow(df) > 0) {
        df$Protein.Ids <- as.character(df$Protein.Ids)
        df$n_peptides <- as.numeric(df$n_peptides)
      }
      return(df)
    })

    p3 <- ggplot() +
      geom_boxplot(
        data = bind_rows(peptides_per_protein_list, .id = "step") |>
          mutate(step = filtering_progress@steps[as.numeric(step)]),
        aes(
          x = factor(step, levels = filtering_progress@steps),
          y = n_peptides
        ),
        fill = "darkred",
        alpha = 0.5,
        outlier.shape = NA
      ) +
      labs(
        title = "Peptides per Protein Distribution",
        x = "Filtering Step",
        y = "Number of Peptides"
      ) +
      theme_minimal() +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid.major.x = element_blank()
      ) +
      coord_cartesian(
        ylim = c(
          0,
          quantile(bind_rows(peptides_per_protein_list)$n_peptides, 0.95)
        )
      )

    # Ensure consistent data types in peptides_per_run list
    peptides_per_run_list <- lapply(filtering_progress@peptides_per_run, function(df) {
      if (nrow(df) > 0) {
        df$Run <- as.character(df$Run)
        df$n_peptides <- as.numeric(df$n_peptides)
      }
      return(df)
    })

    p5 <- bind_rows(peptides_per_run_list, .id = "step") |>
      mutate(step = filtering_progress@steps[as.numeric(step)]) |>
      group_by(Run) |>
      mutate(avg_peptides = mean(n_peptides)) |>
      ungroup() |>
      # Run is already character from our preprocessing
      mutate(Run = fct_reorder(Run, avg_peptides)) |>
      ggplot(aes(
        x = Run, y = n_peptides,
        group = step,
        color = factor(step, levels = filtering_progress@steps)
      )) +
      geom_line() +
      geom_point() +
      labs(
        title = "Peptides per Run",
        x = "Run ID (ordered by average peptide count)",
        y = "Number of Peptides",
        color = "Step"
      ) +
      theme_minimal() +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid.major.x = element_blank()
      ) +
      scale_color_manual(values = get_color_palette(length(filtering_progress@steps), "forestgreen"))
  }

  # Create plot list based on data type
  plot_list <- list(
    proteins_total = p1,
    proteins_per_run = p4,
    peptides_total = p2,
    peptides_per_protein = p3,
    peptides_per_run = p5
  )

  # Save plots if derived_time_dir is valid and save_plots is TRUE
  if (save_plots) {
    message(sprintf("   [updateProteinFiltering] Saving individual plots to %s in formats: %s...", derived_time_dir, paste(formats, collapse = ", ")))
    for (plot_name in names(plot_list)) {
      message(sprintf("      Saving %s...", plot_name))
      for (fmt in formats) {
        filename <- file.path(
          derived_time_dir,
          sprintf("%s_%s.%s", step_name, plot_name, fmt)
        )
        tryCatch(
          {
            ggsave(filename,
              plot = plot_list[[plot_name]],
              width = 10,
              height = 8,
              dpi = 300
            )
          },
          error = function(e) message(sprintf("Warning: Failed to save %s as %s: %s", plot_name, fmt, e$message))
        )
      }
    }
  }

  # Return/display plots based on return_grid parameter
  if (return_grid) {
    message("   [updateProteinFiltering] Generating final grid...")
    message(sprintf("      Memory before grid: %s", format(sum(gc()[, 2]), units = "auto")))

    tryCatch(
      {
        if (!is_protein_quant || !all(is.na(filtering_progress@total_peptides))) {
          # Create full grid with all plots if peptide data exists
          message("      Combining all 5 plots...")
          grid1 <- gridExtra::arrangeGrob(p1, p2, p3, ncol = 3)
          grid2 <- gridExtra::arrangeGrob(p4, ncol = 1)
          grid3 <- gridExtra::arrangeGrob(p5, ncol = 1)

          # Use arrangeGrob to prevent immediate drawing, which might double-render
          grid_plot <- gridExtra::arrangeGrob(
            grid1,
            grid2,
            grid3,
            heights = c(1, 1, 1)
          )
        } else {
          # For protein_quant_table without peptide data, only show protein plots
          message("      Combining protein plots (p1, p4)...")
          grid_plot <- gridExtra::arrangeGrob(
            p1,
            p4,
            ncol = 1,
            heights = c(1, 1)
          )
        }

        message(sprintf("      Memory after grid creation: %s", format(sum(gc()[, 2]), units = "auto")))

        # Save the grid if derived_time_dir is valid and save_plots is TRUE
        if (save_plots) {
          message(sprintf("      Saving combined grid plot in formats: %s...", paste(formats, collapse = ", ")))
          for (fmt in formats) {
            filename <- file.path(
              derived_time_dir,
              sprintf("%s_combined_plots.%s", step_name, fmt)
            )
            ggsave(filename,
              plot = grid_plot,
              width = 15,
              height = if (!is_protein_quant || !all(is.na(filtering_progress@total_peptides))) 18 else 12,
              dpi = 300
            )
          }
        }

        if (!is.null(grid_plot)) {
          gridExtra::grid.arrange(grid_plot)
        }
        return(invisible(grid_plot))
      },
      error = function(e) {
        message(sprintf("ERROR in grid generation: %s", e$message))
        return(NULL)
      }
    )
  } else {
    # Print each plot individually
    message("   [updateProteinFiltering] Printing individual plots...")
    for (plot_obj in plot_list) {
      print(plot_obj)
    }
    # Return the list invisibly
    invisible(plot_list)
  }
}

