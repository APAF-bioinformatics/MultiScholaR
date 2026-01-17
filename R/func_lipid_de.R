# ============================================================================
# func_lipid_de.R
# ============================================================================
# Purpose: Lipidomics differential abundance analysis functions
# 
# This file contains functions for lipidomics differential abundance
# analysis, including limma-based analysis and result formatting. Functions
# in this file are used by lipidomics DE modules and related workflows.
#
# Functions to extract here:
# - differentialAbundanceAnalysis(): S4 method for DE analysis (lipid)
# - differentialAbundanceAnalysisHelper(): Helper for DE analysis
# - getCountsTable(): Gets counts table from object
# - Additional lipidomics DE helper functions
#
# Dependencies:
# - limma, edgeR
# - func_general_plotting.R (for visualization)
# - func_general_helpers.R (for utility functions)
# ============================================================================

# TODO: Extract the following functions from their current locations:

# Function 1: differentialAbundanceAnalysis() (lipid method)
# Current location: R/lipidVsSamplesS4Objects.R
# Type: S4 method (exportMethods)
# Description: Performs differential abundance analysis on lipids
# setMethod(f = "differentialAbundanceAnalysis", signature = "LipidomicsAssayData", ...) {
#   # Extract from R/lipidVsSamplesS4Objects.R
# }

# Function 2: differentialAbundanceAnalysisHelper()
# Current location: R/lipidVsSamplesS4Objects.R
# Type: S4 method (exportMethods)
# Description: Helper function for lipid DE analysis
# setMethod(f = "differentialAbundanceAnalysisHelper", ...) {
#   # Extract from R/lipidVsSamplesS4Objects.R
# }

# Function 3: getCountsTable()
# Current location: R/lipid_de_analysis_wrapper.R, R/lipids_de_analysis_wrapper.R
# Description: Gets counts table from object
# getCountsTable <- function(obj) {
#   # Extract from R/lipid_de_analysis_wrapper.R or R/lipids_de_analysis_wrapper.R
# }


# ----------------------------------------------------------------------------
# getCountsTable
# ----------------------------------------------------------------------------
# Helper function to get counts table
getCountsTable <- function(obj) {
    if (inherits(obj, "LipidomicsAssayData")) {
        message(sprintf("   Getting counts table for object of class: %s", class(obj)[1]))
        message(sprintf("   Returning lipid_data with dimensions: %d rows, %d cols"
            , nrow(obj@lipid_data), ncol(obj@lipid_data)))
        obj@lipid_data
    } else if (inherits(obj, "ProteinQuantitativeData")) {
        message(sprintf("   Returning protein_quant_table with dimensions: %d rows, %d cols"
            , nrow(obj@protein_quant_table), ncol(obj@protein_quant_table)))
        obj@protein_quant_table
    } else {
        message(sprintf("   ERROR: Unsupported object type: %s", class(obj)[1]))
        stop("Unsupported object type")
    }
}


# ============================================================================
# METABOLOMICS DIFFERENTIAL EXPRESSION ANALYSIS FUNCTIONS
# ============================================================================
# These functions provide differential expression analysis for lipidomics
# data stored as multiple assays (e.g., LCMS_Pos, LCMS_Neg) in the
# LipidomicsAssayData S4 object.
# ============================================================================


# ----------------------------------------------------------------------------
# runTestsContrastsMetab
# ----------------------------------------------------------------------------
#' Run limma differential expression analysis for a single assay
#'
#' @description Performs limma-based differential expression analysis on a
#'   single lipidomics assay matrix. This is the core DE engine adapted from
#'   the proteomics `runTestsContrasts()` function.
#'
#' @param data Numeric matrix with lipids as rows and samples as columns.
#'   Rownames should be lipid IDs.
#' @param contrast_strings Character vector of contrast specifications
#'   (e.g., "Treatment-Control").
#' @param design_matrix Data frame containing sample metadata with sample IDs
#'   as a column (will be converted to rownames).
#' @param formula_string Model formula string (e.g., "~ 0 + group").
#' @param sample_id_col Name of the sample ID column in design_matrix.
#' @param treat_lfc_cutoff Optional log fold-change cutoff for TREAT analysis.
#'   Set to NA for standard analysis.
#' @param eBayes_trend Logical, whether to use trend in eBayes.
#' @param eBayes_robust Logical, whether to use robust eBayes.
#'
#' @return A list containing:
#'   - results: Named list of DE result tables (one per contrast)
#'   - fit.eb: The eBayes fit object
#'   - qvalue_warnings: List of contrasts where qvalue() failed
#'
#' @importFrom limma lmFit contrasts.fit eBayes treat topTable topTreat makeContrasts
#' @importFrom qvalue qvalue
#' @importFrom stats model.matrix model.frame as.formula p.adjust
#' @importFrom dplyr mutate rename
#' @importFrom rlang ensym as_name
#' @importFrom purrr map
#' @export
runTestsContrastsMetab <- function(
    data
    , contrast_strings
    , design_matrix
    , formula_string
    , sample_id_col = "Run"
    , treat_lfc_cutoff = NA
    , eBayes_trend = TRUE
    , eBayes_robust = TRUE
) {
    logger::log_info("--- Entering runTestsContrastsMetab ---")
    logger::log_info(sprintf("   data dims = %d x %d, %d contrasts"
        , nrow(data), ncol(data), length(contrast_strings)))
    logger::log_info(sprintf("   contrasts = %s", paste(contrast_strings, collapse = ", ")))

    # Prepare design matrix with sample IDs as rownames
    dm <- design_matrix
    if (sample_id_col %in% colnames(dm)) {
        rownames(dm) <- dm[[sample_id_col]]
    }

    # Create formula and model matrix
    ff <- stats::as.formula(formula_string)
    mod_frame <- stats::model.frame(ff, dm)
    design_m <- stats::model.matrix(ff, mod_frame)
    logger::log_info(sprintf("   design_m dims = %d x %d", nrow(design_m), ncol(design_m)))

    # [D66:START] ─────────────────────────
    d66_log <- function(...) message(sprintf("[D66] %s", paste0(...)))
    d66_log("  runTestsContrastsMetab - Model matrix created:")
    d66_log("    design_m levels (colnames) = ", paste(colnames(design_m), collapse = ", "))
    d66_log("    contrast_strings = ", paste(contrast_strings, collapse = ", "))
    # [D66:END] ───────────────────────────

    # Subset data to match design matrix samples
    common_samples <- intersect(colnames(data), rownames(design_m))
    if (length(common_samples) == 0) {
        stop("No common samples between data matrix and design matrix")
    }
    data_subset <- data[, common_samples, drop = FALSE]
    design_m <- design_m[common_samples, , drop = FALSE]
    logger::log_info(sprintf("   data_subset dims = %d x %d", nrow(data_subset), ncol(data_subset)))

    # [D66:START] ─────────────────────────
    d66_log("  About to call makeContrasts:")
    d66_log("    contrasts = ", paste(contrast_strings, collapse = ", "))
    d66_log("    levels = ", paste(colnames(design_m), collapse = ", "))
    # [D66:END] ───────────────────────────

    # =====================================================================
    # VALIDATION: Check that contrast terms exist in design matrix levels
    # This catches the common bug where contrasts are saved without the
    # 'group' prefix (e.g., "Treatment-Control" instead of "groupTreatment-groupControl")
    # =====================================================================
    available_levels <- colnames(design_m)
    for (cs in contrast_strings) {
        # Extract terms from contrast string (split on - or +)
        terms <- unlist(strsplit(cs, "[-+]"))
        terms <- trimws(terms)
        terms <- terms[nzchar(terms)]  # Remove empty strings

        missing_terms <- setdiff(terms, available_levels)
        if (length(missing_terms) > 0) {
            error_msg <- sprintf(
                paste0(
                    "Contrast '%s' references undefined levels: [%s].\n"
                    , "Available levels from model matrix: [%s].\n\n"
                    , "This usually means the contrast was saved with a different formula than currently used.\n"
                    , "If using formula '~ 0 + group', contrasts must be in format 'groupX-groupY'.\n"
                    , "Check that the design builder formula matches the DE module formula."
                )
                , cs
                , paste(missing_terms, collapse = ", ")
                , paste(available_levels, collapse = ", ")
            )
            d66_log("  VALIDATION FAILED: ", error_msg)
            stop(error_msg)
        }
    }
    d66_log("  VALIDATION PASSED: All contrast terms exist in design matrix levels")

    # Create contrast matrix
    logger::log_info("   Creating contrast matrix...")
    contr.matrix <- limma::makeContrasts(
        contrasts = contrast_strings
        , levels = colnames(design_m)
    )

    # Run limma analysis
    logger::log_info("   Running lmFit...")
    fit <- limma::lmFit(data_subset, design = design_m)

    logger::log_info("   Running contrasts.fit...")
    cfit <- limma::contrasts.fit(fit, contrasts = contr.matrix)

    logger::log_info("   Running eBayes...")
    eb.fit <- limma::eBayes(cfit, trend = eBayes_trend, robust = eBayes_robust)

    # Run treat or standard analysis
    qvalue_failures <- list()

    if (!is.na(treat_lfc_cutoff) && treat_lfc_cutoff > 0) {
        logger::log_info("   Running TREAT analysis with LFC cutoff...")
        t.fit <- limma::treat(eb.fit, lfc = as.double(treat_lfc_cutoff))
        extract_fn <- function(fit, coef) limma::topTreat(fit, coef = coef, n = Inf)
    } else {
        logger::log_info("   Running standard analysis...")
        t.fit <- eb.fit
        extract_fn <- function(fit, coef) limma::topTable(fit, coef = coef, n = Inf)
    }

    # Process each contrast
    result_tables <- purrr::map(contrast_strings, function(contrast) {
        logger::log_info(sprintf("      Processing contrast: %s", contrast))
        qvalue_failed <- FALSE

        tryCatch({
            de_tbl <- extract_fn(t.fit, contrast)
            logger::log_info(sprintf("      Extracted %d rows", nrow(de_tbl)))

            # Safe qvalue computation
            valid_p_idx <- which(!is.na(de_tbl$P.Value) & is.finite(de_tbl$P.Value))
            q_values_all <- rep(NA_real_, nrow(de_tbl))
            fdr_values_all <- rep(NA_real_, nrow(de_tbl))

            if (length(valid_p_idx) > 0) {
                valid_p_values <- de_tbl$P.Value[valid_p_idx]

                # Edge case checks
                all_zeros <- all(valid_p_values == 0)
                all_ones <- all(valid_p_values == 1)
                too_few <- length(valid_p_values) < 3

                use_qvalue <- !all_zeros && !all_ones && !too_few

                if (use_qvalue) {
                    tryCatch({
                        q_values_all[valid_p_idx] <- qvalue::qvalue(valid_p_values)$qvalues
                        logger::log_info("      qvalue() computation successful")
                    }, error = function(e) {
                        qvalue_failed <<- TRUE
                        logger::log_warn(sprintf("      qvalue() failed: %s, using BH", e$message))
                        q_values_all[valid_p_idx] <<- stats::p.adjust(valid_p_values, method = "BH")
                    })
                } else {
                    qvalue_failed <- TRUE
                    q_values_all[valid_p_idx] <- stats::p.adjust(valid_p_values, method = "BH")
                    logger::log_info("      Using p.adjust() due to edge case")
                }

                # Always compute BH-adjusted FDR
                fdr_values_all[valid_p_idx] <- stats::p.adjust(valid_p_values, method = "BH")
            }

            # Add columns to results
            de_tbl$fdr_qvalue <- q_values_all
            de_tbl$fdr_value_bh <- fdr_values_all
            de_tbl$raw_pvalue <- de_tbl$P.Value

            if (qvalue_failed) {
                qvalue_failures[[contrast]] <<- TRUE
            }

            return(de_tbl)

        }, error = function(e) {
            logger::log_error(sprintf("      ERROR in contrast %s: %s", contrast, e$message))
            stop(e)
        })
    })

    names(result_tables) <- contrast_strings
    logger::log_info("--- Exiting runTestsContrastsMetab ---")

    return(list(
        results = result_tables
        , fit.eb = t.fit
        , qvalue_warnings = qvalue_failures
    ))
}


# ----------------------------------------------------------------------------
# runLipidsDE
# ----------------------------------------------------------------------------
#' Run differential expression analysis on all assays in a LipidomicsAssayData object
#'
#' @description Main entry point for lipidomics differential expression analysis.
#'   Loops over all assays in the LipidomicsAssayData object, runs limma DE on each,
#'   and aggregates results with an assay identifier column.
#'
#' @param theObject LipidomicsAssayData S4 object containing lipid data.
#' @param contrasts_tbl Data frame with contrast definitions. Must have a column
#'   named "contrasts" or "contrast_string" containing the contrast formulas.
#' @param formula_string Model formula (default "~ 0 + group").
#' @param de_q_val_thresh Q-value threshold for significance (default 0.05).
#' @param treat_lfc_cutoff Log fold-change cutoff for TREAT (default 0, standard analysis).
#' @param eBayes_trend Logical, use trend in eBayes (default TRUE).
#' @param eBayes_robust Logical, use robust eBayes (default TRUE).
#'
#' @return A list containing:
#'   - theObject: The input S4 object (unchanged)
#'   - contrasts_results: Per-assay limma fit objects
#'   - de_lipids_long: Combined long-format results with assay column
#'   - per_assay_results: List of per-assay result tables
#'   - significant_counts: Summary of significant lipids per assay
#'   - qvalue_warnings: Any qvalue computation warnings
#'
#' @importFrom tibble rownames_to_column
#' @importFrom dplyr mutate bind_rows select left_join
#' @importFrom purrr map map2 set_names
#' @importFrom logger log_info log_error log_warn
#' @export
runLipidsDE <- function(
    theObject
    , contrasts_tbl
    , formula_string = "~ 0 + group"
    , de_q_val_thresh = 0.05
    , treat_lfc_cutoff = 0
    , eBayes_trend = TRUE
    , eBayes_robust = TRUE
) {
    # [D66:START] ─────────────────────────
    d66_log <- function(...) message(sprintf("[D66] %s", paste0(...)))
    d66_log("=== ENTER runLipidsDE ===")
    d66_log("  formula_string = ", formula_string)
    d66_log("  contrasts_tbl columns = ", paste(colnames(contrasts_tbl), collapse = ", "))
    d66_log("  contrasts_tbl$contrasts = ", paste(contrasts_tbl$contrasts, collapse = ", "))
    # [D66:END] ───────────────────────────

    logger::log_info("=== Starting runLipidsDE ===")

    # Validate input
    if (!inherits(theObject, "LipidomicsAssayData")) {
        stop("theObject must be a LipidomicsAssayData S4 object")
    }

    # Extract object slots
    assay_list <- theObject@lipid_data
    design_matrix <- theObject@design_matrix
    sample_id_col <- theObject@sample_id
    group_col <- theObject@group_id
    lipid_id_col <- theObject@lipid_id_column
    annotation_col <- theObject@annotation_id_column

    # [D66:START] ─────────────────────────
    d66_log("  S4 slots extracted:")
    d66_log("    sample_id_col = ", sample_id_col)
    d66_log("    group_col = ", group_col)
    d66_log("    design_matrix columns = ", paste(colnames(design_matrix), collapse = ", "))
    d66_log("    design_matrix nrow = ", nrow(design_matrix))
    d66_log("    assay_list names = ", paste(names(assay_list), collapse = ", "))
    # [D66:END] ───────────────────────────

    # Extract contrast strings
    if ("contrasts" %in% colnames(contrasts_tbl)) {
        contrast_strings <- as.character(contrasts_tbl$contrasts)
    } else if ("contrast_string" %in% colnames(contrasts_tbl)) {
        contrast_strings <- as.character(contrasts_tbl$contrast_string)
    } else {
        stop("contrasts_tbl must have a 'contrasts' or 'contrast_string' column")
    }

    # Get friendly names if available
    friendly_names <- if ("friendly_names" %in% colnames(contrasts_tbl)) {
        as.character(contrasts_tbl$friendly_names)
    } else {
        contrast_strings
    }

    # [D66:START] ─────────────────────────
    d66_log("  contrast_strings = ", paste(contrast_strings, collapse = ", "))
    d66_log("  friendly_names = ", paste(friendly_names, collapse = ", "))
    # [D66:END] ───────────────────────────

    assay_names <- names(assay_list)
    logger::log_info(sprintf("   Processing %d assays: %s"
        , length(assay_names), paste(assay_names, collapse = ", ")))
    logger::log_info(sprintf("   Contrasts: %s", paste(contrast_strings, collapse = ", ")))

    # Initialize results storage
    per_assay_results <- list()
    contrasts_results <- list()
    all_qvalue_warnings <- list()

    # Process each assay
    for (assay_name in assay_names) {
        logger::log_info(sprintf("--- Processing assay: %s ---", assay_name))

        assay_data <- assay_list[[assay_name]]

        # Build expression matrix (lipids as rows, samples as columns)
        # First, identify sample columns (those that match design matrix)
        sample_cols <- intersect(colnames(assay_data), design_matrix[[sample_id_col]])

        if (length(sample_cols) == 0) {
            logger::log_warn(sprintf("   No matching samples for assay: %s, skipping", assay_name))
            next
        }

        # Extract numeric data and set rownames to lipid IDs
        expr_matrix <- as.matrix(assay_data[, sample_cols, drop = FALSE])
        rownames(expr_matrix) <- assay_data[[lipid_id_col]]

        # Get annotation info for later joining
        annotation_info <- data.frame(
            lipid_id = assay_data[[lipid_id_col]]
            , lipid_name = assay_data[[annotation_col]]
            , stringsAsFactors = FALSE
        )

        logger::log_info(sprintf("   Expression matrix: %d lipids x %d samples"
            , nrow(expr_matrix), ncol(expr_matrix)))

        # Run limma DE
        tryCatch({
            limma_results <- runTestsContrastsMetab(
                data = expr_matrix
                , contrast_strings = contrast_strings
                , design_matrix = design_matrix
                , formula_string = formula_string
                , sample_id_col = sample_id_col
                , treat_lfc_cutoff = if (treat_lfc_cutoff > 0) treat_lfc_cutoff else NA
                , eBayes_trend = eBayes_trend
                , eBayes_robust = eBayes_robust
            )

            # Store fit object for Glimma
            contrasts_results[[assay_name]] <- limma_results

            # Collect qvalue warnings
            if (length(limma_results$qvalue_warnings) > 0) {
                all_qvalue_warnings[[assay_name]] <- limma_results$qvalue_warnings
            }

            # Process results for each contrast
            assay_results_list <- purrr::map2(
                limma_results$results
                , names(limma_results$results)
                , function(de_tbl, contrast_name) {
                    # Add lipid_id from rownames
                    de_tbl <- tibble::rownames_to_column(de_tbl, "lipid_id")

                    # Add metadata columns
                    de_tbl$assay <- assay_name
                    de_tbl$comparison <- contrast_name

                    # Add friendly name
                    idx <- match(contrast_name, contrast_strings)
                    de_tbl$friendly_name <- if (!is.na(idx)) friendly_names[idx] else contrast_name

                    # Join annotation info
                    de_tbl <- dplyr::left_join(
                        de_tbl
                        , annotation_info
                        , by = "lipid_id"
                    )

                    # Classify significance
                    de_tbl$significant <- ifelse(
                        de_tbl$fdr_qvalue < de_q_val_thresh & abs(de_tbl$logFC) >= treat_lfc_cutoff
                        , ifelse(de_tbl$logFC > 0, "Up", "Down")
                        , "NS"
                    )

                    # Add sample intensity columns using createMetabDeResultsLongFormat
                    cat(sprintf("[D66] Calling createMetabDeResultsLongFormat for contrast: %s\n", contrast_name))
                    de_tbl <- createMetabDeResultsLongFormat(
                        lfc_qval_tbl = de_tbl
                        , expr_matrix = expr_matrix
                        , design_matrix = design_matrix
                        , sample_id_col = sample_id_col
                        , group_id_col = group_col
                        , lipid_id_col = "lipid_id"
                    )

                    return(de_tbl)
                }
            )

            # Combine all contrasts for this assay
            assay_combined <- dplyr::bind_rows(assay_results_list)
            per_assay_results[[assay_name]] <- assay_combined

            logger::log_info(sprintf("   Completed assay %s: %d total results"
                , assay_name, nrow(assay_combined)))

        }, error = function(e) {
            cat(sprintf("[ERROR] Processing assay %s: %s\n", assay_name, e$message))
            logger::log_error(sprintf("   ERROR processing assay %s: %s", assay_name, e$message))
            per_assay_results[[assay_name]] <<- NULL
        })
    }

    # Combine all assay results
    de_lipids_long <- dplyr::bind_rows(per_assay_results)
    logger::log_info(sprintf("   Combined results: %d total rows", nrow(de_lipids_long)))

    # Calculate significant counts per assay
    significant_counts <- purrr::map(per_assay_results, function(df) {
        if (is.null(df) || nrow(df) == 0) {
            return(list(up = 0, down = 0, ns = 0))
        }
        list(
            up = sum(df$significant == "Up", na.rm = TRUE)
            , down = sum(df$significant == "Down", na.rm = TRUE)
            , ns = sum(df$significant == "NS", na.rm = TRUE)
        )
    }) |> purrr::set_names(names(per_assay_results))

    logger::log_info("=== Completed runLipidsDE ===")

    return(list(
        theObject = theObject
        , contrasts_results = contrasts_results
        , de_lipids_long = de_lipids_long
        , per_assay_results = per_assay_results
        , significant_counts = significant_counts
        , qvalue_warnings = all_qvalue_warnings
        , de_q_val_thresh = de_q_val_thresh
        , treat_lfc_cutoff = treat_lfc_cutoff
    ))
}


# ----------------------------------------------------------------------------
# createMetabDeResultsLongFormat
# ----------------------------------------------------------------------------
#' Create lipidomics DE results in long format with sample intensity columns
#'
#' @description Creates a long-format DE results table that includes individual
#'   sample intensity values, mirroring the proteomics `createDeResultsLongFormat()`.
#'   The output includes columns for each sample's intensity value, named as
#'   `intensity.{sample_id}.{group}`.
#'
#' @param lfc_qval_tbl Data frame with DE statistics (logFC, pvalue, fdr_qvalue).
#'   Must have a lipid_id column and a comparison column containing the
#'   contrast string in format "groupA-groupB".
#' @param expr_matrix Expression matrix with lipids as rows, samples as columns.
#'   Rownames must match lipid_id column in lfc_qval_tbl.
#' @param design_matrix Design matrix with sample and group assignments.
#' @param sample_id_col Name of sample ID column in design_matrix.
#' @param group_id_col Name of group ID column in design_matrix.
#' @param lipid_id_col Name of lipid ID column.
#'
#' @return Data frame with DE results plus intensity columns for each sample.
#'
#' @importFrom dplyr left_join mutate select rename arrange distinct filter pull
#' @importFrom tidyr pivot_longer pivot_wider separate_wider_delim
#' @importFrom tibble rownames_to_column
#' @importFrom rlang sym set_names
#' @importFrom stringr str_replace_all
#' @export
createMetabDeResultsLongFormat <- function(
    lfc_qval_tbl
    , expr_matrix
    , design_matrix
    , sample_id_col = "sample_id"
    , group_id_col = "group"
    , lipid_id_col = "lipid_id"
) {
    logger::log_info("--- Entering createMetabDeResultsLongFormat ---")
    logger::log_info(sprintf("   lfc_qval_tbl: %d rows", nrow(lfc_qval_tbl)))
    logger::log_info(sprintf("   expr_matrix: %d lipids x %d samples"
        , nrow(expr_matrix), ncol(expr_matrix)))

    # Convert expression matrix to long format with intensity values
    intensity_long <- expr_matrix |>
        as.data.frame() |>
        tibble::rownames_to_column(lipid_id_col) |>
        tidyr::pivot_longer(
            cols = -!!rlang::sym(lipid_id_col)
            , names_to = sample_id_col
            , values_to = "intensity"
        ) |>
        dplyr::left_join(design_matrix, by = sample_id_col)

    logger::log_info(sprintf("   intensity_long: %d rows", nrow(intensity_long)))

    # Create intensity columns per group
    # Format: intensity.{sample_id}.{group}
    intensity_wide <- intensity_long |>
        dplyr::mutate(
            col_name = paste0("intensity.", !!rlang::sym(sample_id_col), ".", !!rlang::sym(group_id_col))
        ) |>
        dplyr::select(!!rlang::sym(lipid_id_col), col_name, intensity) |>
        tidyr::pivot_wider(
            id_cols = !!rlang::sym(lipid_id_col)
            , names_from = col_name
            , values_from = intensity
        )

    logger::log_info(sprintf("   intensity_wide: %d lipids x %d columns"
        , nrow(intensity_wide), ncol(intensity_wide)))

    # Parse contrast to get numerator/denominator groups
    # Contrast format is "groupA-groupB" where groupA is numerator, groupB is denominator
    if ("comparison" %in% colnames(lfc_qval_tbl)) {
        # Extract unique comparisons
        contrasts <- unique(lfc_qval_tbl$comparison)
        logger::log_info(sprintf("   Processing %d contrasts: %s"
            , length(contrasts), paste(contrasts, collapse = ", ")))

        # Parse first contrast to get group prefix pattern
        # Assumes format like "groupTreatment-groupControl"
        first_contrast <- contrasts[1]
        parts <- strsplit(first_contrast, "-")[[1]]

        if (length(parts) == 2) {
            # Try to detect group prefix (e.g., "group" from "groupTreatment")
            left_part <- parts[1]
            right_part <- parts[2]

            # Add numerator/denominator columns by parsing comparison
            tryCatch({
                lfc_qval_tbl <- lfc_qval_tbl |>
                    tidyr::separate_wider_delim(
                        comparison
                        , delim = "-"
                        , names = c("numerator", "denominator")
                        , cols_remove = FALSE
                    )
            }, error = function(e) {
                cat(sprintf("[WARN] Could not parse contrast groups: %s\n", e$message))
                lfc_qval_tbl$numerator <<- NA_character_
                lfc_qval_tbl$denominator <<- NA_character_
            })
        }
    }

    # Join DE results with intensity values
    de_results_long <- lfc_qval_tbl |>
        dplyr::left_join(intensity_wide, by = lipid_id_col) |>
        dplyr::arrange(comparison, fdr_qvalue, logFC) |>
        dplyr::distinct()

    logger::log_info(sprintf("   Final de_results_long: %d rows x %d columns"
        , nrow(de_results_long), ncol(de_results_long)))

    logger::log_info("--- Exiting createMetabDeResultsLongFormat ---")

    return(de_results_long)
}


# ----------------------------------------------------------------------------
# getLipidQuantData
# ----------------------------------------------------------------------------
#' Extract quantitative data columns from an assay data frame
#'
#' @description Helper function to separate lipid quantitative data
#'   (sample columns) from metadata columns (ID, annotation, etc.).
#'
#' @param assay_df Data frame containing lipid data for one assay.
#' @param lipid_id_col Name of the lipid ID column.
#' @param annotation_col Name of the annotation column.
#' @param additional_meta_cols Additional columns to exclude from quant data.
#'
#' @return A list with:
#'   - quant_data: Data frame with only sample columns
#'   - meta_data: Data frame with only metadata columns
#'   - sample_cols: Names of sample columns
#'
#' @export
getLipidQuantData <- function(
    assay_df
    , lipid_id_col = "Alignment ID"
    , annotation_col = "Lipid name"
    , additional_meta_cols = NULL
) {
    # Identify metadata columns to exclude
    meta_cols <- c(lipid_id_col, annotation_col)
    if (!is.null(additional_meta_cols)) {
        meta_cols <- c(meta_cols, additional_meta_cols)
    }

    # Get sample columns (everything that's not metadata)
    all_cols <- colnames(assay_df)
    sample_cols <- setdiff(all_cols, meta_cols)

    # Also exclude any obviously non-numeric columns
    sample_cols <- sample_cols[sapply(assay_df[, sample_cols, drop = FALSE], is.numeric)]

    list(
        quant_data = assay_df[, sample_cols, drop = FALSE]
        , meta_data = assay_df[, intersect(meta_cols, all_cols), drop = FALSE]
        , sample_cols = sample_cols
    )
}


# ----------------------------------------------------------------------------
# generateMetabVolcanoPlotGlimma
# ----------------------------------------------------------------------------
#' Generate interactive Glimma volcano plot for lipidomics DE results
#'
#' @description Creates an interactive volcano plot using the Glimma package
#'   for lipidomics differential expression results. Supports per-assay
#'   or combined viewing.
#'
#' @param de_results_list Results list from `runLipidsDE()`.
#' @param selected_contrast Contrast to display (from friendly_name or comparison).
#' @param selected_assay Optional assay to filter (NULL for combined view).
#' @param de_q_val_thresh Q-value threshold for significance marking.
#' @param lipid_id_column Column name for lipid IDs.
#' @param annotation_column Column name for lipid annotations (for labels).
#'
#' @return A Glimma HTML widget or NULL if generation fails.
#'
#' @importFrom Glimma glimmaVolcano
#' @importFrom dplyr filter mutate case_when select
#' @importFrom rlang sym
#' @importFrom stringr str_extract
#' @importFrom logger log_info log_error log_warn
#' @export
generateMetabVolcanoPlotGlimma <- function(
    de_results_list
    , selected_contrast = NULL
    , selected_assay = NULL
    , de_q_val_thresh = 0.05
    , lipid_id_column = "lipid_id"
    , annotation_column = "lipid_name"
) {
    # [D66:START] ─────────────────────────
    d66_log <- function(...) message(sprintf("[D66] %s", paste0(...)))
    d66_log("=== ENTER generateMetabVolcanoPlotGlimma ===")
    d66_log("  selected_contrast = ", if(is.null(selected_contrast)) "NULL" else selected_contrast)
    d66_log("  selected_assay = ", if(is.null(selected_assay)) "NULL" else selected_assay)
    # [D66:END] ───────────────────────────

    logger::log_info("--- Entering generateMetabVolcanoPlotGlimma ---")
    logger::log_info(sprintf("   selected_contrast = %s", selected_contrast))
    logger::log_info(sprintf("   selected_assay = %s", ifelse(is.null(selected_assay), "NULL (combined)", selected_assay)))

    if (is.null(de_results_list) || is.null(de_results_list$de_lipids_long)) {
        # [D66:START]
        d66_log("  ERROR: de_results_list or de_lipids_long is NULL")
        d66_log("    de_results_list is NULL = ", is.null(de_results_list))
        if (!is.null(de_results_list)) {
            d66_log("    de_results_list names = ", paste(names(de_results_list), collapse = ", "))
            d66_log("    de_lipids_long is NULL = ", is.null(de_results_list$de_lipids_long))
        }
        # [D66:END]
        logger::log_warn("   No DE results available")
        return(NULL)
    }

    if (is.null(selected_contrast)) {
        logger::log_warn("   No contrast selected")
        return(NULL)
    }

    # Get data
    de_lipids_long <- de_results_list$de_lipids_long
    contrasts_results <- de_results_list$contrasts_results

    # [D66:START] ─────────────────────────
    d66_log("  de_lipids_long dims = ", nrow(de_lipids_long), " x ", ncol(de_lipids_long))
    d66_log("  de_lipids_long columns = ", paste(colnames(de_lipids_long), collapse = ", "))
    if (nrow(de_lipids_long) > 0 && "comparison" %in% colnames(de_lipids_long)) {
        d66_log("  unique comparisons = ", paste(unique(de_lipids_long$comparison), collapse = ", "))
    }
    if (nrow(de_lipids_long) > 0 && "friendly_name" %in% colnames(de_lipids_long)) {
        d66_log("  unique friendly_names = ", paste(unique(de_lipids_long$friendly_name), collapse = ", "))
    }
    # [D66:END] ───────────────────────────

    # Extract comparison name (handle "=" format from full_format)
    comparison_to_search <- stringr::str_extract(selected_contrast, "^[^=]+")
    if (is.na(comparison_to_search)) {
        comparison_to_search <- selected_contrast
    }

    # [D66:START]
    d66_log("  comparison_to_search = ", comparison_to_search)
    # [D66:END]

    # Glimma requires matching dimensions between fit object and data.
    # For "Combined" view across multiple assays, this is not possible since
    # each assay has its own fit object with different genes.
    # Return NULL with a message - user should use static plot for combined view.
    if (is.null(selected_assay) || selected_assay == "Combined") {
        d66_log("  NOTE: Combined view not supported for Glimma - use static plot")
        logger::log_info("   Glimma not supported for Combined view (dimension mismatch). Use static plot.")
        return(NULL)
    }

    # Filter for selected contrast
    contrast_data <- de_lipids_long |>
        dplyr::filter(comparison == comparison_to_search | friendly_name == comparison_to_search)

    # Filter by selected assay (required for Glimma)
    contrast_data <- contrast_data |>
        dplyr::filter(assay == selected_assay)

    if (nrow(contrast_data) == 0) {
        logger::log_warn(sprintf("   No data found for contrast: %s", comparison_to_search))
        return(NULL)
    }

    logger::log_info(sprintf("   Found %d rows for contrast", nrow(contrast_data)))

    # Determine which assay to use for the fit object
    if (!is.null(selected_assay) && selected_assay != "Combined" && selected_assay %in% names(contrasts_results)) {
        fit_assay <- selected_assay
    } else {
        # Use first available assay for combined view
        fit_assay <- names(contrasts_results)[1]
    }

    if (is.null(contrasts_results[[fit_assay]])) {
        logger::log_warn(sprintf("   No fit object for assay: %s", fit_assay))
        return(NULL)
    }

    fit_obj <- contrasts_results[[fit_assay]]$fit.eb

    # [D66:START] ─────────────────────────
    d66_log("  fit_assay = ", fit_assay)
    d66_log("  fit_obj class = ", class(fit_obj)[1])
    # [D66:END] ───────────────────────────

    # Find coefficient index
    # NOTE: coef_names contains the actual contrast strings (e.g., "groupTreatment-groupControl")
    # but comparison_to_search is the friendly name (e.g., "Treatment_vs_Control")
    # We need to get the actual contrast string from the data
    coef_names <- colnames(fit_obj$coefficients)

    # [D66:START] ─────────────────────────
    d66_log("  Looking for coefficient:")
    d66_log("    comparison_to_search (friendly) = ", comparison_to_search)
    d66_log("    coef_names (actual) = ", paste(coef_names, collapse = ", "))
    # [D66:END] ───────────────────────────

    # First, try to find the actual contrast string from the filtered data
    # The 'comparison' column contains the actual limma contrast string
    actual_contrast_string <- NULL
    if (nrow(contrast_data) > 0 && "comparison" %in% colnames(contrast_data)) {
        actual_contrast_string <- unique(contrast_data$comparison)[1]
        d66_log("    actual_contrast_string from data = ", actual_contrast_string)
    }

    # Try to find coefficient by actual contrast string first
    coef_index <- integer(0)
    if (!is.null(actual_contrast_string)) {
        coef_index <- which(coef_names == actual_contrast_string)
        d66_log("    Match by actual_contrast_string: coef_index = ", paste(coef_index, collapse = ", "))
    }

    # Fallback: try friendly name (for backwards compatibility)
    if (length(coef_index) == 0) {
        coef_index <- which(coef_names == comparison_to_search)
        d66_log("    Match by friendly name: coef_index = ", paste(coef_index, collapse = ", "))
    }

    # Fallback: try pattern match on friendly name
    if (length(coef_index) == 0) {
        coef_index <- grep(paste0("^", comparison_to_search), coef_names)
        d66_log("    Match by pattern: coef_index = ", paste(coef_index, collapse = ", "))
    }

    if (length(coef_index) == 0) {
        d66_log("  ERROR: No coefficient found!")
        logger::log_warn(sprintf("   No coefficient found for: %s", comparison_to_search))
        logger::log_info(sprintf("   Available coefficients: %s", paste(coef_names, collapse = ", ")))
        return(NULL)
    }

    coef_index <- coef_index[1]
    d66_log("  FINAL coef_index = ", coef_index, " (", coef_names[coef_index], ")")
    logger::log_info(sprintf("   Using coefficient index %d: %s", coef_index, coef_names[coef_index]))

    # Prepare volcano plot annotation table
    volcano_tab <- contrast_data |>
        dplyr::mutate(
            lqm = -log10(fdr_qvalue)
            , label = dplyr::case_when(
                abs(logFC) >= 1 & fdr_qvalue >= de_q_val_thresh ~ "Not sig., |logFC| >= 1"
                , abs(logFC) >= 1 & fdr_qvalue < de_q_val_thresh ~ "Sig., |logFC| >= 1"
                , abs(logFC) < 1 & fdr_qvalue < de_q_val_thresh ~ "Sig., |logFC| < 1"
                , TRUE ~ "Not sig."
            )
            , display_name = ifelse(
                !is.na(lipid_name) & lipid_name != ""
                , lipid_name
                , lipid_id
            )
        ) |>
        dplyr::select(
            lipid_id
            , lipid_name
            , display_name
            , assay
            , logFC
            , raw_pvalue
            , fdr_qvalue
            , lqm
            , label
            , significant
        )

    # Get counts matrix for the selected assay
    theObject <- de_results_list$theObject
    assay_data <- theObject@lipid_data[[fit_assay]]
    sample_cols <- intersect(colnames(assay_data), theObject@design_matrix[[theObject@sample_id]])

    counts_mat <- as.matrix(assay_data[, sample_cols, drop = FALSE])
    rownames(counts_mat) <- assay_data[[theObject@lipid_id_column]]

    # Get groups
    dm <- theObject@design_matrix
    rownames(dm) <- dm[[theObject@sample_id]]
    groups <- dm[sample_cols, theObject@group_id]

    # [D66:START] ─────────────────────────
    d66_log("  Preparing Glimma inputs:")
    d66_log("    counts_mat dims = ", nrow(counts_mat), " x ", ncol(counts_mat))
    d66_log("    counts_mat rownames (first 5) = ", paste(head(rownames(counts_mat), 5), collapse = ", "))
    d66_log("    groups = ", paste(groups, collapse = ", "))
    d66_log("    volcano_tab nrow = ", nrow(volcano_tab))
    d66_log("    fit_obj nrow = ", nrow(fit_obj$coefficients))
    d66_log("    coef_index = ", coef_index)
    # [D66:END] ───────────────────────────

    # Generate Glimma widget
    tryCatch({
        logger::log_info("   Generating Glimma widget...")

        # [D66:START]
        d66_log("  Calling Glimma::glimmaVolcano...")
        # [D66:END]

        # Use Glimma's glimmaVolcano function directly
        # NOTE: transform.counts = "none" because lipidomics data is already
        # log-transformed and contains negative values. Glimma's default is to
        # apply cpm(log=TRUE) which fails on negative values.
        glimma_widget <- Glimma::glimmaVolcano(
            x = fit_obj
            , coef = coef_index
            , counts = counts_mat
            , groups = groups
            , anno = data.frame(
                ID = volcano_tab$lipid_id
                , Name = volcano_tab$display_name
                , Assay = volcano_tab$assay
            )
            , status = ifelse(volcano_tab$significant == "Up", 1
                , ifelse(volcano_tab$significant == "Down", -1, 0))
            , main = paste("Volcano Plot:", selected_contrast)
            , transform.counts = "none"
            , html = NULL
        )

        # [D66:START]
        d66_log("  Glimma widget created successfully!")
        d66_log("    widget class = ", class(glimma_widget)[1])
        # [D66:END]

        logger::log_info("--- Exiting generateMetabVolcanoPlotGlimma (success) ---")
        return(glimma_widget)

    }, error = function(e) {
        # [D66:START]
        d66_log("  ERROR in Glimma::glimmaVolcano: ", e$message)
        d66_log("  Error call: ", deparse(e$call))
        # [D66:END]
        logger::log_error(sprintf("   Glimma widget generation failed: %s", e$message))
        return(NULL)
    })
}


# ----------------------------------------------------------------------------
# generateMetabDEHeatmap
# ----------------------------------------------------------------------------
#' Generate heatmap for lipidomics DE results
#'
#' @description Creates a heatmap of top differentially expressed lipids
#'   with customizable clustering and scaling options.
#'
#' @param de_results_list Results list from `runLipidsDE()`.
#' @param selected_contrast Contrast to display.
#' @param selected_assay Optional assay filter (NULL for combined).
#' @param top_n Number of top lipids to include (by |logFC|).
#' @param clustering_method Hierarchical clustering method.
#' @param distance_method Distance metric for clustering.
#' @param cluster_rows Logical, cluster rows.
#' @param cluster_cols Logical, cluster columns.
#' @param scale_data Scaling option: "row", "column", "both", or "none".
#' @param color_scheme Color palette name.
#' @param show_lipid_names Logical, show lipid name labels.
#' @param de_q_val_thresh Q-value threshold for significance.
#'
#' @return A ggplot2 or ComplexHeatmap object, or NULL.
#'
#' @importFrom ComplexHeatmap Heatmap
#' @importFrom circlize colorRamp2
#' @importFrom stats hclust dist cor
#' @importFrom dplyr filter arrange desc slice_head
#' @importFrom logger log_info log_error log_warn
#' @export
generateMetabDEHeatmap <- function(
    de_results_list
    , selected_contrast = NULL
    , selected_assay = NULL
    , top_n = 50
    , clustering_method = "ward.D2"
    , distance_method = "euclidean"
    , cluster_rows = TRUE
    , cluster_cols = TRUE
    , scale_data = "row"
    , color_scheme = "RdBu"
    , show_lipid_names = FALSE
    , de_q_val_thresh = 0.05
) {
    logger::log_info("--- Entering generateMetabDEHeatmap ---")
    logger::log_info(sprintf("   selected_contrast = %s, top_n = %d", selected_contrast, top_n))

    if (is.null(de_results_list) || is.null(de_results_list$de_lipids_long)) {
        logger::log_warn("   No DE results available")
        return(NULL)
    }

    if (is.null(selected_contrast)) {
        logger::log_warn("   No contrast selected")
        return(NULL)
    }

    # Get data
    de_lipids_long <- de_results_list$de_lipids_long

    # Extract comparison name
    comparison_to_search <- stringr::str_extract(selected_contrast, "^[^=]+")
    if (is.na(comparison_to_search)) {
        comparison_to_search <- selected_contrast
    }

    # Filter for significant results in selected contrast
    contrast_data <- de_lipids_long |>
        dplyr::filter(comparison == comparison_to_search | friendly_name == comparison_to_search) |>
        dplyr::filter(fdr_qvalue < de_q_val_thresh)

    # Optionally filter by assay
    if (!is.null(selected_assay) && selected_assay != "Combined") {
        contrast_data <- contrast_data |>
            dplyr::filter(assay == selected_assay)
    }

    if (nrow(contrast_data) == 0) {
        logger::log_warn("   No significant lipids found")
        return(NULL)
    }

    # Get top N by absolute logFC
    top_lipids <- contrast_data |>
        dplyr::arrange(dplyr::desc(abs(logFC))) |>
        dplyr::slice_head(n = top_n)

    logger::log_info(sprintf("   Selected %d top lipids", nrow(top_lipids)))

    # Get expression matrix
    theObject <- de_results_list$theObject

    # Determine which assay(s) to use
    if (!is.null(selected_assay) && selected_assay != "Combined") {
        assays_to_use <- selected_assay
    } else {
        assays_to_use <- unique(top_lipids$assay)
    }

    # Build expression matrix from relevant assays
    expr_list <- lapply(assays_to_use, function(assay_name) {
        assay_data <- theObject@lipid_data[[assay_name]]
        if (is.null(assay_data)) return(NULL)

        # Get lipid IDs for this assay
        assay_lipid_ids <- top_lipids |>
            dplyr::filter(assay == assay_name) |>
            dplyr::pull(lipid_id)

        if (length(assay_lipid_ids) == 0) return(NULL)

        # Filter to selected lipids
        rows_to_keep <- assay_data[[theObject@lipid_id_column]] %in% assay_lipid_ids
        assay_subset <- assay_data[rows_to_keep, , drop = FALSE]

        # Get sample columns
        sample_cols <- intersect(colnames(assay_subset), theObject@design_matrix[[theObject@sample_id]])

        # Build matrix
        mat <- as.matrix(assay_subset[, sample_cols, drop = FALSE])
        rownames(mat) <- assay_subset[[theObject@lipid_id_column]]

        return(mat)
    })

    # Combine matrices (if multiple assays)
    expr_matrix <- do.call(rbind, expr_list[!sapply(expr_list, is.null)])

    if (is.null(expr_matrix) || nrow(expr_matrix) == 0) {
        logger::log_warn("   Could not build expression matrix")
        return(NULL)
    }

    logger::log_info(sprintf("   Expression matrix: %d x %d", nrow(expr_matrix), ncol(expr_matrix)))

    # Apply scaling
    if (scale_data == "row") {
        expr_matrix <- t(scale(t(expr_matrix)))
    } else if (scale_data == "column") {
        expr_matrix <- scale(expr_matrix)
    } else if (scale_data == "both") {
        expr_matrix <- t(scale(t(expr_matrix)))
        expr_matrix <- scale(expr_matrix)
    }

    # Handle NA/Inf from scaling
    expr_matrix[is.na(expr_matrix)] <- 0
    expr_matrix[is.infinite(expr_matrix)] <- 0

    # Build row labels (lipid names if requested)
    if (show_lipid_names) {
        # Map IDs to names
        id_to_name <- stats::setNames(top_lipids$lipid_name, top_lipids$lipid_id)
        row_labels <- id_to_name[rownames(expr_matrix)]
        row_labels[is.na(row_labels)] <- rownames(expr_matrix)[is.na(row_labels)]
    } else {
        row_labels <- rownames(expr_matrix)
    }

    # Calculate clustering
    row_clust <- NULL
    col_clust <- NULL

    if (cluster_rows && nrow(expr_matrix) > 1) {
        if (distance_method %in% c("pearson", "spearman")) {
            row_dist <- stats::as.dist(1 - stats::cor(t(expr_matrix), method = distance_method, use = "pairwise.complete.obs"))
        } else {
            row_dist <- stats::dist(expr_matrix, method = distance_method)
        }
        row_clust <- stats::hclust(row_dist, method = clustering_method)
    }

    if (cluster_cols && ncol(expr_matrix) > 1) {
        if (distance_method %in% c("pearson", "spearman")) {
            col_dist <- stats::as.dist(1 - stats::cor(expr_matrix, method = distance_method, use = "pairwise.complete.obs"))
        } else {
            col_dist <- stats::dist(t(expr_matrix), method = distance_method)
        }
        col_clust <- stats::hclust(col_dist, method = clustering_method)
    }

    # Color palette
    color_fn <- switch(color_scheme
        , "RdBu" = circlize::colorRamp2(c(-2, 0, 2), c("blue", "white", "red"))
        , "RdYlBu" = circlize::colorRamp2(c(-2, 0, 2), c("blue", "yellow", "red"))
        , "coolwarm" = circlize::colorRamp2(c(-2, 0, 2), c("#3B4CC0", "white", "#B40426"))
        , "viridis" = circlize::colorRamp2(seq(-2, 2, length.out = 256), viridisLite::viridis(256))
        , "plasma" = circlize::colorRamp2(seq(-2, 2, length.out = 256), viridisLite::plasma(256))
        , "inferno" = circlize::colorRamp2(seq(-2, 2, length.out = 256), viridisLite::inferno(256))
        , circlize::colorRamp2(c(-2, 0, 2), c("blue", "white", "red"))
    )

    # Get column annotations (groups)
    dm <- theObject@design_matrix
    rownames(dm) <- dm[[theObject@sample_id]]
    col_groups <- dm[colnames(expr_matrix), theObject@group_id]

    # Create heatmap
    tryCatch({
        hm <- ComplexHeatmap::Heatmap(
            expr_matrix
            , name = "Z-score"
            , col = color_fn
            , cluster_rows = if (is.null(row_clust)) cluster_rows else row_clust
            , cluster_columns = if (is.null(col_clust)) cluster_cols else col_clust
            , show_row_names = show_lipid_names
            , row_labels = row_labels
            , show_column_names = TRUE
            , column_title = paste("Top", nrow(expr_matrix), "DE Lipids:", selected_contrast)
            , row_title = "Lipids"
            , top_annotation = ComplexHeatmap::HeatmapAnnotation(
                Group = col_groups
                , col = list(Group = c("Control" = "#4DBBD5", "Treatment" = "#E64B35"))
            )
        )

        logger::log_info("--- Exiting generateMetabDEHeatmap (success) ---")
        return(hm)

    }, error = function(e) {
        logger::log_error(sprintf("   Heatmap generation failed: %s", e$message))
        return(NULL)
    })
}


# ----------------------------------------------------------------------------
# generateMetabVolcanoStatic
# ----------------------------------------------------------------------------
#' Generate static ggplot2 volcano plot for lipidomics DE results
#'
#' @description Creates a static volcano plot using ggplot2 as a fallback
#'   when Glimma is not available or for export purposes.
#'
#' @param de_results_list Results list from `runLipidsDE()`.
#' @param selected_contrast Contrast to display.
#' @param selected_assay Optional assay filter (NULL for combined/faceted).
#' @param de_q_val_thresh Q-value threshold for significance.
#' @param lfc_threshold Log fold-change threshold for significance lines.
#' @param show_labels Logical, label top significant lipids.
#' @param n_labels Number of top lipids to label.
#'
#' @return A ggplot2 object.
#'
#' @importFrom ggplot2 ggplot aes geom_point geom_hline geom_vline labs theme_minimal scale_color_manual facet_wrap
#' @importFrom ggrepel geom_text_repel
#' @importFrom dplyr filter arrange slice_head
#' @export
generateMetabVolcanoStatic <- function(
    de_results_list
    , selected_contrast = NULL
    , selected_assay = NULL
    , de_q_val_thresh = 0.05
    , lfc_threshold = 1
    , show_labels = TRUE
    , n_labels = 10
) {
    if (is.null(de_results_list) || is.null(de_results_list$de_lipids_long)) {
        return(NULL)
    }

    if (is.null(selected_contrast)) {
        return(NULL)
    }

    # Get data
    de_lipids_long <- de_results_list$de_lipids_long

    # Extract comparison name
    comparison_to_search <- stringr::str_extract(selected_contrast, "^[^=]+")
    if (is.na(comparison_to_search)) {
        comparison_to_search <- selected_contrast
    }

    # Filter
    plot_data <- de_lipids_long |>
        dplyr::filter(comparison == comparison_to_search | friendly_name == comparison_to_search)

    if (!is.null(selected_assay) && selected_assay != "Combined") {
        plot_data <- plot_data |>
            dplyr::filter(assay == selected_assay)
    }

    if (nrow(plot_data) == 0) {
        return(NULL)
    }

    # Add plot columns
    plot_data <- plot_data |>
        dplyr::mutate(
            neg_log10_q = -log10(fdr_qvalue)
            , display_name = ifelse(!is.na(lipid_name) & lipid_name != ""
                , lipid_name, lipid_id)
        )

    # Get top lipids for labeling
    if (show_labels) {
        top_to_label <- plot_data |>
            dplyr::filter(significant != "NS") |>
            dplyr::arrange(fdr_qvalue) |>
            dplyr::slice_head(n = n_labels)
    }

    # Build plot
    p <- ggplot2::ggplot(plot_data, ggplot2::aes(
        x = logFC
        , y = neg_log10_q
        , color = significant
    )) +
        ggplot2::geom_point(alpha = 0.6, size = 2) +
        ggplot2::geom_hline(
            yintercept = -log10(de_q_val_thresh)
            , linetype = "dashed"
            , color = "gray50"
        ) +
        ggplot2::geom_vline(
            xintercept = c(-lfc_threshold, lfc_threshold)
            , linetype = "dashed"
            , color = "gray50"
        ) +
        ggplot2::scale_color_manual(
            values = c("Up" = "#E64B35", "Down" = "#4DBBD5", "NS" = "gray70")
            , name = "Significance"
        ) +
        ggplot2::labs(
            title = paste("Volcano Plot:", selected_contrast)
            , x = "Log2 Fold Change"
            , y = "-log10(Q-value)"
        ) +
        ggplot2::theme_minimal() +
        ggplot2::theme(legend.position = "bottom")

    # Add labels
    if (show_labels && nrow(top_to_label) > 0) {
        p <- p + ggrepel::geom_text_repel(
            data = top_to_label
            , ggplot2::aes(label = display_name)
            , size = 3
            , max.overlaps = 20
        )
    }

    # Facet if multiple assays
    if (is.null(selected_assay) || selected_assay == "Combined") {
        if (length(unique(plot_data$assay)) > 1) {
            p <- p + ggplot2::facet_wrap(~assay, scales = "free")
        }
    }

    return(p)
}


# ----------------------------------------------------------------------------
# outputMetabDeResultsAllContrasts
# ----------------------------------------------------------------------------
#' Output all lipidomics DE results to disk
#'
#' @description Writes DE results tables, volcano plots, and heatmaps to disk
#'   for all contrasts in a lipidomics DE analysis. Outputs are split by
#'   assay mode (posmode/negmode) and contrast, matching the proteomics workflow.
#'
#' @details Output filenames follow the pattern:
#'   - `de_posmode_lipids_{contrast}_long_annot.xlsx`
#'   - `de_negmode_lipids_{contrast}_long_annot.xlsx`
#'
#' @param de_results_list Results list from `runLipidsDE()`.
#' @param de_output_dir Directory for DE result tables.
#' @param publication_graphs_dir Directory for publication-quality figures.
#' @param de_q_val_thresh Q-value threshold for significance (default 0.05).
#' @param lfc_threshold Log fold-change threshold for volcano plot lines.
#' @param heatmap_top_n Number of top lipids for heatmap (default 50).
#' @param heatmap_clustering Clustering option: "both", "row", "column", "none".
#' @param heatmap_color_scheme Color scheme for heatmap.
#'
#' @return TRUE if successful, FALSE otherwise.
#'
#' @importFrom vroom vroom_write
#' @importFrom writexl write_xlsx
#' @importFrom ggplot2 ggsave
#' @importFrom grDevices pdf dev.off png
#' @importFrom dplyr filter group_by summarise n mutate select
#' @importFrom purrr walk map
#' @importFrom logger log_info log_error log_warn
#' @export
outputMetabDeResultsAllContrasts <- function(
    de_results_list
    , de_output_dir
    , publication_graphs_dir
    , de_q_val_thresh = 0.05
    , lfc_threshold = 1
    , heatmap_top_n = 50
    , heatmap_clustering = "both"
    , heatmap_color_scheme = "RdBu"
) {
    logger::log_info("--- Entering outputMetabDeResultsAllContrasts ---")
    logger::log_info(sprintf("   de_output_dir = %s", de_output_dir))
    logger::log_info(sprintf("   publication_graphs_dir = %s", publication_graphs_dir))

    # Validate inputs
    if (is.null(de_results_list) || is.null(de_results_list$de_lipids_long)) {
        logger::log_error("   No DE results available")
        return(FALSE)
    }

    # Normalize paths
    if (!is.null(de_output_dir)) {
        de_output_dir <- gsub("//+", "/", de_output_dir)
        de_output_dir <- normalizePath(de_output_dir, winslash = "/", mustWork = FALSE)
    }

    if (!is.null(publication_graphs_dir)) {
        publication_graphs_dir <- gsub("//+", "/", publication_graphs_dir)
        publication_graphs_dir <- normalizePath(publication_graphs_dir, winslash = "/", mustWork = FALSE)
    }

    # Create output directories
    if (!dir.exists(de_output_dir)) {
        dir.create(de_output_dir, recursive = TRUE, showWarnings = FALSE)
        logger::log_info(sprintf("   Created de_output_dir: %s", de_output_dir))
    }

    volcano_dir <- file.path(publication_graphs_dir, "Volcano_Plots")
    if (!dir.exists(volcano_dir)) {
        dir.create(volcano_dir, recursive = TRUE, showWarnings = FALSE)
        logger::log_info(sprintf("   Created Volcano_Plots directory: %s", volcano_dir))
    }

    heatmap_dir <- file.path(publication_graphs_dir, "Heatmaps")
    if (!dir.exists(heatmap_dir)) {
        dir.create(heatmap_dir, recursive = TRUE, showWarnings = FALSE)
        logger::log_info(sprintf("   Created Heatmaps directory: %s", heatmap_dir))
    }

    numsigde_dir <- file.path(publication_graphs_dir, "NumSigDeMolecules")
    if (!dir.exists(numsigde_dir)) {
        dir.create(numsigde_dir, recursive = TRUE, showWarnings = FALSE)
        logger::log_info(sprintf("   Created NumSigDeMolecules directory: %s", numsigde_dir))
    }

    # Get DE results with sample intensity columns (already included from runLipidsDE)
    de_lipids_long <- de_results_list$de_lipids_long

    # Get unique assays and contrasts
    assays <- unique(de_lipids_long$assay)
    contrasts <- unique(de_lipids_long$comparison)

    logger::log_info(sprintf("   Found %d assays: %s", length(assays), paste(assays, collapse = ", ")))
    logger::log_info(sprintf("   Found %d contrasts: %s", length(contrasts), paste(contrasts, collapse = ", ")))

    # Map assay names to mode prefixes
    # LCMS_Pos -> posmode, LCMS_Neg -> negmode
    get_mode_prefix <- function(assay_name) {
        if (grepl("pos", assay_name, ignore.case = TRUE)) {
            return("posmode")
        } else if (grepl("neg", assay_name, ignore.case = TRUE)) {
            return("negmode")
        } else {
            # Fallback: use sanitized assay name
            return(gsub("[^A-Za-z0-9]", "", tolower(assay_name)))
        }
    }

    # Storage for combined PDFs
    all_volcano_plots <- list()
    all_heatmap_plots <- list()
    all_numsig_tables <- list()

    # Process each ASSAY x CONTRAST combination
    for (assay_name in assays) {
        mode_prefix <- get_mode_prefix(assay_name)
        logger::log_info(sprintf("   Processing assay: %s (mode: %s)", assay_name, mode_prefix))

        for (contrast_name in contrasts) {
            logger::log_info(sprintf("      Processing contrast: %s", contrast_name))

            # Filter data for this assay + contrast
            assay_contrast_data <- de_lipids_long |>
                dplyr::filter(assay == assay_name & comparison == contrast_name)

            if (nrow(assay_contrast_data) == 0) {
                logger::log_warn(sprintf("      No data for %s / %s, skipping", assay_name, contrast_name))
                next
            }

            # =====================================================================
            # 1. Write DE results tables (TSV and Excel)
            # Filename format: de_{mode}_lipids_{contrast}_long_annot.xlsx
            # =====================================================================
            tryCatch({
                # Reorder columns: ID, name first, then stats, then intensity columns
                priority_cols <- c("lipid_id", "lipid_name"
                    , "logFC", "raw_pvalue", "fdr_qvalue", "significant"
                    , "comparison", "friendly_name", "numerator", "denominator")
                priority_cols <- intersect(priority_cols, colnames(assay_contrast_data))

                # Get sample intensity columns
                intensity_cols <- grep("^intensity\\.", colnames(assay_contrast_data), value = TRUE)

                # Get any remaining columns (exclude assay since we're splitting by it)
                other_cols <- setdiff(
                    colnames(assay_contrast_data)
                    , c(priority_cols, intensity_cols, "assay")
                )

                # Final column order
                final_col_order <- c(priority_cols, other_cols, sort(intensity_cols))
                output_data <- assay_contrast_data[, final_col_order, drop = FALSE]

                # Build filename: de_{mode}_lipids_{contrast}_long_annot
                file_base <- paste0("de_", mode_prefix, "_lipids_", contrast_name, "_long_annot")

                # Write TSV
                tsv_path <- file.path(de_output_dir, paste0(file_base, ".tsv"))
                vroom::vroom_write(output_data, tsv_path)
                logger::log_info(sprintf("      Wrote TSV: %s (%d rows, %d cols)"
                    , basename(tsv_path), nrow(output_data), ncol(output_data)))

                # Write Excel
                xlsx_path <- file.path(de_output_dir, paste0(file_base, ".xlsx"))
                writexl::write_xlsx(output_data, xlsx_path)
                logger::log_info(sprintf("      Wrote Excel: %s", basename(xlsx_path)))

            }, error = function(e) {
                logger::log_error(sprintf("      Error writing tables for %s / %s: %s"
                    , assay_name, contrast_name, e$message))
            })

            # =================================================================
            # 2. Generate and save volcano plots (per assay, per contrast)
            # =================================================================
            tryCatch({
                # Build a filtered results list for just this assay
                assay_filtered_results <- de_results_list
                assay_filtered_results$de_lipids_long <- assay_contrast_data

                volcano_plot <- generateMetabVolcanoStatic(
                    de_results_list = assay_filtered_results
                    , selected_contrast = contrast_name
                    , selected_assay = assay_name
                    , de_q_val_thresh = de_q_val_thresh
                    , lfc_threshold = lfc_threshold
                    , show_labels = TRUE
                    , n_labels = 15
                )

                if (!is.null(volcano_plot)) {
                    # Unique key for combined PDF
                    plot_key <- paste0(mode_prefix, "_", contrast_name)
                    all_volcano_plots[[plot_key]] <- volcano_plot

                    # Filename: {mode}_{contrast}_volcano
                    volcano_base <- paste0(mode_prefix, "_", contrast_name)

                    # Save PNG
                    volcano_png <- file.path(volcano_dir, paste0(volcano_base, ".png"))
                    ggplot2::ggsave(volcano_png, volcano_plot, width = 8, height = 7, dpi = 300)
                    logger::log_info(sprintf("      Saved volcano PNG: %s", basename(volcano_png)))

                    # Save PDF
                    volcano_pdf <- file.path(volcano_dir, paste0(volcano_base, ".pdf"))
                    ggplot2::ggsave(volcano_pdf, volcano_plot, width = 8, height = 7)
                    logger::log_info(sprintf("      Saved volcano PDF: %s", basename(volcano_pdf)))
                }

            }, error = function(e) {
                logger::log_error(sprintf("      Error generating volcano for %s / %s: %s"
                    , assay_name, contrast_name, e$message))
            })

            # =================================================================
            # 3. Generate and save heatmaps (per assay, per contrast)
            # =================================================================
            tryCatch({
                # Build a filtered results list for just this assay
                assay_filtered_results <- de_results_list
                assay_filtered_results$de_lipids_long <- assay_contrast_data

                heatmap_obj <- generateMetabDEHeatmap(
                    de_results_list = assay_filtered_results
                    , selected_contrast = contrast_name
                    , selected_assay = assay_name
                    , top_n = heatmap_top_n
                    , clustering_method = "ward.D2"
                    , distance_method = "euclidean"
                    , cluster_rows = heatmap_clustering %in% c("both", "row")
                    , cluster_cols = heatmap_clustering %in% c("both", "column")
                    , scale_data = "row"
                    , color_scheme = heatmap_color_scheme
                    , show_lipid_names = FALSE
                    , de_q_val_thresh = de_q_val_thresh
                )

                if (!is.null(heatmap_obj)) {
                    # Unique key for combined PDF
                    plot_key <- paste0(mode_prefix, "_", contrast_name)
                    all_heatmap_plots[[plot_key]] <- heatmap_obj

                    # Filename: {mode}_{contrast}_heatmap
                    heatmap_base <- paste0(mode_prefix, "_", contrast_name)

                    # Save PNG
                    heatmap_png <- file.path(heatmap_dir, paste0(heatmap_base, "_heatmap.png"))
                    grDevices::png(heatmap_png, width = 10, height = 8, units = "in", res = 300)
                    ComplexHeatmap::draw(heatmap_obj)
                    grDevices::dev.off()
                    logger::log_info(sprintf("      Saved heatmap PNG: %s", basename(heatmap_png)))

                    # Save PDF
                    heatmap_pdf <- file.path(heatmap_dir, paste0(heatmap_base, "_heatmap.pdf"))
                    grDevices::pdf(heatmap_pdf, width = 10, height = 8)
                    ComplexHeatmap::draw(heatmap_obj)
                    grDevices::dev.off()
                    logger::log_info(sprintf("      Saved heatmap PDF: %s", basename(heatmap_pdf)))
                } else {
                    logger::log_warn(sprintf("      No significant lipids for heatmap: %s / %s"
                        , assay_name, contrast_name))
                }

            }, error = function(e) {
                logger::log_error(sprintf("      Error generating heatmap for %s / %s: %s"
                    , assay_name, contrast_name, e$message))
            })

            # =================================================================
            # 4. Calculate NumSigDeMolecules for this assay/contrast
            # =================================================================
            tryCatch({
                sig_summary <- assay_contrast_data |>
                    dplyr::summarise(
                        total = dplyr::n()
                        , significant = sum(significant != "NS", na.rm = TRUE)
                        , up_regulated = sum(significant == "Up", na.rm = TRUE)
                        , down_regulated = sum(significant == "Down", na.rm = TRUE)
                        , .groups = "drop"
                    ) |>
                    dplyr::mutate(
                        assay = assay_name
                        , mode = mode_prefix
                        , contrast = contrast_name
                        , q_threshold = de_q_val_thresh
                    )

                table_key <- paste0(mode_prefix, "_", contrast_name)
                all_numsig_tables[[table_key]] <- sig_summary

            }, error = function(e) {
                logger::log_error(sprintf("      Error calculating NumSigDE for %s / %s: %s"
                    , assay_name, contrast_name, e$message))
            })

        } # End contrast loop
    } # End assay loop

    # =========================================================================
    # 5. Create combined volcano plots PDF
    # =========================================================================
    if (length(all_volcano_plots) > 0) {
        tryCatch({
            combined_volcano_pdf <- file.path(volcano_dir, "all_volcano_plots_combined.pdf")
            grDevices::pdf(combined_volcano_pdf, width = 8, height = 7, onefile = TRUE)
            purrr::walk(all_volcano_plots, print)
            grDevices::dev.off()
            logger::log_info(sprintf("   Created combined volcano PDF: %d plots"
                , length(all_volcano_plots)))
        }, error = function(e) {
            logger::log_error(sprintf("   Error creating combined volcano PDF: %s", e$message))
        })
    }

    # =========================================================================
    # 6. Create combined heatmaps PDF
    # =========================================================================
    if (length(all_heatmap_plots) > 0) {
        tryCatch({
            combined_heatmap_pdf <- file.path(heatmap_dir, "all_heatmaps_combined.pdf")
            grDevices::pdf(combined_heatmap_pdf, width = 10, height = 8, onefile = TRUE)
            purrr::walk(all_heatmap_plots, function(hm) {
                ComplexHeatmap::draw(hm)
            })
            grDevices::dev.off()
            logger::log_info(sprintf("   Created combined heatmap PDF: %d plots"
                , length(all_heatmap_plots)))
        }, error = function(e) {
            logger::log_error(sprintf("   Error creating combined heatmap PDF: %s", e$message))
        })
    }

    # =========================================================================
    # 7. Write NumSigDeMolecules summary table
    # =========================================================================
    if (length(all_numsig_tables) > 0) {
        tryCatch({
            combined_numsig <- dplyr::bind_rows(all_numsig_tables)

            # Write TSV
            numsig_tsv <- file.path(numsigde_dir, "lipids_num_sig_de_molecules.tab")
            vroom::vroom_write(combined_numsig, numsig_tsv)
            logger::log_info(sprintf("   Wrote NumSigDE table: %s", basename(numsig_tsv)))

            # Write Excel
            numsig_xlsx <- file.path(numsigde_dir, "lipids_num_sig_de_molecules.xlsx")
            writexl::write_xlsx(combined_numsig, numsig_xlsx)
            logger::log_info(sprintf("   Wrote NumSigDE Excel: %s", basename(numsig_xlsx)))

            # Create bar plot
            numsig_plot <- ggplot2::ggplot(combined_numsig
                , ggplot2::aes(x = contrast, y = significant, fill = mode)) +
                ggplot2::geom_bar(stat = "identity", position = "dodge") +
                ggplot2::labs(
                    title = "Number of Significant DE Lipids"
                    , x = "Contrast"
                    , y = "Number of Significant Lipids"
                    , fill = "Assay"
                ) +
                ggplot2::theme_minimal() +
                ggplot2::theme(
                    axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)
                )

            # Save bar plot
            numsig_png <- file.path(numsigde_dir, "lipids_num_sig_barplot.png")
            ggplot2::ggsave(numsig_png, numsig_plot, width = 10, height = 6, dpi = 300)
            logger::log_info(sprintf("   Saved NumSigDE barplot: %s", basename(numsig_png)))

        }, error = function(e) {
            logger::log_error(sprintf("   Error writing NumSigDE summary: %s", e$message))
        })
    }

    logger::log_info("--- Exiting outputMetabDeResultsAllContrasts ---")
    return(TRUE)
}

