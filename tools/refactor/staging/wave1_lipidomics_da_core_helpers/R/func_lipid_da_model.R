# ----------------------------------------------------------------------------
# runTestsContrastsLipidDA
# ----------------------------------------------------------------------------
#' Run limma differential abundance analysis for a single assay
#'
#' @description Performs limma-based differential abundance analysis on a
#'   single lipidomics assay matrix.
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
runTestsContrastsLipidDA <- function(
  data,
  contrast_strings,
  design_matrix,
  formula_string,
  sample_id_col = "Run",
  treat_lfc_cutoff = NA,
  eBayes_trend = TRUE,
  eBayes_robust = TRUE
) {
    logger::log_info("--- Entering runTestsContrastsLipidDA ---")
    logger::log_info(sprintf(
        "   data dims = %d x %d, %d contrasts",
        nrow(data), ncol(data), length(contrast_strings)
    ))
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

    # [D66:START] -------------------------
    d66_log <- function(...) message(sprintf("[D66] %s", paste0(...)))
    d66_log("  runTestsContrastsLipidDA - Model matrix created:")
    d66_log("    design_m levels (colnames) = ", paste(colnames(design_m), collapse = ", "))
    d66_log("    contrast_strings = ", paste(contrast_strings, collapse = ", "))
    # [D66:END] ---------------------------

    # Subset data to match design matrix samples
    common_samples <- intersect(colnames(data), rownames(design_m))
    if (length(common_samples) == 0) {
        stop("No common samples between data matrix and design matrix")
    }
    data_subset <- data[, common_samples, drop = FALSE]
    design_m <- design_m[common_samples, , drop = FALSE]
    logger::log_info(sprintf("   data_subset dims = %d x %d", nrow(data_subset), ncol(data_subset)))

    # [D66:START] -------------------------
    d66_log("  About to call makeContrasts:")
    d66_log("    contrasts = ", paste(contrast_strings, collapse = ", "))
    d66_log("    levels = ", paste(colnames(design_m), collapse = ", "))
    # [D66:END] ---------------------------

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
        terms <- terms[nzchar(terms)] # Remove empty strings

        missing_terms <- setdiff(terms, available_levels)
        if (length(missing_terms) > 0) {
            error_msg <- sprintf(
                paste0(
                    "Contrast '%s' references undefined levels: [%s].\n",
                    "Available levels from model matrix: [%s].\n\n",
                    "This usually means the contrast was saved with a different formula than currently used.\n",
                    "If using formula '~ 0 + group', contrasts must be in format 'groupX-groupY'.\n",
                    "Check that the design builder formula matches the DA module formula."
                ),
                cs,
                paste(missing_terms, collapse = ", "),
                paste(available_levels, collapse = ", ")
            )
            d66_log("  VALIDATION FAILED: ", error_msg)
            stop(error_msg)
        }
    }
    d66_log("  VALIDATION PASSED: All contrast terms exist in design matrix levels")

    # Create contrast matrix
    logger::log_info("   Creating contrast matrix...")
    contr.matrix <- limma::makeContrasts(
        contrasts = contrast_strings,
        levels = colnames(design_m)
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

        tryCatch(
            {
                da_tbl <- extract_fn(t.fit, contrast)
                logger::log_info(sprintf("      Extracted %d rows", nrow(da_tbl)))

                # Safe qvalue computation
                valid_p_idx <- which(!is.na(da_tbl$P.Value) & is.finite(da_tbl$P.Value))
                q_values_all <- rep(NA_real_, nrow(da_tbl))
                fdr_values_all <- rep(NA_real_, nrow(da_tbl))

                if (length(valid_p_idx) > 0) {
                    valid_p_values <- da_tbl$P.Value[valid_p_idx]

                    # Edge case checks
                    all_zeros <- all(valid_p_values == 0)
                    all_ones <- all(valid_p_values == 1)
                    too_few <- length(valid_p_values) < 3

                    use_qvalue <- !all_zeros && !all_ones && !too_few

                    if (use_qvalue) {
                        tryCatch(
                            {
                                q_values_all[valid_p_idx] <- qvalue::qvalue(valid_p_values)$qvalues
                                logger::log_info("      qvalue() computation successful")
                            },
                            error = function(e) {
                                qvalue_failed <<- TRUE
                                logger::log_warn(sprintf("      qvalue() failed: %s, using BH", e$message))
                                q_values_all[valid_p_idx] <<- stats::p.adjust(valid_p_values, method = "BH")
                            }
                        )
                    } else {
                        qvalue_failed <- TRUE
                        q_values_all[valid_p_idx] <- stats::p.adjust(valid_p_values, method = "BH")
                        logger::log_info("      Using p.adjust() due to edge case")
                    }

                    # Always compute BH-adjusted FDR
                    fdr_values_all[valid_p_idx] <- stats::p.adjust(valid_p_values, method = "BH")
                }

                # Add columns to results
                da_tbl$fdr_qvalue <- q_values_all
                da_tbl$fdr_value_bh <- fdr_values_all
                da_tbl$raw_pvalue <- da_tbl$P.Value

                if (qvalue_failed) {
                    qvalue_failures[[contrast]] <<- TRUE
                }

                return(da_tbl)
            },
            error = function(e) {
                logger::log_error(sprintf("      ERROR in contrast %s: %s", contrast, e$message))
                stop(e)
            }
        )
    })

    names(result_tables) <- contrast_strings
    logger::log_info("--- Exiting runTestsContrastsLipidDA ---")

    return(list(
        results = result_tables,
        fit.eb = t.fit,
        qvalue_warnings = qvalue_failures
    ))
}

# ----------------------------------------------------------------------------
# runLipidsDA
# ----------------------------------------------------------------------------
#' Run differential abundance analysis on all assays in a LipidomicsAssayData object
#'
#' @description Main entry point for lipidomics differential abundance analysis.
#'   Loops over all assays in the LipidomicsAssayData object, runs limma DA on each,
#'   and aggregates results with an assay identifier column.
#'
#' @param theObject LipidomicsAssayData S4 object containing lipid data.
#' @param contrasts_tbl Data frame with contrast definitions. Must have a column
#'   named "contrasts" or "contrast_string" containing the contrast formulas.
#' @param formula_string Model formula (default "~ 0 + group").
#' @param da_q_val_thresh Q-value threshold for significance (default 0.05).
#' @param treat_lfc_cutoff Log fold-change cutoff for TREAT (default 0, standard analysis).
#' @param eBayes_trend Logical, use trend in eBayes (default TRUE).
#' @param eBayes_robust Logical, use robust eBayes (default TRUE).
#'
#' @return A list containing:
#'   - theObject: The input S4 object (unchanged)
#'   - contrasts_results: Per-assay limma fit objects
#'   - da_lipids_long: Combined long-format results with assay column
#'   - per_assay_results: List of per-assay result tables
#'   - significant_counts: Summary of significant lipids per assay
#'   - qvalue_warnings: Any qvalue computation warnings
#'
#' @importFrom tibble rownames_to_column
#' @importFrom dplyr mutate bind_rows select left_join
#' @importFrom purrr map map2 set_names
#' @importFrom logger log_info log_error log_warn
#' @export
runLipidsDA <- function(
  theObject,
  contrasts_tbl,
  formula_string = "~ 0 + group",
  da_q_val_thresh = 0.05,
  treat_lfc_cutoff = 0,
  eBayes_trend = TRUE,
  eBayes_robust = TRUE
) {
    # [D66:START] -------------------------
    d66_log <- function(...) message(sprintf("[D66] %s", paste0(...)))
    d66_log("=== ENTER runLipidsDA ===")
    d66_log("  formula_string = ", formula_string)
    d66_log("  contrasts_tbl columns = ", paste(colnames(contrasts_tbl), collapse = ", "))
    d66_log("  contrasts_tbl$contrasts = ", paste(contrasts_tbl$contrasts, collapse = ", "))
    # [D66:END] ---------------------------

    logger::log_info("=== Starting runLipidsDA ===")

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

    # [D66:START] -------------------------
    d66_log("  S4 slots extracted:")
    d66_log("    sample_id_col = ", sample_id_col)
    d66_log("    group_col = ", group_col)
    d66_log("    design_matrix columns = ", paste(colnames(design_matrix), collapse = ", "))
    d66_log("    design_matrix nrow = ", nrow(design_matrix))
    d66_log("    assay_list names = ", paste(names(assay_list), collapse = ", "))
    # [D66:END] ---------------------------

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

    # [D66:START] -------------------------
    d66_log("  contrast_strings = ", paste(contrast_strings, collapse = ", "))
    d66_log("  friendly_names = ", paste(friendly_names, collapse = ", "))
    # [D66:END] ---------------------------

    assay_names <- names(assay_list)
    logger::log_info(sprintf(
        "   Processing %d assays: %s",
        length(assay_names), paste(assay_names, collapse = ", ")
    ))
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
            lipid_id = assay_data[[lipid_id_col]],
            lipid_name = assay_data[[annotation_col]],
            stringsAsFactors = FALSE
        )

        logger::log_info(sprintf(
            "   Expression matrix: %d lipids x %d samples",
            nrow(expr_matrix), ncol(expr_matrix)
        ))

        # Run limma DE
        tryCatch(
            {
                limma_results <- runTestsContrastsLipidDA(
                    data = expr_matrix,
                    contrast_strings = contrast_strings,
                    design_matrix = design_matrix,
                    formula_string = formula_string,
                    sample_id_col = sample_id_col,
                    treat_lfc_cutoff = if (treat_lfc_cutoff > 0) treat_lfc_cutoff else NA,
                    eBayes_trend = eBayes_trend,
                    eBayes_robust = eBayes_robust
                )

                # Store fit object for Glimma
                contrasts_results[[assay_name]] <- limma_results

                # Collect qvalue warnings
                if (length(limma_results$qvalue_warnings) > 0) {
                    all_qvalue_warnings[[assay_name]] <- limma_results$qvalue_warnings
                }

                # Process results for each contrast
                assay_results_list <- purrr::map2(
                    limma_results$results,
                    names(limma_results$results),
                    function(da_tbl, contrast_name) {
                        # Add lipid_id from rownames
                        da_tbl <- tibble::rownames_to_column(da_tbl, "lipid_id")

                        # Add metadata columns
                        da_tbl$assay <- assay_name
                        da_tbl$comparison <- contrast_name

                        # Add friendly name
                        idx <- match(contrast_name, contrast_strings)
                        da_tbl$friendly_name <- if (!is.na(idx)) friendly_names[idx] else contrast_name

                        # Join annotation info
                        da_tbl <- dplyr::left_join(
                            da_tbl,
                            annotation_info,
                            by = "lipid_id"
                        )

                        # Classify significance
                        da_tbl$significant <- ifelse(
                            da_tbl$fdr_qvalue < da_q_val_thresh & abs(da_tbl$logFC) >= treat_lfc_cutoff,
                            ifelse(da_tbl$logFC > 0, "Up", "Down"),
                            "NS"
                        )

                        # Add sample intensity columns using createLipidDaResultsLongFormat
                        cat(sprintf("[D66] Calling createLipidDaResultsLongFormat for contrast: %s\n", contrast_name))
                        da_tbl <- createLipidDaResultsLongFormat(
                            lfc_qval_tbl = da_tbl,
                            expr_matrix = expr_matrix,
                            design_matrix = design_matrix,
                            sample_id_col = sample_id_col,
                            group_id_col = group_col,
                            lipid_id_col = "lipid_id"
                        )

                        return(da_tbl)
                    }
                )

                # Combine all contrasts for this assay
                assay_combined <- dplyr::bind_rows(assay_results_list)
                per_assay_results[[assay_name]] <- assay_combined

                logger::log_info(sprintf(
                    "   Completed assay %s: %d total results",
                    assay_name, nrow(assay_combined)
                ))
            },
            error = function(e) {
                cat(sprintf("[ERROR] Processing assay %s: %s\n", assay_name, e$message))
                logger::log_error(sprintf("   ERROR processing assay %s: %s", assay_name, e$message))
                per_assay_results[[assay_name]] <<- NULL
            }
        )
    }

    # Combine all assay results
    da_lipids_long <- dplyr::bind_rows(per_assay_results)
    logger::log_info(sprintf("   Combined results: %d total rows", nrow(da_lipids_long)))

    # Calculate significant counts per assay
    significant_counts <- purrr::map(per_assay_results, function(df) {
        if (is.null(df) || nrow(df) == 0) {
            return(list(up = 0, down = 0, ns = 0))
        }
        list(
            up = sum(df$significant == "Up", na.rm = TRUE),
            down = sum(df$significant == "Down", na.rm = TRUE),
            ns = sum(df$significant == "NS", na.rm = TRUE)
        )
    }) |> purrr::set_names(names(per_assay_results))

    logger::log_info("=== Completed runLipidsDA ===")

    return(list(
        theObject = theObject,
        contrasts_results = contrasts_results,
        da_lipids_long = da_lipids_long,
        per_assay_results = per_assay_results,
        significant_counts = significant_counts,
        qvalue_warnings = all_qvalue_warnings,
        da_q_val_thresh = da_q_val_thresh,
        treat_lfc_cutoff = treat_lfc_cutoff
    ))
}

