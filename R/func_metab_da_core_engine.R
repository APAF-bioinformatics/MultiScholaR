# ----------------------------------------------------------------------------
# runTestsContrastsMetabDA
# ----------------------------------------------------------------------------
#' Run limma differential abundance analysis for a single assay
#'
#' @description Performs limma-based differential abundance analysis on a
#'   single metabolomics assay matrix. This is the core DA engine adapted from
#'   the proteomics `runTestsContrastsDA()` function.
#'
#' @param data Numeric matrix with metabolites as rows and samples as columns.
#'   Rownames should be metabolite IDs.
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
runTestsContrastsMetabDA <- function(
  data,
  contrast_strings,
  design_matrix,
  formula_string,
  sample_id_col = "Run",
  treat_lfc_cutoff = NA,
  eBayes_trend = TRUE,
  eBayes_robust = TRUE
) {
    logger::log_info("--- Entering runTestsContrastsMetabDA ---")
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
    d66_log("  runTestsContrastsMetabDA - Model matrix created:")
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
    logger::log_info("--- Exiting runTestsContrastsMetabDA ---")
    return(list(
        results = result_tables,
        fit.eb = t.fit,
        qvalue_warnings = qvalue_failures
    ))
}

