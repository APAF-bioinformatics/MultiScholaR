# KIT-076: the fit half of `runTestsContrastsMetabDA`, moved unchanged: the model matrix on every sample, the samples the data
# shares with it, and `lmFit`. It depends on no contrast, so one saved fit serves every contrast of an assay.
#' @param data,design_matrix,formula_string,sample_id_col As `runTestsContrastsMetabDA`.
#' @param group_id,technical_replicate_id The object's declared columns (KIT-128). Without `technical_replicate_id` the
#'   fit is the unblocked `lmFit`, as before. With it the proteomics rule (`resolveProtTechnicalReplicateBlock`) reads the
#'   design rows of the fitted samples: declared technical replicates are fitted with the biological unit as block
#'   (`duplicateCorrelation`, then `lmFit(block, correlation)`; limma User's Guide 11.3 and 9.7), any other design as before.
#' @return The uncontrasted `MArrayLM`; its `design` is the model matrix on the shared samples. With
#'   `technical_replicate_id` it also carries `block_reason`, `uses_duplicate_correlation` and `consensus_correlation`.
#' @noRd
fitMetabDaAssayModel <- function(data, design_matrix, formula_string, sample_id_col = "Run", group_id = NULL,
                                 technical_replicate_id = NULL) {
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


    # Subset data to match design matrix samples
    common_samples <- intersect(colnames(data), rownames(design_m))
    if (length(common_samples) == 0) {
        stop("No common samples between data matrix and design matrix")
    }
    data_subset <- data[, common_samples, drop = FALSE]
    design_m <- design_m[common_samples, , drop = FALSE]
    logger::log_info(sprintf("   data_subset dims = %d x %d", nrow(data_subset), ncol(data_subset)))

    logger::log_info("   Running lmFit...")
    if (is.null(technical_replicate_id)) return(limma::lmFit(data_subset, design = design_m))
    # KIT-128: the rule proteomics uses, on the design rows of the fitted samples.
    replicate <- resolveProtTechnicalReplicateBlock(as.data.frame(dm)[common_samples, , drop = FALSE], group_id, technical_replicate_id)
    refuseNoBiologicalReplication(replicate, as.data.frame(dm)[common_samples, group_id])
    correlation <- NA_real_
    if (isTRUE(replicate$technical)) {
        correlation <- limma::duplicateCorrelation(data_subset, design = design_m, block = replicate$block)$consensus.correlation
        if (length(correlation) != 1L || !is.finite(correlation)) {
            rlang::abort(sprintf(paste("the correlation between technical replicates of one biological unit could not be",
                "estimated (limma's duplicateCorrelation gave %s), so the model with technical replicates cannot be fitted"),
                paste(format(correlation), collapse = ", ")), class = "omics_da_block_correlation_failed")
        }
        logger::log_info(sprintf("   Technical replicates: %d blocks for %d samples, consensus correlation %.4f",
            replicate$blocks, length(common_samples), correlation))
        fit <- limma::lmFit(data_subset, design = design_m, block = replicate$block, correlation = correlation)
    } else {
        fit <- limma::lmFit(data_subset, design = design_m)
    }
    fit$block_reason <- replicate$reason
    fit$uses_duplicate_correlation <- isTRUE(replicate$technical)
    fit$consensus_correlation <- correlation
    fit
}

# KIT-076: the contrast half of `runTestsContrastsMetabDA`, moved unchanged: the terms checked against the model, then
# `contrasts.fit`, `eBayes`, `treat` when asked, the extraction and the q-values, for one or more contrasts.
#' @param fit The fit half's `MArrayLM`.
#' @param contrast_strings,treat_lfc_cutoff,eBayes_trend,eBayes_robust As `runTestsContrastsMetabDA`.
#' @param contrastsFitFn The contrast step: `limma::contrasts.fit` (today's), or the exact one of the factors-in-the-model
#'   design (`ruvExactContrastsFit`, decision M29).
#' @return `list(results, fit.eb, qvalue_warnings)` for those contrasts.
#' @noRd
applyMetabDaContrasts <- function(fit, contrast_strings, treat_lfc_cutoff = NA, eBayes_trend = TRUE, eBayes_robust = TRUE,
                                  contrastsFitFn = limma::contrasts.fit) {
    design_m <- fit$design


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
            stop(error_msg)
        }
    }

    # Create contrast matrix
    logger::log_info("   Creating contrast matrix...")
    contr.matrix <- limma::makeContrasts(
        contrasts = contrast_strings,
        levels = colnames(design_m)
    )

    # Run limma analysis
    logger::log_info("   Running contrasts.fit...")
    cfit <- contrastsFitFn(fit, contrasts = contr.matrix)

    logger::log_info("   Running eBayes...")
    eb.fit <- limma::eBayes(cfit, trend = eBayes_trend, robust = eBayes_robust)

    # Run treat or standard analysis
    qvalue_failures <- list()

    if (!is.na(treat_lfc_cutoff) && treat_lfc_cutoff > 0) {
        logger::log_info("   Running TREAT analysis with LFC cutoff...")
        # KIT-201 (VAL-127, VAL-131): limma's treat() re-estimates the prior itself (trend and robust FALSE by default);
        # it takes the moderated test's settings, except while a DA made before that change replays.
        moderated <- workflowTreatModerationCurrent()
        t.fit <- limma::treat(eb.fit, lfc = as.double(treat_lfc_cutoff), trend = moderated && isTRUE(eBayes_trend),
            robust = moderated && isTRUE(eBayes_robust))
        extract_fn <- function(fit, coef) limma::topTreat(fit, coef = coef, n = Inf)
    } else {
        logger::log_info("   Running standard analysis...")
        t.fit <- eb.fit
        extract_fn <- function(fit, coef) limma::topTable(fit, coef = coef, number = Inf)
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

    list(
        results = result_tables,
        fit.eb = t.fit,
        qvalue_warnings = qvalue_failures
    )
}

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
#' @param group_id,technical_replicate_id Optional names of the design's group column and of the column that names
#'   the biological unit of each run. When given and the design has technical replicates (several runs of one unit,
#'   not spanned by the group), the model is fitted with the unit as block (`limma::duplicateCorrelation`, then
#'   `limma::lmFit(block, correlation)`). `NULL` (the default) fits every run as an independent sample.
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
  eBayes_robust = TRUE,
  group_id = NULL,
  technical_replicate_id = NULL
) {
    logger::log_info("--- Entering runTestsContrastsMetabDA ---")
    logger::log_info(sprintf(
        "   data dims = %d x %d, %d contrasts",
        nrow(data), ncol(data), length(contrast_strings)
    ))
    logger::log_info(sprintf("   contrasts = %s", paste(contrast_strings, collapse = ", ")))

    # KIT-076: the two halves, in the order they ran.
    fit <- fitMetabDaAssayModel(data, design_matrix, formula_string, sample_id_col, group_id, technical_replicate_id)
    results <- applyMetabDaContrasts(fit, contrast_strings, treat_lfc_cutoff, eBayes_trend, eBayes_robust)
    logger::log_info("--- Exiting runTestsContrastsMetabDA ---")

    return(results)
}
