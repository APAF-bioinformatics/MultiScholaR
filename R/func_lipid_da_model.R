# KIT-076: the fit half of `runTestsContrastsLipidDA`, moved unchanged: the model matrix on every sample, the samples the data
# shares with it, and `lmFit`. It depends on no contrast, so one saved fit serves every contrast of an assay.
#' @param data,design_matrix,formula_string,sample_id_col As `runTestsContrastsLipidDA`.
#' @param group_id,technical_replicate_id The object's declared columns (KIT-128). Without `technical_replicate_id` the
#'   fit is the unblocked `lmFit`, as before. With it the proteomics rule (`resolveProtTechnicalReplicateBlock`) reads the
#'   design rows of the fitted samples: declared technical replicates are fitted with the biological unit as block
#'   (`duplicateCorrelation`, then `lmFit(block, correlation)`; limma User's Guide 11.3 and 9.7), any other design as before.
#' @return The uncontrasted `MArrayLM`; its `design` is the model matrix on the shared samples. With
#'   `technical_replicate_id` it also carries `block_reason`, `uses_duplicate_correlation` and `consensus_correlation`.
#' @noRd
fitLipidDaAssayModel <- function(data, design_matrix, formula_string, sample_id_col = "Run", group_id = NULL,
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

# KIT-076: the contrast half of `runTestsContrastsLipidDA`, moved unchanged: the terms checked against the model, then
# `contrasts.fit`, `eBayes`, `treat` when asked, the extraction and the q-values, for one or more contrasts.
#' @param fit The fit half's `MArrayLM`.
#' @param contrast_strings,treat_lfc_cutoff,eBayes_trend,eBayes_robust As `runTestsContrastsLipidDA`.
#' @param contrastsFitFn The contrast step: `limma::contrasts.fit` (today's), or the exact one of the factors-in-the-model
#'   design (`ruvExactContrastsFit`, decision M29).
#' @return `list(results, fit.eb, qvalue_warnings)` for those contrasts, with `features_without_values` (their ids)
#'   when some feature has no value in any fitted sample (KIT-196: not tested, a row without statistics per table).
#' @noRd
applyLipidDaContrasts <- function(fit, contrast_strings, treat_lfc_cutoff = NA, eBayes_trend = TRUE, eBayes_robust = TRUE,
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

    # KIT-196: a feature with no value in any fitted sample has no mean (Amean NaN), so limma's empirical-Bayes step
    # cannot use it (the trend covariate; limma 3.68.4 fitFDistUnequalDF1.R:17). It carries no information and is not
    # tested: it leaves the moderation here and is returned as a row without statistics in every contrast table.
    without_values <- if (!is.null(cfit$Amean)) rownames(cfit$coefficients)[is.na(cfit$Amean)] else character()
    if (length(without_values) == nrow(cfit$coefficients)) {
        rlang::abort("no feature has a value in any sample, so nothing can be tested", class = "omics_da_no_observed_features")
    }
    if (length(without_values)) {
        logger::log_info(sprintf("   %d feature(s) with no value in any sample are not tested", length(without_values)))
        cfit <- cfit[!is.na(cfit$Amean), ]
    }

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

                # KIT-196: the features without values, as rows without statistics (AveExpr NA), typed as the table's.
                if (length(without_values)) {
                    untested <- da_tbl[rep(NA_integer_, length(without_values)), , drop = FALSE]
                    rownames(untested) <- without_values
                    da_tbl <- rbind(da_tbl, untested)
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

    result <- list(
        results = result_tables,
        fit.eb = t.fit,
        qvalue_warnings = qvalue_failures
    )
    # KIT-196: present only when some feature had no value in any sample, as `failed_assays` is.
    if (length(without_values)) result$features_without_values <- without_values
    result
}

# ----------------------------------------------------------------------------
# runTestsContrastsLipidDA
# ----------------------------------------------------------------------------
#' Run limma differential abundance analysis for a single assay
#'
#' @description Performs limma-based differential abundance analysis on a
#'   single lipidomics assay matrix. A feature with no value in any sample is not tested and
#'   is listed (`features_without_values`).
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
#' @param group_id,technical_replicate_id Optional names of the design's group column and of the column that names
#'   the biological unit of each run. When given and the design has technical replicates (several runs of one unit,
#'   not spanned by the group), the model is fitted with the unit as block (`limma::duplicateCorrelation`, then
#'   `limma::lmFit(block, correlation)`). `NULL` (the default) fits every run as an independent sample.
#'
#' @return A list containing:
#'   - results: Named list of DE result tables (one per contrast)
#'   - fit.eb: The eBayes fit object
#'   - qvalue_warnings: List of contrasts where qvalue() failed
#'   - features_without_values: present only when the assay held features with no value in any sample; their
#'     ids. Such a feature is not tested: it is a row of every table without statistics (`AveExpr` missing)
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
  eBayes_robust = TRUE,
  group_id = NULL,
  technical_replicate_id = NULL
) {
    logger::log_info("--- Entering runTestsContrastsLipidDA ---")
    logger::log_info(sprintf(
        "   data dims = %d x %d, %d contrasts",
        nrow(data), ncol(data), length(contrast_strings)
    ))
    logger::log_info(sprintf("   contrasts = %s", paste(contrast_strings, collapse = ", ")))

    # KIT-076: the two halves, in the order they ran.
    fit <- fitLipidDaAssayModel(data, design_matrix, formula_string, sample_id_col, group_id, technical_replicate_id)
    results <- applyLipidDaContrasts(fit, contrast_strings, treat_lfc_cutoff, eBayes_trend, eBayes_robust)
    logger::log_info("--- Exiting runTestsContrastsLipidDA ---")

    return(results)
}

# KIT-076: the orchestrator's per-assay and per-contrast steps, moved unchanged so that `runLipidsDA` and the
# per-contrast declaration build the same tables.

#' The assay's expression matrix and annotation table.
#' @return `list(expr_matrix, annotation_info)`, or NULL when no sample of the assay is in the design.
#' @noRd
prepareLipidDaAssayInput <- function(assay_data, assay_name, design_matrix, sample_id_col, lipid_id_col, annotation_col) {
    # Build expression matrix (lipids as rows, samples as columns)
    # First, identify sample columns (those that match design matrix)
    sample_cols <- intersect(colnames(assay_data), design_matrix[[sample_id_col]])

    if (length(sample_cols) == 0) {
        logger::log_warn(sprintf("   No matching samples for assay: %s, skipping", assay_name))
        return(NULL)
    }

    # Extract numeric data and set rownames to lipid IDs
    expr_matrix <- as.matrix(assay_data[, sample_cols, drop = FALSE])
    rownames(expr_matrix) <- assay_data[[lipid_id_col]]

    # Get annotation info for later joining
    lipid_names <- if (annotation_col %in% colnames(assay_data)) {
        assay_data[[annotation_col]]
    } else {
        assay_data[[lipid_id_col]]
    }
    annotation_info <- data.frame(
        lipid_id = assay_data[[lipid_id_col]],
        lipid_name = lipid_names,
        stringsAsFactors = FALSE
    )

    logger::log_info(sprintf(
        "   Expression matrix: %d lipids x %d samples",
        nrow(expr_matrix), ncol(expr_matrix)
    ))

    list(expr_matrix = expr_matrix, annotation_info = annotation_info)
}

#' One contrast's table of one assay: the id column, the labels, the annotation, the significance call and the
#' sample intensity columns.
#' @noRd
formatLipidDaContrastTable <- function(da_tbl, contrast_name, assay_name, friendly_name, annotation_info, expr_matrix,
                                        design_matrix, sample_id_col, group_col, da_q_val_thresh, treat_lfc_cutoff) {
    # Add lipid_id from rownames
    da_tbl <- tibble::rownames_to_column(da_tbl, "lipid_id")

    # Add metadata columns
    da_tbl$assay <- assay_name
    da_tbl$comparison <- contrast_name

    # Add friendly name
    da_tbl$friendly_name <- friendly_name

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

#' The up, down and not-significant counts of one results table.
#' @noRd
countLipidDaSignificant <- function(df) {
    if (is.null(df) || nrow(df) == 0) {
        return(list(up = 0, down = 0, ns = 0))
    }
    list(
        up = sum(df$significant == "Up", na.rm = TRUE),
        down = sum(df$significant == "Down", na.rm = TRUE),
        ns = sum(df$significant == "NS", na.rm = TRUE)
    )
}

# ----------------------------------------------------------------------------
# runLipidsDA
# ----------------------------------------------------------------------------
#' Run differential abundance analysis on all assays in a LipidomicsAssayData object
#'
#' @description Main entry point for lipidomics differential abundance analysis.
#'   Loops over all assays in the LipidomicsAssayData object, runs limma DA on each,
#'   and aggregates results with an assay identifier column. A feature with no value in any sample is not
#'   tested and is listed (`features_without_values`).
#'
#' @param theObject LipidomicsAssayData S4 object containing lipid data.
#' @param contrasts_tbl Data frame with contrast definitions. Must have a column
#'   named "contrasts" or "contrast_string" containing the contrast formulas.
#' @param formula_string Model formula (default "~ 0 + group").
#' @param da_q_val_thresh Q-value threshold for significance (default 0.05).
#' @param treat_lfc_cutoff Log fold-change cutoff for TREAT (default 0, standard analysis).
#' @param eBayes_trend Logical, use trend in eBayes (default TRUE).
#' @param eBayes_robust Logical, use robust eBayes (default TRUE).
#' @param ruv_factors NULL (the default): the test reads the object's values with `formula_string`, as before.
#'   Or the per-assay RUV-III factors of the factors-in-the-model design (KIT-127): `theObject` is then the
#'   uncorrected object (the corrected object is refused), and per assay the factors W are added to the design and
#'   to the model and the supplied controls are not tested. A feature is left out only when its observed samples
#'   cannot estimate a W coefficient. A feature with no value in one whole group is tested in the comparisons
#'   between the groups where it is observed (limma's handling of missing values) and listed as not observed in
#'   that group. A model that leaves no residual degree of freedom stops. Fold changes and p-values come from the
#'   uncorrected values with W in the model; the RUV-corrected matrix is for display and export and need not
#'   reproduce them exactly.
#'
#' @return A list containing:
#'   - theObject: The input S4 object (unchanged)
#'   - contrasts_results: Per-assay limma fit objects
#'   - da_lipids_long: Combined long-format results with assay column
#'   - per_assay_results: List of per-assay result tables
#'   - significant_counts: Summary of significant lipids per assay
#'   - qvalue_warnings: Any qvalue computation warnings
#'   - failed_assays: Present only when an assay's analysis stopped: the message per failed assay. The
#'     other assays are kept and a warning (class `omics_da_assay_failed`) names the failed ones; when every
#'     assay failed the function stops (class `omics_da_all_assays_failed`).
#'   - features_without_values: present only when an assay held features with no value in any sample; per such
#'     assay, their ids. Such a feature is not tested: it is a row of every contrast table without statistics
#'     (`AveExpr` missing), and the other features' results are those of the data without it
#'   - ruv_features_not_tested, ruv_features_not_observed: with `ruv_factors` only. Per assay, the ids not tested
#'     (`supplied_controls`, `not_estimable`), and a table (`feature`, `group`) of the tested features with no
#'     value in one whole group
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
  eBayes_robust = TRUE,
  ruv_factors = NULL
) {

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


    preflight <- validateLipidDaAnalysisPreflight(
        currentS4 = theObject,
        contrastsTbl = contrasts_tbl,
        formulaString = formula_string
    )
    if (!isTRUE(preflight$valid)) {
        stop(paste(preflight$errors, collapse = "; "))
    }

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


    assay_names <- names(assay_list)
    logger::log_info(sprintf(
        "   Processing %d assays: %s",
        length(assay_names), paste(assay_names, collapse = ", ")
    ))
    logger::log_info(sprintf("   Contrasts: %s", paste(contrast_strings, collapse = ", ")))

    # KIT-127: an object the factors kernel has corrected is not tested with the factors as well.
    if (!is.null(ruv_factors)) refuseOmicsRuvCorrectedObject(theObject)
    ruv_not_tested <- list()
    ruv_not_observed <- list()

    # Initialize results storage
    per_assay_results <- list()
    contrasts_results <- list()
    all_qvalue_warnings <- list()
    # KIT-128: the assays whose kernel stopped, with the message.
    failed_assays <- list()
    # KIT-196: per assay, the features with no value in any sample (not tested).
    features_without_values <- list()
    # KIT-128 K128-1: a design without biological replication is refused before any assay.
    refuseNoBiologicalReplication(resolveProtTechnicalReplicateBlock(as.data.frame(design_matrix), group_col,
        theObject@technical_replicate_id), design_matrix[[group_col]])

    # Process each assay
    for (assay_name in assay_names) {
        logger::log_info(sprintf("--- Processing assay: %s ---", assay_name))

        assay_data <- assay_list[[assay_name]]

        # KIT-076: the assay's matrix and annotation, the contrast path's too.
        assay_input <- prepareLipidDaAssayInput(assay_data, assay_name, design_matrix, sample_id_col, lipid_id_col, annotation_col)
        if (is.null(assay_input)) next
        expr_matrix <- assay_input$expr_matrix
        annotation_info <- assay_input$annotation_info
        # KIT-127: with RUV factors the test reads these (uncorrected) values with W in the model.
        da_design <- design_matrix
        da_formula <- formula_string
        if (!is.null(ruv_factors)) {
            da_input <- applyOmicsRuvFactorsToDaInput(expr_matrix, design_matrix, formula_string, sample_id_col,
                ruv_factors[[assay_name]], assay_name, group_col, theObject@technical_replicate_id)
            ruv_not_tested[[assay_name]] <- da_input$not_tested
            ruv_not_observed[[assay_name]] <- da_input$not_observed
            expr_matrix <- da_input$expr_matrix
            da_design <- da_input$design_matrix
            da_formula <- da_input$formula_string
        }

        # Run limma DE
        tryCatch(
            {
                # KIT-129 B1: with RUV factors, a technical-replicate correlation that cannot be estimated names k.
                # Decision M29: an assay with W in the model forms each comparison's standard error exactly.
                limma_results <- withOmicsRuvKOnBlockFailure(if (!is.null(ruv_factors[[assay_name]]$W)) {
                    runOmicsRuvFactorsDaKernel(fitLipidDaAssayModel, applyLipidDaContrasts, expr_matrix, contrast_strings, da_design, da_formula,
                        sample_id_col, if (treat_lfc_cutoff > 0) treat_lfc_cutoff else NA, eBayes_trend, eBayes_robust,
                        group_col, theObject@technical_replicate_id)
                } else {
                    runTestsContrastsLipidDA(
                        data = expr_matrix,
                        contrast_strings = contrast_strings,
                        design_matrix = da_design,
                        formula_string = da_formula,
                        sample_id_col = sample_id_col,
                        treat_lfc_cutoff = if (treat_lfc_cutoff > 0) treat_lfc_cutoff else NA,
                        eBayes_trend = eBayes_trend,
                        eBayes_robust = eBayes_robust,
                        # KIT-128: the declared technical replicates, for the block.
                        group_id = group_col,
                        technical_replicate_id = theObject@technical_replicate_id
                    )
                }, if (!is.null(ruv_factors)) ruv_factors[[assay_name]])

                # Collect qvalue warnings
                if (length(limma_results$qvalue_warnings) > 0) {
                    all_qvalue_warnings[[assay_name]] <- limma_results$qvalue_warnings
                }

                # Process results for each contrast
                assay_results_list <- purrr::map2(
                    limma_results$results,
                    names(limma_results$results),
                    function(da_tbl, contrast_name) {
                        idx <- match(contrast_name, contrast_strings)
                        formatLipidDaContrastTable(da_tbl, contrast_name, assay_name,
                            if (!is.na(idx)) friendly_names[idx] else contrast_name, annotation_info, expr_matrix, design_matrix,
                            sample_id_col, group_col, da_q_val_thresh, treat_lfc_cutoff)
                    }
                )

                # Combine all contrasts for this assay
                assay_combined <- dplyr::bind_rows(assay_results_list)
                per_assay_results[[assay_name]] <- assay_combined
                # Store fit object for Glimma (KIT-128 K128-5: only once the assay's tables are formatted).
                contrasts_results[[assay_name]] <- limma_results
                if (length(limma_results$features_without_values)) {
                    features_without_values[[assay_name]] <- limma_results$features_without_values
                }

                logger::log_info(sprintf(
                    "   Completed assay %s: %d total results",
                    assay_name, nrow(assay_combined)
                ))
            },
            error = function(e) {
                cat(sprintf("[ERROR] Processing assay %s: %s\n", assay_name, e$message))
                logger::log_error(sprintf("   ERROR processing assay %s: %s", assay_name, e$message))
                per_assay_results[[assay_name]] <<- NULL
                # KIT-128: the root cause; an error of the formatting step arrives wrapped by purrr (its message is the index).
                cause <- e
                while (inherits(cause$parent, "condition")) cause <- cause$parent
                failed_assays[[assay_name]] <<- conditionMessage(cause)
            }
        )
    }

    # KIT-128 (parked row 24 b): a failed assay is named and warned about and the others are kept; the run stops only
    # when every assay failed (before: printed, logged and dropped, down to an empty table).
    if (length(failed_assays) > 0) {
        failures <- stats::setNames(sprintf("%s: %s", names(failed_assays), unlist(failed_assays)), rep("x", length(failed_assays)))
        if (length(per_assay_results) == 0) {
            rlang::abort(c("differential abundance failed for every assay", failures), class = "omics_da_all_assays_failed")
        }
        rlang::warn(c("differential abundance failed for some assays; the other assays are kept", failures),
            class = "omics_da_assay_failed")
    }

    # Combine all assay results
    da_lipids_long <- dplyr::bind_rows(per_assay_results)
    logger::log_info(sprintf("   Combined results: %d total rows", nrow(da_lipids_long)))

    # Calculate significant counts per assay
    significant_counts <- purrr::map(per_assay_results, countLipidDaSignificant) |> purrr::set_names(names(per_assay_results))

    logger::log_info("=== Completed runLipidsDA ===")

    result <- list(
        theObject = theObject,
        contrasts_results = contrasts_results,
        da_lipids_long = da_lipids_long,
        per_assay_results = per_assay_results,
        significant_counts = significant_counts,
        qvalue_warnings = all_qvalue_warnings,
        da_q_val_thresh = da_q_val_thresh,
        treat_lfc_cutoff = treat_lfc_cutoff
    )
    if (length(failed_assays) > 0) result$failed_assays <- failed_assays
    # KIT-196: present only when some assay had features without values, as `failed_assays` is.
    if (length(features_without_values) > 0) result$features_without_values <- features_without_values
    # KIT-127: with RUV factors, what was not tested and what was not observed in a whole group, by assay.
    if (!is.null(ruv_factors)) {
        result$ruv_features_not_tested <- ruv_not_tested
        result$ruv_features_not_observed <- ruv_not_observed
    }
    return(result)
}
