# KIT-076: the orchestrator's per-assay and per-contrast steps, moved unchanged so that `runMetabolitesDA` and the
# per-contrast declaration build the same tables.

#' The assay's expression matrix and annotation table.
#' @return `list(expr_matrix, annotation_info)`, or NULL when no sample of the assay is in the design.
#' @noRd
prepareMetabDaAssayInput <- function(assay_data, assay_name, design_matrix, sample_id_col, metabolite_id_col, annotation_col) {
    # Build expression matrix (metabolites as rows, samples as columns)
    # First, identify sample columns (those that match design matrix)
    sample_cols <- intersect(colnames(assay_data), design_matrix[[sample_id_col]])

    if (length(sample_cols) == 0) {
        logger::log_warn(sprintf("   No matching samples for assay: %s, skipping", assay_name))
        return(NULL)
    }

    # Extract numeric data and set rownames to metabolite IDs
    expr_matrix <- as.matrix(assay_data[, sample_cols, drop = FALSE])
    rownames(expr_matrix) <- assay_data[[metabolite_id_col]]

    # Get annotation info for later joining. Custom imports may not have a
    # distinct annotation column; in that case keep DA output usable by
    # falling back to the primary metabolite ID.
    annotation_values <- if (
        length(annotation_col) == 1L &&
        !is.na(annotation_col) &&
        nzchar(annotation_col) &&
        annotation_col %in% colnames(assay_data)
    ) {
        assay_data[[annotation_col]]
    } else {
        assay_data[[metabolite_id_col]]
    }

    annotation_info <- data.frame(
        metabolite_id = assay_data[[metabolite_id_col]],
        metabolite_name = annotation_values,
        stringsAsFactors = FALSE
    )

    logger::log_info(sprintf(
        "   Expression matrix: %d metabolites x %d samples",
        nrow(expr_matrix), ncol(expr_matrix)
    ))

    list(expr_matrix = expr_matrix, annotation_info = annotation_info)
}

#' One contrast's table of one assay: the id column, the labels, the annotation, the significance call and the
#' sample intensity columns.
#' @noRd
formatMetabDaContrastTable <- function(da_tbl, contrast_name, assay_name, friendly_name, annotation_info, expr_matrix,
                                        design_matrix, sample_id_col, group_col, da_q_val_thresh, treat_lfc_cutoff) {
    tryCatch(
        {

            # Validate da_tbl
            if (!is.data.frame(da_tbl)) {
                stop(sprintf("da_tbl is not a data frame, it is a %s", class(da_tbl)[1]))
            }

            # Add metabolite_id from rownames
            da_tbl <- tibble::rownames_to_column(da_tbl, "metabolite_id")

            # Add metadata columns
            da_tbl$assay <- assay_name
            da_tbl$comparison <- contrast_name

            # Add friendly name
            da_tbl$friendly_name <- friendly_name

            # Join annotation info
            da_tbl <- dplyr::left_join(
                da_tbl,
                annotation_info,
                by = "metabolite_id"
            )

            # Classify significance
            da_tbl$significant <- ifelse(
                da_tbl$fdr_qvalue < da_q_val_thresh & abs(da_tbl$logFC) >= treat_lfc_cutoff,
                ifelse(da_tbl$logFC > 0, "Up", "Down"),
                "NS"
            )

            # Add sample intensity columns using createMetabDaResultsLongFormat
            da_tbl <- createMetabDaResultsLongFormat(
                lfc_qval_tbl = da_tbl,
                expr_matrix = expr_matrix,
                design_matrix = design_matrix,
                sample_id_col = sample_id_col,
                group_id_col = group_col,
                metabolite_id_col = "metabolite_id"
            )

            return(da_tbl)
        },
        error = function(inner_e) {
            cat(sprintf("[ERROR] Inner loop failure for contrast %s: %s\n", contrast_name, inner_e$message))
            stop(inner_e)
        }
    )
}

#' The up, down and not-significant counts of one results table.
#' @noRd
countMetabDaSignificant <- function(df) {
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
# runMetabolitesDA
# ----------------------------------------------------------------------------
#' Run differential abundance analysis on all assays in a MetaboliteAssayData object
#'
#' @description Main entry point for metabolomics differential abundance analysis.
#'   Loops over all assays in the MetaboliteAssayData object, runs limma DA on each,
#'   and aggregates results with an assay identifier column. A feature with no value in any sample is not
#'   tested and is listed (`features_without_values`).
#'
#' @param theObject MetaboliteAssayData S4 object containing metabolite data.
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
#'   - da_metabolites_long: Combined long-format results with assay column
#'   - per_assay_results: List of per-assay result tables
#'   - significant_counts: Summary of significant metabolites per assay
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
runMetabolitesDA <- function(
  theObject,
  contrasts_tbl,
  formula_string = "~ 0 + group",
  da_q_val_thresh = 0.05,
  treat_lfc_cutoff = 0,
  eBayes_trend = TRUE,
  eBayes_robust = TRUE,
  ruv_factors = NULL
) {

    logger::log_info("=== Starting runMetabolitesDA ===")

    # Validate input
    if (!inherits(theObject, "MetaboliteAssayData")) {
        stop("theObject must be a MetaboliteAssayData S4 object")
    }

    # Extract object slots
    assay_list <- theObject@metabolite_data
    design_matrix <- theObject@design_matrix
    sample_id_col <- theObject@sample_id
    group_col <- theObject@group_id
    metabolite_id_col <- theObject@metabolite_id_column
    annotation_col <- theObject@annotation_id_column


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
        assay_input <- prepareMetabDaAssayInput(assay_data, assay_name, design_matrix, sample_id_col, metabolite_id_col, annotation_col)
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

        # Run limma DA
        tryCatch(
            {
                # KIT-129 B1: with RUV factors, a technical-replicate correlation that cannot be estimated names k.
                # Decision M29: an assay with W in the model forms each comparison's standard error exactly.
                limma_results <- withOmicsRuvKOnBlockFailure(if (!is.null(ruv_factors[[assay_name]]$W)) {
                    runOmicsRuvFactorsDaKernel(fitMetabDaAssayModel, applyMetabDaContrasts, expr_matrix, contrast_strings, da_design, da_formula,
                        sample_id_col, if (treat_lfc_cutoff > 0) treat_lfc_cutoff else NA, eBayes_trend, eBayes_robust,
                        group_col, theObject@technical_replicate_id)
                } else {
                    runTestsContrastsMetabDA(
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
                        formatMetabDaContrastTable(da_tbl, contrast_name, assay_name,
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
    da_metabolites_long <- dplyr::bind_rows(per_assay_results)
    logger::log_info(sprintf("   Combined results: %d total rows", nrow(da_metabolites_long)))

    # Calculate significant counts per assay
    significant_counts <- purrr::map(per_assay_results, countMetabDaSignificant) |> purrr::set_names(names(per_assay_results))

    logger::log_info("=== Completed runMetabolitesDA ===")

    result <- list(
        theObject = theObject,
        contrasts_results = contrasts_results,
        da_metabolites_long = da_metabolites_long,
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
