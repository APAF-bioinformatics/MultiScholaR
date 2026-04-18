# ----------------------------------------------------------------------------
# runMetabolitesDA
# ----------------------------------------------------------------------------
#' Run differential abundance analysis on all assays in a MetaboliteAssayData object
#'
#' @description Main entry point for metabolomics differential abundance analysis.
#'   Loops over all assays in the MetaboliteAssayData object, runs limma DA on each,
#'   and aggregates results with an assay identifier column.
#'
#' @param theObject MetaboliteAssayData S4 object containing metabolite data.
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
#'   - da_metabolites_long: Combined long-format results with assay column
#'   - per_assay_results: List of per-assay result tables
#'   - significant_counts: Summary of significant metabolites per assay
#'   - qvalue_warnings: Any qvalue computation warnings
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
  eBayes_robust = TRUE
) {
    # [D66:START] -------------------------
    d66_log <- function(...) message(sprintf("[D66] %s", paste0(...)))
    cat("[D66] === ENTER runMetabolitesDA [AG-v3] ===\n")
    d66_log("  formula_string = ", formula_string)
    d66_log("  contrasts_tbl columns = ", paste(colnames(contrasts_tbl), collapse = ", "))
    d66_log("  contrasts_tbl$contrasts = ", paste(contrasts_tbl$contrasts, collapse = ", "))
    # [D66:END] ---------------------------

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

        # Build expression matrix (metabolites as rows, samples as columns)
        # First, identify sample columns (those that match design matrix)
        sample_cols <- intersect(colnames(assay_data), design_matrix[[sample_id_col]])

        if (length(sample_cols) == 0) {
            logger::log_warn(sprintf("   No matching samples for assay: %s, skipping", assay_name))
            next
        }

        # Extract numeric data and set rownames to metabolite IDs
        expr_matrix <- as.matrix(assay_data[, sample_cols, drop = FALSE])
        rownames(expr_matrix) <- assay_data[[metabolite_id_col]]

        # Get annotation info for later joining
        annotation_info <- data.frame(
            metabolite_id = assay_data[[metabolite_id_col]],
            metabolite_name = assay_data[[annotation_col]],
            stringsAsFactors = FALSE
        )

        logger::log_info(sprintf(
            "   Expression matrix: %d metabolites x %d samples",
            nrow(expr_matrix), ncol(expr_matrix)
        ))

        # Run limma DA
        tryCatch(
            {
                limma_results <- runTestsContrastsMetabDA(
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

                cat(sprintf("[D66] Assay %s: limma_results$results length = %d\n", assay_name, length(limma_results$results)))
                cat(sprintf("[D66] Assay %s: names(limma_results$results) = [%s]\n", assay_name, paste(names(limma_results$results), collapse = ", ")))

                # Process results for each contrast
                assay_results_list <- purrr::map2(
                    limma_results$results,
                    names(limma_results$results),
                    function(da_tbl, contrast_name) {
                        tryCatch(
                            {
                                cat(sprintf("[D66] Processing contrast: %s (assay: %s)\n", contrast_name, assay_name))

                                # Validate da_tbl
                                if (!is.data.frame(da_tbl)) {
                                    stop(sprintf("da_tbl is not a data frame, it is a %s", class(da_tbl)[1]))
                                }

                                # Add metabolite_id from rownames
                                cat("[D66]   Step: rownames_to_column\n")
                                da_tbl <- tibble::rownames_to_column(da_tbl, "metabolite_id")

                                # Add metadata columns
                                cat("[D66]   Step: adding assay/comparison/friendly_name\n")
                                da_tbl$assay <- assay_name
                                da_tbl$comparison <- contrast_name

                                # Add friendly name
                                idx <- match(contrast_name, contrast_strings)
                                da_tbl$friendly_name <- if (!is.na(idx)) friendly_names[idx] else contrast_name

                                # Join annotation info
                                cat("[D66]   Step: left_join annotation\n")
                                da_tbl <- dplyr::left_join(
                                    da_tbl,
                                    annotation_info,
                                    by = "metabolite_id"
                                )

                                # Classify significance
                                cat("[D66]   Step: classifying significance\n")
                                da_tbl$significant <- ifelse(
                                    da_tbl$fdr_qvalue < da_q_val_thresh & abs(da_tbl$logFC) >= treat_lfc_cutoff,
                                    ifelse(da_tbl$logFC > 0, "Up", "Down"),
                                    "NS"
                                )

                                # Add sample intensity columns using createMetabDaResultsLongFormat
                                cat(sprintf("[D66]   Calling createMetabDaResultsLongFormat for contrast: %s\n", contrast_name))
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
    da_metabolites_long <- dplyr::bind_rows(per_assay_results)
    logger::log_info(sprintf("   Combined results: %d total rows", nrow(da_metabolites_long)))

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

    logger::log_info("=== Completed runMetabolitesDA ===")

    return(list(
        theObject = theObject,
        contrasts_results = contrasts_results,
        da_metabolites_long = da_metabolites_long,
        per_assay_results = per_assay_results,
        significant_counts = significant_counts,
        qvalue_warnings = all_qvalue_warnings,
        da_q_val_thresh = da_q_val_thresh,
        treat_lfc_cutoff = treat_lfc_cutoff
    ))
}

