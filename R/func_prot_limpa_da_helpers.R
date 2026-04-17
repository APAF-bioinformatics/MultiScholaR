# ----------------------------------------------------------------------------
# convertDpcDAToStandardFormat
# ----------------------------------------------------------------------------
#' Convert limpa dpcDE Results to Standard Format
#'
#' This function converts the output from limpa::dpcDE() to the same format as
#' runTestsContrasts() to ensure compatibility with existing DE analysis workflows.
#'
#' @param dpc_fit MArrayLM object from limpa::dpcDE()
#' @param contrast_strings Character vector of contrast strings in "name=expression" format
#' @param design_matrix Design matrix used in the analysis
#' @param eBayes_trend Logical, whether empirical Bayes trend was used
#' @param eBayes_robust Logical, whether robust empirical Bayes was used
#'
#' @return List with same structure as runTestsContrasts output
#' @export
convertDpcDAToStandardFormat <- function(dpc_fit,
                                         contrast_strings,
                                         design_matrix,
                                         eBayes_trend = TRUE,
                                         eBayes_robust = TRUE) {
  if (!requireNamespace("limma", quietly = TRUE)) {
    stop("limma package is required for convertDpcDAToStandardFormat")
  }

  cat("   convertDpcDAToStandardFormat: Converting dpcDE results to standard format\n")
  cat("   convertDpcDAToStandardFormat: Processing", length(contrast_strings), "contrasts\n")

  # Create contrast matrix from contrast strings
  # Extract just the contrast expressions (after the "=" sign)
  contrast_expressions <- sapply(contrast_strings, function(x) {
    parts <- strsplit(x, "=")[[1]]
    if (length(parts) == 2) {
      return(parts[2])
    } else {
      return(x) # Fallback if no "=" found
    }
  })

  # Create contrast matrix using limma
  contrast_matrix <- limma::makeContrasts(
    contrasts = contrast_expressions,
    levels = colnames(design_matrix)
  )

  cat("   convertDpcDAToStandardFormat: Contrast matrix dims:", nrow(contrast_matrix), "x", ncol(contrast_matrix), "\n")

  # Apply contrasts to the dpcDE fit
  contrast_fit <- limma::contrasts.fit(dpc_fit, contrast_matrix)

  # Apply empirical Bayes
  eb_fit <- limma::eBayes(contrast_fit, trend = eBayes_trend, robust = eBayes_robust)

  # Extract results for each contrast
  results_list <- list()

  for (i in seq_along(contrast_strings)) {
    contrast_name <- names(contrast_strings)[i]
    if (is.null(contrast_name) || contrast_name == "") {
      # Extract friendly name from contrast string if no name provided
      parts <- strsplit(contrast_strings[i], "=")[[1]]
      if (length(parts) == 2) {
        contrast_name <- parts[1]
      } else {
        contrast_name <- paste0("Contrast_", i)
      }
    }

    # Get results for this contrast using limma::topTable
    contrast_results <- limma::topTable(
      eb_fit,
      coef = i,
      number = Inf,
      sort.by = "none", # Keep original order
      adjust.method = "BH"
    )

    # Rename columns to match runTestsContrasts output format, and add qvalue
    if (!"P.Value" %in% colnames(contrast_results)) {
      stop("P.Value column not found in topTable results for dpcDE.")
    }

    contrast_results <- contrast_results |>
      dplyr::mutate(
        # Ensure the qvalue library is used for fdr_qvalue, consistent with runTestsContrasts
        fdr_qvalue = qvalue::qvalue(P.Value)$qvalues,
        # Keep the BH adjustment in its own column for full compatibility
        fdr_value_bh_adjustment = adj.P.Val
      ) |>
      dplyr::rename(
        raw_pvalue = P.Value
      ) |>
      dplyr::mutate(
        comparison = contrast_name,
        uniprot_acc = rownames(contrast_results),
        log_intensity = AveExpr
      ) |>
      dplyr::select(
        uniprot_acc,
        comparison,
        logFC,
        log_intensity,
        raw_pvalue,
        fdr_qvalue,
        fdr_value_bh_adjustment,
        everything(),
        -adj.P.Val # remove original adj.P.Val to avoid confusion
      )

    # Store in results list (matching runTestsContrasts structure)
    # The name of the list element MUST be the full contrast string for downstream parsing
    full_contrast_string <- contrast_strings[i]
    results_list[[full_contrast_string]] <- contrast_results

    cat("   convertDpcDAToStandardFormat: Processed contrast", contrast_name, "with", nrow(contrast_results), "proteins\n")
  }

  # Return in same format as runTestsContrasts, which is a list of tables
  return_object <- list(
    results = results_list, # Return the list of data frames
    fit.eb = eb_fit,
    dpc_method_used = TRUE
  )

  cat("   convertDpcDAToStandardFormat: Conversion completed successfully\n")

  return(return_object)
}

