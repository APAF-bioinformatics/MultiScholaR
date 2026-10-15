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
  contrast_expressions <- purrr::map_chr(contrast_strings, function(x) {
    parts <- strsplit(x, "=")[[1]]
    if (length(parts) >= 2) {
      purrr::pluck(parts, 2, .default = x)
    } else {
      x
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

  # Extract results for each contrast using functional purrr::map
  c_names <- names(contrast_strings)
  results_list <- purrr::map(seq_along(contrast_strings), function(i) {
    c_str <- contrast_strings[i]
    contrast_name <- if (!is.null(c_names) && length(c_names) >= i && nzchar(c_names[i])) {
      c_names[i]
    } else {
      parts <- strsplit(c_str, "=")[[1]]
      if (length(parts) >= 2 && nzchar(parts[1])) {
        parts[1]
      } else {
        paste0("Contrast_", i)
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

    computed_qvalues <- tryCatch({
      if (requireNamespace("qvalue", quietly = TRUE) && length(contrast_results$P.Value) >= 3) {
        qvalue::qvalue(contrast_results$P.Value)$qvalues
      } else {
        contrast_results$adj.P.Val
      }
    }, error = function(e) {
      contrast_results$adj.P.Val
    })

    contrast_results <- contrast_results |>
      dplyr::mutate(
        fdr_qvalue = computed_qvalues,
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
        dplyr::everything(),
        -adj.P.Val # remove original adj.P.Val to avoid confusion
      )

    cat("   convertDpcDAToStandardFormat: Processed contrast", contrast_name, "with", nrow(contrast_results), "proteins\n")

    contrast_results
  })

  names(results_list) <- contrast_strings

  # Return in same format as runTestsContrasts, which is a list of tables
  return_object <- list(
    results = results_list, # Return the list of data frames
    fit.eb = eb_fit,
    dpc_method_used = TRUE
  )

  cat("   convertDpcDAToStandardFormat: Conversion completed successfully\n")

  return(return_object)
}

