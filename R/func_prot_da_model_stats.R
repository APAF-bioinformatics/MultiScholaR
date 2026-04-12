# ----------------------------------------------------------------------------
# ebFit
# ----------------------------------------------------------------------------
#' Run the Empircal Bayes Statistics for Differential Expression in the limma package
#' @param ID List of protein accessions / row names.
#' @param design Output from running the function \code{\link{model.matrix}}.
#' @param contr.matrix Output from the function \code{\link{makeContrasts}}.
#' @seealso \code{\link{model.matrix}}
#' @seealso \code{\link{makeContrasts}}
#' @export
ebFit <- function(data, design, contr.matrix) {
  fit <- lmFit(data, design)
  fit.c <- contrasts.fit(fit, contrasts = contr.matrix)

  fit.eb <- suppressWarnings(eBayes(fit.c))

  logFC <- fit.eb$coefficients[, 1]
  df.r <- fit.eb$df.residual
  df.0 <- rep(fit.eb$df.prior, dim(data)[1])
  s2.0 <- rep(fit.eb$s2.prior, dim(data)[1])
  s2 <- (fit.eb$sigma)^2
  s2.post <- fit.eb$s2.post
  t.ord <- fit.eb$coefficients[, 1] /
    fit.eb$sigma /
    fit.eb$stdev.unscaled[, 1]
  t.mod <- fit.eb$t[, 1]
  p.ord <- 2 * pt(-abs(t.ord), fit.eb$df.residual)
  raw_pvalue <- fit.eb$p.value[, 1]
  q.ord <- qvalue(p.ord)$q
  fdr_qvalue <- qvalue(raw_pvalue)$q

  return(list(
    table = data.frame(logFC, t.ord, t.mod, p.ord, raw_pvalue, q.ord, fdr_qvalue, df.r, df.0, s2.0, s2, s2.post),
    fit.eb = fit.eb
  ))
}

# ----------------------------------------------------------------------------
# runTest
# ----------------------------------------------------------------------------
## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#' Analyse one contrast (e.g. compare a pair of experimental groups) and output the q-values per protein.
#' @param ID List of protein accessions / row names.
#' @param A String representing the name of experimental group A for pairwise comparison of B - A.
#' @param B String representing the name of experimental group B for pairwise comparison of B - A.
#' @param group_A Names of all the columns / samples that are in experimental group A.
#' @param group_B Names of all the columns / samples that are in experimental group B.
#' @param design_matrix A data frame with a column containing the sample ID (as per the sample_id param) and the experimental group (as per the group param). Each row as the sample ID as row name in the data frame.
#' @param formula_string A formula string representing the experimental design. e.g. ("~ 0 + group")
#' @param contrast_variable String representing the contrast variable, which is also used in the formula string. (e.g. "group")
#' @param weights Numeric matrix for adjusting each sample and gene.
#' @return A data frame with the following columns:
#' row.names = the protein accessions
#' comparison A string showing log({group B's name}) minus log({group A's name})
#' meanA     mean of the normalised log abundance value of the gene across samples from experimental group A
#' meanB     mean of the normalised log abundance value of the gene across samples from experimental group B
#' logFC     log fold-change
#' tstats    t-test statistics
#' tmod      moderated t-test statistics
#' pval      t-test p-value
#' raw_pvalue      moderated t-test p-value
#' qval      t-test q-value
#' fdr_qvalue     moderated t-test q-value
#' @export
runTest <- function(ID, A, B, group_A, group_B, design_matrix, formula_string,
                    contrast_variable = "group",
                    weights = NA) {
  ff <- as.formula(formula_string)
  mod_frame <- model.frame(ff, design_matrix)
  design_m <- model.matrix(ff, mod_frame)


  # print("My design matrix")
  # print(design_m)
  # print( paste( "nrow(weights)", nrow(weights), "nrow(design_m)", nrow(design_m)))

  if (!is.na(weights)) {
    if (nrow(weights) == nrow(design_m)) {
      design_m <- cbind(design_m, weights)
    } else {
      stop("Stop: nrow(weights) should be equal to nrow(design_m)")
    }
  }

  # print(paste("group_A = ", group_A))
  # print(paste("group_B = ", group_B))

  contr.matrix <- makeContrasts(
    contrasts = paste0(group_B, "vs", group_A, "=", contrast_variable, group_B, "-", contrast_variable, group_A),
    levels = colnames(design_m)
  )

  eb_fit_list <- ebFit(cbind(A, B), design_m, contr.matrix = contr.matrix)

  r <- eb_fit_list$table
  fit.eb <- eb_fit_list$fit.eb

  return(list(
    table = data.frame(
      row.names = row.names(r),
      comparison = paste("log(", group_B, ") minus log(", group_A, ")", sep = ""),
      meanA = rowMeans(A),
      meanB = rowMeans(B),
      logFC = r$logFC,
      tstats = r$t.ord,
      tmod = r$t.mod,
      pval = r$p.ord,
      raw_pvalue = r$raw_pvalue,
      qval = r$q.ord,
      fdr_qvalue = r$fdr_qvalue
    ),
    fit.eb = fit.eb
  ))
}

# ----------------------------------------------------------------------------
# runTests
# ----------------------------------------------------------------------------
#' Compare a pair of experimental groups and output the log fold-change and q-values per protein.
#' @param ID List of protein accessions / row names.
#' @param data Data frame containing the log (base 2) protein abundance values where each column represents a sample and each row represents a protein group, and proteins as rows. The data is preferably median-scaled, with missing values imputed, and batch-effects removed.
#' @param test_pairs Input file with a table listing all the pairs of experimental groups to compare. First column represents group A and second column represents group B. Linear model comparisons (e.g. Contrasts) would be group B minus group A.
#' @param sample_columns A vector of column names (e.g. strings) representing samples which would be used in the statistical tests. Each column contains protein abundance values.
#' @param sample_rows_list A list, the name of each element is the sample ID and each element is a vector containing the protein accessions (e.g. row_id) with enough number of values. It is usually the output from the function \code{get_rows_to_keep_list}.
#' @param type_of_grouping A list where each element name is the name of a treatment group and each element is a vector containing the sample IDs within the treatment group. It is usually the output from the function \code{get_type_of_grouping}.
#' @param design_matrix A data frame with a column containing the sample ID (as per the sample_id param) and the experimental group (as per the group param). Each row as the sample ID as row name in the data frame.
#' @param formula_string A formula string representing the experimental design. e.g. ("~ 0 + group")
#' @param contrast_variable String representing the contrast variable, which is also used in the formula string. (e.g. "group")
#' @param weights Numeric matrix for adjusting each sample and gene.
#' @return A list of data frames, the name of each element represents each pairwise comparison. Each data frame has the following columns:
#' row.names = the protein accessions
#' comparison A string showing log({group B's name}) minus log({group A's name})
#' meanA     mean of the normalised log abundance value of the gene across samples from experimental group A
#' meanB     mean of the normalised log abundance value of the gene across samples from experimental group B
#' logFC     log fold-change
#' tstats    t-test statistics
#' tmod      moderated t-test statistics
#' pval      t-test p-value
#' raw_pvalue      moderated t-test p-value
#' qval      t-test q-value
#' fdr_qvalue     moderated t-test q-value
#' @seealso \code{\link{get_rows_to_keep_list}}
#' @seealso \code{\link{get_type_of_grouping}}
#' @export
runTests <- function(ID, data, test_pairs, sample_columns, sample_rows_list = NA, type_of_grouping, design_matrix, formula_string, contrast_variable = "group", weights = NA) {
  r <- list()
  for (i in 1:nrow(test_pairs)) {
    rows_to_keep <- rownames(data)


    if (length(sample_rows_list) > 0) {
      if (!is.na(sample_rows_list) &
        #  Check that sample group exists as names inside sample_rows_list
        length(which(c(test_pairs[i, "A"], test_pairs[i, "B"]) %in% names(sample_rows_list))) > 0) {
        rows_to_keep <- unique(
          sample_rows_list[[test_pairs[[i, "A"]]]],
          sample_rows_list[[test_pairs[[i, "B"]]]]
        )
      }
    }

    tmp <- data[rows_to_keep, sample_columns]
    rep <- colnames(tmp)

    # print( paste( test_pairs[i,]$A, test_pairs[i,]$B) )
    A <- tmp[, type_of_grouping[test_pairs[i, ]$A][[1]]]
    B <- tmp[, type_of_grouping[test_pairs[i, ]$B][[1]]]

    subset_weights <- NA

    if (!is.na(weights)) {
      subset_weights <- weights[c(colnames(A), colnames(B)), ]
    }

    # print(colnames(A))
    # print(colnames(B))
    tmp <- unname(cbind(A, B))
    Aname <- paste(test_pairs[i, ]$A, 1:max(1, ncol(A)), sep = "_")
    Bname <- paste(test_pairs[i, ]$B, 1:max(1, ncol(B)), sep = "_")
    colnames(tmp) <- c(Aname, Bname)

    selected_sample_ids <- c(type_of_grouping[test_pairs[i, ]$A][[1]], type_of_grouping[test_pairs[i, ]$B][[1]])
    design_matrix_subset <- design_matrix[selected_sample_ids, , drop = FALSE]

    # print("My design matrix 1")
    # print( selected_sample_ids)
    # print(design_matrix)
    # print( dim(design_matrix))

    group_A <- test_pairs[i, ]$A
    group_B <- test_pairs[i, ]$B

    x <- runTest(ID, A, B, group_A, group_B,
      design_matrix = design_matrix_subset,
      formula_string = formula_string, contrast_variable = contrast_variable,
      weights = subset_weights
    )

    comparison <- paste(group_B, " vs ", group_A, sep = "")

    r[[comparison]] <- list(results = x$table, counts = t(cbind(A, B)), fit.eb = x$fit.eb)
  }
  r
}

# ----------------------------------------------------------------------------
# runTestsContrasts
# ----------------------------------------------------------------------------
#' Run the linear model fitting and statistical tests for a set of contrasts, then adjust with Empirical Bayes function
#' @param data Data frame containing the log (base 2) protein abundance values where each column represents a sample and each row represents a protein group, and proteins as rows. The data is preferably median-scaled, with missing values imputed, and batch-effects removed.
#' @param contrast_strings Input file with a table listing all the experimental contrasts to analyse. It will be in the format required for the function \code{makeContrasts} in the limma package.
#' The contrast string consists of variable that each consist of concatenating the column name (e.g. group) and the string representing the group type (e.g. A) in the design matrix.
#' @param design_matrix A data frame with a column containing the sample ID (as per the sample_id param) and the experimental group (as per the group param). Each row as the sample ID as row name in the data frame.
#' @param formula_string A formula string representing the experimental design. e.g. ("~ 0 + group")
#' @param p_value_column The name of the raw p-value column (tidyverse style).
#' @param q_value_column The name of the q-value column (tidyverse style).
#' @param fdr_value_column The name of the fdr-value column (tidyverse style).
#' @return A list containing two elements. $results returns a list of tables containing logFC and q-values. $fit.eb returns the Empiracle Bayes output object.
#' @export
runTestsContrasts <- function(data,
                              contrast_strings,
                              design_matrix,
                              formula_string,
                              p_value_column = raw_pvalue,
                              q_value_column = fdr_qvalue,
                              fdr_value_column = fdr_value_bh_adjustment,
                              weights = NA,
                              treat_lfc_cutoff = NA,
                              eBayes_trend = FALSE,
                              eBayes_robust = FALSE) {
  message("--- Entering runTestsContrasts ---")
  message(sprintf("   runTestsContrasts: data dims = %d x %d, %d contrasts", nrow(data), ncol(data), length(contrast_strings)))
  message(sprintf("   runTestsContrasts: contrasts = %s", paste(contrast_strings, collapse = ", ")))
  message(sprintf("   runTestsContrasts: treat_lfc_cutoff = %s", treat_lfc_cutoff))

  # Create formula and design matrix
  ff <- as.formula(formula_string)
  mod_frame <- model.frame(ff, design_matrix)
  design_m <- model.matrix(ff, mod_frame)
  message(sprintf("   runTestsContrasts: design_m dims = %d x %d", nrow(design_m), ncol(design_m)))

  # Subset data to match design matrix
  data_subset <- data[, rownames(design_m)]
  message(sprintf("   runTestsContrasts: data_subset dims = %d x %d", nrow(data_subset), ncol(data_subset)))

  # Create contrast matrix
  message("   runTestsContrasts: Creating contrast matrix...")
  contr.matrix <- makeContrasts(
    contrasts = contrast_strings,
    levels = colnames(design_m)
  )
  message(sprintf("   runTestsContrasts: contr.matrix dims = %d x %d", nrow(contr.matrix), ncol(contr.matrix)))

  # Check weights
  if (!is.na(weights)) {
    message("   runTestsContrasts: Attaching weights...")
    if (nrow(weights) == nrow(design_m)) {
      design_m <- cbind(design_m, weights)
    } else {
      stop("Stop: nrow(weights) should be equal to nrow(design_m)")
    }
  }

  # Run limma analysis
  message("   runTestsContrasts: Running lmFit...")

  # Check for technical replicates
  # Blocking factor should be constructed from replicate ID (biological replicate)
  # But we need to ensure it maps correctly to the samples in the subset
  # The design_matrix here is already subsetted/ordered match mod_frame in model.matrix construction
  # However, we need the original columns (group, replicates) which might be in the 'design_matrix' argument (which is the params data.frame)
  # BUT 'design_matrix' input argument is strictly the data.frame with metadata

  # Ensure we have the metadata for the subsetted samples
  # The samples in data_subset are columns matching rownames(design_m)
  # design_m was created from design_matrix (the metadata dataframe)
  samples_in_model <- rownames(design_m)

  # Check if we have replicates column
  has_replicates <- "replicates" %in% colnames(design_matrix)

  fit <- NULL

  if (has_replicates) {
    # Extract replicates for the samples in the model, maintaining order
    # We assume design_matrix rownames are the sample IDs (which was set in helper/wrapper functions)
    design_matrix_subset <- design_matrix[samples_in_model, ]

    # Construct blocking factor: Biological Replicate ID (e.g. "Control_1")
    # This identifies the biological unit that is technically replicated
    # Combining group + replicate number gives unique biological ID
    # Use paste0 to ensure character vector
    block <- paste(design_matrix_subset$group, design_matrix_subset$replicates, sep = "_")

    # Check if there are any technical replicates (duplicated blocks)
    if (any(duplicated(block))) {
      message("   runTestsContrasts: Detected technical replicates. Calculating duplicateCorrelation...")
      message(sprintf(
        "   runTestsContrasts: Block defined by group_replicates. %d unique blocks for %d samples.",
        length(unique(block)), length(block)
      ))

      # Calculate consensus correlation
      # Note: duplicateCorrelation can be slow for large datasets
      dup_cor <- duplicateCorrelation(data_subset, design = design_m, block = block)

      message(sprintf("   runTestsContrasts: Consensus correlation = %.4f", dup_cor$consensus.correlation))

      # Run lmFit with correlation and block
      message("   runTestsContrasts: Running lmFit with duplicateCorrelation...")
      fit <- lmFit(data_subset, design = design_m, block = block, correlation = dup_cor$consensus.correlation)
    } else {
      message("   runTestsContrasts: No technical replicates detected (unique blocks). Running standard lmFit...")
      fit <- lmFit(data_subset, design = design_m)
    }
  } else {
    message("   runTestsContrasts: 'replicates' column not found. Running standard lmFit...")
    fit <- lmFit(data_subset, design = design_m)
  }

  message("   runTestsContrasts: Running contrasts.fit...")
  cfit <- contrasts.fit(fit, contrasts = contr.matrix)

  message("   runTestsContrasts: Running eBayes...")
  eb.fit <- eBayes(cfit, trend = eBayes_trend, robust = eBayes_robust)

  # Run treat or standard analysis
  t.fit <- NA
  result_tables <- NA
  if (!is.na(treat_lfc_cutoff)) {
    message("   runTestsContrasts: Running treat analysis...")
    t.fit <- treat(eb.fit, lfc = as.double(treat_lfc_cutoff))

    message(sprintf("   runTestsContrasts: Processing %d contrasts with topTreat...", length(contrast_strings)))
    # Track qvalue failures for user notification
    qvalue_failures <- list()

    result_tables <- purrr::map(
      contrast_strings,
      function(contrast) {
        message(sprintf("      [map] Processing contrast: %s", contrast))
        qvalue_failed <- FALSE

        tryCatch(
          {
            message(sprintf("      About to call topTreat with coef = %s", contrast))
            da_tbl <- topTreat(t.fit, coef = contrast, n = Inf)
            message(sprintf("      [map] topTreat success: %d rows", nrow(da_tbl)))

            message("      Adding qvalue column...")
            # Safe qvalue computation: handle invalid p-values (NA, Inf, NaN)
            valid_p_idx <- which(!is.na(da_tbl$P.Value) & is.finite(da_tbl$P.Value))
            if (length(valid_p_idx) > 0) {
              # Compute q-values only for valid p-values
              valid_p_values <- da_tbl$P.Value[valid_p_idx]

              # Diagnostic: Log p-value distribution statistics
              message(sprintf("      Diagnostic: Valid p-values: %d of %d total", length(valid_p_idx), nrow(da_tbl)))
              message(sprintf("      Diagnostic: P-value range: [%.6f, %.6f]", min(valid_p_values), max(valid_p_values)))
              message(sprintf("      Diagnostic: P-value mean: %.6f, median: %.6f", mean(valid_p_values), median(valid_p_values)))

              # Edge case checks that might cause qvalue() to fail
              all_zeros <- all(valid_p_values == 0)
              all_ones <- all(valid_p_values == 1)
              too_few <- length(valid_p_values) < 3

              if (all_zeros) {
                message("      Warning: All p-values are 0 - qvalue() cannot compute, using p.adjust()")
                use_qvalue <- FALSE
              } else if (all_ones) {
                message("      Warning: All p-values are 1 - qvalue() may fail, using p.adjust()")
                use_qvalue <- FALSE
              } else if (too_few) {
                message(sprintf("      Warning: Too few p-values (%d < 3) for qvalue() estimation, using p.adjust()", length(valid_p_values)))
                use_qvalue <- FALSE
              } else {
                use_qvalue <- TRUE
              }

              q_values_all <- rep(NA_real_, nrow(da_tbl))
              if (use_qvalue) {
                tryCatch(
                  {
                    q_values_valid <- qvalue(valid_p_values)$q
                    q_values_all[valid_p_idx] <- q_values_valid
                    message("      qvalue() computation successful")
                  },
                  error = function(e) {
                    qvalue_failed <<- TRUE
                    message(sprintf("      Warning: qvalue() failed during computation: %s", e$message))
                    message(sprintf("      Diagnostic: P-value distribution may be problematic for qvalue smoothing algorithm"))
                    message(sprintf("      Diagnostic: Falling back to p.adjust() method='BH' (Benjamini-Hochberg FDR)"))
                    # Fallback to p.adjust if qvalue fails
                    q_values_all[valid_p_idx] <- p.adjust(valid_p_values, method = "BH")
                    message(sprintf("      Diagnostic: Assigned %d p.adjust() values to q-value column", length(valid_p_idx)))
                  }
                )
              } else {
                # Use p.adjust directly for edge cases
                qvalue_failed <<- TRUE
                q_values_all[valid_p_idx] <- p.adjust(valid_p_values, method = "BH")
                message("      Using p.adjust() due to edge case detection")
                message(sprintf("      Diagnostic: Assigned %d p.adjust() values to q-value column", length(valid_p_idx)))
              }
            } else {
              # All p-values are invalid, set all q-values to NA
              q_values_all <- rep(NA_real_, nrow(da_tbl))
              message("      Warning: All p-values are invalid (NA, Inf, or NaN), setting q-values to NA")
            }

            # Debug: Verify assignment before mutate
            message(sprintf("      Diagnostic: q_values_all has %d non-NA values before assignment", sum(!is.na(q_values_all))))

            da_tbl <- da_tbl |>
              mutate({{ q_value_column }} := q_values_all)

            # Debug: Verify assignment after mutate
            assigned_col <- da_tbl[[rlang::as_name(rlang::ensym(q_value_column))]]
            message(sprintf("      Diagnostic: Assigned column has %d non-NA values after mutate", sum(!is.na(assigned_col))))
            message("      qvalue column added")

            message("      Adding FDR column...")
            # Use the same safe logic for FDR column - only compute for valid p-values
            fdr_values_all <- rep(NA_real_, nrow(da_tbl))
            if (length(valid_p_idx) > 0) {
              fdr_values_all[valid_p_idx] <- p.adjust(valid_p_values, method = "BH")
            }
            da_tbl <- da_tbl |>
              mutate({{ fdr_value_column }} := fdr_values_all)
            message("      FDR column added")

            message("      Renaming P.Value column...")
            da_tbl <- da_tbl |>
              dplyr::rename({{ p_value_column }} := P.Value)
            message("      P.Value column renamed")

            message(sprintf("   [map] Completed processing contrast: %s", contrast))
            if (qvalue_failed) {
              qvalue_failures[[contrast]] <<- TRUE
            }
            return(da_tbl)
          },
          error = function(e) {
            message(sprintf("      [map] ERROR in contrast %s: %s", contrast, e$message))
            message(sprintf("      [map] ERROR call stack: %s", capture.output(traceback())))
            stop(e)
          }
        )
      }
    )
  } else {
    message("   runTestsContrasts: Running standard analysis...")
    t.fit <- eb.fit

    message(sprintf("   runTestsContrasts: Processing %d contrasts with topTable...", length(contrast_strings)))
    # Track qvalue failures for user notification
    qvalue_failures <- list()

    result_tables <- purrr::map(
      contrast_strings,
      function(contrast) {
        message(sprintf("      [map] Processing contrast: %s", contrast))
        qvalue_failed <- FALSE

        tryCatch(
          {
            message(sprintf("      About to call topTable with coef = %s", contrast))
            da_tbl <- topTable(t.fit, coef = contrast, n = Inf)
            message(sprintf("      [map] topTable success: %d rows", nrow(da_tbl)))

            message("      Adding qvalue column...")
            # Safe qvalue computation: handle invalid p-values (NA, Inf, NaN)
            valid_p_idx <- which(!is.na(da_tbl$P.Value) & is.finite(da_tbl$P.Value))
            if (length(valid_p_idx) > 0) {
              # Compute q-values only for valid p-values
              valid_p_values <- da_tbl$P.Value[valid_p_idx]

              # Diagnostic: Log p-value distribution statistics
              message(sprintf("      Diagnostic: Valid p-values: %d of %d total", length(valid_p_idx), nrow(da_tbl)))
              message(sprintf("      Diagnostic: P-value range: [%.6f, %.6f]", min(valid_p_values), max(valid_p_values)))
              message(sprintf("      Diagnostic: P-value mean: %.6f, median: %.6f", mean(valid_p_values), median(valid_p_values)))

              # Edge case checks that might cause qvalue() to fail
              all_zeros <- all(valid_p_values == 0)
              all_ones <- all(valid_p_values == 1)
              too_few <- length(valid_p_values) < 3

              if (all_zeros) {
                message("      Warning: All p-values are 0 - qvalue() cannot compute, using p.adjust()")
                use_qvalue <- FALSE
              } else if (all_ones) {
                message("      Warning: All p-values are 1 - qvalue() may fail, using p.adjust()")
                use_qvalue <- FALSE
              } else if (too_few) {
                message(sprintf("      Warning: Too few p-values (%d < 3) for qvalue() estimation, using p.adjust()", length(valid_p_values)))
                use_qvalue <- FALSE
              } else {
                use_qvalue <- TRUE
              }

              q_values_all <- rep(NA_real_, nrow(da_tbl))
              if (use_qvalue) {
                tryCatch(
                  {
                    q_values_valid <- qvalue(valid_p_values)$q
                    q_values_all[valid_p_idx] <- q_values_valid
                    message("      qvalue() computation successful")
                  },
                  error = function(e) {
                    qvalue_failed <<- TRUE
                    message(sprintf("      Warning: qvalue() failed during computation: %s", e$message))
                    message(sprintf("      Diagnostic: P-value distribution may be problematic for qvalue smoothing algorithm"))
                    message(sprintf("      Diagnostic: Falling back to p.adjust() method='BH' (Benjamini-Hochberg FDR)"))
                    # Fallback to p.adjust if qvalue fails
                    q_values_all[valid_p_idx] <- p.adjust(valid_p_values, method = "BH")
                    message(sprintf("      Diagnostic: Assigned %d p.adjust() values to q-value column", length(valid_p_idx)))
                  }
                )
              } else {
                # Use p.adjust directly for edge cases
                qvalue_failed <<- TRUE
                q_values_all[valid_p_idx] <- p.adjust(valid_p_values, method = "BH")
                message("      Using p.adjust() due to edge case detection")
                message(sprintf("      Diagnostic: Assigned %d p.adjust() values to q-value column", length(valid_p_idx)))
              }
            } else {
              # All p-values are invalid, set all q-values to NA
              q_values_all <- rep(NA_real_, nrow(da_tbl))
              message("      Warning: All p-values are invalid (NA, Inf, or NaN), setting q-values to NA")
            }

            # Debug: Verify assignment before mutate
            message(sprintf("      Diagnostic: q_values_all has %d non-NA values before assignment", sum(!is.na(q_values_all))))

            da_tbl <- da_tbl |>
              mutate({{ q_value_column }} := q_values_all)

            # Debug: Verify assignment after mutate
            assigned_col <- da_tbl[[rlang::as_name(rlang::ensym(q_value_column))]]
            message(sprintf("      Diagnostic: Assigned column has %d non-NA values after mutate", sum(!is.na(assigned_col))))
            message("      qvalue column added")

            message("      Adding FDR column...")
            # Use the same safe logic for FDR column - only compute for valid p-values
            fdr_values_all <- rep(NA_real_, nrow(da_tbl))
            if (length(valid_p_idx) > 0) {
              fdr_values_all[valid_p_idx] <- p.adjust(valid_p_values, method = "BH")
            }
            da_tbl <- da_tbl |>
              mutate({{ fdr_value_column }} := fdr_values_all)
            message("      FDR column added")

            message("      Renaming P.Value column...")
            da_tbl <- da_tbl |>
              dplyr::rename({{ p_value_column }} := P.Value)
            message("      P.Value column renamed")

            message(sprintf("   [map] Completed processing contrast: %s", contrast))
            if (qvalue_failed) {
              qvalue_failures[[contrast]] <<- TRUE
            }
            return(da_tbl)
          },
          error = function(e) {
            message(sprintf("      [map] ERROR in contrast %s: %s", contrast, e$message))
            message(sprintf("      [map] ERROR call stack: %s", capture.output(traceback())))
            stop(e)
          }
        )
      }
    )
  }

  names(result_tables) <- contrast_strings
  message("--- Exiting runTestsContrasts ---")
  return(list(results = result_tables, fit.eb = t.fit, qvalue_warnings = qvalue_failures))
}

