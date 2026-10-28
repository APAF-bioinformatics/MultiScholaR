# ----------------------------------------------------------------------------
# ebFit
# ----------------------------------------------------------------------------
#' Run the Empircal Bayes Statistics for Differential Expression in the limma package
#' @param design Output from running the function \code{\link{model.matrix}}.
#' @param contr.matrix Output from the function \code{\link{makeContrasts}}.
#' @seealso \code{\link{model.matrix}}
#' @seealso \code{\link{makeContrasts}}
#' @export
#' @param data Runtime inputs used by this function; see the usage section for accepted values.
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
#' comparison A string showing `log(group B) - log(group A)`.
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

  if (!is.null(weights) && !(length(weights) == 1L && is.na(weights))) { # KIT-064 (roar BUG-05): a matrix has no scalar NA
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
#' @param sample_rows_list A list, the name of each element is the sample ID and each element is a vector containing the protein accessions (e.g. row_id) with enough number of values. It is usually the output from [getRowsToKeepList()].
#' @param type_of_grouping A list where each element name is the name of a treatment group and each element is a vector containing the sample IDs within the treatment group. It is usually the output from [getTypeOfGrouping()].
#' @param design_matrix A data frame with a column containing the sample ID (as per the sample_id param) and the experimental group (as per the group param). Each row as the sample ID as row name in the data frame.
#' @param formula_string A formula string representing the experimental design. e.g. ("~ 0 + group")
#' @param contrast_variable String representing the contrast variable, which is also used in the formula string. (e.g. "group")
#' @param weights Numeric matrix for adjusting each sample and gene.
#' @return A list of data frames, the name of each element represents each pairwise comparison. Each data frame has the following columns:
#' row.names = the protein accessions
#' comparison A string showing `log(group B) - log(group A)`.
#' meanA     mean of the normalised log abundance value of the gene across samples from experimental group A
#' meanB     mean of the normalised log abundance value of the gene across samples from experimental group B
#' logFC     log fold-change
#' tstats    t-test statistics
#' tmod      moderated t-test statistics
#' pval      t-test p-value
#' raw_pvalue      moderated t-test p-value
#' qval      t-test q-value
#' fdr_qvalue     moderated t-test q-value
#' @seealso [getRowsToKeepList()]
#' @seealso [getTypeOfGrouping()]
#' @export
runTests <- function(ID, data, test_pairs, sample_columns, sample_rows_list = NA, type_of_grouping, design_matrix, formula_string, contrast_variable = "group", weights = NA) {
  r <- list()
  for (i in seq_len(nrow(test_pairs))) { # KIT-064 (roar BUG-16)
    rows_to_keep <- rownames(data)


    if (length(sample_rows_list) > 0) {
      # KIT-064 (roar BUG-04/05): a scalar guard (the former `!is.na(list)` had the list's length) and the union of the two
      # groups' rows (the former `unique(a, b)` passed the second group as `incomparables`).
      if (!(is.atomic(sample_rows_list) && length(sample_rows_list) == 1L && is.na(sample_rows_list)) &&
        #  Check that sample group exists as names inside sample_rows_list
        length(which(c(test_pairs[i, "A"], test_pairs[i, "B"]) %in% names(sample_rows_list))) > 0) {
        rows_to_keep <- union(
          sample_rows_list[[test_pairs[[i, "A"]]]],
          sample_rows_list[[test_pairs[[i, "B"]]]]
        )
      }
    }

    tmp <- data[rows_to_keep, sample_columns]
    rep <- colnames(tmp)

    # print( paste( test_pairs[i,]$A, test_pairs[i,]$B) )
    A <- tmp[, type_of_grouping[test_pairs[i, ]$A][[1]], drop = FALSE] # KIT-064 (roar BUG-14)
    B <- tmp[, type_of_grouping[test_pairs[i, ]$B][[1]], drop = FALSE]

    subset_weights <- NA

    if (!is.null(weights) && !(length(weights) == 1L && is.na(weights))) { # KIT-064 (roar BUG-05)
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
#' Fit the existing uncontrasted limma model and retain its effective provenance.
#' @param data_subset Ordered feature-by-sample matrix.
#' @param design_m Existing model matrix, including historical weight handling.
#' @param design_matrix Original ordered sample metadata.
#' @param samples_in_model Ordered model samples.
#' @param duplicateCorrelationFn,lmFitFn Existing limma operations.
#' @param counterFn Optional operation counter used only by standalone proof.
#' @param messageFn Existing diagnostic sink.
#' @return Explicit uncontrasted fit artifact; never the public moderated fit.eb.
#' @noRd
fitProtDaBaseModel <- function(data_subset, design_m, design_matrix,
    samples_in_model, duplicateCorrelationFn = limma::duplicateCorrelation,
    lmFitFn = limma::lmFit, counterFn = \(operation) invisible(NULL),
    messageFn = message, group_id = "group", technical_replicate_id = "replicates") {
  # KIT-035: the block comes from the object's declared technical replicate column (the design step sets
  # `tech_rep_group`; the class default is `replicates`), through the rule RUV-III shares. A block spanned by the
  # design is standard lmFit with the reason recorded, never duplicateCorrelation with a zeroed correlation.
  design_subset <- design_matrix[samples_in_model, , drop = FALSE]
  replicate <- resolveProtTechnicalReplicateBlock(design_subset, group_id, technical_replicate_id)
  block <- replicate$block
  duplicate_correlation <- NULL
  uses_duplicate_correlation <- isTRUE(replicate$technical)
  if (uses_duplicate_correlation) {
    messageFn("   runTestsContrasts: Detected technical replicates. Calculating duplicateCorrelation...")
    messageFn(sprintf(
      "   runTestsContrasts: Block defined by %s. %d unique blocks for %d samples.",
      technical_replicate_id, length(unique(block)), length(block)
    ))
    counterFn("duplicateCorrelation")
    duplicate_correlation <- duplicateCorrelationFn(data_subset,
      design = design_m, block = block)
    messageFn(sprintf("   runTestsContrasts: Consensus correlation = %.4f",
      duplicate_correlation$consensus.correlation))
    messageFn("   runTestsContrasts: Running lmFit with duplicateCorrelation...")
    counterFn("lmFit")
    fit <- lmFitFn(data_subset, design = design_m, block = block,
      correlation = duplicate_correlation$consensus.correlation)
  } else {
    messageFn(switch(replicate$reason,
      no_declaration = sprintf("   runTestsContrasts: technical replicate column '%s' not found. Running standard lmFit...", technical_replicate_id),
      spanned_by_design = "   runTestsContrasts: Replicate block is spanned by the design (no technical replicates). Running standard lmFit...",
      missing_values = "   runTestsContrasts: Replicate column has missing values. Running standard lmFit...",
      "   runTestsContrasts: No technical replicates detected (unique blocks). Running standard lmFit..."))
    counterFn("lmFit")
    fit <- lmFitFn(data_subset, design = design_m)
  }
  list(schema = "prot_da_base_fit_v1", fit = fit,
    sample_ids = colnames(data_subset), feature_ids = rownames(data_subset),
    design_columns = colnames(design_m), block = block,
    group_id = group_id, technical_replicate_id = technical_replicate_id, block_reason = replicate$reason,
    uses_duplicate_correlation = uses_duplicate_correlation,
    consensus_correlation = if (is.null(duplicate_correlation)) NA_real_ else
      duplicate_correlation$consensus.correlation)
}

#' Apply the existing contrast, moderation and optional TREAT sequence.
#' @param base_fit_artifact Artifact returned by fitProtDaBaseModel.
#' @param contrast_matrix Existing limma contrast matrix.
#' @param treat_lfc_cutoff Existing TREAT cutoff or NA.
#' @param eBayes_trend,eBayes_robust Existing moderation settings.
#' @param contrastsFitFn,eBayesFn,treatFn Existing limma operations.
#' @param counterFn Optional standalone operation counter.
#' @param messageFn Existing diagnostic sink.
#' @return Explicit contrasted/moderated state and unchanged public fit.eb value.
#' @noRd
applyProtDaContrastModel <- function(base_fit_artifact, contrast_matrix,
    treat_lfc_cutoff = NA, eBayes_trend = FALSE, eBayes_robust = FALSE,
    contrastsFitFn = limma::contrasts.fit, eBayesFn = limma::eBayes,
    treatFn = limma::treat, counterFn = \(operation) invisible(NULL),
    messageFn = message) {
  if (!is.list(base_fit_artifact) ||
      !identical(base_fit_artifact$schema, "prot_da_base_fit_v1") ||
      is.null(base_fit_artifact$fit) || !is.matrix(contrast_matrix)) {
    cli::cli_abort("invalid DA base-fit artifact", class = "saved_entity_invalid")
  }
  messageFn("   runTestsContrasts: Running contrasts.fit...")
  counterFn("contrasts.fit")
  contrasted_fit <- contrastsFitFn(base_fit_artifact$fit, contrasts = contrast_matrix)
  messageFn("   runTestsContrasts: Running eBayes...")
  counterFn("eBayes")
  moderated_fit <- eBayesFn(contrasted_fit, trend = eBayes_trend,
    robust = eBayes_robust)
  treat_applied <- !is.na(treat_lfc_cutoff)
  public_fit <- moderated_fit
  if (treat_applied) {
    messageFn("   runTestsContrasts: Running treat analysis...")
    counterFn("treat")
    public_fit <- treatFn(moderated_fit, lfc = as.double(treat_lfc_cutoff))
  }
  list(schema = "prot_da_contrast_model_v1",
    base_fit_entity = NA_character_, contrast_matrix = contrast_matrix,
    contrasted_fit = contrasted_fit, moderated_fit = moderated_fit,
    public_fit = public_fit, treat_applied = treat_applied,
    treat_lfc_cutoff = treat_lfc_cutoff,
    eBayes_trend = eBayes_trend, eBayes_robust = eBayes_robust)
}

#' Apply the existing qvalue fallback and BH populations to one full table.
#' @param table Full topTable/topTreat output.
#' @param p_value_name,q_value_name,fdr_value_name Output column names.
#' @param qvalueFn,pAdjustFn Existing adjustment functions.
#' @param counterFn Optional standalone operation counter.
#' @param messageFn Existing diagnostic sink.
#' @return Adjusted table and whether qvalue used its historical fallback.
#' @noRd
adjustProtDaResultTable <- function(table, p_value_name, q_value_name,
    fdr_value_name, qvalueFn = qvalue::qvalue, pAdjustFn = stats::p.adjust,
    counterFn = \(operation) invisible(NULL), messageFn = message) {
  if (!is.data.frame(table) || !"P.Value" %in% names(table)) {
    cli::cli_abort("invalid limma result table", class = "saved_entity_invalid")
  }
  for (name in c(p_value_name, q_value_name, fdr_value_name)) {
    if (!is.character(name) || length(name) != 1L || is.na(name) || !nzchar(name)) {
      cli::cli_abort("invalid DA result column", class = "saved_entity_invalid")
    }
  }
  valid <- which(!is.na(table$P.Value) & is.finite(table$P.Value))
  p_values <- table$P.Value[valid]
  q_values <- rep(NA_real_, nrow(table))
  qvalue_failed <- FALSE
  use_qvalue <- length(p_values) >= 3L && !all(p_values == 0) && !all(p_values == 1)
  if (length(p_values) && use_qvalue) {
    q_values[valid] <- tryCatch({
      counterFn("qvalue")
      qvalueFn(p_values)$q
    }, error = function(error) {
      qvalue_failed <<- TRUE
      messageFn(sprintf("      Warning: qvalue() failed during computation: %s",
        conditionMessage(error)))
      counterFn("p.adjust:qvalue_fallback")
      pAdjustFn(p_values, method = "BH")
    })
  } else if (length(p_values)) {
    qvalue_failed <- TRUE
    counterFn("p.adjust:qvalue_fallback")
    q_values[valid] <- pAdjustFn(p_values, method = "BH")
  }
  fdr_values <- rep(NA_real_, nrow(table))
  if (length(p_values)) {
    counterFn("p.adjust:fdr")
    fdr_values[valid] <- pAdjustFn(p_values, method = "BH")
  }
  table[[q_value_name]] <- q_values
  table[[fdr_value_name]] <- fdr_values
  names(table)[names(table) == "P.Value"] <- p_value_name
  list(table = table, qvalue_failed = qvalue_failed,
    valid_p_value_mask = seq_len(nrow(table)) %in% valid)
}

#' Build full-population contrast result tables from a reusable contrast model.
#' @param contrast_artifact Artifact returned by applyProtDaContrastModel.
#' @param contrast_strings Ordered contrast definitions.
#' @param p_value_name,q_value_name,fdr_value_name Output column names.
#' @param topTableFn,topTreatFn,qvalueFn,pAdjustFn Existing operations.
#' @param counterFn Optional standalone operation counter.
#' @param messageFn Existing diagnostic sink.
#' @return Named full result tables, warnings and finite-p-value masks.
#' @noRd
buildProtDaResultTables <- function(contrast_artifact, contrast_strings,
    p_value_name = "raw_pvalue", q_value_name = "fdr_qvalue",
    fdr_value_name = "fdr_value_bh_adjustment",
    topTableFn = limma::topTable, topTreatFn = limma::topTreat,
    qvalueFn = qvalue::qvalue, pAdjustFn = stats::p.adjust,
    counterFn = \(operation) invisible(NULL), messageFn = message) {
  if (!is.list(contrast_artifact) ||
      !identical(contrast_artifact$schema, "prot_da_contrast_model_v1") ||
      is.null(contrast_artifact$public_fit) || !is.character(contrast_strings) ||
      !length(contrast_strings) || anyNA(contrast_strings)) {
    cli::cli_abort("invalid DA contrast artifact", class = "saved_entity_invalid")
  }
  results <- vector("list", length(contrast_strings))
  warnings <- list()
  masks <- vector("list", length(contrast_strings))
  for (index in seq_along(contrast_strings)) {
    contrast <- contrast_strings[[index]]
    if (isTRUE(contrast_artifact$treat_applied)) {
      counterFn("topTreat")
      table <- topTreatFn(contrast_artifact$public_fit, coef = contrast, n = Inf)
    } else {
      counterFn("topTable")
      table <- topTableFn(contrast_artifact$public_fit, coef = contrast, number = Inf)
    }
    adjusted <- adjustProtDaResultTable(table, p_value_name, q_value_name,
      fdr_value_name, qvalueFn, pAdjustFn, counterFn, messageFn)
    results[[index]] <- adjusted$table
    masks[[index]] <- adjusted$valid_p_value_mask
    if (isTRUE(adjusted$qvalue_failed)) warnings[[contrast]] <- TRUE
  }
  names(results) <- contrast_strings
  names(masks) <- contrast_strings
  list(schema = "prot_da_result_tables_v1", results = results,
    qvalue_warnings = warnings, valid_p_value_masks = masks,
    fit.eb = contrast_artifact$public_fit)
}

#' @param p_value_column The name of the raw p-value column (tidyverse style).
#' @param q_value_column The name of the q-value column (tidyverse style).
#' @param fdr_value_column The name of the fdr-value column (tidyverse style).
#' @return A list containing two elements. $results returns a list of tables containing logFC and q-values. $fit.eb returns the Empiracle Bayes output object.
#' @export
#' @param weights,treat_lfc_cutoff,eBayes_trend,eBayes_robust Runtime inputs used by this function; see the usage section for accepted values.
#' @param group_id,technical_replicate_id Column names of the design's group and of the object's declared technical replicate column (KIT-035): the replicate block comes from the declared column, and a block spanned by the design runs standard lmFit with the reason recorded.
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
                              eBayes_robust = FALSE,
                              group_id = "group",
                              technical_replicate_id = "replicates") {
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
  data_subset <- data[, rownames(design_m), drop = FALSE] # KIT-064 (roar BUG-14)
  message(sprintf("   runTestsContrasts: data_subset dims = %d x %d", nrow(data_subset), ncol(data_subset)))

  # Create contrast matrix
  message("   runTestsContrasts: Creating contrast matrix...")
  contr.matrix <- makeContrasts(
    contrasts = contrast_strings,
    levels = colnames(design_m)
  )
  message(sprintf("   runTestsContrasts: contr.matrix dims = %d x %d", nrow(contr.matrix), ncol(contr.matrix)))

  # Check weights
  if (!is.null(weights) && !(length(weights) == 1L && is.na(weights))) { # KIT-064 (roar BUG-05): a matrix has no scalar NA
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

  base_fit_artifact <- fitDIAExecutionBaseModel(data_subset, design_m,
    design_matrix, samples_in_model, formula_string, weights,
    group_id = group_id, technical_replicate_id = technical_replicate_id)

  saved_manager <- getDIAActiveManager()
  saved_identity <- NULL
  if (usesDIAReferenceStates(saved_manager) && isTRUE(saved_manager$.context$da_reference_execution)) {
    saved_identity <- newWorkflowDiaResultIdentity(saved_manager$.context$da_base_fit_reference,
      contrast_strings, treat_lfc_cutoff, eBayes_trend, eBayes_robust,
      paste("full-population-qvalue-BH-v1", rlang::as_name(rlang::ensym(p_value_column)),
        rlang::as_name(rlang::ensym(q_value_column)), rlang::as_name(rlang::ensym(fdr_value_column)), sep = ":"),
      saved_manager$.context$software)
    saved_result <- resolveWorkflowDiaResult(saved_manager$.context$store, saved_identity)
    if (identical(saved_result$status, "hit")) return(saved_result$value$value$public_result)
  }

  contrast_artifact <- applyProtDaContrastModel(base_fit_artifact, contr.matrix,
    treat_lfc_cutoff = treat_lfc_cutoff, eBayes_trend = eBayes_trend,
    eBayes_robust = eBayes_robust, contrastsFitFn = contrasts.fit,
    eBayesFn = eBayes, treatFn = treat)
  eb.fit <- contrast_artifact$moderated_fit

  # Run treat or standard analysis
  t.fit <- contrast_artifact$public_fit
  result_tables <- NA
  if (!is.na(treat_lfc_cutoff)) {
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
                # Use p.adjust directly for edge cases. KIT-064 addendum (roar BUG-21, the reviewer's correction): this
                # branch runs in the contrast function's own frame, where `qvalue_failed` is declared, so a `<<-` here
                # skipped it (the contrast's q-value warning was never recorded and a global leaked); `<-` it is.
                qvalue_failed <- TRUE
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
            da_tbl <- topTable(t.fit, coef = contrast, number = Inf)
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
                # Use p.adjust directly for edge cases. KIT-064 addendum (roar BUG-21, the reviewer's correction): this
                # branch runs in the contrast function's own frame, where `qvalue_failed` is declared, so a `<<-` here
                # skipped it (the contrast's q-value warning was never recorded and a global leaked); `<-` it is.
                qvalue_failed <- TRUE
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
  result <- list(results = result_tables, fit.eb = t.fit, qvalue_warnings = qvalue_failures)
  if (!is.null(saved_identity)) {
    saved_reference <- saveWorkflowDiaResult(saved_manager$.context$store,
      list(schema = "prot_da_result_tables_v1", results = result_tables, public_result = result),
      basename(tempfile("dia-contrasted-fit-")), saved_identity)
    rootDIAProduct(saved_manager, saved_reference, list(saved_manager$.context$da_base_fit_reference))
  }
  return(result)
}
