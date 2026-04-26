library(testthat)

# fidelity-coverage-compare: shared

plotPcaDispatch <- get(
  "plotPcaDispatch",
  envir = asNamespace("MultiScholaR"),
  inherits = FALSE
)

.peptide_calculateAdaptiveMaxK <- get(
  ".peptide_calculateAdaptiveMaxK",
  envir = asNamespace("MultiScholaR"),
  inherits = FALSE
)

.peptide_calculateSeparationScore <- get(
  ".peptide_calculateSeparationScore",
  envir = asNamespace("MultiScholaR"),
  inherits = FALSE
)

.peptide_calculateCompositeScore <- get(
  ".peptide_calculateCompositeScore",
  envir = asNamespace("MultiScholaR"),
  inherits = FALSE
)

newPeptideNormObject <- function(peptide_count = 6, sample_ids = paste0("S", 1:4)) {
  peptide_ids <- paste0("PEP", seq_len(peptide_count))
  protein_ids <- paste0("P", seq_len(peptide_count))
  peptide_data <- expand.grid(
    peptide_index = seq_len(peptide_count),
    Run = sample_ids,
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  peptide_data$Protein.Ids <- protein_ids[peptide_data$peptide_index]
  peptide_data$Stripped.Sequence <- peptide_ids[peptide_data$peptide_index]
  peptide_data$Precursor.Quantity <- 100 + peptide_data$peptide_index * 10 +
    match(peptide_data$Run, sample_ids)
  peptide_data$Precursor.Normalised <- 10 + peptide_data$peptide_index +
    match(peptide_data$Run, sample_ids) / 10
  peptide_data$Q.Value <- 0.001
  peptide_data <- peptide_data[
    c(
      "Protein.Ids",
      "Stripped.Sequence",
      "Run",
      "Precursor.Quantity",
      "Precursor.Normalised",
      "Q.Value"
    )
  ]

  peptide_object <- new(
    "PeptideQuantitativeData",
    peptide_data = peptide_data,
    design_matrix = data.frame(
      Run = sample_ids,
      group = rep(c("G1", "G2"), length.out = length(sample_ids)),
      batch = rep(c("B1", "B2"), length.out = length(sample_ids)),
      replicates = rep(c("R1", "R2"), length.out = length(sample_ids)),
      stringsAsFactors = FALSE
    ),
    sample_id = "Run",
    group_id = "group",
    technical_replicate_id = "replicates",
    protein_id_column = "Protein.Ids",
    peptide_sequence_column = "Stripped.Sequence",
    q_value_column = "Q.Value",
    raw_quantity_column = "Precursor.Quantity",
    norm_quantity_column = "Precursor.Normalised",
    args = list()
  )

  calcPeptideMatrix(peptide_object)
}

newCancorPlot <- function(delta = 0.25) {
  ggplot2::ggplot(
    data.frame(
      featureset = rep(c("Control", "All"), each = 3),
      cc = c(0.1, 0.2, 0.3, 0.1 + delta, 0.2 + delta * 2, 0.3 + delta * 3),
      K = rep(1:3, times = 2)
    ),
    ggplot2::aes(K, cc)
  ) +
    ggplot2::geom_point()
}

localGlobalBinding <- function(name, value, env = parent.frame()) {
  had_binding <- exists(name, envir = .GlobalEnv, inherits = FALSE)
  old_value <- if (had_binding) {
    get(name, envir = .GlobalEnv, inherits = FALSE)
  } else {
    NULL
  }

  assign(name, value, envir = .GlobalEnv)
  withr::defer({
    if (had_binding) {
      assign(name, old_value, envir = .GlobalEnv)
    } else if (exists(name, envir = .GlobalEnv, inherits = FALSE)) {
      rm(list = name, envir = .GlobalEnv)
    }
  }, envir = env)

  invisible(value)
}

multischolarPackageRoot <- function() {
  ns_path <- tryCatch(
    getNamespaceInfo(asNamespace("MultiScholaR"), "path"),
    error = function(...) ""
  )
  if (nzchar(ns_path) && dir.exists(ns_path)) {
    return(ns_path)
  }

  normalizePath(file.path("..", ".."), winslash = "/", mustWork = FALSE)
}

is_refactored_peptide_qc_methods <- function() {
  file.exists(file.path(multischolarPackageRoot(), "R", "func_pept_s4_qc_methods.R"))
}

allow_legacy_peptide_qc_error <- function(expr, expected_message) {
  result <- tryCatch(
    force(expr),
    error = function(err) err
  )
  if (inherits(result, "error")) {
    if (is_refactored_peptide_qc_methods()) {
      stop(conditionMessage(result), call. = FALSE)
    }
    expect_match(conditionMessage(result), expected_message)
    return(NULL)
  }

  result
}

callPeptideS4Method <- function(generic, theObject, ...) {
  method <- methods::selectMethod(generic, "PeptideQuantitativeData")
  method(theObject, ...)
}

withPeptidePearsonMethodBinding <- function(expr) {
  method <- methods::selectMethod("pearsonCorForSamplePairs", "PeptideQuantitativeData")
  local_mocked_bindings(
    pearsonCorForSamplePairs = function(theObject, ...) {
      method(theObject, ...)
    },
    .env = environment(methods::selectMethod("plotPearson", "PeptideQuantitativeData")@.Data)
  )
  force(expr)
}

test_that("PeptideQuantitativeData QC plotting methods cover RLE, PCA, density, and Pearson paths", {
  peptide_object <- newPeptideNormObject(peptide_count = 6, sample_ids = c("S1", "S2", "S3", "S4"))
  peptide_object@peptide_matrix <- matrix(
    c(
      10, 12, 14, 16,
      16, 14, 12, 10,
      10, 14, 16, 12,
      12, 16, 10, 14,
      14, 10, 16, 12,
      16, 12, 10, 14
    ),
    nrow = 6,
    byrow = TRUE,
    dimnames = dimnames(peptide_object@peptide_matrix)
  )
  rle_object <- peptide_object
  rle_object@peptide_matrix[1, 1] <- Inf

  rle_plot <- callPeptideS4Method(
    "plotRle",
    rle_object,
    grouping_variable = "group",
    yaxis_limit = c(-1, 1),
    sample_label = "Run"
  )
  expect_s3_class(rle_plot, "ggplot")

  pca_plot <- allow_legacy_peptide_qc_error(callPeptideS4Method(
    "plotPca",
    peptide_object,
    grouping_variable = "group",
    shape_variable = "batch",
    label_column = NA_character_,
    title = "Peptide PCA",
    font_size = 5,
    cv_percentile = 0
  ), "mixOmics")
  if (!is.null(pca_plot)) {
    expect_s3_class(pca_plot, "ggplot")
    expect_identical(pca_plot$labels$title, "Peptide PCA")
  }

  dispatch_plot <- allow_legacy_peptide_qc_error(plotPcaDispatch(
    peptide_object,
    grouping_variable = "group",
    shape_variable = NULL,
    label_column = "",
    title = "Peptide PCA dispatch",
    font_size = 6,
    cv_percentile = 0
  ), "mixOmics")
  if (!is.null(dispatch_plot)) {
    expect_s3_class(dispatch_plot, "ggplot")
    expect_identical(dispatch_plot$labels$title, "Peptide PCA dispatch")
  }

  density_plot <- callPeptideS4Method(
    "plotDensity",
    peptide_object,
    grouping_variable = "group",
    title = "Peptide density",
    font_size = 7
  )
  expect_s3_class(density_plot, "ggplot")
  expect_identical(density_plot$labels$title, "Peptide density")
  expect_identical(density_plot$labels$x, "Log2 Intensity")

  correlations <- callPeptideS4Method(
    "pearsonCorForSamplePairs",
    peptide_object,
    correlation_group = "replicates",
    exclude_pool_samples = FALSE
  )
  expect_named(correlations, c("Run.x", "Run.y", "pearson_correlation", "replicates"))
  expect_equal(nrow(correlations), 2)
  expect_true(all(is.finite(correlations$pearson_correlation)))

  pearson_plot <- withPeptidePearsonMethodBinding(callPeptideS4Method(
    "plotPearson",
    peptide_object,
    correlation_group = "replicates",
    exclude_pool_samples = FALSE
  ))
  expect_s3_class(pearson_plot, "ggplot")
  expect_identical(pearson_plot$labels$x, "Pearson Correlation")

  empty_group_object <- peptide_object
  empty_group_object@design_matrix$unique_sample <- empty_group_object@design_matrix$Run
  empty_plot <- withPeptidePearsonMethodBinding(callPeptideS4Method(
    "plotPearson",
    empty_group_object,
    correlation_group = "unique_sample",
    exclude_pool_samples = FALSE
  ))
  expect_s3_class(empty_plot, "ggplot")
  expect_identical(empty_plot$labels$title, "No within-group sample pairs found for correlation")
})

test_that("PeptideQuantitativeData QC plotting methods report invalid PCA inputs", {
  peptide_object <- newPeptideNormObject(peptide_count = 4, sample_ids = c("S1", "S2", "S3", "S4"))
  plot_pca_method <- methods::selectMethod("plotPca", "PeptideQuantitativeData")

  expect_error(
    plot_pca_method(
      peptide_object,
      grouping_variable = c("group", "batch"),
      label_column = "",
      title = "bad"
    ),
    "grouping_variable must be a single character string",
    fixed = TRUE
  )
  expect_error(
    plot_pca_method(
      peptide_object,
      grouping_variable = "missing_group",
      label_column = "",
      title = "bad"
    ),
    "grouping_variable 'missing_group' not found in design matrix",
    fixed = TRUE
  )
  expect_error(
    plot_pca_method(
      peptide_object,
      grouping_variable = "group",
      shape_variable = c("batch", "replicates"),
      label_column = "",
      title = "bad"
    ),
    "shape_variable must be NULL or a single character string",
    fixed = TRUE
  )
  expect_error(
    plot_pca_method(
      peptide_object,
      grouping_variable = "group",
      shape_variable = "missing_shape",
      label_column = "",
      title = "bad"
    ),
    "shape_variable 'missing_shape' not found in design matrix",
    fixed = TRUE
  )

  method_env <- environment(methods::selectMethod("plotPca", "PeptideQuantitativeData")@.Data)
  local_mocked_bindings(
    plotPcaHelper = function(...) {
      stop("helper exploded", call. = FALSE)
    },
    .env = method_env
  )

  expect_error(
    plot_pca_method(
      peptide_object,
      grouping_variable = "group",
      label_column = "",
      title = "bad"
    ),
    "Error in plotPcaHelper for PeptideQuantitativeData: helper exploded",
    fixed = TRUE
  )
})

test_that("PeptideQuantitativeData chooseBestProteinAccession and log2Transformation cover cleanup paths", {
  peptide_object <- newPeptideNormObject(peptide_count = 2, sample_ids = c("S1", "S2"))
  duplicate_rows <- peptide_object@peptide_data[
    peptide_object@peptide_data$Stripped.Sequence == "PEP1",
  ]
  duplicate_rows$Protein.Ids <- "P_CANON_ALT"
  duplicate_rows$Precursor.Quantity <- duplicate_rows$Precursor.Quantity + 100
  duplicate_rows$Precursor.Normalised <- duplicate_rows$Precursor.Normalised + 10
  peptide_object@peptide_data <- rbind(peptide_object@peptide_data, duplicate_rows)

  method_env <- environment(methods::selectMethod("chooseBestProteinAccession", "PeptideQuantitativeData")@.Data)
  local_mocked_bindings(
    chooseBestProteinAccessionHelper = function(input_tbl,
                                                acc_detail_tab,
                                                accessions_column,
                                                row_id_column,
                                                group_id,
                                                delim) {
      data.frame(
        Protein.Ids = unique(input_tbl$Protein.Ids),
        uniprot_acc = "P_CANON",
        stringsAsFactors = FALSE
      )
    },
    .env = method_env
  )

  accession_method <- methods::selectMethod("chooseBestProteinAccession", "PeptideQuantitativeData")
  peptide_object@args$chooseBestProteinAccession <- list()
  cleaned <- accession_method(
    peptide_object,
    delim = ";",
    seqinr_obj = list(aa_seq_tbl_final = data.frame(uniprot_acc = "P_CANON", stringsAsFactors = FALSE)),
    seqinr_accession_column = "uniprot_acc",
    aggregation_method = "mean"
  )

  expect_s4_class(cleaned, "PeptideQuantitativeData")
  expect_equal(unique(cleaned@peptide_data$Protein.Ids), "P_CANON")
  expect_true("P_CANON%PEP1" %in% rownames(cleaned@peptide_matrix))

  input_matrix <- matrix(c(0, 100, NA, 400), nrow = 2)
  transformed <- log2Transformation(input_matrix)
  expect_equal(transformed[1, 1], -Inf)
  expect_equal(transformed[2, 1], log2(101))
  expect_true(is.na(transformed[1, 2]))
  expect_equal(transformed[2, 2], log2(401))
})

test_that("PeptideQuantitativeData normalization preserves none mode and log2 transforms matrices", {
  peptide_object <- newPeptideNormObject(peptide_count = 3, sample_ids = c("S1", "S2"))
  original_matrix <- peptide_object@peptide_matrix

  normalized <- normaliseBetweenSamples(peptide_object, normalisation_method = "none")

  expect_equal(normalized@peptide_matrix, original_matrix)
  expect_equal(normalized@peptide_data$Precursor.Normalised, peptide_object@peptide_data$Precursor.Normalised)

  logged <- log2TransformPeptideMatrix(peptide_object)

  expect_true(logged@is_logged_data)
  expect_equal(logged@peptide_matrix, log2(original_matrix))
  expect_equal(
    logged@peptide_data$Precursor.Normalised[
      logged@peptide_data$Run == "S2" & logged@peptide_data$Stripped.Sequence == "PEP3"
    ],
    log2(original_matrix["P3%PEP3", "S2"])
  )
})

test_that("PeptideQuantitativeData log2 transform skips already logged objects", {
  peptide_object <- newPeptideNormObject(peptide_count = 2, sample_ids = c("S1", "S2"))
  peptide_object@is_logged_data <- TRUE

  logged <- suppressWarnings(log2TransformPeptideMatrix(peptide_object))

  expect_identical(logged@peptide_matrix, peptide_object@peptide_matrix)
  expect_true(logged@is_logged_data)
})

test_that("PeptideQuantitativeData negative-control selection forwards matrix and resolved parameters", {
  peptide_object <- newPeptideNormObject(peptide_count = 6, sample_ids = c("S1", "S2"))
  helper_call <- NULL
  method_env <- environment(methods::selectMethod("getNegCtrlProtAnovaPeptides", "PeptideQuantitativeData")@.Data)

  local_mocked_bindings(
    getNegCtrlProtAnovaHelper = function(data_matrix,
                                         design_matrix,
                                         grouping_variable,
                                         percentage_as_neg_ctrl,
                                         num_neg_ctrl,
                                         ruv_qval_cutoff,
                                         ruv_fdr_method) {
      helper_call <<- list(
        data_matrix = data_matrix,
        design_matrix = design_matrix,
        grouping_variable = grouping_variable,
        percentage_as_neg_ctrl = percentage_as_neg_ctrl,
        num_neg_ctrl = num_neg_ctrl,
        ruv_qval_cutoff = ruv_qval_cutoff,
        ruv_fdr_method = ruv_fdr_method
      )
      setNames(c(TRUE, TRUE, FALSE, FALSE, FALSE, FALSE), rownames(data_matrix))
    },
    .env = method_env
  )

  controls <- getNegCtrlProtAnovaPeptides(
    peptide_object,
    ruv_grouping_variable = "batch",
    percentage_as_neg_ctrl = 50,
    ruv_qval_cutoff = 0.2,
    ruv_fdr_method = "BH"
  )

  expect_identical(sum(controls), 2L)
  expect_equal(helper_call$data_matrix, peptide_object@peptide_matrix)
  expect_identical(helper_call$grouping_variable, "batch")
  expect_identical(helper_call$percentage_as_neg_ctrl, 50)
  expect_identical(helper_call$num_neg_ctrl, 3)
  expect_identical(helper_call$ruv_qval_cutoff, 0.2)
  expect_identical(helper_call$ruv_fdr_method, "BH")
})

test_that("PeptideQuantitativeData negative-control optimization returns per-percentage scores", {
  peptide_object <- newPeptideNormObject(peptide_count = 6, sample_ids = c("S1", "S2", "S3", "S4"))
  method_env <- environment(methods::selectMethod("findBestNegCtrlPercentagePeptides", "PeptideQuantitativeData")@.Data)

  local_mocked_bindings(
    getNegCtrlProtAnovaHelper = function(data_matrix,
                                         design_matrix,
                                         grouping_variable,
                                         percentage_as_neg_ctrl,
                                         num_neg_ctrl,
                                         ruv_qval_cutoff,
                                         ruv_fdr_method) {
      setNames(rep(TRUE, nrow(data_matrix)), rownames(data_matrix))
    },
    findBestK = function(cancorplot) {
      2
    },
    .env = method_env
  )

  results <- findBestNegCtrlPercentagePeptides(
    peptide_object,
    percentage_range = c(10, 20),
    num_components_to_impute = 3,
    ruv_grouping_variable = "batch",
    ruv_qval_cutoff = 0.2,
    ruv_fdr_method = "BH",
    separation_metric = "max_difference",
    k_penalty_weight = 0.25,
    max_acceptable_k = 3,
    adaptive_k_penalty = FALSE,
    verbose = FALSE
  )

  expect_true(results$best_percentage %in% c(10, 20))
  expect_identical(results$best_k, 2)
  expect_identical(results$separation_metric_used, "max_difference")
  expect_identical(results$adaptive_k_penalty_used, FALSE)
  expect_equal(results$optimization_results$percentage, c(10, 20))
  expect_identical(names(results$best_control_genes_index), rownames(peptide_object@peptide_matrix))
})

test_that("PeptideQuantitativeData negative-control optimization covers validation and fallback branches", {
  peptide_object <- newPeptideNormObject(peptide_count = 6, sample_ids = c("S1", "S2", "S3", "S4"))
  method_env <- environment(methods::selectMethod("findBestNegCtrlPercentagePeptides", "PeptideQuantitativeData")@.Data)

  expect_error(
    findBestNegCtrlPercentagePeptides(peptide_object, percentage_range = 0, verbose = FALSE),
    "between 0 and 100",
    fixed = TRUE
  )
  expect_error(
    findBestNegCtrlPercentagePeptides(peptide_object, separation_metric = "unsupported", verbose = FALSE),
    "separation_metric",
    fixed = TRUE
  )
  expect_error(
    findBestNegCtrlPercentagePeptides(peptide_object, k_penalty_weight = 2, verbose = FALSE),
    "between 0 and 1",
    fixed = TRUE
  )
  expect_error(
    findBestNegCtrlPercentagePeptides(peptide_object, max_acceptable_k = 0, verbose = FALSE),
    "positive number",
    fixed = TRUE
  )
  expect_error(
    findBestNegCtrlPercentagePeptides(peptide_object, adaptive_k_penalty = TRUE, max_acceptable_k = 2, verbose = FALSE),
    "at least 3",
    fixed = TRUE
  )

  local_mocked_bindings(
    getNegCtrlProtAnovaHelper = function(data_matrix,
                                         design_matrix,
                                         grouping_variable,
                                         percentage_as_neg_ctrl,
                                         num_neg_ctrl,
                                         ruv_qval_cutoff,
                                         ruv_fdr_method) {
      setNames(c(TRUE, TRUE, FALSE, FALSE, FALSE, FALSE), rownames(data_matrix))
    },
    .env = method_env
  )

  expect_error(
    findBestNegCtrlPercentagePeptides(
      peptide_object,
      percentage_range = c(10),
      ruv_grouping_variable = "batch",
      max_acceptable_k = 3,
      adaptive_k_penalty = TRUE,
      verbose = TRUE
    ),
    "No valid percentage",
    fixed = TRUE
  )

  local_mocked_bindings(
    getNegCtrlProtAnovaHelper = function(data_matrix,
                                         design_matrix,
                                         grouping_variable,
                                         percentage_as_neg_ctrl,
                                         num_neg_ctrl,
                                         ruv_qval_cutoff,
                                         ruv_fdr_method) {
      setNames(rep(TRUE, nrow(data_matrix)), rownames(data_matrix))
    },
    findBestK = function(cancorplot) {
      stop("best-k unavailable")
    },
    .env = method_env
  )

  expect_error(
    suppressWarnings(findBestNegCtrlPercentagePeptides(
      peptide_object,
      percentage_range = c(40),
      ruv_grouping_variable = "batch",
      adaptive_k_penalty = FALSE,
      verbose = TRUE
    )),
    "No valid percentage",
    fixed = TRUE
  )

  local_mocked_bindings(
    getNegCtrlProtAnovaHelper = function(data_matrix,
                                         design_matrix,
                                         grouping_variable,
                                         percentage_as_neg_ctrl,
                                         num_neg_ctrl,
                                         ruv_qval_cutoff,
                                         ruv_fdr_method) {
      setNames(rep(TRUE, nrow(data_matrix)), rownames(data_matrix))
    },
    findBestK = function(cancorplot) {
      1
    },
    .env = method_env
  )

  adaptive_results <- suppressWarnings(findBestNegCtrlPercentagePeptides(
    peptide_object,
    percentage_range = c(40),
    ruv_grouping_variable = "batch",
    max_acceptable_k = 3,
    adaptive_k_penalty = TRUE,
    verbose = TRUE
  ))

  expect_identical(adaptive_results$adaptive_k_penalty_used, TRUE)
  expect_identical(adaptive_results$max_acceptable_k, 2L)
  expect_identical(adaptive_results$sample_size, ncol(peptide_object@peptide_matrix))
})

test_that("PeptideQuantitativeData canonical-correlation methods validate and call RUV helpers", {
  peptide_object <- newPeptideNormObject(peptide_count = 6, sample_ids = c("S1", "S2", "S3", "S4"))
  controls <- setNames(rep(TRUE, 6), rownames(peptide_object@peptide_matrix))
  ruv_call <- NULL
  fast_plot <- suppressWarnings(ruvCancorFast(
    peptide_object,
    ctrl = controls,
    num_components_to_impute = 2,
    ruv_grouping_variable = "group",
    simple_imputation_method = "mean"
  ))

  expect_s3_class(fast_plot, "ggplot")
  expect_setequal(fast_plot$data$featureset, c("All", "Control"))
  expect_true(all(c("K", "cc") %in% colnames(fast_plot$data)))

  peptide_object_with_na <- peptide_object
  peptide_object_with_na@peptide_matrix[1, 1] <- NA_real_
  mean_imputed_plot <- suppressWarnings(ruvCancorFast(
    peptide_object_with_na,
    ctrl = controls,
    num_components_to_impute = 2,
    ruv_grouping_variable = "group",
    simple_imputation_method = "mean"
  ))
  median_imputed_plot <- suppressWarnings(ruvCancorFast(
    peptide_object_with_na,
    ctrl = controls,
    num_components_to_impute = 2,
    ruv_grouping_variable = "group",
    simple_imputation_method = "median"
  ))
  min_imputed_plot <- suppressWarnings(ruvCancorFast(
    peptide_object_with_na,
    ctrl = controls,
    num_components_to_impute = 2,
    ruv_grouping_variable = "group",
    simple_imputation_method = "min"
  ))

  expect_s3_class(mean_imputed_plot, "ggplot")
  expect_s3_class(median_imputed_plot, "ggplot")
  expect_s3_class(min_imputed_plot, "ggplot")

  localGlobalBinding("dpc", function(x) {
    list(centered = x)
  })
  localGlobalBinding("dpcImpute", function(x, dpc) {
    list(E = x)
  })

  cancor_plot <- suppressWarnings(ruvCancor(
    peptide_object,
    ctrl = controls,
    num_components_to_impute = 2,
    ruv_grouping_variable = "group"
  ))

  expect_s3_class(cancor_plot, "ggplot")
  expect_setequal(cancor_plot$data$featureset, c("All", "Control"))
  expect_true(all(c("K", "cc") %in% colnames(cancor_plot$data)))

  localGlobalBinding(
    "RUVIII_C_Varying",
    function(k, Y, M, toCorrect, potentialControls) {
      ruv_call <<- list(k = k, Y = Y, M = M, toCorrect = toCorrect, potentialControls = potentialControls)
      Y + 1
    }
  )

  corrected <- ruvIII_C_Varying(
    peptide_object,
    ruv_grouping_variable = "replicates",
    ruv_number_k = 1,
    ctrl = controls
  )

  expect_s4_class(corrected, "PeptideQuantitativeData")
  expect_identical(ruv_call$k, 1)
  expect_setequal(ruv_call$potentialControls, names(controls))
  expect_equal(corrected@peptide_matrix, t(ruv_call$Y + 1))
})

test_that("PeptideQuantitativeData canonical-correlation methods reject invalid inputs", {
  peptide_object <- newPeptideNormObject(peptide_count = 6, sample_ids = c("S1", "S2"))
  controls <- setNames(rep(TRUE, 6), rownames(peptide_object@peptide_matrix))

  expect_error(
    ruvCancorFast(
      peptide_object,
      ctrl = controls,
      num_components_to_impute = 2,
      ruv_grouping_variable = "missing_group"
    ),
    "not a column in the design matrix",
    fixed = TRUE
  )

  expect_error(
    ruvCancor(
      peptide_object,
      ctrl = controls,
      num_components_to_impute = 0,
      ruv_grouping_variable = "group"
    ),
    "value is invalid",
    fixed = TRUE
  )

  expect_error(
    ruvCancor(
      peptide_object,
      ctrl = controls,
      num_components_to_impute = 2,
      ruv_grouping_variable = "missing_group"
    ),
    "not a column in the design matrix",
    fixed = TRUE
  )

  expect_error(
    ruvCancorFast(
      peptide_object,
      ctrl = controls[1:4],
      num_components_to_impute = 2,
      ruv_grouping_variable = "group"
    ),
    "less than 5",
    fixed = TRUE
  )

  expect_error(
    ruvCancor(
      peptide_object,
      ctrl = controls[1:4],
      num_components_to_impute = 2,
      ruv_grouping_variable = "group"
    ),
    "less than 5",
    fixed = TRUE
  )
})

test_that("peptide percentage-scoring helpers preserve current score contracts", {
  cancor_plot <- newCancorPlot(delta = 0.2)

  expect_identical(.peptide_calculateAdaptiveMaxK(10), 2L)
  expect_identical(.peptide_calculateAdaptiveMaxK(20), 3L)
  expect_identical(.peptide_calculateAdaptiveMaxK(50), 4L)
  expect_identical(.peptide_calculateAdaptiveMaxK(100), 5L)
  expect_equal(.peptide_calculateSeparationScore(cancor_plot, "max_difference"), 0.6)
  expect_equal(.peptide_calculateSeparationScore(cancor_plot, "mean_difference"), 0.4)
  expect_equal(.peptide_calculateSeparationScore(cancor_plot, "weighted_difference"), (0.2 + 0.4 * 2 + 0.6 * 3) / 6)
  expect_equal(.peptide_calculateSeparationScore(cancor_plot, "auc"), 0.8)
  expect_true(is.na(.peptide_calculateSeparationScore(list(data = data.frame()), "max_difference")))
  expect_true(is.na(.peptide_calculateCompositeScore(NA_real_, 1, 0.5, 3)))
  expect_identical(.peptide_calculateCompositeScore(-1, 1, 0.5, 3), 0)
  expect_equal(.peptide_calculateCompositeScore(10, 1, 0.5, 3), 10)
  expect_equal(.peptide_calculateCompositeScore(10, 3, 0.5, 3), 5)
  expect_lt(.peptide_calculateCompositeScore(10, 5, 0.5, 3), 5)
})
