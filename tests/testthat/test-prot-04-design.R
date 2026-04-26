# testthat for Proteomics S4 Objects
# Phase 4 of Proteomics GUI Test Strategy

repoRoot <- pkgload::pkg_path()

skipIfMissingProtDesignExtractedFiles <- function() {
  required_paths <- c(
    "R/mod_prot_design_builder_action_helpers.R",
    "R/mod_prot_design_builder_display_helpers.R",
    "R/mod_prot_design_builder_server_helpers.R",
    "R/mod_prot_design_builder_state_helpers.R",
    "R/mod_prot_design_import_helpers.R"
  )
  missing <- required_paths[!file.exists(file.path(repoRoot, required_paths))]
  if (length(missing) > 0) {
    testthat::skip(
      sprintf(
        "Target-only prot design extracted file(s) not present: %s",
        paste(basename(missing), collapse = ", ")
      )
    )
  }
}

skipIfMissingProtDesignBindings <- function(...) {
  missing <- setdiff(c(...), ls(envir = asNamespace("MultiScholaR")))
  if (length(missing) > 0) {
    testthat::skip(
      sprintf(
        "Target-only extracted helper(s) not present: %s",
        paste(missing, collapse = ", ")
      )
    )
  }
}

skipIfMissingProtDesignExtractedFiles()

skipIfMissingProtDesignBindings(
  "resolveProtDesignImportArtifacts",
  "runProtDesignImportConfirmationFlow",
  "runProtDesignBuilderServerEntryShell"
)

test_that("ProteinQuantitativeData constructor works", {
  # Check for captured checkpoint
  cp_file <- test_path("..", "testdata", "sepsis", "proteomics", "cp04_design_matrix.rds")
  
  if (file.exists(cp_file)) {
    first_line <- tryCatch(readLines(cp_file, n = 1, warn = FALSE), error = function(e) "")
    if (length(first_line) > 0 && identical(first_line[[1]], "version https://git-lfs.github.com/spec/v1")) {
      skip("Snapshot cp04 is a Git LFS pointer and the binary artifact is not present")
    }

    pqd <- readRDS(cp_file)
    # If the CP04 was from a DIA run, it might be PeptideQuantitativeData
    expect_true(inherits(pqd, "ProteinQuantitativeData") || inherits(pqd, "PeptideQuantitativeData"))
    expect_true(nrow(pqd@design_matrix) > 0)
  } else {
    # Fallback to mock
    pqd <- ProteinQuantitativeData(
      protein_quant_table = data.frame(
        Protein.Ids = c("P1", "P2"),
        S1 = c(10, 11),
        S2 = c(20, 21),
        stringsAsFactors = FALSE
      ),
      design_matrix = data.frame(
        Run = c("S1", "S2"),
        group = c("G1", "G2"),
        stringsAsFactors = FALSE
      ),
      sample_id = "Run",
      group_id = "group",
      protein_id_column = "Protein.Ids"
    )
    
    expect_s4_class(pqd, "ProteinQuantitativeData")
    expect_equal(pqd@sample_id, "Run")
    expect_equal(nrow(pqd@protein_quant_table), 2)
  }
})

test_that("PeptideQuantitativeData constructor works", {
  # Mock peptide data
  pept_data <- data.frame(
    Protein.Ids = c("P1", "P1"),
    Peptide.Sequence = c("PEPT1", "PEPT2"),
    Run = c("S1", "S1"),
    Precursor.Quantity = c(100, 200),
    Q.Value = c(0.001, 0.001),
    stringsAsFactors = FALSE
  )
  
  pqd <- PeptideQuantitativeData(
    peptide_data = pept_data,
    design_matrix = data.frame(
      Run = "S1",
      group = "G1",
      stringsAsFactors = FALSE
    ),
    sample_id = "Run",
    group_id = "group",
    protein_id_column = "Protein.Ids",
    peptide_sequence_column = "Peptide.Sequence",
    q_value_column = "Q.Value",
    raw_quantity_column = "Precursor.Quantity"
  )
  
  expect_s4_class(pqd, "PeptideQuantitativeData")
  expect_equal(pqd@peptide_sequence_column, "Peptide.Sequence")
})

test_that("calcPeptideMatrix builds the normalized peptide matrix", {
  pept_data <- data.frame(
    Protein.Ids = c("P1", "P1", "P2", "P2"),
    Stripped.Sequence = c("PEP1", "PEP1", "PEP2", "PEP2"),
    Run = c("S1", "S2", "S1", "S2"),
    Precursor.Quantity = c(100, 110, 200, 210),
    Precursor.Normalised = c(10, 11, 20, 21),
    Q.Value = c(0.001, 0.001, 0.001, 0.001),
    stringsAsFactors = FALSE
  )

  pqd <- new(
    "PeptideQuantitativeData",
    peptide_data = pept_data,
    design_matrix = data.frame(
      Run = c("S1", "S2"),
      group = c("G1", "G2"),
      stringsAsFactors = FALSE
    ),
    sample_id = "Run",
    group_id = "group",
    protein_id_column = "Protein.Ids",
    peptide_sequence_column = "Stripped.Sequence",
    q_value_column = "Q.Value",
    raw_quantity_column = "Precursor.Quantity",
    norm_quantity_column = "Precursor.Normalised"
  )

  pqd <- calcPeptideMatrix(pqd)

  expect_equal(dim(pqd@peptide_matrix), c(2, 2))
  expect_equal(rownames(pqd@peptide_matrix), c("P1%PEP1", "P2%PEP2"))
  expect_equal(colnames(pqd@peptide_matrix), c("S1", "S2"))
  expect_equal(unname(pqd@peptide_matrix["P1%PEP1", "S2"]), 11)
})

test_that("peptide normalization helpers preserve values and log2-transform the matrix", {
  pept_data <- data.frame(
    Protein.Ids = c("P1", "P1", "P2", "P2"),
    Stripped.Sequence = c("PEP1", "PEP1", "PEP2", "PEP2"),
    Run = c("S1", "S2", "S1", "S2"),
    Precursor.Quantity = c(100, 110, 200, 210),
    Precursor.Normalised = c(10, 11, 20, 21),
    Q.Value = c(0.001, 0.001, 0.001, 0.001),
    stringsAsFactors = FALSE
  )

  pqd <- new(
    "PeptideQuantitativeData",
    peptide_data = pept_data,
    design_matrix = data.frame(
      Run = c("S1", "S2"),
      group = c("G1", "G2"),
      stringsAsFactors = FALSE
    ),
    sample_id = "Run",
    group_id = "group",
    protein_id_column = "Protein.Ids",
    peptide_sequence_column = "Stripped.Sequence",
    q_value_column = "Q.Value",
    raw_quantity_column = "Precursor.Quantity",
    norm_quantity_column = "Precursor.Normalised"
  )

  pqd <- calcPeptideMatrix(pqd)
  no_op_norm <- normaliseBetweenSamples(pqd, normalisation_method = "none")

  expect_equal(no_op_norm@peptide_matrix, pqd@peptide_matrix)
  expect_equal(no_op_norm@peptide_data$Precursor.Normalised, pqd@peptide_data$Precursor.Normalised)

  logged <- log2TransformPeptideMatrix(pqd)

  expect_true(logged@is_logged_data)
  expect_equal(unname(logged@peptide_matrix["P1%PEP1", "S1"]), log2(10))
  expect_equal(
    logged@peptide_data$Precursor.Normalised[
      logged@peptide_data$Run == "S2" & logged@peptide_data$Stripped.Sequence == "PEP2"
    ],
    log2(21)
  )
})

test_that("plotRle keeps peptide sample ordering and grouping labels", {
  pept_data <- data.frame(
    Protein.Ids = c("P1", "P1", "P2", "P2"),
    Stripped.Sequence = c("PEP1", "PEP1", "PEP2", "PEP2"),
    Run = c("S1", "S2", "S1", "S2"),
    Precursor.Quantity = c(100, 110, 200, 210),
    Precursor.Normalised = c(10, 11, 20, 21),
    Q.Value = c(0.001, 0.001, 0.001, 0.001),
    stringsAsFactors = FALSE
  )

  pqd <- new(
    "PeptideQuantitativeData",
    peptide_data = pept_data,
    design_matrix = data.frame(
      Run = c("S1", "S2"),
      group = c("G1", "G2"),
      stringsAsFactors = FALSE
    ),
    sample_id = "Run",
    group_id = "group",
    protein_id_column = "Protein.Ids",
    peptide_sequence_column = "Stripped.Sequence",
    q_value_column = "Q.Value",
    raw_quantity_column = "Precursor.Quantity",
    norm_quantity_column = "Precursor.Normalised"
  )

  pqd <- calcPeptideMatrix(pqd)
  pqd@peptide_matrix["P1%PEP1", "S1"] <- Inf

  rle_plot <- plotRle(pqd, grouping_variable = "group")

  expect_s3_class(rle_plot, "ggplot")
  expect_identical(as.character(rle_plot$data$rle.x.factor), c("S1", "S2"))
  expect_identical(rle_plot$data$rowinfo, c("G1", "G2"))
})

test_that("plotPca keeps peptide design annotations and coerces shape labels", {
  pept_data <- data.frame(
    Protein.Ids = rep(c("P1", "P2", "P3", "P4"), each = 3),
    Stripped.Sequence = rep(c("PEP1", "PEP2", "PEP3", "PEP4"), each = 3),
    Run = rep(c("S1", "S2", "S3"), times = 4),
    Precursor.Quantity = c(
      100, 200, 300,
      200, 400, 600,
      300, 600, 900,
      400, 800, 1200
    ),
    Precursor.Normalised = c(
      10, 20, 30,
      20, 40, 60,
      30, 60, 90,
      40, 80, 120
    ),
    Q.Value = rep(0.001, 12),
    stringsAsFactors = FALSE
  )

  pqd <- new(
    "PeptideQuantitativeData",
    peptide_data = pept_data,
    design_matrix = data.frame(
      Run = c("S1", "S2", "S3"),
      group = c("G1", "G1", "G2"),
      batch = c("B1", "B2", "B1"),
      label = c("Sample 1", "Sample 2", "Sample 3"),
      stringsAsFactors = FALSE
    ),
    sample_id = "Run",
    group_id = "group",
    protein_id_column = "Protein.Ids",
    peptide_sequence_column = "Stripped.Sequence",
    q_value_column = "Q.Value",
    raw_quantity_column = "Precursor.Quantity",
    norm_quantity_column = "Precursor.Normalised"
  )

  pqd <- calcPeptideMatrix(pqd)

  pca_plot <- plotPca(
    pqd,
    grouping_variable = "group",
    shape_variable = "batch",
    label_column = "label",
    title = "Peptide PCA"
  )

  expect_s3_class(pca_plot, "ggplot")
  expect_identical(as.character(pca_plot$data$Run), c("S1", "S2", "S3"))
  expect_identical(as.character(pca_plot$data$group), c("G1", "G1", "G2"))
  expect_true(is.factor(pca_plot$data$batch))
  expect_identical(as.character(pca_plot$data$label), c("Sample 1", "Sample 2", "Sample 3"))
  expect_identical(pca_plot$labels$title, "Peptide PCA")
})

test_that("plotDensity filters missing peptide intensities and keeps sample colouring data", {
  pept_data <- data.frame(
    Protein.Ids = c("P1", "P1", "P2", "P2"),
    Stripped.Sequence = c("PEP1", "PEP1", "PEP2", "PEP2"),
    Run = c("S1", "S2", "S1", "S2"),
    Precursor.Quantity = c(100, 200, 300, 400),
    Precursor.Normalised = c(10, NA, 30, 40),
    Q.Value = rep(0.001, 4),
    stringsAsFactors = FALSE
  )

  pqd <- new(
    "PeptideQuantitativeData",
    peptide_data = pept_data,
    design_matrix = data.frame(
      Run = c("S1", "S2"),
      group = c("G1", "G2"),
      stringsAsFactors = FALSE
    ),
    sample_id = "Run",
    group_id = "group",
    protein_id_column = "Protein.Ids",
    peptide_sequence_column = "Stripped.Sequence",
    q_value_column = "Q.Value",
    raw_quantity_column = "Precursor.Quantity",
    norm_quantity_column = "Precursor.Normalised"
  )

  pqd <- calcPeptideMatrix(pqd)

  density_plot <- plotDensity(
    pqd,
    grouping_variable = "group",
    title = "Peptide density"
  )

  expect_s3_class(density_plot, "ggplot")
  expect_identical(density_plot$data$Precursor.Normalised, c(10, 30, 40))
  expect_identical(as.character(density_plot$data$Run), c("S1", "S1", "S2"))
  expect_identical(density_plot$labels$title, "Peptide density")
  expect_identical(density_plot$labels$x, "Log2 Intensity")
  expect_identical(density_plot$labels$y, "Density")
})

test_that("peptide correlation helpers keep within-group pairs and drop pool samples by default", {
  pept_data <- data.frame(
    Protein.Ids = c("P1", "P1", "P1", "P1", "P2", "P2", "P2", "P2", "P3", "P3", "P3", "P3"),
    Stripped.Sequence = c("PEP1", "PEP1", "PEP1", "PEP1", "PEP2", "PEP2", "PEP2", "PEP2", "PEP3", "PEP3", "PEP3", "PEP3"),
    Run = c("S1", "S2", "S3", "S4", "S1", "S2", "S3", "S4", "S1", "S2", "S3", "S4"),
    Precursor.Quantity = c(10, 11, 5, 6, 20, 19, 6, 7, 30, 31, 7, 8),
    Precursor.Normalised = c(10, 11, 5, 6, 20, 19, 6, 7, 30, 31, 7, 8),
    Q.Value = rep(0.001, 12),
    stringsAsFactors = FALSE
  )

  pqd <- new(
    "PeptideQuantitativeData",
    peptide_data = pept_data,
    design_matrix = data.frame(
      Run = c("S1", "S2", "S3", "S4"),
      group = c("bio", "bio", "pool", "pool"),
      replicates = c("rep1", "rep1", "pool", "pool"),
      stringsAsFactors = FALSE
    ),
    sample_id = "Run",
    group_id = "group",
    technical_replicate_id = "replicates",
    protein_id_column = "Protein.Ids",
    peptide_sequence_column = "Stripped.Sequence",
    q_value_column = "Q.Value",
    raw_quantity_column = "Precursor.Quantity",
    norm_quantity_column = "Precursor.Normalised"
  )

  pqd <- calcPeptideMatrix(pqd)

  within_group_pairs <- pearsonCorForSamplePairs(pqd, exclude_pool_samples = FALSE)
  filtered_pairs <- pearsonCorForSamplePairs(pqd)
  pearson_plot <- plotPearson(pqd)

  expect_equal(nrow(within_group_pairs), 2)
  expect_setequal(within_group_pairs$Run.x, c("S1", "S3"))
  expect_setequal(within_group_pairs$Run.y, c("S2", "S4"))
  expect_equal(nrow(filtered_pairs), 1)
  expect_identical(filtered_pairs$Run.x, "S1")
  expect_identical(filtered_pairs$Run.y, "S2")
  expect_identical(filtered_pairs$replicates, "rep1")
  expect_true(all(filtered_pairs$pearson_correlation <= 1))
  expect_true(all(filtered_pairs$pearson_correlation >= -1))
  expect_s3_class(pearson_plot, "ggplot")
  expect_equal(pearson_plot$data, filtered_pairs)
  expect_identical(pearson_plot$labels$x, "Pearson Correlation")
  expect_identical(pearson_plot$labels$y, "Counts")
})

test_that("peptide accession cleanup resolves the best accession and aggregates duplicate peptide rows", {
  pept_data <- data.frame(
    Protein.Ids = c("P1;P2", "P2"),
    Stripped.Sequence = c("PEP1", "PEP1"),
    Run = c("S1", "S1"),
    Precursor.Quantity = c(10, 20),
    Precursor.Normalised = c(100, 200),
    Q.Value = c(0.001, 0.001),
    stringsAsFactors = FALSE
  )

  pqd <- new(
    "PeptideQuantitativeData",
    peptide_data = pept_data,
    design_matrix = data.frame(
      Run = "S1",
      group = "G1",
      stringsAsFactors = FALSE
    ),
    sample_id = "Run",
    group_id = "group",
    protein_id_column = "Protein.Ids",
    peptide_sequence_column = "Stripped.Sequence",
    q_value_column = "Q.Value",
    raw_quantity_column = "Precursor.Quantity",
    norm_quantity_column = "Precursor.Normalised"
  )

  accession_tbl <- data.frame(
    uniprot_acc = c("P1", "P2"),
    gene_name = c("GENE1", "GENE1"),
    cleaned_acc = c("P1", "P2"),
    protein_evidence = c(2, 1),
    status = c("reviewed", "reviewed"),
    is_isoform = c("Canonical", "Canonical"),
    isoform_num = c(0, 0),
    seq_length = c(100, 200),
    annotation_score = c(1, 5),
    stringsAsFactors = FALSE
  )

  cleaned <- chooseBestProteinAccession(
    pqd,
    seqinr_obj = accession_tbl,
    aggregation_method = "mean"
  )

  expect_equal(nrow(cleaned@peptide_data), 1)
  expect_identical(cleaned@peptide_data$Protein.Ids, "P2")
  expect_equal(cleaned@peptide_data$Precursor.Quantity, 15)
  expect_equal(cleaned@peptide_data$Precursor.Normalised, 150)
  expect_equal(dim(cleaned@peptide_matrix), c(1, 1))
  expect_identical(rownames(cleaned@peptide_matrix), "P2%PEP1")
  expect_equal(unname(cleaned@peptide_matrix[1, 1]), 150)
})

test_that("peptide negative-control selection forwards the peptide matrix and records resolved arguments", {
  pept_data <- data.frame(
    Protein.Ids = c("P1", "P1", "P2", "P2"),
    Stripped.Sequence = c("PEP1", "PEP1", "PEP2", "PEP2"),
    Run = c("S1", "S2", "S1", "S2"),
    Precursor.Quantity = c(100, 110, 200, 220),
    Precursor.Normalised = c(10, 11, 20, 22),
    Q.Value = rep(0.001, 4),
    stringsAsFactors = FALSE
  )

  pqd <- new(
    "PeptideQuantitativeData",
    peptide_data = pept_data,
    design_matrix = data.frame(
      Run = c("S1", "S2"),
      group = c("G1", "G2"),
      batch = c("B1", "B2"),
      stringsAsFactors = FALSE
    ),
    sample_id = "Run",
    group_id = "group",
    protein_id_column = "Protein.Ids",
    peptide_sequence_column = "Stripped.Sequence",
    q_value_column = "Q.Value",
    raw_quantity_column = "Precursor.Quantity",
    norm_quantity_column = "Precursor.Normalised"
  )

  pqd <- calcPeptideMatrix(pqd)
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
      c(P1 = TRUE, P2 = FALSE)
    },
    .env = method_env
  )

  controls <- getNegCtrlProtAnovaPeptides(
    pqd,
    ruv_grouping_variable = "batch",
    percentage_as_neg_ctrl = 50,
    ruv_qval_cutoff = 0.2,
    ruv_fdr_method = "BH"
  )

  expect_identical(controls, c(P1 = TRUE, P2 = FALSE))
  expect_equal(helper_call$data_matrix, pqd@peptide_matrix)
  expect_equal(
    helper_call$design_matrix,
    data.frame(batch = c("B1", "B2"), row.names = c("S1", "S2"), stringsAsFactors = FALSE)
  )
  expect_identical(helper_call$grouping_variable, "batch")
  expect_identical(helper_call$percentage_as_neg_ctrl, 50)
  expect_identical(helper_call$num_neg_ctrl, 1)
  expect_identical(helper_call$ruv_qval_cutoff, 0.2)
  expect_identical(helper_call$ruv_fdr_method, "BH")
})

test_that("peptide negative-control optimization keeps the best percentage and adaptive k penalty wiring", {
  pept_data <- data.frame(
    Protein.Ids = rep(paste0("P", 1:6), each = 2),
    Stripped.Sequence = rep(paste0("PEP", 1:6), each = 2),
    Run = rep(c("S1", "S2"), times = 6),
    Precursor.Quantity = seq(100, 210, by = 10),
    Precursor.Normalised = seq(10, 21, by = 1),
    Q.Value = rep(0.001, 12),
    stringsAsFactors = FALSE
  )

  pqd <- new(
    "PeptideQuantitativeData",
    peptide_data = pept_data,
    design_matrix = data.frame(
      Run = c("S1", "S2"),
      group = c("G1", "G2"),
      batch = c("B1", "B2"),
      stringsAsFactors = FALSE
    ),
    sample_id = "Run",
    group_id = "group",
    protein_id_column = "Protein.Ids",
    peptide_sequence_column = "Stripped.Sequence",
    q_value_column = "Q.Value",
    raw_quantity_column = "Precursor.Quantity",
    norm_quantity_column = "Precursor.Normalised"
  )

  pqd <- calcPeptideMatrix(pqd)
  neg_ctrl_calls <- list()
  cancor_calls <- list()
  composite_calls <- list()
  adaptive_calls <- integer()
  method_env <- environment(methods::selectMethod("findBestNegCtrlPercentagePeptides", "PeptideQuantitativeData")@.Data)

  local_mocked_bindings(
    getNegCtrlProtAnovaPeptides = function(theObject,
                                           ruv_grouping_variable,
                                           percentage_as_neg_ctrl,
                                           ruv_qval_cutoff,
                                           ruv_fdr_method) {
      neg_ctrl_calls[[length(neg_ctrl_calls) + 1]] <<- list(
        grouping_variable = ruv_grouping_variable,
        percentage_as_neg_ctrl = percentage_as_neg_ctrl,
        ruv_qval_cutoff = ruv_qval_cutoff,
        ruv_fdr_method = ruv_fdr_method
      )

      setNames(rep(TRUE, 5), paste0("pct", percentage_as_neg_ctrl, "_", seq_len(5)))
    },
    ruvCancorFast = function(theObject,
                             ctrl,
                             num_components_to_impute,
                             ruv_grouping_variable,
                             simple_imputation_method) {
      current_percentage <- as.numeric(sub("^pct([0-9]+)_.*$", "\\1", names(ctrl)[[1]]))

      cancor_calls[[length(cancor_calls) + 1]] <<- list(
        current_percentage = current_percentage,
        num_components_to_impute = num_components_to_impute,
        grouping_variable = ruv_grouping_variable,
        simple_imputation_method = simple_imputation_method
      )

      structure(list(current_percentage = current_percentage), class = "mock_cancorplot")
    },
    .peptide_calculateSeparationScore = function(cancorplot, metric) {
      cancorplot$current_percentage / 10
    },
    findBestKElbow = function(cancorplot, epsilon = 0.05, min_effect = 0.05) {
      if (identical(cancorplot$current_percentage, 10)) {
        2L
      } else {
        5L
      }
    },
    .peptide_calculateCompositeScore = function(separation_score,
                                                best_k,
                                                k_penalty_weight,
                                                max_acceptable_k) {
      composite_calls[[length(composite_calls) + 1]] <<- list(
        separation_score = separation_score,
        best_k = best_k,
        k_penalty_weight = k_penalty_weight,
        max_acceptable_k = max_acceptable_k
      )

      separation_score - best_k / 10
    },
    .peptide_calculateAdaptiveMaxK = function(sample_size) {
      adaptive_calls <<- c(adaptive_calls, sample_size)
      2L
    },
    .env = method_env
  )

  results <- findBestNegCtrlPercentagePeptides(
    pqd,
    percentage_range = c(10, 20),
    num_components_to_impute = 3,
    ruv_grouping_variable = "batch",
    ruv_qval_cutoff = 0.2,
    ruv_fdr_method = "BH",
    separation_metric = "weighted_difference",
    adaptive_k_penalty = TRUE,
    verbose = FALSE
  )

  expect_identical(adaptive_calls, 2L)
  expect_length(neg_ctrl_calls, 2)
  expect_identical(vapply(neg_ctrl_calls, `[[`, character(1), "grouping_variable"), c("batch", "batch"))
  expect_identical(vapply(neg_ctrl_calls, `[[`, numeric(1), "percentage_as_neg_ctrl"), c(10, 20))
  expect_identical(vapply(neg_ctrl_calls, `[[`, numeric(1), "ruv_qval_cutoff"), c(0.2, 0.2))
  expect_identical(vapply(neg_ctrl_calls, `[[`, character(1), "ruv_fdr_method"), c("BH", "BH"))
  expect_length(cancor_calls, 2)
  expect_identical(vapply(cancor_calls, `[[`, numeric(1), "current_percentage"), c(10, 20))
  expect_identical(vapply(cancor_calls, `[[`, numeric(1), "num_components_to_impute"), c(3, 3))
  expect_identical(vapply(cancor_calls, `[[`, character(1), "grouping_variable"), c("batch", "batch"))
  expect_identical(vapply(cancor_calls, `[[`, character(1), "simple_imputation_method"), c("mean", "mean"))
  expect_length(composite_calls, 2)
  expect_identical(vapply(composite_calls, `[[`, numeric(1), "max_acceptable_k"), c(2, 2))
  expect_identical(results$best_percentage, 20)
  expect_identical(results$best_k, 5)
  expect_equal(results$sample_size, 2)
  expect_identical(results$separation_metric_used, "weighted_difference")
  expect_identical(results$adaptive_k_penalty_used, TRUE)
  expect_equal(results$optimization_results$percentage, c(10, 20))
  expect_equal(results$optimization_results$composite_score, c(0.8, 1.5))
  expect_true(all(startsWith(names(results$best_control_genes_index), "pct20_")))
})

test_that("ProteinQuantitativeData validation catches mismatches", {
  pqd <- ProteinQuantitativeData(
    protein_quant_table = data.frame(
      Protein.Ids = c("P1", "P2"),
      S1 = c(10, 11),
      S2 = c(20, 21),
      stringsAsFactors = FALSE
    ),
    design_matrix = data.frame(
      Run = factor(c("S1", "S2")),
      group = c("G1", "G2"),
      stringsAsFactors = TRUE
    ),
    sample_id = "Run",
    group_id = "group",
    protein_id_column = "Protein.Ids"
  )

  expect_true(is.character(pqd@design_matrix$Run))
  expect_identical(pqd@design_matrix$Run, c("S1", "S2"))

  expect_error(
    ProteinQuantitativeData(
      protein_quant_table = data.frame(
        Accession = c("P1", "P2"),
        S1 = c(10, 11),
        S2 = c(20, 21),
        stringsAsFactors = FALSE
      ),
      design_matrix = data.frame(
        Run = c("S1", "S2"),
        group = c("G1", "G2"),
        stringsAsFactors = FALSE
      ),
      sample_id = "Run",
      group_id = "group",
      protein_id_column = "Protein.Ids"
    ),
    "Protein ID column must be in the protein data table",
    fixed = TRUE
  )

  expect_error(
    ProteinQuantitativeData(
      protein_quant_table = data.frame(
        Protein.Ids = c("P1", "P2"),
        S1 = c(10, 11),
        S2 = c(20, 21),
        stringsAsFactors = FALSE
      ),
      design_matrix = data.frame(
        Run = c("S1", "S3"),
        group = c("G1", "G2"),
        stringsAsFactors = FALSE
      ),
      sample_id = "Run",
      group_id = "group",
      protein_id_column = "Protein.Ids"
    ),
    "Samples in protein data and design matrix must be the same",
    fixed = TRUE
  )
})

test_that("resolveProtDesignImportArtifacts prefers a selected FASTA path", {
  import_dir <- tempfile("prot-design-import-")
  dir.create(import_dir)
  writeLines("Run\tgroup", file.path(import_dir, "design_matrix.tab"))
  writeLines("Protein.Ids\tRun\tAbundance", file.path(import_dir, "data_cln.tab"))
  writeLines(">auto\nSEQUENCE", file.path(import_dir, "auto.fasta"))

  selected_fasta <- tempfile("prot-design-selected-", fileext = ".fasta")
  writeLines(">selected\nSEQUENCE", selected_fasta)

  artifacts <- resolveProtDesignImportArtifacts(
    importPath = import_dir,
    selectedFastaPath = selected_fasta
  )

  expect_true(artifacts$ok)
  expect_equal(artifacts$fastaPath, selected_fasta)
  expect_equal(artifacts$designFile, file.path(import_dir, "design_matrix.tab"))
  expect_equal(artifacts$dataClnFile, file.path(import_dir, "data_cln.tab"))
  expect_equal(artifacts$contrastFile, file.path(import_dir, "contrast_strings.tab"))
})

test_that("resolveProtDesignImportArtifacts reports missing required files", {
  import_dir <- tempfile("prot-design-import-missing-")
  dir.create(import_dir)

  artifacts <- resolveProtDesignImportArtifacts(importPath = import_dir)

  expect_false(artifacts$ok)
  expect_match(artifacts$errorMessage, "design_matrix.tab")
})

test_that("runProtDesignImportConfirmationFlow orchestrates imported sidecars and checkpoint handoff", {
  call_log <- character()
  workflow_data <- new.env(parent = emptyenv())
  experiment_paths <- list(
    source_dir = tempfile("prot-design-source-"),
    results_dir = tempfile("prot-design-results-")
  )
  import_artifacts <- list(
    importPath = tempfile("prot-design-import-"),
    fastaPath = tempfile("prot-design-fasta-", fileext = ".fasta")
  )
  imported_artifacts <- list(
    importedDesign = data.frame(Run = "S1", group = "G1", stringsAsFactors = FALSE),
    importedDataCln = data.frame(Protein.Ids = "P1", Run = "S1", Abundance = 1, stringsAsFactors = FALSE),
    importedContrasts = data.frame(contrasts = "groupA-groupB", stringsAsFactors = FALSE)
  )
  qc_trigger_calls <- logical()
  qc_trigger <- function(value) {
    qc_trigger_calls <<- c(qc_trigger_calls, value)
  }
  session_obj <- structure(list(user = "test-session"), class = "test_session")

  result <- runProtDesignImportConfirmationFlow(
    workflowData = workflow_data,
    experimentPaths = experiment_paths,
    importArtifacts = import_artifacts,
    importedArtifacts = imported_artifacts,
    taxonId = 9606,
    organismName = "Homo sapiens",
    session = session_obj,
    qcTrigger = qc_trigger,
    hydrateFastaSidecar = function(workflowData, importPath, sourceDir, resultsDir, fastaPath, organismName) {
      expect_identical(workflowData, workflow_data)
      expect_identical(importPath, import_artifacts$importPath)
      expect_identical(sourceDir, experiment_paths$source_dir)
      expect_identical(resultsDir, experiment_paths$results_dir)
      expect_identical(fastaPath, import_artifacts$fastaPath)
      expect_identical(organismName, "Homo sapiens")
      call_log <<- c(call_log, "fasta")
      invisible(NULL)
    },
    hydrateUniprotSidecar = function(workflowData, importPath, sourceDir) {
      expect_identical(workflowData, workflow_data)
      expect_identical(importPath, import_artifacts$importPath)
      expect_identical(sourceDir, experiment_paths$source_dir)
      call_log <<- c(call_log, "uniprot")
      invisible(NULL)
    },
    initializeWorkflowState = function(workflowData, importedDesign, importedDataCln, importedContrasts, taxonId, organismName) {
      expect_identical(workflowData, workflow_data)
      expect_identical(importedDesign, imported_artifacts$importedDesign)
      expect_identical(importedDataCln, imported_artifacts$importedDataCln)
      expect_identical(importedContrasts, imported_artifacts$importedContrasts)
      expect_identical(taxonId, 9606)
      expect_identical(organismName, "Homo sapiens")
      call_log <<- c(call_log, "initialize")
      "TMT"
    },
    buildStateCheckpointFn = function(workflowData, workflowType, actionLabel, validateColumnMapping) {
      expect_identical(workflowData, workflow_data)
      expect_identical(workflowType, "TMT")
      expect_identical(actionLabel, "Import")
      expect_true(validateColumnMapping)
      call_log <<- c(call_log, "checkpoint")
      "protein_s4_initial"
    },
    completePostCheckpointFn = function(workflowData, experimentPaths, session, qcTrigger, successMessage, successNotificationId, debugQcTrigger = FALSE) {
      expect_identical(workflowData, workflow_data)
      expect_identical(experimentPaths, experiment_paths)
      expect_identical(session, session_obj)
      expect_identical(successMessage, "Design imported successfully!")
      expect_identical(successNotificationId, "importing_design")
      expect_false(debugQcTrigger)
      qcTrigger(TRUE)
      call_log <<- c(call_log, "complete")
      invisible(NULL)
    }
  )

  expect_identical(call_log, c("fasta", "uniprot", "initialize", "checkpoint", "complete"))
  expect_identical(qc_trigger_calls, TRUE)
  expect_identical(result$workflowType, "TMT")
  expect_identical(result$stateName, "protein_s4_initial")
})

test_that("runProtDesignImportObserverShell loads imported artifacts and routes the notification shell", {
  call_log <- character()
  notifications <- list()
  workflow_data <- new.env(parent = emptyenv())
  experiment_paths <- list(
    source_dir = tempfile("prot-design-source-"),
    results_dir = tempfile("prot-design-results-")
  )
  session_obj <- structure(list(user = "test-session"), class = "test_session")

  result <- runProtDesignImportObserverShell(
    workflowData = workflow_data,
    experimentPaths = experiment_paths,
    importPath = "/tmp/prot-design-import",
    selectedFastaPath = "/tmp/prot-design-selected.fasta",
    taxonId = 9606,
    organismName = "Homo sapiens",
    session = session_obj,
    resolveImportArtifacts = function(importPath, selectedFastaPath = NULL) {
      expect_identical(importPath, "/tmp/prot-design-import")
      expect_identical(selectedFastaPath, "/tmp/prot-design-selected.fasta")
      call_log <<- c(call_log, "resolve")
      list(
        ok = TRUE,
        importPath = importPath,
        fastaPath = selectedFastaPath,
        designFile = "/tmp/prot-design-import/design_matrix.tab",
        dataClnFile = "/tmp/prot-design-import/data_cln.tab",
        contrastFile = "/tmp/prot-design-import/contrast_strings.tab"
      )
    },
    loadImportedArtifacts = function(workflowData, experimentPaths, designFile, dataClnFile, contrastFile) {
      expect_identical(workflowData, workflow_data)
      expect_identical(experimentPaths, experiment_paths)
      expect_identical(designFile, "/tmp/prot-design-import/design_matrix.tab")
      expect_identical(dataClnFile, "/tmp/prot-design-import/data_cln.tab")
      expect_identical(contrastFile, "/tmp/prot-design-import/contrast_strings.tab")
      call_log <<- c(call_log, "load")
      list(
        importedDesign = data.frame(Run = "S1", group = "G1", stringsAsFactors = FALSE),
        importedDataCln = data.frame(Protein.Ids = "P1", Run = "S1", Abundance = 1, stringsAsFactors = FALSE),
        importedContrasts = NULL
      )
    },
    runImportConfirmationFlow = function(workflowData, experimentPaths, importArtifacts, importedArtifacts, taxonId, organismName, session, qcTrigger = NULL) {
      expect_identical(workflowData, workflow_data)
      expect_identical(experimentPaths, experiment_paths)
      expect_identical(importArtifacts$importPath, "/tmp/prot-design-import")
      expect_identical(importArtifacts$fastaPath, "/tmp/prot-design-selected.fasta")
      expect_identical(importedArtifacts$importedDesign$Run, "S1")
      expect_identical(taxonId, 9606)
      expect_identical(organismName, "Homo sapiens")
      expect_identical(session, session_obj)
      expect_null(qcTrigger)
      call_log <<- c(call_log, "flow")
      list(stateName = "protein_s4_initial")
    },
    showNotification = function(message, type = NULL, duration = NULL, id = NULL) {
      notifications[[length(notifications) + 1]] <<- list(
        message = message,
        type = type,
        duration = duration,
        id = id
      )
      invisible(id)
    },
    removeNotification = function(id) {
      notifications[[length(notifications) + 1]] <<- list(remove = id)
      invisible(id)
    }
  )

  expect_true(result$ok)
  expect_identical(call_log, c("resolve", "load", "flow"))
  expect_length(notifications, 1)
  expect_identical(notifications[[1]]$message, "Importing design files...")
  expect_identical(notifications[[1]]$id, "importing_design")
  expect_null(notifications[[1]]$type)
  expect_null(notifications[[1]]$duration)
  expect_identical(result$flowResult$stateName, "protein_s4_initial")
})

test_that("runProtDesignBuilderObserverShell shows the processing modal and delegates save flow", {
  call_log <- character()
  workflow_data <- new.env(parent = emptyenv())
  experiment_paths <- list(source_dir = tempfile("prot-design-source-"))
  session_obj <- structure(list(user = "test-session"), class = "test_session")
  result_payload <- list(
    design_matrix = data.frame(Run = "S1", group = "G1", stringsAsFactors = FALSE),
    contrasts = data.frame(contrast = "groupA-groupB", stringsAsFactors = FALSE)
  )

  result <- runProtDesignBuilderObserverShell(
    results = result_payload,
    workflowData = workflow_data,
    experimentPaths = experiment_paths,
    session = session_obj,
    qcTrigger = TRUE,
    showProcessingModal = function() {
      call_log <<- c(call_log, "modal")
      invisible(NULL)
    },
    hydrateBuilderResults = function(results, workflowData) {
      expect_identical(results, result_payload)
      expect_identical(workflowData, workflow_data)
      call_log <<- c(call_log, "hydrate")
      invisible(NULL)
    },
    runBuilderSaveFlow = function(results, workflowData, experimentPaths, session, qcTrigger = NULL) {
      expect_identical(results, result_payload)
      expect_identical(workflowData, workflow_data)
      expect_identical(experimentPaths, experiment_paths)
      expect_identical(session, session_obj)
      expect_true(qcTrigger)
      call_log <<- c(call_log, "save")
      "saved"
    },
    logInfo = function(message) {
      expect_match(message, "Received results from design matrix builder module")
      call_log <<- c(call_log, "log")
      invisible(NULL)
    }
  )

  expect_identical(call_log, c("modal", "log", "hydrate", "save"))
  expect_identical(result, "saved")
})

test_that("formatProtDesignTechRepSummary reports when no technical replicates are assigned", {
  dm <- data.frame(
    Run = c("S1", "S2"),
    group = c("GroupA", "GroupA"),
    replicates = c(1L, 2L),
    tech_reps = c(NA_integer_, NA_integer_),
    stringsAsFactors = FALSE
  )

  expect_identical(
    formatProtDesignTechRepSummary(dm),
    "No technical replicates assigned yet."
  )
})

test_that("formatProtDesignTechRepSummary groups assigned technical replicates by group and replicate", {
  dm <- data.frame(
    Run = c("S1", "S2", "S3"),
    group = c("GroupA", "GroupA", "GroupB"),
    replicates = c(1L, 1L, 2L),
    tech_reps = c(1L, 2L, 1L),
    stringsAsFactors = FALSE
  )

  expect_identical(
    formatProtDesignTechRepSummary(dm),
    paste(
      "Group: GroupA, Biological Replicate: 1\n  Samples: S1, S2\n  Technical Replicates: 1, 2",
      "Group: GroupB, Biological Replicate: 2\n  Samples: S3\n  Technical Replicates: 1",
      sep = "\n\n"
    )
  )
})

test_that("registerProtDesignTechRepSummaryOutput wires the tech-replicate summary render shell", {
  call_log <- list()
  output <- new.env(parent = emptyenv())
  current_design_matrix <- data.frame(
    Run = c("S1", "S2"),
    group = c("GroupA", "GroupA"),
    replicates = c(1L, 1L),
    tech_reps = c(1L, 2L),
    stringsAsFactors = FALSE
  )

  render_value <- registerProtDesignTechRepSummaryOutput(
    output = output,
    designMatrix = function() {
      call_log[[length(call_log) + 1]] <<- list(kind = "getDesignMatrix")
      current_design_matrix
    },
    renderTextFn = function(expr) {
      expr_sub <- substitute(expr)
      expr_env <- parent.frame()
      call_log[[length(call_log) + 1]] <<- list(kind = "renderText")
      rendered <- eval(expr_sub, expr_env)
      call_log[[length(call_log) + 1]] <<- list(kind = "rendered", value = rendered)
      "rendered-text"
    },
    reqFn = function(value) {
      call_log[[length(call_log) + 1]] <<- list(kind = "req", value = value)
      invisible(NULL)
    },
    formatTechRepSummaryFn = function(designMatrix) {
      call_log[[length(call_log) + 1]] <<- list(kind = "format", value = designMatrix)
      "formatted-summary"
    }
  )

  expect_identical(render_value, "rendered-text")
  expect_identical(output$tech_rep_summary, "rendered-text")
  expect_identical(call_log, list(
    list(kind = "renderText"),
    list(kind = "getDesignMatrix"),
    list(kind = "req", value = current_design_matrix),
    list(kind = "format", value = current_design_matrix),
    list(kind = "rendered", value = "formatted-summary")
  ))
})

test_that("registerProtDesignRemovedSamplesDisplayOutput wires the removed-samples render shell", {
  call_log <- list()
  output <- new.env(parent = emptyenv())

  render_value <- registerProtDesignRemovedSamplesDisplayOutput(
    output = output,
    removedSamples = function() {
      call_log[[length(call_log) + 1]] <<- list(kind = "getRemovedSamples")
      c("Sample10", "Sample2")
    },
    renderTextFn = function(expr) {
      expr_sub <- substitute(expr)
      expr_env <- parent.frame()
      call_log[[length(call_log) + 1]] <<- list(kind = "renderText")
      rendered <- eval(expr_sub, expr_env)
      call_log[[length(call_log) + 1]] <<- list(kind = "rendered", value = rendered)
      "rendered-text"
    },
    formatRemovedSamplesDisplayFn = function(removedSamples) {
      call_log[[length(call_log) + 1]] <<- list(kind = "format", value = removedSamples)
      "formatted-removed-samples"
    }
  )

  expect_identical(render_value, "rendered-text")
  expect_identical(output$removed_samples_display, "rendered-text")
  expect_identical(call_log, list(
    list(kind = "renderText"),
    list(kind = "getRemovedSamples"),
    list(kind = "format", value = c("Sample10", "Sample2")),
    list(kind = "rendered", value = "formatted-removed-samples")
  ))
})

test_that("registerProtDesignRangePreviewOutput wires the range-preview render shell", {
  call_log <- list()
  output <- new.env(parent = emptyenv())
  input <- list(
    samples_to_transform = c("SampleA_T1_Rep1", "SampleB_T2_Rep1"),
    range_start = 2,
    range_end = 2
  )

  render_value <- registerProtDesignRangePreviewOutput(
    output = output,
    input = input,
    renderTextFn = function(expr) {
      expr_sub <- substitute(expr)
      expr_env <- parent.frame()
      call_log[[length(call_log) + 1]] <<- list(kind = "renderText")
      rendered <- eval(expr_sub, expr_env)
      call_log[[length(call_log) + 1]] <<- list(kind = "rendered", value = rendered)
      "rendered-text"
    },
    reqFn = function(...) {
      call_log[[length(call_log) + 1]] <<- list(kind = "req", value = list(...))
      invisible(NULL)
    },
    formatRangePreviewFn = function(selectedSamples, rangeStart, rangeEnd) {
      call_log[[length(call_log) + 1]] <<- list(
        kind = "format",
        selectedSamples = selectedSamples,
        rangeStart = rangeStart,
        rangeEnd = rangeEnd
      )
      "formatted-range-preview"
    }
  )

  expect_identical(render_value, "rendered-text")
  expect_identical(output$range_preview, "rendered-text")
  expect_identical(call_log, list(
    list(kind = "renderText"),
    list(kind = "req", value = list(input$samples_to_transform)),
    list(kind = "req", value = list(input$range_start, input$range_end)),
    list(
      kind = "format",
      selectedSamples = input$samples_to_transform,
      rangeStart = input$range_start,
      rangeEnd = input$range_end
    ),
    list(kind = "rendered", value = "formatted-range-preview")
  ))
})

test_that("registerProtDesignAvailableFactorsDisplayOutput wires the available-factors render shell", {
  call_log <- list()
  output <- new.env(parent = emptyenv())

  render_value <- registerProtDesignAvailableFactorsDisplayOutput(
    output = output,
    factors = function() {
      call_log[[length(call_log) + 1]] <<- list(kind = "getFactors")
      c("Condition", "Batch")
    },
    renderUIFn = function(expr) {
      expr_sub <- substitute(expr)
      expr_env <- parent.frame()
      call_log[[length(call_log) + 1]] <<- list(kind = "renderUI")
      rendered <- eval(expr_sub, expr_env)
      call_log[[length(call_log) + 1]] <<- list(kind = "rendered", value = rendered)
      "rendered-ui"
    },
    buildAvailableFactorsDisplayFn = function(currentFactors) {
      call_log[[length(call_log) + 1]] <<- list(kind = "build", value = currentFactors)
      "available-factors-ui"
    }
  )

  expect_identical(render_value, "rendered-ui")
  expect_identical(output$available_factors_display, "rendered-ui")
  expect_identical(call_log, list(
    list(kind = "renderUI"),
    list(kind = "getFactors"),
    list(kind = "build", value = c("Condition", "Batch")),
    list(kind = "rendered", value = "available-factors-ui")
  ))
})

test_that("registerProtDesignDefinedContrastsDisplayOutput wires the defined-contrasts render shell", {
  call_log <- list()
  output <- new.env(parent = emptyenv())
  contrast_data <- data.frame(
    contrast_name = "GA.Control.vs.GA.Elevated",
    numerator = "GA_Control",
    denominator = "GA_Elevated",
    stringsAsFactors = FALSE
  )

  render_value <- registerProtDesignDefinedContrastsDisplayOutput(
    output = output,
    contrasts = function() {
      call_log[[length(call_log) + 1]] <<- list(kind = "getContrasts")
      contrast_data
    },
    formulaString = function() {
      call_log[[length(call_log) + 1]] <<- list(kind = "getFormulaString")
      "~ 0 + group"
    },
    renderUIFn = function(expr) {
      expr_sub <- substitute(expr)
      expr_env <- parent.frame()
      call_log[[length(call_log) + 1]] <<- list(kind = "renderUI")
      rendered <- eval(expr_sub, expr_env)
      call_log[[length(call_log) + 1]] <<- list(kind = "rendered", value = rendered)
      "rendered-ui"
    },
    buildDefinedContrastsDisplayFn = function(contrastData, formulaString) {
      call_log[[length(call_log) + 1]] <<- list(
        kind = "build",
        contrastData = contrastData,
        formulaString = formulaString
      )
      "defined-contrasts-ui"
    }
  )

  expect_identical(render_value, "rendered-ui")
  expect_identical(output$defined_contrasts_display, "rendered-ui")
  expect_identical(call_log, list(
    list(kind = "renderUI"),
    list(kind = "getContrasts"),
    list(kind = "getFormulaString"),
    list(
      kind = "build",
      contrastData = contrast_data,
      formulaString = "~ 0 + group"
    ),
    list(kind = "rendered", value = "defined-contrasts-ui")
  ))
})

test_that("registerProtDesignContrastFactorsInfoOutput wires the contrast-info render shell", {
  call_log <- list()
  output <- new.env(parent = emptyenv())

  render_value <- registerProtDesignContrastFactorsInfoOutput(
    output = output,
    formulaString = function() {
      call_log[[length(call_log) + 1]] <<- list(kind = "getFormulaString")
      "~ 0 + group"
    },
    renderTextFn = function(expr) {
      expr_sub <- substitute(expr)
      expr_env <- parent.frame()
      call_log[[length(call_log) + 1]] <<- list(kind = "renderText")
      rendered <- eval(expr_sub, expr_env)
      call_log[[length(call_log) + 1]] <<- list(kind = "rendered", value = rendered)
      "rendered-text"
    },
    formatContrastFactorsInfoFn = function(formulaString) {
      call_log[[length(call_log) + 1]] <<- list(kind = "format", value = formulaString)
      "formatted-contrast-info"
    }
  )

  expect_identical(render_value, "rendered-text")
  expect_identical(output$contrast_factors_info, "rendered-text")
  expect_identical(call_log, list(
    list(kind = "renderText"),
    list(kind = "getFormulaString"),
    list(kind = "format", value = "~ 0 + group"),
    list(kind = "rendered", value = "formatted-contrast-info")
  ))
})

test_that("registerProtDesignDataTableOutput wires the data-table render shell", {
  call_log <- list()
  output <- new.env(parent = emptyenv())
  design_matrix <- data.frame(
    Run = c("Run1", "Run2", "Run3"),
    group = c("G1", "G2", "G3"),
    stringsAsFactors = FALSE
  )
  filtered_matrix <- design_matrix[c(1, 3), , drop = FALSE]

  render_value <- registerProtDesignDataTableOutput(
    output = output,
    designMatrix = function() {
      call_log[[length(call_log) + 1]] <<- list(kind = "getDesignMatrix")
      design_matrix
    },
    removedSamples = function() {
      call_log[[length(call_log) + 1]] <<- list(kind = "getRemovedSamples")
      "Run2"
    },
    renderDTFn = function(expr, selection = NULL, options = NULL) {
      expr_sub <- substitute(expr)
      expr_env <- parent.frame()
      call_log[[length(call_log) + 1]] <<- list(
        kind = "renderDT",
        selection = selection,
        options = options
      )
      rendered <- eval(expr_sub, expr_env)
      call_log[[length(call_log) + 1]] <<- list(kind = "rendered", value = rendered)
      "rendered-dt"
    },
    reqFn = function(value) {
      call_log[[length(call_log) + 1]] <<- list(kind = "req", value = value)
      invisible(value)
    },
    filterDesignMatrixFn = function(currentDesignMatrix, currentlyRemoved) {
      call_log[[length(call_log) + 1]] <<- list(
        kind = "filter",
        designMatrix = currentDesignMatrix,
        removedSamples = currentlyRemoved
      )
      filtered_matrix
    }
  )

  expect_identical(render_value, "rendered-dt")
  expect_identical(output$data_table, "rendered-dt")
  expect_identical(call_log, list(
    list(
      kind = "renderDT",
      selection = "none",
      options = list(pageLength = 10, scrollX = TRUE, server = FALSE)
    ),
    list(kind = "getDesignMatrix"),
    list(kind = "req", value = design_matrix),
    list(kind = "getDesignMatrix"),
    list(kind = "getRemovedSamples"),
    list(
      kind = "filter",
      designMatrix = design_matrix,
      removedSamples = "Run2"
    ),
    list(kind = "rendered", value = filtered_matrix)
  ))
})

test_that("registerProtDesignDataTableProxyRefreshObserver wires the proxy-refresh observer shell", {
  call_log <- list()
  design_matrix <- data.frame(
    Run = c("Run1", "Run2", "Run3"),
    group = c("G1", "G2", "G3"),
    stringsAsFactors = FALSE
  )
  filtered_matrix <- design_matrix[c(1, 3), , drop = FALSE]

  observer <- registerProtDesignDataTableProxyRefreshObserver(
    proxyDataTable = "proxy-handle",
    designMatrix = function() {
      call_log[[length(call_log) + 1]] <<- list(kind = "getDesignMatrix")
      design_matrix
    },
    removedSamples = function() {
      call_log[[length(call_log) + 1]] <<- list(kind = "getRemovedSamples")
      "Run2"
    },
    observeFn = function(expr) {
      expr_sub <- substitute(expr)
      expr_env <- parent.frame()
      call_log[[length(call_log) + 1]] <<- list(kind = "observe")
      eval(expr_sub, expr_env)
      "observer-registered"
    },
    reqFn = function(value) {
      call_log[[length(call_log) + 1]] <<- list(kind = "req", value = value)
      invisible(value)
    },
    filterDesignMatrixFn = function(currentDesignMatrix, currentlyRemoved) {
      call_log[[length(call_log) + 1]] <<- list(
        kind = "filter",
        designMatrix = currentDesignMatrix,
        removedSamples = currentlyRemoved
      )
      filtered_matrix
    },
    replaceDataFn = function(proxyDataTable, data, resetPaging = TRUE) {
      call_log[[length(call_log) + 1]] <<- list(
        kind = "replaceData",
        proxyDataTable = proxyDataTable,
        data = data,
        resetPaging = resetPaging
      )
      invisible(NULL)
    }
  )

  expect_identical(observer, "observer-registered")
  expect_identical(call_log, list(
    list(kind = "observe"),
    list(kind = "req", value = "proxy-handle"),
    list(kind = "getDesignMatrix"),
    list(kind = "req", value = design_matrix),
    list(kind = "getDesignMatrix"),
    list(kind = "getRemovedSamples"),
    list(
      kind = "filter",
      designMatrix = design_matrix,
      removedSamples = "Run2"
    ),
    list(
      kind = "replaceData",
      proxyDataTable = "proxy-handle",
      data = filtered_matrix,
      resetPaging = FALSE
    )
  ))
})

test_that("registerProtDesignSampleSelectionSyncObserver wires the sample-selection sync observer shell", {
  call_log <- list()
  input <- list(
    sample_to_rename = "Run2",
    selected_runs = c("Run1", "Run4"),
    samples_to_transform = c("Run2", "Run3"),
    tech_rep_samples = c("Run1", "Run4"),
    samples_to_remove = c("Run3", "Run4")
  )
  design_matrix <- data.frame(
    Run = c("Run3", "Run1", "Run2"),
    stringsAsFactors = FALSE
  )

  observer <- registerProtDesignSampleSelectionSyncObserver(
    input = input,
    session = "builder-session",
    designMatrix = function() {
      call_log[[length(call_log) + 1]] <<- list(kind = "getDesignMatrix")
      design_matrix
    },
    removedSamples = function() {
      call_log[[length(call_log) + 1]] <<- list(kind = "getRemovedSamples")
      "Run2"
    },
    observeFn = function(expr) {
      expr_sub <- substitute(expr)
      expr_env <- parent.frame()
      call_log[[length(call_log) + 1]] <<- list(kind = "observe")
      eval(expr_sub, expr_env)
      "observer-registered"
    },
    reqFn = function(value) {
      call_log[[length(call_log) + 1]] <<- list(kind = "req", value = value)
      invisible(value)
    },
    sortRunsFn = function(runs) {
      call_log[[length(call_log) + 1]] <<- list(kind = "sortRuns", runs = runs)
      c("Run1", "Run2", "Run3")
    },
    isolateFn = function(expr) {
      expr_sub <- substitute(expr)
      expr_env <- parent.frame()
      call_log[[length(call_log) + 1]] <<- list(kind = "isolate")
      eval(expr_sub, expr_env)
    },
    updateSelectizeInputFn = function(session, inputId, choices, selected) {
      call_log[[length(call_log) + 1]] <<- list(
        kind = "updateSelectizeInput",
        session = session,
        inputId = inputId,
        choices = choices,
        selected = selected
      )
      invisible(NULL)
    }
  )

  expect_identical(observer, "observer-registered")
  expect_identical(call_log, list(
    list(kind = "observe"),
    list(kind = "getDesignMatrix"),
    list(kind = "req", value = design_matrix),
    list(kind = "getDesignMatrix"),
    list(kind = "sortRuns", runs = c("Run3", "Run1", "Run2")),
    list(kind = "getRemovedSamples"),
    list(kind = "isolate"),
    list(
      kind = "updateSelectizeInput",
      session = "builder-session",
      inputId = "sample_to_rename",
      choices = c("Run1", "Run3"),
      selected = ""
    ),
    list(
      kind = "updateSelectizeInput",
      session = "builder-session",
      inputId = "selected_runs",
      choices = c("Run1", "Run3"),
      selected = "Run1"
    ),
    list(
      kind = "updateSelectizeInput",
      session = "builder-session",
      inputId = "samples_to_transform",
      choices = c("Run1", "Run3"),
      selected = "Run3"
    ),
    list(
      kind = "updateSelectizeInput",
      session = "builder-session",
      inputId = "tech_rep_samples",
      choices = c("Run1", "Run3"),
      selected = "Run1"
    ),
    list(
      kind = "updateSelectizeInput",
      session = "builder-session",
      inputId = "samples_to_remove",
      choices = c("Run1", "Run3"),
      selected = "Run3"
    )
  ))
})

test_that("registerProtDesignFactorGroupSyncObserver wires the factor/group sync observer shell", {
  call_log <- list()
  input <- list(
    factor1_select = "Condition",
    factor2_select = "Batch",
    factor3_select = "",
    contrast_group1 = "Control",
    contrast_group2 = "Treatment"
  )

  observer <- registerProtDesignFactorGroupSyncObserver(
    input = input,
    session = "builder-session",
    factors = function() {
      call_log[[length(call_log) + 1]] <<- list(kind = "getFactors")
      c("Condition", "Batch", "Timepoint")
    },
    groups = function() {
      call_log[[length(call_log) + 1]] <<- list(kind = "getGroups")
      c("Control", "Treatment", "Reference")
    },
    observeFn = function(expr) {
      expr_sub <- substitute(expr)
      expr_env <- parent.frame()
      call_log[[length(call_log) + 1]] <<- list(kind = "observe")
      eval(expr_sub, expr_env)
      "observer-registered"
    },
    isolateFn = function(expr) {
      expr_sub <- substitute(expr)
      expr_env <- parent.frame()
      call_log[[length(call_log) + 1]] <<- list(kind = "isolate")
      eval(expr_sub, expr_env)
    },
    updateSelectInputFn = function(session, inputId, choices, selected) {
      call_log[[length(call_log) + 1]] <<- list(
        kind = "updateSelectInput",
        session = session,
        inputId = inputId,
        choices = choices,
        selected = selected
      )
      invisible(NULL)
    }
  )

  expect_identical(observer, "observer-registered")
  expect_identical(call_log, list(
    list(kind = "observe"),
    list(kind = "isolate"),
    list(kind = "getFactors"),
    list(kind = "getGroups"),
    list(
      kind = "updateSelectInput",
      session = "builder-session",
      inputId = "factor1_select",
      choices = c("", "Condition", "Batch", "Timepoint"),
      selected = "Condition"
    ),
    list(
      kind = "updateSelectInput",
      session = "builder-session",
      inputId = "factor2_select",
      choices = c("", "Condition", "Batch", "Timepoint"),
      selected = "Batch"
    ),
    list(
      kind = "updateSelectInput",
      session = "builder-session",
      inputId = "factor3_select",
      choices = c("", "Condition", "Batch", "Timepoint"),
      selected = ""
    ),
    list(
      kind = "updateSelectInput",
      session = "builder-session",
      inputId = "contrast_group1",
      choices = c("", "Control", "Treatment", "Reference"),
      selected = "Control"
    ),
    list(
      kind = "updateSelectInput",
      session = "builder-session",
      inputId = "contrast_group2",
      choices = c("", "Control", "Treatment", "Reference"),
      selected = "Treatment"
    )
  ))
})

test_that("registerProtDesignInitialStateSyncObserver wires the initial-state reset observer shell", {
  call_log <- list()
  initial_state <- list(
    design_matrix = data.frame(
      Run = c("Run2", "Run1"),
      batch = NA_character_,
      stringsAsFactors = FALSE
    ),
    data_cln = data.frame(
      Run = c("Run2", "Run1"),
      Batch = c("B2", "B1"),
      stringsAsFactors = FALSE
    ),
    groups = c("Control", "Treatment"),
    factors = c("Condition"),
    formula = "~ 0 + group",
    contrasts = data.frame(
      contrast_name = "Treatment_vs_Control",
      numerator = "Treatment",
      denominator = "Control",
      stringsAsFactors = FALSE
    )
  )
  merged_design_matrix <- data.frame(
    Run = c("Run2", "Run1"),
    batch = c("B2", "B1"),
    stringsAsFactors = FALSE
  )
  batch_assignments <- data.frame(
    Run = c("Run2", "Run1"),
    Batch = c("B2", "B1"),
    stringsAsFactors = FALSE
  )

  observer <- registerProtDesignInitialStateSyncObserver(
    input = list(reset_changes = 2),
    session = "builder-session",
    initialState = function() {
      call_log[[length(call_log) + 1]] <<- list(kind = "getInitialState")
      initial_state
    },
    dataTbl = function() {
      call_log[[length(call_log) + 1]] <<- list(kind = "getDataTbl")
      "data-table-signal"
    },
    designMatrix = function(value) {
      call_log[[length(call_log) + 1]] <<- list(kind = "setDesignMatrix", value = value)
      invisible(NULL)
    },
    dataClnReactive = function(value) {
      call_log[[length(call_log) + 1]] <<- list(kind = "setDataCln", value = value)
      invisible(NULL)
    },
    groups = function(value) {
      call_log[[length(call_log) + 1]] <<- list(kind = "setGroups", value = value)
      invisible(NULL)
    },
    factors = function(value) {
      call_log[[length(call_log) + 1]] <<- list(kind = "setFactors", value = value)
      invisible(NULL)
    },
    contrasts = function(value) {
      call_log[[length(call_log) + 1]] <<- list(kind = "setContrasts", value = value)
      invisible(NULL)
    },
    removedSamples = function(value) {
      call_log[[length(call_log) + 1]] <<- list(kind = "setRemovedSamples", value = value)
      invisible(NULL)
    },
    observeFn = function(expr) {
      expr_sub <- substitute(expr)
      expr_env <- parent.frame()
      call_log[[length(call_log) + 1]] <<- list(kind = "observe")
      eval(expr_sub, expr_env)
      "observer-registered"
    },
    bindEventFn = function(x, ..., ignoreNULL = TRUE, ignoreInit = FALSE, once = FALSE, label = NULL) {
      event_subs <- as.list(substitute(list(...)))[-1]
      event_env <- parent.frame()
      event_values <- lapply(event_subs, eval, envir = event_env)
      call_log[[length(call_log) + 1]] <<- list(
        kind = "bindEvent",
        observer = x,
        eventValues = event_values
      )
      "bound-observer"
    },
    reqFn = function(value) {
      call_log[[length(call_log) + 1]] <<- list(kind = "req", value = value)
      invisible(value)
    },
    hasBatchColumnFn = function(dataCln) {
      call_log[[length(call_log) + 1]] <<- list(kind = "hasBatchColumn", value = dataCln)
      TRUE
    },
    buildBatchAssignmentsFn = function(dataCln) {
      call_log[[length(call_log) + 1]] <<- list(kind = "buildBatchAssignments", value = dataCln)
      batch_assignments
    },
    mergeBatchAssignmentsFn = function(currentDesignMatrix, currentBatchAssignments) {
      call_log[[length(call_log) + 1]] <<- list(
        kind = "mergeBatchAssignments",
        designMatrix = currentDesignMatrix,
        batchAssignments = currentBatchAssignments
      )
      merged_design_matrix
    },
    updateSelectizeInputFn = function(session, inputId, choices, selected) {
      call_log[[length(call_log) + 1]] <<- list(
        kind = "updateSelectizeInput",
        session = session,
        inputId = inputId,
        choices = choices,
        selected = selected
      )
      invisible(NULL)
    },
    updateTextInputFn = function(session, inputId, value) {
      call_log[[length(call_log) + 1]] <<- list(
        kind = "updateTextInput",
        session = session,
        inputId = inputId,
        value = value
      )
      invisible(NULL)
    },
    logInfoFn = function(message) {
      call_log[[length(call_log) + 1]] <<- list(kind = "logInfo", message = message)
      invisible(NULL)
    }
  )

  expect_identical(observer, "bound-observer")
  expect_identical(call_log, list(
    list(kind = "observe"),
    list(kind = "getInitialState"),
    list(kind = "req", value = initial_state),
    list(kind = "hasBatchColumn", value = initial_state$data_cln),
    list(
      kind = "logInfo",
      message = "Design Matrix: 'Batch' column detected. Auto-assigning to batch column."
    ),
    list(kind = "buildBatchAssignments", value = initial_state$data_cln),
    list(
      kind = "mergeBatchAssignments",
      designMatrix = initial_state$design_matrix,
      batchAssignments = batch_assignments
    ),
    list(
      kind = "logInfo",
      message = "Design Matrix: Pre-populated 'batch' column with batch assignments."
    ),
    list(
      kind = "logInfo",
      message = "Note: Batch is kept separate from experimental factors to enable cross-batch comparisons."
    ),
    list(kind = "setDesignMatrix", value = merged_design_matrix),
    list(kind = "setDataCln", value = initial_state$data_cln),
    list(kind = "setGroups", value = initial_state$groups),
    list(kind = "setFactors", value = initial_state$factors),
    list(kind = "setContrasts", value = initial_state$contrasts),
    list(kind = "setRemovedSamples", value = character(0)),
    list(
      kind = "updateSelectizeInput",
      session = "builder-session",
      inputId = "sample_to_rename",
      choices = c("Run2", "Run1"),
      selected = ""
    ),
    list(
      kind = "updateSelectizeInput",
      session = "builder-session",
      inputId = "selected_runs",
      choices = c("Run2", "Run1"),
      selected = ""
    ),
    list(
      kind = "updateSelectizeInput",
      session = "builder-session",
      inputId = "samples_to_transform",
      choices = c("Run2", "Run1"),
      selected = ""
    ),
    list(
      kind = "updateSelectizeInput",
      session = "builder-session",
      inputId = "tech_rep_samples",
      choices = c("Run2", "Run1"),
      selected = ""
    ),
    list(
      kind = "updateSelectizeInput",
      session = "builder-session",
      inputId = "samples_to_remove",
      choices = c("Run2", "Run1"),
      selected = ""
    ),
    list(
      kind = "updateTextInput",
      session = "builder-session",
      inputId = "formula_string",
      value = "~ 0 + group"
    ),
    list(kind = "getDataTbl"),
    list(
      kind = "bindEvent",
      observer = "observer-registered",
      eventValues = list("data-table-signal", 2)
    )
  ))
})

test_that("buildProtDesignInitialState normalizes quantity columns and seeds builder defaults", {
  initial_state <- buildProtDesignInitialState(
    dataTbl = data.frame(
      Run = c("Run10", "Run2", "Run1", "Run2"),
      norm_qty = c("10.5", "11.5", "12.5", "13.5"),
      raw_qty = c("100", "101", "102", "103"),
      quantity = c("1000", "1001", "1002", "1003"),
      stringsAsFactors = FALSE
    ),
    configList = list(
      deAnalysisParameters = list(
        formula_string = "~ 0 + group"
      )
    ),
    columnMapping = list(
      norm_quantity_col = "norm_qty",
      raw_quantity_col = "raw_qty",
      quantity_col = "quantity"
    )
  )

  expect_identical(initial_state$design_matrix$Run, c("Run1", "Run2", "Run10"))
  expect_true(is.double(initial_state$data_cln$norm_qty))
  expect_true(is.double(initial_state$data_cln$raw_qty))
  expect_true(is.double(initial_state$data_cln$quantity))
  expect_equal(initial_state$data_cln$norm_qty, c(10.5, 11.5, 12.5, 13.5))
  expect_equal(initial_state$data_cln$raw_qty, c(100, 101, 102, 103))
  expect_equal(initial_state$data_cln$quantity, c(1000, 1001, 1002, 1003))
  expect_identical(initial_state$factors, character(0))
  expect_identical(initial_state$groups, character(0))
  expect_identical(initial_state$formula, "~ 0 + group")
  expect_named(
    initial_state$contrasts,
    c("contrast_name", "numerator", "denominator")
  )
  expect_equal(nrow(initial_state$contrasts), 0)
})

test_that("createProtDesignBuilderState seeds mutable builder reactive shells", {
  seeds <- list()
  reactive_val_fn <- function(value = NULL) {
    current_value <- value
    seeds[length(seeds) + 1] <<- list(value)

    function(new_value) {
      if (!missing(new_value)) {
        current_value <<- new_value
      }

      current_value
    }
  }

  builder_state <- createProtDesignBuilderState(
    reactiveValFn = reactive_val_fn
  )

  expect_identical(
    names(builder_state),
    c(
      "resultRv",
      "designMatrix",
      "dataClnReactive",
      "groups",
      "factors",
      "contrasts",
      "removedSamples"
    )
  )
  expect_identical(
    seeds,
    list(NULL, NULL, NULL, NULL, NULL, NULL, character(0))
  )
  expect_null(builder_state$resultRv())
  expect_null(builder_state$designMatrix())
  expect_null(builder_state$dataClnReactive())
  expect_null(builder_state$groups())
  expect_null(builder_state$factors())
  expect_null(builder_state$contrasts())
  expect_identical(builder_state$removedSamples(), character(0))

  updated_design_matrix <- data.frame(Run = "Run1", stringsAsFactors = FALSE)
  builder_state$designMatrix(updated_design_matrix)
  builder_state$removedSamples("Run1")

  expect_identical(builder_state$designMatrix(), updated_design_matrix)
  expect_identical(builder_state$removedSamples(), "Run1")
})

test_that("createProtDesignMutableStateShells exposes the mutable builder state aliases", {
  builder_state <- list(
    resultRv = function(value) value,
    designMatrix = function(value) value,
    dataClnReactive = function(value) value,
    groups = function(value) value,
    factors = function(value) value,
    contrasts = function(value) value,
    removedSamples = function(value) value
  )
  create_calls <- 0

  mutable_state <- createProtDesignMutableStateShells(
    createBuilderStateFn = function() {
      create_calls <<- create_calls + 1
      builder_state
    }
  )

  expect_identical(create_calls, 1)
  expect_identical(
    names(mutable_state),
    c(
      "resultRv",
      "designMatrix",
      "dataClnReactive",
      "groups",
      "factors",
      "contrasts",
      "removedSamples"
    )
  )
  expect_identical(mutable_state$resultRv, builder_state$resultRv)
  expect_identical(mutable_state$designMatrix, builder_state$designMatrix)
  expect_identical(mutable_state$dataClnReactive, builder_state$dataClnReactive)
  expect_identical(mutable_state$groups, builder_state$groups)
  expect_identical(mutable_state$factors, builder_state$factors)
  expect_identical(mutable_state$contrasts, builder_state$contrasts)
  expect_identical(mutable_state$removedSamples, builder_state$removedSamples)
})

test_that("initializeProtDesignBuilderServerState wires the builder-state setup seam in order", {
  call_log <- character()
  data_tbl_calls <- 0L
  config_list_calls <- 0L
  column_mapping_calls <- 0L
  build_initial_state_calls <- 0L
  fake_input <- list(reset_builder = 1)
  fake_session <- list(ns = function(id) paste0("builder-shell-", id))
  initial_state_seed <- NULL

  data_tbl <- function(...) {
    data_tbl_calls <<- data_tbl_calls + 1L
    data.frame(Run = "S1", stringsAsFactors = FALSE)
  }
  config_list <- function(...) {
    config_list_calls <<- config_list_calls + 1L
    list(globalParameters = list(workflow = "design"))
  }
  column_mapping <- function(...) {
    column_mapping_calls <<- column_mapping_calls + 1L
    list(sample_id = "Run")
  }

  mutable_state <- list(
    resultRv = function() "builder-results",
    designMatrix = function() "design-matrix",
    dataClnReactive = function() "clean-data",
    groups = function() c("Control", "Treatment"),
    factors = function() c("Condition", "Batch"),
    contrasts = function() "contrasts",
    removedSamples = function() "Run2"
  )

  builder_state <- initializeProtDesignBuilderServerState(
    input = fake_input,
    session = fake_session,
    dataTbl = data_tbl,
    configList = config_list,
    columnMapping = column_mapping,
    buildInitialState = function(...) {
      build_initial_state_calls <<- build_initial_state_calls + 1L
      list(...)
    },
    createInitialStateReactive = function(dataTbl, configList, columnMapping, buildInitialState) {
      expect_identical(dataTbl, data_tbl)
      expect_identical(configList, config_list)
      expect_identical(columnMapping, column_mapping)
      expect_true(is.function(buildInitialState))
      initial_state_seed <<- function() "initial-state"
      call_log <<- c(call_log, "initialStateReactive")
      initial_state_seed
    },
    createMutableStateShells = function() {
      call_log <<- c(call_log, "mutableState")
      mutable_state
    },
    registerInitialStateSyncObserver = function(input, session, initialState, dataTbl, designMatrix, dataClnReactive, groups, factors, contrasts, removedSamples) {
      expect_identical(input, fake_input)
      expect_identical(session, fake_session)
      expect_identical(initialState, initial_state_seed)
      expect_identical(dataTbl, data_tbl)
      expect_identical(designMatrix, mutable_state$designMatrix)
      expect_identical(dataClnReactive, mutable_state$dataClnReactive)
      expect_identical(groups, mutable_state$groups)
      expect_identical(factors, mutable_state$factors)
      expect_identical(contrasts, mutable_state$contrasts)
      expect_identical(removedSamples, mutable_state$removedSamples)
      call_log <<- c(call_log, "initialObserver")
      invisible(NULL)
    }
  )

  expect_identical(
    names(builder_state),
    c(
      "initialState",
      "resultRv",
      "designMatrix",
      "dataClnReactive",
      "groups",
      "factors",
      "contrasts",
      "removedSamples"
    )
  )
  expect_identical(builder_state$initialState, initial_state_seed)
  expect_identical(builder_state$resultRv, mutable_state$resultRv)
  expect_identical(builder_state$designMatrix, mutable_state$designMatrix)
  expect_identical(builder_state$dataClnReactive, mutable_state$dataClnReactive)
  expect_identical(builder_state$groups, mutable_state$groups)
  expect_identical(builder_state$factors, mutable_state$factors)
  expect_identical(builder_state$contrasts, mutable_state$contrasts)
  expect_identical(builder_state$removedSamples, mutable_state$removedSamples)
  expect_identical(
    call_log,
    c("initialStateReactive", "mutableState", "initialObserver")
  )
  expect_identical(data_tbl_calls, 0L)
  expect_identical(config_list_calls, 0L)
  expect_identical(column_mapping_calls, 0L)
  expect_identical(build_initial_state_calls, 0L)
})

test_that("createProtDesignDataTableProxy seeds the default builder proxy id", {
  proxy_calls <- list()

  proxy <- createProtDesignDataTableProxy(
    dataTableProxyFn = function(outputId) {
      proxy_calls[[length(proxy_calls) + 1]] <<- outputId
      list(proxyId = outputId)
    }
  )

  expect_identical(proxy_calls, list("data_table"))
  expect_identical(proxy, list(proxyId = "data_table"))
})

test_that("registerProtDesignRenderOutputShells wires the builder output registration fan-out", {
  call_log <- list()
  input_values <- list(
    formula_string = "~ 0 + group",
    samples_to_transform = c("Run1", "Run2")
  )
  output_env <- new.env(parent = emptyenv())
  proxy_data_table <- list(proxy = "data-table")
  design_matrix <- function() "design-matrix"
  removed_samples <- function() "Run2"
  factors_rv <- function() c("Condition", "Batch")
  contrasts_rv <- function() data.frame(
    contrast_name = "ConditionA.vs.ConditionB",
    numerator = "ConditionA",
    denominator = "ConditionB",
    stringsAsFactors = FALSE
  )
  session_obj <- list(ns = function(id) paste0("builder-", id))

  result <- registerProtDesignRenderOutputShells(
    input = input_values,
    output = output_env,
    session = session_obj,
    proxyDataTable = proxy_data_table,
    designMatrix = design_matrix,
    removedSamples = removed_samples,
    factors = factors_rv,
    contrasts = contrasts_rv,
    registerDataTableOutput = function(output, designMatrix, removedSamples) {
      call_log[[length(call_log) + 1]] <<- list(
        kind = "dataTable",
        output = identical(output, output_env),
        designMatrix = identical(designMatrix, design_matrix),
        removedSamples = identical(removedSamples, removed_samples)
      )
    },
    registerDataTableProxyRefreshObserver = function(proxyDataTable, designMatrix, removedSamples) {
      call_log[[length(call_log) + 1]] <<- list(
        kind = "proxyRefresh",
        proxyDataTable = identical(proxyDataTable, proxy_data_table),
        designMatrix = identical(designMatrix, design_matrix),
        removedSamples = identical(removedSamples, removed_samples)
      )
    },
    registerAvailableFactorsDisplayOutput = function(output, factors) {
      call_log[[length(call_log) + 1]] <<- list(
        kind = "availableFactors",
        output = identical(output, output_env),
        factors = identical(factors, factors_rv)
      )
    },
    registerDefinedContrastsDisplayOutput = function(output, contrasts, formulaString) {
      call_log[[length(call_log) + 1]] <<- list(
        kind = "definedContrasts",
        output = identical(output, output_env),
        contrasts = identical(contrasts, contrasts_rv),
        formulaString = formulaString()
      )
    },
    registerRangePreviewOutput = function(output, input) {
      call_log[[length(call_log) + 1]] <<- list(
        kind = "rangePreview",
        output = identical(output, output_env),
        input = identical(input, input_values)
      )
    },
    registerTechRepSummaryOutput = function(output, designMatrix) {
      call_log[[length(call_log) + 1]] <<- list(
        kind = "techRepSummary",
        output = identical(output, output_env),
        designMatrix = identical(designMatrix, design_matrix)
      )
    },
    registerRemovedSamplesDisplayOutput = function(output, removedSamples) {
      call_log[[length(call_log) + 1]] <<- list(
        kind = "removedSamples",
        output = identical(output, output_env),
        removedSamples = identical(removedSamples, removed_samples)
      )
    },
    registerReplicateInputsOutput = function(output, input, nsFn) {
      call_log[[length(call_log) + 1]] <<- list(
        kind = "replicateInputs",
        output = identical(output, output_env),
        input = identical(input, input_values),
        namespacedInput = nsFn("replicate_start")
      )
    },
    registerContrastFactorsInfoOutput = function(output, formulaString) {
      call_log[[length(call_log) + 1]] <<- list(
        kind = "contrastInfo",
        output = identical(output, output_env),
        formulaString = formulaString()
      )
    }
  )

  expect_null(result)
  expect_identical(call_log, list(
    list(kind = "dataTable", output = TRUE, designMatrix = TRUE, removedSamples = TRUE),
    list(kind = "proxyRefresh", proxyDataTable = TRUE, designMatrix = TRUE, removedSamples = TRUE),
    list(kind = "availableFactors", output = TRUE, factors = TRUE),
    list(kind = "definedContrasts", output = TRUE, contrasts = TRUE, formulaString = "~ 0 + group"),
    list(kind = "rangePreview", output = TRUE, input = TRUE),
    list(kind = "techRepSummary", output = TRUE, designMatrix = TRUE),
    list(kind = "removedSamples", output = TRUE, removedSamples = TRUE),
    list(kind = "replicateInputs", output = TRUE, input = TRUE, namespacedInput = "builder-replicate_start"),
    list(kind = "contrastInfo", output = TRUE, formulaString = "~ 0 + group")
  ))
})

test_that("registerProtDesignInputSyncObserverShells wires the input sync registration fan-out", {
  call_log <- list()
  input_values <- list(
    sample_to_rename = "S1",
    factor1_select = "Condition"
  )
  session_obj <- list(ns = function(id) paste0("builder-", id))
  design_matrix <- function() "design-matrix"
  removed_samples <- function() "Run2"
  factors_rv <- function() c("Condition", "Batch")
  groups_rv <- function() c("Control", "Treatment")

  result <- registerProtDesignInputSyncObserverShells(
    input = input_values,
    session = session_obj,
    designMatrix = design_matrix,
    removedSamples = removed_samples,
    factors = factors_rv,
    groups = groups_rv,
    registerSampleSelectionSyncObserver = function(input, session, designMatrix, removedSamples) {
      call_log[[length(call_log) + 1]] <<- list(
        kind = "sampleSelection",
        input = identical(input, input_values),
        session = identical(session, session_obj),
        designMatrix = identical(designMatrix, design_matrix),
        removedSamples = identical(removedSamples, removed_samples)
      )
    },
    registerFactorGroupSyncObserver = function(input, session, factors, groups) {
      call_log[[length(call_log) + 1]] <<- list(
        kind = "factorGroup",
        input = identical(input, input_values),
        session = identical(session, session_obj),
        factors = identical(factors, factors_rv),
        groups = identical(groups, groups_rv)
      )
    }
  )

  expect_null(result)
  expect_identical(call_log, list(
    list(
      kind = "sampleSelection",
      input = TRUE,
      session = TRUE,
      designMatrix = TRUE,
      removedSamples = TRUE
    ),
    list(
      kind = "factorGroup",
      input = TRUE,
      session = TRUE,
      factors = TRUE,
      groups = TRUE
    )
  ))
})

test_that("registerProtDesignResetAndSaveObserverShells wires the reset/save registration fan-out", {
  call_log <- list()
  input_values <- list(
    reset_changes = 1,
    confirm_reset = 2,
    save_results = 3
  )
  session_obj <- list(ns = function(id) paste0("builder-", id))
  initial_state <- function() "initial-state"
  design_matrix <- function() "design-matrix"
  data_cln_reactive <- function() "clean-data"
  removed_samples <- function() "Run2"
  factors_rv <- function() c("Condition", "Batch")
  groups_rv <- function() c("Control", "Treatment")
  contrasts_rv <- function() data.frame(
    contrast_name = "ConditionA.vs.ConditionB",
    numerator = "ConditionA",
    denominator = "ConditionB",
    stringsAsFactors = FALSE
  )
  config_list <- function() list(formula = "~ 0 + group")
  result_rv <- function(value) value

  result <- registerProtDesignResetAndSaveObserverShells(
    input = input_values,
    session = session_obj,
    initialState = initial_state,
    designMatrix = design_matrix,
    dataClnReactive = data_cln_reactive,
    removedSamples = removed_samples,
    factors = factors_rv,
    groups = groups_rv,
    contrasts = contrasts_rv,
    configList = config_list,
    resultRv = result_rv,
    registerResetRequestObserver = function(input, session) {
      call_log[[length(call_log) + 1]] <<- list(
        kind = "resetRequest",
        input = identical(input, input_values),
        session = identical(session, session_obj)
      )
    },
    registerResetConfirmationObserver = function(
      input,
      initialState,
      designMatrix,
      dataClnReactive,
      removedSamples,
      factors,
      groups,
      contrasts,
      session
    ) {
      call_log[[length(call_log) + 1]] <<- list(
        kind = "resetConfirm",
        input = identical(input, input_values),
        initialState = identical(initialState, initial_state),
        designMatrix = identical(designMatrix, design_matrix),
        dataClnReactive = identical(dataClnReactive, data_cln_reactive),
        removedSamples = identical(removedSamples, removed_samples),
        factors = identical(factors, factors_rv),
        groups = identical(groups, groups_rv),
        contrasts = identical(contrasts, contrasts_rv),
        session = identical(session, session_obj)
      )
    },
    registerSaveResultsObserver = function(
      input,
      designMatrix,
      removedSamples,
      dataCln,
      contrastData,
      configList,
      resultRv
    ) {
      call_log[[length(call_log) + 1]] <<- list(
        kind = "saveResults",
        input = identical(input, input_values),
        designMatrix = identical(designMatrix, design_matrix),
        removedSamples = identical(removedSamples, removed_samples),
        dataCln = identical(dataCln, data_cln_reactive),
        contrastData = identical(contrastData, contrasts_rv),
        configList = identical(configList, config_list),
        resultRv = identical(resultRv, result_rv)
      )
    }
  )

  expect_null(result)
  expect_identical(call_log, list(
    list(kind = "resetRequest", input = TRUE, session = TRUE),
    list(
      kind = "resetConfirm",
      input = TRUE,
      initialState = TRUE,
      designMatrix = TRUE,
      dataClnReactive = TRUE,
      removedSamples = TRUE,
      factors = TRUE,
      groups = TRUE,
      contrasts = TRUE,
      session = TRUE
    ),
    list(
      kind = "saveResults",
      input = TRUE,
      designMatrix = TRUE,
      removedSamples = TRUE,
      dataCln = TRUE,
      contrastData = TRUE,
      configList = TRUE,
      resultRv = TRUE
    )
  ))
})

test_that("registerProtDesignEventObserverShells wires the event observer registration fan-out", {
  call_log <- list()
  input_values <- list(
    rename_sample = 1,
    add_factor = 2,
    add_contrast = 3,
    save_results = 4
  )
  session_obj <- list(ns = function(id) paste0("builder-", id))
  initial_state <- function() "initial-state"
  design_matrix <- function() "design-matrix"
  data_cln_reactive <- function() "clean-data"
  removed_samples <- function() "Run2"
  factors_rv <- function() c("Condition", "Batch")
  groups_rv <- function() c("Control", "Treatment")
  contrasts_rv <- function() "contrasts"
  config_list <- function() list(formula = "~ 0 + group")
  result_rv <- function(value) value

  result <- registerProtDesignEventObserverShells(
    input = input_values,
    session = session_obj,
    initialState = initial_state,
    designMatrix = design_matrix,
    dataClnReactive = data_cln_reactive,
    removedSamples = removed_samples,
    factors = factors_rv,
    groups = groups_rv,
    contrasts = contrasts_rv,
    configList = config_list,
    resultRv = result_rv,
    registerRenameObserverShells = function(input, designMatrix, dataCln, session) {
      call_log[[length(call_log) + 1]] <<- list(
        kind = "rename",
        input = identical(input, input_values),
        designMatrix = identical(designMatrix, design_matrix),
        dataCln = identical(dataCln, data_cln_reactive),
        session = identical(session, session_obj)
      )
    },
    registerFactorMetadataObserverShells = function(input, factors, designMatrix, groups, session) {
      call_log[[length(call_log) + 1]] <<- list(
        kind = "factorMetadata",
        input = identical(input, input_values),
        factors = identical(factors, factors_rv),
        designMatrix = identical(designMatrix, design_matrix),
        groups = identical(groups, groups_rv),
        session = identical(session, session_obj)
      )
    },
    registerActionObserverShells = function(input, designMatrix, contrasts, removedSamples, session) {
      call_log[[length(call_log) + 1]] <<- list(
        kind = "action",
        input = identical(input, input_values),
        designMatrix = identical(designMatrix, design_matrix),
        contrasts = identical(contrasts, contrasts_rv),
        removedSamples = identical(removedSamples, removed_samples),
        session = identical(session, session_obj)
      )
    },
    registerResetAndSaveObserverShells = function(
      input,
      session,
      initialState,
      designMatrix,
      dataClnReactive,
      removedSamples,
      factors,
      groups,
      contrasts,
      configList,
      resultRv
    ) {
      call_log[[length(call_log) + 1]] <<- list(
        kind = "resetAndSave",
        input = identical(input, input_values),
        session = identical(session, session_obj),
        initialState = identical(initialState, initial_state),
        designMatrix = identical(designMatrix, design_matrix),
        dataClnReactive = identical(dataClnReactive, data_cln_reactive),
        removedSamples = identical(removedSamples, removed_samples),
        factors = identical(factors, factors_rv),
        groups = identical(groups, groups_rv),
        contrasts = identical(contrasts, contrasts_rv),
        configList = identical(configList, config_list),
        resultRv = identical(resultRv, result_rv)
      )
    }
  )

  expect_null(result)
  expect_identical(call_log, list(
    list(
      kind = "rename",
      input = TRUE,
      designMatrix = TRUE,
      dataCln = TRUE,
      session = TRUE
    ),
    list(
      kind = "factorMetadata",
      input = TRUE,
      factors = TRUE,
      designMatrix = TRUE,
      groups = TRUE,
      session = TRUE
    ),
    list(
      kind = "action",
      input = TRUE,
      designMatrix = TRUE,
      contrasts = TRUE,
      removedSamples = TRUE,
      session = TRUE
    ),
    list(
      kind = "resetAndSave",
      input = TRUE,
      session = TRUE,
      initialState = TRUE,
      designMatrix = TRUE,
      dataClnReactive = TRUE,
      removedSamples = TRUE,
      factors = TRUE,
      groups = TRUE,
      contrasts = TRUE,
      configList = TRUE,
      resultRv = TRUE
    )
  ))
})

test_that("registerProtDesignEventObserverShells passes callback seeds opaquely through the registration fan-out", {
  call_log <- character()
  initial_state_calls <- 0L
  design_matrix_calls <- 0L
  data_cln_calls <- 0L
  removed_samples_calls <- 0L
  factors_calls <- 0L
  groups_calls <- 0L
  contrasts_calls <- 0L
  config_list_calls <- 0L
  result_rv_calls <- 0L
  input_values <- list(
    rename_sample = 1,
    add_factor = 2,
    add_contrast = 3,
    save_results = 4
  )
  session_obj <- list(ns = function(id) paste0("builder-", id))
  initial_state <- function(...) {
    initial_state_calls <<- initial_state_calls + 1L
    "initial-state"
  }
  design_matrix <- function(...) {
    design_matrix_calls <<- design_matrix_calls + 1L
    "design-matrix"
  }
  data_cln_reactive <- function(...) {
    data_cln_calls <<- data_cln_calls + 1L
    "clean-data"
  }
  removed_samples <- function(...) {
    removed_samples_calls <<- removed_samples_calls + 1L
    "Run2"
  }
  factors_rv <- function(...) {
    factors_calls <<- factors_calls + 1L
    c("Condition", "Batch")
  }
  groups_rv <- function(...) {
    groups_calls <<- groups_calls + 1L
    c("Control", "Treatment")
  }
  contrasts_rv <- function(...) {
    contrasts_calls <<- contrasts_calls + 1L
    "contrasts"
  }
  config_list <- function(...) {
    config_list_calls <<- config_list_calls + 1L
    list(formula = "~ 0 + group")
  }
  result_rv <- function(...) {
    result_rv_calls <<- result_rv_calls + 1L
    "result-rv"
  }

  result <- registerProtDesignEventObserverShells(
    input = input_values,
    session = session_obj,
    initialState = initial_state,
    designMatrix = design_matrix,
    dataClnReactive = data_cln_reactive,
    removedSamples = removed_samples,
    factors = factors_rv,
    groups = groups_rv,
    contrasts = contrasts_rv,
    configList = config_list,
    resultRv = result_rv,
    registerRenameObserverShells = function(input, designMatrix, dataCln, session) {
      expect_identical(input, input_values)
      expect_identical(designMatrix, design_matrix)
      expect_identical(dataCln, data_cln_reactive)
      expect_identical(session, session_obj)
      call_log <<- c(call_log, "rename")
    },
    registerFactorMetadataObserverShells = function(input, factors, designMatrix, groups, session) {
      expect_identical(input, input_values)
      expect_identical(factors, factors_rv)
      expect_identical(designMatrix, design_matrix)
      expect_identical(groups, groups_rv)
      expect_identical(session, session_obj)
      call_log <<- c(call_log, "factorMetadata")
    },
    registerActionObserverShells = function(input, designMatrix, contrasts, removedSamples, session) {
      expect_identical(input, input_values)
      expect_identical(designMatrix, design_matrix)
      expect_identical(contrasts, contrasts_rv)
      expect_identical(removedSamples, removed_samples)
      expect_identical(session, session_obj)
      call_log <<- c(call_log, "action")
    },
    registerResetAndSaveObserverShells = function(
      input,
      session,
      initialState,
      designMatrix,
      dataClnReactive,
      removedSamples,
      factors,
      groups,
      contrasts,
      configList,
      resultRv
    ) {
      expect_identical(input, input_values)
      expect_identical(session, session_obj)
      expect_identical(initialState, initial_state)
      expect_identical(designMatrix, design_matrix)
      expect_identical(dataClnReactive, data_cln_reactive)
      expect_identical(removedSamples, removed_samples)
      expect_identical(factors, factors_rv)
      expect_identical(groups, groups_rv)
      expect_identical(contrasts, contrasts_rv)
      expect_identical(configList, config_list)
      expect_identical(resultRv, result_rv)
      call_log <<- c(call_log, "resetAndSave")
    }
  )

  expect_null(result)
  expect_identical(
    call_log,
    c("rename", "factorMetadata", "action", "resetAndSave")
  )
  expect_identical(initial_state_calls, 0L)
  expect_identical(design_matrix_calls, 0L)
  expect_identical(data_cln_calls, 0L)
  expect_identical(removed_samples_calls, 0L)
  expect_identical(factors_calls, 0L)
  expect_identical(groups_calls, 0L)
  expect_identical(contrasts_calls, 0L)
  expect_identical(config_list_calls, 0L)
  expect_identical(result_rv_calls, 0L)
})

test_that("registerProtDesignBuilderServerShells wires the top-level builder registration fan-out", {
  call_log <- list()
  input_values <- list(reset_builder = 1)
  output_env <- new.env(parent = emptyenv())
  session_obj <- list(ns = function(id) paste0("builder-", id))
  initial_state <- function() "initial-state"
  proxy_data_table <- list(proxy = "data-table")
  design_matrix <- function() "design-matrix"
  data_cln_reactive <- function() "clean-data"
  removed_samples <- function() "Run2"
  factors_rv <- function() c("Condition", "Batch")
  groups_rv <- function() c("Control", "Treatment")
  contrasts_rv <- function() "contrasts"
  config_list <- function() list(formula = "~ 0 + group")
  result_rv <- function(value) value

  result <- registerProtDesignBuilderServerShells(
    input = input_values,
    output = output_env,
    session = session_obj,
    initialState = initial_state,
    proxyDataTable = proxy_data_table,
    designMatrix = design_matrix,
    dataClnReactive = data_cln_reactive,
    removedSamples = removed_samples,
    factors = factors_rv,
    groups = groups_rv,
    contrasts = contrasts_rv,
    configList = config_list,
    resultRv = result_rv,
    registerInputSyncObserverShells = function(input, session, designMatrix, removedSamples, factors, groups) {
      call_log[[length(call_log) + 1]] <<- list(
        kind = "inputSync",
        input = identical(input, input_values),
        session = identical(session, session_obj),
        designMatrix = identical(designMatrix, design_matrix),
        removedSamples = identical(removedSamples, removed_samples),
        factors = identical(factors, factors_rv),
        groups = identical(groups, groups_rv)
      )
    },
    registerRenderOutputShells = function(input, output, session, proxyDataTable, designMatrix, removedSamples, factors, contrasts) {
      call_log[[length(call_log) + 1]] <<- list(
        kind = "renderOutputs",
        input = identical(input, input_values),
        output = identical(output, output_env),
        session = identical(session, session_obj),
        proxyDataTable = identical(proxyDataTable, proxy_data_table),
        designMatrix = identical(designMatrix, design_matrix),
        removedSamples = identical(removedSamples, removed_samples),
        factors = identical(factors, factors_rv),
        contrasts = identical(contrasts, contrasts_rv)
      )
    },
    registerEventObserverShells = function(input, session, initialState, designMatrix, dataClnReactive, removedSamples, factors, groups, contrasts, configList, resultRv) {
      call_log[[length(call_log) + 1]] <<- list(
        kind = "eventShells",
        input = identical(input, input_values),
        session = identical(session, session_obj),
        initialState = identical(initialState, initial_state),
        designMatrix = identical(designMatrix, design_matrix),
        dataClnReactive = identical(dataClnReactive, data_cln_reactive),
        removedSamples = identical(removedSamples, removed_samples),
        factors = identical(factors, factors_rv),
        groups = identical(groups, groups_rv),
        contrasts = identical(contrasts, contrasts_rv),
        configList = identical(configList, config_list),
        resultRv = identical(resultRv, result_rv)
      )
    }
  )

  expect_null(result)
  expect_identical(call_log, list(
    list(
      kind = "inputSync",
      input = TRUE,
      session = TRUE,
      designMatrix = TRUE,
      removedSamples = TRUE,
      factors = TRUE,
      groups = TRUE
    ),
    list(
      kind = "renderOutputs",
      input = TRUE,
      output = TRUE,
      session = TRUE,
      proxyDataTable = TRUE,
      designMatrix = TRUE,
      removedSamples = TRUE,
      factors = TRUE,
      contrasts = TRUE
    ),
    list(
      kind = "eventShells",
      input = TRUE,
      session = TRUE,
      initialState = TRUE,
      designMatrix = TRUE,
      dataClnReactive = TRUE,
      removedSamples = TRUE,
      factors = TRUE,
      groups = TRUE,
      contrasts = TRUE,
      configList = TRUE,
      resultRv = TRUE
    )
  ))
})

test_that("registerProtDesignBuilderServerShells passes callback seeds opaquely through the top-level registration fan-out", {
  call_log <- character()
  initial_state_calls <- 0L
  design_matrix_calls <- 0L
  data_cln_calls <- 0L
  removed_samples_calls <- 0L
  factors_calls <- 0L
  groups_calls <- 0L
  contrasts_calls <- 0L
  config_list_calls <- 0L
  result_rv_calls <- 0L
  input_values <- list(reset_builder = 1)
  output_env <- new.env(parent = emptyenv())
  session_obj <- list(ns = function(id) paste0("builder-", id))
  proxy_data_table <- list(proxy = "data-table")
  initial_state <- function(...) {
    initial_state_calls <<- initial_state_calls + 1L
    "initial-state"
  }
  design_matrix <- function(...) {
    design_matrix_calls <<- design_matrix_calls + 1L
    "design-matrix"
  }
  data_cln_reactive <- function(...) {
    data_cln_calls <<- data_cln_calls + 1L
    "clean-data"
  }
  removed_samples <- function(...) {
    removed_samples_calls <<- removed_samples_calls + 1L
    "Run2"
  }
  factors_rv <- function(...) {
    factors_calls <<- factors_calls + 1L
    c("Condition", "Batch")
  }
  groups_rv <- function(...) {
    groups_calls <<- groups_calls + 1L
    c("Control", "Treatment")
  }
  contrasts_rv <- function(...) {
    contrasts_calls <<- contrasts_calls + 1L
    "contrasts"
  }
  config_list <- function(...) {
    config_list_calls <<- config_list_calls + 1L
    list(formula = "~ 0 + group")
  }
  result_rv <- function(...) {
    result_rv_calls <<- result_rv_calls + 1L
    "result-rv"
  }

  result <- registerProtDesignBuilderServerShells(
    input = input_values,
    output = output_env,
    session = session_obj,
    initialState = initial_state,
    proxyDataTable = proxy_data_table,
    designMatrix = design_matrix,
    dataClnReactive = data_cln_reactive,
    removedSamples = removed_samples,
    factors = factors_rv,
    groups = groups_rv,
    contrasts = contrasts_rv,
    configList = config_list,
    resultRv = result_rv,
    registerInputSyncObserverShells = function(input, session, designMatrix, removedSamples, factors, groups) {
      expect_identical(input, input_values)
      expect_identical(session, session_obj)
      expect_identical(designMatrix, design_matrix)
      expect_identical(removedSamples, removed_samples)
      expect_identical(factors, factors_rv)
      expect_identical(groups, groups_rv)
      call_log <<- c(call_log, "inputSync")
    },
    registerRenderOutputShells = function(input, output, session, proxyDataTable, designMatrix, removedSamples, factors, contrasts) {
      expect_identical(input, input_values)
      expect_identical(output, output_env)
      expect_identical(session, session_obj)
      expect_identical(proxyDataTable, proxy_data_table)
      expect_identical(designMatrix, design_matrix)
      expect_identical(removedSamples, removed_samples)
      expect_identical(factors, factors_rv)
      expect_identical(contrasts, contrasts_rv)
      call_log <<- c(call_log, "renderOutputs")
    },
    registerEventObserverShells = function(input, session, initialState, designMatrix, dataClnReactive, removedSamples, factors, groups, contrasts, configList, resultRv) {
      expect_identical(input, input_values)
      expect_identical(session, session_obj)
      expect_identical(initialState, initial_state)
      expect_identical(designMatrix, design_matrix)
      expect_identical(dataClnReactive, data_cln_reactive)
      expect_identical(removedSamples, removed_samples)
      expect_identical(factors, factors_rv)
      expect_identical(groups, groups_rv)
      expect_identical(contrasts, contrasts_rv)
      expect_identical(configList, config_list)
      expect_identical(resultRv, result_rv)
      call_log <<- c(call_log, "eventShells")
    }
  )

  expect_null(result)
  expect_identical(call_log, c("inputSync", "renderOutputs", "eventShells"))
  expect_identical(initial_state_calls, 0L)
  expect_identical(design_matrix_calls, 0L)
  expect_identical(data_cln_calls, 0L)
  expect_identical(removed_samples_calls, 0L)
  expect_identical(factors_calls, 0L)
  expect_identical(groups_calls, 0L)
  expect_identical(contrasts_calls, 0L)
  expect_identical(config_list_calls, 0L)
  expect_identical(result_rv_calls, 0L)
})

test_that("createProtDesignInitialStateReactive defers the builder initial-state bootstrap until the reactive is invoked", {
  data_tbl_calls <- 0L
  config_list_calls <- 0L
  column_mapping_calls <- 0L
  build_initial_state_calls <- 0L
  req_calls <- 0L
  reactive_expr <- NULL
  reactive_env <- NULL

  data_tbl <- function(...) {
    data_tbl_calls <<- data_tbl_calls + 1L
    data.frame(Run = "S1", stringsAsFactors = FALSE)
  }
  config_list <- function(...) {
    config_list_calls <<- config_list_calls + 1L
    list(globalParameters = list(workflow = "design"))
  }
  column_mapping <- function(...) {
    column_mapping_calls <<- column_mapping_calls + 1L
    list(sample_id = "Run")
  }

  initial_state <- createProtDesignInitialStateReactive(
    dataTbl = data_tbl,
    configList = config_list,
    columnMapping = column_mapping,
    buildInitialState = function(dataTbl, configList, columnMapping) {
      build_initial_state_calls <<- build_initial_state_calls + 1L
      list(
        dataTbl = dataTbl,
        configList = configList,
        columnMapping = columnMapping
      )
    },
    reactiveFn = function(expr) {
      reactive_expr <<- substitute(expr)
      reactive_env <<- parent.frame()
      function() eval(reactive_expr, envir = reactive_env)
    },
    reqFn = function(value) {
      req_calls <<- req_calls + 1L
      value
    }
  )

  expect_true(is.function(initial_state))
  expect_identical(data_tbl_calls, 0L)
  expect_identical(config_list_calls, 0L)
  expect_identical(column_mapping_calls, 0L)
  expect_identical(build_initial_state_calls, 0L)
  expect_identical(req_calls, 0L)

  result <- initial_state()

  expect_identical(
    result,
    list(
      dataTbl = data.frame(Run = "S1", stringsAsFactors = FALSE),
      configList = list(globalParameters = list(workflow = "design")),
      columnMapping = list(sample_id = "Run")
    )
  )
  expect_identical(data_tbl_calls, 2L)
  expect_identical(config_list_calls, 2L)
  expect_identical(column_mapping_calls, 1L)
  expect_identical(build_initial_state_calls, 1L)
  expect_identical(req_calls, 2L)
})

test_that("completeProtDesignBuilderServerBootstrap wires the remaining proxy bootstrap and registration fan-out", {
  call_log <- character()
  fake_input <- list(reset_builder = 1)
  fake_output <- structure(list(), class = "shinyoutput")
  fake_session <- list(ns = function(id) paste0("builder-shell-", id))
  proxy_seed <- structure(list(id = "data_table"), class = "prot_design_proxy")
  config_list <- function(...) list(globalParameters = list(workflow = "design"))

  builder_state <- list(
    initialState = function() "initial-state",
    resultRv = function() "builder-results",
    designMatrix = function() "design-matrix",
    dataClnReactive = function() "clean-data",
    groups = function() c("Control", "Treatment"),
    factors = function() c("Condition", "Batch"),
    contrasts = function() "contrasts",
    removedSamples = function() "Run2"
  )

  result <- completeProtDesignBuilderServerBootstrap(
    input = fake_input,
    output = fake_output,
    session = fake_session,
    builderState = builder_state,
    configList = config_list,
    createDataTableProxy = function() {
      call_log <<- c(call_log, "proxy")
      proxy_seed
    },
    registerBuilderServerShells = function(input, output, session, initialState, proxyDataTable, designMatrix, dataClnReactive, removedSamples, factors, groups, contrasts, configList, resultRv) {
      expect_identical(input, fake_input)
      expect_identical(output, fake_output)
      expect_identical(session, fake_session)
      expect_identical(initialState, builder_state$initialState)
      expect_identical(proxyDataTable, proxy_seed)
      expect_identical(designMatrix, builder_state$designMatrix)
      expect_identical(dataClnReactive, builder_state$dataClnReactive)
      expect_identical(removedSamples, builder_state$removedSamples)
      expect_identical(factors, builder_state$factors)
      expect_identical(groups, builder_state$groups)
      expect_identical(contrasts, builder_state$contrasts)
      expect_identical(configList, config_list)
      expect_identical(resultRv, builder_state$resultRv)
      call_log <<- c(call_log, "builderServerShells")
      invisible(NULL)
    }
  )

  expect_identical(result, builder_state$resultRv)
  expect_identical(call_log, c("proxy", "builderServerShells"))
})

test_that("runProtDesignBuilderModuleServerShell wires the remaining builder module-server shell in order", {
  call_log <- character()
  fake_input <- list(reset_builder = 1)
  fake_output <- structure(list(), class = "shinyoutput")
  fake_session <- list(ns = function(id) paste0("builder-shell-", id))
  data_tbl <- function(...) data.frame(Run = "S1", stringsAsFactors = FALSE)
  config_list <- function(...) list(globalParameters = list(workflow = "design"))
  column_mapping <- function(...) list(sample_id = "Run")

  builder_state <- list(
    initialState = function() "initial-state",
    resultRv = function() "builder-results",
    designMatrix = function() "design-matrix",
    dataClnReactive = function() "clean-data",
    groups = function() c("Control", "Treatment"),
    factors = function() c("Condition", "Batch"),
    contrasts = function() "contrasts",
    removedSamples = function() "Run2"
  )

  result <- runProtDesignBuilderModuleServerShell(
    input = fake_input,
    output = fake_output,
    session = fake_session,
    dataTbl = data_tbl,
    configList = config_list,
    columnMapping = column_mapping,
    initializeBuilderServerState = function(input, session, dataTbl, configList, columnMapping, buildInitialState, createInitialStateReactive, createMutableStateShells, registerInitialStateSyncObserver) {
      expect_identical(input, fake_input)
      expect_identical(session, fake_session)
      expect_identical(dataTbl, data_tbl)
      expect_identical(configList, config_list)
      expect_identical(columnMapping, column_mapping)
      expect_true(is.function(buildInitialState))
      expect_true(is.function(createInitialStateReactive))
      expect_true(is.function(createMutableStateShells))
      expect_true(is.function(registerInitialStateSyncObserver))
      call_log <<- c(call_log, "builderState")
      builder_state
    },
    completeBuilderServerBootstrap = function(input, output, session, builderState, configList, createDataTableProxy, registerBuilderServerShells) {
      expect_identical(input, fake_input)
      expect_identical(output, fake_output)
      expect_identical(session, fake_session)
      expect_identical(builderState, builder_state)
      expect_identical(configList, config_list)
      expect_true(is.function(createDataTableProxy))
      expect_true(is.function(registerBuilderServerShells))
      call_log <<- c(call_log, "builderBootstrap")
      builder_state$resultRv
    }
  )

  expect_identical(result, builder_state$resultRv)
  expect_identical(
    call_log,
    c("builderState", "builderBootstrap")
  )
})

test_that("runProtDesignBuilderServerEntryShell wires the remaining public builder moduleServer entry shell in order", {
  call_log <- character()
  data_tbl <- function() data.frame(Run = "S1", stringsAsFactors = FALSE)
  config_list <- function() list(globalParameters = list(workflow = "design"))
  column_mapping <- function() list(sample_id = "Run")
  fake_input <- list(reset_builder = 1)
  fake_output <- structure(list(), class = "shinyoutput")
  fake_session <- list(ns = function(id) paste0("builder-entry-", id))
  entry_id <- "builder-entry"

  result <- runProtDesignBuilderServerEntryShell(
    id = entry_id,
    dataTbl = data_tbl,
    configList = config_list,
    columnMapping = column_mapping,
    moduleServer = function(id, module, ...) {
      expect_identical(id, entry_id)
      call_log <<- c(call_log, "moduleServer")
      module(fake_input, fake_output, fake_session)
    },
    runBuilderModuleServerShell = function(input, output, session, dataTbl, configList, columnMapping) {
      expect_identical(input, fake_input)
      expect_identical(output, fake_output)
      expect_identical(session, fake_session)
      expect_identical(dataTbl, data_tbl)
      expect_identical(configList, config_list)
      expect_identical(columnMapping, column_mapping)
      call_log <<- c(call_log, "builderModuleServerShell")
      "builder-results"
    }
  )

  expect_identical(result, "builder-results")
  expect_identical(
    call_log,
    c("moduleServer", "builderModuleServerShell")
  )
})

test_that("mod_prot_design_builder_server wires the public builder entry seam in order", {
  call_log <- character()
  data_tbl <- function() data.frame(Run = "S1", stringsAsFactors = FALSE)
  config_list <- function() list(globalParameters = list(workflow = "design"))
  column_mapping <- function() list(sample_id = "Run")
  entry_id <- "builder-entry"

  local_mocked_bindings(
    runProtDesignBuilderServerEntryShell = function(id, dataTbl, configList, columnMapping, moduleServer = shiny::moduleServer, runBuilderModuleServerShell = runProtDesignBuilderModuleServerShell) {
      expect_identical(id, entry_id)
      expect_identical(dataTbl, data_tbl)
      expect_identical(configList, config_list)
      expect_identical(columnMapping, column_mapping)
      expect_identical(moduleServer, shiny::moduleServer)
      expect_identical(runBuilderModuleServerShell, runProtDesignBuilderModuleServerShell)
      call_log <<- c(call_log, "builderEntryShell")
      "builder-results"
    },
    .env = environment(mod_prot_design_builder_server)
  )

  result <- mod_prot_design_builder_server(
    entry_id,
    data_tbl = data_tbl,
    config_list = config_list,
    column_mapping = column_mapping
  )

  expect_identical(result, "builder-results")
  expect_identical(call_log, "builderEntryShell")
})

test_that("mod_prot_design_builder_server passes callback seeds opaquely through the public builder entry seam", {
  call_log <- character()
  data_tbl_calls <- 0L
  config_list_calls <- 0L
  column_mapping_calls <- 0L
  entry_id <- "builder-entry"

  data_tbl <- function(...) {
    data_tbl_calls <<- data_tbl_calls + 1L
    data.frame(Run = "S1", stringsAsFactors = FALSE)
  }
  config_list <- function(...) {
    config_list_calls <<- config_list_calls + 1L
    list(globalParameters = list(workflow = "design"))
  }
  column_mapping <- function(...) {
    column_mapping_calls <<- column_mapping_calls + 1L
    list(sample_id = "Run")
  }

  local_mocked_bindings(
    runProtDesignBuilderServerEntryShell = function(id, dataTbl, configList, columnMapping, moduleServer = shiny::moduleServer, runBuilderModuleServerShell = runProtDesignBuilderModuleServerShell) {
      expect_identical(id, entry_id)
      expect_identical(dataTbl, data_tbl)
      expect_identical(configList, config_list)
      expect_identical(columnMapping, column_mapping)
      expect_identical(moduleServer, shiny::moduleServer)
      expect_identical(runBuilderModuleServerShell, runProtDesignBuilderModuleServerShell)
      call_log <<- c(call_log, "builderEntryShell")
      "builder-results"
    },
    .env = environment(mod_prot_design_builder_server)
  )

  result <- mod_prot_design_builder_server(
    entry_id,
    data_tbl = data_tbl,
    config_list = config_list,
    column_mapping = column_mapping
  )

  expect_identical(result, "builder-results")
  expect_identical(
    call_log,
    "builderEntryShell"
  )
  expect_identical(data_tbl_calls, 0L)
  expect_identical(config_list_calls, 0L)
  expect_identical(column_mapping_calls, 0L)
})

test_that("registerProtDesignRenameObserverShells wires the rename registration fan-out", {
  call_log <- list()
  input_values <- list(
    rename_sample = 1,
    bulk_rename = 2
  )
  design_matrix <- function() "design-matrix"
  data_cln_reactive <- function() "clean-data"
  session_obj <- list(ns = function(id) paste0("builder-", id))

  result <- registerProtDesignRenameObserverShells(
    input = input_values,
    designMatrix = design_matrix,
    dataCln = data_cln_reactive,
    session = session_obj,
    registerRenameSampleObserver = function(input, designMatrix, dataCln, session) {
      call_log[[length(call_log) + 1]] <<- list(
        kind = "renameSample",
        input = identical(input, input_values),
        designMatrix = identical(designMatrix, design_matrix),
        dataCln = identical(dataCln, data_cln_reactive),
        session = identical(session, session_obj)
      )
    },
    registerBulkRenameObserver = function(input, designMatrix, dataCln) {
      call_log[[length(call_log) + 1]] <<- list(
        kind = "bulkRename",
        input = identical(input, input_values),
        designMatrix = identical(designMatrix, design_matrix),
        dataCln = identical(dataCln, data_cln_reactive)
      )
    }
  )

  expect_null(result)
  expect_identical(call_log, list(
    list(
      kind = "renameSample",
      input = TRUE,
      designMatrix = TRUE,
      dataCln = TRUE,
      session = TRUE
    ),
    list(
      kind = "bulkRename",
      input = TRUE,
      designMatrix = TRUE,
      dataCln = TRUE
    )
  ))
})

test_that("registerProtDesignFactorMetadataObserverShells wires the factor and metadata registration fan-out", {
  call_log <- list()
  input_values <- list(
    add_factor = 1,
    assign_metadata = 2
  )
  factors_rv <- function() c("Condition", "Batch")
  design_matrix <- function() "design-matrix"
  groups_rv <- function(value) invisible(value)
  session_obj <- list(ns = function(id) paste0("builder-", id))

  result <- registerProtDesignFactorMetadataObserverShells(
    input = input_values,
    factors = factors_rv,
    designMatrix = design_matrix,
    groups = groups_rv,
    session = session_obj,
    registerAddFactorObserver = function(input, factors, session) {
      call_log[[length(call_log) + 1]] <<- list(
        kind = "addFactor",
        input = identical(input, input_values),
        factors = identical(factors, factors_rv),
        session = identical(session, session_obj)
      )
    },
    registerAssignMetadataObserver = function(input, designMatrix, groups) {
      call_log[[length(call_log) + 1]] <<- list(
        kind = "assignMetadata",
        input = identical(input, input_values),
        designMatrix = identical(designMatrix, design_matrix),
        groups = identical(groups, groups_rv)
      )
    }
  )

  expect_null(result)
  expect_identical(call_log, list(
    list(
      kind = "addFactor",
      input = TRUE,
      factors = TRUE,
      session = TRUE
    ),
    list(
      kind = "assignMetadata",
      input = TRUE,
      designMatrix = TRUE,
      groups = TRUE
    )
  ))
})

test_that("registerProtDesignActionObserverShells wires the action observer registration fan-out", {
  call_log <- list()
  input_values <- list(
    assign_tech_reps = 1,
    add_contrast = 2,
    remove_samples = 3
  )
  design_matrix <- function() "design-matrix"
  contrasts_rv <- function() "contrasts"
  removed_samples <- function() "Run2"
  session_obj <- list(ns = function(id) paste0("builder-", id))

  result <- registerProtDesignActionObserverShells(
    input = input_values,
    designMatrix = design_matrix,
    contrasts = contrasts_rv,
    removedSamples = removed_samples,
    session = session_obj,
    registerAssignTechRepsObserver = function(input, designMatrix) {
      call_log[[length(call_log) + 1]] <<- list(
        kind = "assignTechReps",
        input = identical(input, input_values),
        designMatrix = identical(designMatrix, design_matrix)
      )
    },
    registerAddContrastObserver = function(input, contrasts) {
      call_log[[length(call_log) + 1]] <<- list(
        kind = "addContrast",
        input = identical(input, input_values),
        contrasts = identical(contrasts, contrasts_rv)
      )
    },
    registerRemoveSamplesObserver = function(input, removedSamples, session) {
      call_log[[length(call_log) + 1]] <<- list(
        kind = "removeSamples",
        input = identical(input, input_values),
        removedSamples = identical(removedSamples, removed_samples),
        session = identical(session, session_obj)
      )
    }
  )

  expect_null(result)
  expect_identical(call_log, list(
    list(
      kind = "assignTechReps",
      input = TRUE,
      designMatrix = TRUE
    ),
    list(
      kind = "addContrast",
      input = TRUE,
      contrasts = TRUE
    ),
    list(
      kind = "removeSamples",
      input = TRUE,
      removedSamples = TRUE,
      session = TRUE
    )
  ))
})

test_that("registerProtDesignReplicateInputsOutput wires the replicate-input render shell", {
  call_log <- list()
  output <- new.env(parent = emptyenv())
  input <- list(selected_runs = c("Run1", "Run2", "Run3"))

  render_value <- registerProtDesignReplicateInputsOutput(
    output = output,
    input = input,
    nsFn = function(id) paste0("builder-", id),
    renderUIFn = function(expr) {
      expr_sub <- substitute(expr)
      expr_env <- parent.frame()
      call_log[[length(call_log) + 1]] <<- list(kind = "renderUI")
      rendered <- eval(expr_sub, expr_env)
      call_log[[length(call_log) + 1]] <<- list(kind = "rendered", value = rendered)
      "rendered-ui"
    },
    reqFn = function(value) {
      call_log[[length(call_log) + 1]] <<- list(kind = "req", value = value)
      invisible(NULL)
    },
    buildReplicateInputsFn = function(selectedRuns, nsFn) {
      call_log[[length(call_log) + 1]] <<- list(
        kind = "build",
        selectedRuns = selectedRuns,
        namespacedInput = nsFn("replicate_start")
      )
      "replicate-ui"
    }
  )

  expect_identical(render_value, "rendered-ui")
  expect_identical(output$replicate_inputs, "rendered-ui")
  expect_identical(call_log, list(
    list(kind = "renderUI"),
    list(kind = "req", value = c("Run1", "Run2", "Run3")),
    list(
      kind = "build",
      selectedRuns = c("Run1", "Run2", "Run3"),
      namespacedInput = "builder-replicate_start"
    ),
    list(kind = "rendered", value = "replicate-ui")
  ))
})

test_that("buildProtDesignDefinedContrastsDisplay reports when no contrasts are defined", {
  display <- buildProtDesignDefinedContrastsDisplay(
    contrastData = NULL,
    formulaString = "~ group"
  )

  expect_equal(
    as.character(htmltools::renderTags(display)$html),
    "<p>No contrasts defined yet.</p>"
  )
})

test_that("buildProtDesignDefinedContrastsDisplay formats friendly names and group-prefixed contrasts", {
  contrast_data <- data.frame(
    contrast_name = c("GA.Control.vs.GA.Elevated", "GB.Control.vs.GB.Elevated"),
    numerator = c("GA_Control", "GB_Control"),
    denominator = c("GA_Elevated", "GB_Elevated"),
    stringsAsFactors = FALSE
  )

  display <- buildProtDesignDefinedContrastsDisplay(
    contrastData = contrast_data,
    formulaString = "~ 0 + group"
  )
  html <- htmltools::renderTags(display)$html

  expect_match(
    html,
    "GA_Control_vs_GA_Elevated=groupGA_Control-groupGA_Elevated",
    fixed = TRUE
  )
  expect_match(
    html,
    "GB_Control_vs_GB_Elevated=groupGB_Control-groupGB_Elevated",
    fixed = TRUE
  )
})

test_that("buildProtDesignAvailableFactorsDisplay reports when no factors are defined", {
  display <- buildProtDesignAvailableFactorsDisplay(character(0))

  expect_equal(
    as.character(htmltools::renderTags(display)$html),
    "<p>No factors defined yet (use the 'Factors' tab).</p>"
  )
})

test_that("buildProtDesignAvailableFactorsDisplay renders a comma-separated factor list", {
  display <- buildProtDesignAvailableFactorsDisplay(c("Condition", "Batch", "Timepoint"))

  expect_equal(
    as.character(htmltools::renderTags(display)$html),
    "<p>Condition, Batch, Timepoint</p>"
  )
})

test_that("applyProtDesignFactorAppendReset appends a trimmed unique factor and clears the input", {
  updated <- applyProtDesignFactorAppendReset(
    currentFactors = c("Condition", "Batch"),
    newFactorInput = "  Timepoint  "
  )

  expect_identical(updated$factors, c("Condition", "Batch", "Timepoint"))
  expect_identical(updated$newFactorValue, "")
})

test_that("applyProtDesignFactorAppendReset ignores blanks and duplicates but still clears the input", {
  blank_update <- applyProtDesignFactorAppendReset(
    currentFactors = c("Condition", "Batch"),
    newFactorInput = "   "
  )
  duplicate_update <- applyProtDesignFactorAppendReset(
    currentFactors = c("Condition", "Batch"),
    newFactorInput = "Batch"
  )

  expect_identical(blank_update$factors, c("Condition", "Batch"))
  expect_identical(blank_update$newFactorValue, "")
  expect_identical(duplicate_update$factors, c("Condition", "Batch"))
  expect_identical(duplicate_update$newFactorValue, "")
})

test_that("registerProtDesignAddFactorObserver registers the add-factor handoff shell", {
  call_log <- list()
  current_factors <- c("Condition", "Batch")
  updated_factors <- c("Condition", "Batch", "Timepoint")
  input <- list(
    add_factor = 1,
    new_factor = "  Timepoint  "
  )

  factors_rv <- function(value) {
    if (missing(value)) {
      call_log[[length(call_log) + 1]] <<- list(kind = "get")
      current_factors
    } else {
      call_log[[length(call_log) + 1]] <<- list(kind = "set", value = value)
      current_factors <<- value
      invisible(NULL)
    }
  }

  observer <- registerProtDesignAddFactorObserver(
    input = input,
    factors = factors_rv,
    session = "builder-session",
    observeEventFn = function(eventExpr, handlerExpr, ignoreNULL = FALSE) {
      event_sub <- substitute(eventExpr)
      handler_sub <- substitute(handlerExpr)
      event_value <- eval(event_sub, parent.frame())

      expect_identical(event_value, 1)
      expect_false(ignoreNULL)
      call_log[[length(call_log) + 1]] <<- list(kind = "observe")

      eval(handler_sub, parent.frame())
      "observer-registered"
    },
    reqFn = function(value) {
      expect_identical(value, "  Timepoint  ")
      call_log[[length(call_log) + 1]] <<- list(kind = "req", value = value)
      invisible(NULL)
    },
    applyFactorAppendResetFn = function(currentFactors, newFactorInput) {
      expect_identical(currentFactors, c("Condition", "Batch"))
      expect_identical(newFactorInput, "  Timepoint  ")
      call_log[[length(call_log) + 1]] <<- list(
        kind = "apply",
        currentFactors = currentFactors,
        newFactorInput = newFactorInput
      )
      list(
        factors = updated_factors,
        newFactorValue = ""
      )
    },
    updateTextInputFn = function(session, inputId, value = NULL) {
      call_log[[length(call_log) + 1]] <<- list(
        kind = "updateTextInput",
        session = session,
        inputId = inputId,
        value = value
      )
      invisible(NULL)
    }
  )

  expect_identical(observer, "observer-registered")
  expect_identical(current_factors, updated_factors)
  expect_identical(call_log, list(
    list(kind = "observe"),
    list(kind = "req", value = "  Timepoint  "),
    list(kind = "get"),
    list(
      kind = "apply",
      currentFactors = c("Condition", "Batch"),
      newFactorInput = "  Timepoint  "
    ),
    list(kind = "set", value = c("Condition", "Batch", "Timepoint")),
    list(
      kind = "updateTextInput",
      session = "builder-session",
      inputId = "new_factor",
      value = ""
    )
  ))
})

test_that("registerProtDesignAssignMetadataObserver registers the metadata-assignment handoff shell", {
  call_log <- list()
  current_design_matrix <- data.frame(
    Run = c("Run1", "Run2", "Run3"),
    factor1 = c(NA_character_, "Control", NA_character_),
    factor2 = c(NA_character_, "BatchA", NA_character_),
    factor3 = c(NA_character_, NA_character_, NA_character_),
    replicates = c(NA_integer_, 2L, NA_integer_),
    group = c(NA_character_, "Control_BatchA", NA_character_),
    stringsAsFactors = FALSE
  )
  input <- list(
    assign_metadata = 1,
    selected_runs = c("Run1", "Run3"),
    factor1_select = "Treatment",
    factor2_select = "Day1",
    factor3_select = "",
    replicate_start = 5L
  )

  design_matrix_rv <- function(value) {
    if (missing(value)) {
      call_log[[length(call_log) + 1]] <<- list(kind = "getDesignMatrix")
      current_design_matrix
    } else {
      call_log[[length(call_log) + 1]] <<- list(kind = "setDesignMatrix", value = value)
      current_design_matrix <<- value
      invisible(NULL)
    }
  }
  groups_rv <- function(value) {
    call_log[[length(call_log) + 1]] <<- list(kind = "setGroups", value = value)
    invisible(NULL)
  }

  observer <- registerProtDesignAssignMetadataObserver(
    input = input,
    designMatrix = design_matrix_rv,
    groups = groups_rv,
    observeEventFn = function(eventExpr, handlerExpr, ignoreNULL = FALSE) {
      event_sub <- substitute(eventExpr)
      handler_sub <- substitute(handlerExpr)
      event_value <- eval(event_sub, parent.frame())

      expect_identical(event_value, 1)
      expect_false(ignoreNULL)
      call_log[[length(call_log) + 1]] <<- list(kind = "observe")

      eval(handler_sub, parent.frame())
      "observer-registered"
    },
    reqFn = function(...) {
      expect_identical(list(...), list(c("Run1", "Run3"), "Treatment"))
      call_log[[length(call_log) + 1]] <<- list(kind = "req", value = list(...))
      invisible(NULL)
    },
    seqFn = function(from, length.out) {
      expect_identical(from, 5L)
      expect_identical(length.out, 2L)
      call_log[[length(call_log) + 1]] <<- list(kind = "seq", from = from, length.out = length.out)
      c(5L, 6L)
    }
  )

  expect_identical(observer, "observer-registered")
  expect_identical(
    current_design_matrix,
    data.frame(
      Run = c("Run1", "Run2", "Run3"),
      factor1 = c("Treatment", "Control", "Treatment"),
      factor2 = c("Day1", "BatchA", "Day1"),
      factor3 = c(NA_character_, NA_character_, NA_character_),
      replicates = c(5L, 2L, 6L),
      group = c("Treatment_Day1", "Control_BatchA", "Treatment_Day1"),
      stringsAsFactors = FALSE
    )
  )
  expect_identical(call_log, list(
    list(kind = "observe"),
    list(kind = "req", value = list(c("Run1", "Run3"), "Treatment")),
    list(kind = "getDesignMatrix"),
    list(kind = "seq", from = 5L, length.out = 2L),
    list(
      kind = "setDesignMatrix",
      value = data.frame(
        Run = c("Run1", "Run2", "Run3"),
        factor1 = c("Treatment", "Control", "Treatment"),
        factor2 = c("Day1", "BatchA", "Day1"),
        factor3 = c(NA_character_, NA_character_, NA_character_),
        replicates = c(5L, 2L, 6L),
        group = c("Treatment_Day1", "Control_BatchA", "Treatment_Day1"),
        stringsAsFactors = FALSE
      )
    ),
    list(kind = "setGroups", value = c("Treatment_Day1", "Control_BatchA"))
  ))
})

test_that("registerProtDesignAssignTechRepsObserver registers the tech-replicate handoff shell", {
  call_log <- list()
  current_design_matrix <- data.frame(
    Run = c("Run1", "Run2", "Run3", "Run4"),
    group = c("GroupA", "GroupA", "GroupA", "GroupB"),
    replicates = c(1L, 2L, 3L, 1L),
    tech_reps = c(NA_integer_, NA_integer_, NA_integer_, NA_integer_),
    stringsAsFactors = FALSE
  )
  input <- list(
    assign_tech_reps = 1,
    tech_rep_samples = c("Run1", "Run2"),
    tech_rep_assignment_mode = "lowest",
    manual_replicate_number = 9L
  )

  design_matrix_rv <- function(value) {
    if (missing(value)) {
      call_log[[length(call_log) + 1]] <<- list(kind = "getDesignMatrix")
      current_design_matrix
    } else {
      call_log[[length(call_log) + 1]] <<- list(kind = "setDesignMatrix", value = value)
      current_design_matrix <<- value
      invisible(NULL)
    }
  }

  observer <- registerProtDesignAssignTechRepsObserver(
    input = input,
    designMatrix = design_matrix_rv,
    observeEventFn = function(eventExpr, handlerExpr, ignoreNULL = FALSE) {
      event_sub <- substitute(eventExpr)
      handler_sub <- substitute(handlerExpr)
      event_value <- eval(event_sub, parent.frame())

      expect_identical(event_value, 1)
      expect_false(ignoreNULL)
      call_log[[length(call_log) + 1]] <<- list(kind = "observe")

      eval(handler_sub, parent.frame())
      "observer-registered"
    },
    reqFn = function(...) {
      expect_identical(list(...), list(c("Run1", "Run2")))
      call_log[[length(call_log) + 1]] <<- list(kind = "req", value = list(...))
      invisible(NULL)
    },
    showNotificationFn = function(message, type) {
      call_log[[length(call_log) + 1]] <<- list(kind = "showNotification", message = message, type = type)
      invisible(NULL)
    }
  )

  expect_identical(observer, "observer-registered")
  expect_identical(
    current_design_matrix,
    data.frame(
      Run = c("Run1", "Run2", "Run3", "Run4"),
      group = c("GroupA", "GroupA", "GroupA", "GroupB"),
      replicates = c(1, 1, 2, 1),
      tech_reps = c(1L, 2L, NA_integer_, NA_integer_),
      stringsAsFactors = FALSE
    )
  )
  expect_identical(call_log, list(
    list(kind = "observe"),
    list(kind = "req", value = list(c("Run1", "Run2"))),
    list(kind = "getDesignMatrix"),
    list(
      kind = "setDesignMatrix",
      value = data.frame(
        Run = c("Run1", "Run2", "Run3", "Run4"),
        group = c("GroupA", "GroupA", "GroupA", "GroupB"),
        replicates = c(1, 1, 2, 1),
        tech_reps = c(1L, 2L, NA_integer_, NA_integer_),
        stringsAsFactors = FALSE
      )
    ),
    list(
      kind = "showNotification",
      message = "Assigned 2 samples as technical replicates with biological replicate number 1 and adjusted subsequent replicate numbers",
      type = "message"
    )
  ))
})

test_that("applyProtDesignContrastAppend adds a new unique contrast row", {
  updated <- applyProtDesignContrastAppend(
    currentContrasts = data.frame(
      contrast_name = "GA_Control.vs.GA_Elevated",
      numerator = "GA_Control",
      denominator = "GA_Elevated",
      stringsAsFactors = FALSE
    ),
    group1 = "GB_Control",
    group2 = "GB_Elevated"
  )

  expect_identical(
    updated$contrast_name,
    c("GA_Control.vs.GA_Elevated", "GB_Control.vs.GB_Elevated")
  )
  expect_identical(updated$numerator, c("GA_Control", "GB_Control"))
  expect_identical(updated$denominator, c("GA_Elevated", "GB_Elevated"))
})

test_that("applyProtDesignContrastAppend ignores blank, duplicate, and self contrasts", {
  existing <- data.frame(
    contrast_name = "GA_Control.vs.GA_Elevated",
    numerator = "GA_Control",
    denominator = "GA_Elevated",
    stringsAsFactors = FALSE
  )

  blank_update <- applyProtDesignContrastAppend(
    currentContrasts = existing,
    group1 = "",
    group2 = "GB_Elevated"
  )
  duplicate_update <- applyProtDesignContrastAppend(
    currentContrasts = existing,
    group1 = "GA_Control",
    group2 = "GA_Elevated"
  )
  self_update <- applyProtDesignContrastAppend(
    currentContrasts = existing,
    group1 = "GA_Control",
    group2 = "GA_Control"
  )

  expect_identical(blank_update, existing)
  expect_identical(duplicate_update, existing)
  expect_identical(self_update, existing)
})

test_that("registerProtDesignAddContrastObserver registers the add-contrast handoff shell", {
  call_log <- character()
  current_contrasts <- data.frame(
    contrast_name = "GA_Control.vs.GA_Elevated",
    numerator = "GA_Control",
    denominator = "GA_Elevated",
    stringsAsFactors = FALSE
  )
  appended_contrasts <- data.frame(
    contrast_name = c("GA_Control.vs.GA_Elevated", "GB_Control.vs.GB_Elevated"),
    numerator = c("GA_Control", "GB_Control"),
    denominator = c("GA_Elevated", "GB_Elevated"),
    stringsAsFactors = FALSE
  )
  input <- list(
    add_contrast = 1,
    contrast_group1 = "GB_Control",
    contrast_group2 = "GB_Elevated"
  )

  contrasts_rv <- function(value) {
    if (missing(value)) {
      current_contrasts
    } else {
      call_log <<- c(call_log, "set")
      current_contrasts <<- value
      invisible(NULL)
    }
  }

  observer <- registerProtDesignAddContrastObserver(
    input = input,
    contrasts = contrasts_rv,
    observeEventFn = function(eventExpr, handlerExpr, ignoreNULL = FALSE) {
      event_sub <- substitute(eventExpr)
      handler_sub <- substitute(handlerExpr)
      event_value <- eval(event_sub, parent.frame())

      expect_identical(event_value, 1)
      expect_false(ignoreNULL)
      call_log <<- c(call_log, "observe")

      eval(handler_sub, parent.frame())
      "observer-registered"
    },
    reqFn = function(...) {
      expect_identical(list(...), list("GB_Control", "GB_Elevated"))
      call_log <<- c(call_log, "req")
      invisible(NULL)
    },
    appendContrastFn = function(currentContrasts, group1, group2) {
      expect_identical(currentContrasts, current_contrasts)
      expect_identical(group1, "GB_Control")
      expect_identical(group2, "GB_Elevated")
      call_log <<- c(call_log, "append")
      appended_contrasts
    }
  )

  expect_identical(observer, "observer-registered")
  expect_identical(current_contrasts, appended_contrasts)
  expect_identical(call_log, c("observe", "req", "append", "set"))
})

test_that("registerProtDesignRemoveSamplesObserver registers the remove-samples handoff shell", {
  call_log <- list()
  current_removed_samples <- c("Run0")
  updated_removed_samples <- c("Run0", "Run2", "Run5")
  input <- list(
    remove_samples = 1,
    samples_to_remove = c("Run2", "Run5")
  )

  removed_samples <- function(value) {
    if (missing(value)) {
      call_log[[length(call_log) + 1]] <<- list(kind = "get")
      current_removed_samples
    } else {
      call_log[[length(call_log) + 1]] <<- list(kind = "set", value = value)
      current_removed_samples <<- value
      invisible(NULL)
    }
  }

  observer <- registerProtDesignRemoveSamplesObserver(
    input = input,
    removedSamples = removed_samples,
    session = "builder-session",
    observeEventFn = function(eventExpr, handlerExpr, ignoreNULL = FALSE) {
      event_sub <- substitute(eventExpr)
      handler_sub <- substitute(handlerExpr)
      event_value <- eval(event_sub, parent.frame())

      expect_identical(event_value, 1)
      expect_false(ignoreNULL)
      call_log[[length(call_log) + 1]] <<- list(kind = "observe")

      eval(handler_sub, parent.frame())
      "observer-registered"
    },
    reqFn = function(value) {
      expect_identical(value, c("Run2", "Run5"))
      call_log[[length(call_log) + 1]] <<- list(kind = "req", value = value)
      invisible(NULL)
    },
    applyRemovedSamplesUpdateFn = function(currentRemovedSamples, selectedSamples) {
      expect_identical(currentRemovedSamples, c("Run0"))
      expect_identical(selectedSamples, c("Run2", "Run5"))
      call_log[[length(call_log) + 1]] <<- list(
        kind = "apply",
        currentRemovedSamples = currentRemovedSamples,
        selectedSamples = selectedSamples
      )
      updated_removed_samples
    },
    updateSelectizeInputFn = function(session, inputId, selected = NULL) {
      call_log[[length(call_log) + 1]] <<- list(
        kind = "updateSelectizeInput",
        session = session,
        inputId = inputId,
        selected = selected
      )
      invisible(NULL)
    },
    showNotificationFn = function(message, type = NULL) {
      call_log[[length(call_log) + 1]] <<- list(
        kind = "showNotification",
        message = message,
        type = type
      )
      invisible(NULL)
    }
  )

  expect_identical(observer, "observer-registered")
  expect_identical(current_removed_samples, updated_removed_samples)
  expect_identical(call_log, list(
    list(kind = "observe"),
    list(kind = "req", value = c("Run2", "Run5")),
    list(kind = "get"),
    list(
      kind = "apply",
      currentRemovedSamples = c("Run0"),
      selectedSamples = c("Run2", "Run5")
    ),
    list(kind = "set", value = c("Run0", "Run2", "Run5")),
    list(
      kind = "updateSelectizeInput",
      session = "builder-session",
      inputId = "samples_to_remove",
      selected = ""
    ),
    list(
      kind = "showNotification",
      message = "Removed 2 sample(s) from analysis.",
      type = "message"
    )
  ))
})

test_that("buildProtDesignSaveResultsContrastsTable returns NULL when no contrasts are defined", {
  expect_null(
    buildProtDesignSaveResultsContrastsTable(
      contrastData = NULL,
      formulaString = "~ group"
    )
  )
})

test_that("buildProtDesignSaveResultsContrastsTable formats friendly names and contrast strings for save results", {
  contrast_data <- data.frame(
    contrast_name = c("GA.Control.vs.GA.Elevated", "GB.Control.vs.GB.Elevated"),
    numerator = c("GA_Control", "GB_Control"),
    denominator = c("GA_Elevated", "GB_Elevated"),
    stringsAsFactors = FALSE
  )

  expect_identical(
    buildProtDesignSaveResultsContrastsTable(
      contrastData = contrast_data,
      formulaString = "~ 0 + group"
    ),
    data.frame(
      contrasts = c(
        "groupGA_Control-groupGA_Elevated",
        "groupGB_Control-groupGB_Elevated"
      ),
      friendly_names = c(
        "GA_Control_vs_GA_Elevated",
        "GB_Control_vs_GB_Elevated"
      ),
      full_format = c(
        "GA_Control_vs_GA_Elevated=groupGA_Control-groupGA_Elevated",
        "GB_Control_vs_GB_Elevated=groupGB_Control-groupGB_Elevated"
      ),
      stringsAsFactors = FALSE
    )
  )
})

test_that("buildProtDesignSaveResultsPayload returns NULL when no assigned samples remain", {
  design_matrix <- data.frame(
    Run = c("S1", "S2"),
    group = c(NA_character_, ""),
    replicates = c(1L, 2L),
    stringsAsFactors = FALSE
  )

  expect_null(
    buildProtDesignSaveResultsPayload(
      designMatrix = design_matrix,
      currentRemovedSamples = character(0),
      dataCln = data.frame(Run = c("S1", "S2"), value = c(10, 20), stringsAsFactors = FALSE),
      contrastData = NULL,
      configList = list(deAnalysisParameters = list(formula_string = "~ old")),
      formulaString = "~ group"
    )
  )
})

test_that("buildProtDesignSaveResultsPayload filters removed samples and updates final save payload", {
  design_matrix <- data.frame(
    Run = c("S1", "S2", "S3"),
    group = c("GA", "GB", "GC"),
    replicates = c(1L, 2L, 3L),
    stringsAsFactors = FALSE
  )
  data_cln <- data.frame(
    Run = c("S1", "S1", "S2", "S3"),
    value = c(10, 11, 20, 30),
    stringsAsFactors = FALSE
  )
  contrast_data <- data.frame(
    contrast_name = "GA.vs.GC",
    numerator = "GA",
    denominator = "GC",
    stringsAsFactors = FALSE
  )
  config_list <- list(
    deAnalysisParameters = list(formula_string = "~ old")
  )

  result <- buildProtDesignSaveResultsPayload(
    designMatrix = design_matrix,
    currentRemovedSamples = "S2",
    dataCln = data_cln,
    contrastData = contrast_data,
    configList = config_list,
    formulaString = "~ 0 + group"
  )

  expect_identical(result$design_matrix$Run, c("S1", "S3"))
  expect_identical(result$design_matrix$tech_rep_group, c("GA_1", "GC_3"))
  expect_identical(result$data_cln$Run, c("S1", "S1", "S3"))
  expect_identical(result$config_list$deAnalysisParameters$formula_string, "~ 0 + group")
  expect_identical(
    result$contrasts_tbl,
    data.frame(
      contrasts = "groupGA-groupGC",
      friendly_names = "GA_vs_GC",
      full_format = "GA_vs_GC=groupGA-groupGC",
      stringsAsFactors = FALSE
    )
  )
})

test_that("runProtDesignSaveResultsObserverShell warns when no assigned samples remain", {
  notifications <- list()
  assigned_result <- "unset"

  result <- runProtDesignSaveResultsObserverShell(
    designMatrix = data.frame(Run = c("S1", "S2"), stringsAsFactors = FALSE),
    currentRemovedSamples = character(0),
    dataCln = data.frame(Run = c("S1", "S2"), stringsAsFactors = FALSE),
    contrastData = NULL,
    configList = list(deAnalysisParameters = list(formula_string = "~ old")),
    formulaString = "~ group",
    resultSetter = function(value) {
      assigned_result <<- value
    },
    buildSaveResultsPayload = function(...) {
      NULL
    },
    showNotification = function(message, type = NULL, duration = NULL) {
      notifications[[length(notifications) + 1]] <<- list(
        message = message,
        type = type,
        duration = duration
      )
      invisible(NULL)
    }
  )

  expect_null(result)
  expect_identical(assigned_result, "unset")
  expect_identical(
    notifications,
    list(list(
      message = "No samples have been assigned to groups. Please assign metadata before saving.",
      type = "warning",
      duration = NULL
    ))
  )
})

test_that("runProtDesignSaveResultsObserverShell stores the result payload and shows success", {
  notifications <- list()
  assigned_result <- NULL
  expected_result <- list(
    design_matrix = data.frame(Run = "S1", group = "GA", stringsAsFactors = FALSE),
    data_cln = data.frame(Run = "S1", value = 10, stringsAsFactors = FALSE)
  )

  result <- runProtDesignSaveResultsObserverShell(
    designMatrix = data.frame(Run = "S1", stringsAsFactors = FALSE),
    currentRemovedSamples = "S2",
    dataCln = data.frame(Run = c("S1", "S2"), stringsAsFactors = FALSE),
    contrastData = NULL,
    configList = list(deAnalysisParameters = list(formula_string = "~ old")),
    formulaString = "~ 0 + group",
    resultSetter = function(value) {
      assigned_result <<- value
    },
    buildSaveResultsPayload = function(
        designMatrix,
        currentRemovedSamples,
        dataCln,
        contrastData,
        configList,
        formulaString
    ) {
      expect_identical(designMatrix$Run, "S1")
      expect_identical(currentRemovedSamples, "S2")
      expect_identical(dataCln$Run, c("S1", "S2"))
      expect_null(contrastData)
      expect_identical(configList$deAnalysisParameters$formula_string, "~ old")
      expect_identical(formulaString, "~ 0 + group")
      expected_result
    },
    showNotification = function(message, type = NULL, duration = NULL) {
      notifications[[length(notifications) + 1]] <<- list(
        message = message,
        type = type,
        duration = duration
      )
      invisible(NULL)
    }
  )

  expect_identical(result, expected_result)
  expect_identical(assigned_result, expected_result)
  expect_identical(
    notifications,
    list(list(
      message = "Design saved successfully. You can close this builder.",
      type = "message",
      duration = 5
    ))
  )
})

test_that("registerProtDesignSaveResultsObserver registers the save-results handoff shell", {
  call_log <- character()
  input <- list(
    save_results = 1,
    formula_string = "~ 0 + group"
  )
  design_matrix <- function() data.frame(Run = "S1", stringsAsFactors = FALSE)
  removed_samples <- function() "S2"
  data_cln <- function() data.frame(Run = c("S1", "S2"), stringsAsFactors = FALSE)
  contrast_data <- function() data.frame(
    contrast_name = "GA.vs.GC",
    numerator = "GA",
    denominator = "GC",
    stringsAsFactors = FALSE
  )
  config_list <- function() list(deAnalysisParameters = list(formula_string = "~ old"))
  result_rv <- function(value) value

  observer <- registerProtDesignSaveResultsObserver(
    input = input,
    designMatrix = design_matrix,
    removedSamples = removed_samples,
    dataCln = data_cln,
    contrastData = contrast_data,
    configList = config_list,
    resultRv = result_rv,
    observeEventFn = function(eventExpr, handlerExpr, ignoreNULL = FALSE) {
      event_sub <- substitute(eventExpr)
      handler_sub <- substitute(handlerExpr)
      event_value <- eval(event_sub, parent.frame())

      expect_identical(event_value, 1)
      expect_false(ignoreNULL)
      call_log <<- c(call_log, "observe")

      eval(handler_sub, parent.frame())
      "observer-registered"
    },
    runSaveResultsObserverShell = function(
        designMatrix,
        currentRemovedSamples = character(0),
        dataCln,
        contrastData,
        configList,
        formulaString,
        resultSetter
    ) {
      expect_identical(designMatrix$Run, "S1")
      expect_identical(currentRemovedSamples, "S2")
      expect_identical(dataCln$Run, c("S1", "S2"))
      expect_identical(contrastData$contrast_name, "GA.vs.GC")
      expect_identical(configList$deAnalysisParameters$formula_string, "~ old")
      expect_identical(formulaString, "~ 0 + group")
      expect_identical(resultSetter, result_rv)
      call_log <<- c(call_log, "shell")
      invisible("handled")
    }
  )

  expect_identical(observer, "observer-registered")
  expect_identical(call_log, c("observe", "shell"))
})

test_that("showProtDesignResetConfirmationModal builds the reset confirmation modal shell", {
  modal <- showProtDesignResetConfirmationModal(
    resetScope = "removed_samples",
    nsFn = function(id) paste0("builder-", id),
    showModalFn = function(modal) modal,
    modalDialogFn = function(title, body, footer, easyClose) {
      list(
        title = title,
        body = body,
        footer = footer,
        easyClose = easyClose
      )
    },
    htmlFn = function(text) text,
    tagListFn = function(...) list(...),
    modalButtonFn = function(label) list(type = "modalButton", label = label),
    actionButtonFn = function(inputId, label, class = NULL) {
      list(
        type = "actionButton",
        inputId = inputId,
        label = label,
        class = class
      )
    }
  )

  expect_identical(modal$title, "Confirm Reset")
  expect_match(modal$body, "removed_samples", fixed = TRUE)
  expect_true(modal$easyClose)
  expect_identical(
    modal$footer[[1]],
    list(type = "modalButton", label = "Cancel")
  )
  expect_identical(
    modal$footer[[2]],
    list(
      type = "actionButton",
      inputId = "builder-confirm_reset",
      label = "Reset",
      class = "btn-danger"
    )
  )
})

test_that("registerProtDesignResetRequestObserver registers the reset-request modal shell", {
  call_log <- character()
  input <- list(
    reset_changes = 1,
    reset_scope = "removed_samples"
  )

  observer <- registerProtDesignResetRequestObserver(
    input = input,
    session = list(ns = function(id) paste0("builder-", id)),
    observeEventFn = function(eventExpr, handlerExpr, ignoreNULL = FALSE) {
      event_sub <- substitute(eventExpr)
      handler_sub <- substitute(handlerExpr)
      event_value <- eval(event_sub, parent.frame())

      expect_identical(event_value, 1)
      expect_false(ignoreNULL)
      call_log <<- c(call_log, "observe")

      eval(handler_sub, parent.frame())
      "observer-registered"
    },
    showResetConfirmationModalFn = function(resetScope, nsFn) {
      expect_identical(resetScope, "removed_samples")
      expect_identical(nsFn("confirm_reset"), "builder-confirm_reset")
      call_log <<- c(call_log, "modal")
      invisible("handled")
    }
  )

  expect_identical(observer, "observer-registered")
  expect_identical(call_log, c("observe", "modal"))
})

test_that("runProtDesignResetConfirmationObserverShell applies the reset flow and closes the modal", {
  call_log <- list()
  design_matrix_value <- data.frame(Run = "ResetSample1", stringsAsFactors = FALSE)

  result <- runProtDesignResetConfirmationObserverShell(
    scope = "formula",
    initialState = list(formula = "~ condition"),
    designMatrix = function(value) {
      if (missing(value)) {
        call_log[[length(call_log) + 1]] <<- list(kind = "designMatrixGetter")
        return(design_matrix_value)
      }

      call_log[[length(call_log) + 1]] <<- list(kind = "designMatrixSetter", value = value)
      design_matrix_value <<- value
      invisible(NULL)
    },
    dataClnReactive = function(value) {
      call_log[[length(call_log) + 1]] <<- list(kind = "dataClnSetter", value = value)
      invisible(NULL)
    },
    removedSamples = function(value) {
      if (missing(value)) {
        call_log[[length(call_log) + 1]] <<- list(kind = "removedSamplesGetter")
        return(c("Drop1", "Drop2"))
      }

      call_log[[length(call_log) + 1]] <<- list(kind = "removedSamplesSetter", value = value)
      invisible(NULL)
    },
    factors = function(value) {
      call_log[[length(call_log) + 1]] <<- list(kind = "factorsSetter", value = value)
      invisible(NULL)
    },
    groups = function(value) {
      call_log[[length(call_log) + 1]] <<- list(kind = "groupsSetter", value = value)
      invisible(NULL)
    },
    contrasts = function(value) {
      call_log[[length(call_log) + 1]] <<- list(kind = "contrastsSetter", value = value)
      invisible(NULL)
    },
    session = "builder-session",
    applyResetStateFn = function(
      scope,
      initialState,
      designMatrixGetter,
      designMatrixSetter,
      dataClnSetter,
      currentRemovedSamples,
      removedSamplesSetter,
      factorsSetter,
      groupsSetter,
      contrastsSetter,
      updateFormulaFn
    ) {
      call_log[[length(call_log) + 1]] <<- list(
        kind = "applyResetState",
        scope = scope,
        initialState = initialState,
        currentRemovedSamples = currentRemovedSamples
      )
      expect_identical(scope, "formula")
      expect_identical(initialState$formula, "~ condition")
      expect_identical(currentRemovedSamples, c("Drop1", "Drop2"))
      expect_identical(designMatrixGetter(), design_matrix_value)
      designMatrixSetter(design_matrix_value)
      dataClnSetter(data.frame(Run = "ResetSample1", stringsAsFactors = FALSE))
      removedSamplesSetter(character(0))
      factorsSetter("condition")
      groupsSetter("GroupA")
      contrastsSetter(data.frame(contrast_name = "A.vs.B", stringsAsFactors = FALSE))
      updateFormulaFn("~ condition")
      invisible(NULL)
    },
    removeModalFn = function() {
      call_log[[length(call_log) + 1]] <<- list(kind = "removeModal")
      invisible(NULL)
    },
    showNotificationFn = function(message, type) {
      call_log[[length(call_log) + 1]] <<- list(kind = "showNotification", message = message, type = type)
      invisible(NULL)
    },
    updateTextInputFn = function(session, inputId, value) {
      call_log[[length(call_log) + 1]] <<- list(
        kind = "updateTextInput",
        session = session,
        inputId = inputId,
        value = value
      )
      invisible(NULL)
    }
  )

  expect_null(result)
  expect_identical(call_log[[1]]$kind, "removedSamplesGetter")
  expect_identical(call_log[[2]]$kind, "applyResetState")
  expect_identical(call_log[[3]]$kind, "designMatrixGetter")
  expect_identical(call_log[[4]]$kind, "designMatrixSetter")
  expect_identical(call_log[[5]]$kind, "dataClnSetter")
  expect_identical(call_log[[6]]$kind, "removedSamplesSetter")
  expect_identical(call_log[[7]]$kind, "factorsSetter")
  expect_identical(call_log[[8]]$kind, "groupsSetter")
  expect_identical(call_log[[9]]$kind, "contrastsSetter")
  expect_identical(call_log[[10]], list(
    kind = "updateTextInput",
    session = "builder-session",
    inputId = "formula_string",
    value = "~ condition"
  ))
  expect_identical(call_log[[11]], list(kind = "removeModal"))
  expect_identical(call_log[[12]], list(
    kind = "showNotification",
    message = "Reset of formula completed.",
    type = "message"
  ))
})

test_that("registerProtDesignResetConfirmationObserver registers the reset-confirmation handoff shell", {
  call_log <- character()
  input <- list(
    confirm_reset = 1,
    reset_scope = "all"
  )
  initial_state <- function() list(formula = "~ 0 + group")
  design_matrix <- function() data.frame(Run = "S1", stringsAsFactors = FALSE)
  data_cln <- function() data.frame(Run = "S1", stringsAsFactors = FALSE)
  removed_samples <- function() c("Drop1", "Drop2")
  factors_setter <- function(value) value
  groups_setter <- function(value) value
  contrasts_setter <- function(value) value

  observer <- registerProtDesignResetConfirmationObserver(
    input = input,
    initialState = initial_state,
    designMatrix = design_matrix,
    dataClnReactive = data_cln,
    removedSamples = removed_samples,
    factors = factors_setter,
    groups = groups_setter,
    contrasts = contrasts_setter,
    session = "builder-session",
    observeEventFn = function(eventExpr, handlerExpr, ignoreNULL = FALSE) {
      event_sub <- substitute(eventExpr)
      handler_sub <- substitute(handlerExpr)
      event_value <- eval(event_sub, parent.frame())

      expect_identical(event_value, 1)
      expect_false(ignoreNULL)
      call_log <<- c(call_log, "observe")

      eval(handler_sub, parent.frame())
      "observer-registered"
    },
    runResetConfirmationObserverShell = function(
      scope,
      initialState,
      designMatrix,
      dataClnReactive,
      removedSamples,
      factors,
      groups,
      contrasts,
      session
    ) {
      expect_identical(scope, "all")
      expect_identical(initialState$formula, "~ 0 + group")
      expect_identical(designMatrix, design_matrix)
      expect_identical(dataClnReactive, data_cln)
      expect_identical(removedSamples, removed_samples)
      expect_identical(factors, factors_setter)
      expect_identical(groups, groups_setter)
      expect_identical(contrasts, contrasts_setter)
      expect_identical(session, "builder-session")
      call_log <<- c(call_log, "shell")
      invisible("handled")
    }
  )

  expect_identical(observer, "observer-registered")
  expect_identical(call_log, c("observe", "shell"))
})

test_that("applyProtDesignResetState resets all scoped builder state through the provided setters", {
  design_matrix_value <- data.frame(
    Run = c("StateSample1", "StateSample2"),
    factor1 = c("treated", "control"),
    factor2 = c("batch1", "batch2"),
    factor3 = c("late", "early"),
    group = c("G1", "G2"),
    tech_reps = c(1L, 2L),
    batch = c("B1", "B2"),
    stringsAsFactors = FALSE
  )
  setter_log <- list()

  applyProtDesignResetState(
    scope = "all",
    initialState = list(
      design_matrix = design_matrix_value,
      data_cln = data.frame(Run = c("StateSample1", "StateSample2"), stringsAsFactors = FALSE),
      factors = c("condition", "batch"),
      groups = c("GroupA", "GroupB"),
      contrasts = data.frame(contrast_name = "A.vs.B", stringsAsFactors = FALSE),
      formula = "~ 0 + group"
    ),
    designMatrixGetter = function() {
      design_matrix_value
    },
    designMatrixSetter = function(value) {
      setter_log[[length(setter_log) + 1]] <<- list(kind = "designMatrix", value = value)
      design_matrix_value <<- value
      invisible(NULL)
    },
    dataClnSetter = function(value) {
      setter_log[[length(setter_log) + 1]] <<- list(kind = "dataCln", value = value)
      invisible(NULL)
    },
    currentRemovedSamples = c("Drop1", "Drop2"),
    removedSamplesSetter = function(value) {
      setter_log[[length(setter_log) + 1]] <<- list(kind = "removedSamples", value = value)
      invisible(NULL)
    },
    factorsSetter = function(value) {
      setter_log[[length(setter_log) + 1]] <<- list(kind = "factors", value = value)
      invisible(NULL)
    },
    groupsSetter = function(value) {
      setter_log[[length(setter_log) + 1]] <<- list(kind = "groups", value = value)
      invisible(NULL)
    },
    contrastsSetter = function(value) {
      setter_log[[length(setter_log) + 1]] <<- list(kind = "contrasts", value = value)
      invisible(NULL)
    },
    updateFormulaFn = function(value) {
      setter_log[[length(setter_log) + 1]] <<- list(kind = "formula", value = value)
      invisible(NULL)
    },
    applyRemovedSamplesUpdateFn = function(currentRemovedSamples, reset = FALSE) {
      expect_identical(currentRemovedSamples, c("Drop1", "Drop2"))
      expect_true(reset)
      "cleared"
    }
  )

  expect_identical(vapply(setter_log, `[[`, character(1), "kind"), c(
    "designMatrix",
    "dataCln",
    "removedSamples",
    "factors",
    "designMatrix",
    "groups",
    "contrasts",
    "formula"
  ))
  expect_identical(setter_log[[1]]$value, data.frame(
    Run = c("StateSample1", "StateSample2"),
    factor1 = c("treated", "control"),
    factor2 = c("batch1", "batch2"),
    factor3 = c("late", "early"),
    group = c("G1", "G2"),
    tech_reps = c(1L, 2L),
    batch = c("B1", "B2"),
    stringsAsFactors = FALSE
  ))
  expect_identical(setter_log[[3]]$value, "cleared")
  expect_identical(setter_log[[4]]$value, c("condition", "batch"))
  expect_identical(setter_log[[6]]$value, c("GroupA", "GroupB"))
  expect_identical(setter_log[[7]]$value$contrast_name, "A.vs.B")
  expect_identical(setter_log[[8]]$value, "~ 0 + group")
  expect_true(all(is.na(setter_log[[5]]$value$factor1)))
  expect_true(all(is.na(setter_log[[5]]$value$factor2)))
  expect_true(all(is.na(setter_log[[5]]$value$factor3)))
  expect_true(all(is.na(setter_log[[5]]$value$group)))
  expect_true(all(is.na(setter_log[[5]]$value$tech_reps)))
  expect_identical(setter_log[[5]]$value$batch, c("B1", "B2"))
})

test_that("applyProtDesignResetState only updates the formula for formula-only resets", {
  setter_log <- character()

  applyProtDesignResetState(
    scope = "formula",
    initialState = list(
      design_matrix = data.frame(Run = "StateSample1", stringsAsFactors = FALSE),
      data_cln = data.frame(Run = "StateSample1", stringsAsFactors = FALSE),
      factors = "condition",
      groups = "GroupA",
      contrasts = data.frame(contrast_name = "A.vs.B", stringsAsFactors = FALSE),
      formula = "~ condition"
    ),
    designMatrixGetter = function() {
      stop("designMatrixGetter should not be called for formula-only reset")
    },
    designMatrixSetter = function(value) {
      setter_log <<- c(setter_log, "designMatrix")
      invisible(NULL)
    },
    dataClnSetter = function(value) {
      setter_log <<- c(setter_log, "dataCln")
      invisible(NULL)
    },
    currentRemovedSamples = "Drop1",
    removedSamplesSetter = function(value) {
      setter_log <<- c(setter_log, "removedSamples")
      invisible(NULL)
    },
    factorsSetter = function(value) {
      setter_log <<- c(setter_log, "factors")
      invisible(NULL)
    },
    groupsSetter = function(value) {
      setter_log <<- c(setter_log, "groups")
      invisible(NULL)
    },
    contrastsSetter = function(value) {
      setter_log <<- c(setter_log, "contrasts")
      invisible(NULL)
    },
    updateFormulaFn = function(value) {
      expect_identical(value, "~ condition")
      setter_log <<- c(setter_log, "formula")
      invisible(NULL)
    },
    applyRemovedSamplesUpdateFn = function(...) {
      stop("applyRemovedSamplesUpdateFn should not be called for formula-only reset")
    }
  )

  expect_identical(setter_log, "formula")
})

test_that("applyProtDesignRemovedSamplesUpdate accumulates unique removed samples in append order", {
  updated <- applyProtDesignRemovedSamplesUpdate(
    currentRemovedSamples = c("Sample1", "Sample3"),
    selectedSamples = c("Sample2", "Sample3", "Sample4")
  )

  expect_identical(updated, c("Sample1", "Sample3", "Sample2", "Sample4"))
})

test_that("applyProtDesignRemovedSamplesUpdate clears removed samples on reset", {
  updated <- applyProtDesignRemovedSamplesUpdate(
    currentRemovedSamples = c("Sample1", "Sample2"),
    reset = TRUE
  )

  expect_identical(updated, character(0))
})

test_that("formatProtDesignRemovedSamplesDisplay reports when no samples are removed", {
  expect_identical(
    formatProtDesignRemovedSamplesDisplay(character(0)),
    "No samples have been removed."
  )
})

test_that("formatProtDesignRemovedSamplesDisplay counts and mixed-sorts removed samples", {
  expect_identical(
    formatProtDesignRemovedSamplesDisplay(c("Sample10", "Sample2", "Sample1")),
    paste(
      "Removed 3 sample(s):",
      paste(c("Sample1", "Sample2", "Sample10"), collapse = "\n"),
      sep = "\n"
    )
  )
})

test_that("formatProtDesignContrastFactorsInfo reports prefixed group contrasts for grouped formulas", {
  expect_identical(
    formatProtDesignContrastFactorsInfo("~ 0 + group"),
    paste(
      "Note: Contrasts will use 'group' prefix (e.g., groupGA_Control-groupGA_Elevated)",
      "based on current formula: ~ 0 + group",
      sep = "\n"
    )
  )
})

test_that("formatProtDesignContrastFactorsInfo reports as-is contrasts for non-grouped formulas", {
  expect_identical(
    formatProtDesignContrastFactorsInfo("~ condition + batch"),
    "Note: Contrasts will use group names as-is (e.g., GA_Control-GA_Elevated)"
  )
})

test_that("buildProtDesignReplicateInputs wires the replicate_start numeric input for selected runs", {
  numeric_input <- buildProtDesignReplicateInputs(
    selectedRuns = c("Sample1", "Sample2", "Sample3"),
    nsFn = function(id) paste0("builder-", id),
    numericInputFn = function(inputId, label, value, min) {
      list(
        inputId = inputId,
        label = label,
        value = value,
        min = min
      )
    }
  )

  expect_identical(
    numeric_input,
    list(
      inputId = "builder-replicate_start",
      label = "Starting replicate number for 3 selected runs:",
      value = 1,
      min = 1
    )
  )
})

test_that("formatProtDesignRangePreview renders the first selected sample preview", {
  preview <- formatProtDesignRangePreview(
    selectedSamples = c("SampleA_T1_Rep1", "SampleB_T2_Rep1"),
    rangeStart = 2,
    rangeEnd = 2,
    extractExperimentFn = function(sampleName, mode, start, end) {
      expect_identical(sampleName, "SampleA_T1_Rep1")
      expect_identical(mode, "range")
      expect_identical(start, 2)
      expect_identical(end, 2)
      "T1"
    }
  )

  expect_identical(preview, "\"SampleA_T1_Rep1\" -> \"T1\"")
})

test_that("formatProtDesignRangePreview reports extraction errors", {
  preview <- formatProtDesignRangePreview(
    selectedSamples = "SampleA_T1_Rep1",
    rangeStart = 4,
    rangeEnd = 2,
    extractExperimentFn = function(...) {
      stop("bad range")
    }
  )

  expect_identical(preview, "Error: bad range")
})

test_that("transformProtDesignSampleNames routes each supported transform mode", {
  extracted_calls <- list()
  req_calls <- list()
  format_optional <- function(value) if (is.null(value)) "" else as.character(value)

  extract_fn <- function(sampleName, mode, start = NULL, end = NULL) {
    extracted_calls[[length(extracted_calls) + 1]] <<- list(
      sampleName = sampleName,
      mode = mode,
      start = start,
      end = end
    )
    paste(sampleName, mode, format_optional(start), format_optional(end))
  }

  req_fn <- function(...) {
    req_calls[[length(req_calls) + 1]] <<- list(...)
    invisible(NULL)
  }

  expect_identical(
    transformProtDesignSampleNames(
      selectedSamples = c("SampleA_T1_Rep1", "SampleB_T2_Rep1"),
      transformMode = "range",
      rangeStart = 2,
      rangeEnd = 2,
      extractExperimentFn = extract_fn,
      reqFn = req_fn
    ),
    c("SampleA_T1_Rep1 range 2 2", "SampleB_T2_Rep1 range 2 2")
  )
  expect_length(req_calls, 2)
  expect_identical(
    transformProtDesignSampleNames(
      selectedSamples = "SampleA_T1_Rep1",
      transformMode = "before_underscore",
      extractExperimentFn = extract_fn
    ),
    "SampleA_T1_Rep1 start  "
  )
  expect_identical(
    transformProtDesignSampleNames(
      selectedSamples = "SampleA_T1_Rep1",
      transformMode = "after_underscore",
      extractExperimentFn = extract_fn
    ),
    "SampleA_T1_Rep1 end  "
  )

  expect_identical(extracted_calls[[1]], list(
    sampleName = "SampleA_T1_Rep1",
    mode = "range",
    start = 2,
    end = 2
  ))
  expect_identical(extracted_calls[[4]], list(
    sampleName = "SampleA_T1_Rep1",
    mode = "end",
    start = NULL,
    end = NULL
  ))
})

test_that("transformProtDesignSampleNames rejects unsupported transform modes", {
  expect_error(
    transformProtDesignSampleNames(
      selectedSamples = "SampleA_T1_Rep1",
      transformMode = "invalid",
      extractExperimentFn = function(...) "ignored"
    ),
    "Unsupported transform mode: invalid"
  )
})

test_that("applyProtDesignBulkRenameUpdates rewrites both design and data tables", {
  design_matrix <- data.frame(
    Run = c("SampleA_T1_Rep1", "SampleB_T2_Rep1", "SampleC_T3_Rep1"),
    group = c("G1", "G2", "G3"),
    stringsAsFactors = FALSE
  )
  data_cln <- data.frame(
    Run = c("SampleA_T1_Rep1", "SampleB_T2_Rep1", "SampleA_T1_Rep1"),
    intensity = c(10, 20, 30),
    stringsAsFactors = FALSE
  )

  updated <- applyProtDesignBulkRenameUpdates(
    designMatrix = design_matrix,
    dataCln = data_cln,
    selectedSamples = c("SampleA_T1_Rep1", "SampleB_T2_Rep1"),
    newNames = c("SampleA", "SampleB")
  )

  expect_identical(
    updated$designMatrix$Run,
    c("SampleA", "SampleB", "SampleC_T3_Rep1")
  )
  expect_identical(
    updated$dataCln$Run,
    c("SampleA", "SampleB", "SampleA")
  )
  expect_identical(updated$dataCln$intensity, c(10, 20, 30))
})

test_that("applyProtDesignSingleRenameUpdate rewrites both design and data tables", {
  design_matrix <- data.frame(
    Run = c("SampleA_T1_Rep1", "SampleB_T2_Rep1", "SampleC_T3_Rep1"),
    group = c("G1", "G2", "G3"),
    stringsAsFactors = FALSE
  )
  data_cln <- data.frame(
    Run = c("SampleA_T1_Rep1", "SampleB_T2_Rep1", "SampleA_T1_Rep1"),
    intensity = c(10, 20, 30),
    stringsAsFactors = FALSE
  )

  updated <- applyProtDesignSingleRenameUpdate(
    designMatrix = design_matrix,
    dataCln = data_cln,
    originalName = "SampleA_T1_Rep1",
    newName = "SampleA"
  )

  expect_identical(
    updated$designMatrix$Run,
    c("SampleA", "SampleB_T2_Rep1", "SampleC_T3_Rep1")
  )
  expect_identical(
    updated$dataCln$Run,
    c("SampleA", "SampleB_T2_Rep1", "SampleA")
  )
  expect_identical(updated$dataCln$intensity, c(10, 20, 30))
})

test_that("registerProtDesignRenameSampleObserver registers the single-rename handoff shell", {
  call_log <- list()
  current_design_matrix <- data.frame(
    Run = c("SampleA_T1_Rep1", "SampleB_T2_Rep1"),
    group = c("G1", "G2"),
    stringsAsFactors = FALSE
  )
  current_data_cln <- data.frame(
    Run = c("SampleA_T1_Rep1", "SampleB_T2_Rep1"),
    intensity = c(10, 20),
    stringsAsFactors = FALSE
  )
  renamed_design_matrix <- transform(current_design_matrix, Run = c("SampleA", "SampleB_T2_Rep1"))
  renamed_data_cln <- transform(current_data_cln, Run = c("SampleA", "SampleB_T2_Rep1"))
  input <- list(
    rename_sample = 1,
    sample_to_rename = "SampleA_T1_Rep1",
    new_sample_name = "SampleA"
  )

  design_matrix_rv <- function(value) {
    if (missing(value)) {
      call_log[[length(call_log) + 1]] <<- list(kind = "getDesignMatrix")
      current_design_matrix
    } else {
      call_log[[length(call_log) + 1]] <<- list(kind = "setDesignMatrix", value = value)
      current_design_matrix <<- value
      invisible(NULL)
    }
  }
  data_cln_rv <- function(value) {
    if (missing(value)) {
      call_log[[length(call_log) + 1]] <<- list(kind = "getDataCln")
      current_data_cln
    } else {
      call_log[[length(call_log) + 1]] <<- list(kind = "setDataCln", value = value)
      current_data_cln <<- value
      invisible(NULL)
    }
  }

  observer <- registerProtDesignRenameSampleObserver(
    input = input,
    designMatrix = design_matrix_rv,
    dataCln = data_cln_rv,
    session = "builder-session",
    observeEventFn = function(eventExpr, handlerExpr, ignoreNULL = FALSE) {
      event_sub <- substitute(eventExpr)
      handler_sub <- substitute(handlerExpr)
      event_value <- eval(event_sub, parent.frame())

      expect_identical(event_value, 1)
      expect_false(ignoreNULL)
      call_log[[length(call_log) + 1]] <<- list(kind = "observe")

      eval(handler_sub, parent.frame())
      "observer-registered"
    },
    reqFn = function(...) {
      expect_identical(list(...), list("SampleA_T1_Rep1", "SampleA"))
      call_log[[length(call_log) + 1]] <<- list(kind = "req", value = list(...))
      invisible(NULL)
    },
    applySingleRenameUpdateFn = function(designMatrix, dataCln, originalName, newName) {
      expect_identical(designMatrix, current_design_matrix)
      expect_identical(dataCln, current_data_cln)
      expect_identical(originalName, "SampleA_T1_Rep1")
      expect_identical(newName, "SampleA")
      call_log[[length(call_log) + 1]] <<- list(
        kind = "apply",
        designMatrix = designMatrix,
        dataCln = dataCln,
        originalName = originalName,
        newName = newName
      )
      list(
        designMatrix = renamed_design_matrix,
        dataCln = renamed_data_cln
      )
    },
    updateTextInputFn = function(session, inputId, value = NULL) {
      call_log[[length(call_log) + 1]] <<- list(
        kind = "updateTextInput",
        session = session,
        inputId = inputId,
        value = value
      )
      invisible(NULL)
    }
  )

  expect_identical(observer, "observer-registered")
  expect_identical(current_design_matrix, renamed_design_matrix)
  expect_identical(current_data_cln, renamed_data_cln)
  expect_identical(call_log, list(
    list(kind = "observe"),
    list(kind = "req", value = list("SampleA_T1_Rep1", "SampleA")),
    list(kind = "getDesignMatrix"),
    list(kind = "getDataCln"),
    list(
      kind = "apply",
      designMatrix = data.frame(
        Run = c("SampleA_T1_Rep1", "SampleB_T2_Rep1"),
        group = c("G1", "G2"),
        stringsAsFactors = FALSE
      ),
      dataCln = data.frame(
        Run = c("SampleA_T1_Rep1", "SampleB_T2_Rep1"),
        intensity = c(10, 20),
        stringsAsFactors = FALSE
      ),
      originalName = "SampleA_T1_Rep1",
      newName = "SampleA"
    ),
    list(kind = "setDesignMatrix", value = renamed_design_matrix),
    list(kind = "setDataCln", value = renamed_data_cln),
    list(
      kind = "updateTextInput",
      session = "builder-session",
      inputId = "new_sample_name",
      value = ""
    )
  ))
})

test_that("registerProtDesignBulkRenameObserver registers the bulk-rename handoff shell", {
  call_log <- list()
  current_design_matrix <- data.frame(
    Run = c("SampleA_T1_Rep1", "SampleB_T2_Rep1", "SampleC_T3_Rep1"),
    group = c("G1", "G2", "G3"),
    stringsAsFactors = FALSE
  )
  current_data_cln <- data.frame(
    Run = c("SampleA_T1_Rep1", "SampleB_T2_Rep1", "SampleC_T3_Rep1"),
    intensity = c(10, 20, 30),
    stringsAsFactors = FALSE
  )
  renamed_design_matrix <- transform(
    current_design_matrix,
    Run = c("SampleA", "SampleB", "SampleC_T3_Rep1")
  )
  renamed_data_cln <- transform(
    current_data_cln,
    Run = c("SampleA", "SampleB", "SampleC_T3_Rep1")
  )
  input <- list(
    bulk_rename = 1,
    samples_to_transform = c("SampleA_T1_Rep1", "SampleB_T2_Rep1"),
    transform_mode = "before_underscore",
    range_start = 1,
    range_end = 1
  )

  design_matrix_rv <- function(value) {
    if (missing(value)) {
      call_log[[length(call_log) + 1]] <<- list(kind = "getDesignMatrix")
      current_design_matrix
    } else {
      call_log[[length(call_log) + 1]] <<- list(kind = "setDesignMatrix", value = value)
      current_design_matrix <<- value
      invisible(NULL)
    }
  }
  data_cln_rv <- function(value) {
    if (missing(value)) {
      call_log[[length(call_log) + 1]] <<- list(kind = "getDataCln")
      current_data_cln
    } else {
      call_log[[length(call_log) + 1]] <<- list(kind = "setDataCln", value = value)
      current_data_cln <<- value
      invisible(NULL)
    }
  }

  observer <- registerProtDesignBulkRenameObserver(
    input = input,
    designMatrix = design_matrix_rv,
    dataCln = data_cln_rv,
    observeEventFn = function(eventExpr, handlerExpr, ignoreNULL = FALSE) {
      event_sub <- substitute(eventExpr)
      handler_sub <- substitute(handlerExpr)
      event_value <- eval(event_sub, parent.frame())

      expect_identical(event_value, 1)
      expect_false(ignoreNULL)
      call_log[[length(call_log) + 1]] <<- list(kind = "observe")

      eval(handler_sub, parent.frame())
      "observer-registered"
    },
    reqFn = function(value) {
      expect_identical(value, c("SampleA_T1_Rep1", "SampleB_T2_Rep1"))
      call_log[[length(call_log) + 1]] <<- list(kind = "req", value = value)
      invisible(NULL)
    },
    transformSampleNamesFn = function(selectedSamples, transformMode, rangeStart, rangeEnd) {
      expect_identical(selectedSamples, c("SampleA_T1_Rep1", "SampleB_T2_Rep1"))
      expect_identical(transformMode, "before_underscore")
      expect_identical(rangeStart, 1)
      expect_identical(rangeEnd, 1)
      call_log[[length(call_log) + 1]] <<- list(
        kind = "transform",
        selectedSamples = selectedSamples,
        transformMode = transformMode,
        rangeStart = rangeStart,
        rangeEnd = rangeEnd
      )
      c("SampleA", "SampleB")
    },
    applyBulkRenameUpdatesFn = function(designMatrix, dataCln, selectedSamples, newNames) {
      expect_identical(designMatrix, current_design_matrix)
      expect_identical(dataCln, current_data_cln)
      expect_identical(selectedSamples, c("SampleA_T1_Rep1", "SampleB_T2_Rep1"))
      expect_identical(newNames, c("SampleA", "SampleB"))
      call_log[[length(call_log) + 1]] <<- list(
        kind = "apply",
        designMatrix = designMatrix,
        dataCln = dataCln,
        selectedSamples = selectedSamples,
        newNames = newNames
      )
      list(
        designMatrix = renamed_design_matrix,
        dataCln = renamed_data_cln
      )
    }
  )

  expect_identical(observer, "observer-registered")
  expect_identical(current_design_matrix, renamed_design_matrix)
  expect_identical(current_data_cln, renamed_data_cln)
  expect_identical(call_log, list(
    list(kind = "observe"),
    list(kind = "req", value = c("SampleA_T1_Rep1", "SampleB_T2_Rep1")),
    list(kind = "getDesignMatrix"),
    list(kind = "getDataCln"),
    list(
      kind = "transform",
      selectedSamples = c("SampleA_T1_Rep1", "SampleB_T2_Rep1"),
      transformMode = "before_underscore",
      rangeStart = 1,
      rangeEnd = 1
    ),
    list(
      kind = "apply",
      designMatrix = data.frame(
        Run = c("SampleA_T1_Rep1", "SampleB_T2_Rep1", "SampleC_T3_Rep1"),
        group = c("G1", "G2", "G3"),
        stringsAsFactors = FALSE
      ),
      dataCln = data.frame(
        Run = c("SampleA_T1_Rep1", "SampleB_T2_Rep1", "SampleC_T3_Rep1"),
        intensity = c(10, 20, 30),
        stringsAsFactors = FALSE
      ),
      selectedSamples = c("SampleA_T1_Rep1", "SampleB_T2_Rep1"),
      newNames = c("SampleA", "SampleB")
    ),
    list(kind = "setDesignMatrix", value = renamed_design_matrix),
    list(kind = "setDataCln", value = renamed_data_cln)
  ))
})

test_that("registerProtDesignPreviewOutputs wires availability flags and preview renders", {
  workflow_data <- new.env(parent = emptyenv())
  workflow_data$data_tbl <- data.frame(Run = "S1", stringsAsFactors = FALSE)
  workflow_data$config_list <- list(globalParameters = list())
  workflow_data$design_matrix <- data.frame(Run = "S1", group = "G1", stringsAsFactors = FALSE)
  workflow_data$contrasts_tbl <- data.frame(contrast = "groupA-groupB", stringsAsFactors = FALSE)

  output <- new.env(parent = emptyenv())
  output_option_calls <- list()

  reactiveFn <- function(expr) {
    expr_sub <- substitute(expr)
    expr_env <- parent.frame()
    function() eval(expr_sub, expr_env)
  }

  renderDT <- function(expr, options = NULL) {
    list(
      expr = substitute(expr),
      env = parent.frame(),
      options = options
    )
  }

  result <- registerProtDesignPreviewOutputs(
    output = output,
    workflowData = workflow_data,
    reactiveFn = reactiveFn,
    outputOptionsFn = function(output, name, suspendWhenHidden = FALSE) {
      output_option_calls[[length(output_option_calls) + 1]] <<- list(
        output = output,
        name = name,
        suspendWhenHidden = suspendWhenHidden
      )
      invisible(NULL)
    },
    renderDT = renderDT,
    reqFn = function(value) {
      if (is.null(value)) {
        stop("required value was NULL")
      }
      value
    }
  )

  expect_identical(result, output)
  expect_true(output$data_available())
  expect_true(output$design_matrix_exists())
  expect_length(output_option_calls, 2)
  expect_identical(output_option_calls[[1]]$name, "data_available")
  expect_false(output_option_calls[[1]]$suspendWhenHidden)
  expect_identical(output_option_calls[[2]]$name, "design_matrix_exists")
  expect_false(output_option_calls[[2]]$suspendWhenHidden)

  expect_identical(
    eval(output$design_matrix_preview$expr, output$design_matrix_preview$env),
    workflow_data$design_matrix
  )
  expect_identical(output$design_matrix_preview$options, list(pageLength = 5, scrollX = TRUE))
  expect_identical(
    eval(output$contrasts_preview$expr, output$contrasts_preview$env),
    workflow_data$contrasts_tbl
  )
  expect_identical(output$contrasts_preview$options, list(pageLength = 5, scrollX = TRUE))
})

test_that("registerProtDesignBuilderModule wires builder inputs and fallback reactiveVal", {
  workflow_data <- new.env(parent = emptyenv())
  workflow_data$data_tbl <- data.frame(Run = "S1", stringsAsFactors = FALSE)
  workflow_data$config_list <- list(globalParameters = list(workflow = "design"))
  workflow_data$column_mapping <- list(sample_id = "Run")

  reactiveFn <- function(expr) {
    expr_sub <- substitute(expr)
    expr_env <- parent.frame()
    function() eval(expr_sub, expr_env)
  }

  fallback <- registerProtDesignBuilderModule(
    workflowData = workflow_data,
    builderServerExists = FALSE,
    reactiveFn = reactiveFn,
    reactiveValFn = function(value) {
      function() value
    }
  )

  expect_null(fallback())

  builder_call <- NULL
  builder_results <- function() "builder-results"
  registered <- registerProtDesignBuilderModule(
    workflowData = workflow_data,
    moduleId = "builder-test",
    builderServerExists = TRUE,
    builderServerFn = function(id, data_tbl, config_list, column_mapping) {
      builder_call <<- list(
        id = id,
        data_tbl = data_tbl,
        config_list = config_list,
        column_mapping = column_mapping
      )
      builder_results
    },
    reactiveFn = reactiveFn
  )

  expect_identical(registered, builder_results)
  expect_identical(builder_call$id, "builder-test")
  expect_identical(builder_call$data_tbl(), workflow_data$data_tbl)
  expect_identical(builder_call$config_list(), workflow_data$config_list)
  expect_identical(builder_call$column_mapping(), workflow_data$column_mapping)
})

test_that("registerProtDesignBuilderModule uses the default builder id and server binding", {
  workflow_data <- new.env(parent = emptyenv())
  workflow_data$data_tbl <- data.frame(Run = "S1", stringsAsFactors = FALSE)
  workflow_data$config_list <- list(globalParameters = list(workflow = "design"))
  workflow_data$column_mapping <- list(sample_id = "Run")

  reactiveFn <- function(expr) {
    expr_sub <- substitute(expr)
    expr_env <- parent.frame()
    function() eval(expr_sub, expr_env)
  }

  builder_call <- NULL
  builder_results <- function() "builder-results"
  builder_helper <- registerProtDesignBuilderModule
  environment(builder_helper) <- list2env(
    list(
      exists = function(x, where = -1, inherits = TRUE) {
        if (identical(x, "mod_prot_design_builder_server")) {
          TRUE
        } else {
          base::exists(x, where = where, inherits = inherits)
        }
      },
      mod_prot_design_builder_server = function(id, data_tbl, config_list, column_mapping) {
        builder_call <<- list(
          id = id,
          data_tbl = data_tbl,
          config_list = config_list,
          column_mapping = column_mapping
        )
        builder_results
      }
    ),
    parent = environment(registerProtDesignBuilderModule)
  )

  registered <- builder_helper(
    workflowData = workflow_data,
    reactiveFn = reactiveFn
  )

  expect_identical(registered, builder_results)
  expect_identical(builder_call$id, "builder")
  expect_identical(builder_call$data_tbl(), workflow_data$data_tbl)
  expect_identical(builder_call$config_list(), workflow_data$config_list)
  expect_identical(builder_call$column_mapping(), workflow_data$column_mapping)
})

test_that("initializeProtDesignImportBootstrap resolves volumes and wires shinyFiles setup", {
  base_dir <- tempfile("prot-design-base-")
  dir.create(base_dir)

  input <- list(import_dir = NULL, import_fasta_file = NULL)
  session_obj <- structure(list(user = "test-session"), class = "test_session")
  experiment_paths <- list(base_dir = base_dir)
  dir_choose_call <- NULL
  file_choose_call <- NULL
  reactive_seed <- new.env(parent = emptyenv())
  log_messages <- character()

  bootstrap <- initializeProtDesignImportBootstrap(
    input = input,
    session = session_obj,
    experimentPaths = experiment_paths,
    volumes = function() c(Home = "/tmp/home-volume"),
    dirChooseFn = function(input, id, roots, session) {
      dir_choose_call <<- list(input = input, id = id, roots = roots, session = session)
      invisible(NULL)
    },
    fileChooseFn = function(input, id, roots, session, filetypes) {
      file_choose_call <<- list(
        input = input,
        id = id,
        roots = roots,
        session = session,
        filetypes = filetypes
      )
      invisible(NULL)
    },
    reactiveValFn = function(value) {
      reactive_seed$value <- value
      function(new_value) {
        if (!missing(new_value)) {
          reactive_seed$value <<- new_value
        }
        reactive_seed$value
      }
    },
    isolateFn = function(expr) eval(substitute(expr), parent.frame()),
    dirExistsFn = function(path) identical(path, base_dir),
    logInfo = function(message) {
      log_messages <<- c(log_messages, message)
      invisible(NULL)
    }
  )

  expect_identical(
    bootstrap$resolvedVolumes,
    c("Project Base Dir" = base_dir, Home = "/tmp/home-volume")
  )
  expect_identical(dir_choose_call$id, "import_dir")
  expect_identical(dir_choose_call$roots, bootstrap$resolvedVolumes)
  expect_identical(dir_choose_call$session, session_obj)
  expect_identical(file_choose_call$id, "import_fasta_file")
  expect_identical(file_choose_call$roots, bootstrap$resolvedVolumes)
  expect_identical(file_choose_call$session, session_obj)
  expect_identical(file_choose_call$filetypes, c("fasta", "fa", "faa"))
  expect_null(bootstrap$importFastaPath())
  bootstrap$importFastaPath("/tmp/prot-design-selected.fasta")
  expect_identical(bootstrap$importFastaPath(), "/tmp/prot-design-selected.fasta")
  expect_length(log_messages, 1)
  expect_match(log_messages[[1]], "Added base_dir to volumes")
})

test_that("registerProtDesignImportConfirmationObserver registers the import handoff shell", {
  call_log <- character()
  workflow_data <- new.env(parent = emptyenv())
  experiment_paths <- list(source_dir = tempfile("prot-design-source-"))
  session_obj <- structure(list(user = "test-session"), class = "test_session")
  input <- list(
    confirm_import = 1,
    import_dir = "import-token",
    import_taxon_id = 9606,
    import_organism_name = "Homo sapiens"
  )
  resolved_volumes <- c("Project Base Dir" = tempdir())
  import_fasta_path <- function() "/tmp/prot-design-selected.fasta"

  observer <- registerProtDesignImportConfirmationObserver(
    input = input,
    resolvedVolumes = resolved_volumes,
    importFastaPath = import_fasta_path,
    workflowData = workflow_data,
    experimentPaths = experiment_paths,
    session = session_obj,
    qcTrigger = TRUE,
    observeEventFn = function(eventExpr, handlerExpr) {
      event_sub <- substitute(eventExpr)
      handler_sub <- substitute(handlerExpr)
      event_value <- eval(event_sub, parent.frame())

      expect_identical(event_value, 1)
      call_log <<- c(call_log, "observe")

      eval(handler_sub, parent.frame())
      "observer-registered"
    },
    reqFn = function(value) {
      call_log <<- c(call_log, paste0("req:", if (identical(value, "import-token")) "dir" else "path"))
      value
    },
    parseDirPathFn = function(roots, selection) {
      expect_identical(roots, resolved_volumes)
      expect_identical(selection, "import-token")
      call_log <<- c(call_log, "parse")
      "/tmp/prot-design-import"
    },
    removeModalFn = function() {
      call_log <<- c(call_log, "remove")
      invisible(NULL)
    },
    runImportObserverShell = function(workflowData, experimentPaths, importPath, selectedFastaPath = NULL, taxonId, organismName, session, qcTrigger = NULL) {
      expect_identical(workflowData, workflow_data)
      expect_identical(experimentPaths, experiment_paths)
      expect_identical(importPath, "/tmp/prot-design-import")
      expect_identical(selectedFastaPath, "/tmp/prot-design-selected.fasta")
      expect_identical(taxonId, 9606)
      expect_identical(organismName, "Homo sapiens")
      expect_identical(session, session_obj)
      expect_true(qcTrigger)
      call_log <<- c(call_log, "shell")
      invisible("handled")
    }
  )

  expect_identical(observer, "observer-registered")
  expect_identical(call_log, c("observe", "req:dir", "parse", "req:path", "remove", "shell"))
})

test_that("registerProtDesignBuilderResultsObserver registers the builder handoff shell", {
  call_log <- character()
  workflow_data <- new.env(parent = emptyenv())
  experiment_paths <- list(source_dir = tempfile("prot-design-source-"))
  session_obj <- structure(list(user = "test-session"), class = "test_session")
  builder_results <- list(
    design_matrix = data.frame(Run = "S1", group = "G1", stringsAsFactors = FALSE)
  )
  builder_results_rv <- function() builder_results

  observer <- registerProtDesignBuilderResultsObserver(
    builderResultsRv = builder_results_rv,
    workflowData = workflow_data,
    experimentPaths = experiment_paths,
    session = session_obj,
    qcTrigger = TRUE,
    observeEventFn = function(eventExpr, handlerExpr, ignoreNULL = FALSE) {
      event_sub <- substitute(eventExpr)
      handler_sub <- substitute(handlerExpr)
      event_value <- eval(event_sub, parent.frame())

      expect_identical(event_value, builder_results)
      expect_true(ignoreNULL)
      call_log <<- c(call_log, "observe")

      eval(handler_sub, parent.frame())
      "observer-registered"
    },
    reqFn = function(value) {
      expect_identical(value, builder_results)
      call_log <<- c(call_log, "req")
      value
    },
    runBuilderObserverShell = function(results, workflowData, experimentPaths, session, qcTrigger = NULL) {
      expect_identical(results, builder_results)
      expect_identical(workflowData, workflow_data)
      expect_identical(experimentPaths, experiment_paths)
      expect_identical(session, session_obj)
      expect_true(qcTrigger)
      call_log <<- c(call_log, "shell")
      invisible("handled")
    }
  )

  expect_identical(observer, "observer-registered")
  expect_identical(call_log, c("observe", "req", "shell"))
})

test_that("registerProtDesignServerShells wires the wrapper shell seam in order", {
  call_log <- character()
  workflow_data <- shiny::reactiveValues(
    data_tbl = data.frame(Run = "S1", stringsAsFactors = FALSE),
    config_list = list(globalParameters = list(workflow = "design")),
    column_mapping = list(sample_id = "Run")
  )
  experiment_paths <- list(source_dir = tempfile("prot-design-source-"))
  fake_input <- list(confirm_import = NULL)
  fake_output <- structure(list(), class = "shinyoutput")
  fake_session <- list(ns = function(id) paste0("design-entry-", id))
  resolved_volumes <- c(Home = tempdir())
  import_fasta_seed <- function() "/tmp/prot-design-selected.fasta"
  builder_results_rv <- function() "builder-results"

  result <- registerProtDesignServerShells(
    input = fake_input,
    output = fake_output,
    session = fake_session,
    workflowData = workflow_data,
    experimentPaths = experiment_paths,
    volumes = "volumes-seed",
    qcTrigger = "qc-seed",
    initializeImportBootstrap = function(input, session, experimentPaths, volumes = NULL) {
      expect_identical(input, fake_input)
      expect_identical(session, fake_session)
      expect_identical(experimentPaths, experiment_paths)
      expect_identical(volumes, "volumes-seed")
      call_log <<- c(call_log, "bootstrap")

      list(
        resolvedVolumes = resolved_volumes,
        importFastaPath = import_fasta_seed
      )
    },
    registerImportModalShell = function(input, output, session, resolvedVolumes, importFastaPath) {
      expect_identical(input, fake_input)
      expect_identical(output, fake_output)
      expect_identical(session, fake_session)
      expect_identical(resolvedVolumes, resolved_volumes)
      expect_identical(importFastaPath, import_fasta_seed)
      call_log <<- c(call_log, "modal")
      invisible(NULL)
    },
    registerImportConfirmationObserver = function(input, resolvedVolumes, importFastaPath, workflowData, experimentPaths, session, qcTrigger = NULL) {
      expect_identical(input, fake_input)
      expect_identical(resolvedVolumes, resolved_volumes)
      expect_identical(importFastaPath, import_fasta_seed)
      expect_identical(workflowData, workflow_data)
      expect_identical(experimentPaths, experiment_paths)
      expect_identical(session, fake_session)
      expect_identical(qcTrigger, "qc-seed")
      call_log <<- c(call_log, "importObserver")
      invisible(NULL)
    },
    registerPreviewOutputs = function(output, workflowData) {
      expect_identical(output, fake_output)
      expect_identical(workflowData, workflow_data)
      call_log <<- c(call_log, "preview")
      invisible(output)
    },
    registerBuilderModule = function(workflowData) {
      expect_identical(workflowData, workflow_data)
      call_log <<- c(call_log, "builderModule")
      builder_results_rv
    },
    registerBuilderResultsObserver = function(builderResultsRv, workflowData, experimentPaths, session, qcTrigger = NULL) {
      expect_identical(builderResultsRv, builder_results_rv)
      expect_identical(workflowData, workflow_data)
      expect_identical(experimentPaths, experiment_paths)
      expect_identical(session, fake_session)
      expect_identical(qcTrigger, "qc-seed")
      call_log <<- c(call_log, "builderObserver")
      invisible(NULL)
    }
  )

  expect_identical(result, builder_results_rv)
  expect_identical(
    call_log,
    c("bootstrap", "modal", "importObserver", "preview", "builderModule", "builderObserver")
  )
})

test_that("registerProtDesignServerShells passes callback seeds opaquely through the wrapper shell", {
  call_log <- character()
  workflow_data <- shiny::reactiveValues(
    data_tbl = data.frame(Run = "S1", stringsAsFactors = FALSE),
    config_list = list(globalParameters = list(workflow = "design")),
    column_mapping = list(sample_id = "Run")
  )
  experiment_paths <- list(source_dir = tempfile("prot-design-source-"))
  fake_input <- list(confirm_import = NULL)
  fake_output <- structure(list(), class = "shinyoutput")
  fake_session <- list(ns = function(id) paste0("design-entry-", id))
  resolved_volumes_seed <- c(
    "Project Base Dir" = tempfile("prot-design-base-"),
    Home = tempfile("prot-design-home-")
  )
  import_fasta_calls <- 0L
  builder_results_calls <- 0L
  qc_trigger_calls <- 0L
  import_fasta_seed <- function(...) {
    import_fasta_calls <<- import_fasta_calls + 1L
    "/tmp/prot-design-selected.fasta"
  }
  builder_results_rv <- function(...) {
    builder_results_calls <<- builder_results_calls + 1L
    "builder-results"
  }
  qc_trigger_seed <- function(...) {
    qc_trigger_calls <<- qc_trigger_calls + 1L
    invisible(NULL)
  }

  result <- registerProtDesignServerShells(
    input = fake_input,
    output = fake_output,
    session = fake_session,
    workflowData = workflow_data,
    experimentPaths = experiment_paths,
    volumes = "volumes-seed",
    qcTrigger = qc_trigger_seed,
    initializeImportBootstrap = function(input, session, experimentPaths, volumes = NULL) {
      expect_identical(input, fake_input)
      expect_identical(session, fake_session)
      expect_identical(experimentPaths, experiment_paths)
      expect_identical(volumes, "volumes-seed")
      call_log <<- c(call_log, "bootstrap")

      list(
        resolvedVolumes = resolved_volumes_seed,
        importFastaPath = import_fasta_seed
      )
    },
    registerImportModalShell = function(input, output, session, resolvedVolumes, importFastaPath) {
      expect_identical(input, fake_input)
      expect_identical(output, fake_output)
      expect_identical(session, fake_session)
      expect_identical(resolvedVolumes, resolved_volumes_seed)
      expect_identical(importFastaPath, import_fasta_seed)
      call_log <<- c(call_log, "modal")
      invisible(NULL)
    },
    registerImportConfirmationObserver = function(input, resolvedVolumes, importFastaPath, workflowData, experimentPaths, session, qcTrigger = NULL) {
      expect_identical(input, fake_input)
      expect_identical(resolvedVolumes, resolved_volumes_seed)
      expect_identical(importFastaPath, import_fasta_seed)
      expect_identical(workflowData, workflow_data)
      expect_identical(experimentPaths, experiment_paths)
      expect_identical(session, fake_session)
      expect_identical(qcTrigger, qc_trigger_seed)
      call_log <<- c(call_log, "importObserver")
      invisible(NULL)
    },
    registerPreviewOutputs = function(output, workflowData) {
      expect_identical(output, fake_output)
      expect_identical(workflowData, workflow_data)
      call_log <<- c(call_log, "preview")
      invisible(output)
    },
    registerBuilderModule = function(workflowData) {
      expect_identical(workflowData, workflow_data)
      call_log <<- c(call_log, "builderModule")
      builder_results_rv
    },
    registerBuilderResultsObserver = function(builderResultsRv, workflowData, experimentPaths, session, qcTrigger = NULL) {
      expect_identical(builderResultsRv, builder_results_rv)
      expect_identical(workflowData, workflow_data)
      expect_identical(experimentPaths, experiment_paths)
      expect_identical(session, fake_session)
      expect_identical(qcTrigger, qc_trigger_seed)
      call_log <<- c(call_log, "builderObserver")
      invisible(NULL)
    }
  )

  expect_identical(result, builder_results_rv)
  expect_identical(
    call_log,
    c("bootstrap", "modal", "importObserver", "preview", "builderModule", "builderObserver")
  )
  expect_identical(import_fasta_calls, 0L)
  expect_identical(builder_results_calls, 0L)
  expect_identical(qc_trigger_calls, 0L)
})

test_that("registerProtDesignServerShells forwards NULL optional inputs through the wrapper shell", {
  call_log <- character()
  workflow_data <- shiny::reactiveValues(
    data_tbl = data.frame(Run = "S1", stringsAsFactors = FALSE),
    config_list = list(globalParameters = list(workflow = "design")),
    column_mapping = list(sample_id = "Run")
  )
  experiment_paths <- list(source_dir = tempfile("prot-design-source-"))
  fake_input <- list(confirm_import = NULL)
  fake_output <- structure(list(), class = "shinyoutput")
  fake_session <- list(ns = function(id) paste0("design-entry-", id))
  builder_results_rv <- function() "builder-results"

  result <- registerProtDesignServerShells(
    input = fake_input,
    output = fake_output,
    session = fake_session,
    workflowData = workflow_data,
    experimentPaths = experiment_paths,
    volumes = NULL,
    qcTrigger = NULL,
    initializeImportBootstrap = function(input, session, experimentPaths, volumes = NULL) {
      expect_identical(input, fake_input)
      expect_identical(session, fake_session)
      expect_identical(experimentPaths, experiment_paths)
      expect_null(volumes)
      call_log <<- c(call_log, "bootstrap")

      list(
        resolvedVolumes = NULL,
        importFastaPath = NULL
      )
    },
    registerImportModalShell = function(input, output, session, resolvedVolumes, importFastaPath) {
      expect_identical(input, fake_input)
      expect_identical(output, fake_output)
      expect_identical(session, fake_session)
      expect_null(resolvedVolumes)
      expect_null(importFastaPath)
      call_log <<- c(call_log, "modal")
      invisible(NULL)
    },
    registerImportConfirmationObserver = function(input, resolvedVolumes, importFastaPath, workflowData, experimentPaths, session, qcTrigger = NULL) {
      expect_identical(input, fake_input)
      expect_null(resolvedVolumes)
      expect_null(importFastaPath)
      expect_identical(workflowData, workflow_data)
      expect_identical(experimentPaths, experiment_paths)
      expect_identical(session, fake_session)
      expect_null(qcTrigger)
      call_log <<- c(call_log, "importObserver")
      invisible(NULL)
    },
    registerPreviewOutputs = function(output, workflowData) {
      expect_identical(output, fake_output)
      expect_identical(workflowData, workflow_data)
      call_log <<- c(call_log, "preview")
      invisible(output)
    },
    registerBuilderModule = function(workflowData) {
      expect_identical(workflowData, workflow_data)
      call_log <<- c(call_log, "builderModule")
      builder_results_rv
    },
    registerBuilderResultsObserver = function(builderResultsRv, workflowData, experimentPaths, session, qcTrigger = NULL) {
      expect_identical(builderResultsRv, builder_results_rv)
      expect_identical(workflowData, workflow_data)
      expect_identical(experimentPaths, experiment_paths)
      expect_identical(session, fake_session)
      expect_null(qcTrigger)
      call_log <<- c(call_log, "builderObserver")
      invisible(NULL)
    }
  )

  expect_identical(result, builder_results_rv)
  expect_identical(
    call_log,
    c("bootstrap", "modal", "importObserver", "preview", "builderModule", "builderObserver")
  )
})

test_that("registerProtDesignServerShells uses the current default callback seeds from the wrapper environment", {
  call_log <- character()
  workflow_data <- shiny::reactiveValues(
    data_tbl = data.frame(Run = "S1", stringsAsFactors = FALSE),
    config_list = list(globalParameters = list(workflow = "design")),
    column_mapping = list(sample_id = "Run")
  )
  experiment_paths <- list(source_dir = tempfile("prot-design-source-"))
  fake_input <- list(confirm_import = NULL)
  fake_output <- structure(list(), class = "shinyoutput")
  fake_session <- list(ns = function(id) paste0("design-entry-", id))
  resolved_volumes <- c(Home = tempdir())
  import_fasta_seed <- function() "/tmp/prot-design-selected.fasta"
  builder_results_rv <- function() "builder-results"
  server_shells <- registerProtDesignServerShells

  environment(server_shells) <- list2env(
    list(
      initializeProtDesignImportBootstrap = function(input, session, experimentPaths, volumes = NULL) {
        expect_identical(input, fake_input)
        expect_identical(session, fake_session)
        expect_identical(experimentPaths, experiment_paths)
        expect_identical(volumes, "volumes-seed")
        call_log <<- c(call_log, "bootstrap")

        list(
          resolvedVolumes = resolved_volumes,
          importFastaPath = import_fasta_seed
        )
      },
      registerProtDesignImportModalShell = function(input, output, session, resolvedVolumes, importFastaPath) {
        expect_identical(input, fake_input)
        expect_identical(output, fake_output)
        expect_identical(session, fake_session)
        expect_identical(resolvedVolumes, resolved_volumes)
        expect_identical(importFastaPath, import_fasta_seed)
        call_log <<- c(call_log, "modal")
        invisible(NULL)
      },
      registerProtDesignImportConfirmationObserver = function(input, resolvedVolumes, importFastaPath, workflowData, experimentPaths, session, qcTrigger = NULL) {
        expect_identical(input, fake_input)
        expect_identical(resolvedVolumes, resolved_volumes)
        expect_identical(importFastaPath, import_fasta_seed)
        expect_identical(workflowData, workflow_data)
        expect_identical(experimentPaths, experiment_paths)
        expect_identical(session, fake_session)
        expect_identical(qcTrigger, "qc-seed")
        call_log <<- c(call_log, "importObserver")
        invisible(NULL)
      },
      registerProtDesignPreviewOutputs = function(output, workflowData) {
        expect_identical(output, fake_output)
        expect_identical(workflowData, workflow_data)
        call_log <<- c(call_log, "preview")
        invisible(output)
      },
      registerProtDesignBuilderModule = function(workflowData) {
        expect_identical(workflowData, workflow_data)
        call_log <<- c(call_log, "builderModule")
        builder_results_rv
      },
      registerProtDesignBuilderResultsObserver = function(builderResultsRv, workflowData, experimentPaths, session, qcTrigger = NULL) {
        expect_identical(builderResultsRv, builder_results_rv)
        expect_identical(workflowData, workflow_data)
        expect_identical(experimentPaths, experiment_paths)
        expect_identical(session, fake_session)
        expect_identical(qcTrigger, "qc-seed")
        call_log <<- c(call_log, "builderObserver")
        invisible(NULL)
      }
    ),
    parent = environment(registerProtDesignServerShells)
  )

  result <- server_shells(
    input = fake_input,
    output = fake_output,
    session = fake_session,
    workflowData = workflow_data,
    experimentPaths = experiment_paths,
    volumes = "volumes-seed",
    qcTrigger = "qc-seed"
  )

  expect_identical(result, builder_results_rv)
  expect_identical(
    call_log,
    c("bootstrap", "modal", "importObserver", "preview", "builderModule", "builderObserver")
  )
})

test_that("mod_prot_design_server wires the wrapper helper seams in order", {
  call_log <- character()
  workflow_data <- shiny::reactiveValues(
    data_tbl = data.frame(Run = "S1", stringsAsFactors = FALSE),
    config_list = list(globalParameters = list(workflow = "design")),
    column_mapping = list(sample_id = "Run")
  )
  experiment_paths <- list(source_dir = tempfile("prot-design-source-"))
  builder_results_rv <- function() "builder-results"

  local_mocked_bindings(
    initializeProtDesignImportBootstrap = function(input, session, experimentPaths, volumes = NULL) {
      expect_identical(experimentPaths, experiment_paths)
      expect_identical(volumes, "volumes-seed")
      call_log <<- c(call_log, "bootstrap")

      list(
        resolvedVolumes = c(Home = tempdir()),
        importFastaPath = function() "/tmp/prot-design-selected.fasta"
      )
    },
    registerProtDesignImportModalShell = function(input, output, session, resolvedVolumes, importFastaPath) {
      expect_identical(resolvedVolumes, c(Home = tempdir()))
      expect_identical(importFastaPath(), "/tmp/prot-design-selected.fasta")
      call_log <<- c(call_log, "modal")
      invisible(NULL)
    },
    registerProtDesignImportConfirmationObserver = function(input, resolvedVolumes, importFastaPath, workflowData, experimentPaths, session, qcTrigger = NULL) {
      expect_identical(resolvedVolumes, c(Home = tempdir()))
      expect_identical(importFastaPath(), "/tmp/prot-design-selected.fasta")
      expect_identical(workflowData, workflow_data)
      expect_identical(experimentPaths, experiment_paths)
      expect_true(qcTrigger)
      call_log <<- c(call_log, "importObserver")
      invisible(NULL)
    },
    registerProtDesignPreviewOutputs = function(output, workflowData) {
      expect_identical(workflowData, workflow_data)
      call_log <<- c(call_log, "preview")
      invisible(output)
    },
    registerProtDesignBuilderModule = function(workflowData) {
      expect_identical(workflowData, workflow_data)
      call_log <<- c(call_log, "builderModule")
      builder_results_rv
    },
    registerProtDesignBuilderResultsObserver = function(builderResultsRv, workflowData, experimentPaths, session, qcTrigger = NULL) {
      expect_identical(builderResultsRv, builder_results_rv)
      expect_identical(workflowData, workflow_data)
      expect_identical(experimentPaths, experiment_paths)
      expect_true(qcTrigger)
      call_log <<- c(call_log, "builderObserver")
      invisible(NULL)
    },
    .env = environment(mod_prot_design_server)
  )

  testServer(
    mod_prot_design_server,
    args = list(
      workflow_data = workflow_data,
      experiment_paths = experiment_paths,
      volumes = "volumes-seed",
      qc_trigger = TRUE
    ),
    {
      expect_identical(
        call_log,
        c("bootstrap", "modal", "importObserver", "preview", "builderModule", "builderObserver")
      )
    }
  )
})

test_that("mod_prot_design_server forwards the wrapper entry id into moduleServer and returns its result", {
  call_log <- character()
  workflow_data <- shiny::reactiveValues(
    data_tbl = data.frame(Run = "S1", stringsAsFactors = FALSE),
    config_list = list(globalParameters = list(workflow = "design")),
    column_mapping = list(sample_id = "Run")
  )
  experiment_paths <- list(source_dir = tempfile("prot-design-source-"))
  builder_results_rv <- function() "builder-results"
  fake_input <- list(confirm_import = NULL)
  fake_output <- structure(list(), class = "shinyoutput")
  fake_session <- list(ns = function(id) paste0("design-entry-", id))
  entry_id <- "design-entry"

  local_mocked_bindings(
    initializeProtDesignImportBootstrap = function(input, session, experimentPaths, volumes = NULL) {
      expect_identical(input, fake_input)
      expect_identical(session, fake_session)
      expect_identical(experimentPaths, experiment_paths)
      expect_identical(volumes, "volumes-seed")
      call_log <<- c(call_log, "bootstrap")

      list(
        resolvedVolumes = c(Home = tempdir()),
        importFastaPath = function() "/tmp/prot-design-selected.fasta"
      )
    },
    registerProtDesignImportModalShell = function(input, output, session, resolvedVolumes, importFastaPath) {
      expect_identical(input, fake_input)
      expect_identical(output, fake_output)
      expect_identical(session, fake_session)
      expect_identical(resolvedVolumes, c(Home = tempdir()))
      expect_identical(importFastaPath(), "/tmp/prot-design-selected.fasta")
      call_log <<- c(call_log, "modal")
      invisible(NULL)
    },
    registerProtDesignImportConfirmationObserver = function(input, resolvedVolumes, importFastaPath, workflowData, experimentPaths, session, qcTrigger = NULL) {
      expect_identical(input, fake_input)
      expect_identical(resolvedVolumes, c(Home = tempdir()))
      expect_identical(importFastaPath(), "/tmp/prot-design-selected.fasta")
      expect_identical(workflowData, workflow_data)
      expect_identical(experimentPaths, experiment_paths)
      expect_identical(session, fake_session)
      expect_identical(qcTrigger, "qc-seed")
      call_log <<- c(call_log, "importObserver")
      invisible(NULL)
    },
    registerProtDesignPreviewOutputs = function(output, workflowData) {
      expect_identical(output, fake_output)
      expect_identical(workflowData, workflow_data)
      call_log <<- c(call_log, "preview")
      invisible(output)
    },
    registerProtDesignBuilderModule = function(workflowData) {
      expect_identical(workflowData, workflow_data)
      call_log <<- c(call_log, "builderModule")
      builder_results_rv
    },
    registerProtDesignBuilderResultsObserver = function(builderResultsRv, workflowData, experimentPaths, session, qcTrigger = NULL) {
      expect_identical(builderResultsRv, builder_results_rv)
      expect_identical(workflowData, workflow_data)
      expect_identical(experimentPaths, experiment_paths)
      expect_identical(session, fake_session)
      expect_identical(qcTrigger, "qc-seed")
      call_log <<- c(call_log, "builderObserver")
      invisible(NULL)
    },
    .env = environment(mod_prot_design_server)
  )
  local_mocked_bindings(
    moduleServer = function(id, module, ...) {
      expect_identical(id, entry_id)
      call_log <<- c(call_log, "moduleServer")
      module(fake_input, fake_output, fake_session)
      "moduleServer-result"
    },
    .package = "shiny"
  )

  result <- mod_prot_design_server(
    entry_id,
    workflow_data = workflow_data,
    experiment_paths = experiment_paths,
    volumes = "volumes-seed",
    qc_trigger = "qc-seed"
  )

  expect_identical(result, "moduleServer-result")
  expect_identical(
    call_log,
    c("moduleServer", "bootstrap", "modal", "importObserver", "preview", "builderModule", "builderObserver")
  )
})

test_that("mod_prot_design_server forwards NULL optional inputs through the wrapper seams", {
  call_log <- character()
  workflow_data <- shiny::reactiveValues(
    data_tbl = data.frame(Run = "S1", stringsAsFactors = FALSE),
    config_list = list(globalParameters = list(workflow = "design")),
    column_mapping = list(sample_id = "Run")
  )
  experiment_paths <- list(source_dir = tempfile("prot-design-source-"))

  local_mocked_bindings(
    initializeProtDesignImportBootstrap = function(input, session, experimentPaths, volumes = NULL) {
      expect_identical(experimentPaths, experiment_paths)
      expect_null(volumes)
      call_log <<- c(call_log, "bootstrap")

      list(
        resolvedVolumes = c(Home = tempdir()),
        importFastaPath = function() NULL
      )
    },
    registerProtDesignImportModalShell = function(input, output, session, resolvedVolumes, importFastaPath) {
      expect_identical(resolvedVolumes, c(Home = tempdir()))
      expect_null(importFastaPath())
      call_log <<- c(call_log, "modal")
      invisible(NULL)
    },
    registerProtDesignImportConfirmationObserver = function(input, resolvedVolumes, importFastaPath, workflowData, experimentPaths, session, qcTrigger = NULL) {
      expect_identical(resolvedVolumes, c(Home = tempdir()))
      expect_null(importFastaPath())
      expect_identical(workflowData, workflow_data)
      expect_identical(experimentPaths, experiment_paths)
      expect_null(qcTrigger)
      call_log <<- c(call_log, "importObserver")
      invisible(NULL)
    },
    registerProtDesignPreviewOutputs = function(output, workflowData) {
      expect_identical(workflowData, workflow_data)
      call_log <<- c(call_log, "preview")
      invisible(output)
    },
    registerProtDesignBuilderModule = function(workflowData) {
      expect_identical(workflowData, workflow_data)
      call_log <<- c(call_log, "builderModule")
      function() NULL
    },
    registerProtDesignBuilderResultsObserver = function(builderResultsRv, workflowData, experimentPaths, session, qcTrigger = NULL) {
      expect_true(is.function(builderResultsRv))
      expect_identical(workflowData, workflow_data)
      expect_identical(experimentPaths, experiment_paths)
      expect_null(qcTrigger)
      call_log <<- c(call_log, "builderObserver")
      invisible(NULL)
    },
    .env = environment(mod_prot_design_server)
  )

  testServer(
    mod_prot_design_server,
    args = list(
      workflow_data = workflow_data,
      experiment_paths = experiment_paths
    ),
    {
      expect_identical(
        call_log,
        c("bootstrap", "modal", "importObserver", "preview", "builderModule", "builderObserver")
      )
    }
  )
})

test_that("mod_prot_design_server preserves bootstrap import handoff objects across wrapper seams", {
  call_log <- character()
  workflow_data <- shiny::reactiveValues(
    data_tbl = data.frame(Run = "S1", stringsAsFactors = FALSE),
    config_list = list(globalParameters = list(workflow = "design")),
    column_mapping = list(sample_id = "Run")
  )
  experiment_paths <- list(source_dir = tempfile("prot-design-source-"))
  resolved_volumes_seed <- c(
    "Project Base Dir" = tempfile("prot-design-base-"),
    Home = tempfile("prot-design-home-")
  )
  import_fasta_seed <- local({
    selected_path <- "/tmp/prot-design-selected.fasta"
    function(new_value) {
      if (!missing(new_value)) {
        selected_path <<- new_value
      }
      selected_path
    }
  })

  local_mocked_bindings(
    initializeProtDesignImportBootstrap = function(input, session, experimentPaths, volumes = NULL) {
      expect_identical(experimentPaths, experiment_paths)
      expect_identical(volumes, "volumes-seed")
      call_log <<- c(call_log, "bootstrap")

      list(
        resolvedVolumes = resolved_volumes_seed,
        importFastaPath = import_fasta_seed
      )
    },
    registerProtDesignImportModalShell = function(input, output, session, resolvedVolumes, importFastaPath) {
      expect_identical(resolvedVolumes, resolved_volumes_seed)
      expect_identical(importFastaPath, import_fasta_seed)
      call_log <<- c(call_log, "modal")
      invisible(NULL)
    },
    registerProtDesignImportConfirmationObserver = function(input, resolvedVolumes, importFastaPath, workflowData, experimentPaths, session, qcTrigger = NULL) {
      expect_identical(resolvedVolumes, resolved_volumes_seed)
      expect_identical(importFastaPath, import_fasta_seed)
      expect_identical(workflowData, workflow_data)
      expect_identical(experimentPaths, experiment_paths)
      expect_true(qcTrigger)
      call_log <<- c(call_log, "importObserver")
      invisible(NULL)
    },
    registerProtDesignPreviewOutputs = function(output, workflowData) {
      expect_identical(workflowData, workflow_data)
      call_log <<- c(call_log, "preview")
      invisible(output)
    },
    registerProtDesignBuilderModule = function(workflowData) {
      expect_identical(workflowData, workflow_data)
      call_log <<- c(call_log, "builderModule")
      function() NULL
    },
    registerProtDesignBuilderResultsObserver = function(builderResultsRv, workflowData, experimentPaths, session, qcTrigger = NULL) {
      expect_true(is.function(builderResultsRv))
      expect_identical(workflowData, workflow_data)
      expect_identical(experimentPaths, experiment_paths)
      expect_true(qcTrigger)
      call_log <<- c(call_log, "builderObserver")
      invisible(NULL)
    },
    .env = environment(mod_prot_design_server)
  )

  testServer(
    mod_prot_design_server,
    args = list(
      workflow_data = workflow_data,
      experiment_paths = experiment_paths,
      volumes = "volumes-seed",
      qc_trigger = TRUE
    ),
    {
      expect_identical(
        call_log,
        c("bootstrap", "modal", "importObserver", "preview", "builderModule", "builderObserver")
      )
    }
  )
})

test_that("mod_prot_design_server reuses moduleServer IO objects across wrapper seams", {
  call_log <- character()
  workflow_data <- shiny::reactiveValues(
    data_tbl = data.frame(Run = "S1", stringsAsFactors = FALSE),
    config_list = list(globalParameters = list(workflow = "design")),
    column_mapping = list(sample_id = "Run")
  )
  experiment_paths <- list(source_dir = tempfile("prot-design-source-"))
  input_ref <- NULL
  output_ref <- NULL
  session_ref <- NULL

  local_mocked_bindings(
    initializeProtDesignImportBootstrap = function(input, session, experimentPaths, volumes = NULL) {
      expect_identical(experimentPaths, experiment_paths)
      expect_identical(volumes, "volumes-seed")
      input_ref <<- input
      session_ref <<- session
      call_log <<- c(call_log, "bootstrap")

      list(
        resolvedVolumes = c(Home = tempdir()),
        importFastaPath = function() NULL
      )
    },
    registerProtDesignImportModalShell = function(input, output, session, resolvedVolumes, importFastaPath) {
      expect_identical(input, input_ref)
      expect_identical(session, session_ref)
      expect_identical(resolvedVolumes, c(Home = tempdir()))
      expect_null(importFastaPath())
      output_ref <<- output
      call_log <<- c(call_log, "modal")
      invisible(NULL)
    },
    registerProtDesignImportConfirmationObserver = function(input, resolvedVolumes, importFastaPath, workflowData, experimentPaths, session, qcTrigger = NULL) {
      expect_identical(input, input_ref)
      expect_identical(session, session_ref)
      expect_identical(resolvedVolumes, c(Home = tempdir()))
      expect_null(importFastaPath())
      expect_identical(workflowData, workflow_data)
      expect_identical(experimentPaths, experiment_paths)
      expect_true(qcTrigger)
      call_log <<- c(call_log, "importObserver")
      invisible(NULL)
    },
    registerProtDesignPreviewOutputs = function(output, workflowData) {
      expect_identical(output, output_ref)
      expect_identical(workflowData, workflow_data)
      call_log <<- c(call_log, "preview")
      invisible(output)
    },
    registerProtDesignBuilderModule = function(workflowData) {
      expect_identical(workflowData, workflow_data)
      call_log <<- c(call_log, "builderModule")
      function() NULL
    },
    registerProtDesignBuilderResultsObserver = function(builderResultsRv, workflowData, experimentPaths, session, qcTrigger = NULL) {
      expect_true(is.function(builderResultsRv))
      expect_identical(workflowData, workflow_data)
      expect_identical(experimentPaths, experiment_paths)
      expect_identical(session, session_ref)
      expect_true(qcTrigger)
      call_log <<- c(call_log, "builderObserver")
      invisible(NULL)
    },
    .env = environment(mod_prot_design_server)
  )

  testServer(
    mod_prot_design_server,
    args = list(
      workflow_data = workflow_data,
      experiment_paths = experiment_paths,
      volumes = "volumes-seed",
      qc_trigger = TRUE
    ),
    {
      expect_identical(
        call_log,
        c("bootstrap", "modal", "importObserver", "preview", "builderModule", "builderObserver")
      )
    }
  )
})

test_that("mod_prot_design_server forwards the exact qc trigger callback through observer seams", {
  call_log <- character()
  workflow_data <- shiny::reactiveValues(
    data_tbl = data.frame(Run = "S1", stringsAsFactors = FALSE),
    config_list = list(globalParameters = list(workflow = "design")),
    column_mapping = list(sample_id = "Run")
  )
  experiment_paths <- list(source_dir = tempfile("prot-design-source-"))
  qc_trigger_seed <- local({
    trigger_count <- 0L
    function(new_value) {
      if (!missing(new_value)) {
        trigger_count <<- trigger_count + 1L
      }
      trigger_count
    }
  })
  builder_results_rv <- function() NULL

  local_mocked_bindings(
    initializeProtDesignImportBootstrap = function(input, session, experimentPaths, volumes = NULL) {
      expect_identical(experimentPaths, experiment_paths)
      expect_identical(volumes, "volumes-seed")
      call_log <<- c(call_log, "bootstrap")

      list(
        resolvedVolumes = c(Home = tempdir()),
        importFastaPath = function() NULL
      )
    },
    registerProtDesignImportModalShell = function(input, output, session, resolvedVolumes, importFastaPath) {
      expect_identical(resolvedVolumes, c(Home = tempdir()))
      expect_null(importFastaPath())
      call_log <<- c(call_log, "modal")
      invisible(NULL)
    },
    registerProtDesignImportConfirmationObserver = function(input, resolvedVolumes, importFastaPath, workflowData, experimentPaths, session, qcTrigger = NULL) {
      expect_identical(resolvedVolumes, c(Home = tempdir()))
      expect_null(importFastaPath())
      expect_identical(workflowData, workflow_data)
      expect_identical(experimentPaths, experiment_paths)
      expect_identical(qcTrigger, qc_trigger_seed)
      call_log <<- c(call_log, "importObserver")
      invisible(NULL)
    },
    registerProtDesignPreviewOutputs = function(output, workflowData) {
      expect_identical(workflowData, workflow_data)
      call_log <<- c(call_log, "preview")
      invisible(output)
    },
    registerProtDesignBuilderModule = function(workflowData) {
      expect_identical(workflowData, workflow_data)
      call_log <<- c(call_log, "builderModule")
      builder_results_rv
    },
    registerProtDesignBuilderResultsObserver = function(builderResultsRv, workflowData, experimentPaths, session, qcTrigger = NULL) {
      expect_identical(builderResultsRv, builder_results_rv)
      expect_identical(workflowData, workflow_data)
      expect_identical(experimentPaths, experiment_paths)
      expect_identical(qcTrigger, qc_trigger_seed)
      call_log <<- c(call_log, "builderObserver")
      invisible(NULL)
    },
    .env = environment(mod_prot_design_server)
  )

  testServer(
    mod_prot_design_server,
    args = list(
      workflow_data = workflow_data,
      experiment_paths = experiment_paths,
      volumes = "volumes-seed",
      qc_trigger = qc_trigger_seed
    ),
    {
      expect_identical(
        call_log,
        c("bootstrap", "modal", "importObserver", "preview", "builderModule", "builderObserver")
      )
    }
  )
})

test_that("mod_prot_design_server passes callback seeds opaquely through wrapper seams", {
  call_log <- character()
  workflow_data <- shiny::reactiveValues(
    data_tbl = data.frame(Run = "S1", stringsAsFactors = FALSE),
    config_list = list(globalParameters = list(workflow = "design")),
    column_mapping = list(sample_id = "Run")
  )
  experiment_paths <- list(source_dir = tempfile("prot-design-source-"))
  import_fasta_calls <- 0L
  builder_results_calls <- 0L
  qc_trigger_calls <- 0L
  import_fasta_seed <- function(...) {
    import_fasta_calls <<- import_fasta_calls + 1L
    "/tmp/prot-design-selected.fasta"
  }
  builder_results_rv <- function(...) {
    builder_results_calls <<- builder_results_calls + 1L
    "builder-results"
  }
  qc_trigger_seed <- function(...) {
    qc_trigger_calls <<- qc_trigger_calls + 1L
    invisible(NULL)
  }

  local_mocked_bindings(
    initializeProtDesignImportBootstrap = function(input, session, experimentPaths, volumes = NULL) {
      expect_identical(experimentPaths, experiment_paths)
      expect_identical(volumes, "volumes-seed")
      call_log <<- c(call_log, "bootstrap")

      list(
        resolvedVolumes = c(Home = tempdir()),
        importFastaPath = import_fasta_seed
      )
    },
    registerProtDesignImportModalShell = function(input, output, session, resolvedVolumes, importFastaPath) {
      expect_identical(resolvedVolumes, c(Home = tempdir()))
      expect_identical(importFastaPath, import_fasta_seed)
      call_log <<- c(call_log, "modal")
      invisible(NULL)
    },
    registerProtDesignImportConfirmationObserver = function(input, resolvedVolumes, importFastaPath, workflowData, experimentPaths, session, qcTrigger = NULL) {
      expect_identical(resolvedVolumes, c(Home = tempdir()))
      expect_identical(importFastaPath, import_fasta_seed)
      expect_identical(workflowData, workflow_data)
      expect_identical(experimentPaths, experiment_paths)
      expect_identical(qcTrigger, qc_trigger_seed)
      call_log <<- c(call_log, "importObserver")
      invisible(NULL)
    },
    registerProtDesignPreviewOutputs = function(output, workflowData) {
      expect_identical(workflowData, workflow_data)
      call_log <<- c(call_log, "preview")
      invisible(output)
    },
    registerProtDesignBuilderModule = function(workflowData) {
      expect_identical(workflowData, workflow_data)
      call_log <<- c(call_log, "builderModule")
      builder_results_rv
    },
    registerProtDesignBuilderResultsObserver = function(builderResultsRv, workflowData, experimentPaths, session, qcTrigger = NULL) {
      expect_identical(builderResultsRv, builder_results_rv)
      expect_identical(workflowData, workflow_data)
      expect_identical(experimentPaths, experiment_paths)
      expect_identical(qcTrigger, qc_trigger_seed)
      call_log <<- c(call_log, "builderObserver")
      invisible(NULL)
    },
    .env = environment(mod_prot_design_server)
  )

  testServer(
    mod_prot_design_server,
    args = list(
      workflow_data = workflow_data,
      experiment_paths = experiment_paths,
      volumes = "volumes-seed",
      qc_trigger = qc_trigger_seed
    ),
    {
      expect_identical(
        call_log,
        c("bootstrap", "modal", "importObserver", "preview", "builderModule", "builderObserver")
      )
    }
  )

  expect_identical(import_fasta_calls, 0L)
  expect_identical(builder_results_calls, 0L)
  expect_identical(qc_trigger_calls, 0L)
})

test_that("mod_prot_design_server uses the current default wrapper shell helper from its environment", {
  call_log <- character()
  workflow_data <- shiny::reactiveValues(
    data_tbl = data.frame(Run = "S1", stringsAsFactors = FALSE),
    config_list = list(globalParameters = list(workflow = "design")),
    column_mapping = list(sample_id = "Run")
  )
  experiment_paths <- list(source_dir = tempfile("prot-design-source-"))
  fake_input <- list(confirm_import = NULL)
  fake_output <- structure(list(), class = "shinyoutput")
  fake_session <- list(ns = function(id) paste0("design-entry-", id))
  entry_id <- "design-entry"
  server_fn <- mod_prot_design_server

  environment(server_fn) <- list2env(
    list(
      registerProtDesignServerShells = function(input, output, session, workflowData, experimentPaths, volumes = NULL, qcTrigger = NULL) {
        expect_identical(input, fake_input)
        expect_identical(output, fake_output)
        expect_identical(session, fake_session)
        expect_identical(workflowData, workflow_data)
        expect_identical(experimentPaths, experiment_paths)
        expect_identical(volumes, "volumes-seed")
        expect_identical(qcTrigger, "qc-seed")
        call_log <<- c(call_log, "serverShells")
        "wrapper-shell-result"
      }
    ),
    parent = environment(mod_prot_design_server)
  )

  local_mocked_bindings(
    moduleServer = function(id, module, ...) {
      expect_identical(id, entry_id)
      call_log <<- c(call_log, "moduleServer")
      expect_identical(
        module(fake_input, fake_output, fake_session),
        "wrapper-shell-result"
      )
      "moduleServer-result"
    },
    .package = "shiny"
  )

  result <- server_fn(
    entry_id,
    workflow_data = workflow_data,
    experiment_paths = experiment_paths,
    volumes = "volumes-seed",
    qc_trigger = "qc-seed"
  )

  expect_identical(result, "moduleServer-result")
  expect_identical(call_log, c("moduleServer", "serverShells"))
})

test_that("mod_prot_design_server keeps the wrapper entry logging ordered around the module shell", {
  call_log <- character()
  workflow_data <- shiny::reactiveValues(
    data_tbl = data.frame(Run = "S1", stringsAsFactors = FALSE),
    config_list = list(globalParameters = list(workflow = "design")),
    column_mapping = list(sample_id = "Run")
  )
  experiment_paths <- list(source_dir = tempfile("prot-design-source-"))
  fake_input <- list(confirm_import = NULL)
  fake_output <- structure(list(), class = "shinyoutput")
  fake_session <- list(ns = function(id) paste0("design-entry-", id))
  entry_id <- "design-entry"
  builder_results_rv <- function() NULL

  local_mocked_bindings(
    initializeProtDesignImportBootstrap = function(input, session, experimentPaths, volumes = NULL) {
      expect_identical(input, fake_input)
      expect_identical(session, fake_session)
      expect_identical(experimentPaths, experiment_paths)
      expect_identical(volumes, "volumes-seed")
      call_log <<- c(call_log, "bootstrap")

      list(
        resolvedVolumes = c(Home = tempdir()),
        importFastaPath = function() NULL
      )
    },
    registerProtDesignImportModalShell = function(input, output, session, resolvedVolumes, importFastaPath) {
      expect_identical(input, fake_input)
      expect_identical(output, fake_output)
      expect_identical(session, fake_session)
      expect_identical(resolvedVolumes, c(Home = tempdir()))
      expect_null(importFastaPath())
      call_log <<- c(call_log, "modal")
      invisible(NULL)
    },
    registerProtDesignImportConfirmationObserver = function(input, resolvedVolumes, importFastaPath, workflowData, experimentPaths, session, qcTrigger = NULL) {
      expect_identical(input, fake_input)
      expect_identical(resolvedVolumes, c(Home = tempdir()))
      expect_null(importFastaPath())
      expect_identical(workflowData, workflow_data)
      expect_identical(experimentPaths, experiment_paths)
      expect_identical(session, fake_session)
      expect_identical(qcTrigger, "qc-seed")
      call_log <<- c(call_log, "importObserver")
      invisible(NULL)
    },
    registerProtDesignPreviewOutputs = function(output, workflowData) {
      expect_identical(output, fake_output)
      expect_identical(workflowData, workflow_data)
      call_log <<- c(call_log, "preview")
      invisible(output)
    },
    registerProtDesignBuilderModule = function(workflowData) {
      expect_identical(workflowData, workflow_data)
      call_log <<- c(call_log, "builderModule")
      builder_results_rv
    },
    registerProtDesignBuilderResultsObserver = function(builderResultsRv, workflowData, experimentPaths, session, qcTrigger = NULL) {
      expect_identical(builderResultsRv, builder_results_rv)
      expect_identical(workflowData, workflow_data)
      expect_identical(experimentPaths, experiment_paths)
      expect_identical(session, fake_session)
      expect_identical(qcTrigger, "qc-seed")
      call_log <<- c(call_log, "builderObserver")
      invisible(NULL)
    },
    .env = environment(mod_prot_design_server)
  )
  local_mocked_bindings(
    moduleServer = function(id, module, ...) {
      expect_identical(id, entry_id)
      call_log <<- c(call_log, "moduleServer")
      module(fake_input, fake_output, fake_session)
      "moduleServer-result"
    },
    .package = "shiny"
  )

  result <- withCallingHandlers(
    mod_prot_design_server(
      entry_id,
      workflow_data = workflow_data,
      experiment_paths = experiment_paths,
      volumes = "volumes-seed",
      qc_trigger = "qc-seed"
    ),
    message = function(cnd) {
      call_log <<- c(
        call_log,
        paste0("message:", sub("\n$", "", conditionMessage(cnd)))
      )
      invokeRestart("muffleMessage")
    }
  )

  expect_identical(result, "moduleServer-result")
  expect_identical(
    call_log,
    c(
      "message:--- Entering mod_prot_design_server ---",
      "message:   mod_prot_design_server Arg: id = design-entry",
      "moduleServer",
      "message:   mod_prot_design_server Step: Inside moduleServer function",
      "bootstrap",
      "modal",
      "importObserver",
      "preview",
      "builderModule",
      "builderObserver"
    )
  )
})

test_that("mod_prot_design_ui renders the builder module when it is available", {
  ui_with_builder <- mod_prot_design_ui
  environment(ui_with_builder) <- list2env(
    list(
      exists = function(x, where = -1, inherits = TRUE) {
        if (identical(x, "mod_prot_design_builder_ui")) {
          TRUE
        } else {
          base::exists(x, where = where, inherits = inherits)
        }
      },
      mod_prot_design_builder_ui = function(id) {
        shiny::div(id = paste0(id, "-builder-ready"), "Builder ready")
      }
    ),
    parent = environment(mod_prot_design_ui)
  )

  html <- htmltools::renderTags(ui_with_builder("design"))$html

  expect_match(html, "Design Matrix Builder", fixed = TRUE)
  expect_match(html, "design-builder-builder-ready", fixed = TRUE)
  expect_match(html, "Builder ready", fixed = TRUE)
  expect_false(grepl("Design builder module not loaded", html, fixed = TRUE))
})

test_that("mod_prot_design_ui falls back when the builder module is unavailable", {
  ui_without_builder <- mod_prot_design_ui
  environment(ui_without_builder) <- list2env(
    list(
      exists = function(x, where = -1, inherits = TRUE) {
        if (identical(x, "mod_prot_design_builder_ui")) {
          FALSE
        } else {
          base::exists(x, where = where, inherits = inherits)
        }
      }
    ),
    parent = environment(mod_prot_design_ui)
  )

  html <- htmltools::renderTags(ui_without_builder("design"))$html

  expect_match(html, "Design Matrix Builder", fixed = TRUE)
  expect_match(html, "Design builder module not loaded", fixed = TRUE)
})

test_that("mod_prot_design_ui preserves the import and preview wrapper scaffold", {
  ui_with_builder <- mod_prot_design_ui
  environment(ui_with_builder) <- list2env(
    list(
      exists = function(x, where = -1, inherits = TRUE) {
        if (identical(x, "mod_prot_design_builder_ui")) {
          TRUE
        } else {
          base::exists(x, where = where, inherits = inherits)
        }
      },
      mod_prot_design_builder_ui = function(id) {
        shiny::div(id = paste0(id, "-builder-ready"), "Builder ready")
      }
    ),
    parent = environment(mod_prot_design_ui)
  )

  html <- htmltools::renderTags(ui_with_builder("design"))$html

  expect_match(html, "id=\"design-show_import_modal\"", fixed = TRUE)
  expect_match(html, "Import Existing Design", fixed = TRUE)
  expect_match(html, "Saved Results Preview", fixed = TRUE)
  expect_match(html, "Current Design Matrix", fixed = TRUE)
  expect_match(html, "Defined Contrasts", fixed = TRUE)
  expect_match(html, "design-design_matrix_preview", fixed = TRUE)
  expect_match(html, "design-contrasts_preview", fixed = TRUE)
  expect_match(html, "design-data_available", fixed = TRUE)
  expect_match(html, "design-design_matrix_exists", fixed = TRUE)
})

test_that("mod_prot_design_ui preserves namespaced wrapper bindings for non-default ids", {
  ui_with_builder <- mod_prot_design_ui
  environment(ui_with_builder) <- list2env(
    list(
      exists = function(x, where = -1, inherits = TRUE) {
        if (identical(x, "mod_prot_design_builder_ui")) {
          TRUE
        } else {
          base::exists(x, where = where, inherits = inherits)
        }
      },
      mod_prot_design_builder_ui = function(id) {
        shiny::div(id = paste0(id, "-builder-ready"), "Builder ready")
      }
    ),
    parent = environment(mod_prot_design_ui)
  )

  html <- htmltools::renderTags(ui_with_builder("prot-design"))$html

  expect_match(html, "id=\"prot-design-show_import_modal\"", fixed = TRUE)
  expect_match(html, "prot-design-builder-builder-ready", fixed = TRUE)
  expect_match(html, "prot-design-design_matrix_preview", fixed = TRUE)
  expect_match(html, "prot-design-contrasts_preview", fixed = TRUE)
  expect_match(
    html,
    "output\\[&#39;prot-design-data_available&#39;\\]",
    perl = TRUE
  )
  expect_match(
    html,
    "output\\[&#39;prot-design-design_matrix_exists&#39;\\]",
    perl = TRUE
  )
})

test_that("mod_prot_design_ui keeps the fallback scaffold namespaced for non-default ids", {
  ui_without_builder <- mod_prot_design_ui
  environment(ui_without_builder) <- list2env(
    list(
      exists = function(x, where = -1, inherits = TRUE) {
        if (identical(x, "mod_prot_design_builder_ui")) {
          FALSE
        } else {
          base::exists(x, where = where, inherits = inherits)
        }
      }
    ),
    parent = environment(mod_prot_design_ui)
  )

  html <- htmltools::renderTags(ui_without_builder("prot-design"))$html

  expect_match(html, "id=\"prot-design-show_import_modal\"", fixed = TRUE)
  expect_match(html, "Design builder module not loaded", fixed = TRUE)
  expect_match(html, "prot-design-design_matrix_preview", fixed = TRUE)
  expect_match(html, "prot-design-contrasts_preview", fixed = TRUE)
  expect_match(
    html,
    "output\\[&#39;prot-design-design_matrix_exists&#39;\\]",
    perl = TRUE
  )
  expect_match(
    html,
    "data-display-if=\"!output\\[&#39;prot-design-data_available&#39;\\]\"",
    perl = TRUE
  )
})

test_that("mod_prot_design_ui keeps the default fallback scaffold single-instance and namespaced", {
  ui_without_builder <- mod_prot_design_ui
  environment(ui_without_builder) <- list2env(
    list(
      exists = function(x, where = -1, inherits = TRUE) {
        if (identical(x, "mod_prot_design_builder_ui")) {
          FALSE
        } else {
          base::exists(x, where = where, inherits = inherits)
        }
      }
    ),
    parent = environment(mod_prot_design_ui)
  )

  html <- htmltools::renderTags(ui_without_builder("design"))$html

  count_occurrences <- function(pattern) {
    matches <- gregexpr(pattern, html, fixed = TRUE)[[1]]
    if (identical(matches[1], -1L)) {
      return(0L)
    }

    length(matches)
  }

  expect_identical(count_occurrences("Design builder module not loaded"), 1L)
  expect_identical(count_occurrences("id=\"design-show_import_modal\""), 1L)
  expect_identical(count_occurrences("id=\"design-design_matrix_preview\""), 1L)
  expect_identical(count_occurrences("id=\"design-contrasts_preview\""), 1L)
  expect_identical(
    count_occurrences("data-display-if=\"output[&#39;design-data_available&#39;]\""),
    1L
  )
  expect_identical(
    count_occurrences("data-display-if=\"output[&#39;design-design_matrix_exists&#39;]\""),
    1L
  )
  expect_identical(
    count_occurrences("data-display-if=\"!output[&#39;design-data_available&#39;]\""),
    1L
  )
  expect_false(grepl("id=\"show_import_modal\"", html, fixed = TRUE))
  expect_false(grepl("id=\"design_matrix_preview\"", html, fixed = TRUE))
  expect_false(grepl("id=\"contrasts_preview\"", html, fixed = TRUE))
})

test_that("mod_prot_design_ui keeps non-default wrapper ids fully namespaced", {
  ui_with_builder <- mod_prot_design_ui
  environment(ui_with_builder) <- list2env(
    list(
      exists = function(x, where = -1, inherits = TRUE) {
        if (identical(x, "mod_prot_design_builder_ui")) {
          TRUE
        } else {
          base::exists(x, where = where, inherits = inherits)
        }
      },
      mod_prot_design_builder_ui = function(id) {
        shiny::div(id = paste0(id, "-builder-ready"), "Builder ready")
      }
    ),
    parent = environment(mod_prot_design_ui)
  )

  html <- htmltools::renderTags(ui_with_builder("prot-design"))$html

  expect_match(html, "id=\"prot-design-show_import_modal\"", fixed = TRUE)
  expect_match(html, "prot-design-builder-builder-ready", fixed = TRUE)
  expect_match(html, "prot-design-design_matrix_preview", fixed = TRUE)
  expect_match(html, "prot-design-contrasts_preview", fixed = TRUE)
  expect_false(grepl("id=\"show_import_modal\"", html, fixed = TRUE))
  expect_false(grepl("id=\"builder-builder-ready\"", html, fixed = TRUE))
  expect_false(grepl("id=\"design_matrix_preview\"", html, fixed = TRUE))
  expect_false(grepl("id=\"contrasts_preview\"", html, fixed = TRUE))
  expect_false(grepl("output[&#39;data_available&#39;]", html, fixed = TRUE))
  expect_false(grepl("output[&#39;design_matrix_exists&#39;]", html, fixed = TRUE))
})

test_that("mod_prot_design_ui keeps progress-handler hooks global while wrapper ids stay namespaced", {
  html <- htmltools::renderTags(mod_prot_design_ui("prot-design"))$html

  expect_match(html, "id=\"prot-design-show_import_modal\"", fixed = TRUE)
  expect_match(html, "prot-design-design_matrix_preview", fixed = TRUE)
  expect_match(html, "prot-design-contrasts_preview", fixed = TRUE)
  expect_match(html, "uniprot_progress_bar", fixed = TRUE)
  expect_match(html, "uniprot_progress_text", fixed = TRUE)
  expect_match(html, "$('#uniprot_progress_bar').css('width', message.percent + '%');", fixed = TRUE)
  expect_match(html, "$('#uniprot_progress_text').text(message.text);", fixed = TRUE)
  expect_false(grepl("prot-design-uniprot_progress_bar", html, fixed = TRUE))
  expect_false(grepl("prot-design-uniprot_progress_text", html, fixed = TRUE))
})

test_that("mod_prot_design_ui preserves conditional-panel output bindings", {
  ui_with_builder <- mod_prot_design_ui
  environment(ui_with_builder) <- list2env(
    list(
      exists = function(x, where = -1, inherits = TRUE) {
        if (identical(x, "mod_prot_design_builder_ui")) {
          TRUE
        } else {
          base::exists(x, where = where, inherits = inherits)
        }
      },
      mod_prot_design_builder_ui = function(id) {
        shiny::div(id = paste0(id, "-builder-ready"), "Builder ready")
      }
    ),
    parent = environment(mod_prot_design_ui)
  )

  html <- htmltools::renderTags(ui_with_builder("design"))$html

  expect_match(
    html,
    "data-display-if=\"output\\[&#39;design-data_available&#39;\\]\"",
    perl = TRUE
  )
  expect_match(
    html,
    "data-display-if=\"output\\[&#39;design-design_matrix_exists&#39;\\]\"",
    perl = TRUE
  )
  expect_match(
    html,
    "data-display-if=\"!output\\[&#39;design-data_available&#39;\\]\"",
    perl = TRUE
  )
  expect_equal(
    lengths(regmatches(html, gregexpr("data-ns-prefix=\"\"", html, fixed = TRUE))),
    3L
  )
})

test_that("mod_prot_design_ui preserves the progress handler and empty-state guidance scaffold", {
  html <- htmltools::renderTags(mod_prot_design_ui("design"))$html

  expect_match(html, "Shiny.addCustomMessageHandler('updateUniprotProgress'", fixed = TRUE)
  expect_match(html, "$('#uniprot_progress_bar').css('width', message.percent + '%');", fixed = TRUE)
  expect_match(html, "$('#uniprot_progress_text').text(message.text);", fixed = TRUE)
  expect_match(html, "uniprot_progress_bar", fixed = TRUE)
  expect_match(html, "uniprot_progress_text", fixed = TRUE)
  expect_match(html, "btn-info pull-right", fixed = TRUE)
  expect_match(html, "fa-folder-open", fixed = TRUE)
  expect_match(html, "!output[&#39;design-data_available&#39;]", fixed = TRUE)
  expect_match(html, "alert alert-info", fixed = TRUE)
  expect_match(
    html,
    "Please complete the 'Setup &amp; Import' step first. The builder will appear here once data is available.",
    fixed = TRUE
  )
})

test_that("mod_prot_design_ui keeps the default empty-state guidance single-instance and ordered inside the alert shell", {
  ui_with_builder <- mod_prot_design_ui
  environment(ui_with_builder) <- list2env(
    list(
      exists = function(x, where = -1, inherits = TRUE) {
        if (identical(x, "mod_prot_design_builder_ui")) {
          TRUE
        } else {
          base::exists(x, where = where, inherits = inherits)
        }
      },
      mod_prot_design_builder_ui = function(id) {
        shiny::div(id = paste0(id, "-builder-ready"), "Builder ready")
      }
    ),
    parent = environment(mod_prot_design_ui)
  )

  html <- htmltools::renderTags(ui_with_builder("design"))$html

  count_occurrences <- function(pattern) {
    matches <- gregexpr(pattern, html, fixed = TRUE)[[1]]
    if (identical(matches[1], -1L)) {
      return(0L)
    }

    length(matches)
  }

  empty_state_binding_pos <- as.integer(regexpr(
    "data-display-if=\"!output[&#39;design-data_available&#39;]\"",
    html,
    fixed = TRUE
  ))
  alert_pos <- as.integer(regexpr(
    "alert alert-info",
    html,
    fixed = TRUE
  ))
  guidance_pos <- as.integer(regexpr(
    "Please complete the 'Setup &amp; Import' step first. The builder will appear here once data is available.",
    html,
    fixed = TRUE
  ))

  expect_identical(
    count_occurrences("Please complete the 'Setup &amp; Import' step first. The builder will appear here once data is available."),
    1L
  )
  expect_identical(count_occurrences("alert alert-info"), 1L)
  expect_gt(empty_state_binding_pos, 0L)
  expect_gt(alert_pos, empty_state_binding_pos)
  expect_gt(guidance_pos, alert_pos)
  expect_false(grepl("!output[&#39;data_available&#39;]", html, fixed = TRUE))
  expect_false(grepl("Design builder module not loaded", html, fixed = TRUE))
})

test_that("mod_prot_design_ui keeps the non-default empty-state guidance single-instance and ordered inside the alert shell", {
  ui_with_builder <- mod_prot_design_ui
  environment(ui_with_builder) <- list2env(
    list(
      exists = function(x, where = -1, inherits = TRUE) {
        if (identical(x, "mod_prot_design_builder_ui")) {
          TRUE
        } else {
          base::exists(x, where = where, inherits = inherits)
        }
      },
      mod_prot_design_builder_ui = function(id) {
        shiny::div(id = paste0(id, "-builder-ready"), "Builder ready")
      }
    ),
    parent = environment(mod_prot_design_ui)
  )

  html <- htmltools::renderTags(ui_with_builder("prot-design"))$html

  count_occurrences <- function(pattern) {
    matches <- gregexpr(pattern, html, fixed = TRUE)[[1]]
    if (identical(matches[1], -1L)) {
      return(0L)
    }

    length(matches)
  }

  empty_state_binding_pos <- as.integer(regexpr(
    "data-display-if=\"!output[&#39;prot-design-data_available&#39;]\"",
    html,
    fixed = TRUE
  ))
  alert_pos <- as.integer(regexpr(
    "alert alert-info",
    html,
    fixed = TRUE
  ))
  guidance_pos <- as.integer(regexpr(
    "Please complete the 'Setup &amp; Import' step first. The builder will appear here once data is available.",
    html,
    fixed = TRUE
  ))

  expect_identical(
    count_occurrences("Please complete the 'Setup &amp; Import' step first. The builder will appear here once data is available."),
    1L
  )
  expect_identical(count_occurrences("alert alert-info"), 1L)
  expect_gt(empty_state_binding_pos, 0L)
  expect_gt(alert_pos, empty_state_binding_pos)
  expect_gt(guidance_pos, alert_pos)
  expect_false(grepl("!output[&#39;data_available&#39;]", html, fixed = TRUE))
  expect_false(grepl("Design builder module not loaded", html, fixed = TRUE))
})

test_that("mod_prot_design_ui registers the UniProt progress handler only once", {
  html <- htmltools::renderTags(mod_prot_design_ui("prot-design"))$html

  count_occurrences <- function(pattern) {
    matches <- gregexpr(pattern, html, fixed = TRUE)[[1]]
    if (identical(matches[1], -1L)) {
      return(0L)
    }

    length(matches)
  }

  expect_identical(
    count_occurrences("Shiny.addCustomMessageHandler('updateUniprotProgress'"),
    1L
  )
  expect_identical(
    count_occurrences("$('#uniprot_progress_bar').css('width', message.percent + '%');"),
    1L
  )
  expect_identical(
    count_occurrences("$('#uniprot_progress_text').text(message.text);"),
    1L
  )
})

test_that("mod_prot_design_ui keeps the wrapper shell single-instance for non-default ids", {
  html <- htmltools::renderTags(mod_prot_design_ui("prot-design"))$html

  count_occurrences <- function(pattern) {
    matches <- gregexpr(pattern, html, fixed = TRUE)[[1]]
    if (identical(matches[1], -1L)) {
      return(0L)
    }

    length(matches)
  }

  expect_identical(count_occurrences("id=\"prot-design-show_import_modal\""), 1L)
  expect_identical(count_occurrences("Saved Results Preview"), 1L)
  expect_identical(count_occurrences("id=\"prot-design-design_matrix_preview\""), 1L)
  expect_identical(count_occurrences("id=\"prot-design-contrasts_preview\""), 1L)
  expect_identical(
    count_occurrences("data-display-if=\"output[&#39;prot-design-data_available&#39;]\""),
    1L
  )
  expect_identical(
    count_occurrences("data-display-if=\"output[&#39;prot-design-design_matrix_exists&#39;]\""),
    1L
  )
  expect_identical(
    count_occurrences("data-display-if=\"!output[&#39;prot-design-data_available&#39;]\""),
    1L
  )
})

test_that("mod_prot_design_ui keeps the wrapper shell single-instance for the default wrapper id", {
  html <- htmltools::renderTags(mod_prot_design_ui("design"))$html

  count_occurrences <- function(pattern) {
    matches <- gregexpr(pattern, html, fixed = TRUE)[[1]]
    if (identical(matches[1], -1L)) {
      return(0L)
    }

    length(matches)
  }

  expect_identical(count_occurrences("id=\"design-show_import_modal\""), 1L)
  expect_identical(count_occurrences("Saved Results Preview"), 1L)
  expect_identical(count_occurrences("id=\"design-design_matrix_preview\""), 1L)
  expect_identical(count_occurrences("id=\"design-contrasts_preview\""), 1L)
  expect_identical(
    count_occurrences("data-display-if=\"output[&#39;design-data_available&#39;]\""),
    1L
  )
  expect_identical(
    count_occurrences("data-display-if=\"output[&#39;design-design_matrix_exists&#39;]\""),
    1L
  )
  expect_identical(
    count_occurrences("data-display-if=\"!output[&#39;design-data_available&#39;]\""),
    1L
  )
  expect_false(grepl("id=\"show_import_modal\"", html, fixed = TRUE))
  expect_false(grepl("id=\"design_matrix_preview\"", html, fixed = TRUE))
  expect_false(grepl("id=\"contrasts_preview\"", html, fixed = TRUE))
})

test_that("mod_prot_design_ui keeps the saved-results preview headings ordered for the default wrapper id", {
  ui_with_builder <- mod_prot_design_ui
  environment(ui_with_builder) <- list2env(
    list(
      exists = function(x, where = -1, inherits = TRUE) {
        if (identical(x, "mod_prot_design_builder_ui")) {
          TRUE
        } else {
          base::exists(x, where = where, inherits = inherits)
        }
      },
      mod_prot_design_builder_ui = function(id) {
        shiny::div(id = paste0(id, "-builder-ready"), "Builder ready")
      }
    ),
    parent = environment(mod_prot_design_ui)
  )

  html <- htmltools::renderTags(ui_with_builder("design"))$html

  count_occurrences <- function(pattern) {
    matches <- gregexpr(pattern, html, fixed = TRUE)[[1]]
    if (identical(matches[1], -1L)) {
      return(0L)
    }

    length(matches)
  }

  current_heading_pos <- as.integer(regexpr("Current Design Matrix", html, fixed = TRUE))
  design_preview_pos <- as.integer(regexpr("id=\"design-design_matrix_preview\"", html, fixed = TRUE))
  contrasts_heading_pos <- as.integer(regexpr("Defined Contrasts", html, fixed = TRUE))
  contrasts_preview_pos <- as.integer(regexpr("id=\"design-contrasts_preview\"", html, fixed = TRUE))

  expect_identical(count_occurrences("Current Design Matrix"), 1L)
  expect_identical(count_occurrences("Defined Contrasts"), 1L)
  expect_identical(count_occurrences("id=\"design-design_matrix_preview\""), 1L)
  expect_identical(count_occurrences("id=\"design-contrasts_preview\""), 1L)
  expect_gt(design_preview_pos, current_heading_pos)
  expect_gt(contrasts_heading_pos, design_preview_pos)
  expect_gt(contrasts_preview_pos, contrasts_heading_pos)
})

test_that("mod_prot_design_server passes the exact volumes callback seed opaquely into bootstrap", {
  call_log <- character()
  workflow_data <- shiny::reactiveValues(
    data_tbl = data.frame(Run = "S1", stringsAsFactors = FALSE),
    config_list = list(globalParameters = list(workflow = "design")),
    column_mapping = list(sample_id = "Run")
  )
  experiment_paths <- list(source_dir = tempfile("prot-design-source-"))
  volumes_calls <- 0L
  volumes_seed <- function(...) {
    volumes_calls <<- volumes_calls + 1L
    c(Home = "/tmp/should-not-be-used")
  }

  local_mocked_bindings(
    initializeProtDesignImportBootstrap = function(input, session, experimentPaths, volumes = NULL) {
      expect_identical(experimentPaths, experiment_paths)
      expect_identical(volumes, volumes_seed)
      call_log <<- c(call_log, "bootstrap")

      list(
        resolvedVolumes = c(Home = tempdir()),
        importFastaPath = function() NULL
      )
    },
    registerProtDesignImportModalShell = function(input, output, session, resolvedVolumes, importFastaPath) {
      expect_identical(resolvedVolumes, c(Home = tempdir()))
      expect_null(importFastaPath())
      call_log <<- c(call_log, "modal")
      invisible(NULL)
    },
    registerProtDesignImportConfirmationObserver = function(input, resolvedVolumes, importFastaPath, workflowData, experimentPaths, session, qcTrigger = NULL) {
      expect_identical(resolvedVolumes, c(Home = tempdir()))
      expect_null(importFastaPath())
      expect_identical(workflowData, workflow_data)
      expect_identical(experimentPaths, experiment_paths)
      expect_null(qcTrigger)
      call_log <<- c(call_log, "importObserver")
      invisible(NULL)
    },
    registerProtDesignPreviewOutputs = function(output, workflowData) {
      expect_identical(workflowData, workflow_data)
      call_log <<- c(call_log, "preview")
      invisible(output)
    },
    registerProtDesignBuilderModule = function(workflowData) {
      expect_identical(workflowData, workflow_data)
      call_log <<- c(call_log, "builderModule")
      function() NULL
    },
    registerProtDesignBuilderResultsObserver = function(builderResultsRv, workflowData, experimentPaths, session, qcTrigger = NULL) {
      expect_true(is.function(builderResultsRv))
      expect_identical(workflowData, workflow_data)
      expect_identical(experimentPaths, experiment_paths)
      expect_null(qcTrigger)
      call_log <<- c(call_log, "builderObserver")
      invisible(NULL)
    },
    .env = environment(mod_prot_design_server)
  )

  testServer(
    mod_prot_design_server,
    args = list(
      workflow_data = workflow_data,
      experiment_paths = experiment_paths,
      volumes = volumes_seed
    ),
    {
      expect_identical(
        call_log,
        c("bootstrap", "modal", "importObserver", "preview", "builderModule", "builderObserver")
      )
    }
  )

  expect_identical(volumes_calls, 0L)
})

test_that("mod_prot_design_ui keeps the embedded builder binding single-instance for non-default ids", {
  ui_with_builder <- mod_prot_design_ui
  environment(ui_with_builder) <- list2env(
    list(
      exists = function(x, where = -1, inherits = TRUE) {
        if (identical(x, "mod_prot_design_builder_ui")) {
          TRUE
        } else {
          base::exists(x, where = where, inherits = inherits)
        }
      },
      mod_prot_design_builder_ui = function(id) {
        shiny::div(
          id = paste0(id, "-builder-ready"),
          `data-builder-id` = id,
          "Builder ready"
        )
      }
    ),
    parent = environment(mod_prot_design_ui)
  )

  html <- htmltools::renderTags(ui_with_builder("prot-design"))$html

  count_occurrences <- function(pattern) {
    matches <- gregexpr(pattern, html, fixed = TRUE)[[1]]
    if (identical(matches[1], -1L)) {
      return(0L)
    }

    length(matches)
  }

  expect_identical(
    count_occurrences("id=\"prot-design-builder-builder-ready\""),
    1L
  )
  expect_identical(
    count_occurrences("data-builder-id=\"prot-design-builder\""),
    1L
  )
  expect_false(grepl("id=\"builder-builder-ready\"", html, fixed = TRUE))
  expect_false(grepl("data-builder-id=\"builder\"", html, fixed = TRUE))
})

test_that("mod_prot_design_ui forwards the exact namespaced builder id into the embedded module", {
  observed_builder_ids <- character()
  ui_with_builder <- mod_prot_design_ui
  environment(ui_with_builder) <- list2env(
    list(
      exists = function(x, where = -1, inherits = TRUE) {
        if (identical(x, "mod_prot_design_builder_ui")) {
          TRUE
        } else {
          base::exists(x, where = where, inherits = inherits)
        }
      },
      mod_prot_design_builder_ui = function(id) {
        observed_builder_ids <<- c(observed_builder_ids, id)
        shiny::div(
          id = paste0(id, "-builder-ready"),
          `data-builder-id` = id,
          "Builder ready"
        )
      }
    ),
    parent = environment(mod_prot_design_ui)
  )

  html <- htmltools::renderTags(ui_with_builder("design"))$html

  expect_identical(observed_builder_ids, "design-builder")
  expect_match(html, "id=\"design-builder-builder-ready\"", fixed = TRUE)
  expect_match(html, "data-builder-id=\"design-builder\"", fixed = TRUE)
})

test_that("mod_prot_design_ui keeps the embedded builder binding single-instance for the default wrapper id", {
  ui_with_builder <- mod_prot_design_ui
  environment(ui_with_builder) <- list2env(
    list(
      exists = function(x, where = -1, inherits = TRUE) {
        if (identical(x, "mod_prot_design_builder_ui")) {
          TRUE
        } else {
          base::exists(x, where = where, inherits = inherits)
        }
      },
      mod_prot_design_builder_ui = function(id) {
        shiny::div(
          id = paste0(id, "-builder-ready"),
          `data-builder-id` = id,
          "Builder ready"
        )
      }
    ),
    parent = environment(mod_prot_design_ui)
  )

  html <- htmltools::renderTags(ui_with_builder("design"))$html

  count_occurrences <- function(pattern) {
    matches <- gregexpr(pattern, html, fixed = TRUE)[[1]]
    if (identical(matches[1], -1L)) {
      return(0L)
    }

    length(matches)
  }

  expect_identical(
    count_occurrences("id=\"design-builder-builder-ready\""),
    1L
  )
  expect_identical(
    count_occurrences("data-builder-id=\"design-builder\""),
    1L
  )
  expect_false(grepl("Design builder module not loaded", html, fixed = TRUE))
  expect_false(grepl("id=\"builder-builder-ready\"", html, fixed = TRUE))
  expect_false(grepl("data-builder-id=\"builder\"", html, fixed = TRUE))
})

test_that("mod_prot_design_ui keeps default builder-enabled conditional bindings namespaced", {
  ui_with_builder <- mod_prot_design_ui
  environment(ui_with_builder) <- list2env(
    list(
      exists = function(x, where = -1, inherits = TRUE) {
        if (identical(x, "mod_prot_design_builder_ui")) {
          TRUE
        } else {
          base::exists(x, where = where, inherits = inherits)
        }
      },
      mod_prot_design_builder_ui = function(id) {
        shiny::div(
          id = paste0(id, "-builder-ready"),
          `data-builder-id` = id,
          "Builder ready"
        )
      }
    ),
    parent = environment(mod_prot_design_ui)
  )

  html <- htmltools::renderTags(ui_with_builder("design"))$html

  count_occurrences <- function(pattern) {
    matches <- gregexpr(pattern, html, fixed = TRUE)[[1]]
    if (identical(matches[1], -1L)) {
      return(0L)
    }

    length(matches)
  }

  expect_identical(
    count_occurrences("id=\"design-builder-builder-ready\""),
    1L
  )
  expect_identical(
    count_occurrences("data-builder-id=\"design-builder\""),
    1L
  )
  expect_identical(
    count_occurrences("data-display-if=\"output[&#39;design-data_available&#39;]\""),
    1L
  )
  expect_identical(
    count_occurrences("data-display-if=\"output[&#39;design-design_matrix_exists&#39;]\""),
    1L
  )
  expect_identical(
    count_occurrences("data-display-if=\"!output[&#39;design-data_available&#39;]\""),
    1L
  )
  expect_false(grepl("id=\"builder-builder-ready\"", html, fixed = TRUE))
  expect_false(grepl("data-builder-id=\"builder\"", html, fixed = TRUE))
  expect_false(grepl("output[&#39;data_available&#39;]", html, fixed = TRUE))
  expect_false(grepl("output[&#39;design_matrix_exists&#39;]", html, fixed = TRUE))
})

test_that("mod_prot_design_ui keeps non-default builder-enabled conditional bindings on the shared empty namespace prefix", {
  ui_with_builder <- mod_prot_design_ui
  environment(ui_with_builder) <- list2env(
    list(
      exists = function(x, where = -1, inherits = TRUE) {
        if (identical(x, "mod_prot_design_builder_ui")) {
          TRUE
        } else {
          base::exists(x, where = where, inherits = inherits)
        }
      },
      mod_prot_design_builder_ui = function(id) {
        shiny::div(
          id = paste0(id, "-builder-ready"),
          `data-builder-id` = id,
          "Builder ready"
        )
      }
    ),
    parent = environment(mod_prot_design_ui)
  )

  html <- htmltools::renderTags(ui_with_builder("prot-design"))$html

  count_occurrences <- function(pattern) {
    matches <- gregexpr(pattern, html, fixed = TRUE)[[1]]
    if (identical(matches[1], -1L)) {
      return(0L)
    }

    length(matches)
  }

  expect_identical(
    count_occurrences("id=\"prot-design-builder-builder-ready\""),
    1L
  )
  expect_identical(
    count_occurrences("data-builder-id=\"prot-design-builder\""),
    1L
  )
  expect_identical(count_occurrences("data-ns-prefix=\"\""), 3L)
  expect_identical(
    count_occurrences("data-display-if=\"output[&#39;prot-design-data_available&#39;]\""),
    1L
  )
  expect_identical(
    count_occurrences("data-display-if=\"output[&#39;prot-design-design_matrix_exists&#39;]\""),
    1L
  )
  expect_identical(
    count_occurrences("data-display-if=\"!output[&#39;prot-design-data_available&#39;]\""),
    1L
  )
  expect_false(grepl("id=\"builder-builder-ready\"", html, fixed = TRUE))
  expect_false(grepl("data-builder-id=\"builder\"", html, fixed = TRUE))
  expect_false(grepl("output[&#39;data_available&#39;]", html, fixed = TRUE))
  expect_false(grepl("output[&#39;design_matrix_exists&#39;]", html, fixed = TRUE))
})

test_that("mod_prot_design_ui keeps default builder-missing conditional bindings on the shared empty namespace prefix", {
  ui_without_builder <- mod_prot_design_ui
  environment(ui_without_builder) <- list2env(
    list(
      exists = function(x, where = -1, inherits = TRUE) {
        if (identical(x, "mod_prot_design_builder_ui")) {
          FALSE
        } else {
          base::exists(x, where = where, inherits = inherits)
        }
      }
    ),
    parent = environment(mod_prot_design_ui)
  )

  html <- htmltools::renderTags(ui_without_builder("design"))$html

  count_occurrences <- function(pattern) {
    matches <- gregexpr(pattern, html, fixed = TRUE)[[1]]
    if (identical(matches[1], -1L)) {
      return(0L)
    }

    length(matches)
  }

  expect_identical(count_occurrences("Design builder module not loaded"), 1L)
  expect_identical(count_occurrences("data-ns-prefix=\"\""), 3L)
  expect_identical(
    count_occurrences("data-display-if=\"output[&#39;design-data_available&#39;]\""),
    1L
  )
  expect_identical(
    count_occurrences("data-display-if=\"output[&#39;design-design_matrix_exists&#39;]\""),
    1L
  )
  expect_identical(
    count_occurrences("data-display-if=\"!output[&#39;design-data_available&#39;]\""),
    1L
  )
  expect_false(grepl("design-builder-builder-ready", html, fixed = TRUE))
  expect_false(grepl("data-builder-id=\"design-builder\"", html, fixed = TRUE))
  expect_false(grepl("id=\"builder-builder-ready\"", html, fixed = TRUE))
  expect_false(grepl("data-builder-id=\"builder\"", html, fixed = TRUE))
  expect_false(grepl("output[&#39;data_available&#39;]", html, fixed = TRUE))
  expect_false(grepl("output[&#39;design_matrix_exists&#39;]", html, fixed = TRUE))
})

test_that("mod_prot_design_ui keeps non-default builder-missing conditional bindings on the shared empty namespace prefix", {
  ui_without_builder <- mod_prot_design_ui
  environment(ui_without_builder) <- list2env(
    list(
      exists = function(x, where = -1, inherits = TRUE) {
        if (identical(x, "mod_prot_design_builder_ui")) {
          FALSE
        } else {
          base::exists(x, where = where, inherits = inherits)
        }
      }
    ),
    parent = environment(mod_prot_design_ui)
  )

  html <- htmltools::renderTags(ui_without_builder("prot-design"))$html

  count_occurrences <- function(pattern) {
    matches <- gregexpr(pattern, html, fixed = TRUE)[[1]]
    if (identical(matches[1], -1L)) {
      return(0L)
    }

    length(matches)
  }

  expect_identical(count_occurrences("Design builder module not loaded"), 1L)
  expect_identical(count_occurrences("data-ns-prefix=\"\""), 3L)
  expect_identical(
    count_occurrences("data-display-if=\"output[&#39;prot-design-data_available&#39;]\""),
    1L
  )
  expect_identical(
    count_occurrences("data-display-if=\"output[&#39;prot-design-design_matrix_exists&#39;]\""),
    1L
  )
  expect_identical(
    count_occurrences("data-display-if=\"!output[&#39;prot-design-data_available&#39;]\""),
    1L
  )
  expect_false(grepl("prot-design-builder-builder-ready", html, fixed = TRUE))
  expect_false(grepl("data-builder-id=\"prot-design-builder\"", html, fixed = TRUE))
  expect_false(grepl("id=\"builder-builder-ready\"", html, fixed = TRUE))
  expect_false(grepl("data-builder-id=\"builder\"", html, fixed = TRUE))
  expect_false(grepl("output[&#39;data_available&#39;]", html, fixed = TRUE))
  expect_false(grepl("output[&#39;design_matrix_exists&#39;]", html, fixed = TRUE))
})

test_that("mod_prot_design_ui keeps the shared progress-percent text hook single-instance and unnamespaced", {
  html <- htmltools::renderTags(mod_prot_design_ui("prot-design"))$html

  count_occurrences <- function(pattern) {
    matches <- gregexpr(pattern, html, fixed = TRUE)[[1]]
    if (identical(matches[1], -1L)) {
      return(0L)
    }

    length(matches)
  }

  expect_identical(
    count_occurrences("$('#uniprot_progress_bar').text(message.percent + '%');"),
    1L
  )
  expect_false(grepl("$('#prot-design-uniprot_progress_bar').text(message.percent + '%');", html, fixed = TRUE))
  expect_false(grepl("prot-design-uniprot_progress_bar", html, fixed = TRUE))
})

test_that("mod_prot_design_ui keeps the shared progress-message text hook single-instance and unnamespaced", {
  html <- htmltools::renderTags(mod_prot_design_ui("prot-design"))$html

  count_occurrences <- function(pattern) {
    matches <- gregexpr(pattern, html, fixed = TRUE)[[1]]
    if (identical(matches[1], -1L)) {
      return(0L)
    }

    length(matches)
  }

  expect_identical(
    count_occurrences("$('#uniprot_progress_text').text(message.text);"),
    1L
  )
  expect_false(grepl("$('#prot-design-uniprot_progress_text').text(message.text);", html, fixed = TRUE))
  expect_false(grepl("prot-design-uniprot_progress_text", html, fixed = TRUE))
})

test_that("mod_prot_design_ui keeps the shared progress-width css hook single-instance and unnamespaced", {
  html <- htmltools::renderTags(mod_prot_design_ui("prot-design"))$html

  count_occurrences <- function(pattern) {
    matches <- gregexpr(pattern, html, fixed = TRUE)[[1]]
    if (identical(matches[1], -1L)) {
      return(0L)
    }

    length(matches)
  }

  expect_identical(
    count_occurrences("$('#uniprot_progress_bar').css('width', message.percent + '%');"),
    1L
  )
  expect_false(grepl("$('#prot-design-uniprot_progress_bar').css('width', message.percent + '%');", html, fixed = TRUE))
  expect_false(grepl("prot-design-uniprot_progress_bar", html, fixed = TRUE))
})

test_that("mod_prot_design_ui keeps the shared progress handler name single-instance and unnamespaced", {
  html <- htmltools::renderTags(mod_prot_design_ui("prot-design"))$html

  count_occurrences <- function(pattern) {
    matches <- gregexpr(pattern, html, fixed = TRUE)[[1]]
    if (identical(matches[1], -1L)) {
      return(0L)
    }

    length(matches)
  }

  expect_identical(
    count_occurrences("Shiny.addCustomMessageHandler('updateUniprotProgress'"),
    1L
  )
  expect_false(grepl("Shiny.addCustomMessageHandler('prot-design-updateUniprotProgress'", html, fixed = TRUE))
  expect_false(grepl("Shiny.addCustomMessageHandler('design-updateUniprotProgress'", html, fixed = TRUE))
  expect_false(grepl("prot-design-updateUniprotProgress", html, fixed = TRUE))
  expect_false(grepl("design-updateUniprotProgress", html, fixed = TRUE))
})

test_that("mod_prot_design_ui keeps the shared UniProt progress script body unnamespaced for the default wrapper id", {
  html <- htmltools::renderTags(mod_prot_design_ui("design"))$html

  count_occurrences <- function(pattern) {
    matches <- gregexpr(pattern, html, fixed = TRUE)[[1]]
    if (identical(matches[1], -1L)) {
      return(0L)
    }

    length(matches)
  }

  expect_identical(
    count_occurrences("Shiny.addCustomMessageHandler('updateUniprotProgress'"),
    1L
  )
  expect_identical(
    count_occurrences("$('#uniprot_progress_bar').css('width', message.percent + '%');"),
    1L
  )
  expect_identical(
    count_occurrences("$('#uniprot_progress_bar').text(message.percent + '%');"),
    1L
  )
  expect_identical(
    count_occurrences("$('#uniprot_progress_text').text(message.text);"),
    1L
  )
  expect_false(grepl("Shiny.addCustomMessageHandler('design-updateUniprotProgress'", html, fixed = TRUE))
  expect_false(grepl("$('#design-uniprot_progress_bar').css('width', message.percent + '%');", html, fixed = TRUE))
  expect_false(grepl("$('#design-uniprot_progress_bar').text(message.percent + '%');", html, fixed = TRUE))
  expect_false(grepl("$('#design-uniprot_progress_text').text(message.text);", html, fixed = TRUE))
  expect_false(grepl("design-updateUniprotProgress", html, fixed = TRUE))
  expect_false(grepl("design-uniprot_progress_bar", html, fixed = TRUE))
  expect_false(grepl("design-uniprot_progress_text", html, fixed = TRUE))
})

test_that("mod_prot_design_ui keeps the shared UniProt progress update order stable for the default wrapper id", {
  html <- htmltools::renderTags(mod_prot_design_ui("design"))$html

  count_occurrences <- function(pattern) {
    matches <- gregexpr(pattern, html, fixed = TRUE)[[1]]
    if (identical(matches[1], -1L)) {
      return(0L)
    }

    length(matches)
  }

  width_update_pos <- as.integer(regexpr(
    "$('#uniprot_progress_bar').css('width', message.percent + '%');",
    html,
    fixed = TRUE
  ))
  percent_text_pos <- as.integer(regexpr(
    "$('#uniprot_progress_bar').text(message.percent + '%');",
    html,
    fixed = TRUE
  ))
  message_text_pos <- as.integer(regexpr(
    "$('#uniprot_progress_text').text(message.text);",
    html,
    fixed = TRUE
  ))

  expect_identical(
    count_occurrences("$('#uniprot_progress_bar').css('width', message.percent + '%');"),
    1L
  )
  expect_identical(
    count_occurrences("$('#uniprot_progress_bar').text(message.percent + '%');"),
    1L
  )
  expect_identical(
    count_occurrences("$('#uniprot_progress_text').text(message.text);"),
    1L
  )
  expect_gt(width_update_pos, 0L)
  expect_gt(percent_text_pos, width_update_pos)
  expect_gt(message_text_pos, percent_text_pos)
  expect_false(grepl("$('#design-uniprot_progress_bar').css('width', message.percent + '%');", html, fixed = TRUE))
  expect_false(grepl("$('#design-uniprot_progress_bar').text(message.percent + '%');", html, fixed = TRUE))
  expect_false(grepl("$('#design-uniprot_progress_text').text(message.text);", html, fixed = TRUE))
})

test_that("mod_prot_design_ui keeps the shared UniProt progress update order stable for the non-default wrapper id", {
  html <- htmltools::renderTags(mod_prot_design_ui("prot-design"))$html

  count_occurrences <- function(pattern) {
    matches <- gregexpr(pattern, html, fixed = TRUE)[[1]]
    if (identical(matches[1], -1L)) {
      return(0L)
    }

    length(matches)
  }

  width_update_pos <- as.integer(regexpr(
    "$('#uniprot_progress_bar').css('width', message.percent + '%');",
    html,
    fixed = TRUE
  ))
  percent_text_pos <- as.integer(regexpr(
    "$('#uniprot_progress_bar').text(message.percent + '%');",
    html,
    fixed = TRUE
  ))
  message_text_pos <- as.integer(regexpr(
    "$('#uniprot_progress_text').text(message.text);",
    html,
    fixed = TRUE
  ))

  expect_identical(
    count_occurrences("$('#uniprot_progress_bar').css('width', message.percent + '%');"),
    1L
  )
  expect_identical(
    count_occurrences("$('#uniprot_progress_bar').text(message.percent + '%');"),
    1L
  )
  expect_identical(
    count_occurrences("$('#uniprot_progress_text').text(message.text);"),
    1L
  )
  expect_gt(width_update_pos, 0L)
  expect_gt(percent_text_pos, width_update_pos)
  expect_gt(message_text_pos, percent_text_pos)
  expect_false(grepl("$('#prot-design-uniprot_progress_bar').css('width', message.percent + '%');", html, fixed = TRUE))
  expect_false(grepl("$('#prot-design-uniprot_progress_bar').text(message.percent + '%');", html, fixed = TRUE))
  expect_false(grepl("$('#prot-design-uniprot_progress_text').text(message.text);", html, fixed = TRUE))
})

test_that("mod_prot_design_ui emits the shared progress script ahead of the default wrapper shell", {
  html <- htmltools::renderTags(mod_prot_design_ui("design"))$html

  script_pos <- as.integer(regexpr(
    "Shiny.addCustomMessageHandler('updateUniprotProgress'",
    html,
    fixed = TRUE
  ))
  import_button_pos <- as.integer(regexpr(
    "id=\"design-show_import_modal\"",
    html,
    fixed = TRUE
  ))
  builder_panel_pos <- as.integer(regexpr(
    "data-display-if=\"output[&#39;design-data_available&#39;]\"",
    html,
    fixed = TRUE
  ))
  empty_state_pos <- as.integer(regexpr(
    "data-display-if=\"!output[&#39;design-data_available&#39;]\"",
    html,
    fixed = TRUE
  ))

  expect_gt(script_pos, 0L)
  expect_gt(import_button_pos, script_pos)
  expect_gt(builder_panel_pos, import_button_pos)
  expect_gt(empty_state_pos, builder_panel_pos)
  expect_false(grepl(
    "Shiny.addCustomMessageHandler('design-updateUniprotProgress'",
    html,
    fixed = TRUE
  ))
})

test_that("mod_prot_design_ui emits the shared progress script ahead of the non-default wrapper shell", {
  html <- htmltools::renderTags(mod_prot_design_ui("prot-design"))$html

  script_pos <- as.integer(regexpr(
    "Shiny.addCustomMessageHandler('updateUniprotProgress'",
    html,
    fixed = TRUE
  ))
  import_button_pos <- as.integer(regexpr(
    "id=\"prot-design-show_import_modal\"",
    html,
    fixed = TRUE
  ))
  builder_panel_pos <- as.integer(regexpr(
    "data-display-if=\"output[&#39;prot-design-data_available&#39;]\"",
    html,
    fixed = TRUE
  ))
  empty_state_pos <- as.integer(regexpr(
    "data-display-if=\"!output[&#39;prot-design-data_available&#39;]\"",
    html,
    fixed = TRUE
  ))

  expect_gt(script_pos, 0L)
  expect_gt(import_button_pos, script_pos)
  expect_gt(builder_panel_pos, import_button_pos)
  expect_gt(empty_state_pos, builder_panel_pos)
  expect_false(grepl(
    "Shiny.addCustomMessageHandler('prot-design-updateUniprotProgress'",
    html,
    fixed = TRUE
  ))
  expect_false(grepl("id=\"show_import_modal\"", html, fixed = TRUE))
  expect_false(grepl("output[&#39;data_available&#39;]", html, fixed = TRUE))
})

test_that("mod_prot_design_ui emits the shared progress script ahead of the default fallback shell", {
  ui_without_builder <- mod_prot_design_ui
  environment(ui_without_builder) <- list2env(
    list(
      exists = function(x, where = -1, inherits = TRUE) {
        if (identical(x, "mod_prot_design_builder_ui")) {
          FALSE
        } else {
          base::exists(x, where = where, inherits = inherits)
        }
      }
    ),
    parent = environment(mod_prot_design_ui)
  )

  html <- htmltools::renderTags(ui_without_builder("design"))$html

  script_pos <- as.integer(regexpr(
    "Shiny.addCustomMessageHandler('updateUniprotProgress'",
    html,
    fixed = TRUE
  ))
  import_button_pos <- as.integer(regexpr(
    "id=\"design-show_import_modal\"",
    html,
    fixed = TRUE
  ))
  fallback_pos <- as.integer(regexpr(
    "Design builder module not loaded",
    html,
    fixed = TRUE
  ))
  preview_heading_pos <- as.integer(regexpr(
    "Saved Results Preview",
    html,
    fixed = TRUE
  ))
  empty_state_pos <- as.integer(regexpr(
    "data-display-if=\"!output[&#39;design-data_available&#39;]\"",
    html,
    fixed = TRUE
  ))

  expect_gt(script_pos, 0L)
  expect_gt(import_button_pos, script_pos)
  expect_gt(fallback_pos, import_button_pos)
  expect_gt(preview_heading_pos, fallback_pos)
  expect_gt(empty_state_pos, preview_heading_pos)
  expect_false(grepl(
    "Shiny.addCustomMessageHandler('design-updateUniprotProgress'",
    html,
    fixed = TRUE
  ))
  expect_false(grepl("design-builder-builder-ready", html, fixed = TRUE))
})

test_that("mod_prot_design_ui emits the shared progress script ahead of the non-default fallback shell", {
  ui_without_builder <- mod_prot_design_ui
  environment(ui_without_builder) <- list2env(
    list(
      exists = function(x, where = -1, inherits = TRUE) {
        if (identical(x, "mod_prot_design_builder_ui")) {
          FALSE
        } else {
          base::exists(x, where = where, inherits = inherits)
        }
      }
    ),
    parent = environment(mod_prot_design_ui)
  )

  html <- htmltools::renderTags(ui_without_builder("prot-design"))$html

  script_pos <- as.integer(regexpr(
    "Shiny.addCustomMessageHandler('updateUniprotProgress'",
    html,
    fixed = TRUE
  ))
  import_button_pos <- as.integer(regexpr(
    "id=\"prot-design-show_import_modal\"",
    html,
    fixed = TRUE
  ))
  fallback_pos <- as.integer(regexpr(
    "Design builder module not loaded",
    html,
    fixed = TRUE
  ))
  preview_heading_pos <- as.integer(regexpr(
    "Saved Results Preview",
    html,
    fixed = TRUE
  ))
  empty_state_pos <- as.integer(regexpr(
    "data-display-if=\"!output[&#39;prot-design-data_available&#39;]\"",
    html,
    fixed = TRUE
  ))

  expect_gt(script_pos, 0L)
  expect_gt(import_button_pos, script_pos)
  expect_gt(fallback_pos, import_button_pos)
  expect_gt(preview_heading_pos, fallback_pos)
  expect_gt(empty_state_pos, preview_heading_pos)
  expect_false(grepl(
    "Shiny.addCustomMessageHandler('prot-design-updateUniprotProgress'",
    html,
    fixed = TRUE
  ))
  expect_false(grepl("id=\"show_import_modal\"", html, fixed = TRUE))
  expect_false(grepl("output[&#39;data_available&#39;]", html, fixed = TRUE))
  expect_false(grepl("prot-design-builder-builder-ready", html, fixed = TRUE))
})

test_that("mod_prot_design_ui emits the shared progress script ahead of the default builder shell", {
  ui_with_builder <- mod_prot_design_ui
  environment(ui_with_builder) <- list2env(
    list(
      exists = function(x, where = -1, inherits = TRUE) {
        if (identical(x, "mod_prot_design_builder_ui")) {
          TRUE
        } else {
          base::exists(x, where = where, inherits = inherits)
        }
      },
      mod_prot_design_builder_ui = function(id) {
        shiny::div(
          id = paste0(id, "-builder-ready"),
          `data-builder-id` = id,
          "Builder ready"
        )
      }
    ),
    parent = environment(mod_prot_design_ui)
  )

  html <- htmltools::renderTags(ui_with_builder("design"))$html

  script_pos <- as.integer(regexpr(
    "Shiny.addCustomMessageHandler('updateUniprotProgress'",
    html,
    fixed = TRUE
  ))
  import_button_pos <- as.integer(regexpr(
    "id=\"design-show_import_modal\"",
    html,
    fixed = TRUE
  ))
  builder_pos <- as.integer(regexpr(
    "id=\"design-builder-builder-ready\"",
    html,
    fixed = TRUE
  ))
  preview_heading_pos <- as.integer(regexpr(
    "Saved Results Preview",
    html,
    fixed = TRUE
  ))
  empty_state_pos <- as.integer(regexpr(
    "data-display-if=\"!output[&#39;design-data_available&#39;]\"",
    html,
    fixed = TRUE
  ))

  expect_gt(script_pos, 0L)
  expect_gt(import_button_pos, script_pos)
  expect_gt(builder_pos, import_button_pos)
  expect_gt(preview_heading_pos, builder_pos)
  expect_gt(empty_state_pos, preview_heading_pos)
  expect_false(grepl(
    "Shiny.addCustomMessageHandler('design-updateUniprotProgress'",
    html,
    fixed = TRUE
  ))
  expect_false(grepl("Design builder module not loaded", html, fixed = TRUE))
  expect_false(grepl("data-builder-id=\"builder\"", html, fixed = TRUE))
})

test_that("mod_prot_design_ui keeps the default builder guidance paragraph ordered ahead of the embedded builder shell", {
  ui_with_builder <- mod_prot_design_ui
  environment(ui_with_builder) <- list2env(
    list(
      exists = function(x, where = -1, inherits = TRUE) {
        if (identical(x, "mod_prot_design_builder_ui")) {
          TRUE
        } else {
          base::exists(x, where = where, inherits = inherits)
        }
      },
      mod_prot_design_builder_ui = function(id) {
        shiny::div(
          id = paste0(id, "-builder-ready"),
          `data-builder-id` = id,
          "Builder ready"
        )
      }
    ),
    parent = environment(mod_prot_design_ui)
  )

  html <- htmltools::renderTags(ui_with_builder("design"))$html

  count_occurrences <- function(pattern) {
    matches <- gregexpr(pattern, html, fixed = TRUE)[[1]]
    if (identical(matches[1], -1L)) {
      return(0L)
    }

    length(matches)
  }

  import_button_pos <- as.integer(regexpr(
    "id=\"design-show_import_modal\"",
    html,
    fixed = TRUE
  ))
  guidance_pos <- as.integer(regexpr(
    "Use the tools below to define your experimental groups and the contrasts for differential analysis. Alternatively, import an existing design from a previous analysis.",
    html,
    fixed = TRUE
  ))
  builder_pos <- as.integer(regexpr(
    "id=\"design-builder-builder-ready\"",
    html,
    fixed = TRUE
  ))
  preview_heading_pos <- as.integer(regexpr(
    "Saved Results Preview",
    html,
    fixed = TRUE
  ))

  expect_identical(
    count_occurrences("Use the tools below to define your experimental groups and the contrasts for differential analysis. Alternatively, import an existing design from a previous analysis."),
    1L
  )
  expect_gt(import_button_pos, 0L)
  expect_gt(guidance_pos, import_button_pos)
  expect_gt(builder_pos, guidance_pos)
  expect_gt(preview_heading_pos, builder_pos)
  expect_false(grepl("Design builder module not loaded", html, fixed = TRUE))
  expect_false(grepl("data-builder-id=\"builder\"", html, fixed = TRUE))
})

test_that("mod_prot_design_ui keeps the non-default builder guidance paragraph ordered ahead of the embedded builder shell", {
  ui_with_builder <- mod_prot_design_ui
  environment(ui_with_builder) <- list2env(
    list(
      exists = function(x, where = -1, inherits = TRUE) {
        if (identical(x, "mod_prot_design_builder_ui")) {
          TRUE
        } else {
          base::exists(x, where = where, inherits = inherits)
        }
      },
      mod_prot_design_builder_ui = function(id) {
        shiny::div(
          id = paste0(id, "-builder-ready"),
          `data-builder-id` = id,
          "Builder ready"
        )
      }
    ),
    parent = environment(mod_prot_design_ui)
  )

  html <- htmltools::renderTags(ui_with_builder("prot-design"))$html

  count_occurrences <- function(pattern) {
    matches <- gregexpr(pattern, html, fixed = TRUE)[[1]]
    if (identical(matches[1], -1L)) {
      return(0L)
    }

    length(matches)
  }

  import_button_pos <- as.integer(regexpr(
    "id=\"prot-design-show_import_modal\"",
    html,
    fixed = TRUE
  ))
  guidance_pos <- as.integer(regexpr(
    "Use the tools below to define your experimental groups and the contrasts for differential analysis. Alternatively, import an existing design from a previous analysis.",
    html,
    fixed = TRUE
  ))
  builder_pos <- as.integer(regexpr(
    "id=\"prot-design-builder-builder-ready\"",
    html,
    fixed = TRUE
  ))
  preview_heading_pos <- as.integer(regexpr(
    "Saved Results Preview",
    html,
    fixed = TRUE
  ))

  expect_identical(
    count_occurrences("Use the tools below to define your experimental groups and the contrasts for differential analysis. Alternatively, import an existing design from a previous analysis."),
    1L
  )
  expect_gt(import_button_pos, 0L)
  expect_gt(guidance_pos, import_button_pos)
  expect_gt(builder_pos, guidance_pos)
  expect_gt(preview_heading_pos, builder_pos)
  expect_false(grepl("Design builder module not loaded", html, fixed = TRUE))
  expect_false(grepl("id=\"builder-builder-ready\"", html, fixed = TRUE))
  expect_false(grepl("data-builder-id=\"builder\"", html, fixed = TRUE))
})

test_that("mod_prot_design_ui keeps the non-default fallback builder guidance paragraph ordered ahead of the fallback shell", {
  ui_without_builder <- mod_prot_design_ui
  environment(ui_without_builder) <- list2env(
    list(
      exists = function(x, where = -1, inherits = TRUE) {
        if (identical(x, "mod_prot_design_builder_ui")) {
          FALSE
        } else {
          base::exists(x, where = where, inherits = inherits)
        }
      }
    ),
    parent = environment(mod_prot_design_ui)
  )

  html <- htmltools::renderTags(ui_without_builder("prot-design"))$html

  count_occurrences <- function(pattern) {
    matches <- gregexpr(pattern, html, fixed = TRUE)[[1]]
    if (identical(matches[1], -1L)) {
      return(0L)
    }

    length(matches)
  }

  import_button_pos <- as.integer(regexpr(
    "id=\"prot-design-show_import_modal\"",
    html,
    fixed = TRUE
  ))
  guidance_pos <- as.integer(regexpr(
    "Use the tools below to define your experimental groups and the contrasts for differential analysis. Alternatively, import an existing design from a previous analysis.",
    html,
    fixed = TRUE
  ))
  fallback_pos <- as.integer(regexpr(
    "Design builder module not loaded",
    html,
    fixed = TRUE
  ))
  preview_heading_pos <- as.integer(regexpr(
    "Saved Results Preview",
    html,
    fixed = TRUE
  ))

  expect_identical(
    count_occurrences("Use the tools below to define your experimental groups and the contrasts for differential analysis. Alternatively, import an existing design from a previous analysis."),
    1L
  )
  expect_identical(count_occurrences("Design builder module not loaded"), 1L)
  expect_gt(import_button_pos, 0L)
  expect_gt(guidance_pos, import_button_pos)
  expect_gt(fallback_pos, guidance_pos)
  expect_gt(preview_heading_pos, fallback_pos)
  expect_false(grepl("id=\"show_import_modal\"", html, fixed = TRUE))
  expect_false(grepl("prot-design-builder-builder-ready", html, fixed = TRUE))
  expect_false(grepl("data-builder-id=\"prot-design-builder\"", html, fixed = TRUE))
})

test_that("mod_prot_design_ui keeps the default fallback builder guidance paragraph ordered ahead of the fallback shell", {
  ui_without_builder <- mod_prot_design_ui
  environment(ui_without_builder) <- list2env(
    list(
      exists = function(x, where = -1, inherits = TRUE) {
        if (identical(x, "mod_prot_design_builder_ui")) {
          FALSE
        } else {
          base::exists(x, where = where, inherits = inherits)
        }
      }
    ),
    parent = environment(mod_prot_design_ui)
  )

  html <- htmltools::renderTags(ui_without_builder("design"))$html

  count_occurrences <- function(pattern) {
    matches <- gregexpr(pattern, html, fixed = TRUE)[[1]]
    if (identical(matches[1], -1L)) {
      return(0L)
    }

    length(matches)
  }

  import_button_pos <- as.integer(regexpr(
    "id=\"design-show_import_modal\"",
    html,
    fixed = TRUE
  ))
  guidance_pos <- as.integer(regexpr(
    "Use the tools below to define your experimental groups and the contrasts for differential analysis. Alternatively, import an existing design from a previous analysis.",
    html,
    fixed = TRUE
  ))
  fallback_pos <- as.integer(regexpr(
    "Design builder module not loaded",
    html,
    fixed = TRUE
  ))
  preview_heading_pos <- as.integer(regexpr(
    "Saved Results Preview",
    html,
    fixed = TRUE
  ))

  expect_identical(
    count_occurrences("Use the tools below to define your experimental groups and the contrasts for differential analysis. Alternatively, import an existing design from a previous analysis."),
    1L
  )
  expect_identical(count_occurrences("Design builder module not loaded"), 1L)
  expect_gt(import_button_pos, 0L)
  expect_gt(guidance_pos, import_button_pos)
  expect_gt(fallback_pos, guidance_pos)
  expect_gt(preview_heading_pos, fallback_pos)
  expect_false(grepl("id=\"show_import_modal\"", html, fixed = TRUE))
  expect_false(grepl("design-builder-builder-ready", html, fixed = TRUE))
  expect_false(grepl("data-builder-id=\"design-builder\"", html, fixed = TRUE))
})

test_that("mod_prot_design_ui keeps the default preview guidance paragraph ordered ahead of the saved previews", {
  ui_with_builder <- mod_prot_design_ui
  environment(ui_with_builder) <- list2env(
    list(
      exists = function(x, where = -1, inherits = TRUE) {
        if (identical(x, "mod_prot_design_builder_ui")) {
          TRUE
        } else {
          base::exists(x, where = where, inherits = inherits)
        }
      },
      mod_prot_design_builder_ui = function(id) {
        shiny::div(
          id = paste0(id, "-builder-ready"),
          `data-builder-id` = id,
          "Builder ready"
        )
      }
    ),
    parent = environment(mod_prot_design_ui)
  )

  html <- htmltools::renderTags(ui_with_builder("design"))$html

  count_occurrences <- function(pattern) {
    matches <- gregexpr(pattern, html, fixed = TRUE)[[1]]
    if (identical(matches[1], -1L)) {
      return(0L)
    }

    length(matches)
  }

  preview_heading_pos <- as.integer(regexpr(
    "Saved Results Preview",
    html,
    fixed = TRUE
  ))
  preview_guidance_pos <- as.integer(regexpr(
    "This section shows the design matrix and contrasts that have been saved to the workflow.",
    html,
    fixed = TRUE
  ))
  design_preview_pos <- as.integer(regexpr(
    "id=\"design-design_matrix_preview\"",
    html,
    fixed = TRUE
  ))
  contrasts_preview_pos <- as.integer(regexpr(
    "id=\"design-contrasts_preview\"",
    html,
    fixed = TRUE
  ))

  expect_identical(
    count_occurrences("This section shows the design matrix and contrasts that have been saved to the workflow."),
    1L
  )
  expect_gt(preview_heading_pos, 0L)
  expect_gt(preview_guidance_pos, preview_heading_pos)
  expect_gt(design_preview_pos, preview_guidance_pos)
  expect_gt(contrasts_preview_pos, design_preview_pos)
  expect_false(grepl("Design builder module not loaded", html, fixed = TRUE))
  expect_false(grepl("id=\"builder-builder-ready\"", html, fixed = TRUE))
  expect_false(grepl("data-builder-id=\"builder\"", html, fixed = TRUE))
})

test_that("mod_prot_design_ui keeps the non-default preview guidance paragraph ordered ahead of the saved previews", {
  ui_with_builder <- mod_prot_design_ui
  environment(ui_with_builder) <- list2env(
    list(
      exists = function(x, where = -1, inherits = TRUE) {
        if (identical(x, "mod_prot_design_builder_ui")) {
          TRUE
        } else {
          base::exists(x, where = where, inherits = inherits)
        }
      },
      mod_prot_design_builder_ui = function(id) {
        shiny::div(
          id = paste0(id, "-builder-ready"),
          `data-builder-id` = id,
          "Builder ready"
        )
      }
    ),
    parent = environment(mod_prot_design_ui)
  )

  html <- htmltools::renderTags(ui_with_builder("prot-design"))$html

  count_occurrences <- function(pattern) {
    matches <- gregexpr(pattern, html, fixed = TRUE)[[1]]
    if (identical(matches[1], -1L)) {
      return(0L)
    }

    length(matches)
  }

  preview_heading_pos <- as.integer(regexpr(
    "Saved Results Preview",
    html,
    fixed = TRUE
  ))
  preview_guidance_pos <- as.integer(regexpr(
    "This section shows the design matrix and contrasts that have been saved to the workflow.",
    html,
    fixed = TRUE
  ))
  design_preview_pos <- as.integer(regexpr(
    "id=\"prot-design-design_matrix_preview\"",
    html,
    fixed = TRUE
  ))
  contrasts_preview_pos <- as.integer(regexpr(
    "id=\"prot-design-contrasts_preview\"",
    html,
    fixed = TRUE
  ))

  expect_identical(
    count_occurrences("This section shows the design matrix and contrasts that have been saved to the workflow."),
    1L
  )
  expect_gt(preview_heading_pos, 0L)
  expect_gt(preview_guidance_pos, preview_heading_pos)
  expect_gt(design_preview_pos, preview_guidance_pos)
  expect_gt(contrasts_preview_pos, design_preview_pos)
  expect_false(grepl("Design builder module not loaded", html, fixed = TRUE))
  expect_false(grepl("id=\"builder-builder-ready\"", html, fixed = TRUE))
  expect_false(grepl("data-builder-id=\"builder\"", html, fixed = TRUE))
})

test_that("mod_prot_design_ui keeps the non-default fallback preview guidance paragraph ordered ahead of the saved previews", {
  ui_without_builder <- mod_prot_design_ui
  environment(ui_without_builder) <- list2env(
    list(
      exists = function(x, where = -1, inherits = TRUE) {
        if (identical(x, "mod_prot_design_builder_ui")) {
          FALSE
        } else {
          base::exists(x, where = where, inherits = inherits)
        }
      }
    ),
    parent = environment(mod_prot_design_ui)
  )

  html <- htmltools::renderTags(ui_without_builder("prot-design"))$html

  count_occurrences <- function(pattern) {
    matches <- gregexpr(pattern, html, fixed = TRUE)[[1]]
    if (identical(matches[1], -1L)) {
      return(0L)
    }

    length(matches)
  }

  fallback_pos <- as.integer(regexpr(
    "Design builder module not loaded",
    html,
    fixed = TRUE
  ))
  preview_heading_pos <- as.integer(regexpr(
    "Saved Results Preview",
    html,
    fixed = TRUE
  ))
  preview_guidance_pos <- as.integer(regexpr(
    "This section shows the design matrix and contrasts that have been saved to the workflow.",
    html,
    fixed = TRUE
  ))
  design_preview_pos <- as.integer(regexpr(
    "id=\"prot-design-design_matrix_preview\"",
    html,
    fixed = TRUE
  ))
  contrasts_preview_pos <- as.integer(regexpr(
    "id=\"prot-design-contrasts_preview\"",
    html,
    fixed = TRUE
  ))

  expect_identical(count_occurrences("Design builder module not loaded"), 1L)
  expect_identical(
    count_occurrences("This section shows the design matrix and contrasts that have been saved to the workflow."),
    1L
  )
  expect_gt(fallback_pos, 0L)
  expect_gt(preview_heading_pos, fallback_pos)
  expect_gt(preview_guidance_pos, preview_heading_pos)
  expect_gt(design_preview_pos, preview_guidance_pos)
  expect_gt(contrasts_preview_pos, design_preview_pos)
  expect_false(grepl("prot-design-builder-builder-ready", html, fixed = TRUE))
  expect_false(grepl("data-builder-id=\"prot-design-builder\"", html, fixed = TRUE))
  expect_false(grepl("id=\"builder-builder-ready\"", html, fixed = TRUE))
  expect_false(grepl("data-builder-id=\"builder\"", html, fixed = TRUE))
})

test_that("mod_prot_design_ui keeps the default fallback preview guidance paragraph ordered ahead of the saved previews", {
  ui_without_builder <- mod_prot_design_ui
  environment(ui_without_builder) <- list2env(
    list(
      exists = function(x, where = -1, inherits = TRUE) {
        if (identical(x, "mod_prot_design_builder_ui")) {
          FALSE
        } else {
          base::exists(x, where = where, inherits = inherits)
        }
      }
    ),
    parent = environment(mod_prot_design_ui)
  )

  html <- htmltools::renderTags(ui_without_builder("design"))$html

  count_occurrences <- function(pattern) {
    matches <- gregexpr(pattern, html, fixed = TRUE)[[1]]
    if (identical(matches[1], -1L)) {
      return(0L)
    }

    length(matches)
  }

  fallback_pos <- as.integer(regexpr(
    "Design builder module not loaded",
    html,
    fixed = TRUE
  ))
  preview_heading_pos <- as.integer(regexpr(
    "Saved Results Preview",
    html,
    fixed = TRUE
  ))
  preview_guidance_pos <- as.integer(regexpr(
    "This section shows the design matrix and contrasts that have been saved to the workflow.",
    html,
    fixed = TRUE
  ))
  design_preview_pos <- as.integer(regexpr(
    "id=\"design-design_matrix_preview\"",
    html,
    fixed = TRUE
  ))
  contrasts_preview_pos <- as.integer(regexpr(
    "id=\"design-contrasts_preview\"",
    html,
    fixed = TRUE
  ))

  expect_identical(count_occurrences("Design builder module not loaded"), 1L)
  expect_identical(
    count_occurrences("This section shows the design matrix and contrasts that have been saved to the workflow."),
    1L
  )
  expect_gt(fallback_pos, 0L)
  expect_gt(preview_heading_pos, fallback_pos)
  expect_gt(preview_guidance_pos, preview_heading_pos)
  expect_gt(design_preview_pos, preview_guidance_pos)
  expect_gt(contrasts_preview_pos, design_preview_pos)
  expect_false(grepl("id=\"design_matrix_preview\"", html, fixed = TRUE))
  expect_false(grepl("id=\"contrasts_preview\"", html, fixed = TRUE))
  expect_false(grepl("design-builder-builder-ready", html, fixed = TRUE))
  expect_false(grepl("data-builder-id=\"design-builder\"", html, fixed = TRUE))
  expect_false(grepl("id=\"builder-builder-ready\"", html, fixed = TRUE))
  expect_false(grepl("data-builder-id=\"builder\"", html, fixed = TRUE))
})

test_that("mod_prot_design_ui keeps the default wrapper heading single-instance and ordered ahead of the fallback body", {
  ui_without_builder <- mod_prot_design_ui
  environment(ui_without_builder) <- list2env(
    list(
      exists = function(x, where = -1, inherits = TRUE) {
        if (identical(x, "mod_prot_design_builder_ui")) {
          FALSE
        } else {
          base::exists(x, where = where, inherits = inherits)
        }
      }
    ),
    parent = environment(mod_prot_design_ui)
  )

  html <- htmltools::renderTags(ui_without_builder("design"))$html

  count_occurrences <- function(pattern) {
    matches <- gregexpr(pattern, html, fixed = TRUE)[[1]]
    if (identical(matches[1], -1L)) {
      return(0L)
    }

    length(matches)
  }

  heading_pos <- as.integer(regexpr(
    "Design Matrix Builder",
    html,
    fixed = TRUE
  ))
  import_button_pos <- as.integer(regexpr(
    "id=\"design-show_import_modal\"",
    html,
    fixed = TRUE
  ))
  guidance_pos <- as.integer(regexpr(
    "Use the tools below to define your experimental groups and the contrasts for differential analysis. Alternatively, import an existing design from a previous analysis.",
    html,
    fixed = TRUE
  ))
  fallback_pos <- as.integer(regexpr(
    "Design builder module not loaded",
    html,
    fixed = TRUE
  ))

  expect_identical(count_occurrences("Design Matrix Builder"), 1L)
  expect_gt(heading_pos, 0L)
  expect_gt(import_button_pos, heading_pos)
  expect_gt(guidance_pos, import_button_pos)
  expect_gt(fallback_pos, guidance_pos)
  expect_false(grepl("id=\"show_import_modal\"", html, fixed = TRUE))
  expect_false(grepl("design-builder-builder-ready", html, fixed = TRUE))
})

test_that("mod_prot_design_ui keeps the default wrapper heading single-instance and ordered ahead of the builder body", {
  ui_with_builder <- mod_prot_design_ui
  environment(ui_with_builder) <- list2env(
    list(
      exists = function(x, where = -1, inherits = TRUE) {
        if (identical(x, "mod_prot_design_builder_ui")) {
          TRUE
        } else {
          base::exists(x, where = where, inherits = inherits)
        }
      },
      mod_prot_design_builder_ui = function(id) {
        shiny::div(
          id = paste0(id, "-builder-ready"),
          `data-builder-id` = id,
          "Builder ready"
        )
      }
    ),
    parent = environment(mod_prot_design_ui)
  )

  html <- htmltools::renderTags(ui_with_builder("design"))$html

  count_occurrences <- function(pattern) {
    matches <- gregexpr(pattern, html, fixed = TRUE)[[1]]
    if (identical(matches[1], -1L)) {
      return(0L)
    }

    length(matches)
  }

  heading_pos <- as.integer(regexpr(
    "Design Matrix Builder",
    html,
    fixed = TRUE
  ))
  import_button_pos <- as.integer(regexpr(
    "id=\"design-show_import_modal\"",
    html,
    fixed = TRUE
  ))
  guidance_pos <- as.integer(regexpr(
    "Use the tools below to define your experimental groups and the contrasts for differential analysis. Alternatively, import an existing design from a previous analysis.",
    html,
    fixed = TRUE
  ))
  builder_pos <- as.integer(regexpr(
    "id=\"design-builder-builder-ready\"",
    html,
    fixed = TRUE
  ))
  preview_heading_pos <- as.integer(regexpr(
    "Saved Results Preview",
    html,
    fixed = TRUE
  ))

  expect_identical(count_occurrences("Design Matrix Builder"), 1L)
  expect_gt(heading_pos, 0L)
  expect_gt(import_button_pos, heading_pos)
  expect_gt(guidance_pos, import_button_pos)
  expect_gt(builder_pos, guidance_pos)
  expect_gt(preview_heading_pos, builder_pos)
  expect_false(grepl("Design builder module not loaded", html, fixed = TRUE))
  expect_false(grepl("id=\"builder-builder-ready\"", html, fixed = TRUE))
  expect_false(grepl("data-builder-id=\"builder\"", html, fixed = TRUE))
})

test_that("mod_prot_design_ui keeps the non-default wrapper heading single-instance and ordered ahead of the builder body", {
  ui_with_builder <- mod_prot_design_ui
  environment(ui_with_builder) <- list2env(
    list(
      exists = function(x, where = -1, inherits = TRUE) {
        if (identical(x, "mod_prot_design_builder_ui")) {
          TRUE
        } else {
          base::exists(x, where = where, inherits = inherits)
        }
      },
      mod_prot_design_builder_ui = function(id) {
        shiny::div(
          id = paste0(id, "-builder-ready"),
          `data-builder-id` = id,
          "Builder ready"
        )
      }
    ),
    parent = environment(mod_prot_design_ui)
  )

  html <- htmltools::renderTags(ui_with_builder("prot-design"))$html

  count_occurrences <- function(pattern) {
    matches <- gregexpr(pattern, html, fixed = TRUE)[[1]]
    if (identical(matches[1], -1L)) {
      return(0L)
    }

    length(matches)
  }

  heading_pos <- as.integer(regexpr(
    "Design Matrix Builder",
    html,
    fixed = TRUE
  ))
  import_button_pos <- as.integer(regexpr(
    "id=\"prot-design-show_import_modal\"",
    html,
    fixed = TRUE
  ))
  guidance_pos <- as.integer(regexpr(
    "Use the tools below to define your experimental groups and the contrasts for differential analysis. Alternatively, import an existing design from a previous analysis.",
    html,
    fixed = TRUE
  ))
  builder_pos <- as.integer(regexpr(
    "id=\"prot-design-builder-builder-ready\"",
    html,
    fixed = TRUE
  ))
  preview_heading_pos <- as.integer(regexpr(
    "Saved Results Preview",
    html,
    fixed = TRUE
  ))

  expect_identical(count_occurrences("Design Matrix Builder"), 1L)
  expect_gt(heading_pos, 0L)
  expect_gt(import_button_pos, heading_pos)
  expect_gt(guidance_pos, import_button_pos)
  expect_gt(builder_pos, guidance_pos)
  expect_gt(preview_heading_pos, builder_pos)
  expect_false(grepl("Design builder module not loaded", html, fixed = TRUE))
  expect_false(grepl("id=\"builder-builder-ready\"", html, fixed = TRUE))
  expect_false(grepl("data-builder-id=\"builder\"", html, fixed = TRUE))
})

test_that("mod_prot_design_ui keeps the non-default wrapper heading single-instance and ordered ahead of the fallback body", {
  ui_without_builder <- mod_prot_design_ui
  environment(ui_without_builder) <- list2env(
    list(
      exists = function(x, where = -1, inherits = TRUE) {
        if (identical(x, "mod_prot_design_builder_ui")) {
          FALSE
        } else {
          base::exists(x, where = where, inherits = inherits)
        }
      }
    ),
    parent = environment(mod_prot_design_ui)
  )

  html <- htmltools::renderTags(ui_without_builder("prot-design"))$html

  count_occurrences <- function(pattern) {
    matches <- gregexpr(pattern, html, fixed = TRUE)[[1]]
    if (identical(matches[1], -1L)) {
      return(0L)
    }

    length(matches)
  }

  heading_pos <- as.integer(regexpr(
    "Design Matrix Builder",
    html,
    fixed = TRUE
  ))
  import_button_pos <- as.integer(regexpr(
    "id=\"prot-design-show_import_modal\"",
    html,
    fixed = TRUE
  ))
  guidance_pos <- as.integer(regexpr(
    "Use the tools below to define your experimental groups and the contrasts for differential analysis. Alternatively, import an existing design from a previous analysis.",
    html,
    fixed = TRUE
  ))
  fallback_pos <- as.integer(regexpr(
    "Design builder module not loaded",
    html,
    fixed = TRUE
  ))

  expect_identical(count_occurrences("Design Matrix Builder"), 1L)
  expect_gt(heading_pos, 0L)
  expect_gt(import_button_pos, heading_pos)
  expect_gt(guidance_pos, import_button_pos)
  expect_gt(fallback_pos, guidance_pos)
  expect_false(grepl("id=\"show_import_modal\"", html, fixed = TRUE))
  expect_false(grepl("id=\"design-show_import_modal\"", html, fixed = TRUE))
  expect_false(grepl("builder-builder-ready", html, fixed = TRUE))
})

test_that("mod_prot_design_ui keeps the non-default wrapper heading single-instance and ordered ahead of the saved-preview guidance", {
  ui_with_builder <- mod_prot_design_ui
  environment(ui_with_builder) <- list2env(
    list(
      exists = function(x, where = -1, inherits = TRUE) {
        if (identical(x, "mod_prot_design_builder_ui")) {
          TRUE
        } else {
          base::exists(x, where = where, inherits = inherits)
        }
      },
      mod_prot_design_builder_ui = function(id) {
        shiny::div(
          id = paste0(id, "-builder-ready"),
          `data-builder-id` = id,
          "Builder ready"
        )
      }
    ),
    parent = environment(mod_prot_design_ui)
  )

  html <- htmltools::renderTags(ui_with_builder("prot-design"))$html

  count_occurrences <- function(pattern) {
    matches <- gregexpr(pattern, html, fixed = TRUE)[[1]]
    if (identical(matches[1], -1L)) {
      return(0L)
    }

    length(matches)
  }

  heading_pos <- as.integer(regexpr(
    "Design Matrix Builder",
    html,
    fixed = TRUE
  ))
  import_button_pos <- as.integer(regexpr(
    "id=\"prot-design-show_import_modal\"",
    html,
    fixed = TRUE
  ))
  top_guidance_pos <- as.integer(regexpr(
    "Use the tools below to define your experimental groups and the contrasts for differential analysis. Alternatively, import an existing design from a previous analysis.",
    html,
    fixed = TRUE
  ))
  builder_pos <- as.integer(regexpr(
    "id=\"prot-design-builder-builder-ready\"",
    html,
    fixed = TRUE
  ))
  preview_heading_pos <- as.integer(regexpr(
    "Saved Results Preview",
    html,
    fixed = TRUE
  ))
  preview_guidance_pos <- as.integer(regexpr(
    "This section shows the design matrix and contrasts that have been saved to the workflow.",
    html,
    fixed = TRUE
  ))

  expect_identical(count_occurrences("Design Matrix Builder"), 1L)
  expect_gt(heading_pos, 0L)
  expect_gt(import_button_pos, heading_pos)
  expect_gt(top_guidance_pos, import_button_pos)
  expect_gt(builder_pos, top_guidance_pos)
  expect_gt(preview_heading_pos, builder_pos)
  expect_gt(preview_guidance_pos, preview_heading_pos)
  expect_false(grepl("Design builder module not loaded", html, fixed = TRUE))
  expect_false(grepl("id=\"builder-builder-ready\"", html, fixed = TRUE))
  expect_false(grepl("data-builder-id=\"builder\"", html, fixed = TRUE))
})

test_that("mod_prot_design_ui keeps the non-default saved-results preview headings ordered ahead of the preview tables", {
  ui_with_builder <- mod_prot_design_ui
  environment(ui_with_builder) <- list2env(
    list(
      exists = function(x, where = -1, inherits = TRUE) {
        if (identical(x, "mod_prot_design_builder_ui")) {
          TRUE
        } else {
          base::exists(x, where = where, inherits = inherits)
        }
      },
      mod_prot_design_builder_ui = function(id) {
        shiny::div(
          id = paste0(id, "-builder-ready"),
          `data-builder-id` = id,
          "Builder ready"
        )
      }
    ),
    parent = environment(mod_prot_design_ui)
  )

  html <- htmltools::renderTags(ui_with_builder("prot-design"))$html

  count_occurrences <- function(pattern) {
    matches <- gregexpr(pattern, html, fixed = TRUE)[[1]]
    if (identical(matches[1], -1L)) {
      return(0L)
    }

    length(matches)
  }

  preview_heading_pos <- as.integer(regexpr(
    "Saved Results Preview",
    html,
    fixed = TRUE
  ))
  current_heading_pos <- as.integer(regexpr(
    "Current Design Matrix",
    html,
    fixed = TRUE
  ))
  design_preview_pos <- as.integer(regexpr(
    "id=\"prot-design-design_matrix_preview\"",
    html,
    fixed = TRUE
  ))
  contrasts_heading_pos <- as.integer(regexpr(
    "Defined Contrasts",
    html,
    fixed = TRUE
  ))
  contrasts_preview_pos <- as.integer(regexpr(
    "id=\"prot-design-contrasts_preview\"",
    html,
    fixed = TRUE
  ))

  expect_identical(count_occurrences("Saved Results Preview"), 1L)
  expect_identical(count_occurrences("Current Design Matrix"), 1L)
  expect_identical(count_occurrences("Defined Contrasts"), 1L)
  expect_identical(count_occurrences("id=\"prot-design-design_matrix_preview\""), 1L)
  expect_identical(count_occurrences("id=\"prot-design-contrasts_preview\""), 1L)
  expect_gt(current_heading_pos, preview_heading_pos)
  expect_gt(design_preview_pos, current_heading_pos)
  expect_gt(contrasts_heading_pos, design_preview_pos)
  expect_gt(contrasts_preview_pos, contrasts_heading_pos)
  expect_false(grepl("id=\"design-design_matrix_preview\"", html, fixed = TRUE))
  expect_false(grepl("id=\"design-contrasts_preview\"", html, fixed = TRUE))
})

test_that("mod_prot_design_ui keeps the default saved-results preview headings ordered ahead of the preview tables", {
  ui_with_builder <- mod_prot_design_ui
  environment(ui_with_builder) <- list2env(
    list(
      exists = function(x, where = -1, inherits = TRUE) {
        if (identical(x, "mod_prot_design_builder_ui")) {
          TRUE
        } else {
          base::exists(x, where = where, inherits = inherits)
        }
      },
      mod_prot_design_builder_ui = function(id) {
        shiny::div(
          id = paste0(id, "-builder-ready"),
          `data-builder-id` = id,
          "Builder ready"
        )
      }
    ),
    parent = environment(mod_prot_design_ui)
  )

  html <- htmltools::renderTags(ui_with_builder("design"))$html

  count_occurrences <- function(pattern) {
    matches <- gregexpr(pattern, html, fixed = TRUE)[[1]]
    if (identical(matches[1], -1L)) {
      return(0L)
    }

    length(matches)
  }

  preview_heading_pos <- as.integer(regexpr(
    "Saved Results Preview",
    html,
    fixed = TRUE
  ))
  current_heading_pos <- as.integer(regexpr(
    "Current Design Matrix",
    html,
    fixed = TRUE
  ))
  design_preview_pos <- as.integer(regexpr(
    "id=\"design-design_matrix_preview\"",
    html,
    fixed = TRUE
  ))
  contrasts_heading_pos <- as.integer(regexpr(
    "Defined Contrasts",
    html,
    fixed = TRUE
  ))
  contrasts_preview_pos <- as.integer(regexpr(
    "id=\"design-contrasts_preview\"",
    html,
    fixed = TRUE
  ))

  expect_identical(count_occurrences("Saved Results Preview"), 1L)
  expect_identical(count_occurrences("Current Design Matrix"), 1L)
  expect_identical(count_occurrences("Defined Contrasts"), 1L)
  expect_identical(count_occurrences("id=\"design-design_matrix_preview\""), 1L)
  expect_identical(count_occurrences("id=\"design-contrasts_preview\""), 1L)
  expect_gt(current_heading_pos, preview_heading_pos)
  expect_gt(design_preview_pos, current_heading_pos)
  expect_gt(contrasts_heading_pos, design_preview_pos)
  expect_gt(contrasts_preview_pos, contrasts_heading_pos)
  expect_false(grepl("id=\"prot-design-design_matrix_preview\"", html, fixed = TRUE))
  expect_false(grepl("id=\"prot-design-contrasts_preview\"", html, fixed = TRUE))
})

test_that("mod_prot_design_ui keeps the default wrapper heading single-instance and ordered ahead of the saved-preview guidance", {
  ui_with_builder <- mod_prot_design_ui
  environment(ui_with_builder) <- list2env(
    list(
      exists = function(x, where = -1, inherits = TRUE) {
        if (identical(x, "mod_prot_design_builder_ui")) {
          TRUE
        } else {
          base::exists(x, where = where, inherits = inherits)
        }
      },
      mod_prot_design_builder_ui = function(id) {
        shiny::div(
          id = paste0(id, "-builder-ready"),
          `data-builder-id` = id,
          "Builder ready"
        )
      }
    ),
    parent = environment(mod_prot_design_ui)
  )

  html <- htmltools::renderTags(ui_with_builder("design"))$html

  count_occurrences <- function(pattern) {
    matches <- gregexpr(pattern, html, fixed = TRUE)[[1]]
    if (identical(matches[1], -1L)) {
      return(0L)
    }

    length(matches)
  }

  heading_pos <- as.integer(regexpr(
    "Design Matrix Builder",
    html,
    fixed = TRUE
  ))
  import_button_pos <- as.integer(regexpr(
    "id=\"design-show_import_modal\"",
    html,
    fixed = TRUE
  ))
  top_guidance_pos <- as.integer(regexpr(
    "Use the tools below to define your experimental groups and the contrasts for differential analysis. Alternatively, import an existing design from a previous analysis.",
    html,
    fixed = TRUE
  ))
  builder_pos <- as.integer(regexpr(
    "id=\"design-builder-builder-ready\"",
    html,
    fixed = TRUE
  ))
  preview_heading_pos <- as.integer(regexpr(
    "Saved Results Preview",
    html,
    fixed = TRUE
  ))
  preview_guidance_pos <- as.integer(regexpr(
    "This section shows the design matrix and contrasts that have been saved to the workflow.",
    html,
    fixed = TRUE
  ))

  expect_identical(count_occurrences("Design Matrix Builder"), 1L)
  expect_gt(heading_pos, 0L)
  expect_gt(import_button_pos, heading_pos)
  expect_gt(top_guidance_pos, import_button_pos)
  expect_gt(builder_pos, top_guidance_pos)
  expect_gt(preview_heading_pos, builder_pos)
  expect_gt(preview_guidance_pos, preview_heading_pos)
  expect_false(grepl("Design builder module not loaded", html, fixed = TRUE))
  expect_false(grepl("id=\"prot-design-show_import_modal\"", html, fixed = TRUE))
  expect_false(grepl("id=\"builder-builder-ready\"", html, fixed = TRUE))
  expect_false(grepl("data-builder-id=\"builder\"", html, fixed = TRUE))
  expect_false(grepl("id=\"prot-design-design_matrix_preview\"", html, fixed = TRUE))
  expect_false(grepl("id=\"prot-design-contrasts_preview\"", html, fixed = TRUE))
})

test_that("mod_prot_design_ui keeps the default wrapper heading ordered ahead of the fallback preview guidance", {
  ui_without_builder <- mod_prot_design_ui
  environment(ui_without_builder) <- list2env(
    list(
      exists = function(x, where = -1, inherits = TRUE) {
        if (identical(x, "mod_prot_design_builder_ui")) {
          FALSE
        } else {
          base::exists(x, where = where, inherits = inherits)
        }
      }
    ),
    parent = environment(mod_prot_design_ui)
  )

  html <- htmltools::renderTags(ui_without_builder("design"))$html

  count_occurrences <- function(pattern) {
    matches <- gregexpr(pattern, html, fixed = TRUE)[[1]]
    if (identical(matches[1], -1L)) {
      return(0L)
    }

    length(matches)
  }

  heading_pos <- as.integer(regexpr(
    "Design Matrix Builder",
    html,
    fixed = TRUE
  ))
  import_button_pos <- as.integer(regexpr(
    "id=\"design-show_import_modal\"",
    html,
    fixed = TRUE
  ))
  top_guidance_pos <- as.integer(regexpr(
    "Use the tools below to define your experimental groups and the contrasts for differential analysis. Alternatively, import an existing design from a previous analysis.",
    html,
    fixed = TRUE
  ))
  fallback_pos <- as.integer(regexpr(
    "Design builder module not loaded",
    html,
    fixed = TRUE
  ))
  preview_heading_pos <- as.integer(regexpr(
    "Saved Results Preview",
    html,
    fixed = TRUE
  ))
  preview_guidance_pos <- as.integer(regexpr(
    "This section shows the design matrix and contrasts that have been saved to the workflow.",
    html,
    fixed = TRUE
  ))
  design_preview_pos <- as.integer(regexpr(
    "id=\"design-design_matrix_preview\"",
    html,
    fixed = TRUE
  ))
  contrasts_preview_pos <- as.integer(regexpr(
    "id=\"design-contrasts_preview\"",
    html,
    fixed = TRUE
  ))

  expect_identical(count_occurrences("Design Matrix Builder"), 1L)
  expect_identical(
    count_occurrences("This section shows the design matrix and contrasts that have been saved to the workflow."),
    1L
  )
  expect_identical(count_occurrences("Design builder module not loaded"), 1L)
  expect_gt(heading_pos, 0L)
  expect_gt(import_button_pos, heading_pos)
  expect_gt(top_guidance_pos, import_button_pos)
  expect_gt(fallback_pos, top_guidance_pos)
  expect_gt(preview_heading_pos, fallback_pos)
  expect_gt(preview_guidance_pos, preview_heading_pos)
  expect_gt(design_preview_pos, preview_guidance_pos)
  expect_gt(contrasts_preview_pos, design_preview_pos)
  expect_false(grepl("id=\"show_import_modal\"", html, fixed = TRUE))
  expect_false(grepl("design-builder-builder-ready", html, fixed = TRUE))
  expect_false(grepl("data-builder-id=\"design-builder\"", html, fixed = TRUE))
  expect_false(grepl("id=\"prot-design-design_matrix_preview\"", html, fixed = TRUE))
  expect_false(grepl("id=\"prot-design-contrasts_preview\"", html, fixed = TRUE))
})

test_that("mod_prot_design_ui keeps the non-default wrapper heading ordered ahead of the fallback preview guidance", {
  ui_without_builder <- mod_prot_design_ui
  environment(ui_without_builder) <- list2env(
    list(
      exists = function(x, where = -1, inherits = TRUE) {
        if (identical(x, "mod_prot_design_builder_ui")) {
          FALSE
        } else {
          base::exists(x, where = where, inherits = inherits)
        }
      }
    ),
    parent = environment(mod_prot_design_ui)
  )

  html <- htmltools::renderTags(ui_without_builder("prot-design"))$html

  count_occurrences <- function(pattern) {
    matches <- gregexpr(pattern, html, fixed = TRUE)[[1]]
    if (identical(matches[1], -1L)) {
      return(0L)
    }

    length(matches)
  }

  heading_pos <- as.integer(regexpr(
    "Design Matrix Builder",
    html,
    fixed = TRUE
  ))
  import_button_pos <- as.integer(regexpr(
    "id=\"prot-design-show_import_modal\"",
    html,
    fixed = TRUE
  ))
  top_guidance_pos <- as.integer(regexpr(
    "Use the tools below to define your experimental groups and the contrasts for differential analysis. Alternatively, import an existing design from a previous analysis.",
    html,
    fixed = TRUE
  ))
  fallback_pos <- as.integer(regexpr(
    "Design builder module not loaded",
    html,
    fixed = TRUE
  ))
  preview_heading_pos <- as.integer(regexpr(
    "Saved Results Preview",
    html,
    fixed = TRUE
  ))
  preview_guidance_pos <- as.integer(regexpr(
    "This section shows the design matrix and contrasts that have been saved to the workflow.",
    html,
    fixed = TRUE
  ))
  design_preview_pos <- as.integer(regexpr(
    "id=\"prot-design-design_matrix_preview\"",
    html,
    fixed = TRUE
  ))
  contrasts_preview_pos <- as.integer(regexpr(
    "id=\"prot-design-contrasts_preview\"",
    html,
    fixed = TRUE
  ))

  expect_identical(count_occurrences("Design Matrix Builder"), 1L)
  expect_identical(
    count_occurrences("This section shows the design matrix and contrasts that have been saved to the workflow."),
    1L
  )
  expect_identical(count_occurrences("Design builder module not loaded"), 1L)
  expect_gt(heading_pos, 0L)
  expect_gt(import_button_pos, heading_pos)
  expect_gt(top_guidance_pos, import_button_pos)
  expect_gt(fallback_pos, top_guidance_pos)
  expect_gt(preview_heading_pos, fallback_pos)
  expect_gt(preview_guidance_pos, preview_heading_pos)
  expect_gt(design_preview_pos, preview_guidance_pos)
  expect_gt(contrasts_preview_pos, design_preview_pos)
  expect_false(grepl("id=\"show_import_modal\"", html, fixed = TRUE))
  expect_false(grepl("id=\"design-show_import_modal\"", html, fixed = TRUE))
  expect_false(grepl("design-builder-builder-ready", html, fixed = TRUE))
  expect_false(grepl("data-builder-id=\"design-builder\"", html, fixed = TRUE))
  expect_false(grepl("id=\"design-design_matrix_preview\"", html, fixed = TRUE))
  expect_false(grepl("id=\"design-contrasts_preview\"", html, fixed = TRUE))
})

test_that("mod_prot_design_ui keeps the non-default saved-preview scaffold ordered ahead of the empty-state alert", {
  ui_with_builder <- mod_prot_design_ui
  environment(ui_with_builder) <- list2env(
    list(
      exists = function(x, where = -1, inherits = TRUE) {
        if (identical(x, "mod_prot_design_builder_ui")) {
          TRUE
        } else {
          base::exists(x, where = where, inherits = inherits)
        }
      },
      mod_prot_design_builder_ui = function(id) {
        shiny::div(
          id = paste0(id, "-builder-ready"),
          `data-builder-id` = id,
          "Builder ready"
        )
      }
    ),
    parent = environment(mod_prot_design_ui)
  )

  html <- htmltools::renderTags(ui_with_builder("prot-design"))$html

  count_occurrences <- function(pattern) {
    matches <- gregexpr(pattern, html, fixed = TRUE)[[1]]
    if (identical(matches[1], -1L)) {
      return(0L)
    }

    length(matches)
  }

  builder_binding_pos <- as.integer(regexpr(
    "data-display-if=\"output[&#39;prot-design-data_available&#39;]\"",
    html,
    fixed = TRUE
  ))
  builder_pos <- as.integer(regexpr(
    "id=\"prot-design-builder-builder-ready\"",
    html,
    fixed = TRUE
  ))
  preview_heading_pos <- as.integer(regexpr(
    "Saved Results Preview",
    html,
    fixed = TRUE
  ))
  preview_guidance_pos <- as.integer(regexpr(
    "This section shows the design matrix and contrasts that have been saved to the workflow.",
    html,
    fixed = TRUE
  ))
  current_heading_pos <- as.integer(regexpr(
    "Current Design Matrix",
    html,
    fixed = TRUE
  ))
  design_preview_pos <- as.integer(regexpr(
    "id=\"prot-design-design_matrix_preview\"",
    html,
    fixed = TRUE
  ))
  contrasts_heading_pos <- as.integer(regexpr(
    "Defined Contrasts",
    html,
    fixed = TRUE
  ))
  contrasts_preview_pos <- as.integer(regexpr(
    "id=\"prot-design-contrasts_preview\"",
    html,
    fixed = TRUE
  ))
  empty_state_binding_pos <- as.integer(regexpr(
    "data-display-if=\"!output[&#39;prot-design-data_available&#39;]\"",
    html,
    fixed = TRUE
  ))
  alert_pos <- as.integer(regexpr(
    "alert alert-info",
    html,
    fixed = TRUE
  ))
  empty_state_guidance_pos <- as.integer(regexpr(
    "Please complete the 'Setup &amp; Import' step first. The builder will appear here once data is available.",
    html,
    fixed = TRUE
  ))

  expect_identical(count_occurrences("Saved Results Preview"), 1L)
  expect_identical(count_occurrences("Current Design Matrix"), 1L)
  expect_identical(count_occurrences("Defined Contrasts"), 1L)
  expect_identical(count_occurrences("alert alert-info"), 1L)
  expect_gt(builder_binding_pos, 0L)
  expect_gt(builder_pos, builder_binding_pos)
  expect_gt(preview_heading_pos, builder_pos)
  expect_gt(preview_guidance_pos, preview_heading_pos)
  expect_gt(current_heading_pos, preview_guidance_pos)
  expect_gt(design_preview_pos, current_heading_pos)
  expect_gt(contrasts_heading_pos, design_preview_pos)
  expect_gt(contrasts_preview_pos, contrasts_heading_pos)
  expect_gt(empty_state_binding_pos, contrasts_preview_pos)
  expect_gt(alert_pos, empty_state_binding_pos)
  expect_gt(empty_state_guidance_pos, alert_pos)
  expect_false(grepl("id=\"design-design_matrix_preview\"", html, fixed = TRUE))
  expect_false(grepl("id=\"design-contrasts_preview\"", html, fixed = TRUE))
  expect_false(grepl("id=\"design-show_import_modal\"", html, fixed = TRUE))
  expect_false(grepl("Design builder module not loaded", html, fixed = TRUE))
})

test_that("mod_prot_design_ui keeps the default saved-preview scaffold ordered ahead of the empty-state alert", {
  ui_with_builder <- mod_prot_design_ui
  environment(ui_with_builder) <- list2env(
    list(
      exists = function(x, where = -1, inherits = TRUE) {
        if (identical(x, "mod_prot_design_builder_ui")) {
          TRUE
        } else {
          base::exists(x, where = where, inherits = inherits)
        }
      },
      mod_prot_design_builder_ui = function(id) {
        shiny::div(
          id = paste0(id, "-builder-ready"),
          `data-builder-id` = id,
          "Builder ready"
        )
      }
    ),
    parent = environment(mod_prot_design_ui)
  )

  html <- htmltools::renderTags(ui_with_builder("design"))$html

  count_occurrences <- function(pattern) {
    matches <- gregexpr(pattern, html, fixed = TRUE)[[1]]
    if (identical(matches[1], -1L)) {
      return(0L)
    }

    length(matches)
  }

  builder_binding_pos <- as.integer(regexpr(
    "data-display-if=\"output[&#39;design-data_available&#39;]\"",
    html,
    fixed = TRUE
  ))
  builder_pos <- as.integer(regexpr(
    "id=\"design-builder-builder-ready\"",
    html,
    fixed = TRUE
  ))
  preview_heading_pos <- as.integer(regexpr(
    "Saved Results Preview",
    html,
    fixed = TRUE
  ))
  preview_guidance_pos <- as.integer(regexpr(
    "This section shows the design matrix and contrasts that have been saved to the workflow.",
    html,
    fixed = TRUE
  ))
  current_heading_pos <- as.integer(regexpr(
    "Current Design Matrix",
    html,
    fixed = TRUE
  ))
  design_preview_pos <- as.integer(regexpr(
    "id=\"design-design_matrix_preview\"",
    html,
    fixed = TRUE
  ))
  contrasts_heading_pos <- as.integer(regexpr(
    "Defined Contrasts",
    html,
    fixed = TRUE
  ))
  contrasts_preview_pos <- as.integer(regexpr(
    "id=\"design-contrasts_preview\"",
    html,
    fixed = TRUE
  ))
  empty_state_binding_pos <- as.integer(regexpr(
    "data-display-if=\"!output[&#39;design-data_available&#39;]\"",
    html,
    fixed = TRUE
  ))
  alert_pos <- as.integer(regexpr(
    "alert alert-info",
    html,
    fixed = TRUE
  ))
  empty_state_guidance_pos <- as.integer(regexpr(
    "Please complete the 'Setup &amp; Import' step first. The builder will appear here once data is available.",
    html,
    fixed = TRUE
  ))

  expect_identical(count_occurrences("Saved Results Preview"), 1L)
  expect_identical(count_occurrences("Current Design Matrix"), 1L)
  expect_identical(count_occurrences("Defined Contrasts"), 1L)
  expect_identical(count_occurrences("alert alert-info"), 1L)
  expect_gt(builder_binding_pos, 0L)
  expect_gt(builder_pos, builder_binding_pos)
  expect_gt(preview_heading_pos, builder_pos)
  expect_gt(preview_guidance_pos, preview_heading_pos)
  expect_gt(current_heading_pos, preview_guidance_pos)
  expect_gt(design_preview_pos, current_heading_pos)
  expect_gt(contrasts_heading_pos, design_preview_pos)
  expect_gt(contrasts_preview_pos, contrasts_heading_pos)
  expect_gt(empty_state_binding_pos, contrasts_preview_pos)
  expect_gt(alert_pos, empty_state_binding_pos)
  expect_gt(empty_state_guidance_pos, alert_pos)
  expect_false(grepl("id=\"prot-design-design_matrix_preview\"", html, fixed = TRUE))
  expect_false(grepl("id=\"prot-design-contrasts_preview\"", html, fixed = TRUE))
  expect_false(grepl("id=\"prot-design-show_import_modal\"", html, fixed = TRUE))
  expect_false(grepl("Design builder module not loaded", html, fixed = TRUE))
})

test_that("mod_prot_design_ui keeps the non-default fallback saved-preview scaffold ordered ahead of the empty-state alert", {
  ui_without_builder <- mod_prot_design_ui
  environment(ui_without_builder) <- list2env(
    list(
      exists = function(x, where = -1, inherits = TRUE) {
        if (identical(x, "mod_prot_design_builder_ui")) {
          FALSE
        } else {
          base::exists(x, where = where, inherits = inherits)
        }
      }
    ),
    parent = environment(mod_prot_design_ui)
  )

  html <- htmltools::renderTags(ui_without_builder("prot-design"))$html

  count_occurrences <- function(pattern) {
    matches <- gregexpr(pattern, html, fixed = TRUE)[[1]]
    if (identical(matches[1], -1L)) {
      return(0L)
    }

    length(matches)
  }

  builder_binding_pos <- as.integer(regexpr(
    "data-display-if=\"output[&#39;prot-design-data_available&#39;]\"",
    html,
    fixed = TRUE
  ))
  fallback_pos <- as.integer(regexpr(
    "Design builder module not loaded",
    html,
    fixed = TRUE
  ))
  preview_heading_pos <- as.integer(regexpr(
    "Saved Results Preview",
    html,
    fixed = TRUE
  ))
  preview_guidance_pos <- as.integer(regexpr(
    "This section shows the design matrix and contrasts that have been saved to the workflow.",
    html,
    fixed = TRUE
  ))
  current_heading_pos <- as.integer(regexpr(
    "Current Design Matrix",
    html,
    fixed = TRUE
  ))
  design_preview_pos <- as.integer(regexpr(
    "id=\"prot-design-design_matrix_preview\"",
    html,
    fixed = TRUE
  ))
  contrasts_heading_pos <- as.integer(regexpr(
    "Defined Contrasts",
    html,
    fixed = TRUE
  ))
  contrasts_preview_pos <- as.integer(regexpr(
    "id=\"prot-design-contrasts_preview\"",
    html,
    fixed = TRUE
  ))
  empty_state_binding_pos <- as.integer(regexpr(
    "data-display-if=\"!output[&#39;prot-design-data_available&#39;]\"",
    html,
    fixed = TRUE
  ))
  alert_pos <- as.integer(regexpr(
    "alert alert-info",
    html,
    fixed = TRUE
  ))
  empty_state_guidance_pos <- as.integer(regexpr(
    "Please complete the 'Setup &amp; Import' step first. The builder will appear here once data is available.",
    html,
    fixed = TRUE
  ))

  expect_identical(count_occurrences("Saved Results Preview"), 1L)
  expect_identical(count_occurrences("Current Design Matrix"), 1L)
  expect_identical(count_occurrences("Defined Contrasts"), 1L)
  expect_identical(count_occurrences("alert alert-info"), 1L)
  expect_identical(count_occurrences("Design builder module not loaded"), 1L)
  expect_gt(builder_binding_pos, 0L)
  expect_gt(fallback_pos, builder_binding_pos)
  expect_gt(preview_heading_pos, fallback_pos)
  expect_gt(preview_guidance_pos, preview_heading_pos)
  expect_gt(current_heading_pos, preview_guidance_pos)
  expect_gt(design_preview_pos, current_heading_pos)
  expect_gt(contrasts_heading_pos, design_preview_pos)
  expect_gt(contrasts_preview_pos, contrasts_heading_pos)
  expect_gt(empty_state_binding_pos, contrasts_preview_pos)
  expect_gt(alert_pos, empty_state_binding_pos)
  expect_gt(empty_state_guidance_pos, alert_pos)
  expect_false(grepl("id=\"design-design_matrix_preview\"", html, fixed = TRUE))
  expect_false(grepl("id=\"design-contrasts_preview\"", html, fixed = TRUE))
  expect_false(grepl("id=\"design-show_import_modal\"", html, fixed = TRUE))
  expect_false(grepl("prot-design-builder-builder-ready", html, fixed = TRUE))
})

test_that("mod_prot_design_ui keeps the default preview divider between the builder shell and preview heading", {
  ui_with_builder <- mod_prot_design_ui
  environment(ui_with_builder) <- list2env(
    list(
      exists = function(x, where = -1, inherits = TRUE) {
        if (identical(x, "mod_prot_design_builder_ui")) {
          TRUE
        } else {
          base::exists(x, where = where, inherits = inherits)
        }
      },
      mod_prot_design_builder_ui = function(id) {
        shiny::div(
          id = paste0(id, "-builder-ready"),
          `data-builder-id` = id,
          "Builder ready"
        )
      }
    ),
    parent = environment(mod_prot_design_ui)
  )

  html <- htmltools::renderTags(ui_with_builder("design"))$html

  count_occurrences <- function(pattern) {
    matches <- gregexpr(pattern, html, fixed = TRUE)[[1]]
    if (identical(matches[1], -1L)) {
      return(0L)
    }

    length(matches)
  }

  builder_binding_pos <- as.integer(regexpr(
    "data-display-if=\"output[&#39;design-data_available&#39;]\"",
    html,
    fixed = TRUE
  ))
  builder_pos <- as.integer(regexpr(
    "id=\"design-builder-builder-ready\"",
    html,
    fixed = TRUE
  ))
  divider_pos <- as.integer(regexpr(
    "<hr/>",
    html,
    fixed = TRUE
  ))
  preview_heading_pos <- as.integer(regexpr(
    "Saved Results Preview",
    html,
    fixed = TRUE
  ))
  preview_guidance_pos <- as.integer(regexpr(
    "This section shows the design matrix and contrasts that have been saved to the workflow.",
    html,
    fixed = TRUE
  ))

  expect_identical(count_occurrences("<hr/>"), 1L)
  expect_identical(count_occurrences("Saved Results Preview"), 1L)
  expect_gt(builder_binding_pos, 0L)
  expect_gt(builder_pos, builder_binding_pos)
  expect_gt(divider_pos, builder_pos)
  expect_gt(preview_heading_pos, divider_pos)
  expect_gt(preview_guidance_pos, preview_heading_pos)
  expect_false(grepl("id=\"prot-design-show_import_modal\"", html, fixed = TRUE))
  expect_false(grepl("id=\"prot-design-design_matrix_preview\"", html, fixed = TRUE))
  expect_false(grepl("id=\"prot-design-contrasts_preview\"", html, fixed = TRUE))
  expect_false(grepl("Design builder module not loaded", html, fixed = TRUE))
})

test_that("mod_prot_design_ui keeps the non-default preview divider between the builder shell and preview heading", {
  ui_with_builder <- mod_prot_design_ui
  environment(ui_with_builder) <- list2env(
    list(
      exists = function(x, where = -1, inherits = TRUE) {
        if (identical(x, "mod_prot_design_builder_ui")) {
          TRUE
        } else {
          base::exists(x, where = where, inherits = inherits)
        }
      },
      mod_prot_design_builder_ui = function(id) {
        shiny::div(
          id = paste0(id, "-builder-ready"),
          `data-builder-id` = id,
          "Builder ready"
        )
      }
    ),
    parent = environment(mod_prot_design_ui)
  )

  html <- htmltools::renderTags(ui_with_builder("prot-design"))$html

  count_occurrences <- function(pattern) {
    matches <- gregexpr(pattern, html, fixed = TRUE)[[1]]
    if (identical(matches[1], -1L)) {
      return(0L)
    }

    length(matches)
  }

  builder_binding_pos <- as.integer(regexpr(
    "data-display-if=\"output[&#39;prot-design-data_available&#39;]\"",
    html,
    fixed = TRUE
  ))
  builder_pos <- as.integer(regexpr(
    "id=\"prot-design-builder-builder-ready\"",
    html,
    fixed = TRUE
  ))
  divider_pos <- as.integer(regexpr(
    "<hr/>",
    html,
    fixed = TRUE
  ))
  preview_heading_pos <- as.integer(regexpr(
    "Saved Results Preview",
    html,
    fixed = TRUE
  ))
  preview_guidance_pos <- as.integer(regexpr(
    "This section shows the design matrix and contrasts that have been saved to the workflow.",
    html,
    fixed = TRUE
  ))

  expect_identical(count_occurrences("<hr/>"), 1L)
  expect_identical(count_occurrences("Saved Results Preview"), 1L)
  expect_gt(builder_binding_pos, 0L)
  expect_gt(builder_pos, builder_binding_pos)
  expect_gt(divider_pos, builder_pos)
  expect_gt(preview_heading_pos, divider_pos)
  expect_gt(preview_guidance_pos, preview_heading_pos)
  expect_false(grepl("id=\"design-show_import_modal\"", html, fixed = TRUE))
  expect_false(grepl("id=\"design-design_matrix_preview\"", html, fixed = TRUE))
  expect_false(grepl("id=\"design-contrasts_preview\"", html, fixed = TRUE))
  expect_false(grepl("Design builder module not loaded", html, fixed = TRUE))
})

test_that("mod_prot_design_ui keeps the default fallback preview divider between the fallback shell and preview heading", {
  ui_without_builder <- mod_prot_design_ui
  environment(ui_without_builder) <- list2env(
    list(
      exists = function(x, where = -1, inherits = TRUE) {
        if (identical(x, "mod_prot_design_builder_ui")) {
          FALSE
        } else {
          base::exists(x, where = where, inherits = inherits)
        }
      }
    ),
    parent = environment(mod_prot_design_ui)
  )

  html <- htmltools::renderTags(ui_without_builder("design"))$html

  count_occurrences <- function(pattern) {
    matches <- gregexpr(pattern, html, fixed = TRUE)[[1]]
    if (identical(matches[1], -1L)) {
      return(0L)
    }

    length(matches)
  }

  builder_binding_pos <- as.integer(regexpr(
    "data-display-if=\"output[&#39;design-data_available&#39;]\"",
    html,
    fixed = TRUE
  ))
  fallback_pos <- as.integer(regexpr(
    "Design builder module not loaded",
    html,
    fixed = TRUE
  ))
  divider_pos <- as.integer(regexpr(
    "<hr/>",
    html,
    fixed = TRUE
  ))
  preview_heading_pos <- as.integer(regexpr(
    "Saved Results Preview",
    html,
    fixed = TRUE
  ))
  preview_guidance_pos <- as.integer(regexpr(
    "This section shows the design matrix and contrasts that have been saved to the workflow.",
    html,
    fixed = TRUE
  ))

  expect_identical(count_occurrences("<hr/>"), 1L)
  expect_identical(count_occurrences("Saved Results Preview"), 1L)
  expect_identical(count_occurrences("Design builder module not loaded"), 1L)
  expect_gt(builder_binding_pos, 0L)
  expect_gt(fallback_pos, builder_binding_pos)
  expect_gt(divider_pos, fallback_pos)
  expect_gt(preview_heading_pos, divider_pos)
  expect_gt(preview_guidance_pos, preview_heading_pos)
  expect_false(grepl("id=\"prot-design-show_import_modal\"", html, fixed = TRUE))
  expect_false(grepl("id=\"prot-design-design_matrix_preview\"", html, fixed = TRUE))
  expect_false(grepl("id=\"prot-design-contrasts_preview\"", html, fixed = TRUE))
  expect_false(grepl("id=\"design-builder-builder-ready\"", html, fixed = TRUE))
})

test_that("mod_prot_design_ui keeps the default saved-preview content behind the design-matrix-exists binding", {
  ui_with_builder <- mod_prot_design_ui
  environment(ui_with_builder) <- list2env(
    list(
      exists = function(x, where = -1, inherits = TRUE) {
        if (identical(x, "mod_prot_design_builder_ui")) {
          TRUE
        } else {
          base::exists(x, where = where, inherits = inherits)
        }
      },
      mod_prot_design_builder_ui = function(id) {
        shiny::div(
          id = paste0(id, "-builder-ready"),
          `data-builder-id` = id,
          "Builder ready"
        )
      }
    ),
    parent = environment(mod_prot_design_ui)
  )

  html <- htmltools::renderTags(ui_with_builder("design"))$html

  count_occurrences <- function(pattern) {
    matches <- gregexpr(pattern, html, fixed = TRUE)[[1]]
    if (identical(matches[1], -1L)) {
      return(0L)
    }

    length(matches)
  }

  preview_guidance_pos <- as.integer(regexpr(
    "This section shows the design matrix and contrasts that have been saved to the workflow.",
    html,
    fixed = TRUE
  ))
  design_matrix_exists_binding_pos <- as.integer(regexpr(
    "data-display-if=\"output[&#39;design-design_matrix_exists&#39;]\"",
    html,
    fixed = TRUE
  ))
  current_heading_pos <- as.integer(regexpr(
    "Current Design Matrix",
    html,
    fixed = TRUE
  ))
  design_preview_pos <- as.integer(regexpr(
    "id=\"design-design_matrix_preview\"",
    html,
    fixed = TRUE
  ))
  contrasts_heading_pos <- as.integer(regexpr(
    "Defined Contrasts",
    html,
    fixed = TRUE
  ))
  contrasts_preview_pos <- as.integer(regexpr(
    "id=\"design-contrasts_preview\"",
    html,
    fixed = TRUE
  ))

  expect_identical(
    count_occurrences("data-display-if=\"output[&#39;design-design_matrix_exists&#39;]\""),
    1L
  )
  expect_identical(count_occurrences("Current Design Matrix"), 1L)
  expect_identical(count_occurrences("Defined Contrasts"), 1L)
  expect_identical(count_occurrences("id=\"design-design_matrix_preview\""), 1L)
  expect_identical(count_occurrences("id=\"design-contrasts_preview\""), 1L)
  expect_gt(preview_guidance_pos, 0L)
  expect_gt(design_matrix_exists_binding_pos, preview_guidance_pos)
  expect_gt(current_heading_pos, design_matrix_exists_binding_pos)
  expect_gt(design_preview_pos, current_heading_pos)
  expect_gt(contrasts_heading_pos, design_preview_pos)
  expect_gt(contrasts_preview_pos, contrasts_heading_pos)
  expect_false(grepl("id=\"prot-design-design_matrix_preview\"", html, fixed = TRUE))
  expect_false(grepl("id=\"prot-design-contrasts_preview\"", html, fixed = TRUE))
})

test_that("mod_prot_design_ui keeps the non-default saved-preview content behind the design-matrix-exists binding", {
  ui_with_builder <- mod_prot_design_ui
  environment(ui_with_builder) <- list2env(
    list(
      exists = function(x, where = -1, inherits = TRUE) {
        if (identical(x, "mod_prot_design_builder_ui")) {
          TRUE
        } else {
          base::exists(x, where = where, inherits = inherits)
        }
      },
      mod_prot_design_builder_ui = function(id) {
        shiny::div(
          id = paste0(id, "-builder-ready"),
          `data-builder-id` = id,
          "Builder ready"
        )
      }
    ),
    parent = environment(mod_prot_design_ui)
  )

  html <- htmltools::renderTags(ui_with_builder("prot-design"))$html

  count_occurrences <- function(pattern) {
    matches <- gregexpr(pattern, html, fixed = TRUE)[[1]]
    if (identical(matches[1], -1L)) {
      return(0L)
    }

    length(matches)
  }

  preview_guidance_pos <- as.integer(regexpr(
    "This section shows the design matrix and contrasts that have been saved to the workflow.",
    html,
    fixed = TRUE
  ))
  design_matrix_exists_binding_pos <- as.integer(regexpr(
    "data-display-if=\"output[&#39;prot-design-design_matrix_exists&#39;]\"",
    html,
    fixed = TRUE
  ))
  current_heading_pos <- as.integer(regexpr(
    "Current Design Matrix",
    html,
    fixed = TRUE
  ))
  design_preview_pos <- as.integer(regexpr(
    "id=\"prot-design-design_matrix_preview\"",
    html,
    fixed = TRUE
  ))
  contrasts_heading_pos <- as.integer(regexpr(
    "Defined Contrasts",
    html,
    fixed = TRUE
  ))
  contrasts_preview_pos <- as.integer(regexpr(
    "id=\"prot-design-contrasts_preview\"",
    html,
    fixed = TRUE
  ))

  expect_identical(
    count_occurrences("data-display-if=\"output[&#39;prot-design-design_matrix_exists&#39;]\""),
    1L
  )
  expect_identical(count_occurrences("Current Design Matrix"), 1L)
  expect_identical(count_occurrences("Defined Contrasts"), 1L)
  expect_identical(count_occurrences("id=\"prot-design-design_matrix_preview\""), 1L)
  expect_identical(count_occurrences("id=\"prot-design-contrasts_preview\""), 1L)
  expect_gt(preview_guidance_pos, 0L)
  expect_gt(design_matrix_exists_binding_pos, preview_guidance_pos)
  expect_gt(current_heading_pos, design_matrix_exists_binding_pos)
  expect_gt(design_preview_pos, current_heading_pos)
  expect_gt(contrasts_heading_pos, design_preview_pos)
  expect_gt(contrasts_preview_pos, contrasts_heading_pos)
  expect_false(grepl("id=\"design-design_matrix_preview\"", html, fixed = TRUE))
  expect_false(grepl("id=\"design-contrasts_preview\"", html, fixed = TRUE))
})

test_that("mod_prot_design_ui keeps the non-default saved-preview well shell inside the design-matrix-exists binding", {
  ui_with_builder <- mod_prot_design_ui
  environment(ui_with_builder) <- list2env(
    list(
      exists = function(x, where = -1, inherits = TRUE) {
        if (identical(x, "mod_prot_design_builder_ui")) {
          TRUE
        } else {
          base::exists(x, where = where, inherits = inherits)
        }
      },
      mod_prot_design_builder_ui = function(id) {
        shiny::div(
          id = paste0(id, "-builder-ready"),
          `data-builder-id` = id,
          "Builder ready"
        )
      }
    ),
    parent = environment(mod_prot_design_ui)
  )

  html <- htmltools::renderTags(ui_with_builder("prot-design"))$html

  preview_binding_pos <- as.integer(regexpr(
    "data-display-if=\"output[&#39;prot-design-design_matrix_exists&#39;]\"",
    html,
    fixed = TRUE
  ))
  preview_subtree <- substr(html, preview_binding_pos, nchar(html))

  count_occurrences <- function(text, pattern) {
    matches <- gregexpr(pattern, text, fixed = TRUE)[[1]]
    if (identical(matches[1], -1L)) {
      return(0L)
    }

    length(matches)
  }

  preview_well_pos <- as.integer(regexpr(
    "<div class=\"well\">",
    preview_subtree,
    fixed = TRUE
  ))
  current_heading_pos <- as.integer(regexpr(
    "Current Design Matrix",
    preview_subtree,
    fixed = TRUE
  ))
  design_preview_pos <- as.integer(regexpr(
    "id=\"prot-design-design_matrix_preview\"",
    preview_subtree,
    fixed = TRUE
  ))
  contrasts_heading_pos <- as.integer(regexpr(
    "Defined Contrasts",
    preview_subtree,
    fixed = TRUE
  ))
  contrasts_preview_pos <- as.integer(regexpr(
    "id=\"prot-design-contrasts_preview\"",
    preview_subtree,
    fixed = TRUE
  ))

  expect_gt(preview_binding_pos, 0L)
  expect_identical(count_occurrences(preview_subtree, "<div class=\"well\">"), 1L)
  expect_gt(preview_well_pos, 0L)
  expect_gt(current_heading_pos, preview_well_pos)
  expect_gt(design_preview_pos, current_heading_pos)
  expect_gt(contrasts_heading_pos, design_preview_pos)
  expect_gt(contrasts_preview_pos, contrasts_heading_pos)
  expect_false(grepl("id=\"design-design_matrix_preview\"", preview_subtree, fixed = TRUE))
  expect_false(grepl("id=\"design-contrasts_preview\"", preview_subtree, fixed = TRUE))
})

test_that("mod_prot_design_ui keeps the default saved-preview well shell inside the design-matrix-exists binding", {
  ui_with_builder <- mod_prot_design_ui
  environment(ui_with_builder) <- list2env(
    list(
      exists = function(x, where = -1, inherits = TRUE) {
        if (identical(x, "mod_prot_design_builder_ui")) {
          TRUE
        } else {
          base::exists(x, where = where, inherits = inherits)
        }
      },
      mod_prot_design_builder_ui = function(id) {
        shiny::div(
          id = paste0(id, "-builder-ready"),
          `data-builder-id` = id,
          "Builder ready"
        )
      }
    ),
    parent = environment(mod_prot_design_ui)
  )

  html <- htmltools::renderTags(ui_with_builder("design"))$html

  preview_binding_pos <- as.integer(regexpr(
    "data-display-if=\"output[&#39;design-design_matrix_exists&#39;]\"",
    html,
    fixed = TRUE
  ))
  preview_subtree <- substr(html, preview_binding_pos, nchar(html))

  count_occurrences <- function(text, pattern) {
    matches <- gregexpr(pattern, text, fixed = TRUE)[[1]]
    if (identical(matches[1], -1L)) {
      return(0L)
    }

    length(matches)
  }

  preview_well_pos <- as.integer(regexpr(
    "<div class=\"well\">",
    preview_subtree,
    fixed = TRUE
  ))
  current_heading_pos <- as.integer(regexpr(
    "Current Design Matrix",
    preview_subtree,
    fixed = TRUE
  ))
  design_preview_pos <- as.integer(regexpr(
    "id=\"design-design_matrix_preview\"",
    preview_subtree,
    fixed = TRUE
  ))
  contrasts_heading_pos <- as.integer(regexpr(
    "Defined Contrasts",
    preview_subtree,
    fixed = TRUE
  ))
  contrasts_preview_pos <- as.integer(regexpr(
    "id=\"design-contrasts_preview\"",
    preview_subtree,
    fixed = TRUE
  ))

  expect_gt(preview_binding_pos, 0L)
  expect_identical(count_occurrences(preview_subtree, "<div class=\"well\">"), 1L)
  expect_gt(preview_well_pos, 0L)
  expect_gt(current_heading_pos, preview_well_pos)
  expect_gt(design_preview_pos, current_heading_pos)
  expect_gt(contrasts_heading_pos, design_preview_pos)
  expect_gt(contrasts_preview_pos, contrasts_heading_pos)
  expect_false(grepl("id=\"prot-design-design_matrix_preview\"", preview_subtree, fixed = TRUE))
  expect_false(grepl("id=\"prot-design-contrasts_preview\"", preview_subtree, fixed = TRUE))
})

test_that("mod_prot_design_ui keeps the default fallback saved-preview well shell inside the design-matrix-exists binding", {
  ui_without_builder <- mod_prot_design_ui
  environment(ui_without_builder) <- list2env(
    list(
      exists = function(x, where = -1, inherits = TRUE) {
        if (identical(x, "mod_prot_design_builder_ui")) {
          FALSE
        } else {
          base::exists(x, where = where, inherits = inherits)
        }
      }
    ),
    parent = environment(mod_prot_design_ui)
  )

  html <- htmltools::renderTags(ui_without_builder("design"))$html

  preview_binding_pos <- as.integer(regexpr(
    "data-display-if=\"output[&#39;design-design_matrix_exists&#39;]\"",
    html,
    fixed = TRUE
  ))
  preview_subtree <- substr(html, preview_binding_pos, nchar(html))

  count_occurrences <- function(text, pattern) {
    matches <- gregexpr(pattern, text, fixed = TRUE)[[1]]
    if (identical(matches[1], -1L)) {
      return(0L)
    }

    length(matches)
  }

  preview_well_pos <- as.integer(regexpr(
    "<div class=\"well\">",
    preview_subtree,
    fixed = TRUE
  ))
  current_heading_pos <- as.integer(regexpr(
    "Current Design Matrix",
    preview_subtree,
    fixed = TRUE
  ))
  design_preview_pos <- as.integer(regexpr(
    "id=\"design-design_matrix_preview\"",
    preview_subtree,
    fixed = TRUE
  ))
  contrasts_heading_pos <- as.integer(regexpr(
    "Defined Contrasts",
    preview_subtree,
    fixed = TRUE
  ))
  contrasts_preview_pos <- as.integer(regexpr(
    "id=\"design-contrasts_preview\"",
    preview_subtree,
    fixed = TRUE
  ))

  expect_gt(preview_binding_pos, 0L)
  expect_identical(count_occurrences(preview_subtree, "<div class=\"well\">"), 1L)
  expect_gt(preview_well_pos, 0L)
  expect_gt(current_heading_pos, preview_well_pos)
  expect_gt(design_preview_pos, current_heading_pos)
  expect_gt(contrasts_heading_pos, design_preview_pos)
  expect_gt(contrasts_preview_pos, contrasts_heading_pos)
  expect_false(grepl("design-builder-builder-ready", preview_subtree, fixed = TRUE))
  expect_false(grepl("id=\"prot-design-design_matrix_preview\"", preview_subtree, fixed = TRUE))
  expect_false(grepl("id=\"prot-design-contrasts_preview\"", preview_subtree, fixed = TRUE))
})

test_that("mod_prot_design_ui keeps the non-default fallback saved-preview well shell inside the design-matrix-exists binding", {
  ui_without_builder <- mod_prot_design_ui
  environment(ui_without_builder) <- list2env(
    list(
      exists = function(x, where = -1, inherits = TRUE) {
        if (identical(x, "mod_prot_design_builder_ui")) {
          FALSE
        } else {
          base::exists(x, where = where, inherits = inherits)
        }
      }
    ),
    parent = environment(mod_prot_design_ui)
  )

  html <- htmltools::renderTags(ui_without_builder("prot-design"))$html

  preview_binding_pos <- as.integer(regexpr(
    "data-display-if=\"output[&#39;prot-design-design_matrix_exists&#39;]\"",
    html,
    fixed = TRUE
  ))
  preview_subtree <- substr(html, preview_binding_pos, nchar(html))

  count_occurrences <- function(text, pattern) {
    matches <- gregexpr(pattern, text, fixed = TRUE)[[1]]
    if (identical(matches[1], -1L)) {
      return(0L)
    }

    length(matches)
  }

  preview_well_pos <- as.integer(regexpr(
    "<div class=\"well\">",
    preview_subtree,
    fixed = TRUE
  ))
  current_heading_pos <- as.integer(regexpr(
    "Current Design Matrix",
    preview_subtree,
    fixed = TRUE
  ))
  design_preview_pos <- as.integer(regexpr(
    "id=\"prot-design-design_matrix_preview\"",
    preview_subtree,
    fixed = TRUE
  ))
  contrasts_heading_pos <- as.integer(regexpr(
    "Defined Contrasts",
    preview_subtree,
    fixed = TRUE
  ))
  contrasts_preview_pos <- as.integer(regexpr(
    "id=\"prot-design-contrasts_preview\"",
    preview_subtree,
    fixed = TRUE
  ))

  expect_gt(preview_binding_pos, 0L)
  expect_identical(count_occurrences(preview_subtree, "<div class=\"well\">"), 1L)
  expect_gt(preview_well_pos, 0L)
  expect_gt(current_heading_pos, preview_well_pos)
  expect_gt(design_preview_pos, current_heading_pos)
  expect_gt(contrasts_heading_pos, design_preview_pos)
  expect_gt(contrasts_preview_pos, contrasts_heading_pos)
  expect_false(grepl("prot-design-builder-builder-ready", preview_subtree, fixed = TRUE))
  expect_false(grepl("id=\"design-design_matrix_preview\"", preview_subtree, fixed = TRUE))
  expect_false(grepl("id=\"design-contrasts_preview\"", preview_subtree, fixed = TRUE))
})

test_that("mod_prot_design_ui keeps the default saved-preview spacer between the preview tables", {
  ui_with_builder <- mod_prot_design_ui
  environment(ui_with_builder) <- list2env(
    list(
      exists = function(x, where = -1, inherits = TRUE) {
        if (identical(x, "mod_prot_design_builder_ui")) {
          TRUE
        } else {
          base::exists(x, where = where, inherits = inherits)
        }
      },
      mod_prot_design_builder_ui = function(id) {
        shiny::div(
          id = paste0(id, "-builder-ready"),
          `data-builder-id` = id,
          "Builder ready"
        )
      }
    ),
    parent = environment(mod_prot_design_ui)
  )

  html <- htmltools::renderTags(ui_with_builder("design"))$html

  count_occurrences <- function(text, pattern) {
    matches <- gregexpr(pattern, text, fixed = TRUE)[[1]]
    if (identical(matches[1], -1L)) {
      return(0L)
    }

    length(matches)
  }

  preview_binding_pos <- as.integer(regexpr(
    "data-display-if=\"output[&#39;design-design_matrix_exists&#39;]\"",
    html,
    fixed = TRUE
  ))
  preview_subtree <- substr(html, preview_binding_pos, nchar(html))
  design_preview_pos <- as.integer(regexpr(
    "id=\"design-design_matrix_preview\"",
    preview_subtree,
    fixed = TRUE
  ))
  spacer_pos <- as.integer(regexpr(
    "<br/>",
    preview_subtree,
    fixed = TRUE
  ))
  contrasts_heading_pos <- as.integer(regexpr(
    "Defined Contrasts",
    preview_subtree,
    fixed = TRUE
  ))
  contrasts_preview_pos <- as.integer(regexpr(
    "id=\"design-contrasts_preview\"",
    preview_subtree,
    fixed = TRUE
  ))

  expect_gt(preview_binding_pos, 0L)
  expect_identical(count_occurrences(preview_subtree, "<br/>"), 1L)
  expect_gt(design_preview_pos, 0L)
  expect_gt(spacer_pos, design_preview_pos)
  expect_gt(contrasts_heading_pos, spacer_pos)
  expect_gt(contrasts_preview_pos, contrasts_heading_pos)
  expect_false(grepl("id=\"prot-design-design_matrix_preview\"", preview_subtree, fixed = TRUE))
  expect_false(grepl("id=\"prot-design-contrasts_preview\"", preview_subtree, fixed = TRUE))
})

test_that("mod_prot_design_ui keeps the non-default saved-preview spacer between the preview tables", {
  ui_with_builder <- mod_prot_design_ui
  environment(ui_with_builder) <- list2env(
    list(
      exists = function(x, where = -1, inherits = TRUE) {
        if (identical(x, "mod_prot_design_builder_ui")) {
          TRUE
        } else {
          base::exists(x, where = where, inherits = inherits)
        }
      },
      mod_prot_design_builder_ui = function(id) {
        shiny::div(
          id = paste0(id, "-builder-ready"),
          `data-builder-id` = id,
          "Builder ready"
        )
      }
    ),
    parent = environment(mod_prot_design_ui)
  )

  html <- htmltools::renderTags(ui_with_builder("prot-design"))$html

  count_occurrences <- function(text, pattern) {
    matches <- gregexpr(pattern, text, fixed = TRUE)[[1]]
    if (identical(matches[1], -1L)) {
      return(0L)
    }

    length(matches)
  }

  preview_binding_pos <- as.integer(regexpr(
    "data-display-if=\"output[&#39;prot-design-design_matrix_exists&#39;]\"",
    html,
    fixed = TRUE
  ))
  preview_subtree <- substr(html, preview_binding_pos, nchar(html))
  design_preview_pos <- as.integer(regexpr(
    "id=\"prot-design-design_matrix_preview\"",
    preview_subtree,
    fixed = TRUE
  ))
  spacer_pos <- as.integer(regexpr(
    "<br/>",
    preview_subtree,
    fixed = TRUE
  ))
  contrasts_heading_pos <- as.integer(regexpr(
    "Defined Contrasts",
    preview_subtree,
    fixed = TRUE
  ))
  contrasts_preview_pos <- as.integer(regexpr(
    "id=\"prot-design-contrasts_preview\"",
    preview_subtree,
    fixed = TRUE
  ))

  expect_gt(preview_binding_pos, 0L)
  expect_identical(count_occurrences(preview_subtree, "<br/>"), 1L)
  expect_gt(design_preview_pos, 0L)
  expect_gt(spacer_pos, design_preview_pos)
  expect_gt(contrasts_heading_pos, spacer_pos)
  expect_gt(contrasts_preview_pos, contrasts_heading_pos)
  expect_false(grepl("id=\"design-design_matrix_preview\"", preview_subtree, fixed = TRUE))
  expect_false(grepl("id=\"design-contrasts_preview\"", preview_subtree, fixed = TRUE))
})

test_that("mod_prot_design_ui keeps the non-default fallback saved-preview spacer between the preview tables", {
  ui_without_builder <- mod_prot_design_ui
  environment(ui_without_builder) <- list2env(
    list(
      exists = function(x, where = -1, inherits = TRUE) {
        if (identical(x, "mod_prot_design_builder_ui")) {
          FALSE
        } else {
          base::exists(x, where = where, inherits = inherits)
        }
      }
    ),
    parent = environment(mod_prot_design_ui)
  )

  html <- htmltools::renderTags(ui_without_builder("prot-design"))$html

  count_occurrences <- function(text, pattern) {
    matches <- gregexpr(pattern, text, fixed = TRUE)[[1]]
    if (identical(matches[1], -1L)) {
      return(0L)
    }

    length(matches)
  }

  preview_binding_pos <- as.integer(regexpr(
    "data-display-if=\"output[&#39;prot-design-design_matrix_exists&#39;]\"",
    html,
    fixed = TRUE
  ))
  preview_subtree <- substr(html, preview_binding_pos, nchar(html))
  design_preview_pos <- as.integer(regexpr(
    "id=\"prot-design-design_matrix_preview\"",
    preview_subtree,
    fixed = TRUE
  ))
  spacer_pos <- as.integer(regexpr(
    "<br/>",
    preview_subtree,
    fixed = TRUE
  ))
  contrasts_heading_pos <- as.integer(regexpr(
    "Defined Contrasts",
    preview_subtree,
    fixed = TRUE
  ))
  contrasts_preview_pos <- as.integer(regexpr(
    "id=\"prot-design-contrasts_preview\"",
    preview_subtree,
    fixed = TRUE
  ))

  expect_gt(preview_binding_pos, 0L)
  expect_identical(count_occurrences(preview_subtree, "<br/>"), 1L)
  expect_gt(design_preview_pos, 0L)
  expect_gt(spacer_pos, design_preview_pos)
  expect_gt(contrasts_heading_pos, spacer_pos)
  expect_gt(contrasts_preview_pos, contrasts_heading_pos)
  expect_false(grepl("prot-design-builder-builder-ready", preview_subtree, fixed = TRUE))
  expect_false(grepl("id=\"design-design_matrix_preview\"", preview_subtree, fixed = TRUE))
  expect_false(grepl("id=\"design-contrasts_preview\"", preview_subtree, fixed = TRUE))
})

test_that("mod_prot_design_ui keeps the default fallback saved-preview spacer between the preview tables", {
  ui_without_builder <- mod_prot_design_ui
  environment(ui_without_builder) <- list2env(
    list(
      exists = function(x, where = -1, inherits = TRUE) {
        if (identical(x, "mod_prot_design_builder_ui")) {
          FALSE
        } else {
          base::exists(x, where = where, inherits = inherits)
        }
      }
    ),
    parent = environment(mod_prot_design_ui)
  )

  html <- htmltools::renderTags(ui_without_builder("design"))$html

  count_occurrences <- function(text, pattern) {
    matches <- gregexpr(pattern, text, fixed = TRUE)[[1]]
    if (identical(matches[1], -1L)) {
      return(0L)
    }

    length(matches)
  }

  preview_binding_pos <- as.integer(regexpr(
    "data-display-if=\"output[&#39;design-design_matrix_exists&#39;]\"",
    html,
    fixed = TRUE
  ))
  preview_subtree <- substr(html, preview_binding_pos, nchar(html))
  design_preview_pos <- as.integer(regexpr(
    "id=\"design-design_matrix_preview\"",
    preview_subtree,
    fixed = TRUE
  ))
  spacer_pos <- as.integer(regexpr(
    "<br/>",
    preview_subtree,
    fixed = TRUE
  ))
  contrasts_heading_pos <- as.integer(regexpr(
    "Defined Contrasts",
    preview_subtree,
    fixed = TRUE
  ))
  contrasts_preview_pos <- as.integer(regexpr(
    "id=\"design-contrasts_preview\"",
    preview_subtree,
    fixed = TRUE
  ))

  expect_gt(preview_binding_pos, 0L)
  expect_identical(count_occurrences(preview_subtree, "<br/>"), 1L)
  expect_gt(design_preview_pos, 0L)
  expect_gt(spacer_pos, design_preview_pos)
  expect_gt(contrasts_heading_pos, spacer_pos)
  expect_gt(contrasts_preview_pos, contrasts_heading_pos)
  expect_false(grepl("design-builder-builder-ready", preview_subtree, fixed = TRUE))
  expect_false(grepl("id=\"prot-design-design_matrix_preview\"", preview_subtree, fixed = TRUE))
  expect_false(grepl("id=\"prot-design-contrasts_preview\"", preview_subtree, fixed = TRUE))
})

# APAF Bioinformatics | test-prot-04-design.R | Approved | 2026-03-13
