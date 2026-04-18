library(testthat)

test_that("getTypeOfGrouping returns named sample groups", {
  design_matrix <- data.frame(
    Run = c("S1", "S2", "S3", "S4"),
    group = c("A", "A", "B", "B"),
    stringsAsFactors = FALSE
  )

  result <- getTypeOfGrouping(
    design_matrix = design_matrix,
    group_id = "group",
    sample_id = "Run"
  )

  expect_type(result, "list")
  expect_named(result, c("A", "B"))
  expect_equal(result$A, c("S1", "S2"))
  expect_equal(result$B, c("S3", "S4"))
})

test_that("extractResults preserves names and pulls nested results", {
  results_list <- list(
    contrast_a = list(results = data.frame(logFC = 1, stringsAsFactors = FALSE)),
    contrast_b = list(results = data.frame(logFC = -1, stringsAsFactors = FALSE))
  )

  extracted <- extractResults(results_list)

  expect_named(extracted, c("contrast_a", "contrast_b"))
  expect_s3_class(extracted$contrast_a, "data.frame")
  expect_equal(extracted$contrast_a$logFC, 1)
  expect_equal(extracted$contrast_b$logFC, -1)
})

test_that("countStatDaGenesHelper handles named limma-style result lists", {
  da_table <- list(
    "A_vs_B=groupA-groupB" = data.frame(
      logFC = c(1.2, -1.4, 0.2),
      fdr_qvalue = c(0.01, 0.02, 0.8),
      stringsAsFactors = FALSE
    )
  )

  result <- countStatDaGenesHelper(
    da_table = da_table,
    description = "RUV applied",
    comparison_column = "comparison",
    expression_column = "expression"
  )

  expect_s3_class(result, "data.frame")
  expect_true(all(c("comparison", "expression", "status", "counts", "analysis_type") %in% names(result)))
  expect_equal(unique(result$comparison), "A_vs_B")
  expect_equal(unique(result$expression), "groupA-groupB")
  expect_equal(sum(result$counts), 3)
})

test_that("getSignificantData adds fallback expression column and significance colours", {
  list_of_da_tables <- list(
    list(
      "A_vs_B" = data.frame(
        uniprot_acc = c("P1", "P2", "P3"),
        logFC = c(1.5, 0.2, -2.0),
        raw_pvalue = c(0.001, 0.1, 0.002),
        fdr_qvalue = c(0.01, 0.8, 0.02),
        stringsAsFactors = FALSE
      )
    )
  )

  result <- getSignificantData(
    list_of_da_tables = list_of_da_tables,
    list_of_descriptions = list("RUV applied"),
    row_id = uniprot_acc,
    p_value_column = raw_pvalue,
    q_value_column = fdr_qvalue,
    log_q_value_column = lqm,
    log_fc_column = logFC,
    comparison_column = "comparison",
    expression_column = "expression",
    facet_column = analysis_type,
    q_val_thresh = 0.05
  )

  expect_s3_class(result, "data.frame")
  expect_true("expression" %in% names(result))
  expect_true(all(is.na(result$expression)))
  expect_equal(as.character(result$colour), c("purple", "black", "purple"))
  expect_equal(result$comparison, rep("A_vs_B", 3))
})

test_that("createDaResultsLongFormat renames replicate columns for a single contrast", {
  lfc_qval_tbl <- data.frame(
    uniprot_acc = "P1",
    comparison = "A_vs_B",
    log_intensity = "groupA-groupB",
    analysis_type = "RUV applied",
    lqm = 2,
    colour = "purple",
    fdr_qvalue = 0.01,
    raw_pvalue = 0.001,
    log2FC = 1.25,
    stringsAsFactors = FALSE
  )

  norm_counts_input_tbl <- matrix(
    c(10, 11, 20, 21),
    nrow = 1,
    dimnames = list("P1", c("S1", "S2", "S3", "S4"))
  )
  raw_counts_input_tbl <- matrix(
    c(100, 110, 200, 210),
    nrow = 1,
    dimnames = list("P1", c("S1", "S2", "S3", "S4"))
  )

  design_matrix <- data.frame(
    Run = c("S1", "S2", "S3", "S4"),
    group = c("A", "A", "B", "B"),
    stringsAsFactors = FALSE
  )

  protein_id_table <- data.frame(
    uniprot_acc = "P1",
    gene_names = "GENE1",
    stringsAsFactors = FALSE
  )

  result <- createDaResultsLongFormat(
    lfc_qval_tbl = lfc_qval_tbl,
    norm_counts_input_tbl = norm_counts_input_tbl,
    raw_counts_input_tbl = raw_counts_input_tbl,
    row_id = "uniprot_acc",
    sample_id = "Run",
    group_id = "group",
    group_pattern = "S[0-9]+",
    design_matrix_norm = design_matrix,
    design_matrix_raw = design_matrix,
    protein_id_table = protein_id_table
  )

  expect_s3_class(result, "data.frame")
  expect_equal(result$numerator, "A")
  expect_equal(result$denominator, "B")
  expect_true(all(c(
    "log2norm.S1.A", "log2norm.S2.A", "log2norm.S3.B", "log2norm.S4.B",
    "raw.S1.A", "raw.S2.A", "raw.S3.B", "raw.S4.B"
  ) %in% names(result)))
  expect_equal(result$gene_names, "GENE1")
})
