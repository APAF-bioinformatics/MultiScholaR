# fidelity-coverage-compare: shared
library(testthat)

test_that("protein accession ranking helpers preserve grouped accession ordering and fallback metadata", {
  input_tbl <- data.frame(
    group_id = c("G1", "G2"),
    accession_bundle = c("P11111:P22222-2:REV__P33333", "CON__Q11111:Q22222"),
    stringsAsFactors = FALSE
  )

  acc_detail_tab <- data.frame(
    uniprot_acc = c("P11111", "P22222", "Q22222"),
    gene_name = c("GENE1", "GENE2", "GENE3"),
    cleaned_acc = c("P11111", "P22222", "Q22222"),
    protein_evidence = c(2, 1, 1),
    status = c("reviewed", "reviewed", "reviewed"),
    is_isoform = c("Canonical", "Isoform", "Canonical"),
    isoform_num = c(0, 2, 0),
    seq_length = c(100, 250, 150),
    annotation_score = c(3, 10, 8),
    stringsAsFactors = FALSE
  )

  chosen <- chooseBestProteinAccessionHelper(
    input_tbl = input_tbl,
    acc_detail_tab = acc_detail_tab,
    accessions_column = accession_bundle,
    row_id_column = "uniprot_acc",
    group_id = group_id,
    delim = ":"
  )

  ranked <- rankProteinAccessionHelper(
    input_tbl = input_tbl,
    acc_detail_tab = acc_detail_tab,
    accessions_column = accession_bundle,
    row_id_column = "uniprot_acc",
    group_id = group_id,
    delim = ";"
  )

  fallback <- chooseBestProteinAccessionHelper(
    input_tbl = data.frame(
      group_id = "G3",
      accession_bundle = "R11111:R22222",
      stringsAsFactors = FALSE
    ),
    acc_detail_tab = data.frame(
      uniprot_acc = c("R11111", "R22222"),
      annotation_score = c(0, 0),
      stringsAsFactors = FALSE
    ),
    accessions_column = accession_bundle,
    row_id_column = "uniprot_acc",
    group_id = group_id,
    delim = ":"
  )

  expect_identical(chosen$group_id, c("G1", "G2"))
  expect_true(all(c("num_gene_names", "gene_names", "uniprot_acc", "is_unique") %in% names(chosen)))
  expect_true(grepl("P22222", chosen$uniprot_acc[chosen$group_id == "G1"], fixed = TRUE))
  expect_true(grepl("P11111", chosen$uniprot_acc[chosen$group_id == "G1"], fixed = TRUE))
  expect_identical(chosen$is_unique[chosen$group_id == "G2"], "Unique")
  expect_identical(chosen$uniprot_acc[chosen$group_id == "G2"], "Q22222")

  expect_identical(ranked$group_id, c("G1", "G2"))
  expect_true(all(c("num_gene_names", "gene_names", "uniprot_acc", "is_unique") %in% names(ranked)))
  expect_true(grepl("P22222", ranked$uniprot_acc[ranked$group_id == "G1"], fixed = TRUE))
  expect_true(grepl("Q22222", ranked$uniprot_acc[ranked$group_id == "G2"], fixed = TRUE))

  expect_identical(fallback$group_id, "G3")
  expect_true(all(c("num_gene_names", "gene_names", "uniprot_acc", "is_unique") %in% names(fallback)))
  expect_true(grepl("R11111", fallback$uniprot_acc[[1]], fixed = TRUE))
  expect_identical(fallback$num_gene_names[[1]], 1L)
})
