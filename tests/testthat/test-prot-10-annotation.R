# fidelity-coverage-compare: shared
# testthat for Proteomics Annotation
# Phase 4 of Proteomics GUI Test Strategy
library(testthat)

if (!methods::isClass("mockProtAnnotationProteinData")) {
  methods::setClass(
    "mockProtAnnotationProteinData",
    slots = c(
      protein_quant_table = "data.frame",
      protein_id_table = "data.frame",
      protein_id_column = "character"
    )
  )
}

test_that("cleanIsoformNumber removes suffix correctly", {
  expect_equal(cleanIsoformNumber("Q8K4R4-2"), "Q8K4R4")
  expect_equal(cleanIsoformNumber("P12345"), "P12345")
  expect_equal(cleanIsoformNumber("P12345-10"), "P12345")
  # Should not remove if not at the end or not following a dash
  expect_equal(cleanIsoformNumber("P12345-PROT"), "P12345-PROT")
})

test_that("goIdToTerm works with mock GO data", {
  mock_goterms <- c("GO:0008150" = "biological_process", "GO:0003674" = "molecular_function", "GO:0005575" = "cellular_component")
  mock_gotypes <- c("GO:0008150" = "BP", "GO:0003674" = "MF", "GO:0005575" = "CC")
  
  result <- goIdToTerm("GO:0008150; GO:0003674; GO:0005575", goterms = mock_goterms, gotypes = mock_gotypes)
  
  expect_true("go_biological_process" %in% names(result))
  expect_equal(result$go_biological_process, "biological_process")
  expect_equal(result$go_molecular_function, "molecular_function")
  expect_equal(result$go_cellular_compartment, "cellular_component")
})

test_that("matchAnnotations handles mock data", {
  # Create mock da_results_for_enrichment
  da_obj <- new("da_results_for_enrichment")
  da_obj@da_data <- list(
    "T_vs_C" = data.frame(
      uniprot_acc = c("P12345", "P67890"),
      log2FC = c(2, -2),
      fdr_qvalue = c(0.01, 0.01),
      stringsAsFactors = FALSE
    )
  )
  
  # Mock UniProt annotations
  uniprot_ann <- data.frame(
    Entry = c("P12345", "P67890"),
    gene_names = c("GENE1", "GENE2"),
    Protein.names = c("Protein 1", "Protein 2"),
    stringsAsFactors = FALSE
  )
  
  result <- matchAnnotations(
    da_results_s4 = da_obj,
    uniprot_annotations = uniprot_ann,
    protein_id_column = "uniprot_acc"
  )
  
  expect_type(result, "list")
  expect_named(result, c("annotated_da_results", "match_statistics", "unmatched_proteins"))
  expect_equal(result$match_statistics$match_rate, 100)
  expect_true("gene_names" %in% names(result$annotated_da_results))
})

test_that("matchAnnotations keeps unmatched proteins while using version-based fuzzy matches", {
  da_obj <- new("da_results_for_enrichment")
  da_obj@da_data <- list(
    "T_vs_C" = data.frame(
      uniprot_acc = c("P12345.2", "NOT_A_MATCH.1"),
      log2FC = c(1, -1),
      fdr_qvalue = c(0.01, 0.02),
      stringsAsFactors = FALSE
    )
  )

  uniprot_ann <- data.frame(
    Entry = "P12345",
    gene_names = "GENE1",
    stringsAsFactors = FALSE
  )

  result <- matchAnnotations(
    da_results_s4 = da_obj,
    uniprot_annotations = uniprot_ann,
    protein_id_column = "uniprot_acc"
  )

  matched_row <- result$annotated_da_results |>
    dplyr::filter(.data$original_id == "P12345.2")

  expect_equal(nrow(matched_row), 1)
  expect_equal(matched_row$Entry, "P12345")
  expect_equal(result$match_statistics$exact_matches, 0)
  expect_equal(result$match_statistics$fuzzy_matches, 1)
  expect_equal(result$match_statistics$unmatched_proteins, 1)
  expect_equal(result$match_statistics$match_rate, 50)
  expect_identical(result$unmatched_proteins, "NOT_A_MATCH.1")
})

test_that("matchAnnotations reads non-enrichment S4 protein_quant_table IDs", {
  da_obj <- methods::new(
    "mockProtAnnotationProteinData",
    protein_quant_table = data.frame(
      Protein.Ids = c("P12345-2", "P99999;SECOND"),
      Sample_1 = c(10, 20),
      stringsAsFactors = FALSE
    ),
    protein_id_table = data.frame(),
    protein_id_column = "Protein.Ids"
  )

  uniprot_ann <- data.frame(
    Entry = c("P12345", "P99999"),
    gene_names = c("", "GENE9;ALT"),
    stringsAsFactors = FALSE
  )

  result <- matchAnnotations(
    da_results_s4 = da_obj,
    uniprot_annotations = uniprot_ann
  )

  expect_equal(result$match_statistics$total_proteins, 2)
  expect_equal(result$match_statistics$matched_proteins, 2)
  expect_equal(result$match_statistics$exact_matches, 2)
  expect_equal(result$match_statistics$match_rate, 100)
  expect_true(is.na(result$annotated_da_results$gene_name[[1]]))
  expect_identical(result$annotated_da_results$gene_name[[2]], "GENE9")
  expect_identical(result$unmatched_proteins, character(0))
})

test_that("matchAnnotations falls back to protein_id_table and validates inputs", {
  da_obj <- methods::new(
    "mockProtAnnotationProteinData",
    protein_quant_table = data.frame(
      OtherId = c("ignored-1", "ignored-2"),
      Sample_1 = c(10, 20),
      stringsAsFactors = FALSE
    ),
    protein_id_table = data.frame(
      Accession = c("Q11111", "Q22222.1"),
      stringsAsFactors = FALSE
    ),
    protein_id_column = "Accession"
  )

  result <- matchAnnotations(
    da_results_s4 = da_obj,
    uniprot_annotations = data.frame(
      Entry = c("Q11111", "Q22222"),
      stringsAsFactors = FALSE
    ),
    protein_id_column = "Accession"
  )

  expect_equal(result$match_statistics$total_proteins, 2)
  expect_equal(result$match_statistics$matched_proteins, 2)
  expect_equal(result$match_statistics$unmatched_proteins, 0)
  expect_equal(result$match_statistics$fuzzy_matches, 1)
  expect_false("gene_name" %in% names(result$annotated_da_results))
  expect_identical(result$unmatched_proteins, character(0))

  expect_error(
    matchAnnotations(NULL, data.frame(Entry = "P1")),
    "da_results_s4 cannot be NULL",
    fixed = TRUE
  )
  expect_error(
    matchAnnotations(da_obj, data.frame(), protein_id_column = "Accession"),
    "uniprot_annotations cannot be NULL or empty",
    fixed = TRUE
  )
  expect_error(
    matchAnnotations(
      da_obj,
      data.frame(NotEntry = "Q11111"),
      protein_id_column = "Accession"
    ),
    "UniProt ID column 'Entry' not found in annotations",
    fixed = TRUE
  )
  expect_error(
    matchAnnotations(
      da_obj,
      data.frame(Entry = "Q11111"),
      protein_id_column = "MissingId"
    ),
    "Protein ID column 'MissingId' not found",
    fixed = TRUE
  )
})

# APAF Bioinformatics | test-prot-10-annotation.R | Approved | 2026-03-13
