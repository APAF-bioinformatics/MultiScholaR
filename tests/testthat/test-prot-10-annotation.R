# testthat for Proteomics Annotation
# Phase 4 of Proteomics GUI Test Strategy

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

# APAF Bioinformatics | test-prot-10-annotation.R | Approved | 2026-03-13
