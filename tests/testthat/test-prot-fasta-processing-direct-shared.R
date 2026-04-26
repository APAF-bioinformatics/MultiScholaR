# fidelity-coverage-compare: shared
library(testthat)

processFastaFile <- get(
  "processFastaFile",
  envir = asNamespace("MultiScholaR"),
  inherits = FALSE
)

test_that("processFastaFile preserves standard UniProt FASTA parsing and metadata capture", {
  skip_if_not_installed("seqinr")
  skip_if_not_installed("vroom")

  fixture_dir <- tempfile("prot-fasta-standard-")
  dir.create(fixture_dir, recursive = TRUE)
  withr::defer(unlink(fixture_dir, recursive = TRUE, force = TRUE))
  withr::local_dir(fixture_dir)

  fasta_path <- file.path(fixture_dir, "standard.fasta")
  meta_path <- file.path(fixture_dir, "standard_meta.rds")

  writeLines(
    c(
      ">sp|P11111|PROTA_HUMAN OS=Homo sapiens OX=9606 GN=GENEA PE=1 SV=1",
      "MPEPTIDESEQ",
      ">tr|Q22222|PROTB_HUMAN OS=Homo sapiens OX=9606 GN=GENEB PE=3 SV=2",
      "MSEQUENCEAA"
    ),
    fasta_path
  )

  result <- processFastaFile(
    fasta_file_path = fasta_path,
    fasta_meta_file = meta_path,
    organism_name = "Homo sapiens"
  )

  expect_named(result, c("aa_seq_tbl_final", "fasta_metadata"))
  expect_identical(result$fasta_metadata$fasta_format, "standard_uniprot")
  expect_identical(result$fasta_metadata$num_sequences, 2L)
  expect_true(result$fasta_metadata$has_protein_evidence)
  expect_true(result$fasta_metadata$has_gene_names)
  expect_true(result$fasta_metadata$has_status_info)
  expect_identical(result$aa_seq_tbl_final$accession, c("P11111", "Q22222"))
  expect_identical(result$aa_seq_tbl_final$cleaned_acc, c("P11111", "Q22222"))
  expect_identical(result$aa_seq_tbl_final$gene_name, c("GENEA", "GENEB"))
  expect_identical(result$aa_seq_tbl_final$status, c("unreviewed", "unreviewed"))
  expect_equal(result$aa_seq_tbl_final$protein_evidence, c(1L, 3L))
  expect_true(file.exists(meta_path))

  saved <- readRDS(meta_path)
  expect_identical(saved$accession, c("P11111", "Q22222"))
})

test_that("processFastaFile preserves non-standard FASTA parsing and external mapping reconciliation", {
  skip_if_not_installed("seqinr")
  skip_if_not_installed("vroom")

  fixture_dir <- tempfile("prot-fasta-nonstandard-")
  dir.create(fixture_dir, recursive = TRUE)
  withr::defer(unlink(fixture_dir, recursive = TRUE, force = TRUE))
  withr::local_dir(fixture_dir)

  fasta_path <- file.path(fixture_dir, "nonstandard.fasta")
  meta_path <- file.path(fixture_dir, "nonstandard_meta.rds")

  writeLines(
    c(
      ">ref|ABC123 [locus_tag=LT1] [protein=Protein alpha] [protein_id=WP_000001]",
      "MALPHASEQ",
      ">ref|XYZ789 [locus_tag=LT2] [protein=Protein beta] [protein_id=WP_000002]",
      "MBETASEQ"
    ),
    fasta_path
  )

  uniprot_search_results <- data.frame(
    Organism = "Bacillus subtilis",
    ncbi_refseq = "WP_000001",
    uniprot_id = "UPI000001",
    stringsAsFactors = FALSE
  )
  uniparc_search_results <- data.frame(
    ncbi_refseq = c("WP_000001", "WP_000002"),
    uniprot_id = c("UPARC_A", "UPARC_B"),
    stringsAsFactors = FALSE
  )

  result <- processFastaFile(
    fasta_file_path = fasta_path,
    uniprot_search_results = uniprot_search_results,
    uniparc_search_results = uniparc_search_results,
    fasta_meta_file = meta_path,
    organism_name = "Bacillus subtilis"
  )

  expect_named(result, c("aa_seq_tbl_final", "fasta_metadata"))
  expect_identical(result$fasta_metadata$fasta_format, "non_standard")
  expect_identical(result$fasta_metadata$num_sequences, 2L)
  expect_false(result$fasta_metadata$has_protein_evidence)
  expect_false(result$fasta_metadata$has_gene_names)
  expect_identical(result$aa_seq_tbl_final$accession, c("ABC123", "XYZ789"))
  expect_identical(result$aa_seq_tbl_final$protein_id, c("LT1", "LT2"))
  expect_identical(result$aa_seq_tbl_final$protein, c("Protein alpha", "Protein beta"))
  expect_identical(result$aa_seq_tbl_final$database_id, c("UPI000001", "UPARC_B"))
  expect_true(file.exists(file.path(fixture_dir, "aa_seq_tbl.tsv")))
  expect_true(file.exists(meta_path))

  saved <- readRDS(meta_path)
  expect_identical(saved$database_id, c("UPI000001", "UPARC_B"))
})
