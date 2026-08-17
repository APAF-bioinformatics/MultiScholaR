library(testthat)

makePeptideExclusionFixture <- function() {
  data.frame(
    .source_row_id = 1:7,
    Run = paste0("S", 1:7),
    Precursor.Id = paste0("prec", 1:7),
    Protein.Group = paste0("GROUP", 1:7),
    Protein.Ids = c(
      "P_TARGET",
      "P_DECOY_COLUMN",
      "CON__P_TAGGED",
      "P02769",
      "P02769;P_TARGET",
      "REV__P_DECOY;P_TARGET",
      "XP02769_SUFFIX"
    ),
    Stripped.Sequence = paste0("PEP", 1:7),
    Modified.Sequence = paste0("PEP", 1:7),
    Decoy = c(0, 1, 0, 0, 0, 0, 0),
    stringsAsFactors = FALSE
  )
}

test_that("DIA-NN decoys and contaminants produce an immutable analysis split", {
  input <- makePeptideExclusionFixture()
  result <- classifyPeptideBiologicalExclusions(
    input,
    contaminant_manifest = c("P02769")
  )

  expect_identical(result$raw_data, input)
  expect_identical(
    result$exclusion_ledger$.source_row_id,
    c(2L, 3L, 4L, 6L)
  )
  expect_identical(
    result$analysis_data$.source_row_id,
    c(1L, 5L, 7L)
  )
  expect_true(result$classified_data$contaminant_manifest_partial_match[[5]])
  expect_false(result$classified_data$is_contaminant[[5]])
  expect_false(result$classified_data$is_contaminant[[7]])
  expect_identical(
    result$classified_data$exclusion_reason[[5]],
    "mixed_group_partial_contaminant_manifest_match"
  )
  expect_identical(result$summary$decoy_rows, 2L)
  expect_identical(result$summary$contaminant_rows, 2L)
  expect_identical(result$summary$excluded_rows, 4L)
})

test_that("contaminant inclusion is explicit without including decoys", {
  result <- classifyPeptideBiologicalExclusions(
    makePeptideExclusionFixture(),
    contaminant_manifest = c("P02769"),
    exclude_contaminants = FALSE
  )

  expect_identical(result$exclusion_ledger$.source_row_id, c(2L, 6L))
  expect_true(all(c(3L, 4L) %in% result$analysis_data$.source_row_id))
  expect_true(result$classified_data$is_contaminant[[3]])
  expect_true(result$classified_data$is_contaminant[[4]])
})

test_that("standard reports without exclusion signals record an explicit no-op", {
  input <- data.frame(
    Protein.Group = "P1",
    Protein.Ids = "P1",
    stringsAsFactors = FALSE
  )
  result <- classifyPeptideBiologicalExclusions(input)

  expect_identical(result$analysis_data$Protein.Group, "P1")
  expect_identical(nrow(result$exclusion_ledger), 0L)
  expect_identical(
    result$summary$classification_status,
    "not_classified_no_reliable_signal"
  )
})

test_that("user manifest files are normalised and checksummed", {
  manifest_file <- tempfile(fileext = ".tsv")
  writeLines(c("accession\tversion", "sp|P02769|ALBU_BOVIN\t2012.01.01"), manifest_file)

  manifest <- readPeptideContaminantManifest(manifest_file)

  expect_identical(manifest$accessions, "P02769")
  expect_identical(manifest$version, "2012.01.01")
  expect_true(nzchar(manifest$checksum))
})

test_that("biological exclusions occur before frozen identification evidence", {
  input <- data.frame(
    Run = c("S1", "S2", "S1", "S2"),
    Precursor.Id = paste0("prec", 1:4),
    Protein.Group = c("P_TARGET", "P_TARGET", "P_CONT", "P_CONT"),
    Protein.Ids = c("P_TARGET", "P_TARGET", "P02769", "P02769"),
    Stripped.Sequence = c("PEP_A", "PEP_B", "CONT_A", "CONT_B"),
    Modified.Sequence = c("PEP_A", "PEP_B", "CONT_A", "CONT_B"),
    Precursor.Charge = 2L,
    Precursor.Quantity = 100,
    Precursor.Normalised = 10,
    Q.Value = 0.001,
    Global.Q.Value = 0.001,
    Global.PG.Q.Value = 0.001,
    Proteotypic = 1,
    stringsAsFactors = FALSE
  )

  filtered <- srlQvalueProteotypicPeptideCleanHelper(
    input,
    input_matrix_column_ids = c(
      "Run", "Precursor.Id", "Protein.Group", "Protein.Ids",
      "Stripped.Sequence", "Modified.Sequence", "Precursor.Charge",
      "Precursor.Quantity", "Precursor.Normalised"
    ),
    protein_id_column = Protein.Group,
    contaminant_manifest = "P02769",
    return_exclusion_result = TRUE
  )

  expect_identical(unique(filtered$data$Protein.Group), "P_TARGET")
  expect_identical(unique(filtered$data$identification_peptide_count), 2L)
  expect_identical(unique(filtered$data$identification_peptidoform_count), 2L)
  expect_identical(unique(filtered$exclusion_ledger$Protein.Group), "P_CONT")
  expect_identical(filtered$exclusion_summary$excluded_rows, 2L)
})

test_that("DIA-NN S4 filtering stores the exclusion ledger and manifest provenance", {
  manifest_file <- tempfile(fileext = ".tsv")
  writeLines(
    c(
      paste(
        "accession", "manifest_schema_version", "version", "source_name",
        sep = "\t"
      ),
      paste("P02769", "1.0.0", "test-v1", "synthetic test manifest", sep = "\t")
    ),
    manifest_file
  )
  input <- data.frame(
    Run = c("S1", "S2"),
    Precursor.Id = c("prec1", "prec2"),
    Protein.Group = c("P_TARGET", "P_CONT"),
    Protein.Ids = c("P_TARGET", "P02769"),
    Stripped.Sequence = c("PEP_A", "CONT_A"),
    Modified.Sequence = c("PEP_A", "CONT_A"),
    Precursor.Charge = 2L,
    Precursor.Quantity = 100,
    Precursor.Normalised = 10,
    Q.Value = 0.001,
    Global.Q.Value = 0.001,
    Global.PG.Q.Value = 0.001,
    Proteotypic = 1,
    stringsAsFactors = FALSE
  )
  design <- data.frame(
    Run = c("S1", "S2"),
    group = "A",
    replicates = c("R1", "R2"),
    stringsAsFactors = FALSE
  )
  object <- PeptideQuantitativeDataDiann(input, design, args = list())

  filtered <- srlQvalueProteotypicPeptideClean(
    object,
    input_matrix_column_ids = c(
      "Run", "Precursor.Id", "Protein.Group", "Protein.Ids",
      "Stripped.Sequence", "Modified.Sequence", "Precursor.Charge",
      "Precursor.Quantity", "Precursor.Normalised"
    ),
    contaminant_manifest = manifest_file
  )
  metadata <- filtered@args$srlQvalueProteotypicPeptideClean
  local_root <- normalizePath(tempdir(), winslash = "/", mustWork = TRUE)

  expect_identical(filtered@peptide_data$Protein.Group, "P_TARGET")
  expect_identical(metadata$biological_exclusion_summary$excluded_rows, 1L)
  expect_identical(metadata$biological_exclusion_ledger$Protein.Group, "P_CONT")
  expect_identical(metadata$contaminant_manifest_provenance$accessions, "P02769")
  expect_identical(
    metadata$contaminant_manifest_provenance$validation_status,
    "valid_versioned_manifest"
  )
  expect_s3_class(metadata$contaminant_manifest, "data.frame")
  expect_false(any(grepl(local_root, unlist(metadata$contaminant_manifest), fixed = TRUE)))
  expect_false(any(grepl(
    local_root,
    unlist(metadata$contaminant_manifest_provenance),
    fixed = TRUE
  )))
})

test_that("an explicit contaminant column takes precedence over missing manifests", {
  input <- data.frame(
    Protein.Group = c("P1", "P2"),
    Protein.Ids = c("P1", "P2"),
    Potential.contaminant = c("", "+"),
    stringsAsFactors = FALSE
  )
  result <- classifyPeptideBiologicalExclusions(input)

  expect_identical(result$analysis_data$Protein.Group, "P1")
  expect_identical(
    result$exclusion_ledger$exclusion_reason,
    "explicit_contaminant_column"
  )
})
