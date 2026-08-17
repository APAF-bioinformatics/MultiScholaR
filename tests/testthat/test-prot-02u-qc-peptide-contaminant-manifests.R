library(testthat)

syntheticContaminantManifestPath <- function() {
  testthat::test_path(
    "fixtures",
    "contaminants",
    "synthetic-contaminant-manifest.tsv"
  )
}

test_that("synthetic manifest proves exact token and whole-group classification", {
  manifest_path <- syntheticContaminantManifestPath()
  input <- data.frame(
    .source_row_id = 1:7,
    Protein.Group = paste0("GROUP", 1:7),
    Protein.Ids = c(
      "TARGET_A",
      "SYNTH_CONT_A",
      "SYNTH_CONT_A;TARGET_B",
      "sp|SYNTH_CONT_B|SYNTHETIC_B",
      "XSYNTH_CONT_A",
      "CON__TAGGED_CONTAMINANT",
      "TARGET_KERATIN_DESCRIPTION"
    ),
    Description = c(rep("target protein", 6), "common contaminant keratin"),
    stringsAsFactors = FALSE
  )

  result <- classifyPeptideBiologicalExclusions(
    input,
    contaminant_manifest = manifest_path
  )

  expect_identical(result$raw_data, input)
  expect_identical(result$exclusion_ledger$.source_row_id, c(2L, 4L, 6L))
  expect_identical(result$analysis_data$.source_row_id, c(1L, 3L, 5L, 7L))
  expect_true(result$classified_data$contaminant_manifest_partial_match[[3]])
  expect_false(result$classified_data$is_contaminant[[3]])
  expect_false(result$classified_data$is_contaminant[[5]])
  expect_false(result$classified_data$is_contaminant[[7]])
  expect_identical(
    result$classified_data$exclusion_reason[[3]],
    "mixed_group_partial_contaminant_manifest_match"
  )
  expect_identical(
    result$classified_data$classification_source[[2]],
    "user_file:synthetic-contaminant-manifest.tsv"
  )
})

test_that("synthetic manifest provenance is versioned and exactly fingerprinted", {
  manifest <- readPeptideContaminantManifest(syntheticContaminantManifestPath())

  expect_identical(manifest$accessions, c("SYNTH_CONT_A", "SYNTH_CONT_B"))
  expect_identical(manifest$schema_version, "1.0.0")
  expect_identical(manifest$declared_schema_version, "1.0.0")
  expect_identical(manifest$version, "1.0.0")
  expect_identical(
    manifest$source_name,
    "MultiScholaR synthetic contaminant fixture"
  )
  expect_identical(manifest$license, "LGPL-3.0-or-later")
  expect_identical(manifest$contract, "versioned_manifest_v1")
  expect_false(manifest$legacy_adapter)
  expect_identical(manifest$fingerprint_algorithm, "sha256")
  expect_identical(
    manifest$fingerprint,
    "ab222a67227f8367eb02aea26ef97c0bd494798ee42a8f3c2ce9cc5fbf65952a"
  )
  expect_identical(manifest$checksum, manifest$fingerprint)
  expect_identical(
    manifest$source_file_fingerprint,
    "054c1bac6295a41a196fd94093135bc9745007e4256861f3ea30f59b796e37fb"
  )
})

test_that("portable manifest provenance omits local paths and round trips", {
  manifest_file <- tempfile(fileext = ".tsv")
  writeLines(
    c(
      paste(
        "accession", "manifest_schema_version", "version", "source_name",
        "source_uri", "license",
        sep = "\t"
      ),
      paste(
        "USER_CONT_A", "1.0.0", "2026.08", "APAF user contaminant set",
        "urn:apaf:user-contaminants:2026.08", "user-supplied",
        sep = "\t"
      )
    ),
    manifest_file
  )

  manifest <- readPeptideContaminantManifest(manifest_file)
  portable <- .portablePeptideContaminantManifest(manifest)
  round_trip <- readPeptideContaminantManifest(portable)
  local_root <- normalizePath(tempdir(), winslash = "/", mustWork = TRUE)

  expect_match(manifest$fingerprint, "^[0-9a-f]{64}$")
  expect_match(manifest$source_file_fingerprint, "^[0-9a-f]{64}$")
  expect_false(any(grepl(local_root, unlist(manifest), fixed = TRUE)))
  expect_false(any(grepl(local_root, unlist(portable), fixed = TRUE)))
  expect_identical(round_trip$accessions, manifest$accessions)
  expect_identical(round_trip$fingerprint, manifest$fingerprint)
  expect_identical(round_trip$input_source, manifest$input_source)
  expect_identical(round_trip$validation_status, "valid_versioned_manifest")
})

test_that("valid legacy manifests use an explicit fingerprinted adapter", {
  manifest <- readPeptideContaminantManifest(
    c("sp|P02769|ALBU_BOVIN", "P02769-2")
  )

  expect_identical(manifest$accessions, c("P02769", "P02769-2"))
  expect_true(manifest$legacy_adapter)
  expect_identical(manifest$contract, "legacy_adapter")
  expect_identical(manifest$validation_status, "valid_legacy_adapter")
  expect_identical(manifest$checksum_algorithm, "sha256")
  expect_match(manifest$checksum, "^[0-9a-f]{64}$")
})

test_that("missing and malformed manifests fail without heuristic fallback", {
  missing_path <- file.path(tempdir(), "missing-contaminants.tsv")
  expect_error(
    readPeptideContaminantManifest(missing_path),
    "manifest file does not exist"
  )
  expect_error(
    readPeptideContaminantManifest(data.frame(
      accession = "SYNTH_A",
      manifest_schema_version = "2.0.0",
      version = "1",
      source_name = "synthetic",
      stringsAsFactors = FALSE
    )),
    "Unsupported contaminant manifest schema version"
  )
  expect_error(
    readPeptideContaminantManifest(data.frame(
      accession = c("SYNTH_A", "SYNTH_B"),
      manifest_schema_version = "1.0.0",
      version = c("1", "2"),
      source_name = "synthetic",
      stringsAsFactors = FALSE
    )),
    "metadata 'version' must contain one non-empty value"
  )
  expect_error(
    readPeptideContaminantManifest(data.frame(
      accession = "SYNTH_A",
      manifest_schema_version = "1.0.0",
      version = "1",
      source_name = normalizePath(tempdir(), winslash = "/", mustWork = TRUE),
      stringsAsFactors = FALSE
    )),
    "must not contain an absolute local path"
  )
})
