# fidelity-coverage-compare: shared
# testthat for Proteomics Rollup
# Phase 4 of Proteomics GUI Test Strategy

test_that("rolled-up protein snapshot is valid", {
  cp_file <- test_path("..", "testdata", "sepsis", "proteomics", "cp03_rolled_up_protein.rds")
  
  if (file.exists(cp_file)) {
    first_line <- tryCatch(readLines(cp_file, n = 1, warn = FALSE), error = function(e) "")
    if (length(first_line) > 0 && identical(first_line[[1]], "version https://git-lfs.github.com/spec/v1")) {
      skip("Snapshot cp03 is a Git LFS pointer and the binary artifact is not present")
    }

    obj <- readRDS(cp_file)
    expect_s4_class(obj, "ProteinQuantitativeData")
    expect_true(nrow(obj@protein_quant_table) > 0)
  } else {
    skip("Snapshot cp03 not found")
  }
})

test_that("calcPeptidesPerProtein works with mock data", {
  pept_data <- data.frame(
    Protein.Ids = c("P1", "P1", "P2"),
    Stripped.Sequence = c("PEP1", "PEP2", "PEP3"),
    Run = rep("S1", 3),
    Q.Value = rep(0.01, 3),
    Precursor.Quantity = rep(100, 3),
    stringsAsFactors = FALSE
  )
  
  # Create mock S4 with peptide_data and design_matrix
  obj <- new("PeptideQuantitativeData",
             peptide_data = pept_data,
             design_matrix = data.frame(Run = "S1", group = "G1", stringsAsFactors = FALSE),
             sample_id = "Run",
             group_id = "group",
             protein_id_column = "Protein.Ids",
             peptide_sequence_column = "Stripped.Sequence",
             q_value_column = "Q.Value",
             raw_quantity_column = "Precursor.Quantity")
  
  result <- calcPeptidesPerProtein(obj)
  
  expect_equal(nrow(result), 2)
  expect_equal(result$n_peptides[result$Protein.Ids == "P1"], 2)
  expect_equal(result$n_peptides[result$Protein.Ids == "P2"], 1)
})

test_that("rollUpPrecursorToPeptideHelper sums quantities correctly", {
  # Mock input table
  input_table <- data.frame(
    Protein.Ids = rep("P1", 4),
    Stripped.Sequence = rep("PEP1", 4),
    Modified.Sequence = c("PEP1_mod1", "PEP1_mod1", "PEP1_mod2", "PEP1_mod2"),
    Precursor.Id = c("PEP1_mod1_2", "PEP1_mod1_2", "PEP1_mod2_3", "PEP1_mod2_3"),
    Precursor.Charge = c(2L, 2L, 3L, 3L),
    Run = c("S1", "S2", "S1", "S2"),
    Precursor.Quantity = c(10, 20, 30, 40),
    Precursor.Normalised = c(10, 20, 30, 40),
    identification_peptide_count = 1L,
    identification_peptidoform_count = 2L,
    stringsAsFactors = FALSE
  )
  
  # Run helper with core_utilisation = 1 (No parallel)
  result <- rollUpPrecursorToPeptideHelper(
    input_table = input_table,
    sample_id_column = Run,
    protein_id_column = Protein.Ids,
    peptide_sequence_column = Stripped.Sequence,
    modified_peptide_sequence_column = Modified.Sequence,
    precursor_quantity_column = Precursor.Quantity,
    precursor_normalised_column = Precursor.Normalised,
    core_utilisation = 1
  )
  
  # result should be a data frame with quantities summed per sample and peptide
  expect_s3_class(result, "data.frame")
  expect_equal(result$Peptide.RawQuantity[result$Run == "S1"], 40)
  expect_equal(result$Peptide.RawQuantity[result$Run == "S2"], 60)
  expect_identical(unique(result$identification_peptide_count), 1L)
  expect_identical(unique(result$identification_peptidoform_count), 2L)
  expect_false("peptidoform_ids" %in% names(result))
})

makeRollupIntegrityFixture <- function() {
  data.frame(
    Run = rep("S1", 6),
    Protein.Group = rep("PG1", 6),
    Protein.Ids = c("P1;P2", "P1", "P2", "P1", "P2", "P1;P2"),
    Stripped.Sequence = rep("PEPTIDE", 6),
    Modified.Sequence = rep(c("PEPTIDE", "PEP[MOD1]TIDE", "PEP[MOD2]TIDE"), each = 2),
    Precursor.Charge = rep(c(2L, 3L), 3),
    Precursor.Id = paste0("precursor_", seq_len(6)),
    Precursor.Quantity = c(1, 2, 4, 8, 16, 32),
    Precursor.Normalised = c(2, 4, 8, 16, 32, 64),
    identification_peptide_count = rep(1L, 6),
    identification_peptidoform_count = rep(3L, 6),
    stringsAsFactors = FALSE
  )
}

test_that("rollup intentionally merges modified forms and charge states on linear scale", {
  result <- rollUpPrecursorToPeptideHelper(
    makeRollupIntegrityFixture(),
    protein_id_column = Protein.Group,
    core_utilisation = NA,
    return_rollup_result = TRUE
  )

  expect_equal(nrow(result$data), 1L)
  expect_equal(result$data$Peptide.RawQuantity, 63)
  expect_equal(result$data$Peptide.Normalised, 126)
  expect_equal(result$data$peptidoform_count, 3L)
  expect_identical(result$data$Protein.Group, "PG1")
  expect_identical(result$data$Protein.Ids, "P1;P2")
  expect_identical(result$data$identification_peptide_count, 1L)
  expect_identical(result$data$identification_peptidoform_count, 3L)
  expect_identical(
    result$summary$aggregation_rule,
    "linear_sum_modified_unmodified_and_charge_states"
  )
  expect_identical(result$summary$active_protein_key, "Protein.Group")
  expect_equal(result$summary$input_precursor_rows, 6L)
  expect_equal(result$summary$unique_precursor_identities, 6L)
  expect_equal(result$summary$output_stripped_peptide_identities, 1L)
})

test_that("duplicate precursor identities fail with bounded key evidence", {
  input <- makeRollupIntegrityFixture()
  input$Precursor.Id[2] <- input$Precursor.Id[1]

  expect_error(
    rollUpPrecursorToPeptideHelper(
      input,
      protein_id_column = Protein.Group,
      core_utilisation = NA
    ),
    "duplicate precursor identities.*Run=S1, Protein.Group=PG1, Precursor.Id=precursor_1"
  )
})

test_that("partial missing components sum observed values and all-missing stays NA", {
  input <- makeRollupIntegrityFixture()[1:2, ]
  input$Precursor.Quantity <- c(NA_real_, 5)
  input$Precursor.Normalised <- c(NA_real_, NA_real_)

  result <- rollUpPrecursorToPeptideHelper(
    input,
    protein_id_column = Protein.Group,
    core_utilisation = NA,
    return_rollup_result = TRUE
  )

  expect_equal(result$data$Peptide.RawQuantity, 5)
  expect_true(is.na(result$data$Peptide.Normalised))
  expect_equal(result$summary$raw_partial_missing_outputs, 1L)
  expect_equal(result$summary$raw_all_missing_outputs, 0L)
  expect_equal(result$summary$normalised_partial_missing_outputs, 0L)
  expect_equal(result$summary$normalised_all_missing_outputs, 1L)
})

test_that("negative, infinite, NaN, and logged abundances are rejected", {
  input <- makeRollupIntegrityFixture()
  negative <- input
  negative$Precursor.Quantity[1] <- -1
  infinite <- input
  infinite$Precursor.Normalised[2] <- Inf
  nan_input <- input
  nan_input$Precursor.Quantity[3] <- NaN

  expect_error(
    rollUpPrecursorToPeptideHelper(negative, protein_id_column = Protein.Group, core_utilisation = NA),
    "negative or non-finite.*row\\(s\\) 1"
  )
  expect_error(
    rollUpPrecursorToPeptideHelper(infinite, protein_id_column = Protein.Group, core_utilisation = NA),
    "negative or non-finite.*row\\(s\\) 2"
  )
  expect_error(
    rollUpPrecursorToPeptideHelper(nan_input, protein_id_column = Protein.Group, core_utilisation = NA),
    "negative or non-finite.*row\\(s\\) 3"
  )
  expect_error(
    rollUpPrecursorToPeptideHelper(
      input,
      protein_id_column = Protein.Group,
      core_utilisation = NA,
      is_logged_data = TRUE
    ),
    "logged abundance cannot be summed"
  )
})

test_that("legacy inputs derive identity explicitly and still reject duplicates", {
  input <- makeRollupIntegrityFixture()
  input$Precursor.Id <- NULL

  expect_warning(
    result <- rollUpPrecursorToPeptideHelper(
      input,
      protein_id_column = Protein.Group,
      core_utilisation = NA,
      return_rollup_result = TRUE
    ),
    "declared precursor identity"
  )
  expect_identical(
    result$summary$precursor_identity_mode,
    "legacy_derived_modified_sequence_and_charge"
  )

  no_evidence <- input
  no_evidence$identification_peptide_count <- NULL
  no_evidence$identification_peptidoform_count <- NULL
  expect_warning(
    no_evidence_result <- rollUpPrecursorToPeptideHelper(
      no_evidence,
      protein_id_column = Protein.Group,
      core_utilisation = NA,
      return_rollup_result = TRUE
    ),
    "declared precursor identity"
  )
  expect_identical(
    no_evidence_result$summary$frozen_identification_evidence,
    "absent"
  )

  duplicate <- rbind(input, input[1, ])
  expect_error(
    suppressWarnings(rollUpPrecursorToPeptideHelper(
      duplicate,
      protein_id_column = Protein.Group,
      core_utilisation = NA
    )),
    "duplicate precursor identities"
  )
})

test_that("sequential and multidplyr rollup paths are identical", {
  skip_if_not_installed("multidplyr")
  input <- rbind(
    makeRollupIntegrityFixture(),
    transform(
      makeRollupIntegrityFixture(),
      Run = "S2",
      Protein.Group = "PG2",
      Protein.Ids = "P3",
      Precursor.Id = paste0("second_", seq_len(6))
    )
  )
  sequential <- rollUpPrecursorToPeptideHelper(
    input,
    protein_id_column = Protein.Group,
    core_utilisation = NA
  )

  cluster <- multidplyr::new_cluster(2)
  withr::defer(
    invisible(lapply(unclass(cluster), function(worker) worker$kill())),
    envir = environment()
  )
  parallel <- rollUpPrecursorToPeptideHelper(
    input,
    protein_id_column = Protein.Group,
    core_utilisation = cluster
  )

  expect_identical(parallel, sequential)
})

# APAF Bioinformatics | test-prot-03-rollup.R | Approved | 2026-03-13
