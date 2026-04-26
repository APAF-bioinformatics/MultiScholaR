# fidelity-coverage-compare: shared
library(testthat)

localBinding <- function(env, name, value, .local_envir = parent.frame()) {
  had_binding <- exists(name, envir = env, inherits = FALSE)
  old_value <- if (had_binding) get(name, envir = env, inherits = FALSE) else NULL
  was_locked <- had_binding && bindingIsLocked(name, env)

  if (was_locked) {
    unlockBinding(name, env)
  }
  assign(name, value, envir = env)
  if (was_locked) {
    lockBinding(name, env)
  }

  withr::defer({
    if (exists(name, envir = env, inherits = FALSE) && bindingIsLocked(name, env)) {
      unlockBinding(name, env)
    }
    if (had_binding) {
      assign(name, old_value, envir = env)
    } else if (exists(name, envir = env, inherits = FALSE)) {
      rm(list = name, envir = env)
    }
    if (was_locked && exists(name, envir = env, inherits = FALSE)) {
      lockBinding(name, env)
    }
  }, envir = .local_envir)
}

test_that("protein QC support helpers preserve data summaries and replicate transforms", {
  peptide_tbl <- data.frame(
    Protein.Ids = c("P1", "P1", "P1", "P2", "P2"),
    Run = c("S1", "S1", "S2", "S2", "S3"),
    Stripped.Sequence = c("PEP1", "PEP2", "PEP1", "PEP3", "PEP3"),
    Peptide.Normalised = c(10, NA, 20, 30, NA),
    Log2.Protein.Imputed = c(5, 5, 7, 8, NA),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  metadata_tbl <- data.frame(
    Run = c("S1", "S2", "S3"),
    collaborator_patient_id = c("Patient1", "Patient1", "Patient2"),
    condition = c("A", "A", "B"),
    batch = c("B1", "B2", "B1"),
    stringsAsFactors = FALSE
  )

  averaged <- avgReplicateProteinIntensity(
    peptide_tbl,
    metadata_tbl,
    protein_id_column = Protein.Ids,
    input_table_sample_id_column = Run,
    sample_id_tbl_sample_id_column = Run,
    replicate_group_column = collaborator_patient_id,
    quantity_column = Log2.Protein.Imputed,
    avg_quantity_column = Avg.Log2.Protein.Imputed
  )

  peptide_missing <- calculatePercentMissingPeptidePerReplicate(
    peptide_tbl,
    metadata_tbl,
    protein_id_column = Protein.Ids,
    intensity_column = Peptide.Normalised,
    replicate_id_column = Run,
    peptide_sequence_column = Stripped.Sequence
  )

  protein_missing <- calculatePercentMissingProteinPerReplicate(
    peptide_tbl,
    metadata_tbl,
    protein_id_column = Protein.Ids,
    intensity_column = Log2.Protein.Imputed,
    replicate_id_column = Run
  )

  intensity_wide <- data.frame(
    uniprot_acc = c("P1", "P2"),
    A1 = c(1, NA),
    A2 = c(2, NA),
    B1 = c(NA, 3),
    B2 = c(4, NaN),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  experimental_design <- data.frame(
    sample_collaborator_sample_id = c("A1", "A2", "B1", "B2"),
    condition = c("A", "A", "B", "B"),
    batch = c("X", "Y", "X", "Y"),
    stringsAsFactors = FALSE
  )
  localBinding(
    parent.env(environment(calculatePercentMissingPerProtein)),
    "is_missing",
    rlang::sym("is_missing")
  )

  missing_by_protein <- calculatePercentMissingPerProtein(
    intensity_wide_table = intensity_wide,
    experimental_design_table = experimental_design
  )

  localBinding(
    parent.env(environment(calculateMissingValuesPerProteinFishersTest)),
    "plan",
    function(...) invisible(NULL)
  )
  fisher_tbl <- calculateMissingValuesPerProteinFishersTest(
    contrasts_table = data.frame(
      contrasts = "contrast=conditionA-conditionB",
      stringsAsFactors = FALSE
    ),
    missing_value_per_category = missing_by_protein
  )

  rows_to_keep <- getRowsToKeepList(
    input_table = intensity_wide,
    cols = !matches("uniprot_acc"),
    design_matrix = experimental_design,
    sample_id = sample_collaborator_sample_id,
    row_id = uniprot_acc,
    grouping_variable = condition,
    min_num_samples_per_group = 1,
    abundance_threshold = 1.5
  )

  averaged_matrix <- averageValuesFromReplicates(
    input_table = structure(
      c(1, 2, 3, 4, 5, 6),
      dim = c(2, 3),
      dimnames = list(c("P1", "P2"), c("Rep1_A", "Rep2_A", "Rep1_B"))
    ),
    design_matrix = data.frame(
      sample = c("Rep1_A", "Rep2_A", "Rep1_B"),
      avg_id = c("A", "A", "B"),
      stringsAsFactors = FALSE
    ),
    group_pattern = "Rep",
    row_id = "Protein.Ids",
    sample_id = "sample",
    average_replicates_id = "avg_id"
  )

  tech_rep_corr <- proteinTechRepCorrelationHelper(
    design_matrix_tech_rep = data.frame(
      Sample_ID = c("S1", "S2", "Pool"),
      replicates = c("Bio1", "Bio1", "poolQC"),
      tech_rep_num = c(1, 2, 1),
      stringsAsFactors = FALSE
    ),
    data_matrix = structure(
      c(1, 2, 2, 4, 5, 6),
      dim = c(2, 3),
      dimnames = list(c("P1", "P2"), c("S1", "S2", "Pool"))
    )
  )

  expect_equal(
    averaged$Avg.Log2.Protein.Imputed,
    c(17 / 3, 8, NaN),
    tolerance = 1e-8
  )
  expect_equal(peptide_missing$Run, c("S1", "S2"))
  expect_equal(peptide_missing$percent_missing, c(66.66667, 33.33333), tolerance = 1e-5)
  expect_equal(protein_missing$Run, c("S1", "S2"))
  expect_equal(protein_missing$percent_missing, c(0, 0))
  expect_true(all(c("num_missing", "num_present", "perc_missing", "compare_column") %in% names(missing_by_protein)))
  expect_true(all(c("contrast_name", "fisher_test", "fdr") %in% names(fisher_tbl)))
  expect_identical(sort(names(rows_to_keep)), c("A", "B"))
  expect_true("P1" %in% rows_to_keep$A)
  expect_equal(unname(averaged_matrix[, "A"]), c(2, 3))
  expect_equal(unname(averaged_matrix[, "B"]), c(5, 6))
  expect_true(all(is.na(tech_rep_corr$pearson)))
  expect_true(all(is.na(tech_rep_corr$spearman)))
})
