# fidelity-coverage-compare: shared
library(testthat)

makeFunctionWithOverrides <- function(fun, replacements) {
  fun_override <- fun
  environment(fun_override) <- list2env(replacements, parent = environment(fun))
  fun_override
}

test_that("peptide replicate filter helpers preserve direct and partitioned branches", {
  input_table <- data.frame(
    Run = c("S1", "S2", "S3", "S4", "S1", "S2", "S3"),
    Protein.Ids = c("P1", "P1", "P1", "P1", "P2", "P2", "P3"),
    Stripped.Sequence = c("pep1", "pep1", "pep1", "pep2", "pep3", "pep3", "pep4"),
    peptidoform_count = c(2, 2, 1, 1, 2, 2, 1),
    Q.Value = c(0.001, 0.002, 0.001, 0.020, 0.003, 0.004, 0.001),
    Global.Q.Value = c(0.002, 0.002, 0.002, 0.020, 0.004, 0.004, 0.001),
    Proteotypic = c(1, 1, 0, 1, 1, 1, 1),
    Precursor.Id = paste0("prec", 1:7),
    Modified.Sequence = paste0("mod", 1:7),
    Precursor.Charge = c(2, 2, 2, 3, 2, 2, 2),
    Precursor.Quantity = c(100, 110, 90, 80, 120, 130, 140),
    Precursor.Normalised = c(10, 11, 9, 8, 12, 13, 14),
    stringsAsFactors = FALSE
  )

  samples_id_tbl <- data.frame(
    ms_filename = c("S1", "S2", "S3", "S4"),
    general_sample_info = c("PatientA", "PatientA", "PatientB", "PatientB"),
    stringsAsFactors = FALSE
  )

  remove_parallel <- makeFunctionWithOverrides(
    removePeptidesWithOnlyOneReplicateHelper,
    list(
      partition = function(x, cores) x,
      collect = function(x) x
    )
  )
  protein_parallel <- makeFunctionWithOverrides(
    filterMinNumPeptidesPerProteinHelper,
    list(
      partition = function(x, cores) x,
      collect = function(x) x
    )
  )
  sample_parallel <- makeFunctionWithOverrides(
    filterMinNumPeptidesPerSampleHelper,
    list(
      partition = function(x, cores) x,
      collect = function(x) x
    )
  )

  direct_removed <- removePeptidesWithOnlyOneReplicateHelper(
    input_table = input_table,
    samples_id_tbl = samples_id_tbl,
    input_table_sample_id_column = Run,
    sample_id_tbl_sample_id_column = ms_filename,
    replicate_group_column = general_sample_info,
    protein_id_column = Protein.Ids,
    peptide_sequence_column = Stripped.Sequence,
    core_utilisation = NA
  )
  parallel_removed <- remove_parallel(
    input_table = input_table,
    samples_id_tbl = samples_id_tbl,
    input_table_sample_id_column = Run,
    sample_id_tbl_sample_id_column = ms_filename,
    replicate_group_column = general_sample_info,
    protein_id_column = Protein.Ids,
    peptide_sequence_column = Stripped.Sequence,
    core_utilisation = 2
  )

  expect_identical(sort(unique(direct_removed$Protein.Ids)), c("P1", "P2"))
  expect_equal(nrow(direct_removed), nrow(parallel_removed))

  direct_protein <- filterMinNumPeptidesPerProteinHelper(
    input_table = input_table,
    num_peptides_per_protein_thresh = 2,
    num_peptidoforms_per_protein_thresh = 2,
    protein_id_column = Protein.Ids,
    core_utilisation = NA
  )
  parallel_protein <- protein_parallel(
    input_table = input_table,
    num_peptides_per_protein_thresh = 2,
    num_peptidoforms_per_protein_thresh = 2,
    protein_id_column = Protein.Ids,
    core_utilisation = 2
  )

  expect_setequal(unique(direct_protein$Protein.Ids), c("P1", "P2"))
  expect_equal(nrow(direct_protein), nrow(parallel_protein))
  expect_error(
    filterMinNumPeptidesPerProteinHelper(
      input_table = input_table,
      num_peptides_per_protein_thresh = NA_real_,
      num_peptidoforms_per_protein_thresh = 2,
      protein_id_column = Protein.Ids,
      core_utilisation = NA
    ),
    "must be provided",
    fixed = TRUE
  )

  direct_sample <- filterMinNumPeptidesPerSampleHelper(
    input_table = input_table,
    peptides_per_sample_cutoff = 2,
    sample_id_column = Run,
    core_utilisation = NA,
    inclusion_list = "S4"
  )
  parallel_sample <- sample_parallel(
    input_table = input_table,
    peptides_per_sample_cutoff = 2,
    sample_id_column = Run,
    core_utilisation = 2,
    inclusion_list = "S4"
  )

  expect_setequal(unique(direct_sample$Run), c("S1", "S2", "S3", "S4"))
  expect_equal(nrow(direct_sample), nrow(parallel_sample))
})

test_that("peptide q-value cleanup helper preserves validation and filtered output", {
  input_table <- data.frame(
    Run = c("S1", "S2", "S3"),
    Precursor.Id = c("prec1", "prec2", "prec3"),
    Protein.Ids = c("P1", "P1", "P2"),
    Stripped.Sequence = c("pep1", "pep2", "pep3"),
    Modified.Sequence = c("mod1", "mod2", "mod3"),
    Precursor.Charge = c(2, 2, 3),
    Precursor.Quantity = c(100, 120, 90),
    Precursor.Normalised = c(10, 12, 9),
    Q.Value = c(0.001, 0.020, 0.003),
    Global.Q.Value = c(0.002, 0.030, 0.004),
    Proteotypic = c(1, 1, 0),
    stringsAsFactors = FALSE
  )

  cleaned <- srlQvalueProteotypicPeptideCleanHelper(
    input_table = input_table,
    qvalue_threshold = 0.01,
    global_qvalue_threshold = 0.01,
    choose_only_proteotypic_peptide = 1
  )

  expect_equal(nrow(cleaned), 1)
  expect_identical(cleaned$Run, "S1")

  expect_error(
    srlQvalueProteotypicPeptideCleanHelper(
      input_table = input_table[, setdiff(names(input_table), "Precursor.Id")],
      input_matrix_column_ids = c("Run", "Precursor.Id", "Protein.Ids")
    ),
    "Required output columns not found",
    fixed = TRUE
  )

  expect_error(
    srlQvalueProteotypicPeptideCleanHelper(
      input_table = input_table[, setdiff(names(input_table), "Q.Value")],
      input_matrix_column_ids = c("Run", "Protein.Ids")
    ),
    "Required filter columns not found",
    fixed = TRUE
  )
})
