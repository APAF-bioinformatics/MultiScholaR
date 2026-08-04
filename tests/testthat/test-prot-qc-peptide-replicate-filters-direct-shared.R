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
    Global.PG.Q.Value = c(0.002, 0.002, 0.002, 0.020, 0.004, 0.004, 0.001),
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

  expect_identical(unique(direct_protein$Protein.Ids), "P1")
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
    protein_id_column = Protein.Ids,
    peptide_sequence_column = Stripped.Sequence,
    core_utilisation = NA,
    inclusion_list = "S4"
  )
  parallel_sample <- sample_parallel(
    input_table = input_table,
    peptides_per_sample_cutoff = 2,
    sample_id_column = Run,
    protein_id_column = Protein.Ids,
    peptide_sequence_column = Stripped.Sequence,
    core_utilisation = 2,
    inclusion_list = "S4"
  )

  expect_setequal(unique(direct_sample$Run), c("S1", "S2", "S3", "S4"))
  expect_equal(nrow(direct_sample), nrow(parallel_sample))
})

makeDistinctIdentityFilterFixture <- function() {
  design <- data.frame(
    Run = c("A1", "A2", "B1", "B2", "B3"),
    replicate_group = c("A", "A", "B", "B", "B"),
    stringsAsFactors = FALSE
  )
  input <- data.frame(
    Run = c("A1", "A2", "B1", "A1", "B1", "A1", "A2"),
    Protein.Group = c("PG1", "PG1", "PG1", "PG2", "PG2", "PG3", "PG3"),
    Protein.Ids = c("P1", "P1", "P1", "P2", "P2", "P3", "P3"),
    Stripped.Sequence = c(
      "SHARED", "SHARED", "SHARED", "SHARED", "SHARED", "EQUAL", "EQUAL"
    ),
    Peptide.Normalised = c(10, 10, 10, 20, 20, 7, 7),
    stringsAsFactors = FALSE
  )
  list(input = input, design = design)
}

test_that("replicate support counts distinct runs and retains singleton rows in other groups", {
  fixture <- makeDistinctIdentityFilterFixture()
  result <- removePeptidesWithOnlyOneReplicateHelper(
    input_table = fixture$input,
    samples_id_tbl = fixture$design,
    input_table_sample_id_column = "Run",
    sample_id_tbl_sample_id_column = "Run",
    replicate_group_column = "replicate_group",
    protein_id_column = "Protein.Group",
    peptide_sequence_column = "Stripped.Sequence",
    core_utilisation = NA,
    return_filter_result = TRUE
  )

  expect_setequal(unique(result$data$Protein.Group), c("PG1", "PG3"))
  expect_true(any(result$data$Protein.Group == "PG1" & result$data$Run == "B1"))
  pg1_support <- result$support_table[
    result$support_table$Protein.Group == "PG1",
  ]
  expect_equal(pg1_support$distinct_run_count, c(2L, 1L))
  expect_identical(pg1_support$group_supports_feature, c(TRUE, FALSE))
  expect_identical(
    result$summary$retention_policy,
    "supported_in_any_group"
  )
  expect_identical(result$summary$identity_definition, "Protein.Group + Stripped.Sequence")
  expect_identical(result$removal_ledger$Protein.Group, "PG2")

  duplicated <- removePeptidesWithOnlyOneReplicateHelper(
    input_table = rbind(fixture$input, fixture$input),
    samples_id_tbl = fixture$design,
    input_table_sample_id_column = "Run",
    sample_id_tbl_sample_id_column = "Run",
    replicate_group_column = "replicate_group",
    protein_id_column = "Protein.Group",
    peptide_sequence_column = "Stripped.Sequence",
    core_utilisation = NA,
    return_filter_result = TRUE
  )
  expect_identical(duplicated$support_table, result$support_table)
  expect_identical(duplicated$removal_ledger, result$removal_ledger)
  expect_setequal(unique(duplicated$data$Protein.Group), unique(result$data$Protein.Group))
})

test_that("sample support counts distinct Protein.Group and stripped-peptide pairs", {
  input <- data.frame(
    Run = c("S1", "S1", "S1", "S2"),
    Protein.Group = c("PG1", "PG1", "PG2", "PG1"),
    Protein.Ids = c("P1", "P1", "P2", "P1"),
    Stripped.Sequence = "SHARED",
    stringsAsFactors = FALSE
  )
  design <- data.frame(Run = c("S1", "S2"), stringsAsFactors = FALSE)
  result <- filterMinNumPeptidesPerSampleHelper(
    input_table = input,
    peptides_per_sample_cutoff = 2,
    sample_id_column = "Run",
    core_utilisation = NA,
    protein_id_column = "Protein.Group",
    peptide_sequence_column = "Stripped.Sequence",
    design_matrix = design,
    design_sample_id_column = "Run",
    return_filter_result = TRUE
  )

  expect_identical(unique(result$data$Run), "S1")
  expect_equal(result$support_table$distinct_peptide_identity_count, c(2L, 1L))
  expect_identical(
    result$summary$count_definition,
    "distinct_protein_group_stripped_peptide_pairs_per_run"
  )
  duplicated <- filterMinNumPeptidesPerSampleHelper(
    input_table = rbind(input, input),
    peptides_per_sample_cutoff = 2,
    sample_id_column = "Run",
    core_utilisation = NA,
    protein_id_column = "Protein.Group",
    peptide_sequence_column = "Stripped.Sequence",
    design_matrix = design,
    design_sample_id_column = "Run",
    return_filter_result = TRUE
  )
  expect_identical(duplicated$support_table, result$support_table)

  included <- filterMinNumPeptidesPerSampleHelper(
    input_table = input,
    peptides_per_sample_cutoff = 2,
    sample_id_column = "Run",
    core_utilisation = NA,
    inclusion_list = c("S2", "NOT_IN_DATA"),
    protein_id_column = "Protein.Group",
    peptide_sequence_column = "Stripped.Sequence",
    design_matrix = design,
    design_sample_id_column = "Run",
    return_filter_result = TRUE
  )
  expect_setequal(unique(included$data$Run), c("S1", "S2"))
  expect_false("NOT_IN_DATA" %in% included$support_table$Run)
  expect_identical(
    included$summary$inclusion_entries_absent_from_data,
    "NOT_IN_DATA"
  )
})

test_that("sample and replicate filters reject broken design mappings", {
  fixture <- makeDistinctIdentityFilterFixture()
  duplicate_design <- rbind(fixture$design, fixture$design[1, ])
  expect_error(
    removePeptidesWithOnlyOneReplicateHelper(
      fixture$input,
      duplicate_design,
      input_table_sample_id_column = "Run",
      sample_id_tbl_sample_id_column = "Run",
      replicate_group_column = "replicate_group",
      protein_id_column = "Protein.Group",
      peptide_sequence_column = "Stripped.Sequence",
      core_utilisation = NA
    ),
    "duplicate run ID\\(s\\): A1"
  )

  unmatched <- fixture$input
  unmatched$Run[1] <- "UNMAPPED"
  expect_error(
    removePeptidesWithOnlyOneReplicateHelper(
      unmatched,
      fixture$design,
      input_table_sample_id_column = "Run",
      sample_id_tbl_sample_id_column = "Run",
      replicate_group_column = "replicate_group",
      protein_id_column = "Protein.Group",
      peptide_sequence_column = "Stripped.Sequence",
      core_utilisation = NA
    ),
    "Unmapped input run\\(s\\): UNMAPPED"
  )

  sample_input <- data.frame(
    Run = c("S1", "S2"),
    Protein.Group = c("PG1", "PG2"),
    Stripped.Sequence = c("A", "B"),
    stringsAsFactors = FALSE
  )
  expect_error(
    filterMinNumPeptidesPerSampleHelper(
      sample_input,
      peptides_per_sample_cutoff = 1,
      sample_id_column = "Run",
      core_utilisation = NA,
      protein_id_column = "Protein.Group",
      peptide_sequence_column = "Stripped.Sequence",
      design_matrix = data.frame(Run = c("S1", "S1")),
      design_sample_id_column = "Run"
    ),
    "duplicate run ID\\(s\\): S1"
  )
})

test_that("packaged sample threshold is 500", {
  config_lines <- readLines(
    test_path("..", "..", "inst", "config", "config.ini"),
    warn = FALSE
  )
  section_start <- match("[filterMinNumPeptidesPerSample]", config_lines)
  expect_false(is.na(section_start))
  expect_identical(
    config_lines[[section_start + 1L]],
    "peptides_per_sample_cutoff=500"
  )
})

test_that("S4 sample threshold resolves to 500 and preserves explicit config override", {
  peptide_data <- data.frame(
    Run = "S1",
    Protein.Group = paste0("PG", seq_len(500)),
    Protein.Ids = paste0("P", seq_len(500)),
    Stripped.Sequence = paste0("PEP", seq_len(500)),
    Modified.Sequence = paste0("PEP", seq_len(500)),
    Precursor.Id = paste0("PREC", seq_len(500)),
    Precursor.Charge = 2L,
    Precursor.Quantity = 100,
    Precursor.Normalised = 100,
    Q.Value = 0.001,
    Global.Q.Value = 0.001,
    Global.PG.Q.Value = 0.001,
    Proteotypic = 1,
    stringsAsFactors = FALSE
  )
  design <- data.frame(
    Run = "S1",
    group = "A",
    replicates = 1L,
    stringsAsFactors = FALSE
  )
  default_object <- PeptideQuantitativeDataDiann(peptide_data, design, args = list())
  default_filtered <- filterMinNumPeptidesPerSample(default_object)
  expect_identical(
    default_filtered@args$filterMinNumPeptidesPerSample$peptides_per_sample_cutoff,
    500
  )
  expect_identical(
    default_filtered@args$filterMinNumPeptidesPerSample$filter_summary$peptides_per_sample_cutoff,
    500L
  )

  configured_object <- PeptideQuantitativeDataDiann(peptide_data, design, args = list())
  configured_object@args$filterMinNumPeptidesPerSample <- list(
    peptides_per_sample_cutoff = 499,
    inclusion_list = character(),
    core_utilisation = NA
  )
  configured_filtered <- filterMinNumPeptidesPerSample(configured_object)
  expect_identical(
    configured_filtered@args$filterMinNumPeptidesPerSample$peptides_per_sample_cutoff,
    499
  )
  expect_identical(
    configured_filtered@args$filterMinNumPeptidesPerSample$filter_summary$peptides_per_sample_cutoff,
    499L
  )
})

test_that("replicate S4 method accepts historical grouping_variable config alias", {
  fixture <- makeDistinctIdentityFilterFixture()
  peptide_data <- fixture$input |>
    dplyr::mutate(
      Modified.Sequence = Stripped.Sequence,
      Precursor.Id = paste0("PREC", dplyr::row_number()),
      Precursor.Charge = 2L,
      Precursor.Quantity = Peptide.Normalised,
      Q.Value = 0.001,
      Global.Q.Value = 0.001,
      Global.PG.Q.Value = 0.001,
      Proteotypic = 1
    )
  design <- fixture$design |>
    dplyr::rename(group = replicate_group) |>
    dplyr::mutate(replicates = seq_len(dplyr::n()))
  object <- PeptideQuantitativeDataDiann(peptide_data, design, args = list())
  object@args$removePeptidesWithOnlyOneReplicate <- list(
    grouping_variable = "group",
    core_utilisation = NA
  )

  filtered <- removePeptidesWithOnlyOneReplicate(object)
  expect_identical(
    filtered@args$removePeptidesWithOnlyOneReplicate$replicate_group_column,
    "group"
  )
  expect_identical(
    filtered@args$removePeptidesWithOnlyOneReplicate$retention_policy,
    "supported_in_any_group"
  )
  expect_setequal(unique(filtered@peptide_data$Protein.Group), c("PG1", "PG3"))
})

test_that("protein peptide evidence counts identities rather than repeated run rows", {
  runs <- paste0("S", seq_len(12))
  repeated_single_peptide <- data.frame(
    Run = runs,
    Protein.Ids = "P_REPEAT",
    Stripped.Sequence = "ONLY_ONE",
    peptidoform_count = 1L,
    peptidoform_ids = I(rep(list("ONLY_ONE"), length(runs))),
    stringsAsFactors = FALSE
  )
  two_unique_peptides <- data.frame(
    Run = rep(runs, each = 2L),
    Protein.Ids = "P_TWO",
    Stripped.Sequence = rep(c("PEPTIDE_A", "PEPTIDE_B"), times = length(runs)),
    peptidoform_count = 1L,
    peptidoform_ids = I(rep(list("PEPTIDE_A", "PEPTIDE_B"), times = length(runs))),
    stringsAsFactors = FALSE
  )

  expect_warning(
    filtered <- filterMinNumPeptidesPerProteinHelper(
      input_table = rbind(repeated_single_peptide, two_unique_peptides),
      num_peptides_per_protein_thresh = 2,
      num_peptidoforms_per_protein_thresh = 2,
      protein_id_column = Protein.Ids,
      peptide_sequence_column = Stripped.Sequence,
      core_utilisation = NA
    ),
    "legacy peptidoform_ids",
    fixed = TRUE
  )

  expect_identical(unique(filtered$Protein.Ids), "P_TWO")
  expect_true(all(filtered$peptides_for_protein_count == 2L))
  expect_true(all(filtered$peptidoforms_for_protein_count == 2L))

  refiltered <- suppressWarnings(filterMinNumPeptidesPerProteinHelper(
    input_table = filtered,
    num_peptides_per_protein_thresh = 2,
    num_peptidoforms_per_protein_thresh = 2,
    protein_id_column = Protein.Ids,
    peptide_sequence_column = Stripped.Sequence,
    core_utilisation = NA
  ))
  expect_equal(refiltered, filtered)
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
    Global.PG.Q.Value = c(0.002, 0.030, 0.004),
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

test_that("DIA-NN confidence fields are distinct fixed gates", {
  input_table <- data.frame(
    Run = paste0("S", 1:5),
    Precursor.Id = paste0("prec", 1:5),
    Protein.Ids = paste0("P", 1:5),
    Stripped.Sequence = paste0("PEP", 1:5),
    Modified.Sequence = paste0("PEP", 1:5),
    Precursor.Charge = 2L,
    Precursor.Quantity = 100,
    Precursor.Normalised = 10,
    Q.Value = c(0.001, 0.02, 0.001, 0.001, 0.001),
    Global.Q.Value = c(0.001, 0.001, 0.02, 0.001, 0.001),
    Global.PG.Q.Value = c(0.001, 0.001, 0.001, 0.02, 0.001),
    Proteotypic = c(1, 1, 1, 1, 0),
    PG.Q.Value = 0.001,
    stringsAsFactors = FALSE
  )

  filtered <- srlQvalueProteotypicPeptideCleanHelper(input_table)

  expect_identical(filtered$Run, "S1")
  expect_error(
    srlQvalueProteotypicPeptideCleanHelper(
      input_table[, setdiff(names(input_table), "Global.Q.Value")]
    ),
    "Global.Q.Value",
    fixed = TRUE
  )
  expect_error(
    srlQvalueProteotypicPeptideCleanHelper(
      input_table,
      global_qvalue_threshold = 0.05
    ),
    "fixed at 0.01",
    fixed = TRUE
  )
  expect_error(
    srlQvalueProteotypicPeptideCleanHelper(
      input_table,
      choose_only_proteotypic_peptide = 0
    ),
    "Proteotypic == 1",
    fixed = TRUE
  )
})

test_that("upstream-prefiltered confidence mode records unavailable metrics", {
  input_table <- data.frame(
    Run = "S1",
    Precursor.Id = "prec1",
    Protein.Ids = "P1",
    Stripped.Sequence = "PEP1",
    Modified.Sequence = "PEP1",
    Precursor.Charge = 2L,
    Precursor.Quantity = 100,
    Precursor.Normalised = 10,
    stringsAsFactors = FALSE
  )

  expect_warning(
    filtered <- srlQvalueProteotypicPeptideCleanHelper(
      input_table,
      confidence_provenance_mode = "upstream_prefiltered"
    ),
    "were not verified",
    fixed = TRUE
  )
  expect_identical(nrow(filtered), 1L)
  expect_false(any(c(
    "Q.Value",
    "Global.Q.Value",
    "Global.PG.Q.Value",
    "Proteotypic"
  ) %in% names(filtered)))
})
