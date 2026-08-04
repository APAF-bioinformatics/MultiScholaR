library(testthat)

test_that("MCI-006.1 peptide intensity filter covers zero, boundary, above-max, missing, all-zero, and non-finite cases", {
  tbl <- module_ci_prot_peptide_table()
  design <- module_ci_prot_peptide_design()

  zero_threshold <- module_ci_prot_intensity_filter(tbl, design, threshold = 0)
  expect_setequal(
    module_ci_prot_peptide_ids(zero_threshold),
    c("P_PASS%PEP_PASS", "P_BOUNDARY%PEP_BOUNDARY", "P_MISSING%PEP_MISSING", "P_ZERO%PEP_ZERO")
  )
  module_ci_prot_assert_peptide_alignment(zero_threshold, design)

  boundary_threshold <- module_ci_prot_intensity_filter(tbl, design, threshold = 10)
  expect_setequal(
    module_ci_prot_peptide_ids(boundary_threshold),
    c("P_PASS%PEP_PASS", "P_BOUNDARY%PEP_BOUNDARY", "P_MISSING%PEP_MISSING")
  )
  expect_false("P_ZERO%PEP_ZERO" %in% module_ci_prot_peptide_ids(boundary_threshold))
  expect_false("P_INF%PEP_INF" %in% module_ci_prot_peptide_ids(boundary_threshold))

  strict_missing <- module_ci_prot_intensity_filter(
    tbl,
    design,
    threshold = 10,
    group_cutoff = 49,
    max_groups_cutoff = 49
  )
  expect_false("P_MISSING%PEP_MISSING" %in% module_ci_prot_peptide_ids(strict_missing))

  above_max <- module_ci_prot_intensity_filter(tbl, design, threshold = 1000)
  expect_equal(nrow(above_max), 0L)
})

test_that("MCI-006.2 peptide q-value filter covers valid, inclusive-boundary, missing, malformed, all-pass, and all-fail cases", {
  tbl <- module_ci_prot_peptide_table()

  filtered <- module_ci_prot_qvalue_filter(tbl, q_threshold = 0.01, global_threshold = 0.01)
  expect_setequal(
    module_ci_prot_peptide_ids(filtered),
    c("P_PASS%PEP_PASS", "P_BOUNDARY%PEP_BOUNDARY")
  )
  expect_true(all(filtered$Q.Value <= 0.01))
  expect_true(all(filtered$Global.Q.Value <= 0.01))
  expect_false(any(is.na(filtered$Q.Value) | is.na(filtered$Global.Q.Value)))

  all_pass_tbl <- module_ci_prot_peptide_table(
    q_values = rep(0.001, 5),
    global_q_values = rep(0.001, 5),
    proteotypic = rep(1, 5)
  )
  expect_setequal(
    module_ci_prot_peptide_ids(module_ci_prot_qvalue_filter(all_pass_tbl)),
    module_ci_prot_peptide_ids(all_pass_tbl)
  )
  expect_error(
    module_ci_prot_qvalue_filter(all_pass_tbl, q_threshold = 0, global_threshold = 0),
    "fixed at 0.01",
    fixed = TRUE
  )

  expect_error(
    module_ci_prot_qvalue_filter(tbl[, setdiff(names(tbl), "Q.Value")]),
    "Required filter columns not found",
    fixed = TRUE
  )
  expect_error(
    module_ci_prot_qvalue_filter(tbl[, setdiff(names(tbl), "Precursor.Id")]),
    "Required output columns not found",
    fixed = TRUE
  )
})

test_that("MCI-006.3 peptide sample and replicate filters cover balanced, unbalanced, duplicate, and missing-design cases", {
  tbl <- module_ci_prot_peptide_table()
  design <- module_ci_prot_peptide_design()

  sample_tbl <- rbind(
    tbl[tbl$Stripped.Sequence %in% c("PEP_PASS", "PEP_BOUNDARY"), ],
    tbl[tbl$Run == "B2" & tbl$Stripped.Sequence == "PEP_MISSING", ]
  )
  sample_filtered <- suppressMessages(filterMinNumPeptidesPerSampleHelper(
    sample_tbl,
    peptides_per_sample_cutoff = 2,
    sample_id_column = Run,
    core_utilisation = NA,
    inclusion_list = character()
  ))
  expect_setequal(unique(sample_filtered$Run), c("A1", "A2", "B1", "B2"))

  sample_filtered_strict <- suppressMessages(filterMinNumPeptidesPerSampleHelper(
    sample_tbl,
    peptides_per_sample_cutoff = 3,
    sample_id_column = Run,
    core_utilisation = NA,
    inclusion_list = "B2"
  ))
  expect_setequal(unique(sample_filtered_strict$Run), c("B2"))

  replicate_tbl <- rbind(
    tbl[tbl$Stripped.Sequence == "PEP_PASS", ],
    tbl[tbl$Stripped.Sequence == "PEP_MISSING" & tbl$Run %in% c("A1", "B1"), ],
    tbl[tbl$Stripped.Sequence == "PEP_PASS" & tbl$Run == "A1", ]
  )
  replicate_filtered <- suppressMessages(removePeptidesWithOnlyOneReplicateHelper(
    input_table = replicate_tbl,
    samples_id_tbl = design,
    input_table_sample_id_column = Run,
    sample_id_tbl_sample_id_column = Run,
    replicate_group_column = replicate_group,
    protein_id_column = Protein.Ids,
    peptide_sequence_column = Stripped.Sequence,
    core_utilisation = NA
  ))
  expect_identical(module_ci_prot_peptide_ids(replicate_filtered), "P_PASS%PEP_PASS")
  expect_true(nrow(replicate_filtered) >= 4L)

  expect_error(
    module_ci_prot_validate_peptide_design(tbl, design[design$Run != "B2", , drop = FALSE]),
    "peptide runs and design runs must match exactly once",
    fixed = TRUE
  )
})

test_that("MCI-006.4 peptide imputation matrix covers deterministic imputation, no-imputation, all-missing, and minimum-observation guards", {
  fixture <- module_ci_prot_peptide_imputation_fixture()

  imputed <- module_ci_prot_impute(fixture, proportion_missing_values = 0.75)
  p1_a2 <- imputed$Peptide.Imputed[imputed$Protein.Ids == "P1" & imputed$Stripped.Sequence == "PEP1" & imputed$Run == "A2"]
  p1_b2 <- imputed$Peptide.Imputed[imputed$Protein.Ids == "P1" & imputed$Stripped.Sequence == "PEP1" & imputed$Run == "B2"]
  all_missing <- imputed$Peptide.Imputed[imputed$Stripped.Sequence == "ALL_MISSING"]

  expect_equal(p1_a2, 10)
  expect_equal(p1_b2, 30)
  expect_true(all(is.na(all_missing)))
  expect_true(any(imputed$is_imputed))

  boundary_impute <- module_ci_prot_impute(fixture, proportion_missing_values = 0.50)
  expect_true(any(boundary_impute$is_imputed, na.rm = TRUE))
  expect_true(all(is.na(boundary_impute$Peptide.Imputed[boundary_impute$Stripped.Sequence == "ALL_MISSING"])))

  empty_fixture <- fixture
  empty_fixture$data <- fixture$data[0, , drop = FALSE]
  empty_impute <- module_ci_prot_impute(
    empty_fixture,
    proportion_missing_values = 0.75,
    core_utilisation = structure(list(), class = "multidplyr_cluster")
  )
  expect_identical(nrow(empty_impute), 0L)
  expect_true(all(c("Peptide.Imputed", "is_imputed") %in% names(empty_impute)))
  expect_identical(empty_impute$is_imputed, logical(0))
})

test_that("MCI-006.5 peptide QC plots and summaries survive normal and degenerate small fixtures", {
  finite_tbl <- module_ci_prot_peptide_table(values = matrix(
    c(
      10, 11, 12, 13,
      20, 21, 22, 23,
      30, 31, 32, 33,
      40, 41, 42, 43,
      50, 51, 52, 53
    ),
    nrow = 5,
    byrow = TRUE
  ))
  finite_obj <- module_ci_prot_peptide_object(finite_tbl)

  density_plot <- plotDensity(finite_obj, grouping_variable = "group", title = "density")
  pca_plot <- plotPca(finite_obj, grouping_variable = "group", label_column = NULL, title = "pca")
  rle_plot <- plotRle(finite_obj, grouping_variable = "group")
  summary_tbl <- summarisePeptideObject(finite_obj)

  expect_s3_class(density_plot, "ggplot")
  expect_s3_class(pca_plot, "ggplot")
  expect_s3_class(rle_plot, "ggplot")
  expect_identical(summary_tbl$num_peptides, 5L)
  expect_identical(summary_tbl$num_samples, 4L)

  degenerate_tbl <- module_ci_prot_peptide_table(values = matrix(1, nrow = 5, ncol = 4))
  degenerate_obj <- module_ci_prot_peptide_object(degenerate_tbl)
  expect_s3_class(plotPca(degenerate_obj, grouping_variable = "group", label_column = NULL, title = "degenerate pca"), "ggplot")
  expect_s3_class(plotRle(degenerate_obj, grouping_variable = "group"), "ggplot")
})

test_that("MCI-006.6 peptide QC state transitions are recorded once and protein QC is gated on peptide completion", {
  workflow <- module_ci_prot_workflow()
  now <- as.POSIXct("2026-05-05 00:00:00", tz = "UTC")
  no_log <- function(...) invisible(NULL)
  identity_filter <- function(theObject, ...) theObject

  runPeptideQvalueApplyStep(
    workflow,
    qvalueThreshold = 0.01,
    globalQvalueThreshold = 0.01,
    proteotypicOnly = TRUE,
    updateConfigParameterFn = module_ci_prot_update_config,
    qvalueFilterFn = identity_filter,
    logInfoFn = no_log,
    logWarnFn = no_log,
    nowFn = function() now
  )
  runPeptideRollupApplyStep(
    workflow,
    rollupFn = identity_filter,
    logInfoFn = no_log,
    nowFn = function() now
  )
  runPeptideIntensityApplyStep(
    workflow,
    useStrictMode = FALSE,
    minRepsPerGroup = 2,
    minGroups = 2,
    intensityCutoffPercentile = 1,
    updateConfigParameterFn = module_ci_prot_update_config,
    updateMissingValueParametersFn = module_ci_prot_update_missing,
    peptideIntensityFilteringFn = identity_filter,
    logInfoFn = no_log,
    nowFn = function() now
  )
  runProteinPeptideApplyStep(
    workflow,
    minPeptidesPerProtein = 2,
    minPeptidoformsPerProtein = 2,
    updateConfigParameterFn = module_ci_prot_update_config,
    filterMinNumPeptidesPerProteinFn = identity_filter,
    logInfoFn = no_log,
    nowFn = function() now
  )
  runPeptideSampleApplyStep(
    workflow,
    minPeptidesPerSample = 2,
    updateConfigParameterFn = module_ci_prot_update_config,
    filterMinNumPeptidesPerSampleFn = identity_filter,
    logInfoFn = no_log,
    nowFn = function() now
  )
  runPeptideReplicateApplyStep(
    workflow,
    replicateGroupColumn = "replicate_group",
    removePeptidesWithOnlyOneReplicateFn = identity_filter,
    logInfoFn = no_log,
    nowFn = function() now
  )

  testthat::local_mocked_bindings(
    updateConfigParameter = module_ci_prot_update_config,
    peptideMissingValueImputation = identity_filter,
    .capture_checkpoint = function(...) invisible(TRUE),
    .env = environment(runPeptideImputationStep)
  )
  runPeptideImputationStep(workflow, proportionMissingValues = 0.75)

  expected_states <- c(
    "qvalue_filtered",
    "precursor_rollup",
    "intensity_filtered",
    "protein_peptide_filtered",
    "sample_filtered",
    "replicate_filtered",
    "imputed"
  )
  expect_identical(workflow$state_manager$getHistory(), expected_states)
  expect_false(anyDuplicated(workflow$state_manager$getHistory()) > 0L)

  blocked_workflow <- module_ci_prot_workflow()
  expect_error(module_ci_prot_validate_peptide_prereqs(blocked_workflow), "peptide state 'imputed' is required", fixed = TRUE)
  expect_silent(module_ci_prot_validate_peptide_prereqs(workflow))
})
