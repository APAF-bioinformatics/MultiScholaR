# fidelity-coverage-compare: shared
library(testthat)

fakeBangBangSym <- function(x) {
  structure(x, class = c("fake_bangbang_sym", "character"))
}

`!.fake_bangbang_sym` <- function(e1) {
  structure(unclass(e1), class = c("fake_bangbang_unquote", "character"))
}

`!.fake_bangbang_unquote` <- function(e1) {
  unclass(e1)
}

newPeptideQcMethodsObject <- function(raw_values, norm_values) {
  peptide_data <- data.frame(
    Protein.Ids = c("P1", "P1", "P2", "P2"),
    Stripped.Sequence = c("SEQ1", "SEQ2", "SEQ1", "SEQ2"),
    Proteotypic.Sequence = c("SEQ1", "SEQ2", "SEQ1", "SEQ2"),
    Run = c("S1", "S2", "S1", "S2"),
    Precursor.Quantity = raw_values,
    Precursor.Normalised = norm_values,
    Q.Value = c(0.001, 0.002, 0.003, 0.004),
    Global.Q.Value = c(0.001, 0.002, 0.003, 0.004),
    Global.PG.Q.Value = c(0.001, 0.002, 0.003, 0.004),
    stringsAsFactors = FALSE
  )

  methods::new(
    "PeptideQuantitativeData",
    peptide_data = peptide_data,
    peptide_matrix = matrix(
      seq_len(4),
      nrow = 2,
      dimnames = list(c("SEQ1", "SEQ2"), c("S1", "S2"))
    ),
    protein_id_column = "Protein.Ids",
    peptide_sequence_column = "Stripped.Sequence",
    q_value_column = "Q.Value",
    global_q_value_column = "Global.Q.Value",
    proteotypic_peptide_sequence_column = "Proteotypic.Sequence",
    raw_quantity_column = "Precursor.Quantity",
    norm_quantity_column = "Precursor.Normalised",
    is_logged_data = FALSE,
    design_matrix = data.frame(
      Run = c("S1", "S2"),
      group = c("A", "B"),
      replicates = c("R1", "R1"),
      stringsAsFactors = FALSE
    ),
    sample_id = "Run",
    group_id = "group",
    technical_replicate_id = "replicates",
    args = list()
  )
}

appendMethodCall <- function(log_env, key, value) {
  existing <- if (exists(key, envir = log_env, inherits = FALSE)) {
    get(key, envir = log_env, inherits = FALSE)
  } else {
    list()
  }
  assign(key, c(existing, list(value)), envir = log_env)
}

localPeptideQcMethodBindings <- function(log_env, .local_envir = parent.frame()) {
  local_mocked_bindings(
    checkParamsObjectFunctionSimplify = function(theObject, function_param, default = NULL) {
      switch(
        function_param,
        grouping_variable = "group",
        groupwise_percentage_cutoff = 25,
        max_groups_percentage_cutoff = 50,
        peptides_intensity_cutoff_percentile = 50,
        core_utilisation = NA,
        num_peptides_per_protein_thresh = 2,
        num_peptidoforms_per_protein_thresh = 3,
        peptides_per_sample_cutoff = 100,
        qvalue_threshold = 0.01,
        global_qvalue_threshold = 0.01,
        global_pg_qvalue_threshold = 0.01,
        choose_only_proteotypic_peptide = 1,
        confidence_provenance_mode = "diann_main_report",
        input_matrix_column_ids = c("Run", "CustomColumn"),
        default
      )
    },
    checkParamsObjectFunctionSimplifyAcceptNull = function(theObject, function_param, default = NULL) {
      switch(
        function_param,
        replicate_group_column = "replicates",
        inclusion_list = c("S1"),
        default
      )
    },
    sym = fakeBangBangSym,
    updateParamInObject = function(theObject, function_param) {
      appendMethodCall(log_env, "updated_params", function_param)
      theObject
    },
    cleanDesignMatrixPeptide = function(theObject) {
      appendMethodCall(log_env, "clean_calls", TRUE)
      theObject
    },
    peptideIntensityFilteringHelper = function(...) {
      args <- list(...)
      appendMethodCall(log_env, "intensity_calls", args)
      args$input_table[1:2, , drop = FALSE]
    },
    removePeptidesWithMissingValuesPercentHelper = function(...) {
      args <- list(...)
      appendMethodCall(log_env, "missing_calls", args)
      args$input_table[1, , drop = FALSE]
    },
    removePeptidesWithOnlyOneReplicateHelper = function(...) {
      args <- list(...)
      appendMethodCall(log_env, "replicate_calls", args)
      args$input_table[1, , drop = FALSE]
    },
    filterMinNumPeptidesPerProteinHelper = function(...) {
      args <- list(...)
      appendMethodCall(log_env, "protein_calls", args)
      args$input_table[1, , drop = FALSE]
    },
    filterMinNumPeptidesPerSampleHelper = function(input_table,
                                                   peptides_per_sample_cutoff,
                                                   sample_id_column,
                                                   core_utilisation,
                                                   inclusion_list = NULL,
                                                   ...) {
      appendMethodCall(
        log_env,
        "sample_calls",
        list(
          input_table = input_table,
          peptides_per_sample_cutoff = peptides_per_sample_cutoff,
          sample_id_column = sample_id_column,
          core_utilisation = core_utilisation,
          inclusion_list = inclusion_list,
          extra_args = list(...)
        )
      )
      input_table[1, , drop = FALSE]
    },
    srlQvalueProteotypicPeptideCleanHelper = function(...) {
      args <- list(...)
      appendMethodCall(log_env, "srl_calls", args)
      args$input_table[1, , drop = FALSE]
    },
    .env = .local_envir
  )
}

test_that("peptide QC S4 methods preserve threshold resolution and invalid-value fallback branches", {
  valid_object <- newPeptideQcMethodsObject(
    raw_values = c(10, 20, 30, 40),
    norm_values = c(1, 2, 3, 4)
  )
  invalid_object <- newPeptideQcMethodsObject(
    raw_values = c(NA_real_, Inf, NaN, NA_real_),
    norm_values = c(NA_real_, Inf, NaN, NA_real_)
  )

  log_env <- new.env(parent = emptyenv())
  localPeptideQcMethodBindings(log_env)

  valid_intensity <- peptideIntensityFiltering(valid_object)
  invalid_intensity <- peptideIntensityFiltering(invalid_object)
  valid_missing <- removePeptidesWithMissingValuesPercent(valid_object)
  invalid_missing <- removePeptidesWithMissingValuesPercent(invalid_object)

  expect_s4_class(valid_intensity, "PeptideQuantitativeData")
  expect_s4_class(invalid_intensity, "PeptideQuantitativeData")
  expect_s4_class(valid_missing, "PeptideQuantitativeData")
  expect_s4_class(invalid_missing, "PeptideQuantitativeData")

  expect_equal(log_env$intensity_calls[[1]]$min_peptide_intensity_threshold, 2)
  expect_identical(
    log_env$intensity_calls[[1]]$peptide_quantity_column,
    "Precursor.Normalised"
  )
  expect_equal(log_env$intensity_calls[[2]]$min_peptide_intensity_threshold, 0)
  expect_equal(log_env$missing_calls[[1]]$abundance_threshold, 3)
  expect_equal(log_env$missing_calls[[2]]$abundance_threshold, 0)
  expect_equal(nrow(valid_intensity@peptide_data), 2)
  expect_equal(nrow(valid_missing@peptide_data), 1)
  expect_true(length(log_env$clean_calls) >= 4)
})

test_that("peptide intensity helpers resolve S4 slot column expressions", {
  peptide_object <- newPeptideQcMethodsObject(
    raw_values = c(10, 20, 30, 40),
    norm_values = c(1, 2, 3, 4)
  )

  filtered <- peptideIntensityFilteringHelper(
    input_table = peptide_object@peptide_data,
    design_matrix = peptide_object@design_matrix,
    min_peptide_intensity_threshold = 10,
    sample_id_column = peptide_object@sample_id,
    grouping_variable = peptide_object@group_id,
    groupwise_percentage_cutoff = 100,
    max_groups_percentage_cutoff = 100,
    protein_id_column = peptide_object@protein_id_column,
    peptide_sequence_column = peptide_object@peptide_sequence_column,
    peptide_quantity_column = peptide_object@raw_quantity_column,
    core_utilisation = NA,
    min_reps_per_group = 1,
    min_groups = 1
  )

  missing_filtered <- removePeptidesWithMissingValuesPercentHelper(
    input_table = peptide_object@peptide_data,
    design_matrix = peptide_object@design_matrix,
    sample_id = peptide_object@sample_id,
    protein_id_column = peptide_object@protein_id_column,
    peptide_sequence_column = peptide_object@peptide_sequence_column,
    grouping_variable = peptide_object@group_id,
    groupwise_percentage_cutoff = 100,
    max_groups_percentage_cutoff = 100,
    abundance_threshold = 0,
    abundance_column = peptide_object@norm_quantity_column
  )

  expect_equal(nrow(filtered), nrow(peptide_object@peptide_data))
  expect_equal(nrow(missing_filtered), nrow(peptide_object@peptide_data))
  expect_identical(names(filtered), names(peptide_object@peptide_data))
  expect_identical(names(missing_filtered), names(peptide_object@peptide_data))
})

test_that("peptide intensity method defaults to normalised quantity and permits raw override", {
  peptide_object <- newPeptideQcMethodsObject(
    raw_values = c(10, 20, 30, 40),
    norm_values = c(1, 2, 3, 4)
  )

  normalised <- peptideIntensityFiltering(
    peptide_object,
    min_reps_per_group = 1,
    min_groups = 1,
    peptides_intensity_cutoff_percentile = 50
  )
  raw <- peptideIntensityFiltering(
    peptide_object,
    min_reps_per_group = 1,
    min_groups = 1,
    peptides_intensity_cutoff_percentile = 50,
    peptide_quantity_column = "Precursor.Quantity"
  )

  expect_identical(
    normalised@args$peptideIntensityFiltering$peptide_quantity_column,
    "Precursor.Normalised"
  )
  expect_equal(
    normalised@args$peptideIntensityFiltering$min_peptide_intensity_threshold,
    3
  )
  expect_identical(
    raw@args$peptideIntensityFiltering$peptide_quantity_column,
    "Precursor.Quantity"
  )
  expect_equal(
    raw@args$peptideIntensityFiltering$min_peptide_intensity_threshold,
    25
  )
})

test_that("config-only direct peptide intensity counts execute without percentage fallback", {
  peptide_object <- newPeptideQcMethodsObject(
    raw_values = c(10, 20, 30, 40),
    norm_values = c(1, 2, 3, 4)
  )
  peptide_object@args$peptideIntensityFiltering <- list(
    grouping_variable = "group",
    min_reps_per_group = 1,
    min_groups = 1,
    strict_mode = FALSE,
    peptides_intensity_cutoff_percentile = 50
  )

  expect_no_warning(filtered <- peptideIntensityFiltering(peptide_object))
  expect_identical(
    filtered@args$peptideIntensityFiltering$filter_summary$rule_mode,
    "direct_counts"
  )
  expect_identical(
    filtered@args$peptideIntensityFiltering$filter_summary$min_reps_per_group,
    1L
  )
  expect_s3_class(
    filtered@args$peptideIntensityFiltering$support_table,
    "data.frame"
  )
})

test_that("saved percentage-only peptide intensity configs use the explicit adapter", {
  peptide_object <- newPeptideQcMethodsObject(
    raw_values = c(10, 20, 30, 40),
    norm_values = c(1, 2, 3, 4)
  )
  peptide_object@args$peptideIntensityFiltering <- list(
    grouping_variable = "group",
    groupwise_percentage_cutoff = 50,
    max_groups_percentage_cutoff = 50,
    peptides_intensity_cutoff_percentile = 50
  )

  expect_warning(
    filtered <- peptideIntensityFiltering(peptide_object),
    "percentage-only missingness arguments are deprecated"
  )
  expect_identical(
    filtered@args$peptideIntensityFiltering$filter_summary$rule_mode,
    "legacy_percentage_adapter"
  )
  expect_identical(
    filtered@args$peptideIntensityFiltering$filter_summary$min_groups,
    1L
  )
})

test_that("peptide sample filtering honors object args from module updates", {
  peptide_object <- newPeptideQcMethodsObject(
    raw_values = c(10, 20, 30, 40),
    norm_values = c(1, 2, 3, 4)
  )
  peptide_object@args$filterMinNumPeptidesPerSample <- list(
    peptides_per_sample_cutoff = 1,
    inclusion_list = character(0),
    core_utilisation = NA
  )

  filtered <- filterMinNumPeptidesPerSample(peptide_object)

  expect_s4_class(filtered, "PeptideQuantitativeData")
  expect_equal(nrow(filtered@peptide_data), nrow(peptide_object@peptide_data))
  expect_identical(sort(unique(filtered@peptide_data$Run)), c("S1", "S2"))
})

test_that("peptide QC S4 filtering and DIA cleanup methods preserve helper delegation and slot updates", {
  peptide_object <- newPeptideQcMethodsObject(
    raw_values = c(11, 21, 31, 41),
    norm_values = c(2, 4, 6, 8)
  )

  log_env <- new.env(parent = emptyenv())
  localPeptideQcMethodBindings(log_env)

  replicate_filtered <- removePeptidesWithOnlyOneReplicate(
    peptide_object,
    replicate_group_column = "replicates"
  )
  protein_filtered <- filterMinNumPeptidesPerProtein(peptide_object)
  sample_filtered <- filterMinNumPeptidesPerSample(peptide_object)
  cleaned <- srlQvalueProteotypicPeptideClean(peptide_object)

  expect_s4_class(replicate_filtered, "PeptideQuantitativeData")
  expect_s4_class(protein_filtered, "PeptideQuantitativeData")
  expect_s4_class(sample_filtered, "PeptideQuantitativeData")
  expect_s4_class(cleaned, "PeptideQuantitativeData")

  expect_equal(nrow(replicate_filtered@peptide_data), 1)
  expect_equal(nrow(protein_filtered@peptide_data), 1)
  expect_equal(nrow(sample_filtered@peptide_data), 1)
  expect_equal(nrow(cleaned@peptide_data), 1)

  expect_identical(as.character(log_env$replicate_calls[[1]]$replicate_group_column), "replicates")
  expect_identical(log_env$protein_calls[[1]]$num_peptides_per_protein_thresh, 2)
  expect_identical(log_env$protein_calls[[1]]$num_peptidoforms_per_protein_thresh, 3)
  expect_identical(log_env$sample_calls[[1]]$peptides_per_sample_cutoff, 100)
  expect_identical(log_env$sample_calls[[1]]$inclusion_list, "S1")
  expect_equal(
    log_env$srl_calls[[1]]$input_matrix_column_ids,
    c(
      "Run",
      "CustomColumn",
      "Protein.Ids",
      "Stripped.Sequence",
      "Modified.Sequence"
    )
  )
  expect_identical(log_env$srl_calls[[1]]$global_qvalue_threshold, 0.01)
  expect_identical(log_env$srl_calls[[1]]$global_pg_qvalue_threshold, 0.01)
  expect_identical(log_env$srl_calls[[1]]$qvalue_threshold, 0.01)
  expect_identical(log_env$srl_calls[[1]]$choose_only_proteotypic_peptide, 1)
  expect_identical(
    cleaned@args$srlQvalueProteotypicPeptideClean$qvalue_threshold,
    0.01
  )
  expect_identical(
    cleaned@args$srlQvalueProteotypicPeptideClean$global_qvalue_threshold,
    0.01
  )
  expect_identical(
    cleaned@args$srlQvalueProteotypicPeptideClean$global_pg_qvalue_threshold,
    0.01
  )
  expect_identical(
    cleaned@args$srlQvalueProteotypicPeptideClean$choose_only_proteotypic_peptide,
    1
  )
  expect_identical(
    cleaned@args$srlQvalueProteotypicPeptideClean$input_matrix_column_ids,
    c("Run", "CustomColumn")
  )
  expect_identical(
    protein_filtered@args$filterMinNumPeptidesPerProtein$num_peptides_per_protein_thresh,
    2
  )
  expect_identical(
    protein_filtered@args$filterMinNumPeptidesPerProtein$num_peptidoforms_per_protein_thresh,
    3
  )
  expect_identical(
    replicate_filtered@args$removePeptidesWithOnlyOneReplicate$replicate_group_column,
    "replicates"
  )
  expect_identical(
    replicate_filtered@args$removePeptidesWithOnlyOneReplicate$minimum_distinct_runs,
    2L
  )
  expect_identical(
    replicate_filtered@args$removePeptidesWithOnlyOneReplicate$retention_policy,
    "supported_in_any_group"
  )
  expect_identical(
    sample_filtered@args$filterMinNumPeptidesPerSample$peptides_per_sample_cutoff,
    100
  )
  expect_identical(
    sample_filtered@args$filterMinNumPeptidesPerSample$inclusion_list,
    "S1"
  )
})

makePeptideImputationContractFixture <- function() {
  design <- data.frame(
    Run = c("A1", "A2", "B1", "B2", "B3", "S1", "HEK_subject_1", "H2", "QC1", "QC2"),
    technical_group = c("A", "A", "B", "B", "B", "SINGLE", "H", "H", "QC", "QC"),
    exclude_imputation = c(FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE),
    stringsAsFactors = FALSE
  )
  data <- data.frame(
    Run = c(
      "A1", "A2",
      "A1", "A2",
      "A1", "A2", "B1", "B2", "B3",
      "B1", "B2", "B3",
      "B1", "B2", "B3",
      "S1",
      "HEK_subject_1", "H2",
      "QC1", "QC2"
    ),
    Protein.Group = c(
      "P_EXACT", "P_EXACT",
      "P_ALL", "P_ALL",
      rep("P_MEAN", 5),
      rep("P_EQUAL", 3),
      rep("P_OVER", 3),
      "P_SINGLE",
      "P_HEK", "P_HEK",
      "P_QC", "P_QC"
    ),
    Stripped.Sequence = c(
      "EXACT", "EXACT",
      "ALL", "ALL",
      rep("MEAN", 5),
      rep("EQUAL", 3),
      rep("OVER", 3),
      "SINGLE",
      "HEK_OK", "HEK_OK",
      "QC", "QC"
    ),
    Precursor.Quantity = c(
      10, NA,
      NA, NA,
      100, 100, 10, 20, NA,
      5, 5, NA,
      3, NA, NA,
      7,
      8, NA,
      9, NA
    ),
    stringsAsFactors = FALSE
  )
  list(data = data, design = design)
}

runPeptideImputationContractFixture <- function(fixture = makePeptideImputationContractFixture(),
                                                proportion_missing_values = 0.5,
                                                exclusion_column = "exclude_imputation",
                                                hek_string = NULL,
                                                return_imputation_result = TRUE) {
  peptideMissingValueImputationHelper(
    input_table = fixture$data,
    metadata_table = fixture$design,
    input_table_sample_id_column = Run,
    sample_id_tbl_sample_id_column = Run,
    replicate_group_column = technical_group,
    protein_id_column = Protein.Group,
    peptide_sequence_column = Stripped.Sequence,
    quantity_to_impute_column = Precursor.Quantity,
    imputed_value_column = Peptide.Imputed,
    hek_string = hek_string,
    proportion_missing_values = proportion_missing_values,
    core_utilisation = NA,
    exclusion_column = exclusion_column,
    return_imputation_result = return_imputation_result
  )
}

test_that("technical-replicate imputation uses distinct runs and an inclusive missing boundary", {
  expect_output(
    result <- runPeptideImputationContractFixture(),
    NA
  )

  value_for <- function(protein, run) {
    result$data$Peptide.Imputed[
      result$data$Protein.Group == protein & result$data$Run %in% run
    ]
  }
  flag_for <- function(protein, run) {
    result$data$is_imputed[
      result$data$Protein.Group == protein & result$data$Run %in% run
    ]
  }

  expect_equal(value_for("P_EXACT", "A2"), 10)
  expect_true(flag_for("P_EXACT", "A2"))
  expect_equal(value_for("P_MEAN", "B3"), 15)
  expect_equal(value_for("P_EQUAL", "B3"), 5)
  expect_equal(value_for("P_HEK", "H2"), 8)
  expect_true(all(is.na(value_for("P_ALL", c("A1", "A2")))))
  expect_true(is.na(value_for("P_OVER", "B2")))
  expect_true(is.na(value_for("P_OVER", "B3")))
  expect_equal(value_for("P_SINGLE", "S1"), 7)
  expect_true(is.na(value_for("P_QC", "QC2")))
  expect_false(any(flag_for("P_QC", c("QC1", "QC2"))))

  equal_support <- result$support_table[
    result$support_table$technical_replicate_group == "B" &
      result$support_table$protein_identity == "P_EQUAL",
    ,
    drop = FALSE
  ]
  expect_identical(equal_support$observed_distinct_runs, 2L)
  expect_equal(equal_support$mean_observed_raw_quantity, 5)
  expect_equal(equal_support$missing_fraction, 1 / 3)
  expect_true(equal_support$eligible_for_imputation)

  exact_support <- result$support_table[
    result$support_table$technical_replicate_group == "A" &
      result$support_table$protein_identity == "P_EXACT",
    ,
    drop = FALSE
  ]
  expect_equal(exact_support$missing_fraction, 0.5)
  expect_true(exact_support$eligible_for_imputation)
  expect_identical(result$summary$eligibility_operator, "<=")
  expect_identical(result$summary$quantity_column, "Precursor.Quantity")
  expect_identical(result$summary$exclusion_source, "design_column:exclude_imputation")
  expect_setequal(result$summary$excluded_runs, c("QC1", "QC2"))
  expect_identical(result$summary$imputed_rows, 4L)
  expect_identical(result$summary$observed_rows, 11L)
  expect_identical(result$summary$missing_not_imputed_rows, 5L)

  observed <- dplyr::inner_join(
    makePeptideImputationContractFixture()$data,
    result$data,
    by = c("Run", "Protein.Group", "Stripped.Sequence"),
    suffix = c(".before", ".after")
  )
  expect_equal(observed$Precursor.Quantity.before, observed$Precursor.Quantity.after)
  expect_equal(
    observed$Peptide.Imputed[!is.na(observed$Precursor.Quantity.before)],
    observed$Precursor.Quantity.before[!is.na(observed$Precursor.Quantity.before)]
  )
  expect_false(any(observed$is_imputed[!is.na(observed$Precursor.Quantity.before)]))
})

test_that("technical-replicate imputation rejects duplicate quantitative identities", {
  fixture <- makePeptideImputationContractFixture()
  fixture$data <- rbind(fixture$data, fixture$data[1, , drop = FALSE])

  expect_error(
    runPeptideImputationContractFixture(fixture),
    "A1/P_EXACT/EXACT",
    fixed = TRUE
  )

  duplicate_design <- makePeptideImputationContractFixture()
  duplicate_design$design <- rbind(
    duplicate_design$design,
    duplicate_design$design[1, , drop = FALSE]
  )
  expect_error(
    runPeptideImputationContractFixture(duplicate_design),
    "duplicate run ID(s): A1",
    fixed = TRUE
  )

  unmatched <- makePeptideImputationContractFixture()
  unmatched$data$Run[1] <- "UNMAPPED"
  expect_error(
    runPeptideImputationContractFixture(unmatched),
    "unmapped input run ID(s): UNMAPPED",
    fixed = TRUE
  )
})

test_that("technical-replicate imputation has explicit legacy exclusion and safe skip semantics", {
  fixture <- makePeptideImputationContractFixture()
  default_result <- runPeptideImputationContractFixture(
    fixture,
    exclusion_column = NULL
  )
  expect_equal(
    default_result$data$Peptide.Imputed[
      default_result$data$Protein.Group == "P_HEK" & default_result$data$Run == "H2"
    ],
    8
  )

  expect_warning(
    legacy_result <- runPeptideImputationContractFixture(
      fixture,
      exclusion_column = NULL,
      hek_string = "HEK"
    ),
    "`hek_string` is deprecated",
    fixed = TRUE
  )
  expect_true(is.na(legacy_result$data$Peptide.Imputed[
    legacy_result$data$Protein.Group == "P_HEK" & legacy_result$data$Run == "H2"
  ]))
  expect_identical(legacy_result$summary$exclusion_source, "deprecated_hek_string")

  singleton <- list(
    data = data.frame(
      Run = "ONLY",
      Protein.Group = "P1",
      Stripped.Sequence = "PEP",
      Precursor.Quantity = NA_real_,
      stringsAsFactors = FALSE
    ),
    design = data.frame(
      Run = "ONLY",
      technical_group = "SINGLE",
      exclude_imputation = FALSE,
      stringsAsFactors = FALSE
    )
  )
  skipped <- runPeptideImputationContractFixture(singleton)
  expect_identical(skipped$summary$status, "skipped")
  expect_identical(
    skipped$summary$skip_reason,
    "no_eligible_technical_replicate_groups"
  )
  expect_identical(nrow(skipped$data), 1L)
  expect_true(is.na(skipped$data$Peptide.Imputed))
  expect_false(skipped$data$is_imputed)
})

test_that("peptide imputation S4 method records raw-scale provenance and no-replicate skips", {
  peptide_object <- newPeptideQcMethodsObject(
    raw_values = c(10, 20, 30, 40),
    norm_values = c(1, 2, 3, 4)
  )
  imputed <- peptideMissingValueImputation(
    peptide_object,
    proportion_missing_values = 0.5
  )

  summary <- imputed@args$peptideMissingValueImputation$imputation_summary
  expect_identical(summary$status, "applied")
  expect_identical(summary$quantity_column, "Precursor.Quantity")
  expect_identical(summary$maximum_missing_fraction, 0.5)
  expect_identical(summary$eligibility_operator, "<=")
  expect_identical(sum(imputed@peptide_data$is_imputed), 4L)

  peptide_object@technical_replicate_id <- NA_character_
  skipped <- peptideMissingValueImputation(peptide_object)
  skipped_summary <- skipped@args$peptideMissingValueImputation$imputation_summary
  expect_identical(skipped_summary$status, "skipped")
  expect_identical(
    skipped_summary$skip_reason,
    "technical_replicate_column_not_declared"
  )
  expect_identical(nrow(skipped@peptide_data), nrow(peptide_object@peptide_data))
  expect_equal(skipped@peptide_data$Peptide.Imputed, peptide_object@peptide_data$Precursor.Quantity)
  expect_false(any(skipped@peptide_data$is_imputed))
})

test_that("peptide imputation S4 config and explicit thresholds execute exactly", {
  peptide_object <- newPeptideQcMethodsObject(
    raw_values = c(10, 20, 30, 40),
    norm_values = c(1, 2, 3, 4)
  )
  peptide_object@args$peptideMissingValueImputation <- list(
    imputed_value_column = "Configured.Imputed",
    proportion_missing_values = 0.25,
    exclusion_column = NULL
  )

  configured <- peptideMissingValueImputation(peptide_object)
  configured_summary <- configured@args$peptideMissingValueImputation$imputation_summary
  expect_identical(configured_summary$maximum_missing_fraction, 0.25)
  expect_identical(configured_summary$imputed_value_column, "Configured.Imputed")
  expect_false(any(configured@peptide_data$is_imputed))

  explicit <- peptideMissingValueImputation(
    peptide_object,
    imputed_value_column = "Explicit.Imputed",
    proportion_missing_values = 0.5
  )
  explicit_summary <- explicit@args$peptideMissingValueImputation$imputation_summary
  expect_identical(explicit_summary$maximum_missing_fraction, 0.5)
  expect_identical(explicit_summary$imputed_value_column, "Explicit.Imputed")
  expect_identical(sum(explicit@peptide_data$is_imputed), 4L)
})

test_that("packaged peptide imputation config has no implicit sample-name exclusion", {
  config_lines <- readLines(testthat::test_path("..", "..", "inst", "config", "config.ini"))
  section_start <- match("[peptideMissingValueImputation]", trimws(config_lines))
  next_section <- which(seq_along(config_lines) > section_start & grepl("^\\[", trimws(config_lines)))[1]
  section <- trimws(config_lines[section_start:(next_section - 1L)])

  expect_identical(section[grepl("^proportion_missing_values=", section)], "proportion_missing_values=0.5")
  expect_identical(section[grepl("^exclusion_column=", section)], "exclusion_column=")
  expect_false(any(grepl("hek_string", section, fixed = TRUE)))
})
