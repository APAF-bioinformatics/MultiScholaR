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
    Global.Q.Value = c(0.011, 0.012, 0.013, 0.014),
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
  method_env <- environment(methods::selectMethod("peptideIntensityFiltering", "PeptideQuantitativeData")@.Data)

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
        global_qvalue_threshold = 0.05,
        choose_only_proteotypic_peptide = 1,
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
                                                   inclusion_list = NULL) {
      appendMethodCall(
        log_env,
        "sample_calls",
        list(
          input_table = input_table,
          peptides_per_sample_cutoff = peptides_per_sample_cutoff,
          sample_id_column = sample_id_column,
          core_utilisation = core_utilisation,
          inclusion_list = inclusion_list
        )
      )
      input_table[1, , drop = FALSE]
    },
    srlQvalueProteotypicPeptideCleanHelper = function(...) {
      args <- list(...)
      appendMethodCall(log_env, "srl_calls", args)
      args$input_table[1, , drop = FALSE]
    },
    .env = method_env
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

  expect_equal(log_env$intensity_calls[[1]]$min_peptide_intensity_threshold, 25)
  expect_equal(log_env$intensity_calls[[2]]$min_peptide_intensity_threshold, 0)
  expect_equal(log_env$missing_calls[[1]]$abundance_threshold, 3)
  expect_equal(log_env$missing_calls[[2]]$abundance_threshold, 0)
  expect_equal(nrow(valid_intensity@peptide_data), 2)
  expect_equal(nrow(valid_missing@peptide_data), 1)
  expect_true(length(log_env$clean_calls) >= 4)
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
    c("Run", "CustomColumn", "Protein.Ids", "Stripped.Sequence")
  )
  expect_identical(log_env$srl_calls[[1]]$global_qvalue_threshold, 0.05)
  expect_identical(log_env$srl_calls[[1]]$qvalue_threshold, 0.01)
  expect_identical(log_env$srl_calls[[1]]$choose_only_proteotypic_peptide, 1)
  expect_true(length(log_env$updated_params) >= 10)
})
