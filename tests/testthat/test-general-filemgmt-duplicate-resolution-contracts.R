library(testthat)

repoRoot <- pkgload::pkg_path()

readTopLevelFunction <- function(path, name, parent = baseenv()) {
  exprs <- parse(file = file.path(repoRoot, path), keep.source = TRUE)

  for (expr in exprs) {
    if (is.call(expr) &&
        identical(expr[[1]], as.name("<-")) &&
        identical(as.character(expr[[2]]), name)) {
      env <- new.env(parent = parent)
      eval(expr, envir = env)
      return(get(name, envir = env, inherits = FALSE))
    }
  }

  NULL
}

readTopLevelFunctionsFromFiles <- function(paths, names, parent = baseenv()) {
  env <- new.env(parent = parent)

  for (path in paths) {
    exprs <- parse(file = file.path(repoRoot, path), keep.source = TRUE)

    for (expr in exprs) {
      if (is.call(expr) &&
          identical(expr[[1]], as.name("<-")) &&
          identical(length(expr), 3L)) {
        expr_name <- as.character(expr[[2]])
        if (expr_name %in% names) {
          eval(expr, envir = env)
        }
      }
    }
  }

  stats::setNames(
    lapply(names, function(name) {
      if (exists(name, envir = env, inherits = FALSE)) {
        get(name, envir = env, inherits = FALSE)
      } else {
        NULL
      }
    }),
    names
  )
}

captureWarningResult <- function(fun, ...) {
  warning_message <- NULL
  result <- withCallingHandlers(
    fun(...),
    warning = function(w) {
      warning_message <<- conditionMessage(w)
      invokeRestart("muffleWarning")
    }
  )

  list(result = result, warning = warning_message)
}

captureErrorMessage <- function(fun, ...) {
  tryCatch(fun(...), error = conditionMessage)
}

captureMessageResult <- function(fun, ...) {
  messages <- character()
  value <- withCallingHandlers(
    fun(...),
    message = function(m) {
      messages <<- c(messages, conditionMessage(m))
      invokeRestart("muffleMessage")
    }
  )

  list(value = value, messages = messages)
}

capturePrintedResult <- function(fun, ...) {
  output <- capture.output(
    value <- fun(...),
    type = "output"
  )

  list(value = value, output = output)
}

runUpdateRuvParameters <- function(fun,
                                   config_list,
                                   best_k,
                                   control_genes_index,
                                   percentage_as_neg_ctrl) {
  withVisible(fun(
    config_list = config_list,
    best_k = best_k,
    control_genes_index = control_genes_index,
    percentage_as_neg_ctrl = percentage_as_neg_ctrl
  ))
}

if (!methods::isClass("GeneralFilemgmtDuplicateResolutionObject")) {
  methods::setClass(
    "GeneralFilemgmtDuplicateResolutionObject",
    slots = c(design_matrix = "data.frame", args = "list")
  )
}

makeDuplicateResolutionObject <- function(design_matrix) {
  methods::new(
    "GeneralFilemgmtDuplicateResolutionObject",
    design_matrix = design_matrix,
    args = list()
  )
}

runUpdateMissingValueParameters <- function(fun, design_matrix, config_list = list(), ...) {
  contract_env <- new.env(parent = emptyenv())
  contract_env$config_list <- config_list
  captured <- captureMessageResult(
    fun,
    makeDuplicateResolutionObject(design_matrix),
    config_list_name = "config_list",
    env = contract_env,
    ...
  )

  list(
    object = captured$value,
    messages = captured$messages,
    config_list = contract_env$config_list
  )
}

runChooseBestProteinAccession <- function(fun,
                                          data_tbl,
                                          seqinr_obj,
                                          protein_id_column = "Protein.Group",
                                          ...) {
  output <- capture.output(
    value <- fun(
      data_tbl = data_tbl,
      protein_id_column = protein_id_column,
      seqinr_obj = seqinr_obj,
      ...
    ),
    type = "output"
  )

  list(value = value, output = output)
}

makePeptideNAContractObject <- function() {
  peptide_data <- data.frame(
    Protein.Ids = c("P1", "P1", "P1", "P1", "P2", "P2", "P2", "P2"),
    Peptide.Sequence = c("PEP1", "PEP1", "PEP1", "PEP1", "PEP2", "PEP2", "PEP2", "PEP2"),
    Run = c("S1", "S2", "S3", "S4", "S1", "S2", "S3", "S4"),
    Precursor.Quantity = c(10, NA, 11, NA, 20, 21, NA, 24),
    Precursor.Normalised = c(10, NA, 11, NA, 20, 21, NA, 24),
    Q.Value = rep(0.001, 8),
    stringsAsFactors = FALSE
  )

  peptide_obj <- PeptideQuantitativeData(
    peptide_data = peptide_data,
    design_matrix = data.frame(
      Run = c("S1", "S2", "S3", "S4"),
      group = c("G1", "G1", "G2", "G2"),
      stringsAsFactors = FALSE
    ),
    sample_id = "Run",
    group_id = "group",
    protein_id_column = "Protein.Ids",
    peptide_sequence_column = "Peptide.Sequence",
    q_value_column = "Q.Value",
    raw_quantity_column = "Precursor.Quantity"
  )

  calcPeptideMatrix(peptide_obj)
}

makeProteinNAContractObject <- function() {
  ProteinQuantitativeData(
    protein_quant_table = data.frame(
      Protein.Ids = c("P1", "P2", "P3", "P4"),
      S1 = c(10, 11, NA, 13),
      S2 = c(NA, 12, 13, 14),
      S3 = c(20, NA, 22, 23),
      S4 = c(NA, 32, 33, 34),
      stringsAsFactors = FALSE
    ),
    design_matrix = data.frame(
      Run = c("S1", "S2", "S3", "S4"),
      group = c("G1", "G1", "G2", "G2"),
      stringsAsFactors = FALSE
    ),
    sample_id = "Run",
    group_id = "group",
    protein_id_column = "Protein.Ids"
  )
}

test_that("DESCRIPTION keeps general filemgmt ahead of general helpers for the canonical helper load order", {
  collate_lines <- strsplit(
    read.dcf(file.path(repoRoot, "DESCRIPTION"), fields = "Collate")[1, 1],
    "\n"
  )[[1]]
  collate_entries <- trimws(gsub("'", "", collate_lines, fixed = TRUE))

  expect_lt(
    match("func_general_filemgmt.R", collate_entries),
    match("func_general_helpers.R", collate_entries)
  )
})

test_that("extract_experiment duplicate implementation is removed from general filemgmt", {
  filemgmt_impl <- readTopLevelFunction("R/func_general_filemgmt.R", "extract_experiment")
  helper_impl <- readTopLevelFunction("R/func_general_helpers.R", "extract_experiment")

  expect_null(filemgmt_impl)
  expect_true(is.function(helper_impl))
})

test_that("package-level extract_experiment matches the canonical helper implementation", {
  helper_impl <- readTopLevelFunction("R/func_general_helpers.R", "extract_experiment")
  inputs <- c(
    "20140602_ffs_expt1_r1_junk",
    "20140603_ffs_expt2_r2_test"
  )
  expect_identical(
    extract_experiment(inputs, mode = "range", start = 2, end = 3),
    helper_impl(inputs, mode = "range", start = 2, end = 3)
  )
  expect_identical(
    unname(extract_experiment(inputs, mode = "range", start = 2, end = 3)),
    c("ffs_expt1", "ffs_expt2")
  )

  expect_identical(
    unname(extract_experiment(inputs, mode = "start")),
    c("20140602", "20140603")
  )
  expect_identical(
    unname(extract_experiment(inputs, mode = "end")),
    c("junk", "test")
  )

  pkg_error <- captureErrorMessage(extract_experiment, "a_b", mode = "range")
  helper_error <- captureErrorMessage(helper_impl, "a_b", mode = "range")

  expect_identical(pkg_error, "End position required for range mode")
  expect_identical(pkg_error, helper_error)

  pkg_warning <- captureWarningResult(
    extract_experiment,
    "a_b",
    mode = "range",
    start = 2,
    end = 5
  )
  helper_warning <- captureWarningResult(
    helper_impl,
    "a_b",
    mode = "range",
    start = 2,
    end = 5
  )

  expect_identical(pkg_warning, helper_warning)
  expect_identical(pkg_warning$warning, "Position out of bounds for string: a_b")
  expect_identical(unname(pkg_warning$result), NA_character_)
})

test_that("updateMissingValueParameters now resolves directly to the canonical helper signature", {
  filemgmt_impl <- readTopLevelFunction("R/func_general_filemgmt.R", "updateMissingValueParameters")
  helper_impl <- readTopLevelFunction("R/func_general_helpers.R", "updateMissingValueParameters")
  helper_formals <- c(
    "theObject",
    "min_reps_per_group",
    "min_groups",
    "function_name",
    "grouping_variable",
    "config_list_name",
    "env"
  )

  expect_null(filemgmt_impl)
  expect_true(is.function(helper_impl))
  expect_identical(names(formals(helper_impl)), helper_formals)
  expect_identical(formals(updateMissingValueParameters), formals(helper_impl))
})

test_that("package-level updateMissingValueParameters matches the canonical helper on default and helper-only paths", {
  helper_impl <- readTopLevelFunction(
    "R/func_general_helpers.R",
    "updateMissingValueParameters",
    parent = environment(updateMissingValueParameters)
  )
  design_matrix <- data.frame(
    sample = paste0("sample_", seq_len(5)),
    group = c("A", "A", "B", "B", "B"),
    condition = c("X", "X", "Y", "Y", "Y"),
    stringsAsFactors = FALSE
  )

  package_default <- runUpdateMissingValueParameters(
    updateMissingValueParameters,
    design_matrix = design_matrix,
    config_list = list(removeRowsWithMissingValuesPercent = list(existing = TRUE)),
    min_reps_per_group = 2,
    min_groups = 2
  )
  helper_default <- runUpdateMissingValueParameters(
    helper_impl,
    design_matrix = design_matrix,
    config_list = list(removeRowsWithMissingValuesPercent = list(existing = TRUE)),
    min_reps_per_group = 2,
    min_groups = 2
  )

  expect_identical(package_default$object@args, helper_default$object@args)
  expect_identical(package_default$config_list, helper_default$config_list)
  expect_identical(package_default$messages, helper_default$messages)
  expect_identical(
    package_default$config_list$removeRowsWithMissingValuesPercent$proteins_intensity_cutoff_percentile,
    1
  )

  package_custom <- runUpdateMissingValueParameters(
    updateMissingValueParameters,
    design_matrix = design_matrix,
    min_reps_per_group = 2,
    min_groups = 1,
    function_name = "peptideIntensityFiltering",
    grouping_variable = "condition"
  )
  helper_custom <- runUpdateMissingValueParameters(
    helper_impl,
    design_matrix = design_matrix,
    min_reps_per_group = 2,
    min_groups = 1,
    function_name = "peptideIntensityFiltering",
    grouping_variable = "condition"
  )

  expect_identical(package_custom$object@args, helper_custom$object@args)
  expect_identical(package_custom$config_list, helper_custom$config_list)
  expect_identical(package_custom$messages, helper_custom$messages)
  expect_identical(
    package_custom$config_list$peptideIntensityFiltering$peptides_intensity_cutoff_percentile,
    1
  )
})

test_that("chooseBestProteinAccession_s3 duplicate implementation is removed from general filemgmt", {
  filemgmt_impl <- readTopLevelFunction(
    "R/func_general_filemgmt.R",
    "chooseBestProteinAccession_s3",
    parent = environment(chooseBestProteinAccession_s3)
  )
  helper_impl <- readTopLevelFunction(
    "R/func_general_helpers.R",
    "chooseBestProteinAccession_s3",
    parent = environment(chooseBestProteinAccession_s3)
  )

  expect_null(filemgmt_impl)
  expect_true(is.function(helper_impl))
  expect_identical(formals(chooseBestProteinAccession_s3), formals(helper_impl))
})

test_that("package-level chooseBestProteinAccession_s3 matches the canonical helper on current fallback paths", {
  helper_impl <- readTopLevelFunction(
    "R/func_general_helpers.R",
    "chooseBestProteinAccession_s3",
    parent = environment(chooseBestProteinAccession_s3)
  )

  default_data <- data.frame(
    Protein.Group = c("P2;P1", "P2", "P3;P4", "P4"),
    sampleA = c(10, 6, 0, 5),
    sampleB = c(0, 3, 4, 6),
    annotation = c("grp", "single", "grp34", "single4"),
    stringsAsFactors = FALSE
  )
  default_seq <- data.frame(
    uniprot_acc = c("P1", "P2", "P3", "P4"),
    length = c(100, 120, 90, 90),
    stringsAsFactors = FALSE
  )
  expected_default <- default_data
  expected_default$sampleA[expected_default$sampleA == 0] <- NA_real_
  expected_default$sampleB[expected_default$sampleB == 0] <- NA_real_

  package_default <- runChooseBestProteinAccession(
    chooseBestProteinAccession_s3,
    data_tbl = default_data,
    seqinr_obj = default_seq
  )
  filemgmt_default <- runChooseBestProteinAccession(
    helper_impl,
    data_tbl = default_data,
    seqinr_obj = default_seq
  )
  helper_default <- runChooseBestProteinAccession(
    helper_impl,
    data_tbl = default_data,
    seqinr_obj = default_seq
  )

  expect_identical(package_default, helper_default)
  expect_identical(package_default, filemgmt_default)
  expect_identical(package_default$value, expected_default)
  expect_true(any(grepl("S3 CLEANUP FATAL ERROR: Error combining results", package_default$output, fixed = TRUE)))
  expect_true(any(grepl("Returning original data", package_default$output, fixed = TRUE)))

  fallback_data <- data.frame(
    Protein.Group = c("X2|X1", "X3|X4"),
    sampleA = c(1, 2),
    sampleB = c(4, 5),
    stringsAsFactors = FALSE
  )
  fallback_seq <- data.frame(
    uniprot_acc = c("X1", "X2", "X3", "X4"),
    stringsAsFactors = FALSE
  )

  package_fallback <- runChooseBestProteinAccession(
    chooseBestProteinAccession_s3,
    data_tbl = fallback_data,
    seqinr_obj = fallback_seq,
    delim = "\\|",
    replace_zero_with_na = FALSE,
    aggregation_method = "sum"
  )
  filemgmt_fallback <- runChooseBestProteinAccession(
    helper_impl,
    data_tbl = fallback_data,
    seqinr_obj = fallback_seq,
    delim = "\\|",
    replace_zero_with_na = FALSE,
    aggregation_method = "sum"
  )
  helper_fallback <- runChooseBestProteinAccession(
    helper_impl,
    data_tbl = fallback_data,
    seqinr_obj = fallback_seq,
    delim = "\\|",
    replace_zero_with_na = FALSE,
    aggregation_method = "sum"
  )

  expect_identical(package_fallback, helper_fallback)
  expect_identical(package_fallback, filemgmt_fallback)
  expect_identical(package_fallback$value, fallback_data)
  expect_true(any(grepl("Creating default length column", package_fallback$output, fixed = TRUE)))
  expect_true(any(grepl("Returning original data", package_fallback$output, fixed = TRUE)))
})

test_that("updateRuvParameters duplicate implementation is removed from general filemgmt while the normalization helper remains canonical", {
  collate_lines <- strsplit(
    read.dcf(file.path(repoRoot, "DESCRIPTION"), fields = "Collate")[1, 1],
    "\n"
  )[[1]]
  collate_entries <- trimws(gsub("'", "", collate_lines, fixed = TRUE))

  filemgmt_impl <- readTopLevelFunction("R/func_general_filemgmt.R", "updateRuvParameters")
  norm_impl <- readTopLevelFunction("R/func_prot_norm_optimization_helpers.R", "updateRuvParameters")

  expect_lt(
    match("func_general_filemgmt.R", collate_entries),
    match("func_prot_norm_optimization_helpers.R", collate_entries)
  )
  expect_null(filemgmt_impl)
  expect_true(is.function(norm_impl))
  expect_identical(formals(updateRuvParameters), formals(norm_impl))
})

test_that("package-level updateRuvParameters matches the canonical normalization helper on current mutation paths", {
  norm_impl <- readTopLevelFunction("R/func_prot_norm_optimization_helpers.R", "updateRuvParameters")

  existing_config <- list(
    ruvParameters = list(
      best_k = 1L,
      num_neg_ctrl = 99L,
      percentage_as_neg_ctrl = 50,
      preserved = "keep"
    ),
    metadata = list(tag = "existing")
  )
  package_existing <- runUpdateRuvParameters(
    updateRuvParameters,
    config_list = existing_config,
    best_k = 4L,
    control_genes_index = c(2L, 4L, 6L),
    percentage_as_neg_ctrl = 12.5
  )
  norm_existing <- runUpdateRuvParameters(
    norm_impl,
    config_list = existing_config,
    best_k = 4L,
    control_genes_index = c(2L, 4L, 6L),
    percentage_as_neg_ctrl = 12.5
  )

  expect_identical(package_existing, norm_existing)
  expect_identical(
    package_existing$value$ruvParameters,
    list(
      best_k = 4L,
      num_neg_ctrl = 3L,
      percentage_as_neg_ctrl = 12.5,
      preserved = "keep"
    )
  )
  expect_identical(package_existing$value$metadata, existing_config$metadata)
  expect_true(package_existing$visible)

  empty_ctrl_config <- list(
    ruvParameters = list(
      best_k = NULL,
      num_neg_ctrl = NULL,
      percentage_as_neg_ctrl = NULL
    ),
    metadata = list(tag = "empty")
  )
  package_empty <- runUpdateRuvParameters(
    updateRuvParameters,
    config_list = empty_ctrl_config,
    best_k = 0L,
    control_genes_index = integer(),
    percentage_as_neg_ctrl = 0
  )
  norm_empty <- runUpdateRuvParameters(
    norm_impl,
    config_list = empty_ctrl_config,
    best_k = 0L,
    control_genes_index = integer(),
    percentage_as_neg_ctrl = 0
  )

  expect_identical(package_empty, norm_empty)
  expect_identical(package_empty$value$ruvParameters$best_k, 0L)
  expect_identical(package_empty$value$ruvParameters$num_neg_ctrl, 0L)
  expect_identical(package_empty$value$ruvParameters$percentage_as_neg_ctrl, 0)
  expect_identical(package_empty$value$metadata, empty_ctrl_config$metadata)
  expect_true(package_empty$visible)
})

test_that("DESCRIPTION keeps general filemgmt ahead of the NA-validation canonical helpers", {
  collate_lines <- strsplit(
    read.dcf(file.path(repoRoot, "DESCRIPTION"), fields = "Collate")[1, 1],
    "\n"
  )[[1]]
  collate_entries <- trimws(gsub("'", "", collate_lines, fixed = TRUE))

  expect_lt(
    match("func_general_filemgmt.R", collate_entries),
    match("func_prot_qc_peptide_support.R", collate_entries)
  )
  expect_lt(
    match("func_general_filemgmt.R", collate_entries),
    match("func_peptide_qc_imputation.R", collate_entries)
  )
  expect_lt(
    match("func_general_filemgmt.R", collate_entries),
    match("func_prot_qc_reporting_helpers.R", collate_entries)
  )
})

test_that("NA-validation duplicates are removed from general filemgmt while canonical package signatures stay aligned", {
  namespace_env <- environment(checkPeptideNAPercentages)
  filemgmt_impls <- readTopLevelFunctionsFromFiles(
    "R/func_general_filemgmt.R",
    c(
      "checkPeptideNAPercentages",
      "validatePostImputationData",
      "getProteinNARecommendations",
      "checkProteinNAPercentages",
      "validatePostImputationProteinData"
    ),
    parent = namespace_env
  )
  peptide_impls <- readTopLevelFunctionsFromFiles(
    c("R/func_prot_qc_peptide_support.R", "R/func_peptide_qc_imputation.R"),
    c("checkPeptideNAPercentages", "validatePostImputationData"),
    parent = namespace_env
  )
  protein_impls <- readTopLevelFunctionsFromFiles(
    "R/func_prot_qc_reporting_helpers.R",
    c(
      "getProteinNARecommendations",
      "checkProteinNAPercentages",
      "validatePostImputationProteinData"
    ),
    parent = namespace_env
  )

  for (name in names(filemgmt_impls)) {
    expect_null(filemgmt_impls[[name]])
  }

  expect_true(is.function(peptide_impls$checkPeptideNAPercentages))
  expect_true(is.function(peptide_impls$validatePostImputationData))
  expect_true(is.function(protein_impls$getProteinNARecommendations))
  expect_true(is.function(protein_impls$checkProteinNAPercentages))
  expect_true(is.function(protein_impls$validatePostImputationProteinData))

  expect_identical(
    formals(checkPeptideNAPercentages),
    formals(peptide_impls$checkPeptideNAPercentages)
  )
  expect_identical(
    formals(validatePostImputationData),
    formals(peptide_impls$validatePostImputationData)
  )
  expect_identical(
    formals(getProteinNARecommendations),
    formals(protein_impls$getProteinNARecommendations)
  )
  expect_identical(
    formals(checkProteinNAPercentages),
    formals(protein_impls$checkProteinNAPercentages)
  )
  expect_identical(
    formals(validatePostImputationProteinData),
    formals(protein_impls$validatePostImputationProteinData)
  )
})

test_that("package-level peptide NA helpers match the canonical helper implementations on the current analysis and validation paths", {
  namespace_env <- environment(checkPeptideNAPercentages)
  canonical_impls <- readTopLevelFunctionsFromFiles(
    c("R/func_prot_qc_peptide_support.R", "R/func_peptide_qc_imputation.R"),
    c("checkPeptideNAPercentages", "validatePostImputationData"),
    parent = namespace_env
  )
  peptide_obj <- makePeptideNAContractObject()

  package_check <- capturePrintedResult(checkPeptideNAPercentages, peptide_obj, verbose = TRUE)
  canonical_check <- capturePrintedResult(canonical_impls$checkPeptideNAPercentages, peptide_obj, verbose = TRUE)

  expect_identical(package_check, canonical_check)
  expect_identical(package_check$value$total_na_percent, 37.5)
  expect_true(any(grepl("Peptide Data Missing Value Analysis", package_check$output, fixed = TRUE)))

  package_validate <- capturePrintedResult(
    validatePostImputationData,
    peptide_obj,
    expected_na_percent = 0,
    tolerance = 0.1
  )
  canonical_validate <- capturePrintedResult(
    canonical_impls$validatePostImputationData,
    peptide_obj,
    expected_na_percent = 0,
    tolerance = 0.1
  )

  expect_identical(package_validate, canonical_validate)
  expect_false(package_validate$value$is_valid)
  expect_identical(package_validate$value$actual_na_percent, 37.5)
  expect_true(any(grepl("VALIDATION FAILED", package_validate$output, fixed = TRUE)))
  expect_true(any(grepl("High NA percentage suggests imputation problems", package_validate$output, fixed = TRUE)))
})

test_that("package-level protein NA helpers match the canonical helper implementations on the current recommendation and validation paths", {
  namespace_env <- environment(checkProteinNAPercentages)
  canonical_impls <- readTopLevelFunctionsFromFiles(
    "R/func_prot_qc_reporting_helpers.R",
    c(
      "getProteinNARecommendations",
      "checkProteinNAPercentages",
      "validatePostImputationProteinData"
    ),
    parent = namespace_env
  )
  protein_obj <- makeProteinNAContractObject()

  package_check <- capturePrintedResult(checkProteinNAPercentages, protein_obj, verbose = TRUE)
  canonical_check <- capturePrintedResult(canonical_impls$checkProteinNAPercentages, protein_obj, verbose = TRUE)

  expect_identical(package_check, canonical_check)
  expect_identical(package_check$value$total_na_percent, 25)
  expect_true(any(grepl("Protein Data Missing Value Analysis", package_check$output, fixed = TRUE)))

  package_recommend <- capturePrintedResult(
    getProteinNARecommendations,
    protein_obj,
    include_code = TRUE
  )
  canonical_recommend <- capturePrintedResult(
    canonical_impls$getProteinNARecommendations,
    protein_obj,
    include_code = TRUE
  )

  expect_identical(package_recommend, canonical_recommend)
  expect_identical(package_recommend$value$primary_recommendation, "imputation_or_filtering")
  expect_true(any(grepl("Consider Protein-Level Imputation", package_recommend$output, fixed = TRUE)))

  package_validate <- capturePrintedResult(
    validatePostImputationProteinData,
    protein_obj,
    expected_na_percent = NULL,
    tolerance = 10
  )
  canonical_validate <- capturePrintedResult(
    canonical_impls$validatePostImputationProteinData,
    protein_obj,
    expected_na_percent = NULL,
    tolerance = 10
  )

  expect_identical(package_validate, canonical_validate)
  expect_true(package_validate$value$is_valid)
  expect_identical(package_validate$value$actual_na_percent, 25)
  expect_identical(package_validate$value$expected_na_percent, 35)
  expect_true(any(grepl("Using default expected NA% of 35.0% for protein data", package_validate$output, fixed = TRUE)))
})
