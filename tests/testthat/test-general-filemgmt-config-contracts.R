# fidelity-coverage-compare: shared
library(testthat)

makeFunctionWithOverrides <- function(fun, replacements) {
  fun_override <- fun
  environment(fun_override) <- list2env(replacements, parent = environment(fun))
  fun_override
}

if (!methods::isClass("mockGeneralFileMgmtConfigCarrier")) {
  methods::setClass("mockGeneralFileMgmtConfigCarrier", slots = c(args = "list"))
}

test_that("readConfigFile normalizes parsed config entries and forwards file_type", {
  parser_call <- NULL
  new_cluster_arg <- NULL
  cluster_library_call <- NULL
  cluster_token <- structure(list(id = "cluster-1"), class = "mock_cluster")

  local_mocked_bindings(
    read.config = function(file, file.type) {
      parser_call <<- list(file = file, file_type = file.type)

      list(
        globalParameters = list(
          number_of_cpus = "4",
          plots_format = "pdf,png"
        ),
        srlQvalueProteotypicPeptideClean = list(
          input_matrix_column_ids = "Run, Precursor.Id, , Intensity",
          qvalue_threshold = "0.01",
          global_qvalue_threshold = "0.05",
          choose_only_proteotypic_peptide = "1"
        ),
        peptideIntensityFiltering = list(
          peptides_intensity_cutoff_percentile = "0.2",
          peptides_proportion_of_samples_below_cutoff = "0.6"
        ),
        filterMinNumPeptidesPerProtein = list(
          peptides_per_protein_cutoff = "2",
          peptidoforms_per_protein_cutoff = "3"
        ),
        filterMinNumPeptidesPerSample = list(
          peptides_per_sample_cutoff = "3",
          inclusion_list = "control,treated"
        ),
        peptideMissingValueImputation = list(
          proportion_missing_values = "0.15"
        ),
        removeRowsWithMissingValuesPercent = list(
          groupwise_percentage_cutoff = "0.25",
          max_groups_percentage_cutoff = "0.5",
          proteins_intensity_cutoff_percentile = "0.1"
        ),
        ruvIII_C_Varying = list(
          ruv_number_k = "7"
        ),
        plotRle = list(
          yaxis_limit = "1.5,2.5"
        ),
        deAnalysisParameters = list(
          plots_format = "pdf,png",
          formula_string = "~ 0 + group",
          da_q_val_thresh = "0.2",
          eBayes_trend = "TRUE",
          eBayes_robust = "false"
        )
      )
    },
    new_cluster = function(n) {
      new_cluster_arg <<- n
      cluster_token
    },
    cluster_library = function(cluster, packages) {
      cluster_library_call <<- list(cluster = cluster, packages = packages)
      invisible(NULL)
    },
    .env = environment(readConfigFile)
  )

  supports_file_type <- "file_type" %in% names(formals(readConfigFile))
  read_args <- list(file = "/tmp/mock-config.ini")
  if (supports_file_type) {
    read_args$file_type <- "cfg"
  }

  config_list <- suppressMessages(do.call(readConfigFile, read_args))

  expected_file_type <- if (supports_file_type) "cfg" else "ini"
  expect_identical(parser_call, list(file = "/tmp/mock-config.ini", file_type = expected_file_type))
  expect_identical(new_cluster_arg, "4")
  expect_identical(cluster_library_call$cluster, cluster_token)
  expect_identical(
    cluster_library_call$packages,
    c("tidyverse", "glue", "rlang", "lazyeval")
  )

  expect_identical(config_list$globalParameters$plots_format, c("pdf", "png"))
  expect_identical(
    config_list$rollUpPrecursorToPeptide$core_utilisation,
    cluster_token
  )
  expect_identical(
    config_list$removeProteinsWithOnlyOneReplicate$core_utilisation,
    cluster_token
  )
  expect_identical(
    config_list$srlQvalueProteotypicPeptideClean$input_matrix_column_ids,
    c("Run", "Precursor.Id", "Intensity")
  )
  expect_equal(config_list$srlQvalueProteotypicPeptideClean$qvalue_threshold, 0.01)
  expect_equal(config_list$srlQvalueProteotypicPeptideClean$global_qvalue_threshold, 0.05)
  expect_equal(config_list$srlQvalueProteotypicPeptideClean$choose_only_proteotypic_peptide, 1)
  expect_equal(config_list$peptideIntensityFiltering$peptides_intensity_cutoff_percentile, 0.2)
  expect_equal(config_list$peptideIntensityFiltering$peptides_proportion_of_samples_below_cutoff, 0.6)
  expect_equal(config_list$filterMinNumPeptidesPerProtein$peptides_per_protein_cutoff, 2)
  expect_equal(config_list$filterMinNumPeptidesPerProtein$peptidoforms_per_protein_cutoff, 3)
  expect_equal(config_list$filterMinNumPeptidesPerSample$peptides_per_sample_cutoff, 3)
  expect_identical(
    config_list$filterMinNumPeptidesPerSample$inclusion_list,
    c("control", "treated")
  )
  expect_equal(config_list$peptideMissingValueImputation$proportion_missing_values, 0.15)
  expect_equal(config_list$removeRowsWithMissingValuesPercent$groupwise_percentage_cutoff, 0.25)
  expect_equal(config_list$removeRowsWithMissingValuesPercent$max_groups_percentage_cutoff, 0.5)
  expect_equal(config_list$removeRowsWithMissingValuesPercent$proteins_intensity_cutoff_percentile, 0.1)
  expect_equal(config_list$ruvIII_C_Varying$ruv_number_k, 7)
  expect_equal(config_list$plotRle$yaxis_limit, c(1.5, 2.5))
  expect_identical(config_list$deAnalysisParameters$plots_format, c("pdf", "png"))
  expect_false(config_list$deAnalysisParameters$lfc_cutoff)
  expect_equal(config_list$deAnalysisParameters$treat_lfc_cutoff, 0)
  expect_equal(config_list$deAnalysisParameters$da_q_val_thresh, 0.2)
  expect_true(config_list$deAnalysisParameters$eBayes_trend)
  expect_false(config_list$deAnalysisParameters$eBayes_robust)
})

test_that("readConfigFileSection updates either a full section or one parameter", {
  read_config_calls <- list()

  section_reader <- makeFunctionWithOverrides(
    readConfigFileSection,
    list(
      readConfigFile = function(file, file_type = "ini") {
        read_config_calls[[length(read_config_calls) + 1]] <<- list(
          file = file,
          file_type = file_type
        )

        list(
          deAnalysisParameters = list(
            formula_string = "~ 0 + group",
            eBayes_trend = TRUE
          ),
          plotRle = list(
            yaxis_limit = c(1, 2)
          )
        )
      }
    )
  )

  carrier <- methods::new(
    "mockGeneralFileMgmtConfigCarrier",
    args = list(
      deAnalysisParameters = list(
        formula_string = "~ old",
        keep_me = "present"
      ),
      plotRle = list(
        yaxis_limit = c(0, 1),
        keep_me = "still-here"
      )
    )
  )

  updated_section <- section_reader(
    carrier,
    file = "/tmp/mock-config.ini",
    function_name = "deAnalysisParameters"
  )
  updated_parameter <- section_reader(
    carrier,
    file = "/tmp/mock-config.ini",
    function_name = "plotRle",
    parameter_name = "yaxis_limit"
  )

  expect_identical(
    updated_section@args$deAnalysisParameters,
    list(formula_string = "~ 0 + group", eBayes_trend = TRUE)
  )
  expect_identical(updated_section@args$plotRle, carrier@args$plotRle)

  expect_identical(updated_parameter@args$plotRle$yaxis_limit, c(1, 2))
  expect_identical(updated_parameter@args$plotRle$keep_me, "still-here")
  expect_identical(
    updated_parameter@args$deAnalysisParameters,
    carrier@args$deAnalysisParameters
  )

  expect_identical(
    read_config_calls,
    list(
      list(file = "/tmp/mock-config.ini", file_type = "ini"),
      list(file = "/tmp/mock-config.ini", file_type = "ini")
    )
  )
})
