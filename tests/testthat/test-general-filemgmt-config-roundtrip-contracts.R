library(testthat)

restoreGlobalBinding <- function(name, had_value, old_value = NULL) {
  if (had_value) {
    assign(name, old_value, envir = .GlobalEnv)
  } else if (exists(name, envir = .GlobalEnv, inherits = FALSE)) {
    rm(list = name, envir = .GlobalEnv)
  }
}

if (!methods::isClass("mockGeneralFileMgmtRoundtripCarrier")) {
  methods::setClass("mockGeneralFileMgmtRoundtripCarrier", slots = c(args = "list"))
}

test_that("formatConfigList omits internal sources and preserves nested display structure", {
  config_list <- list(
    peptideIntensityFiltering = list(
      peptides_proportion_of_samples_below_cutoff = 0.7,
      keep_me = "yes"
    ),
    nested = list(items = list("first", "second")),
    internal_workflow_source_dir = "/tmp/skip-me",
    seqinr_obj = data.frame(id = 1:2)
  )

  formatted <- suppressMessages(formatConfigList(config_list))

  expect_true(any(grepl("^Nested:$", formatted)))
  expect_true(any(grepl("Peptides Proportion of Samples Below Cutoff: 0.7", formatted, fixed = TRUE)))
  expect_true(any(grepl("Keep Me: yes", formatted, fixed = TRUE)))
  expect_true(any(grepl("- first", formatted, fixed = TRUE)))
  expect_true(any(grepl("- second", formatted, fixed = TRUE)))
  expect_false(any(grepl("internal_workflow_source_dir", formatted, fixed = TRUE)))
  expect_false(any(grepl("skip-me", formatted, fixed = TRUE)))
  expect_false(any(grepl("seqinr_obj", formatted, fixed = TRUE)))
})

test_that("config round-trip helpers keep updated values aligned across env and file writers", {
  config_env <- new.env(parent = emptyenv())
  config_env$config_list <- list(
    peptideIntensityFiltering = list(
      peptides_proportion_of_samples_below_cutoff = 0.2,
      keep_me = "yes"
    ),
    internal_workflow_source_dir = "/tmp/skip-me"
  )

  carrier <- methods::new(
    "mockGeneralFileMgmtRoundtripCarrier",
    args = list(
      peptideIntensityFiltering = list(
        peptides_proportion_of_samples_below_cutoff = 0.2,
        keep_me = "yes"
      )
    )
  )

  updated <- suppressMessages(updateConfigParameter(
    theObject = carrier,
    function_name = "peptideIntensityFiltering",
    parameter_name = "peptides_proportion_of_samples_below_cutoff",
    new_value = 0.7,
    config_list_name = "config_list",
    env = config_env
  ))

  expect_equal(
    updated@args$peptideIntensityFiltering$peptides_proportion_of_samples_below_cutoff,
    0.7
  )
  expect_equal(
    config_env$config_list$peptideIntensityFiltering$peptides_proportion_of_samples_below_cutoff,
    0.7
  )

  had_global_config <- exists("config_list", envir = .GlobalEnv, inherits = FALSE)
  old_global_config <- if (had_global_config) get("config_list", envir = .GlobalEnv)
  on.exit(
    restoreGlobalBinding("config_list", had_global_config, old_global_config),
    add = TRUE
  )

  assign("config_list", get("config_list", envir = config_env), envir = .GlobalEnv)

  study_dir <- tempfile("general-filemgmt-study-")
  dir.create(study_dir)
  study_file <- suppressMessages(suppressWarnings(createStudyParametersFile(
    workflow_name = "proteomics",
    source_dir_path = study_dir,
    contrasts_tbl = data.frame(contrasts = "A-B", stringsAsFactors = FALSE),
    config_list_name = "config_list",
    env = config_env
  )))
  study_lines <- readLines(study_file)

  expect_identical(basename(study_file), "study_parameters.txt")
  expect_true(any(grepl("^Configuration Parameters:$", study_lines)))
  expect_true(any(grepl("Peptides Proportion of Samples Below Cutoff: 0.7", study_lines, fixed = TRUE)))
  expect_true(any(grepl("Keep Me: yes", study_lines, fixed = TRUE)))
  expect_true(any(grepl("A-B", study_lines, fixed = TRUE)))
  expect_false(any(grepl("internal_workflow_source_dir", study_lines, fixed = TRUE)))

  workflow_dir <- tempfile("general-filemgmt-workflow-")
  dir.create(workflow_dir)
  workflow_file <- suppressMessages(createWorkflowArgsFromConfig(
    workflow_name = "proteomics",
    source_dir_path = workflow_dir,
    final_s4_object = updated,
    contrasts_tbl = data.frame(contrasts = "A-B", stringsAsFactors = FALSE)
  ))
  workflow_lines <- readLines(workflow_file)

  expect_identical(basename(workflow_file), "study_parameters.txt")
  expect_true(any(grepl("^Version Information:$", workflow_lines)))
  expect_true(any(grepl("^Additional Configuration Parameters:$", workflow_lines)))
  expect_true(any(grepl("^\\[peptideIntensityFiltering\\]$", workflow_lines)))
  expect_true(any(grepl("peptides_proportion_of_samples_below_cutoff = 0.7", workflow_lines, fixed = TRUE)))
  expect_true(any(grepl("Peptides Proportion of Samples Below Cutoff: 0.7", workflow_lines, fixed = TRUE)))
  expect_true(any(grepl("A-B", workflow_lines, fixed = TRUE)))
  expect_false(any(grepl("internal_workflow_source_dir", workflow_lines, fixed = TRUE)))
})
