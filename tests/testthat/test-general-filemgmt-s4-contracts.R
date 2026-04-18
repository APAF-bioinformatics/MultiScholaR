library(testthat)

restoreGlobalBinding <- function(name, had_value, old_value = NULL) {
  if (had_value) {
    assign(name, old_value, envir = .GlobalEnv)
  } else if (exists(name, envir = .GlobalEnv, inherits = FALSE)) {
    rm(list = name, envir = .GlobalEnv)
  }
}

makeMetabolomicsWorkflowObject <- function() {
  design_matrix <- data.frame(
    Sample_ID = c("S1", "S2", "S3"),
    group = c("A", "A", "B"),
    stringsAsFactors = FALSE
  )

  assay_pos <- data.frame(
    database_identifier = c("M1", "M2"),
    metabolite_identification = c("Met1", "Met2"),
    S1 = c(10, 20),
    S2 = c(11, 21),
    S3 = c(12, 22),
    check.names = FALSE
  )
  assay_neg <- data.frame(
    database_identifier = c("M3", "M4", "M5"),
    metabolite_identification = c("Met3", "Met4", "Met5"),
    S1 = c(30, 40, 50),
    S2 = c(31, 41, 51),
    S3 = c(32, 42, 52),
    check.names = FALSE
  )

  createMetaboliteAssayData(
    metabolite_data = list(assay_pos = assay_pos, assay_neg = assay_neg),
    design_matrix = design_matrix,
    metabolite_id_column = "database_identifier",
    annotation_id_column = "metabolite_identification",
    sample_id = "Sample_ID",
    group_id = "group",
    database_identifier_type = "HMDB",
    internal_standard_regex = "^ISTD",
    args = list(
      ITSDNormalization = list(
        applied = TRUE,
        method_type = "median",
        itsd_aggregation = "sum",
        itsd_pattern_columns = c("ISTD_1", "ISTD_2"),
        removed_itsd = TRUE,
        timestamp = as.POSIXct("2026-04-18 09:00:00", tz = "UTC"),
        itsd_features_per_assay = list(
          assay_pos = c("ISTD_A1", "ISTD_A2"),
          assay_neg = c("ISTD_B1")
        ),
        itsd_counts_per_assay = list(
          assay_pos = 2L,
          assay_neg = 1L
        )
      ),
      log_transformed = TRUE,
      log_transform_offset = 1,
      normalisation_method = "PQN",
      ruv_grouping_variable = "group",
      ruv_number_k = list(assay_pos = 2L, assay_neg = 3L),
      ctrl = list(
        assay_pos = c(TRUE, FALSE),
        assay_neg = c(TRUE, TRUE, FALSE)
      )
    )
  )
}

makeLipidomicsWorkflowObject <- function() {
  design_matrix <- data.frame(
    Sample_ID = c("S1", "S2", "S3"),
    group = c("A", "A", "B"),
    stringsAsFactors = FALSE
  )

  assay_data <- data.frame(
    database_identifier = c("L1", "L2", "L3"),
    lipid_identification = c("Lipid1", "Lipid2", "Lipid3"),
    S1 = c(100, 200, 300),
    S2 = c(101, 201, 301),
    S3 = c(102, 202, 302),
    check.names = FALSE
  )

  createLipidomicsAssayData(
    lipid_data = list(assay_main = assay_data),
    design_matrix = design_matrix,
    lipid_id_column = "database_identifier",
    annotation_id_column = "lipid_identification",
    sample_id = "Sample_ID",
    group_id = "group",
    database_identifier_type = "LipidMaps",
    internal_standard_regex = "^LIPID_ISTD",
    args = list(
      ITSDNormalization = list(applied = FALSE),
      log_transformed = FALSE,
      normalisation_method = "median",
      ruv_grouping_variable = "group",
      ruv_number_k = 4L,
      ctrl = c(TRUE, FALSE, TRUE)
    )
  )
}

test_that("createWorkflowArgsFromConfig serializes metabolomics S4 helper output", {
  had_global_config <- exists("config_list", envir = .GlobalEnv, inherits = FALSE)
  old_global_config <- if (had_global_config) get("config_list", envir = .GlobalEnv)
  on.exit(
    restoreGlobalBinding("config_list", had_global_config, old_global_config),
    add = TRUE
  )

  assign("config_list", list(), envir = .GlobalEnv)

  workflow_dir <- tempfile("general-filemgmt-metab-workflow-")
  dir.create(workflow_dir)
  workflow_file <- NULL

  invisible(capture.output(
    workflow_file <- createWorkflowArgsFromConfig(
      workflow_name = "metabolomics",
      source_dir_path = workflow_dir,
      final_s4_object = makeMetabolomicsWorkflowObject(),
      contrasts_tbl = data.frame(contrast = "A_vs_B", stringsAsFactors = FALSE)
    ),
    type = "output"
  ))
  workflow_lines <- readLines(workflow_file)

  expect_identical(basename(workflow_file), "study_parameters.txt")
  expect_true(any(grepl("^Assay Information:$", workflow_lines)))
  expect_true(any(grepl("* Number of Assays: 2", workflow_lines, fixed = TRUE)))
  expect_true(any(grepl("* Total Metabolites: 5", workflow_lines, fixed = TRUE)))
  expect_true(any(grepl("assay_pos: 2", workflow_lines, fixed = TRUE)))
  expect_true(any(grepl("assay_neg: 3", workflow_lines, fixed = TRUE)))
  expect_true(any(grepl("^ITSD Normalization:$", workflow_lines)))
  expect_true(any(grepl("* Applied: Yes", workflow_lines, fixed = TRUE)))
  expect_true(any(grepl("assay_pos (2 features): ISTD_A1, ISTD_A2", workflow_lines, fixed = TRUE)))
  expect_true(any(grepl("^Between-Sample Normalization:$", workflow_lines)))
  expect_true(any(grepl("* Method: PQN", workflow_lines, fixed = TRUE)))
  expect_true(any(grepl("^Metabolomics Metadata:$", workflow_lines)))
  expect_true(any(grepl("* Number of Groups: 2", workflow_lines, fixed = TRUE)))
  expect_true(any(grepl("^RUV-III Batch Correction:$", workflow_lines)))
  expect_true(any(grepl("assay_pos: k=2, controls=1", workflow_lines, fixed = TRUE)))
  expect_true(any(grepl("assay_neg: k=3, controls=2", workflow_lines, fixed = TRUE)))
})

test_that("createWorkflowArgsFromConfig serializes lipidomics S4 helper output", {
  had_global_config <- exists("config_list", envir = .GlobalEnv, inherits = FALSE)
  old_global_config <- if (had_global_config) get("config_list", envir = .GlobalEnv)
  on.exit(
    restoreGlobalBinding("config_list", had_global_config, old_global_config),
    add = TRUE
  )

  assign("config_list", list(), envir = .GlobalEnv)

  workflow_dir <- tempfile("general-filemgmt-lipid-workflow-")
  dir.create(workflow_dir)
  workflow_file <- NULL

  invisible(capture.output(
    workflow_file <- createWorkflowArgsFromConfig(
      workflow_name = "lipidomics",
      source_dir_path = workflow_dir,
      final_s4_object = makeLipidomicsWorkflowObject(),
      contrasts_tbl = data.frame(contrast = "A_vs_B", stringsAsFactors = FALSE)
    ),
    type = "output"
  ))
  workflow_lines <- readLines(workflow_file)

  expect_identical(basename(workflow_file), "study_parameters.txt")
  expect_true(any(grepl("^Assay Information:$", workflow_lines)))
  expect_true(any(grepl("* Number of Assays: 1", workflow_lines, fixed = TRUE)))
  expect_true(any(grepl("* Total Lipids: 3", workflow_lines, fixed = TRUE)))
  expect_true(any(grepl("assay_main: 3", workflow_lines, fixed = TRUE)))
  expect_true(any(grepl("^ITSD Normalization:$", workflow_lines)))
  expect_true(any(grepl("* Applied: No", workflow_lines, fixed = TRUE)))
  expect_true(any(grepl("^Log Transformation:$", workflow_lines)))
  expect_true(any(grepl("^Lipidomics Metadata:$", workflow_lines)))
  expect_true(any(grepl("* Database Identifier Type: LipidMaps", workflow_lines, fixed = TRUE)))
  expect_true(any(grepl("* Internal Standard Regex: ^LIPID_ISTD", workflow_lines, fixed = TRUE)))
  expect_true(any(grepl("^RUV-III Batch Correction:$", workflow_lines)))
  expect_true(any(grepl("* k value: 4", workflow_lines, fixed = TRUE)))
  expect_true(any(grepl("* Control features: 2", workflow_lines, fixed = TRUE)))
})
