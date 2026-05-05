library(testthat)
library(shiny)

test_that("MCI-004.1 DIA-NN matrix preserves peptide/protein/run contracts and rejects malformed schemas", {
  canonical <- module_ci_prot_make_diann(variant = "canonical")
  imported <- suppressMessages(importDIANNData(canonical))
  module_ci_prot_assert_import_contract(
    imported,
    expected_data_type = "peptide",
    expected_samples = c("WT_1", "WT_2", "KO_1", "KO_2"),
    expected_features = c("P001", "P002", "P003"),
    expected_quantity_col = "Precursor.Normalised"
  )
  expect_true("Extra.Metadata" %in% names(imported$data))
  expect_identical(imported$column_mapping$peptide_col, "Stripped.Sequence")
  expect_setequal(unique(imported$data$Stripped.Sequence), c("PEP1", "PEP2", "PEP3", "PEP4"))

  duplicate_proteins <- module_ci_prot_make_diann(variant = "duplicate_proteins")
  duplicate_import <- suppressMessages(importDIANNData(duplicate_proteins))
  expect_setequal(module_ci_prot_feature_values(duplicate_import), c("P001", "P002"))
  expect_true(any(duplicated(duplicate_import$data$Protein.Group)))

  zero_na <- module_ci_prot_make_diann(variant = "zero_na")
  zero_import <- suppressMessages(importDIANNData(zero_na, use_precursor_norm = FALSE))
  expect_identical(zero_import$column_mapping$quantity_col, "Precursor.Quantity")
  expect_true(any(zero_import$data$Precursor.Quantity == 0, na.rm = TRUE))
  expect_true(any(is.na(zero_import$data$Precursor.Quantity)))

  bad_names <- module_ci_prot_make_diann(variant = "bad_run_names")
  bad_import <- suppressMessages(importDIANNData(bad_names))
  workflow_data <- module_ci_prot_workflow_data()
  applyProtImportResultToWorkflow(
    workflowData = workflow_data,
    dataImportResult = bad_import,
    format = "diann",
    fastaPath = NULL,
    sanitizeNames = TRUE,
    logInfo = function(...) NULL,
    showNotification = function(...) NULL
  )
  expected_clean <- janitor::make_clean_names(c("123 Sample-A", "WT sample 2", "KO/sample 1", "KO sample 2"))
  module_ci_assert_sample_identity(expected_clean, unique(workflow_data$data_tbl$Run))

  missing_required <- module_ci_prot_make_diann(variant = "missing_required")
  expect_error(
    suppressMessages(importDIANNData(missing_required)),
    "Missing required DIA-NN columns: Run",
    fixed = TRUE
  )

  no_fasta_state <- module_ci_prot_apply_import_state(imported, "diann", fasta_path = NULL)
  module_ci_prot_assert_workflow_state(no_fasta_state$workflow_data, "diann")
  expect_true("cp01:raw_imported" %in% no_fasta_state$checkpoints)
})

test_that("MCI-004.2 Proteome Discoverer TMT matrix covers real names, simplified names, collisions, shared peptides, and all-NA channels", {
  canonical <- module_ci_prot_make_tmt(variant = "canonical")
  imported <- suppressMessages(suppressWarnings(importProteomeDiscovererTMTData(canonical)))
  module_ci_prot_assert_import_contract(
    imported,
    expected_data_type = "protein",
    expected_samples = c("126_WT_1", "127N_WT_2", "128C_KO_1", "129N_KO_2"),
    expected_features = c("P001", "P002", "P001;P002"),
    expected_quantity_col = "Abundance"
  )
  expect_true("Annotated.Sequence" %in% names(imported$data))
  expect_true("P001;P002" %in% imported$data$Protein.Ids)

  simplified <- module_ci_prot_make_tmt(variant = "simplified")
  simplified_import <- suppressMessages(suppressWarnings(importProteomeDiscovererTMTData(simplified)))
  module_ci_assert_sample_identity(
    c("126", "127N", "128C", "129N"),
    module_ci_prot_run_values(simplified_import)
  )

  all_na <- module_ci_prot_make_tmt(variant = "all_na_channel")
  all_na_import <- suppressMessages(suppressWarnings(importProteomeDiscovererTMTData(all_na)))
  expect_true(all(is.na(all_na_import$data$Abundance[all_na_import$data$Run == "129N_KO_2"])))

  collision <- module_ci_prot_make_tmt(variant = "collision")
  expect_error(
    suppressMessages(suppressWarnings(importProteomeDiscovererTMTData(collision))),
    "Duplicate TMT reporter channel/sample names after normalization",
    fixed = TRUE
  )

  missing_accession <- module_ci_prot_make_tmt(variant = "missing_accession")
  expect_error(
    suppressMessages(suppressWarnings(importProteomeDiscovererTMTData(missing_accession))),
    "Required column 'Accession' not found",
    fixed = TRUE
  )

  state <- module_ci_prot_apply_import_state(imported, "pd_tmt")
  module_ci_prot_assert_workflow_state(state$workflow_data, "pd_tmt")
})

test_that("MCI-004.3 MaxQuant LFQ matrix covers LFQ detection, flags, missing genes, duplicate accessions, and zero coercion", {
  canonical <- module_ci_prot_make_maxquant(variant = "canonical")
  imported <- suppressWarnings(importMaxQuantData(canonical))
  module_ci_prot_assert_import_contract(
    imported,
    expected_data_type = "protein",
    expected_samples = c("WT_1", "WT_2", "KO_1", "KO_2"),
    expected_features = c("P001", "P002;P003"),
    expected_quantity_col = "Intensity"
  )
  expect_false(any(grepl("^CON__|^REV__", imported$data$Protein.Ids)))
  expect_true(any(imported$data$Intensity == 0))
  expect_true("Gene.names" %in% names(imported$data))

  duplicates <- module_ci_prot_make_maxquant(variant = "duplicate_accessions")
  duplicate_import <- suppressWarnings(importMaxQuantData(duplicates, filter_contaminants = FALSE))
  expect_true(any(duplicated(duplicate_import$data$Protein.Ids)))

  no_intensity <- module_ci_prot_make_maxquant(variant = "missing_intensity")
  expect_error(
    importMaxQuantData(no_intensity),
    "No MaxQuant intensity columns found",
    fixed = TRUE
  )

  no_protein <- module_ci_prot_make_maxquant(variant = "missing_protein")
  expect_error(
    importMaxQuantData(no_protein),
    "Required MaxQuant protein identifier column not found",
    fixed = TRUE
  )

  state <- module_ci_prot_apply_import_state(imported, "maxquant")
  module_ci_prot_assert_workflow_state(state$workflow_data, "maxquant")
})

test_that("MCI-004.4 FragPipe LFQ matrix covers MaxLFQ columns, alternate names, invalid schemas, and LFQ contract equivalence", {
  canonical <- module_ci_prot_make_fragpipe(variant = "canonical")
  imported <- suppressMessages(importFragPipeData(canonical))
  module_ci_prot_assert_import_contract(
    imported,
    expected_data_type = "protein",
    expected_samples = c("WT A", "WT-B", "KO/sample 1", "KO sample 2"),
    expected_features = c("sp|P001|PROT1", "sp|P002|PROT2"),
    expected_quantity_col = "Intensity"
  )

  regular <- module_ci_prot_make_fragpipe(variant = "regular_intensity")
  regular_import <- suppressMessages(importFragPipeData(regular, use_maxlfq = FALSE))
  expect_setequal(module_ci_prot_run_values(regular_import), c("WT A", "WT-B", "KO/sample 1", "KO sample 2"))
  expect_identical(regular_import$data_type, "protein")
  expect_identical(regular_import$column_mapping$protein_col, "Protein.Ids")
  expect_identical(regular_import$column_mapping$run_col, "Run")
  expect_identical(regular_import$column_mapping$quantity_col, "Intensity")

  no_intensity <- module_ci_prot_make_fragpipe(variant = "missing_intensity")
  expect_error(
    suppressMessages(importFragPipeData(no_intensity)),
    "No intensity columns found",
    fixed = TRUE
  )

  no_protein <- module_ci_prot_make_fragpipe(variant = "missing_protein")
  expect_error(
    suppressMessages(importFragPipeData(no_protein)),
    "Required 'Protein ID' column not found",
    fixed = TRUE
  )

  state <- module_ci_prot_apply_import_state(imported, "fragpipe")
  module_ci_prot_assert_workflow_state(state$workflow_data, "fragpipe")
})

test_that("MCI-004.5 proteomics import UI exposes test-mode uploads and production shinyFiles controls", {
  old_option <- getOption("multischolar.test_mode", FALSE)
  on.exit(options(multischolar.test_mode = old_option), add = TRUE)

  options(multischolar.test_mode = TRUE)
  test_mode_html <- htmltools::renderTags(mod_prot_import_ui("proteomics_workflow-setup_import"))$html
  expect_match(test_mode_html, "search_results_standard", fixed = TRUE)
  expect_match(test_mode_html, "fasta_file_standard", fixed = TRUE)
  expect_no_match(test_mode_html, "search_results_path", fixed = TRUE)

  options(multischolar.test_mode = FALSE)
  production_html <- htmltools::renderTags(mod_prot_import_ui("proteomics_workflow-setup_import"))$html
  if (requireNamespace("shinyFiles", quietly = TRUE)) {
    expect_match(production_html, "search_results", fixed = TRUE)
    expect_match(production_html, "fasta_file", fixed = TRUE)
    expect_match(production_html, "search_results_path", fixed = TRUE)
    expect_no_match(production_html, "search_results_standard", fixed = TRUE)
  }
})

test_that("MCI-004.6 import artifact and state assertions gate downstream design unlock", {
  diann <- suppressMessages(importDIANNData(module_ci_prot_make_diann()))
  state <- module_ci_prot_apply_import_state(diann, "diann")

  expect_true("cp01:raw_imported" %in% state$checkpoints)
  expect_identical(state$workflow_data$processing_log$setup_import$search_file, "diann_fixture.tsv")
  expect_identical(state$workflow_data$processing_log$setup_import$fasta_file, basename(state$workflow_data$fasta_file_path))
  expect_identical(state$workflow_data$processing_log$setup_import$detected_format, "diann")
  expect_identical(state$workflow_data$tab_status$setup_import, "complete")
  expect_identical(state$workflow_data$tab_status$design_matrix, "pending")
  expect_identical(module_ci_prot_expected_report_template("diann"), "DIANN_report.rmd")
  expect_true(dir.exists(file.path(state$artifact_dir, "results")))
  expect_true(dir.exists(file.path(state$artifact_dir, "scripts")))
})
