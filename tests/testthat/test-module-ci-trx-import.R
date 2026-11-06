# KIT-096 MCI-026: the transcriptomics import, the lipid import module under the transcriptomics profile. One format (the
# notebook's count table), the id and name columns chosen by the user, no sanitisation, KIT-093's count policy at import with
# its refusal shown by class.

test_that("MCI-026.1 the preview offers the count table's columns, the sample columns in file order", {
  truth <- module_ci_trx_truth()
  preview <- trxImportPreview(file.path(module_ci_trx_corpus(), "counts.tsv"))
  expect_identical(preview$detectedFormat, "trx_counts")
  expect_identical(preview$updates$lipidId$selected, "GeneID")
  expect_identical(preview$updates$annotation$selected, "GeneName")
  expect_identical(preview$updates$lipidId$choices, c("GeneID", "GeneName", truth$column_order))
  expect_identical(preview$importResult$sample_columns, truth$column_order)
  expect_null(preview$updates$isPattern)
})

test_that("MCI-026.1 the import accepts the screen corpus: names and order as in the file, the owner bound to TRANSCRIPTOMICS", {
  skip_if_not(workflowOwnerCapable())
  truth <- module_ci_trx_truth()
  for (dir in c(module_ci_trx_corpus(), module_ci_trx_fixture("counts_whole_double"))) {
    fixture <- module_ci_trx_fixture_owner()
    run <- module_ci_trx_import(fixture, file.path(dir, "counts.tsv"))
    workflow <- fixture$workflow
    expect_identical(run$result$status, "success", info = dir)
    expect_length(run$errors, 0L)
    expect_identical(workflow$data_format, "trx_counts")
    expect_identical(workflow$data_type, "transcriptomics")
    expect_identical(workflow$state_manager$workflow_type, "TRANSCRIPTOMICS")
    expect_identical(workflow$state_manager$.context$store$identity$workflow_id, "TRANSCRIPTOMICS")
    # KIT-075: the import's table is held by reference.
    expect_true(isSavedWorkflowReference(workflow$data_tbl))
    tables <- getWorkflowOmicsBusTables(workflow, "data_tbl")
    expect_identical(names(tables), "counts")
    expect_identical(names(tables$counts), c("GeneID", "GeneName", truth$column_order))
    expect_identical(nrow(tables$counts), 300L)
    expect_identical(workflow$column_mapping$lipid_id_col, "GeneID")
    expect_identical(workflow$column_mapping$annotation_col, "GeneName")
    expect_identical(workflow$column_mapping$sample_columns, truth$column_order)
    expect_identical(workflow$tab_status$setup_import, "complete")
    expect_identical(workflow$tab_status$design_matrix, "pending")
    for (artifact in c("data_cln_counts.tab", "column_mapping.json", "transcriptomics_import_summary.tsv")) {
      expect_true(file.exists(file.path(fixture$paths$source_dir, artifact)), info = artifact)
    }
    expect_false(file.exists(file.path(fixture$paths$source_dir, "lipidomics_import_summary.tsv")))
    # The source copy the design's provenance records (trx_counts_1) is the file, byte for byte.
    copies <- list.files(file.path(fixture$paths$source_dir, "import_sources"), full.names = TRUE)
    expect_length(copies, 1L)
    expect_identical(unname(tools::md5sum(copies)), unname(tools::md5sum(file.path(dir, "counts.tsv"))))
  }
})

test_that("MCI-026.2 the import refuses each KIT-040 invalid fixture and shows KIT-093's condition class", {
  skip_if_not(workflowOwnerCapable())
  refused <- c(counts_na = "trx_counts_na", counts_negative = "trx_counts_negative", counts_fractional = "trx_counts_fractional")
  for (name in names(refused)) {
    fixture <- module_ci_trx_fixture_owner()
    run <- module_ci_trx_import(fixture, file.path(module_ci_trx_fixture(name), "counts.tsv"))
    expect_null(run$result, info = name)
    shown <- vapply(Filter(function(note) identical(note$type, "error"), run$notes), `[[`, "", "message")
    expect_length(shown, 1L)
    expect_match(shown, refused[[name]], fixed = TRUE, info = name)
    expect_match(run$errors, refused[[name]], fixed = TRUE, info = name)
    # Nothing was held or bound, and the tab did not complete.
    expect_null(fixture$workflow$data_tbl)
    expect_identical(fixture$workflow$tab_status$setup_import, "pending")
    expect_identical(fixture$workflow$state_manager$.context$store$identity$workflow_id, "PENDING")
    # The check itself refuses by the class, with no rounding.
    table <- utils::read.delim(file.path(module_ci_trx_fixture(name), "counts.tsv"), check.names = FALSE, stringsAsFactors = FALSE)
    expect_error(checkTrxImportCounts(list(counts = table), "GeneID", "GeneName"), class = refused[[name]])
  }
})

test_that("MCI-026.2 a count table without its chosen name column, or with a non-numeric sample column, is refused by class", {
  table <- utils::read.delim(file.path(module_ci_trx_corpus(), "counts.tsv"), check.names = FALSE, stringsAsFactors = FALSE)
  expect_error(checkTrxImportCounts(list(counts = table), "GeneID", ""), class = "trx_counts_invalid")
  table$note <- "x"
  expect_error(checkTrxImportCounts(list(counts = table), "GeneID", "GeneName"), class = "trx_counts_invalid")
})

test_that("MCI-026.3 the import tab under the transcriptomics profile: one format, gene columns, no sanitisation, the policy stated", {
  html <- htmltools::renderTags(mod_lipid_import_ui("import", profile = trxOmicProfile()))$html
  expect_match(html, "Transcriptomics Data Import", fixed = TRUE)
  expect_match(html, "Gene ID Column", fixed = TRUE)
  expect_match(html, "Gene Name Column", fixed = TRUE)
  expect_match(html, "whole number", fixed = TRUE)
  expect_match(html, "value=\"trx_counts\"", fixed = TRUE)
  for (vendor in c("msdial", "lipidsearch", "progenesis", "custom")) expect_false(grepl(sprintf("value=\"%s\"", vendor), html, fixed = TRUE))
  for (id in c("import-sanitize_names", "import-is_pattern", "import-exclude_norm", "import-assay2_name")) {
    expect_false(grepl(id, html, fixed = TRUE), info = id)
  }
  for (id in c("trx-import-format", "trx-import-assay1-file", "trx-import-id-col", "trx-import-annotation-col", "trx-import-process")) {
    expect_match(html, sprintf("data-testid=\"%s\"", id), fixed = TRUE)
  }
  expect_false(grepl("data-testid=\"lipid-import", html, fixed = TRUE))
  # The lipid tab is unchanged.
  lipid <- htmltools::renderTags(mod_lipid_import_ui("import"))$html
  expect_match(lipid, "Lipidomics Data Import", fixed = TRUE)
  expect_match(lipid, "import-sanitize_names", fixed = TRUE)
  expect_match(lipid, "data-testid=\"lipid-import-process\"", fixed = TRUE)
})
