# KIT-097 MCI-032: the transcriptomics Summary & Export tab, the lipid summary module under the transcriptomics profile:
# the RUVIII export after DA with its M (group, or the technical replicates when the design has them), the integration
# export (KIT-095's RUVIII matrix), the analysis record, the results summary (the output map) and the transcriptomics
# report.

module_ci_trx_ruviii_value <- function(owner) {
  export <- owner$.context$snapshot$exports$trx_ruviii_export
  if (is.null(export)) return(NULL)
  hydrateSavedWorkflowEntity(owner$.context$store, export$reference)$value
}

#' The summary tab's save, copy and report, in order.
module_ci_trx_summary_steps <- function(session, output, out) {
  session$setInputs(`summary-experiment_label` = "kit097", `summary-description` = "KIT-097 synthetic transcriptomics")
  session$setInputs(`summary-save_workflow_args` = 1L)
  session$flushReact()
  session$setInputs(`summary-copy_to_publication` = 1L)
  session$flushReact()
  session$setInputs(`summary-generate_report` = 1L)
  session$flushReact()
  out$session_summary <- output$`summary-session_summary`
}

#' The required and forbidden statements of the rendered report (F05).
module_ci_trx_report_required <- c("edgeR quasi-likelihood F-test", "robust", "~ 0 + group + W_1", "Benjamini-Hochberg",
  "false discovery rate", "RUVs factors W_1", "model covariates", "display only", "inferred empirical negative control",
  "least significant", "selection mode", "separation score", "composite score", "RUVIII", "integration export",
  "replicate structure M")
module_ci_trx_report_forbidden <- c("limma", "eBayes", "TREAT", "q-value", "qvalue", "RUVIII-C", "RUV-III-C", "stable", "housekeeping",
  "GSE152840", "Sepsis")

#' The text of a rendered Word report.
module_ci_trx_docx_text <- function(path) {
  xml <- utils::unzip(path, "word/document.xml", exdir = tempfile())
  gsub("<[^>]+>", " ", paste(readLines(xml, warn = FALSE), collapse = " "))
}

test_that("MCI-032.1 after DA the RUVIII export records M = group, leaves the DA unchanged and its newY is the reference's; the summary writes the record, the integration export, the results summary and the report", {
  skipWithoutTrxDeclReference()
  reference <- trxDeclReference("trx_edge_12")
  fixture <- module_ci_trx_project(reference$dir)
  assign("project_dirs", fixture$project_dirs, envir = globalenv())
  withr::defer(if (exists("project_dirs", envir = globalenv())) rm("project_dirs", envir = globalenv()))
  owner <- fixture$owner
  out <- module_ci_trx_tabs_journey(fixture, ruviii = "group", summary = module_ci_trx_summary_steps)
  expect_identical(out$errors, character())
  record <- buildWorkflowAnalysisRecord(owner)
  expect_identical(record$analyses$trx_ruviii_export$outcomes$ruviii_m_source, "group")
  expect_identical(record$analyses$differential_abundance$source_state, "trx_ruvs")
  expect_identical(owner$current_state, "trx_ruvs")
  value <- module_ci_trx_ruviii_value(owner)
  expectTrxDeclProducts(list(ruviii_newY = value$newY), reference, "the summary tab's RUVIII export")
  expect_match(out$ruviii_status, "group", fixed = TRUE)
  expect_identical(owner$.context$snapshot$da_results$entity_id, out$da_reference)
  # The record beside study_parameters.txt, and the integration export: the RUVIII matrix as the notebook writes it.
  paths <- fixture$paths
  expect_true(file.exists(file.path(paths$source_dir, "analysis_record.json")))
  written <- jsonlite::read_json(file.path(paths$source_dir, "analysis_record.json"), simplifyVector = FALSE)
  expect_identical(written$analyses$trx_ruviii_export$outcomes$ruviii_m_source, "group")
  tsv <- file.path(paths$integration_dir, "transcriptomics_normalized_data.tsv")
  expect_true(file.exists(tsv))
  # The export's matrix is the reference newY (above); the file is the notebook's (vroom_write of that matrix, as rmd:714-724),
  # byte-identical to KIT-040's stable output.
  stable <- reference$store$stable_outputs
  expect_identical(digest::digest(file = tsv, algo = "sha256"),
    stable$sha256[stable$path == "results/transcriptomics/transcriptomics_normalized_data.tsv"])
  exported <- as.data.frame(vroom::vroom(tsv, show_col_types = FALSE))
  expect_identical(exported$Run, rownames(value$newY))
  expect_identical(names(exported)[-1L], colnames(value$newY))
  # The results summary: the output map's files.
  summary_dir <- paths$results_summary_dir
  expected <- c("Publication_tables/canonical_correlation_optimization_summary.xlsx", "Publication_figures/canonical_correlation_panels.png",
    "Publication_figures/canonical_correlation_panels.pdf", "Publication_tables/ruv_normalised_log_values.tsv",
    "QC_figures/normalization_comparison_grid.png", "QC_figures/normalization_comparison_grid.pdf",
    "Publication_tables/de_results_Sera_vs_RPMI.csv", "Publication_tables/DA_results_transcriptomics.xlsx",
    "Publication_tables/top_de_genes.tsv", "Publication_figures/Volcano_Plots/volcano_plot_Sera_vs_RPMI.png",
    "Publication_figures/NumSigDeMolecules/summary_barplot_Sera_vs_RPMI.png", "QC_figures/ruviii_before_and_after.png",
    "QC_figures/ruviii_before_and_after.pdf", "Study_report/study_parameters.txt")
  present <- file.exists(file.path(summary_dir, expected))
  expect_true(all(present), info = paste(expected[!present], collapse = ", "))
  workbook <- openxlsx::getSheetNames(file.path(summary_dir, "Publication_tables", "DA_results_transcriptomics.xlsx"))
  expect_true(length(workbook) >= 2L)
  # The report: rendered, every required statement and none of the forbidden ones.
  report <- list.files(summary_dir, pattern = "^transcriptomics_report.*\\.docx$", full.names = TRUE)
  expect_length(report, 1L)
  text <- module_ci_trx_docx_text(report[[1L]])
  for (statement in module_ci_trx_report_required) expect_true(grepl(statement, text, fixed = TRUE), info = statement)
  for (statement in module_ci_trx_report_forbidden) expect_false(grepl(statement, text, ignore.case = TRUE), info = statement)
})

test_that("MCI-032.2 a design with technical replicates offers M = technical_replicate; its export is recorded as such and the DA does not change", {
  skipWithoutTrxDeclReference()
  fixture <- module_ci_trx_project(module_ci_trx_corpus(), technical = TRUE)
  owner <- fixture$owner
  out <- module_ci_trx_tabs_journey(fixture, ruviii = c("group", "technical_replicate"))
  expect_identical(out$errors, character())
  record <- buildWorkflowAnalysisRecord(owner)
  expect_identical(record$analyses$trx_ruviii_export$outcomes$ruviii_m_source, "technical_replicate")
  value <- module_ci_trx_ruviii_value(owner)
  expect_identical(value$m_source, "technical_replicate")
  # The DA results are those the tab accepted: the manifest on the head is unchanged by either export.
  expect_false(is.null(out$da_reference))
  expect_identical(owner$.context$snapshot$da_results$entity_id, out$da_reference)
  expect_identical(record$analyses$differential_abundance$source_state, "trx_ruvs")
})

test_that("MCI-032.3 a design without technical replicates disables the option, and a request for it shows the classed refusal", {
  skipWithoutTrxDeclReference()
  fixture <- module_ci_trx_project(module_ci_trx_corpus(), technical = FALSE)
  out <- module_ci_trx_tabs_journey(fixture, ruviii = "technical_replicate")
  expect_null(fixture$owner$.context$snapshot$exports$trx_ruviii_export)
  shown <- vapply(out$notes, `[[`, "", "message")
  expect_true(any(grepl("trx_ruviii_replicates_invalid", shown, fixed = TRUE)), info = paste(shown, collapse = " | "))
  expect_identical(trxOmicProfile()$ruviii_m_choices(FALSE), c(`Group (the notebook)` = "group"))
  expect_true("technical_replicate" %in% trxOmicProfile()$ruviii_m_choices(TRUE))
})

test_that("MCI-032.4 under a managed session the export, its figure and the integration file are the worker's and equal the in-process ones", {
  skipWithoutTrxDeclReference()
  in_process <- module_ci_trx_project(module_ci_trx_corpus())
  module_ci_trx_tabs_journey(in_process, ruviii = "group")
  fixture <- module_ci_trx_project(module_ci_trx_corpus())
  probe <- new.env()
  omicsManagedChainMocks(fixture, probe)
  out <- module_ci_trx_tabs_journey(fixture, managed = TRUE, ruviii = "group")
  expect_identical(out$errors, character())
  expect_true(all(c("trx_ruviii_export", "trx_composite_figure", "trx_qc_figures") %in% probe$operations), info = paste(probe$operations, collapse = " "))
  expect_true(all(probe$statuses == "delivered"))
  expect_identical(module_ci_trx_ruviii_value(fixture$owner)$newY, module_ci_trx_ruviii_value(in_process$owner)$newY)
  expect_true(file.exists(file.path(fixture$paths$publication_graphs_dir, "ruviii_before_and_after.png")))
})

test_that("MCI-032.5 the summary tab under the transcriptomics profile: the RUVIII export and its M, the transcriptomics template", {
  html <- htmltools::renderTags(mod_lipid_summary_ui("summary", profile = trxOmicProfile()))$html
  # The M choice is rendered by the server (its options follow the design's technical replicates).
  for (text in c("RUVIII", "summary-ruviii_m_ui", "trx-summary-ruviii-export", "trx-summary-generate-report", "trx-summary-save-args",
    "trx-summary-copy-publication", "trx-summary-export-state")) expect_match(html, text, fixed = TRUE)
  expect_false(grepl("lipid", html, ignore.case = TRUE))
  lipid <- htmltools::renderTags(mod_lipid_summary_ui("summary"))$html
  expect_match(lipid, "lipid-summary-generate-report", fixed = TRUE)
  expect_false(grepl("ruviii", lipid, ignore.case = TRUE))
})
