# KIT-047: module contracts of the DIA_LIMPA peptide stage's three tabs (log2 matrix, peptide accession cleanup,
# DPC quantification): the UI carries the inputs the servers read; the in-memory apply steps (the path without an
# owner) save the route's states with their configuration; the DPC tab offers the installed engines only.

limpa_tab_inputs <- function(ui) {
  html <- as.character(ui)
  regmatches(html, gregexpr('id="[^"]+"', html))[[1L]]
}

test_that("the three DIA_LIMPA tabs are the route's QC specs and carry the inputs their servers read", {
  peptide <- getProtQcPeptideModuleSpecs("DIA_LIMPA")
  expect_identical(vapply(peptide, `[[`, character(1), "id"), c("qvalue_filter", "rollup", "log2_matrix", "accession_cleanup"))
  expect_true(all(vapply(peptide, function(spec) is.null(spec$placeholder), logical(1))))
  protein <- getProtQcProteinModuleSpecs("DIA_LIMPA")
  expect_identical(protein[[1L]]$id, "dpc_quant")
  expect_null(protein[[1L]]$placeholder)
  expect_true(isTRUE(protein[[1L]]$needsExperimentPaths))
  for (spec in c(peptide[3:4], protein[1L])) {
    expect_true(exists(spec$uiFn, mode = "function"), info = spec$id)
    expect_true(exists(spec$serverFn, mode = "function"), info = spec$id)
  }
  log2 <- limpa_tab_inputs(mod_prot_qc_peptide_log2_ui("x"))
  expect_true(all(c('id="x-apply_log2_matrix"', 'id="x-revert_log2_matrix"', 'id="x-log2_results"', 'id="x-log2_plot"') %in% log2))
  accession <- limpa_tab_inputs(mod_prot_qc_peptide_accession_ui("x"))
  expect_true(all(c('id="x-delimiter"', 'id="x-aggregation_method"', 'id="x-apply_peptide_accession"',
    'id="x-peptide_accession_results"') %in% accession))
  dpc <- limpa_tab_inputs(mod_prot_qc_protein_dpc_quant_ui("x"))
  expect_true(all(c('id="x-dpc_slope"', 'id="x-chunk"', 'id="x-apply_dpc_quant"', 'id="x-dpc_quant_results"') %in% dpc))
  if (requireNamespace("limpa", quietly = TRUE)) expect_true('id="x-engine"' %in% dpc)
  # The engine choices follow the installed namespaces, the shipped default first.
  expect_identical(unname(limpaDpcEngineChoices(available = function(p) p %in% c("fastlimpadpc", "limpa"))), c("limpa", "fastlimpadpc"))
  expect_identical(unname(limpaDpcEngineChoices(available = function(p) p == "limpa")), "limpa")
  expect_length(limpaDpcEngineChoices(available = function(p) FALSE), 0L)
})

test_that("the in-memory apply steps save the route's states: log2 matrix, peptide accession (skipped without FASTA), DPC quantification", {
  now <- as.POSIXct("2026-10-21 00:00:00", tz = "UTC")
  no_log <- function(...) invisible(NULL)
  rolled <- module_ci_prot_peptide_rollup_object()
  workflow <- module_ci_prot_workflow(initial_state = rolled)
  workflow$config_list <- list(globalParameters = list(workflow_type = "DIA_LIMPA"))

  log2 <- runPeptideLog2ApplyStep(workflow, logInfoFn = no_log, nowFn = function() now)
  expect_true(log2$transformedS4@is_logged_data)
  expect_identical(dim(log2$transformedS4@peptide_matrix), c(4L, 4L))
  expect_match(log2$resultText, "State saved as: 'log2_matrix'", fixed = TRUE)
  expect_identical(workflow$state_manager$getHistory(), "log2_matrix")
  expect_identical(workflow$qc_params$peptide_qc$log2_matrix$peptides, 4L)
  expect_identical(workflow$qc_params$peptide_qc$log2_matrix$timestamp, now)

  accession <- runPeptideAccessionApplyStep(workflow, delimiter = ";", aggregationMethod = "mean",
    logInfoFn = no_log, nowFn = function() now)
  expect_false(accession$cleanupApplied)
  expect_match(accession$resultText, "Peptide Accession Cleanup Skipped", fixed = TRUE)
  expect_identical(workflow$state_manager$getHistory(), c("log2_matrix", "peptide_accession_cleaned"))
  expect_identical(workflow$qc_params$peptide_qc$peptide_accession_cleanup$aggregation_method, "mean")
  expect_identical(workflow$accession_cleanup_results$cleanup_applied, FALSE)

  fake_engine <- function(theObject, dpc_slope, verbose, chunk, engine) {
    ids <- unique(theObject@peptide_data[[theObject@protein_id_column]])
    samples <- theObject@design_matrix[[theObject@sample_id]]
    E <- matrix(seq_len(length(ids) * length(samples)), nrow = length(ids), dimnames = list(ids, samples))
    table <- data.frame(Protein.Ids = ids, as.data.frame(E, check.names = FALSE), check.names = FALSE, stringsAsFactors = FALSE)
    object <- ProteinQuantitativeData(protein_quant_table = table, protein_id_column = theObject@protein_id_column,
      protein_id_table = data.frame(Protein.Ids = ids, stringsAsFactors = FALSE), design_matrix = theObject@design_matrix,
      sample_id = theObject@sample_id, group_id = theObject@group_id,
      technical_replicate_id = theObject@technical_replicate_id, args = theObject@args)
    object@args$limpa_dpc_quant_results <- list(engine = engine, engine_version = "9.9.9", dpc_parameters_used = c(-5, 0.7),
      slope_interpretation = "strong intensity-dependent missing", total_proteins_quantified = length(ids),
      total_peptides_used = nrow(theObject@peptide_matrix), missing_percentage_peptides = 12.5,
      quantified_elist = list(E = E), standard_errors = E * 0,
      peptide_counts_per_protein = data.frame(Protein.Ids = ids, peptide_count = 1L, peptidoform_count = 1L))
    object
  }
  dpc <- runProteinDpcQuantStep(workflow, engine = "limpa", dpcSlope = 0.8, chunk = 1000,
    proteinMissingValueImputationLimpaFn = fake_engine, logInfoFn = no_log, nowFn = function() now)
  expect_s4_class(dpc$quantifiedS4, "ProteinQuantitativeData")
  expect_match(dpc$resultText, "Engine: limpa 9.9.9", fixed = TRUE)
  expect_match(dpc$resultText, "State saved as: 'protein_dpc_quant'", fixed = TRUE)
  expect_identical(workflow$state_manager$getHistory(), c("log2_matrix", "peptide_accession_cleaned", "protein_dpc_quant"))
  saved <- workflow$state_manager$saved()
  expect_identical(saved[[3L]]$config_object$engine, "limpa")
  expect_identical(saved[[3L]]$config_object$dpc_intercept, -5)
  expect_identical(saved[[3L]]$s4_data_object@args$limpaDpcQuant$engine_version, "9.9.9")
  expect_identical(workflow$qc_params$protein_qc$dpc_quant$slope_interpretation, "strong intensity-dependent missing")
})
