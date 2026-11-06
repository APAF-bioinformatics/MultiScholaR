# KIT-097 MCI-030: the transcriptomics normalisation tab, the lipid norm module under the transcriptomics profile: TMM and
# logCPM, the RUV grid (automatic, the notebook's; or manual, k within the cap) with its CCA panels and grid table, RUVs (the
# DA head), the PCA and RLE at raw, TMM and RUVs. In process through the wrappers, managed through the worker.

module_ci_trx_ruvs_outcomes <- function(owner) {
  record <- buildWorkflowAnalysisRecord(owner)
  steps <- Filter(function(step) identical(step$operation, "trx_ruvs"), record$steps)
  if (length(steps)) steps[[length(steps)]]$outcomes else list()
}

test_that("MCI-030.1 the automatic run commits trx_ruvs and records the mode, grid, k, percentage, controls and scores; k, percentage and controls are the reference's", {
  skipWithoutTrxDeclReference()
  reference <- trxDeclReference("trx_edge_12")
  fixture <- module_ci_trx_project(reference$dir)
  owner <- fixture$owner
  out <- module_ci_trx_tabs_journey(fixture, da = FALSE)
  expect_identical(out$errors, character())
  expect_identical(owner$current_state, "trx_ruvs")
  expect_identical(fixture$workflow$tab_status$normalization, "complete")
  outcomes <- module_ci_trx_ruvs_outcomes(owner)
  expect_identical(outcomes$ruv_mode, "automatic")
  expect_identical(unlist(outcomes[c("percentage_min", "percentage_max", "percentage_step", "k_penalty_weight")]),
    c(percentage_min = 5, percentage_max = 30, percentage_step = 5, k_penalty_weight = 0.05))
  for (name in c("k", "negative_control_percentage", "negative_control_count", "separation_score", "composite_score", "max_acceptable_k")) {
    expect_false(is.null(outcomes[[name]]), info = name)
  }
  expect_identical(outcomes$controls, "inferred")
  # The grid the RUVs state used: its winner's k, percentage and control set are the KIT-040 reference's.
  grid <- withWorkflowStateScope(owner, function() trxProductValue(owner, trxStateRuvs(owner, getWorkflowStateReference(owner))$grid))
  expectTrxDeclProducts(list(best_control_genes = grid$control_genes,
    best_res = grid$best[c("percentage", "n_controls", "best_k", "separation_score", "k_penalty", "composite_score", "ctl_idx", "control_genes")]),
    reference, "the norm tab's automatic grid")
  expect_identical(outcomes$k, grid$chosen$k)
  expect_identical(outcomes$negative_control_percentage, grid$chosen$percentage)
  # What the tab shows: the grid's summary with the inferred controls, the log, the figures at the three stages.
  expect_match(out$ruv_summary, sprintf("Best k: %d", grid$chosen$k), fixed = TRUE)
  expect_match(out$ruv_summary, sprintf("Controls: %d", length(grid$control_genes)), fixed = TRUE)
  expect_match(out$norm_log, "RUVs", fixed = TRUE)
  for (stage in c("raw", "tmm", "ruvs")) for (kind in c("pca", "rle")) {
    expect_true(file.exists(file.path(fixture$paths$gene_qc_dir, sprintf("trx_%s_%s.png", stage, kind))), info = paste(stage, kind))
  }
  expect_true(file.exists(file.path(fixture$paths$publication_graphs_dir, "canonical_correlation_panels.png")))
  expect_true(file.exists(file.path(fixture$paths$publication_graphs_dir, "canonical_correlation_optimization_summary.xlsx")))
  expect_true(file.exists(file.path(fixture$paths$publication_graphs_dir, "normalization_comparison_grid.png")))
  expect_true(file.exists(file.path(fixture$paths$gene_qc_dir, "ruv_normalised_log_values.tsv")))
  expect_true(validateWorkflowAnalysisRecord(buildWorkflowAnalysisRecord(owner))$valid)
})

test_that("MCI-030.2 the manual run records its mode, k and percentage; a k above the cap is refused with its class and the head kept", {
  skipWithoutTrxDeclReference()
  fixture <- module_ci_trx_project(module_ci_trx_corpus())
  out <- module_ci_trx_tabs_journey(fixture, da = FALSE, norm = list(ruv_mode = "manual", ruv_k = 2L, ruv_percentage = 10))
  expect_identical(fixture$owner$current_state, "trx_ruvs")
  outcomes <- module_ci_trx_ruvs_outcomes(fixture$owner)
  expect_identical(outcomes$ruv_mode, "manual")
  expect_identical(outcomes$k, 2L)
  expect_identical(outcomes$negative_control_percentage, 10)
  above <- module_ci_trx_project(module_ci_trx_corpus())
  refused <- module_ci_trx_tabs_journey(above, da = FALSE, norm = list(ruv_mode = "manual", ruv_k = 4L, ruv_percentage = 10))
  expect_identical(above$owner$current_state, "trx_qc_complete")
  shown <- vapply(refused$notes, `[[`, "", "message")
  expect_true(any(grepl("trx_ruv_k_above_cap", shown, fixed = TRUE)), info = paste(shown, collapse = " | "))
})

test_that("MCI-030.3 under a managed session the norm tab runs on the worker and commits what the wrappers commit", {
  skipWithoutTrxDeclReference()
  in_process <- module_ci_trx_project(module_ci_trx_corpus())
  module_ci_trx_tabs_journey(in_process, da = FALSE)
  fixture <- module_ci_trx_project(module_ci_trx_corpus())
  probe <- new.env()
  omicsManagedChainMocks(fixture, probe)
  out <- module_ci_trx_tabs_journey(fixture, managed = TRUE, da = FALSE)
  expect_identical(out$errors, character())
  expect_true(all(c("trx_tmm_logcpm", "trx_ruv_grid", "trx_ruvs", "trx_qc_figures", "trx_ruv_publication", "trx_composite_figure") %in%
    probe$operations), info = paste(probe$operations, collapse = " "))
  expect_true(all(probe$statuses == "delivered"))
  expect_identical(fixture$owner$getHistory(), in_process$owner$getHistory())
  expect_identical(module_ci_trx_ruvs_outcomes(fixture$owner), module_ci_trx_ruvs_outcomes(in_process$owner))
  expect_true(file.exists(file.path(fixture$paths$gene_qc_dir, "trx_ruvs_pca.png")))
})

test_that("MCI-030.4 the norm tab under the transcriptomics profile: automatic or manual RUVs, inferred controls, no skip, the stage figures", {
  html <- htmltools::renderTags(mod_lipid_norm_ui("norm", profile = trxOmicProfile()))$html
  expect_match(html, "inferred negative controls (the least significant genes in a preliminary group GLM)", fixed = TRUE)
  expect_match(html, "trx-norm-run", fixed = TRUE)
  expect_match(html, "trx-norm-ruv-mode", fixed = TRUE)
  expect_match(html, "norm-ruv_k", fixed = TRUE)
  expect_match(html, "norm-ruv_percentage", fixed = TRUE)
  for (label in c("RUV QC", "PCA", "RLE", "Raw logCPM", "TMM logCPM", "RUVs", "norm-ruv_qc_ui", "norm-pca_ruv_corrected_assay1")) {
    expect_match(html, label, fixed = TRUE)
  }
  expect_false(grepl("value=\"skip\"", html, fixed = TRUE))
  expect_false(grepl("stable|housekeeping|ITSD|Density|Correlation Filtering|Cyclic Loess", html))
  # The lipid tab is unchanged.
  lipid <- htmltools::renderTags(mod_lipid_norm_ui("norm"))$html
  expect_match(lipid, "ITSD Selection", fixed = TRUE)
  expect_match(lipid, "value=\"skip\"", fixed = TRUE)
})
