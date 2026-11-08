# KIT-097 MCI-031: the transcriptomics DA tab, the lipid DA module under the transcriptomics profile, presented as the
# edgeR QL analysis KIT-095 records: robust glmQLFit, glmQLFTest, BH FDR (labelled "FDR (BH)", never a q-value), the
# realised formula with its W terms, no TREAT. One pairwise result per contrast line on one fit (KIT-093 revision 1.2).

test_that("MCI-031.1 the DA tab runs on the RUVs head; the record names the edgeR QL model; counts and the FDR column are the reference's", {
  skipWithoutTrxDeclReference()
  reference <- trxDeclReference("trx_edge_12")
  fixture <- module_ci_trx_project(reference$dir)
  out <- module_ci_trx_tabs_journey(fixture)
  expect_identical(out$errors, character())
  record <- buildWorkflowAnalysisRecord(fixture$owner)
  da <- record$analyses$differential_abundance
  expect_identical(da$source_state, "trx_ruvs")
  expect_identical(da$model$model, "edgeR_QL")
  expect_identical(da$model$fit, "glmQLFit(robust = TRUE)")
  expect_identical(da$model$test, "QL F-test (glmQLFTest)")
  expect_identical(da$model$adjust, "BH")
  expect_false(isTRUE(da$model$treat_applied))
  expect_identical(da$model$treat_applied, FALSE)
  expect_match(da$model$formula_used, "^~ 0 \\+ group \\+ W_1")
  expect_gt(length(da$model$w_columns), 0L)
  expect_identical(fixture$workflow$tab_status$differential_analysis, "complete")
  # The first contrast's topTags table is the reference's (sha256 of the saved product), and the tab's rows are that table:
  # its FDR column is the reference's FDR.
  owner <- fixture$owner
  manifest <- withWorkflowStateScope(owner, function() resolveWorkflowStateValue(owner$.context, owner$.context$snapshot$da_results))
  product <- trxDeclValue(owner, manifest$results[[1L]]$products$da_contrast)$assays[[1L]]$results$top_table
  expectTrxDeclProducts(list(top_genes = product), reference, "the DA tab's first contrast")
  first <- module_ci_trx_first_contrast(out$da)
  top <- module_ci_trx_top_table(first)
  expect_identical(top, product)
  up <- sum(top$FDR < 0.05 & top$logFC > 0)
  down <- sum(top$FDR < 0.05 & top$logFC < 0)
  expect_identical(c(sum(first$significant == "Up"), sum(first$significant == "Down")), c(up, down))
  expect_match(out$summary_stats, sprintf("Up-regulated: %d", up), fixed = TRUE)
  expect_match(out$summary_stats, sprintf("Down-regulated: %d", down), fixed = TRUE)
  expect_match(out$summary_stats, "FDR (BH)", fixed = TRUE)
  expect_false(grepl("Q <|q-value", out$summary_stats))
  # The displayed table: the FDR column, labelled FDR (BH), identical to the reference's FDR.
  # KIT-106 W3: the shared panel library's table (mod_omic_da_panels.R) with trx's configured columns.
  trx <- omicGuiConfig("transcriptomics")$da
  widget <- omicDaTableWidget(out$da$da_genes_long, unique(first$comparison), "All", "all", 0.05, 1, 100000,
    columns = trx$columns)
  shown <- widget$x$data
  # DT keeps the data's column names and shows the label in the table's header.
  expect_true(grepl(">FDR (BH)<", paste(as.character(widget$x$container), collapse = ""), fixed = TRUE))
  expect_identical(shown[["FDR"]][match(rownames(top), shown[["GeneID"]])], product$FDR)
  expect_false(any(grepl("qvalue|q_value|q-value", names(shown))))
})

test_that("MCI-031.2 several pairwise contrasts run on one fit and are offered by the tab; the managed tab equals the in-process one", {
  skipWithoutTrxDeclReference()
  entries <- c("Sera_vs_RPMI=groupSera-groupRPMI", "RPMI_vs_Sera=groupRPMI-groupSera")
  in_process <- module_ci_trx_project(module_ci_trx_corpus(), contrasts = entries)
  expected <- module_ci_trx_tabs_journey(in_process)
  record <- buildWorkflowAnalysisRecord(in_process$owner)
  expect_setequal(vapply(record$analyses$differential_abundance$results, `[[`, "", "contrast"), entries)
  expect_setequal(unique(expected$da$da_genes_long$comparison), c("Sera_vs_RPMI", "RPMI_vs_Sera"))
  fixture <- module_ci_trx_project(module_ci_trx_corpus(), contrasts = entries)
  probe <- new.env()
  omicsManagedChainMocks(fixture, probe)
  out <- module_ci_trx_tabs_journey(fixture, managed = TRUE)
  expect_identical(out$errors, character())
  expect_true(all(c("trx_da_contrast", "trx_da_tab_results", "trx_da_publication") %in% probe$operations), info = paste(probe$operations, collapse = " "))
  expect_true(all(probe$statuses == "delivered"))
  expect_identical(out$da$da_genes_long, expected$da$da_genes_long)
  # The worker's publication reported success (no "Could not save" warning) and wrote the notebook's files.
  expect_false(any(grepl("Could not save", vapply(out$notes, `[[`, "", "message"), fixed = TRUE)))
  expect_true(all(file.exists(file.path(fixture$paths$da_output_dir, c("de_results_Sera_vs_RPMI.csv", "de_results_RPMI_vs_Sera.csv", "top_de_genes.tsv")))))
})

test_that("MCI-031.3 the volcanoes: the static one is the notebook's on FDR (BH); the interactive status is FDR (BH) < 0.05 with no lfc test", {
  skipWithoutTrxDeclReference()
  fixture <- module_ci_trx_project(module_ci_trx_corpus())
  out <- module_ci_trx_tabs_journey(fixture)
  # KIT-106 W3: the static volcano from trx's configuration; the interactive one is the binding's widget, given the
  # contrast's entry as the page action gives it (D3).
  trx <- omicGuiConfig("transcriptomics")$da
  contrast <- unique(out$da$da_genes_long$comparison)[[1L]]
  plot <- trx$views$static_volcano$plot(out$da, contrast, "counts", 0.05, 1)
  expect_s3_class(plot, "ggplot")
  expect_match(rlang::as_label(plot$mapping$y), "FDR", fixed = TRUE)
  results <- withWorkflowStateScope(fixture$owner, function() {
    manifest <- resolveWorkflowStateValue(fixture$owner$.context, fixture$owner$.context$snapshot$da_results)
    r <- readWorkflowOmicsDaResults(transcriptomicsDeclarationBinding, fixture$owner, manifest)
    r$theObject <- readTrxWorkflowState(fixture$owner, manifest$source)
    r
  })
  widget <- transcriptomicsDeclarationBinding$volcano_widget(results, trx$contrast_entry(fixture$workflow, contrast), "counts", 0.05)
  expect_s3_class(widget, "htmlwidget")
  table <- widget$x$data$table
  first <- module_ci_trx_first_contrast(out$da)
  status <- table$status[match(first$GeneID, rownames(table))]
  expect_false(anyNA(status))
  expect_identical(status != "nonDE", first$FDR < 0.05)
})

test_that("MCI-031.4 the DA tab under the transcriptomics profile names the model, labels FDR (BH), keeps the fold change display-only", {
  html <- htmltools::renderTags(mod_lipid_da_ui("de", profile = trxOmicProfile()))$html
  for (text in c("edgeR", "glmQLFit(robust = TRUE)", "glmQLFTest", "FDR (BH) threshold", "display only", "trx-da-run", "trx-da-download",
    "trx-da-results-tabs", "de-formula_display")) expect_match(html, text, fixed = TRUE)
  expect_false(grepl("Q-value|q-value|TREAT|\\blimma\\b|Heatmap|Load Filtered Session", html))
  lipid <- htmltools::renderTags(mod_lipid_da_ui("de"))$html
  for (text in c("Q-value Threshold:", "Heatmap", "lipid-da-run", "Load Filtered Session")) expect_match(lipid, text, fixed = TRUE)
})
