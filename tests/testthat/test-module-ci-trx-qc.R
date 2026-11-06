# KIT-096 MCI-028: the transcriptomics QC tab: the expression-filter sub-tab (KIT-095's trx_expression_filter, its counts
# metrics, the filtering figure, the raw-logCPM PCA and RLE) and the lipid QC finalise under the transcriptomics profile
# (trx_qc_finalize). In process through the wrappers, managed through the worker.

#' The QC tab's sub-modules, as the QC orchestrator starts them, with the runner when managed. The filter sub-tab's
#' displayed values are kept on the session for the test.
trxQcTabServer <- function(fixture, managed) {
  force(fixture)
  force(managed)
  function(input, output, session) {
    if (managed) {
      runner <- installWorkflowOmicsActions(session, fixture$workflow, "transcriptomics", "kit096", fixture$paths)
      runner$progress <- function(message) list(set = function(...) invisible(NULL), inc = function(...) invisible(NULL), close = function() invisible(NULL))
    }
    session$userData$shown <- mod_trx_qc_filter_server("filter", fixture$workflow, "transcriptomics", "kit096")
    mod_lipid_qc_s4_server("finalize", fixture$workflow, "transcriptomics", "kit096", profile = trxOmicProfile())
  }
}

#' Apply the filter, then finalise; returns what the filter sub-tab shows and the finalise's text.
trxQcTabJourney <- function(fixture, managed) {
  out <- new.env(parent = emptyenv())
  suppressWarnings(suppressMessages(utils::capture.output(shiny::testServer(trxQcTabServer(fixture, managed), {
    session$setInputs(`filter-apply_filter` = 1L)
    session$flushReact()
    out$shown <- shiny::isolate(session$userData$shown())
    out$filter_text <- output$`filter-filter_results`
    session$setInputs(`finalize-finalize_qc` = 1L)
    session$flushReact()
    out$finalize_text <- output$`finalize-finalize_results`
  }))))
  as.list(out)
}

test_that("MCI-028.1 in process the filter commits trx_expression_filtered and the finalise trx_qc_complete; the record is valid", {
  skip_if_not(workflowOwnerCapable())
  skip_if_not_installed("edgeR")
  fixture <- trxDeclarationsOwner(withr::local_tempdir(), module_ci_trx_corpus())
  withr::defer(closeDIAStateAdapter(fixture$owner))
  owner <- fixture$owner
  before <- trxAssayInputs(fixture$object)
  kernel <- trxExpressionFilter(before$counts, before$group)
  out <- trxQcTabJourney(fixture, managed = FALSE)
  expect_identical(owner$getHistory(), c("initial", "trx_raw_data_s4", "trx_expression_filtered", "trx_qc_complete"))
  expect_identical(fixture$workflow$tab_status$quality_control, "complete")
  # The counts at each step, the zero-variance step included, are the kernel's.
  expect_identical(out$shown$counts, list(before = nrow(before$counts), after_filter_by_expr = sum(kernel$keep_expr),
    after_zero_variance = sum(kernel$keep_var)))
  expect_match(out$filter_text, "Zero-variance", fixed = TRUE)
  filtered <- trxAssayInputs(withWorkflowStateScope(owner, function() readTrxWorkflowState(owner, "trx_expression_filtered")))
  expect_identical(out$shown$metrics, trxFilteringMetrics(withWorkflowStateScope(owner, function() readTrxWorkflowState(owner, "trx_expression_filtered"))))
  expect_identical(names(out$shown$metrics$library_size), colnames(filtered$counts))
  raw <- trxTmmLogCpm(filtered$counts, filtered$design$group)$log_cpm_raw
  expect_identical(out$shown$stats$pca, trxPcaStats(raw, trxFigureDesign(filtered)))
  expect_identical(out$shown$stats$rle, trxRleStats(raw, trxFigureDesign(filtered)))
  expect_s3_class(trxQcFilteringFigure(out$shown$counts, out$shown$metrics), "ggplot")
  expect_match(out$finalize_text, "trx_qc_complete", fixed = TRUE)
  expect_true(validateWorkflowAnalysisRecord(buildWorkflowAnalysisRecord(owner))$valid)
})

test_that("MCI-028.2 the QC tab's filter mask on a GUI-built design is identical to the KIT-040 reference mask", {
  skipWithoutTrxDeclReference()
  reference <- trxDeclReference("trx_edge_12")
  fixture <- module_ci_trx_fixture_owner()
  module_ci_trx_designed(fixture, reference$dir)
  owner <- fixture$workflow$state_manager
  out <- trxQcTabJourney(list(workflow = fixture$workflow, paths = fixture$paths), managed = FALSE)
  expect_identical(owner$current_state, "trx_qc_complete")
  filtered <- trxAssayInputs(withWorkflowStateScope(owner, function() readTrxWorkflowState(owner, "trx_expression_filtered")))$counts
  # The GUI's design keeps the counts file's run order; the reference is in the design file's order, so the comparison
  # puts the runs in that order by name. The genes kept, their order and every count are the reference's.
  design_order <- utils::read.delim(file.path(reference$dir, "design_matrix.tab"), stringsAsFactors = FALSE)$Run
  expectTrxDeclProducts(list(count_matrix_filtered = filtered[, design_order]), reference, "the QC tab's filter")
  prints <- reference$store$fingerprints
  expect_identical(out$shown$counts$after_filter_by_expr, as.integer(prints$dims[prints$product == "gene_vars"]))
  expect_identical(out$shown$counts$after_zero_variance, nrow(filtered))
})

test_that("MCI-028.3 under a managed session the metrics and figures are the worker's, and equal the in-process ones", {
  skip_if_not(workflowOwnerCapable())
  skip_if_not_installed("edgeR")
  in_process <- trxDeclarationsOwner(file.path(withr::local_tempdir(), "in-process"), module_ci_trx_corpus())
  withr::defer(closeDIAStateAdapter(in_process$owner))
  expected <- trxQcTabJourney(in_process, managed = FALSE)
  fixture <- trxDeclarationsOwner(file.path(withr::local_tempdir(), "managed"), module_ci_trx_corpus())
  withr::defer(closeDIAStateAdapter(fixture$owner))
  probe <- new.env()
  omicsManagedChainMocks(fixture, probe)
  out <- trxQcTabJourney(fixture, managed = TRUE)
  expect_true(all(c("trx_expression_filter", "trx_qc_figures", "trx_filtering_metrics") %in% probe$operations), info = paste(probe$operations, collapse = " "))
  expect_true(all(probe$statuses == "delivered"))
  expect_identical(fixture$owner$getHistory(), in_process$owner$getHistory())
  expect_identical(out$shown, expected$shown)
  expect_identical(out$filter_text, expected$filter_text)
  # The worker staged the raw figures for the QC directory.
  expect_true(all(file.exists(file.path(fixture$paths$gene_qc_dir, c("trx_raw_pca.png", "trx_raw_rle.png")))))
})

test_that("MCI-028.4 the QC tab under the transcriptomics profile: the filter and finalise sub-tabs, the orchestrator's class", {
  html <- htmltools::renderTags(mod_lipid_qc_ui("qc", profile = trxOmicProfile()))$html
  expect_match(html, "Transcriptomics Quality Control", fixed = TRUE)
  tabs <- htmltools::renderTags(trxOmicProfile()$qc_tabs(shiny::NS("qc")))$html
  expect_match(tabs, "Expression Filter", fixed = TRUE)
  expect_match(tabs, "Finalize QC", fixed = TRUE)
  expect_match(tabs, "qc-filter-apply_filter", fixed = TRUE)
  expect_match(tabs, "qc-s4_finalize-finalize_qc", fixed = TRUE)
  expect_false(grepl("Intensity Filter|Internal Standard|Duplicate", tabs))
  # The finalise shows no per-assay table for the one count table.
  expect_false(grepl("assay_stats_table", tabs, fixed = TRUE))
})

test_that("KIT-101 F-11: the expression filter before any imported state refuses with a message naming what is missing", {
  skip_if_not(workflowOwnerCapable())
  notes <- character()
  testthat::local_mocked_bindings(showNotification = function(ui, ..., type = "default") {
    notes <<- c(notes, paste0("[", type, "] ", ui))
    invisible(NULL)
  }, .package = "shiny")
  paths <- list(source_dir = withr::local_tempdir())
  suppressMessages(shiny::testServer(function(input, output, session) {
    # The session's owner as the transcriptomics module creates it at start, before the import.
    workflow_data <- shiny::reactiveValues(state_manager = NULL, config_list = list(), tab_status = list())
    workflow_data$state_manager <- createWorkflowSessionOwner(workflow_data, paths, session, omic_type = "transcriptomics")
    session$userData$owner <- workflow_data$state_manager
    mod_trx_qc_filter_server("filter", workflow_data, "transcriptomics", "kit101")
  }, {
    session$setInputs(`filter-apply_filter` = 1L)
    expect_identical(diaStateHistory(session$userData$owner), "initial")
  }))
  expect_identical(notes, "[warning] Please import the count data before applying the expression filter.")
})
