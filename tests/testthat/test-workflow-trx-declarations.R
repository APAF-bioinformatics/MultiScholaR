# KIT-095: the transcriptomics binding on the shared omics core. Every state and product of the in-process journey and of
# the worker is identical to the KIT-040 reference (and so to KIT-093's kernels, which are); the record is valid and replays
# exactly, the realised DA formula included; the manual mode, the refusals, the technical-replicate RUVIII, the partitioned
# variants, every worker step and publication operation, the registry shape and the owner refusals.
library(testthat)

trxDeclMocks <- function(env = parent.frame()) testthat::local_mocked_bindings(workflowExecutionBuildIdentity = function() strrep("f", 64),
  finishWorkflowExecutionJob = function(owner, status) { owner$job <- NULL; list(status = status) }, .package = "MultiScholaR", .env = env)

#' A replay project: a copy of the fixture's store cut to the design state.
trxDeclReplayProject <- function(fixture, harness, api) {
  copy <- withr::local_tempdir(.local_envir = parent.frame())
  stopifnot(system2("cp", c("-a", file.path(fixture$source_dir, "."), copy)) == 0L)
  store_root <- dirname(dirname(list.files(fixture$source_dir, "^store[.]json$", recursive = TRUE, all.files = TRUE, full.names = TRUE)[[1L]]))
  copy_root <- file.path(copy, substring(store_root, nchar(fixture$source_dir) + 2L))
  list(prepared = harness$kitPrepareReplayProject(copy_root, c("initial", "trx_raw_data_s4"), api), origin = store_root)
}

#' Open a replayed project's store for reading its products.
trxDeclOpen <- function(root) {
  identity <- readSavedStoreIdentity(file.path(root, ".workflow-entities", "store.json"))
  newWorkflowStateOwner(newWorkflowEntityStore(root, identity$project_id, identity$workflow_id), reopen = TRUE)
}

test_that("KIT-095: the journey's states and products are the KIT-040 reference's; the record is valid and replays exactly", {
  skipWithoutTrxDeclReference()
  reference <- trxDeclReference("trx_edge_12")
  fixture <- trxDeclarationsOwner(withr::local_tempdir(), reference$dir)
  withr::defer(closeDIAStateAdapter(fixture$owner))
  owner <- fixture$owner
  out <- suppressWarnings(suppressMessages(trxDeclarationsJourney(fixture)))
  history <- unlist(owner$getHistory())
  expect_identical(history, c("initial", "trx_raw_data_s4", "trx_expression_filtered", "trx_qc_complete", "trx_ruvs"))
  products <- trxDeclJourneyProducts(fixture, out)
  expectTrxDeclProducts(products, reference, "journey")
  # The filter's counts are the kernel's.
  x <- trxAssayInputs(fixture$object)
  kernel <- trxExpressionFilter(x$counts, x$design$group)
  expect_identical(out$filter$value$counts, list(before = nrow(x$counts), after_filter_by_expr = sum(kernel$keep_expr),
    after_zero_variance = sum(kernel$keep_var)))
  # RUVs: the counts unchanged, the design with W_1..W_k; the DA manifest's model is the edgeR QL one, its formula the realised one.
  expect_identical(products$ruvs_counts, products$count_matrix_filtered)
  manifest <- out$da$manifest
  expect_identical(manifest$schema, "trx.da_results.v2")
  k <- out$ruvs$value$k
  expect_identical(manifest$parameters$formula, paste("~ 0 + group", paste(paste0("W_", seq_len(k)), collapse = " + "), sep = " + "))
  expect_identical(manifest$model[c("model", "dispersion", "fit", "test", "adjust", "formula_used", "w_columns", "treat_applied")],
    list(model = "edgeR_QL", dispersion = "estimateDisp(robust = TRUE)", fit = "glmQLFit(robust = TRUE)", test = "QL F-test (glmQLFTest)",
      adjust = "BH", formula_used = manifest$parameters$formula, w_columns = paste0("W_", seq_len(k)), treat_applied = FALSE))
  expect_true(isTrxResultManifest(manifest))
  # The record: valid, the filter's counts, the RUV outcomes, the DA model and the RUVIII export with its M.
  record <- buildWorkflowAnalysisRecord(owner)
  expect_true(validateWorkflowAnalysisRecord(record)$valid)
  step <- function(name) Filter(function(s) identical(s$state, name), record$steps)[[1L]]
  expect_identical(step("trx_expression_filtered")$outcomes[c("genes_before", "genes_after_filter_by_expr", "genes_after_zero_variance")],
    list(genes_before = out$filter$value$counts$before, genes_after_filter_by_expr = out$filter$value$counts$after_filter_by_expr,
      genes_after_zero_variance = out$filter$value$counts$after_zero_variance))
  ruv <- step("trx_ruvs")$outcomes
  expect_identical(ruv$ruv_mode, "automatic")
  expect_identical(ruv$k, k)
  expect_identical(ruv$controls, "inferred")
  expect_true(ruv$ruv_applied)
  expect_identical(ruv$negative_control_count, products$best_res$n_controls)
  expect_identical(vapply(step("trx_ruvs")$preceded_by, `[[`, "", "operation"), c("trx_tmm_logcpm", "trx_ruv_grid"))
  expect_identical(record$analyses$differential_abundance$operation, "trx_da_contrast")
  expect_identical(record$analyses$differential_abundance$model$model, "edgeR_QL")
  export <- record$analyses$trx_ruviii_export
  expect_identical(export$operation, "trx_ruviii_export")
  expect_identical(export$outcomes$ruviii_m_source, "group")
  expect_identical(export$outcomes$k, k)
  # The report overlay reads the names KIT-097's template prints.
  json <- file.path(withr::local_tempdir(), "analysis_record.json")
  writeWorkflowAnalysisRecord(owner, json)
  overlay <- workflowRecordReportParameters(json, list())
  expect_identical(overlay[c("ruv_mode", "ruviii_m_source", "da_model", "s4_formula_string")],
    list(ruv_mode = "automatic", ruviii_m_source = "group", da_model = "edgeR_QL", s4_formula_string = manifest$parameters$formula))
  expect_identical(overlay$best_k, k)
  expect_true(overlay$ruv_applied)
  # The replay on a copy cut to the design state: every state hash-identical, the DA and the export re-run to the same products.
  harness <- new.env(parent = globalenv())
  sys.source(file.path(trxDeclRoot(), "tools", "profiling", "workflow_kit.R"), harness)
  api <- function(name) get(name, asNamespace("MultiScholaR"))
  inputs <- c(trx_counts_1 = file.path(fixture$source_dir, "import_sources", "counts.tsv"), design_file = file.path(fixture$source_dir, "design_matrix.tab"))
  record_read <- jsonlite::read_json(json, simplifyVector = FALSE)
  project <- trxDeclReplayProject(fixture, harness, api)
  replay <- suppressWarnings(suppressMessages(harness$kitReplayJourney(record_read, project$prepared, project$origin, inputs, fixture$paths,
    withr::local_tempdir(), api)))
  expect_identical(replay$states, history)
  expect_identical(as.character(replay$mismatches), character())
  replayed <- trxDeclOpen(project$prepared)
  withr::defer(closeWorkflowStateOwner(replayed))
  again <- resolveWorkflowStateValue(replayed$.context, replayed$.context$snapshot$da_results)
  expect_identical(again$parameters, manifest$parameters)
  expect_identical(trxDeclValue(replayed, again$results[[1L]]$products$da_contrast)$assays[[1L]]$results$top_table, products$top_genes)
  expect_identical(trxDeclValue(replayed, replayed$.context$snapshot$exports$trx_ruviii_export$reference)$newY, products$ruviii_newY)
  # A record whose DA formula is not the replayed head's is refused.
  tampered <- record_read
  for (i in seq_along(tampered$analyses$differential_abundance$parameters)) {
    if (identical(tampered$analyses$differential_abundance$parameters[[i]]$name, "formula"))
      tampered$analyses$differential_abundance$parameters[[i]]$value <- encodeWorkflowRecordValue("~ 0 + group + W_1")
  }
  project <- trxDeclReplayProject(fixture, harness, api)
  expect_error(suppressWarnings(suppressMessages(replayWorkflowAnalysisRecord(tampered, project$prepared, inputs, fixture$paths,
    withr::local_tempdir()))), class = "workflow_replay_formula_mismatch")
})

test_that("KIT-095: the manual mode records its mode and outcomes; skip and a manual k above the cap are refused", {
  skipWithoutTrxDeclReference()
  reference <- trxDeclReference("trx_edge_12")
  fixture <- trxDeclarationsOwner(withr::local_tempdir(), reference$dir)
  withr::defer(closeDIAStateAdapter(fixture$owner))
  owner <- fixture$owner
  suppressWarnings(suppressMessages(runTrxReferenceExpressionFilter(owner, list())))
  commitWorkflowOmicsStoreStep(owner, "trx_qc_finalize", "trx_qc_complete", list(), "QC processing complete - ready for normalization")
  expect_error(suppressWarnings(suppressMessages(runTrxReferenceRuvs(owner, trxRuvGridParameters(ruv_mode = "skip"), list()))),
    class = "trx_ruv_mode_unsupported")
  # trx_edge_12 has 12 samples: the cap is 3.
  expect_error(suppressWarnings(suppressMessages(runTrxReferenceRuvs(owner, trxRuvGridParameters(ruv_mode = "manual",
    manual_percentage = 10, manual_k = 4L), list()))), class = "trx_ruv_k_above_cap")
  expect_identical(owner$current_state, "trx_qc_complete")
  manual <- suppressWarnings(suppressMessages(runTrxReferenceRuvs(owner, trxRuvGridParameters(ruv_mode = "manual",
    manual_percentage = 10, manual_k = 2L), list())))
  expect_identical(manual$grid$value$chosen$k, 2L)
  expect_identical(vapply(manual$grid$value$grid_members, `[[`, numeric(1), "percentage"), 10)
  expect_identical(manual$value$k, 2L)
  expect_identical(colnames(manual$value$W), c("W_1", "W_2"))
  outcomes <- Filter(function(s) identical(s$state, "trx_ruvs"), buildWorkflowAnalysisRecord(owner)$steps)[[1L]]$outcomes
  expect_identical(outcomes[c("ruv_mode", "k", "negative_control_percentage", "controls", "ruv_applied")],
    list(ruv_mode = "manual", k = 2L, negative_control_percentage = 10, controls = "inferred", ruv_applied = TRUE))
})

test_that("KIT-095: RUVIII with technical replicates is a distinct product with its M, recorded as such", {
  skipWithoutTrxDeclReference()
  reference <- trxDeclReference("trx_edge_12")
  design <- trxDeclTechnicalDesign(reference$dir)
  fixture <- trxDeclarationsOwner(withr::local_tempdir(), reference$dir, design = design, technical_replicate_id = "tech_rep_group")
  withr::defer(closeDIAStateAdapter(fixture$owner))
  owner <- fixture$owner
  suppressWarnings(suppressMessages({
    runTrxReferenceExpressionFilter(owner, list())
    commitWorkflowOmicsStoreStep(owner, "trx_qc_finalize", "trx_qc_complete", list(), "QC processing complete - ready for normalization")
    ruvs <- runTrxReferenceRuvs(owner, trxRuvGridParameters(), list())
    group <- runTrxReferenceRuviiiExport(owner, "group")
    tech <- runTrxReferenceRuviiiExport(owner, "technical_replicate")
  }))
  expect_false(identical(group$reference$entity_id, tech$reference$entity_id))
  log_cpm <- withWorkflowStateScope(owner, function() trxSavedMatrix(owner, ruvs$logcpm$value$log_cpm))
  M <- ruv::replicate.matrix(design$tech_rep_group)
  expect_identical(tech$value$M, M)
  expect_identical(tech$value$newY, ruv::RUVIII(Y = t(log_cpm), M = M, ctl = group$value$ctl_mask, k = ruvs$value$k, return.info = TRUE)$newY)
  expect_false(identical(tech$value$newY, group$value$newY))
  expect_identical(c(group$value$m_source, tech$value$m_source), c("group", "technical_replicate"))
  export <- buildWorkflowAnalysisRecord(owner)$analyses$trx_ruviii_export
  expect_identical(export$outcomes$ruviii_m_source, "technical_replicate")
  expect_identical(decodeWorkflowRecordValue(Filter(function(p) identical(p$name, "ruviii_m_source"), export$parameters)[[1L]]$value),
    "technical_replicate")
  expect_identical(export$product_entity_id, tech$reference$entity_id)
  # Without technical replicates the export is refused by class.
  plain <- trxDeclarationsOwner(withr::local_tempdir(), reference$dir)
  withr::defer(closeDIAStateAdapter(plain$owner))
  suppressWarnings(suppressMessages({
    runTrxReferenceExpressionFilter(plain$owner, list())
    runTrxReferenceRuvs(plain$owner, trxRuvGridParameters(), list())
  }))
  expect_error(suppressWarnings(suppressMessages(runTrxReferenceRuviiiExport(plain$owner, "technical_replicate"))),
    class = "trx_ruviii_replicates_invalid")
})

test_that("KIT-095: the partitioned filter and logCPM equal whole on a corpus whose TMM factors are not 1; the template never calls edgeR", {
  skipWithoutTrxDeclReference()
  reference <- trxDeclReference("trx_edge_16")
  root <- withr::local_tempdir()
  whole <- trxDeclarationsOwner(file.path(root, "whole"), reference$dir)
  part <- trxDeclarationsOwner(file.path(root, "partitioned"), reference$dir)
  withr::defer({ closeDIAStateAdapter(whole$owner); closeDIAStateAdapter(part$owner) })
  calls <- new.env()
  real <- edgeR::DGEList
  local_mocked_bindings(DGEList = function(counts, ...) {
    calls$rows <- c(calls$rows, nrow(counts))
    real(counts, ...)
  }, .package = "edgeR")
  w <- suppressWarnings(suppressMessages(runTrxReferenceExpressionFilter(whole$owner, list())))
  calls$rows <- NULL
  p <- suppressWarnings(suppressMessages(runTrxReferenceExpressionFilter(part$owner, list(), variant = "partitioned")))
  # The pre-pass is the one DGEList call, on every gene; the template path makes none.
  expect_identical(calls$rows, nrow(whole$object@count_data[[1L]]))
  expect_identical(part$owner$.context$snapshot$step_provenance$trx_expression_filtered$variant, "partitioned")
  expect_identical(p$value[setdiff(names(p$value), "final")], w$value[setdiff(names(w$value), "final")])
  expect_false(isTRUE(all.equal(w$value$effective_lib_size, colSums(trxAssayInputs(whole$object)$counts))))
  read <- function(fixture) withWorkflowStateScope(fixture$owner, function() readTrxWorkflowState(fixture$owner))
  expect_identical(trxAssayInputs(read(part)), trxAssayInputs(read(whole)))
  tmm <- function(fixture, variant) {
    entry <- workflowOperationEntry("trx_tmm_logcpm")
    withWorkflowStateScope(fixture$owner, function() {
      run <- runSavedOperation(fixture$owner, entry$declaration, getWorkflowStateReference(fixture$owner), list(variant = variant))
      list(samples = run$value$samples, log_cpm = trxSavedMatrix(fixture$owner, run$value$log_cpm))
    })
  }
  tw <- suppressWarnings(suppressMessages(tmm(whole, "whole")))
  tp <- suppressWarnings(suppressMessages(tmm(part, "partitioned")))
  expect_true(any(tw$samples$norm.factors != 1))
  expect_identical(tp, tw)
  expect_identical(trxDeclSha(tw$log_cpm), reference$store$fingerprints$sha256[reference$store$fingerprints$product == "log_cpm_tmm"])
  expect_identical(workflowOperationEntry("trx_expression_filter")$variants, c("whole", "partitioned"))
  expect_identical(workflowOperationEntry("trx_tmm_logcpm")$variants, c("whole", "partitioned"))
})

#' A runner that dispatches each action to the worker harness; `staged` receives the session file the export staged.
trxDeclSyncRunner <- function(fixture, controller, staged) {
  list(run = function(actions, done = NULL, failed = NULL) {
    for (action in actions) {
      result <- omicsWorkerAction(fixture, controller, action$operation, action$parameters, action$context)
      staged$session <- readRDS(file.path(workflowExportPaths(action$parameters$root, action$parameters$request_id)$request_dir, "session.rds"))
      action$deliver(result)
    }
    if (is.function(done)) done()
    invisible(NULL)
  }, retrying = function() FALSE, progress = function(message) list(close = function() invisible(NULL)))
}

test_that("KIT-095: every trx worker step and publication operation completes through the worker harness, as in process", {
  skipWithoutTrxDeclReference()
  trxDeclMocks()
  reference <- trxDeclReference("trx_edge_12")
  fixture <- trxDeclarationsOwner(withr::local_tempdir(), reference$dir)
  withr::defer(closeDIAStateAdapter(fixture$owner))
  owner <- fixture$owner
  binding <- transcriptomicsDeclarationBinding
  controller <- omicsWorkerController(fixture)
  act <- function(operation, parameters = list(), context = NULL) {
    done <- suppressWarnings(suppressMessages(omicsWorkerAction(fixture, controller, operation, parameters, context)))
    expect_identical(done$acceptance, "complete", info = operation)
    done
  }
  commit <- function(operation, parameters, state_name) act(operation, parameters, list(state_name = state_name, config = list(), description = operation))
  filter <- commit("trx_expression_filter", list(), "trx_expression_filtered")
  expect_true(isSavedWorkflowReference(filter$display$filtering_metrics))
  commitWorkflowOmicsStoreStep(owner, "trx_qc_finalize", "trx_qc_complete", list(), "QC processing complete - ready for normalization")
  logcpm <- act("trx_tmm_logcpm")
  grid <- act("trx_ruv_grid", c(list(logcpm = logcpm$reference$entity_id), trxRuvGridParameters()))
  ruvs <- commit("trx_ruvs", list(logcpm = logcpm$reference$entity_id, grid = grid$reference$entity_id), "trx_ruvs")
  source <- getWorkflowStateReference(owner)
  settings <- list(formula_string = trxRealisedFormula(owner, source), da_q_val_thresh = 0.05, treat_lfc_cutoff = 1)
  values <- list()
  for (step in omicsDaContrastSteps(binding, settings, fixture$workflow$contrasts_tbl)) values[[step$contrast]] <- act(step$operation, step$parameters)$value
  da <- withWorkflowStateScope(owner, function() acceptWorkflowOmicsDaManifest(binding, owner, source, settings, fixture$workflow$contrasts_tbl, values))
  ruviii <- act("trx_ruviii_export", list(logcpm = logcpm$reference$entity_id, grid = grid$reference$entity_id, ruviii_m_source = "group"))
  acceptTrxRuviiiExport(owner, source, ruviii)
  out <- list(filter = filter, ruvs = list(value = ruvs$value, logcpm = logcpm, grid = grid), da = da, ruviii = ruviii)
  expectTrxDeclProducts(trxDeclJourneyProducts(fixture, out), reference, "worker")
  expect_identical(unlist(owner$getHistory()), c("initial", "trx_raw_data_s4", "trx_expression_filtered", "trx_qc_complete", "trx_ruvs"))
  expect_true(validateWorkflowAnalysisRecord(buildWorkflowAnalysisRecord(owner))$valid)
  # Every compute declaration on the omics role was run above; the store steps refuse dispatch.
  ran <- c("trx_expression_filter", "trx_tmm_logcpm", "trx_ruv_grid", "trx_ruvs", "trx_da_contrast", "trx_ruviii_export")
  expect_setequal(Filter(function(name) identical(workflowActionSpec(name)$omic_type, "transcriptomics") && identical(workflowActionSpec(name)$compute, "omics") &&
    !is.null(workflowActionSpec(name)$declaration), workflowOperationNames()), ran)
  for (name in c("trx_qc_finalize", "trx_reset")) expect_error(assertWorkflowOperationDispatchable(name), class = "workflow_execution_unsupported")

  # The publication operations. The notebook's CSV and TSV outputs equal KIT-040's stable hashes.
  stable <- reference$store$stable_outputs
  hash_of <- function(path) stable$sha256[stable$path == path]
  metrics <- act("trx_filtering_metrics")
  expect_identical(withWorkflowStateScope(owner, function() resolveWorkflowStateValue(owner$.context, metrics$display$filtering_metrics)),
    binding$filtering_metrics(readTrxWorkflowState(owner)))
  expect_true(act("trx_ruv_publication")$display$saved)
  expect_identical(digest::digest(file = file.path(fixture$paths$gene_qc_dir, "ruv_normalised_log_values.tsv"), algo = "sha256"),
    hash_of("results/transcriptomics/ruv_normalised_log_values.tsv"))
  expect_true(all(file.exists(file.path(fixture$paths$publication_graphs_dir, c("canonical_correlation_panels.png", "canonical_correlation_panels.pdf",
    "canonical_correlation_optimization_summary.xlsx")))))
  workbook <- lapply(c("Optimal_Parameters", "Grid_Optimization"), function(sheet) openxlsx::read.xlsx(file.path(fixture$paths$publication_graphs_dir,
    "canonical_correlation_optimization_summary.xlsx"), sheet = sheet))
  expect_equal(workbook[[2L]], as.data.frame(grid$value$grid_summary))
  prints <- reference$store$fingerprints
  for (stage in c("raw", "tmm", "ruvs", "ruviii")) {
    done <- act("trx_qc_figures", list(stage = stage))
    stats <- withWorkflowStateScope(owner, function() resolveWorkflowStateValue(owner$.context, done$display$stats))
    name <- if (identical(stage, "ruviii")) "post_ruviii" else stage
    expect_identical(trxDeclSha(stats$pca$pca_scores), prints$sha256[prints$product == paste0("pca_", name, "_scores")], info = stage)
    expect_identical(trxDeclSha(stats$rle$rle_whiskers), prints$sha256[prints$product == paste0("rle_", name, "_whiskers")], info = stage)
  }
  for (figure in c("normalization_comparison", "ruviii_before_and_after")) expect_true(act("trx_composite_figure", list(figure = figure))$display$saved, info = figure)
  expect_true(all(file.exists(file.path(fixture$paths$publication_graphs_dir, c("normalization_comparison_grid.png", "ruviii_before_and_after.pdf")))))
  expect_true(act("trx_da_publication", list(dataset = "kit095 trx_edge_12"))$display$saved)
  name <- da$manifest$results[[1L]]$comparison
  expect_identical(digest::digest(file = file.path(fixture$paths$da_output_dir, paste0("de_results_", name, ".csv")), algo = "sha256"),
    hash_of(sprintf("results/transcriptomics/de_genes/de_results_%s.csv", name)))
  expect_identical(digest::digest(file = file.path(fixture$paths$da_output_dir, "top_de_genes.tsv"), algo = "sha256"),
    hash_of("results/transcriptomics/de_genes/top_de_genes.tsv"))
  summary <- readLines(file.path(fixture$paths$da_output_dir, "ruv_optimization_summary.txt"))
  expect_identical(summary[[4L]], "Dataset: kit095 trx_edge_12")
  expect_identical(summary[[10L]], paste("Dynamic edgeR Model Formula:", settings$formula_string))
  # The rest is the notebook's text: the timestamp excluded, the dataset line the declared difference (the project's).
  text <- summary[!startsWith(summary, "Timestamp: ")]
  text[text == "Dataset: kit095 trx_edge_12"] <- "Dataset: GSE152840 (Sepsis Transcriptomics Workshop)"
  expect_identical(trxDeclSha(text), prints$sha256[prints$product == "summary_text"])
  expect_true(all(file.exists(file.path(fixture$paths$publication_graphs_dir, c(sprintf("Volcano_Plots/volcano_plot_%s.png", name),
    sprintf("NumSigDeMolecules/summary_barplot_%s.pdf", name))))))
  tab <- withWorkflowStateScope(owner, function() resolveWorkflowStateValue(owner$.context, act("trx_da_tab_results")$display$tab))
  expect_identical(tab$significant_counts, da$results$significant_counts)
  pca <- act("trx_final_pca", list(state = owner$current_state, grouping_variable = "group", shape_variable = NULL))
  expect_null(pca$display$message)
  corrected <- withWorkflowStateScope(owner, function() trxSavedMatrix(owner, ruvs$value$corrected))
  expect_identical(trxStateCorrectedLogCpm(trxAssayInputs(readTrxWorkflowState(owner))), corrected)
  # The Glimma page: saved, or its failure recorded with the reason.
  volcano <- act("trx_volcano_page", list(contrast = names(da$manifest$results)[[1L]], assay = "counts", qvalue = 0.05))
  expect_true(isSavedWorkflowReference(volcano$display$page) || (is.character(volcano$display$message) && nzchar(volcano$display$message)))
  expect_true(isSavedWorkflowReference(volcano$display$page), info = volcano$display$message)
  for (path in file.path(c(fixture$paths$integration_dir, file.path(fixture$paths$integration_dir, "Transcriptomics")), "transcriptomics_normalized_data.tsv")) {
    dir.create(dirname(path), showWarnings = FALSE)
    act("trx_integration_export", list(path = path, state_name = owner$current_state))
    expect_identical(digest::digest(file = path, algo = "sha256"), hash_of("results/transcriptomics/transcriptomics_normalized_data.tsv"))
  }
  staged <- new.env()
  collect <- function(workflowData, normData, input) list(r6_current_state_name = owner$current_state, current_s4_object = NULL, omic_type = "transcriptomics")
  startWorkflowManagedSessionExport(owner, trxDeclSyncRunner(fixture, controller, staged), fixture$workflow, NULL, NULL, fixture$source_dir, collect, "trx_",
    publish = function() invisible(NULL), failed = function(error) stop(error))
  expect_equal(staged$session$current_s4_object, owner$getState(owner$current_state))
  # None of the publications committed a state.
  expect_identical(unlist(owner$getHistory()), c("initial", "trx_raw_data_s4", "trx_expression_filtered", "trx_qc_complete", "trx_ruvs"))
})

test_that("KIT-095: the registry shape of the transcriptomics binding, and every wrapper refuses a manager that is not the owner", {
  binding <- transcriptomicsDeclarationBinding
  expect_identical(omicsDeclarationBindings()$transcriptomics, binding)
  steps <- c("trx_expression_filter", "trx_qc_finalize", "trx_tmm_logcpm", "trx_ruv_grid", "trx_ruvs", "trx_da_contrast", "trx_ruviii_export", "trx_reset")
  publications <- paste0("trx_", c("filtering_metrics", "ruv_publication", "qc_figures", "composite_figure", "da_publication", "da_tab_results",
    "final_pca", "volcano_page", "session_export", "integration_export"))
  names <- Filter(function(name) identical(workflowOperationEntry(name)$omic_type, "transcriptomics"), workflowOperationNames())
  expect_setequal(names, c(steps, publications))
  worker <- c("trx_expression_filter", "trx_tmm_logcpm", "trx_ruv_grid", "trx_ruvs", "trx_da_contrast", "trx_ruviii_export")
  for (name in names) {
    entry <- workflowOperationEntry(name)
    expect_identical(entry$input_class, "TranscriptomicsAssayData", info = name)
    expect_identical(entry$default_variant, "whole", info = name)
    expect_identical(entry$compute, if (name %in% setdiff(steps, worker)) "in_process" else "omics", info = name)
    # No transcriptomics kernel writes an `@args` section: every worker step and publication declares `section = NULL`.
    if (!identical(entry$compute, "in_process")) expect_true("section" %in% names(entry) && is.null(entry$section), info = name)
    expect_setequal(names(entry$defaults) %||% character(), entry$parameters)
    expect_setequal(names(workflowParameterLabels[[name]]) %||% character(), entry$parameters)
  }
  expect_identical(vapply(steps, function(name) workflowOperationEntry(name)$step_label %||% NA_character_, ""),
    c(trx_expression_filter = "2_Expression_Filtered", trx_qc_finalize = "3_QC_Complete", trx_tmm_logcpm = NA, trx_ruv_grid = NA,
      trx_ruvs = "4_RUVs", trx_da_contrast = NA, trx_ruviii_export = NA, trx_reset = "0_Reset"))
  expect_identical(workflowOperationEntry("trx_ruv_grid")$parameters, c("logcpm", "ruv_mode", "percentage_min", "percentage_max",
    "percentage_step", "k_penalty_weight", "manual_k", "manual_percentage"))
  expect_identical(workflowOperationEntry("trx_ruviii_export")$parameters, c("logcpm", "grid", "ruviii_m_source"))
  expect_identical(workflowOperationEntry("trx_da_contrast")$declaration$fixed$model, "edgeR_QL")
  expect_null(workflowOperationEntry("trx_da_contrast")$declaration$fixed$eBayes_trend)
  expect_identical(binding[c("dependency_parameters", "native_steps", "filtering_steps", "design_input_prefix", "qc_dir_role")],
    list(dependency_parameters = c("logcpm", "grid"), native_steps = "trx_da_contrast", filtering_steps = "trx_expression_filter",
      design_input_prefix = "trx_counts_", qc_dir_role = "gene_qc_dir"))
  output <- workflowOutputDeclaration("trx.da_results.v2")
  expect_identical(output[c("operation", "long_table", "identifies", "version")],
    list(operation = "trx_da_contrast", long_table = "da_genes_long", identifies = "isTrxResultManifest", version = 2L))
  expect_identical(workflowRouteForType("TRANSCRIPTOMICS"), "transcriptomics")
  expect_true("ruv_mode" %in% names(workflowRecordReportSources) && "ruviii_m_source" %in% names(workflowRecordReportSources))
  # Every wrapper refuses a manager that is not the session's owner.
  not_owner <- list(.context = list())
  object <- methods::new("TranscriptomicsAssayData")
  refusals <- list(
    function() runTrxReferenceExpressionFilter(not_owner, list()),
    function() runTrxReferenceRuvs(not_owner, trxRuvGridParameters(), list()),
    function() runTrxReferenceDaAnalysis(not_owner, NULL, list(), data.frame()),
    function() runTrxReferenceRuviiiExport(not_owner, "group"),
    function() acceptTrxRuviiiExport(not_owner, NULL, list()),
    function() recordTrxDesignProvenance(not_owner, object),
    function() readTrxWorkflowState(not_owner),
    function() trxRealisedFormula(not_owner, NULL),
    function() commitWorkflowOmicsStoreStep(not_owner, "trx_reset", "trx_reset", list(), "reset"))
  for (i in seq_along(refusals)) expect_error(refusals[[i]](), class = "saved_entity_scope", info = i)
})
