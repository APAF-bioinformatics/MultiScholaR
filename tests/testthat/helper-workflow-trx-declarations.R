# KIT-095: the transcriptomics declarations on a session owner, against KIT-040's references
# (docs/qa/workflow-execution/kit/kit-040/reference-<corpus>.rds) and KIT-093's kernels. Used by
# test-workflow-trx-declarations.R. The reference helpers repeat those of test-trx-kernels.R under their own names (each
# test file runs on its own).

TRX_DECL_STORE_SHA256 <- c(
  trx_edge_12 = "f460b376008219c667fdca32bd1184813ff1f5b8ba8c53d31f45c63b728cb12c",
  trx_edge_16 = "ac118c0ada3411ad5510efdd385e0ae08a4b6b3cecfb6c9a4ba436c75a7ecf59",
  trx_screen = "dd058d7bd344e606899f95dac67786ef8a60df473ffede1e71918bdca0d66840")
TRX_DECL_PER_GROUP <- c(trx_edge_12 = 3L, trx_edge_16 = 4L, trx_screen = 6L)

trxDeclRoot <- function() normalizePath(testthat::test_path("..", ".."), mustWork = FALSE)

skipWithoutTrxDeclReference <- function() {
  root <- trxDeclRoot()
  testthat::skip_if_not(file.exists(file.path(root, "docs/qa/workflow-execution/kit/kit-040/reference-trx_screen.rds")) &&
    file.exists(file.path(root, "tools/profiling/kit_omics_corpus.R")), "the KIT-040 references are not shipped (R CMD check)")
  for (package in c("edgeR", "RUVSeq", "EDASeq", "ruv")) testthat::skip_if_not_installed(package)
  testthat::skip_if_not(workflowOwnerCapable())
}

trxDeclSha <- function(x) digest::digest(serialize(x, NULL, version = 3L), algo = "sha256", serialize = FALSE)

#' The KIT-040 store (hash-checked) and its corpus regenerated (manifest hash-checked).
trxDeclReference <- function(corpus, env = parent.frame()) {
  path <- file.path(trxDeclRoot(), "docs/qa/workflow-execution/kit/kit-040", sprintf("reference-%s.rds", corpus))
  testthat::expect_identical(digest::digest(file = path, algo = "sha256"), TRX_DECL_STORE_SHA256[[corpus]])
  store <- readRDS(path)
  h <- new.env(parent = globalenv())
  sys.source(file.path(trxDeclRoot(), "tools/profiling/kit_omics_corpus.R"), h)
  dir <- withr::local_tempdir(.local_envir = env)
  h$writeKitTrxCorpus(dir, h$kitOmicsRecipe("trx_screen", TRX_DECL_PER_GROUP[[corpus]]))
  testthat::expect_identical(digest::digest(file = file.path(dir, "manifest.json"), algo = "sha256"), store$corpus_manifest_sha256)
  list(store = store, dir = dir, corpus = corpus)
}

#' Every named product's sha256 equals its KIT-040 fingerprint.
expectTrxDeclProducts <- function(products, reference, what = "") {
  prints <- reference$store$fingerprints
  for (name in names(products)) {
    testthat::expect_identical(trxDeclSha(products[[name]]), prints$sha256[prints$product == name], info = paste(reference$corpus, what, name))
  }
}

#' A transcriptomics session owner on a corpus: the counts and the design read as the notebook reads them (rmd:121, 136),
#' the import's source copy and the design file on disk, the design state committed with its provenance and inputs.
#' @param design A design to use instead of the corpus's `design_matrix.tab` (written to the design file).
trxDeclarationsOwner <- function(root, dir, design = NULL, technical_replicate_id = NA_character_) {
  source_dir <- file.path(root, "source")
  dir.create(file.path(source_dir, "import_sources"), recursive = TRUE, showWarnings = FALSE)
  file.copy(file.path(dir, "counts.tsv"), file.path(source_dir, "import_sources", "counts.tsv"))
  counts <- vroom::vroom(file.path(dir, "counts.tsv"), show_col_types = FALSE)
  if (is.null(design)) {
    file.copy(file.path(dir, "design_matrix.tab"), file.path(source_dir, "design_matrix.tab"))
    design <- vroom::vroom(file.path(dir, "design_matrix.tab"), show_col_types = FALSE)
  } else {
    utils::write.table(design, file.path(source_dir, "design_matrix.tab"), sep = "\t", quote = FALSE, row.names = FALSE)
  }
  paths <- list(source_dir = source_dir, gene_qc_dir = file.path(root, "gene_qc"), da_output_dir = file.path(root, "da_genes"),
    publication_graphs_dir = file.path(root, "publication_graphs"), integration_dir = file.path(root, "integration"))
  for (p in paths[-1]) dir.create(p, recursive = TRUE, showWarnings = FALSE)
  session <- omicsOwnerSession()
  contrasts <- trxContrastsTable(trxParseContrastFile(file.path(dir, "contrast_strings.tab")))
  workflow <- list2env(list(state_update_trigger = NULL, config_list = list(), contrasts_tbl = contrasts,
    tab_status = list(quality_control = "pending", normalization = "pending")))
  workflow$state_manager <- createDIASessionOwner(workflow, paths, session, omic_type = "transcriptomics")
  shiny::withReactiveDomain(session, bindWorkflowOwnerTypeAtDesign(workflow, "TRANSCRIPTOMICS", "trx_counts"))
  object <- createTranscriptomicsAssayData(list(counts = counts), design, c("RPMI", "Sera"),
    roles = list(technical_replicate_id = technical_replicate_id))
  owner <- workflow$state_manager
  recordTrxDesignProvenance(owner, object, source_dir = source_dir)
  owner$saveState("trx_raw_data_s4", object, list(), "Initial TranscriptomicsAssayData S4 object created after design matrix")
  list(owner = owner, workflow = workflow, session = session, paths = paths, source_dir = source_dir, object = object, dir = dir)
}

#' The corpus's design with its technical-replicate sets (design.tsv's tech_rep_group, joined by run).
trxDeclTechnicalDesign <- function(dir) {
  design <- vroom::vroom(file.path(dir, "design_matrix.tab"), show_col_types = FALSE)
  tech <- utils::read.delim(file.path(dir, "design.tsv"))
  design$tech_rep_group <- tech$tech_rep_group[match(design$Run, tech$sample)]
  design
}

#' The whole journey on the owner through the wrappers (design, filter, QC finalise, TMM/logCPM, grid, RUVs, DA, RUVIII).
trxDeclarationsJourney <- function(fixture, grid = trxRuvGridParameters(), m_source = "group") {
  owner <- fixture$owner
  config <- fixture$workflow$config_list
  out <- list()
  out$filter <- runTrxReferenceExpressionFilter(owner, config)
  commitWorkflowOmicsStoreStep(owner, "trx_qc_finalize", "trx_qc_complete", config, "QC processing complete - ready for normalization")
  out$ruvs <- runTrxReferenceRuvs(owner, grid, config)
  out$da <- runTrxReferenceDaAnalysis(owner, getWorkflowStateReference(owner), list(da_q_val_thresh = 0.05, treat_lfc_cutoff = 1),
    fixture$workflow$contrasts_tbl)
  out$ruviii <- runTrxReferenceRuviiiExport(owner, m_source)
  out
}

#' A saved product's value, from the owner's store.
trxDeclValue <- function(owner, reference) hydrateSavedWorkflowEntity(owner$.context$store, reference)$value

#' The driver's grid_members shape (test-trx-kernels.R's trxTestGridMembers).
trxDeclGridMembers <- function(grid, percentage_range = seq(5, 30, by = 5)) {
  ok <- stats::setNames(grid$grid_members, vapply(grid$grid_members, function(r) as.character(r$percentage), ""))
  do.call(rbind, lapply(percentage_range, function(pct) {
    r <- ok[[as.character(pct)]]
    if (is.null(r)) return(data.frame(percentage = pct, status = "failed", n_controls = NA_real_, best_k = NA_integer_,
      separation_score = NA_real_, k_penalty = NA_real_, composite_score = NA_real_))
    data.frame(percentage = r$percentage, status = "complete", n_controls = r$n_controls, best_k = r$best_k,
      separation_score = r$separation_score, k_penalty = r$k_penalty, composite_score = r$composite_score)
  }))
}

#' The journey's states and products under their KIT-040 reference names (kit-040/product-map.tsv).
trxDeclJourneyProducts <- function(fixture, out) {
  owner <- fixture$owner
  read <- function(name) trxAssayInputs(withWorkflowStateScope(owner, function() readTrxWorkflowState(owner, name)))
  logcpm <- out$ruvs$logcpm$value
  grid <- out$ruvs$grid$value
  ruvs <- out$ruvs$value
  table <- function(reference) withWorkflowStateScope(owner, function() trxSavedMatrix(owner, reference))
  contrast <- out$da$manifest$results[[1L]]
  applied <- trxDeclValue(owner, contrast$products$da_contrast)$assays[[1L]]$results
  # A DA fit hydrates as its envelope.
  fit <- hydrateSavedWorkflowEntity(owner$.context$store, contrast$assays[[1L]]$fit)
  members <- grid$grid_members
  pct <- vapply(members, function(r) paste0("pct_", r$percentage), "")
  ruviii <- out$ruviii$value
  list(filter_effective_lib_size = out$filter$value$effective_lib_size, filter_min_sample_size = out$filter$value$min_sample_size,
    filter_median_lib_size = out$filter$value$median_lib_size, count_matrix_filtered = read("trx_expression_filtered")$counts,
    dge_filtered_samples = logcpm$samples, log_cpm_tmm = table(logcpm$log_cpm),
    ordered_features_idx = grid$prelim_order, max_acceptable_k = grid$max_acceptable_k, grid_members = trxDeclGridMembers(grid),
    grid_cancor_data = stats::setNames(lapply(members, `[[`, "cancor_data"), pct),
    grid_control_genes = stats::setNames(lapply(members, `[[`, "control_genes"), pct), df_summary = grid$grid_summary,
    best_res = grid$best[c("percentage", "n_controls", "best_k", "separation_score", "k_penalty", "composite_score", "ctl_idx", "control_genes")],
    best_control_genes = grid$control_genes, w_matrix = ruvs$W, ruvs_counts = read("trx_ruvs")$counts,
    log_cpm_corrected = table(ruvs$corrected), design_final = fit$fit$design, design_formula_str = out$da$manifest$parameters$formula,
    contrast_entry = names(out$da$manifest$results)[[1L]], contrasts_matrix = applied$contrasts_matrix, fit_final = fit$fit,
    qlf = applied$qlf, qlf_table = applied$qlf$table, top_genes = applied$top_table,
    ruviii_newY = ruviii$newY, ruviii_M = ruviii$M, ruviii_ctl_mask = ruviii$ctl_mask, best_k = ruviii$k)
}
