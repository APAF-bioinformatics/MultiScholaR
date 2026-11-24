# KIT-186 (VAL-59): the transcriptomics DA tab on a normalised state without RUV factors (the RUV design `none`, or a
# declared design whose k resolves to 0). Such a state is ready for the DA when its design has W_ columns (as before) or
# its record carries the RUVs step's product (`trx_ruvs_science`); a state never normalised is still refused. The DA is
# edgeR's quasi-likelihood analysis on ~ 0 + group, checked here against edgeR's own functions.
# Every test that fails before the fix computes what derives from the DA result only when that result exists, so that it
# fails at an expectation, not at an error (ticket section 9).

trxNoRuvNotes <- function(out, type = "error") {
    shown <- Filter(function(note) identical(note$type, type), out$notes)
    vapply(shown, function(note) note$message, character(1))
}

trxNoRuvRuvsOutcomes <- function(record) {
    steps <- Filter(function(step) identical(step$operation, "trx_ruvs"), record$steps)
    if (length(steps)) steps[[length(steps)]]$outcomes else list()
}

trxNoRuvDesignColumns <- function(owner) {
    names(readWorkflowOmicsMetadata(transcriptomicsDeclarationBinding, owner)$design_matrix)
}

trxNoRuvTop <- function(out) {
    first <- if (!is.null(out$da)) module_ci_trx_first_contrast(out$da)
    if (!is.null(first)) module_ci_trx_top_table(first)
}

trxNoRuvNotNormalised <- "Run the normalisation (RUVs) first; the formula is realised from its design."
trxNoRuvRefusalW1 <- "Run the normalisation first: the differential analysis runs on the normalised data."

#' The first contrast's table computed with edgeR's own functions on the counts and design of `trx_qc_complete` (no
#' package kernel: neither trxDaFit nor trxDaApply).
trxNoRuvEdgeRTable <- function(owner) {
    qc <- withWorkflowStateScope(owner, function() readTrxWorkflowState(owner, "trx_qc_complete"))
    counts <- trxAssayInputs(qc)$counts
    pdata <- trxDaDesign(qc)
    y <- edgeR::calcNormFactors(edgeR::DGEList(counts, group = pdata$group), method = "TMM")
    X <- stats::model.matrix(~ 0 + group, data = pdata)
    fit <- edgeR::glmQLFit(edgeR::estimateDisp(y, X, robust = TRUE), X, robust = TRUE)
    edgeR::topTags(edgeR::glmQLFTest(fit,
        contrast = as.numeric(colnames(X) == "groupSera") - as.numeric(colnames(X) == "groupRPMI")),
    n = Inf, sort.by = "PValue")$table
}

test_that("KIT-186 T1 T2: design none, in process and managed: the DA runs on ~ 0 + group and is edgeR QL", {
    skipWithoutTrxDeclReference()
    withr::local_options(MultiScholaR.ruv_design_panel = TRUE)
    # T1: in process.
    fixture <- module_ci_trx_project(module_ci_trx_corpus())
    owner <- fixture$owner
    out <- module_ci_trx_tabs_journey(fixture, norm = list(ruv_mode = "automatic", ruv_design = "none"))
    expect_identical(out$errors, character())
    expect_identical(owner$current_state, "trx_ruvs")
    expect_false(any(startsWith(trxNoRuvDesignColumns(owner), "W_")))
    record <- buildWorkflowAnalysisRecord(owner)
    ruv <- trxNoRuvRuvsOutcomes(record)
    expect_identical(ruv[c("ruv_design", "ruv_applied", "ruv_k", "ruv_da_method")],
        list(ruv_design = "none", ruv_applied = FALSE, ruv_k = 0L, ruv_da_method = "edger_qlf_none"))
    da <- record$analyses$differential_abundance
    expect_false(is.null(da))
    expect_identical(da$model$formula_used, "~ 0 + group")
    expect_identical(length(da$model$w_columns), 0L)
    expect_identical(da$source_state, "trx_ruvs")
    expect_identical(fixture$workflow$tab_status$differential_analysis, "complete")
    expect_identical(trxNoRuvNotes(out), character())
    expect_true(startsWith(trxDaModelText(fixture$workflow), "Formula: ~ 0 + group\n"))
    expect_false(is.null(out$da))
    top <- trxNoRuvTop(out)
    expected <- trxNoRuvEdgeRTable(owner)
    columns <- c("logFC", "logCPM", "F", "PValue", "FDR")
    expect_identical(rownames(top), rownames(expected))
    expect_equal(top[columns], expected[columns], tolerance = 1e-10)

    # T2: managed, on a second fixture; the same table, the worker's publications delivered and written.
    managed <- module_ci_trx_project(module_ci_trx_corpus())
    probe <- new.env()
    omicsManagedChainMocks(managed, probe)
    worker <- module_ci_trx_tabs_journey(managed, managed = TRUE, norm = list(ruv_mode = "automatic", ruv_design = "none"))
    expect_identical(worker$errors, character())
    expect_identical(trxNoRuvNotes(worker), character())
    expect_true(all(c("trx_da_contrast", "trx_da_tab_results", "trx_da_publication") %in% probe$operations),
        info = paste(probe$operations, collapse = " "))
    expect_true(all(probe$statuses == "delivered"))
    expect_false(is.null(worker$da))
    expect_identical(worker$da$da_genes_long, out$da$da_genes_long)
    expect_true(all(file.exists(file.path(managed$paths$da_output_dir, c("de_results_Sera_vs_RPMI.csv", "top_de_genes.tsv")))))
})

test_that("KIT-186 T3: a declared design whose estimated k is 0 runs the same DA as the design none", {
    skipWithoutTrxDeclReference()
    withr::local_options(MultiScholaR.ruv_design_panel = TRUE)
    reference <- module_ci_trx_project(module_ci_trx_corpus())
    none <- module_ci_trx_tabs_journey(reference, norm = list(ruv_mode = "automatic", ruv_design = "none"))
    fixture <- module_ci_trx_project(module_ci_trx_corpus())
    # The estimated count forced to 0 and nothing else.
    real <- estimateRuvFactorCount
    testthat::local_mocked_bindings(estimateRuvFactorCount = function(...) {
        e <- real(...)
        e$k <- 0L
        e
    }, .package = "MultiScholaR")
    out <- module_ci_trx_tabs_journey(fixture, norm = list(ruv_mode = "automatic", ruv_design = "groups_all_features",
        ruv_k_rule = "num_sv_be"))
    expect_identical(out$errors, character())
    expect_identical(trxNoRuvNotes(out), character())
    record <- buildWorkflowAnalysisRecord(fixture$owner)
    ruv <- trxNoRuvRuvsOutcomes(record)
    expect_identical(ruv[c("ruv_design", "ruv_k_rule", "ruv_k", "ruv_applied")],
        list(ruv_design = "groups_all_features", ruv_k_rule = "num_sv_be", ruv_k = 0L, ruv_applied = FALSE))
    expect_identical(record$analyses$differential_abundance$model$formula_used, "~ 0 + group")
    expect_false(is.null(out$da))
    expect_false(is.null(none$da))
    expect_identical(out$da$da_genes_long, none$da$da_genes_long)
})

test_that("KIT-186 T4: a state never normalised is refused; a legacy-style state with W_ columns is accepted", {
    skipWithoutTrxDeclReference()
    withr::local_options(MultiScholaR.ruv_design_panel = TRUE)
    # (a) at trx_qc_complete
    fixture <- module_ci_trx_project(module_ci_trx_corpus())
    expect_identical(fixture$owner$current_state, "trx_qc_complete")
    expect_null(trxDaReadMetadata(fixture$workflow))
    expect_identical(trxDaModelText(fixture$workflow), trxNoRuvNotNormalised)
    # (b) after a none run, the tab's reset to the state the run started from
    reset <- module_ci_trx_project(module_ci_trx_corpus())
    owner <- reset$owner
    out <- module_ci_trx_tabs_journey(reset, norm = list(ruv_mode = "automatic", ruv_design = "none"), da = FALSE)
    expect_identical(out$errors, character())
    expect_identical(owner$current_state, "trx_ruvs")
    resetTrxNormalisation(owner,
        normData = list2env(list(post_filter_obj = getWorkflowStateReference(owner, "trx_qc_complete"))),
        config = list(), addLog = function(...) invisible(NULL), notify = function(...) invisible(NULL),
        setStatus = function(...) invisible(NULL))
    expect_identical(owner$current_state, "trx_reset")
    expect_null(trxDaReadMetadata(reset$workflow))
    expect_identical(trxDaModelText(reset$workflow), trxNoRuvNotNormalised)
    # (c) a design with a W_ column on a state that carries no RUVs product (a legacy-style state): accepted as today
    legacy <- module_ci_trx_project(module_ci_trx_corpus())
    lowner <- legacy$owner
    object <- withWorkflowStateScope(lowner, function() readTrxWorkflowState(lowner, "trx_qc_complete"))
    object@design_matrix$W_1 <- as.numeric(seq_len(nrow(object@design_matrix)))
    lowner$saveState("trx_legacy_ruvs", object, list(), "A legacy-style state: W_1 in the design, no RUVs product")
    expect_identical(lowner$current_state, "trx_legacy_ruvs")
    product <- withWorkflowStateScope(lowner, function() {
        readWorkflowStateMetadata(lowner, getWorkflowStateReference(lowner))$data$products$trx_ruvs_science
    })
    expect_null(product)
    meta <- trxDaReadMetadata(legacy$workflow)
    expect_false(is.null(meta))
    expect_true("W_1" %in% names(meta$design_matrix))
    # (d) the refusal's text
    expect_identical(omicGuiConfig("transcriptomics")$da$source$no_data, trxNoRuvRefusalW1)
})

test_that("KIT-186 T5: a RUV head (today's grid, gate off) takes today's branch", {
    skipWithoutTrxDeclReference()
    withr::local_options(MultiScholaR.ruv_design_panel = FALSE)
    fixture <- module_ci_trx_project(module_ci_trx_corpus())
    owner <- fixture$owner
    out <- module_ci_trx_tabs_journey(fixture, da = FALSE)
    expect_identical(out$errors, character())
    expect_identical(owner$current_state, "trx_ruvs")
    expect_identical(trxDaReadMetadata(fixture$workflow), readWorkflowOmicsMetadata(transcriptomicsDeclarationBinding, owner))
    expect_true(startsWith(trxDaModelText(fixture$workflow), "Formula: ~ 0 + group + W_1"))
})

test_that("KIT-186 T6: the scripted path and replay of a none run: the formula check passes and the DA replays exactly", {
    skipWithoutTrxDeclReference()
    reference <- trxDeclReference("trx_edge_12")
    fixture <- trxDeclarationsOwner(withr::local_tempdir(), reference$dir)
    withr::defer(closeDIAStateAdapter(fixture$owner))
    owner <- fixture$owner
    testthat::local_mocked_bindings(workflowExecutionBuildIdentity = function() strrep("f", 64),
        finishWorkflowExecutionJob = function(owner, status) {
            owner$job <- NULL
            list(status = status)
        }, .package = "MultiScholaR")
    suppressWarnings(suppressMessages(runTrxReferenceExpressionFilter(owner, list())))
    commitWorkflowOmicsStoreStep(owner, "trx_qc_finalize", "trx_qc_complete", list(), "QC processing complete - ready for normalization")
    run <- suppressWarnings(suppressMessages(runTrxReferenceRuvs(owner, NULL, list(),
        design_parameters = trxRuvDesignParameters("none"))))
    expect_identical(run$value$k, 0L)
    da <- suppressWarnings(suppressMessages(runTrxReferenceDaAnalysis(owner, getWorkflowStateReference(owner),
        list(da_q_val_thresh = 0.05, treat_lfc_cutoff = 1), fixture$workflow$contrasts_tbl)))
    manifest <- da$manifest
    expect_identical(manifest$parameters$formula, "~ 0 + group")
    recorded <- trxDeclValue(owner, manifest$results[[1L]]$products$da_contrast)$assays[[1L]]$results$top_table
    # The formula check the replay makes, on the head.
    expect_null(withWorkflowStateScope(owner, function() trxDaReplayCheck(owner, list(formula_string = "~ 0 + group"))))
    history <- unlist(owner$getHistory())
    json <- file.path(withr::local_tempdir(), "analysis_record.json")
    writeWorkflowAnalysisRecord(owner, json)
    record <- jsonlite::read_json(json, simplifyVector = FALSE)
    # The replay on a copy cut to the design state (KIT-095's way).
    harness <- new.env(parent = globalenv())
    sys.source(file.path(trxDeclRoot(), "tools", "profiling", "workflow_kit.R"), harness)
    api <- function(name) get(name, asNamespace("MultiScholaR"))
    copy <- withr::local_tempdir()
    stopifnot(system2("cp", c("-a", file.path(fixture$source_dir, "."), copy)) == 0L)
    store_root <- dirname(dirname(list.files(fixture$source_dir, "^store[.]json$", recursive = TRUE, all.files = TRUE,
        full.names = TRUE)[[1L]]))
    copy_root <- file.path(copy, substring(store_root, nchar(fixture$source_dir) + 2L))
    prepared <- harness$kitPrepareReplayProject(copy_root, c("initial", "trx_raw_data_s4"), api)
    inputs <- c(trx_counts_1 = file.path(fixture$source_dir, "import_sources", "counts.tsv"),
        design_file = file.path(fixture$source_dir, "design_matrix.tab"))
    replay <- suppressWarnings(suppressMessages(harness$kitReplayJourney(record, prepared, store_root, inputs, fixture$paths,
        withr::local_tempdir(), api)))
    expect_identical(replay$states, history)
    expect_identical(as.character(replay$mismatches), character())
    identity <- readSavedStoreIdentity(file.path(prepared, ".workflow-entities", "store.json"))
    replayed <- newWorkflowStateOwner(newWorkflowEntityStore(prepared, identity$project_id, identity$workflow_id), reopen = TRUE)
    withr::defer(closeWorkflowStateOwner(replayed))
    again <- resolveWorkflowStateValue(replayed$.context, replayed$.context$snapshot$da_results)
    expect_identical(again$parameters, manifest$parameters)
    expect_identical(trxDeclValue(replayed, again$results[[1L]]$products$da_contrast)$assays[[1L]]$results$top_table, recorded)
})
