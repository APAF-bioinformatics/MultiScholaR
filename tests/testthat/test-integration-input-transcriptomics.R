# KIT-136: transcriptomics writes the declared integration input (KIT-061's contract, utils_workflow_integration_input.R)
# beside today's `transcriptomics_normalized_data.tsv`: the RUVIII export where RUV was applied (INT-1), the TMM logCPM
# where it was not (design `none`, KIT-129), in process at the Summary hook and on the worker (`trx_integration_export`).
# Fixture: the KIT-040 corpus trx_edge_12 on a session owner (helper-workflow-trx-declarations.R).

integrationTrxQuiet <- function(expr) {
    value <- NULL
    utils::capture.output(value <- suppressWarnings(suppressMessages(expr)))
    value
}

integrationTrxOwner <- function(env = parent.frame()) {
    skipWithoutTrxDeclReference()
    reference <- trxDeclReference("trx_edge_12", env = env)
    fixture <- trxDeclarationsOwner(withr::local_tempdir(.local_envir = env), reference$dir)
    withr::defer(closeDIAStateAdapter(fixture$owner), envir = env)
    integrationTrxQuiet({
        runTrxReferenceExpressionFilter(fixture$owner, list())
        commitWorkflowOmicsStoreStep(fixture$owner, "trx_qc_finalize", "trx_qc_complete", list(), "QC processing complete")
    })
    fixture
}

integrationTrxSave <- function(fixture, root, label = "kit136") {
    paths <- list(source_dir = file.path(root, "scripts", "transcriptomics"), integration_dir = file.path(root, "integration"),
        base_dir = root)
    dir.create(paths$source_dir, recursive = TRUE, showWarnings = FALSE)
    shown <- character()
    saved <- NULL
    failure <- NULL
    integrationTrxQuiet(saveWorkflowSummaryParameters(omicGuiConfig("transcriptomics")$summary, fixture$workflow,
        projectDirs = list(transcriptomics = paths), omicType = "transcriptomics", label = label, description = "d",
        notify = function(message, ...) shown <<- c(shown, message), done = function(value) saved <<- value,
        failed = function(error) failure <<- error, createWorkflowArgsFn = function(...) "study_parameters.txt",
        writeAnalysisRecordFn = function(...) NULL, catFn = function(...) NULL))
    if (!is.null(failure)) stop(failure)
    list(shown = shown, saved = saved, dir = file.path(paths$integration_dir, "transcriptomics"),
        integration_dir = paths$integration_dir, tsv = file.path(paths$integration_dir, "transcriptomics_normalized_data.tsv"))
}

#' A genes-by-runs matrix as the input's matrix table, in the input's sample order.
integrationTrxExpected <- function(values, runs) {
    data.frame(feature_id = rownames(values), values[, runs, drop = FALSE], check.names = FALSE, stringsAsFactors = FALSE,
        row.names = NULL)
}

test_that("KIT-136: with RUV applied the input is the accepted RUVIII export (INT-1), corrected, beside today's TSV", {
    fixture <- integrationTrxOwner()
    owner <- fixture$owner
    root <- withr::local_tempdir()
    ruvs <- integrationTrxQuiet(runTrxReferenceRuvs(owner, trxRuvGridParameters(), list()))
    # Before the export the save refuses the integration write, as today; no input is written.
    before <- integrationTrxSave(fixture, root)
    expect_identical(before$saved$integrationFile, "not written (no integration export yet)")
    expect_false(dir.exists(before$dir))
    export <- integrationTrxQuiet(runTrxReferenceRuviiiExport(owner, "group"))
    run <- integrationTrxSave(fixture, root)
    expect_identical(run$saved$integrationFile, "transcriptomics_normalized_data.tsv")
    expect_true(file.exists(run$tsv))
    expect_setequal(list.files(run$integration_dir), c("transcriptomics", basename(run$tsv)))
    input <- readIntegrationInputWithMatrices(run$dir)
    manifest <- input$manifest
    expect_identical(names(input$views), "transcriptomics")
    # KIT-153: the role of the data set's matrix file (the normalised one); the corrected matrix beside it.
    expect_identical(manifest$matrix_role, "normalised")
    expect_true(input$views$transcriptomics$corrected)
    expect_identical(manifest$normalisation$source$entity_id, ruvs$logcpm$reference$entity_id)
    expect_identical(manifest$scale, "log2_cpm")
    expect_identical(manifest$source_state$name, "trx_ruvs")
    expect_identical(manifest$source_state$entity_id, getWorkflowStateReference(owner, "trx_ruvs")$entity_id)
    expect_identical(manifest$source_state$product[c("type", "entity_id")], list(type = "trx.ruviii", entity_id = export$reference$entity_id))
    correction <- manifest$correction
    expect_identical(correction[c("state", "id", "method", "m_source", "k", "control_count", "replicate_structure",
        "grouping_variable", "matrix", "da_tests")], list(state = "trx_ruvs", id = "C3b", method = "RUVIII", m_source = "group",
        k = as.integer(ruvs$value$k), control_count = as.integer(sum(export$value$ctl_mask)), replicate_structure = "group",
        grouping_variable = "group", matrix = "ruviii_export_values", da_tests = "counts_with_ruvs_factors_in_model"))
    expect_identical(unlist(correction$grouping_levels), unique(as.character(input$samples$group)))
    # Today's grid names no new-design field.
    expect_false("ruv_design" %in% names(correction))
    expect_equal(input$views$transcriptomics$corrected_matrix, integrationTrxExpected(t(export$value$newY), input$samples$Run),
        tolerance = 1e-12, ignore_attr = TRUE)
    # KIT-153: the data set's own matrix is the TMM logCPM the export's RUV-III read; RUV-III with the export's
    # replicate matrix, controls and k, run again on it, gives the corrected matrix.
    log_cpm <- withWorkflowStateScope(owner, function() trxSavedMatrix(owner, ruvs$logcpm$value$log_cpm))
    expect_equal(input$views$transcriptomics$matrix, integrationTrxExpected(log_cpm, input$samples$Run), tolerance = 0,
        ignore_attr = TRUE)
    uncorrected <- as.matrix(input$views$transcriptomics$matrix[-1L])
    rownames(uncorrected) <- input$views$transcriptomics$matrix$feature_id
    uncorrected <- uncorrected[rownames(log_cpm), rownames(export$value$newY), drop = FALSE]
    again <- ruv::RUVIII(Y = t(uncorrected), M = export$value$M, ctl = export$value$ctl_mask, k = export$value$k,
        return.info = TRUE)$newY
    expect_equal(unname(again), unname(export$value$newY), tolerance = 1e-10)
    # KIT-184 K4 (gate G1): the written factors are what RUV-III removed: they span RUV-III's own W, recomputed here from
    # its `fullalpha` (the RUVSeq factors of the state's design are the DA model's, another method's).
    factors <- expectIntegrationRemovedFactors(run$dir, "transcriptomics", export$value$k, input = input)
    info <- ruv::RUVIII(Y = t(uncorrected), M = export$value$M, ctl = export$value$ctl_mask, k = export$value$k,
        return.info = TRUE)
    alpha <- info$fullalpha[seq_len(export$value$k), , drop = FALSE]
    ac <- alpha[, export$value$ctl_mask, drop = FALSE]
    w <- t(uncorrected)[, export$value$ctl_mask, drop = FALSE] %*% t(ac) %*% solve(ac %*% t(ac))
    w <- scale(w[input$samples$Run, , drop = FALSE], scale = FALSE)
    expect_lt(max(abs(w - factors %*% crossprod(factors, w))) / max(abs(w)), 1e-8)
    # The samples: the state's design without the RUVs factors (the DA model's, not the export's).
    expect_false(any(startsWith(names(input$samples), "W_")))
    expect_true(manifest$biological_sample_column %in% names(input$samples))
    # The features: the gene id and its name from the count table.
    counts <- fixture$object@count_data[[1L]]
    features <- input$views$transcriptomics$features
    expect_identical(features$GeneName, as.character(counts$GeneName)[match(features$feature_id, as.character(counts$GeneID))])
    expect_identical(features$display_name, ifelse(is.na(features$GeneName) | !nzchar(features$GeneName), features$feature_id,
        features$GeneName))
})

test_that("KIT-136: a new design names its RUV record fields; an export from an earlier RUV step is refused until run again", {
    fixture <- integrationTrxOwner()
    owner <- fixture$owner
    root <- withr::local_tempdir()
    integrationTrxQuiet(runTrxReferenceRuvs(owner, trxRuvGridParameters(), list()))
    integrationTrxQuiet(runTrxReferenceRuviiiExport(owner, "group"))
    # The RUV step run again (a new design): the accepted export is the earlier state's.
    integrationTrxQuiet(runTrxReferenceRuvs(owner, NULL, list(), design_parameters = trxRuvDesignParameters("groups_all_features", 2L, "user")))
    expect_error(integrationTrxSave(fixture, root), class = "integration_input_invalid", regexp = "earlier RUV step")
    export <- integrationTrxQuiet(runTrxReferenceRuviiiExport(owner, "group"))
    trx_dir <- integrationTrxSave(fixture, root)$dir
    input <- readIntegrationInputWithMatrices(trx_dir)
    # KIT-184 K4 (gate G1): the declared design's written factors are what RUV-III removed.
    expectIntegrationRemovedFactors(trx_dir, "transcriptomics", 2L, input = input)
    correction <- input$manifest$correction
    expect_identical(input$manifest$matrix_role, "normalised")
    expect_true(all(RUV_DESIGN_RECORD_FIELDS %in% names(correction)))
    expect_identical(correction[c("k", "ruv_design", "ruv_k_rule", "ruv_k", "ruv_da_method")], list(k = 2L,
        ruv_design = "groups_all_features", ruv_k_rule = "user", ruv_k = 2L, ruv_da_method = "edger_qlf_factors_in_model"))
    expect_equal(input$views$transcriptomics$corrected_matrix, integrationTrxExpected(t(export$value$newY), input$samples$Run),
        tolerance = 1e-12, ignore_attr = TRUE)
})

test_that("KIT-136: with RUV off (design none) the TMM logCPM is written, normalised, without a RUVIII export; the worker writes the same", {
    fixture <- integrationTrxOwner()
    owner <- fixture$owner
    root <- withr::local_tempdir()
    ruvs <- integrationTrxQuiet(runTrxReferenceRuvs(owner, NULL, list(), design_parameters = trxRuvDesignParameters("none")))
    expect_null(owner$.context$snapshot$exports$trx_ruviii_export)
    expect_null(trxSummaryIntegrationRequires(fixture$workflow))
    run <- integrationTrxSave(fixture, root)
    expect_identical(run$saved$integrationFile, "transcriptomics_normalized_data.tsv")
    log_cpm <- withWorkflowStateScope(owner, function() trxSavedMatrix(owner, ruvs$logcpm$value$log_cpm))
    # Today's file, which no RUV-off route wrote before, holds the TMM logCPM in its layout (runs by genes, `Run`).
    tsv <- utils::read.delim(run$tsv, check.names = FALSE)
    expect_identical(tsv$Run, colnames(log_cpm))
    expect_equal(unname(as.matrix(tsv[-1L])), unname(t(log_cpm)), tolerance = 1e-6)
    input <- readIntegrationInputWithMatrices(run$dir)
    manifest <- input$manifest
    expect_identical(manifest$matrix_role, "normalised")
    expect_null(manifest$correction)
    expect_false(input$views$transcriptomics$corrected)
    # KIT-184 K4: design none removes nothing: no factors.
    expect_null(manifest$views[[1L]]$removed_factors)
    expect_false("removed_transcriptomics.parquet" %in% names(manifest$files))
    expect_identical(manifest$source_state$product[c("type", "entity_id")],
        list(type = "trx.tmm_logcpm", entity_id = ruvs$logcpm$reference$entity_id))
    expect_equal(input$views$transcriptomics$matrix, integrationTrxExpected(log_cpm, input$samples$Run), tolerance = 1e-12,
        ignore_attr = TRUE)

    # The worker's integration export writes the same input beside its file (the in-test worker under load_all, as
    # test-workflow-omics-application-reads.R's `omicsReadMocks`).
    testthat::local_mocked_bindings(workflowExecutionBuildIdentity = function() strrep("f", 64),
        finishWorkflowExecutionJob = function(owner, status) { owner$job <- NULL; list(status = status) }, .package = "MultiScholaR")
    integration <- file.path(withr::local_tempdir(), "integration")
    dir.create(integration, recursive = TRUE)
    path <- file.path(integration, "transcriptomics_normalized_data.tsv")
    done <- integrationTrxQuiet(omicsWorkerAction(fixture, omicsWorkerController(fixture), "trx_integration_export",
        list(path = path, state_name = owner$current_state)))
    expect_identical(done$acceptance, "complete")
    worker <- readIntegrationInputWithMatrices(file.path(integration, "transcriptomics"))
    drop <- c("written_at", "software")
    expect_identical(worker$manifest[setdiff(names(worker$manifest), drop)], manifest[setdiff(names(manifest), drop)])
    expect_identical(worker$views, input$views)
    expect_identical(worker$samples, input$samples)
})
