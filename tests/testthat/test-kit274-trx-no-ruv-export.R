# KIT-274 (VAL-318, VAL-332): in transcriptomics, "Run RUVIII Export" is refused when the head's RUVs state applied no
# correction for unwanted variation (the design `none`, or a k rule that found 0); an export made in that state by an
# earlier build is not listed in the analysis record (the record's own last `trx_ruvs` step decides) or in the report (the
# template's guard on the same rule); a replay of a record that lists one skips it with a statement and a warning; the
# integration input in that state is the TMM logCPM, unchanged.
# Written before the fix (a tests-only commit; nothing under R/ or inst/). The texts below are the owner's decided ones
# (VAL-332), drafts for the wording pass: each is one constant here, so a wording pass is a text edit.
#
# PREDICTIONS (written before the first run; ticket r1.2 section 7). Base = this commit without the fix.
#   T1  design none, then the export in process: an error of class trx_ruviii_no_ruv with the decided message; no export
#       in the snapshot; history and head unchanged.                                      base: no error, export accepted (fails by value)
#   T2  a declared design whose estimated k is forced to 0, then the export in process: the same refusal.   base: accepted (fails)
#   T3  the design none through the mounted tabs in process: the error notification "RUVIII export error: <message>
#       (trx_ruviii_no_ruv)", the status "Not run yet.", no ruviii_before_and_after file, no export.   base: the status reads
#       "RUVIII export complete ... k = 0 ...", no notification (fails)
#   T4  the same, managed (the journey helper's runner on the chain mocks): one error notification equal to T3's, the status
#       "Not run yet.", nothing escaped the observer, no export, the products unchanged by the press.   base: complete (fails)
#   T5  a planted older export (the state of T1 under the old acceptance): the record does not list it.   base: listed (fails)
#   T6  the report chunks on the record an old build wrote (the analysis injected): "No RUVIII integration export was recorded
#       for this analysis." and "Figure not available: ruviii_before_and_after.png"; no "was applied".   base: the applied
#       sentence with k = 0 and the image (fails)
#   T7  replay of that record on its project: completes; the result carries the statement and a warning of class
#       workflow_replay_export_skipped; the replayed snapshot holds no export.   base: the export is repeated, no statement (fails)
#   T8  guard: the planted older export, then Save: succeeds; the integration input is the TMM logCPM, normalised, product
#       trx.tmm_logcpm, not corrected; the file's digest is printed (A5: equal on head).   base: same; head: same
#   T9  guard: after a legacy run the export is accepted, listed with the run's k, described, and the integration input is
#       the export (corrected, product trx.ruviii).   base: same; head: same
#   T10 guard: KIT-264's routes (run k = 2, export, run k = 1; run k = 2, export, run with the design none): no entry.
#       base: same; head: same. R1b (VAL-54 instrument re-send 1): one shared `nd` per project and a precondition that the
#       second run committed a new trx_ruvs entity; with it met T10 passes on base. T1 to T9, T11, T12 use no k274Run: as in R1.
#   T11 guard: a manager that is not the session's owner: saved_entity_scope from both wrappers.   base: same; head: same
#   T12 the call sites: trxRuviiiExportParameters and acceptTrxRuviiiExport on a no-RUV head each raise trx_ruviii_no_ruv.
#       base: no error (fail)
# R1 (this commit): T1 to T7 and T12 fail by value; T8 to T11 pass; X1 of test-kit264-trx-export-current.R prints the base line.

library(testthat)

K274Reason <- "the design none: no correction for unwanted variation"
K274Refusal <- function(reason) paste0("no RUVIII export is made: this analysis applied no correction for unwanted variation (k = 0; ",
    reason, "), so RUVIII would return the TMM logCPM unchanged; the integration uses the TMM logCPM.")
K274Statement <- paste("This record lists a RUVIII export made when no correction for unwanted variation had been applied (k = 0).",
    "The replay does not repeat it: such an export is no longer made, and in this state the integration uses the TMM logCPM,",
    "which the replay writes as before.")
K274None <- "No RUVIII integration export was recorded for this analysis."
K274Missing <- "Figure not available: ruviii_before_and_after.png "
K274Class <- "trx_ruviii_no_ruv"
K274Design <- list(ruv_mode = "automatic", ruv_design = "none")
K274Legacy2 <- list(ruv_mode = "manual", ruv_k = 2L, ruv_percentage = 10)
K274Legacy1 <- list(ruv_mode = "manual", ruv_k = 1L, ruv_percentage = 10)

k274Row <- function(row, ...) cat(sprintf("KIT-274 %s %s\n", row, paste(..., collapse = " ")))

# ------------------------------------------------------------------------------------------------
# Helpers (file-local copies of test-kit264-trx-export-current.R's and of test-integration-input-transcriptomics.R's).
# ------------------------------------------------------------------------------------------------

k274NoneRun <- function(owner) kit172Quiet(runTrxReferenceRuvs(owner, NULL, list(), design_parameters = trxRuvDesignParameters("none")))

#' The export's declaration run on the head, without any acceptance (what `runTrxReferenceRuviiiExport` did before its
#' refusal: its body without the check).
k274Operate <- function(owner, m = "group") kit172Quiet(withWorkflowStateScope(owner, function() {
    source <- getWorkflowStateReference(owner)
    ruvs <- trxStateRuvs(owner, source)
    runSavedOperation(owner, workflowOperationEntry("trx_ruviii_export")$declaration, source,
        list(logcpm = ruvs$logcpm$entity_id, grid = ruvs$grid$entity_id, ruviii_m_source = m), list(ruvs$logcpm, ruvs$grid))
}))

#' `acceptTrxRuviiiExport`'s body without the check: the state an earlier build leaves.
k274Accept <- function(owner, source, run) kit172Quiet(withWorkflowStateScope(owner, function() {
    context <- owner$.context
    clearWorkflowStepProvenance(owner)
    if (is.null(context$snapshot$products[[run$reference$entity_id]])) rootWorkflowProduct(owner, run$reference, list(source))
    snapshot <- context$snapshot
    snapshot$exports$trx_ruviii_export <- list(operation = "trx_ruviii_export", source = source, reference = run$reference,
        parameters = list(ruviii_m_source = run$value$m_source))
    publishWorkflowStateHead(context, snapshot)
    invisible(run)
}))

k274Plant <- function(owner, m = "group") {
    run <- k274Operate(owner, m)
    k274Accept(owner, withWorkflowStateScope(owner, function() getWorkflowStateReference(owner)), run)
    invisible(run)
}

k274Export <- function(owner, m = "group") kit172Quiet(runTrxReferenceRuviiiExport(owner, m))

k274Snapshot <- function(owner) owner$.context$snapshot

#' What an in-process export did: the condition it raised (NULL when none) and its message ("<no error>" when none).
k274Attempt <- function(owner, m = "group") {
    error <- tryCatch({
        k274Export(owner, m)
        NULL
    }, error = function(e) e)
    list(error = error, message = if (inherits(error, "condition")) conditionMessage(error) else "<no error>")
}

k274Record <- function(owner) {
    path <- tempfile(fileext = ".json")
    on.exit(unlink(path))
    kit172Quiet(writeWorkflowAnalysisRecord(owner, path))
    jsonlite::read_json(path, simplifyVector = FALSE)
}

#' The record an older build wrote for the planted state: the builder's record with the export's analysis injected as the
#' builder before KIT-274 listed it, through the same JSON as Save's.
k274OldRecord <- function(owner) {
    snapshot <- k274Snapshot(owner)
    analysis <- workflowRecordExportAnalysis(owner$.context, snapshot, snapshot$exports$trx_ruviii_export, readWorkflowShippedDefaults())
    record <- kit172Quiet(buildWorkflowAnalysisRecord(owner))
    record$analyses$trx_ruviii_export <- analysis
    path <- tempfile(fileext = ".json")
    on.exit(unlink(path))
    jsonlite::write_json(record, path, auto_unbox = TRUE, digits = I(17), null = "null", na = "null", pretty = TRUE)
    jsonlite::read_json(path, simplifyVector = FALSE)
}

k274Listed <- function(record) "trx_ruviii_export" %in% names(record$analyses)

k274LastRuvs <- function(record) {
    steps <- Filter(function(step) identical(step$operation, "trx_ruvs"), record$steps)
    if (length(steps)) steps[[length(steps)]]$outcomes else list()
}

k274Exprs <- function(label) {
    path <- system.file("reports", "transcriptomics", "transcriptomics_report.rmd", package = "MultiScholaR")
    lines <- readLines(path, warn = FALSE, encoding = "UTF-8")
    start <- grep(paste0("^```\\{r ", label, "[,}]"), lines)
    stopifnot(length(start) == 1L)
    end <- start + which(lines[-seq_len(start)] == "```")[[1L]]
    chunk <- tempfile(fileext = ".Rmd")
    writeLines(lines[start:end], chunk, useBytes = TRUE)
    code <- tempfile(fileext = ".R")
    suppressMessages(knitr::purl(chunk, output = code, quiet = TRUE, documentation = 0))
    parse(code, keep.source = FALSE, encoding = "UTF-8")
}

#' An environment holding what the head template's "Project Setup" chunk defines from the record: its functions and its
#' own `ruv`, `no_ruv` and `export` assignments (so the template's own guard decides what `export` is).
k274Env <- function(record, results_dir) {
    env <- new.env(parent = globalenv())
    env$report_results_summary_dir <- results_dir
    env$record <- record
    keep <- c("ruv", "no_ruv", "export")
    for (e in k274Exprs("Project Setup")) {
        if (!is.call(e) || !identical(e[[1L]], as.name("<-"))) next
        target <- if (is.name(e[[2L]])) as.character(e[[2L]]) else ""
        is_function <- is.call(e[[3L]]) && identical(e[[3L]][[1L]], as.name("function"))
        if (is_function || target %in% keep) eval(e, env)
    }
    env
}

k274Chunk <- function(env, label) {
    result <- NULL
    out <- utils::capture.output(for (e in k274Exprs(label)) result <- eval(e, env))
    list(out = paste(out, collapse = "\n"), value = result)
}

k274Figure <- function(results_dir) {
    file <- file.path(results_dir, "QC_figures", "ruviii_before_and_after.png")
    dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
    writeBin(charToRaw("kit274"), file)
    file
}

#' Copy of `k264Run`: trx's "Run Normalization Pipeline" in process, on one `nd` shared by a project's runs (R1b, VAL-54).
k274Run <- function(project, inputs, nd) {
    testthat::local_mocked_bindings(workflowOmicsManagedActions = function(omic_type) NULL, .package = "MultiScholaR")
    notify <- function(ui, ..., type = "default", duration = 5) invisible(NULL)
    set_status <- function(tab, value) {
        status <- project$workflow$tab_status
        status[[tab]] <- value
        project$workflow$tab_status <- status
    }
    withr::local_seed(20261124L)
    kit172Quiet(utils::capture.output(shiny::isolate(runTrxNormalisation(list(managed = runTrxNormManaged,
        in_process = runTrxNormInProcess), transcriptomicsDeclarationBinding, project$workflow$state_manager,
        trxNormGridParameters(inputs), nd, project$workflow$config_list, project$paths,
        omicNormAddLog(nd), notify, set_status))))
    invisible(NULL)
}

#' The module-CI project (the MCI-030 corpus) with the RUV design option on.
k274Project <- function(env = parent.frame()) {
    skipWithoutTrxDeclReference()
    withr::local_options(MultiScholaR.ruv_design_panel = TRUE, .local_envir = env)
    module_ci_trx_project(module_ci_trx_corpus(), env = env)
}

#' The KIT-040 corpus trx_edge_12 on a session owner at the QC state (`integrationTrxOwner`).
k274EdgeOwner <- function(env = parent.frame()) {
    skipWithoutTrxDeclReference()
    reference <- trxDeclReference("trx_edge_12", env = env)
    fixture <- trxDeclarationsOwner(withr::local_tempdir(.local_envir = env), reference$dir)
    withr::defer(closeDIAStateAdapter(fixture$owner), envir = env)
    kit172Quiet({
        runTrxReferenceExpressionFilter(fixture$owner, list())
        commitWorkflowOmicsStoreStep(fixture$owner, "trx_qc_finalize", "trx_qc_complete", list(), "QC processing complete")
    })
    fixture
}

#' Copy of `integrationTrxSave`: Save on the summary's own function, the record and the study-parameter writer replaced.
k274Save <- function(fixture, root, label = "kit274") {
    paths <- list(source_dir = file.path(root, "scripts", "transcriptomics"), integration_dir = file.path(root, "integration"),
        base_dir = root)
    dir.create(paths$source_dir, recursive = TRUE, showWarnings = FALSE)
    saved <- NULL
    failure <- NULL
    kit172Quiet(utils::capture.output(saveWorkflowSummaryParameters(omicGuiConfig("transcriptomics")$summary, fixture$workflow,
        projectDirs = list(transcriptomics = paths), omicType = "transcriptomics", label = label, description = "d",
        notify = function(message, ...) NULL, done = function(value) saved <<- value,
        failed = function(error) failure <<- error, createWorkflowArgsFn = function(...) "study_parameters.txt",
        writeAnalysisRecordFn = function(...) NULL, catFn = function(...) NULL)))
    list(saved = saved, failure = failure, dir = file.path(paths$integration_dir, "transcriptomics"),
        tsv = file.path(paths$integration_dir, "transcriptomics_normalized_data.tsv"))
}

# ------------------------------------------------------------------------------------------------
# T1, T2: in process.
# ------------------------------------------------------------------------------------------------

test_that("KIT-274 T1: the design none, then the export in process: refused, nothing accepted", {
    skip_if_not_installed("edgeR")
    project <- k274Project()
    owner <- project$owner
    k274NoneRun(owner)
    expect_identical(owner$current_state, "trx_ruvs", info = "run precondition")
    history <- unlist(owner$getHistory())
    attempt <- k274Attempt(owner)
    k274Row("T1", "class:", paste(class(attempt$error), collapse = "/"), "message:", attempt$message,
        "export in snapshot:", !is.null(k274Snapshot(owner)$exports$trx_ruviii_export))
    expect_s3_class(attempt$error, K274Class)
    expect_identical(attempt$message, K274Refusal(K274Reason))
    expect_null(k274Snapshot(owner)$exports$trx_ruviii_export)
    expect_identical(unlist(owner$getHistory()), history)
    expect_identical(owner$current_state, "trx_ruvs")
})

test_that("KIT-274 T2: a declared design whose estimated k is forced to 0, then the export in process: the same refusal", {
    skip_if_not_installed("edgeR")
    project <- k274Project()
    owner <- project$owner
    real <- estimateRuvFactorCount
    testthat::local_mocked_bindings(estimateRuvFactorCount = function(...) {
        e <- real(...)
        e$k <- 0L
        e
    }, .package = "MultiScholaR")
    kit172Quiet(runTrxReferenceRuvs(owner, NULL, list(),
        design_parameters = trxRuvDesignParameters("groups_all_features", ruv_k_rule = "num_sv_be")))
    expect_identical(owner$current_state, "trx_ruvs", info = "run precondition")
    outcomes <- k274LastRuvs(k274Record(owner))
    expect_identical(outcomes$ruv_applied, FALSE, info = "precondition: the k rule found no factor")
    expect_equal(outcomes$k, 0, info = "precondition: k is 0")
    reason <- outcomes$ruv_not_applied_reason
    expect_true(is.character(reason) && length(reason) == 1L && nzchar(reason), info = "precondition: a reason is recorded")
    attempt <- k274Attempt(owner)
    k274Row("T2", "class:", paste(class(attempt$error), collapse = "/"), "message:", attempt$message)
    expect_s3_class(attempt$error, K274Class)
    expect_identical(attempt$message, K274Refusal(reason))
    expect_null(k274Snapshot(owner)$exports$trx_ruviii_export)
})

# ------------------------------------------------------------------------------------------------
# T3, T4: the mounted tabs, in process and managed.
# ------------------------------------------------------------------------------------------------

k274Notes <- function(out) Filter(function(note) startsWith(note$message, "RUVIII export error:"), out$notes)

k274Files <- function(fixture) list.files(c(fixture$paths$publication_graphs_dir, fixture$paths$results_summary_dir),
    pattern = "^ruviii_before_and_after", recursive = TRUE)

test_that("KIT-274 T3: the design none through the mounted tabs, in process: the refusal is shown, nothing is run or drawn", {
    fixture <- k274Project()
    out <- module_ci_trx_tabs_journey(fixture, managed = FALSE, norm = K274Design, ruviii = "group")
    expect_identical(out$errors, character(), info = "precondition: no step escaped")
    notes <- k274Notes(out)
    k274Row("T3", "notes:", length(notes), "status:", out$ruviii_status, "figure files:", length(k274Files(fixture)),
        "export in snapshot:", !is.null(k274Snapshot(fixture$owner)$exports$trx_ruviii_export))
    expect_length(notes, 1L)
    if (length(notes)) {
        expect_identical(notes[[1L]]$type, "error")
        expect_identical(notes[[1L]]$message, paste("RUVIII export error:", K274Refusal(K274Reason), paste0("(", K274Class, ")")))
    }
    expect_identical(out$ruviii_status, "Not run yet.")
    expect_length(k274Files(fixture), 0L)
    expect_null(k274Snapshot(fixture$owner)$exports$trx_ruviii_export)
})

test_that("KIT-274 T4: the same, managed (the journey's runner on the chain mocks): one notification, nothing escaped, nothing accepted", {
    fixture <- k274Project()
    probe <- new.env()
    omicsManagedChainMocks(fixture, probe)
    # The press is made in the summary hook (after the DA), so that the products are read before and after it.
    press <- function(session, output, out) {
        out$products_before <- names(fixture$owner$.context$snapshot$products)
        session$setInputs(`summary-ruviii_m_source` = "group")
        session$setInputs(`summary-run_ruviii_export` = 1L)
        session$flushReact()
        out$status <- output$`summary-ruviii_status`
        out$products_after <- names(fixture$owner$.context$snapshot$products)
    }
    out <- module_ci_trx_tabs_journey(fixture, managed = TRUE, norm = K274Design, summary = press)
    expect_identical(out$errors, character(), info = "nothing escaped the observer")
    notes <- k274Notes(out)
    k274Row("T4", "notes:", length(notes), "status:", out$status, "products added:",
        length(setdiff(out$products_after, out$products_before)), "export in snapshot:",
        !is.null(k274Snapshot(fixture$owner)$exports$trx_ruviii_export))
    expect_length(notes, 1L)
    if (length(notes)) {
        expect_identical(notes[[1L]]$type, "error")
        expect_identical(notes[[1L]]$message, paste("RUVIII export error:", K274Refusal(K274Reason), paste0("(", K274Class, ")")))
    }
    expect_identical(out$status, "Not run yet.")
    expect_null(k274Snapshot(fixture$owner)$exports$trx_ruviii_export)
    expect_identical(out$products_after, out$products_before)
})

# ------------------------------------------------------------------------------------------------
# T5, T6, T7: an export made by an earlier build (planted): the record, the report, the replay.
# ------------------------------------------------------------------------------------------------

test_that("KIT-274 T5: a planted older export is not listed in the analysis record", {
    skip_if_not_installed("edgeR")
    project <- k274Project()
    owner <- project$owner
    k274NoneRun(owner)
    k274Plant(owner)
    expect_false(is.null(k274Snapshot(owner)$exports$trx_ruviii_export), info = "precondition: the older export is in the state")
    record <- k274Record(owner)
    k274Row("T5", "ruv_applied:", format(k274LastRuvs(record)$ruv_applied), "listed:", k274Listed(record), "k:",
        format(record$analyses$trx_ruviii_export$outcomes$k %||% NA))
    expect_identical(k274LastRuvs(record)$ruv_applied, FALSE)
    expect_false(k274Listed(record))
})

test_that("KIT-274 T6: the report of the record an older build wrote says that no RUVIII export was recorded", {
    skip_if_not_installed("edgeR")
    project <- k274Project()
    owner <- project$owner
    k274NoneRun(owner)
    k274Plant(owner)
    record <- k274OldRecord(owner)
    expect_true(k274Listed(record), info = "precondition: the older record lists the export")
    expect_identical(k274LastRuvs(record)$ruv_applied, FALSE, info = "precondition: its last trx_ruvs step applied no RUV")
    results <- withr::local_tempdir()
    figure <- k274Figure(results)
    env <- k274Env(record, results)
    text <- k274Chunk(env, "ruviii-text")$out
    shown <- k274Chunk(env, "ruviii-figure")
    k274Row("T6", "text:", text, "figure file:", file.exists(figure), "figure output:", shown$out, "figure class:",
        paste(class(shown$value), collapse = ","))
    expect_identical(text, K274None)
    expect_false(grepl("was applied", text, fixed = TRUE))
    expect_identical(shown$out, K274Missing)
    expect_false(inherits(shown$value, "knit_image_paths"))
})

test_that("KIT-274 T7: the replay of a record that lists such an export skips it, says so and warns", {
    fixture <- k274EdgeOwner()
    owner <- fixture$owner
    k274NoneRun(owner)
    kit172Quiet(runTrxReferenceDaAnalysis(owner, getWorkflowStateReference(owner), list(da_q_val_thresh = 0.05, treat_lfc_cutoff = 1),
        fixture$workflow$contrasts_tbl))
    k274Plant(owner)
    record <- k274OldRecord(owner)
    expect_true(k274Listed(record), info = "precondition: the record lists the export")
    harness <- new.env(parent = globalenv())
    sys.source(file.path(trxDeclRoot(), "tools", "profiling", "workflow_kit.R"), harness)
    api <- function(name) get(name, asNamespace("MultiScholaR"))
    copy <- withr::local_tempdir()
    stopifnot(system2("cp", c("-a", file.path(fixture$source_dir, "."), copy)) == 0L)
    store_root <- dirname(dirname(list.files(fixture$source_dir, "^store[.]json$", recursive = TRUE, all.files = TRUE, full.names = TRUE)[[1L]]))
    prepared <- harness$kitPrepareReplayProject(file.path(copy, substring(store_root, nchar(fixture$source_dir) + 2L)),
        c("initial", "trx_raw_data_s4"), api)
    inputs <- c(trx_counts_1 = file.path(fixture$source_dir, "import_sources", "counts.tsv"),
        design_file = file.path(fixture$source_dir, "design_matrix.tab"))
    classes <- character()
    replayed <- NULL
    utils::capture.output(replayed <- withCallingHandlers(suppressMessages(replayWorkflowAnalysisRecord(record, prepared, inputs,
        fixture$paths, withr::local_tempdir())), warning = function(w) {
        classes <<- c(classes, class(w))
        invokeRestart("muffleWarning")
    }))
    identity <- readSavedStoreIdentity(file.path(prepared, ".workflow-entities", "store.json"))
    again <- newWorkflowStateOwner(newWorkflowEntityStore(prepared, identity$project_id, identity$workflow_id), reopen = TRUE)
    withr::defer(closeWorkflowStateOwner(again))
    k274Row("T7", "statement:", replayed$export_statement %||% "<none>", "warning classes:", paste(unique(classes), collapse = ","),
        "export in replayed snapshot:", !is.null(again$.context$snapshot$exports$trx_ruviii_export))
    expect_identical(replayed$export_statement, K274Statement)
    expect_true("workflow_replay_export_skipped" %in% classes)
    expect_null(again$.context$snapshot$exports$trx_ruviii_export)
})

# ------------------------------------------------------------------------------------------------
# T8, T9: the integration input, with and without RUV (guards).
# ------------------------------------------------------------------------------------------------

test_that("KIT-274 T8 (guard): the planted older export, then Save: the integration input is the TMM logCPM, normalised", {
    fixture <- k274EdgeOwner()
    owner <- fixture$owner
    root <- withr::local_tempdir()
    ruvs <- k274NoneRun(owner)
    k274Plant(owner)
    expect_false(is.null(k274Snapshot(owner)$exports$trx_ruviii_export), info = "precondition: the older export is in the state")
    run <- k274Save(fixture, root)
    expect_null(run$failure)
    expect_identical(run$saved$integrationFile, "transcriptomics_normalized_data.tsv")
    input <- readIntegrationInputWithMatrices(run$dir)
    manifest <- input$manifest
    k274Row("T8", "matrix_role:", manifest$matrix_role, "corrected:", input$views$transcriptomics$corrected, "product:",
        manifest$source_state$product$type, "digest of the written file:", digest::digest(file = run$tsv, algo = "sha256"))
    expect_identical(manifest$matrix_role, "normalised")
    expect_false(input$views$transcriptomics$corrected)
    expect_identical(manifest$source_state$product[c("type", "entity_id")],
        list(type = "trx.tmm_logcpm", entity_id = ruvs$logcpm$reference$entity_id))
    log_cpm <- withWorkflowStateScope(owner, function() trxSavedMatrix(owner, ruvs$logcpm$value$log_cpm))
    tsv <- utils::read.delim(run$tsv, check.names = FALSE)
    expect_identical(tsv$Run, colnames(log_cpm))
    expect_equal(unname(as.matrix(tsv[-1L])), unname(t(log_cpm)), tolerance = 1e-6)
})

test_that("KIT-274 T9 (guard): after a legacy run the export is accepted, listed with the run's k, described, and is the integration input", {
    fixture <- k274EdgeOwner()
    owner <- fixture$owner
    root <- withr::local_tempdir()
    ruvs <- kit172Quiet(runTrxReferenceRuvs(owner, trxRuvGridParameters(), list()))
    export <- k274Export(owner)
    expect_false(is.null(k274Snapshot(owner)$exports$trx_ruviii_export))
    record <- k274Record(owner)
    k <- as.integer(ruvs$value$k)
    expect_true(k274Listed(record))
    expect_equal(record$analyses$trx_ruviii_export$outcomes$k, k)
    expect_true(k >= 1L, info = "precondition: the legacy grid's k is at least 1")
    results <- withr::local_tempdir()
    figure <- k274Figure(results)
    env <- k274Env(record, results)
    text <- k274Chunk(env, "ruviii-text")$out
    shown <- k274Chunk(env, "ruviii-figure")
    expect_identical(text, paste0("RUVIII (ruv::RUVIII) was applied to the TMM logCPM, with the RUVs inferred negative controls and k = ", k,
        ", for the integration export only, with the replicate structure M taken from the group. The differential expression ",
        "results do not use it."))
    expect_true(inherits(shown$value, "knit_image_paths"))
    run <- k274Save(fixture, root)
    expect_null(run$failure)
    input <- readIntegrationInputWithMatrices(run$dir)
    k274Row("T9", "listed:", k274Listed(record), "k:", k, "corrected:", input$views$transcriptomics$corrected, "product:",
        input$manifest$source_state$product$type)
    expect_true(input$views$transcriptomics$corrected)
    expect_identical(input$manifest$source_state$product[c("type", "entity_id")],
        list(type = "trx.ruviii", entity_id = export$reference$entity_id))
})

# ------------------------------------------------------------------------------------------------
# T10, T11: guards.
# ------------------------------------------------------------------------------------------------

test_that("KIT-274 T10 (guard): KIT-264's routes still list no earlier export (run k = 2, export, run k = 1; and with the design none)", {
    skip_if_not_installed("edgeR")
    ruvs_entity <- function(owner) k274Snapshot(owner)$states[["trx_ruvs"]]$entity_id
    one <- k274Project()
    nd_one <- omicNormReactiveState()
    k274Run(one, K274Legacy2, nd_one)
    expect_identical(one$owner$current_state, "trx_ruvs", info = "route A run 1")
    first_one <- ruvs_entity(one$owner)
    k274Export(one$owner)
    expect_true(k274Listed(k274Record(one$owner)), info = "the true export is listed")
    k274Run(one, K274Legacy1, nd_one)
    expect_false(identical(first_one, ruvs_entity(one$owner)), info = "precondition: route A's second run committed a new trx_ruvs entity")
    two <- k274Project()
    nd_two <- omicNormReactiveState()
    k274Run(two, K274Legacy2, nd_two)
    first_two <- ruvs_entity(two$owner)
    k274Export(two$owner)
    k274Run(two, K274Design, nd_two)
    expect_identical(two$owner$current_state, "trx_ruvs", info = "route B run 2")
    expect_false(identical(first_two, ruvs_entity(two$owner)), info = "precondition: route B's second run committed a new trx_ruvs entity")
    a <- k274Listed(k274Record(one$owner))
    b <- k274Listed(k274Record(two$owner))
    k274Row("T10", "route A listed:", a, "route B listed:", b)
    expect_false(a)
    expect_false(b)
})

# ------------------------------------------------------------------------------------------------
# T12: the call sites.
# ------------------------------------------------------------------------------------------------

test_that("KIT-274 T12: the parameters of the managed request and the acceptance each refuse a head that applied no RUV", {
    skip_if_not_installed("edgeR")
    project <- k274Project()
    owner <- project$owner
    k274NoneRun(owner)
    source <- getWorkflowStateReference(owner)
    parameters <- tryCatch(trxRuviiiExportParameters(owner, source, "group"), error = function(e) e)
    k274Row("T12", "parameters:", if (inherits(parameters, "condition")) paste(class(parameters), collapse = "/") else "<no error>")
    expect_s3_class(parameters, K274Class)
    run <- k274Operate(owner)
    accepted <- tryCatch(acceptTrxRuviiiExport(owner, source, run), error = function(e) e)
    k274Row("T12", "acceptance:", if (inherits(accepted, "condition")) paste(class(accepted), collapse = "/") else "<no error>")
    expect_s3_class(accepted, K274Class)
    expect_null(k274Snapshot(owner)$exports$trx_ruviii_export)
})
