# KIT-110 W2: the transcriptomics normalisation core (R/utils_workflow_trx_norm.R) and the generic
# in-process publication runner (`runWorkflowOmicsPublicationInProcess`, R/utils_workflow_omics_actions.R),
# moved from R/mod_trx_norm.R (D4, D5). Written before the move; the first two blocks fail on the
# start commit (d5244ffe), the third errors there (the handlers' new signatures)
# (docs/qa/workflow-execution/kit/kit-110/w2/test-first.txt). The handlers are called outside a
# session with recording `notify` and `setStatus`, on MCI-030's fixture; their log, notifications,
# flags and status are W1 (b-trx)'s in-process captures (test-trx-norm-module.R).

trxNormCoreMoved <- c("trxNormGridParameters", "trxNormRuvQcResult", "trxNormRuvsStep", "runTrxNormPublicationsInProcess",
    "runTrxNormManaged", "runTrxNormInProcess", "runTrxNormalisation", "resetTrxNormalisation")

test_that("the moved functions are the namespace's, defined in utils_workflow_trx_norm.R, and call no Shiny UI function", {
    ns <- asNamespace("MultiScholaR")
    files <- c(stats::setNames(rep("utils_workflow_trx_norm.R", length(trxNormCoreMoved)), trxNormCoreMoved),
        runWorkflowOmicsPublicationInProcess = "utils_workflow_omics_actions.R")
    for (name in names(files)) {
        fn <- get0(name, envir = ns, inherits = FALSE)
        expect_true(is.function(fn), info = name)
        if (!is.function(fn)) next
        expect_identical(environment(fn), ns, info = name)
        # Source references are kept under load_all only.
        if (!is.null(utils::getSrcref(fn))) {
            expect_identical(basename(utils::getSrcFilename(fn, full.names = TRUE)), files[[name]], info = name)
        }
        expect_false(any(c("showNotification", "req", "withProgress", "incProgress") %in% all.names(body(fn))), info = name)
    }
    # Made generic (D5): the trx-only runner is gone.
    expect_false(exists("runTrxPublicationInProcess", envir = ns, inherits = FALSE))
    expect_identical(names(formals(runWorkflowOmicsPublicationInProcess)),
        c("binding", "manager", "operation", "parameters", "experiment_paths"))
    expect_identical(names(formals(runTrxNormManaged)),
        c("binding", "runner", "manager", "settings", "normData", "config", "addLog", "finish", "failed"))
    expect_identical(names(formals(runTrxNormInProcess)),
        c("binding", "manager", "settings", "normData", "config", "addLog", "experimentPaths", "finish", "failed"))
    expect_identical(names(formals(runTrxNormalisation)),
        c("run", "binding", "manager", "settings", "normData", "config", "experimentPaths", "addLog", "notify", "setStatus"))
    expect_identical(names(formals(resetTrxNormalisation)), c("manager", "normData", "config", "addLog", "notify", "setStatus"))
})

test_that("runWorkflowOmicsPublicationInProcess resolves through binding$publications, and refuses an unknown operation by class", {
    root <- withr::local_tempdir()
    paths <- list(gene_qc_dir = file.path(root, "gene_qc"), lipid_qc_dir = file.path(root, "lipid_qc"))
    seen <- new.env(parent = emptyenv())
    seen$calls <- list()
    # Each handler records its request and stages one file under its omic's QC role.
    handler <- function(role) {
        force(role)
        function(manager, request, context) {
            seen$calls[[length(seen$calls) + 1L]] <- list(operation = request$operation, parameters = request$parameters,
                state = request$input_references$state, context = context, scoped = isTRUE(seen$scoped))
            target <- file.path(request$directory, "omics-publication", role)
            dir.create(target, recursive = TRUE)
            writeLines(request$operation, file.path(target, paste0(request$operation, ".txt")))
            list(display = list(saved = TRUE, operation = request$operation))
        }
    }
    testthat::local_mocked_bindings(
        computeWorkflowTrxQcFiguresAction = handler("gene_qc_dir"),
        computeWorkflowOmicsQcFiguresAction = handler("lipid_qc_dir"),
        getWorkflowStateReference = function(manager) list(entity_id = paste0("head-of-", manager)),
        withWorkflowStateScope = function(manager, fn) {
            seen$scoped <- TRUE
            on.exit(seen$scoped <- FALSE)
            fn()
        },
        .package = "MultiScholaR")
    display <- runWorkflowOmicsPublicationInProcess(transcriptomicsDeclarationBinding, "trx-owner", "trx_qc_figures",
        list(stage = "raw"), paths)
    expect_identical(display, list(saved = TRUE, operation = "trx_qc_figures"))
    expect_identical(seen$calls[[1L]], list(operation = "trx_qc_figures", parameters = list(stage = "raw"),
        state = list(entity_id = "head-of-trx-owner"), context = NULL, scoped = TRUE))
    expect_identical(list.files(paths$gene_qc_dir), "trx_qc_figures.txt")
    display <- runWorkflowOmicsPublicationInProcess(lipidomicsDeclarationBinding, "lipid-owner", "lipid_qc_figures",
        list(stage = "post_norm"), paths)
    expect_identical(display$operation, "lipid_qc_figures")
    expect_identical(seen$calls[[2L]]$state, list(entity_id = "head-of-lipid-owner"))
    expect_identical(list.files(paths$lipid_qc_dir), "lipid_qc_figures.txt")
    # An operation the binding does not publish, including another omic's, is refused before anything runs.
    expect_error(runWorkflowOmicsPublicationInProcess(transcriptomicsDeclarationBinding, "trx-owner", "trx_unknown", list(), paths),
        "unknown transcriptomics publication", class = "workflow_operation_parameters")
    expect_error(runWorkflowOmicsPublicationInProcess(transcriptomicsDeclarationBinding, "trx-owner", "lipid_qc_figures", list(), paths),
        class = "workflow_operation_parameters")
    expect_error(runWorkflowOmicsPublicationInProcess(lipidomicsDeclarationBinding, "lipid-owner", "lipid_ruv_publication", list(), paths),
        "unknown lipidomics publication", class = "workflow_operation_parameters")
    expect_error(runWorkflowOmicsPublicationInProcess(lipidomicsDeclarationBinding, "lipid-owner", "trx_qc_figures", list(), paths),
        class = "workflow_operation_parameters")
    expect_length(seen$calls, 2L)
})

#' Drive the handlers outside a session on MCI-030's fixture: `steps` are named functions of the
#' context, each followed by a capture (the new log lines, the notifications, the flags, the history
#' and the normalisation status).
trxNormCoreDrive <- function(steps) {
    fixture <- module_ci_trx_project(module_ci_trx_corpus())
    owner <- fixture$owner
    nd <- omicNormReactiveState()
    add_log <- omicNormAddLog(nd)
    seen <- new.env(parent = emptyenv())
    seen$notes <- character()
    seen$log <- 0L
    notify <- function(ui, ..., type = "default", duration = 5) {
        seen$notes <- c(seen$notes, sprintf("[%s %s] %s", type, format(duration), paste(as.character(ui), collapse = " ")))
        invisible(NULL)
    }
    set_status <- function(tab, value) {
        status <- fixture$workflow$tab_status
        status[[tab]] <- value
        fixture$workflow$tab_status <- status
    }
    run_forms <- list(managed = runTrxNormManaged, in_process = runTrxNormInProcess)
    ctx <- list(
        run = function(inputs) runTrxNormalisation(run_forms, transcriptomicsDeclarationBinding, owner, trxNormGridParameters(inputs),
            nd, fixture$workflow$config_list, fixture$paths, add_log, notify, set_status),
        reset = function() resetTrxNormalisation(owner, nd, fixture$workflow$config_list, add_log, notify, set_status))
    captures <- list()
    withr::local_seed(20261108L)
    for (step in names(steps)) {
        suppressWarnings(suppressMessages(utils::capture.output(shiny::isolate(steps[[step]](ctx)))))
        values <- shiny::isolate(shiny::reactiveValuesToList(nd))
        lines <- values$normalization_log
        state_of <- function(value) {
            if (is.null(value)) return("<NULL>")
            named <- names(Filter(function(ref) identical(ref$entity_id, value$entity_id), owner$.context$snapshot$states))
            paste(named, collapse = "|")
        }
        captures[[step]] <- list(
            log = sub("^\\[[0-9:]+\\] ", "", lines[seq_along(lines) > seen$log]),
            notes = seen$notes,
            norm = c(normalization_complete = as.character(isTRUE(values$normalization_complete)),
                ruv_complete = as.character(isTRUE(values$ruv_complete)), plot_refresh_trigger = format(values$plot_refresh_trigger),
                post_filter_obj = state_of(values$post_filter_obj), ruv_corrected_obj = state_of(values$ruv_corrected_obj)),
            history = as.character(owner$getHistory()),
            normalization = fixture$workflow$tab_status$normalization)
        seen$log <- length(lines)
        seen$notes <- character()
    }
    captures
}

# W1 (b-trx)'s in-process captures at d5244ffe (test-trx-norm-module.R), the fields the handlers write.
trxNormCoreHistory <- c("initial", "trx_raw_data_s4", "trx_expression_filtered", "trx_qc_complete")
trxNormCoreRunLog <- function(mode) {
    c(sprintf("Running TMM, the RUV grid (%s) and RUVs", mode), "RUVs applied: W_1..W_2 added to the design (state 'trx_ruvs')",
        "PCA and RLE (raw) generated", "PCA and RLE (tmm) generated", "PCA and RLE (ruvs) generated",
        "CCA panels, grid workbook and the normalisation comparison grid saved",
        "Normalization complete - ready for differential expression")
}
trxNormCoreCompleted <- c(normalization_complete = "TRUE", ruv_complete = "TRUE", plot_refresh_trigger = "1",
    post_filter_obj = "trx_qc_complete", ruv_corrected_obj = "trx_ruvs")

test_that("the handlers outside a session: an automatic run, then a reset, as W1 (b-trx) in process", {
    skip_if_not_installed("edgeR")
    seen <- trxNormCoreDrive(list(
        run = function(ctx) ctx$run(list(ruv_mode = "automatic")),
        reset = function(ctx) ctx$reset()))
    expect_identical(seen$run, list(log = trxNormCoreRunLog("automatic"),
        notes = "[message 5] Normalization complete! Proceeding to DE analysis.", norm = trxNormCoreCompleted,
        history = c(trxNormCoreHistory, "trx_ruvs"), normalization = "complete"))
    expect_identical(seen$reset, list(log = "Reset to pre-normalization state", notes = "[message 5] Reset to pre-normalization state",
        norm = c(normalization_complete = "FALSE", ruv_complete = "FALSE", plot_refresh_trigger = "1",
            post_filter_obj = "trx_qc_complete", ruv_corrected_obj = "<NULL>"),
        history = c(trxNormCoreHistory, "trx_ruvs", "trx_reset"), normalization = "pending"))
})

test_that("the handlers outside a session: a manual run (k 2, 10 %), as W1 (b-trx) in process", {
    skip_if_not_installed("edgeR")
    seen <- trxNormCoreDrive(list(run = function(ctx) ctx$run(list(ruv_mode = "manual", ruv_k = 2L, ruv_percentage = 10))))
    expect_identical(seen$run, list(log = trxNormCoreRunLog("manual"),
        notes = "[message 5] Normalization complete! Proceeding to DE analysis.", norm = trxNormCoreCompleted,
        history = c(trxNormCoreHistory, "trx_ruvs"), normalization = "complete"))
})

test_that("the handlers outside a session: a reset before any run, then a manual k above the cap (refused), as W1 (b-trx) in process", {
    skip_if_not_installed("edgeR")
    seen <- trxNormCoreDrive(list(
        reset_before_run = function(ctx) ctx$reset(),
        refused = function(ctx) ctx$run(list(ruv_mode = "manual", ruv_k = 4L, ruv_percentage = 10))))
    expect_identical(seen$reset_before_run, list(log = "Reset to pre-normalization state",
        notes = "[message 5] Reset to pre-normalization state",
        norm = c(normalization_complete = "FALSE", ruv_complete = "FALSE", plot_refresh_trigger = "0",
            post_filter_obj = "<NULL>", ruv_corrected_obj = "<NULL>"),
        history = trxNormCoreHistory, normalization = "pending"))
    refusal <- "Normalization error: manual k 4 is above the cap of 3 for 12 samples (trx_ruv_k_above_cap)"
    expect_identical(seen$refused, list(log = c("Running TMM, the RUV grid (manual) and RUVs", refusal),
        notes = paste("[error 15]", refusal),
        norm = c(normalization_complete = "FALSE", ruv_complete = "FALSE", plot_refresh_trigger = "0",
            post_filter_obj = "trx_qc_complete", ruv_corrected_obj = "<NULL>"),
        history = trxNormCoreHistory, normalization = "pending"))
})
