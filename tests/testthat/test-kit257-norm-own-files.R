# KIT-257 (VAL-250): each normalisation-stage run first removes its own earlier files that this run would not write
# again, by exact name and in the omic's own folders, after the run's own refusals and before anything is computed,
# saved or published; a file that cannot be removed refuses the run. Written before the fix: on the tests-only commit
# E1 to E9, G3 and G5 fail at the values the ticket states (section 10); G1 and G2 pass on both commits. The functions
# the fix adds are reached only through `get0`, so the tests-only commit runs the same tests and fails at a value.
# Lines starting `KIT-257` are the evidence of docs/qa/workflow-execution/kit/kit-257/.

# ------------------------------------------------------------------------------------------------
# The harness.
# ------------------------------------------------------------------------------------------------

k257Fn <- function(name) get0(name, envir = asNamespace("MultiScholaR"), inherits = FALSE)

#' Sentinel bytes at each path (an earlier run's file, or a user's).
k257Plant <- function(files, bytes = "earlier") {
    for (file in files) {
        dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
        writeBin(charToRaw(bytes), file)
    }
    invisible(files)
}

#' Whether a regular file holds exactly `bytes`.
k257Holds <- function(file, bytes = "earlier") {
    file.exists(file) && !dir.exists(file) && identical(readBin(file, "raw", 4096L), charToRaw(bytes))
}

k257Present <- function(files) stats::setNames(file.exists(files), basename(files))

k257Log <- function(case, ...) cat(sprintf("KIT-257 %s %s\n", case, paste(..., collapse = " ")))

#' A runner that records and runs nothing (B04's form): `step$active()` is NULL while idle; a press while `busy` is
#' queued and `start()` starts it; actions given as a function are resolved when they start (as the omics runner,
#' `R/utils_workflow_omics_actions.R:43`). `present` is which of `files` existed when the first actions started.
k257Recorder <- function(files, busy = FALSE) {
    rec <- new.env(parent = emptyenv())
    rec$present <- NULL
    rec$dispatched <- 0L
    rec$queued <- list()
    rec$busy <- busy
    execute <- function(actions) {
        resolved <- if (is.function(actions)) actions() else actions
        if (is.null(rec$present)) rec$present <- k257Present(files)
        rec$dispatched <- rec$dispatched + length(resolved)
        invisible(NULL)
    }
    runner <- new.env(parent = emptyenv())
    runner$step <- list(active = function() if (isTRUE(rec$busy)) "a running action")
    runner$retrying <- function() FALSE
    runner$progress <- function(...) list(set = function(...) invisible(NULL), inc = function(...) invisible(NULL),
        close = function(...) invisible(NULL))
    runner$run <- function(actions, done = NULL, failed = NULL, dropped = NULL) {
        if (isTRUE(rec$busy)) {
            rec$queued[[length(rec$queued) + 1L]] <- list(actions = actions, failed = failed)
            return(invisible(NULL))
        }
        execute(actions)
    }
    rec$start <- function() {
        rec$busy <- FALSE
        queued <- rec$queued
        rec$queued <- list()
        for (waiting in queued) tryCatch(execute(waiting$actions), error = function(error) {
            if (is.function(waiting$failed)) waiting$failed(error) else stop(error)
        })
        invisible(NULL)
    }
    rec$runner <- runner
    rec
}

#' Copy (`copyToResultsSummary`) on `paths`: its `COPY:` lines, its failures and what reached the summary folders.
k257Copy <- function(paths, omic) {
    for (path in unlist(paths)) dir.create(path, recursive = TRUE, showWarnings = FALSE)
    testthat::local_mocked_bindings(getProjectPaths = function(...) paths, .package = "MultiScholaR")
    failures <- NULL
    console <- utils::capture.output(failures <- suppressWarnings(suppressMessages(copyToResultsSummary(omic, "k257",
        contrasts_tbl = data.frame(contrasts = "groupB-groupA"), design_matrix = data.frame(Run = "S1", group = "A"),
        force = TRUE))))
    summary <- paths$results_summary_dir
    list(console = grep("^COPY:", console, value = TRUE),
        failed = vapply(failures, function(failure) basename(as.character(failure$source %||% "")), ""),
        tables = sort(list.files(file.path(summary, "Publication_tables"))),
        qc = sort(list.files(file.path(summary, "QC_figures"))))
}

#' Copy's paths beside an omic's own: `feature_qc_dir` is the omic's QC folder; the rest under `root`.
k257CopyPaths <- function(paths, root, feature_qc_dir) {
    defaults <- list(results_dir = file.path(root, "copy", "results"),
        results_summary_dir = file.path(root, "copy", "results_summary"), publication_graphs_dir = file.path(root, "copy", "pub"),
        time_dir = file.path(root, "copy", "time"), qc_dir = file.path(root, "copy", "qc"),
        da_output_dir = file.path(root, "copy", "da"), pathway_dir = file.path(root, "copy", "pathway"),
        source_dir = file.path(root, "copy", "source"))
    for (name in names(defaults)) if (is.null(paths[[name]])) paths[[name]] <- defaults[[name]]
    paths$feature_qc_dir <- feature_qc_dir
    paths
}

# --- Proteomics: the native DIA owner of helper-prot-norm-core.R (KIT-206's clicks, copied) ---

k257ProtSettings <- function(mode = "skip", ...) utils::modifyList(list(norm_method = "quantile", ruv_mode = mode,
    auto_percentage_min = NULL, auto_percentage_max = NULL, separation_metric = NULL, k_penalty_weight = NULL,
    adaptive_k_penalty = NULL, ruv_percentage = 50, ruv_k = 1, color_variable = "group", shape_variable = NULL), list(...))

k257ProtAesthetics <- function() list(color_var = "group", shape_var = NULL)

k257ProtHead <- function(fixture) getDIAStateReference(fixture$manager)

k257ProtRuvFigures <- c(sprintf("ruv_corrected_%s.png", c("pca", "rle", "density", "correlation")), "ruv_corrected_cancor.png")

k257ProtTables <- c("ruv_normalised_results_cln_with_replicates.tsv", "normalised_results_cln_with_replicates.tsv",
    "ruv_normalised_results_cln_with_replicates.RDS", "normalised_results_cln_with_replicates.RDS")

#' The tab's in-process normalisation click; `ruv_qc` is the tab's RUV figure generator (none by default).
k257ProtClick <- function(fixture, norm, mode = "skip", ruv_qc = function(reference) invisible(NULL)) {
    protNormCoreQuiet(runDIAReferenceNormalization(k257ProtSettings(mode), fixture$workflow, norm, fixture$paths,
        "proteomics", "core", function(reference) invisible(NULL), ruv_qc, function() "group", pre_qc = NULL))
}

#' The correlation filter in process (Apply, or Skip with `skipped`).
k257ProtCorrelation <- function(fixture, norm, threshold = 0.5, skipped = FALSE) {
    final <- protNormCoreQuiet(runDIAReferenceCorrelation(fixture$workflow, norm, fixture$paths, threshold, skipped,
        "proteomics", "core", if (skipped) function() NULL else function() "group", k257ProtAesthetics))
    norm$correlation_filtering_complete <- TRUE
    final
}

#' "Run Normalization & RUV": the tab's entry (`runProtNormNormalizationWorkflow`), in process (`runner` NULL) or
#' managed, its notices and errors recorded.
k257ProtTab <- function(fixture, norm, mode = "skip", runner = NULL, settings = k257ProtSettings(mode)) {
    testthat::local_mocked_bindings(workflowOmicsManagedActions = function(omic_type) runner, .package = "MultiScholaR")
    seen <- new.env(parent = emptyenv())
    seen$notes <- character()
    seen$errors <- list()
    note <- function(ui, type = "default", duration = 5, ...) {
        seen$notes <- c(seen$notes, sprintf("[%s] %s", type, paste(as.character(ui), collapse = " ")))
        invisible(NULL)
    }
    protNormCoreQuiet(shiny::isolate(runProtNormNormalizationWorkflow(input = settings,
        workflowData = fixture$workflow, normData = norm, experimentPaths = fixture$paths, omicType = "proteomics",
        experimentLabel = "core", generatePostNormalizationQcFn = function(reference) invisible(NULL),
        generateRuvCorrectedQcFn = function(reference) invisible(NULL), getRuvGroupingVariableFn = function() "group",
        showNotificationFn = note, handleErrorFn = function(error) {
            seen$errors <- c(seen$errors, list(error))
            invisible(FALSE)
        }, messageFn = function(...) NULL)))
    as.list(seen)
}

#' The DA's re-commit of the head's arguments (KIT-206's `k206Da`); an error is returned, not raised.
k257ProtDa <- function(fixture) tryCatch({
    protNormCoreQuiet(withDIAStateScope(fixture$manager, function() commitDIADaArguments(fixture$manager,
        list(da_q_val_thresh = 0.05, treat_lfc_cutoff = 0, formula_string = "~ 0 + group"), fixture$workflow,
        data.frame(contrasts = "groupB-groupA", full_format = "B_vs_A=groupB-groupA", stringsAsFactors = FALSE),
        "Protein.Ids")))
    NULL
}, error = function(e) e)

# --- Lipidomics and metabolomics: F172 at QC complete (KIT-194's fixture, copied) ---

k257Omic <- function(class) c(lipid = "lipidomics", metab = "metabolomics")[[class]]

k257OmicFixture <- function(root, class, env = parent.frame()) {
    fixture <- kit172Owner(root, kit172Object(kit172Values(172001L)$values, class), class)
    withr::defer(try(closeDIAStateAdapter(fixture$workflow$state_manager), silent = TRUE), envir = env)
    binding <- kit172Binding(class)
    kit172Quiet({
        runWorkflowOmicsReferenceIntensityFilter(binding, fixture$owner, 1, 0.5, fixture$workflow$config_list)
        commitWorkflowOmicsStoreStep(fixture$owner, paste0(binding$prefix, "_qc_finalize"), paste0(binding$prefix, "_qc_complete"),
            fixture$workflow$config_list, "QC complete")
    })
    fixture$class <- class
    fixture$omic <- k257Omic(class)
    fixture$binding <- binding
    config <- omicGuiConfig(fixture$omic)
    config$binding$binding$qc_plots <- function(...) invisible(NULL)
    fixture$config <- config
    fixture$qc_dir <- fixture$paths[[config$binding$binding$qc_dir_role]]
    log <- new.env(parent = emptyenv())
    log$lines <- character()
    fixture$log <- log
    fixture$addLog <- function(entry) log$lines <- c(log$lines, entry)
    fixture
}

k257OmicSnapshot <- function(fixture) list(history = as.character(unlist(fixture$workflow$state_manager$getHistory())),
    head = fixture$workflow$state_manager$current_state)

#' The head's assay names, read as the tab reads them on entry (`omicNormDetectAssays`).
k257OmicAssays <- function(fixture) tryCatch(fixture$config$binding$read_metadata(fixture$workflow$state_manager)$assay_names,
    error = function(e) NULL)

k257OmicRuvFigures <- function(dir, assays) {
    safe <- gsub("[^A-Za-z0-9]", "_", tolower(assays))
    file.path(dir, sprintf("%s_ruv_corrected_%s.png", rep(safe, each = 4L), c("pca", "density", "rle", "correlation")))
}

k257OmicInput <- function(ruv = NULL) {
    list(apply_itsd = FALSE, itsd_aggregation = "median", log_offset = 0, norm_method = "quantile",
        ruv_mode = if (is.null(ruv)) "skip" else ruv$ruv_mode, auto_percentage_min = 1, auto_percentage_max = 20,
        ruv_grouping_variable = "group", separation_metric = "max_difference", k_penalty_weight = 0.5,
        adaptive_k_penalty = TRUE, ruv_k = ruv$manual_k %||% 1L, ruv_percentage = ruv$manual_percentage %||% 20,
        min_pearson_correlation_threshold = 0.5)
}

#' Shiny's notifications recorded as "[type] text"; progress and selector updates quiet.
k257Notes <- function(env = parent.frame()) {
    seen <- new.env(parent = emptyenv())
    seen$notes <- character()
    testthat::local_mocked_bindings(
        showNotification = function(ui, ..., type = "default") {
            seen$notes <- c(seen$notes, paste0("[", type, "] ", paste(as.character(ui), collapse = " ")))
            invisible(NULL)
        },
        removeNotification = function(...) invisible(NULL),
        withProgress = function(expr, ...) expr,
        incProgress = function(...) invisible(NULL),
        updateSelectInput = function(...) invisible(NULL),
        .package = "shiny", .env = env)
    seen
}

#' "Run Normalization Pipeline": the frame's `runOmicNormPipeline`, in process (`runner` NULL) or managed, on `paths`.
k257OmicRun <- function(fixture, nd, input = k257OmicInput(), paths = fixture$paths, runner = NULL) {
    testthat::local_mocked_bindings(workflowOmicsManagedActions = function(omic_type) runner, .package = "MultiScholaR")
    kit172Quiet(shiny::isolate(runOmicNormPipeline(input, fixture$workflow, paths, nd, fixture$addLog,
        function(...) list(color_var = "group", shape_var = NULL), fixture$config)))
}

# --- Transcriptomics: MCI-030's project (KIT-194's run, copied) ---

k257TrxCca <- function(paths) file.path(paths$publication_graphs_dir,
    c("canonical_correlation_panels.png", "canonical_correlation_panels.pdf", "canonical_correlation_optimization_summary.xlsx"))

k257TrxRun <- function(project, nd, inputs, seen, runner = NULL) {
    testthat::local_mocked_bindings(workflowOmicsManagedActions = function(omic_type) runner, .package = "MultiScholaR")
    notify <- function(ui, ..., type = "default", duration = 5) {
        seen$notes <- c(seen$notes, paste0("[", type, "] ", paste(as.character(ui), collapse = " ")))
        invisible(NULL)
    }
    set_status <- function(tab, value) {
        status <- project$workflow$tab_status
        status[[tab]] <- value
        project$workflow$tab_status <- status
    }
    withr::local_seed(20261124L)
    kit172Quiet(utils::capture.output(shiny::isolate(runTrxNormalisation(list(managed = runTrxNormManaged,
        in_process = runTrxNormInProcess), transcriptomicsDeclarationBinding, project$workflow$state_manager,
        trxNormGridParameters(inputs), nd, project$workflow$config_list, project$paths, omicNormAddLog(nd), notify,
        set_status))))
    invisible(NULL)
}

k257TrxLegacy <- list(ruv_mode = "manual", ruv_k = 2L, ruv_percentage = 10)
k257TrxDeclared <- list(ruv_mode = "automatic", ruv_design = "none")

k257TrxErrors <- function(seen) grep("[[]error", seen$notes, value = TRUE)

k257Unremoved <- "could not be removed"

# ------------------------------------------------------------------------------------------------
# Transcriptomics: E1, E2, G1 (with G2's user files), G3.
# ------------------------------------------------------------------------------------------------

test_that("KIT-257 E1: trx in process, a legacy grid run then a declared design run: the CCA files are gone (fails before)", {
    skip_if_not_installed("edgeR")
    withr::local_options(MultiScholaR.ruv_design_panel = TRUE)
    seen <- new.env(parent = emptyenv())
    seen$notes <- character()
    project <- module_ci_trx_project(module_ci_trx_corpus())
    nd <- omicNormReactiveState()
    cca <- k257TrxCca(project$paths)
    k257TrxRun(project, nd, k257TrxLegacy, seen)
    expect_true(all(file.exists(cca)), info = paste("run 1 precondition:", paste(seen$notes, collapse = " | ")))
    k257TrxRun(project, nd, k257TrxDeclared, seen)
    info <- paste("E1", paste(seen$notes, collapse = " | "))
    expect_identical(k257TrxErrors(seen), character(), info = info)
    expect_identical(project$owner$current_state, "trx_ruvs", info = info)
    sources <- vapply(transcriptomicsResultsSummaryFiles(project$paths), function(spec) spec$source, "")
    k257Log("E1", "after the declared run present:", paste(names(k257Present(cca))[k257Present(cca)], collapse = ","))
    expect_identical(unname(k257Present(cca)), c(FALSE, FALSE, FALSE), info = info)
    expect_false(any(file.exists(sources[basename(sources) %in% basename(cca)])), info = info)
})

test_that("KIT-257 E2: trx managed, the run's start: the CCA files are removed before the first action is dispatched (fails before)", {
    skip_if_not_installed("edgeR")
    seen <- new.env(parent = emptyenv())
    seen$notes <- character()
    project <- module_ci_trx_project(module_ci_trx_corpus())
    cca <- k257Plant(k257TrxCca(project$paths))
    rec <- k257Recorder(cca)
    k257TrxRun(project, omicNormReactiveState(), k257TrxLegacy, seen, runner = rec$runner)
    k257Log("E2", "dispatched:", rec$dispatched, "present at the first dispatch:", paste(rec$present, collapse = ","))
    expect_gt(rec$dispatched, 0L)
    expect_identical(unname(rec$present), c(FALSE, FALSE, FALSE), info = paste(seen$notes, collapse = " | "))
})

test_that("KIT-257 G1, G2 (trx): a declared run then a legacy grid run writes new CCA files; users' and other folders' files are untouched", {
    skip_if_not_installed("edgeR")
    withr::local_options(MultiScholaR.ruv_design_panel = TRUE)
    seen <- new.env(parent = emptyenv())
    seen$notes <- character()
    project <- module_ci_trx_project(module_ci_trx_corpus())
    nd <- omicNormReactiveState()
    graphs <- project$paths$publication_graphs_dir
    other <- file.path(dirname(graphs), "kit257_other_omic")
    user <- k257Plant(c(file.path(graphs, c("my_notes.png", "canonical_correlation_panels_v2.png")),
        file.path(other, basename(k257TrxCca(project$paths)))), "user")
    k257TrxRun(project, nd, k257TrxDeclared, seen)
    cca <- k257Plant(k257TrxCca(project$paths))
    k257TrxRun(project, nd, k257TrxLegacy, seen)
    info <- paste("G1 trx", paste(seen$notes, collapse = " | "))
    expect_identical(k257TrxErrors(seen), character(), info = info)
    expect_true(all(file.exists(cca)), info = info)
    expect_false(any(vapply(cca, k257Holds, NA)), info = info)
    expect_true(all(vapply(user, k257Holds, NA, bytes = "user")), info = info)
})

test_that("KIT-257 G3: trx, a directory named canonical_correlation_panels.png refuses the run with workflow_figure_unwritten (fails before)", {
    skip_if_not_installed("edgeR")
    seen <- new.env(parent = emptyenv())
    seen$notes <- character()
    project <- module_ci_trx_project(module_ci_trx_corpus())
    blocked <- k257TrxCca(project$paths)[[1L]]
    dir.create(blocked, recursive = TRUE)
    before <- as.character(unlist(project$owner$getHistory()))
    k257TrxRun(project, omicNormReactiveState(), k257TrxLegacy, seen)
    after <- as.character(unlist(project$owner$getHistory()))
    k257Log("G3", "history moved:", !identical(after, before), "notes:", paste(seen$notes, collapse = " | "))
    expect_identical(after, before)
    expect_true(any(grepl(k257Unremoved, seen$notes, fixed = TRUE) & grepl(blocked, seen$notes, fixed = TRUE)),
        info = paste(seen$notes, collapse = " | "))
    remove <- k257Fn("removeWorkflowOwnFiles")
    expect_false(is.null(remove))
    if (!is.null(remove)) {
        error <- expect_error(remove(blocked, "not removed: %s"), class = "workflow_figure_unwritten")
        expect_match(conditionMessage(error), blocked, fixed = TRUE)
        expect_true(dir.exists(blocked))
    }
})

# ------------------------------------------------------------------------------------------------
# Proteomics: E3, E4, E7, E9, G1 (with G2's user files), G5.
# ------------------------------------------------------------------------------------------------

test_that("KIT-257 E3: a correlation after a skipped-RUV normalisation removes the earlier RUV table; Copy delivers the current one (fails before)", {
    fixture <- protNormCoreOwner()
    norm <- new.env(parent = emptyenv())
    qc <- fixture$paths$protein_qc_dir
    k257ProtClick(fixture, norm, "skip")
    k257Plant(file.path(qc, k257ProtTables[c(1L, 3L)]))
    k257ProtCorrelation(fixture, norm, 0.5)
    expect_identical(fixture$manager$current_state, "correlation_filtered")
    present <- k257Present(file.path(qc, k257ProtTables))
    copy <- k257Copy(k257CopyPaths(fixture$paths, fixture$root, qc), "proteomics")
    k257Log("E3", "present:", paste(names(present)[present], collapse = ","), "copy:", paste(copy$console, collapse = " | "),
        "tables:", paste(copy$tables, collapse = ","))
    expect_identical(unname(present), c(FALSE, TRUE, FALSE, FALSE))
    expect_true("normalised_results.tsv" %in% copy$tables)
    expect_false(any(c("RUV_normalised_results.tsv", "ruv_normalised_results.RDS") %in% copy$tables))
    # Outside native sessions the .RDS is written by the S4 compatibility download only: Copy names it missing.
    expect_true("normalised_results_cln_with_replicates.RDS" %in% copy$failed)
})

test_that("KIT-257 E4: a normalisation with RUV skipped removes the earlier RUV-corrected figures (fails before)", {
    fixture <- protNormCoreOwner()
    norm <- new.env(parent = emptyenv())
    figures <- k257Plant(file.path(fixture$paths$protein_qc_dir, k257ProtRuvFigures))
    k257ProtClick(fixture, norm, "skip")
    expect_identical(fixture$manager$current_state, "normalized")
    k257Log("E4", "present:", paste(k257Present(figures), collapse = ","))
    expect_identical(unname(k257Present(figures)), rep(FALSE, length(figures)))
})

test_that("KIT-257 E7: RUV, correlation, then RUV skipped and the DA without a correlation: no earlier table, Copy delivers none (fails before)", {
    fixture <- protNormCoreOwner()
    norm <- new.env(parent = emptyenv())
    qc <- fixture$paths$protein_qc_dir
    k257ProtClick(fixture, norm, "manual")
    k257ProtCorrelation(fixture, norm, 0.5)
    first <- k257Present(file.path(qc, k257ProtTables))
    k257ProtClick(fixture, norm, "skip")
    expect_identical(fixture$manager$current_state, "normalized")
    da <- k257ProtDa(fixture)
    k257Log("E7", "after the first correlation present:", paste(names(first)[first], collapse = ","),
        "DA at", fixture$manager$current_state, "=", if (is.null(da)) "ran" else paste("refused:", conditionMessage(da)))
    present <- k257Present(file.path(qc, k257ProtTables))
    copy <- k257Copy(k257CopyPaths(fixture$paths, fixture$root, qc), "proteomics")
    k257Log("E7", "present:", paste(names(present)[present], collapse = ","), "copy:", paste(copy$console, collapse = " | "),
        "tables:", paste(copy$tables, collapse = ","))
    expect_true(first[["ruv_normalised_results_cln_with_replicates.tsv"]])
    expect_identical(unname(present), rep(FALSE, 4L))
    expect_false(any(c("RUV_normalised_results.tsv", "normalised_results.tsv") %in% copy$tables))
    expect_true(any(grepl("neither RUV nor normalized files found", copy$console, fixed = TRUE)))
})

test_that("KIT-257 E9: managed, a normalisation pressed while a correlation's action runs removes that action's table when it starts (fails before)", {
    fixture <- protNormCoreOwner()
    norm <- new.env(parent = emptyenv())
    qc <- fixture$paths$protein_qc_dir
    k257ProtClick(fixture, norm, "skip")
    table <- file.path(qc, "normalised_results_cln_with_replicates.tsv")
    rec <- k257Recorder(table, busy = TRUE)
    protNormCoreQuiet(dispatchProtNormNormalization(rec$runner, k257ProtSettings("skip"), fixture$workflow, norm,
        fixture$paths, "proteomics", "core", function(reference) invisible(NULL), function(reference) invisible(NULL),
        function() "group", function(completed) completed, function(error) stop(error)))
    expect_length(rec$queued, 1L)
    # The running correlation's publication writes its table after the click.
    k257Plant(table)
    rec$start()
    k257Log("E9", "present when the queued normalisation started:", paste(rec$present, collapse = ","))
    expect_gt(rec$dispatched, 0L)
    expect_identical(unname(rec$present), FALSE)
})

test_that("KIT-257 G1, G2 (proteomics): a skip then an RUV run writes new RUV figures; users' and other folders' files are untouched", {
    fixture <- protNormCoreOwner()
    norm <- new.env(parent = emptyenv())
    qc <- fixture$paths$protein_qc_dir
    user <- k257Plant(c(file.path(qc, c("my_notes.png", "ruv_corrected_pca_v2.png")),
        file.path(fixture$root, "results", "lipidomics", "lipid_qc", k257ProtRuvFigures)), "user")
    k257ProtClick(fixture, norm, "skip")
    figures <- k257Plant(file.path(qc, k257ProtRuvFigures[1:4]))
    draw <- function(reference) for (name in k257ProtRuvFigures[1:4]) {
        writeBin(charToRaw("new"), file.path(norm$publication_paths$protein_qc_dir, name))
    }
    k257ProtClick(fixture, norm, "manual", ruv_qc = draw)
    expect_identical(fixture$manager$current_state, "ruv_corrected")
    expect_true(all(vapply(figures, k257Holds, NA, bytes = "new")))
    expect_true(all(vapply(user, k257Holds, NA, bytes = "user")))
})

test_that("KIT-257 G5 (proteomics): a file that cannot be removed refuses the normalisation (in process and managed) and the correlation (in process and managed) with nothing committed (fails before)", {
    fixture <- protNormCoreOwner()
    norm <- new.env(parent = emptyenv())
    qc <- fixture$paths$protein_qc_dir
    k257ProtClick(fixture, norm, "skip")
    other <- k257Plant(file.path(qc, "pre_norm_pca.png"))
    blocked <- file.path(qc, "ruv_normalised_results_cln_with_replicates.tsv")
    dir.create(blocked)
    head <- k257ProtHead(fixture)
    check <- function(case, error, notes = character()) {
        k257Log("G5", case, "head:", fixture$manager$current_state, "error:",
            if (is.null(error)) "none" else paste(class(error)[[1L]], conditionMessage(error)), "notes:", paste(notes, collapse = " | "))
        expect_identical(k257ProtHead(fixture), head, info = case)
        expect_s3_class(error, "workflow_figure_unwritten")
        if (!is.null(error)) expect_match(conditionMessage(error), k257Unremoved, fixed = TRUE)
        if (!is.null(error)) expect_match(conditionMessage(error), blocked, fixed = TRUE)
        expect_true(k257Holds(other), info = case)
    }
    tab <- k257ProtTab(fixture, norm, "skip")
    check("normalisation in process", if (length(tab$errors)) tab$errors[[1L]], tab$notes)
    rec <- k257Recorder(blocked)
    tab <- k257ProtTab(fixture, norm, "skip", runner = rec$runner)
    check("normalisation managed", if (length(tab$errors)) tab$errors[[1L]], tab$notes)
    error <- tryCatch({
        k257ProtCorrelation(fixture, norm, 0.5)
        NULL
    }, error = function(e) e)
    check("correlation in process", error)
    rec <- k257Recorder(blocked)
    error <- tryCatch({
        protNormCoreQuiet(dispatchProtNormCorrelation(rec$runner, fixture$workflow, norm, fixture$paths, 0.5, "proteomics",
            "core", function() "group", k257ProtAesthetics, function(final) final, function(error) stop(error)))
        NULL
    }, error = function(e) e)
    check("correlation managed", error)
})

# ------------------------------------------------------------------------------------------------
# DIA with limpa: E8, G5.
# ------------------------------------------------------------------------------------------------

test_that("KIT-257 E8: DIA with limpa, the pass-through pressed again removes both result-table names (fails before)", {
    skip_if_not_installed("limpa")
    skip_if_not_installed("mixOmics")
    fixture <- protNormW1LimpaOwner()
    norm <- new.env(parent = emptyenv())
    qc <- fixture$paths$protein_qc_dir
    designated <- fixture$manager$current_state
    first <- k257ProtTab(fixture, norm, "covariates")
    expect_length(first$errors, 0L)
    k257ProtCorrelation(fixture, norm, 0.5)
    written <- k257Present(file.path(qc, k257ProtTables))
    k257Plant(file.path(qc, k257ProtTables)[!written])
    protNormCoreQuiet(withDIAStateScope(fixture$manager, function() selectDIAStateReference(fixture$manager, designated)))
    again <- k257ProtTab(fixture, norm, "covariates")
    present <- k257Present(file.path(qc, k257ProtTables))
    k257Log("E8", "the correlation wrote:", paste(names(written)[written], collapse = ","), "after the second press present:",
        paste(names(present)[present], collapse = ","), "errors:", length(again$errors))
    expect_length(again$errors, 0L)
    expect_identical(unname(present), rep(FALSE, 4L))
})

test_that("KIT-257 G5 (DIA with limpa): a file that cannot be removed refuses the pass-through, in process and managed (fails before)", {
    skip_if_not_installed("limpa")
    skip_if_not_installed("mixOmics")
    fixture <- protNormW1LimpaOwner()
    norm <- new.env(parent = emptyenv())
    qc <- fixture$paths$protein_qc_dir
    other <- k257Plant(file.path(qc, "pre_norm_pca.png"))
    blocked <- file.path(qc, "normalised_results_cln_with_replicates.tsv")
    dir.create(blocked, recursive = TRUE)
    head <- k257ProtHead(fixture)
    for (case in c("in process", "managed")) {
        rec <- k257Recorder(blocked)
        tab <- k257ProtTab(fixture, norm, "covariates", runner = if (identical(case, "managed")) rec$runner)
        error <- if (length(tab$errors)) tab$errors[[1L]]
        k257Log("G5 limpa", case, "error:", if (is.null(error)) "none" else conditionMessage(error), "dispatched:", rec$dispatched,
            "designated:", !is.null(norm$normalization_input_reference))
        expect_identical(k257ProtHead(fixture), head, info = case)
        expect_s3_class(error, "workflow_figure_unwritten")
        if (!is.null(error)) expect_match(conditionMessage(error), k257Unremoved, fixed = TRUE)
        expect_identical(rec$dispatched, 0L, info = case)
        expect_true(k257Holds(other), info = case)
    }
    expect_false(isTRUE(norm$normalization_complete))
})

# ------------------------------------------------------------------------------------------------
# Lipidomics and metabolomics: E5, E6, G5 (lipid and metab); G1, G2 (lipid: KIT-194 runs RUV on F172 for lipid only).
# ------------------------------------------------------------------------------------------------

for (class in c("lipid", "metab")) {
    test_that(paste("KIT-257 E5:", class, "in process, RUV skipped removes the earlier RUV-corrected figures; Copy lists none (fails before)"), {
        skip_if_not(workflowOwnerCapable())
        root <- withr::local_tempdir()
        seen <- k257Notes()
        fx <- k257OmicFixture(file.path(root, "a"), class)
        nd <- omicNormReactiveState()
        assays <- k257OmicAssays(fx)
        expect_gt(length(assays), 0L)
        figures <- k257Plant(k257OmicRuvFigures(fx$qc_dir, assays))
        expect_gt(length(figures), 0L)
        expect_null(shiny::isolate(nd$assay_names))
        k257OmicRun(fx, nd)
        info <- paste("E5", class, paste(seen$notes, collapse = " | "), paste(fx$log$lines, collapse = " | "))
        expect_identical(k257OmicSnapshot(fx)$head, paste0(fx$binding$prefix, "_normalized"), info = info)
        copy <- k257Copy(k257CopyPaths(fx$paths, root, fx$qc_dir), fx$omic)
        k257Log("E5", class, "present:", paste(k257Present(figures), collapse = ","), "QC_figures:", paste(copy$qc, collapse = ","))
        expect_identical(unname(k257Present(figures)), rep(FALSE, length(figures)), info = info)
        expect_false(any(grepl("_ruv_corrected_", copy$qc, fixed = TRUE)))
    })

    test_that(paste("KIT-257 E6:", class, "managed, the run's start removes the earlier RUV-corrected figures before the first action (fails before)"), {
        skip_if_not(workflowOwnerCapable())
        root <- withr::local_tempdir()
        seen <- k257Notes()
        fx <- k257OmicFixture(file.path(root, "a"), class)
        nd <- omicNormReactiveState()
        assays <- k257OmicAssays(fx)
        expect_gt(length(assays), 0L)
        figures <- k257Plant(k257OmicRuvFigures(fx$qc_dir, assays))
        expect_gt(length(figures), 0L)
        rec <- k257Recorder(figures)
        k257OmicRun(fx, nd, runner = rec$runner)
        k257Log("E6", class, "dispatched:", rec$dispatched, "present at the first dispatch:", paste(rec$present, collapse = ","))
        expect_gt(rec$dispatched, 0L)
        expect_identical(unname(rec$present), rep(FALSE, length(figures)), info = paste(seen$notes, collapse = " | "))
    })

    test_that(paste("KIT-257 G5:", class, "a file that cannot be removed refuses the run with nothing committed (fails before)"), {
        skip_if_not(workflowOwnerCapable())
        root <- withr::local_tempdir()
        seen <- k257Notes()
        fx <- k257OmicFixture(file.path(root, "a"), class)
        assays <- k257OmicAssays(fx)
        blocked <- k257OmicRuvFigures(fx$qc_dir, assays)[[1L]]
        dir.create(blocked, recursive = TRUE)
        other <- k257Plant(file.path(fx$qc_dir, sprintf("%s_pre_norm_pca.png", gsub("[^A-Za-z0-9]", "_", tolower(assays[[1L]])))))
        before <- k257OmicSnapshot(fx)
        k257OmicRun(fx, omicNormReactiveState())
        info <- paste("G5", class, paste(seen$notes, collapse = " | "))
        k257Log("G5", class, "history moved:", !identical(k257OmicSnapshot(fx), before), "notes:", paste(seen$notes, collapse = " | "))
        expect_identical(k257OmicSnapshot(fx), before, info = info)
        expect_true(any(grepl(k257Unremoved, seen$notes, fixed = TRUE) & grepl(blocked, seen$notes, fixed = TRUE)), info = info)
        expect_true(k257Holds(other), info = info)
    })
}

test_that("KIT-257 G1 (lipid): a skip then an RUV run writes new RUV-corrected figures", {
    skip_if_not(workflowOwnerCapable())
    root <- withr::local_tempdir()
    seen <- k257Notes()
    fx <- k257OmicFixture(file.path(root, "a"), "lipid")
    data_slot <- fx$config$binding$binding$data_slot
    qc_dir <- fx$qc_dir
    fx$config$binding$binding$qc_plots <- function(theObject, experiment_paths, stage, ...) {
        if (!identical(stage, "ruv_corrected")) return(invisible(NULL))
        for (file in k257OmicRuvFigures(qc_dir, names(methods::slot(theObject, data_slot)))) writeBin(charToRaw("new"), file)
        invisible(NULL)
    }
    nd <- omicNormReactiveState()
    k257OmicRun(fx, nd)
    figures <- k257Plant(k257OmicRuvFigures(qc_dir, k257OmicAssays(fx)))
    k257OmicRun(fx, nd, k257OmicInput(ruv = kit172LegacyRuv("manual")))
    info <- paste("G1 lipid", paste(seen$notes, collapse = " | "), paste(fx$log$lines, collapse = " | "))
    expect_identical(k257OmicSnapshot(fx)$head, "lipid_ruv_corrected", info = info)
    expect_gt(length(figures), 0L)
    expect_true(all(vapply(figures, k257Holds, NA, bytes = "new")), info = info)
})

test_that("KIT-257 G2 (lipid): users' files, another omic's folder and a sibling project matched by a bracket in the path are untouched", {
    skip_if_not(workflowOwnerCapable())
    root <- withr::local_tempdir()
    seen <- k257Notes()
    fx <- k257OmicFixture(file.path(root, "a"), "lipid")
    assays <- k257OmicAssays(fx)
    paths <- fx$paths
    paths$lipid_qc_dir <- file.path(root, "Study [pilot]", "lipid_qc")
    own <- k257Plant(k257OmicRuvFigures(paths$lipid_qc_dir, assays)[[1L]])
    name <- basename(own)
    user <- k257Plant(c(file.path(paths$lipid_qc_dir, c("my_notes.png", sub("[.]png$", "_v2.png", name))),
        file.path(root, "Study p", "lipid_qc", name), file.path(root, "Study [pilot]", "metabolite_qc", name)), "user")
    k257OmicRun(fx, omicNormReactiveState(), paths = paths)
    info <- paste("G2", paste(seen$notes, collapse = " | "))
    k257Log("G2", "own file present:", file.exists(own), "untouched:", paste(vapply(user, k257Holds, NA, bytes = "user"), collapse = ","))
    expect_identical(k257OmicSnapshot(fx)$head, "lipid_normalized", info = info)
    expect_true(all(vapply(user, k257Holds, NA, bytes = "user")), info = info)
    # Where the removal exists (the fix), the project's own file is removed through the bracket.
    if (!is.null(k257Fn("removeWorkflowOwnFiles"))) expect_false(file.exists(own))
})

# ------------------------------------------------------------------------------------------------
# Correction pass (AUDIT-KIT-257-impl-20261126): E8b (I01), G2 with a `~` path (I02), E10 (I05). Written before the
# repair; their base is the branch merged with merge train 8, before the repair (R2').
# ------------------------------------------------------------------------------------------------

k257ProtRuvCovariates <- function(fixture, norm, ruv_qc = function(reference) invisible(NULL)) {
    protNormCoreQuiet(runDIAReferenceRuvCovariates(list(ruv_mode = "covariates", color_variable = "group",
        shape_variable = NULL), fixture$workflow, norm, fixture$paths, "proteomics", "core", ruv_qc, function() "group"))
}

k257Md5 <- function(files) unname(tools::md5sum(files))

test_that("KIT-257 E8b (I01): DIA with limpa, a press again at correlation_filtered and at ruv_covariates keeps the current result's table and RUV figures, in process and managed (fails before)", {
    skip_if_not_installed("limpa")
    skip_if_not_installed("mixOmics")
    for (case in c("in process", "managed")) {
        fixture <- protNormW1LimpaOwner()
        norm <- new.env(parent = emptyenv())
        qc <- fixture$paths$protein_qc_dir
        first <- k257ProtTab(fixture, norm, "covariates")
        expect_length(first$errors, 0L)
        k257ProtCorrelation(fixture, norm, 0.5)
        tables <- file.path(qc, k257ProtTables[1:2])
        table <- tables[file.exists(tables)]
        expect_length(table, 1L)
        bytes <- k257Md5(table)
        rec <- k257Recorder(table)
        runner <- if (identical(case, "managed")) rec$runner
        again <- k257ProtTab(fixture, norm, "covariates", runner = runner)
        k257Log("E8b", case, "at correlation_filtered: errors", length(again$errors), "table present", file.exists(table),
            "head", fixture$manager$current_state)
        expect_length(again$errors, 0L)
        expect_identical(fixture$manager$current_state, "correlation_filtered", info = case)
        expect_true(file.exists(table), info = case)
        if (file.exists(table)) expect_identical(k257Md5(table), bytes, info = case)
        draw <- function(reference) for (name in k257ProtRuvFigures[1:4]) {
            writeBin(charToRaw("covariates"), file.path(norm$publication_paths$protein_qc_dir, name))
        }
        k257ProtRuvCovariates(fixture, norm, draw)
        expect_identical(fixture$manager$current_state, "ruv_covariates", info = case)
        figures <- file.path(qc, k257ProtRuvFigures[1:4])
        expect_true(all(file.exists(figures)), info = case)
        third <- k257ProtTab(fixture, norm, "covariates", runner = runner)
        k257Log("E8b", case, "at ruv_covariates: errors", length(third$errors), "table present", file.exists(table),
            "figures", paste(k257Present(figures), collapse = ","))
        expect_length(third$errors, 0L)
        expect_true(file.exists(table), info = case)
        expect_true(all(vapply(figures, k257Holds, NA, bytes = "covariates")), info = case)
    }
})

test_that("KIT-257 G2 (I02, lipid): a project path beginning with `~` removes the project's own file and the run proceeds (fails before)", {
    skip_if_not(workflowOwnerCapable())
    root <- withr::local_tempdir()
    home <- file.path(root, "home")
    dir.create(home)
    fx <- k257OmicFixture(file.path(root, "a"), "lipid")
    assays <- k257OmicAssays(fx)
    # HOME is a temporary folder for this test only (restored on exit): the real home directory is never touched.
    withr::local_envvar(HOME = home)
    expanded <- normalizePath(path.expand("~"), mustWork = FALSE)
    k257Log("G2 ~", "path.expand('~') follows HOME:", identical(expanded, normalizePath(home)))
    skip_if_not(identical(expanded, normalizePath(home)), "this R does not read HOME again after start-up")
    seen <- k257Notes()
    paths <- fx$paths
    paths$lipid_qc_dir <- "~/Study/lipid_qc"
    own <- k257Plant(k257OmicRuvFigures(file.path(home, "Study", "lipid_qc"), assays)[[1L]])
    user <- k257Plant(file.path(home, "Study", "lipid_qc", "my_notes.png"), "user")
    k257OmicRun(fx, omicNormReactiveState(), paths = paths)
    info <- paste("G2 ~", paste(seen$notes, collapse = " | "))
    k257Log("G2 ~", "own file present:", file.exists(own), "head:", k257OmicSnapshot(fx)$head, "notes:", paste(seen$notes, collapse = " | "))
    expect_identical(k257OmicSnapshot(fx)$head, "lipid_normalized", info = info)
    expect_false(file.exists(own), info = info)
    expect_true(k257Holds(user, "user"), info = info)
})

test_that("KIT-257 E10 (I05): a proteomics normalisation that stops after its start tells the user its files were removed (fails before)", {
    fixture <- protNormCoreOwner()
    norm <- new.env(parent = emptyenv())
    qc <- fixture$paths$protein_qc_dir
    k257ProtClick(fixture, norm, "skip")
    k257ProtCorrelation(fixture, norm, 0.5)
    expect_true(any(file.exists(file.path(qc, k257ProtTables[1:2]))))
    head <- k257ProtHead(fixture)
    # The tab's error handler shows the error's message ("Error in normalization: ..."); the recorder keeps the error.
    tab <- k257ProtTab(fixture, norm, "manual", settings = k257ProtSettings("manual", ruv_k = 50))
    shown <- if (length(tab$errors)) tab$errors[[1L]]$message else NA_character_
    k257Log("E10", "errors:", length(tab$errors), "shown:", shown)
    expect_length(tab$errors, 1L)
    expect_identical(k257ProtHead(fixture), head)
    expect_match(shown, "were removed when this run started", fixed = TRUE)
})

# ------------------------------------------------------------------------------------------------
# Merge obligation M1 (AUDIT-KIT-257-impl): KIT-243's DA removal routed through removeWorkflowOwnFiles. Written before
# the routing: on the merged tree KIT-243's core does not expand a leading `~`, so the removal refuses (fails before).
# ------------------------------------------------------------------------------------------------

test_that("KIT-257 M1 (I02 at the DA): the DA's own earlier figure is removed under a project path beginning with `~` (fails before)", {
    root <- withr::local_tempdir()
    home <- file.path(root, "home")
    dir.create(home)
    # HOME is a temporary folder for this test only (restored on exit): the real home directory is never touched.
    withr::local_envvar(HOME = home)
    skip_if_not(identical(normalizePath(path.expand("~"), mustWork = FALSE), normalizePath(home)),
        "this R does not read HOME again after start-up")
    graphs <- "~/Study/publication_graphs"
    own <- file.path(graphs, "Heatmaps", "posmode_groupB-groupA_heatmap.png")
    planted <- k257Plant(file.path(home, "Study", "publication_graphs", "Heatmaps", basename(own)))
    user <- k257Plant(file.path(home, "Study", "publication_graphs", "Heatmaps", "my_notes.png"), "user")
    # KIT-243 E16's form: an owner whose accepted manifest's names are `own`.
    testthat::local_mocked_bindings(usesWorkflowReferenceStates = function(manager) TRUE,
        isSavedWorkflowReference = function(reference) TRUE, workflowDaMadeFiles = function(manager, paths, omicType) own,
        .package = "MultiScholaR")
    manager <- list(.context = list(snapshot = list(da_results = list(entity_id = "kit257"))))
    error <- tryCatch({
        removeWorkflowDaOwnFigures(manager, list(publication_graphs_dir = graphs), "lipidomics")
        NULL
    }, error = function(e) e)
    k257Log("M1 ~", "error:", if (is.null(error)) "none" else conditionMessage(error), "own present:", file.exists(planted))
    expect_null(error)
    expect_false(file.exists(planted))
    expect_true(k257Holds(user, "user"))
})
