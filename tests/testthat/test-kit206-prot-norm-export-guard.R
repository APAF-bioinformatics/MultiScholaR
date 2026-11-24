# KIT-206 (VAL-134 (1); parity with VAL-105 (2)): the proteomics session exports refuse when the current table is not
# the normalisation tab's last result ("Export Filtered Data Session", the Summary's archive refresh, the S4
# compatibility download), and the native export names its results file and RUV statement from the run (part (b)).
# On the native DIA owner of helper-prot-norm-core.R (the "small" recipe, committed as `protein_replicate_filtered`)
# and its DIA_LIMPA owner. The runs, Reverts and the QC change are KIT-182's (test-prot-norm-input-current.R), copied.
# "Export" is the tab's observer inside the owner's state scope, its notices, messages and errors recorded; on the
# native owner the export workflow is the real one (its progress display replaced), elsewhere a recorder. The new
# function is called directly only where it exists, after the observer's assertions, so the tests-only commit runs the
# same tests and fails at a value. Lines starting `KIT-206` are the evidence of docs/qa/workflow-execution/kit/kit-206/.

k206QcParameters <- list(strict_mode = FALSE, min_reps_per_group = 2, min_groups = 2,
    proteins_intensity_cutoff_percentile = 10)

k206Settings <- function(mode = "skip", ...) utils::modifyList(list(norm_method = "quantile", ruv_mode = mode,
    auto_percentage_min = NULL, auto_percentage_max = NULL, separation_metric = NULL, k_penalty_weight = NULL,
    adaptive_k_penalty = NULL, ruv_percentage = 50, ruv_k = 1, color_variable = "group", shape_variable = NULL), list(...))

k206Aesthetics <- function() list(color_var = "group", shape_var = NULL)

#' A call with only the arguments `fn` has.
k206Call <- function(fn, args) do.call(fn, args[intersect(names(args), names(formals(fn)))])

k206Scope <- function(fixture, code) withDIAStateScope(fixture$manager, code)

k206Head <- function(fixture) getDIAStateReference(fixture$manager)

k206Revert <- function(fixture, state) k206Scope(fixture, function() selectDIAStateReference(fixture$manager, state))

k206Error <- function(code) tryCatch({
    force(code)
    NULL
}, error = function(e) e)

#' The native DIA owner with the contrasts table the design step sets (KIT-182 F1 (h), VAL-138).
k206Owner <- function(envir = parent.frame()) {
    fixture <- protNormCoreOwner(envir)
    fixture$workflow$contrasts_tbl <- data.frame(contrasts = "grouptreated-groupcontrol")
    fixture
}

#' The tab's in-process normalisation click.
k206Click <- function(fixture, norm, mode = "skip") {
    protNormCoreQuiet(k206Call(runDIAReferenceNormalization, list(settings = k206Settings(mode),
        workflow = fixture$workflow, norm_data = norm, paths = fixture$paths, omic_type = "proteomics",
        experiment_label = "core", post_qc = function(reference) invisible(NULL),
        ruv_qc = function(reference) invisible(NULL), grouping = function() "group", pre_qc = NULL)))
}

#' The correlation filter in process (Apply, or Skip with `skipped`), and the flag the tab then sets.
k206Correlation <- function(fixture, norm, threshold = 0.5, skipped = FALSE) {
    final <- protNormCoreQuiet(runDIAReferenceCorrelation(fixture$workflow, norm, fixture$paths, threshold, skipped,
        "proteomics", "core", if (skipped) function() NULL else function() "group", k206Aesthetics))
    norm$correlation_filtering_complete <- TRUE
    final
}

#' KIT-182's QC change: the Protein Intensity Filter on the head through its kernel, committed as the step commits it,
#' its display published; should it remove no protein, the head's object less its first protein row. The new head.
k206QcChange <- function(fixture) {
    manager <- fixture$manager
    rlang::local_bindings(filtering_progress = methods::new("FilteringProgress"), .env = .GlobalEnv)
    grDevices::pdf(NULL)
    device <- grDevices::dev.cur()
    on.exit(if (device %in% grDevices::dev.list()) grDevices::dev.off(device), add = TRUE)
    before <- k206Scope(fixture, function() readDIAStateMetadata(manager, k206Head(fixture))$data$dimensions$rows)
    protNormCoreQuiet(k206Scope(fixture, function() {
        value <- computeDIAReferenceProteinIntensityFilter(manager, getDIAStateReference(manager), k206QcParameters)
        if (identical(as.integer(value$protein_count), as.integer(before))) {
            object <- manager$getState()
            object@protein_quant_table <- object@protein_quant_table[-1L, , drop = FALSE]
            manager$saveState("protein_intensity_filtered", object, list(), "qc change")
            return(invisible(NULL))
        }
        result <- finishDIAReferenceProteinStep(fixture$workflow, value, publish = FALSE)
        publishDIAProteinProgress(result, "11_protein_intensity_filtered", "proteomics", "core", updateProteinFiltering)
        invisible(NULL)
    }))
    expect_identical(manager$current_state, "protein_intensity_filtered")
    k206Head(fixture)
}

#' The DA's re-commit of the head's arguments (KIT-182 F1 (f)).
k206Da <- function(fixture) protNormCoreQuiet(k206Scope(fixture, function() commitDIADaArguments(fixture$manager,
    list(da_q_val_thresh = 0.05, treat_lfc_cutoff = 0, formula_string = "~ 0 + group"), fixture$workflow,
    data.frame(contrasts = "groupB-groupA", full_format = "B_vs_A=groupB-groupA", stringsAsFactors = FALSE),
    "Protein.Ids")))

#' A completed enrichment on the head, as the enrichment module commits its result class (test-workflow-integration-
#' input.R's fixture); an error is returned, not raised.
k206Enrichment <- function(fixture) k206Error({
    results <- methods::new("EnrichmentResults", contrasts = tibble::tibble(contrasts = "groupB-groupA"),
        enrichment_data = list(`groupB-groupA` = list(up = data.frame(term_id = "GO:0008150", p_value = 0.01))),
        enrichment_plots = list(), enrichment_plotly = list(), enrichment_summaries = list())
    protNormCoreQuiet(fixture$manager$saveState("enrichment_completed", results, list(q_cutoff = 0.05),
        "fixture enrichment"))
})

k206PassThrough <- function(fixture, norm) protNormCoreQuiet(runDIAReferenceLimpaNormalisationPassThrough(
    k206Settings("covariates"), fixture$workflow, norm, function(reference) invisible(NULL)))

k206RuvCovariates <- function(fixture, norm, mode) protNormCoreQuiet(runDIAReferenceRuvCovariates(list(ruv_mode = mode,
    color_variable = "group", shape_variable = NULL), fixture$workflow, norm, fixture$paths, "proteomics", "core",
    function(reference) invisible(NULL), function() "group"))

#' The ticket's 7 text (VAL-41, proposed) after the export's prefix.
k206Text <- paste0("the current table is not the normalisation tab's last result (from the correlation filter, its ",
    "Skip, or on DIA with limpa the RUV covariates step): a step was applied or reverted since that result was made, ",
    "so the exported session would not describe the current table. Press \"Run Normalization & RUV\" again, apply or ",
    "skip the correlation filter (on DIA with limpa, then the RUV covariates step), then export.")
k206ExportNote <- paste("[error 10] Error exporting session:", k206Text)

#' An export workflow that only records its calls (the non-native owners; the guard's decision alone).
k206Recorder <- function() {
    seen <- new.env(parent = emptyenv())
    seen$calls <- 0L
    seen$run <- function(workflowData, normData, input, sourceDir) {
        seen$calls <- seen$calls + 1L
        list(exportArtifacts = list(sessionFilename = "recorded"))
    }
    seen
}

#' "Export Filtered Data Session": the tab's observer inside the owner's scope with the run's RUV mode; the real export
#' workflow unless `recorder` is given. Its notices (`[type duration] text`) and errors.
k206Export <- function(fixture, norm, mode = "skip", recorder = NULL) {
    log <- new.env(parent = emptyenv())
    log$notes <- character()
    log$errors <- list()
    note <- function(ui, type = "default", duration = 5, ...) {
        log$notes <- c(log$notes, sprintf("[%s %s] %s", type, format(duration), paste(as.character(ui), collapse = " ")))
        invisible(NULL)
    }
    quiet <- function(...) invisible(NULL)
    run <- if (is.null(recorder)) function(workflowData, normData, input, sourceDir) runProtNormExportSessionWorkflow(
        workflowData, normData, input, sourceDir, withProgressFn = function(message, value, expr) force(expr),
        incProgressFn = function(...) NULL) else recorder$run
    log$value <- protNormCoreQuiet(withDIAStateScope(fixture$manager, function() runProtNormExportObserver(
        list(norm_method = "quantile", ruv_mode = mode), fixture$workflow, norm, fixture$paths,
        runExportSessionWorkflowFn = run, showNotificationFn = note, messageFn = quiet,
        handleExportErrorFn = function(e) {
            log$errors <- c(log$errors, list(e))
            handleProtNormExportError(e, showNotificationFn = note, messageFn = quiet)
        })))
    as.list(log)
}

k206Archives <- function(paths) sort(list.files(paths$source_dir, pattern = "^filtered_session_.*[.]msdia[.]tar$",
    full.names = TRUE))

#' What an export leaves in the source folder: the archives, the latest pointer and the session summary.
k206Snapshot <- function(paths) {
    read <- function(name) {
        path <- file.path(paths$source_dir, name)
        if (file.exists(path)) readLines(path, warn = FALSE) else NULL
    }
    list(archives = k206Archives(paths), latest = read("filtered_session_latest.json"),
        summary = read("filtered_session_summary.txt"))
}

#' One archive, opened into a private folder: its recorded current state, members (destinations), correlation result
#' and RUV statement.
k206Archive <- function(path) {
    root <- tempfile("k206-archive-")
    dir.create(root)
    imported <- protNormCoreQuiet(openDIANativeSession(path, root))
    on.exit(closeDIAStateAdapter(imported$manager), add = TRUE)
    list(state = imported$metadata$r6_current_state_name,
        destinations = vapply(imported$manifest$dependencies, function(dependency) dependency$destination,
            character(1), USE.NAMES = FALSE),
        correlation = imported$extra_references$correlation_filtered$entity_id,
        ruv_mode = imported$metadata$ruv_mode, ruv_applied = imported$metadata$ruv_applied)
}

k206States <- function(paths) paste(vapply(paths, function(path) k206Archive(path)$state, character(1)), collapse = ",")

#' The archive opened into a new owner and project by the package's own path (KIT-182 F1 (h): `openDIANativeSession`,
#' `adoptDIANativeSession`, as the DA tab's Load Filtered Session does), and the tab's restore there.
k206Reopen <- function(archive, envir = parent.frame()) {
    root <- withr::local_tempdir(.local_envir = envir)
    paths <- list(source_dir = file.path(root, "scripts", "proteomics"),
        protein_qc_dir = file.path(root, "results", "proteomics", "protein_qc"),
        peptide_qc_dir = file.path(root, "results", "proteomics", "peptide_qc"))
    paths$publication_graphs_dir <- file.path(root, "results", "proteomics", "publication_graphs")
    paths$time_dir <- file.path(paths$publication_graphs_dir, "filtering_qc")
    for (path in unlist(paths)) dir.create(path, recursive = TRUE, showWarnings = FALSE)
    session <- new.env(parent = emptyenv())
    session$userData <- new.env(parent = emptyenv())
    session$onFlush <- session$onFlushed <- session$onSessionEnded <- function(...) function() invisible(NULL)
    workflow <- new.env(parent = emptyenv())
    workflow$state_manager <- WorkflowState$new(audit_enabled = FALSE)
    workflow$qc_params <- list()
    imported <- protNormCoreQuiet(openDIANativeSession(archive, paths$source_dir))
    protNormCoreQuiet(adoptDIANativeSession(imported, workflow, session, paths$source_dir, paths))
    withr::defer(closeDIAStateAdapter(workflow$state_manager), envir = envir)
    norm <- new.env(parent = emptyenv())
    testthat::local_mocked_bindings(updateSelectInput = function(...) invisible(NULL), .package = "shiny")
    protNormCoreQuiet(restoreProtNormTabState(NULL, workflow, norm))
    list(fixture = list(manager = workflow$state_manager, workflow = workflow, paths = paths), norm = norm)
}

#' A refused export: the notice, nothing written or moved in the source folder, nothing pending; then, where the
#' function exists, the guard's class on the same head.
k206ExpectRefused <- function(case, log, fixture, norm, before) {
    after <- k206Snapshot(fixture$paths)
    added <- setdiff(after$archives, before$archives)
    cat(sprintf("KIT-206 %s head=%s notes=%s new_archives=%d archive_states=%s pending=%s\n", case,
        fixture$manager$current_state, paste(log$notes, collapse = " | "), length(added), k206States(added),
        !is.null(fixture$manager$.context$export_pending)))
    expect_identical(log$notes, k206ExportNote, info = case)
    expect_identical(after, before, info = case)
    expect_null(fixture$manager$.context$export_pending, info = case)
    if (exists("assertProtNormExportHeadCurrent", mode = "function")) {
        error <- k206Error(assertProtNormExportHeadCurrent(fixture$manager, norm$correlation_filtered_obj))
        expect_s3_class(error, "prot_norm_correlation_input_stale")
    }
}

#' The archive an export wrote, checked for its recorded current state.
k206ExpectExported <- function(case, log, paths, before, state) {
    added <- setdiff(k206Archives(paths), before$archives)
    cat(sprintf("KIT-206 %s notes=%s new_archives=%d archive_states=%s\n", case, paste(log$notes, collapse = " | "),
        length(added), k206States(added)))
    expect_length(log$errors, 0L)
    expect_length(added, 1L)
    if (length(added) == 1L) {
        archive <- k206Archive(added)
        expect_identical(archive$state, state, info = case)
        return(invisible(archive))
    }
    invisible(NULL)
}

# ------------------------------------------------------------------------------------------------
# Part (a): the export refuses a stale table (E1 to E3, E9, E11, E12 fail before; E4 to E6, E8, E10 pass on both)
# ------------------------------------------------------------------------------------------------

test_that("KIT-206 E1: after a Revert past normalisation (S1) the export is refused and nothing is written (fails before)", {
    fixture <- k206Owner()
    norm <- new.env(parent = emptyenv())
    k206Click(fixture, norm, "skip")
    k206Correlation(fixture, norm, 0.5)
    k206Revert(fixture, "protein_replicate_filtered")
    before <- k206Snapshot(fixture$paths)
    expect_identical(before, list(archives = character(), latest = NULL, summary = NULL))
    log <- k206Export(fixture, norm, "skip")
    k206ExpectRefused("E1", log, fixture, norm, before)
})

test_that("KIT-206 E2: after a Revert and a QC step applied again (S2) the export is refused (fails before)", {
    fixture <- k206Owner()
    norm <- new.env(parent = emptyenv())
    k206Click(fixture, norm, "skip")
    k206Correlation(fixture, norm, 0.5)
    k206Revert(fixture, "protein_replicate_filtered")
    k206QcChange(fixture)
    before <- k206Snapshot(fixture$paths)
    log <- k206Export(fixture, norm, "skip")
    k206ExpectRefused("E2", log, fixture, norm, before)
})

test_that("KIT-206 E2b: S2, then a second normalisation without the correlation filter: the export is refused (fails before)", {
    fixture <- k206Owner()
    norm <- new.env(parent = emptyenv())
    k206Click(fixture, norm, "skip")
    k206Correlation(fixture, norm, 0.5)
    first <- norm$correlation_filtered_obj$entity_id
    k206Revert(fixture, "protein_replicate_filtered")
    k206QcChange(fixture)
    k206Click(fixture, norm, "skip")
    expect_identical(fixture$manager$current_state, "normalized")
    before <- k206Snapshot(fixture$paths)
    log <- k206Export(fixture, norm, "skip")
    added <- setdiff(k206Archives(fixture$paths), before$archives)
    if (length(added) == 1L) cat(sprintf("KIT-206 E2b archive correlation is the first run's=%s\n",
        identical(k206Archive(added)$correlation, first)))
    k206ExpectRefused("E2b", log, fixture, norm, before)
})

test_that("KIT-206 E3: after RUV and two Reverts (head `normalized`, before RUV) the export is refused (fails before)", {
    fixture <- k206Owner()
    norm <- new.env(parent = emptyenv())
    k206Click(fixture, norm, "manual")
    k206Correlation(fixture, norm, 0.5)
    k206Revert(fixture, "normalized")
    before <- k206Snapshot(fixture$paths)
    log <- k206Export(fixture, norm, "manual")
    k206ExpectRefused("E3", log, fixture, norm, before)
})

test_that("KIT-206 E4: at the correlation result the first export of a fresh project writes the archive (guard)", {
    fixture <- k206Owner()
    norm <- new.env(parent = emptyenv())
    k206Click(fixture, norm, "skip")
    k206Correlation(fixture, norm, 0.5)
    before <- k206Snapshot(fixture$paths)
    expect_null(before$summary)
    log <- k206Export(fixture, norm, "skip")
    archive <- k206ExpectExported("E4", log, fixture$paths, before, "correlation_filtered")
    expect_true(all(c("protein_qc/normalised_results_cln_with_replicates.tsv", "source/filtered_session_summary.txt") %in%
        archive$destinations))
    expect_null(fixture$manager$.context$export_pending)
})

test_that("KIT-206 E5: after Skip the export writes; on DIA with limpa after the RUV covariates it proceeds (guard)", {
    fixture <- k206Owner()
    norm <- new.env(parent = emptyenv())
    k206Click(fixture, norm, "skip")
    k206Correlation(fixture, norm, 0, skipped = TRUE)
    before <- k206Snapshot(fixture$paths)
    log <- k206Export(fixture, norm, "skip")
    k206ExpectExported("E5 skip", log, fixture$paths, before, "correlation_filtered")
    skip_if_not_installed("limpa")
    skip_if_not_installed("mixOmics")
    fixture <- protNormW1LimpaOwner()
    norm <- new.env(parent = emptyenv())
    k206PassThrough(fixture, norm)
    k206Correlation(fixture, norm, 0.5)
    k206RuvCovariates(fixture, norm, "covariates")
    expect_identical(fixture$manager$current_state, "ruv_covariates")
    recorder <- k206Recorder()
    log <- k206Export(fixture, norm, "covariates", recorder)
    cat(sprintf("KIT-206 E5 limpa notes=%s calls=%d\n", paste(log$notes, collapse = " | "), recorder$calls))
    expect_identical(recorder$calls, 1L)
    expect_length(log$errors, 0L)
})

test_that("KIT-206 E6: after the DA's re-commit the export writes; KIT-182's enrichment stand-in proceeds (guard)", {
    fixture <- k206Owner()
    norm <- new.env(parent = emptyenv())
    k206Click(fixture, norm, "skip")
    k206Correlation(fixture, norm, 0.5)
    k206Da(fixture)
    expect_identical(fixture$manager$current_state, "correlation_filtered")
    before <- k206Snapshot(fixture$paths)
    log <- k206Export(fixture, norm, "skip")
    k206ExpectExported("E6 da", log, fixture$paths, before, "correlation_filtered")
    # KIT-182 T9 (a)'s stand-in: the correlation's record committed as `enrichment_completed` (the components branch
    # only: the stand-in keeps the protein record); the guard's decision through the observer, with a recorder.
    k206Scope(fixture, function() {
        record <- readDIAStateMetadata(fixture$manager)
        commitDIAReferenceState(fixture$manager, "enrichment_completed", record$data, record$config, "fixture")
    })
    recorder <- k206Recorder()
    log <- k206Export(fixture, norm, "skip", recorder)
    cat(sprintf("KIT-206 E6 stand-in notes=%s calls=%d\n", paste(log$notes, collapse = " | "), recorder$calls))
    expect_identical(recorder$calls, 1L)
})

test_that("KIT-206 E7: the S4 compatibility download is refused after S1, before the export, with a notice (fails before)", {
    fixture <- k206Owner()
    norm <- new.env(parent = emptyenv())
    k206Click(fixture, norm, "skip")
    k206Correlation(fixture, norm, 0.5)
    k206Revert(fixture, "protein_replicate_filtered")
    seen <- new.env(parent = emptyenv())
    seen$states <- character()
    seen$notes <- character()
    testthat::local_mocked_bindings(
        downloadHandler = function(filename, content, ...) list(filename = filename, content = content),
        showNotification = function(ui, action = NULL, duration = 5, closeButton = TRUE, id = NULL,
            type = c("default", "message", "warning", "error"), session = NULL) {
            seen$notes <- c(seen$notes, sprintf("[%s %s] %s", type[[1L]], format(duration),
                paste(as.character(ui), collapse = " ")))
            invisible(id %||% "notification")
        }, .package = "shiny")
    # The cgroup guard of the diagnostic runner cannot hold on a runner: the export itself is recorded.
    testthat::local_mocked_bindings(exportDIAS4Compatibility = function(manager, path, state_name = NULL) {
        seen$states <- c(seen$states, manager$current_state)
        saveRDS("compatibility object", path)
        invisible(path)
    }, .package = "MultiScholaR")
    handler <- protNormS4DownloadHandler(fixture$workflow, norm, fixture$paths, "proteomics", "core")
    error <- k206Error(protNormCoreQuiet(handler$content(file.path(withr::local_tempdir(), "download.RDS"))))
    cat(sprintf("KIT-206 E7 exported_states=%s notes=%s error=%s\n", paste(seen$states, collapse = ","),
        paste(seen$notes, collapse = " | "), if (is.null(error)) "none" else class(error)[[1L]]))
    expect_s3_class(error, "prot_norm_correlation_input_stale")
    expect_identical(seen$states, character())
    expect_identical(seen$notes, k206ExportNote)
})

test_that("KIT-206 E8: an owner without reference states and an owned-payload result pass to the export (guard)", {
    # (i) An in-memory owner and a correlation object that is not a saved reference.
    root <- withr::local_tempdir()
    plain <- list(manager = WorkflowState$new(audit_enabled = FALSE), paths = list(source_dir = root))
    plain$workflow <- list2env(list(state_manager = plain$manager))
    norm <- new.env(parent = emptyenv())
    norm$correlation_filtering_complete <- TRUE
    norm$correlation_filtered_obj <- "filtered_obj"
    recorder <- k206Recorder()
    log <- k206Export(plain, norm, "skip", recorder)
    expect_identical(recorder$calls, 1L)
    expect_length(log$errors, 0L)
    # (ii) A reference owner and a list `normData` whose result is not a saved reference.
    fixture <- k206Owner()
    owned <- new.env(parent = emptyenv())
    k206Click(fixture, owned, "skip")
    k206Correlation(fixture, owned, 0.5)
    k206Revert(fixture, "protein_replicate_filtered")
    listed <- list(correlation_filtering_complete = TRUE, correlation_filtered_obj = "filtered_obj")
    recorder <- k206Recorder()
    log <- k206Export(fixture, listed, "skip", recorder)
    expect_identical(recorder$calls, 1L)
    expect_length(log$errors, 0L)
})

test_that("KIT-206 E9: an exported project reopened into a new owner exports again; after a Revert there it is refused (second part fails before)", {
    fixture <- k206Owner()
    norm <- new.env(parent = emptyenv())
    k206Click(fixture, norm, "skip")
    k206Correlation(fixture, norm, 0.5)
    before <- k206Snapshot(fixture$paths)
    log <- k206Export(fixture, norm, "skip")
    k206ExpectExported("E9 first", log, fixture$paths, before, "correlation_filtered")
    reopened <- k206Reopen(k206Archives(fixture$paths))
    expect_false(identical(reopened$fixture$manager$.context$store$directory, fixture$manager$.context$store$directory))
    expect_identical(reopened$fixture$manager$current_state, "correlation_filtered")
    before <- k206Snapshot(reopened$fixture$paths)
    log <- k206Export(reopened$fixture, reopened$norm, "skip")
    k206ExpectExported("E9 reopened", log, reopened$fixture$paths, before, "correlation_filtered")
    k206Revert(reopened$fixture, "protein_replicate_filtered")
    before <- k206Snapshot(reopened$fixture$paths)
    log <- k206Export(reopened$fixture, reopened$norm, "skip")
    k206ExpectRefused("E9 reverted", log, reopened$fixture, reopened$norm, before)
})

test_that("KIT-206 E10: DIA with limpa, the tab restored after the RUV covariates, exports; also after the DA and an enrichment (guard)", {
    skip_if_not_installed("limpa")
    skip_if_not_installed("mixOmics")
    fixture <- protNormW1LimpaOwner()
    norm <- new.env(parent = emptyenv())
    k206PassThrough(fixture, norm)
    k206Correlation(fixture, norm, 0.5)
    k206RuvCovariates(fixture, norm, "covariates")
    expect_identical(fixture$manager$current_state, "ruv_covariates")
    testthat::local_mocked_bindings(updateSelectInput = function(...) invisible(NULL), .package = "shiny")
    protNormCoreQuiet(restoreProtNormTabState(NULL, fixture$workflow, norm))
    expect_identical(norm$correlation_filtered_obj$entity_id,
        getDIAStateReference(fixture$manager, "correlation_filtered")$entity_id)
    differ <- k206Scope(fixture, function() !identical(readDIALineageRecord(fixture$manager, k206Head(fixture))$components,
        readDIALineageRecord(fixture$manager, norm$correlation_filtered_obj)$components))
    cat(sprintf("KIT-206 E10 ruv_covariates and correlation_filtered components differ=%s\n", differ))
    recorder <- k206Recorder()
    k206Export(fixture, norm, "covariates", recorder)
    expect_identical(recorder$calls, 1L)
    k206Da(fixture)
    expect_identical(fixture$manager$current_state, "ruv_covariates")
    k206Export(fixture, norm, "covariates", recorder)
    expect_identical(recorder$calls, 2L)
    error <- k206Enrichment(fixture)
    cat(sprintf("KIT-206 E10 enrichment head built=%s\n", is.null(error)))
    if (is.null(error)) {
        k206Export(fixture, norm, "covariates", recorder)
        expect_identical(recorder$calls, 3L)
    }
})

test_that("KIT-206 E11: the Summary's archive refresh is refused after a Revert past normalisation; without it, it writes (fails before)", {
    fixture <- k206Owner()
    norm <- new.env(parent = emptyenv())
    k206Click(fixture, norm, "skip")
    k206Correlation(fixture, norm, 0.5)
    k206ExpectExported("E11 export", k206Export(fixture, norm, "skip"), fixture$paths, k206Snapshot(fixture$paths),
        "correlation_filtered")
    k206Da(fixture)
    k206Revert(fixture, "protein_replicate_filtered")
    before <- k206Snapshot(fixture$paths)
    error <- k206Error(protNormCoreQuiet(refreshDIANativeSummaryArchive(fixture$workflow, fixture$paths$source_dir)))
    after <- k206Snapshot(fixture$paths)
    added <- setdiff(after$archives, before$archives)
    cat(sprintf("KIT-206 E11 error=%s new_archives=%d archive_states=%s latest_moved=%s\n",
        if (is.null(error)) "none" else class(error)[[1L]], length(added), k206States(added),
        !identical(after$latest, before$latest)))
    expect_s3_class(error, "prot_norm_correlation_input_stale")
    expect_identical(after$archives, before$archives)
    expect_identical(after$latest, before$latest)
    # Without the Revert the refresh writes a new archive at the DA's re-commit and moves the latest pointer to it.
    fixture <- k206Owner()
    norm <- new.env(parent = emptyenv())
    k206Click(fixture, norm, "skip")
    k206Correlation(fixture, norm, 0.5)
    k206Export(fixture, norm, "skip")
    k206Da(fixture)
    before <- k206Snapshot(fixture$paths)
    error <- k206Error(protNormCoreQuiet(refreshDIANativeSummaryArchive(fixture$workflow, fixture$paths$source_dir)))
    after <- k206Snapshot(fixture$paths)
    added <- setdiff(after$archives, before$archives)
    cat(sprintf("KIT-206 E11 guard error=%s new_archives=%d archive_states=%s\n",
        if (is.null(error)) "none" else conditionMessage(error), length(added), k206States(added)))
    expect_null(error)
    expect_length(added, 1L)
    expect_false(identical(after$latest, before$latest))
})

test_that("KIT-206 E12: S3, the DA, then an enrichment result on the head: the export judges the table it was made from (fails before)", {
    fixture <- k206Owner()
    norm <- new.env(parent = emptyenv())
    k206Click(fixture, norm, "manual")
    k206Correlation(fixture, norm, 0.5)
    k206Revert(fixture, "normalized")
    k206Da(fixture)
    error <- k206Enrichment(fixture)
    if (!is.null(error)) {
        cat(sprintf("KIT-206 E12 enrichment head not built: %s\n", conditionMessage(error)))
        skip("an enrichment result could not be committed on the fixture")
    }
    record <- k206Scope(fixture, function() tryCatch(readComponentMetadata(fixture$manager, k206Head(fixture)),
        error = function(e) NULL))
    cat(sprintf("KIT-206 E12 head=%s class=%s previous_state=%s\n", fixture$manager$current_state,
        format(record$data$class), format(record$previous_state)))
    before <- k206Snapshot(fixture$paths)
    log <- k206Export(fixture, norm, "manual")
    k206ExpectRefused("E12", log, fixture, norm, before)
})

# ------------------------------------------------------------------------------------------------
# Part (b): the native export follows the run, and a missing file leaves nothing pending (E13, fails before)
# ------------------------------------------------------------------------------------------------

test_that("KIT-206 E13: the export takes its results file and RUV statement from the run; a missing file is refused with nothing pending (fails before)", {
    fixture <- k206Owner()
    norm <- new.env(parent = emptyenv())
    k206Click(fixture, norm, "skip")
    k206Correlation(fixture, norm, 0.5)
    before <- k206Snapshot(fixture$paths)
    # The RUV setting changed to Automatic after a run that skipped RUV.
    log <- k206Export(fixture, norm, "automatic")
    added <- setdiff(k206Archives(fixture$paths), before$archives)
    archive <- if (length(added) == 1L) k206Archive(added) else NULL
    applied <- k206Error(k206Correlation(fixture, norm, 0.5))
    cat(sprintf("KIT-206 E13 radio notes=%s ruv_mode=%s ruv_applied=%s pending=%s then_correlation=%s\n",
        paste(log$notes, collapse = " | "), format(archive$ruv_mode), format(archive$ruv_applied),
        !is.null(fixture$manager$.context$export_pending), if (is.null(applied)) "ran" else conditionMessage(applied)))
    expect_length(log$errors, 0L)
    expect_length(added, 1L)
    expect_true("protein_qc/normalised_results_cln_with_replicates.tsv" %in% archive$destinations)
    expect_identical(archive$ruv_applied, FALSE)
    expect_null(applied)
    # A file the export needs is missing: refused before the request exists, nothing pending, the next step runs.
    fixture <- k206Owner()
    norm <- new.env(parent = emptyenv())
    k206Click(fixture, norm, "skip")
    k206Correlation(fixture, norm, 0.5)
    results <- file.path(fixture$paths$protein_qc_dir, "normalised_results_cln_with_replicates.tsv")
    expect_true(file.exists(results))
    unlink(results)
    log <- k206Export(fixture, norm, "skip")
    pending <- !is.null(fixture$manager$.context$export_pending)
    applied <- k206Error(k206Correlation(fixture, norm, 0.5))
    cat(sprintf("KIT-206 E13 missing notes=%s class=%s pending=%s then_correlation=%s\n",
        paste(log$notes, collapse = " | "), if (length(log$errors)) class(log$errors[[1L]])[[1L]] else "none", pending,
        if (is.null(applied)) "ran" else conditionMessage(applied)))
    expect_length(log$notes, 1L)
    expect_match(log$notes, paste("[error 10] Error exporting session: the file normalised_results_cln_with_replicates.tsv",
        "that this export needs is missing from"), fixed = TRUE)
    expect_match(log$notes, "; run the step that writes it, then export again.", fixed = TRUE)
    expect_length(log$errors, 1L)
    if (length(log$errors) == 1L) expect_s3_class(log$errors[[1L]], "workflow_export_dependency_missing")
    expect_false(pending)
    expect_null(applied)
})
