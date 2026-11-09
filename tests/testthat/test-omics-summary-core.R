# KIT-108 W3 (D3, D4): the Summary & Export tab's core (R/utils_workflow_summary.R), its entries tested directly.
#   - The 20 lipid, metab and trx functions D3's map moves are in the core file, unchanged by name; the core calls no
#     Shiny UI function.
#   - Core side (live): each entry, called directly in W1 (b)'s session on W1 (b)'s lipid, metab and trx fixtures in
#     place of the tab's button, gives W1 (b)'s frozen notifications, progress, arguments (the recorders' calls), the
#     `project_dirs` global and the files of that step: `saveWorkflowSummaryParameters` (in process, managed, the record
#     failing, trx's gate), `copyWorkflowSummaryResults` (complete, an optional artifact missing, required ones missing),
#     `renderWorkflowSummaryReport` (the template copied from the package, in the project, the GitHub fallback failing
#     offline; the record bundle failing), `writeWorkflowSummarySessionState` and `runWorkflowSummaryExportStep` (in
#     process and managed). F108-03: the save before trx's RUVIII export shows no "Saved Integration" notice.
#   - Comparison side: the start commit's handlers (the old modules' buttons, the same W1 (b) cases) give the same fields
#     as the core side, apart from F108-03's notice. One block, quarantined when W3 deletes the handlers (Q-175).
#   - The copy's classification by `copy$required`, and the template's ensure order (the project's copy, the package's,
#     the GitHub download; a failed retrieval stops the report unless `report$render_after_failed_template`).
#   - W5 (D8): the copy builds the combined DA workbook one sheet at a time; its R heap's peak stays below holding
#     both sheets' read data at once, and nothing larger than one sheet survives it.

# W1 (b)'s helpers and frozen captures: every top-level assignment of test-omic-summary-module.R (its blocks are not run).
for (omicsSummaryCoreW1Expression in parse(testthat::test_path("test-omic-summary-module.R"), keep.source = FALSE)) {
    if (is.call(omicsSummaryCoreW1Expression) && identical(omicsSummaryCoreW1Expression[[1L]], as.name("<-"))) {
        eval(omicsSummaryCoreW1Expression, environment())
    }
}
rm(omicsSummaryCoreW1Expression)

test_that("(W3) the moved lipid, metab and trx functions are in the core file; the core calls no Shiny UI function", {
    moved <- c("buildLipidSummarySessionState", "lipidSummaryFirstNonNull", "getLipidSummaryWorkflowValue",
        "getLipidSummaryStateObject", "getLipidSummaryObjectArgs", "buildLipidSummaryParameterPayload",
        "prepareLipidSummarySessionState", "metabSummaryFirstNonNull", "getMetabSummaryWorkflowValue",
        "getMetabSummaryStateObject", "getMetabSummaryObjectArgs", "buildMetabSummaryParameterPayload",
        "prepareMetabSummarySessionState", "trxRuviiiMChoices", "trxHeadHasTechnicalReplicates", "trxRuviiiStatusText",
        "trxSummaryContext", "trxSummaryStateObject", "trxSummaryIntegrationWrite", "trxSummaryIntegrationReady")
    entries <- c("workflowSummaryFinalState", "saveWorkflowSummaryParameters", "startWorkflowSummaryIntegrationExport",
        "copyWorkflowSummaryResults", "renderWorkflowSummaryReport", "writeWorkflowSummarySessionState",
        "workflowSummaryTemplateStatus", "runWorkflowSummaryExportStep", "trxRuviiiExportParameters",
        "trxSummarySessionState")
    expect_length(moved, 20L)
    ns <- asNamespace("MultiScholaR")
    for (name in c(moved, entries)) {
        expect_true(exists(name, envir = ns, inherits = FALSE), info = name)
        if (!exists(name, envir = ns, inherits = FALSE)) next
        fn <- get(name, envir = ns, inherits = FALSE)
        # Source files are known under load_all only; the installed package keeps no srcref (KIT-104 ca35c462's rule).
        if (!is.null(utils::getSrcref(fn))) expect_identical(utils::getSrcFilename(fn), "utils_workflow_summary.R", info = name)
    }
    if (!is.null(utils::getSrcref(get("saveWorkflowSummaryParameters", envir = ns)))) {
        source <- readLines(testthat::test_path("..", "..", "R", "utils_workflow_summary.R"))
        code <- sub("#.*$", "", source)
        expect_false(any(grepl("shiny::", code, fixed = TRUE)))
    }
})

# ------------------------------------------------------------------------------------------------
# Core side: the entries in W1 (b)'s session, in place of the tab's buttons.
# ------------------------------------------------------------------------------------------------

#' The fields an entry's step is compared on (the frame renders the outputs and holds the reactive values).
omicsSummaryCoreFields <- c("notes", "progress", "calls", "global", "files", "study_parameters", "integration",
    "session_state")

#' The cases, and the steps each replaces with its entry.
omicsSummaryCoreCases <- c("lipid_main", "lipid_variants", "lipid_required_missing", "lipid_github_fallback",
    "lipid_managed", "metab_main", "metab_variants", "metab_required_missing", "metab_github_fallback", "metab_managed",
    "trx_main", "trx_managed")

#' A W1 (b) case with its save, copy, report, export and RUVIII steps calling the core entries directly. The tab's
#' values the frame would hold (saved, copied, the report) are kept in `held`, as the frame keeps them.
omicsSummaryCoreSteps <- function(case, steps, setup) {
    omic <- c(lipid = "lipidomics", metab = "metabolomics", trx = "transcriptomics", prot = "proteomics")[[setup$omic]]
    held <- new.env(parent = emptyenv())
    held$saved <- FALSE
    held$copied <- FALSE
    held$generated <- FALSE
    label <- "kit108"
    config <- function() omicGuiConfig(omic)$summary
    save <- function(ctx) shiny::isolate(saveWorkflowSummaryParameters(config(), ctx$bus(), setup$project_dirs, omic, label,
        setup$description, notify = omicSummaryNotify, done = function(saved) held$saved <- TRUE,
        failed = function(error) held$saved <- TRUE))
    copy <- function(ctx) {
        copied <- shiny::isolate(shiny::withProgress(message = "Copying files to publication directory...",
            copyWorkflowSummaryResults(config(), ctx$bus(), setup$project_dirs, omic, label, setup$description,
                argsSaved = held$saved, notify = omicSummaryNotify)))
        if (isTRUE(copied$fallbackWritten)) held$saved <- TRUE
        if (is.null(copied$error)) held$copied <- TRUE
    }
    report <- function(ctx) {
        report <- shiny::isolate(renderWorkflowSummaryReport(config(), ctx$bus(), setup$project_dirs, omic, label,
            notify = omicSummaryNotify, progress = omicSummaryProgress,
            within = function(message, code) shiny::withProgress(message = message, code)))
        if (!is.null(report$path)) {
            held$generated <- TRUE
            held$path <- report$path
        }
    }
    export <- function(ctx) shiny::isolate(writeWorkflowSummarySessionState(config(), ctx$bus(), setup$project_dirs, omic,
        label, setup$description, values = list(workflow_args_saved = held$saved, files_copied = held$copied,
            report_generated = held$generated, report_path = held$path), notify = omicSummaryNotify))
    # The head's export as W1 (b) keeps it (`omicSummaryW1Ruviii`), read when the step is done: on a managed session
    # the runner delivers at the step's flush.
    accepted <- function(ctx) {
        export <- ctx$setup$owner$.context$snapshot$exports$trx_ruviii_export
        value <- if (!is.null(export)) hydrateSavedWorkflowEntity(ctx$setup$owner$.context$store, export$reference)$value
        ctx$keep("ruviii", if (is.null(value)) "<none>" else c(m_source = value$m_source, k = as.character(value$k),
            controls = as.character(value$control_count %||% sum(value$ctl_mask)),
            newY = sprintf("%dx%d %s", nrow(value$newY), ncol(value$newY), omicSummaryW1Short(value$newY))))
    }
    ruviii <- function(m_source) function(ctx) {
        step <- config()$export_steps[[1L]]
        accepted(ctx)
        shiny::isolate(runWorkflowSummaryExportStep(step, ctx$bus(), setup$project_dirs[[omic]], m_source,
            done = function(shown) accepted(ctx), binding = transcriptomicsDeclarationBinding,
            failed = function(error) shiny::showNotification(paste(step$error$prefix, trxConditionText(error)),
                type = step$error$type, duration = step$error$duration)))
    }
    fault <- function(name, value, then) function(ctx) {
        ctx$fault(name, value)
        then(ctx)
    }
    entries <- list(save = save, save_before_export = save, save_record_fails = fault("record", TRUE, save),
        copy = if (grepl("_variants$", case)) fault("record", FALSE, copy) else copy, report = report,
        report_fallback = report, report_bundle_fails = fault("bundle", TRUE, report), export = export,
        ruviii = ruviii("group"), ruviii_unavailable = ruviii("technical_replicate"))
    # A report refused at the tab's gate (the copy failed) is the frame's, not the entry's.
    if (grepl("_required_missing$", case)) entries$report <- NULL
    # KIT-108 W4: proteomics' owner variants, as W1 (b) drives them: the copy refused by the export gate (an export
    # pending on the owner), then accepted once the export is no longer pending.
    if (identical(case, "prot_owner_variants")) {
        pending <- function(ctx, value) {
            manager <- ctx$setup$manager()
            manager$.context$export_pending <- value
        }
        entries$copy_gate_refuses <- function(ctx) {
            ctx$fault("record", FALSE)
            pending(ctx, list(label = "session-export"))
            copy(ctx)
        }
        entries$copy <- function(ctx) {
            pending(ctx, NULL)
            copy(ctx)
        }
    }
    replaced <- intersect(names(steps), names(entries))
    steps[replaced] <- entries[replaced]
    list(steps = steps, replaced = replaced)
}

#' A case driven on the core: the captures of its entry steps.
omicsSummaryCoreRun <- function(case) {
    env <- environment()
    spec <- omicSummaryW1CaseSpec(case, env)
    core <- omicsSummaryCoreSteps(case, spec$steps, spec$setup)
    captures <- omicSummaryW1Drive(spec$setup, core$steps)
    list(captures = captures, replaced = core$replaced)
}

#' F108-03 (the lead's decision): the save before trx's RUVIII export shows the gate's warning and no "Saved
#' Integration" notice (W1 (b)'s capture re-frozen once at W3; this holds before and after that re-freeze).
omicsSummaryCoreExpected <- function(case, step, field) {
    frozen <- omicSummaryW1Frozen[[case]][[step]][[field]]
    if (identical(c(case, step, field), c("trx_main", "save_before_export", "notes"))) {
        frozen <- setdiff(frozen, "[message 5 -] Saved Integration S4 Object")
    }
    frozen
}

for (omicsSummaryCoreCase in omicsSummaryCoreCases) {
    local({
        case <- omicsSummaryCoreCase
        test_that(sprintf("(W3, core side) %s: each entry gives W1 (b)'s notifications, progress, arguments, global and files", case), {
            run <- omicsSummaryCoreRun(case)
            expect_gt(length(run$replaced), 0L)
            for (step in run$replaced) {
                for (field in omicsSummaryCoreFields) {
                    expect_identical(run$captures[[step]][[field]], omicsSummaryCoreExpected(case, step, field),
                        info = paste(case, step, field))
                }
                if (!is.null(omicSummaryW1Frozen[[case]][[step]]$ruviii)) {
                    expect_identical(run$captures[[step]]$ruviii, omicSummaryW1Frozen[[case]][[step]]$ruviii,
                        info = paste(case, step))
                }
            }
        })
    })
}
rm(omicsSummaryCoreCase)

# Comparison side: quarantined when W3 deletes the start commit's lipid, metab and trx handlers, Q-175.
test_that("(W3, comparison side) the start commit's handlers give the core entries' notifications, progress, arguments, global and files", {
    testthat::skip("quarantined: Q-175 (KIT-108 W3) the comparison side drives the start commit's lipid, metab and trx handlers (mod_lipid_summary_*_helpers.R, mod_metab_summary_*_helpers.R, mod_trx_summary.R, deleted by W3) through their modules' buttons; the core side's blocks above stay live and compare each entry with W1 (b)'s frozen steps, and W1 (b) drives the tab on the core")
    for (case in omicsSummaryCoreCases) {
        core <- omicsSummaryCoreRun(case)
        old <- omicSummaryW1Case(case)
        for (step in core$replaced) {
            for (field in omicsSummaryCoreFields) {
                expected <- old[[step]][[field]]
                # F108-03: the start commit's handler showed the notice after a refused write.
                if (identical(c(case, step, field), c("trx_main", "save_before_export", "notes"))) {
                    expect_true("[message 5 -] Saved Integration S4 Object" %in% expected)
                    expected <- setdiff(expected, "[message 5 -] Saved Integration S4 Object")
                }
                expect_identical(core$captures[[step]][[field]], expected, info = paste(case, step, field))
            }
        }
    }
})

# ------------------------------------------------------------------------------------------------
# KIT-108 W4: proteomics on the core (its route pieces, utils_workflow_dia_summary.R).
# ------------------------------------------------------------------------------------------------

#' W1 (b)'s proteomics cases with an entry step (the refusals are the frame's).
omicsSummaryCoreProtCases <- c("prot_route_dia", "prot_route_dia_limpa", "prot_route_tmt", "prot_route_lfq", "prot_stale",
    "prot_owner", "prot_owner_variants", "prot_native", "prot_github_fallback")

test_that("(W4) the seven moved proteomics functions and the route's plugins are in utils_workflow_dia_summary.R", {
    moved <- c("resolveProtSummaryFinalS4State", "buildProtSummaryTemplateStatus", "resolveProtSummaryExpectedTemplate",
        "resolveProtSummaryTemplateSelection", "resolveProtSummaryReportTemplate", "summariseProtSummaryWorkflowMetadata",
        "prepareProtSummarySessionStateExport")
    plugins <- c("protSummaryCopyGate", "protSummaryIntegrationFilename", "protSummaryIntegrationWrite",
        "protSummarySessionState")
    expect_length(moved, 7L)
    ns <- asNamespace("MultiScholaR")
    for (name in c(moved, plugins)) {
        expect_true(exists(name, envir = ns, inherits = FALSE), info = name)
        if (!exists(name, envir = ns, inherits = FALSE)) next
        fn <- get(name, envir = ns, inherits = FALSE)
        if (!is.null(utils::getSrcref(fn))) expect_identical(utils::getSrcFilename(fn), "utils_workflow_dia_summary.R", info = name)
    }
    if (!is.null(utils::getSrcref(get("protSummaryCopyGate", envir = ns)))) {
        code <- sub("#.*$", "", readLines(testthat::test_path("..", "..", "R", "utils_workflow_dia_summary.R")))
        expect_false(any(grepl("shiny::", code, fixed = TRUE)))
    }
})

for (omicsSummaryCoreCase in omicsSummaryCoreProtCases) {
    local({
        case <- omicsSummaryCoreCase
        test_that(sprintf("(W4, core side) %s: each entry gives W1 (b)'s notifications, progress, arguments, global and files", case), {
            run <- omicsSummaryCoreRun(case)
            expect_gt(length(run$replaced), 0L)
            for (step in run$replaced) {
                for (field in omicsSummaryCoreFields) {
                    expect_identical(run$captures[[step]][[field]], omicSummaryW1Frozen[[case]][[step]][[field]],
                        info = paste(case, step, field))
                }
            }
        })
    })
}
rm(omicsSummaryCoreCase)

test_that("(W4, comparison side) the start commit's proteomics handlers give the core entries' notifications, progress, arguments, global and files", {
    testthat::skip("quarantined: Q-183 (KIT-108 W4) the comparison side drives the start commit's proteomics handlers (mod_prot_summary_{server,publication,session,report}_helpers.R, deleted by W4) through mod_prot_summary_server's buttons; the core side's blocks above stay live and compare each entry with W1 (b)'s frozen steps, and W1 (b) drives the tab on the core")
    for (case in omicsSummaryCoreProtCases) {
        core <- omicsSummaryCoreRun(case)
        old <- omicSummaryW1Case(case)
        for (step in core$replaced) {
            for (field in omicsSummaryCoreFields) {
                expect_identical(core$captures[[step]][[field]], old[[step]][[field]], info = paste(case, step, field))
            }
        }
    }
})

# ------------------------------------------------------------------------------------------------
# The classification, the template's ensure order and the integration gate.
# ------------------------------------------------------------------------------------------------

#' A project's summary directories under a temporary base.
omicsSummaryCoreDirs <- function(omic = "lipidomics", env = parent.frame()) {
    base <- withr::local_tempdir(.local_envir = env)
    paths <- list(base_dir = base, source_dir = file.path(base, "scripts", omic), integration_dir = file.path(base, "integration"),
        results_summary_dir = file.path(base, "results_summary"))
    for (dir in paths) dir.create(dir, recursive = TRUE, showWarnings = FALSE)
    stats::setNames(list(paths), omic)
}

omicsSummaryCoreNotes <- function() {
    notes <- new.env(parent = emptyenv())
    notes$shown <- character()
    notes$progress <- character()
    notes$notify <- function(message, type, duration) {
        notes$shown <- c(notes$shown, sprintf("[%s %s] %s", type, if (is.null(duration)) "NULL" else format(duration), message))
    }
    notes$step <- function(value, detail) notes$progress <- c(notes$progress, paste(value, detail))
    notes
}

test_that("(W3) the copy classifies copyToResultsSummary's failures by copy$required: a required one stops the copy; NULL leaves them unclassified", {
    dirs <- omicsSummaryCoreDirs()
    withr::local_options(list())
    had <- exists("project_dirs", envir = .GlobalEnv, inherits = FALSE)
    old <- if (had) get("project_dirs", envir = .GlobalEnv)
    withr::defer(if (had) assign("project_dirs", old, envir = .GlobalEnv) else if (exists("project_dirs", envir = .GlobalEnv,
        inherits = FALSE)) rm("project_dirs", envir = .GlobalEnv))
    config <- omicGuiConfig("lipidomics")$summary
    failures <- list(list(display_name = "Design Matrix"), list(display_name = "Heatmaps"))
    copy <- function(...) failures
    notes <- omicsSummaryCoreNotes()
    out <- copyWorkflowSummaryResults(config, list(), dirs, "lipidomics", "kit108", "d", argsSaved = TRUE, notify = notes$notify,
        copyFn = copy)
    expect_identical(out$error, "Required publication artifacts missing: Design Matrix")
    expect_identical(notes$shown, "[error 10] Copy error: Required publication artifacts missing: Design Matrix")
    failures <- list(list(display_name = "Heatmaps"))
    notes <- omicsSummaryCoreNotes()
    out <- copyWorkflowSummaryResults(config, list(), dirs, "lipidomics", "kit108", "d", argsSaved = TRUE, notify = notes$notify,
        copyFn = copy)
    expect_null(out$error)
    expect_identical(out$optional, failures)
    expect_identical(notes$shown, "[warning 5] Publication files copied with optional missing artifacts")
    # Unclassified (proteomics' copy$required = NULL): copied whatever is missing.
    config$copy$required <- NULL
    failures <- list(list(display_name = "Design Matrix"))
    notes <- omicsSummaryCoreNotes()
    out <- copyWorkflowSummaryResults(config, list(), dirs, "lipidomics", "kit108", "d", argsSaved = TRUE, notify = notes$notify,
        copyFn = copy)
    expect_null(out$error)
    expect_identical(out$optional, list())
    expect_identical(notes$shown, "[message 5] Publication files copied")
    # Not saved: the fallback study parameters are written once, and reported.
    out <- copyWorkflowSummaryResults(config, list(), dirs, "lipidomics", "kit108", "d", argsSaved = FALSE,
        notify = notes$notify, copyFn = copy)
    expect_true(out$fallbackWritten)
    expect_match(readLines(file.path(dirs$lipidomics$source_dir, "study_parameters.txt"))[[4L]], "Workflow Name: kit108")
    expect_false(copyWorkflowSummaryResults(config, list(), dirs, "lipidomics", "kit108", "d", argsSaved = FALSE,
        notify = notes$notify, copyFn = copy)$fallbackWritten)
})

test_that("(W3) the report ensures its template in order: the project's copy, the package's, the GitHub download; a failed retrieval stops it unless render_after_failed_template", {
    config <- omicGuiConfig("lipidomics")$summary
    rendered <- character()
    render <- function(omic_type, experiment_label, rmd_filename) {
        rendered <<- c(rendered, rmd_filename)
        path <- file.path(tempdir(), "kit108-core-report.docx")
        writeLines("report", path)
        path
    }
    # The package's copy.
    dirs <- omicsSummaryCoreDirs()
    notes <- omicsSummaryCoreNotes()
    out <- renderWorkflowSummaryReport(config, NULL, dirs, "lipidomics", "kit108", notes$notify, notes$step, renderFn = render)
    expect_identical(out$templateSource, "package")
    expect_true(file.exists(file.path(dirs$lipidomics$source_dir, "lipidomics_report.rmd")))
    expect_identical(notes$progress, c("0.2 Checking for lipidomics_report.rmd template...", "0.5 Rendering report..."))
    expect_identical(notes$shown, c("[message 5] lipidomics_report.rmd template copied from package",
        "[message 5] Report generated successfully!"))
    # The project's copy: no check.
    notes <- omicsSummaryCoreNotes()
    out <- renderWorkflowSummaryReport(config, NULL, dirs, "lipidomics", "kit108", notes$notify, notes$step, renderFn = render)
    expect_identical(out$templateSource, "existing")
    expect_identical(notes$progress, "0.5 Rendering report...")
    # The check's progress always (proteomics' flag).
    always <- config
    always$report$checking_progress_always <- TRUE
    notes <- omicsSummaryCoreNotes()
    renderWorkflowSummaryReport(always, NULL, dirs, "lipidomics", "kit108", notes$notify, notes$step, renderFn = render)
    expect_identical(notes$progress, c("0.2 Checking for lipidomics_report.rmd template...", "0.5 Rendering report..."))
    # The GitHub download, when the package has no template (a directory it does not ship).
    missing <- config
    missing$report$dir <- "kit108-absent"
    downloads <- character()
    download <- function(url, destfile, quiet) {
        downloads <<- c(downloads, url)
        writeLines("template", destfile)
    }
    dirs <- omicsSummaryCoreDirs()
    notes <- omicsSummaryCoreNotes()
    out <- renderWorkflowSummaryReport(missing, NULL, dirs, "lipidomics", "kit108", notes$notify, notes$step,
        renderFn = render, downloadFileFn = download)
    expect_identical(out$templateSource, "github")
    expect_identical(downloads,
        "https://raw.githubusercontent.com/APAF-bioinformatics/MultiScholaR/main/inst/reports/kit108-absent/lipidomics_report.rmd")
    # A failed retrieval stops the report, or (metab's flag) renders anyway.
    offline <- function(url, destfile, quiet) stop("offline")
    for (render_after in c(FALSE, TRUE)) {
        missing$report$render_after_failed_template <- render_after
        dirs <- omicsSummaryCoreDirs()
        notes <- omicsSummaryCoreNotes()
        rendered <- character()
        out <- renderWorkflowSummaryReport(missing, NULL, dirs, "lipidomics", "kit108", notes$notify, notes$step,
            renderFn = render, downloadFileFn = offline)
        expect_null(out$templateSource)
        expect_identical(notes$shown[[1L]], "[error 10] Template retrieval failed: offline")
        expect_identical(length(rendered), as.integer(render_after), info = render_after)
    }
    # Project directories not initialised: refused before the progress display.
    notes <- omicsSummaryCoreNotes()
    started <- FALSE
    out <- renderWorkflowSummaryReport(config, NULL, list(), "lipidomics", "kit108", notes$notify, notes$step,
        within = function(message, code) started <<- TRUE, renderFn = render)
    expect_null(out$path)
    expect_false(started)
    expect_identical(notes$shown, "[error 10] Error: Project directories not properly initialized")
})

test_that("(W3, F108-03) a save the integration gate refuses shows the gate's warning and no \"Saved Integration\" notice; an accepted export writes and says so", {
    config <- omicGuiConfig("transcriptomics")$summary
    dirs <- omicsSummaryCoreDirs("transcriptomics")
    paths <- dirs$transcriptomics
    written <- character()
    config$final_state$context <- function(workflowData, catFn) list(finalS4Object = NULL, contrastsTbl = NULL,
        dataStateUsed = "trx_ruvs")
    config$integration$write <- function(workflowData, path) {
        written <<- c(written, path)
        writeLines("integration", path)
    }
    args <- function(...) file.path(paths$source_dir, "study_parameters.txt")
    record <- function(manager, path) invisible(path)
    for (accepted in c(FALSE, TRUE)) {
        manager <- new.env()
        manager$.context <- list(snapshot = list(exports = if (accepted) list(trx_ruviii_export = list(reference = "r"))))
        notes <- omicsSummaryCoreNotes()
        saved <- NULL
        saveWorkflowSummaryParameters(config, list(state_manager = manager), dirs, "transcriptomics", "kit108", "d",
            notify = notes$notify, done = function(value) saved <<- value, failed = function(error) stop(error),
            createWorkflowArgsFn = args, writeAnalysisRecordFn = record)
        if (accepted) {
            expect_identical(notes$shown, c("[message 5] Saved Integration S4 Object", "[message 5] Study parameters saved successfully"))
            expect_identical(saved$integrationFile, "transcriptomics_normalized_data.tsv")
            expect_identical(written, file.path(paths$integration_dir, "transcriptomics_normalized_data.tsv"))
        } else {
            expect_identical(notes$shown, c(paste("[warning 10] The integration export needs the RUVIII export: run it in the",
                "panel above, then save again."), "[message 5] Study parameters saved successfully"))
            expect_identical(saved$integrationFile, "not written (no integration export yet)")
            expect_identical(written, character())
        }
    }
})

# ------------------------------------------------------------------------------------------------
# KIT-108 W4: F104-01 (D5), test 1 (the core): shown failing on the W3 merge 1cb2a79b (kit-108/w4/test-first.txt).
# ------------------------------------------------------------------------------------------------

#' A native owner of a route (`test-workflow-native-session.R`'s fixture, :5-21, as the LFQ block at :500-512 collects
#' it): the store and the bus typed by the route, its final state, the native format.
omicsSummaryCoreNativeOwner <- function(type, env = parent.frame()) {
    manager <- newDIAStateAdapter(newWorkflowEntityStore(withr::local_tempdir(.local_envir = env), "kit108-native", type),
        WorkflowState$new(audit_enabled = FALSE))
    withr::defer(closeDIAStateAdapter(manager), envir = env)
    manager$setWorkflowType(type)
    final <- makeWorkflowCheckpointProtein(readWorkflowCheckpointRecipe("small"))
    manager$saveState("correlation_filtered", final, list(), "final")
    manager$.context$session_format <- "native_v1"
    workflow <- new.env(parent = emptyenv())
    workflow$state_manager <- manager
    workflow$config_list <- list(globalParameters = list(workflow_type = type))
    workflow$contrasts_tbl <- data.frame(contrasts = "B-A", friendly_names = "B_vs_A")
    workflow$design_matrix <- final@design_matrix
    workflow
}

test_that("(W4, F104-01) the native export's metadata names the route's template as the report step resolves it", {
    expected <- c(DIA = "DIANN_report.rmd", DIA_LIMPA = "DIANN_limpa_report.rmd", TMT = "TMT_report.rmd",
        LFQ = "LFQ_report.rmd")
    for (type in names(expected)) {
        workflow <- omicsSummaryCoreNativeOwner(type)
        norm <- new.env(parent = emptyenv())
        norm$best_k <- NA
        norm$correlation_threshold <- 0.5
        collected <- withDIAStateScope(workflow$state_manager, function()
            collectDIANativeSessionMetadata(workflow, norm, list(norm_method = "cyclicloess", ruv_mode = "skip")))
        resolved <- resolveProtSummaryReportTemplate(workflow, catFn = function(...) invisible(NULL))$templateFilename
        expect_identical(resolved, expected[[type]], info = type)
        expect_identical(collected$metadata$report_template, resolved, info = type)
    }
})

test_that("(W4, F104-01) the summary's archive refresh fills a missing stored template with the route's, and keeps a set one", {
    refresh <- function(stored, type) {
        source_dir <- withr::local_tempdir(.local_envir = parent.frame())
        metadata <- list(final_protein_count = 120L, final_sample_count = 12L,
            contrasts_tbl = data.frame(friendly_names = "B_vs_A"), normalization_method = "cyclicloess",
            ruv_mode = "skip", ruv_k = NA, correlation_threshold = 0.5, workflow_type = type, report_template = stored)
        manager <- new.env(parent = emptyenv())
        manager$.context <- list(owned_export = list(metadata_reference = list(entity_id = "metadata"),
            generation = "earlier", extra_references = list()), head = list(entity_id = "current"))
        manager$current_state <- "correlation_filtered"
        manager$state_history <- list("correlation_filtered")
        workflow <- list(state_manager = manager, config_list = list(globalParameters = list(workflow_type = type)))
        exported <- NULL
        testthat::with_mocked_bindings(
            suppressMessages(refreshDIANativeSummaryArchive(workflow, source_dir)),
            assertDIAExportReady = function(...) invisible(TRUE),
            resolveDIAStateValue = function(context, reference) metadata,
            collectDIANativeGuiFiles = function(...) list(),
            runDIANativeSessionExport = function(manager, metadata, source_dir, extra_references, dependencies,
                prepare_outputs) {
                exported <<- metadata
                prepare_outputs(metadata, list(archive_name = "filtered_session_data_kit108.zip"))
                list(status = "complete")
            })
        list(metadata = exported$report_template,
            text = grep("Report Template", readLines(file.path(source_dir, "filtered_session_summary.txt")), value = TRUE))
    }
    filled <- refresh(NULL, "TMT")
    expect_identical(filled$metadata, "TMT_report.rmd")
    expect_identical(filled$text, "- Report Template: TMT_report.rmd")
    kept <- refresh("default", "TMT")
    expect_identical(kept$metadata, "default")
    expect_identical(kept$text, "- Report Template: default")
})

test_that("(W5, D8) the copy's DA workbook is built one sheet at a time: its heap peak stays below both sheets' read data held at once", {
    # A synthetic DA directory in the shape the copy reads: 2 assays x 2,000 features x 200 samples, as xlsx.
    root <- withr::local_tempdir()
    da_dir <- file.path(root, "da")
    dir.create(da_dir)
    set.seed(108L)
    for (assay in c("gcms", "posmode")) {
        n <- 2000L
        sheet <- data.frame(metabolite_id = as.character(seq_len(n)),
            metabolite_name = sprintf("Metabolite_%s_%d", assay, seq_len(n)), logFC = stats::rnorm(n),
            raw_pvalue = stats::runif(n), fdr_qvalue = stats::runif(n),
            significant = sample(c("Up", "Down", "NS"), n, replace = TRUE), assay = assay,
            comparison = "groupWT-groupKO", stringsAsFactors = FALSE)
        intensities <- matrix(stats::rnorm(n * 200L, 20), n, 200L,
            dimnames = list(NULL, sprintf("intensity.S%03d", seq_len(200L))))
        openxlsx::write.xlsx(cbind(sheet, as.data.frame(intensities)),
            file.path(da_dir, sprintf("de_%s_metabolites_groupWT-groupKO_long_annot.xlsx", assay)))
    }
    # The heap is measured in a fresh R whose collection triggers are set high (R_VSIZE, R_NSIZE), so a collection runs
    # only where the code asks for one: gc's max used then counts all the garbage alive at once, independent of R's
    # adaptive triggers and of this process's history. A first copy warms the namespaces and the byte-compiled closures.
    package_path <- getNamespaceInfo("MultiScholaR", "path")
    measured <- callr::r(function(package_path, root) {
        # An installed package also has DESCRIPTION and R/; only an installed one has Meta/.
        if (!dir.exists(file.path(package_path, "Meta"))) {
            suppressMessages(pkgload::load_all(package_path, quiet = TRUE, helpers = FALSE))
        } else {
            loadNamespace("MultiScholaR", lib.loc = dirname(package_path))
        }
        copy_function <- get("copyToResultsSummary", envir = asNamespace("MultiScholaR"))
        da_files <- list.files(file.path(root, "da"), full.names = TRUE)
        used_mb <- function(g) sum(g[, 2])
        max_used_mb <- function(g) sum(g[, 6])
        copy_once <- function(run) {
            run_dir <- file.path(root, run)
            paths <- list(results_dir = file.path(run_dir, "results"),
                results_summary_dir = file.path(run_dir, "results_summary", "metabolomics_W5"),
                publication_graphs_dir = file.path(run_dir, "publication_graphs"), time_dir = file.path(run_dir, "time"),
                qc_dir = file.path(run_dir, "qc"), da_output_dir = file.path(root, "da"),
                pathway_dir = file.path(run_dir, "pathway"), source_dir = file.path(run_dir, "source"),
                feature_qc_dir = file.path(run_dir, "feature_qc"))
            for (path in paths) dir.create(path, recursive = TRUE, showWarnings = FALSE)
            copy <- copy_function
            environment(copy) <- list2env(list(getProjectPaths = function(...) paths), parent = environment(copy_function))
            utils::capture.output(suppressWarnings(copy("metabolomics", "W5", data.frame(contrasts = "B-A"),
                data.frame(Sample_ID = "S1"), force = TRUE)))
            file.path(paths$results_summary_dir, "Publication_tables", "DA_results_metabolomics.xlsx")
        }
        copy_once("warm")
        before <- gc(reset = TRUE)
        workbook <- copy_once("measured")
        after <- gc()
        copy <- c(peak = max_used_mb(after) - used_mb(before), residue = used_mb(after) - used_mb(before))
        before <- gc(reset = TRUE)
        sheets <- lapply(da_files, openxlsx::read.xlsx)
        after <- gc()
        list(copy = copy, both = max_used_mb(after) - used_mb(before),
            one_sheet = as.numeric(utils::object.size(sheets[[1L]])) / 2^20,
            sheet_names = openxlsx::getSheetNames(workbook),
            index = openxlsx::read.xlsx(workbook, sheet = "DA_Results_Index"))
    }, args = list(package_path, root), libpath = .libPaths(),
        env = c(callr::rcmd_safe_env(), R_VSIZE = "3G", R_NSIZE = "40M"))
    expect_identical(measured$sheet_names, c("DA_Results_Index", "DA_Sheet1", "DA_Sheet2"))
    expect_identical(measured$index, data.frame(Sheet = c("DA_Sheet1", "DA_Sheet2"),
        Description = c("gcms_metabolites_groupWT-groupKO", "posmode_metabolites_groupWT-groupKO")))
    expect_lt(measured$copy[["peak"]], measured$both,
        label = sprintf("the copy's gc max used (%.1f MB)", measured$copy[["peak"]]),
        expected.label = sprintf("both sheets' read data held at once (%.1f MB)", measured$both))
    expect_lt(measured$copy[["residue"]], measured$one_sheet,
        label = sprintf("the heap left after the copy (%.1f MB)", measured$copy[["residue"]]),
        expected.label = sprintf("one sheet's data (%.1f MB)", measured$one_sheet))
})
