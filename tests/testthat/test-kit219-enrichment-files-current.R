# KIT-219 (VAL-179 SD-5, VAL-187): Copy and the Report present only the enrichment files the current enrichment made;
# the lipidomics and metabolomics report templates print their enrichment section only when the publication folder holds
# an enrichment file. The tests-only commit (no R/ change): each E row is written to fail on the base at the value the
# ticket states (section 10) and to pass on the fix; the G rows hold on both. The new function
# (`workflowEnrichmentCurrentFiles`) and the new record (`enrichment_not_current`) are looked at only where they exist,
# after the base assertions. Lines starting `KIT-219` are the evidence of docs/qa/workflow-execution/kit/kit-219/.
#
# The proteomics rows run on the native DIA owner of helper-prot-norm-core.R (the "small" recipe): normalisation,
# correlation filter and export as the Summary finds them, a real DA manifest, an enrichment committed through the
# application's commit (`saveWorkflowEnrichmentState`), the real `copyWorkflowSummaryResults` and
# `renderWorkflowSummaryReport` (a recording render), with the real `copyToResultsSummary` inside the owner's scope.
# The enrichment's files are written as the test kernel writes them (`writeProtEnrichDeterministicResultTable`) and, for
# the plots, by the real kernel's naming; the real kernel's names are the subject of G0 (kit-219/g0.R).

library(testthat)

# --- Helpers copied from tests/testthat/test-kit216-prot-session-after-enrichment.R (lines named) ---------------------

# :17-19.
k219Settings <- function(mode = "skip", norm_method = "quantile") list(norm_method = norm_method, ruv_mode = mode,
    auto_percentage_min = NULL, auto_percentage_max = NULL, separation_metric = NULL, k_penalty_weight = NULL,
    adaptive_k_penalty = NULL, ruv_percentage = 50, ruv_k = 1, color_variable = "group", shape_variable = NULL)

k219Aesthetics <- function() list(color_var = "group", shape_var = NULL) # :21

k219Call <- function(fn, args) do.call(fn, args[intersect(names(args), names(formals(fn)))]) # :23

k219Scope <- function(fixture, code) withDIAStateScope(fixture$manager, code) # :25

k219Text <- function(x) if (is.null(x) || !length(x)) "<none>" else paste(format(x), collapse = ",") # :127

# :37-39.
k219Contrasts <- data.frame(contrasts = "grouptreated-groupcontrol",
    full_format = "treated_vs_control=grouptreated-groupcontrol", friendly_names = "treated_vs_control",
    stringsAsFactors = FALSE)

# :42-45.
k219Owner <- function(envir = parent.frame()) {
    fixture <- protNormCoreOwner(envir)
    fixture$workflow$contrasts_tbl <- k219Contrasts
    fixture
}

# :49-55.
k219Click <- function(fixture, norm, mode = "skip", norm_method = "quantile") {
    protNormCoreQuiet(k219Call(runDIAReferenceNormalization, list(settings = k219Settings(mode, norm_method),
        workflow = fixture$workflow, norm_data = norm, paths = fixture$paths, omic_type = "proteomics",
        experiment_label = "core", post_qc = function(reference) invisible(NULL),
        ruv_qc = function(reference) invisible(NULL), grouping = function() "group", pre_qc = NULL)))
}

# :57-63.
k219Correlation <- function(fixture, norm, threshold = 0.5, skipped = FALSE) {
    final <- protNormCoreQuiet(runDIAReferenceCorrelation(fixture$workflow, norm, fixture$paths, threshold, skipped,
        "proteomics", "core", if (skipped) function() NULL else function() "group", k219Aesthetics))
    norm$correlation_filtering_complete <- TRUE
    final
}

# :90-110, without the recorder: "Export Filtered Data Session" in process (the archive Copy reads).
k219Export <- function(fixture, norm, mode = "skip", norm_method = "quantile") {
    errors <- list()
    note <- function(ui, type = "default", duration = 5, ...) invisible(NULL)
    quiet <- function(...) invisible(NULL)
    run <- function(workflowData, normData, input, sourceDir) runProtNormExportSessionWorkflow(
        workflowData, normData, input, sourceDir, withProgressFn = function(message, value, expr) force(expr),
        incProgressFn = function(...) NULL)
    protNormCoreQuiet(withDIAStateScope(fixture$manager, function() runProtNormExportObserver(
        list(norm_method = norm_method, ruv_mode = mode), fixture$workflow, norm, fixture$paths,
        runExportSessionWorkflowFn = run, showNotificationFn = note, messageFn = quiet,
        handleExportErrorFn = function(e) {
            errors <<- c(errors, list(e))
            handleProtNormExportError(e, showNotificationFn = note, messageFn = quiet)
        })))
    expect_length(errors, 0L)
    invisible(errors)
}

# :131-137: a DA on the head (the manifest returned; no publication).
k219RealDa <- function(fixture, q = 0.05) {
    da_data <- new.env(parent = emptyenv())
    da_data$current_s4_object <- k219Scope(fixture, function() getDIAStateReference(fixture$manager))
    da_data$contrasts_available <- k219Contrasts$contrasts
    protNormCoreQuiet(runDIAReferenceDaAnalysis(list(formula_string = "~ 0 + group", da_q_val_thresh = q,
        treat_lfc_cutoff = 0), da_data, fixture$workflow, fixture$paths, contrasts = k219Contrasts, publish = FALSE))
}

# --- New helpers ---------------------------------------------------------------------------------------------------

k219Texts <- list(
    t1_lead = "Publication files copied.",
    t1 = "were not copied, because the current enrichment did not make them",
    t1b = "this application has no enrichment analysis for",
    t2 = "holds enrichment results that the current enrichment did not make")

#' The Summary's folders for the proteomics fixture: its paths, the results folders `copyToResultsSummary` reads, and an
#' empty `results_summary_dir`.
k219Dirs <- function(fixture) {
    root <- fixture$root
    results <- file.path(root, "results", "proteomics")
    paths <- c(fixture$paths, list(results_dir = results, results_summary_dir = file.path(root, "results_summary"),
        qc_dir = file.path(results, "qc"), da_output_dir = file.path(results, "da_output"),
        pathway_dir = file.path(results, "pathway_enrichment"), feature_qc_dir = fixture$paths$protein_qc_dir,
        integration_dir = file.path(root, "integration"), base_dir = root))
    for (path in unlist(paths)) dir.create(path, recursive = TRUE, showWarnings = FALSE)
    writeLines("KIT-219 study parameters", file.path(paths$source_dir, "study_parameters.txt"))
    list(proteomics = paths)
}

#' The default path to the Summary: normalised (RUV skipped), correlation filter 0.5, exported, the DA run on the head;
#' the contrast is the manifest's first result name (printed by the rows).
k219Journey <- function(envir = parent.frame()) {
    fixture <- k219Owner(envir)
    norm <- new.env(parent = emptyenv())
    k219Click(fixture, norm)
    k219Correlation(fixture, norm)
    k219Export(fixture, norm)
    fixture$norm <- norm
    fixture$dirs <- k219Dirs(fixture)
    fixture$manifest <- k219RealDa(fixture)
    fixture$contrast <- (names(fixture$manifest$results) %||% "grouptreated-groupcontrol")[[1L]]
    fixture$pathway <- fixture$dirs$proteomics$pathway_dir
    fixture
}

#' The names one enrichment run writes for a contrast and its directions (the test kernel's and the real kernel's naming:
#' the table with the contrast as written, the plots with the sanitised contrast).
k219Names <- function(contrast, directions = c("up", "down")) {
    safe <- gsub("[^A-Za-z0-9_.-]", "_", contrast)
    as.character(unlist(lapply(directions, function(direction) c(paste0(contrast, "_", direction, "_enrichment_results.tsv"),
        paste0(safe, "_", direction, "_enrichment_plot.", c("html", "png"))))))
}

#' An `EnrichmentResults` with a result for each of `directions` of one contrast.
k219Results <- function(contrast, directions) methods::new("EnrichmentResults",
    contrasts = tibble::tibble(contrasts = contrast),
    enrichment_data = stats::setNames(list(stats::setNames(lapply(directions, function(direction)
        data.frame(term_id = "GO:0008150", p_value = 0.01)), directions)), contrast),
    enrichment_plots = list(), enrichment_plotly = list(), enrichment_summaries = list())

#' One enrichment run: the product committed through the application's commit (`saveWorkflowEnrichmentState`, the
#' manifest's source passed) and the run's files written into `pathway_dir` (`content` marks the run).
k219Enrich <- function(fixture, content = "run1", directions = c("up", "down"), contrast = fixture$contrast) {
    parameters <- list(selected_contrast = "treated_vs_control", up_cutoff = 0, down_cutoff = 0, q_cutoff = 0.05,
        organism_taxid = "9606")
    declared <- list(taxon_id = "9606", up_cutoff = 0, down_cutoff = 0, q_cutoff = 0.05, correction_method = "gSCS")
    saved <- protNormCoreQuiet(k219Scope(fixture, function() k219Call(saveWorkflowEnrichmentState, list(
        workflow = fixture$workflow, config = omicGuiConfig("proteomics")$enrichment,
        results = k219Results(contrast, directions), declared = declared,
        method = list(method = "gprofiler", supported = TRUE), parameters = parameters, pathway_dir = fixture$pathway,
        source = fixture$manifest$source))))
    expect_true(isTRUE(saved$saved), info = paste("enrichment commit:", k219Text(saved$warning)))
    safe <- gsub("[^A-Za-z0-9_.-]", "_", contrast)
    for (direction in directions) {
        writeProtEnrichDeterministicResultTable(data.frame(term_id = "GO:0008150", description = content), fixture$pathway,
            contrast, direction)
        for (extension in c("html", "png")) {
            writeLines(paste("KIT-219", content), file.path(fixture$pathway, paste0(safe, "_", direction,
                "_enrichment_plot.", extension)))
        }
    }
    invisible(saved)
}

#' Files placed in `pathway_dir` outside the run (an earlier run's, or another tool's).
k219Plant <- function(fixture, names, content = "planted") {
    for (name in names) writeLines(paste("KIT-219", content, name), file.path(fixture$pathway, name))
    invisible(names)
}

#' The `project_dirs` global as it was before the calling test (Copy assigns it when absent; k212KeepGlobal's form).
k219KeepGlobal <- function(env = parent.frame()) {
    had <- exists("project_dirs", envir = .GlobalEnv, inherits = FALSE)
    old <- if (had) get("project_dirs", envir = .GlobalEnv, inherits = FALSE)
    withr::defer(if (had) assign("project_dirs", old, envir = .GlobalEnv) else
        if (exists("project_dirs", envir = .GlobalEnv, inherits = FALSE)) rm("project_dirs", envir = .GlobalEnv), envir = env)
}

#' "Copy to Publication Directory": the core's `copyWorkflowSummaryResults` under the omic's group with the real
#' `copyToResultsSummary` (inside `scope`, the owner's scope, for proteomics), recording the notices ("[type] text"), the
#' arguments of the copy call and its return value.
k219CopyCore <- function(workflow, omic, dirs, scope = NULL, env = parent.frame()) {
    rec <- new.env(parent = emptyenv())
    rec$notes <- character()
    rec$args <- NULL
    rec$returned <- NULL
    rec$value <- NULL
    k219KeepGlobal(env)
    assign("project_dirs", dirs, envir = .GlobalEnv)
    copy <- function(...) {
        args <- list(...)
        rec$args <- args
        rec$returned <- if (is.null(scope)) do.call(copyToResultsSummary, args) else
            withDIAStateScope(scope, function() do.call(copyToResultsSummary, args))
        invisible(rec$returned)
    }
    testthat::local_mocked_bindings(workflowOmicsManagedActions = function(omic_type) NULL, .package = "MultiScholaR")
    suppressWarnings(suppressMessages(utils::capture.output(rec$value <- copyWorkflowSummaryResults(
        omicGuiConfig(omic)$summary, workflow, dirs, omic, "k219", "d", argsSaved = TRUE,
        notify = function(message, type, duration) rec$notes <- c(rec$notes, paste0("[", type, "] ", message)),
        copyFn = copy, catFn = function(...) NULL))))
    rec
}

k219Copy <- function(fixture, env = parent.frame()) k219CopyCore(fixture$workflow, "proteomics", fixture$dirs,
    scope = fixture$manager, env = env)

#' `copyToResultsSummary` as the workbooks call it: no `enrichment_files`, on the fixture's project.
k219CopyDirect <- function(fixture, env = parent.frame()) {
    k219KeepGlobal(env)
    assign("project_dirs", fixture$dirs, envir = .GlobalEnv)
    returned <- NULL
    suppressWarnings(suppressMessages(utils::capture.output(returned <- withDIAStateScope(fixture$manager, function()
        copyToResultsSummary(omic_type = "proteomics", experiment_label = "k219",
            contrasts_tbl = fixture$workflow$contrasts_tbl, design_matrix = fixture$workflow$design_matrix, force = TRUE)))))
    returned
}

#' "Generate Report" in process with a recorded render and a template the test supplies (k216Report's form,
#' test-kit216-prot-session-after-enrichment.R:270-286, with the Summary's own folders).
k219Report <- function(fixture) {
    calls <- new.env(parent = emptyenv())
    calls$render <- 0L
    calls$bundle <- 0L
    calls$notes <- character()
    calls$report <- file.path(fixture$root, "k219-report.html")
    notify <- function(message, type, duration) calls$notes <- c(calls$notes, paste0("[", type, "] ", message))
    protNormCoreQuiet(utils::capture.output(renderWorkflowSummaryReport(omicGuiConfig("proteomics")$summary,
        fixture$workflow, fixture$dirs, "proteomics", "k219", notify = notify, progress = function(...) NULL,
        renderFn = function(omic_type, experiment_label, rmd_filename) {
            calls$render <- calls$render + 1L
            writeLines("KIT-219 report", calls$report)
            calls$report
        }, exportRecordFn = function(...) {
            calls$bundle <- calls$bundle + 1L
            invisible(NULL)
        }, downloadFileFn = function(url, destfile, ...) {
            writeLines("KIT-219 template", destfile)
            invisible(0L)
        }, catFn = function(...) NULL)))
    as.list(calls)
}

k219Publication <- function(dirs, omic) file.path(dirs[[omic]]$results_summary_dir, "Publication_figures", "Enrichment_Plots")

#' The files the publication folder holds, relative to it.
k219Held <- function(dirs, omic) {
    folder <- k219Publication(dirs, omic)
    if (!dir.exists(folder)) character() else sort(list.files(folder, recursive = TRUE))
}

#' The enrichment workbook's sheets: the number of result sheets and the index's (contrast, direction) pairs.
k219Workbook <- function(dirs, omic) {
    path <- file.path(dirs[[omic]]$results_summary_dir, "Publication_tables", paste0("Pathway_enrichment_results_", omic,
        ".xlsx"))
    if (!file.exists(path)) return(list(present = FALSE, n = NA_integer_, pairs = character()))
    sheets <- openxlsx::getSheetNames(path)
    index <- openxlsx::read.xlsx(path, "Enrichment_Index")
    index <- index[grepl("^Enrich_Sheet", as.character(index$Sheet)), , drop = FALSE]
    list(present = TRUE, n = sum(grepl("^Enrich_Sheet", sheets)), pairs = paste(index$Contrast, index$Direction))
}

#' The record `copyToResultsSummary` returned for enrichment files left out, or NULL.
k219Record <- function(returned) {
    found <- Filter(function(failure) is.list(failure) && identical(failure$type, "enrichment_not_current"), returned)
    if (length(found)) found[[1L]]
}

k219Note <- function(rec, ...) any(vapply(rec$notes, function(note) all(vapply(c(...), function(part)
    grepl(part, note, fixed = TRUE), logical(1))), logical(1)))

k219Line <- function(row, rec, extra = "") cat(sprintf("KIT-219 %s\tnotes=%s%s\n", row, paste(rec$notes, collapse = " | "), extra))

# ------------------------------------------------------------------------------------------------
# E rows: failing first (ticket section 10).
# ------------------------------------------------------------------------------------------------

test_that("KIT-219 E1 (S2): a planted earlier run's files are left out of the workbook and the folder, and Copy says so", {
    fixture <- k219Journey()
    k219Enrich(fixture)
    product <- k219Names(fixture$contrast)
    planted <- c("C_vs_A_up_enrichment_results.tsv", "C_vs_A_up_enrichment_plot.png")
    k219Plant(fixture, planted)
    rec <- k219Copy(fixture)
    workbook <- k219Workbook(fixture$dirs, "proteomics")
    held <- k219Held(fixture$dirs, "proteomics")
    record <- k219Record(rec$returned)
    k219Line("E1", rec, sprintf("\tcontrast=%s\tsheets=%s\tpairs=%s\theld=%s\trecord=%s", fixture$contrast, workbook$n,
        paste(workbook$pairs, collapse = ","), paste(held, collapse = ","), k219Text(record$files)))
    expect_null(rec$value$error)
    expect_true(workbook$present)
    expect_identical(workbook$n, 2L)
    expect_setequal(workbook$pairs, paste(fixture$contrast, c("up", "down")))
    expect_setequal(held, product)
    expect_false(is.null(record))
    expect_setequal(record$files, planted)
    expect_identical(record$display_name, "Enrichment files not made by the current enrichment")
    expect_true(k219Note(rec, "[warning]", k219Texts$t1_lead, k219Texts$t1, planted[[1L]]), info = paste(rec$notes, collapse = " | "))
    expect_false("[message] Publication files copied" %in% rec$notes)
    expect_true(all(file.exists(file.path(fixture$pathway, c(planted, product)))))
    if (exists("workflowEnrichmentCurrentFiles", mode = "function")) {
        expect_setequal(workflowEnrichmentCurrentFiles(fixture$workflow, "proteomics"), product)
        expect_identical(workflowEnrichmentCurrentFiles(fixture$workflow, "metabolomics"), character())
        fake <- list2env(list(state_manager = WorkflowState$new(audit_enabled = FALSE)))
        expect_null(workflowEnrichmentCurrentFiles(fake, "proteomics"))
    }
})

test_that("KIT-219 E2 (S1): after the DA is run again the earlier enrichment is not copied and Copy names every file", {
    fixture <- k219Journey()
    k219Enrich(fixture)
    expect_false(is.null(diaEnrichmentProduct(fixture$manager)), info = "precondition: a current product")
    k219RealDa(fixture, q = 0.01)
    expect_null(diaEnrichmentProduct(fixture$manager), info = "precondition: the product ended")
    on_disk <- sort(list.files(fixture$pathway))
    rec <- k219Copy(fixture)
    workbook <- k219Workbook(fixture$dirs, "proteomics")
    held <- k219Held(fixture$dirs, "proteomics")
    record <- k219Record(rec$returned)
    k219Line("E2", rec, sprintf("\tsheets=%s\theld=%d\trecord=%s", workbook$n, length(held), k219Text(record$files)))
    expect_null(rec$value$error)
    expect_identical(workbook$n, 0L)
    expect_length(held, 0L)
    expect_false(is.null(record))
    expect_setequal(record$files, on_disk)
    expect_true(k219Note(rec, "[warning]", k219Texts$t1_lead, k219Texts$t1), info = paste(rec$notes, collapse = " | "))
})

test_that("KIT-219 E3 (S3): with no enrichment in the store, tables left in the enrichment folder are not copied", {
    fixture <- k219Journey()
    tables <- c("X_up_enrichment_results.tsv", "X_down_enrichment_results.tsv")
    k219Plant(fixture, tables)
    rec <- k219Copy(fixture)
    workbook <- k219Workbook(fixture$dirs, "proteomics")
    held <- k219Held(fixture$dirs, "proteomics")
    record <- k219Record(rec$returned)
    k219Line("E3", rec, sprintf("\tsheets=%s\theld=%d\trecord=%s", workbook$n, length(held), k219Text(record$files)))
    expect_null(rec$value$error)
    expect_identical(workbook$n, 0L)
    expect_length(held, 0L)
    expect_false(is.null(record))
    expect_setequal(record$files, tables)
})

test_that("KIT-219 E4 (S1): Copy, then the DA run again: the Report is refused before anything is rendered", {
    fixture <- k219Journey()
    k219Enrich(fixture)
    rec <- k219Copy(fixture)
    expect_null(rec$value$error)
    expect_setequal(k219Held(fixture$dirs, "proteomics"), k219Names(fixture$contrast))
    k219RealDa(fixture, q = 0.01)
    report <- k219Report(fixture)
    held <- k219Held(fixture$dirs, "proteomics")
    cat(sprintf("KIT-219 E4\trender=%d\tnotes=%s\theld=%s\tfile=%s\n", report$render, paste(report$notes, collapse = " | "),
        paste(held, collapse = ","), file.exists(report$report)))
    expect_identical(report$render, 0L)
    expect_identical(report$bundle, 0L)
    expect_true(any(startsWith(report$notes, "[error] Report generation failed:") & grepl(k219Texts$t2, report$notes,
        fixed = TRUE)), info = paste(report$notes, collapse = " | "))
    expect_false(file.exists(report$report))
})

test_that("KIT-219 E6 (S6): Copy, then the enrichment run again with other content, then the Report: refused", {
    fixture <- k219Journey()
    k219Enrich(fixture, content = "run1")
    rec <- k219Copy(fixture)
    expect_null(rec$value$error)
    k219Enrich(fixture, content = "run2")
    expect_setequal(k219Held(fixture$dirs, "proteomics"), k219Names(fixture$contrast))
    report <- k219Report(fixture)
    cat(sprintf("KIT-219 E6\trender=%d\tnotes=%s\n", report$render, paste(report$notes, collapse = " | ")))
    expect_identical(report$render, 0L)
    expect_true(any(startsWith(report$notes, "[error] Report generation failed:") & grepl(k219Texts$t2, report$notes,
        fixed = TRUE)), info = paste(report$notes, collapse = " | "))
    expect_false(file.exists(report$report))
})

#' A lipid or metab Summary project of the module-CI helpers whose folders hold no DA file (KIT-212's guard then lets
#' Copy through): the inputs written, the DA tables and the lipid table of significant counts removed.
k219OmicProject <- function(omic, pathway = TRUE) {
    metab <- identical(omic, "metabolomics")
    dirs <- if (metab) module_ci_metab_summary_paths() else module_ci_lipid_summary_paths()
    workflow <- if (metab) module_ci_metab_summary_workflow() else module_ci_lipid_summary_workflow()
    if (metab) module_ci_metab_summary_write_publication_inputs(dirs, include_optional = TRUE) else
        module_ci_lipid_summary_write_publication_inputs(dirs, include_optional = TRUE)
    paths <- dirs[[omic]]
    unlink(list.files(paths$da_output_dir, full.names = TRUE))
    unlink(file.path(paths$publication_graphs_dir, "NumSigDeMolecules", "lipids_num_sig_de_molecules.tab"))
    if (!pathway) unlink(list.files(paths$pathway_dir, full.names = TRUE))
    withr::defer(unlink(paths$base_dir, recursive = TRUE), envir = parent.frame())
    list(dirs = dirs, workflow = workflow)
}

test_that("KIT-219 E5 (OD-2 (i)): metabolomics has no enrichment analysis, so tables in its enrichment folder are not copied", {
    project <- k219OmicProject("metabolomics")
    pathway <- project$dirs$metabolomics$pathway_dir
    tables <- sort(list.files(pathway))
    expect_length(tables, 2L)
    rec <- k219CopyCore(project$workflow, "metabolomics", project$dirs)
    workbook <- k219Workbook(project$dirs, "metabolomics")
    held <- k219Held(project$dirs, "metabolomics")
    record <- k219Record(rec$returned)
    k219Line("E5", rec, sprintf("\tsheets=%s\theld=%d\trecord=%s", workbook$n, length(held), k219Text(record$files)))
    expect_null(rec$value$error)
    expect_identical(workbook$n, 0L)
    expect_length(held, 0L)
    expect_false(is.null(record))
    expect_setequal(record$files, tables)
    expect_true(k219Note(rec, "[warning]", k219Texts$t1_lead, k219Texts$t1b, "metabolomics"), info = paste(rec$notes, collapse = " | "))
    expect_identical(sort(list.files(pathway)), tables)
})

test_that("KIT-219 E9 (K219-I03): the decided loop end to end: refused, Copy again leaves the file out and says so, the Report renders", {
    fixture <- k219Journey()
    k219Enrich(fixture)
    product <- k219Names(fixture$contrast)
    stale <- "C_vs_A_up_enrichment_results.tsv"
    k219Plant(fixture, stale)
    # What an earlier Copy (an earlier build, or the folder copy) left in the publication folder: the stale file beside the current ones.
    folder <- k219Publication(fixture$dirs, "proteomics")
    dir.create(folder, recursive = TRUE)
    file.copy(file.path(fixture$pathway, c(stale, product)), folder)
    t2 <- paste("The report was not generated: the publication folder holds enrichment results that the current",
        "enrichment did not make. Press 'Copy to Publication Directory' again, then generate the report.")
    if (exists("assertWorkflowSummaryEnrichmentCurrent", mode = "function")) {
        refusal <- tryCatch(assertWorkflowSummaryEnrichmentCurrent(fixture$workflow, "proteomics", fixture$dirs$proteomics),
            error = function(e) e)
        expect_s3_class(refusal, "workflow_summary_enrichment_not_current")
        expect_identical(conditionMessage(refusal), t2)
    }
    refused <- k219Report(fixture)
    cat(sprintf("KIT-219 E9\trefused_render=%d\tnotes=%s\n", refused$render, paste(refused$notes, collapse = " | ")))
    expect_identical(refused$render, 0L)
    expect_identical(refused$bundle, 0L)
    expect_identical(refused$notes, paste("[error] Report generation failed:", t2))
    expect_false(file.exists(refused$report))
    rec <- k219Copy(fixture)
    expect_null(rec$value$error)
    expect_setequal(k219Held(fixture$dirs, "proteomics"), product)
    expect_true(file.exists(file.path(fixture$pathway, stale)))
    expect_true(k219Note(rec, "[warning]", k219Texts$t1_lead, "1 enrichment file(s)", k219Texts$t1, stale),
        info = paste(rec$notes, collapse = " | "))
    rendered <- k219Report(fixture)
    expect_identical(rendered$render, 1L)
    expect_true("[message] Report generated successfully!" %in% rendered$notes, info = paste(rendered$notes, collapse = " | "))
    expect_true(file.exists(rendered$report))
})

# ------------------------------------------------------------------------------------------------
# G rows: allowed, and unchanged.
# ------------------------------------------------------------------------------------------------

#' The publication folder as a comparison sees it: file names, sha256 of every file but the workbooks (their bytes carry
#' the time they were written), and each workbook's sheets with their contents.
k219Snapshot <- function(dirs, omic) {
    root <- dirs[[omic]]$results_summary_dir
    files <- sort(list.files(root, recursive = TRUE))
    hashes <- vapply(files, function(file) if (grepl("[.]xlsx$", file)) NA_character_ else
        digest::digest(file = file.path(root, file), algo = "sha256"), "")
    workbooks <- files[grepl("[.]xlsx$", files)]
    sheets <- lapply(stats::setNames(workbooks, workbooks), function(file) {
        names <- openxlsx::getSheetNames(file.path(root, file))
        stats::setNames(lapply(names, function(sheet) openxlsx::read.xlsx(file.path(root, file), sheet)), names)
    })
    list(files = files, hashes = hashes, sheets = sheets)
}

test_that("KIT-219 G1: with the current product's files only, Copy with the list equals the unchanged copy", {
    fixture <- k219Journey()
    k219Enrich(fixture)
    rec <- k219Copy(fixture)
    expect_null(rec$value$error)
    with_list <- k219Snapshot(fixture$dirs, "proteomics")
    k219CopyDirect(fixture)
    without <- k219Snapshot(fixture$dirs, "proteomics")
    cat(sprintf("KIT-219 G1\tfiles=%d\tnotes=%s\n", length(with_list$files), paste(rec$notes, collapse = " | ")))
    expect_true(length(with_list$files) > 0L)
    expect_identical(with_list$files, without$files)
    expect_identical(with_list$hashes, without$hashes)
    expect_identical(with_list$sheets, without$sheets)
    expect_null(k219Record(rec$returned))
    expect_true("[message] Publication files copied" %in% rec$notes, info = paste(rec$notes, collapse = " | "))
})

test_that("KIT-219 G2: copyToResultsSummary called as the workbooks call it copies every file of the enrichment folder", {
    fixture <- k219Journey()
    k219Enrich(fixture)
    planted <- c("C_vs_A_up_enrichment_results.tsv", "C_vs_A_up_enrichment_plot.png")
    k219Plant(fixture, planted)
    returned <- k219CopyDirect(fixture)
    workbook <- k219Workbook(fixture$dirs, "proteomics")
    held <- k219Held(fixture$dirs, "proteomics")
    cat(sprintf("KIT-219 G2\tsheets=%s\theld=%d\n", workbook$n, length(held)))
    expect_identical(workbook$n, 3L)
    expect_setequal(held, c(k219Names(fixture$contrast), planted))
    expect_null(k219Record(returned))
})

test_that("KIT-219 G3: lipidomics and metabolomics with an empty enrichment folder: no record, the plain notice, the workbook", {
    for (omic in c("lipidomics", "metabolomics")) {
        project <- k219OmicProject(omic, pathway = FALSE)
        expect_length(list.files(project$dirs[[omic]]$pathway_dir), 0L)
        rec <- k219CopyCore(project$workflow, omic, project$dirs)
        workbook <- k219Workbook(project$dirs, omic)
        k219Line(paste("G3", omic), rec, sprintf("\tsheets=%s", workbook$n))
        expect_null(rec$value$error, info = omic)
        expect_null(k219Record(rec$returned), info = omic)
        expect_true("[message] Publication files copied" %in% rec$notes, info = paste(omic, paste(rec$notes, collapse = " | ")))
        expect_true(workbook$present, info = omic)
        expect_identical(workbook$n, 0L, info = omic)
    }
})

# ------------------------------------------------------------------------------------------------
# Report templates (3.2, 3.1 item 5, G5): the sources, and the chunk that decides, run on a folder.
# ------------------------------------------------------------------------------------------------

k219TemplateDir <- function(omic) {
    dir <- system.file("reports", omic, package = "MultiScholaR")
    if (!nzchar(dir)) dir <- file.path(rprojroot::find_package_root_file(), "inst", "reports", omic)
    dir
}

k219Lines <- function(omic, template) readLines(file.path(k219TemplateDir(omic), template), warn = FALSE)

#' The code lines of the chunk labelled `label` (empty when there is none or more than one).
k219Chunk <- function(lines, label) {
    start <- grep(paste0("^```[{]r ", label, "[,}]"), lines)
    if (length(start) != 1L || start == length(lines)) return(character())
    close <- which(grepl("^```[[:space:]]*$", lines[(start + 1L):length(lines)]))
    if (!length(close)) return(character())
    lines[(start + 1L):(start + close[[1L]] - 1L)]
}

#' The header line of the chunk labelled `label`.
k219Header <- function(lines, label) {
    at <- grep(paste0("^```[{]r ", label, "[,}]"), lines)
    if (length(at) == 1L) lines[[at]] else ""
}

#' A chunk's code evaluated in a new environment; its printed text and `enrichment_ran`.
k219Run <- function(code, folder, text, pathway_dir = NULL) {
    env <- new.env(parent = globalenv())
    env$report_publication_figures_dir <- folder
    env$current_project_paths <- list(pathway_dir = pathway_dir)
    env$results_data <- list(content = list(pathway_enrichment = text, go_enrichment = text))
    printed <- utils::capture.output(eval(parse(text = code), envir = env))
    list(out = paste(printed, collapse = "\n"), ran = env$enrichment_ran)
}

test_that("KIT-219 E7a (OD-5 (i)): the lipidomics and metabolomics templates print the enrichment section only when files exist", {
    for (case in list(list("lipidomics", "lipidomics_report.rmd"), list("metabolomics", "metabolomics_report.rmd"))) {
        lines <- k219Lines(case[[1L]], case[[2L]])
        plots <- k219Chunk(lines, "Enrichment Plots")
        info <- case[[2L]]
        expect_true(any(grepl("enrichment_ran <- ", lines, fixed = TRUE)), info = info)
        expect_true(grepl("eval=enrichment_ran", k219Header(lines, "Enrichment Plots"), fixed = TRUE), info = info)
        expect_true(grepl("eval=enrichment_ran", k219Header(lines, "Enrichment Caption"), fixed = TRUE), info = info)
        expect_true(any(grepl("No enrichment analysis was run for this project", lines, fixed = TRUE)), info = info)
        expect_true(length(plots) > 0L && !any(grepl('cat(glue::glue("Warning:', plots, fixed = TRUE)), info = info)
    }
})

test_that("KIT-219 E7b (OD-5 (i)): the enrichment_presence chunk of each template, run on three publication folders", {
    for (case in list(list("lipidomics", "lipidomics_report.rmd"), list("metabolomics", "metabolomics_report.rmd"))) {
        code <- k219Chunk(k219Lines(case[[1L]], case[[2L]]), "enrichment_presence")
        info <- case[[2L]]
        expect_true(length(code) > 0L, info = paste(info, "chunk enrichment_presence exists"))
        if (!length(code)) next
        root <- withr::local_tempdir()
        absent <- k219Run(code, file.path(root, "none"), "PARAGRAPH")
        dir.create(file.path(root, "empty", "Enrichment_Plots"), recursive = TRUE)
        empty <- k219Run(code, file.path(root, "empty"), "PARAGRAPH")
        dir.create(file.path(root, "full", "Enrichment_Plots"), recursive = TRUE)
        writeLines("x", file.path(root, "full", "Enrichment_Plots", "x_up_enrichment_plot.png"))
        full <- k219Run(code, file.path(root, "full"), "PARAGRAPH")
        for (none in list(absent, empty)) {
            expect_false(none$ran, info = info)
            expect_true(grepl("# Pathway Enrichment Analysis", none$out, fixed = TRUE), info = info)
            expect_true(grepl("No enrichment analysis was run for this project.", none$out, fixed = TRUE), info = info)
            expect_false(grepl("PARAGRAPH", none$out, fixed = TRUE), info = info)
        }
        expect_true(full$ran, info = info)
        expect_true(grepl("PARAGRAPH", full$out, fixed = TRUE), info = info)
    }
})

test_that("KIT-219 E8 (3.1 item 5): the proteomics templates' flag reads the publication folder only", {
    for (template in c("LFQ_report.rmd", "TMT_report.rmd", "DIANN_report.rmd", "DIANN_limpa_report.rmd")) {
        code <- k219Chunk(k219Lines("proteomics", template), "enrichment_presence")
        expect_true(length(code) > 0L, info = template)
        if (!length(code)) next
        root <- withr::local_tempdir()
        pathway <- file.path(root, "pathway_enrichment")
        dir.create(pathway)
        writeLines("x", file.path(pathway, "x_up_enrichment_results.tsv"))
        dir.create(file.path(root, "empty", "Enrichment_Plots"), recursive = TRUE)
        empty <- k219Run(code, file.path(root, "empty"), "PARAGRAPH", pathway_dir = pathway)
        dir.create(file.path(root, "full", "Enrichment_Plots"), recursive = TRUE)
        writeLines("x", file.path(root, "full", "Enrichment_Plots", "x_up_enrichment_results.tsv"))
        full <- k219Run(code, file.path(root, "full"), "PARAGRAPH", pathway_dir = pathway)
        expect_false(empty$ran, info = template)
        expect_true(grepl("No enrichment analysis was run for this project.", empty$out, fixed = TRUE), info = template)
        expect_true(full$ran, info = template)
        expect_true(grepl("PARAGRAPH", full$out, fixed = TRUE), info = template)
    }
})

test_that("KIT-219 G5 (P07): the limpa template's STRING chunks read the publication folder, not the enrichment folder", {
    lines <- k219Lines("proteomics", "DIANN_limpa_report.rmd")
    for (label in c("stringdb_presence", "STRING DB Enrichment Plots")) {
        code <- k219Chunk(lines, label)
        expect_true(length(code) > 0L, info = label)
        defined <- code[grepl("stringdb_plots_dir <- ", code, fixed = TRUE)]
        expect_true(length(defined) > 0L, info = label)
        expect_true(all(grepl("report_publication_figures_dir", defined, fixed = TRUE)), info = label)
        expect_false(any(grepl("pathway_dir", code, fixed = TRUE)), info = label)
    }
})
