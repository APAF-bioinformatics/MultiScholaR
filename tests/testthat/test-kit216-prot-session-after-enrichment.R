# KIT-216 (VAL-175, VAL-179): proteomics session files after enrichment. The tests-only commit (no R/ change): each
# row is written to fail on the base at the value the ticket states (sections 3.0 and 10; the gate ids G1 to G10
# become E1 to E10) and to pass on the fix, or, for the G rows, to hold on both as stated. On the native DIA owner of
# helper-prot-norm-core.R (the "small" recipe, committed as `protein_replicate_filtered`) and its DIA_LIMPA owner.
# Every enrichment is committed through the application's commit, `saveWorkflowEnrichmentState` with the proteomics
# enrichment configuration (ticket P10); its new argument `source` is passed only where the function has it, and a new
# function (`diaEnrichmentProduct`) is called only where it exists, after the base assertions. Lines starting
# `KIT-216` are the evidence of docs/qa/workflow-execution/kit/kit-216/.
#
# Written on staging 16112917 (KIT-215 and KIT-212 merged). E9 (the base-made archives of section 4) is added in a
# second tests-only commit once R1 has made them; with `KIT216_ARCHIVE_OUT` set, the row "R1 archives" writes them.
# G12 runs on DIA only (no LFQ or TMT native owner with a DA among the fixtures); G13 is test-go-propagation.R's.

# --- Helpers copied from test-kit206-prot-norm-export-guard.R (source lines named) --------------------------------

# :14-16, with the normalisation method as an argument.
k216Settings <- function(mode = "skip", norm_method = "quantile") list(norm_method = norm_method, ruv_mode = mode,
    auto_percentage_min = NULL, auto_percentage_max = NULL, separation_metric = NULL, k_penalty_weight = NULL,
    adaptive_k_penalty = NULL, ruv_percentage = 50, ruv_k = 1, color_variable = "group", shape_variable = NULL)

k216Aesthetics <- function() list(color_var = "group", shape_var = NULL) # :18

k216Call <- function(fn, args) do.call(fn, args[intersect(names(args), names(formals(fn)))]) # :21

k216Scope <- function(fixture, code) withDIAStateScope(fixture$manager, code) # :23

k216Head <- function(fixture) getDIAStateReference(fixture$manager) # :25

k216Revert <- function(fixture, state) k216Scope(fixture, function() selectDIAStateReference(fixture$manager, state)) # :27

k216Error <- function(code) tryCatch({ # :29-32
    force(code)
    NULL
}, error = function(e) e)

#' The design's contrast on the recipe's groups (control, treated).
k216Contrasts <- data.frame(contrasts = "grouptreated-groupcontrol",
    full_format = "treated_vs_control=grouptreated-groupcontrol", friendly_names = "treated_vs_control",
    stringsAsFactors = FALSE)

#' :35-39, with the full contrasts table.
k216Owner <- function(envir = parent.frame()) {
    fixture <- protNormCoreOwner(envir)
    fixture$workflow$contrasts_tbl <- k216Contrasts
    fixture
}

#' :42-47, with the normalisation method as an argument.
k216Click <- function(fixture, norm, mode = "skip", norm_method = "quantile") {
    protNormCoreQuiet(k216Call(runDIAReferenceNormalization, list(settings = k216Settings(mode, norm_method),
        workflow = fixture$workflow, norm_data = norm, paths = fixture$paths, omic_type = "proteomics",
        experiment_label = "core", post_qc = function(reference) invisible(NULL),
        ruv_qc = function(reference) invisible(NULL), grouping = function() "group", pre_qc = NULL)))
}

#' :50-55.
k216Correlation <- function(fixture, norm, threshold = 0.5, skipped = FALSE) {
    final <- protNormCoreQuiet(runDIAReferenceCorrelation(fixture$workflow, norm, fixture$paths, threshold, skipped,
        "proteomics", "core", if (skipped) function() NULL else function() "group", k216Aesthetics))
    norm$correlation_filtering_complete <- TRUE
    final
}

#' :83-86 (the DA's re-commit of the head's arguments only; no manifest).
k216Da <- function(fixture) protNormCoreQuiet(k216Scope(fixture, function() commitDIADaArguments(fixture$manager,
    list(da_q_val_thresh = 0.05, treat_lfc_cutoff = 0, formula_string = "~ 0 + group"), fixture$workflow,
    data.frame(contrasts = "groupB-groupA", full_format = "B_vs_A=groupB-groupA", stringsAsFactors = FALSE),
    "Protein.Ids")))

#' :98-103.
k216PassThrough <- function(fixture, norm) protNormCoreQuiet(runDIAReferenceLimpaNormalisationPassThrough(
    k216Settings("covariates"), fixture$workflow, norm, function(reference) invisible(NULL)))

k216RuvCovariates <- function(fixture, norm, mode) protNormCoreQuiet(runDIAReferenceRuvCovariates(list(ruv_mode = mode,
    color_variable = "group", shape_variable = NULL), fixture$workflow, norm, fixture$paths, "proteomics", "core",
    function(reference) invisible(NULL), function() "group"))

#' :113-121.
k216Recorder <- function() {
    seen <- new.env(parent = emptyenv())
    seen$calls <- 0L
    seen$run <- function(workflowData, normData, input, sourceDir) {
        seen$calls <- seen$calls + 1L
        list(exportArtifacts = list(sessionFilename = "recorded"))
    }
    seen
}

#' :125-145, with the tab's normalisation selector as an argument.
k216Export <- function(fixture, norm, mode = "skip", norm_method = "quantile", recorder = NULL) {
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
        list(norm_method = norm_method, ruv_mode = mode), fixture$workflow, norm, fixture$paths,
        runExportSessionWorkflowFn = run, showNotificationFn = note, messageFn = quiet,
        handleExportErrorFn = function(e) {
            log$errors <- c(log$errors, list(e))
            handleProtNormExportError(e, showNotificationFn = note, messageFn = quiet)
        })))
    as.list(log)
}

k216Archives <- function(paths) sort(list.files(paths$source_dir, pattern = "^filtered_session_.*[.]msdia[.]tar$",
    full.names = TRUE)) # :147-148

#' :151-158.
k216Snapshot <- function(paths) {
    read <- function(name) {
        path <- file.path(paths$source_dir, name)
        if (file.exists(path)) readLines(path, warn = FALSE) else NULL
    }
    list(archives = k216Archives(paths), latest = read("filtered_session_latest.json"),
        summary = read("filtered_session_summary.txt"))
}

# --- New helpers ---------------------------------------------------------------------------------------------------

k216Text <- function(x) if (is.null(x) || !length(x)) "<none>" else paste(format(x), collapse = ",")

#' A DA on the head (`runDIAReferenceDaAnalysis`, as helper-go-propagation.R:151-156, without its publication): the
#' head's arguments re-committed and the accepted manifest on the head (`snapshot$da_results`).
k216RealDa <- function(fixture, q = 0.05) {
    da_data <- new.env(parent = emptyenv())
    da_data$current_s4_object <- k216Scope(fixture, function() getDIAStateReference(fixture$manager))
    da_data$contrasts_available <- k216Contrasts$contrasts
    protNormCoreQuiet(runDIAReferenceDaAnalysis(list(formula_string = "~ 0 + group", da_q_val_thresh = q,
        treat_lfc_cutoff = 0), da_data, fixture$workflow, fixture$paths, contrasts = k216Contrasts, publish = FALSE))
}

#' KIT-206's enrichment object (:91-93) on the recipe's contrast.
k216EnrichmentResults <- function() methods::new("EnrichmentResults",
    contrasts = tibble::tibble(contrasts = "grouptreated-groupcontrol"),
    enrichment_data = list(`grouptreated-groupcontrol` = list(up = data.frame(term_id = "GO:0008150", p_value = 0.01))),
    enrichment_plots = list(), enrichment_plotly = list(), enrichment_summaries = list())

#' KIT-215's `not_tested` row shape (`diaEnrichmentNoIdentifiers`), committed in the record's config.
k216NotTested <- data.frame(comparison = "treated_vs_control", contrast = "grouptreated-groupcontrol", direction = "up",
    proteins = 3L, route = "gprofiler", reason = "no_gene_name", stringsAsFactors = FALSE)

#' A completed enrichment through the application's commit (`saveWorkflowEnrichmentState`, the proteomics enrichment
#' configuration), its DA manifest's source passed where the function takes it. The commit's own result.
k216Enrich <- function(fixture, manifest = NULL, q_cutoff = 0.05, not_tested = NULL) {
    parameters <- list(selected_contrast = "treated_vs_control", up_cutoff = 0, down_cutoff = 0, q_cutoff = q_cutoff,
        organism_taxid = "9606")
    declared <- list(taxon_id = "9606", up_cutoff = 0, down_cutoff = 0, q_cutoff = q_cutoff, correction_method = "gSCS")
    saved <- protNormCoreQuiet(k216Scope(fixture, function() k216Call(saveWorkflowEnrichmentState, list(
        workflow = fixture$workflow, config = omicGuiConfig("proteomics")$enrichment, results = k216EnrichmentResults(),
        declared = declared, method = list(method = "gprofiler", supported = TRUE), parameters = parameters,
        pathway_dir = file.path(fixture$root, "pathway"), source = manifest$source, not_tested = not_tested))))
    expect_true(isTRUE(saved$saved), info = paste("enrichment commit:", k216Text(saved$warning)))
    invisible(saved)
}

#' One archive opened into a private folder: its recorded current state and statements, the stored state names and the
#' enrichment marker; or the error Open raised.
k216Open <- function(path) {
    root <- tempfile("k216-archive-")
    dir.create(root)
    imported <- tryCatch(protNormCoreQuiet(openDIANativeSession(path, root)), error = function(e) e)
    if (inherits(imported, "error")) return(list(error = imported))
    on.exit(closeDIAStateAdapter(imported$manager), add = TRUE)
    snapshot <- imported$manager$.context$snapshot
    list(error = NULL, state = imported$metadata$r6_current_state_name,
        normalization_method = imported$metadata$normalization_method, ruv_mode = imported$metadata$ruv_mode,
        ruv_applied = imported$metadata$ruv_applied, checkpoint_kind = imported$metadata$checkpoint_kind,
        output_policy = imported$metadata$output_policy,
        correlation = imported$extra_references$correlation_filtered$entity_id, stored = names(snapshot$states),
        marker = snapshot$enrichment_results)
}

#' KIT-206's `k206Reopen` (:178-199): the archive opened into a new owner and project by the package's own path
#' (`openDIANativeSession`, `adoptDIANativeSession`, as Load Filtered Session does), then the tab's restore with the
#' selector's updates recorded. The stage that raised, if any.
k216Reopen <- function(archive, envir = parent.frame()) {
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
    out <- list(error = NULL, stage = "none", selected = character(), root = root)
    imported <- tryCatch(protNormCoreQuiet(openDIANativeSession(archive, paths$source_dir)), error = function(e) e)
    if (inherits(imported, "error")) return(utils::modifyList(out, list(error = imported, stage = "open")))
    view <- NULL
    adopted <- k216Error(view <- protNormCoreQuiet(adoptDIANativeSession(imported, workflow, session, paths$source_dir,
        paths)))
    if (!is.null(adopted)) {
        try(closeDIAStateAdapter(imported$manager), silent = TRUE)
        return(utils::modifyList(out, list(error = adopted, stage = "adopt")))
    }
    withr::defer(closeDIAStateAdapter(workflow$state_manager), envir = envir)
    norm <- new.env(parent = emptyenv())
    updates <- character()
    testthat::local_mocked_bindings(updateSelectInput = function(session, inputId, label = NULL, choices = NULL,
        selected = NULL) {
        updates <<- c(updates, paste0(inputId, "=", k216Text(selected)))
        invisible(NULL)
    }, .package = "shiny")
    restored <- k216Error(protNormCoreQuiet(restoreProtNormTabState(NULL, workflow, norm)))
    out$error <- restored
    out$stage <- if (is.null(restored)) "none" else "restore"
    out$selected <- updates
    out$view <- view
    out$fixture <- list(manager = workflow$state_manager, workflow = workflow, paths = paths, root = root)
    out$norm <- norm
    out
}

k216Dirs <- function(fixture) list(proteomics = c(fixture$paths, list(integration_dir = file.path(fixture$root,
    "integration"), base_dir = fixture$root)))

k216Input <- function(fixture) file.path(fixture$root, "integration", "proteomics")

#' "Save Workflow Arguments" in process (KIT-206 E11b's harness, :637-653): the notices with their type and duration,
#' `done`, `failed`, and the parameters writer and the record writer recorded (each writes its file).
k216Save <- function(fixture, dirs = k216Dirs(fixture)) {
    seen <- new.env(parent = emptyenv())
    seen$notes <- character()
    seen$types <- character()
    seen$persistent <- logical()
    seen$failed <- list()
    seen$done <- list()
    seen$args <- list()
    seen$records <- 0L
    protNormCoreQuiet(saveWorkflowSummaryParameters(omicGuiConfig("proteomics")$summary, fixture$workflow, dirs,
        "proteomics", "k216", "d", notify = function(message, type = NULL, duration = NULL) {
            seen$notes <- c(seen$notes, as.character(message))
            seen$types <- c(seen$types, k216Text(type))
            seen$persistent <- c(seen$persistent, is.null(duration))
            invisible(NULL)
        }, done = function(saved) seen$done <- c(seen$done, list(saved)),
        failed = function(error) seen$failed <- c(seen$failed, list(error)),
        createWorkflowArgsFn = function(...) {
            args <- list(...)
            seen$args <- c(seen$args, list(args))
            path <- file.path(args$source_dir_path, "study_parameters.txt")
            writeLines("KIT-216 study parameters", path)
            path
        }, writeAnalysisRecordFn = function(manager, path) {
            seen$records <- seen$records + 1L
            writeLines("{}", path)
            invisible(path)
        }, catFn = function(...) NULL))
    as.list(seen)
}

k216SaveLine <- function(case, fixture, seen, extra = "") cat(sprintf(
    "KIT-216 %s head=%s notes=%s types=%s done=%d failed=%s args=%d records=%d input=%s%s\n", case,
    fixture$manager$current_state, paste(seen$notes, collapse = " | "), paste(seen$types, collapse = ","),
    length(seen$done), if (length(seen$failed)) class(seen$failed[[1L]])[[1L]] else "none", length(seen$args),
    seen$records, file.exists(file.path(k216Input(fixture), "integration_input.json")), extra))

#' "Generate Report" in process with a recorded render (KIT-207 E8's harness, test-kit207-summary-head-guard.R:509-526).
k216Report <- function(fixture, dirs = k216Dirs(fixture)) {
    calls <- new.env(parent = emptyenv())
    calls$render <- 0L
    calls$notes <- character()
    notify <- function(message, type, duration) calls$notes <- c(calls$notes, paste0("[", type, "] ", message))
    protNormCoreQuiet(utils::capture.output(renderWorkflowSummaryReport(omicGuiConfig("proteomics")$summary,
        fixture$workflow, dirs, "proteomics", "k216", notify = notify, progress = function(...) NULL,
        renderFn = function(omic_type, experiment_label, rmd_filename) {
            calls$render <- calls$render + 1L
            path <- file.path(fixture$root, "k216-report.html")
            writeLines("KIT-216 report", path)
            path
        }, exportRecordFn = function(...) invisible(NULL), downloadFileFn = function(url, destfile, ...) {
            writeLines("KIT-216 template", destfile)
            invisible(0L)
        }, catFn = function(...) NULL)))
    as.list(calls)
}

#' "Export Session State" in process (KIT-207 E8's harness, :530-536): the notices and the files written.
k216SessionState <- function(fixture, dirs = k216Dirs(fixture)) {
    notes <- character()
    protNormCoreQuiet(utils::capture.output(writeWorkflowSummarySessionState(omicGuiConfig("proteomics")$summary,
        fixture$workflow, dirs, "proteomics", "k216", "", values = list(workflow_args_saved = FALSE,
            files_copied = FALSE, report_generated = FALSE, report_path = NULL),
        notify = function(message, type, duration) notes <<- c(notes, paste0("[", type, "] ", message)))))
    list(notes = notes, files = list.files(fixture$paths$source_dir, "^session_state_.*[.]RDS$", full.names = TRUE))
}

#' The proteomics Summary guard's condition on the current head, or NULL when it passes.
k216Guard <- function(fixture) k216Error(assertWorkflowSummaryHeadFinal(fixture$workflow, "proteomics"))

#' The record's enrichment steps.
k216EnrichmentSteps <- function(fixture) {
    record <- buildWorkflowAnalysisRecord(fixture$manager)
    Filter(function(step) identical(step$operation, "enrichment"), record$steps)
}

#' Normalise (RUV skip), apply the correlation filter (0.5): the default path's table.
k216Normalised <- function(envir = parent.frame()) {
    fixture <- k216Owner(envir)
    norm <- new.env(parent = emptyenv())
    k216Click(fixture, norm, "skip")
    k216Correlation(fixture, norm, 0.5)
    list(fixture = fixture, norm = norm)
}

#' Where the fix's `diaEnrichmentProduct` exists: the head has no current enrichment product.
k216ExpectNoProduct <- function(fixture, case) {
    if (exists("diaEnrichmentProduct", mode = "function")) expect_null(diaEnrichmentProduct(fixture$manager), info = case)
}

# ------------------------------------------------------------------------------------------------
# Step 0's gates as failing-first rows (G1 to G10 become E1 to E10)
# ------------------------------------------------------------------------------------------------

test_that("KIT-216 E1 (G1, D1 (a)): after the DA and an enrichment the export writes the table with the product (fails before)", {
    q <- k216Normalised()
    fixture <- q$fixture
    manifest <- k216RealDa(fixture)
    k216Enrich(fixture, manifest)
    before <- k216Snapshot(fixture$paths)
    log <- k216Export(fixture, q$norm, "skip")
    added <- setdiff(k216Archives(fixture$paths), before$archives)
    archive <- if (length(added) == 1L) k216Open(added) else list()
    cat(sprintf("KIT-216 E1 head=%s notes=%s errors=%s new_archives=%d archive_state=%s stored=%s marker=%s\n",
        fixture$manager$current_state, paste(log$notes, collapse = " | "),
        k216Text(vapply(log$errors, function(e) class(e)[[1L]], "")), length(added), k216Text(archive$state),
        k216Text(archive$stored), !is.null(archive$marker)))
    expect_identical(fixture$manager$current_state, "correlation_filtered")
    expect_length(log$errors, 0L)
    expect_length(added, 1L)
    expect_null(archive$error)
    expect_identical(archive$state, "correlation_filtered")
    expect_true("enrichment_completed" %in% archive$stored)
    expect_false(is.null(archive$marker))
})

test_that("KIT-216 E2 (G2, D1 (b), SD-2): with no export in the session the Save writes the integration input and skips the session file (fails before)", {
    q <- k216Normalised()
    fixture <- q$fixture
    manifest <- k216RealDa(fixture)
    k216Enrich(fixture, manifest)
    input <- k216Input(fixture)
    dir.create(input, recursive = TRUE)
    writeLines("KIT-216 planted", file.path(input, "integration_input.json"))
    before <- k216Snapshot(fixture$paths)
    seen <- k216Save(fixture)
    planted <- identical(readLines(file.path(input, "integration_input.json"), warn = FALSE), "KIT-216 planted")
    added <- setdiff(k216Archives(fixture$paths), before$archives)
    k216SaveLine("E2", fixture, seen, sprintf(" planted_kept=%s new_archives=%d", planted, length(added)))
    expect_length(seen$failed, 0L)
    expect_length(seen$done, 1L)
    expect_false("Study parameters saved with warnings" %in% seen$notes)
    expect_true(file.exists(file.path(input, "integration_input.json")))
    expect_false(planted)
    expect_length(added, 0L)
    # T2 (VAL-41 draft): a persistent warning after the success notice.
    expect_true(any(seen$types == "warning" & seen$persistent))
})

test_that("KIT-216 E3 (G3, D1 (c)): the archive the Summary writes after an enrichment opens and is adopted at the table (fails before)", {
    q <- k216Normalised()
    fixture <- q$fixture
    exported <- k216Export(fixture, q$norm, "skip")
    expect_length(exported$errors, 0L)
    manifest <- k216RealDa(fixture)
    k216Enrich(fixture, manifest, not_tested = k216NotTested)
    before <- k216Snapshot(fixture$paths)
    seen <- k216Save(fixture)
    added <- setdiff(k216Archives(fixture$paths), before$archives)
    k216SaveLine("E3 save", fixture, seen, sprintf(" new_archives=%d", length(added)))
    expect_length(added, 1L)
    skip_if(length(added) != 1L, "no archive written by the Save")
    reopened <- k216Reopen(added)
    restored <- if (is.null(reopened$error)) restoreWorkflowEnrichment(reopened$fixture$manager,
        omicGuiConfig("proteomics")$enrichment)
    cat(sprintf("KIT-216 E3 reopen stage=%s error=%s head=%s restored_contrasts=%s\n", reopened$stage,
        if (is.null(reopened$error)) "none" else paste(class(reopened$error)[[1L]], conditionMessage(reopened$error)),
        k216Text(reopened$fixture$manager$current_state), k216Text(names(restored$all_enrichment_results))))
    expect_null(reopened$error)
    expect_identical(reopened$fixture$manager$current_state, "correlation_filtered")
    expect_false(is.null(restored))
    expect_true(length(restored$all_enrichment_results) > 0L)
    # KIT-215's `not_tested` from the product path (3.1 item 4).
    expect_identical(restored$not_tested, k216NotTested)
})

test_that("KIT-216 E4 (G4, D1 (e)): an archive exported after a Revert from an enrichment opens (fails before)", {
    q <- k216Normalised()
    fixture <- q$fixture
    manifest <- k216RealDa(fixture)
    k216Enrich(fixture, manifest)
    k216Revert(fixture, "correlation_filtered")
    before <- k216Snapshot(fixture$paths)
    log <- k216Export(fixture, q$norm, "skip")
    added <- setdiff(k216Archives(fixture$paths), before$archives)
    expect_length(log$errors, 0L)
    expect_length(added, 1L)
    skip_if(length(added) != 1L, "no archive written by the export")
    reopened <- k216Reopen(added)
    cat(sprintf("KIT-216 E4 export_notes=%s reopen stage=%s error=%s head=%s\n", paste(log$notes, collapse = " | "),
        reopened$stage, if (is.null(reopened$error)) "none" else paste(class(reopened$error)[[1L]],
            conditionMessage(reopened$error)), k216Text(reopened$fixture$manager$current_state)))
    expect_null(reopened$error)
    expect_identical(reopened$fixture$manager$current_state, "correlation_filtered")
})

test_that("KIT-216 E11 (G6, I05): DIA with limpa, the refreshed archive states the RUV mode of the last covariates run (fails before)", {
    skip_if_not_installed("limpa")
    skip_if_not_installed("mixOmics")
    fixture <- protNormW1LimpaOwner()
    fixture$manager$.context$session_format <- "native_v1"
    fixture$manager$.context$experiment_paths <- fixture$paths
    norm <- new.env(parent = emptyenv())
    k216PassThrough(fixture, norm)
    k216Correlation(fixture, norm, 0.5)
    k216RuvCovariates(fixture, norm, "covariates")
    expect_identical(fixture$manager$current_state, "ruv_covariates")
    before <- k216Snapshot(fixture$paths)
    log <- k216Export(fixture, norm, "covariates")
    added <- setdiff(k216Archives(fixture$paths), before$archives)
    feasible <- length(log$errors) == 0L && length(added) == 1L
    cat(sprintf("KIT-216 E11 feasibility export_notes=%s new_archives=%d feasible=%s\n", paste(log$notes,
        collapse = " | "), length(added), feasible))
    # Stop rule 4: the native limpa store needs a package change; recorded under "Not checked", no claim.
    skip_if(!feasible, "G6 feasibility: the limpa store cannot export natively without a package change")
    k216RuvCovariates(fixture, norm, "correction")
    k216Da(fixture)
    before <- k216Snapshot(fixture$paths)
    error <- k216Error(protNormCoreQuiet(refreshDIANativeSummaryArchive(fixture$workflow, fixture$paths$source_dir)))
    added <- setdiff(k216Archives(fixture$paths), before$archives)
    archive <- if (length(added) == 1L) k216Open(added) else list()
    cat(sprintf("KIT-216 E11 head=%s refresh_error=%s new_archives=%d ruv_mode=%s ruv_applied=%s\n",
        fixture$manager$current_state, if (is.null(error)) "none" else conditionMessage(error), length(added),
        k216Text(archive$ruv_mode), k216Text(archive$ruv_applied)))
    expect_null(error)
    expect_length(added, 1L)
    expect_identical(archive$ruv_mode, "correction")
    expect_identical(archive$ruv_applied, TRUE)
})

test_that("KIT-216 E7 (G7, I04): the archive and the reopened tab state the normalisation method that was run (fails before)", {
    q <- k216Normalised()
    fixture <- q$fixture
    before <- k216Snapshot(fixture$paths)
    log <- k216Export(fixture, q$norm, "skip", norm_method = "cyclicloess")
    added <- setdiff(k216Archives(fixture$paths), before$archives)
    expect_length(added, 1L)
    skip_if(length(added) != 1L, "no archive written by the export")
    archive <- k216Open(added)
    reopened <- k216Reopen(added)
    cat(sprintf("KIT-216 E7 normalization_method=%s reopen_error=%s selector=%s\n", k216Text(archive$normalization_method),
        if (is.null(reopened$error)) "none" else conditionMessage(reopened$error), k216Text(reopened$selected)))
    expect_identical(archive$normalization_method, "quantile")
    expect_null(reopened$error)
    expect_true("norm_method=quantile" %in% reopened$selected)
    expect_false("norm_method=cyclicloess" %in% reopened$selected)
})

test_that("KIT-216 E7b (I04, P08): after two runs the archive states the later run's method (fails before)", {
    fixture <- k216Owner()
    norm <- new.env(parent = emptyenv())
    k216Click(fixture, norm, "manual", "quantile")
    k216Click(fixture, norm, "skip", "cyclicloess")
    k216Correlation(fixture, norm, 0.5)
    cat(sprintf("KIT-216 E7b history=%s\n", k216Text(utils::tail(diaStateHistory(fixture$manager), 3L))))
    before <- k216Snapshot(fixture$paths)
    log <- k216Export(fixture, norm, "skip", norm_method = "quantile")
    added <- setdiff(k216Archives(fixture$paths), before$archives)
    expect_length(added, 1L)
    skip_if(length(added) != 1L, "no archive written by the export")
    archive <- k216Open(added)
    cat(sprintf("KIT-216 E7b normalization_method=%s\n", k216Text(archive$normalization_method)))
    expect_identical(archive$normalization_method, "cyclicloess")
})

test_that("KIT-216 E8 (G8, D2): after a second enrichment run the Save writes the integration input and the record's input is the table (fails before)", {
    q <- k216Normalised()
    fixture <- q$fixture
    exported <- k216Export(fixture, q$norm, "skip")
    expect_length(exported$errors, 0L)
    manifest <- k216RealDa(fixture)
    k216Enrich(fixture, manifest)
    k216Enrich(fixture, manifest)
    seen <- k216Save(fixture)
    steps <- k216EnrichmentSteps(fixture)
    inputs <- lapply(steps, function(step) unlist(step$input_entity_ids))
    own <- fixture$manager$.context$snapshot$states$enrichment_completed$entity_id
    k216SaveLine("E8", fixture, seen, sprintf(" enrichment_steps=%d step_input=%s own_entity=%s table=%s",
        length(steps), k216Text(unlist(inputs)), k216Text(own), k216Text(manifest$source$entity_id)))
    expect_true(file.exists(file.path(k216Input(fixture), "integration_input.json")))
    expect_length(steps, 1L)
    expect_identical(unlist(inputs), manifest$source$entity_id)
})

# E9 (G9) and E10 (G10): the stale head of KIT-206 E11 (an export at the correlation result, the DA's re-commit, a
# Revert to `protein_replicate_filtered`).
k216StaleHead <- function(envir = parent.frame()) {
    q <- k216Normalised(envir)
    exported <- k216Export(q$fixture, q$norm, "skip")
    expect_length(exported$errors, 0L)
    k216Da(q$fixture)
    k216Revert(q$fixture, "protein_replicate_filtered")
    q
}

test_that("KIT-216 E9 (G9, I10): at a stale head Generate Report and Export Session State are refused before anything is written (fails before)", {
    q <- k216StaleHead()
    fixture <- q$fixture
    guard <- k216Guard(fixture)
    report <- k216Report(fixture)
    state <- k216SessionState(fixture)
    cat(sprintf("KIT-216 E9 head=%s guard=%s report_render=%d report_notes=%s state_files=%d state_notes=%s\n",
        fixture$manager$current_state, if (is.null(guard)) "pass" else class(guard)[[1L]], report$render,
        paste(report$notes, collapse = " | "), length(state$files), paste(state$notes, collapse = " | ")))
    expect_s3_class(guard, "omic_norm_correlation_input_stale")
    expect_identical(report$render, 0L)
    expect_length(state$files, 0L)
    if (!is.null(guard)) {
        expect_true(paste("[error] Report generation failed:", cli::ansi_strip(conditionMessage(guard))) %in% report$notes)
        expect_true(paste("[error] Export failed:", conditionMessage(guard)) %in% state$notes)
    }
})

test_that("KIT-216 E10 (G10, I10): at a stale head the Save is refused before the parameters and the record are written (fails before)", {
    q <- k216StaleHead()
    fixture <- q$fixture
    input <- k216Input(fixture)
    dir.create(input, recursive = TRUE)
    writeLines("KIT-216 planted", file.path(input, "integration_input.json"))
    seen <- k216Save(fixture)
    written <- file.exists(file.path(fixture$paths$source_dir, c("study_parameters.txt", "analysis_record.json")))
    k216SaveLine("E10", fixture, seen, sprintf(" written=%s", k216Text(written)))
    expect_length(seen$failed, 1L)
    if (length(seen$failed)) expect_s3_class(seen$failed[[1L]], "omic_norm_correlation_input_stale")
    expect_length(seen$args, 0L)
    expect_identical(seen$records, 0L)
    expect_identical(written, c(FALSE, FALSE))
    expect_false(dir.exists(input))
})

# ------------------------------------------------------------------------------------------------
# Section 10's further rows
# ------------------------------------------------------------------------------------------------

test_that("KIT-216 E5: a DA run again with other arguments ends the product (head only; base printed)", {
    q <- k216Normalised()
    fixture <- q$fixture
    manifest <- k216RealDa(fixture)
    k216Enrich(fixture, manifest)
    again <- k216Error(k216RealDa(fixture, q = 0.01))
    steps <- k216EnrichmentSteps(fixture)
    cat(sprintf("KIT-216 E5 head=%s second_da=%s enrichment_steps=%d marker=%s\n", fixture$manager$current_state,
        if (is.null(again)) "ran" else conditionMessage(again), length(steps),
        !is.null(fixture$manager$.context$snapshot$enrichment_results)))
    expect_null(again)
    expect_length(steps, 0L)
    k216ExpectNoProduct(fixture, "E5")
})

test_that("KIT-216 E5b (P03): an enrichment of an earlier table's DA after a Revert is not a product of the head (fails before)", {
    q <- k216Normalised()
    fixture <- q$fixture
    manifest <- k216RealDa(fixture)
    k216Revert(fixture, "normalized")
    k216Enrich(fixture, manifest)
    steps <- k216EnrichmentSteps(fixture)
    cat(sprintf("KIT-216 E5b head=%s enrichment_steps=%d\n", fixture$manager$current_state, length(steps)))
    expect_identical(fixture$manager$current_state, "normalized")
    expect_length(steps, 0L)
    k216ExpectNoProduct(fixture, "E5b")
})

test_that("KIT-216 E6: a Revert to an earlier table after the DA and an enrichment ends the product (head only; base printed)", {
    q <- k216Normalised()
    fixture <- q$fixture
    manifest <- k216RealDa(fixture)
    k216Enrich(fixture, manifest)
    reverted <- k216Error(k216Revert(fixture, "normalized"))
    cat(sprintf("KIT-216 E6 revert=%s head=%s marker=%s\n", if (is.null(reverted)) "ran" else conditionMessage(reverted),
        fixture$manager$current_state, !is.null(fixture$manager$.context$snapshot$enrichment_results)))
    expect_null(reverted)
    expect_null(fixture$manager$.context$snapshot$enrichment_results)
    k216ExpectNoProduct(fixture, "E6")
})

test_that("KIT-216 E10c (SD-1): a project never normalised is refused Save, Report and Session State (fails before)", {
    fixture <- k216Owner()
    expect_identical(fixture$manager$current_state, "protein_replicate_filtered")
    guard <- k216Guard(fixture)
    seen <- k216Save(fixture)
    report <- k216Report(fixture)
    state <- k216SessionState(fixture)
    k216SaveLine("E10c", fixture, seen, sprintf(" guard=%s report_render=%d state_files=%d",
        if (is.null(guard)) "pass" else conditionMessage(guard), report$render, length(state$files)))
    expect_s3_class(guard, "omic_norm_correlation_input_stale")
    expect_length(seen$failed, 1L)
    expect_length(seen$args, 0L)
    expect_identical(report$render, 0L)
    expect_length(state$files, 0L)
})

test_that("KIT-216 E14 (P06): at `normalized` the Save names the `normalized` table as final (fails before)", {
    fixture <- k216Owner()
    norm <- new.env(parent = emptyenv())
    k216Click(fixture, norm, "skip")
    expect_identical(fixture$manager$current_state, "normalized")
    k216RealDa(fixture)
    seen <- k216Save(fixture)
    final <- if (length(seen$args)) seen$args[[1L]]$final_s4_object$entity_id
    ids <- vapply(c("normalized", "protein_replicate_filtered"), function(state)
        k216Text(getDIAStateReference(fixture$manager, state)$entity_id), "")
    k216SaveLine("E14", fixture, seen, sprintf(" final=%s normalized=%s protein_replicate_filtered=%s", k216Text(final),
        ids[["normalized"]], ids[["protein_replicate_filtered"]]))
    expect_identical(final, getDIAStateReference(fixture$manager, "normalized")$entity_id)
})

test_that("KIT-216 E9b: an old-shape archive whose current state is the enrichment opens at the DA's table (fails before)", {
    q <- k216Normalised()
    fixture <- q$fixture
    exported <- k216Export(fixture, q$norm, "skip")
    expect_length(exported$errors, 0L)
    manifest <- k216RealDa(fixture)
    # The store shape before the fix, built on purpose: the enrichment committed as the head state (KIT-206 :94-95).
    built <- k216Error(protNormCoreQuiet(fixture$manager$saveState("enrichment_completed", k216EnrichmentResults(),
        list(q_cutoff = 0.05), "fixture enrichment")))
    expect_null(built)
    expect_identical(fixture$manager$current_state, "enrichment_completed")
    before <- k216Snapshot(fixture$paths)
    error <- k216Error(protNormCoreQuiet(refreshDIANativeSummaryArchive(fixture$workflow, fixture$paths$source_dir)))
    added <- setdiff(k216Archives(fixture$paths), before$archives)
    expect_null(error)
    expect_length(added, 1L)
    skip_if(length(added) != 1L, "no archive written by the refresh")
    reopened <- k216Reopen(added)
    cat(sprintf("KIT-216 E9b reopen stage=%s error=%s head=%s marker=%s\n", reopened$stage,
        if (is.null(reopened$error)) "none" else paste(class(reopened$error)[[1L]], conditionMessage(reopened$error)),
        k216Text(reopened$fixture$manager$current_state),
        !is.null(reopened$fixture$manager$.context$snapshot$enrichment_results)))
    expect_null(reopened$error)
    expect_identical(reopened$fixture$manager$current_state, "correlation_filtered")
    expect_false(is.null(reopened$fixture$manager$.context$snapshot$enrichment_results))
    if (exists("diaEnrichmentProduct", mode = "function")) {
        expect_false(is.null(diaEnrichmentProduct(reopened$fixture$manager)))
    }
})

# ------------------------------------------------------------------------------------------------
# Guard rows (G15, G16, G18, G19, G20)
# ------------------------------------------------------------------------------------------------

test_that("KIT-216 G15: the Summary saves after the DA, after an enrichment, and on DIA with limpa after the covariates (guard)", {
    q <- k216Normalised()
    fixture <- q$fixture
    exported <- k216Export(fixture, q$norm, "skip")
    expect_length(exported$errors, 0L)
    manifest <- k216RealDa(fixture)
    seen <- k216Save(fixture)
    k216SaveLine("G15 da", fixture, seen)
    expect_length(seen$done, 1L)
    expect_length(seen$failed, 0L)
    k216Enrich(fixture, manifest)
    seen <- k216Save(fixture)
    k216SaveLine("G15 enrichment", fixture, seen)
    expect_length(seen$done, 1L)
    expect_length(seen$failed, 0L)
    skip_if_not_installed("limpa")
    skip_if_not_installed("mixOmics")
    limpa <- protNormW1LimpaOwner()
    norm <- new.env(parent = emptyenv())
    k216PassThrough(limpa, norm)
    k216Correlation(limpa, norm, 0.5)
    k216RuvCovariates(limpa, norm, "covariates")
    seen <- k216Save(limpa)
    k216SaveLine("G15 limpa", limpa, seen)
    expect_length(seen$done, 1L)
    expect_length(seen$failed, 0L)
})

test_that("KIT-216 G16: a second enrichment run with other cut-offs replaces the first; the record lists one step with the second run's parameters (guard)", {
    q <- k216Normalised()
    fixture <- q$fixture
    manifest <- k216RealDa(fixture)
    k216Enrich(fixture, manifest, q_cutoff = 0.05)
    k216Enrich(fixture, manifest, q_cutoff = 0.01)
    steps <- k216EnrichmentSteps(fixture)
    # R1 instrument repair (VAL-54; run 45588146575): a record value is encoded (`encodeWorkflowRecordValue`); it is
    # decoded before it is compared.
    values <- if (length(steps)) stats::setNames(lapply(steps[[1L]]$parameters, function(p)
        decodeWorkflowRecordValue(p$value)), vapply(steps[[1L]]$parameters, `[[`, "", "name"))
    cat(sprintf("KIT-216 G16 enrichment_steps=%d q_cutoff=%s\n", length(steps), k216Text(values$q_cutoff)))
    expect_length(steps, 1L)
    expect_equal(values$q_cutoff, 0.01)
})

test_that("KIT-216 G18 (SD-4): after an export and a new normalisation the Report runs and the Save writes all but the session file (fails before at the Save)", {
    q <- k216Normalised()
    fixture <- q$fixture
    exported <- k216Export(fixture, q$norm, "skip")
    expect_length(exported$errors, 0L)
    k216Click(fixture, q$norm, "skip", "cyclicloess")
    k216Correlation(fixture, q$norm, 0.5)
    k216RealDa(fixture)
    report <- k216Report(fixture)
    before <- k216Snapshot(fixture$paths)
    seen <- k216Save(fixture)
    added <- setdiff(k216Archives(fixture$paths), before$archives)
    k216SaveLine("G18", fixture, seen, sprintf(" report_render=%d new_archives=%d", report$render, length(added)))
    expect_identical(report$render, 1L)
    expect_length(seen$failed, 0L)
    expect_length(seen$done, 1L)
    expect_true(file.exists(file.path(k216Input(fixture), "integration_input.json")))
    expect_length(added, 0L)
    expect_true(any(seen$types == "warning" & seen$persistent))
})

test_that("KIT-216 G19 (P05, R1-03): after a Revert to `normalized` the Save completes; with an export, all but the session file (base printed)", {
    q <- k216Normalised()
    fixture <- q$fixture
    exported <- k216Export(fixture, q$norm, "skip")
    expect_length(exported$errors, 0L)
    k216Revert(fixture, "normalized")
    # R1 (run 45588146575): on this recipe the correlation filter at 0.5 removes no protein, so the export may hold the
    # same table as `normalized`; T2b is expected exactly when KIT-206's guard finds that it does not (R1-03).
    stale <- !is.null(k216Error(assertProtNormExportHeadCurrent(fixture$manager,
        fixture$manager$.context$owned_export$extra_references$correlation_filtered)))
    seen <- k216Save(fixture)
    k216SaveLine("G19", fixture, seen, sprintf(" export_stale=%s", stale))
    expect_length(seen$failed, 0L)
    expect_length(seen$done, 1L)
    expect_true(file.exists(file.path(k216Input(fixture), "integration_input.json")))
    expect_identical(any(seen$types == "warning" & seen$persistent), stale)
})

test_that("KIT-216 G20 (R1-03): S3, a Revert into an earlier run's `normalized` with an export: the Save writes all but the session file (fails before)", {
    fixture <- k216Owner()
    norm <- new.env(parent = emptyenv())
    k216Click(fixture, norm, "skip")
    k216Click(fixture, norm, "manual")
    k216Correlation(fixture, norm, 0.5)
    exported <- k216Export(fixture, norm, "manual")
    expect_length(exported$errors, 0L)
    k216Revert(fixture, "normalized")
    seen <- k216Save(fixture)
    k216SaveLine("G20", fixture, seen)
    expect_length(seen$failed, 0L)
    expect_length(seen$done, 1L)
    expect_true(file.exists(file.path(k216Input(fixture), "integration_input.json")))
    expect_true(any(seen$types == "warning" & seen$persistent))
})

# ------------------------------------------------------------------------------------------------
# E15 (3.1 item 13; KIT-212 F3): the enrichment Run on a DA that is no longer current
# ------------------------------------------------------------------------------------------------

test_that("KIT-216 E15 (item 13): the enrichment Run is refused when the DA tab has cleared its result (fails before)", {
    q <- k216Normalised()
    fixture <- q$fixture
    manifest <- k216RealDa(fixture)
    tab <- shiny::reactiveValues(da_results_data = manifest, current_s4_object = manifest$source)
    fixture$workflow$da_analysis_results_list <- manifest
    k216Revert(fixture, "normalized")
    da_data <- new.env(parent = emptyenv())
    da_data$da_results_list <- manifest
    k216Scope(fixture, function() invalidateDIAVisibleResults(da_data, fixture$workflow))
    # R1 instrument repair (VAL-54; run 45588146575): on this recipe the correlation filter at 0.5 removes no protein, so
    # the DA tab's check finds the same table and keeps the bus's manifest; the DA tab's clearing is applied as it
    # applies it (ticket E15: "or set to NULL as it does").
    cleared_by_tab <- is.null(fixture$workflow$da_analysis_results_list)
    fixture$workflow$da_analysis_results_list <- NULL
    cleared <- is.null(fixture$workflow$da_analysis_results_list)
    calls <- 0L
    notes <- character()
    testthat::local_mocked_bindings(runWorkflowEnrichment = function(...) {
        calls <<- calls + 1L
        list(saved = list(saved = TRUE), method = list(method = "gprofiler", supported = TRUE))
    }, holdOmicEnrichRun = function(tab, run) invisible(run), .package = "MultiScholaR")
    testthat::local_mocked_bindings(showNotification = function(ui, action = NULL, duration = 5, closeButton = TRUE,
        id = NULL, type = c("default", "message", "warning", "error"), session = NULL) {
        notes <<- c(notes, sprintf("[%s] %s", type[[1L]], paste(as.character(ui), collapse = " ")))
        invisible(id %||% "notification")
    }, withProgress = function(expr, ...) force(expr), incProgress = function(...) invisible(NULL),
        removeNotification = function(...) invisible(NULL), .package = "shiny")
    input <- list(selected_contrast = "treated_vs_control", organism_taxid = "9606", up_cutoff = 0, down_cutoff = 0,
        q_cutoff = 0.05, correction_method = "gSCS", enable_organism_filter = FALSE)
    out <- k216Error(protNormCoreQuiet(shiny::isolate(runOmicEnrichAnalysis(input, tab, fixture$workflow,
        fixture$paths, omicGuiConfig("proteomics")$enrichment))))
    cat(sprintf("KIT-216 E15 head=%s cleared_by_tab=%s bus_cleared=%s run_calls=%d error=%s notes=%s\n",
        fixture$manager$current_state, cleared_by_tab, cleared, calls, if (is.null(out)) "none" else conditionMessage(out), paste(notes, collapse = " | ")))
    expect_true(cleared)
    expect_identical(calls, 0L)
    expect_true(any(startsWith(notes, "[warning]")))
})

# ------------------------------------------------------------------------------------------------
# E7c, G11, G13b, and the base-made archives of section 4 (R1)
# ------------------------------------------------------------------------------------------------

test_that("KIT-216 E7c (I04, P08): DIA with limpa, the archive states the protein normalisation step's method (fails before)", {
    skip_if_not_installed("limpa")
    skip_if_not_installed("mixOmics")
    fixture <- protNormW1LimpaOwner()
    fixture$manager$.context$session_format <- "native_v1"
    norm <- new.env(parent = emptyenv())
    k216PassThrough(fixture, norm)
    k216Correlation(fixture, norm, 0.5)
    k216RuvCovariates(fixture, norm, "covariates")
    before <- k216Snapshot(fixture$paths)
    log <- k216Export(fixture, norm, "covariates", norm_method = "quantile")
    added <- setdiff(k216Archives(fixture$paths), before$archives)
    archive <- if (length(added) == 1L) k216Open(added) else list()
    cat(sprintf("KIT-216 E7c export_notes=%s new_archives=%d normalization_method=%s\n", paste(log$notes,
        collapse = " | "), length(added), k216Text(archive$normalization_method)))
    skip_if(length(added) != 1L, "G6 feasibility: the limpa store cannot export natively without a package change")
    expect_identical(archive$normalization_method, "none")
})

#' The record with store-local ids, timestamps and builds masked, one line per step: what G11 compares base against
#' head (the R2 logs).
k216RecordLines <- function(record) {
    mask <- function(x) {
        if (is.character(x)) {
            x[!is.na(x) & grepl("^[a-f0-9]{16,64}$", x)] <- "<id>"
            x[!is.na(x) & grepl("^[0-9]{4}-[0-9]{2}-[0-9]{2}T", x)] <- "<time>"
            return(x)
        }
        if (is.list(x)) x[] <- lapply(x, mask)
        x
    }
    vapply(seq_along(record$steps), function(i) {
        step <- mask(record$steps[[i]])
        sprintf("step %d %s %s fields=%s sha=%s", i, step$operation, step$state, paste(names(step), collapse = ","),
            digest::digest(step[setdiff(names(step), c("committed_at", "software"))], algo = "sha256"))
    }, "")
}

test_that("KIT-216 G11: the default path's record keeps its steps and fields (printed for base against head)", {
    q <- k216Normalised()
    fixture <- q$fixture
    exported <- k216Export(fixture, q$norm, "skip")
    expect_length(exported$errors, 0L)
    manifest <- k216RealDa(fixture)
    k216Enrich(fixture, manifest)
    record <- buildWorkflowAnalysisRecord(fixture$manager)
    lines <- k216RecordLines(record)
    cat(paste0("KIT-216 G11 ", lines, "\n"), sep = "")
    cat(sprintf("KIT-216 G11 top-level fields=%s design_state=%s\n", paste(names(record), collapse = ","),
        k216Text(record$design$state)))
    steps <- Filter(function(step) identical(step$operation, "enrichment"), record$steps)
    expect_length(steps, 1L)
    expect_identical(record$steps[[length(record$steps)]]$operation, "enrichment")
    expect_true(all(c("committed_at", "software", "defaults_sha256") %in% names(steps[[1L]])))
})

#' A DIA project as test-workflow-analysis-record.R:979-1016 builds it (the module CI's DA fixture, the DA, its
#' publication, then `runs` enrichment runs as the tab runs them on the test double), with its record written.
k216ReplayProject <- function(root, runs = 1L) {
    manager <- newDIAStateAdapter(newWorkflowEntityStore(root, "k216-replay", "DIA"), WorkflowState$new(audit_enabled = FALSE))
    on.exit(closeDIAStateAdapter(manager), add = TRUE)
    manager$.context$reference_states <- TRUE
    manager$setWorkflowType("DIA")
    design <- module_ci_prot_da_design()
    rownames(design) <- NULL
    quant <- module_ci_prot_da_table()
    rownames(quant) <- NULL
    setDIAStepProvenance(manager, "design_builder")
    manager$saveState("correlation_filtered", module_ci_prot_da_object(design = design, protein_quant_table = quant), list(),
        "fixture")
    paths <- list(source_dir = file.path(root, "source"), results_dir = file.path(root, "results"),
        da_output_dir = file.path(root, "results", "da_proteins"), pathway_dir = file.path(root, "results", "pathway_enrichment"),
        publication_graphs_dir = file.path(root, "publication_graphs"))
    for (p in paths) dir.create(p, recursive = TRUE, showWarnings = FALSE)
    manager$.context$experiment_paths <- paths
    contrasts <- module_ci_prot_da_contrasts("two_group")
    workflow <- list2env(list(state_manager = manager, contrasts_tbl = contrasts, tab_status = list(),
        qc_params = list(protein_qc = list()), config_list = list()))
    withDIAStateScope(manager, function() {
        da_data <- new.env(parent = emptyenv())
        da_data$current_s4_object <- getDIAStateReference(manager)
        da_data$contrasts_available <- contrasts$contrasts
        manifest <- protNormCoreQuiet(runDIAReferenceDaAnalysis(list(formula_string = "~ 0 + group",
            da_q_val_thresh = 0.05, treat_lfc_cutoff = 0), da_data, workflow, paths, contrasts = contrasts, publish = FALSE))
        protNormCoreQuiet(produceDIADaPublication(manager, manifest, data.frame(), 0.05,
            paths[c("da_output_dir", "publication_graphs_dir")], save_assets = FALSE))
        for (run in seq_len(runs)) utils::capture.output(protNormCoreQuiet(runWorkflowEnrichment(workflow,
            omicGuiConfig("proteomics")$enrichment, manifest, list(selected_contrast = "B_vs_A", organism_taxid = "9606",
                up_cutoff = 0, down_cutoff = 0, q_cutoff = c(0.05, 0.01)[[run]], correction_method = "gSCS",
                enable_organism_filter = FALSE), paths)))
    })
    json <- file.path(root, "analysis_record.json")
    writeWorkflowAnalysisRecord(manager, json)
    list(record = jsonlite::read_json(json, simplifyVector = FALSE), paths = paths)
}

test_that("KIT-216 G13b (P01): records made with one and with two enrichment runs replay to completion", {
    withr::local_envvar(MULTISCHOLAR_TEST_MODE = "true")
    for (runs in 1:2) {
        root <- withr::local_tempdir()
        project <- k216ReplayProject(root, runs)
        copy <- withr::local_tempdir()
        stopifnot(system2("cp", c("-a", file.path(root, "."), copy)) == 0L)
        h <- new.env(parent = globalenv())
        sys.source(file.path(testthat::test_path("..", ".."), "tools", "profiling", "workflow_kit.R"), h)
        api <- function(name) get(name, asNamespace("MultiScholaR"))
        prepared <- h$kitPrepareReplayProject(copy, c("initial", "correlation_filtered"), api)
        product <- "<not reached>"
        replay <- k216Error(utils::capture.output(protNormCoreQuiet(replayWorkflowAnalysisRecord(project$record, prepared,
            character(), project$paths, withr::local_tempdir(), on_step = function(step, manager) {
                if (identical(step$operation, "enrichment")) product <<- if (exists("diaEnrichmentProduct",
                    mode = "function")) !is.null(diaEnrichmentProduct(manager)) else "<function absent>"
            }))))
        cat(sprintf("KIT-216 G13b runs=%d enrichment_steps=%d replay=%s product=%s\n", runs,
            sum(vapply(project$record$steps, function(step) identical(step$operation, "enrichment"), logical(1))),
            if (is.null(replay)) "completed" else paste(class(replay)[[1L]], conditionMessage(replay)), format(product)))
        expect_null(replay, info = paste("runs", runs))
        if (exists("diaEnrichmentProduct", mode = "function")) expect_identical(product, TRUE, info = paste("runs", runs))
    }
})

test_that("KIT-216 R1 archives (section 4): the latest archive the Summary writes after one and two enrichment runs", {
    out <- Sys.getenv("KIT216_ARCHIVE_OUT")
    skip_if(!nzchar(out), "made only in R1, with KIT216_ARCHIVE_OUT set")
    dir.create(out, recursive = TRUE, showWarnings = FALSE)
    for (runs in 1:2) {
        q <- k216Normalised()
        fixture <- q$fixture
        exported <- k216Export(fixture, q$norm, "skip")
        expect_length(exported$errors, 0L)
        manifest <- k216RealDa(fixture)
        for (run in seq_len(runs)) k216Enrich(fixture, manifest, q_cutoff = c(0.05, 0.01)[[run]])
        before <- k216Snapshot(fixture$paths)
        error <- k216Error(protNormCoreQuiet(refreshDIANativeSummaryArchive(fixture$workflow, fixture$paths$source_dir)))
        added <- setdiff(k216Archives(fixture$paths), before$archives)
        expect_null(error)
        expect_length(added, 1L)
        if (length(added) == 1L) {
            target <- file.path(out, sprintf("base-enrichment-%d-run%s.msdia.tar", runs, if (runs > 1L) "s" else ""))
            file.copy(added, target, overwrite = TRUE)
            cat(sprintf("KIT-216 R1 archive runs=%d head=%s bytes=%d sha256=%s da_source=%s\n", runs,
                fixture$manager$current_state, file.size(target), digest::digest(file = target, algo = "sha256"),
                manifest$source$entity_id))
        }
    }
})

# ------------------------------------------------------------------------------------------------
# E9 (section 4, P07): the archives the base wrote in R1 (run 45588146575), opened on this commit
# ------------------------------------------------------------------------------------------------

test_that("KIT-216 E9 (P07): archives written before the fix with the enrichment as current state open at the DA's table (fails before)", {
    for (name in c("base-enrichment-1-run.msdia.tar", "base-enrichment-2-runs.msdia.tar")) {
        archive <- file.path(withr::local_tempdir(), name)
        file.copy(testthat::test_path("..", "testdata", "kit216", name), archive)
        reopened <- k216Reopen(archive)
        manager <- reopened$fixture$manager
        steps <- if (is.null(reopened$error)) Filter(function(step) identical(step$operation, "enrichment"),
            buildWorkflowAnalysisRecord(manager)$steps)
        table <- if (!is.null(manager)) getDIACurrentReference(manager)$entity_id
        cat(sprintf("KIT-216 E9 %s stage=%s error=%s head=%s marker=%s enrichment_steps=%d input_is_table=%s\n", name,
            reopened$stage, if (is.null(reopened$error)) "none" else paste(class(reopened$error)[[1L]],
                conditionMessage(reopened$error)), k216Text(manager$current_state),
            !is.null(manager$.context$snapshot$enrichment_results), length(steps),
            length(steps) == 1L && identical(unlist(steps[[1L]]$input_entity_ids), table)))
        expect_null(reopened$error, info = name)
        expect_true(is.character(reopened$view$enrichment_conversion), info = paste(name, "T1"))
        expect_identical(manager$current_state, "correlation_filtered", info = name)
        expect_false(is.null(manager$.context$snapshot$enrichment_results), info = name)
        if (exists("diaEnrichmentProduct", mode = "function")) expect_false(is.null(diaEnrichmentProduct(manager)), info = name)
        expect_length(steps, 1L)
        if (length(steps) == 1L) expect_identical(unlist(steps[[1L]]$input_entity_ids), table, info = name)
        # r1.4 (audit F02, (c)): the enrichment tab's restore on the converted owner returns the enrichment.
        restored <- if (is.null(reopened$error)) restoreWorkflowEnrichment(manager, omicGuiConfig("proteomics")$enrichment)
        cat(sprintf("KIT-216 E9 %s restored_contrasts=%s\n", name, k216Text(names(restored$all_enrichment_results))))
        expect_true(length(restored$all_enrichment_results) > 0L, info = paste(name, "restore"))
    }
})

# ------------------------------------------------------------------------------------------------
# G21 (r1.4; audit F02, F03): an enrichment below a table head is never restored as current
# ------------------------------------------------------------------------------------------------

test_that("KIT-216 G21 (F02): an old-shape store with the enrichment below a table head restores no enrichment (guard)", {
    q <- k216Normalised()
    fixture <- q$fixture
    manifest <- k216RealDa(fixture)
    # The store shape of a base project: the enrichment committed as the head state (KIT-206 :94-95), then the
    # normalisation and the correlation filter committed again on top of it (their names are already in the history,
    # so the history still ends with `enrichment_completed`).
    built <- k216Error(protNormCoreQuiet(fixture$manager$saveState("enrichment_completed", k216EnrichmentResults(),
        list(q_cutoff = 0.05), "fixture enrichment")))
    expect_null(built)
    for (state in c("normalized", "correlation_filtered")) protNormCoreQuiet(k216Scope(fixture, function() {
        record <- readDIAStateMetadata(fixture$manager, getDIAStateReference(fixture$manager, state))
        commitDIAReferenceState(fixture$manager, state, record$data, record$config, "fixture re-run")
    }))
    history <- diaStateHistory(fixture$manager)
    cat(sprintf("KIT-216 G21 head=%s history=%s\n", fixture$manager$current_state, k216Text(history)))
    expect_identical(fixture$manager$current_state, "correlation_filtered")
    expect_true("enrichment_completed" %in% history)
    before <- k216Snapshot(fixture$paths)
    log <- k216Export(fixture, q$norm, "skip")
    added <- setdiff(k216Archives(fixture$paths), before$archives)
    expect_length(added, 1L)
    skip_if(length(added) != 1L, "no archive written by the export")
    reopened <- k216Reopen(added)
    manager <- reopened$fixture$manager
    restored <- if (is.null(reopened$error)) restoreWorkflowEnrichment(manager, omicGuiConfig("proteomics")$enrichment)
    steps <- if (is.null(reopened$error)) Filter(function(step) identical(step$operation, "enrichment"),
        buildWorkflowAnalysisRecord(manager)$steps)
    cat(sprintf("KIT-216 G21 reopen stage=%s error=%s head=%s restored=%s enrichment_steps=%d inputs=%s\n", reopened$stage,
        if (is.null(reopened$error)) "none" else conditionMessage(reopened$error), k216Text(manager$current_state),
        !is.null(restored), length(steps), k216Text(unlist(lapply(steps, function(step) step$input_entity_ids)))))
    expect_null(restored)
})
