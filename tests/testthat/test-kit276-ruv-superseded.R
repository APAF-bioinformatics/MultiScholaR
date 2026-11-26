# KIT-276 (VAL-334, VAL-346, VAL-359; AUDIT-KIT-276-plan-20261126; revision notes 1.4 and 1.5; rulings of 2026-11-27 in
# docs/qa/workflow-execution/kit/kit-276/step0.md): DIA with limpa, what the tab, the analysis record, the session archive
# and the report say about RUV agrees with what the differential analysis uses. After a correlation re-apply that follows a
# RUV run the history, the record, the report and the archive carry no RUV step, and "Skip RUV" after a covariates or
# correction run takes the RUV result out of the analysis (the head goes back to `correlation_filtered`).
# The fixture and the presses are test-kit273-limpa-ruv-delivery.R's and test-kit244-limpa-pass-through.R's (the limpa owner
# on protNormW1LimpaOwner's recipe; the correlation and the RUV in process and through the worker), the replay is
# test-kit248-correlation-skip-record.R's, the report chunks test-kit227-report-ruv-sentences.R's, the enrichment recipe
# test-kit213-study-parameters-labels.R's; all copied under the prefix k276 and never sourced. The texts are drafts for the
# owner's wording pass (VAL-41); the tests assert them as drafts. Each row that fails on the old code fails BY VALUE.

k276Quiet <- function(code) suppressWarnings(suppressMessages(code))
k276Noop <- function(...) invisible(NULL)
k276Aesthetics <- function() list(color_var = "group", shape_var = NULL)
k276Paths <- function() c("in process", if (workflowOwnerCapable()) "worker")
k276Need <- function() {
    skip_if_not_installed("limpa")
    skip_if_not_installed("mixOmics")
}

# ---- the texts (drafts, VAL-41) ----

k276TextSkip <- paste0("RUV-III skipped. The earlier RUV-III result is set aside: the differential analysis will use the ",
    "correlation-filtered table without RUV factors or correction.")
k276TextRunDa <- " Run the differential analysis again."
k276TextAgain <- paste0("The correlation filter was applied again, so the table no longer carries the RUV-III result. ",
    "Compute the RUV-III covariates again to use them.")
k276TextBusy <- "A step is still running. Press \"Compute RUV-III Covariates\" again when it has finished."
k276TextOldSkip <- "RUV-III covariates skipped by the user; the DA model carries no W terms."

# ---- fixtures ----

k276Cache <- new.env(parent = emptyenv())
k276Quantified <- function() {
    if (is.null(k276Cache$quantified)) {
        precursors <- makeWorkflowCheckpointPeptide(readWorkflowCheckpointRecipe("small"))
        rolled <- rollUpPrecursorToPeptide(precursors, core_utilisation = NA)
        k276Cache$quantified <- k276Quiet(proteinMissingValueImputationLimpa(
            log2TransformPeptideMatrix(calcPeptideMatrix(rolled)), dpc_slope = 0.8, verbose = FALSE, chunk = 1000,
            engine = "limpa"))
    }
    k276Cache$quantified
}

k276HeadRef <- function(fixture) withDIAStateScope(fixture$manager, function() getDIAStateReference(fixture$manager))

#' A DIA with limpa owner (protNormW1LimpaOwner's commits), with the DA's and the enrichment's folders and the tab's DA state.
k276Owner <- function(envir = parent.frame()) {
    withr::local_options(multischolar.dia_native_sessions = TRUE, multischolar.capture_test_checkpoints = FALSE,
        .local_envir = envir)
    root <- withr::local_tempdir(.local_envir = envir)
    paths <- list(source_dir = file.path(root, "scripts", "proteomics"),
        protein_qc_dir = file.path(root, "results", "proteomics", "protein_qc"),
        peptide_qc_dir = file.path(root, "results", "proteomics", "peptide_qc"),
        da_output_dir = file.path(root, "results", "proteomics", "da_proteins"),
        pathway_dir = file.path(root, "results", "proteomics", "pathway_enrichment"))
    paths$publication_graphs_dir <- file.path(root, "results", "proteomics", "publication_graphs")
    paths$time_dir <- file.path(paths$publication_graphs_dir, "filtering_qc")
    for (path in unlist(paths)) dir.create(path, recursive = TRUE, showWarnings = FALSE)
    manager <- newDIAStateAdapter(newWorkflowEntityStore(paths$source_dir, "k276-limpa", "DIA"),
        WorkflowState$new(audit_enabled = FALSE))
    withr::defer(closeDIAStateAdapter(manager), envir = envir)
    manager$.context$reference_states <- TRUE
    manager$setWorkflowType("DIA_LIMPA")
    manager$.context$experiment_paths <- paths
    workflow <- list2env(list(state_manager = manager, qc_params = list(protein_qc = list()), config_list = list(),
        tab_status = list()))
    fixture <- list(root = root, paths = paths, workflow = workflow, manager = manager,
        da = new.env(parent = emptyenv()), seen = new.env(parent = emptyenv()))
    setDIAStepProvenance(manager, "design_builder")
    k276Quiet(manager$saveState("protein_dpc_quant", k276Quantified(), list(), "k276"))
    k276Quiet(withDIAStateScope(manager, function() {
        v1 <- computeDIAReferenceProteinNormalisation(manager, getDIAStateReference(manager), list(norm_method = "none"))
        setDIAStepProvenance(manager, "protein_between_sample_normalisation", "whole", "norm_method", NULL,
            values = list(norm_method = "none"))
        finishDIAReferenceProteinStep(workflow, v1, publish = FALSE)
        v2 <- computeDIAReferenceProteinPeptideCountFilter(manager, getDIAStateReference(manager),
            list(num_peptides_per_protein_thresh = 1, num_peptidoforms_per_protein_thresh = 2))
        setDIAStepProvenance(manager, "protein_peptide_count_filter", "whole",
            c("num_peptides_per_protein_thresh", "num_peptidoforms_per_protein_thresh"), NULL, values = v2$config[1:2])
        finishDIAReferenceProteinStep(workflow, v2, publish = FALSE)
    }))
    fixture
}

k276W2Mocks <- function(env = parent.frame()) {
    testthat::local_mocked_bindings(workflowExecutionBuildIdentity = function() strrep("f", 64),
        finishWorkflowExecutionJob = function(owner, status) { owner$job <- NULL; list(status = status) },
        .package = "MultiScholaR", .env = env)
}

#' The limpa owner on a path; on the worker path the W2 mocks in the test's frame and a controller.
k276Start <- function(path, envir = parent.frame()) {
    fixture <- k276Owner(envir = envir)
    if (identical(path, "worker")) {
        k276W2Mocks(env = envir)
        fixture$seen$controller <- omicsWorkerController(list(owner = fixture$manager))
    }
    fixture
}

k276Error <- function(code) tryCatch({
    k276Quiet(code)
    NA_character_
}, error = function(e) paste0(class(e)[[1L]], ": ", conditionMessage(e)))

k276Object <- function(fixture, reference) {
    manager <- fixture$manager
    withDIAStateScope(manager, function() {
        object <- resolveDIAStateValue(manager$.context, reference)$data
        if (!isS4(object)) object <- rebuildDIACheckpointObject(manager$.context, object, TRUE, list())
        object
    })
}

# ---- the tab's actions (kit273's and kit244's, copied) ----

k276Press <- function(fixture, norm) {
    settings <- list(norm_method = "quantile", ruv_mode = "covariates", auto_percentage_min = NULL,
        auto_percentage_max = NULL, separation_metric = NULL, k_penalty_weight = NULL, adaptive_k_penalty = NULL,
        ruv_percentage = 50, ruv_k = 1, color_variable = "group", shape_variable = NULL)
    k276Error(runDIAReferenceLimpaNormalisationPassThrough(settings, fixture$workflow, norm, k276Noop))
}

k276Dispatch <- function(fixture) function(operation, parameters, input, state_name = NULL) {
    context <- if (!is.null(state_name)) list(state_name = state_name, config = list(), description = "dispatch")
    protNormW2Dispatch(fixture, fixture$seen$controller, operation, parameters, input, context)
}

k276WorkerCorrelation <- function(fixture, norm, threshold) {
    manager <- fixture$manager
    source <- withDIAStateScope(manager, function() {
        assertDIAExportReady(manager, allow_action = "correlation")
        source <- diaCorrelationSource(manager, norm)
        readDIAStateMetadata(manager, source)
        source
    })
    dispatch <- k276Dispatch(fixture)
    pairs <- dispatch("pearson-pairs", diaPearsonPairsParameters("group"), source)
    selection <- dispatch("pearson-selection", list(pairs = pairs$reference$entity_id, threshold = threshold), source,
        "correlation_filtered")
    withDIAStateScope(manager, function() publishDIAReferenceCorrelation(manager, list(pairs = pairs$reference,
        final = selection$value), fixture$workflow, norm, fixture$paths, threshold, FALSE, "proteomics", "core",
        k276Aesthetics))
}

#' Correlation Apply at `threshold`, below the tab (the presses of the journeys' history).
k276Correlate <- function(fixture, norm, path, threshold = 0.5) k276Error(if (identical(path, "worker"))
    k276WorkerCorrelation(fixture, norm, threshold) else runDIAReferenceCorrelation(fixture$workflow, norm, fixture$paths,
    threshold, FALSE, "proteomics", "core", function() "group", k276Aesthetics))

k276RuvCovariates <- function(fixture, norm, path, mode = "covariates") {
    manager <- fixture$manager
    settings <- list(ruv_mode = mode, color_variable = "group", shape_variable = NULL)
    if (!identical(path, "worker")) return(k276Error(runDIAReferenceRuvCovariates(settings, fixture$workflow, norm,
        fixture$paths, "proteomics", "core", k276Noop, function() "group")))
    k276Error({
        request <- withDIAStateScope(manager, function() {
            assertDIAExportReady(manager, allow_action = "ruv-covariates")
            list(source = diaRuvCovariatesSource(manager, norm),
                parameters = resolveDIARuvCovariatesParameters(fixture$workflow$config_list, "group", mode))
        })
        result <- k276Dispatch(fixture)("protein_ruv_covariates", request$parameters, request$source,
            if (identical(mode, "correction")) "ruv_corrected" else "ruv_covariates")
        withDIAStateScope(manager, function() publishDIAReferenceRuvCovariates(manager, settings,
            result[c("value", "reference")], fixture$workflow, norm, fixture$paths, k276Noop))
    })
}

#' The presses as letters: P pass-through, C correlation, K covariates, R correction (RUV-III correction mode).
k276Journey <- function(fixture, norm, path, presses) {
    for (press in presses) {
        error <- switch(press, P = k276Press(fixture, norm), C = k276Correlate(fixture, norm, path),
            K = k276RuvCovariates(fixture, norm, path, "covariates"), R = k276RuvCovariates(fixture, norm, path, "correction"))
        expect_identical(error, NA_character_, info = paste("press", press, "of", paste(presses, collapse = "")))
    }
    invisible(fixture)
}

#' "Skip RUV" through the tab's observer; its summary text, its notices and the in-application builds of Figure 1.
k276RuvSkip <- function(fixture, norm, env = parent.frame()) {
    seen <- new.env(parent = emptyenv())
    seen$notes <- character()
    seen$builds <- list()
    note <- function(ui, type = "default", duration = 5, ...) {
        seen$notes <- c(seen$notes, paste(as.character(ui), collapse = " "))
        invisible(NULL)
    }
    testthat::local_mocked_bindings(buildDiaLimpaCombinedQcFigure = function(paths, omic_type, ruv_mode, manager)
        seen$builds <- c(seen$builds, list(ruv_mode)), .package = "MultiScholaR", .env = env)
    output <- new.env(parent = emptyenv())
    norm$correlation_filtering_complete <- TRUE
    k276Quiet(runProtNormRuvCovariatesObserver(list(ruv_mode = "skip", color_variable = "group", shape_variable = NULL),
        output, fixture$workflow, norm, fixture$paths, "proteomics", "core", k276Noop, function() "group",
        renderTextFn = function(x) x, showNotificationFn = note, messageFn = k276Noop))
    list(text = output$ruv_covariates_summary, notes = seen$notes, builds = seen$builds)
}

#' The completion as the tab's observers call it, with the summary text and the notices readable.
k276Complete <- function(output) function(...) completeProtNormCorrelationWorkflow(..., showNotificationFn = k276Noop,
    renderTextFn = function(x) x, messageFn = k276Noop)

#' Correlation Apply or Skip through the tab's in-process observers, inside the owner's scope.
k276Reapply <- function(fixture, norm, skipped = FALSE, env = parent.frame()) {
    seen <- new.env(parent = emptyenv())
    seen$builds <- list()
    testthat::local_mocked_bindings(buildDiaLimpaCombinedQcFigure = function(paths, omic_type, ruv_mode, manager)
        seen$builds <- c(seen$builds, list(ruv_mode)), .package = "MultiScholaR", .env = env)
    output <- new.env(parent = emptyenv())
    fail <- function(error, ...) stop(error)
    k276Quiet(withDIAStateScope(fixture$manager, function() if (skipped) runProtNormSkipCorrelationObserver(output,
        fixture$workflow, norm, fixture$paths, "proteomics", "core", k276Aesthetics, k276Noop,
        completeCorrelationWorkflowFn = k276Complete(output), handleCorrelationErrorFn = fail) else
        runProtNormApplyCorrelationObserver(list(min_pearson_correlation_threshold = 0.5), output, fixture$workflow, norm,
            fixture$paths, "proteomics", "core", function() "group", k276Aesthetics, k276Noop,
            completeCorrelationWorkflowFn = k276Complete(output), handleCorrelationErrorFn = fail)))
    list(output = output, builds = seen$builds)
}

# ---- the DA, the record, the history ----

k276Contrasts <- data.frame(contrasts = "treated.vs.control", full_format = "treated.vs.control=grouptreated-groupcontrol",
    friendly_names = "treated.vs.control", stringsAsFactors = FALSE)

#' The tab's DA click on the head, in process, no publication (kit248Da's).
k276Da <- function(fixture, publish = FALSE) {
    fixture$workflow$contrasts_tbl <- k276Contrasts
    fixture$da$contrasts_available <- k276Contrasts$contrasts
    fixture$da$current_s4_object <- k276HeadRef(fixture)
    k276Quiet(withDIAStateScope(fixture$manager, function() runDIAReferenceDaAnalysis(list(formula_string = "~ 0 + group",
        da_q_val_thresh = 0.05, treat_lfc_cutoff = 0), fixture$da, fixture$workflow, fixture$paths,
        contrasts = k276Contrasts, publish = publish)))
}

#' The DA's own model record (W terms and RUV mode it used).
k276Model <- function(fixture) {
    manager <- fixture$manager
    withDIAStateScope(manager, function() resolveWorkflowStateValue(manager$.context,
        manager$.context$snapshot$da_results))$model
}

k276History <- function(fixture) diaStateHistory(fixture$manager)

#' The analysis record, written and read back.
k276Record <- function(fixture) {
    json <- file.path(fixture$root, paste0("analysis_record_", length(list.files(fixture$root, "^analysis_record_")), ".json"))
    k276Quiet(writeWorkflowAnalysisRecord(fixture$manager, json))
    jsonlite::read_json(json, simplifyVector = FALSE)
}
k276StepStates <- function(record) vapply(record$steps, function(step) step$state, "")
k276StepOperations <- function(record) vapply(record$steps, function(step) step$operation, "")

# ---- the report's chunks (kit227's, for the limpa template only) ----

k276Template <- function() system.file("reports", "proteomics", "DIANN_limpa_report.rmd", package = "MultiScholaR")
k276Exprs <- function(path, label) {
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
k276ReportEnv <- function() {
    env <- new.env(parent = globalenv())
    env$params <- list(workflow_name = "maxquant", omic_type = "proteomics", experiment_label = "kit276")
    env$d66_log <- function(...) invisible(NULL)
    env$tibble <- tibble::tibble
    for (label in c("New Hierarchical Parser", "Method Section", "Result Section", "Sample Analysis Functions",
        "Main Execution Function")) {
        for (e in k276Exprs(k276Template(), label)) {
            if (is.call(e) && identical(e[[1L]], as.name("<-")) && is.call(e[[3L]]) &&
                identical(e[[3L]][[1L]], as.name("function"))) eval(e, env)
        }
    }
    env
}
k276ResultsCall <- function(env) {
    find <- function(expr) {
        if (!is.call(expr)) return(NULL)
        if (!is.null(names(expr)) && all(c("results", "method_details") %in% names(expr))) return(expr[["results"]])
        for (i in seq_along(expr)) if (is.call(expr[[i]])) {
            hit <- find(expr[[i]])
            if (!is.null(hit)) return(hit)
        }
        NULL
    }
    call <- find(body(env$Main))
    stopifnot(!is.null(call))
    call
}
k276Param <- function(name, value) list(name = name, value = list(type = typeof(value), value = as.list(value),
    length = length(value)))
k276Step <- function(operation, outcomes = list()) list(operation = operation, parameters = list(), outcomes = outcomes)

#' The report's Methods and Results text on the RUV outcomes of a record (the dpc and RUV steps' outcomes; the template's own
#' parser, AnalyzeSamples and its Results call).
k276ReportText <- function(record) {
    kept <- Filter(function(step) step$operation %in% c("protein_dpc_quant", "protein_ruv_covariates"), record$steps)
    steps <- lapply(kept, function(step) k276Step(step$operation, step$outcomes))
    env <- k276ReportEnv()
    dir <- withr::local_tempdir(pattern = "kit276-")
    env$report_source_dir <- dir
    params_path <- file.path(dir, "study_parameters.txt")
    writeLines(c("Study Parameters", "================", "", "Basic Information:", "-----------------",
        "* Workflow Name: kit276", ""), params_path, useBytes = TRUE)
    jsonlite::write_json(list(steps = steps), file.path(dir, "analysis_record.json"), auto_unbox = TRUE, null = "null",
        digits = NA)
    study <- env$ReadStudyParamsNew(params_path)
    numbers <- env$AnalyzeSamples(NULL, study)$numbers
    stopifnot(length(numbers) > 0L)
    results_env <- new.env(parent = env)
    results_env$analysis_results <- list(numbers = numbers, table = data.frame())
    results_env$correlation_data <- list()
    list(method = paste(as.character(env$FormatMethodDetails(numbers)), collapse = ""),
        results = paste(as.character(eval(k276ResultsCall(env), results_env)), collapse = ""))
}
k276Has <- function(text, pattern) grepl(pattern, text, fixed = TRUE)

# ---- the replay (kit248's, with the replay journey's own comparison) ----

k276Kit <- function() {
    if (is.null(k276Cache$kit)) {
        kit <- new.env(parent = globalenv())
        sys.source(file.path(testthat::test_path("..", ".."), "tools", "profiling", "workflow_kit.R"), kit)
        k276Cache$kit <- kit
    }
    k276Cache$kit
}

#' Replay the fixture's record on a copy of its store cut to the record's import and design states. `list(error, states,
#' mismatches, listed)`: the replay's states, the replay journey's own comparison over every state of the ORIGIN store
#' (`mismatches`) and the same comparison over the states of the origin's history only (`listed`).
k276Replay <- function(fixture, record, envir = parent.frame()) {
    copy <- withr::local_tempdir(.local_envir = envir)
    stopifnot(system2("cp", c("-a", file.path(fixture$paths$source_dir, "."), copy)) == 0L)
    states <- vapply(record$steps, function(step) step$state, "")
    operations <- vapply(record$steps, function(step) step$operation, "")
    preconditions <- states[seq_len(sum(cumprod(operations %in% c("import", "design_builder"))))]
    kit <- k276Kit()
    api <- function(name) get(name, asNamespace("MultiScholaR"))
    prepared <- kit$kitPrepareReplayProject(copy, preconditions, api)
    root <- withr::local_tempdir(.local_envir = envir)
    paths <- list(source_dir = file.path(root, "scripts", "proteomics"),
        protein_qc_dir = file.path(root, "results", "proteomics", "protein_qc"),
        peptide_qc_dir = file.path(root, "results", "proteomics", "peptide_qc"),
        da_output_dir = file.path(root, "results", "proteomics", "da"))
    paths$publication_graphs_dir <- file.path(root, "results", "proteomics", "publication_graphs")
    paths$time_dir <- file.path(paths$publication_graphs_dir, "filtering_qc")
    for (path in unlist(paths)) dir.create(path, recursive = TRUE, showWarnings = FALSE)
    journey <- tryCatch(k276Quiet(kit$kitReplayJourney(record, prepared, fixture$paths$source_dir, character(), paths,
        withr::local_tempdir(.local_envir = envir), api)), error = function(e) e)
    if (inherits(journey, "error")) return(list(error = journey, states = NULL, mismatches = NULL, listed = NULL))
    history <- k276History(fixture)
    populations <- function(fingerprints) lapply(fingerprints$states, function(x) list(fingerprint = x))
    origin <- populations(journey$origin)
    listed <- as.character(kit$kitReproduceMismatches(origin[intersect(history, names(origin))],
        populations(journey$replay)))
    list(error = NULL, states = journey$states, mismatches = as.character(journey$mismatches), listed = listed)
}
k276ErrorClass <- function(error) if (is.null(error)) NA_character_ else class(error)[[1L]]

# ---- rows ----

test_that("KIT-276 T1: pass-through, correlation, RUV covariates, correlation again leaves no RUV step; the DA used none", {
    k276Need()
    for (path in k276Paths()) {
        fixture <- k276Start(path)
        norm <- new.env(parent = emptyenv())
        k276Journey(fixture, norm, path, c("P", "C", "K", "C"))
        history <- k276History(fixture)
        expect_identical(utils::tail(history, 1L), "correlation_filtered", info = path)
        expect_false("ruv_covariates" %in% history, info = path)
        record <- k276Record(fixture)
        expect_false("ruv_covariates" %in% k276StepStates(record), info = path)
        expect_false("protein_ruv_covariates" %in% k276StepOperations(record), info = path)
        k276Da(fixture)
        model <- k276Model(fixture)
        expect_identical(model$w_terms, character(0), info = path)
        expect_identical(model$ruv_mode, "none", info = path)
    }
})

test_that("KIT-276 T2: the report's Methods and Results on T1's record state no RUV (the record's steps are T1's)", {
    k276Need()
    skip_if_not_installed("knitr")
    fixture <- k276Start("in process")
    norm <- new.env(parent = emptyenv())
    k276Journey(fixture, norm, "in process", c("P", "C", "K", "C"))
    text <- k276ReportText(k276Record(fixture))
    expect_false(k276Has(text$method, "appended to the design as covariates"))
    expect_false(k276Has(text$results, "carried the RUV-III factors as covariates"))
    expect_true(k276Has(text$results, "Figure 1 shows four quality-control views at two stages of processing"))
})

test_that("KIT-276 T3: covariates, correction, covariates leaves the covariates step as the only RUV step", {
    k276Need()
    fixture <- k276Start("in process")
    norm <- new.env(parent = emptyenv())
    k276Journey(fixture, norm, "in process", c("P", "C", "K", "R", "K"))
    history <- k276History(fixture)
    expect_identical(utils::tail(history, 2L), c("correlation_filtered", "ruv_covariates"))
    expect_false("ruv_corrected" %in% history)
    record <- k276Record(fixture)
    ruv <- Filter(function(step) identical(step$operation, "protein_ruv_covariates"), record$steps)
    expect_length(ruv, 1L)
    expect_identical(unlist(ruv[[1L]]$outcomes$ruv_mode), "covariates")
    expect_identical(withDIAStateScope(fixture$manager, function() resolveDIARuvRunDescription(fixture$manager))$ruv_mode,
        "covariates")
    k276Da(fixture)
    w_terms <- k276Model(fixture)$w_terms
    expect_gt(length(w_terms), 0L)
    expect_true(all(grepl("^W[0-9]+$", w_terms)))
})

test_that("KIT-276 T4: Skip RUV after covariates takes the RUV result out: head, history, record, DA, texts, Figure 1", {
    k276Need()
    fixture <- k276Start("in process")
    norm <- new.env(parent = emptyenv())
    k276Journey(fixture, norm, "in process", c("P", "C", "K"))
    expect_identical(fixture$manager$current_state, "ruv_covariates")
    skipped <- k276RuvSkip(fixture, norm)
    expect_identical(fixture$manager$current_state, "correlation_filtered")
    expect_identical(utils::tail(k276History(fixture), 1L), "correlation_filtered")
    expect_false("protein_ruv_covariates" %in% k276StepOperations(k276Record(fixture)))
    expect_identical(skipped$text, k276TextSkip)
    expect_identical(skipped$notes, k276TextSkip)
    expect_identical(skipped$builds, list("skip"))
    expect_null(shiny::isolate(norm$best_k))
    k276Da(fixture)
    model <- k276Model(fixture)
    expect_identical(model$w_terms, character(0))
    expect_identical(model$ruv_mode, "none")
})

test_that("KIT-276 T5: Skip RUV after the correction mode: the DA reads the uncorrected table", {
    k276Need()
    fixture <- k276Start("in process")
    norm <- new.env(parent = emptyenv())
    k276Journey(fixture, norm, "in process", c("P", "C", "R"))
    expect_identical(fixture$manager$current_state, "ruv_corrected")
    k276RuvSkip(fixture, norm)
    expect_identical(fixture$manager$current_state, "correlation_filtered")
    k276Da(fixture)
    expect_identical(k276Model(fixture)$ruv_mode, "none")
    head <- k276Object(fixture, k276HeadRef(fixture))
    filtered <- k276Object(fixture, withDIAStateScope(fixture$manager, function()
        getDIAStateReference(fixture$manager, "correlation_filtered")))
    expect_identical(head@protein_quant_table, filtered@protein_quant_table)
    corrected <- k276Object(fixture, withDIAStateScope(fixture$manager, function()
        getDIAStateReference(fixture$manager, "ruv_corrected")))
    expect_false(identical(head@protein_quant_table, corrected@protein_quant_table))
})

test_that("KIT-276 T6: Skip RUV after a DA and an enrichment on the RUV table sets both not current", {
    k276Need()
    skip_if_not_installed("Glimma")
    withr::local_envvar(MULTISCHOLAR_TEST_MODE = "true")
    enrichment <- omicGuiProteomicsEnrichment
    testthat::local_mocked_bindings(omicGuiProteomicsEnrichment = function() {
        config <- enrichment()
        config[["organism"]] <- config[["organism"]] %||% list(supported = buildProtEnrichSupportedOrganisms,
            method = resolveProtEnrichAnalysisMethod)
        config$annotation <- config$annotation %||% list(resolve = resolveProtEnrichUniprotAnnotations,
            match = resolveProtEnrichAnnotationMatching)
        config$kernel <- config$kernel %||% list(select = resolveProtEnrichProcessEnrichmentsFn,
            test_double = runProtEnrichDeterministicProcessEnrichments)
        config$organism_filter$mapping <- config$organism_filter$mapping %||% resolveProtEnrichOrganismMapping
        config$organism_filter$apply <- config$organism_filter$apply %||% applyProtEnrichOrganismFilter
        config
    }, .package = "MultiScholaR")
    fixture <- k276Start("in process")
    norm <- new.env(parent = emptyenv())
    k276Journey(fixture, norm, "in process", c("P", "C", "K"))
    manager <- fixture$manager
    manifest <- k276Da(fixture)
    withDIAStateScope(manager, function() {
        k276Quiet(produceDIADaPublication(manager, manifest, data.frame(), 0.05,
            fixture$paths[c("da_output_dir", "publication_graphs_dir")], save_assets = FALSE))
        utils::capture.output(k276Quiet(runWorkflowEnrichment(fixture$workflow, omicGuiConfig("proteomics")$enrichment,
            manifest, list(selected_contrast = "treated.vs.control", organism_taxid = "9606", up_cutoff = 0, down_cutoff = 0,
                q_cutoff = 0.05, correction_method = "gSCS", enable_organism_filter = FALSE), fixture$paths)))
    })
    expect_identical(workflowDaResultStatus(manager), "current")
    expect_false(is.null(withDIAStateScope(manager, function() diaEnrichmentProduct(manager))), info = "precondition")
    # the DA tab holds the result as the tab does
    da_data <- new.env(parent = emptyenv())
    da_data$da_results_list <- withDIAStateScope(manager, function() resolveWorkflowStateValue(manager$.context,
        manager$.context$snapshot$da_results))
    skipped <- k276RuvSkip(fixture, norm)
    expect_identical(skipped$notes, paste0(k276TextSkip, k276TextRunDa))
    expect_identical(workflowDaResultStatus(manager), "stale")
    expect_null(withDIAStateScope(manager, function() diaEnrichmentProduct(manager)))
    k276Quiet(invalidateDIAVisibleResults(da_data, fixture$workflow))
    expect_null(da_data$da_results_list)
})

test_that("KIT-276 T7: Skip RUV is refused while the managed runner has a step active; nothing is selected", {
    k276Need()
    fixture <- k276Start("in process")
    norm <- new.env(parent = emptyenv())
    k276Journey(fixture, norm, "in process", c("P", "C", "K"))
    runner <- list(step = list(active = function() "running"), run = function(...) stop("a refused Skip must not dispatch"))
    testthat::local_mocked_bindings(workflowOmicsManagedActions = function(omic_type) runner, .package = "MultiScholaR")
    skipped <- k276RuvSkip(fixture, norm)
    expect_identical(skipped$notes, k276TextBusy)
    expect_identical(fixture$manager$current_state, "ruv_covariates")
})

test_that("KIT-276 T7b: Skip RUV at a RUV head on a managed session dispatches the composite once and builds nothing in the application", {
    k276Need()
    fixture <- k276Start("in process")
    norm <- new.env(parent = emptyenv())
    k276Journey(fixture, norm, "in process", c("P", "C", "K"))
    dispatched <- list()
    runner <- list(step = list(active = function() NULL), run = function(...) stop("only the composite is dispatched"))
    testthat::local_mocked_bindings(workflowOmicsManagedActions = function(omic_type) runner,
        dispatchProtNormLimpaComposite = function(runner, workflowData, normData, experimentPaths, omicType, aesthetics, failed) {
            dispatched <<- c(dispatched, list(list(omic = omicType, head = workflowData$state_manager$current_state)))
            invisible(TRUE)
        }, .package = "MultiScholaR")
    skipped <- k276RuvSkip(fixture, norm)
    expect_identical(skipped$builds, list())
    expect_identical(dispatched, list(list(omic = "proteomics", head = "correlation_filtered")))
})

test_that("KIT-276 T8: the correlation re-applied after RUV clears the tab's RUV holders and says so; the first one does not", {
    k276Need()
    fixture <- k276Start("in process")
    norm <- new.env(parent = emptyenv())
    k276Journey(fixture, norm, "in process", c("P"))
    first <- k276Reapply(fixture, norm)
    # T12 (guard, P05): after the FIRST correlation the RUV summary is not assigned
    expect_null(first$output$ruv_covariates_summary)
    k276Journey(fixture, norm, "in process", c("K"))
    expect_false(is.null(shiny::isolate(norm$best_k)))
    expect_false(is.null(norm$ruv_optimization_result))
    fixture$workflow$ruv_optimization_result <- norm$ruv_optimization_result
    again <- k276Reapply(fixture, norm)
    expect_identical(again$output$ruv_covariates_summary, k276TextAgain)
    for (holder in c("best_k", "control_genes_index", "ruv_optimization_result", "ruv_optimization_table")) {
        expect_null(shiny::isolate(norm[[holder]]), info = holder)
    }
    expect_null(fixture$workflow$ruv_optimization_result)
    expect_identical(again$builds, list("skip"))
    # the correlation Skip at a RUV head says the same
    k276Journey(fixture, norm, "in process", c("K"))
    skipped <- k276Reapply(fixture, norm, skipped = TRUE)
    expect_identical(skipped$output$ruv_covariates_summary, k276TextAgain)
})

test_that("KIT-276 T8m: the managed Apply's finish, called as the dispatcher calls it at a RUV head, does the same", {
    k276Need()
    fixture <- k276Start("in process")
    norm <- new.env(parent = emptyenv())
    k276Journey(fixture, norm, "in process", c("P", "C", "K"))
    planted <- file.path(fixture$paths$source_dir, "ruv_optimization_results.RDS")
    saveRDS(norm$ruv_optimization_result, planted)
    finish <- NULL
    runner <- list(step = list(active = function() NULL))
    testthat::local_mocked_bindings(workflowOmicsManagedActions = function(omic_type) runner,
        dispatchProtNormCorrelation = function(runner, workflowData, normData, experimentPaths, threshold, omicType,
            experimentLabel, grouping, aesthetics, finish_fn, failed, composite = NULL) {
            finish <<- finish_fn
            invisible(TRUE)
        }, .package = "MultiScholaR")
    output <- new.env(parent = emptyenv())
    k276Quiet(runProtNormApplyCorrelationObserver(list(min_pearson_correlation_threshold = 0.5), output, fixture$workflow,
        norm, fixture$paths, "proteomics", "core", function() "group", k276Aesthetics, k276Noop,
        completeCorrelationWorkflowFn = k276Complete(output)))
    expect_true(is.function(finish))
    final <- k276Quiet(runDIAReferenceCorrelation(fixture$workflow, norm, fixture$paths, 0.5, FALSE, "proteomics", "core",
        function() "group", k276Aesthetics))
    expect_identical(k276Error(finish(final, ruv_removed = TRUE)), NA_character_)
    expect_identical(output$ruv_covariates_summary, k276TextAgain)
    expect_null(shiny::isolate(norm$best_k))
    expect_false(file.exists(planted))
})

test_that("KIT-276 T19: the managed dispatcher hands the head-before flag to finish, by value", {
    k276Need()
    run_dispatch <- function(presses) {
        fixture <- k276Start("in process")
        norm <- new.env(parent = emptyenv())
        k276Journey(fixture, norm, "in process", presses)
        recorded <- list()
        queue <- list()
        runner <- list(run = function(actions, done = NULL, failed = NULL) {
            if (is.function(actions)) actions <- actions()
            queue[[length(queue) + 1L]] <<- list(actions = actions, done = done)
            invisible(TRUE)
        })
        testthat::local_mocked_bindings(
            protNormFiguresAction = function(source, aesthetics, figures = "stage", deliver = function(result) invisible(NULL),
                composite = NULL) list(deliver = deliver),
            proteomicsDiaBinding = list(norm = list(correlation = list(pairs = "p", selection = "s",
                publish = function(...) "FINAL"))), .package = "MultiScholaR")
        k276Quiet(withDIAStateScope(fixture$manager, function() dispatchProtNormCorrelation(runner, fixture$workflow, norm,
            fixture$paths, 0.5, "proteomics", "core", function() "group", k276Aesthetics,
            function(final, ...) recorded <<- c(recorded, list(c(list(final = final), list(...)))), function(e) stop(e))))
        queue[[1L]]$actions[[1L]]$deliver(list(reference = list(entity_id = "PAIRS")))
        queue[[1L]]$done()
        queue[[2L]]$actions[[1L]]$deliver(list(value = "SELECTED"))
        queue[[2L]]$done()
        k276Quiet(queue[[3L]]$actions[[1L]]$deliver(list()))
        recorded
    }
    expect_identical(run_dispatch(c("P", "C", "K")), list(list(final = "FINAL", ruv_removed = TRUE)))
    expect_identical(run_dispatch(c("P", "C")), list(list(final = "FINAL")))
})

test_that("KIT-276 T8b: the archive's RUV fields after the re-apply say no RUV", {
    k276Need()
    fixture <- k276Start("in process")
    norm <- new.env(parent = emptyenv())
    k276Journey(fixture, norm, "in process", c("P", "C", "K"))
    k276Reapply(fixture, norm)
    fixture$manager$.context$session_format <- "native_v1"
    collected <- k276Quiet(withDIAStateScope(fixture$manager, function() collectDIANativeSessionMetadata(fixture$workflow,
        norm, list(norm_method = "quantile", ruv_mode = "covariates"))))
    expect_identical(collected$metadata$ruv_mode, "none")
    expect_identical(collected$metadata$ruv_applied, FALSE)
    expect_null(collected$metadata$ruv_k)
})

#' The leftover project file of an earlier export (the RUV run's own), planted, and the study parameters' section.
k276Leftover <- function(fixture, norm) {
    planted <- file.path(fixture$paths$source_dir, "ruv_optimization_results.RDS")
    saveRDS(norm$ruv_optimization_result, planted)
    planted
}
k276StudySection <- function(fixture) {
    k276Quiet(createWorkflowArgsFromConfig("k276", source_dir_path = fixture$paths$source_dir,
        final_s4_object = k276HeadRef(fixture), workflow_data = list(state_manager = fixture$manager, config_list = list())))
    any(grepl("Automatic RUV Optimization Results", readLines(file.path(fixture$paths$source_dir, "study_parameters.txt")),
        fixed = TRUE))
}

test_that("KIT-276 T8c: the leftover ruv_optimization_results.RDS is removed at the re-apply, the correlation Skip and the RUV Skip", {
    k276Need()
    for (how in c("apply", "correlation skip", "ruv skip")) {
        fixture <- k276Start("in process")
        norm <- new.env(parent = emptyenv())
        k276Journey(fixture, norm, "in process", c("P", "C", "K"))
        planted <- k276Leftover(fixture, norm)
        switch(how, apply = k276Reapply(fixture, norm), `correlation skip` = k276Reapply(fixture, norm, skipped = TRUE),
            `ruv skip` = k276RuvSkip(fixture, norm))
        expect_false(file.exists(planted), info = how)
        expect_false(k276StudySection(fixture), info = how)
    }
})

test_that("KIT-276 T9: Reset clears the RUV summary on the limpa route and leaves it on the other routes", {
    k276Need()
    clear <- function(output, type) {
        manager <- if (identical(type, "DIA_LIMPA")) list(workflow_type = "DIA_LIMPA") else list(workflow_type = type)
        workflow <- new.env(parent = emptyenv())
        workflow$state_manager <- manager
        output$ruv_covariates_summary <- "RUV-III covariates appended to the design: k = 2"
        testthat::local_mocked_bindings(renderText = function(expr, ...) expr, .package = "shiny")
        k276Quiet(runProtNormResetWorkflow(workflow, new.env(parent = emptyenv()), output,
            revertStateManagerFn = function(...) list(previousState = NULL), resetReactiveStateFn = k276Noop,
            resetOutputsFn = function(...) NULL, buildNotificationMessageFn = function(...) "reset",
            showNotificationFn = k276Noop, messageFn = k276Noop))
        output$ruv_covariates_summary
    }
    expect_identical(clear(new.env(parent = emptyenv()), "DIA_LIMPA"), "")
    for (type in c("DIA", "LFQ", "TMT")) {
        expect_identical(clear(new.env(parent = emptyenv()), type), "RUV-III covariates appended to the design: k = 2",
            info = type)
    }
})

test_that("KIT-276 T12, T13: the in-order journeys and the Skip with no RUV before are as before", {
    k276Need()
    for (mode in c("K", "R")) {
        fixture <- k276Start("in process")
        norm <- new.env(parent = emptyenv())
        k276Journey(fixture, norm, "in process", c("P", "C", mode))
        expect_identical(utils::tail(k276History(fixture), 2L),
            c("correlation_filtered", if (identical(mode, "K")) "ruv_covariates" else "ruv_corrected"), info = mode)
        expect_length(Filter(function(step) identical(step$operation, "protein_ruv_covariates"),
            k276Record(fixture)$steps), 1L)
    }
    fixture <- k276Start("in process")
    norm <- new.env(parent = emptyenv())
    k276Journey(fixture, norm, "in process", c("P", "C"))
    history <- k276History(fixture)
    skipped <- k276RuvSkip(fixture, norm)
    expect_identical(k276History(fixture), history)
    expect_identical(fixture$manager$current_state, "correlation_filtered")
    expect_identical(skipped$text, k276TextOldSkip)
    expect_identical(skipped$builds, list())
})

test_that("KIT-276 T14: a history written by an earlier build opens and records as before and heals at the next DA that commits and the next correlation", {
    k276Need()
    plant <- function(env = parent.frame()) {
        fixture <- k276Start("in process", envir = env)
        norm <- new.env(parent = emptyenv())
        k276Journey(fixture, norm, "in process", c("P", "C", "K"))
        context <- fixture$manager$.context
        snapshot <- context$snapshot
        expect_identical(utils::tail(unlist(snapshot$state_history), 2L), c("correlation_filtered", "ruv_covariates"))
        snapshot$current_state <- "correlation_filtered"
        withDIAStateScope(fixture$manager, function() publishDIAStateHead(context, snapshot))
        list(fixture = fixture, norm = norm)
    }
    old <- plant()
    expect_identical(utils::tail(k276History(old$fixture), 2L), c("correlation_filtered", "ruv_covariates"))
    expect_true("ruv_covariates" %in% k276StepStates(k276Record(old$fixture)))
    k276Da(old$fixture)
    expect_identical(utils::tail(k276History(old$fixture), 1L), "correlation_filtered")
    expect_false("ruv_covariates" %in% k276History(old$fixture))
    again <- plant()
    # a threshold other than the planted head's, so that the correlation commits
    expect_identical(k276Correlate(again$fixture, again$norm, "in process", 0.4), NA_character_)
    expect_false("ruv_covariates" %in% k276History(again$fixture))
})

test_that("KIT-276 T15: covariates then correction lists only the correction step; the other mode is refused as outside the history and pressing it again re-commits it", {
    k276Need()
    fixture <- k276Start("in process")
    norm <- new.env(parent = emptyenv())
    k276Journey(fixture, norm, "in process", c("P", "C", "K", "R"))
    expect_identical(utils::tail(k276History(fixture), 2L), c("correlation_filtered", "ruv_corrected"))
    expect_false("ruv_covariates" %in% k276History(fixture))
    expect_error(withDIAStateScope(fixture$manager, function() selectDIAStateReference(fixture$manager, "ruv_covariates")),
        class = "saved_entity_scope")
    k276Journey(fixture, norm, "in process", c("K"))
    expect_identical(utils::tail(k276History(fixture), 2L), c("correlation_filtered", "ruv_covariates"))
})

test_that("KIT-276 T16: the classic route's history is as before (RUV skipped, correlation, normalisation again)", {
    fixture <- protNormCoreOwner()
    norm <- new.env(parent = emptyenv())
    settings <- list(norm_method = "quantile", ruv_mode = "skip", auto_percentage_min = NULL, auto_percentage_max = NULL,
        separation_metric = NULL, k_penalty_weight = NULL, adaptive_k_penalty = NULL, ruv_percentage = 50, ruv_k = 1,
        color_variable = "group", shape_variable = NULL)
    normalise <- function() protNormCoreQuiet(runDIAReferenceNormalization(settings, fixture$workflow, norm, fixture$paths,
        "proteomics", "core", k276Noop, k276Noop, function() "group"))
    normalise()
    protNormCoreQuiet(runDIAReferenceCorrelation(fixture$workflow, norm, fixture$paths, 0.5, FALSE, "proteomics", "core",
        function() "group", k276Aesthetics))
    normalise()
    expect_identical(utils::tail(diaStateHistory(fixture$manager), 2L), c("normalized", "correlation_filtered"))
    expect_identical(fixture$manager$current_state, "normalized")
})

test_that("KIT-276 T17: the record built after T1 replays to correlation_filtered and its DA; the replay's own comparison", {
    k276Need()
    fixture <- k276Start("in process")
    norm <- new.env(parent = emptyenv())
    k276Journey(fixture, norm, "in process", c("P", "C", "K", "C"))
    k276Da(fixture)
    replay <- k276Replay(fixture, k276Record(fixture))
    expect_identical(k276ErrorClass(replay$error), NA_character_)
    expect_identical(replay$states, k276History(fixture))
    # The origin store keeps the superseded `ruv_covariates` record by design (ticket 4.1, VAL-346); the record, the history
    # and the replay do not list it, so the comparison over every state of the origin names it, and only it.
    expect_identical(replay$mismatches, "ruv_covariates/missing")
    # Over the states the history lists the replay agrees with the origin.
    expect_identical(replay$listed, character())
})

test_that("KIT-276 T18: the replay of covariates, correction, covariates computes the RUV step from correlation_filtered", {
    k276Need()
    fixture <- k276Start("in process")
    norm <- new.env(parent = emptyenv())
    k276Journey(fixture, norm, "in process", c("P", "C", "K", "R", "K"))
    k276Da(fixture)
    replay <- k276Replay(fixture, k276Record(fixture))
    expect_identical(k276ErrorClass(replay$error), NA_character_)
    expect_identical(replay$states, k276History(fixture))
    # (i) over the states the history lists: the replay agrees with the origin.
    expect_identical(replay$listed, character())
    # (ii) over every state of the origin store, as the replay journey reports it: the origin keeps the record of the
    # superseded `ruv_corrected` state by design (ticket 4.1, VAL-346), which the record, the history and the replay do not
    # list; that state and no other (no `ruv_covariates/*` entry) is what the comparison names.
    expect_identical(replay$mismatches, "ruv_corrected/missing")
})

test_that("KIT-276 note 1.5: a declared-rule run with no sweep leaves no canonical-correlation plot from an earlier run", {
    paths <- list(protein_qc_dir = withr::local_tempdir())
    cancor <- file.path(paths$protein_qc_dir, "ruv_covariates_cancor.png")
    marker <- paste0(cancor, ".placeholder")
    writeLines("x", cancor)
    writeLines("x", marker)
    stages <- file.path(paths$protein_qc_dir, paste0(rep(c("pre_norm", "post_norm", "ruv_corrected"), 4L), "_",
        rep(c("pca", "density", "rle", "correlation"), each = 3L), ".png"))
    for (stage in stages) writeLines("x", stage)
    removeDiaLimpaRuvMarkers(paths)
    expect_false(file.exists(cancor))
    expect_false(file.exists(marker))
    # panel m is not drawn from the earlier run: the layout has the twelve stage panels only
    plot_files <- resolveDiaLimpaCompositeFigureInputs("covariates", paths$protein_qc_dir)$plotFiles
    expect_length(plot_files, 12L)
    expect_false(cancor %in% plot_files)
})
