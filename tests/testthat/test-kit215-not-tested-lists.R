# KIT-215 (VAL-170): a proteomics enrichment list that was not tested is shown and reported as "Not tested", with its
# route's reason, never as a list with no enriched term; a direction with no list says "no list"; a GO comparison that
# cannot be tested (no annotated tested protein, or no term of 10 to 500 of them) and a test that raised an error refuse
# the run. Tests T1 to T12 and guards G1, G3, G5 of the ticket (section 10; r1.3 section 18). No service is called:
# `gprofiler2::gost` is replaced by a recorder that answers NULL (KIT-208's, test-kit208-gprofiler-identifiers.R), the
# proxy points at a closed port, and the clusterProfiler route runs the installed `enricher` on KIT-173's fixture
# (helper-go-propagation.R) or on a synthetic direct-annotation table. Arguments a function on the base commit does not
# take are dropped (`k215Call`), so each failing-first row fails at its stated base value.

library(testthat)

k215Offline <- function(env = parent.frame()) {
    withr::local_envvar(http_proxy = "http://127.0.0.1:9", https_proxy = "http://127.0.0.1:9",
        HTTP_PROXY = "http://127.0.0.1:9", HTTPS_PROXY = "http://127.0.0.1:9", no_proxy = "", NO_PROXY = "",
        MULTISCHOLAR_TEST_MODE = "false", .local_envir = env)
    withr::local_options(multischolar.test_mode = FALSE, .local_envir = env)
}

#' KIT-208's gost recorder (k208RecordGost), with `raise_size`: a query of that many distinct names raises instead
k215RecordGost <- function(raise_size = NULL, env = parent.frame()) {
    log <- new.env(parent = emptyenv())
    log$calls <- character()
    restore <- goPropReplaceBinding("gprofiler2", "gost", function(query, ..., custom_bg = NULL) {
        size <- length(unique(unlist(query, use.names = FALSE)))
        log$calls <- c(log$calls, sprintf("%d,%d", size, length(unique(custom_bg))))
        if (!is.null(raise_size) && identical(size, as.integer(raise_size))) stop("k215 forced service error")
        NULL
    })
    withr::defer(restore(), envir = env)
    log
}

#' The installed `clusterProfiler::enricher`, wrapped: each call's list and whether it returned NULL are logged; a list
#' for which `raise_when(gene)` is TRUE raises instead
k215WrapEnricher <- function(raise_when = NULL, env = parent.frame()) {
    log <- new.env(parent = emptyenv())
    log$calls <- list()
    real <- get("enricher", envir = asNamespace("clusterProfiler"))
    restore <- goPropReplaceBinding("clusterProfiler", "enricher", function(gene, ...) {
        if (!is.null(raise_when) && isTRUE(raise_when(gene))) stop("k215 forced enricher error")
        result <- real(gene = gene, ...)
        log$calls[[length(log$calls) + 1L]] <- list(gene = sort(unique(gene)), null = is.null(result))
        result
    })
    withr::defer(restore(), envir = env)
    log
}

#' Run `expr` with console output and messages muffled; the error (if any) and the warnings' messages are returned
k215Capture <- function(expr) {
    error <- NULL
    warnings <- character()
    utils::capture.output(tryCatch(withCallingHandlers(expr,
        warning = function(w) {
            warnings <<- c(warnings, conditionMessage(w))
            invokeRestart("muffleWarning")
        },
        message = function(m) invokeRestart("muffleMessage")),
        error = function(e) {
            error <<- e
            NULL
        }))
    list(error = error, warnings = warnings)
}

#' Call `fn` with the arguments it takes (the base commit's display functions lack the new ones)
k215Call <- function(fn, args) do.call(fn, args[names(args) %in% names(formals(fn))])

k215Roles <- function() goPropReadTsv("t7-fixture-entries.tsv")

#' KIT-173's E. coli fixture annotated through `getUniprotAnnotations` (no service; 70 rows), as KIT-208's
k215Annotations <- function() goPropAnnotate(goPropReadTsv("t7-annotation-ecoli.tsv"), withr::local_tempdir())

#' KIT-208's T7 annotation (test-kit208-gprofiler-identifiers.R:227): the p_term proteins (the up list) removed
k215LacksUp <- function(full) {
    roles <- k215Roles()
    full[!full$Entry %in% roles$entry[roles$role == "p_term"], , drop = FALSE]
}

#' The enrichment record's metadata, whether KIT-216 has merged or not (ticket section 4)
k215EnrichmentRecord <- function(manager) {
    snapshot <- manager$.context$snapshot
    reference <- snapshot$enrichment_results$record %||% snapshot$states$enrichment_completed
    if (is.null(reference)) NULL else readComponentMetadata(manager, reference)
}

#' A copy of KIT-208's `k208Project` (test-kit208-gprofiler-identifiers.R, itself KIT-175's `k175Project`), on the
#' clusterProfiler route (taxon 999999) by default, with the cutoffs as arguments, `enricher` wrapped (`raise_when`),
#' the analysis record written when the run completes, and `after(manager, workflow, run)` run in the owner's scope
#' before the owner is closed. `gost` is always replaced by the recorder, so no service is called.
k215Project <- function(root, annotations, taxon = "999999", id_gene_names = FALSE, up_cutoff = 0, down_cutoff = 0,
    q_cutoff = 0.05, raise_when = NULL, raise_size = NULL, after = NULL) {
    k215Offline()
    gost <- k215RecordGost(raise_size)
    enricher <- k215WrapEnricher(raise_when)
    roles <- k215Roles()
    paths <- goPropPaths(root)
    for (p in paths) dir.create(p, recursive = TRUE, showWarnings = FALSE)
    saveRDS(annotations, file.path(paths$source_dir, "uniprot_dat_cln.RDS"))
    manager <- newDIAStateAdapter(newWorkflowEntityStore(root, "k215-project", "DIA"),
        WorkflowState$new(audit_enabled = FALSE))
    on.exit(closeDIAStateAdapter(manager), add = TRUE)
    manager$.context$reference_states <- TRUE
    manager$setWorkflowType("DIA")
    entries <- sort(roles$entry)
    design <- data.frame(Run = paste0("S", 1:6), group = rep(c("A", "B"), each = 3L),
        batch = rep(c("b1", "b2", "b3"), times = 2L), stringsAsFactors = FALSE)
    quant <- goPropQuant(entries, roles$entry[roles$role == "p_term"], roles$entry[roles$role == "d_term"])
    id_table <- data.frame(uniprot_acc = entries, stringsAsFactors = FALSE)
    if (id_gene_names) id_table$gene_names <- paste0("GENE", seq_along(entries))
    object <- methods::new("ProteinQuantitativeData", protein_quant_table = quant, design_matrix = design,
        sample_id = "Run", group_id = "group", protein_id_column = "uniprot_acc", protein_id_table = id_table,
        args = list(globalParameters = list(workflow_type = "DIA", use_limpa = FALSE, report_template = "DIA_report"),
            deAnalysisParameters = list(formula_string = "~ 0 + group")))
    setDIAStepProvenance(manager, "design_builder")
    manager$saveState("correlation_filtered", object, list(), "fixture")
    manager$.context$experiment_paths <- paths
    contrasts <- data.frame(contrasts = "groupB-groupA", full_format = "B_vs_A=groupB-groupA",
        friendly_names = "B_vs_A", stringsAsFactors = FALSE)
    workflow <- list2env(list(state_manager = manager, contrasts_tbl = contrasts, tab_status = list(),
        qc_params = list(protein_qc = list()), config_list = list()))
    parameters <- list(selected_contrast = "B_vs_A", organism_taxid = taxon, up_cutoff = up_cutoff,
        down_cutoff = down_cutoff, q_cutoff = q_cutoff, correction_method = "gSCS", enable_organism_filter = FALSE)
    notes <- list()
    run <- NULL
    captured <- k215Capture(withDIAStateScope(manager, function() {
        da_data <- new.env(parent = emptyenv())
        da_data$current_s4_object <- getDIAStateReference(manager)
        da_data$contrasts_available <- contrasts$contrasts
        manifest <- runDIAReferenceDaAnalysis(list(formula_string = "~ 0 + group", da_q_val_thresh = 0.05,
            treat_lfc_cutoff = 0), da_data, workflow, paths, contrasts = contrasts, publish = FALSE)
        produceDIADaPublication(manager, manifest, annotations, 0.05,
            paths[c("da_output_dir", "publication_graphs_dir")], save_assets = FALSE)
        run <<- runWorkflowEnrichment(workflow, omicGuiConfig("proteomics")$enrichment, manifest, parameters, paths,
            process_fn = function(...) processEnrichments(...), notify = function(message, type, duration) {
                notes[[length(notes) + 1L]] <<- list(text = message, type = type, duration = duration)
            })
    }))
    record <- NULL
    if (!is.null(run)) {
        record <- file.path(root, "analysis_record.json")
        writeWorkflowAnalysisRecord(manager, record)
    }
    extra <- if (!is.null(after)) withDIAStateScope(manager, function() after(manager, workflow, run))
    metadata <- k215EnrichmentRecord(manager)
    list(run = run, error = captured$error, warnings = captured$warnings,
        notes = vapply(notes, function(note) note$text, ""), types = vapply(notes, function(note) note$type, ""),
        gost = gost$calls, enricher = enricher$calls, config = metadata$config, outcomes = metadata$config$outcomes,
        record = record, root = root, extra = extra)
}

#' A project made once per file (the T1, G1, T5 and g:Profiler-complete projects are read by several tests)
k215Cached <- local({
    cache <- list()
    function(key, make) {
        if (is.null(cache[[key]])) cache[[key]] <<- make()
        cache[[key]]
    }
})

#' The tab after a simulated reopen of the project's committed state (KIT-024's measure, as W1's `restore` case)
k215Reopen <- function(manager, workflow, run) {
    holder <- shiny::reactiveValues(all_enrichment_results = NULL, enrichment_results_full = NULL,
        analysis_complete = FALSE)
    restored <- restoreProtEnrichResultsFromState(holder, workflow)
    list(restored = restored, results = shiny::isolate(holder$enrichment_results_full),
        not_tested = shiny::isolate(holder$enrichment_not_tested), no_list = shiny::isolate(holder$enrichment_no_list))
}

k215T1 <- function() k215Cached("t1", function() k215Project(tempfile("k215-t1-"), k215LacksUp(k215Annotations()),
    after = k215Reopen))
k215G1 <- function() k215Cached("g1", function() k215Project(tempfile("k215-g1-"), k215Annotations()))
k215T5 <- function() k215Cached("t5", function() k215Project(tempfile("k215-t5-"), k215LacksUp(k215Annotations()),
    taxon = "9606"))
k215GpComplete <- function() k215Cached("gp", function() k215Project(tempfile("k215-gp-"), k215Annotations(),
    taxon = "9606", id_gene_names = TRUE))

#' The text a plot renderer shows (its placeholder's text; a plot's class when it draws one)
k215Plot <- function(renderer, results, direction, not_tested = NULL) {
    shown <- k215Call(renderer, list(analysisComplete = TRUE, enrichmentResultsFull = results,
        rawContrast = names(results@enrichment_data)[[1L]], directionFilter = direction, notTested = not_tested,
        renderPlotlyFn = function(expr) expr, plotLyFn = function(...) list(),
        addTextFn = function(p, x, y, text, showlegend) text))
    if (is.character(shown)) shown else paste(class(shown), collapse = "/")
}

k215Summary <- function(formatter, results, run, filter = "all") {
    k215Call(formatter, list(gprofilerResults = results, clusterprofilerResults = results, directionFilter = filter,
        notTested = run$not_tested, noList = run$no_list, selectedContrast = "B_vs_A"))
}

#' The download archive's summary text (no file is written)
k215Download <- function(run) {
    lines <- character()
    k215Call(writeProtEnrichResultsDownloadArchive, list(file = tempfile(fileext = ".zip"), selectedContrast = "B_vs_A",
        methodInfo = run$method, organismTaxid = "999999", upCutoff = 0, downCutoff = 0, qCutoff = 0.05,
        gprofilerResults = run$selected$gprofilerResults, clusterprofilerResults = run$selected$clusterprofilerResults,
        tempDir = tempdir(), writeTsvFn = function(...) invisible(NULL),
        writeLinesFn = function(text, path) lines <<- c(lines, text), zipFn = function(...) invisible(NULL),
        sysTimeFn = function() "<time>", significanceRule = "bh_whole_family", notTested = run$not_tested,
        noList = run$no_list))
    paste(lines, collapse = "\n")
}

#' A template's `FormatGOEnrichment`, extracted with knitr::purl and evaluated (KIT-205's k205Formatter, KIT-173 T11's
#' way)
k215Formatter <- function(template) {
    lines <- readLines(system.file("reports", "proteomics", template, package = "MultiScholaR"), warn = FALSE,
        encoding = "UTF-8")
    start <- grep("^```\\{r GO Enrichment Section", lines)
    stopifnot(length(start) == 1L)
    end <- start + which(lines[-seq_len(start)] == "```")[[1L]]
    chunk <- tempfile(fileext = ".Rmd")
    writeLines(lines[start:end], chunk, useBytes = TRUE)
    code <- tempfile(fileext = ".R")
    suppressMessages(knitr::purl(chunk, output = code, quiet = TRUE, documentation = 0))
    exprs <- parse(code, keep.source = FALSE, encoding = "UTF-8")
    at <- which(vapply(exprs, function(e) is.call(e) && identical(e[[1L]], as.name("<-")) &&
        identical(e[[2L]], as.name("FormatGOEnrichment")), logical(1)))
    env <- new.env(parent = globalenv())
    eval(exprs[[at]], env)
    env$FormatGOEnrichment
}

k215Templates <- c("LFQ_report.rmd", "TMT_report.rmd", "DIANN_report.rmd", "DIANN_limpa_report.rmd")

#' The report numbers of a record (NULL: no record), as a report reads them
k215Numbers <- function(record, source) {
    numbers <- list(ui_enrichment_organism_selected = if (identical(source, "gprofiler2")) "9606" else "999999",
        ui_enrichment_up_lfc = 0, ui_enrichment_down_lfc = 0, ui_enrichment_q_value_cutoff = 0.05,
        ui_enrichment_database_source = source)
    if (is.null(record)) return(numbers)
    params <- workflowRecordReportParameters(record, list(ui_enrichment_database_source = source))
    c(numbers, params[grepl("^enrichment_", names(params)) | names(params) == "record_present"])
}

k215Limits <- paste("On g:Profiler, a list for which no term is significant cannot be told apart from a list in which",
    "the service recognised none of the gene names.")

k215GoReason <- "none of its 25 proteins has a Gene Ontology annotation that this analysis tests"
k215GpReason <- "none of its 25 proteins has a gene name, and g:Profiler tests gene names"
k215Placeholder <- function(reason) paste0("Not tested: ", reason,
    ". An empty result here does not mean that no term is enriched.")

k215NotTestedGo <- paste("Comparison B_vs_A, up: none of its 25 proteins has a Gene Ontology annotation that this",
    "analysis tests, so this list was not tested (its empty result does not mean that no term is enriched).")

# KIT-208's notice (test-kit208-gprofiler-identifiers.R), unchanged byte for byte (E-7)
k215NotTestedGp <- paste("Comparison B_vs_A, up: none of its 25 proteins has a gene name, so this list was not sent to",
    "g:Profiler and was not tested (its empty result does not mean that no term is enriched).")

k215Bullet <- function(reason) paste0("• Not tested: B_vs_A (up-regulated): ", reason, ". A list that was not ",
    "tested has no result; its absence from the results does not mean that no term is enriched.")

# --- the clusterProfiler route: a list not tested (C2) ---------------------------------------------------------------

test_that("KIT-215 T1: a GO list without an annotated protein is reported once, recorded and committed", {
    k215Offline()
    t1 <- k215T1()
    expect_null(t1$error) # base: completed
    expect_identical(sum(t1$notes == k215NotTestedGo), 1L) # base: no notice containing "not tested"
    expect_identical(sum(grepl("not tested", t1$notes, fixed = TRUE)), 1L)
    expect_identical(t1$types[match(k215NotTestedGo, t1$notes)], "warning")
    expect_identical(t1$outcomes$enrichment_lists_not_tested, "B_vs_A: up") # base: NULL
    expect_identical(NROW(t1$config$not_tested), 1L) # base: 0
    expect_identical(t1$config$not_tested$contrast, "groupB-groupA")
    expect_identical(t1$config$not_tested$direction, "up")
    expect_identical(t1$config$not_tested$route, "clusterprofiler")
    expect_null(t1$config$no_list)
    data <- t1$run$results@enrichment_data[[1L]]
    expect_null(data$up)
    expect_false(is.null(data$down)) # the down list is tested, as on base
})

test_that("KIT-215 T2: the tab shows the GO list as not tested, with its reason", {
    k215Offline()
    run <- k215T1()$run
    expect_identical(k215Plot(renderProtEnrichClusterProfilerPlot, run$results, "up", run$not_tested),
        k215Placeholder(k215GoReason)) # base: "No up-regulated enrichment data"
    summary <- k215Summary(formatProtEnrichClusterProfilerSummaryText, run$selected$clusterprofilerResults, run)
    expect_true(grepl(sprintf("Up-regulated: not tested (%s)", k215GoReason), summary, fixed = TRUE),
        info = summary) # base: "Up-regulated: 0"
    # K215-I01: the summary filtered to the direction not tested says so, not a count of 0.
    filtered <- k215Summary(formatProtEnrichClusterProfilerSummaryText, run$selected$clusterprofilerResults, run,
        filter = "up")
    expect_true(grepl(sprintf("Up-regulated: not tested (%s)", k215GoReason), filtered, fixed = TRUE),
        info = filtered) # base: "Showing 0 up-regulated GO terms"
})

test_that("KIT-215 T2b: after a reopen the tab shows the same texts", {
    k215Offline()
    t1 <- k215T1()
    reopened <- t1$extra
    expect_true(reopened$restored)
    expect_identical(NROW(reopened$not_tested), 1L) # base: 0
    expect_identical(k215Plot(renderProtEnrichClusterProfilerPlot, reopened$results, "up", reopened$not_tested),
        k215Placeholder(k215GoReason))
    summary <- k215Summary(formatProtEnrichClusterProfilerSummaryText, t1$run$selected$clusterprofilerResults,
        list(not_tested = reopened$not_tested, no_list = reopened$no_list))
    expect_true(grepl(sprintf("Up-regulated: not tested (%s)", k215GoReason), summary, fixed = TRUE), info = summary)
})

test_that("KIT-215 T3: the download summary has a section of the lists not tested", {
    k215Offline()
    text <- k215Download(k215T1()$run)
    expect_true(grepl("## Not tested:", text, fixed = TRUE), info = text) # base: no such section
    expect_true(grepl(sprintf("B_vs_A, up-regulated: %s", k215GoReason), text, fixed = TRUE), info = text)
})

# --- the clusterProfiler route: comparisons that cannot be tested (C3, C4) ------------------------------------------

k215C3 <- paste("None of the tested proteins of comparison B_vs_A has a Gene Ontology annotation that this analysis",
    "tests (0 of 70 tested proteins annotated), so nothing would be tested. Check that the protein annotation was",
    "retrieved for this organism (save the design again). This happens when the annotation was not retrieved for",
    "these proteins, or when their identifiers are not UniProt accessions.")

k215Unmatched <- function() {
    unmatched <- k215Annotations()
    unmatched$Entry <- sprintf("Z%05d", seq_len(nrow(unmatched)))
    unmatched
}

test_that("KIT-215 T4: a GO comparison with lists and no annotated tested protein refuses the run (C3)", {
    k215Offline()
    project <- k215Project(withr::local_tempdir(), k215Unmatched())
    expect_true(inherits(project$error, "workflow_enrichment_no_identifiers")) # base: completed, every list NULL
    expect_identical(conditionMessage(project$error), k215C3)
    expect_length(project$enricher, 2L)
    expect_true(all(vapply(project$enricher, function(call) call$null, NA)))
    expect_null(project$config) # nothing committed
})

test_that("KIT-215 T4b (head): a record replayed where the comparison has no annotated tested protein is refused", {
    k215Offline()
    origin <- k215G1()
    expect_null(origin$error)
    copy <- withr::local_tempdir()
    stopifnot(system2("cp", c("-a", file.path(origin$root, "."), copy)) == 0L)
    saveRDS(k215Unmatched(), file.path(copy, "source", "uniprot_dat_cln.RDS"))
    h <- new.env(parent = globalenv())
    sys.source(file.path(testthat::test_path("..", ".."), "tools", "profiling", "workflow_kit.R"), h)
    api <- function(name) get(name, asNamespace("MultiScholaR"))
    prepared <- h$kitPrepareReplayProject(copy, c("initial", "correlation_filtered"), api)
    record <- jsonlite::read_json(origin$record, simplifyVector = FALSE)
    replay <- k215Capture(replayWorkflowAnalysisRecord(record, prepared, character(), goPropPaths(copy),
        withr::local_tempdir()))
    expect_true(inherits(replay$error, "workflow_enrichment_no_identifiers"), info = conditionMessage(replay$error %||%
        simpleError("no error")))
})

test_that("KIT-215 T10: a GO comparison with lists and no term of 10 to 500 annotated tested proteins refuses (C4)", {
    k215Offline()
    roles <- k215Roles()
    full <- k215Annotations()
    nine <- full[full$Entry %in% sort(roles$entry[roles$role == "p_term"])[seq_len(9L)], , drop = FALSE]
    expect_identical(nrow(nine), 9L)
    project <- k215Project(withr::local_tempdir(), nine)
    expect_true(inherits(project$error, "workflow_enrichment_no_identifiers")) # base: completed, every list NULL
    expect_match(conditionMessage(project$error), paste0("^In comparison B_vs_A, [1-9] of the tested proteins have a ",
        "Gene Ontology annotation, and no Gene Ontology term holds between 10 and 500 of them, so no term could be ",
        "tested\\.$"))
    expect_true(length(project$enricher) > 0L)
    expect_true(all(vapply(project$enricher, function(call) call$null, NA)))
    expect_null(project$config)
})

# --- g:Profiler: KIT-208's list not tested (B3), shown --------------------------------------------------------------

test_that("KIT-215 T5: KIT-208's g:Profiler list without a gene name is shown as not tested", {
    k215Offline()
    t5 <- k215T5()
    expect_null(t5$error)
    expect_identical(sum(t5$notes == k215NotTestedGp), 1L) # KIT-208's notice, unchanged
    expect_identical(sum(grepl("was not tested", t5$notes, fixed = TRUE)), 1L)
    expect_identical(t5$outcomes$enrichment_lists_not_tested, "B_vs_A: up")
    expect_identical(t5$config$not_tested$route, "gprofiler2") # base: no such column
    expect_identical(k215Plot(renderProtEnrichGprofilerPlot, t5$run$results, "up", t5$run$not_tested),
        k215Placeholder(k215GpReason)) # base: "No up-regulated enrichment data"
    summary <- k215Summary(formatProtEnrichGprofilerSummaryText, t5$run$selected$gprofilerResults, t5$run)
    expect_true(grepl(sprintf("Up-regulated pathways: not tested (%s)", k215GpReason), summary, fixed = TRUE),
        info = summary) # base: "No gprofiler2 results available."
})

# --- the reports -----------------------------------------------------------------------------------------------------

test_that("KIT-215 T6: each template states the list not tested, with the route's reason; g:Profiler its limit", {
    k215Offline()
    go <- k215Numbers(k215T1()$record, "clusterprofiler")
    gp <- k215Numbers(k215T5()$record, "gprofiler2")
    expect_identical(go$enrichment_lists_not_tested, "B_vs_A: up") # base: NULL (not in the report map)
    for (template in k215Templates) {
        format <- k215Formatter(template)
        go_text <- as.character(format(go))
        gp_text <- as.character(format(gp))
        expect_true(grepl(k215Bullet("none of its proteins has a Gene Ontology annotation that this analysis tests"),
            go_text, fixed = TRUE), info = template)
        expect_true(grepl(k215Bullet("none of its proteins has a gene name, and g:Profiler tests gene names"),
            gp_text, fixed = TRUE), info = template)
        expect_true(grepl(k215Limits, gp_text, fixed = TRUE), info = template)
        expect_false(grepl(k215Limits, go_text, fixed = TRUE), info = template)
    }
})

test_that("KIT-215 G5 (guard): no 'Not tested' without the outcome; texts printed for base against head", {
    k215Offline()
    cases <- list(go_record = k215Numbers(k215G1()$record, "clusterprofiler"),
        gp_record = k215Numbers(k215GpComplete()$record, "gprofiler2"),
        go_none = k215Numbers(NULL, "clusterprofiler"), gp_none = k215Numbers(NULL, "gprofiler2"))
    for (template in k215Templates) {
        format <- k215Formatter(template)
        for (case in names(cases)) {
            text <- as.character(format(cases[[case]]))
            expect_false(grepl("Not tested", text, fixed = TRUE), info = paste(template, case))
            without <- sub(paste0(" ", k215Limits), "", text, fixed = TRUE)
            cat(sprintf("K215-G5 %s %s %s\n", template, case, digest::digest(without, algo = "sha256")))
        }
    }
})

# --- errors refuse the run (VAL-170 (4)) -----------------------------------------------------------------------------

test_that("KIT-215 T7: an error inside a GO list's test refuses the run", {
    k215Offline()
    roles <- k215Roles()
    down <- sort(roles$entry[roles$role == "d_term"])
    project <- k215Project(withr::local_tempdir(), k215Annotations(),
        raise_when = function(gene) identical(sort(unique(as.character(gene))), down))
    expect_true(inherits(project$error, "workflow_enrichment_test_error")) # base: completed, down NULL, warning
    expect_identical(conditionMessage(project$error), paste("The enrichment test of comparison B_vs_A, down stopped",
        "with an error, so the run was not completed: k215 forced enricher error."))
    expect_null(project$config)
})

test_that("KIT-215 T7b: an error in a g:Profiler comparison's processing refuses the run", {
    k215Offline()
    testthat::local_mocked_bindings(stripEnrichmentIdSuffix = function(ids) stop("k215 forced strip error"))
    project <- k215Project(withr::local_tempdir(), k215Annotations(), taxon = "9606", id_gene_names = TRUE)
    expect_true(inherits(project$error, "workflow_enrichment_test_error")) # base: completed, the contrast NULL
    expect_identical(conditionMessage(project$error), paste("The enrichment test of comparison B_vs_A stopped with an",
        "error, so the run was not completed: k215 forced strip error."))
    expect_null(project$config)
})

test_that("KIT-215 T7c: an error inside a g:Profiler list's call refuses the run", {
    k215Offline()
    real <- perform_enrichment
    testthat::local_mocked_bindings(perform_enrichment = function(data_subset, ...) {
        if (all(data_subset$log2FC < 0)) stop("k215 forced list error")
        real(data_subset = data_subset, ...)
    })
    project <- k215Project(withr::local_tempdir(), k215Annotations(), taxon = "9606", id_gene_names = TRUE)
    expect_true(inherits(project$error, "workflow_enrichment_test_error")) # base: completed, down NULL, warning only
    expect_identical(conditionMessage(project$error), paste("The enrichment test of comparison B_vs_A, down stopped",
        "with an error, so the run was not completed: k215 forced list error."))
    expect_true("Error processing down-regulated genes: k215 forced list error" %in% project$warnings)
    expect_null(project$config)
})

test_that("KIT-215 T7d (guard): a service that never answers is refused as a service failure first", {
    k215Offline()
    real <- perform_enrichment
    testthat::local_mocked_bindings(perform_enrichment = function(...) real(..., wait_time = 0))
    project <- k215Project(withr::local_tempdir(), k215Annotations(), taxon = "9606", id_gene_names = TRUE,
        raise_size = 10L)
    expect_true(inherits(project$error, "workflow_enrichment_service"))
    expect_null(project$config)
})

# --- the kernel on a synthetic direct-annotation table (T8, G3) ------------------------------------------------------

#' Direct GO annotations (`go_propagation = "none"`): `sets` names each term's proteins; the first protein also has
#' one molecular-function and one cellular-component term of its own (so each aspect has a row)
k215Synthetic <- function(sets) {
    pairs <- data.frame(term = rep(names(sets), lengths(sets)), entry = unlist(sets, use.names = FALSE),
        stringsAsFactors = FALSE)
    entries <- unique(pairs$entry)
    ids <- vapply(entries, function(e) paste(pairs$term[pairs$entry == e], collapse = "; "), "", USE.NAMES = FALSE)
    out <- data.frame(Entry = entries, go_id_go_biological_process = ids,
        go_term_go_biological_process = gsub("GO:", "k215 term ", ids, fixed = TRUE),
        go_id_go_molecular_function = NA_character_, go_term_go_molecular_function = NA_character_,
        go_id_go_cellular_compartment = NA_character_, go_term_go_cellular_compartment = NA_character_,
        stringsAsFactors = FALSE)
    out[1L, c("go_id_go_molecular_function", "go_term_go_molecular_function")] <- c("GO:0099911", "k215 term 0099911")
    out[1L, c("go_id_go_cellular_compartment", "go_term_go_cellular_compartment")] <- c("GO:0099921",
        "k215 term 0099921")
    out
}

#' One comparison's DA table: `up` and `down` past the cutoffs (q 0.001), every other protein not (q 0.6)
k215Table <- function(proteins, up = character(), down = character()) {
    data.frame(Protein.IDs = proteins, log2FC = ifelse(proteins %in% up, 2, ifelse(proteins %in% down, -2, 0.1)),
        raw_pvalue = ifelse(proteins %in% c(up, down), 0.0001, 0.5),
        fdr_qvalue = ifelse(proteins %in% c(up, down), 0.001, 0.6), stringsAsFactors = FALSE)
}

#' The enrichment declaration's kernel on `tables` (the clusterProfiler route, direct annotations), as KIT-208's T10
k215Kernel <- function(tables, annotations, env = parent.frame()) {
    testthat::local_mocked_bindings(saveSelfContainedWidget = function(...) invisible(NULL), .env = env)
    testthat::local_mocked_bindings(ggsave = function(...) invisible(NULL), .package = "ggplot2", .env = env)
    da <- methods::new("da_results_for_enrichment", contrasts = tibble::tibble(contrasts = names(tables),
        friendly_names = sub("^group([A-Z])-group([A-Z])$", "\\1_vs_\\2", names(tables))), da_data = tables,
        design_matrix = data.frame(Run = "S1", stringsAsFactors = FALSE))
    parameters <- list(da_results = da, taxon_id = "999999", up_cutoff = 0, down_cutoff = 0, q_cutoff = 0.05,
        pathway_dir = withr::local_tempdir(.local_envir = env), go_annotations = annotations, exclude_iea = FALSE,
        protein_id_column = "Protein.IDs", contrast_names = names(tables), correction_method = "BH",
        go_propagation = "none", process_fn = processEnrichments)
    value <- NULL
    captured <- k215Capture(value <- diaEnrichmentOperation()$kernel(NULL, parameters, NULL))
    list(value = value, error = captured$error)
}

k215Members <- sprintf("P%02d", 1:12)

test_that("KIT-215 T8 (guard, C5): a list whose annotated proteins are in no tested term is tested, with no term", {
    skip_if_not_installed("clusterProfiler")
    annotations <- k215Synthetic(list(`GO:0099901` = k215Members, `GO:0099902` = c("L1", "L2", "L3")))
    out <- k215Kernel(list(`groupB-groupA` = k215Table(c(k215Members, "L1", "L2", "L3", "U1", "U2"),
        up = c("L1", "L2", "L3"), down = c("P01", "P02"))), annotations)
    expect_null(out$error)
    expect_null(out$value$not_tested)
    expect_null(out$value$results@enrichment_data[[1L]]$up) # m' = 0: enricher's NULL, a tested list
})

test_that("KIT-215 G3 (guard): a partly annotated comparison and one with no list and no annotation complete", {
    skip_if_not_installed("clusterProfiler")
    annotations <- k215Synthetic(list(`GO:0099901` = k215Members))
    tables <- list(`groupB-groupA` = k215Table(c(k215Members, sprintf("U%d", 1:6)), up = k215Members[1:5],
        down = k215Members[6:8]), `groupC-groupA` = k215Table(sprintf("V%02d", 1:10)))
    out <- k215Kernel(tables, annotations)
    expect_null(out$error)
    expect_null(out$value$not_tested)
    expect_false(is.null(out$value$results@enrichment_data[[1L]]$up))
})

# --- a direction with no list (VAL-170 (5)) --------------------------------------------------------------------------

test_that("KIT-215 T11: a GO direction with no list says so in the summary and the download", {
    k215Offline()
    project <- k215Project(withr::local_tempdir(), k215Annotations(), down_cutoff = 5)
    expect_null(project$error)
    expect_identical(NROW(project$config$no_list), 1L) # base: no `no_list`
    expect_identical(project$config$no_list$direction, "down")
    expect_identical(project$config$no_list$comparison, "B_vs_A")
    expect_null(project$outcomes$enrichment_lists_not_tested)
    summary <- k215Summary(formatProtEnrichClusterProfilerSummaryText, project$run$selected$clusterprofilerResults,
        project$run)
    expect_true(grepl("Down-regulated: no list (no protein passed the cutoffs)", summary, fixed = TRUE),
        info = summary) # base: "Down-regulated: 0"
    # K215-I01: the summary filtered to the direction with no list says so, not a count of 0.
    filtered <- k215Summary(formatProtEnrichClusterProfilerSummaryText, project$run$selected$clusterprofilerResults,
        project$run, filter = "down")
    expect_true(grepl("Down-regulated: no list (no protein passed the cutoffs)", filtered, fixed = TRUE),
        info = filtered) # base: "Showing 0 down-regulated GO terms"
    expect_true(grepl("B_vs_A, down-regulated: no list (no protein passed the cutoffs)", k215Download(project$run),
        fixed = TRUE))
    expect_identical(k215Plot(renderProtEnrichClusterProfilerPlot, project$run$results, "down", project$run$not_tested),
        "No down-regulated enrichment data")
})

# --- the tab's fields (R1-02) ----------------------------------------------------------------------------------------

test_that("KIT-215 T12: a later run with every list tested clears the tab's not-tested field", {
    k215Offline()
    tab <- shiny::reactiveValues(analysis_complete = FALSE, all_enrichment_results = list())
    holdOmicEnrichRun(tab, k215T1()$run)
    first <- shiny::isolate(tab$enrichment_not_tested)
    expect_identical(NROW(first), 1L) # base: NULL after both runs
    expect_identical(first$direction, "up")
    holdOmicEnrichRun(tab, k215G1()$run)
    expect_null(shiny::isolate(tab$enrichment_not_tested))
    expect_null(shiny::isolate(tab$enrichment_no_list))
})

test_that("KIT-215 G1 (guard): with the complete annotation nothing is signalled and the results are as before", {
    k215Offline()
    g1 <- k215G1()
    expect_null(g1$error)
    expect_null(g1$run$not_tested)
    expect_null(g1$run$no_list)
    expect_null(g1$outcomes$enrichment_lists_not_tested)
    expect_false(any(c("not_tested", "no_list") %in% names(g1$config)))
    expect_false(any(grepl("not tested", g1$notes, fixed = TRUE)))
    tables <- lapply(g1$run$results@enrichment_data, function(contrast) lapply(contrast, function(result)
        if (is.null(result)) NULL else result@result))
    cat(sprintf("K215-G1 results %s\n", digest::digest(tables, algo = "sha256")))
    expect_true(all(vapply(g1$enricher, function(call) !call$null, NA)))
})
