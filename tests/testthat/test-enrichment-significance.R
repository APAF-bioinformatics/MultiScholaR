# KIT-203: the significance rule of the proteomics clusterProfiler route (S1 to S4, S6 to S12, S14 of the ticket,
# section 10.1; S5 and S13 are the scripts replay-k203.R and options-k203.R). Fixtures: tests/testdata/
# enrichment-significance/ (README.md: what the base commit wrote, and SHA-256) and KIT-173's go-propagation fixtures.
# No service is called. Helpers: helper-go-propagation.R.
# S11 removed by KIT-234 (the workbooks are archived).

k203Fixture <- function(...) testthat::test_path("..", "testdata", "enrichment-significance", ...)

k203SkipUnlessKernel <- function() {
    testthat::skip_if_not_installed("GO.db")
    testthat::skip_if_not_installed("AnnotationDbi")
    testthat::skip_if_not_installed("clusterProfiler")
}

#' A stand-in result with the slots the family value reads
if (!methods::isClass("k203EnrichResult")) {
    methods::setClass("k203EnrichResult", slots = c(result = "data.frame", geneSets = "list"))
}
if (!methods::isClass("k203EnrichResultNoSets")) {
    methods::setClass("k203EnrichResultNoSets", slots = c(result = "data.frame"))
}
k203Result <- function(pvalue, sets) {
    methods::new("k203EnrichResult", result = data.frame(ID = paste0("GO:", seq_along(pvalue))[seq_along(pvalue)], pvalue = pvalue,
        stringsAsFactors = FALSE), geneSets = stats::setNames(rep(list("P1"), sets), paste0("GO:", seq_len(sets))))
}

#' Silence the kernel's plot files (no assertion reads them), for the caller's lifetime
k203NoPlotFiles <- function(env = parent.frame()) {
    testthat::local_mocked_bindings(saveSelfContainedWidget = function(...) invisible(NULL), .env = env)
    testthat::local_mocked_bindings(ggsave = function(...) invisible(NULL), .package = "ggplot2", .env = env)
}

#' The S3 list (KIT-173's null draw 5) through the kernel on the full E. coli fixture, go_propagation "none"; the same
#' table as capture-k203.R's `s3Run` (KIT-189's runRoute). `...` reaches processEnrichments.
k203S3Run <- function(pathway_dir, ...) {
    list_ids <- utils::read.delim(k203Fixture("e0-accessions.tsv"), colClasses = "character")$entry
    raw <- goPropReadTsv("uniprot-ecoli-83333.tsv.gz")
    annotations <- goPropConvert(raw)
    entries <- sort(unique(raw$Entry), method = "radix")
    table <- data.frame(Protein.IDs = entries, fdr_qvalue = ifelse(entries %in% list_ids, 0.01, 0.5),
        log2FC = ifelse(entries %in% list_ids, 1, 0), stringsAsFactors = FALSE)
    da <- methods::new("da_results_for_enrichment", contrasts = tibble::tibble(contrast = "effect"),
        da_data = list(effect = table), design_matrix = data.frame(sample = "S1", stringsAsFactors = FALSE))
    results <- NULL
    utils::capture.output(results <- suppressMessages(suppressWarnings(processEnrichments(da, "999999", q_cutoff = 0.05,
        pathway_dir = pathway_dir, go_annotations = annotations, exclude_iea = FALSE, protein_id_column = "Protein.IDs",
        contrast_names = "effect", go_propagation = "none", ...))))
    list(up = results@enrichment_data$effect$up, annotations = annotations, entries = entries,
        tsv = file.path(pathway_dir, "effect_up_enrichment_results.tsv"))
}

#' KIT-173's T7 annotation fixture with T12's table; the same as capture-k203.R's `t7Run`
k203T7Run <- function(value, pathway_dir, ...) {
    roles <- goPropReadTsv("t7-fixture-entries.tsv")
    direct <- goPropConvert(goPropReadTsv("t7-annotation-ecoli.tsv"))
    entries <- sort(direct$Entry)
    table <- data.frame(Protein.IDs = entries, fdr_qvalue = 0.01,
        log2FC = ifelse(entries %in% roles$entry[roles$role == "p_term"], 3,
            ifelse(entries %in% roles$entry[roles$role == "d_term"], -3, 0)), stringsAsFactors = FALSE)
    da <- methods::new("da_results_for_enrichment", contrasts = tibble::tibble(contrast = "t7"),
        da_data = list(t7 = table), design_matrix = data.frame(sample = "S1", stringsAsFactors = FALSE))
    results <- NULL
    utils::capture.output(results <- suppressMessages(suppressWarnings(processEnrichments(da, "999999", q_cutoff = 0.05,
        pathway_dir = pathway_dir, go_annotations = direct, exclude_iea = FALSE, protein_id_column = "Protein.IDs",
        contrast_names = "t7", go_propagation = value, ...))))
    results@enrichment_data$t7
}

k203Read <- function(path) utils::read.delim(path, colClasses = "character", na.strings = character(), quote = "")

# --- S1, S2: the family value ---------------------------------------------------------------------------------------

test_that("KIT-203 S1: on KIT-189's E4 vector the family value is BH with denominator m and flags strictly", {
    got <- addEnrichmentFamilyAdjustment(k203Result(c(0.001, 0.01, 0.2), 10L), "bh_whole_family")@result
    expect_identical(got$p_adjust_family, c(0.01, 0.05, 0.6666666666666667))
    expect_identical(got$family_size, rep(10L, 3L)) # m on every row of the result
    expect_identical(got$p_adjust_family < 0.05, c(TRUE, FALSE, FALSE))
    # BH over the kept terms (enrichit's p.adjust) on the same vector, for contrast
    kept <- stats::p.adjust(c(0.001, 0.01, 0.2), "BH")
    expect_identical(kept, c(0.003, 0.015, 0.2))
    expect_identical(kept < 0.05, c(TRUE, TRUE, FALSE))
    expect_identical(names(got), c("ID", "pvalue", "p_adjust_family", "family_size"))
})

test_that("KIT-203 S2: results the family value cannot be formed for are refused; NULL, empty and q-value pass", {
    expect_error(addEnrichmentFamilyAdjustment(k203Result(c(0.01, 0.02, 0.03), 2L), "bh_whole_family"),
        class = "enrichment_family_unavailable")
    no_p <- k203Result(0.01, 10L)
    no_p@result$pvalue <- NULL
    expect_error(addEnrichmentFamilyAdjustment(no_p, "bh_whole_family"), class = "enrichment_family_unavailable")
    expect_error(addEnrichmentFamilyAdjustment(k203Result(c(0.01, NA), 10L), "bh_whole_family"),
        class = "enrichment_family_unavailable")
    no_sets <- methods::new("k203EnrichResultNoSets", result = data.frame(ID = "GO:1", pvalue = 0.01))
    expect_error(addEnrichmentFamilyAdjustment(no_sets, "bh_whole_family"), class = "enrichment_family_unavailable")
    expect_null(addEnrichmentFamilyAdjustment(NULL, "bh_whole_family"))
    empty <- k203Result(numeric(), 10L)
    expect_identical(addEnrichmentFamilyAdjustment(empty, "bh_whole_family"), empty)
    expect_identical(addEnrichmentFamilyAdjustment(no_sets, "qvalue"), no_sets)
    expect_error(checkEnrichmentSignificanceRule("bh"), class = "enrichment_significance_rule_invalid")
    expect_error(checkEnrichmentSignificanceRule(c("qvalue", "bh_whole_family")),
        class = "enrichment_significance_rule_invalid")
    expect_identical(checkEnrichmentSignificanceRule("qvalue"), "qvalue")
    da <- methods::new("da_results_for_enrichment", contrasts = tibble::tibble(contrast = "x"),
        da_data = list(x = data.frame(Protein.IDs = "P1", fdr_qvalue = 0.01, log2FC = 1)),
        design_matrix = data.frame(sample = "S1", stringsAsFactors = FALSE))
    expect_error(suppressMessages(processEnrichments(da, "999999", pathway_dir = withr::local_tempdir(),
        go_annotations = data.frame(Entry = "P1"), protein_id_column = "Protein.IDs", contrast_names = "x",
        significance_rule = "bh")), class = "enrichment_significance_rule_invalid")
})

# --- S3, S4: the route on committed fixtures ----------------------------------------------------------------------

#' The S3 runs, made once per file: the default and the q-value option
k203S3 <- local({
    cached <- list()
    function(rule) {
        if (is.null(cached[[rule]])) {
            k203NoPlotFiles()
            dir <- tempfile("k203-s3-")
            dir.create(dir)
            cached[[rule]] <<- if (identical(rule, "default")) k203S3Run(dir) else k203S3Run(dir, significance_rule = rule)
        }
        cached[[rule]]
    }
})

test_that("KIT-203 S3: KIT-173's null draw 5 calls nothing by default and 120 terms by the q-value; m counted apart", {
    k203SkipUnlessKernel()
    run <- k203S3("default")
    tsv <- k203Read(run$tsv)
    expect_identical(nrow(tsv), nrow(run$up@result))
    expect_identical(sum(tsv$significant %in% "TRUE"), 0L)
    expect_identical(unique(tsv$significance_rule), "bh_whole_family")
    p <- run$up@result$pvalue
    m <- length(run$up@geneSets)
    expect_identical(run$up@result$p_adjust_family, stats::p.adjust(c(p, rep(1, m - length(p))), "BH")[seq_along(p)])
    # m counted apart: the sets the route builds, within the universe enrichit takes (annotated proteins of the
    # background), 10 to 500 members
    sets <- buildUniprotGoGeneSets(run$annotations, "none")$term2gene
    pairs <- unique(data.frame(term = sets[[1L]], gene = sets[[2L]], stringsAsFactors = FALSE))
    universe <- intersect(run$entries, pairs$gene)
    sizes <- table(pairs$term[pairs$gene %in% universe])
    expect_identical(unique(run$up@result$family_size), as.integer(sum(sizes >= 10L & sizes <= 500L)))
    expect_identical(m, as.integer(sum(sizes >= 10L & sizes <= 500L)))
    # nesting (5.1): whole-family calls within kept-family BH calls within q-value calls
    family <- run$up@result$p_adjust_family < 0.05
    kept <- run$up@result$p.adjust < 0.05
    expect_true(all(!family | kept))
    expect_true(all(!kept | (!is.na(run$up@result$qvalue) & run$up@result$qvalue < 0.05)))
    q <- k203Read(k203S3("qvalue")$tsv)
    expect_identical(sum(q$significant %in% "TRUE"), 120L)
})

test_that("KIT-203 S4: the q-value option writes and stores what the base commit did, byte for byte", {
    k203SkipUnlessKernel()
    check <- function(result, tsv, stem) {
        base <- readRDS(k203Fixture(paste0(stem, ".rds")))
        expect_identical(digest::digest(file = tsv, algo = "sha256"),
            digest::digest(file = k203Fixture(paste0(stem, ".tsv")), algo = "sha256"), info = stem)
        expect_identical(result@result, base$result, info = stem)
        expect_identical(attr(result, "enrichment_plot_inputs", exact = TRUE), base$plot_inputs, info = stem)
    }
    run <- k203S3("qvalue")
    check(run$up, run$tsv, "s3-none-up")
    k203NoPlotFiles()
    for (value in c("none", "is_a_part_of_regulates")) {
        dir <- withr::local_tempdir()
        results <- k203T7Run(value, dir, significance_rule = "qvalue")
        short <- if (identical(value, "none")) "none" else "iapor"
        for (direction in c("up", "down")) {
            check(results[[direction]], file.path(dir, paste0("t7_", direction, "_enrichment_results.tsv")),
                paste("t7", short, direction, sep = "-"))
        }
    }
})

# --- S6, S7: declaration, identity, outcome, notice -------------------------------------------------------------

#' KIT-173's T7 project (helper goPropRunProject's steps) with the run's inputs overridden by `extra` and the notices
#' kept: `list(run, record, outcomes, notes)`.
k203RunProject <- function(root, annotations, extra = list(), process_fn = processEnrichments) {
    roles <- goPropReadTsv("t7-fixture-entries.tsv")
    paths <- goPropPaths(root)
    for (p in paths) dir.create(p, recursive = TRUE, showWarnings = FALSE)
    saveRDS(annotations, file.path(paths$source_dir, "uniprot_dat_cln.RDS"))
    manager <- newDIAStateAdapter(newWorkflowEntityStore(root, "go-propagation", "DIA"),
        WorkflowState$new(audit_enabled = FALSE))
    on.exit(closeDIAStateAdapter(manager), add = TRUE)
    manager$.context$reference_states <- TRUE
    manager$setWorkflowType("DIA")
    entries <- sort(roles$entry)
    design <- data.frame(Run = paste0("S", 1:6), group = rep(c("A", "B"), each = 3L),
        batch = rep(c("b1", "b2", "b3"), times = 2L), stringsAsFactors = FALSE)
    quant <- goPropQuant(entries, roles$entry[roles$role == "p_term"], roles$entry[roles$role == "d_term"])
    object <- methods::new("ProteinQuantitativeData", protein_quant_table = quant, design_matrix = design,
        sample_id = "Run", group_id = "group", protein_id_column = "uniprot_acc",
        protein_id_table = data.frame(uniprot_acc = entries, gene_names = paste0("GENE", seq_along(entries)),
            stringsAsFactors = FALSE),
        args = list(globalParameters = list(workflow_type = "DIA", use_limpa = FALSE, report_template = "DIA_report"),
            deAnalysisParameters = list(formula_string = "~ 0 + group")))
    setDIAStepProvenance(manager, "design_builder")
    manager$saveState("correlation_filtered", object, list(), "fixture")
    manager$.context$experiment_paths <- paths
    contrasts <- data.frame(contrasts = "groupB-groupA", full_format = "B_vs_A=groupB-groupA",
        friendly_names = "B_vs_A", stringsAsFactors = FALSE)
    workflow <- list2env(list(state_manager = manager, contrasts_tbl = contrasts, tab_status = list(),
        qc_params = list(protein_qc = list()), config_list = list()))
    parameters <- utils::modifyList(list(selected_contrast = "B_vs_A", organism_taxid = "999999", up_cutoff = 0,
        down_cutoff = 0, q_cutoff = 0.05, correction_method = "gSCS", enable_organism_filter = FALSE), extra)
    notes <- character()
    run <- NULL
    withDIAStateScope(manager, function() {
        da_data <- new.env(parent = emptyenv())
        da_data$current_s4_object <- getDIAStateReference(manager)
        da_data$contrasts_available <- contrasts$contrasts
        manifest <- runDIAReferenceDaAnalysis(list(formula_string = "~ 0 + group", da_q_val_thresh = 0.05,
            treat_lfc_cutoff = 0), da_data, workflow, paths, contrasts = contrasts, publish = FALSE)
        produceDIADaPublication(manager, manifest, annotations, 0.05, paths[c("da_output_dir", "publication_graphs_dir")],
            save_assets = FALSE)
        run <<- runWorkflowEnrichment(workflow, omicGuiConfig("proteomics")$enrichment, manifest, parameters, paths,
            process_fn = process_fn, notify = function(message, type, duration)
                notes <<- c(notes, paste0("[", type, "] ", message)))
    })
    record <- file.path(root, "analysis_record.json")
    writeWorkflowAnalysisRecord(manager, record)
    record_list <- jsonlite::read_json(record, simplifyVector = FALSE)
    step <- Filter(function(s) identical(s$operation, "enrichment"), record_list$steps)[[1L]]
    list(run = run, record = record, step = step, notes = notes)
}

#' The projects, made once per file, keyed by their overrides
k203Project <- local({
    cached <- list()
    function(key, extra = list(), process_fn = processEnrichments) {
        if (is.null(cached[[key]])) {
            work <- tempfile(paste0("k203-", key, "-"))
            dir.create(work)
            annotations <- goPropAnnotate(goPropT7Raw(), file.path(work, "cache"))
            result <- NULL
            withr::with_envvar(c(MULTISCHOLAR_TEST_MODE = "false"), utils::capture.output(result <- suppressMessages(
                suppressWarnings(k203RunProject(file.path(work, "project"), annotations, extra, process_fn)))))
            cached[[key]] <<- result
        }
        cached[[key]]
    }
})

test_that("KIT-203 S6: the default declares its rule; the q-value declares KIT-173's ten as the base did; outcome", {
    k203SkipUnlessKernel()
    names_of <- function(step) vapply(step$parameters, `[[`, "", "name")
    ten <- c(utils::head(diaEnrichmentOperation()$parameters, 9L), "go_propagation")
    default <- k203Project("default")
    # KIT-205: every record step lists the fixed background rule marker last.
    expect_identical(names_of(default$step), c(ten, "significance_rule", "enrichment_background"))
    value <- Filter(function(p) identical(p$name, "significance_rule"), default$step$parameters)[[1L]]$value
    expect_identical(decodeWorkflowRecordValue(value), "bh_whole_family")
    expect_identical(utils::tail(names(default$step$outcomes), 1L), "significance_rule")
    expect_identical(default$step$outcomes$significance_rule, "bh_whole_family")
    # The q-value run: the base run's declared names and values (its identity is formed from them with the store's
    # entity ids, which are store-local, so the values are compared without `results`).
    qvalue <- k203Project("qvalue", list(significance_rule = "qvalue"))
    base <- readRDS(k203Fixture("s6-identity.rds"))
    # KIT-205: the marker last; the names before it are the base run's (the stored fixture is not re-recorded).
    expect_identical(names_of(qvalue$step), c(ten, "enrichment_background"))
    expect_identical(setdiff(names_of(qvalue$step), "enrichment_background"), base$declared_names)
    without_results <- function(x) x[setdiff(names(x), "results")]
    expect_identical(without_results(qvalue$run$parameters), without_results(base$parameters))
    expect_identical(names(qvalue$run$parameters), names(base$parameters))
    expect_identical(qvalue$step$outcomes$significance_rule, "qvalue")
    expect_identical(setdiff(names(qvalue$step$outcomes), "significance_rule"), base$outcome_names)
    # The g:Profiler route declares and passes nothing new, and records no rule.
    expect_identical(diaEnrichmentDeclaredNames(diaEnrichmentOperation()$parameters, list()),
        utils::head(diaEnrichmentOperation()$parameters, 9L))
    expect_identical(diaEnrichmentDeclaredNames(diaEnrichmentOperation()$parameters,
        list(significance_rule = "qvalue")), utils::head(diaEnrichmentOperation()$parameters, 9L))
    gprofiler <- k203Project("gprofiler", list(organism_taxid = "9606"), runProtEnrichDeterministicProcessEnrichments)
    expect_identical(names_of(gprofiler$step), c(utils::head(diaEnrichmentOperation()$parameters, 9L),
        "enrichment_background")) # KIT-205
    expect_false("significance_rule" %in% names(gprofiler$step$outcomes))
    expect_null(buildProtEnrichProcessEnrichmentsArgs(NULL, "9606", 0, 0, 0.05, "p", NULL, "uniprot_acc", "B_vs_A",
        "gSCS")$processArgs$significance_rule)
})

test_that("KIT-203 S7: every clusterProfiler run is told its rule; the correction sentence follows only for gSCS", {
    k203SkipUnlessKernel()
    rules <- ENRICHMENT_SIGNIFICANCE_RULES
    unused <- " The selected correction, gSCS, applies to g:Profiler only and was not used."
    double <- runProtEnrichDeterministicProcessEnrichments
    for (rule in rules$rule) {
        notice <- rules$notice[rules$rule == rule]
        gscs <- k203Project(paste0("notice-gscs-", rule), list(significance_rule = rule), double)
        expect_identical(gscs$notes, paste0("[warning] ", notice, unused), info = rule)
        fdr <- k203Project(paste0("notice-fdr-", rule), list(significance_rule = rule, correction_method = "fdr"), double)
        expect_identical(fdr$notes, paste0("[message] ", notice), info = rule)
    }
    gprofiler <- k203Project("gprofiler", list(organism_taxid = "9606"), double)
    expect_false(any(grepl("called significant", gprofiler$notes, fixed = TRUE)))
})

# --- S8: the plot ---------------------------------------------------------------------------------------------------

test_that("KIT-203 S8: the default draws the family value; a result without the rule draws the q-value as before", {
    result <- methods::new("k203EnrichResult", result = data.frame(ID = paste0("GO:", 1:7),
        Description = paste("Term", 1:7), qvalue = c(0.001, 0.002, 0.003, 0.004, 0.01, 0.02, 0.6),
        p_adjust_family = c(0.002, 0.03, 0.04, 0.045, 0.049, 0.2, 0.9), Count = 2:8, GeneRatio = "2/10",
        BgRatio = "5/100", stringsAsFactors = FALSE), geneSets = list())
    go <- data.frame(ID = paste0("GO:", 1:7), source = "GO:BP", term = paste("Name", 1:7))
    family <- setEnrichmentPlotInputs(result, go, 0.05, significance_rule = "bh_whole_family")
    plot <- suppressWarnings(buildClusterProfilerEnrichmentGgplot(family, "c", "up"))
    expect_equal(plot$data$neg_log10_q, -log10(result@result$p_adjust_family))
    expect_identical(ggplot2::layer_data(plot, 3)$label, paste("Name", c(1, 2, 3, 4, 5)))
    expect_equal(ggplot2::layer_data(plot, 1)$yintercept, -log10(0.05))
    expect_identical(plot$labels$y, "-log10(adjusted p-value, all terms of 10 to 500 proteins)")
    tooltips <- rlang::eval_tidy(plot$mapping$text, plot$data)
    expect_true(all(startsWith(sub("^.*\n", "", tooltips), "Adjusted p (all terms): ")))
    # A stored result without the rule (every result stored before KIT-203, and every q-value run)
    plain <- setEnrichmentPlotInputs(result, go, 0.05)
    expect_identical(names(attr(plain, "enrichment_plot_inputs")), c("go_terms", "q_cutoff"))
    old <- suppressWarnings(buildClusterProfilerEnrichmentGgplot(plain, "c", "up"))
    expect_equal(old$data$neg_log10_q, -log10(result@result$qvalue))
    expect_identical(ggplot2::layer_data(old, 3)$label, paste("Name", 1:5))
    expect_identical(old$labels$y, "-log10(q-value)")
    # The base commit's stored S3 result is drawn exactly as the q-value run's stored result on this commit.
    k203SkipUnlessKernel()
    base <- readRDS(k203Fixture("s3-none-up.rds"))
    stored <- methods::new("k203EnrichResult", result = base$result, geneSets = list())
    attr(stored, "enrichment_plot_inputs") <- base$plot_inputs
    display <- function(x) {
        plots <- suppressWarnings(buildEnrichmentPlot(x, "effect", "up"))
        built <- plotly::plotly_build(plots$interactive)
        list(layers = lapply(seq_along(plots$static$layers), function(i) ggplot2::layer_data(plots$static, i)),
            traces = built$x$data, layout = built$x$layout)
    }
    expect_identical(display(stored), display(k203S3("qvalue")$up))
})

# --- S9: the report paragraph ---------------------------------------------------------------------------------------

#' A template's `FormatGOEnrichment` (as test-go-propagation.R's goPropFormatter)
k203Formatter <- function(template) {
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

test_that("KIT-203 S9: the four reports state the recorded rule; earlier records keep their paragraph and are told", {
    rules <- ENRICHMENT_SIGNIFICANCE_RULES
    method_of <- function(rule) sub("{cutoff}", "0.05", rules$methods[rules$rule == rule], fixed = TRUE)
    numbers <- function(source, extra = list()) c(list(ui_enrichment_organism_selected = "999999",
        ui_enrichment_up_lfc = 0, ui_enrichment_down_lfc = 0, ui_enrichment_q_value_cutoff = 0.05),
        if (!is.null(source)) list(ui_enrichment_database_source = source), extra)
    # KIT-205 (SD-1 (b)): "tested"; these numbers carry no record (the background rule's no-record branch).
    family_clause <- paste0("terms with 10 to 500 annotated tested proteins, whether or not a listed protein falls ",
        "in them (terms without one enter with p = 1).")
    old_flag <- paste0("A term was called significant when its q-value (Storey, 2002) was below 0.05; ",
        "Benjamini-Hochberg adjusted p-values are also given in the results table.")
    for (template in c("LFQ_report.rmd", "TMT_report.rmd", "DIANN_report.rmd", "DIANN_limpa_report.rmd")) {
        format <- k203Formatter(template)
        text <- as.character(format(numbers("clusterprofiler", list(enrichment_significance_rule = "bh_whole_family"))))
        expect_true(grepl(family_clause, text, fixed = TRUE), info = template)
        expect_true(grepl(method_of("bh_whole_family"), text, fixed = TRUE), info = template)
        expect_false(grepl(old_flag, text, fixed = TRUE), info = template)
        expect_false(grepl(method_of("qvalue"), text, fixed = TRUE), info = template)
        for (earlier in list(list(), list(enrichment_significance_rule = "qvalue"),
            list(enrichment_significance_rule = "not recorded"))) {
            text <- as.character(format(numbers("clusterprofiler", earlier)))
            expect_true(grepl(paste0("and at least one protein in the list. ", old_flag, " The background was"), text,
                fixed = TRUE), info = template)
            expect_true(grepl(paste0("annotation. ", method_of("qvalue")), text, fixed = TRUE), info = template)
            expect_false(grepl(family_clause, text, fixed = TRUE), info = template)
        }
        for (source in list("gprofiler2", "N/A", NULL)) {
            text <- as.character(format(numbers(source, list(enrichment_significance_rule = "bh_whole_family"))))
            expect_false(grepl(method_of("bh_whole_family"), text, fixed = TRUE), info = template)
            expect_false(grepl(method_of("qvalue"), text, fixed = TRUE), info = template)
        }
    }
})

# --- S10: the tab -----------------------------------------------------------------------------------------------------

#' Where an input sits in a UI tree: "outside" every conditionalPanel, "inside" one, or NULL when absent
k203Where <- function(tag, id, inside = FALSE) {
    if (inherits(tag, "shiny.tag")) {
        inside <- inside || !is.null(tag$attribs[["data-display-if"]])
        if (identical(tag$attribs$id, id)) return(if (inside) "inside" else "outside")
        tag <- tag$children
    }
    if (is.list(tag)) for (child in tag) {
        found <- k203Where(child, id, inside)
        if (!is.null(found)) return(found)
    }
    NULL
}

#' The output ids of a UI tree (shiny's text outputs and the HTML widgets' outputs)
k203Outputs <- function(tag) {
    if (inherits(tag, "shiny.tag")) {
        class <- paste(unlist(tag$attribs[names(tag$attribs) == "class"]), collapse = " ")
        here <- if (grepl("shiny-[a-z-]*output|html-widget-output", class)) tag$attribs$id
        return(c(here, k203Outputs(tag$children)))
    }
    if (is.list(tag)) return(unlist(lapply(tag, k203Outputs)))
    NULL
}

#' Every conditionalPanel of a UI tree: its condition and its text
k203Panels <- function(tag) {
    if (inherits(tag, "shiny.tag")) {
        condition <- tag$attribs[["data-display-if"]]
        here <- if (is.null(condition)) list() else list(list(condition = condition, text = as.character(tag)))
        return(c(here, k203Panels(tag$children)))
    }
    if (is.list(tag)) return(do.call(c, c(list(list()), lapply(tag, k203Panels))))
    list()
}

test_that("KIT-203 S10: the tab offers both rules before Run, each with its note; the run's rule is shown", {
    config <- omicGuiConfig("proteomics")$enrichment
    ui <- mod_omic_enrich_ui("enrichment_analysis", config)
    id <- "enrichment_analysis-significance_rule"
    expect_identical(k203Where(ui, id), "outside")
    html <- as.character(ui)
    expect_true(grepl('<option value="bh_whole_family" selected>', html, fixed = TRUE))
    expect_true(grepl('<option value="qvalue">', html, fixed = TRUE))
    expect_identical(lengths(regmatches(html, gregexpr("<option value=\"(bh_whole_family|qvalue)\"", html))), 2L)
    rules <- ENRICHMENT_SIGNIFICANCE_RULES
    panels <- k203Panels(ui)
    for (rule in rules$rule) {
        condition <- sprintf("input['%s'] == '%s'", id, rule)
        hit <- Filter(function(p) identical(p$condition, condition), panels)
        expect_length(hit, 1L)
        expect_true(grepl(htmltools::htmlEscape(rules$note[rules$rule == rule]), hit[[1L]]$text, fixed = TRUE),
            info = rule)
        expect_false(grepl(htmltools::htmlEscape(rules$note[rules$rule != rule]), hit[[1L]]$text, fixed = TRUE),
            info = rule)
    }
    expect_true(grepl("p_adjust_family: Benjamini-Hochberg over every term of 10 to 500 tested proteins.", html,
        fixed = TRUE)) # KIT-205 (SD-1 (b))
    expect_true(grepl("Significance cutoff:", html, fixed = TRUE))
    # No output is added: the module's output ids are the base commit's.
    expect_identical(sort(unique(k203Outputs(ui))), c(
        "enrichment_analysis-analysis_method_display", "enrichment_analysis-clusterprofiler_plot",
        "enrichment_analysis-clusterprofiler_results_table", "enrichment_analysis-clusterprofiler_summary_stats",
        "enrichment_analysis-contrasts_display", "enrichment_analysis-enrichment_status",
        "enrichment_analysis-gprofiler_plot", "enrichment_analysis-gprofiler_results_table",
        "enrichment_analysis-gprofiler_summary_stats", "enrichment_analysis-stringdb_plot",
        "enrichment_analysis-stringdb_results_table", "enrichment_analysis-stringdb_summary_stats"))
    # An omic whose services lack clusterProfiler has no selector.
    no_cp <- config
    no_cp$services <- setdiff(config$services, "clusterprofiler")
    expect_null(k203Where(mod_omic_enrich_ui("enrichment_analysis", no_cp), id))

    # The table rounds the family value.
    rounded <- NULL
    renderProtEnrichClusterProfilerResultsTable(
        clusterprofilerResults = data.frame(ID = "GO:1", pvalue = 0.01, p.adjust = 0.02, qvalue = 0.02,
            p_adjust_family = 0.03),
        renderDtFn = function(expr, ...) expr, datatableFn = function(data, ...) data,
        formatRoundFn = function(table, columns, digits) rounded <<- columns)
    expect_identical(rounded, c("pvalue", "p.adjust", "qvalue", "p_adjust_family"))

    # The completed status names the rule on the clusterProfiler method only.
    status <- function(method, rule) formatProtEnrichStatusText(TRUE, list(method = method, species_name = "E. coli"),
        "B_vs_A", 0, 0, 0.05, significanceRule = rule)
    label <- function(rule) rules$label[rules$rule == rule]
    expect_true(grepl(paste0("Significance: ", label("bh_whole_family")), status("clusterprofiler", "bh_whole_family"),
        fixed = TRUE))
    expect_identical(status("gprofiler2", "bh_whole_family"), formatProtEnrichStatusText(TRUE,
        list(method = "gprofiler2", species_name = "E. coli"), "B_vs_A", 0, 0, 0.05))
    expect_false(grepl("Significance:", status("gprofiler2", "bh_whole_family"), fixed = TRUE))

    # The download summary: the run's rule on the clusterProfiler method; a run recorded before the rule is the q-value.
    summary_of <- function(method, rule) {
        written <- NULL
        writeProtEnrichResultsDownloadArchive(file = tempfile(fileext = ".zip"), selectedContrast = "B_vs_A",
            methodInfo = list(method = method, species_name = "E. coli"), organismTaxid = "83333", upCutoff = 0,
            downCutoff = 0, qCutoff = 0.05, tempDir = withr::local_tempdir(.local_envir = parent.frame()),
            writeLinesFn = function(text, path) if (grepl("summary", path)) written <<- text,
            zipFn = function(...) invisible(NULL), significanceRule = rule)
        written
    }
    expect_true(grepl(paste0("Significance cutoff: 0.05\nSignificance: ", label("bh_whole_family")),
        summary_of("clusterprofiler", "bh_whole_family"), fixed = TRUE))
    expect_true(grepl(paste0("Significance: ", label("qvalue")), summary_of("clusterprofiler", NULL), fixed = TRUE))
    # KIT-226 (VAL-196): the g:Profiler cutoff applies to the service's adjusted p-value.
    expect_true(grepl("Adjusted p-value cutoff: 0.05", summary_of("gprofiler2", NULL), fixed = TRUE))
    expect_false(grepl("Significance", summary_of("gprofiler2", NULL), fixed = TRUE))

    # In the server: the status reads the run's rule, not the selector; a restored run without a rule is the q-value;
    # the Run passes the selector's value.
    shiny::testServer(function(input, output, session) {
        tab <- shiny::reactiveValues(analysis_complete = TRUE, gprofiler_results = NULL, clusterprofiler_results = NULL,
            stringdb_results = NULL, contrasts_available = "B_vs_A")
        ran <- shiny::reactiveVal(list(up_cutoff = 0, down_cutoff = 0, q_cutoff = 0.05, significance_rule = "qvalue"))
        method <- shiny::reactive(list(method = "clusterprofiler", species_name = "E. coli"))
        registerOmicEnrichStatusOutputs(output, input, tab, method, ran)
    }, {
        session$setInputs(significance_rule = "bh_whole_family", selected_contrast = "B_vs_A")
        expect_true(grepl(paste0("Significance: ", label("qvalue")), output$enrichment_status, fixed = TRUE))
    })
    testthat::local_mocked_bindings(restoreProtEnrichResultsFromState = function(...) TRUE)
    shiny::testServer(function(input, output, session) {
        tab <- shiny::reactiveValues(analysis_complete = FALSE, all_enrichment_results = NULL)
        ran <- shiny::reactiveVal(NULL)
        manager <- list(.context = list(snapshot = list(step_provenance = list(enrichment_completed = list(
            values = list(up_cutoff = 0, down_cutoff = 0, q_cutoff = 0.05, correction_method = "BH"))))))
        workflow <- shiny::reactiveValues(state_manager = manager)
        registerOmicEnrichSelectedContrastObserver(input, tab, workflow, list(state_name = "enrichment_completed"), ran)
    }, {
        session$setInputs(selected_contrast = "B_vs_A", significance_rule = "bh_whole_family")
        expect_identical(ran()$significance_rule, "qvalue")
        expect_identical(ran()$q_cutoff, 0.05)
    })
    passed <- NULL
    testthat::local_mocked_bindings(withProgress = function(expr, ...) expr, .package = "shiny")
    testthat::local_mocked_bindings(runWorkflowEnrichment = function(workflow, config, manifest, parameters, ...) {
        passed <<- parameters
        list(saved = list(saved = TRUE), selected = list(found = FALSE), method = list(method = "clusterprofiler"),
            contrasts = list(), results = NULL, service = character(), reference = NULL)
    }, omicEnrichManifestCurrent = function(...) TRUE) # KIT-216 (r1.4, F05): item 13's DA check is E15's subject
    shiny::testServer(function(input, output, session) {
        tab <- shiny::reactiveValues(da_results_data = list(source = "s"), analysis_complete = FALSE)
        ran <- shiny::reactiveVal(NULL)
    }, {
        session$setInputs(selected_contrast = "B_vs_A", organism_taxid = "83333", up_cutoff = 0, down_cutoff = 0,
            q_cutoff = 0.05, correction_method = "gSCS", enable_organism_filter = FALSE,
            significance_rule = "qvalue")
        runOmicEnrichAnalysis(input, tab, shiny::reactiveValues(), list(), config, ran)
        expect_identical(passed$significance_rule, "qvalue")
        expect_identical(ran()$significance_rule, "qvalue")
    })
})

# --- S12: record ---

test_that("KIT-203 S12: the record labels the rule and the report reads it, declared or as the outcome", {
    expect_identical(workflowParameterLabels$enrichment$significance_rule,
        "Significance rule for Gene Ontology terms (clusterProfiler route)")
    expect_identical(workflowParameterLabels$enrichment$q_cutoff,
        "Enrichment cutoff (differential-abundance q-value of the listed proteins; significance value of a term)")
    write_record <- function(parameters, outcomes) {
        path <- tempfile(fileext = ".json")
        jsonlite::write_json(list(steps = list(list(operation = "enrichment",
            parameters = lapply(names(parameters), function(name) list(name = name,
                value = encodeWorkflowRecordValue(parameters[[name]]))), outcomes = outcomes))), path,
            auto_unbox = TRUE, null = "null")
        path
    }
    declared <- workflowRecordReportParameters(write_record(list(significance_rule = "bh_whole_family"),
        list(significance_rule = "bh_whole_family")), list())
    expect_identical(declared$enrichment_significance_rule, "bh_whole_family")
    outcome <- workflowRecordReportParameters(write_record(list(q_cutoff = 0.05), list(significance_rule = "qvalue")),
        list())
    expect_identical(outcome$enrichment_significance_rule, "qvalue")
    neither <- workflowRecordReportParameters(write_record(list(q_cutoff = 0.05), list(rows = 1)), list())
    expect_identical(neither$enrichment_significance_rule, "not recorded")
})

# --- S14: the replay default ----------------------------------------------------------------------------------------

test_that("KIT-203 S14: a record without the rule replays with the q-value; a recorded rule replays as recorded", {
    as_recorded <- function(value) decodeWorkflowRecordValue(jsonlite::fromJSON(jsonlite::toJSON(
        encodeWorkflowRecordValue(value), auto_unbox = TRUE, digits = I(17), null = "null", na = "null"),
        simplifyVector = FALSE))
    recorded <- lapply(list(taxon_id = 999999, up_cutoff = 0, down_cutoff = 0, q_cutoff = 0.05, exclude_iea = FALSE,
        correction_method = "BH", protein_id_column = "uniprot_acc", contrast_names = "B_vs_A",
        results = c(B_vs_A = "e")), as_recorded)
    passed <- NULL
    testthat::local_mocked_bindings(
        produceDIADaPublication = function(...) invisible(NULL),
        runWorkflowEnrichment = function(workflow, config, manifest, parameters, paths, ...) {
            passed <<- parameters
            produced <- recorded
            if (!is.null(parameters$significance_rule) && !identical(parameters$significance_rule, "qvalue")) {
                produced$significance_rule <- parameters$significance_rule
            }
            list(parameters = produced, service = character(), organism_filter = list())
        })
    replay <- function(parameters) replayWorkflowEnrichmentStep(NULL, list(state = "enrichment_completed",
        outcomes = list()), parameters, list(parameters = list(qvalue = 0.05)),
        data.frame(friendly_names = "B_vs_A", stringsAsFactors = FALSE), list(), withr::local_tempdir())
    replay(recorded)
    expect_identical(passed$significance_rule, "qvalue")
    replay(c(recorded, list(significance_rule = "bh_whole_family")))
    expect_identical(passed$significance_rule, "bh_whole_family")
})
