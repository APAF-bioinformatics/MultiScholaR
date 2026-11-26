# KIT-205 (VAL-129 (2)): the enrichment background (universe) of each comparison is the proteins tested in that
# comparison (a finite `fdr_qvalue`, and a finite `raw_pvalue` where the table has the column), on both routes; a record
# made before the rule replays with every row (replay rule marker `enrichment_background`). B1 to B9 of the ticket,
# section 10.1. No service is called: the g:Profiler route's `perform_enrichment` is replaced by a recorder; the
# clusterProfiler route runs `clusterProfiler::enricher` on KIT-173's committed E. coli fixture (helper-go-propagation.R).

k205SkipUnlessKernel <- function() {
    testthat::skip_if_not_installed("GO.db")
    testthat::skip_if_not_installed("AnnotationDbi")
    testthat::skip_if_not_installed("clusterProfiler")
}

#' Silence the kernel's plot files (no assertion reads them), for the caller's lifetime
k205NoPlotFiles <- function(env = parent.frame()) {
    testthat::local_mocked_bindings(saveSelfContainedWidget = function(...) invisible(NULL), .env = env)
    testthat::local_mocked_bindings(ggsave = function(...) invisible(NULL), .package = "ggplot2", .env = env)
}

#' The fixture's proteins, made once per file from KIT-173's full E. coli annotation (UniProt 2026_04): a term `term`
#' with exactly 10 of the 120 proteins (`members`), and 110 proteins outside it (`others`): an up list of 15 that share
#' the pool's most frequent term, a down list of 10, and 85 fillers. `pairs` are the direct annotations (term, gene).
k205Ids <- local({
    cached <- NULL
    function() {
        if (is.null(cached)) {
            annotations <- goPropConvert(goPropReadTsv("uniprot-ecoli-83333.tsv.gz"))
            sets <- buildUniprotGoGeneSets(annotations, "none")$term2gene
            pairs <- unique(data.frame(term = as.character(sets[[1L]]), gene = as.character(sets[[2L]]),
                stringsAsFactors = FALSE))
            sizes <- table(pairs$term)
            term <- sort(names(sizes)[sizes >= 10L & sizes <= 40L], method = "radix")[[1L]]
            in_term <- sort(unique(pairs$gene[pairs$term == term]), method = "radix")
            pool <- sort(setdiff(unique(pairs$gene), in_term), method = "radix")
            pool_pairs <- pairs[pairs$gene %in% pool, , drop = FALSE]
            counts <- table(pool_pairs$term)
            big <- names(counts)[order(-as.integer(counts), names(counts), method = "radix")][[1L]]
            up <- sort(unique(pool_pairs$gene[pool_pairs$term == big]), method = "radix")[seq_len(15L)]
            rest <- setdiff(pool, up)
            down <- rest[seq(1L, by = 5L, length.out = 10L)]
            fill <- setdiff(rest[seq(3L, by = 5L, length.out = 120L)], down)[seq_len(85L)]
            cached <<- list(annotations = annotations, pairs = pairs, term = term, members = in_term[seq_len(10L)],
                up = up, down = down, fill = fill, others = c(up, down, fill))
        }
        cached
    }
})

#' The 12 untested proteins of the main variant: 6 of the 10-member term and 6 fillers spread through the table
k205Untested <- function() {
    ids <- k205Ids()
    c(ids$members[seq_len(6L)], ids$fill[c(5L, 19L, 33L, 47L, 61L, 75L)])
}

#' One comparison's DA table on the 120 proteins: `raw_pvalue` uniform (seeded), the planted lists' p-values small,
#' `fdr_qvalue` by Benjamini-Hochberg over the finite p-values (as the package assigns it), `log2FC` +1 for the up
#' list, -1 for the down list, 0 otherwise. The `untested` proteins have p, q and log2FC NA.
k205Table <- function(untested = character(), raw_p = TRUE, seed = 205L) {
    ids <- k205Ids()
    accessions <- c(ids$members, ids$others)
    p <- withr::with_seed(seed, stats::runif(length(accessions)))
    planted <- c(ids$up, ids$down)
    p[match(planted, accessions)] <- seq(1e-6, 1e-5, length.out = length(planted))
    log2fc <- ifelse(accessions %in% ids$up, 1, ifelse(accessions %in% ids$down, -1, 0))
    p[accessions %in% untested] <- NA_real_
    log2fc[accessions %in% untested] <- NA_real_
    q <- rep(NA_real_, length(p))
    q[is.finite(p)] <- stats::p.adjust(p[is.finite(p)], method = "BH")
    table <- data.frame(Protein.IDs = accessions, raw_pvalue = p, fdr_qvalue = q, log2FC = log2fc,
        stringsAsFactors = FALSE)
    if (!raw_p) table$raw_pvalue <- NULL
    table
}

k205Da <- function(tables) {
    methods::new("da_results_for_enrichment", contrasts = tibble::tibble(contrast = names(tables)), da_data = tables,
        design_matrix = data.frame(sample = "S1", stringsAsFactors = FALSE))
}

#' The clusterProfiler route (taxon 999999, direct annotations) on `tables`; `list(results, dir)`
k205Cluster <- function(tables, env = parent.frame()) {
    k205NoPlotFiles(env)
    dir <- withr::local_tempdir(.local_envir = env)
    results <- NULL
    utils::capture.output(results <- suppressMessages(suppressWarnings(processEnrichments(k205Da(tables), "999999",
        q_cutoff = 0.05, pathway_dir = dir, go_annotations = k205Ids()$annotations, exclude_iea = FALSE,
        protein_id_column = "Protein.IDs", contrast_names = names(tables), go_propagation = "none"))))
    list(results = results, dir = dir)
}

#' The g:Profiler route (taxon 9606) on `tables`, all comparisons in one call (so a mask pooled over comparisons would
#' show), with `perform_enrichment` recording each call: the ids sent and the custom background. The route takes the
#' comparisons in order, up then down, and every table here has both lists. Returns the calls keyed
#' "<contrast>/<direction>".
k205Gprofiler <- function(tables, id_column = "Protein.IDs") {
    acc <- new.env()
    acc$calls <- list()
    testthat::local_mocked_bindings(perform_enrichment = function(data_subset, custom_bg, protein_id_column, ...) {
        direction <- if (all(data_subset$log2FC > 0)) "up" else "down"
        acc$calls[[length(acc$calls) + 1L]] <- list(direction = direction,
            ids = as.character(data_subset[[protein_id_column]]), custom_bg = custom_bg)
        NULL
    })
    utils::capture.output(suppressMessages(suppressWarnings(processEnrichments(k205Da(tables), "9606",
        q_cutoff = 0.05, pathway_dir = withr::local_tempdir(), exclude_iea = FALSE, protein_id_column = id_column,
        contrast_names = names(tables)))))
    stopifnot(length(acc$calls) == 2L * length(tables))
    contrasts <- rep(names(tables), each = 2L)
    stats::setNames(acc$calls, vapply(seq_along(acc$calls), function(i) paste0(contrasts[[i]], "/",
        acc$calls[[i]]$direction), ""))
}

#' The tested rows of a table by the rule's definition (VAL-129 (2)), for the expectations
k205Tested <- function(table) {
    tested <- is.finite(table$fdr_qvalue)
    if ("raw_pvalue" %in% names(table)) tested <- tested & is.finite(table$raw_pvalue)
    table[tested, , drop = FALSE]
}

#' Terms of 10 to 500 members among `ids`' annotated proteins, counted apart from the kernel (KIT-203 S3's count)
k205FamilySize <- function(ids) {
    pairs <- k205Ids()$pairs
    universe <- intersect(ids, pairs$gene)
    sizes <- table(pairs$term[pairs$gene %in% universe])
    as.integer(sum(sizes >= 10L & sizes <= 500L))
}

# --- B1, B2: the background is the tested rows, per comparison, on both routes ---------------------------------------

test_that("KIT-205 B1: the clusterProfiler universe excludes the untested proteins; the family counts tested proteins", {
    k205SkipUnlessKernel()
    untested <- k205Untested()
    table <- k205Table(untested)
    up <- k205Cluster(list(k205 = table))$results@enrichment_data$k205$up
    expect_false(is.null(up))
    tested_annotated <- intersect(k205Tested(table)$Protein.IDs, k205Ids()$pairs$gene)
    # Base: the universe holds the 12 annotated untested proteins (length 108 + 12).
    expect_identical(sort(intersect(untested, up@universe)), character(0))
    expect_identical(length(up@universe), length(tested_annotated))
    expect_identical(sort(up@universe, method = "radix"), sort(tested_annotated, method = "radix"))
    # KIT-203's family under the tested universe: the 10-member term holds 4 tested proteins and leaves the family.
    expect_identical(unique(up@result$family_size), k205FamilySize(tested_annotated))
    expect_lt(k205FamilySize(tested_annotated), k205FamilySize(table$Protein.IDs))
})

test_that("KIT-205 B2: the g:Profiler background is each comparison's tested ids, for up and down; a shared name stays", {
    ids <- k205Ids()
    first <- k205Untested()
    second <- ids$fill[c(2L, 16L, 30L, 44L, 58L, 72L, 7L, 21L, 35L, 49L, 63L, 77L)]
    tables <- list(c1 = k205Table(first), c2 = k205Table(second, seed = 206L))
    calls <- k205Gprofiler(tables)
    # (a), (b): per contrast and direction, the stripped unique ids of that contrast's tested rows.
    for (contrast in names(tables)) {
        expected <- unique(stripEnrichmentIdSuffix(k205Tested(tables[[contrast]])$Protein.IDs))
        for (direction in c("up", "down")) {
            key <- paste0(contrast, "/", direction)
            expect_false(is.null(calls[[key]]), info = key)
            expect_identical(calls[[key]]$custom_bg, expected, info = key)
        }
    }
    expect_false(any(first %in% calls[["c1/up"]]$custom_bg))
    expect_true(all(setdiff(second, first) %in% calls[["c1/up"]]$custom_bg))
    expect_false(any(second %in% calls[["c2/down"]]$custom_bg))
    # (c): a gene name shared by one tested and one untested row stays in the background (rows are filtered, not ids).
    named <- k205Table(first)
    named$gene_name <- paste0("g", seq_len(nrow(named)))
    shared_tested <- match(ids$fill[[1L]], named$Protein.IDs)
    shared_untested <- match(first[[1L]], named$Protein.IDs)
    named$gene_name[c(shared_tested, shared_untested)] <- "sharedA"
    by_name <- k205Gprofiler(list(c1 = named), id_column = "gene_name")
    expect_true("sharedA" %in% by_name[["c1/up"]]$custom_bg)
    expect_false(named$gene_name[[match(first[[2L]], named$Protein.IDs)]] %in% by_name[["c1/up"]]$custom_bg)
})

# --- B3, B4, B5: guards and the edge rows -------------------------------------------------------------------------------

test_that("KIT-205 B3: a table with no untested row gives the same universe, background, TSV and results under both rules", {
    k205SkipUnlessKernel()
    table <- k205Table()
    current <- k205Cluster(list(k205 = table))
    earlier <- withr::with_options(list(MultiScholaR.enrichment_background = "all_rows"), k205Cluster(list(k205 = table)))
    for (direction in c("up", "down")) {
        a <- current$results@enrichment_data$k205[[direction]]
        b <- earlier$results@enrichment_data$k205[[direction]]
        expect_false(is.null(a), info = direction)
        expect_identical(a@universe, b@universe, info = direction)
        expect_identical(a@result, b@result, info = direction)
        expect_identical(attr(a, "enrichment_plot_inputs", exact = TRUE), attr(b, "enrichment_plot_inputs", exact = TRUE),
            info = direction)
        tsv <- paste0("k205_", direction, "_enrichment_results.tsv")
        expect_identical(digest::digest(file = file.path(current$dir, tsv), algo = "sha256"),
            digest::digest(file = file.path(earlier$dir, tsv), algo = "sha256"), info = direction)
    }
    expect_identical(length(current$results@enrichment_data$k205$up@universe), 120L)
    gp_current <- k205Gprofiler(list(k205 = table))
    gp_earlier <- withr::with_options(list(MultiScholaR.enrichment_background = "all_rows"),
        k205Gprofiler(list(k205 = table)))
    expect_identical(gp_current, gp_earlier)
    expect_identical(gp_current[["k205/up"]]$custom_bg, table$Protein.IDs)
})

test_that("KIT-205 B4: the lists lie inside the background; a row with a q-value but no p-value is in neither", {
    k205SkipUnlessKernel()
    table <- k205Table(k205Untested())
    calls <- k205Gprofiler(list(k205 = table))
    for (direction in c("up", "down")) {
        call <- calls[[paste0("k205/", direction)]]
        expect_true(length(call$ids) > 0L && all(call$ids %in% call$custom_bg), info = direction)
    }
    cluster <- k205Cluster(list(k205 = table))$results@enrichment_data$k205
    for (direction in c("up", "down")) {
        expect_true(all(intersect(cluster[[direction]]@gene, k205Ids()$pairs$gene) %in% cluster[[direction]]@universe),
            info = direction)
    }
    # Second part: one row with a finite q-value (in the up list's range) and no p-value.
    odd <- k205Table()
    row <- match(k205Ids()$fill[[10L]], odd$Protein.IDs)
    odd$raw_pvalue[[row]] <- NA_real_
    odd$fdr_qvalue[[row]] <- 0.001
    odd$log2FC[[row]] <- 1
    odd_id <- odd$Protein.IDs[[row]]
    odd_calls <- k205Gprofiler(list(k205 = odd))
    expect_false(odd_id %in% odd_calls[["k205/up"]]$ids)
    expect_false(odd_id %in% odd_calls[["k205/up"]]$custom_bg)
    odd_up <- k205Cluster(list(k205 = odd))$results@enrichment_data$k205$up
    expect_false(odd_id %in% odd_up@gene)
    expect_false(odd_id %in% odd_up@universe)
})

test_that("KIT-205 B5: a table without raw_pvalue (KIT-189's form) is reduced by its q-value alone", {
    k205SkipUnlessKernel()
    table <- k205Table(raw_p = FALSE)
    expect_false("raw_pvalue" %in% names(table))
    complete <- k205Cluster(list(k205 = table))$results@enrichment_data$k205$up
    expect_identical(sort(complete@universe, method = "radix"), sort(table$Protein.IDs, method = "radix"))
    # One row without a q-value (its fold change finite): left out of the universe.
    row <- match(k205Ids()$fill[[12L]], table$Protein.IDs)
    table$fdr_qvalue[[row]] <- NA_real_
    table$log2FC[[row]] <- 0.5
    reduced <- k205Cluster(list(k205 = table))$results@enrichment_data$k205$up
    expect_false(table$Protein.IDs[[row]] %in% reduced@universe)
    expect_identical(length(reduced@universe), 119L)
})

# --- B6: replay -------------------------------------------------------------------------------------------------------

test_that("KIT-205 B6 (a): under the earlier rule the background is every row", {
    table <- k205Table(k205Untested())
    calls <- withr::with_options(list(MultiScholaR.enrichment_background = "all_rows"),
        k205Gprofiler(list(k205 = table)))
    expect_identical(calls[["k205/up"]]$custom_bg, table$Protein.IDs)
    expect_identical(calls[["k205/down"]]$custom_bg, table$Protein.IDs)
})

#' The proteomics `enrichment` group with its function fields (as test-workflow-analysis-record.R's
#' `recordEnrichmentConfig`, KIT-107 D1)
k205EnrichmentConfig <- function(config) {
    config[["organism"]] <- config[["organism"]] %||% list(supported = buildProtEnrichSupportedOrganisms,
        method = resolveProtEnrichAnalysisMethod)
    config$annotation <- config$annotation %||% list(resolve = resolveProtEnrichUniprotAnnotations,
        match = resolveProtEnrichAnnotationMatching)
    config$kernel <- config$kernel %||% list(select = resolveProtEnrichProcessEnrichmentsFn,
        test_double = runProtEnrichDeterministicProcessEnrichments)
    config$organism_filter$mapping <- config$organism_filter$mapping %||% resolveProtEnrichOrganismMapping
    config$organism_filter$apply <- config$organism_filter$apply %||% applyProtEnrichOrganismFilter
    config
}

test_that("KIT-205 B6 (b, c): the record replay runs the enrichment step under its own background rule; a fresh run records it", {
    withr::local_envvar(MULTISCHOLAR_TEST_MODE = "true")
    enrichment <- omicGuiProteomicsEnrichment
    testthat::local_mocked_bindings(omicGuiProteomicsEnrichment = function() k205EnrichmentConfig(enrichment()))
    # The enrichment record of test-workflow-analysis-record.R (KIT-107 F-16): a DIA store, the DA, then the enrichment
    # as the tab runs it on the deterministic double (taxon 9606).
    root <- withr::local_tempdir()
    manager <- newDIAStateAdapter(newWorkflowEntityStore(root, "enrichment-record", "DIA"), WorkflowState$new(audit_enabled = FALSE))
    withr::defer(closeDIAStateAdapter(manager))
    manager$.context$reference_states <- TRUE
    manager$setWorkflowType("DIA")
    design <- module_ci_prot_da_design()
    rownames(design) <- NULL
    quant <- module_ci_prot_da_table()
    rownames(quant) <- NULL
    setDIAStepProvenance(manager, "design_builder")
    manager$saveState("correlation_filtered", module_ci_prot_da_object(design = design, protein_quant_table = quant), list(), "fixture")
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
        manifest <- suppressMessages(suppressWarnings(runDIAReferenceDaAnalysis(list(formula_string = "~ 0 + group",
            da_q_val_thresh = 0.05, treat_lfc_cutoff = 0), da_data, workflow, paths, contrasts = contrasts, publish = FALSE)))
        suppressMessages(suppressWarnings(produceDIADaPublication(manager, manifest, data.frame(), 0.05,
            paths[c("da_output_dir", "publication_graphs_dir")], save_assets = FALSE)))
        utils::capture.output(suppressMessages(suppressWarnings(runWorkflowEnrichment(workflow,
            omicGuiConfig("proteomics")$enrichment, manifest, list(selected_contrast = "B_vs_A", organism_taxid = "9606",
                up_cutoff = 0, down_cutoff = 0, q_cutoff = 0.05, correction_method = "gSCS", enable_organism_filter = FALSE),
            paths))))
    })
    expect_identical(manager$current_state, "correlation_filtered") # KIT-216 (VAL-175): the head stays the table
    # (c) the fresh run's step provenance holds the marker.
    expect_identical(manager$.context$snapshot$step_provenance$enrichment_completed$fixed,
        list(enrichment_background = "tested"))
    json <- file.path(root, "analysis_record.json")
    writeWorkflowAnalysisRecord(manager, json)
    record <- jsonlite::read_json(json, simplifyVector = FALSE)
    step_at <- which(vapply(record$steps, function(s) identical(s$operation, "enrichment"), logical(1)))
    expect_length(step_at, 1L)
    # (b) the kernel the replay selects (the double, test mode) is wrapped by a recorder of the option.
    seen <- character()
    original <- runProtEnrichDeterministicProcessEnrichments
    testthat::local_mocked_bindings(runProtEnrichDeterministicProcessEnrichments = function(...) {
        seen <<- c(seen, getOption("MultiScholaR.enrichment_background", "<unset>"))
        original(...)
    })
    h <- new.env(parent = globalenv())
    sys.source(file.path(testthat::test_path("..", ".."), "tools", "profiling", "workflow_kit.R"), h)
    api <- function(name) get(name, asNamespace("MultiScholaR"))
    replayWith <- function(record_used) {
        copy <- withr::local_tempdir(.local_envir = parent.frame())
        stopifnot(system2("cp", c("-a", file.path(root, "."), copy)) == 0L)
        prepared <- h$kitPrepareReplayProject(copy, c("initial", "correlation_filtered"), api)
        seen <<- character()
        utils::capture.output(replay <- suppressMessages(suppressWarnings(replayWorkflowAnalysisRecord(record_used,
            prepared, character(), paths, withr::local_tempdir()))))
        list(states = replay$states, seen = seen)
    }
    as_made <- replayWith(record)
    expect_identical(as_made$states, c("initial", "correlation_filtered")) # KIT-216 (VAL-175)
    expect_gt(length(as_made$seen), 0L)
    expect_identical(unique(as_made$seen), "tested")
    expect_null(getOption("MultiScholaR.enrichment_background"))
    earlier <- record
    earlier$steps[[step_at]]$parameters <- Filter(function(p) !identical(p$name, "enrichment_background"),
        earlier$steps[[step_at]]$parameters)
    old <- replayWith(earlier)
    expect_identical(old$states, c("initial", "correlation_filtered")) # KIT-216 (VAL-175)
    expect_gt(length(old$seen), 0L)
    expect_identical(unique(old$seen), "all_rows")
    expect_null(getOption("MultiScholaR.enrichment_background"))
})

# --- B7: the reports ---------------------------------------------------------------------------------------------------

#' A template's `FormatGOEnrichment`, extracted with knitr::purl and evaluated (KIT-173 T11's way of rendering)
k205Formatter <- function(template) {
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

test_that("KIT-205 B7: the report reads the marker; each template states the background and family of its rule", {
    write_record <- function(parameters) {
        path <- tempfile(fileext = ".json")
        jsonlite::write_json(list(steps = list(list(operation = "enrichment",
            parameters = lapply(names(parameters), function(name) list(name = name, source = "fixed",
                value = encodeWorkflowRecordValue(parameters[[name]]))), outcomes = list()))), path,
            auto_unbox = TRUE, null = "null")
        path
    }
    read <- workflowRecordReportParameters(write_record(list(enrichment_background = "tested")), list())
    expect_identical(read$enrichment_background, "tested")
    absent <- workflowRecordReportParameters(write_record(list(q_cutoff = 0.05)), list())
    expect_null(absent$enrichment_background)
    expect_true(isTRUE(absent$record_present))
    family <- c(tested = "terms with 10 to 500 annotated tested proteins",
        absent = "terms with 10 to 500 annotated proteins of that comparison's differential-abundance table")
    tails <- c(bh_whole_family = ", whether or not a listed protein falls in them (terms without one enter with p = 1).",
        qvalue = " and at least one protein in the list.")
    background <- c(tested = paste("The background was the proteins tested in that comparison (those with a p-value)",
        "with at least one Gene Ontology annotation."),
        absent = paste("The background was every protein listed in that comparison's differential-abundance table,",
            "including any without a p-value, with at least one Gene Ontology annotation."))
    clause <- paste("the first gene name of each protein, from the protein annotation or, where the annotation gives",
        "none, from the identifier table; proteins without a gene name were left out of the lists and the background")
    gprofiler <- c(tested = paste0("• The background was the identifiers of the proteins tested in each comparison ",
        "(those with a p-value in its differential-abundance table; ", clause, ") that g:Profiler recognises (a custom ",
        "background)."),
        absent = paste0("• The background was the identifiers of each comparison's differential-abundance table (",
            clause, ") that g:Profiler recognises (a custom background)."))
    cases <- list(tested = list(expect = "tested", extra = list(enrichment_background = "tested", record_present = TRUE)),
        absent = list(expect = "absent", extra = list(record_present = TRUE)),
        no_record = list(expect = "tested", extra = list()))
    numbers <- function(source, extra) c(list(ui_enrichment_organism_selected = "999999", ui_enrichment_up_lfc = 0,
        ui_enrichment_down_lfc = 0, ui_enrichment_q_value_cutoff = 0.05, ui_enrichment_database_source = source), extra)
    for (template in c("LFQ_report.rmd", "TMT_report.rmd", "DIANN_report.rmd", "DIANN_limpa_report.rmd")) {
        format <- k205Formatter(template)
        for (case in names(cases)) {
            want <- cases[[case]]$expect
            other <- setdiff(names(family), want)
            for (rule in names(tails)) {
                info <- paste(template, case, rule)
                text <- as.character(format(numbers("clusterprofiler", c(cases[[case]]$extra,
                    list(enrichment_significance_rule = rule)))))
                expect_true(grepl(paste0(family[[want]], tails[[rule]]), text, fixed = TRUE), info = info)
                expect_true(grepl(background[[want]], text, fixed = TRUE), info = info)
                expect_false(grepl(family[[other]], text, fixed = TRUE), info = info)
                expect_false(grepl(background[[other]], text, fixed = TRUE), info = info)
                expect_false(grepl("measured proteins", text, fixed = TRUE), info = info)
            }
            text <- as.character(format(numbers("gprofiler2", cases[[case]]$extra)))
            expect_true(grepl(gprofiler[[want]], text, fixed = TRUE), info = paste(template, case, "gprofiler2"))
            expect_false(grepl(gprofiler[[other]], text, fixed = TRUE), info = paste(template, case, "gprofiler2"))
            expect_false(grepl("otherwise protein accessions", text, fixed = TRUE), info = paste(template, case))
        }
    }
})

# --- B8: the organism filter --------------------------------------------------------------------------------------------

test_that("KIT-205 B8: the organism filter commutes with the reduction to tested rows", {
    table <- k205Table(k205Untested())
    names(table)[names(table) == "Protein.IDs"] <- "uniprot_acc"
    mapping <- data.frame(uniprot_acc = table$uniprot_acc, taxon_id = ifelse(seq_len(nrow(table)) %% 3L == 0L,
        "562", "83333"), stringsAsFactors = FALSE)
    filter <- function(t) applyProtEnrichOrganismFilter(k205Da(list(k205 = t)), mapping, "83333",
        catFn = function(...) invisible(NULL))$daResultsForEnrichment@da_data$k205
    one <- k205Tested(filter(table))
    two <- filter(k205Tested(table))
    expect_identical(sort(one$uniprot_acc), sort(two$uniprot_acc))
    expect_true(nrow(one) > 0L && nrow(one) < nrow(table))
})

# --- B9: KIT-203's texts (SD-1 (b)) -------------------------------------------------------------------------------------

#' The output ids of a UI tree (KIT-203 S10's walk)
k205Outputs <- function(tag) {
    if (inherits(tag, "shiny.tag")) {
        class <- paste(unlist(tag$attribs[names(tag$attribs) == "class"]), collapse = " ")
        here <- if (grepl("shiny-[a-z-]*output|html-widget-output", class)) tag$attribs$id
        return(c(here, k205Outputs(tag$children)))
    }
    if (is.list(tag)) return(unlist(lapply(tag, k205Outputs)))
    NULL
}

test_that("KIT-205 B9: KIT-203's label, note, notice and legend say tested proteins; the module adds no output", {
    rules <- ENRICHMENT_SIGNIFICANCE_RULES
    default <- rules[rules$rule == "bh_whole_family", , drop = FALSE]
    expect_identical(default$label, "Benjamini-Hochberg over every term of 10 to 500 tested proteins (default)")
    expect_true(grepl("that has 10 to 500 of your tested proteins, including terms with none", default$note, fixed = TRUE))
    expect_identical(default$notice, paste("For this organism, terms are called significant by the Benjamini-Hochberg",
        "adjusted p-value over every term of 10 to 500 tested proteins."))
    for (column in c("label", "note", "notice", "methods")) {
        expect_false(any(grepl("measured proteins", rules[[column]], fixed = TRUE)), info = column)
    }
    ui <- mod_omic_enrich_ui("enrichment_analysis", omicGuiConfig("proteomics")$enrichment)
    html <- as.character(ui)
    expect_true(grepl("p_adjust_family: Benjamini-Hochberg over every term of 10 to 500 tested proteins.", html,
        fixed = TRUE))
    expect_false(grepl("measured proteins", html, fixed = TRUE))
    expect_identical(sort(unique(k205Outputs(ui))), c(
        "enrichment_analysis-analysis_method_display", "enrichment_analysis-clusterprofiler_plot",
        "enrichment_analysis-clusterprofiler_results_table", "enrichment_analysis-clusterprofiler_summary_stats",
        "enrichment_analysis-contrasts_display", "enrichment_analysis-enrichment_status",
        "enrichment_analysis-gprofiler_plot", "enrichment_analysis-gprofiler_results_table",
        "enrichment_analysis-gprofiler_summary_stats", "enrichment_analysis-stringdb_plot",
        "enrichment_analysis-stringdb_results_table", "enrichment_analysis-stringdb_summary_stats"))
})
