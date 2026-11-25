# KIT-208 (VAL-107 as decided by VAL-147): the DA publication's gene-name join keeps one `gene_names` column (the
# annotation's name, the identifier table's filling blanks, SD-2), and the g:Profiler route never sends an empty list:
# a comparison whose background has no named protein is refused (class `workflow_enrichment_no_identifiers`), a list
# whose proteins are all unnamed is reported as not tested (SD-1 (a)), and the unnamed proteins per comparison are
# recorded (SD-3 (c)). Tests T1 to T9 of the ticket (section 10). No service is called: `gprofiler2::gost` is replaced
# by a recorder that returns NULL (as KIT-175's F1, docs/qa/workflow-execution/kit/kit-175/f1-consumers.R), the proxy
# points at a closed port, and the clusterProfiler route runs the installed `enricher` on KIT-173's fixture.

library(testthat)

k208Offline <- function(env = parent.frame()) {
    withr::local_envvar(http_proxy = "http://127.0.0.1:9", https_proxy = "http://127.0.0.1:9",
        HTTP_PROXY = "http://127.0.0.1:9", HTTPS_PROXY = "http://127.0.0.1:9", no_proxy = "", NO_PROXY = "",
        MULTISCHOLAR_TEST_MODE = "false", .local_envir = env)
    withr::local_options(multischolar.test_mode = FALSE, .local_envir = env)
}

#' Replace `gprofiler2::gost` by a recorder returning NULL for the caller's lifetime; each call is logged as
#' "<query size>,<background size>" (distinct identifiers, as F1 counted them)
k208RecordGost <- function(env = parent.frame()) {
    log <- new.env(parent = emptyenv())
    log$calls <- character()
    restore <- goPropReplaceBinding("gprofiler2", "gost", function(query, ..., custom_bg = NULL) {
        log$calls <- c(log$calls, sprintf("%d,%d", length(unique(unlist(query, use.names = FALSE))),
            length(unique(custom_bg))))
        NULL
    })
    withr::defer(restore(), envir = env)
    log
}

#' Wrap the installed `clusterProfiler::enricher` for the caller's lifetime; each call's list and universe are logged
k208RecordEnricher <- function(env = parent.frame()) {
    log <- new.env(parent = emptyenv())
    log$calls <- list()
    real <- get("enricher", envir = asNamespace("clusterProfiler"))
    restore <- goPropReplaceBinding("clusterProfiler", "enricher", function(gene, universe = NULL, ...) {
        log$calls[[length(log$calls) + 1L]] <- list(gene = sort(unique(gene)), universe = sort(unique(universe)))
        real(gene = gene, universe = universe, ...)
    })
    withr::defer(restore(), envir = env)
    log
}

#' Run `expr` with console output, messages and warnings muffled; the error, if any, is returned
k208Capture <- function(expr) {
    error <- NULL
    utils::capture.output(tryCatch(withCallingHandlers(expr,
        warning = function(w) invokeRestart("muffleWarning"),
        message = function(m) invokeRestart("muffleMessage")),
        error = function(e) {
            error <<- e
            NULL
        }))
    error
}

#' KIT-173's E. coli fixture annotated through `getUniprotAnnotations` (no service; 70 rows)
k208Annotations <- function() goPropAnnotate(goPropReadTsv("t7-annotation-ecoli.tsv"), withr::local_tempdir())

k208Roles <- function() goPropReadTsv("t7-fixture-entries.tsv")

#' The first word of each annotation's `gene_names`, as the DA writer takes it (split on space and colon)
k208FirstName <- function(annotations, entries) {
    first <- vapply(strsplit(annotations$gene_names, " |:"), function(x) if (length(x)) x[[1L]] else "", "")
    first[match(entries, annotations$Entry)]
}

#' A copy of KIT-175's `k175Project` (test-prot-annotation-uniprot-empty.R) on the g:Profiler route (taxon 9606) by
#' default, with `id_gene_names`: TRUE gives the identifier table `gene_names = GENE<i>` (KIT-173's fixture), FALSE
#' the identifier table of an application import (`uniprot_acc` only). `annotations` is written to
#' `uniprot_dat_cln.RDS`, `bus` put on `workflow$uniprot_dat_cln`, `da_annotations` publishes the DA (default: the
#' annotations, or an empty table). `gost` is always replaced by the recorder, so no service is called.
k208Project <- function(root, annotations = NULL, id_gene_names = TRUE, bus = NULL, da_annotations = NULL,
    taxon = "9606", record_enricher = FALSE, q_cutoff = 0.05) {
    gost <- k208RecordGost()
    enricher <- if (record_enricher) k208RecordEnricher()
    roles <- k208Roles()
    paths <- goPropPaths(root)
    for (p in paths) dir.create(p, recursive = TRUE, showWarnings = FALSE)
    if (!is.null(annotations)) saveRDS(annotations, file.path(paths$source_dir, "uniprot_dat_cln.RDS"))
    manager <- newDIAStateAdapter(newWorkflowEntityStore(root, "k208-project", "DIA"),
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
    if (!is.null(bus)) workflow$uniprot_dat_cln <- bus
    parameters <- list(selected_contrast = "B_vs_A", organism_taxid = taxon, up_cutoff = 0, down_cutoff = 0,
        q_cutoff = q_cutoff, correction_method = "gSCS", enable_organism_filter = FALSE)
    published <- da_annotations %||% (if (is.data.frame(annotations)) annotations else createEmptyUniprotTable())
    notes <- list()
    kernel_calls <- 0L
    kernel <- function(...) {
        kernel_calls <<- kernel_calls + 1L
        processEnrichments(...)
    }
    run <- NULL
    error <- k208Capture(withDIAStateScope(manager, function() {
        da_data <- new.env(parent = emptyenv())
        da_data$current_s4_object <- getDIAStateReference(manager)
        da_data$contrasts_available <- contrasts$contrasts
        manifest <- runDIAReferenceDaAnalysis(list(formula_string = "~ 0 + group", da_q_val_thresh = 0.05,
            treat_lfc_cutoff = 0), da_data, workflow, paths, contrasts = contrasts, publish = FALSE)
        produceDIADaPublication(manager, manifest, published, 0.05,
            paths[c("da_output_dir", "publication_graphs_dir")], save_assets = FALSE)
        run <<- runWorkflowEnrichment(workflow, omicGuiConfig("proteomics")$enrichment, manifest, parameters, paths,
            process_fn = kernel, notify = function(message, type, duration) {
                notes[[length(notes) + 1L]] <<- list(text = message, type = type, duration = duration)
            })
    }))
    da_file <- list.files(paths$da_output_dir, "_long_annot\\.tsv$", full.names = TRUE)
    state_ref <- manager$.context$snapshot$states$enrichment_completed
    list(run = run, error = error, notes = notes, kernel_calls = kernel_calls, gost = gost$calls,
        enricher = if (record_enricher) enricher$calls,
        da = if (length(da_file)) utils::read.delim(da_file[[1L]], stringsAsFactors = FALSE),
        outcomes = if (!is.null(state_ref)) readComponentMetadata(manager, state_ref)$config$outcomes)
}

k208Blank <- function(x) sum(is.na(x) | !nzchar(as.character(x)))

k208Refusal <- paste("g:Profiler tests gene names, and none of the tested proteins of comparison B_vs_A has one (0 of",
    "70 tested proteins named), so nothing would be tested. Check that the protein annotation was retrieved for this organism",
    "(save the design again). This happens when the annotation was not retrieved for these proteins, or when their",
    "identifiers are not UniProt accessions.")

k208NotTestedUp <- paste("Comparison B_vs_A, up: none of its 25 proteins has a gene name, so this list was not sent to",
    "g:Profiler and was not tested (its empty result does not mean that no term is enriched).")

# --- (a) the join --------------------------------------------------------------------------------------------------

test_that("KIT-208 T1: both tables carry gene names; the DA table keeps the annotation's", {
    k208Offline()
    full <- k208Annotations()
    project <- k208Project(withr::local_tempdir(), annotations = full)
    da <- project$da
    expect_identical(k208Blank(da$gene_name), 0L) # base: 70 (KIT-175 F1)
    expect_identical(sum(names(da) == "gene_names"), 1L)
    expect_false(any(grepl("^gene_names\\.(x|y|input)$", names(da))))
    expect_identical(as.character(da$gene_name), k208FirstName(full, da$uniprot_acc))
})

test_that("KIT-208 T2: an annotation lacking 25 accessions; the identifier table fills those 25 (SD-2)", {
    k208Offline()
    full <- k208Annotations()
    lacks25 <- full[-seq_len(25L), , drop = FALSE]
    project <- k208Project(withr::local_tempdir(), annotations = lacks25)
    da <- project$da
    expect_identical(k208Blank(da$gene_name), 0L) # base: 70 (KIT-175 F1)
    entries <- sort(k208Roles()$entry)
    from_annotation <- da$uniprot_acc %in% lacks25$Entry
    expect_identical(sum(!from_annotation), 25L)
    expect_identical(as.character(da$gene_name[from_annotation]),
        k208FirstName(lacks25, da$uniprot_acc[from_annotation]))
    expect_identical(as.character(da$gene_name[!from_annotation]),
        paste0("GENE", match(da$uniprot_acc[!from_annotation], entries)))
})

test_that("KIT-208 T3 (guard): an application-shaped table takes the annotation's first name, as before", {
    k208Offline()
    full <- k208Annotations()
    project <- k208Project(withr::local_tempdir(), annotations = full, id_gene_names = FALSE)
    da <- project$da
    expect_null(project$error)
    expect_identical(sum(names(da) == "gene_names"), 1L)
    expect_identical(as.character(da$gene_name), k208FirstName(full, da$uniprot_acc))
    expect_identical(k208Blank(da$gene_name), 0L)
    expect_null(project$outcomes$enrichment_lists_not_tested)
    expect_null(project$outcomes$enrichment_unnamed_proteins)
})

test_that("KIT-208 T4 (guard): without annotation at the DA, the identifier table's names, as before", {
    k208Offline()
    full <- k208Annotations()
    project <- k208Project(withr::local_tempdir(), bus = full, da_annotations = createEmptyUniprotTable())
    da <- project$da
    entries <- sort(k208Roles()$entry)
    expect_identical(as.character(da$gene_name), paste0("GENE", match(da$uniprot_acc, entries)))
    expect_null(project$outcomes$enrichment_lists_not_tested)
    expect_null(project$outcomes$enrichment_unnamed_proteins)
})

# --- (b) g:Profiler ------------------------------------------------------------------------------------------------

test_that("KIT-208 T5: T1's project through g:Profiler sends the names", {
    k208Offline()
    project <- k208Project(withr::local_tempdir(), annotations = k208Annotations())
    # base: 10 calls, each "0,0" (KIT-175 F1 gp-complete); the recorder's NULL is asked five times per list
    expect_identical(project$gost, c(rep("25,70", 5L), rep("10,70", 5L)))
    expect_null(project$error)
    expect_null(project$outcomes$enrichment_lists_not_tested)
    expect_null(project$outcomes$enrichment_unnamed_proteins)
})

test_that("KIT-208 T6: no DA protein has a gene name (case 5); the run is refused before g:Profiler is asked", {
    k208Offline()
    unmatched <- k208Annotations()
    unmatched$Entry <- sprintf("Z%05d", seq_len(nrow(unmatched)))
    project <- k208Project(withr::local_tempdir(), annotations = unmatched, id_gene_names = FALSE)
    expect_identical(project$gost, character()) # base: 10 calls, each "0,0", and the run completed
    expect_true(inherits(project$error, "workflow_enrichment_no_identifiers"))
    expect_identical(conditionMessage(project$error), k208Refusal)
    expect_identical(project$kernel_calls, 1L)
    expect_null(project$outcomes)
})

test_that("KIT-208 T7: the up list has no named protein (case 6); it is not tested, reported and recorded", {
    k208Offline()
    full <- k208Annotations()
    roles <- k208Roles()
    lacks_up <- full[!full$Entry %in% roles$entry[roles$role == "p_term"], , drop = FALSE]
    expect_identical(nrow(lacks_up), 45L)
    project <- k208Project(withr::local_tempdir(), annotations = lacks_up, id_gene_names = FALSE)
    expect_identical(project$gost, rep("10,45", 5L)) # base: five "0,45", then five "10,45"
    expect_null(project$error)
    texts <- vapply(project$notes, function(note) note$text, "")
    expect_identical(sum(texts == k208NotTestedUp), 1L)
    expect_identical(project$notes[[match(k208NotTestedUp, texts)]]$type, "warning")
    expect_identical(sum(grepl("was not tested", texts, fixed = TRUE)), 1L)
    expect_null(project$run$results@enrichment_data[["groupB-groupA"]]$up)
    expect_identical(project$outcomes$enrichment_lists_not_tested, "B_vs_A: up")
    expect_identical(project$outcomes$enrichment_unnamed_proteins, "B_vs_A: 25 of 70")
})

test_that("KIT-208 T8 (guard): the clusterProfiler route sends the same accession lists", {
    k208Offline()
    roles <- k208Roles()
    project <- k208Project(withr::local_tempdir(), annotations = k208Annotations(), taxon = "999999",
        record_enricher = TRUE)
    expect_null(project$error)
    expect_identical(project$gost, character())
    expect_length(project$enricher, 2L)
    expect_identical(project$enricher[[1L]]$gene, sort(roles$entry[roles$role == "p_term"]))
    expect_identical(project$enricher[[2L]]$gene, sort(roles$entry[roles$role == "d_term"]))
    expect_identical(project$enricher[[1L]]$universe, sort(roles$entry))
    expect_identical(project$enricher[[2L]]$universe, sort(roles$entry))
    expect_null(project$outcomes$enrichment_lists_not_tested)
    expect_null(project$outcomes$enrichment_unnamed_proteins)
})

test_that("KIT-208 T9: perform_enrichment called directly with no identifier does not call g:Profiler", {
    log <- k208RecordGost()
    subset <- data.frame(gene_name = rep(NA_character_, 3L), stringsAsFactors = FALSE)
    call <- function() perform_enrichment(data_subset = subset, species = "hsapiens", threshold = 0.05,
        sources = "GO:BP", domain_scope = "custom", custom_bg = character(), max_retries = 5, wait_time = 0,
        protein_id_column = "gene_name")
    value <- "unset"
    suppressMessages(suppressWarnings(value <- call()))
    expect_identical(log$calls, character()) # base: 5 calls, each "0,0"
    expect_null(value)
    caught <- NULL
    suppressMessages(suppressWarnings(withCallingHandlers(call(),
        multischolar_enrichment_no_identifiers = function(condition) caught <<- condition)))
    expect_true(inherits(caught, "multischolar_enrichment_no_identifiers"))
    expect_identical(caught$list_size, 3L)
    expect_identical(caught$query_ids, 0L)
})

# --- At the merge after KIT-205 (implementation audit K208-I08, I01, I07): the counts follow the tested rows -------

#' One comparison's DA table in a `da_results_for_enrichment`, as the run reads it (`groupB-groupA`, shown as B_vs_A)
k208DaResults <- function(table) {
    methods::new("da_results_for_enrichment",
        contrasts = tibble::tibble(contrasts = "groupB-groupA", friendly_names = "B_vs_A"),
        da_data = list(`groupB-groupA` = table), design_matrix = data.frame(Run = "S1", stringsAsFactors = FALSE))
}

test_that("KIT-208 T10: the refusal's n counts the tested proteins only", {
    # Four tested proteins without a name (two past the cutoffs, up), two untested proteins with a name.
    table <- data.frame(gene_name = c(NA, NA, NA, NA, "NAMED1", "NAMED2"),
        log2FC = c(2, 2, 0.1, -0.1, 3, 3), raw_pvalue = c(0.001, 0.001, 0.5, 0.5, NA, NA),
        fdr_qvalue = c(0.01, 0.01, 0.6, 0.6, NA, NA), stringsAsFactors = FALSE)
    log <- k208RecordGost()
    parameters <- list(da_results = k208DaResults(table), taxon_id = "9606", up_cutoff = 0, down_cutoff = 0,
        q_cutoff = 0.05, pathway_dir = withr::local_tempdir(), go_annotations = NULL, exclude_iea = FALSE,
        protein_id_column = "gene_name", contrast_names = "groupB-groupA", correction_method = "gSCS",
        process_fn = processEnrichments)
    error <- k208Capture(diaEnrichmentOperation()$kernel(NULL, parameters, NULL))
    expect_identical(conditionMessage(error), paste("g:Profiler tests gene names, and none of the tested proteins of",
        "comparison B_vs_A has one (0 of 4 tested proteins named), so nothing would be tested. Check that the protein",
        "annotation was retrieved for this organism (save the design again). This happens when the annotation was not",
        "retrieved for these proteins, or when their identifiers are not UniProt accessions."))
    expect_true(inherits(error, "workflow_enrichment_no_identifiers"))
    expect_identical(log$calls, character())
})

test_that("KIT-208 T11: an untested protein without a name is not counted as unnamed", {
    # Three tested proteins with a name, one tested without, one untested without.
    table <- data.frame(gene_name = c("A1", "B1", "C1", NA, NA), log2FC = c(2, -2, 0, 0, 1),
        raw_pvalue = c(0.001, 0.001, 0.5, 0.5, NA), fdr_qvalue = c(0.01, 0.01, 0.6, 0.6, NA), stringsAsFactors = FALSE)
    da <- k208DaResults(table)
    expect_identical(workflowEnrichmentUnnamedProteins(da, da@contrasts, "gene_name"), "B_vs_A: 1 of 4")
})

test_that("KIT-208 T12 (guard, I01): no protein past the cutoffs and none named: the run completes, nothing sent", {
    k208Offline()
    unmatched <- k208Annotations()
    unmatched$Entry <- sprintf("Z%05d", seq_len(nrow(unmatched)))
    project <- k208Project(withr::local_tempdir(), annotations = unmatched, id_gene_names = FALSE, q_cutoff = 0)
    expect_null(project$error)
    expect_identical(project$gost, character())
    expect_null(project$outcomes$enrichment_lists_not_tested)
    expect_identical(project$outcomes$enrichment_unnamed_proteins, "B_vs_A: 70 of 70")
})

test_that("KIT-208 T13 (guard, I07): an annotation entry with an empty name takes the identifier table's", {
    k208Offline()
    full <- k208Annotations()
    entries <- sort(k208Roles()$entry)
    emptied <- full$Entry[1:5]
    full$gene_names[1:5] <- ""
    project <- k208Project(withr::local_tempdir(), annotations = full)
    da <- project$da
    at <- da$uniprot_acc %in% emptied
    expect_identical(sum(at), 5L)
    expect_identical(as.character(da$gene_name[at]), paste0("GENE", match(da$uniprot_acc[at], entries)))
    expect_identical(as.character(da$gene_name[!at]), k208FirstName(full, da$uniprot_acc[!at]))
})
