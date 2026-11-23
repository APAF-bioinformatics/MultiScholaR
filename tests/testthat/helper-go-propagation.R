# KIT-173: fixtures and builders shared by test-go-propagation.R and the ticket's evidence scripts under
# docs/qa/workflow-execution/kit/kit-173/, which source this file into an environment whose parent is the MultiScholaR
# namespace. Everything here also runs on the base commit (T7(a)'s capture runs there), so it calls no function that
# KIT-173 adds.

#' The committed fixtures (README.md gives their sources, releases and SHA-256)
goPropFixtureDir <- function(...) testthat::test_path("..", "testdata", "go-propagation", ...)

#' A committed UniProt TSV (plain or gzip), every column as text, empty fields as ""
goPropReadTsv <- function(name) {
    utils::read.delim(goPropFixtureDir(name), sep = "\t", quote = "", stringsAsFactors = FALSE,
        colClasses = "character", na.strings = character())
}

#' The route's conversion of a raw UniProt table (`uniprotGoIdToTerm`, as `getUniprotAnnotations` calls it)
goPropConvert <- function(raw) {
    uniprotGoIdToTerm(raw, uniprot_id_column = Entry, go_id_column = Gene.Ontology.IDs, sep = "; ")
}

#' A table's content hash: column names and every value as text, in row order (class attributes ignored)
goPropTableHash <- function(x) {
    if (is.null(x)) return(NA_character_)
    x <- as.data.frame(x, stringsAsFactors = FALSE)
    digest::digest(list(names(x), lapply(x, function(column) as.character(unname(column)))), algo = "sha256")
}

#' Replace one namespace binding for the caller's lifetime; returns a function restoring it
goPropReplaceBinding <- function(package, name, value) {
    ns <- asNamespace(package)
    original <- get(name, envir = ns, inherits = FALSE)
    locked <- bindingIsLocked(name, ns)
    if (locked) unlockBinding(name, ns)
    assign(name, value, envir = ns)
    if (locked) lockBinding(name, ns)
    function() {
        if (locked) unlockBinding(name, ns)
        assign(name, original, envir = ns)
        if (locked) lockBinding(name, ns)
    }
}

#' A recording wrapper around the real `clusterProfiler::enricher`: each call's TERM2GENE and TERM2NAME, and the
#' hashes of both, are appended to `acc` (an environment). Returns the restore function.
goPropRecordEnricher <- function(acc) {
    real <- get("enricher", envir = asNamespace("clusterProfiler"))
    acc$calls <- list()
    goPropReplaceBinding("clusterProfiler", "enricher", function(gene, universe = NULL, TERM2GENE, TERM2NAME = NA, ...) {
        acc$calls[[length(acc$calls) + 1L]] <- list(term2gene = TERM2GENE, term2name = TERM2NAME,
            term2gene_sha256 = goPropTableHash(TERM2GENE), term2name_sha256 = goPropTableHash(TERM2NAME))
        real(gene = gene, universe = universe, TERM2GENE = TERM2GENE, TERM2NAME = TERM2NAME, ...)
    })
}

#' The T7 project's annotation table: the committed E. coli fixture (F1's selection), with three GO ids GO.db 3.23.1
#' does not hold appended to the first `fill` protein (T3(b): GO:0000002 obsolete, GO:0000021 secondary to GO:0007131,
#' GO:9999999 unknown), through `getUniprotAnnotations` with the download stubbed (no service call).
goPropT7Raw <- function() {
    raw <- goPropReadTsv("t7-annotation-ecoli.tsv")
    roles <- goPropReadTsv("t7-fixture-entries.tsv")
    target <- sort(roles$entry[roles$role == "fill"])[[1L]]
    at <- match(target, raw$Entry)
    raw$Gene.Ontology.IDs[[at]] <- paste(c(raw$Gene.Ontology.IDs[[at]], "GO:0000002", "GO:0000021", "GO:9999999"),
        collapse = "; ")
    raw
}

goPropAnnotate <- function(raw, cache_dir) {
    restore <- goPropReplaceBinding("MultiScholaR", "directUniprotDownload", function(input_tbl, output_path, ...) raw)
    on.exit(restore(), add = TRUE)
    suppressMessages(getUniprotAnnotations(data.frame(Protein.Ids = raw$Entry, stringsAsFactors = FALSE),
        cache_dir = cache_dir, taxon_id = 83333))
}

#' The project's quantities: three samples per group; the `up` proteins 3 log2 units higher in B, the `down` proteins
#' 3 lower, every other protein the same in both groups (equal group means, the same within-group spread).
goPropQuant <- function(entries, up, down) {
    shift <- (seq_along(entries) %% 5) / 10
    base <- 20 + (seq_along(entries) %% 7) + shift
    effect <- 3 * (entries %in% up) - 3 * (entries %in% down)
    values <- cbind(base - 0.1, base, base + 0.1, base + effect + 0.1, base + effect - 0.1, base + effect)
    colnames(values) <- paste0("S", 1:6)
    data.frame(uniprot_acc = entries, values, check.names = FALSE, stringsAsFactors = FALSE)
}

goPropPaths <- function(root) {
    list(source_dir = file.path(root, "source"), results_dir = file.path(root, "results"),
        da_output_dir = file.path(root, "results", "da_proteins"),
        pathway_dir = file.path(root, "results", "pathway_enrichment"),
        publication_graphs_dir = file.path(root, "publication_graphs"))
}

#' The enrichment product a store's head roots (schema `dia_enrichment_product_v1`), hydrated
goPropProduct <- function(manager) {
    for (reference in manager$.context$snapshot$products) {
        value <- tryCatch(hydrateSavedWorkflowEntity(manager$.context$store, reference)$value, error = function(e) NULL)
        if (is.list(value) && identical(value$schema, "dia_enrichment_product_v1")) return(value)
    }
    NULL
}

#' What T7 compares of one enrichment: the stored result tables, the plots' data and the pathway files
goPropRunEvidence <- function(state, pathway_dir) {
    tables <- lapply(state@enrichment_data, function(contrast) lapply(contrast, function(result)
        if (is.null(result)) NULL else result@result))
    plot_inputs <- lapply(state@enrichment_data, function(contrast) lapply(contrast, function(result)
        if (is.null(result)) NA_character_ else digest::digest(attr(result, "enrichment_plot_inputs"), algo = "sha256")))
    files <- sort(list.files(pathway_dir))
    file_hashes <- data.frame(file = files, sha256 = vapply(files, function(f) if (grepl("\\.tsv$", f))
        digest::digest(file = file.path(pathway_dir, f), algo = "sha256") else "<exists>", ""), stringsAsFactors = FALSE)
    rownames(file_hashes) <- NULL
    list(tables = tables, plot_inputs = plot_inputs, files = file_hashes)
}

#' One proteomics project on the clusterProfiler route (taxon 999999) with the real kernel: the T7 fixture's
#' proteins, the DA fitted, published with the annotation table, then the enrichment as the tab runs it, with the
#' annotation table the project's `uniprot_dat_cln.RDS` (as the replay reads it). Writes the analysis record.
#' @param go_propagation NULL (the tab's case: the run takes its default) or a value passed in the run's inputs.
goPropRunProject <- function(root, annotations, go_propagation = NULL, process_fn = processEnrichments) {
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
    parameters <- c(list(selected_contrast = "B_vs_A", organism_taxid = "999999", up_cutoff = 0, down_cutoff = 0,
        q_cutoff = 0.05, correction_method = "gSCS", enable_organism_filter = FALSE),
        if (!is.null(go_propagation)) list(go_propagation = go_propagation))
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
            process_fn = process_fn)
    })
    record <- file.path(root, "analysis_record.json")
    writeWorkflowAnalysisRecord(manager, record)
    state <- withDIAStateScope(manager, function() manager$getState("enrichment_completed"))
    list(run = run, record = record, paths = paths, state = state, product = goPropProduct(manager),
        provenance = manager$.context$snapshot$step_provenance$enrichment_completed,
        outcomes = readComponentMetadata(manager, manager$.context$snapshot$states$enrichment_completed)$config$outcomes)
}

#' Replay a committed or fresh project's record on a prepared copy (cut to its import and design states), the
#' enrichment with the kernel the GUI selects outside test mode (`processEnrichments`).
goPropReplay <- function(project_root, record, scratch) {
    copy <- file.path(scratch, "project")
    dir.create(copy, recursive = TRUE)
    stopifnot(system2("cp", c("-a", file.path(project_root, "."), copy)) == 0L)
    h <- new.env(parent = globalenv())
    sys.source(file.path(testthat::test_path("..", ".."), "tools", "profiling", "workflow_kit.R"), h)
    api <- function(name) get(name, asNamespace("MultiScholaR"))
    prepared <- h$kitPrepareReplayProject(copy, c("initial", "correlation_filtered"), api)
    replayed <- list()
    out <- file.path(scratch, "replay")
    dir.create(out)
    result <- replayWorkflowAnalysisRecord(record, prepared, character(), goPropPaths(copy), out,
        on_step = function(step, manager) {
            if (identical(step$operation, "enrichment")) {
                replayed$state <<- manager$getState("enrichment_completed")
                replayed$product <<- goPropProduct(manager)
                replayed$provenance <<- manager$.context$snapshot$step_provenance$enrichment_completed
            }
        })
    replayed$pathway_dir <- file.path(out, "enrichment", "pathway")
    replayed$states <- result$states
    replayed
}

# --- T5 (the model-organism check, section 5.2) --------------------------------------------------------------------

#' UniProt entry to ordered locus name(s)
goPropOln <- function(raw) {
    names <- strsplit(raw$Gene.Names..ordered.locus., "[ ;,/]+")
    map <- data.frame(Entry = rep(raw$Entry, lengths(names)), GENE = unlist(names), stringsAsFactors = FALSE)
    unique(map[nzchar(map$GENE), , drop = FALSE])
}

#' The committed g:Profiler GO GMT as (TERM, GENE)
goPropReadGmt <- function() {
    lines <- readLines(goPropFixtureDir("gprofiler-scerevisiae-go.ENSG.gmt.gz"), warn = FALSE)
    parts <- strsplit(lines, "\t", fixed = TRUE)
    data.frame(TERM = rep(vapply(parts, `[[`, "", 1L), lengths(parts) - 2L),
        GENE = unlist(lapply(parts, function(p) p[-(1:2)])), stringsAsFactors = FALSE)
}

#' TERM2GENE as `processEnrichments` passes it to `enricher` for the yeast fixture (one contrast; the 50 first entries
#' called up). `go_propagation` NULL calls the kernel without the argument (its default).
goPropYeastTerm2Gene <- function(annotations, go_propagation = NULL, pathway_dir) {
    entries <- sort(annotations$Entry)
    table <- data.frame(Protein.IDs = entries, fdr_qvalue = ifelse(seq_along(entries) <= 50L, 0.01, 0.5),
        log2FC = ifelse(seq_along(entries) <= 50L, 2, 0), stringsAsFactors = FALSE)
    da <- methods::new("da_results_for_enrichment", contrasts = tibble::tibble(contrast = "yeast"),
        da_data = list(yeast = table), design_matrix = data.frame(sample = "S1", stringsAsFactors = FALSE))
    acc <- new.env()
    restore <- goPropRecordEnricher(acc)
    on.exit(restore(), add = TRUE)
    args <- list(da_results = da, taxon_id = "999999", q_cutoff = 0.05, pathway_dir = pathway_dir,
        go_annotations = annotations, exclude_iea = FALSE, protein_id_column = "Protein.IDs", contrast_names = "yeast")
    if (!is.null(go_propagation)) args$go_propagation <- go_propagation
    utils::capture.output(suppressMessages(suppressWarnings(do.call(processEnrichments, args))))
    acc$calls[[1L]]$term2gene
}

#' Per-term agreement of the route's sets with g:Profiler's (containment C and precision P over genes in both sources)
#' @param sets Named list of TERM2GENE tables (columns TERM, then the UniProt entry), one per value; the scope is set by
#'   `scope_value`: terms in that value's sets and in the GMT, 10 to 500 g:Profiler genes, held by GO.db.
goPropAgreement <- function(sets, raw, gmt, scope_value) {
    oln <- goPropOln(raw)
    universe <- intersect(unique(gmt$GENE), oln$GENE)
    gp <- unique(gmt[gmt$GENE %in% universe, , drop = FALSE])
    gp_size <- table(gp$TERM)
    aspects <- AnnotationDbi::Ontology(GO.db::GOTERM)
    asGenes <- function(t2g) {
        t2g <- data.frame(TERM = t2g[[1L]], Entry = t2g[[2L]], stringsAsFactors = FALSE)
        mapped <- merge(t2g, oln, by = "Entry")
        unique(mapped[mapped$GENE %in% universe, c("TERM", "GENE")])
    }
    genes <- lapply(sets, asGenes)
    scope <- names(gp_size)[gp_size >= 10L & gp_size <= 500L]
    scope <- sort(intersect(intersect(scope, names(aspects)), unique(genes[[scope_value]]$TERM)))
    ref <- split(gp$GENE, gp$TERM)[scope]
    do.call(rbind, lapply(names(genes), function(value) {
        ours <- split(genes[[value]]$GENE, genes[[value]]$TERM)
        do.call(rbind, lapply(scope, function(term) {
            o <- ours[[term]] %||% character()
            r <- ref[[term]]
            both <- length(intersect(o, r))
            data.frame(value = value, term = term, aspect = unname(aspects[[term]]), ours = length(o),
                reference = length(r), both = both, C = both / length(r), P = if (length(o)) both / length(o) else NA_real_,
                stringsAsFactors = FALSE)
        }))
    }))
}
