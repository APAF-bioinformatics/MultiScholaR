# KIT-175 F1 (R0): what each consumer of the UniProt annotation does without it, on the code under test.
#
# Usage, from the package root on a run-tests runner after the commit is installed (R1: the tests-only commit, whose R/
# equals the base):
#   Rscript docs/qa/workflow-execution/kit/kit-175/f1-consumers.R <out dir>
#
# One proteomics project (KIT-173's committed E. coli fixture, 70 proteins; tests/testthat/helper-go-propagation.R),
# DA fitted and published, then the enrichment as the tab runs it, test mode off, no service call: the clusterProfiler
# route (taxon 999999) with the real `clusterProfiler::enricher` wrapped to record its arguments, and the g:Profiler
# route (taxon 9606) with `gprofiler2::gost` replaced by a recorder that returns no result. Three annotation tables:
# `complete` (the fixture's 70 rows), `empty` (zero rows, the table an empty retrieval leaves) and `lacks25` (the first
# 25 entries removed: one chunk of 25 not looked up). Also the organism filter (on, complete and empty) and the
# integration input's feature table. Writes `f1-consumers.tsv` (counts and messages only; no identifier) and
# `f1-session.txt`.

args <- commandArgs(trailingOnly = TRUE)
out <- if (length(args) >= 1L) args[[1L]] else "k175-f1-out"
dir.create(out, showWarnings = FALSE, recursive = TRUE)
out <- normalizePath(out)
suppressPackageStartupMessages(library(MultiScholaR))
Sys.setenv(MULTISCHOLAR_TEST_MODE = "false", http_proxy = "http://127.0.0.1:9", https_proxy = "http://127.0.0.1:9",
    HTTP_PROXY = "http://127.0.0.1:9", HTTPS_PROXY = "http://127.0.0.1:9")
options(multischolar.test_mode = FALSE)
ns <- asNamespace("MultiScholaR")
h <- new.env(parent = ns)
sys.source(file.path("tests", "testthat", "helper-go-propagation.R"), envir = h)
api <- function(name) get(name, envir = ns)

work <- file.path(tempdir(), "k175-f1")
dir.create(work, recursive = TRUE)
complete <- h$goPropAnnotate(h$goPropReadTsv("t7-annotation-ecoli.tsv"), file.path(work, "cache"))
tables <- list(complete = complete, empty = api("createEmptyUniprotTable")(),
    lacks25 = complete[-seq_len(25L), , drop = FALSE])

calls <- new.env()
calls$enricher <- list()
calls$gost <- list()
real_enricher <- get("enricher", envir = asNamespace("clusterProfiler"))
restore_enricher <- h$goPropReplaceBinding("clusterProfiler", "enricher", function(gene, universe = NULL, ...) {
    calls$enricher[[length(calls$enricher) + 1L]] <- c(gene = length(unique(gene)), universe = length(unique(universe)))
    real_enricher(gene = gene, universe = universe, ...)
})
restore_gost <- h$goPropReplaceBinding("gprofiler2", "gost", function(query, ..., custom_bg = NULL) {
    sent <- unique(unlist(query, use.names = FALSE))
    calls$gost[[length(calls$gost) + 1L]] <- c(gene = length(sent), universe = length(unique(custom_bg)))
    NULL
})

run_case <- function(name, annotation, taxon, organism_filter = FALSE) {
    calls$enricher <- list()
    calls$gost <- list()
    root <- file.path(work, name)
    roles <- h$goPropReadTsv("t7-fixture-entries.tsv")
    paths <- h$goPropPaths(root)
    for (p in paths) dir.create(p, recursive = TRUE, showWarnings = FALSE)
    if (nrow(annotation) > 0L) saveRDS(annotation, file.path(paths$source_dir, "uniprot_dat_cln.RDS"))
    manager <- api("newDIAStateAdapter")(api("newWorkflowEntityStore")(root, "k175-f1", "DIA"),
        api("WorkflowState")$new(audit_enabled = FALSE))
    on.exit(api("closeDIAStateAdapter")(manager), add = TRUE)
    manager$.context$reference_states <- TRUE
    manager$setWorkflowType("DIA")
    entries <- sort(roles$entry)
    design <- data.frame(Run = paste0("S", 1:6), group = rep(c("A", "B"), each = 3L),
        batch = rep(c("b1", "b2", "b3"), times = 2L), stringsAsFactors = FALSE)
    quant <- h$goPropQuant(entries, roles$entry[roles$role == "p_term"], roles$entry[roles$role == "d_term"])
    object <- methods::new("ProteinQuantitativeData", protein_quant_table = quant, design_matrix = design,
        sample_id = "Run", group_id = "group", protein_id_column = "uniprot_acc",
        protein_id_table = data.frame(uniprot_acc = entries, gene_names = paste0("GENE", seq_along(entries)),
            stringsAsFactors = FALSE),
        args = list(globalParameters = list(workflow_type = "DIA", use_limpa = FALSE, report_template = "DIA_report"),
            deAnalysisParameters = list(formula_string = "~ 0 + group")))
    api("setDIAStepProvenance")(manager, "design_builder")
    manager$saveState("correlation_filtered", object, list(), "fixture")
    manager$.context$experiment_paths <- paths
    contrasts <- data.frame(contrasts = "groupB-groupA", full_format = "B_vs_A=groupB-groupA",
        friendly_names = "B_vs_A", stringsAsFactors = FALSE)
    workflow <- list2env(list(state_manager = manager, contrasts_tbl = contrasts, tab_status = list(),
        qc_params = list(protein_qc = list()), config_list = list(), uniprot_dat_cln = annotation))
    parameters <- list(selected_contrast = "B_vs_A", organism_taxid = taxon, up_cutoff = 0, down_cutoff = 0,
        q_cutoff = 0.05, correction_method = "gSCS", enable_organism_filter = organism_filter)
    notes <- character()
    outcome <- "completed"
    message_text <- ""
    blank_gene_names <- NA_integer_
    da_rows <- NA_integer_
    utils::capture.output(tryCatch(suppressMessages(suppressWarnings(api("withDIAStateScope")(manager, function() {
        da_data <- new.env(parent = emptyenv())
        da_data$current_s4_object <- api("getDIAStateReference")(manager)
        da_data$contrasts_available <- contrasts$contrasts
        manifest <- api("runDIAReferenceDaAnalysis")(list(formula_string = "~ 0 + group", da_q_val_thresh = 0.05,
            treat_lfc_cutoff = 0), da_data, workflow, paths, contrasts = contrasts, publish = FALSE)
        api("produceDIADaPublication")(manager, manifest, annotation, 0.05,
            paths[c("da_output_dir", "publication_graphs_dir")], save_assets = FALSE)
        da <- utils::read.delim(list.files(paths$da_output_dir, "_long_annot\\.tsv$", full.names = TRUE)[[1L]],
            stringsAsFactors = FALSE)
        da_rows <<- nrow(da)
        if ("gene_name" %in% names(da)) blank_gene_names <<- sum(is.na(da$gene_name) | !nzchar(da$gene_name))
        api("runWorkflowEnrichment")(workflow, api("omicGuiConfig")("proteomics")$enrichment, manifest, parameters,
            paths, process_fn = api("processEnrichments"), notify = function(message, type, duration) {
                notes <<- c(notes, sprintf("[%s] %s", type, message))
            })
    }))), error = function(e) {
        outcome <<- paste0("error (", paste(class(e), collapse = ","), ")")
        message_text <<- conditionMessage(e)
    }))
    route <- if (identical(taxon, "999999")) "clusterProfiler" else "g:Profiler"
    service <- if (identical(route, "clusterProfiler")) calls$enricher else calls$gost
    data.frame(case = name, route = route, annotation_rows = nrow(annotation), organism_filter = organism_filter,
        da_rows = da_rows, da_blank_gene_name = blank_gene_names, outcome = outcome,
        service_calls = length(service),
        identifiers_sent = if (length(service)) paste(vapply(service, function(x) x[["gene"]], 0), collapse = ";") else "",
        background = if (length(service)) paste(vapply(service, function(x) x[["universe"]], 0), collapse = ";") else "",
        notices = paste(notes, collapse = " | "), message = gsub("[\r\n\t]+", " ", message_text),
        stringsAsFactors = FALSE)
}

rows <- list()
for (name in names(tables)) {
    rows[[length(rows) + 1L]] <- run_case(paste0("cp-", name), tables[[name]], "999999")
    rows[[length(rows) + 1L]] <- run_case(paste0("gp-", name), tables[[name]], "9606")
}
rows[[length(rows) + 1L]] <- run_case("cp-empty-filter", tables$empty, "999999", organism_filter = TRUE)
rows[[length(rows) + 1L]] <- run_case("cp-complete-filter", tables$complete, "999999", organism_filter = TRUE)
restore_enricher()
restore_gost()
result <- do.call(rbind, rows)

ids <- sort(h$goPropReadTsv("t7-fixture-entries.tsv")$entry)
features <- lapply(tables, function(table) api("protIntegrationInputFeatures")(ids, table))
integration <- data.frame(case = paste0("integration-", names(features)), route = "integration input",
    annotation_rows = vapply(tables, nrow, 0L), organism_filter = FALSE, da_rows = length(ids), da_blank_gene_name = NA,
    outcome = "completed", service_calls = 0L,
    identifiers_sent = vapply(features, function(f) as.character(sum(is.na(f$gene_names))), ""),
    background = vapply(features, function(f) as.character(sum(f$display_name == f$feature_id)), ""),
    notices = "", message = "identifiers_sent = features with NA gene_names; background = display names that are the id",
    stringsAsFactors = FALSE)
result <- rbind(result, integration)
utils::write.table(result, file.path(out, "f1-consumers.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
writeLines(c(paste("commit:", tryCatch(system2("git", c("rev-parse", "HEAD"), stdout = TRUE), error = function(e) NA)),
    paste("R:", R.version.string), paste("clusterProfiler:", utils::packageVersion("clusterProfiler")),
    paste("gprofiler2:", utils::packageVersion("gprofiler2")), paste("DOSE:", utils::packageVersion("DOSE"))),
    file.path(out, "f1-session.txt"))
print(result[, c("case", "outcome", "service_calls", "identifiers_sent", "background")])
