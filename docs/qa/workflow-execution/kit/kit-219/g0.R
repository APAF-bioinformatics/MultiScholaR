# KIT-219 G0 (ticket section 10, items (i) and (ii)): the files the real enrichment kernel writes on each route, against
# the names the rule derives from the stored result. Run on the base commit, from the package root, after it is installed.
#
#   Rscript docs/qa/workflow-execution/kit/kit-219/g0.R
#
# The project, the DA and the enrichment follow KIT-208's `k208Project` (tests/testthat/test-kit208-gprofiler-identifiers.R:
# 86-142) and KIT-173's `goPropRunProject` (tests/testthat/helper-go-propagation.R:118-166): the T7 fixture, the DA fitted
# and published, then `runWorkflowEnrichment` with the real kernel `processEnrichments` as the tab runs it. The g:Profiler
# route (taxon 9606) has `gprofiler2::gost` replaced by a fixed answer (no service call; the proxy points at a closed
# port), the GO route (taxon 999999) runs offline. Each route is run twice into one `pathway_dir`: the second run's down
# list gives no result (g:Profiler: `gost` answers NULL; GO: a down-regulation cut-off that no protein reaches).
#
# The rule is written here as a function (`g0Derive`): the package's function does not exist on the base commit.
#
# Output: docs/qa/workflow-execution/kit/kit-219/g0-out/g0.tsv (route, run, kind, value) and g0-out/g0-summary.txt.

out_dir <- file.path("docs", "qa", "workflow-execution", "kit", "kit-219", "g0-out")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
suppressPackageStartupMessages(library(MultiScholaR))
Sys.setenv(MULTISCHOLAR_TEST_MODE = "false", http_proxy = "http://127.0.0.1:9", https_proxy = "http://127.0.0.1:9",
    HTTP_PROXY = "http://127.0.0.1:9", HTTPS_PROXY = "http://127.0.0.1:9", no_proxy = "", NO_PROXY = "")
options(multischolar.test_mode = FALSE)
ns <- asNamespace("MultiScholaR")
h <- new.env(parent = ns)
sys.source(file.path("tests", "testthat", "helper-go-propagation.R"), envir = h)

g0Rows <- list()
g0Add <- function(route, run, kind, value) {
    value <- as.character(value)
    if (!length(value)) value <- "<none>"
    g0Rows[[length(g0Rows) + 1L]] <<- data.frame(route = route, run = run, kind = kind, value = value, stringsAsFactors = FALSE)
}

# The rule of 3.1 item 2, as a function of the stored result.
g0Derive <- function(state) {
    data <- methods::slot(state, "enrichment_data")
    as.character(unlist(lapply(names(data), function(contrast) {
        entry <- data[[contrast]]
        directions <- intersect(c("up", "down"), names(Filter(Negate(is.null), entry)))
        safe <- gsub("[^A-Za-z0-9_.-]", "_", contrast)
        unlist(lapply(directions, function(direction) c(paste0(contrast, "_", direction, "_enrichment_results.tsv"),
            paste0(safe, "_", direction, "_enrichment_plot.", c("html", "png")))))
    })))
}

work <- tempfile("kit219-g0-")
dir.create(work)
g0Annotations <- h$goPropAnnotate(h$goPropT7Raw(), file.path(work, "cache"))
g0Roles <- h$goPropReadTsv("t7-fixture-entries.tsv")
g0Pathway <- file.path(work, "pathway_enrichment")
dir.create(g0Pathway)

#' One run of the project into `g0Pathway`; returns what the kernel was given and what the store holds.
g0RunProject <- function(route, run, root, taxon, gost = NULL, down_cutoff = 0) {
    seen <- new.env()
    seen$da_names <- NULL
    seen$contrast_names <- NULL
    restore <- if (!is.null(gost)) h$goPropReplaceBinding("gprofiler2", "gost", gost)
    if (!is.null(restore)) on.exit(restore(), add = TRUE)
    paths <- h$goPropPaths(root)
    paths$pathway_dir <- g0Pathway
    for (p in paths) dir.create(p, recursive = TRUE, showWarnings = FALSE)
    saveRDS(g0Annotations, file.path(paths$source_dir, "uniprot_dat_cln.RDS"))
    manager <- newDIAStateAdapter(newWorkflowEntityStore(root, "kit219-g0", "DIA"),
        WorkflowState$new(audit_enabled = FALSE))
    on.exit(closeDIAStateAdapter(manager), add = TRUE)
    manager$.context$reference_states <- TRUE
    manager$setWorkflowType("DIA")
    entries <- sort(g0Roles$entry)
    design <- data.frame(Run = paste0("S", 1:6), group = rep(c("A", "B"), each = 3L),
        batch = rep(c("b1", "b2", "b3"), times = 2L), stringsAsFactors = FALSE)
    quant <- h$goPropQuant(entries, g0Roles$entry[g0Roles$role == "p_term"], g0Roles$entry[g0Roles$role == "d_term"])
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
    parameters <- list(selected_contrast = "B_vs_A", organism_taxid = taxon, up_cutoff = 0, down_cutoff = down_cutoff,
        q_cutoff = 0.05, correction_method = "gSCS", enable_organism_filter = FALSE)
    kernel <- function(da_results, ...) {
        seen$da_names <- names(da_results@da_data)
        seen$contrast_names <- list(...)$contrast_names
        processEnrichments(da_results, ...)
    }
    state <- NULL
    package_current <- NULL
    manifest_names <- NULL
    error <- NULL
    utils::capture.output(tryCatch(suppressMessages(suppressWarnings(withDIAStateScope(manager, function() {
        da_data <- new.env(parent = emptyenv())
        da_data$current_s4_object <- getDIAStateReference(manager)
        da_data$contrasts_available <- contrasts$contrasts
        manifest <- runDIAReferenceDaAnalysis(list(formula_string = "~ 0 + group", da_q_val_thresh = 0.05,
            treat_lfc_cutoff = 0), da_data, workflow, paths, contrasts = contrasts, publish = FALSE)
        produceDIADaPublication(manager, manifest, g0Annotations, 0.05, paths[c("da_output_dir", "publication_graphs_dir")],
            save_assets = FALSE)
        manifest_names <<- names(manifest$results)
        runWorkflowEnrichment(workflow, omicGuiConfig("proteomics")$enrichment, manifest, parameters, paths,
            process_fn = kernel)
        state <<- manager$getState("enrichment_completed")
        # On a commit that has the package's rule (G4's check on head): its list for this store.
        if (exists("workflowEnrichmentCurrentFiles", envir = ns, inherits = FALSE)) {
            package_current <<- get("workflowEnrichmentCurrentFiles", envir = ns)(workflow, "proteomics")
        }
    }))), error = function(e) error <<- conditionMessage(e)))
    g0Add(route, run, "package_current_files", if (exists("workflowEnrichmentCurrentFiles", envir = ns, inherits = FALSE))
        sort(package_current) else "<function absent>")
    g0Add(route, run, "package_equals_derived", if (is.null(state) || !exists("workflowEnrichmentCurrentFiles",
        envir = ns, inherits = FALSE)) NA else setequal(package_current, g0Derive(state)))
    g0Add(route, run, "error", if (is.null(error)) "<none>" else error)
    g0Add(route, run, "source_da_data_names", seen$da_names)
    g0Add(route, run, "source_contrast_names_argument", seen$contrast_names)
    g0Add(route, run, "source_manifest_results_names", manifest_names)
    state
}

# The fixed g:Profiler answer: the first call of a contrast is the up list, the later ones the down list.
g0MakeGost <- function(down_empty) {
    calls <- 0L
    function(query, ..., custom_bg = NULL) {
        calls <<- calls + 1L
        # The kernel retries a NULL answer (up to five attempts), so once the up list has its answer every later call is NULL.
        if (down_empty && calls >= 2L) return(NULL)
        buildProtEnrichDeterministicGprofilerResult(if (calls %% 2L == 1L) "positive" else "negative")
    }
}

# The package's internals are visible from the helper environment (its parent is the namespace).
environment(g0RunProject) <- h
environment(g0MakeGost) <- h

g0Report <- function(route, run, state) {
    written <- sort(list.files(g0Pathway, recursive = TRUE))
    derived <- if (is.null(state)) character() else g0Derive(state)
    g0Add(route, run, "written", written)
    g0Add(route, run, "derived", sort(derived))
    g0Add(route, run, "written_not_derived", sort(setdiff(written, derived)))
    g0Add(route, run, "derived_not_written", sort(setdiff(derived, written)))
    data <- if (is.null(state)) list() else methods::slot(state, "enrichment_data")
    g0Add(route, run, "enrichment_data_names", names(data))
    g0Add(route, run, "directions_with_result", unlist(lapply(names(data), function(contrast) paste0(contrast, ":",
        paste(names(Filter(Negate(is.null), data[[contrast]])), collapse = "+")))))
}

for (route in c("gprofiler", "go")) {
    taxon <- if (identical(route, "gprofiler")) "9606" else "999999"
    unlink(list.files(g0Pathway, full.names = TRUE, recursive = TRUE))
    state1 <- g0RunProject(route, "run1", file.path(work, paste0(route, "-1")), taxon,
        gost = if (identical(route, "gprofiler")) g0MakeGost(FALSE))
    g0Report(route, "run1", state1)
    first_files <- sort(list.files(g0Pathway, recursive = TRUE))
    state2 <- g0RunProject(route, "run2", file.path(work, paste0(route, "-2")), taxon,
        gost = if (identical(route, "gprofiler")) g0MakeGost(TRUE), down_cutoff = if (identical(route, "go")) 100 else 0)
    g0Report(route, "run2", state2)
    # (ii): the files of the first run that the second run's names exclude, and are still on disk.
    left <- setdiff(first_files, if (is.null(state2)) character() else g0Derive(state2))
    g0Add(route, "run2", "first_run_files_not_in_second_run_names", sort(left))
    g0Add(route, "run2", "of_which_still_on_disk", sort(intersect(left, list.files(g0Pathway, recursive = TRUE))))
}

table <- do.call(rbind, g0Rows)
utils::write.table(table, file.path(out_dir, "g0.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
summary <- vapply(split(table, list(table$route, table$run), drop = TRUE), function(part) {
    pick <- function(kind) part$value[part$kind == kind]
    sprintf("%s %s: error=%s written=%d derived=%d written_not_derived=%s derived_not_written=%s", part$route[[1L]],
        part$run[[1L]], paste(pick("error"), collapse = ","), length(setdiff(pick("written"), "<none>")),
        length(setdiff(pick("derived"), "<none>")), paste(pick("written_not_derived"), collapse = ","),
        paste(pick("derived_not_written"), collapse = ","))
}, "")
writeLines(summary, file.path(out_dir, "g0-summary.txt"))
cat(summary, sep = "\n")
cat("g0 done\n")
