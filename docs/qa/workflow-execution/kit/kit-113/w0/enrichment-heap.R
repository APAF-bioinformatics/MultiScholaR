# KIT-113 W0 (M0p): the enrichment reachability probe (D8). Scratch, not a test; one counted attempt.
#   Rscript docs/qa/workflow-execution/kit/kit-113/w0/enrichment-heap.R <package source> <report file> [plan]
# `plan` is "counted" (default: D8's sizes, the decision case first) or "dry" (a throw-away debug run at a tiny size).
#
# Method: KIT-112 W1's (kit-112/w1/request-half-hydration.R). A synthetic DIA owner at 300 samples on an on-disk store,
# an in-process DA with its publication as the replay runs it (replayWorkflowEnrichmentStep,
# utils_workflow_enrichment.R:1102-1134: runDIAReferenceDaAnalysis(publish = FALSE), then produceDIADaPublication(...,
# save_assets = FALSE) into the run's DA directory). Then each statement group of D1's table, with `gc(reset = TRUE)`
# before the group and the heap's "max used" after it, the group's seconds and the size of the value it produces.
#   - In process: runWorkflowEnrichment's body (:822-925), copied verbatim, group by group, in today's order.
#   - Managed controller: the "stays" groups in the managed order (D3), emulated with today's statements and functions.
#     First the request half (the context, the annotation resolve, the mapping, the method, the kernel's name, the
#     annotation table saved as a transfer value), with no DA table read. Then the worker's part runs in this process,
#     outside the controller's figure: the source, the DA tables, the match, the run and product, the portable results
#     saved as an `rds` state value. Then the moved values are removed and collected. Last, the hydration of the
#     results from that value and the accept half (the records, the bus's filter metadata, the display tables, the
#     payload, the stamp, the UI parameters, the state commit, the status).
# Each sequence starts from a collected heap. Its figure is its highest "max used" minus the live heap at its start.
# Each sequence runs twice: pass 1 (the first Run in the process: "cold") and pass 2 (a second Run: "warm").
# Heap bytes are Ncells x 56 + Vcells x 8. GiB = 2^30 bytes.
# Network: the lane's closed-port proxy (helper-gui-entry.R:275-276). The e2e doubles (helper-e2e-enrichment-doubles.R,
# sourced) replace gprofiler2::gost and clusterProfiler::enricher. Tracers count every call into the HTTP entry points.
# None may be reached.

args <- commandArgs(TRUE)
tree <- normalizePath(args[[1L]])
report_path <- args[[2L]]
plan <- if (length(args) >= 3L) args[[3L]] else "counted"
started <- proc.time()[["elapsed"]]
ceiling_seconds <- 1800
Sys.setenv(http_proxy = "http://127.0.0.1:9", https_proxy = "http://127.0.0.1:9", HTTP_PROXY = "http://127.0.0.1:9",
    HTTPS_PROXY = "http://127.0.0.1:9", no_proxy = "", NO_PROXY = "")
setwd(tree)
suppressMessages(pkgload::load_all(tree, quiet = TRUE, export_all = TRUE, helpers = FALSE))
suppressPackageStartupMessages(library(testthat))
sys.source(file.path(tree, "tests", "testthat", "helper-e2e-enrichment-doubles.R"), envir = globalenv())
ns <- asNamespace("MultiScholaR")

report <- file(report_path, open = "w")
say <- function(...) {
    cat(sprintf(...), "\n", sep = "", file = report)
    flush(report)
}
say("KIT-113 W0 M0p enrichment-heap probe; plan %s; R %s; tree %s; HEAD %s", plan, getRversion(), tree,
    system2("git", c("-C", tree, "rev-parse", "--short", "HEAD"), stdout = TRUE))
say("started %s; TMPDIR %s", format(Sys.time()), tempdir())

# ---- the network guard ------------------------------------------------------------------------------------------
net_calls <- new.env(parent = emptyenv())
.probe_net <- function(name) assign(name, (if (exists(name, envir = net_calls)) get(name, envir = net_calls) else 0L) + 1L,
    envir = net_calls)
traceCount <- function(what, where, prefix = "") {
    if (!exists(what, envir = where, inherits = FALSE)) return(invisible(FALSE))
    tracer <- bquote(base::get(".probe_net", envir = .GlobalEnv)(.(paste0(prefix, environmentName(where), "::", what))))
    suppressMessages(trace(what, tracer = tracer, print = FALSE, where = where))
    invisible(TRUE)
}
for (package in c("curl", "RCurl", "httr", "httr2")) {
    if (!requireNamespace(package, quietly = TRUE)) next
    where <- asNamespace(package)
    for (what in c("curl_fetch_memory", "curl_fetch_disk", "curl_fetch_stream", "curl_fetch_multi", "curl_download",
        "curl", "getURL", "getURLContent", "postForm", "curlPerform", "req_perform")) traceCount(what, where)
}
traceCount("download.file", asNamespace("utils"))
traceCount("url", baseenv())
traceCount("getUniprotAnnotations", ns)
install_enrichment_doubles(.local_envir = globalenv())
# The doubles are in place; count the calls that reach them.
traceCount("gost", asNamespace("gprofiler2"), prefix = "offline double ")
traceCount("enricher", asNamespace("clusterProfiler"), prefix = "offline double ")
# Every counter except "offline double ..." is an HTTP entry point (or the UniProt fallback): each must stay at zero.
networkCalls <- function() {
    names <- ls(net_calls)
    if (!length(names)) return("none")
    paste(sprintf("%s=%d", names, vapply(names, function(n) get(n, envir = net_calls), 1L)), collapse = ", ")
}

# ---- the heap measures --------------------------------------------------------------------------------------------
heapBytes <- function(g, column) g[1L, column] * 56 + g[2L, column] * 8
mib <- function(bytes) bytes / 2^20
valueSize <- function(names, env) {
    names <- intersect(names, ls(env, all.names = TRUE))
    if (!length(names)) return(NA_real_)
    sum(vapply(names, function(n) as.numeric(lobstr::obj_size(get(n, envir = env))), 1))
}
newSequence <- function(label, initial) {
    invisible(gc(full = TRUE))
    g <- gc(reset = TRUE)
    seq <- new.env(parent = emptyenv())
    seq$label <- label
    seq$start <- heapBytes(g, 1L)
    seq$max <- seq$start
    seq$env <- list2env(initial, parent = ns)
    seq$rows <- list()
    say("  [%s] start: live heap %.1f MiB", label, mib(seq$start))
    seq
}
#' One statement group: run `code` in the sequence's frame, between `gc(reset = TRUE)` and `gc()`.
#' `counted = FALSE` (the emulated worker) keeps the group out of the controller's figure.
group <- function(seq, id, label, code, values = character(), counted = TRUE) {
    g0 <- gc(reset = TRUE)
    before <- heapBytes(g0, 1L)
    t0 <- proc.time()[["elapsed"]]
    eval(code, seq$env)
    seconds <- proc.time()[["elapsed"]] - t0
    g1 <- gc()
    max_used <- heapBytes(g1, which(colnames(g1) == "max used")[[1L]])
    if (counted) seq$max <- max(seq$max, max_used)
    size <- valueSize(values, seq$env)
    row <- data.frame(sequence = seq$label, id = id, group = label, counted = counted,
        peak_over_start_mib = mib(max_used - seq$start), peak_over_before_mib = mib(max_used - before),
        live_after_over_start_mib = mib(heapBytes(g1, 1L) - seq$start), value_mib = mib(size), seconds = seconds)
    seq$rows[[length(seq$rows) + 1L]] <- row
    say("    %-4s %-62s %s peak +%8.1f MiB over start (+%7.1f over before), live after +%8.1f, value %8.2f MiB, %6.1f s",
        id, label, if (counted) " " else "w", row$peak_over_start_mib, row$peak_over_before_mib,
        row$live_after_over_start_mib, row$value_mib, seconds)
    invisible(row)
}
figure <- function(seq) {
    value <- seq$max - seq$start
    say("  [%s] FIGURE: highest max used - live heap at start = %.1f MiB = %.4f GiB", seq$label, mib(value), value / 2^30)
    value
}

# ---- the run's statements, copied from runWorkflowEnrichment (utils_workflow_enrichment.R:822-925) ----------------
statements <- list(
    context = quote({ # :822-845
        notify <- notify %||% function(message, type, duration) invisible(NULL)
        progress <- progress %||% function(value, detail) invisible(NULL)
        if (!isTRUE(config$da_source$manifest(manifest))) {
            rlang::abort("Enrichment needs an accepted DA result manifest.", class = "workflow_enrichment_input")
        }
        manager <- workflow$state_manager
        source_object <- readDIAExecutionComponents(manager, manifest$source, fields = config$da_source$components)
        selected <- parameters$selected_contrast
        progress(0.2, "Transforming DE data...")
        progress(0.3, "Creating DE results S4 object...")
        contrasts_tbl <- resolveProtEnrichContrastsTable(workflow)
        selection <- resolveProtEnrichSelectedDaResults(selected, manifest$results, contrasts_tbl)
        if (is.null(selection$selectedDaResults)) {
            rlang::abort(sprintf("Could not find DE results for contrast '%s'. Available contrasts: %s", selected,
                paste(selection$availableKeys, collapse = ", ")), class = "workflow_enrichment_input")
        }
        dependencies <- resolveProtEnrichRunDependencies(source_object, manifest$results, workflow, contrasts_tbl)
        contrasts_tbl <- dependencies$contrastsTbl
        if (is.null(contrasts_tbl) || is.null(dependencies$designMatrix)) {
            rlang::abort(paste("Missing contrasts_tbl or design_matrix required for enrichment analysis.",
                "Please ensure DE analysis was completed successfully."), class = "workflow_enrichment_input")
        }
        directories <- resolveProtEnrichOutputDirectories(paths)
    }),
    da_tables = quote({ # :846-852
        da_results <- config$da_source$read(contrasts_tbl = contrasts_tbl, design_matrix = dependencies$designMatrix,
            da_output_dir = directories$daOutputDir)
        missing <- names(Filter(is.null, da_results@da_data))
        if (isTRUE(refuse_missing_da) && length(missing)) {
            rlang::abort(sprintf("The DA results table is missing for contrast %s.", paste(missing, collapse = ", ")),
                class = "workflow_enrichment_input")
        }
    }),
    resolve = quote({ # :853-854
        annotations <- config$annotation$resolve(workflowData = workflow, experimentPaths = paths,
            currentS4Object = source_object, organismTaxid = parameters$organism_taxid)$uniprotDatCln
    }),
    match = quote({ # :855-856
        matching <- config$annotation$match(uniprotDatCln = annotations, daResultsForEnrichment = da_results,
            currentS4Object = source_object)
    }),
    filter = quote({ # :857-876
        filter_applied <- FALSE
        filter_stats <- list(proteins_before = 0, proteins_after = 0, proteins_removed = 0)
        if (isTRUE(parameters$enable_organism_filter) && !is.null(config$organism_filter)) {
            target <- as.character(parameters$organism_taxid)
            mapping <- config$organism_filter$mapping(workflowData = workflow, uniprotDatCln = annotations,
                targetTaxon = target)$organismMapping
            if (!is.null(mapping) && nrow(mapping) > 0) {
                filtered <- config$organism_filter$apply(daResultsForEnrichment = da_results, organismMapping = mapping,
                    targetTaxon = target, currentS4Object = source_object)
                da_results <- filtered$daResultsForEnrichment
                filter_applied <- filtered$filterApplied
                filter_stats <- filtered$filterStats
            } else {
                notify("Multi-species filtering enabled but no organism mapping available. Proceeding without filtering.",
                    "warning", 8)
            }
        }
        organism_filter <- persistProtEnrichOrganismFilterMetadata(workflowData = workflow,
            organismFilterEnabled = parameters$enable_organism_filter, organismFilterApplied = filter_applied,
            targetTaxonId = parameters$organism_taxid, filterStats = filter_stats)
    }),
    method = quote({ # :877-880
        progress(0.5, "Running enrichment analysis...")
        method <- config[["organism"]]$method(organismTaxid = parameters$organism_taxid,
            supportedOrganisms = config[["organism"]]$supported())$methodInfo
    }),
    kernel_inputs = quote({ # :881-889
        columns <- resolveProtEnrichAnalysisInputColumns(methodInfo = method, daResultsForEnrichment = da_results,
            currentS4Object = source_object)
        arguments <- buildProtEnrichProcessEnrichmentsArgs(daResultsForEnrichment = da_results,
            organismTaxid = parameters$organism_taxid, upCutoff = parameters$up_cutoff, downCutoff = parameters$down_cutoff,
            qCutoff = parameters$q_cutoff, pathwayDir = directories$pathwayDir, goAnnotations = annotations,
            proteinIdColumn = columns$idColumn, contrastNames = names(manifest$results),
            correctionMethod = parameters$correction_method, excludeIea = config$exclude_iea)
        kernel <- process_fn %||% config$kernel$select()
        .capture_checkpoint(arguments$checkpointArgs, "cp10", "enrichment_input")
    }),
    run = quote({ # :890
        declared <- config$run(workflow, manifest, arguments$processArgs, process_fn = kernel)
    }),
    display_tables = quote({ # :891-900
        results <- declared$results
        contrast_results <- buildProtEnrichAllContrastResults(enrichmentResults = results, methodInfo = method)
        for (contrast in names(contrast_results)) {
            contrast_results[[contrast]]$gprofiler_results_table <- saveDIADisplayTable(manager,
                contrast_results[[contrast]]$gprofiler_results)
            contrast_results[[contrast]]$clusterprofiler_results_table <- saveDIADisplayTable(manager,
                contrast_results[[contrast]]$clusterprofiler_results)
        }
        shown <- resolveProtEnrichSelectedContrastResults(selected, contrast_results, contrasts_tbl)
    }),
    payload = quote({ # :901-907
        progress(0.8, "Storing results...")
        workflow$enrichment_analysis_results <- buildProtEnrichAnalysisResultsPayload(
            gprofilerResults = if (shown$found) shown$gprofilerResults,
            clusterprofilerResults = if (shown$found) shown$clusterprofilerResults,
            stringdbResults = if (shown$found) shown$stringdbResults, fullEnrichmentResults = results,
            selectedContrast = selected, analysisMethod = method$method, upCutoff = parameters$up_cutoff,
            downCutoff = parameters$down_cutoff, qCutoff = parameters$q_cutoff, organismTaxid = parameters$organism_taxid)
    }),
    stamp = quote({ # :908-911
        stamped <- propagateProtEnrichResultsArgs(enrichmentResults = results, currentS4Object = source_object,
            selectedContrast = selected, methodInfo = method, upCutoff = parameters$up_cutoff,
            downCutoff = parameters$down_cutoff, qCutoff = parameters$q_cutoff, organismTaxid = parameters$organism_taxid,
            pathwayDir = directories$pathwayDir)$enrichmentResults
    }),
    ui_parameters = quote({ # :912-915
        ui_parameters <- propagateProtEnrichUiParams(currentS4Object = source_object, workflowData = workflow,
            selectedContrast = selected, methodInfo = method, upCutoff = parameters$up_cutoff,
            downCutoff = parameters$down_cutoff, qCutoff = parameters$q_cutoff, organismTaxid = parameters$organism_taxid)
        updateProtEnrichStateManagerUiParams(workflowData = workflow, storedUiParams = ui_parameters$storedUiParams)
    }),
    commit = quote({ # :916-917
        saved <- saveWorkflowEnrichmentState(workflow, config, stamped, declared$parameters, method, parameters,
            directories$pathwayDir)
    }),
    status = quote({ # :918-925
        status <- workflow$tab_status
        status[[config$status_key]] <- "complete"
        workflow$tab_status <- status
        progress(1.0, "Complete!")
        out <- list(results = results, method = method, contrasts = contrast_results, selected = shown,
            raw_contrast = selection$rawContrastName, contrasts_tbl = contrasts_tbl, pathway_dir = directories$pathwayDir,
            annotation_match = matching$annotationMatchResults, organism_filter = organism_filter, service = declared$service,
            reference = declared$reference, parameters = declared$parameters, saved = saved)
    })
)

# ---- the managed emulation's own statements (D3's groupings over today's code) ------------------------------------
managed <- list(
    mapping = quote({ # D3 request: the mapping (bus) and its warning when the filter is on (:859-872, without the application)
        mapping <- NULL
        if (isTRUE(parameters$enable_organism_filter) && !is.null(config$organism_filter)) {
            mapping <- config$organism_filter$mapping(workflowData = workflow, uniprotDatCln = annotations,
                targetTaxon = as.character(parameters$organism_taxid))$organismMapping
            if (is.null(mapping) || nrow(mapping) == 0)
                notify("Multi-species filtering enabled but no organism mapping available. Proceeding without filtering.",
                    "warning", 8)
        }
    }),
    kernel_name = quote({ # D3 request: the kernel's name, looked up by name when called
        kernel <- process_fn %||% config$kernel$select()
        kernel_name <- if (identical(kernel, get("runProtEnrichDeterministicProcessEnrichments", envir = ns)))
            "runProtEnrichDeterministicProcessEnrichments" else if (identical(kernel, get("processEnrichments", envir = ns)))
            "processEnrichments" else stop("unknown kernel")
    }),
    transfer = quote({ # D3 the transfer values: the annotation table (and the mapping) as saved display tables
        annotations_reference <- if (is.data.frame(annotations)) saveDIADisplayTable(manager, annotations)
        mapping_reference <- if (is.data.frame(mapping)) saveDIADisplayTable(manager, mapping)
        request <- list(manifest = manifest, contrasts_tbl = contrasts_tbl, parameters = parameters,
            directories = directories, annotations = annotations_reference, mapping = mapping_reference, method = method,
            kernel = kernel_name, design = dependencies$designMatrix)
    }),
    # The emulated worker (not counted): the compute half of D3.
    w_source = quote({
        w_source <- readDIAExecutionComponents(manager, request$manifest$source, fields = config$da_source$components)
        w_annotations <- if (!is.null(request$annotations)) materializeWorkflowTable(manager$.context, request$annotations)
    }),
    w_da = quote({
        w_da <- config$da_source$read(contrasts_tbl = request$contrasts_tbl, design_matrix = request$design,
            da_output_dir = request$directories$daOutputDir)
        if (isTRUE(refuse_missing_da) && length(names(Filter(is.null, w_da@da_data))))
            rlang::abort("The DA results table is missing", class = "workflow_enrichment_input")
    }),
    w_match = quote({
        w_matching <- config$annotation$match(uniprotDatCln = w_annotations, daResultsForEnrichment = w_da,
            currentS4Object = w_source)
        w_filter <- list(applied = FALSE, stats = list(proteins_before = 0, proteins_after = 0, proteins_removed = 0))
    }),
    w_run = quote({ # workflowEnrichmentRun with computeDIAReferenceEnrichment (utils_workflow_dia_da.R:657-676 at HEAD)
        w_columns <- resolveProtEnrichAnalysisInputColumns(methodInfo = request$method, daResultsForEnrichment = w_da,
            currentS4Object = w_source)
        w_arguments <- buildProtEnrichProcessEnrichmentsArgs(daResultsForEnrichment = w_da,
            organismTaxid = request$parameters$organism_taxid, upCutoff = request$parameters$up_cutoff,
            downCutoff = request$parameters$down_cutoff, qCutoff = request$parameters$q_cutoff,
            pathwayDir = request$directories$pathwayDir, goAnnotations = w_annotations,
            proteinIdColumn = w_columns$idColumn, contrastNames = names(request$manifest$results),
            correctionMethod = request$parameters$correction_method, excludeIea = config$exclude_iea)
        w_kernel <- get(request$kernel, envir = ns)
        w_manifest <- request$manifest
        w_declaration <- diaEnrichmentOperation()
        w_results_ids <- vapply(names(w_manifest$results), function(name)
            w_manifest$results[[name]]$tables$da_proteins_long$entity_id %||% NA_character_, character(1))
        w_parameters <- c(w_arguments$processArgs[setdiff(w_declaration$parameters, "results")],
            list(results = w_results_ids, da_results = w_arguments$processArgs$da_results,
                pathway_dir = w_arguments$processArgs$pathway_dir, go_annotations = w_arguments$processArgs$go_annotations,
                process_fn = w_kernel))
        w_value <- runProjected(manager, w_declaration, w_manifest$source, w_parameters)
        w_identity <- diaComponentScientificIdentity(manager, w_manifest$source, "enrichment",
            c(w_parameters[w_declaration$parameters], list(service = w_value$service)))
        w_product <- list(schema = "dia_enrichment_product_v1", service = w_value$service,
            parameters = w_parameters[setdiff(w_declaration$parameters, "results")], results = w_results_ids,
            contrasts = names(w_value$results@enrichment_data), summaries = w_value$results@enrichment_summaries)
        w_reference <- saveWorkflowConsumerEntity(manager$.context$store,
            newSavedWorkflowScientificProduct(w_product, w_declaration$product_type,
                provenance = list(service = w_value$service, cacheable = FALSE)),
            basename(tempfile("dia-enrichment-")), w_identity, "scientific_product", factors = w_declaration$factors)
        w_tables <- lapply(w_manifest$results, function(result) result$tables$da_proteins_long)
        rootDIAProduct(manager, w_reference, c(list(w_manifest$source), w_tables))
        w_computed <- list(results = w_value$results, service = w_value$service, reference = w_reference,
            parameters = w_parameters[w_declaration$parameters])
    }),
    w_save = quote({ # the results made portable (F107-05) and saved as an rds state value
        w_portable <- portableWorkflowEnrichmentResults(w_computed$results)
        portable_serialized <- length(serialize(w_portable, NULL))
        results_reference <- saveWorkflowStateValue(manager$.context, w_portable, "enrichment-results", codec = "rds")
        delivered <- list(results = results_reference, service = w_computed$service, reference = w_computed$reference,
            parameters = w_computed$parameters, matching = w_matching$annotationMatchResults, filter = w_filter)
    }),
    # The accept half (counted): hydration, then the controller's statements.
    hydrate = quote({
        results_hydrated <- resolveWorkflowStateValue(manager$.context, delivered$results)
        declared <- list(results = results_hydrated, service = delivered$service, reference = delivered$reference,
            parameters = delivered$parameters)
        matching <- list(annotationMatchResults = delivered$matching)
    }),
    records = quote({ # clearWorkflowStepProvenance, then recordDIAEnrichmentResources (dia_da.R:677-683 at HEAD)
        clearWorkflowStepProvenance(manager)
        for (service in setdiff(declared$service, "unavailable")) {
            parts <- strsplit(service, ":", fixed = TRUE)[[1L]]
            recordWorkflowExternalResource(manager, parts[[1L]], paste(parts[-1L], collapse = ":"),
                if (identical(parts[[1L]], "clusterProfiler")) "package_version" else "database_version",
                query_summary = sprintf("enrichment of %d contrasts, taxon %s", length(declared$parameters$results),
                    declared$parameters$taxon_id %||% NA))
        }
    }),
    bus_filter = quote({ # :874-876 with the worker's outcome
        organism_filter <- persistProtEnrichOrganismFilterMetadata(workflowData = workflow,
            organismFilterEnabled = parameters$enable_organism_filter, organismFilterApplied = delivered$filter$applied,
            targetTaxonId = parameters$organism_taxid, filterStats = delivered$filter$stats)
    })
)

# ---- the fixture --------------------------------------------------------------------------------------------------
buildCase <- function(proteins, samples, n_contrasts, root) {
    set.seed(1L)
    levels <- if (n_contrasts == 1L) c("A", "B") else c("A", "B", "C")
    runs <- sprintf("run_%04d", seq_len(samples))
    groups <- rep(levels, each = samples / length(levels))
    design <- data.frame(Run = runs, group = factor(groups, levels = levels),
        tech_rep_group = runs, stringsAsFactors = FALSE)
    # One technical-replicate group per sample: no duplicateCorrelation. With KIT-112 W1's pairs, its per-protein REML
    # added about 10 s at 400 proteins (dry runs 2 and 3: DA 155 s against 145 s), so several minutes at 16,000, inside
    # the 30-minute ceiling. The DA is the probe's setup, not measured, and enrichment reads only its published tables.
    ids <- sprintf("PG_%06d", seq_len(proteins))
    base <- stats::rnorm(proteins, 23, 2)
    # 10% of the proteins up and 10% down in B (and a different 20% in C), so each direction has a significant set.
    effect <- list(A = rep(0, proteins), B = c(rep(1.5, proteins %/% 10), rep(-1.5, proteins %/% 10),
        rep(0, proteins - 2L * (proteins %/% 10))))
    effect$C <- rev(effect$B)
    quant <- as.data.frame(c(list(Protein.Group = ids), stats::setNames(lapply(seq_len(samples), function(i)
        base + effect[[groups[[i]]]] + stats::rnorm(proteins, 0, 1)), runs)), check.names = FALSE,
        stringsAsFactors = FALSE)
    object <- methods::new("ProteinQuantitativeData", protein_quant_table = quant, protein_id_column = "Protein.Group",
        protein_id_table = data.frame(Protein.Group = ids, gene_names = sprintf("GENE%06d", seq_len(proteins)),
            stringsAsFactors = FALSE), design_matrix = design,
        sample_id = "Run", group_id = "group", technical_replicate_id = "tech_rep_group",
        args = list(globalParameters = list(workflow_type = "DIA")))
    annotations <- data.frame(Entry = ids, gene_names = sprintf("GENE%06d", seq_len(proteins)),
        Organism = "Homo sapiens", go_id_go_biological_process = "GO:0006915; GO:0007049",
        go_term_go_biological_process = "apoptotic process; cell cycle", go_id_go_molecular_function = "GO:0005515",
        go_term_go_molecular_function = "protein binding", go_id_go_cellular_compartment = "GO:0005739",
        go_term_go_cellular_compartment = "mitochondrion", stringsAsFactors = FALSE)
    pairs <- if (n_contrasts == 1L) list(c("B", "A")) else list(c("B", "A"), c("C", "A"), c("C", "B"))
    contrasts <- data.frame(contrasts = vapply(pairs, function(p) sprintf("group%s-group%s", p[[1L]], p[[2L]]), ""),
        full_format = vapply(pairs, function(p) sprintf("%s_vs_%s=group%s-group%s", p[[1L]], p[[2L]], p[[1L]], p[[2L]]), ""),
        friendly_names = vapply(pairs, function(p) sprintf("%s_vs_%s", p[[1L]], p[[2L]]), ""), stringsAsFactors = FALSE)
    results <- file.path(root, "results")
    paths <- list(source_dir = file.path(root, "source"), results_dir = results,
        da_output_dir = file.path(results, "da_proteins"), pathway_dir = file.path(results, "pathway_enrichment"),
        publication_graphs_dir = file.path(root, "publication_graphs"))
    for (path in paths) dir.create(path, recursive = TRUE, showWarnings = FALSE)
    manager <- newDIAStateAdapter(newWorkflowEntityStore(root, "k113-enrichment", "DIA"),
        WorkflowState$new(audit_enabled = FALSE))
    manager$.context$reference_states <- TRUE
    manager$.context$experiment_paths <- paths
    manager$setWorkflowType("DIA")
    manager$saveState("correlation_filtered", object, list(), "KIT-113 W0 synthetic state")
    say("  fixture: %d proteins x %d samples, %d groups, %d contrasts; protein table %.1f MiB in R", proteins, samples,
        length(levels), nrow(contrasts), mib(as.numeric(utils::object.size(quant))))
    rm(object, quant)
    # The DA as the replay runs it (utils_workflow_analysis_record.R:633-634), then its publication (:1110).
    t0 <- proc.time()[["elapsed"]]
    da_data <- new.env(parent = emptyenv())
    da_data$current_s4_object <- getWorkflowStateReference(manager)
    da_data$contrasts_available <- contrasts$contrasts
    da_workflow <- new.env(parent = emptyenv())
    da_workflow$state_manager <- manager
    da_workflow$contrasts_tbl <- contrasts
    manifest <- runDIAReferenceDaAnalysis(list(formula_string = "~ 0 + group", da_q_val_thresh = 0.05,
        treat_lfc_cutoff = 0), da_data, da_workflow, paths, contrasts = contrasts, annotations = annotations,
        publish = FALSE)
    t1 <- proc.time()[["elapsed"]]
    withWorkflowStateScope(manager, function() produceDIADaPublication(manager, manifest, annotations,
        manifest$parameters$qvalue, paths[c("da_output_dir", "publication_graphs_dir")], save_assets = FALSE))
    t2 <- proc.time()[["elapsed"]]
    files <- list.files(paths$da_output_dir, pattern = "_long_annot[.]tsv$", full.names = TRUE)
    say("  DA %.1f s, publication %.1f s; DA tables: %s", t1 - t0, t2 - t1,
        paste(sprintf("%s %.1f MB", basename(files), file.size(files) / 1e6), collapse = "; "))
    list(manager = manager, manifest = manifest, contrasts = contrasts, annotations = annotations, paths = paths)
}

#' The bus as the proteomics shell keeps it after the DA (the fields the run reads and writes), fresh per sequence.
newBus <- function(case) {
    workflow <- new.env(parent = emptyenv())
    workflow$state_manager <- case$manager
    workflow$contrasts_tbl <- case$contrasts
    workflow$tab_status <- list(differential_expression = "complete", enrichment_analysis = "pending")
    workflow$uniprot_dat_cln <- case$annotations
    workflow$da_analysis_results_list <- case$manifest
    workflow$enrichment_analysis_results <- NULL
    workflow$enrichment_ui_params <- NULL
    workflow$enrichment_organism_filter <- NULL
    workflow
}

initialFrame <- function(case, kernel) list(workflow = newBus(case), config = omicGuiConfig("proteomics")$enrichment,
    manifest = case$manifest, parameters = list(selected_contrast = case$contrasts$friendly_names[[1L]],
        organism_taxid = "9606", up_cutoff = 0, down_cutoff = 0, q_cutoff = 0.05, correction_method = "gSCS",
        enable_organism_filter = FALSE), paths = case$paths, process_fn = kernel, notify = NULL, progress = NULL,
    refuse_missing_da = TRUE, ns = ns)

inProcessSequence <- function(case, kernel, label) {
    seq <- newSequence(label, initialFrame(case, kernel))
    group(seq, "G1", "context (:822-845)", statements$context, c("source_object", "dependencies", "contrasts_tbl"))
    group(seq, "G2", "DA tables + refusal (:846-852) [moves]", statements$da_tables, "da_results")
    group(seq, "G3", "annotation resolve (:853-854)", statements$resolve, "annotations")
    group(seq, "G4", "annotation match (:855-856) [moves]", statements$match, "matching")
    group(seq, "G5", "organism filter, bus metadata (:857-876) [filter off]", statements$filter, "organism_filter")
    group(seq, "G6", "progress 0.5, method (:877-880)", statements$method, "method")
    group(seq, "G7", "id column, kernel args, kernel, cp10 (:881-889)", statements$kernel_inputs, "arguments")
    group(seq, "G8", "config$run: runDIAReferenceEnrichment (:890) [moves]", statements$run, "declared")
    group(seq, "G9a", "display tables (:891-900)", statements$display_tables, "contrast_results")
    group(seq, "G9b", "bus payload (:901-907)", statements$payload)
    group(seq, "G9c", "results' @args stamped (:908-911)", statements$stamp, "stamped")
    group(seq, "G9d", "UI parameters (:912-915)", statements$ui_parameters, "ui_parameters")
    group(seq, "G9e", "state commit (:916-917)", statements$commit, "saved")
    group(seq, "G9f", "tab_status, progress, the returned list (:918-925)", statements$status, "out")
    value <- figure(seq)
    portable <- portableWorkflowEnrichmentResults(seq$env$stamped)
    say("  [%s] committed: %s (saved = %s); portable results serialized %.2f MiB",
        label, case$manager$current_state, seq$env$saved$saved, mib(length(serialize(portable, NULL))))
    list(figure = value, rows = do.call(rbind, seq$rows))
}

managedSequence <- function(case, kernel, label) {
    seq <- newSequence(label, initialFrame(case, kernel))
    say("    -- the request half (controller)")
    group(seq, "R1", "context (:822-845)", statements$context, c("source_object", "dependencies", "contrasts_tbl"))
    group(seq, "R2", "annotation resolve (:853-854)", statements$resolve, "annotations")
    group(seq, "R3", "mapping and its warning [filter off]", managed$mapping, "mapping")
    group(seq, "R4", "progress 0.5, method (:877-880)", statements$method, "method")
    group(seq, "R5", "the kernel's name", managed$kernel_name, "kernel_name")
    group(seq, "R6", "transfer values: annotations as a saved display table", managed$transfer, "request")
    say("    -- the worker's compute half (emulated in this process; not counted)")
    group(seq, "W1", "source components, annotations materialized", managed$w_source, c("w_source", "w_annotations"),
        counted = FALSE)
    group(seq, "W2", "DA tables + refusal", managed$w_da, "w_da", counted = FALSE)
    group(seq, "W3", "annotation match", managed$w_match, "w_matching", counted = FALSE)
    group(seq, "W4", "kernel inputs, run, product, roots (computeDIAReferenceEnrichment)", managed$w_run,
        "w_computed", counted = FALSE)
    group(seq, "W5", "portable results saved as an rds state value", managed$w_save, c("w_portable", "delivered"),
        counted = FALSE)
    portable_serialized <- seq$env$portable_serialized
    say("    -- the moved values removed and collected")
    rm(list = grep("^w_", ls(seq$env), value = TRUE), envir = seq$env)
    invisible(gc(full = TRUE))
    say("    live heap after removal +%.1f MiB over start; delivered value %.3f MiB",
        mib(heapBytes(gc(), 1L) - seq$start), mib(as.numeric(lobstr::obj_size(seq$env$delivered))))
    say("    -- the accept half (controller)")
    group(seq, "A1", "results hydrated from the rds state value", managed$hydrate, "results_hydrated")
    group(seq, "A2", "provenance cleared; external-resource records", managed$records)
    group(seq, "A3", "bus filter metadata (:874-876)", managed$bus_filter, "organism_filter")
    group(seq, "A4", "display tables (:891-900)", statements$display_tables, "contrast_results")
    group(seq, "A5", "bus payload (:901-907)", statements$payload)
    group(seq, "A6", "results' @args stamped (:908-911)", statements$stamp, "stamped")
    group(seq, "A7", "UI parameters (:912-915)", statements$ui_parameters, "ui_parameters")
    group(seq, "A8", "state commit (:916-917)", statements$commit, "saved")
    group(seq, "A9", "tab_status, progress, the returned list (:918-925)", statements$status, "out")
    value <- figure(seq)
    say("  [%s] committed: %s (saved = %s); portable results serialized %.2f MiB", label, case$manager$current_state,
        seq$env$saved$saved, mib(portable_serialized))
    list(figure = value, rows = do.call(rbind, seq$rows), portable = portable_serialized)
}

# ---- the plan -----------------------------------------------------------------------------------------------------
kernels <- list(live = get("processEnrichments", envir = ns),
    double = get("runProtEnrichDeterministicProcessEnrichments", envir = ns))
cases <- if (identical(plan, "dry")) list(list(proteins = 400L, contrasts = 3L, kernels = "live")) else list(
    list(proteins = 16000L, contrasts = 3L, kernels = c("live", "double")),
    list(proteins = 4000L, contrasts = 3L, kernels = c("live", "double")),
    list(proteins = 16000L, contrasts = 1L, kernels = c("live", "double")),
    list(proteins = 4000L, contrasts = 1L, kernels = c("live", "double")))
samples <- 300L
all_rows <- list()
case_summary <- list()
elapsed <- function() proc.time()[["elapsed"]] - started
for (spec in cases) {
    if (elapsed() > ceiling_seconds - 60) {
        say("CEILING: %.0f s elapsed; case %d x %d x %d not started", elapsed(), spec$proteins, samples, spec$contrasts)
        break
    }
    say("")
    say("== case %d proteins x %d samples x %d contrasts (%.0f s elapsed) ==", spec$proteins, samples, spec$contrasts,
        elapsed())
    root <- tempfile("k113-enrichment-heap-")
    dir.create(root)
    case <- buildCase(spec$proteins, samples, spec$contrasts, root)
    for (kernel_name in spec$kernels) {
        if (elapsed() > ceiling_seconds - 60) {
            say("CEILING: %.0f s elapsed; kernel %s not started", elapsed(), kernel_name)
            break
        }
        tag <- sprintf("%d x %d x %d, %s", spec$proteins, samples, spec$contrasts, kernel_name)
        figures <- list()
        for (pass in c("pass 1 (cold)", "pass 2 (warm)")) {
            say("-- in process, %s, %s", tag, pass)
            ip <- inProcessSequence(case, kernels[[kernel_name]], sprintf("in process %s %s", tag, pass))
            figures[[paste("in_process", pass)]] <- ip$figure
            all_rows[[length(all_rows) + 1L]] <- cbind(case = tag, pass = pass, ip$rows)
        }
        for (pass in c("pass 1 (cold)", "pass 2 (warm)")) {
            say("-- managed controller, %s, %s", tag, pass)
            mc <- managedSequence(case, kernels[[kernel_name]], sprintf("managed %s %s", tag, pass))
            figures[[paste("managed", pass)]] <- mc$figure
            figures$portable <- mc$portable
            all_rows[[length(all_rows) + 1L]] <- cbind(case = tag, pass = pass, mc$rows)
        }
        difference <- c(cold = figures[["in_process pass 1 (cold)"]] - figures[["managed pass 1 (cold)"]],
            warm = figures[["in_process pass 2 (warm)"]] - figures[["managed pass 2 (warm)"]])
        case_summary[[tag]] <- c(in_process_cold = figures[["in_process pass 1 (cold)"]],
            in_process_warm = figures[["in_process pass 2 (warm)"]], managed_cold = figures[["managed pass 1 (cold)"]],
            managed_warm = figures[["managed pass 2 (warm)"]], difference_cold = difference[["cold"]],
            difference_warm = difference[["warm"]], portable = figures$portable)
        say("== %s: in process %.4f / %.4f GiB, managed controller %.4f / %.4f GiB, difference %.4f / %.4f GiB (cold / warm); portable %.2f MiB; network calls: %s",
            tag, figures[["in_process pass 1 (cold)"]] / 2^30, figures[["in_process pass 2 (warm)"]] / 2^30,
            figures[["managed pass 1 (cold)"]] / 2^30, figures[["managed pass 2 (warm)"]] / 2^30,
            difference[["cold"]] / 2^30, difference[["warm"]] / 2^30, mib(figures$portable), networkCalls())
    }
    closeDIAStateAdapter(case$manager)
    rm(case)
    unlink(root, recursive = TRUE)
    invisible(gc(full = TRUE))
}

say("")
say("== summary (GiB = 2^30 bytes) ==")
say("(%d cases)", length(case_summary))
for (tag in names(case_summary)) {
    s <- case_summary[[tag]]
    say("%-28s in process %.4f / %.4f, managed %.4f / %.4f, difference %.4f / %.4f (cold / warm), portable %.2f MiB",
        tag, s[["in_process_cold"]] / 2^30, s[["in_process_warm"]] / 2^30, s[["managed_cold"]] / 2^30,
        s[["managed_warm"]] / 2^30, s[["difference_cold"]] / 2^30, s[["difference_warm"]] / 2^30, mib(s[["portable"]]))
}
rows <- do.call(rbind, all_rows)
utils::write.table(rows, sub("[.]txt$", "-groups.tsv", report_path), sep = "\t", quote = FALSE, row.names = FALSE)
say("network calls over the whole probe: %s", networkCalls())
say("finished %s; %.0f s elapsed", format(Sys.time()), elapsed())
close(report)
